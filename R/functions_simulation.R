# =============================================================================
#  Simulation study: data-generating mechanism, fit and performance measures
#  ADEMP framework of Morris, White and Crowther (2019), Stat Med 38:2074-2102.
#
#  Requires, in the calling environment (defined in scripts/03_simulation.R):
#    beta_true, beta_spline1..3, kappa_svc, t_vec, T0, T_WEEKS, I1..I3,
#    B1..B3, J_sim, E_sim, age65_sim, homb_porc_sim, tipo_suelo_sim,
#    W_sim_sparse, Q_leroux_sim, adj_sim, V_icar, lam_icar, hyper_sim,
#    control_inla_sim, N_POST_SAMPLE
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a


# -----------------------------------------------------------------------------
# Step indicators 1(t >= t_0^(r))
# -----------------------------------------------------------------------------
make_indicators <- function(t_vec, t0) {
  list(I1 = as.integer(t_vec >= t0[1]),
       I2 = as.integer(t_vec >= t0[2]),
       I3 = as.integer(t_vec >= t0[3]))
}


# -----------------------------------------------------------------------------
# Connected subgraph of n areas, grown from a random seed area
# -----------------------------------------------------------------------------
select_connected_areas <- function(adj, n, seed = 42) {
  valid <- which(rowSums(adj) > 0)
  set.seed(seed)
  selected <- sample(valid, 1)
  frontier <- setdiff(valid[adj[selected, valid] != 0], selected)
  while (length(selected) < n && length(frontier) > 0) {
    new      <- sample(frontier, 1)
    selected <- c(selected, new)
    frontier <- setdiff(union(frontier, valid[adj[new, valid] != 0]), selected)
  }
  if (length(selected) < n)
    stop(sprintf("Only %d connected areas could be selected.", length(selected)))
  sort(selected)
}


# -----------------------------------------------------------------------------
# Eigen-decomposition of the ICAR structure, scaled so that the geometric mean
# of the marginal variances equals 1 (as with scale.model = TRUE in INLA).
# -----------------------------------------------------------------------------
scaled_icar_eigen <- function(W) {
  R   <- as.matrix(Diagonal(x = rowSums(W)) - W)
  eig <- eigen(R, symmetric = TRUE)
  keep <- eig$values > 1e-10 * max(eig$values)
  if (sum(!keep) != 1L)
    warning(sprintf("ICAR null space has dimension %d (expected 1).", sum(!keep)))
  V     <- eig$vectors[, keep, drop = FALSE]
  lam   <- eig$values[keep]
  s     <- exp(mean(log(diag(V %*% (t(V) / lam)))))
  lam_s <- lam * s
  check <- exp(mean(log(diag(V %*% (t(V) / lam_s)))))
  stopifnot(abs(check - 1) < 1e-6)
  list(V = V, lam = lam_s, scale = s)
}


# -----------------------------------------------------------------------------
# Draws from the latent components
# -----------------------------------------------------------------------------
draw_icar <- function(kappa) {
  z <- rnorm(length(lam_icar))
  x <- as.numeric(V_icar %*% (z / sqrt(kappa * lam_icar)))
  x - mean(x)                                        # sum-to-zero
}

draw_leroux <- function(W, Q_leroux, kappa_phi, lambda) {
  Q <- kappa_phi * (lambda * Q_leroux + (1 - lambda) * Diagonal(nrow(W)))
  as.numeric(solve(chol(Q), rnorm(nrow(W))))
}

draw_ar1 <- function(TT, rho, sigma_eps) {
  eps    <- numeric(TT)
  eps[1] <- rnorm(1, 0, sigma_eps / sqrt(1 - rho^2))
  for (t in seq(2, TT)) eps[t] <- rho * eps[t - 1] + rnorm(1, 0, sigma_eps)
  eps
}


# -----------------------------------------------------------------------------
# Data generation. `effect` scales the intervention effects and SVCs
# (NULL = observed effect sizes).
# -----------------------------------------------------------------------------
generate_dataset <- function(params, effect = NULL) {

  beta_factor   <- effect$beta_factor   %||% 1
  spline_factor <- effect$spline_factor %||% 1
  svc_factor    <- effect$svc_factor    %||% 1

  beta_sim <- beta_true
  beta_sim[c("int1", "int2", "int3")] <- beta_true[c("int1", "int2", "int3")] *
    beta_factor
  eta_splines <- as.numeric(B1 %*% (beta_spline1 * spline_factor) +
                            B2 %*% (beta_spline2 * spline_factor) +
                            B3 %*% (beta_spline3 * spline_factor))

  phi      <- draw_leroux(W_sim_sparse, Q_leroux_sim,
                          params$kappa_phi, params$lambda)
  phi_int1 <- draw_icar(kappa_svc$int1) * svc_factor
  phi_int2 <- draw_icar(kappa_svc$int2) * svc_factor
  phi_int3 <- draw_icar(kappa_svc$int3) * svc_factor
  gamma1   <- sapply(1:3, function(l) draw_icar(kappa_svc$post1[l]))
  gamma2   <- sapply(1:3, function(l) draw_icar(kappa_svc$post2[l]))
  gamma3   <- sapply(1:3, function(l) draw_icar(kappa_svc$post3[l]))
  f1 <- B1 %*% t(gamma1) * svc_factor
  f2 <- B2 %*% t(gamma2) * svc_factor
  f3 <- B3 %*% t(gamma3) * svc_factor

  eps <- draw_ar1(T_WEEKS, params$rho, params$sigma_eps)   # common to all areas

  Y <- matrix(NA_integer_, nrow = J_sim, ncol = T_WEEKS)
  for (j in seq_len(J_sim)) {
    eta_fixed <- beta_sim["intercept"] +
      beta_sim["int1"] * I1 + beta_sim["int2"] * I2 + beta_sim["int3"] * I3 +
      beta_sim["age65"]      * age65_sim[j] +
      beta_sim["homb_porc"]  * homb_porc_sim[j] +
      beta_sim["tipo_suelo"] * tipo_suelo_sim[j]
    eta_svc <- phi_int1[j] * I1 + phi_int2[j] * I2 + phi_int3[j] * I3 +
      f1[, j] + f2[, j] + f3[, j]
    eta    <- eta_fixed + eta_splines + eta_svc + phi[j] + eps
    Y[j, ] <- rpois(T_WEEKS, E_sim[j, ] * exp(eta))
  }

  list(Y = Y, phi = phi,
       phi_int1 = phi_int1, phi_int2 = phi_int2, phi_int3 = phi_int3,
       gamma1 = gamma1, gamma2 = gamma2, gamma3 = gamma3,
       f1 = f1, f2 = f2, f3 = f3, eps = eps,
       params = params, beta_sim = beta_sim)
}


# -----------------------------------------------------------------------------
# Two-stage INLA fit of one synthetic data set
# -----------------------------------------------------------------------------
fit_inla_sim <- function(dat) {
  TT  <- ncol(dat$Y)
  ind <- make_indicators(seq_len(TT), T0)
  df_long <- expand.grid(j = seq_len(J_sim), t = seq_len(TT)) |>
    dplyr::mutate(
      y          = as.vector(dat$Y),
      E_offset   = as.vector(E_sim),
      int1       = ind$I1[t], int2 = ind$I2[t], int3 = ind$I3[t],
      bs_post1_1 = B1[t, 1], bs_post1_2 = B1[t, 2], bs_post1_3 = B1[t, 3],
      bs_post2_1 = B2[t, 1], bs_post2_2 = B2[t, 2], bs_post2_3 = B2[t, 3],
      bs_post3_1 = B3[t, 1], bs_post3_2 = B3[t, 2], bs_post3_3 = B3[t, 3],
      age65      = age65_sim[j],
      homb_porc  = homb_porc_sim[j],
      tipo_suelo = tipo_suelo_sim[j],
      id_zona    = j, semana = t,
      dep_int1    = j, dep_int2 = j, dep_int3 = j,
      dep_post1_1 = j, dep_post1_2 = j, dep_post1_3 = j,
      dep_post2_1 = j, dep_post2_2 = j, dep_post2_3 = j,
      dep_post3_1 = j, dep_post3_2 = j, dep_post3_3 = j
    )

  h <- hyper_sim
  formula_sim <- y ~
    int1 + int2 + int3 +
    bs_post1_1 + bs_post1_2 + bs_post1_3 +
    bs_post2_1 + bs_post2_2 + bs_post2_3 +
    bs_post3_1 + bs_post3_2 + bs_post3_3 +
    age65 + homb_porc + tipo_suelo +
    f(id_zona, model = "besagproper2", graph = adj_sim, hyper = h$leroux) +
    f(semana, model = "ar1", hyper = h$ar1) +
    f(dep_int1, int1, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_int2, int2, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_int3, int3, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_post1_1, bs_post1_1, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_post1_2, bs_post1_2, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_post1_3, bs_post1_3, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_post2_1, bs_post2_1, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_post2_2, bs_post2_2, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_post2_3, bs_post2_3, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_post3_1, bs_post3_1, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_post3_2, bs_post3_2, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc) +
    f(dep_post3_3, bs_post3_3, model = "besag", graph = adj_sim,
      scale.model = TRUE, constr = TRUE, hyper = h$svc)

  # Stage 1
  err1 <- NA_character_
  fit1 <- tryCatch(
    inla(formula_sim, family = "poisson", data = df_long, E = E_offset,
         control.compute   = list(dic = FALSE, waic = FALSE, cpo = FALSE,
                                  return.marginals = FALSE, config = FALSE),
         control.predictor = list(compute = FALSE, link = 1),
         control.inla      = control_inla_sim,
         inla.mode         = "experimental", verbose = FALSE),
    error = function(e) { err1 <<- conditionMessage(e); NULL }
  )
  if (is.null(fit1)) return(list(fit = NULL, error = err1, stage = 1L))
  theta_ref <- fit1$mode$theta
  rm(fit1); gc()

  # Stage 2
  err2 <- NA_character_
  fit2 <- tryCatch(
    inla(formula_sim, family = "poisson", data = df_long, E = E_offset,
         control.compute   = list(dic = TRUE, waic = TRUE, cpo = FALSE,
                                  return.marginals = FALSE, config = TRUE),
         control.predictor = list(compute = FALSE, link = 1),
         control.inla      = control_inla_sim,
         control.mode      = list(theta = theta_ref, restart = FALSE),
         inla.mode         = "classic", verbose = FALSE),
    error = function(e) { err2 <<- conditionMessage(e); NULL }
  )
  list(fit = fit2, error = err2, stage = if (is.null(fit2)) 2L else NA_integer_)
}


# -----------------------------------------------------------------------------
# Quantities kept from each replicate
# -----------------------------------------------------------------------------
SVC_NAMES <- c("dep_int1", "dep_int2", "dep_int3",
               paste0("dep_post", rep(1:3, each = 3), "_", rep(1:3, 3)))

extract_essentials <- function(fit, dat) {
  if (is.null(fit)) return(NULL)
  vcov_splines <- tryCatch({
    sel <- setNames(as.list(rep(1, 9)),
                    paste0("bs_post", rep(1:3, each = 3), "_", rep(1:3, 3)))
    samps    <- inla.posterior.sample(N_POST_SAMPLE, fit, selection = sel)
    beta_mat <- sapply(samps, function(s) s$latent[, 1])
    lapply(1:3, function(k) cov(t(beta_mat[(3 * k - 2):(3 * k), , drop = FALSE])))
  }, error = function(e) NULL)

  cols <- c("mean", "sd", "0.025quant", "0.975quant")
  list(
    fixed        = fit$summary.fixed[, cols],
    hyper        = fit$summary.hyperpar[, cols],
    vcov_splines = vcov_splines,
    svc = setNames(lapply(SVC_NAMES, function(nm) {
      if (nm %in% names(fit$summary.random)) {
        fit$summary.random[[nm]][seq_len(J_sim), cols]
      } else {
        matrix(NA_real_, J_sim, length(cols), dimnames = list(NULL, cols))
      }
    }), SVC_NAMES),
    true_phi_int1 = dat$phi_int1, true_phi_int2 = dat$phi_int2,
    true_phi_int3 = dat$phi_int3,
    true_f1 = dat$f1, true_f2 = dat$f2, true_f3 = dat$f3,
    true_gamma1 = dat$gamma1, true_gamma2 = dat$gamma2,
    true_gamma3 = dat$gamma3,
    beta_sim = dat$beta_sim
  )
}


# -----------------------------------------------------------------------------
# Performance measures (Morris et al., 2019, Table 6)
# -----------------------------------------------------------------------------
compute_scalar_metrics <- function(hat, lo, hi, sd_post, true_val) {
  ok <- !is.na(hat)
  n  <- sum(ok)
  na_out <- c(Bias = NA, Bias_MCSE = NA, RBias = NA, RBias_MCSE = NA,
              EmpSE = NA, EmpSE_MCSE = NA, ModSE = NA, RelErrModSE = NA,
              Cov95 = NA, Cov95_MCSE = NA, RMSE = NA, MSE_MCSE = NA,
              n_used = n)
  if (n < 2L) return(na_out)
  hat <- hat[ok]; lo <- lo[ok]; hi <- hi[ok]
  true_val <- as.numeric(true_val)

  bias       <- mean(hat) - true_val
  empse      <- stats::sd(hat)
  bias_mcse  <- empse / sqrt(n)
  empse_mcse <- empse / sqrt(2 * (n - 1))
  if (abs(true_val) > .Machine$double.eps^0.5) {
    rbias      <- 100 * bias / true_val
    rbias_mcse <- 100 * bias_mcse / abs(true_val)
  } else {
    rbias <- NA_real_; rbias_mcse <- NA_real_
  }
  sp    <- sd_post[ok]
  modse <- if (all(is.na(sp))) NA_real_ else sqrt(mean(sp^2, na.rm = TRUE))
  rel_err_modse <- if (is.na(modse) || empse <= 0) NA_real_ else
    100 * (modse / empse - 1)
  coverage <- 100 * mean(lo <= true_val & true_val <= hi)
  cov_mcse <- sqrt(coverage * (100 - coverage) / n)
  sq       <- (hat - true_val)^2

  c(Bias = bias, Bias_MCSE = bias_mcse, RBias = rbias, RBias_MCSE = rbias_mcse,
    EmpSE = empse, EmpSE_MCSE = empse_mcse, ModSE = modse,
    RelErrModSE = rel_err_modse, Cov95 = coverage, Cov95_MCSE = cov_mcse,
    RMSE = sqrt(mean(sq)), MSE_MCSE = stats::sd(sq) / sqrt(n), n_used = n)
}


# Area-level SVC recovery: bias, MSE, coverage and model SE averaged over
# areas within each replicate.
svc_summary_by_rep <- function(res, nm, true_getter) {
  t(vapply(seq_along(res), function(i) {
    s   <- res[[i]]$svc[[nm]]
    hat <- s[, "mean"]; lo <- s[, "0.025quant"]; hi <- s[, "0.975quant"]
    sdp <- s[, "sd"]
    tru <- true_getter(res[[i]])
    v   <- !is.na(hat) & !is.na(tru)
    if (!any(v)) return(c(bias = NA, mse = NA, cov = NA, modse2 = NA))
    c(bias   = mean(hat[v] - tru[v]),
      mse    = mean((hat[v] - tru[v])^2),
      cov    = 100 * mean(lo[v] <= tru[v] & tru[v] <= hi[v]),
      modse2 = mean(sdp[v]^2, na.rm = TRUE))
  }, FUN.VALUE = numeric(4)))
}

svc_metrics_from_reps <- function(m) {
  m <- m[stats::complete.cases(m), , drop = FALSE]
  n <- nrow(m)
  if (n < 2L)
    return(c(Bias = NA, Bias_MCSE = NA, RMSE = NA, MSE_MCSE = NA,
             Cov95 = NA, Cov95_MCSE = NA, ModSE = NA, ModSE_ratio = NA,
             n_rep = n))
  mse <- mean(m[, "mse"])
  c(Bias        = mean(m[, "bias"]),
    Bias_MCSE   = stats::sd(m[, "bias"]) / sqrt(n),
    RMSE        = sqrt(mse),
    MSE_MCSE    = stats::sd(m[, "mse"]) / sqrt(n),
    Cov95       = mean(m[, "cov"]),
    Cov95_MCSE  = stats::sd(m[, "cov"]) / sqrt(n),
    ModSE       = sqrt(mean(m[, "modse2"])),
    ModSE_ratio = sqrt(mean(m[, "modse2"])) / sqrt(mse),
    n_rep       = n)
}

validate_svcs_level <- function(res) {
  out <- t(sapply(1:3, function(r) {
    tf <- paste0("true_phi_int", r)
    svc_metrics_from_reps(
      svc_summary_by_rep(res, paste0("dep_int", r), function(x) x[[tf]]))
  }))
  rownames(out) <- paste0("Int_", 1:3)
  out
}

validate_svcs_splines <- function(res, k) {
  out <- t(sapply(1:3, function(l) {
    svc_metrics_from_reps(
      svc_summary_by_rep(res, paste0("dep_post", k, "_", l),
                         function(x) x[[paste0("true_gamma", k)]][, l]))
  }))
  rownames(out) <- paste0("Basis_", 1:3)
  out
}


# Mean integrated squared error of the population-average trajectory f^(r)(t)
# over its active window (x 1000), and its normalized version.
compute_mise <- function(res, k, beta_spl) {
  Bk     <- list(B1, B2, B3)[[k]]
  t_end  <- c(T0[2] - 1, T0[3] - 1, T_WEEKS)[k]
  t_post <- which(t_vec >= T0[k] & t_vec <= t_end)
  nm     <- paste0("bs_post", k, "_", 1:3)
  f_true <- as.numeric(Bk %*% beta_spl)

  ise <- vapply(res, function(r) {
    f_hat <- as.numeric(Bk %*% r$fixed[nm, "mean"])
    sum((f_hat[t_post] - f_true[t_post])^2)
  }, numeric(1))

  n     <- sum(!is.na(ise))
  mise  <- mean(ise, na.rm = TRUE) * 1e3
  mcse  <- stats::sd(ise, na.rm = TRUE) / sqrt(n) * 1e3
  denom <- sum(f_true[t_post]^2)
  if (denom < 1e-10) {
    mise_n <- NA_real_; mcse_n <- NA_real_
  } else {
    mise_n <- mise / (denom * 1e3); mcse_n <- mcse / (denom * 1e3)
  }
  c(MISE = mise, MISE_MCSE = mcse, MISE_norm = mise_n,
    MISE_norm_MCSE = mcse_n, n_rep = n)
}


# All performance measures for one scenario.
extract_metrics <- function(results, true_vals) {
  ok  <- !vapply(results, is.null, logical(1))
  res <- results[ok]
  message(sprintf("%d / %d replicates converged.", sum(ok), length(results)))
  if (length(res) == 0) return(NULL)

  # Fixed effects
  beta_names <- names(true_vals$beta)
  get_fixed <- function(col) {
    sapply(beta_names, function(nm) {
      inla_nm <- if (nm == "intercept") "(Intercept)" else nm
      vapply(res, function(r) r$fixed[inla_nm, col], numeric(1))
    })
  }
  hat <- get_fixed("mean"); lo <- get_fixed("0.025quant")
  hi  <- get_fixed("0.975quant"); sdp <- get_fixed("sd")
  beta_metrics <- t(sapply(beta_names, function(nm)
    compute_scalar_metrics(hat[, nm], lo[, nm], hi[, nm], sdp[, nm],
                           true_vals$beta[nm])))

  # Hyperparameters
  get_hyper <- function(r, pattern) {
    h   <- r$hyper
    idx <- grep(pattern, rownames(h), fixed = TRUE)
    if (length(idx) == 0) return(c(mean = NA, sd = NA, lo = NA, hi = NA))
    c(mean = h[idx[1], "mean"], sd = h[idx[1], "sd"],
      lo = h[idx[1], "0.025quant"], hi = h[idx[1], "0.975quant"])
  }
  svc_prec <- c(int1 = kappa_svc$int1, int2 = kappa_svc$int2,
                int3 = kappa_svc$int3,
                setNames(kappa_svc$post1, paste0("post1_", 1:3)),
                setNames(kappa_svc$post2, paste0("post2_", 1:3)),
                setNames(kappa_svc$post3, paste0("post3_", 1:3)))
  hyper_params <- c(
    list(kappa_phi = list(true = true_vals$kappa_phi,
                          pattern = "Precision for id_zona"),
         lambda    = list(true = true_vals$lambda,
                          pattern = "Lambda for id_zona"),
         rho       = list(true = true_vals$rho,
                          pattern = "Rho for semana"),
         prec_ar1  = list(true = (1 - true_vals$rho^2) / true_vals$sigma_eps^2,
                          pattern = "Precision for semana")),
    lapply(setNames(names(svc_prec), paste0("prec_", names(svc_prec))),
           function(nm) list(true = svc_prec[[nm]],
                             pattern = paste0("Precision for dep_", nm)))
  )
  hyper_metrics <- do.call(rbind, lapply(hyper_params, function(p) {
    v <- t(sapply(res, get_hyper, pattern = p$pattern))
    compute_scalar_metrics(v[, "mean"], v[, "lo"], v[, "hi"], v[, "sd"], p$true)
  }))

  # Trajectories and SVCs
  mise_metrics <- t(sapply(1:3, function(k)
    compute_mise(res, k, true_vals$beta[paste0("bs_post", k, "_", 1:3)])))
  rownames(mise_metrics) <- paste0("Int_", 1:3)

  list(beta_metrics      = beta_metrics,
       hyper_metrics     = hyper_metrics,
       mise_metrics      = mise_metrics,
       svc_level_metrics = validate_svcs_level(res),
       svc_metrics       = lapply(1:3, function(k) validate_svcs_splines(res, k)),
       n_converged       = length(res),
       n_replicates      = length(results))
}


# -----------------------------------------------------------------------------
# Execution: one file per replicate, so an interrupted run can be resumed.
# -----------------------------------------------------------------------------
run_scenario <- function(label, S, seed_base, gen_fun, output_dir,
                         n_workers = 1L, n_cores = 1L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  seeds   <- seed_base + seq_len(S)
  use_par <- n_workers > 1L && .Platform$OS.type != "windows"
  threads <- max(1L, floor(n_cores / max(1L, n_workers)))

  run_one <- function(i) {
    out_file <- file.path(output_dir, sprintf("rep_%s_%03d.rds", label, i))
    if (file.exists(out_file)) return(NULL)
    set.seed(seeds[i])
    if (use_par) INLA::inla.setOption(num.threads = threads)
    t0  <- Sys.time()
    dat <- gen_fun()
    fr  <- fit_inla_sim(dat)
    ess <- extract_essentials(fr$fit, dat)
    saveRDS(ess, out_file)
    data.frame(replicate = i, seed = seeds[i], converged = !is.null(ess),
               error = if (is.null(ess)) fr$error %||% "extract failed" else NA,
               minutes = round(as.numeric(difftime(Sys.time(), t0,
                                                   units = "mins")), 2))
  }

  message(sprintf("=== Scenario %s: %d replicates ===", label, S))
  logs <- if (use_par) {
    parallel::mclapply(seq_len(S), run_one, mc.cores = n_workers,
                       mc.preschedule = FALSE)
  } else {
    lapply(seq_len(S), run_one)
  }
  logs <- do.call(rbind, Filter(is.data.frame, logs))
  if (!is.null(logs))
    utils::write.csv(logs, file.path(output_dir, sprintf("log_%s.csv", label)),
                     row.names = FALSE)
  invisible(logs)
}

load_scenario <- function(label, output_dir) {
  files <- sort(list.files(output_dir, pattern = sprintf("^rep_%s_\\d+\\.rds$", label),
                           full.names = TRUE))
  lapply(files, readRDS)
}
