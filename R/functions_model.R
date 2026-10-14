# =============================================================================
#  Model functions: covariate construction, priors, formula, two-stage fit,
#  linear combinations and model comparison.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Covariates
# -----------------------------------------------------------------------------

# Index columns for the twelve spatially varying coefficients (SVCs).
# Each SVC needs its own copy of the area index in INLA.
prepare_spatial_indices <- function(datos) {
  for (r in 1:3) {
    datos[[paste0("dep_idx_int", r)]] <- datos$id_zona
    for (l in 1:3) {
      datos[[paste0("dep_idx_postint", r, "_", l)]] <- datos$id_zona
    }
  }
  datos
}


# Natural cubic spline bases B_l^(r)(t) = N_l{P_r(t)} for the three
# post-intervention trajectories.
#
# P_r(t) is the elapsed-time counter (1 in the implementation week), set to 0
# outside the active window [t_0^(r), t_0^(r+1)). The bases are computed on the
# full vector {P_r(t)}_{t=1}^T, so they vanish outside the window.
#
# Knots
#   r = 1: ns(df = 3) places the two internal knots at the 1/3 and 2/3
#          quantiles of the FULL vector, zeros included. With 116 zeros out of
#          209 weeks the first quantile is 0, which coincides with the lower
#          boundary knot, and ns() moves it inside. The resulting knots are
#          2.96 and 23.67 (elapsed weeks), with boundary knots 0 and 93.
#   r = 2: internal knots at 6 and 12 weeks.
#   r = 3: internal knots at 16 and 31 weeks.
#
# Input columns: semana, int2, int3 (step indicators 1(t >= t_0^(r))) and
# postint1, postint2, postint3 (elapsed-time counters without truncation).
prepare_spline_basis <- function(datos, df = 3) {

  datos <- datos %>%
    mutate(
      postint1_clean = ifelse(int2 == 1, 0, postint1),
      postint2_clean = ifelse(int3 == 1, 0, postint2),
      postint3_clean = postint3
    )

  df_semana <- datos %>%
    dplyr::select(semana, postint1_clean, postint2_clean, postint3_clean) %>%
    distinct(semana, .keep_all = TRUE) %>%
    arrange(semana)

  bases <- list(
    postint1 = ns(df_semana$postint1_clean, df = df),
    postint2 = ns(df_semana$postint2_clean, knots = c(6, 12)),
    postint3 = ns(df_semana$postint3_clean, knots = c(16, 31))
  )

  row_idx <- match(datos$semana, df_semana$semana)
  for (r in 1:3) {
    B <- bases[[paste0("postint", r)]][row_idx, , drop = FALSE]
    for (l in seq_len(ncol(B))) {
      datos[[paste0("bs_postint", r, "_", l)]] <- B[, l]
    }
  }

  datos
}


# -----------------------------------------------------------------------------
# 2. Hyperpriors
# -----------------------------------------------------------------------------
# PC prior on a precision kappa: P(kappa^{-1/2} > u) = alpha.
# PC prior on the AR(1) correlation (pc.cor0): P(|rho| > u) = alpha.
# Leroux mixing parameter: prior on logit(lambda).
define_hyperpriors <- function(
    leroux_prec_u   = 1,
    leroux_prec_a   = 0.6,
    prior_lambda    = "gaussian",
    lambda_params   = c(0, 2),
    ar1_prec_u      = 0.2,
    ar1_prec_a      = 0.01,
    ar1_rho_u       = 0.7,
    ar1_rho_a       = 0.05,
    svc_prec_u      = 0.1,
    svc_prec_a      = 0.01,
    family          = "poisson",
    nbinomial_u     = 1,
    nbinomial_a     = 0.6
) {
  list(
    space_leroux = list(
      prec   = list(prior = "pc.prec", param = c(leroux_prec_u, leroux_prec_a)),
      lambda = list(prior = prior_lambda, param = lambda_params)
    ),
    temporal_ar1 = list(
      prec = list(prior = "pc.prec", param = c(ar1_prec_u, ar1_prec_a)),
      rho  = list(prior = "pc.cor0", param = c(ar1_rho_u, ar1_rho_a))
    ),
    svc = list(
      prec = list(prior = "pc.prec", param = c(svc_prec_u, svc_prec_a))
    ),
    family = switch(
      family,
      poisson   = list(),
      nbinomial = list(theta = list(prior = "pc.prec",
                                    param = c(nbinomial_u, nbinomial_a))),
      stop("family must be 'poisson' or 'nbinomial'")
    ),
    family_name = family
  )
}


# -----------------------------------------------------------------------------
# 3. Formulas
# -----------------------------------------------------------------------------
FIXED_TERMS <- c(
  "int1", "int2", "int3",
  paste0("bs_postint", rep(1:3, each = 3), "_", rep(1:3, 3)),
  "mayor65_porc", "homb_porc", "tipo_suelo"
)

# Fixed effects only (used for the overdispersion comparison).
build_formula_fixed <- function() {
  as.formula(paste("suma_zona ~", paste(FIXED_TERMS, collapse = " + ")))
}

# Full model of Equation (predictor) in the article:
#   Leroux spatial effect u_j, AR(1) temporal effect gamma_t, and twelve
#   Besag SVCs (three level shifts and nine spline coefficients).
# `graph` is the adjacency matrix of the areas. The SVC weights (int_r and
# bs_postint_r_l) are read from the data.
build_formula_full <- function(hypers, graph) {
  suma_zona ~ int1 + int2 + int3 +
    bs_postint1_1 + bs_postint1_2 + bs_postint1_3 +
    bs_postint2_1 + bs_postint2_2 + bs_postint2_3 +
    bs_postint3_1 + bs_postint3_2 + bs_postint3_3 +
    mayor65_porc + homb_porc + tipo_suelo +
    # Baseline spatial effect (Leroux)
    f(id_zona, model = "besagproper2", graph = graph,
      hyper = hypers$space_leroux) +
    # Temporal effect (AR(1))
    f(semana, model = "ar1", hyper = hypers$temporal_ar1) +
    # SVCs of the immediate level shifts
    f(dep_idx_int1, int1, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    f(dep_idx_int2, int2, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    f(dep_idx_int3, int3, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    # SVCs of the spline coefficients
    f(dep_idx_postint1_1, bs_postint1_1, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    f(dep_idx_postint1_2, bs_postint1_2, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    f(dep_idx_postint1_3, bs_postint1_3, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    f(dep_idx_postint2_1, bs_postint2_1, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    f(dep_idx_postint2_2, bs_postint2_2, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    f(dep_idx_postint2_3, bs_postint2_3, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    f(dep_idx_postint3_1, bs_postint3_1, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    f(dep_idx_postint3_2, bs_postint3_2, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc) +
    f(dep_idx_postint3_3, bs_postint3_3, model = "besag", graph = graph,
      scale.model = TRUE, constr = TRUE, hyper = hypers$svc)
}


# -----------------------------------------------------------------------------
# 4. INLA controls
# -----------------------------------------------------------------------------
# int.strategy = "eb": hyperparameters fixed at their posterior mode.
# cpo = TRUE is needed for the PIT/CPO diagnostics (scripts/04_diagnostics.R).
define_inla_controls <- function(strategy = "simplified.laplace", cpo = TRUE) {
  list(
    compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = cpo,
                   return.marginals.predictor = TRUE),
    predictor = list(compute = TRUE, link = 1),
    inla = list(
      strategy               = strategy,
      int.strategy           = "eb",
      stupid.search.max.iter = 1500,
      h                      = 0.01,
      restart                = 0
    )
  )
}


# -----------------------------------------------------------------------------
# 5. Two-stage fit (Krainski et al., 2019)
# -----------------------------------------------------------------------------
# Stage 1: experimental (compact) mode to locate the hyperparameter modes.
# Stage 2: classic mode started at those modes.
fit_inla_twostage <- function(formula, datos, family = "poisson",
                              offset_var = NULL, hypers = NULL,
                              controls = NULL, lincomb = NULL,
                              verbose = FALSE) {

  if (is.null(hypers))   hypers   <- define_hyperpriors(family = family)
  if (is.null(controls)) controls <- define_inla_controls()
  offset <- if (!is.null(offset_var)) log(datos[[offset_var]]) else NULL

  message(sprintf("Two-stage INLA fit | family = %s | strategy = %s",
                  family, controls$inla$strategy))

  # Stage 1
  time_stage1 <- system.time({
    fit_stage1 <- INLA::inla(
      formula, family = family, data = datos, offset = offset,
      control.compute   = controls$compute,
      control.predictor = controls$predictor,
      control.inla      = controls$inla,
      control.family    = list(hyper = hypers$family),
      inla.mode         = "experimental",
      verbose           = verbose
    )
  })
  theta_mode <- fit_stage1$mode$theta
  message(sprintf("  Stage 1: %.1f s", time_stage1[3]))

  # Stage 2
  control_inla_2 <- controls$inla
  control_inla_2$mode.known <- TRUE

  time_stage2 <- system.time({
    fit_stage2 <- INLA::inla(
      formula, family = family, data = datos, offset = offset,
      control.compute   = controls$compute,
      control.predictor = controls$predictor,
      control.inla      = control_inla_2,
      control.family    = list(hyper = hypers$family),
      control.mode      = list(theta = theta_mode, restart = FALSE),
      inla.mode         = "classic",
      lincomb           = lincomb,
      verbose           = verbose
    )
  })
  message(sprintf("  Stage 2: %.1f s", time_stage2[3]))

  list(fit         = fit_stage2,
       fit_stage1  = fit_stage1,
       theta_mode  = theta_mode,
       timing      = list(stage1 = time_stage1, stage2 = time_stage2),
       family_used = family)
}


# -----------------------------------------------------------------------------
# 6. Linear combinations: psi_j^(r) = beta^(int,r) + b_j^(int,r)
# -----------------------------------------------------------------------------
create_intervention_lincomb <- function(n_areas = N_AREAS) {
  I_n <- diag(n_areas)
  lc  <- lapply(1:3, function(r) {
    out <- inla.make.lincombs(
      setNames(list(rep(1, n_areas), I_n),
               c(paste0("int", r), paste0("dep_idx_int", r))))
    names(out) <- paste0("psi_int", r, "_area_", seq_len(n_areas))
    out
  })
  do.call(c, lc)
}


# -----------------------------------------------------------------------------
# 7. Model comparison (Supplementary Material)
# -----------------------------------------------------------------------------
compare_models <- function(models, labels = names(models)) {
  if (is.null(labels)) labels <- paste0("Model_", seq_along(models))

  get_hyper_mean <- function(fit, pattern) {
    idx <- grep(pattern, rownames(fit$summary.hyperpar), fixed = TRUE)
    if (length(idx) == 0) NA_real_ else fit$summary.hyperpar[idx[1], "mean"]
  }

  out <- do.call(rbind, lapply(seq_along(models), function(i) {
    m   <- models[[i]]
    fit <- m$fit
    data.frame(
      Model     = labels[i],
      Family    = m$family_used,
      DIC       = round(fit$dic$dic, 2),
      WAIC      = round(fit$waic$waic, 2),
      Intercept = round(fit$summary.fixed["(Intercept)", "mean"], 4),
      Lambda    = round(get_hyper_mean(fit, "Lambda for id_zona"), 4),
      Rho       = round(get_hyper_mean(fit, "Rho for semana"), 4),
      Time_s    = round(sum(sapply(m$timing, function(x) x[3])), 2),
      stringsAsFactors = FALSE
    )
  }))
  out$DIC_diff <- round(out$DIC - min(out$DIC, na.rm = TRUE), 2)
  out
}
