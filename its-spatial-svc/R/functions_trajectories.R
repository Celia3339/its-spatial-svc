# =============================================================================
#  Post-intervention trajectories from joint posterior samples
#
#    f^(r)(t)       = B^(r)(t)' beta^(postint,r)                (population)
#    g_j^(r)(t)     = B^(r)(t)' b_j^(postint,r)                 (area deviation)
#    zeta_j^(r)(t)  = f^(r)(t) + g_j^(r)(t)                     (area-specific)
# =============================================================================


# Spline bases on the active window of intervention r (one row per week).
get_basis_grid <- function(datos, r, t0 = T0, t_max = T_WEEKS) {
  t_end <- c(t0[-1] - 1L, t_max)[r]
  cols  <- paste0("bs_postint", r, "_", 1:3)
  grid  <- datos %>%
    filter(semana >= t0[r], semana <= t_end) %>%
    distinct(semana, .keep_all = TRUE) %>%
    arrange(semana)
  list(week = grid$semana, B = as.matrix(grid[, cols]))
}


# Indices of a latent component in the sampled latent vector.
latent_index <- function(latent_names, prefix) {
  grep(paste0("^", prefix, ":"), latent_names)
}


# Posterior samples of the fixed spline coefficients (n_samp x 3).
get_beta_samples <- function(samples, r) {
  latent_names <- rownames(samples[[1]]$latent)
  idx <- sapply(paste0("bs_postint", r, "_", 1:3),
                function(p) latent_index(latent_names, p))
  t(sapply(samples, function(s) s$latent[idx, 1]))
}


# Posterior samples of the area-level spline SVCs: list over samples of
# 3 x J matrices.
get_svc_samples <- function(samples, r) {
  latent_names <- rownames(samples[[1]]$latent)
  idx <- lapply(1:3, function(l)
    latent_index(latent_names, paste0("dep_idx_postint", r, "_", l)))
  lapply(samples, function(s) {
    lat <- s$latent[, 1]
    rbind(lat[idx[[1]]], lat[idx[[2]]], lat[idx[[3]]])
  })
}


# Population-average trajectory f^(r)(t): mean, 95% band, P(f < 0).
summarise_population_trajectory <- function(grid, beta_samples, r) {
  eff <- grid$B %*% t(beta_samples)          # n_t x n_samp
  tibble(
    intervention = r,
    week         = grid$week,
    mean         = rowMeans(eff),
    q025         = apply(eff, 1, quantile, probs = 0.025),
    q975         = apply(eff, 1, quantile, probs = 0.975),
    prob_neg     = rowMeans(eff < 0)
  )
}


# Area-level trajectories: posterior mean and P(. < 0) for each area and week.
# coef_samples: list over samples of 3 x J coefficient matrices.
summarise_area_trajectories <- function(grid, coef_samples, r, value = "zeta") {
  n_samp    <- length(coef_samples)
  sum_eta   <- 0
  count_neg <- 0
  for (s in seq_len(n_samp)) {
    eta_s     <- grid$B %*% coef_samples[[s]]    # n_t x J
    sum_eta   <- sum_eta + eta_s
    count_neg <- count_neg + (eta_s < 0)
  }
  J <- ncol(coef_samples[[1]])
  tibble(
    intervention = r,
    week         = rep(grid$week, times = J),
    id_zona      = rep(seq_len(J), each = length(grid$week)),
    quantity     = value,
    mean         = as.vector(sum_eta / n_samp),
    prob_neg     = as.vector(count_neg / n_samp)
  )
}
