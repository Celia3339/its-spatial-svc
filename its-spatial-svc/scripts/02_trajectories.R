# =============================================================================
#  02. Post-intervention trajectories
#
#  From joint posterior samples of the selected model (fit_1) computes
#    - f^(r)(t): population-average trajectory, 95% band and P(f < 0)
#      (Figure: reconstructed post-intervention trajectories)
#    - zeta_j^(r)(t) = f^(r)(t) + g_j^(r)(t) and g_j^(r)(t) for every area
#      (Figures: area-specific trajectories and trajectory maps)
#
#  Output: output/trajectories_population.rds
#          output/trajectories_area.rds
# =============================================================================

source("R/00_setup.R")
source("R/functions_model.R")
source("R/functions_trajectories.R")
set_threads()

N_SAMPLES <- 2000L
SEED      <- 123L

# -----------------------------------------------------------------------------
# Data and fitted model
# -----------------------------------------------------------------------------
datos_fit <- readRDS(FILE_DATA)
datos_fit <- prepare_spline_basis(datos_fit, df = 3)
modelo    <- readRDS(file.path(PATH_OUTPUT, "fit_1_poisson.rds"))$fit

grids <- lapply(1:3, function(r) get_basis_grid(datos_fit, r))

# -----------------------------------------------------------------------------
# Joint posterior samples (requires control.compute = list(config = TRUE)).
# For bit-for-bit reproducibility also set num.threads = "1:1".
# -----------------------------------------------------------------------------
set.seed(SEED)
samples <- inla.posterior.sample(n = N_SAMPLES, result = modelo, seed = SEED)

# -----------------------------------------------------------------------------
# Population-average trajectories f^(r)(t)
# -----------------------------------------------------------------------------
traj_population <- bind_rows(lapply(1:3, function(r)
  summarise_population_trajectory(grids[[r]], get_beta_samples(samples, r), r)
))
saveRDS(traj_population, file.path(PATH_OUTPUT, "trajectories_population.rds"))

# -----------------------------------------------------------------------------
# Area-level trajectories zeta_j^(r)(t) and deviations g_j^(r)(t)
# -----------------------------------------------------------------------------
traj_area <- bind_rows(lapply(1:3, function(r) {
  beta <- get_beta_samples(samples, r)            # n_samp x 3
  b    <- get_svc_samples(samples, r)             # list of 3 x J
  zeta <- lapply(seq_along(b), function(s) b[[s]] + beta[s, ])
  bind_rows(
    summarise_area_trajectories(grids[[r]], zeta, r, value = "zeta"),
    summarise_area_trajectories(grids[[r]], b,    r, value = "g")
  )
}))
saveRDS(traj_area, file.path(PATH_OUTPUT, "trajectories_area.rds"))
