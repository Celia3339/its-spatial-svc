# =============================================================================
#  03. Simulation study (Section 5 of the article)
#
#  Part 1: dependence scenarios A, B and C (S = 200 each).
#  Part 2 (optional): effect-size sensitivity under Scenario B, scaling all
#          intervention effects and SVCs by alpha in {0, 0.5, 2}
#          (Supplementary Material, S = 100).
#
#  Every replicate is saved as output/simulation/rep_<scenario>_<i>.rds, so an
#  interrupted run resumes where it stopped.
#
#  Inputs: the adjacency matrix and the area-level covariates and population
#  of the real data (FILE_ADJACENCY, FILE_DATA in R/00_setup.R).
#
#  Authors: Celia Robles, Adina Iftimi
# =============================================================================

source("R/00_setup.R")
source("R/functions_simulation.R")

# -----------------------------------------------------------------------------
# Run settings
# -----------------------------------------------------------------------------
TEST_MODE        <- FALSE    # TRUE: a few replicates to check the pipeline
S_MAIN           <- if (TEST_MODE) 5L else 200L
S_EFFECT         <- if (TEST_MODE) 5L else 100L
RUN_EFFECT_SIZE  <- FALSE
N_WORKERS        <- 3L       # parallel replicates (Linux/macOS only)
N_POST_SAMPLE    <- 500L     # posterior samples for the spline covariance
J_sim            <- 100L     # areas in the simulated subgraph
T_WEEKS          <- 205L     # the simulation uses the first 205 weeks

OUTPUT_SIM <- file.path(PATH_OUTPUT, if (TEST_MODE) "simulation_test" else "simulation")
dir.create(OUTPUT_SIM, recursive = TRUE, showWarnings = FALSE)

n_cores <- set_threads()
inla.setOption(num.threads = n_cores)

# -----------------------------------------------------------------------------
# True values: posterior means of the model fitted to the real data
# -----------------------------------------------------------------------------
beta_true <- c(
  intercept  = -6.709,
  int1       = -0.094, int2 = -0.045, int3 = -0.307,
  bs_post1_1 =  0.047, bs_post1_2 = -0.014, bs_post1_3 =  0.014,
  bs_post2_1 = -0.318, bs_post2_2 = -0.430, bs_post2_3 = -0.299,
  bs_post3_1 = -0.096, bs_post3_2 =  0.023, bs_post3_3 =  0.189,
  age65      =  0.554, homb_porc  = -2.875, tipo_suelo =  0.013
)
beta_spline1 <- beta_true[paste0("bs_post1_", 1:3)]
beta_spline2 <- beta_true[paste0("bs_post2_", 1:3)]
beta_spline3 <- beta_true[paste0("bs_post3_", 1:3)]

# SVC precisions, identical across scenarios
kappa_svc <- list(
  int1  = 1.915, int2 = 0.541, int3 = 3.476,
  post1 = c(0.287, 0.431, 0.338),
  post2 = c(1.963, 2.145, 1.681),
  post3 = c(0.762, 0.354, 0.735)
)

# Dependence scenarios of the baseline spatial and temporal effects
scenarios <- list(
  A = list(kappa_phi =  5.0, lambda = 0.85, rho = 0.40, sigma_eps = 0.25),
  B = list(kappa_phi = 15.0, lambda = 0.85, rho = 0.60, sigma_eps = 0.20),
  C = list(kappa_phi = 50.0, lambda = 0.85, rho = 0.80, sigma_eps = 0.15)
)

effect_scenarios <- list(
  null   = list(label = "Null",   beta_factor = 0,   spline_factor = 0,   svc_factor = 1),
  half   = list(label = "Half",   beta_factor = 0.5, spline_factor = 0.5, svc_factor = 1),
  double = list(label = "Double", beta_factor = 2,   spline_factor = 2,   svc_factor = 1)
)

# -----------------------------------------------------------------------------
# Priors and INLA settings used to fit each synthetic data set.
# NOTE: these are the values used to produce the results in the article.
# The SVC and Leroux precisions use pc.prec(1, 0.01), which differs from the
# real-data analysis (SVC: pc.prec(0.1, 0.01); Leroux: pc.prec(1, 0.6)).
# -----------------------------------------------------------------------------
hyper_sim <- list(
  leroux = list(prec   = list(prior = "pc.prec",  param = c(1, 0.01)),
                lambda = list(prior = "gaussian", param = c(0, 2))),
  ar1    = list(prec = list(prior = "pc.prec", param = c(0.2, 0.01)),
                rho  = list(prior = "pc.cor0", param = c(0.7, 0.05))),
  svc    = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
)

control_inla_sim <- list(
  strategy               = "simplified.laplace",
  int.strategy           = "eb",
  stupid.search.max.iter = 1500,
  h                      = 0.01
)

# -----------------------------------------------------------------------------
# Design: time, indicators and spline bases
# -----------------------------------------------------------------------------
t_vec <- seq_len(T_WEEKS)
ind   <- make_indicators(t_vec, T0)
I1 <- ind$I1; I2 <- ind$I2; I3 <- ind$I3

counter <- function(r) ifelse(t_vec >= T0[r], t_vec - T0[r] + 1L, 0L)
P1 <- ifelse(I2 == 1, 0L, counter(1))
P2 <- ifelse(I3 == 1, 0L, counter(2))
P3 <- counter(3)
B1 <- ns(P1, df = 3)
B2 <- ns(P2, knots = c(6, 12))
B3 <- ns(P3, knots = c(16, 31))

# -----------------------------------------------------------------------------
# Spatial structure: connected subgraph of J_sim areas
# -----------------------------------------------------------------------------
load(FILE_ADJACENCY)                                   # adj_dep
adj_dep   <- Matrix::Matrix(adj_dep, sparse = TRUE)
datos_fit <- readRDS(FILE_DATA)

areas_sub <- select_connected_areas(adj_dep, J_sim, seed = 42)
stopifnot(components(graph_from_adjacency_matrix(adj_dep[areas_sub, areas_sub],
                                                 mode = "max"))$no == 1L)

W_sim_sparse <- adj_dep[areas_sub, areas_sub]
Q_leroux_sim <- Diagonal(x = rowSums(W_sim_sparse)) - W_sim_sparse

graph_file <- tempfile(fileext = ".graph")
nb_sim <- mat2listw(as.matrix(W_sim_sparse), style = "B",
                    zero.policy = TRUE)$neighbours
class(nb_sim) <- "nb"
nb2INLA(graph_file, nb_sim)
adj_sim <- inla.read.graph(graph_file)

icar     <- scaled_icar_eigen(W_sim_sparse)
V_icar   <- icar$V
lam_icar <- icar$lam
message(sprintf("ICAR scaling constant: %.4f", icar$scale))

# Area-level covariates and population (offset)
age65_sim      <- tapply(datos_fit$mayor65_porc, datos_fit$id_zona, mean)[areas_sub]
homb_porc_sim  <- tapply(datos_fit$homb_porc,    datos_fit$id_zona, mean)[areas_sub]
tipo_suelo_sim <- tapply(datos_fit$tipo_suelo,   datos_fit$id_zona,
                         function(x) as.integer(names(which.max(table(x)))))[areas_sub]
E_sim <- unname(tapply(datos_fit$poblacion,
                       list(datos_fit$id_zona, datos_fit$semana),
                       mean)[areas_sub, seq_len(T_WEEKS)])
rm(datos_fit); gc()

# -----------------------------------------------------------------------------
# Part 1: dependence scenarios
# -----------------------------------------------------------------------------
seed_base <- c(A = 1000L, B = 2000L, C = 3000L)
for (sc in names(scenarios)) {
  run_scenario(sc, S_MAIN, seed_base[[sc]],
               gen_fun    = function() generate_dataset(scenarios[[sc]]),
               output_dir = OUTPUT_SIM, n_workers = N_WORKERS, n_cores = n_cores)
}

metrics_main <- lapply(setNames(names(scenarios), names(scenarios)), function(sc) {
  true_vals <- c(list(beta = beta_true), scenarios[[sc]])
  extract_metrics(load_scenario(sc, OUTPUT_SIM), true_vals)
})
saveRDS(metrics_main, file.path(OUTPUT_SIM, "metrics_main.rds"))

# EmpSE and relative error of ModSE of the non-spline fixed effects
fixed_non_spline <- c("intercept", "int1", "int2", "int3",
                      "age65", "homb_porc", "tipo_suelo")
tab_empse <- do.call(cbind, lapply(names(metrics_main), function(sc) {
  m <- metrics_main[[sc]]$beta_metrics[fixed_non_spline, c("EmpSE", "RelErrModSE")]
  colnames(m) <- paste(colnames(m), sc, sep = "_")
  m
}))
print(round(tab_empse, 3))

# -----------------------------------------------------------------------------
# Part 2 (optional): effect-size sensitivity under Scenario B
# -----------------------------------------------------------------------------
if (RUN_EFFECT_SIZE) {
  seed_base_eff <- c(null = 5000L, half = 6000L, double = 7000L)
  for (en in names(effect_scenarios)) {
    e <- effect_scenarios[[en]]
    run_scenario(paste0("eff_", en), S_EFFECT, seed_base_eff[[en]],
                 gen_fun    = function() generate_dataset(scenarios$B, effect = e),
                 output_dir = OUTPUT_SIM, n_workers = N_WORKERS, n_cores = n_cores)
  }

  metrics_effect <- lapply(setNames(names(effect_scenarios),
                                    names(effect_scenarios)), function(en) {
    e    <- effect_scenarios[[en]]
    beta <- beta_true
    beta[c("int1", "int2", "int3")] <- beta_true[c("int1", "int2", "int3")] *
      e$beta_factor
    spl <- grep("^bs_post", names(beta))
    beta[spl] <- beta_true[spl] * e$spline_factor
    extract_metrics(load_scenario(paste0("eff_", en), OUTPUT_SIM),
                    c(list(beta = beta), scenarios$B))
  })
  metrics_effect$observed <- metrics_main$B
  saveRDS(metrics_effect, file.path(OUTPUT_SIM, "metrics_effect.rds"))
}
