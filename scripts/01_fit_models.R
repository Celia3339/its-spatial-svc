# =============================================================================
#  01. Fit the spatio-temporal ITS model
#
#  Fits the full model under the Poisson and negative binomial likelihoods,
#  both with a Gaussian(0, 2) prior on logit(lambda), and builds the model
#  comparison table of the Supplementary Material. The model retained in the
#  article is fit_1 (Poisson).
#
#  Output: output/fit_1_poisson.rds, output/fit_2_nbinomial.rds
#          output/model_comparison.csv
# =============================================================================

source("R/00_setup.R")
source("R/functions_model.R")
set_threads()
set.seed(1234)

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
datos_fit <- readRDS(FILE_DATA)
load(FILE_ADJACENCY)                          # adj_dep: N_AREAS x N_AREAS
adj_dep <- Matrix::Matrix(adj_dep, sparse = TRUE)

datos_fit <- prepare_spatial_indices(datos_fit)
datos_fit <- prepare_spline_basis(datos_fit, df = 3)

# -----------------------------------------------------------------------------
# Specifications
# -----------------------------------------------------------------------------
specs <- list(
  fit_1_poisson   = list(family = "poisson",   label = "Poisson"),
  fit_2_nbinomial = list(family = "nbinomial", label = "Negative Binomial")
)

controls <- define_inla_controls(strategy = "simplified.laplace")
lincomb  <- create_intervention_lincomb(N_AREAS)

# -----------------------------------------------------------------------------
# Fit
# -----------------------------------------------------------------------------
fits <- list()
for (nm in names(specs)) {
  s      <- specs[[nm]]
  hypers <- define_hyperpriors(family = s$family)
  fits[[nm]] <- fit_inla_twostage(
    formula    = build_formula_full(hypers, adj_dep),
    datos      = datos_fit,
    family     = s$family,
    offset_var = "poblacion",
    hypers     = hypers,
    controls   = controls,
    lincomb    = lincomb
  )
  saveRDS(fits[[nm]], file.path(PATH_OUTPUT, paste0(nm, ".rds")))
}

# -----------------------------------------------------------------------------
# Model comparison (Supplementary Material)
# -----------------------------------------------------------------------------
comparison <- compare_models(fits, sapply(specs, `[[`, "label"))
print(comparison)
write.csv(comparison, file.path(PATH_OUTPUT, "model_comparison.csv"),
          row.names = FALSE)
