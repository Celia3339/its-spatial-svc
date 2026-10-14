# =============================================================================
#  05. Robustness checks of the model specification (Section 2)
#
#  (a) Overdispersion of the raw counts and of a Poisson GLM with fixed effects
#      (variance-to-mean ratio d and Pearson dispersion phi).
#  (b) Poisson versus negative binomial under two specifications, fixed
#      effects only and fixed effects + AR(1) + Leroux
#      (overdispersion table of the Supplementary Material).
#  (c) Temporal specification: the full model augmented with a secular trend
#      (RW2 on t) and a yearly harmonic, fitted alongside the AR(1) process,
#      compared with the base model (Delta DIC and change in rho).
#
#  Checks (b) and (c) refit the model and take several hours.
#  Output: output/robustness_overdispersion.csv, output/robustness_trend.csv
# =============================================================================

source("R/00_setup.R")
source("R/functions_model.R")
set_threads()

datos_fit <- readRDS(FILE_DATA)
load(FILE_ADJACENCY)                                   # adj_dep
adj_dep <- Matrix::Matrix(adj_dep, sparse = TRUE)
datos_fit <- prepare_spatial_indices(datos_fit)
datos_fit <- prepare_spline_basis(datos_fit, df = 3)

controls <- define_inla_controls(strategy = "simplified.laplace", cpo = FALSE)
offset   <- log(datos_fit$poblacion)

summarise_fit <- function(fit, label, family) {
  data.frame(Specification = label, Family = family,
             DIC  = round(fit$dic$dic, 1),
             pD   = round(fit$dic$p.eff, 1),
             WAIC = round(fit$waic$waic, 1))
}

# -----------------------------------------------------------------------------
# (a) Overdispersion measures
# -----------------------------------------------------------------------------
y <- datos_fit$suma_zona
d_raw <- var(y) / mean(y)

glm_pois <- glm(build_formula_fixed(), family = poisson,
                offset = log(poblacion), data = datos_fit)
phi_glm <- sum(residuals(glm_pois, type = "pearson")^2) / df.residual(glm_pois)

message(sprintf("Variance-to-mean ratio d = %.2f | Pearson dispersion phi = %.2f",
                d_raw, phi_glm))

# -----------------------------------------------------------------------------
# (b) Poisson versus negative binomial
# -----------------------------------------------------------------------------
hyp <- define_hyperpriors()

formula_re <- update(
  build_formula_fixed(),
  . ~ . + f(id_zona, model = "besagproper2", graph = adj_dep,
            hyper = hyp$space_leroux) +
    f(semana, model = "ar1", hyper = hyp$temporal_ar1)
)

specs_b <- list(
  list(label = "Fixed effects only",            formula = build_formula_fixed()),
  list(label = "Random effects (AR(1) + Leroux)", formula = formula_re)
)

tab_overdisp <- do.call(rbind, lapply(specs_b, function(sp) {
  do.call(rbind, lapply(c("poisson", "nbinomial"), function(fam) {
    fit <- INLA::inla(
      sp$formula, family = fam, data = datos_fit, offset = offset,
      control.compute   = controls$compute,
      control.predictor = controls$predictor,
      control.inla      = controls$inla,
      control.family    = list(hyper = define_hyperpriors(family = fam)$family)
    )
    summarise_fit(fit, sp$label, fam)
  }))
}))
tab_overdisp <- tab_overdisp %>%
  group_by(Specification) %>%
  mutate(Delta_DIC = DIC - min(DIC)) %>%
  ungroup()
print(tab_overdisp)
write.csv(tab_overdisp, file.path(PATH_OUTPUT, "robustness_overdispersion.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# (c) AR(1) versus AR(1) + secular trend + yearly harmonic
# -----------------------------------------------------------------------------
# NOTE: set RW2_PRIOR to the prior used in the original analysis.
RW2_PRIOR <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
PERIOD    <- 52.18

datos_fit <- datos_fit %>%
  mutate(semana_rw2 = semana,
         harm_sin   = sin(2 * pi * semana / PERIOD),
         harm_cos   = cos(2 * pi * semana / PERIOD))

formula_base  <- build_formula_full(hyp, adj_dep)
# RW2_PRIOR is found in the global environment when INLA evaluates the formula.
formula_trend <- update(
  formula_base,
  . ~ . + f(semana_rw2, model = "rw2", scale.model = TRUE, hyper = RW2_PRIOR) +
    harm_sin + harm_cos
)

fit_base  <- readRDS(file.path(PATH_OUTPUT, "fit_1_poisson.rds"))$fit
fit_trend <- fit_inla_twostage(formula_trend, datos_fit, family = "poisson",
                               offset_var = "poblacion", hypers = hyp,
                               controls = controls)$fit
saveRDS(fit_trend, file.path(PATH_OUTPUT, "fit_trend_rw2_harmonic.rds"))

get_rho <- function(fit)
  fit$summary.hyperpar[grep("Rho for semana", rownames(fit$summary.hyperpar),
                            fixed = TRUE), "mean"]

tab_trend <- data.frame(
  Model     = c("AR(1)", "AR(1) + RW2 + harmonic"),
  DIC       = round(c(fit_base$dic$dic, fit_trend$dic$dic), 1),
  WAIC      = round(c(fit_base$waic$waic, fit_trend$waic$waic), 1),
  rho       = round(c(get_rho(fit_base), get_rho(fit_trend)), 3)
)
tab_trend$Delta_DIC <- tab_trend$DIC - tab_trend$DIC[1]
print(tab_trend)
write.csv(tab_trend, file.path(PATH_OUTPUT, "robustness_trend.csv"),
          row.names = FALSE)
