# =============================================================================
#  04. Posterior predictive diagnostics of the selected model (Section 4.7)
#
#    - mid-PIT histogram (PIT corrected for discreteness)
#    - ACF of the Pearson residuals averaged over areas by week
#    - Moran's I of the Pearson residuals averaged over weeks by area
#
#  Output: output/fig_pit_histogram.pdf, output/fig_acf_residuals.pdf,
#          output/fig_resid_week.pdf
# =============================================================================

source("R/00_setup.R")
source("R/functions_model.R")

datos_fit <- readRDS(FILE_DATA)
load(FILE_ADJACENCY)                                  # adj_dep
fit <- readRDS(file.path(PATH_OUTPUT, "fit_1_poisson.rds"))$fit

# -----------------------------------------------------------------------------
# PIT. INLA returns the non-randomized PIT P(Y <= y | y_-i), which for counts
# is shifted towards one. The mid-PIT subtracts half the mass at y (the CPO).
# -----------------------------------------------------------------------------
stopifnot(!is.null(fit$cpo$pit))
mid_pit <- fit$cpo$pit - 0.5 * fit$cpo$cpo
message(sprintf("CPO failures: %d", sum(fit$cpo$failure > 0, na.rm = TRUE)))

pdf(file.path(PATH_OUTPUT, "fig_pit_histogram.pdf"), width = 6, height = 5)
hist(mid_pit, breaks = 20, main = "", xlab = "mid-PIT", ylab = "Frequency",
     col = "grey80", border = "white")
dev.off()

# -----------------------------------------------------------------------------
# Pearson residuals
# -----------------------------------------------------------------------------
mu_hat        <- fit$summary.fitted.values$mean
resid_pearson <- (datos_fit$suma_zona - mu_hat) / sqrt(mu_hat)

# Temporal: weekly mean over areas
resid_week <- tapply(resid_pearson, datos_fit$semana, mean)

pdf(file.path(PATH_OUTPUT, "fig_acf_residuals.pdf"), width = 6, height = 5)
acf(resid_week, main = "", lag.max = 25)
dev.off()

weeks <- as.numeric(names(resid_week))
pdf(file.path(PATH_OUTPUT, "fig_resid_week.pdf"), width = 10, height = 3.5)
plot(weeks, resid_week, type = "l", xlab = "Week", ylab = "Mean Pearson residual")
abline(h = 0, lty = 2)
abline(v = T0, lty = 3)
lines(lowess(weeks, resid_week), lwd = 2)
dev.off()

# Spatial: area mean over weeks. zero.policy = TRUE because one area has no
# neighbours in the adjacency matrix.
resid_area <- tapply(resid_pearson, datos_fit$id_zona, mean)
listw      <- spdep::mat2listw(as.matrix(adj_dep), style = "W",
                               zero.policy = TRUE)
print(spdep::moran.test(resid_area, listw, zero.policy = TRUE))
