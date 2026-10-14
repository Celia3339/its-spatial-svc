# =============================================================================
#  Shared setup: packages, paths and study constants
#  Sourced at the top of every script in scripts/.
# =============================================================================

suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
  library(splines)
  library(dplyr)
  library(tidyr)
  library(spdep)
  library(igraph)
})

# -----------------------------------------------------------------------------
# Paths (edit to your local layout). The data are NOT distributed with this
# repository; see README.md for the expected structure of the input files.
# -----------------------------------------------------------------------------
PATH_DATA   <- "data"      # input data (not included)
PATH_OUTPUT <- "output"    # fitted models, posterior summaries, simulation

FILE_DATA      <- file.path(PATH_DATA, "analysis_data.rds")   # area-week data
FILE_ADJACENCY <- file.path(PATH_DATA, "adjacency.RData")     # object adj_dep

dir.create(PATH_OUTPUT, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Study constants
# -----------------------------------------------------------------------------
N_AREAS <- 241L                    # basic health areas
T_WEEKS <- 209L                    # study weeks (week 1 = 1 January 2015)
T0      <- c(53L, 146L, 164L)      # implementation weeks t_0^(r), r = 1, 2, 3

# -----------------------------------------------------------------------------
# Threads
# -----------------------------------------------------------------------------
set_threads <- function(n_threads = NULL) {
  if (is.null(n_threads)) n_threads <- max(1L, parallel::detectCores() - 1L)
  Sys.setenv(OMP_NUM_THREADS        = n_threads,
             OPENBLAS_NUM_THREADS   = n_threads,
             MKL_NUM_THREADS        = n_threads,
             VECLIB_MAXIMUM_THREADS = n_threads)
  invisible(n_threads)
}
