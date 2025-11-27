## sim_main.R
## Top-level control script for IRT-M continuous simulations

library(MASS)
library(mvtnorm)
library(coda)
library(nimble)
library(pscl)
library(sn)
library(parallel)
library(tibble)
library(dplyr)

## Note: before IRTM continuous resubmitted to CRAN load locally via:
devtools::load_all("~/Dropbox/IRTM-Pkg/IRTM")

source("monitornew.R")
source("helpers.R")    # defines mse, coverage, geweke_fcn, ess_fcn, etc.
source("sim_utils.R")
source("run_model.R")
source("sim_one_call.R")

## paths:
simpath  <- "../simulations/"
plotpath <- "../simulations/results/figures/"

## random number generator for parallelism
RNGkind("L'Ecuyer-CMRG")
set.seed(112625)

##------- Simulation grid helper -------##
make_sim_grid <- function(N_vals,
                          K_vals,
                          d_vals,
                          dist_vals,
                          reps = 5) {
  grid <- expand.grid(
    N         = N_vals,
    K         = K_vals,
    d         = d_vals,
    dist_type = dist_vals,
    sim       = seq_len(reps)
  )
  grid$sim_id <- seq_len(nrow(grid))
  grid
}

## Grid
medium_sim <- make_sim_grid(
  N_vals   = c(50, 100, 500),
  K_vals   = c(10, 50),
  d_vals   = c(2, 5),
  dist_vals = c(
    "normal",
    "heavy_t",
    "skewed",
    "mixture",
    "bounded",
    "log_normal",
    "power_law"
  ),
  reps = 10
)

param_grid <- medium_sim

##------- Global simulation control -------##
nsamp   <- 1000
nburn   <- 2000
thin    <- 2
nchains <- 1

lambda_zero_pct <- 0  # fraction of lambda entries zeroed

## Full method universe:
# extend by adding new methods to run_iter_cont()
## and updating below:
all_methods <- c(
  "benchmark",
  "unconstrained",
  "irtM",
  "irtMNoVar",
  "irtAnchors",
  "irtMAnchors"
)

## Subset for run:
sim_methods <- c("irtM", "irtMAnchors", "benchmark")
# or:
# sim_methods <- all_methods

sim_control <- list(
  nsamp = nsamp,
  nburn = nburn,
  thin  = thin,
  nchains = nchains,
  lambda_zero_pct = lambda_zero_pct,
  methods = sim_methods
)

##------- Convert param_grid to list-of-lists for mclapply -------##
param_list <- lapply(seq_len(nrow(param_grid)), function(i) {
  list(
    N         = param_grid$N[i],
    K         = param_grid$K[i],
    d         = param_grid$d[i],
    sim       = param_grid$sim[i],
    dist_type = param_grid$dist_type[i]
  )
})

##------- Parallel execution -------##
n_cores <- detectCores() - 1
n_cores <- max(1, n_cores)

all_res <- mclapply(
  X           = param_list,
  FUN         = run_one_sim,
  control     = sim_control,
  mc.cores    = n_cores,
  mc.set.seed = TRUE
)

##------- Save raw results -------##
save(all_res,
     file = file.path(simpath, "irtm_cont_parallel.Rds"))

## Tidydf form:
results_df <- irtm_cont_results_to_df(all_res)
saveRDS(results_df,
        file = file.path(simpath, "irtm_cont_parallel_df.rds"))
