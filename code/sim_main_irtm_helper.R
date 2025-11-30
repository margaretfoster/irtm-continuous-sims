## sim_irtm_helper.R
## Control script to compare the error in the 
## continuous case for the irt_m() helper 
## without anchor points
## from the cases we already know

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

##!! TODO: Change the helpers to the irtm helper version
source("monitornew.R")
source("helpers.R")    # defines mse, coverage, geweke_fcn, ess_fcn, etc.
source("sim_one_call_irtmhelper.R") #generates Y according to passed-in params and calls IRT-M for 
source("run_model_irtmhelper.R") #calls IRT-M with params passed in by sim_one_call.
source("sim_utils.R") # takes irtm MCMC results specified in run model, handles processing for metrics
source("irtm_cont_results_to_df.R") ## flatten results list to df:

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

## list of supported methods:
# extend by adding new methods to run_iter_cont()
## and updating below:
sim_methods <- c(
  "irtM",              # Old: M_constrained_irt_continuous with Yall+anchors
  "irtM_wrapper",      # New: irt_m() wrapper, no anchors
  "irtMAnchors",       # Old: M + fixed anchors
  "benchmark"
)

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
  FUN         = run_one_sim, #insim_one_call.R
  control     = sim_control,
  mc.cores    = n_cores,
  mc.set.seed = TRUE
)

## How many runs errored out:
num_errors <- sum(sapply(all_res, function(x) "error" %in% names(x)))
num_errors 

##------- Save raw results -------##
save(all_res,
     file = file.path(simpath, "irtm_cont_no_yfake.Rds"))

## Tidydf form:
results_df <- irtm_cont_results_to_df(all_res)
saveRDS(results_df,
        file = file.path(simpath, "irtm_cont_df_no_yfake.rds"))
