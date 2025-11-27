library(MASS)
library(mvtnorm)
library(coda)
library(nimble)
library(pscl)
library(sn) #skew on Y

library(parallel)
n_cores <- detectCores() - 1  # leave 1 core free

## Note: before IRTM continuous resubmitted to CRAN load locally via:
devtools::load_all("~/Dropbox/IRTM-Pkg/IRTM")
#ls("package:IRTM")# to verify that there is a continuous function

source('monitornew.R')
source('helpers.R') #defines mse, sbc, traceplots, etc

## paths:
simpath = "./simulations/"

## Set of possible methods supported by sim engine
all_methods <- c(
  "benchmark",
  "unconstrained",
  "irtM",
  "irtMNoVar",
  "irtAnchors",
  "irtMAnchors"
)

## subset for this run:
sim_methods <- c("irtM", 
                 "irtMAnchors", 
                 "benchmark")



##---- sims grid ------#####
make_sim_grid <- function(N_vals,
                          K_vals,
                          d_vals,
                          dist_vals,
                          reps = 5) {
  
  grid <- expand.grid(
    N   = N_vals,
    K   = K_vals,
    d   = d_vals,
    dist_type =  dist_vals,
    rep = seq_len(reps)
  )
  
  grid$sim_id <- seq_len(nrow(grid))
  return(grid)
}

## Initial diagnostic:
medium_sim <- make_sim_grid(
  N_vals = c(50, 100, 500),   # small vs moderate N
  K_vals = c(10, 50),       # light vs moderate item
  d_vals = c(2, 5),         # low vs moderate dimensions
  reps   = 10,               # 10 replicates per cell
  dist_vals = c("normal", 
                 "heavy_t",
                 "skewed",
                 "mixture",
                 "bounded",
                 "log_normal",
                 "power_law") # types of distributions
)

## Declare param grid that goes into the script:

param_grid = medium_sim

## Parameters:
## types of distribution

nsamp <- 1000
nburn <- 2000
thin <- 2
nchains <- 1 #number of independent MCMC chains to run
## note that nchains > 1 will require adding a wrapper to run and average multiple chains
lambda_zero_pct = 0

##basic simulation:
source("quick_cont.R")

## Save output:
saveRDS(all_res, 
        file = paste0(simpath, "small_cont.rds"))

