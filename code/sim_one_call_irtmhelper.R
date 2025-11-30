## sim_one_call_irtmhelper.R:

## modifies previous irtm simulation code
## by removing fake points and Yall
## because the irt_m() wrapper doesn't allow
## anchor points in the continuous family
## Generate data for one simulation cell
## Then run all requested results

run_one_sim <- function(param, control) {
  ## parallel worker scope begins
  tryCatch({
    
    ## reload packages for parallel workers:
    library(MASS)
    library(mvtnorm)
    library(coda)
    library(nimble)
    library(pscl)
    library(sn)
    
    ## TODO!: Change scripts to irtm helper verions!
    ## source script for workers:
    source("monitornew.R")
    source("helpers.R")    # defines mse, coverage, geweke_fcn, ess_fcn, etc.
    source("run_model_irtmhelper.R") #calls IRT-M with params passed in by sim_one_call.
    source("sim_utils.R") # generate_y(), post-MCMC processing
    source("irtm_cont_results_to_df.R") ## flatten results list to df:
    # param: list with N, K, d, dist_type, sim
    # control: list with nsamp, nburn, thin, nchains, lambda_zero_pct, methods
    
    N         <- param$N
    K         <- param$K
    d         <- param$d
    sim       <- param$sim
    dist_type <- param$dist_type
    dist_type <- as.character(dist_type) # type validate
    
    nsamp           <- control$nsamp
    nburn           <- control$nburn
    thin            <- control$thin
    nchains         <- control$nchains
    lambda_zero_pct <- control$lambda_zero_pct
    methods         <- control$methods
    
    start_time <- Sys.time()
    
    ## latent covariance (Sigma)
    sigmacorr <- matrix(runif(d * d, -1, 1), d, d)
    sigmacorr[lower.tri(sigmacorr)] <- t(sigmacorr)[lower.tri(sigmacorr)]
    true_Sigma <- cov2cor(sigmacorr %*% sigmacorr)
    mu         <- rep(0, d)
    
    ## true theta, lambda, b
    true_theta  <- mvtnorm::rmvnorm(n = N, mean = mu, sigma = true_Sigma)
    true_lambda <- mvtnorm::rmvnorm(n = K, mean = mu, sigma = diag(d))
    true_b      <- rnorm(n = K, mean = 0, sd = 1)
    
    # zero some lambdas, if needed for comparisons:
    num_zero_lambda <- floor(lambda_zero_pct * K * d)
    if (num_zero_lambda > 0) {
      lambda_ind <- sample(1:(d * K),
                           num_zero_lambda,
                           replace = FALSE)
      true_lambda[lambda_ind] <- 0
    }
    
    # M matrix
    M <- array(NA_real_, c(d, d, K))
    for (k in 1:K) {
      M[, , k] <- diag(d) * sign(true_lambda[k, ])
    }
    
    # directions
    direction <- sign(true_lambda)
    
   
    # constraints on thetas (x)
    A <- 10
    x <- list()
    j <- 1
    for (k in 1:d) {
      a  <- list()
      v1 <- rep(0, d)
      v1[k] <- A
      a[[as.character(j)]] <- v1
      j <- j + 1
      a[[as.character(j)]] <- -v1
      j <- j + 1
      x <- c(x, a)
    }
    ind <- 1:(2 * d)
    
    # benchmark constraints (z)
    z <- list()
    j <- 1
    for (k in 1:d) {
      a <- list()
      a[[as.character(j)]] <- true_theta[j, ]
      j <- j + 1
      a[[as.character(j)]] <- true_theta[j, ]
      j <- j + 1
      z <- c(z, a)
    }
    
    # full theta
    true_theta_all <- rbind(
      matrix(unlist(x), length(x), d, byrow = TRUE),
      true_theta
    )
    
    # create Y
    true_Z <- matrix(NA_real_, N, K)
    Y      <- matrix(NA_real_, N, K)
    for (i in 1:N) {
      for (k in 1:K) {
        mu_ik        <- t(true_lambda[k, ]) %*% true_theta[i, ] - true_b[k]
        true_Z[i, k] <- rnorm(1, mean = mu_ik, sd = 1)
        Y[i, k]      <- generate_Y(mu_ik, dist_type)
      }
    }
    
    # run model(s)
    run_results <- run_iter_cont(
      Y = Y,
      Yall = Yall,
      true_theta = true_theta,
      true_theta_all = true_theta_all,
      true_lambda = true_lambda,
      M = M,
      x = x,
      z = z,
      ind = ind,
      d = d,
      nburn = nburn,
      nsamp = nsamp,
      thin = thin,
      nchains = nchains,
      progress = FALSE,
      methods = methods
    )
    
    gc()
    
    end_time        <- Sys.time()
    elapsed_seconds <- as.numeric(end_time - start_time, units = "secs")
    
    iter_res <- list(
      N = N,
      K = K,
      d = d,
      sim = sim,
      dist_type = dist_type,
      lambda_zero_pct = lambda_zero_pct,
      n_fake = 2 * d,
      nsamp = nsamp,
      nburn = nburn,
      thin = thin,
      nchains = nchains,
      methods = methods,
      seed = .Random.seed,
      results = run_results,
      runtime = elapsed_seconds
    )
    
    iter_res
  },
  
  error = function(e) {
    list(
      error = TRUE,
      message = e$message,
      call = deparse(e$call),
      param = param,
      control = control,
      seed = .Random.seed,
      runtime = NA_real_
    )
  })
}
