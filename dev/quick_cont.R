## Quick script for continuous simulation
## using the architecture of the binary simulation

# change seed algorithm for one that handles parallism
RNGkind("L'Ecuyer-CMRG") 
set.seed(112625)

## helper for computing metrics:

#------------------------------------------------------------------
# Helper to update metrics for a given method
#------------------------------------------------------------------
update_method_metrics <- function(method_name,
                                  theta_avg, lambda_avg,
                                  true_theta, true_lambda,
                                  irt_res,
                                  theta_mse, lambda_mse,
                                  theta_coverage, lambda_coverage,
                                  theta_gew, theta_ess, theta_rhat,
                                  lambda_corr) {
  
  theta_mse[method_name]        <- mse(theta_avg,  true_theta, TRUE)
  lambda_mse[method_name]       <- mse(lambda_avg, true_lambda, TRUE)
  
  theta_coverage[method_name]   <- coverage(theta_avg,  true_theta)
  lambda_coverage[method_name]  <- coverage(lambda_avg, true_lambda)
  
  theta_gew[method_name]        <- geweke_fcn(irt_res$theta)
  theta_ess[method_name]        <- ess_fcn(irt_res$theta)
  theta_rhat[method_name]       <- NA_real_
  
  # correlation across all entries
  lambda_corr[method_name] <- cor(
    as.numeric(lambda_avg),
    as.numeric(true_lambda),
    use = "complete.obs"
  )
  
  return(list(
    theta_mse       = theta_mse,
    lambda_mse      = lambda_mse,
    theta_coverage  = theta_coverage,
    lambda_coverage = lambda_coverage,
    theta_gew       = theta_gew,
    theta_ess       = theta_ess,
    theta_rhat      = theta_rhat,
    lambda_corr     = lambda_corr
  ))
}


## Simulation:

run_iter_cont <- function(
    Y, Yall,
    true_theta, true_theta_all, true_lambda,
    M, x, z, ind, d,
    progress = TRUE,
    methods = c('benchmark','unconstrained','irtM','irtMNoVar','irtAnchors','irtMAnchors')
){
  
  nm <- methods
  
  # Preallocate
  theta_mse       <- setNames(rep(NA_real_, length(nm)), nm)
  lambda_mse      <- setNames(rep(NA_real_, length(nm)), nm)
  theta_coverage  <- setNames(rep(NA_real_, length(nm)), nm)
  lambda_coverage <- setNames(rep(NA_real_, length(nm)), nm)
  theta_gew       <- setNames(rep(NA_real_, length(nm)), nm)
  theta_ess       <- setNames(rep(NA_real_, length(nm)), nm)
  theta_rhat      <- setNames(rep(NA_real_, length(nm)), nm)
  lambda_corr     <- setNames(rep(NA_real_, length(nm)), nm)
  
  #--------------------------------------#
  # Benchmark
  #--------------------------------------#
  if ("benchmark" %in% methods) {
    print("Running benchmark...")
    irt_res <- M_constrained_irt_continuous(
      Y, d, M=NULL,
      theta_fix = z,
      which_fix = ind,
      nburn = nburn, nsamp = nsamp, thin = 1,
      learn_Sigma = FALSE, learn_Omega = FALSE,
      display_progress = progress
    )
    
    updated <- update_method_metrics(
      method_name = "benchmark",
      theta_avg   = irt_res$theta,
      lambda_avg  = irt_res$lambda,
      true_theta  = true_theta,
      true_lambda = true_lambda,
      irt_res     = irt_res,
      theta_mse, lambda_mse,
      theta_coverage, lambda_coverage,
      theta_gew, theta_ess, theta_rhat,
      lambda_corr
    )
    
    list2env(updated, envir = environment())
    rm(irt_res)
  }
  
  #--------------------------------------#
  # Unconstrained
  #--------------------------------------#
  if ("unconstrained" %in% methods) {
    print("Running unconstrained...")
    irt_res <- M_constrained_irt_continuous(
      Y, d, M=NULL, theta_fix=NULL, which_fix=NULL,
      nburn = nburn, nsamp = nsamp, thin = 1,
      learn_Sigma = FALSE, learn_Omega = FALSE,
      display_progress = progress
    )
    
    updated <- update_method_metrics(
      "unconstrained",
      theta_avg   = irt_res$theta,
      lambda_avg  = irt_res$lambda,
      true_theta  = true_theta,
      true_lambda = true_lambda,
      irt_res     = irt_res,
      theta_mse, lambda_mse,
      theta_coverage, lambda_coverage,
      theta_gew, theta_ess, theta_rhat,
      lambda_corr
    )
    
    list2env(updated, envir = environment())
    rm(irt_res)
  }
  
  #--------------------------------------#
  # irtM (full M constraint)
  #--------------------------------------#
  if ("irtM" %in% methods) {
    print("Running M-constrained IRT...")
    irt_res <- M_constrained_irt_continuous(
      Yall, d, M=M,
      theta_fix=NULL, which_fix=NULL,
      nburn = nburn, nsamp = nsamp, thin = 1,
      learn_Sigma = TRUE, display_progress = progress
    )
    
    updated <- update_method_metrics(
      "irtM",
      theta_avg   = irt_res$theta,
      lambda_avg  = irt_res$lambda,
      true_theta  = true_theta_all,
      true_lambda = true_lambda,
      irt_res     = irt_res,
      theta_mse, lambda_mse,
      theta_coverage, lambda_coverage,
      theta_gew, theta_ess, theta_rhat,
      lambda_corr
    )
    
    list2env(updated, envir = environment())
    rm(irt_res)
  }
  
  #--------------------------------------#
  # irtM No Variance Prior
  #--------------------------------------#
  if ("irtMNoVar" %in% methods) {
    print("Running M-constrained IRT without variance prior...")
    irt_res <- M_constrained_irt_continuous(
      Yall, d, M=M,
      theta_fix=x, which_fix=NULL,
      nburn = nburn, nsamp = nsamp,
      thin = 1, learn_Sigma = FALSE,
      learn_Omega = FALSE,
      display_progress = progress
    )
    
    updated <- update_method_metrics(
      "irtMNoVar",
      theta_avg   = irt_res$theta,
      lambda_avg  = irt_res$lambda,
      true_theta  = true_theta_all,
      true_lambda = true_lambda,
      irt_res     = irt_res,
      theta_mse, lambda_mse,
      theta_coverage, lambda_coverage,
      theta_gew, theta_ess, theta_rhat,
      lambda_corr
    )
    
    list2env(updated, envir = environment())
    rm(irt_res)
  }
  
  #--------------------------------------#
  # Anchors only
  #--------------------------------------#
  if ("irtAnchors" %in% methods) {
    print("Running IRT with anchors...")
    irt_res <- M_constrained_irt_continuous(
      Yall, d, M=NULL,
      theta_fix=x, which_fix=ind,
      nburn = nburn, nsamp = nsamp,
      thin = 1, learn_Sigma = TRUE,
      display_progress = progress
    )
    
    updated <- update_method_metrics(
      "irtAnchors",
      theta_avg   = irt_res$theta,
      lambda_avg  = irt_res$lambda,
      true_theta  = true_theta_all,
      true_lambda = true_lambda,
      irt_res     = irt_res,
      theta_mse, lambda_mse,
      theta_coverage, lambda_coverage,
      theta_gew, theta_ess, theta_rhat,
      lambda_corr
    )
    
    list2env(updated, envir = environment())
    rm(irt_res)
  }
  
  #--------------------------------------#
  # IRTM with anchors + M
  #--------------------------------------#
  if ("irtMAnchors" %in% methods) {
    print("Running IRT with M and Anchors...")
    irt_res <- M_constrained_irt_continuous(
      Yall, d, M=M,
      theta_fix=x, which_fix=ind,
      nburn = nburn, nsamp = nsamp,
      thin = 1, learn_Sigma = TRUE,
      display_progress = progress
    )
    
    updated <- update_method_metrics(
      "irtMAnchors",
      theta_avg   = irt_res$theta,
      lambda_avg  = irt_res$lambda,
      true_theta  = true_theta_all,
      true_lambda = true_lambda,
      irt_res     = irt_res,
      theta_mse, lambda_mse,
      theta_coverage, lambda_coverage,
      theta_gew, theta_ess, theta_rhat,
      lambda_corr
    )
    
    list2env(updated, envir = environment())
    rm(irt_res)
  }
  
  return(list(
    theta_mse       = theta_mse,
    lambda_mse      = lambda_mse,
    theta_coverage  = theta_coverage,
    lambda_coverage = lambda_coverage,
    theta_gew       = theta_gew,
    theta_ess       = theta_ess,
    theta_rhat      = theta_rhat,
    lambda_corr     = lambda_corr
  ))
}


##----------------- Simulation loop ------------------- ##
all_res <- list()

## parallize the for-loop
## nclusters declared in sims control script (c_sims.R)
all_res <- mclapply(
  X = 1:nrow(param_grid),
  mc.cores = n_cores,
  mc.set.seed = TRUE,
  FUN = function(r) {
    try({
    
    start_time <- Sys.time() ## start timer
    ## ---- Extract simulation parameters ----
    N         <- param_grid$N[r]
    K         <- param_grid$K[r]
    d         <- param_grid$d[r]
    sim       <- param_grid$rep[r]
    dist_type <- param_grid$dist_type[r]
    
    ## report out sim:
    message(sprintf("Simulation %d of %d || N=%d K=%d d=%d dist=%s",
                    r, nrow(param_grid), N, K, d, dist_type))

  ## Generate latent covariance (Sigma)
  sigmacorr <- matrix(runif(d*d, -1, 1), d, d)
  sigmacorr[lower.tri(sigmacorr)] = t(sigmacorr)[lower.tri(sigmacorr)]
  true_Sigma =  cov2cor(sigmacorr%*%sigmacorr)
  mu = rep(0, d)

  ## Generate true theta, lambda, beta
  true_theta  <- rmvnorm(n = N, mean = mu, sigma = true_Sigma)
  true_lambda <- rmvnorm(n = K, mean = mu, sigma = diag(d))
  true_b      <- rnorm(n = K, mean = 0,  sd  = 1)

  # random set some lambdas to be zero
  num_zero_lambda <- floor(lambda_zero_pct * K * d)
  lambda_ind <- sample(1:(d*K), num_zero_lambda, replace = F)
  true_lambda[lambda_ind] <- 0

  # M matrix
  M = array(NA, c(d,d,K))
  for (k in 1:K) {
    M[,,k] = diag(d) * sign(true_lambda[k,])
  }

  # directions
  direction <- sign(true_lambda)

  # Create anchor points for Y and theta
  Flags  <- c(1, -1)
  points <- c()
  for (i in 1:d) {
    for (flag in Flags) {
      points <- c(points, ifelse( direction[ , i] == flag, 1, 0 ))
    }
  }
  fake_points <- data.frame(
    matrix(points, nrow = 2*d, ncol = K, byrow = TRUE),
    row.names = paste0("fake_", 1:(2*d))
  )

  # Create corresponding constraints on thetas
  A <- 10
  x <- list()
  j = 1
  for(k in 1:d){
    a = list()
    v1 = rep(0, d)
    v1[k] = A
    a[[as.character(j)]] = v1
    j = j + 1
    a[[as.character(j)]] = -v1
    j = j + 1
    x = c(x, a)
  }
  ind <- 1:(2*d)

  # Create benchmark constraints
  z <- list()
  j = 1
  for(k in 1:d){
    a = list()
    a[[as.character(j)]] = true_theta[j, ]
    j = j + 1
    a[[as.character(j)]] = true_theta[j, ]
    j = j + 1
    z = c(z, a)
  }

  true_theta_all <- rbind(matrix(unlist(x),
                                 length(x),
                                 d,
                                 byrow=TRUE),
                          true_theta)

  ##### ------------------------ Generate Y ----------------------------- #####
  ## code to produce different forms of continuous Y:
  generate_Y <- function(mu, dist_type) {
    switch(dist_type,
           
           "normal" = mu + rnorm(1, 0, 1),
           
           "heavy_t" = mu + rt(1, df = 3),
           
           "skewed" = {
             # skew-normal via sn package:
             sn::rsn(1, xi = mu, omega = 1, alpha = 5)
           },
           
           "bounded" = {
             val <- mu + rnorm(1,0,1)
             pmax(pmin(val, 3), -3)
           },
           
           "mixture" = {
             if(runif(1) < 0.9) {
               mu + rnorm(1, 0, 1)
             } else {
               mu + rnorm(1, 0, 5)
             }
           },
           
           "power_law" = {
             # Shifted Pareto-like heavy tail
             ## commonly found in social processes
             xm <- 1
             alpha <- 2.5
             x <- xm * (runif(1))^(-1/alpha)
             mu + x
           },
           "log_normal" = {
             ## better behaved than pareto
             ## also common in behavioral data
             sigma <- 1.5
             x <- rlnorm(1, meanlog = 0, sdlog = sigma)
             mu + x
           },
           
           stop("Unknown distribution type")
    )
  }
  
  
  true_Z = matrix(NA, N, K)
  Y = matrix(NA, N, K)
  for (i in 1:N) {
    for (k in 1:K) {
      mu_ik = t(true_lambda[k, ]) %*% true_theta[i, ] - true_b[k]
      true_Z[i, k] = rnorm(1, mean = mu_ik, sd = 1)
      Y[i, k] = generate_Y(mu_ik, dist_type)

    }
  }
  Yall <- rbind( as.matrix(fake_points), Y)

 sim_methods =c('irtM','irtMAnchors', 'benchmark')

  ## run sim:
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
   progress = FALSE,
   methods = c("irtM", "irtMAnchors", "benchmark")
 )
 gc() ## garbage collecting to free up memory
 
 end_time <- Sys.time() #end timer
 elapsed_seconds <- as.numeric(end_time - start_time, 
                               units="secs")
 
 #return param, results, time 
 iter_res <- list(
   N = N,
   K = K,
   d = d,
   sim = sim,
   dist_type = dist_type,
   lambda_zero_pct = lambda_zero_pct, #from control script
   n_fake = 2*d,
   nsamp = nsamp,
   nburn = nburn,
   thin = thin,
   nchains = nchains,
   methods = sim_methods,
   seed = .Random.seed, #diagnose failed runs
   results = run_results,
   runtime = elapsed_seconds
 )
 
 return(iter_res)
   }, silent = FALSE)
} #inside-parallel function
) ## mclapply for all_res


##--------------- Save results ------------------##
save(all_res, 
     file = paste0(simpath, "irtm_cont_parallel.Rds"))

## save time info
#tname = paste0("simulations/irtm_cont_N100K10_time.rds")
#save(model_times, file=tname)

