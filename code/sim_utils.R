
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

update_method_metrics <- function(
    method_name,
    irt_res,#irt-m sampler ouput
    true_theta, true_lambda,
    theta_mse, lambda_mse,
    theta_coverage, lambda_coverage,
    theta_gew, theta_ess, theta_rhat,
    lambda_corr
) {
  
  ## pull MCMC sampler out of the resluts:
  theta_draws  <- irt_res$theta
  lambda_draws <- irt_res$lambda
  
  ## Convert MCMC draw to point estimates
  ## via the mean:
  theta_point <- if (length(dim(theta_draws)) == 3) {
    apply(theta_draws, c(1, 2), mean)
  } else theta_draws
  
  
  lambda_point <- if (length(dim(lambda_draws)) == 3) {
    apply(lambda_draws, c(1, 2), mean)
  } else lambda_draws
  
  ## Verify lambda shape to K × d
  if (!all(dim(lambda_point) == dim(true_lambda))) {
    
    if (length(dim(lambda_point)) == 2 &&
        all(dim(lambda_point) == rev(dim(true_lambda)))) {
      lambda_point <- t(lambda_point)
    } else {
      warning(sprintf(
        "lambda dimension mismatch in method %s: %s vs %s",
        method_name,
        paste(dim(lambda_point), collapse="x"),
        paste(dim(true_lambda), collapse="x")
      ))
      lambda_corr[method_name] <- NA_real_
    }
  }
  
  ## MSE using point estimates
  theta_mse[method_name]  <- mse(theta_point,  true_theta,  TRUE)
  lambda_mse[method_name] <- mse(lambda_point, true_lambda, TRUE)
  
  ## Coverage from full draws
  theta_coverage[method_name]  <- coverage(theta_draws,  true_theta)
  lambda_coverage[method_name] <- coverage(lambda_draws, true_lambda)
  
  ## Convergence diagnostics
  theta_gew[method_name]  <- geweke_fcn(theta_draws)
  theta_ess[method_name]  <- ess_fcn(theta_draws)
  theta_rhat[method_name] <- NA_real_  # nchains = 1
  
  ## Correlation only if lengths match
  if (all(dim(lambda_point) == dim(true_lambda))) {
    la <- as.numeric(lambda_point)
    tl <- as.numeric(true_lambda)
    
    if (length(la) == length(tl)) {
      lambda_corr[method_name] <- suppressWarnings(
        cor(la, tl, use = "complete.obs")
      )
    } else {
      lambda_corr[method_name] <- NA_real_
    }
  }
  
  list(
    theta_mse       = theta_mse,
    lambda_mse      = lambda_mse,
    theta_coverage  = theta_coverage,
    lambda_coverage = lambda_coverage,
    theta_gew       = theta_gew,
    theta_ess       = theta_ess,
    theta_rhat      = theta_rhat,
    lambda_corr     = lambda_corr
  )
}
