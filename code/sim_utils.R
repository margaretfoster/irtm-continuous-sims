## sim_utils.R
## Shared helpers for continuous IRT-M simulations
## Update metrics for method
## generates Y
#------------------------------------------------------------------
# Helper: update metrics for a given method
#------------------------------------------------------------------
update_method_metrics <- function(method_name,
                                  theta_avg, lambda_avg,
                                  true_theta, true_lambda,
                                  irt_res,
                                  theta_mse, lambda_mse,
                                  theta_coverage, lambda_coverage,
                                  theta_gew, theta_ess, theta_rhat,
                                  lambda_corr) {
  
  theta_mse[method_name]       <- mse(theta_avg,  true_theta, TRUE)
  lambda_mse[method_name]      <- mse(lambda_avg, true_lambda, TRUE)
  
  theta_coverage[method_name]  <- coverage(theta_avg,  true_theta)
  lambda_coverage[method_name] <- coverage(lambda_avg, true_lambda)
  
  theta_gew[method_name]       <- geweke_fcn(irt_res$theta)
  theta_ess[method_name]       <- ess_fcn(irt_res$theta)
  theta_rhat[method_name]      <- NA_real_
  
  # correlation across all entries
  lambda_corr[method_name] <- suppressWarnings(
    cor(
      as.numeric(lambda_avg),
      as.numeric(true_lambda),
      use = "complete.obs"
    )
  )
  
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

#------------------------------------------------------------------
# Helper: generate Y given mu and dist_type
#------------------------------------------------------------------
generate_Y <- function(mu, dist_type) {
  switch(
    dist_type,
    
    "normal" = mu + rnorm(1, 0, 1),
    
    "heavy_t" = mu + rt(1, df = 3),
    
    "skewed" = {
      # skew-normal via sn package
      sn::rsn(1, xi = mu, omega = 1, alpha = 5)
    },
    
    "bounded" = {
      val <- mu + rnorm(1, 0, 1)
      pmax(pmin(val, 3), -3)
    },
    
    "mixture" = {
      if (runif(1) < 0.9) {
        mu + rnorm(1, 0, 1)
      } else {
        mu + rnorm(1, 0, 5)
      }
    },
    
    "power_law" = {
      # Shifted Pareto-like heavy tail
      xm    <- 1
      alpha <- 2.5
      x     <- xm * (runif(1))^(-1 / alpha)
      mu + x
    },
    
    "log_normal" = {
      sigma <- 1.5
      x     <- rlnorm(1, meanlog = 0, sdlog = sigma)
      mu + x
    },
    
    stop("Unknown distribution type: ", dist_type)
  )
}

#------------------------------------------------------------------
# Helper: convert all_res list to a tidy tibble
#------------------------------------------------------------------
irtm_cont_results_to_df <- function(all_res) {
  out <- lapply(seq_along(all_res), function(i) {
    res <- all_res[[i]]
    r   <- res$results
    
    tibble::tibble(
      sim_id     = i,
      N          = res$N,
      K          = res$K,
      d          = res$d,
      sim        = res$sim,
      dist_type  = res$dist_type,
      lambda_zero_pct = res$lambda_zero_pct,
      n_fake     = res$n_fake,
      nsamp      = res$nsamp,
      nburn      = res$nburn,
      thin       = res$thin,
      nchains    = res$nchains,
      runtime    = res$runtime,
      model      = names(r$theta_mse),
      theta_mse  = as.numeric(r$theta_mse),
      lambda_mse = as.numeric(r$lambda_mse),
      theta_cov  = as.numeric(r$theta_coverage),
      lambda_cov = as.numeric(r$lambda_coverage),
      theta_ess  = as.numeric(r$theta_ess),
      theta_gew  = as.numeric(r$theta_gew),
      theta_rhat = as.numeric(r$theta_rhat),
      lambda_corr = as.numeric(r$lambda_corr)
    )
  })
  
  dplyr::bind_rows(out)
}
