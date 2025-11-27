## run_model.R
## Run a continuous IRT-M under various identification strategies

run_iter_cont <- function(
    Y, Yall,
    true_theta, 
    true_theta_all, 
    true_lambda,
    M, x, z, ind, d,
    nburn, nsamp, thin, nchains,
    progress = TRUE,
    methods
) {
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
    message("Running benchmark...")
    irt_res <- M_constrained_irt_continuous(
      Y, d, M = NULL,
      theta_fix = z, which_fix = ind,
      nburn = nburn, nsamp = nsamp, thin = thin,
      learn_Sigma = FALSE, learn_Omega = FALSE,
      display_progress = progress
    )
    
    updated <- update_method_metrics(
      method_name      = "benchmark",
      theta_avg        = irt_res$theta,
      lambda_avg       = irt_res$lambda,
      true_theta       = true_theta,
      true_lambda      = true_lambda,
      irt_res          = irt_res,
      theta_mse        = theta_mse,
      lambda_mse       = lambda_mse,
      theta_coverage   = theta_coverage,
      lambda_coverage  = lambda_coverage,
      theta_gew        = theta_gew,
      theta_ess        = theta_ess,
      theta_rhat       = theta_rhat,
      lambda_corr      = lambda_corr
    )
    
    theta_mse       <- updated$theta_mse
    lambda_mse      <- updated$lambda_mse
    theta_coverage  <- updated$theta_coverage
    lambda_coverage <- updated$lambda_coverage
    theta_gew       <- updated$theta_gew
    theta_ess       <- updated$theta_ess
    theta_rhat      <- updated$theta_rhat
    lambda_corr     <- updated$lambda_corr
    
    rm(irt_res)
  }
  
  #--------------------------------------#
  # Unconstrained
  #--------------------------------------#
  if ("unconstrained" %in% methods) {
    message("Running unconstrained...")
    irt_res <- M_constrained_irt_continuous(
      Y, d, M = NULL,
      theta_fix = NULL, which_fix = NULL,
      nburn = nburn, nsamp = nsamp, thin = thin,
      learn_Sigma = FALSE, learn_Omega = FALSE,
      display_progress = progress
    )
    
    updated <- update_method_metrics(
      method_name      = "unconstrained",
      theta_avg        = irt_res$theta,
      lambda_avg       = irt_res$lambda,
      true_theta       = true_theta,
      true_lambda      = true_lambda,
      irt_res          = irt_res,
      theta_mse        = theta_mse,
      lambda_mse       = lambda_mse,
      theta_coverage   = theta_coverage,
      lambda_coverage  = lambda_coverage,
      theta_gew        = theta_gew,
      theta_ess        = theta_ess,
      theta_rhat       = theta_rhat,
      lambda_corr      = lambda_corr
    )
    
    theta_mse       <- updated$theta_mse
    lambda_mse      <- updated$lambda_mse
    theta_coverage  <- updated$theta_coverage
    lambda_coverage <- updated$lambda_coverage
    theta_gew       <- updated$theta_gew
    theta_ess       <- updated$theta_ess
    theta_rhat      <- updated$theta_rhat
    lambda_corr     <- updated$lambda_corr
    
    rm(irt_res)
  }
  
  #--------------------------------------#
  # irtM (full M constraint)
  #--------------------------------------#
  if ("irtM" %in% methods) {
    message("Running M-constrained IRT...")
    irt_res <- M_constrained_irt_continuous(
      Yall, d, M = M,
      theta_fix = NULL, which_fix = NULL,
      nburn = nburn, nsamp = nsamp, thin = thin,
      learn_Sigma = TRUE, display_progress = progress
    )
    
    updated <- update_method_metrics(
      method_name      = "irtM",
      theta_avg        = irt_res$theta,
      lambda_avg       = irt_res$lambda,
      true_theta       = true_theta_all,
      true_lambda      = true_lambda,
      irt_res          = irt_res,
      theta_mse        = theta_mse,
      lambda_mse       = lambda_mse,
      theta_coverage   = theta_coverage,
      lambda_coverage  = lambda_coverage,
      theta_gew        = theta_gew,
      theta_ess        = theta_ess,
      theta_rhat       = theta_rhat,
      lambda_corr      = lambda_corr
    )
    
    theta_mse       <- updated$theta_mse
    lambda_mse      <- updated$lambda_mse
    theta_coverage  <- updated$theta_coverage
    lambda_coverage <- updated$lambda_coverage
    theta_gew       <- updated$theta_gew
    theta_ess       <- updated$theta_ess
    theta_rhat      <- updated$theta_rhat
    lambda_corr     <- updated$lambda_corr
    
    rm(irt_res)
  }
  
  #--------------------------------------#
  # irtMNoVar
  #--------------------------------------#
  if ("irtMNoVar" %in% methods) {
    message("Running M-constrained IRT without variance prior...")
    irt_res <- M_constrained_irt_continuous(
      Yall, d, M = M,
      theta_fix = x, which_fix = NULL,
      nburn = nburn, nsamp = nsamp, thin = thin,
      learn_Sigma = FALSE, learn_Omega = FALSE,
      display_progress = progress
    )
    
    updated <- update_method_metrics(
      method_name      = "irtMNoVar",
      theta_avg        = irt_res$theta,
      lambda_avg       = irt_res$lambda,
      true_theta       = true_theta_all,
      true_lambda      = true_lambda,
      irt_res          = irt_res,
      theta_mse        = theta_mse,
      lambda_mse       = lambda_mse,
      theta_coverage   = theta_coverage,
      lambda_coverage  = lambda_coverage,
      theta_gew        = theta_gew,
      theta_ess        = theta_ess,
      theta_rhat       = theta_rhat,
      lambda_corr      = lambda_corr
    )
    
    theta_mse       <- updated$theta_mse
    lambda_mse      <- updated$lambda_mse
    theta_coverage  <- updated$theta_coverage
    lambda_coverage <- updated$lambda_coverage
    theta_gew       <- updated$theta_gew
    theta_ess       <- updated$theta_ess
    theta_rhat      <- updated$theta_rhat
    lambda_corr     <- updated$lambda_corr
    
    rm(irt_res)
  }
  
  #--------------------------------------#
  # Anchors only
  #--------------------------------------#
  if ("irtAnchors" %in% methods) {
    message("Running IRT with anchors...")
    irt_res <- M_constrained_irt_continuous(
      Yall, d, M = NULL,
      theta_fix = x, which_fix = ind,
      nburn = nburn, nsamp = nsamp, thin = thin,
      learn_Sigma = TRUE,
      display_progress = progress
    )
    
    updated <- update_method_metrics(
      method_name      = "irtAnchors",
      theta_avg        = irt_res$theta,
      lambda_avg       = irt_res$lambda,
      true_theta       = true_theta_all,
      true_lambda      = true_lambda,
      irt_res          = irt_res,
      theta_mse        = theta_mse,
      lambda_mse       = lambda_mse,
      theta_coverage   = theta_coverage,
      lambda_coverage  = lambda_coverage,
      theta_gew        = theta_gew,
      theta_ess        = theta_ess,
      theta_rhat       = theta_rhat,
      lambda_corr      = lambda_corr
    )
    
    theta_mse       <- updated$theta_mse
    lambda_mse      <- updated$lambda_mse
    theta_coverage  <- updated$theta_coverage
    lambda_coverage <- updated$lambda_coverage
    theta_gew       <- updated$theta_gew
    theta_ess       <- updated$theta_ess
    theta_rhat      <- updated$theta_rhat
    lambda_corr     <- updated$lambda_corr
    
    rm(irt_res)
  }
  
  #--------------------------------------#
  # IRTM with anchors + M
  #--------------------------------------#
  if ("irtMAnchors" %in% methods) {
    message("Running IRT with M and anchors...")
    irt_res <- M_constrained_irt_continuous(
      Yall, d, M = M,
      theta_fix = x, which_fix = ind,
      nburn = nburn, nsamp = nsamp, thin = thin,
      learn_Sigma = TRUE,
      display_progress = progress
    )
    
    updated <- update_method_metrics(
      method_name      = "irtMAnchors",
      theta_avg        = irt_res$theta,
      lambda_avg       = irt_res$lambda,
      true_theta       = true_theta_all,
      true_lambda      = true_lambda,
      irt_res          = irt_res,
      theta_mse        = theta_mse,
      lambda_mse       = lambda_mse,
      theta_coverage   = theta_coverage,
      lambda_coverage  = lambda_coverage,
      theta_gew        = theta_gew,
      theta_ess        = theta_ess,
      theta_rhat       = theta_rhat,
      lambda_corr      = lambda_corr
    )
    
    theta_mse       <- updated$theta_mse
    lambda_mse      <- updated$lambda_mse
    theta_coverage  <- updated$theta_coverage
    lambda_coverage <- updated$lambda_coverage
    theta_gew       <- updated$theta_gew
    theta_ess       <- updated$theta_ess
    theta_rhat      <- updated$theta_rhat
    lambda_corr     <- updated$lambda_corr
    
    rm(irt_res)
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
