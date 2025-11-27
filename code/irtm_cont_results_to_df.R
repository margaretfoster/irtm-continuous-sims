irtm_cont_results_to_df <- function(all_res) {
  
  out <- lapply(seq_along(all_res), function(i) {
    res <- all_res[[i]]
    
    tibble::tibble(
      sim_id     = i,
      N          = res$N,
      K          = res$K,
      d          = res$d,
      dist_type  = res$dist_type,
      model      = rep(names(res$results$theta_mse), 1),
      theta_mse  = as.numeric(res$results$theta_mse),
      lambda_mse = as.numeric(res$results$lambda_mse),
      theta_cov  = as.numeric(res$results$theta_cov),
      lambda_cov = as.numeric(res$results$lambda_cov),
      theta_ess  = as.numeric(res$results$theta_ess),
      theta_gew  = as.numeric(res$results$theta_gew),
      theta_rhat = as.numeric(res$results$theta_rhat)
    )
  })
  
  dplyr::bind_rows(out)
}