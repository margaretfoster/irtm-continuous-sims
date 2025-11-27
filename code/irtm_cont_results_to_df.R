irtm_cont_results_to_df <- function(all_res) {
  
  out <- lapply(seq_along(all_res), function(i) {
    res <- all_res[[i]]
    
    tibble::tibble(
      
      # run metadata:
      sim_id     = i,
      N          = res$N,
      K          = res$K,
      d          = res$d,
      dist_type  = res$dist_type,
      model      = names(res$results$theta_mse),
      
      ## run results:
      theta_mse       = as.numeric(res$results$theta_mse),
      lambda_mse      = as.numeric(res$results$lambda_mse),
      theta_coverage  = as.numeric(res$results$theta_coverage),
      lambda_coverage = as.numeric(res$results$lambda_coverage),
      theta_gew       = as.numeric(res$results$theta_gew),
      theta_ess       = as.numeric(res$results$theta_ess),
      theta_rhat      = as.numeric(res$results$theta_rhat),
      lambda_corr     = as.numeric(res$results$lambda_corr)
    )
  })
  
  dplyr::bind_rows(out)
}