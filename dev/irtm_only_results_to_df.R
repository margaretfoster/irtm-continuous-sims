
## load a results dataframe
## load the param scaffold

irtm_results_to_df <- function(results_list, params_df, mode= 1){

  ## Add a catch that mode needs to be one of modelNA or modelResults:
  if(!(mode %in% c(1, 2, 3))){
    print("The mode parameter should be 1 (to format results), 2 (to format time),
          or 3 (data from all runs)!")
  }
  
  
  results_irtm <- data.frame()
  
  ## seq_along() creates an integer sequence of the length of the input:
  for (i in seq_along(results_list)){ ## all results is df of theta mse
    print(i)
    
    if(is.null(results_list[[i]])){
      next
    }
    
    round_res <- results_list[[i]]
    round_params <- params_df[i,]
    
    # Create a data frame for the current list element
    round_df <- data.frame(
      iteration = i, ## track which iteration
      irtm = round_res$irtm,
      #bsem = round_res$bsem,
      #bcfa = round_res$bcfa,
      pass_N = round_params$Var1,
      pass_K = round_params$Var2,
      pass_d = round_params$Var3
    )
    
    # Bind the current data frame to the result data frame
    results_irtm <- rbind(results_irtm, round_df)
  }
  
  if(mode==3){
    return(results_irtm) ## if I want all run results
  }
  
  ## Mode: NA analysis:
  if(mode == 2){
    na_count_by_group <- results_irtm %>%
      # Count NAs for each combination of group_column and another column
      group_by(pass_d) %>%
      summarise(
        irtm = sum(is.na(irtm))/50,  # Count NAs in value_column1
        #bcfa = sum(is.na(bcfa))/50,   # Count NAs in value_column2
        #bsem = sum(is.na(bsem))/50   # Count NAs in value_column2
      ) %>%
      pivot_longer(!pass_d,
                   names_to = "model",
                   values_to = "pct_missing")
    
    return(na_count_by_group)
  } ## close mode == 2
  ## Mode: plot values:
  if(mode == 1){
    ## Process the results into the format that
    ## we want for plotting
    
    ## Group by dimensions (data is 50 iterations for each dimension)
    ## take min, max, mean for each model for each of the varying params:
    agg_results_to_plot <- results_irtm %>%
      group_by(pass_d) %>% ## in this group, dimension is the only thing that varies
      summarise(
        ## irtm
        avg_irtm  = mean(irtm, na.rm = TRUE),
        min_irtm  = min(irtm, na.rm = TRUE),
        max_irtm  = max(irtm, na.rm = TRUE)
      )
    
    ## Pivot longer to make plotting easier:
    ## 
    ggdat <- agg_results_to_plot %>%
      pivot_longer(
        cols = starts_with("avg_") | starts_with("min_") | starts_with("max_"),  # Select the columns that start with avg, min, max
        names_to = c("statistic", "model"),  # Create columns for model and statistic
        names_pattern = "(.*)_(.*)",  # Split column names into model and statistic
        values_to = "value"  # Name of the new column for the MSE values
      ) %>%
      pivot_wider(
        names_from = statistic,  # new column names
        values_from = value,  # new column values
      )
    
    return(ggdat)
  } #close mode == 1 
}

## save template:
#fname = paste0("simulations/irtm_bcfa_bsem_comp", Sys.Date(), ".rds")
#save(results_irtm, file=fname) 
