## Two helper functions: 

## First--  Fast and dirty for-loop
## Takes a list of simulation results
## returns a dataframe with index locations of NA results

## second: for-loop that goes through the index of NA results 
## and goes into list 2 & converts entries into NA. 

## Designed as two standalone functions so that the 
## NA index list can be audited before dropping data.

na_indicies <- function(results_list){
  na_index1 <- c()
  na_index2 <- c()
  
  for(i in 1:length(results_list)){
    #for(i in 1:length(results_list)){
    if(is.null(results_list[[i]])){
      next ##
    }
    for(j in c(1:3)){
      if(is.na(results_list[[i]][j])==TRUE){
        na_index1 <- c(na_index1, i)
        na_index2 <- c(na_index2, j)
      }
    }
  }
  
  na_places <- data.frame(cbind(na_index1, na_index2))
  
  return(na_places)
}


## Second function to NA out time for models that failed:
## not hugely efficient, but N is small:
na_times <- function(time_list, na_index_df){
  if(!class(na_ind) == "data.frame" |!class(time_list) == "list" ){
    print("check input classes, no data changed")
    return(time_list)
  }
  ## want to go through nrows in na_index_df
  ## for each row, go into the time list
  ## for failed models and NA out the time
  if(nrow(na_ind)==0){
    print("NA index dataframe is empty, returning original list")
    return(time_list)
  } #close if no rows
  else{ 
    print("making changes")
    for(r in 1:nrow(na_ind)){
      print(r)
      i = na_ind[r,]$na_index1 ## first level
      j = na_ind[r,]$na_index2 ## second level
      time_list[[i]][j] <-NA 
      } #close for
    } # close else
  
    return(time_list)
} #close function

