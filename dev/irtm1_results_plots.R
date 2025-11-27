## Results:
## takes dataframe of results
## produces graphs

library(dplyr)
library(tidyr)
library(ggplot2)

source('code/irtm_only_results_to_df.R') ## script to format list of results as single df
source('code/na_indicies.R') # helper to find NA values in model comparisons

## Load data:
## params:
load("simulations/irtm_small_sim_params.rds")
## all_results:
load("simulations/irtm_only_N100K10.rds")
## model_times:
load("simulations/irtm_only_N100K10_time.rds")

## %%%%%%%%%%%%% Data Prep
## Standardize list names for the time list:

all_results <- lapply(all_results, function(x) {
  if ("irtM" %in% names(x)) {
    names(x)[names(x) == "irtM"] <- "irtm"
  }
  return(x)
})


all_results <- lapply(all_results, function(x) {
  if ("CFA" %in% names(x)) {
    names(x)[names(x) == "CFA"] <- "bcfa"
  }
  return(x)
})

## get index of NA values in the sim results
## NA out the time runs for those models
## (to not distort estimated time)

na_ind <-na_indicies(all_results)

## update time to NA out times for failed models:
updated_time <- na_times(model_times, na_ind)

##%%%%%%%%%%%%%%%%%
## Process the results from a list to a df:

## Data for Fig 1: MSE by model:
sim_results <- irtm_results_to_df(all_results,
                                  small_sim, 
                                  mode = 1) # model stats

sim_results$statistic <- "mse" # for tracking

## Data for Fig 2: Runtime by model
time_results <- irtm_results_to_df(updated_time,
                                   small_sim,
                                   mode = 1) # model stats

time_results$statistic <- "runtime"

## Data for Fig 3:  count NA by models:
na_count_by_group <- irtm_results_to_df(all_results,
                                         small_sim, 
                                         mode = 2) # NA count

##%%%%%%%%%%%%%%%%%%%%
## Plots:

p_results <- ggplot(sim_results, aes(x= pass_d, 
                       y=avg,
                      ymin=min, 
                       ymax=max,
                      pch=as.factor(model), 
                      linetype=as.factor(model))) + 
  geom_point(size=2.5) +
  geom_line() + 
  geom_ribbon(alpha=0.2, linewidth =0) + 
  labs(linetype='Model', pch='Model') + 
  xlab('Dimensions') + 
  ylab('MSE') + 
  scale_x_continuous(breaks=c(2,3, 5, 8))+
  theme_bw() +
  theme(legend.position="bottom", 
        legend.background = element_rect(color='black'), 
        text=element_text(size=20, family="Times")) 

p_results


ggsave(filename='simulations/results/figures/irt-only.png', 
       plot=p_results,
       dpi = 300,  width = 12, height=6)


## Time comparison

p_time <- ggplot(time_results, aes(x= pass_d, 
                                     y=avg,
                                     ymin=min, 
                                     ymax=max,
                                     pch=as.factor(model), 
                                     linetype=as.factor(model))) + 
  geom_point(size=2.5) +
  geom_line() + 
  geom_ribbon(alpha=0.2, linewidth =0) + 
  labs(linetype='Model', pch='Model') + 
  xlab('Dimensions') + 
  ylab('Runtime (Seconds)') + 
  scale_x_continuous(breaks=c(2,3, 5, 8))+
  theme_bw() +
  theme(legend.position="bottom", 
        legend.background = element_rect(color='black'), 
        text=element_text(size=20, family="Times")) 

p_time


ggsave(filename='simulations/results/figures/runtime_irtm_only.png', 
       plot=p_time,
       dpi = 300,  width = 12, height=6)

## Percent failed models:
## first count NA results by model and d specification

## Summarize the number of failed models
## per dimension-model pairing:


## Boxplot for runtime, to get visibility into the 
## time spike at D=2

## For Fig 4: boxplot of run times

runtime_all <-irtm_results_to_df(updated_time,
                                        small_sim, 
                                        mode = 3) # runtime

runtime_all$statistic <- "runtime" ## to keep track


runtime_all <- runtime_all[, c('irtm','pass_d')] %>% 
  pivot_longer(
  cols = c('irtm'),  # Select the columns that start with avg, min, max
  names_to = c("model"),  # Create columns for model and statistic
  values_to = "value"  # Name of the new column for the MSE values
) 

## Individual plots for runtime box plots; b/c the 
## scales
grtbox <- ggplot(runtime_all, 
       aes(x=as.factor(pass_d),
           y=value)) + 
  geom_boxplot(fill="lightgray",
               alpha=0.2, 
               outlier.colour="darkgray",
               outlier.fill="darkgray",
               outlier.size=3) + 
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  xlab("Dimensions") + 
  ylab("Runtime (seconds) ") + 
  theme_bw() +
  facet_wrap(~model, ncol=1, scales="free")


grtbox

ggsave(filename='simulations/results/figures/irtm_only_runtime_boxplot.png',
       plot=grtbox,
       dpi = 300, 
       width = 12, 
       height=6)

