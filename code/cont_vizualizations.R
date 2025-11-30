##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## Visualizations

library(tidyverse)
library(ggplot2)
source("irtm_cont_results_to_df.R")

simpath = "../simulations/"
plotpath = "../simulations/results/figures/"

## Load:

results_df<- readRDS(paste0(simpath, "irtm_cont_parallel_df_benchmarkvariant.rds"))

## recover the parameters:
head(results_df)
print(unique(results_df$model)) # models run
#"benchmark"     "irtM"    "irtMAnchors"  

print(unique(results_df$dist_type))
#normal     heavy_t    skewed     mixture    bounded    log_normal
#power_law 
print(unique(results_df$K)) #10, 50
print(unique(results_df$N)) #50, 100,500

## Plotting consistencies:
models = c("irtM", "irtMAnchors", "benchmark", "benchmark_learnsigma")
model_plotnames = c("IRT-M",
                    "IRT-M + Anchored Thetas", 
                    "Theta Oracle", 
                    "Theta Oracle + Learn Sigma")

cb_palette <- c(
  "IRT-M" = "#E69F00",
  "IRT-M + Anchored Thetas" = "#0072B2", 
  "Theta Oracle" = "#009E73",
  "Theta Oracle + Learn Sigma" = "#CC79A7"
)

### Thetas:
p_theta <- results_df %>%
  ggplot(aes(x = d,
             y = theta_mse,
             color = model,
             group = model)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  facet_wrap(~ dist_type, scales = "free_y") +
  labs(x = "Dimensions (d)",
       y = "Theta MSE",
       color = "Model",
       title = "Theta MSE by Dimension and Distribution") +
  theme_bw(base_size = 14)

p_theta


p_lambda <- results_df %>%
  ggplot(aes(x = d,
             y = lambda_mse,
             color = model,
             group = model)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  facet_wrap(~ dist_type, scales = "free_y") +
  labs(x = "Dimensions (d)",
       y = "Lambda MSE",
       color = "Model",
       title = "Lambda MSE by Dimension and Distribution") +
  theme_bw(base_size = 14)


p_lambda


# Boxplot of theta:

results_df_clean <- results_df %>%
  filter(model %in% models) %>%
  mutate(
    model = factor(
      model,
      levels = models,
      labels = model_plotnames
    ),
    dist_type = factor(dist_type, 
                       levels = c("normal", "heavy_t",
                                  "skewed", "mixture",
                                  "bounded", "log_normal",
                                  "power_law"),
                       labels = c("Normal", "Heavy-t",
                                  "Skewed", "Mixture",
                                  "Bounded", "Log-Normal",
                                  "Power-Law")),
    d = factor(d, levels = sort(unique(d))),
    K = factor(K, levels = sort(unique(K))),
    N = factor(N, levels = sort(unique(N)))
  )


p_box_theta_mse <- ggplot(results_df_clean,
                aes(x = d,
                    y = theta_mse,
                    fill = model)) +
  geom_boxplot(alpha = 0.7,
               outlier.shape = 19,
               outlier.size = 1,
               linewidth = 0.4) +
  facet_grid(N ~ dist_type,
             scales = "fixed",
             switch = "y") +
  ylim(0, 5) + #NOTE: bounded b/c high benchmark outliers
  scale_fill_manual(values = cb_palette,
                    name = "Model",
                    labels = model_plotnames) +
  labs(x = "Dimension (d)",
       y = expression(Theta~MSE),
       title = "Estimation Error Under Distributional Forms",
       subtitle = "IRT-M Continuous Extension",
       caption = "Lines at theta MSE = 1 and Theta MSE = 0.3") +
  theme_bw(base_size = 14) +
  geom_hline(yintercept = 1, #not great 
             linetype = "dashed",
             color = "lightgray", 
             size = 0.5) +
  geom_hline(yintercept = 0.3, ## recovering signal
             linetype = "dashed",
             color = "lightgray",
             size = 0.5) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white",
                                    color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90",
                                    linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              face = "bold", size = 15)
  )

p_box_theta_mse

ggsave(p_box_theta_mse, 
       width = 16,      # wide enough for 7 columns
       height = 9,      # room for title + legend
       dpi = 320,
       units = "in",
       file= paste0(plotpath,
                    "initial_sims_theta_mse.png"))

## Elegant lambda MSE plot:
p_box_lambda_mse <- ggplot(results_df_clean,
                       aes(x = d,
                           y = lambda_mse,
                           fill = model)) +
  geom_boxplot(alpha = 0.7,
               outlier.shape = 19,
               outlier.size = 1,
               linewidth = 0.4) +
  ylim(0, 5) + #NOTE: bounded b/c high benchmark outliers
  facet_grid(N ~ dist_type,
             scales = "fixed",
             switch = "y") +
  scale_fill_manual(values = cb_palette,
                    name = "Model",
                    labels = model_plotnames) +
  labs(x = "Dimension (d)",
       y = expression(lambda~MSE),
       title = "Lambda MSE",
       subtitle = "IRT-M Continuous Extension",
       caption = "Plot bounded at Y = 5") +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white",
                                    color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90",
                                    linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              face = "bold", size = 15)
  )

p_box_lambda_mse

ggsave(p_box_lambda_mse, 
       width = 16,      # wide enough for 7 columns
       height = 9,      # room for title + legend
       dpi = 320,
       units = "in",
       file= paste0(plotpath,
                    "initial_sims_lambda_mse.png"))


## Lambda correlation:

## Elegant lambda Corr plot:
p_box_lambda_corr <- ggplot(results_df_clean,
                       aes(x = d,
                           y = lambda_corr,
                           fill = model)) +
  geom_boxplot(alpha = 0.7,
               outlier.shape = 19,
               outlier.size = 1,
               linewidth = 0.4) +
  ylim(-.1, 1) + # bound Y b/c corr should be [0, 1]
  facet_grid(N ~ dist_type,
             scales = "fixed",
             switch = "y") +
  scale_fill_manual(values = cb_palette,
                    name = "Model",
                    labels = model_plotnames) +
  labs(x = "Dimension (d)",
       y = expression(lambda~Corr),
       title = "Correlation Estimated vs True Lambda Under Distributional Forms",
       subtitle = "IRT-M Continuous Extension",
       caption = "Y bounded [-0.1, 1]") +
  theme_bw(base_size = 14) +
  geom_hline(yintercept = .75, #High cov 
             linetype = "dashed",
             color = "lightgray", 
             size = 0.5) +
  geom_hline(yintercept = 0, ## no covariance
             linetype = "dashed",
             color = "red",
             alpha =.25,
             size = 0.5) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white",
                                    color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90",
                                    linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              face = "bold", size = 15)
  )

p_box_lambda_corr

ggsave(p_box_lambda_corr, 
       width = 16,      # wide enough for 7 columns
       height = 9,      # room for title + legend
       dpi = 320,
       units = "in",
       file= paste0(plotpath,
                    "initial_sims_lambda_corr.png"))

