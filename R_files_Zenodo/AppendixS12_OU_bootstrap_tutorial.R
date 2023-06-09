## ----setup, message = F, include = FALSE} ------
# You can ignore this chunk of code if simply running through the R code for the tutorial (only relevant for compiling the R Markdown document)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache = FALSE) # Change to 'FALSE' when debugging
library(formatR)
knitr::opts_chunk$set(tidy.opts = list(width.cutoff = 65), tidy = TRUE) 


## ----load_data, message = F} ------
library(OUwie) # Our principal package for OU analysis. Version 2.6. 
library(ouch) # An alternative package for OU analysis. Version 2.17.
library(pmc) # For simple parametric bootstrapping with ouch. Version 1.0.4.
library(parallel) # For parallel analyses. Version 4.1.0.
library(phytools) # Used here for stochastic mapping. Version 0.7-70.
library(tidyverse) # Suite of packages for data manipulation. Key included packages for
# this tutorial are dplyr, tidyr, ggplot2, and purrr. Version 1.3.0. 
source('AppendixS6_tutorial_functions.R') # Custom functions we wrote for this tutorial
load('AppendixS11_OU_bootstrap_data.RData') # R workspace with objects from previous tutorial

## ----95CI_OUwie, warn = F, message = F, results = 'hide'} ------
mod <- OU2_OUwie # Changing the name simplifies later code.
set.seed(1978) # To ensure your results are identical to ours 
nsims <- 100
CTmin_boot_OUwie <- OUwie.boot(phy = tree_ou2, data = dat_OU2, model =  mod$model, nboot = nsims, alpha = mod$solution[1,], sigma.sq = mod$solution[2,], theta = mod$theta[,1], theta0 = mod$theta[2,1], scaleHeight = mod$scaleHeight, root.station = mod$root.station, get.root.theta = mod$get.root.theta, mserr = 'known', algorithm = mod$algorithm, warn = FALSE)


## ----extract_95CI_OUwie} ------
# Calculate 95% CI
apply(CTmin_boot_OUwie, 2, quantile, probs=c(0.025,0.975)) 


## ----dependence_plot, tidy = F, echo = T, fig.height = 3, fig.width = 3.5, fig.cap = 'Plot of $\\\\alpha$ and $\\\\sigma^2$ estimated across parametric bootstrap replicates. Note that we only use this plot for heuristic reasons; the simulated data are all different, so technically this is **not** showing how similar likelihoods for a given dataset can be obtained at different values of $\\\\alpha$ and $\\\\sigma^2$. However, it should be very similar to what we want, given that all the simulations were based on the same parameter values. Both axes are in original units on a logged scale, as their variation seems to fit a lognormal distribution The line has an intercept of 0 and slope of 1.0 to show that as one parameter goes up, the other is estimated similarly higher. Also note that some of our replicates hit the upper bound of 100 in *OUwie* for these two parameters, so this upper bound may be worth increasing in future analyses.'} ------
class(CTmin_boot_OUwie) <- 'matrix'
CTmin_boot_OUwie %>% 
  as_tibble() %>%
  ggplot(aes(alpha_temperate, sigma.sq_temperate)) +
    geom_point(size = 3, alpha = 0.7) + 
    geom_abline(aes(slope = 1, intercept = 0)) + 
    scale_x_log10(name = quote(alpha)) + 
    scale_y_log10(name = quote(sigma^2)) + 
    theme_linedraw() + 
    coord_fixed()


## ----95CI_ouch, warn = F, message = F} ------
CTmin_boot_ouch <- bootstrap(OU2_ouch, nboot = nsims, seed = 1978)
apply(CTmin_boot_ouch, 2, quantile, probs=c(0.025,0.975)) 


## ----OUwie_bootstrap, message = F, warn = F} ------
nsims <- 100
OUwie_args <- list(scaleHeight = T, root.station = TRUE, algorithm = 'invert', diagn = F, quiet = T, warn = F)
# (a) Compare BM (second-most supported model) to OU2 (most supported model)
BMvOU2 <- paraboot_OUwie(trees = list(tree_ou2, tree_ou2), dat_list = list(dat_OU2, dat_OU2), mods = c('BM1', 'OUM'), nsim = nsims, seed = 1978, OUwie_args = OUwie_args, ncores = NULL)
# (b) Compare OU2 (most supported model) to OU3 (more complex model)
ou2v3 <- paraboot_OUwie(trees = list(tree_ou2, tree_ou3), dat_list = list(dat_OU2, dat_OU3), mods = c('OUM', 'OUM'), nsim = nsims, seed = 1978, OUwie_args = OUwie_args, ncores = NULL)


## ----plot_OUwie_boot, tidy = F, fig.width = 6.5, fig.height = 3, fig.cap = 'Results of parametric bootstrapping with *OUwie*. "null" distributions indicate likelihood ratios calculated between the two models when data were simulated under the simpler model. "test" distributions indicate the same, but when data were simulated under the more complex model. The darkest grey color indicates areas of overlap. Dashed vertical lines indicate the critical value for a standard hypothesis test (i.e. critical level of 0.05). Vertical solid red lines indicate the empirical likelihood-ratio test statistic.'} ------
results <- bind_rows(data.frame(comparison = "BMvOU2", 
                                null = BMvOU2$LR_sim[,1], 
                                test = BMvOU2$LR_sim[,2], 
                                lr = BMvOU2$LR_emp, 
                                sig = BMvOU2$LR_sig), 
                     data.frame(comparison = "ou2v3", 
                                null = ou2v3$LR_sim[,1], 
                                test = ou2v3$LR_sim[,2], 
                                lr = ou2v3$LR_emp, 
                                sig = ou2v3$LR_sig)) %>%
      pivot_longer(c(null, test), names_to = 'test_type', values_to = 'LR_ratio') 
# Using the flexibility of ggplot2 you can customize the plot to visualize your observed
# test comparison ('results$lr', red lines) to what you would expect for a comparison when
# the data are simulated under the two models
    ggplot(results, aes(LR_ratio, fill = test_type)) + 
      geom_density(alpha = 0.75) +
      geom_vline(aes(xintercept = lr), color = 'red', linetype = 'solid') +
      geom_vline(aes(xintercept = sig), color = 'black', linetype = 'dashed') +
      scale_fill_manual(values = c('grey60', 'black')) + 
      scale_y_continuous(name = 'Density', expand = c(0,0)) +
      scale_x_continuous(name = 'Likelihood ratio', expand = c(0,0)) +
      theme_classic() +
      theme(
        axis.text = element_text(color = 'black'),
        strip.text = element_text(size = 14, hjust = 0.5, vjust = 2, face = 'bold'),
        strip.background = element_blank()
      ) +
      facet_wrap(~comparison, nrow = 1, scales = 'free')    