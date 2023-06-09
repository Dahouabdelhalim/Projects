##### 'Svalbard reindeer winter diets interannual variation and long-term changes
#####  based on serum d13C and d15N values using source data from winters 2013 and 2019'

# Load all the packages required to load excel files
library(simmr)
library(readxl) 
library(httr) 
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpmisc)
library(viridis)

# Find out what the sheet names are and load in all of them
sheet_names = excel_sheets(path = 'Reindeer_data_schubert.xlsx')
all = lapply(sheet_names,
             read_excel, path = 'Reindeer_data_schubert.xlsx')

##Load the data
# Extract out the different pieces
mix = all[[1]]
# Sources are in different sheets: Source 1 [2], Source 2 [3], Source 3 [4], Source 4 [5]
source = all[[2]]#Source 1
TDF = all[[6]]

str(mix)
str(source)
str(TDF)
names(mix)
names(source)
names(TDF)

## Extract out the different pieces
#Check out the years to be analysed :1995 - 2012 but remember 2003 and 2010 are missing.
table(mix[,'year'])

# Get the data into simmr for all reindeer grouped by year
simmr_reindeer_years_sider = simmr_load(mixtures=as.matrix(mix[,1:2]),
                                        source_names=as.matrix(source[,1]),
                                        source_means=as.matrix(source[,2:3]),
                                        source_sds=as.matrix(source[,4:5]),
                                        correction_means=as.matrix(TDF[,2:3]),
                                        correction_sds=as.matrix(TDF[,4:5]),
                                        group =as.matrix(mix[,3]))
# Plot the iso-space plot
isospace.plot.years.lw<- plot(simmr_reindeer_years_sider,group = 1:16,
                              xlab=expression(paste(delta^13, "C (\\u2030)",sep="")), 
                              ylab=expression(paste(delta^15, "N (\\u2030)",sep="")), 
                              title='Isospace plot of Svalbard Reindeer Data 1995 to 2012',
                              mix_name='')

# Set the priors: late winter

#Late Winter Priors
#prior=simmr_elicit(n_sources = 5,
#                   c(0.05,0.04,0.36,0.29,0.26),
#                 c(0.04,0.03,0.09,0.09,0.07))
#priorlwexample = prior
#save(priorlw, file ="ellicitedlatewinterprior.rda")
load(file ="ellicitedlatewinterprior.rda") # For reproducibility the best ellicited priors were selected and save for future use.

# Run the model
#ITERATIONS = 130 000, BURN = 5 000, THIN = 50, CHAINS = 4

simmr_reindeer_years_out_sider = simmr_mcmc(simmr_reindeer_years_sider, 
                                            prior_control=list(means=priorlw$mean,
                                                               sd=priorlw$sd),
                                            mcmc_control = list(iter = 130000,
                                                                burn = 5000,
                                                                thin = 50,
                                                                n.chain = 4))
# Check convergence
summary(simmr_reindeer_years_out_sider, type = 'diagnostics', group = 1:16) 
#Models Converged

#Check the posterior predicive
post_pred_LW_1995 = posterior_predictive(simmr_reindeer_years_out_sider, group = 1, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_1996 = posterior_predictive(simmr_reindeer_years_out_sider, group = 2, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_1997 = posterior_predictive(simmr_reindeer_years_out_sider, group = 3, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_1998 = posterior_predictive(simmr_reindeer_years_out_sider, group = 4, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_1999 = posterior_predictive(simmr_reindeer_years_out_sider, group = 5, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2000 = posterior_predictive(simmr_reindeer_years_out_sider, group = 6, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2001 = posterior_predictive(simmr_reindeer_years_out_sider, group = 7, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2002 = posterior_predictive(simmr_reindeer_years_out_sider, group = 8, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2004 = posterior_predictive(simmr_reindeer_years_out_sider, group = 9, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2005 = posterior_predictive(simmr_reindeer_years_out_sider, group = 10, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2006 = posterior_predictive(simmr_reindeer_years_out_sider, group = 11, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2007 = posterior_predictive(simmr_reindeer_years_out_sider, group = 12, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2008 = posterior_predictive(simmr_reindeer_years_out_sider, group = 13, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2009 = posterior_predictive(simmr_reindeer_years_out_sider, group = 14, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2011 = posterior_predictive(simmr_reindeer_years_out_sider, group = 15, prob = 0.5, plot_ppc = TRUE)
post_pred_LW_2012 = posterior_predictive(simmr_reindeer_years_out_sider, group = 16, prob = 0.5, plot_ppc = TRUE)

#Proportion outside posterior predictive interval of 50%
PP1995 <- post_pred_LW_1995$prop_outside
PP1996 <- post_pred_LW_1996$prop_outside
PP1997 <- post_pred_LW_1997$prop_outside
PP1998 <- post_pred_LW_1998$prop_outside
PP1999 <- post_pred_LW_1999$prop_outside
PP2000 <- post_pred_LW_2000$prop_outside
PP2001 <- post_pred_LW_2001$prop_outside
PP2002 <- post_pred_LW_2002$prop_outside
PP2004 <- post_pred_LW_2004$prop_outside
PP2005 <- post_pred_LW_2005$prop_outside
PP2006 <- post_pred_LW_2006$prop_outside
PP2007 <- post_pred_LW_2007$prop_outside
PP2008 <- post_pred_LW_2008$prop_outside
PP2009 <- post_pred_LW_2009$prop_outside
PP2011 <- post_pred_LW_2011$prop_outside
PP2012 <- post_pred_LW_2012$prop_outside

##Box plots
plot(simmr_reindeer_years_out_sider,type='boxplot',group=1:16, ggargs = ylim(0,1))

##Density plots
plot(simmr_reindeer_years_out_sider, type = 'density', group = 1,title='Density plot of Svalbard Reindeer Diets late winter 1995')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 2,title='Density plot of Svalbard Reindeer Diets late winter 1996')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 3,title='Density plot of Svalbard Reindeer Diets late winter 1997')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 4,title='Density plot of Svalbard Reindeer Diets late winter 1998')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 5,title='Density plot of Svalbard Reindeer Diets late winter 1999')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 6,title='Density plot of Svalbard Reindeer Diets late winter 2000')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 7,title='Density plot of Svalbard Reindeer Diets late winter 2001')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 8,title='Density plot of Svalbard Reindeer Diets late winter 2002')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 9,title='Density plot of Svalbard Reindeer Diets late winter 2004')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 10,title='Density plot of Svalbard Reindeer Diets late winter 2005')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 11,title='Density plot of Svalbard Reindeer Diets late winter 2006')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 12,title='Density plot of Svalbard Reindeer Diets late winter 2007')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 13,title='Density plot of Svalbard Reindeer Diets late winter 2008')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 14,title='Density plot of Svalbard Reindeer Diets late winter 2009')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 15,title='Density plot of Svalbard Reindeer Diets late winter 2011')+ facet_wrap("~ Source", scales = 'free_y')
plot(simmr_reindeer_years_out_sider, type = 'density', group = 16,title='Density plot of Svalbard Reindeer Diets late winter 2012')+ facet_wrap("~ Source", scales = 'free_y')
#To compare source across the years
table(source[,'Source'])
compare_groups(simmr_reindeer_years_out_sider,source='Graminoids',groups=1:16)
compare_groups(simmr_reindeer_years_out_sider,source='Mosses',groups=1:16)
compare_groups(simmr_reindeer_years_out_sider,source='Salix polaris',groups=1:16)
compare_groups(simmr_reindeer_years_out_sider,source='Forbs',groups=1:16)
compare_groups(simmr_reindeer_years_out_sider,source='Dryas octopetala',groups=1:16)

# Matrix plots
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 1,title='Matrix plot of Svalbard Reindeer Diets late winter 1995')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 2,title='Matrix plot of Svalbard Reindeer Diets late winter 1996')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 3,title='Matrix plot of Svalbard Reindeer Diets late winter 1997')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 4,title='Matrix plot of Svalbard Reindeer Diets late winter 1998')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 5,title='Matrix plot of Svalbard Reindeer Diets late winter 1999')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 6,title='Matrix plot of Svalbard Reindeer Diets late winter 2000')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 7,title='Matrix plot of Svalbard Reindeer Diets late winter 2001')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 8,title='Matrix plot of Svalbard Reindeer Diets late winter 2002')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 9,title='Matrix plot of Svalbard Reindeer Diets late winter 2004')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 10,title='Matrix plot of Svalbard Reindeer Diets late winter 2005')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 11,title='Matrix plot of Svalbard Reindeer Diets late winter 2006')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 12,title='Matrix plot of Svalbard Reindeer Diets late winter 2007')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 13,title='Matrix plot of Svalbard Reindeer Diets late winter 2008')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 14,title='Matrix plot of Svalbard Reindeer Diets late winter 2009')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 15,title='Matrix plot of Svalbard Reindeer Diets late winter 2011')
plot(simmr_reindeer_years_out_sider, type = 'matrix', group = 16,title='Matrix plot of Svalbard Reindeer Diets late winter 2012')

##plot the Prior and posterior distributions
# Prior and posterior distributions
prior_viz(simmr_reindeer_years_out_sider, group = 1, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 2, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 3, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 4, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 5, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 6, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 7, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 8, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 9, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 10, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 11, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 12, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 13, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 14, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 15, plot = TRUE, include_posterior = TRUE, n_sims = 10000)
prior_viz(simmr_reindeer_years_out_sider, group = 16, plot = TRUE, include_posterior = TRUE, n_sims = 10000)

#MODEL COMPARISON
#Lower DIC value better model, can compare between different model runs and groups.
DIC1995 <- simmr_reindeer_years_out_sider$output[[1]]$BUGSoutput$DIC
DIC1996 <-simmr_reindeer_years_out_sider$output[[2]]$BUGSoutput$DIC 
DIC1997 <-simmr_reindeer_years_out_sider$output[[3]]$BUGSoutput$DIC 
DIC1998 <-simmr_reindeer_years_out_sider$output[[4]]$BUGSoutput$DIC 
DIC1999 <-simmr_reindeer_years_out_sider$output[[5]]$BUGSoutput$DIC 
DIC2000 <-simmr_reindeer_years_out_sider$output[[6]]$BUGSoutput$DIC 
DIC2001 <-simmr_reindeer_years_out_sider$output[[7]]$BUGSoutput$DIC 
DIC2002 <-simmr_reindeer_years_out_sider$output[[8]]$BUGSoutput$DIC 
DIC2004 <-simmr_reindeer_years_out_sider$output[[9]]$BUGSoutput$DIC 
DIC2005 <-simmr_reindeer_years_out_sider$output[[10]]$BUGSoutput$DIC 
DIC2006 <-simmr_reindeer_years_out_sider$output[[11]]$BUGSoutput$DIC 
DIC2007 <-simmr_reindeer_years_out_sider$output[[12]]$BUGSoutput$DIC 
DIC2008 <-simmr_reindeer_years_out_sider$output[[13]]$BUGSoutput$DIC 
DIC2009 <-simmr_reindeer_years_out_sider$output[[14]]$BUGSoutput$DIC 
DIC2011 <-simmr_reindeer_years_out_sider$output[[15]]$BUGSoutput$DIC 
DIC2012 <-simmr_reindeer_years_out_sider$output[[16]]$BUGSoutput$DIC 

#Save the posterior distribution
save(simmr_reindeer_years_out_sider , file ="simmr_years_source1.rda")


