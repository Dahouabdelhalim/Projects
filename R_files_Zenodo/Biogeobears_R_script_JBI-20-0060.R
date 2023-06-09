library(optimx) 
library(FD)
library(snow)
library(parallel)
library(BioGeoBEARS)
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)


#######################################################
# WORKING DIRECTORY
#######################################################
wd = np("~/Dropbox/Research/Postdoc/2-Dictyotales/1-Dictyotales_biogeography/Analyses/Ancestral range reconstr")
setwd(wd)
getwd()

#######################################################
#TREE FILE AND GEOGRAPHY FILE
#######################################################

# Phylogeny file
trfn = np("Time tree JBI-20-0060.tre")
# Raw Newick file:
moref(trfn)
# Phylogeny:
tr = read.tree(trfn)
tr
plot(tr, cex=0.3)
title("Dictyotales phylogeny")
axisPhylo() # plots timescale

#######################################################
# Geography file
geogfn = ("Longitudinal_regions.txt") # where A: Atlantic, B: Central Indo-Pacific, C: Eastern Indo-Pacific, D: Western Indo-Pacific, E: Eastern Pacific
# Raw geography text file:
moref(geogfn)

# Geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Maximum number of areas any species may occupy
max_range_size = 3

####################################################

# Empty tables to hold the statistical results
restable = NULL
teststable = NULL

#######################################################
# DEC AND DEC+J ANALYSIS
#######################################################

#######################################################
# Run DEC
#######################################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$speedup=TRUE        # seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run

# Set up a time-stratified analysis 
# (un-comment to use; see example files in extdata_dir, 
#  and BioGeoBEARS google group posts for further hints)
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing if desired
BioGeoBEARS_run_object$num_cores_to_use=1
# (use more cores to speed it up; this requires
# library(parallel), which is default on Macs in
# R 3.0+, but apparently still has to be typed
# on Windows machines. Note: apparently parallel
# works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Read in the input files. Read the error messages if you get them!
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata (stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "Dictyotales_DEC_unconstrained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

#######################################################
# Run DEC+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE        # seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run

# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table

# Set up DEC+J model
# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Dictyotales_DEC+J_unconstrained_v1.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}

#######################################################
# PDF plots
#######################################################
pdffn = "Dictyotales_DEC_vs_DEC+J_unconstrained_v1.pdf"
pdf(pdffn, width=17, height=22)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Dictyotales unconstrained"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.3, statecex=0.3, splitcex=0.3, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.2, statecex=0.3, splitcex=0.3, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on Dictyotales unconstrained"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.3, statecex=0.3, splitcex=0.3, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.3, statecex=0.3, splitcex=0.3, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Statistics
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
if (BioGeoBEARS_run_object$use_optimx == TRUE)
{
  # Using optimx() results
  if (packageVersion("optimx") < 2013)
  {
    # optimx 2012
    LnL_2 = as.numeric(resDEC$optim_result$fvalues)
    LnL_1 = as.numeric(resDECj$optim_result$fvalues)
  } else {
    # optimx 2013
    LnL_2 = as.numeric(resDEC$optim_result$value)
    LnL_1 = as.numeric(resDECj$optim_result$value)
  } # end optimx 2012 vs. 2013
} else {
  # Using optim() results
  LnL_2 = as.numeric(resDEC$optim_result$value)
  LnL_1 = as.numeric(resDECj$optim_result$value)
} # end optim vs. optimx

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

res2    # DEC, null model for Likelihood Ratio Test (LRT)
res1    # DEC+J, alternative model for Likelihood Ratio Test (LRT)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
# confer the same likelihood on the data. See: Brian O'Meara's webpage:
# http://www.brianomeara.info/tutorials/aic
# ...for an intro to LRT, AIC, and AICc

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# DIVALIKE AND DIVALIKE+J ANALYSIS
#######################################################

# NOTE: The BioGeoBEARS "DIVA" model is not identical with 
# Ronquist (1997)'s parsimony DIVA. It is a likelihood
# interpretation of DIVA, constructed by modelling DIVA's
# processes the way DEC does, but only allowing the 
# processes DIVA allows (widespread vicariance: yes; subset
# sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
#
# I thus now call the model "DIVALIKE", and you should also. ;-)
# 

#######################################################
# Run DIVALIKE
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE        # seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run

# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Dictyotales_DIVALIKE_unconstrained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
}

#######################################################
# Run DIVALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE        # seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run

# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table

# Set up DIVALIKE+J model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Dictyotales_DIVALIKE+J_unconstrained_v1.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKEj = res
}

pdffn = "Dictyotales_DIVALIKE_vs_DIVALIKE+J_unconstrained_v1.pdf"
pdf(pdffn, width=8.5, height=11)

#######################################################
# Plot ancestral states - DIVALIKE
#######################################################
analysis_titletxt ="BioGeoBEARS DIVALIKE on Dictyotales unconstrained"

# Setup
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - DIVALIKE+J
#######################################################
analysis_titletxt ="BioGeoBEARS DIVALIKE+J on Dictyotales unconstrained"

# Setup
results_object = resDIVALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Stats
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
if (BioGeoBEARS_run_object$use_optimx == TRUE)
{
  # Using optimx() results
  if (packageVersion("optimx") < 2013)
  {
    # optimx 2012
    LnL_2 = as.numeric(resDIVALIKE$optim_result$fvalues)
    LnL_1 = as.numeric(resDIVALIKEj$optim_result$fvalues)
  } else {
    # optimx 2013
    LnL_2 = as.numeric(resDIVALIKE$optim_result$value)
    LnL_1 = as.numeric(resDIVALIKEj$optim_result$value)
  } # end optimx 2012 vs. 2013
} else {
  # Using optim() results
  LnL_2 = as.numeric(resDIVALIKE$optim_result$value)
  LnL_1 = as.numeric(resDIVALIKEj$optim_result$value)
} # end optim vs. optimx

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

res2    # DIVALIKE basic results
res1    # DIVALIKE+J basic results

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# BAYAREALIKE AND BAYAREALIKE+J ANALYSIS
#######################################################

# NOTE: As with DIVA, the BioGeoBEARS BayArea-like model is 
# not identical with the full Bayesian model implemented 
# in the "BayArea" program of Landis et al. (2013). 
#
# Instead, this is a simplified likelihood interpretation
# of the model.  Basically, "d" and "e" work like they 
# do in the DEC model of Lagrange (and BioGeoBEARS), and 
# then BayArea's cladogenesis assumption (which is that 
# nothing in particular happens at cladogenesis) is 
# replicated by BioGeoBEARS.
#
# This leaves out 3 important things that are in BayArea:
# 1. Distance dependence (you can add this with a distances 
#    matrix in BioGeoBEARS, however)
# 2. A correction for disallowing "e" events that drive
#    a species extinct (a null geographic range)
# 3. The neat Bayesian sampling of histories, which allows
#    analyses on large numbers of areas.
#
# The main purpose of having a "BAYAREALIKE" model is 
# to test the importance of the cladogenesis model on 
# particular datasets. Does it help or hurt the data 
# likelihood if there is no special cladogenesis process?
# 
# I thus now call the model "BAYAREALIKE", and you should also. ;-)
# 

#######################################################
# Run BAYAREALIKE
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE        # seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run

# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Dictyotales_BAYAREALIKE_unconstrained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
}

#######################################################
# Run BAYAREALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size

# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use=1
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table

# Set up BAYAREALIKE+J model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Dictyotales_BAYAREALIKE+J_unconstrained_v1.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resBAYAREALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKEj = res
}

pdffn = "Dictyotales_BAYAREALIKE_vs_BAYAREALIKE+J_unconstrained_v1.pdf"
pdf(pdffn, width=8.5, height=11)

#######################################################
# Plot ancestral states - BAYAREALIKE
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Dictyotales unconstrained"

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - BAYAREALIKE+J
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J on Dictyotales unconstrained"

# Setup
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Stats
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
if (BioGeoBEARS_run_object$use_optimx == TRUE)
{
  # Using optimx() results
  if (packageVersion("optimx") < 2013)
  {
    # optimx 2012
    LnL_2 = as.numeric(resBAYAREALIKE$optim_result$fvalues)
    LnL_1 = as.numeric(resBAYAREALIKEj$optim_result$fvalues)
  } else {
    # optimx 2013
    LnL_2 = as.numeric(resBAYAREALIKE$optim_result$value)
    LnL_1 = as.numeric(resBAYAREALIKEj$optim_result$value)
  } # end optimx 2012 vs. 2013
} else {
  # Using optim() results
  LnL_2 = as.numeric(resBAYAREALIKE$optim_result$value)
  LnL_1 = as.numeric(resBAYAREALIKEj$optim_result$value)
} # end optim vs. optimx

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

res2    # BAYAREALIKE results
res1    # BAYAREALIKE+J results

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#########################################################################
# RESULTS: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
#########################################################################
teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")

# Look at the results!!
restable
teststable

#######################################################
# The p-value of the LRT (Likelihood Ratio Test) tells you whether or not you can reject the
# null hypothesis that DEC and DEC+J confer equal likelihoods on the data
#
# AIC and AIC model weights are also shown, giving a sense of the relative probability of the two models.
#
# (One could easily do model weights between all 6 models, but this is not done here; see Brian O'Meara's AIC webpage)
#######################################################