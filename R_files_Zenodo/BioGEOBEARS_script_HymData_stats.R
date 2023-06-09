
#####install and load required packages
#install.packages("rexpokit") 
#install.packages("cladoRcpp") 
#install.packages("devtools") 
#library(devtools) 
#devtools::install_github(repo="nmatzke/BioGeoBEARS", force = TRUE) 
#install.packages("optimx") install.packages("GenSA") 
#install.packages("FD") 
#install.packages("snow")
#install.packages("parallel")

library(rexpokit)
library(cladoRcpp)
# library(optimx)
library(GenSA)     # GenSA is better than optimx (although somewhat slower)
library(FD)
library(snow)
library(parallel)
library(BioGeoBEARS)


# set working directory
setwd("/path/to/BioGeoBEARS/")


trfn<-np(paste("BioGeoBEARS_input_tree_pruned.newick", sep=""))
tr = read.tree(trfn)
# moref(trfn)       # have a look at tree

geogfn = np(paste("BioGeoBEARS_input_distribution.txt", sep=""))
# moref(geogfn)       # have a look at data

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
# tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

max_range_size=5 



runslow = FALSE

# Load DEC M0
resfn = "BioGeoBEARS_out_DECM0.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDEC = res
} else {
  load(resfn)
  resDEC = res
}

# Load DEC M1
resfn = "BioGeoBEARS_out_DECM1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDECM1 = res
} else {
  load(resfn)
  resDECM1 = res
}

# Load DEC + J M0
resfn = "BioGeoBEARS_out_DECjM0.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDECj = res
} else {
  load(resfn)
  resDECj = res
}

# Load DEC + J M1
resfn = "BioGeoBEARS_out_DECjM1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDECjM1 = res
} else {
  load(resfn)
  resDECjM1 = res
}

# Load DIVALIKE M0
resfn = "BioGeoBEARS_out_DIVALIKEM0.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  load(resfn)
  resDIVALIKE = res
}

# Load DIVALIKE M1
resfn = "BioGeoBEARS_out_DIVALIKEM1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDIVALIKEM1 = res
} else {
  load(resfn)
  resDIVALIKEM1 = res
}

# Load DIVALIKE + J M0
resfn = "BioGeoBEARS_out_DIVALIKEjM0.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDIVALIKEj = res
} else {
  load(resfn)
  resDIVALIKEj = res
}

# Load DIVALIKE + J M1
resfn = "BioGeoBEARS_out_DIVALIKEjM1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDIVALIKEjM1 = res
} else {
  load(resfn)
  resDIVALIKEjM1 = res
}

# Load BAYAREALIKE M0
resfn = "BioGeoBEARS_out_BAYAREALIKEM0.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
  load(resfn)
  resBAYAREALIKE = res
}

# Load BAYAREALIKE M1
resfn = "BioGeoBEARS_out_BAYAREALIKEM1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resBAYAREALIKEM1 = res
} else {
  load(resfn)
  resBAYAREALIKEM1 = res
}

# Load BAYAREALIKE + J M0
resfn = "BioGeoBEARS_out_BAYAREALIKEjM0.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resBAYAREALIKEj = res
} else {
  load(resfn)
  resBAYAREALIKEj = res
}

# Load BAYAREALIKE + J M1
resfn = "BioGeoBEARS_out_BAYAREALIKEjM1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resBAYAREALIKEjM1 = res
} else {
  load(resfn)
  resBAYAREALIKEjM1 = res
}


#######################################################
# Likelihood for each model 
#######################################################
# LnL - DEC M0
get_LnL_from_BioGeoBEARS_results_object(resDEC)

# LnL - DEC M1
get_LnL_from_BioGeoBEARS_results_object(resDECM1)

# LnL - DEC+J M0
get_LnL_from_BioGeoBEARS_results_object(resDECj)

# LnL - DEC+J M1
get_LnL_from_BioGeoBEARS_results_object(resDECjM1)

# LnL - DIVALIKE M0
get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)

# LnL - DIVALIKE M1
get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEM1)

# LnL - DIVALIKE+J M0
get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

# LnL - DIVALIKE+J M1
get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEjM1)

# LnL - BAYAREALIKE M0
get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)

# LnL - BAYAREALIKE M1
get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEM1)

# LnL - BAYAREALIKE+J M0
get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

# LnL - BAYAREALIKE+J M1
get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEjM1)


#########################################################################
# 
# CALCULATE SUMMARY STATISTICS TO COMPARE
# DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J (time-stratified/dispersal multipliers vs. not)
# 
#########################################################################
# REQUIRED READING:
#
# Practical advice / notes / basic principles on statistical model 
#    comparison in general, and in BioGeoBEARS:
# http://phylo.wikidot.com/advice-on-statistical-model-comparison-in-biogeobears
#########################################################################

# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

#######################################################
# Statistics -- DEC vs. DEC+J M0
#######################################################
# We have to extract the log-likelihood differently, depending on the version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models confer the same likelihood on the data. See: Brian O'Meara's webpage: http://www.brianomeara.info/tutorials/aic  ...for an intro to LRT, AIC, and AICc
rbind(res2, res1)
tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+J M0
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
rbind(res2, res1)
tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+J M0
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
rbind(res2, res1)
tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DEC vs. DEC+J M1
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDECM1)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECjM1)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECM1, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECjM1, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
rbind(res2, res1)
tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)


#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+J M1
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEM1)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEjM1)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEM1, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEjM1, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
rbind(res2, res1)
conditional_format_table(stats)
tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)



#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+J M1
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEM1)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEjM1)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEM1, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEjM1, returnwhat="table", addl_params=c("j","w"), paramsstr_digits=4)
rbind(res2, res1)
conditional_format_table(stats)
tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)



# ASSEMBLE RESULTS TABLES
############################
teststable$alt = c("DEC+J M0","DIVALIKE+J M0", "BAYAREALIKE+J M0", "DEC+J M1", "DIVALIKE+J M1","BAYAREALIKE+J M1")
teststable$null = c("DEC M0","DIVALIKE M0", "BAYAREALIKE M0", "DEC M1", "DIVALIKE M1", "BAYAREALIKE M1")
row.names(restable) = c("DEC M0", "DEC+J M0","DIVALIKE M0", "DIVALIKE+J M0", "BAYAREALIKE M0", "BAYAREALIKE+J M0", "DEC M1", "DEC+J M1", "DIVALIKE M1", "DIVALIKE+J M1", "BAYAREALIKE M1", "BAYAREALIKE+J M1")
restable = put_jcol_after_ecol(restable)

# Look at the results!!
restable
teststable

#######################################################
# Model weights of all six models
#######################################################
restable2 = restable

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike

# With AICcs -- factors in sample size
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike

#######################################################
# Save the results tables for later -- check for e.g.
# convergence issues
#######################################################

# Loads to "restable"
save(restable, file="restable_v1.Rdata")
load(file="restable_v1.Rdata")

# Loads to "teststable"
save(teststable, file="teststable_v1.Rdata")
load(file="teststable_v1.Rdata")

# Also save to text files
write.table(restable, file="restable.txt", quote=FALSE, sep="\\t")
write.table(unlist_df(teststable), file="teststable.txt", quote=FALSE, sep="\\t")

# Also save to text files
write.table(restable_AIC_rellike, file="restable_AIC_rellike.txt", quote=FALSE, sep="\\t")
write.table(restable_AICc_rellike, file="restable_AICc_rellike.txt", quote=FALSE, sep="\\t")

# Save with nice conditional formatting
write.table(conditional_format_table(restable_AIC_rellike), file="restable_AIC_rellike_formatted.txt", quote=FALSE, sep="\\t")
write.table(conditional_format_table(restable_AICc_rellike), file="restable_AICc_rellike_formatted.txt", quote=FALSE, sep="\\t")

