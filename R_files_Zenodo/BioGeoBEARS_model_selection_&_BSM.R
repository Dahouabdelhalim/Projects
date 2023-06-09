#######################################################
# This is an adaptation of example scripts in BioGeoBEARS PhyloWiki for 'hyptree.tre' (Boana gr. pulchella)
# Running DEC with model selection in the R package BioGeoBEARS v.1.1.2 by Nick Matzke (downloaded on 20 Jul 2019)
# Using Biogeographical Stochastic Mapping (BSM) in the selected model M5
# by: Mariana M. Vasconcellos
#######################################################

#######################################################
# Installing BioGeoBEARS
#######################################################

#library(devtools)
#devtools::install_github(repo="nmatzke/BioGeoBEARS", INSTALL_opts="--byte-compile", dependencies=FALSE)

#######################################################
# SETUP -- libraries/BioGeoBEARS updates
#######################################################

# Load the package (after installation, see above).
library(GenSA)    # GenSA is better than optimx (although somewhat slower)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
# I have been putting the updates on CRAN/GitHub
# You should use:
# rexpokit version 0.26.6 from CRAN
# cladoRcpp version 0.15 from CRAN
packageVersion("BioGeoBEARS")
packageVersion("rexpokit")
packageVersion("cladoRcpp")

#######################################################
# SETUP: YOUR WORKING DIRECTORY
#######################################################
# You will need to set your working directory to match your local system

# Change this as you like
wd = np("~/GoogleDrive/MARIANA/Manuscripts/Meus/20_Boana_gr_pulchella/ANALYSES/BioGeoBEARS")
setwd(wd)

# Double-check your working directory with getwd()
getwd()

#######################################################
# SETUP: Extension data directory
#######################################################
# When R packages contain extra files, they are stored in the "extdata" directory 
# inside the installed package.
#
# BioGeoBEARS contains various example files and scripts in its extdata directory.
# 
# Each computer operating system might install BioGeoBEARS in a different place, 
# depending on your OS and settings. 
# 
# However, you can find the extdata directory like this:
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

# "system.file" looks in the directory of a specified package (in this case BioGeoBEARS)
# The function "np" is just a shortcut for normalizePath(), which converts the 
# path to the format appropriate for your system (e.g., Mac/Linux use "/", but 
# Windows uses "\\\\", if memory serves).

# Even when using your own data files, you should KEEP these commands in your 
# script, since the plot_BioGeoBEARS_results function needs a script from the 
# extdata directory to calculate the positions of "corners" on the plot. This cannot
# be made into a straight up BioGeoBEARS function because it uses C routines 
# from the package APE which do not pass R CMD check for some reason.

#######################################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE
#######################################################
# To run your own data,
# make the below lines point to your own files
trfn = "./input_files/hyptree_R.tre"
geogfn = "./input_files/input_BioGeoBEARS.txt"

#######################################################
# Phylogeny file
# Notes: 
# 1. Must be binary/bifurcating: no polytomies
# 2. No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
# 3. Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
# 4. You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the 
#    top of the tree as fossils!  This is only a good idea if you actually do have fossils in your tree,
#    as in e.g. Wood, Matzke et al. (2013), Systematic Biology.
# 5. The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of 
#    millions of years, and the tree is 1-1000 units tall. If you have a tree with a total height of
#    e.g. 0.00001, you will need to adjust e.g. the max values of d and e, or (simpler) multiply all
#    your branchlengths to get them into reasonable units.
# 6. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#######################################################

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny (plots to a PDF, which avoids issues with multiple graphics in same window):
pdffn = "output_files/B_pulchella_tree.pdf"
pdf(file=pdffn, width=9, height=12)
tr = read.tree(trfn)
tr
plot(tr)
title("Boana gr. pulchella phylogeny")
axisPhylo() # plots timescale
nodelabels()

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

#######################################################
# Geography file
# Notes:
# 1. This is a PHYLIP-formatted file. This means that in the 
#    first line, 
#    - the 1st number equals the number of rows (species)
#    - the 2nd number equals the number of columns (number of areas)
#    - after a tab, put the areas in parentheses, with spaces: (A B C D)
#
# 1.5. Example first line:
#    10    4    (A B C D)
# 
# 2. The second line, and subsequent lines:
#    speciesA    0110
#    speciesB    0111
#    speciesC    0001
#         ...
# 
# 2.5a. This means a TAB between the species name and the area 0/1s
# 2.5b. This also means NO SPACE AND NO TAB between the area 0/1s.
# 
# 3. See example files at:
#    http://phylo.wikidot.com/biogeobears#files
# 
# 4. Make you understand what a PLAIN-TEXT EDITOR is:
#    http://phylo.wikidot.com/biogeobears#texteditors
#
# 3. The PHYLIP format is the same format used for C++ LAGRANGE geography files.
#
# 4. All names in the geography file must match names in the phylogeny file.
#
# 5. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#
# 6. Operational taxonomic units (OTUs) should ideally be phylogenetic lineages, 
#    i.e. genetically isolated populations.  These may or may not be identical 
#    with species.  You would NOT want to just use specimens, as each specimen 
#    automatically can only live in 1 area, which will typically favor DEC+J 
#    models.  This is fine if the species/lineages really do live in single areas,
#    but you wouldn't want to assume this without thinking about it at least. 
#    In summary, you should collapse multiple specimens into species/lineages if 
#    data indicates they are the same genetic population.
######################################################

# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 3

####################################################
# KEY HINT: The number of states (= number of different possible geographic ranges)
# depends on (a) the number of areas and (b) max_range_size.
# If you have more than about 500-600 states, the calculations will get REALLY slow,
# since the program has to exponentiate a matrix of e.g. 600x600.  Often the computer
# will just sit there and crunch, and never get through the calculation of the first
# likelihood.
# 
# (this is also what is usually happening when LAGRANGE hangs: you have too many states!)
#
# To check the number of states for a given number of ranges, try:
numstates_from_numareas(numareas=6, maxareas=4, include_null_range=TRUE)
numstates_from_numareas(numareas=6, maxareas=4, include_null_range=FALSE)
numstates_from_numareas(numareas=6, maxareas=3, include_null_range=TRUE)
numstates_from_numareas(numareas=6, maxareas=2, include_null_range=TRUE)

# Large numbers of areas have problems:
numstates_from_numareas(numareas=10, maxareas=10, include_null_range=TRUE)

# ...unless you limit the max_range_size:
numstates_from_numareas(numareas=10, maxareas=2, include_null_range=TRUE)
####################################################

#######################################################
# DEC ANALYSES (6 MODEL COMPARISON)
#######################################################

#######################################################
# Run DEC - M0 Model
#######################################################
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 12
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

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
resfn = "output_files/Boana_DEC_M0_null_model.Rdata"
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
# PDF plots - M0 Model
#######################################################
pdffn = "output_files/Boana_DEC_M0_null_model.pdf"
pdf(pdffn, width=8.5, height=11)

#######################################################
# Plot ancestral states - M0 Model
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Boana M0 null model"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Run DEC - M1 Model
#######################################################
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "./input_files/manual_dispersal_multipliers_M1.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 12
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

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

# Fix node 36 (root) (edits Mariana Vasconcellos M1 model)
BioGeoBEARS_run_object$fixnode = 36
# Fix root for E state only (+ null state)
BioGeoBEARS_run_object$fixlikes = c(1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "output_files/Boana_DEC_M1_center_of_origin_E.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC_M1 = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC_M1 = res
}

#######################################################
# PDF plots - M1 Model
#######################################################
pdffn = "output_files/Boana_DEC_M1_center_of_orgin_E.pdf"
pdf(pdffn, width=8.5, height=11)

#######################################################
# Plot ancestral states - M1 Model
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Boana M1 center of origin E"

# Setup
results_object = resDEC_M1
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Run DEC - M2 Model
#######################################################
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "./input_files/manual_dispersal_multipliers_M2.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 12
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

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

# Fix node 36 (root) (edits Mariana Vasconcellos M2 model)
BioGeoBEARS_run_object$fixnode = 36
# Fix root for C state only (+ null state)
BioGeoBEARS_run_object$fixlikes = c(1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "output_files/Boana_DEC_M2_center_of_origin_C.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC_M2 = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC_M2 = res
}

#######################################################
# PDF plots - M2 Model
#######################################################
pdffn = "output_files/Boana_DEC_M2_center_of_orgin_C.pdf"
pdf(pdffn, width=8.5, height=11)

#######################################################
# Plot ancestral states - M2 Model
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Boana M2 center of origin C"

# Setup
results_object = resDEC_M2
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res3 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Run DEC - M3 Model
#######################################################
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "./input_files/manual_dispersal_multipliers_M3.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 12
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

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

# Fix node 36 (root) (edits Mariana Vasconcellos M3 model)
BioGeoBEARS_run_object$fixnode = 36
# Fix root for A state only (+ null state)
BioGeoBEARS_run_object$fixlikes = c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "output_files/Boana_DEC_M3_center_of_origin_A.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC_M3 = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC_M3 = res
}

#######################################################
# PDF plots - M3 Model
#######################################################
pdffn = "output_files/Boana_DEC_M3_center_of_orgin_A.pdf"
pdf(pdffn, width=8.5, height=11)

#######################################################
# Plot ancestral states - M3 Model
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Boana M3 center of origin A"

# Setup
results_object = resDEC_M3
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res4 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Run DEC - M4 Model
#######################################################
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input_files/manual_dispersal_multipliers_M4.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 12
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

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

########################
# Modify the states_list of allowed ranges manually for 3 max_range_size  (edits Mariana Vasconcellos)
BioGeoBEARS_run_object$states_list
areas = names(tipranges@df)
states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
states_list
TF = rep(TRUE, 42) ## Total number of ranges for 6 areas and 3 max_range_size = 42
# Removing 10 areas in your model (non-neighboring areas in the ancestral range)
TF[10] = FALSE
TF[12] = FALSE
TF[14] = FALSE
TF[19] = FALSE
TF[21] = FALSE
TF[24] = FALSE
TF[29] = FALSE
TF[31] = FALSE
TF[37] = FALSE
TF[40] = FALSE
states_list = states_list[TF]
#######################

# Put the new states_list into the BioGeoBEARS_run_object
BioGeoBEARS_run_object$states_list = states_list

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "output_files/Boana_DEC_M4_eco-similarity.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC_M4 = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC_M4 = res
}

#######################################################
# PDF plots - M4 Model
#######################################################
pdffn = "output_files/Boana_DEC_M4_eco-similarity.pdf"
pdf(pdffn, width=8.5, height=11)

#######################################################
# Plot ancestral states - M4 Model
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Boana M4 eco-similarity"

# Setup
results_object = resDEC_M4
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res5 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Run DEC - M5 Model
#######################################################
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input_files/manual_dispersal_multipliers_M5.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 12
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

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

########################
# Modify the states_list of allowed ranges manually for 3 max_range_size  (edits Mariana Vasconcellos)
BioGeoBEARS_run_object$states_list
areas = names(tipranges@df)
states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
states_list
TF = rep(TRUE, 42) ## Total number of ranges for 6 areas and 3 max_range_size = 42
# Removing 10 areas in your model (non-neighboring areas in the ancestral range)
TF[10] = FALSE
TF[12] = FALSE
TF[14] = FALSE
TF[19] = FALSE
TF[21] = FALSE
TF[24] = FALSE
TF[29] = FALSE
TF[31] = FALSE
TF[37] = FALSE
TF[40] = FALSE
states_list = states_list[TF]
#######################

# Put the new states_list into the BioGeoBEARS_run_object
BioGeoBEARS_run_object$states_list = states_list

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "output_files/Boana_DEC_M5_adjacency.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC_M5 = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC_M5 = res
}

#######################################################
# PDF plots - M5 Model
#######################################################
pdffn = "output_files/Boana_DEC_M5_adjacency.pdf"
pdf(pdffn, width=8.5, height=11)

#######################################################
# Plot ancestral states - M5 Model
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Boana M5 adjacency"

# Setup
results_object = resDEC_M5
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res6 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#########################################################################
#########################################################################
# 
# CALCULATE SUMMARY STATISTICS TO COMPARE MODELS
# 
#########################################################################
#########################################################################

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
# Statistics -- DEC
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC_M1)
LnL_3 = get_LnL_from_BioGeoBEARS_results_object(resDEC_M2)
LnL_4 = get_LnL_from_BioGeoBEARS_results_object(resDEC_M3)
LnL_5 = get_LnL_from_BioGeoBEARS_results_object(resDEC_M4)
LnL_6 = get_LnL_from_BioGeoBEARS_results_object(resDEC_M5)

numparams1 = 2
numparams2 = 2
numparams3 = 2
numparams4 = 2
numparams5 = 2
numparams6 = 2

stats = AICstats_2models(LnL_1, LnL_2, LnL_3, LnL_4, LnL_5, LnL_6, numparams1, numparams2, numparams3, numparams4, numparams5, numparams6)
stats

# DEC_M0, null model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC_M1, alternative model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC_M1, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC_M2, alternative model for Likelihood Ratio Test (LRT)
res3 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC_M2, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC_M3, alternative model for Likelihood Ratio Test (LRT)
res4 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC_M3, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC_M4, alternative model for Likelihood Ratio Test (LRT)
res5 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC_M4, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC_M5, alternative model for Likelihood Ratio Test (LRT)
res6 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC_M5, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
# confer the same likelihood on the data. See: Brian O'Meara's webpage:
# http://www.brianomeara.info/tutorials/aic
# ...for an intro to LRT, AIC, and AICc

rbind(res1, res2, res3, res4, res5, res6)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res1, res2, res3, res4, res5, res6)
teststable = rbind(teststable, tmp_tests)

#########################################################################
# ASSEMBLE RESULTS TABLES: DEC: M0 to M5
#########################################################################
teststable$null = c("DEC")
teststable$alt = c("DEC_M1")
teststable$alt = c("DEC_M2")
teststable$alt = c("DEC_M3")
teststable$alt = c("DEC_M4")
teststable$alt = c("DEC_M5")
row.names(restable) = c("DEC", "DEC_M1", "DEC_M2", "DEC_M3", "DEC_M4", "DEC_M5")
restable = put_jcol_after_ecol(restable)

# Look at the results!!
restable
teststable

#######################################################
# Model weights of all models
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

# Save with nice conditional formatting
write.table(conditional_format_table(restable_AIC_rellike), file="output_files/restable_AIC_rellike_formatted.txt", quote=FALSE, sep="\\t")
write.table(conditional_format_table(restable_AICc_rellike), file="output_files/restable_AICc_rellike_formatted.txt", quote=FALSE, sep="\\t")

# For more info on model selection see Brian O'Meara's AIC webpage: http://www.brianomeara.info/tutorials/aic )

###################################################################################
################### Biogeographical Stochastic Maps   #############################
#
# Non-time-stratified script for Biogeographical Stochastic Mapping (BSM)
# Adapted from http://phylo.wikidot.com/biogeographical-stochastic-mapping-example-script
#
###################################################################################
###################################################################################

model_name <- "DEC_M5"
outputdir <- "output_files_BSM/"
#######################################################
# Stochastic mapping
#######################################################
clado_events_tables <- NULL
ana_events_tables <- NULL
lnum <- 0

#######################################################
# From Matzke original script:
# Get the inputs for Biogeographical Stochastic Mapping
# Note: this can be slow for large state spaces and trees, since
# the independent likelihoods for each branch are being pre-calculated
# E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
# for each branch.  On a tree with ~800 tips and thus ~1600 branches,
# this was about 1.6 gigs
# for storage of "BSM_inputs_file.Rdata".
# Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
# the same settings will be used for get_inputs_for_stochastic_mapping().
#######################################################
BSM_inputs_fn <- paste0(outputdir, "BSM_inputs_file.Rdata")

runInputsSlow <- TRUE
if (runInputsSlow)
{
  stochastic_mapping_inputs_list <- get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load(BSM_inputs_fn)
} # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

# Generate stochastic maps
runBSMslow <- TRUE
if (runBSMslow == TRUE)
{
  # Saves to: RES_clado_events_tables.Rdata
  # Saves to: RES_ana_events_tables.Rdata
  BSM_output <- runBSM(res, stochastic_mapping_inputs_list =
                         stochastic_mapping_inputs_list,
                       maxnum_maps_to_try=100, nummaps_goal=100,
                       maxtries_per_branch=40000,
                       save_after_every_try=TRUE,
                       savedir = outputdir,
                       seedval=12345,
                       wait_before_save=0.01)
  RES_clado_events_tables <- BSM_output$RES_clado_events_tables
  RES_ana_events_tables <- BSM_output$RES_ana_events_tables
} else {
  # Load previously saved...
  # Loads to: RES_clado_events_tables
  load(file= paste0(outputdir, "RES_clado_events_tables.Rdata"))
  # Loads to: RES_ana_events_tables
  load(file= paste0(outputdir, "RES_ana_events_tables.Rdata"))
  BSM_output <- NULL
  BSM_output$RES_clado_events_tables <- RES_clado_events_tables
  BSM_output$RES_ana_events_tables <- RES_ana_events_tables
} # END if (runBSMslow == TRUE)

# Extract BSM output
clado_events_tables <- BSM_output$RES_clado_events_tables
ana_events_tables <- BSM_output$RES_ana_events_tables
# Just check data
# head(clado_events_tables[[1]])
# head(ana_events_tables[[1]])
# length(clado_events_tables)
# length(ana_events_tables)

include_null_range <- TRUE
areanames <- names(tipranges@df)
areas <- areanames
max_range_size <- 3

# Note: If you did something to change the states_list from the
# default given the number of areas, you would
# have to manually make that change here as well!
# (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based <-
  rcpp_areas_list_to_states_list(areas=areas,
                                 maxareas=max_range_size,
                                 include_null_range=include_null_range)

colors_list_for_states <-
  get_colors_for_states_list_0based(areanames=areanames,
                                    states_list_0based=states_list_0based,
                                    max_range_size=max_range_size,
                                    plot_null_range=TRUE)

########################
# Next lines used to remove non-adjacent ranges
TF <- rep(TRUE, 42) ## Total number of ranges for 6 areas and 3 max_range_size
## (for max_range_size = 4, Total number of ranges = 57)
# # Removing 10 areas in your model (removing non-adjacent areas in the ancestral range)
TF[10] <- FALSE
TF[12] <- FALSE
TF[14] <- FALSE
TF[19] <- FALSE
TF[21] <- FALSE
TF[24] <- FALSE
TF[29] <- FALSE
TF[31] <- FALSE
TF[37] <- FALSE
TF[40] <- FALSE
states_list_0based <- states_list_0based[TF]
colors_list_for_states <- colors_list_for_states[TF]
#########################

############################################
# Setup for painting a single stochastic map
############################################
scriptdir <- np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified <- FALSE
clado_events_table <- clado_events_tables[[1]]
ana_events_table <- ana_events_tables[[1]]

# cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
# colnums = match(cols_to_get, names(ana_events_table))
# ana_events_table_cols_to_add = ana_events_table[,colnums]
# anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
# ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add,
#      anagenetic_events_txt_below_node)
# rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
# master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)

############################################
# Open a PDF
############################################
pdffn = paste0(outputdir, model_name,
               "_single_stochastic_map_n1.pdf")
pdf(file=pdffn, width=8.5, height=11)

# Convert the BSM into a modified res object
master_table_cladogenetic_events <- clado_events_tables[[1]]
resmod <-
  stochastic_map_states_into_res(res=res,
                                 master_table_cladogenetic_events=
                                   master_table_cladogenetic_events,
                                 stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod,
                         analysis_titletxt= paste0(model_name, " Stochastic map"),
                         addl_params=list("j"),
                         label.offset=0.5, plotwhat="text",
                         cornercoords_loc=scriptdir,
                         root.edge=TRUE,
                         colors_list_for_states=colors_list_for_states,
                         skiptree=FALSE, show.tip.label=TRUE)

# Paint on the branch states
paint_stochastic_map_branches(res=resmod,
                              master_table_cladogenetic_events =
                                master_table_cladogenetic_events,
                              colors_list_for_states=colors_list_for_states,
                              lwd=5, lty=par("lty"), root.edge=TRUE,
                              stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod,
                         analysis_titletxt= paste0(model_name, " Stochastic map"),
                         addl_params=list("j"),
                         plotwhat="text", cornercoords_loc=scriptdir,
                         root.edge=TRUE,
                         colors_list_for_states=colors_list_for_states,
                         skiptree=TRUE, show.tip.label=TRUE)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Plot all 100 stochastic maps to PDF
#######################################################
# Setup
include_null_range <- include_null_range
areanames <-          areanames
areas <-              areanames
max_range_size <-     max_range_size

# Note: If you did something to change the states_list
# from the default given the number of areas, you would
# have to manually make that change here as well!
# (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based <-
  rcpp_areas_list_to_states_list(areas=areas,
                                 maxareas=max_range_size,
                                 include_null_range=include_null_range)
colors_list_for_states <-
  get_colors_for_states_list_0based(areanames= areanames,
                                    states_list_0based=
                                      states_list_0based,
                                    max_range_size=max_range_size,
                                    plot_null_range=TRUE)

########################
# Next lines used to remove non-adjacent ranges
TF <- rep(TRUE, 42) ## Total number of ranges for 6 areas and 3 max_range_size (for max_range_size = 4)
# # Removing 10 areas in your model (removing non-adjacent areas in the ancestral range)
TF[10] <- FALSE
TF[12] <- FALSE
TF[14] <- FALSE
TF[19] <- FALSE
TF[21] <- FALSE
TF[24] <- FALSE
TF[29] <- FALSE
TF[31] <- FALSE
TF[37] <- FALSE
TF[40] <- FALSE
states_list_0based <- states_list_0based[TF]
colors_list_for_states <- colors_list_for_states[TF]
#########################

scriptdir <- np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified <- stratified

# Loop through the maps and plot to PDF
pdffn <- paste0(outputdir,
                model_name, "_",
                length(clado_events_tables),
                "_StochasticMaps.pdf")
pdf(file=pdffn, width=8.5, height=11)

nummaps_goal <- 100
for (i in 1:nummaps_goal)
{
  clado_events_table <- clado_events_tables[[i]]
  analysis_titletxt <- paste0(model_name, " - Stochastic Map #",
                              i, "/", nummaps_goal)
  plot_BSM(results_object=res,
           clado_events_table=clado_events_table,
           stratified=stratified,
           analysis_titletxt=analysis_titletxt,
           addl_params=list("j"),
           label.offset=0.5, plotwhat="text",
           cornercoords_loc=scriptdir,
           root.edge=TRUE,
           colors_list_for_states=colors_list_for_states,
           show.tip.label=TRUE,
           include_null_range=include_null_range)
} # END for (i in 1:nummaps_goal)

dev.off()
cmdstr <- paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Summarize stochastic map tables
#######################################################
# length(clado_events_tables)
# length(ana_events_tables)

# head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

# head(ana_events_tables[[1]])
# tail(ana_events_tables[[1]])

areanames <- names(tipranges@df)
actual_names <- areanames
actual_names

# Get the dmat and times (if any)
dmat_times <- get_dmat_times_from_res(res=res, numstates=NULL)
dmat_times

# Extract BSM output
clado_events_tables <- BSM_output$RES_clado_events_tables
ana_events_tables <- BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas <-
  simulate_source_areas_ana_clado(res,
                                  clado_events_tables,
                                  ana_events_tables,
                                  areanames)
clado_events_tables <- BSMs_w_sourceAreas$clado_events_tables
ana_events_tables <- BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
# NOTE: THE DATA HERE ARE VEY USEFUL TO HISTOGRAMS!!
counts_list <- count_ana_clado_events(clado_events_tables,
                                      ana_events_tables,
                                      areanames, actual_names)

summary_counts_BSMs <- counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

# Histogram of event counts
hist_event_counts(counts_list,
                  pdffn=paste0(outputdir, model_name,
                               "_histograms_of_event_counts.pdf"))

#######################################################
# Print counts to files
#######################################################
# The info in counts_list is important for plotting histograms,
# Especially in the one-Dim lists, which are not output!! --DCC
tmpnames <- names(counts_list)
cat("\\n\\nWriting tables* of counts to tab-delimited text files:
    \\n(* = Tables have dimension=2 (rows and columns).
    Cubes (dimension 3) and lists (dimension 1) will not
    be printed to text files.) \\n\\n")
for (i in 1:length(tmpnames))
{
  cmdtxt <- paste0("item = counts_list$", tmpnames[i])
  eval(parse(text=cmdtxt))
  
  # Skip cubes
  if (length(dim(item)) != 2)
  {
    next()
  }
  
  outfn <- paste0(outputdir, tmpnames[i], "_", model_name, ".txt")
  
  if (length(item) == 0)
  {
    cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
    cat("\\n")
  } else {
    cat(outfn)
    cat("\\n")
    write.table(conditional_format_table(item),
                file=outfn, quote=FALSE,
                sep="\\t", col.names=TRUE,
                row.names=TRUE)
  } # END if (length(item) == 0)
} # END for (i in 1:length(tmpnames))
cat("...done.\\n")

#######################################################
# PLOT
# Check that ML ancestral state/range probabilities and
# the mean of the BSMs approximately line up
#######################################################
library(MultinomialCI)    # For 95% CIs on BSM counts
check_ML_vs_BSM(res, clado_events_tables,
                model_name, tr=NULL,
                plot_each_node=FALSE,
                linreg_plot=TRUE, MultinomialCI=TRUE)

# Extract extra info from tables
# Get data on type of cladogenetic event for each node,
# for each of 100 BSM realizations, and plot on tree
# using piecharts

library("tidyverse")
library("dplyr")
library("ggplot2")
myData <- dplyr::as_tibble(matrix(NA, ncol = 34, nrow = 100))
for (i in 1:100) {
  for (j in 1:34) {
    colnames(myData)[j] <- paste0("node", j+35)
    myData[i, j] <- clado_events_tables[[i]]$clado_event_type[j+35]
  }
}

myTidyData <- pivot_longer(myData, everything(), names_to = "node", values_to = "rangeshift")
# make a two-way table to use as source ofr piechart data on tree
myTable <- table(myTidyData$node, myTidyData$rangeshift)
myTable <- as.matrix(myTable)
# calculate proportions from raw numbers
myTable <- myTable/100
tr = read.tree("./input_files/hyptree_R.tre")
# truncate tiplabels to fit page
tr$tip.label <- str_sub(tr$tip.label, 1, 20)

# Set colors for piecharts and subtitle of plot
if (model_name == "DECJ")  {
  plot_subtitle <-  "red=jump, white=subset, gray=sympatry, green=vicariance"
  myColors <- c("red", "white", "gray", "green")
}  else {   # DEC
  plot_subtitle <-  "white=subset, gray=sympatry, green=vicariance"
  myColors <- c("white", "gray", "green")
}

# Open a PDF
pdffn = paste0(outputdir, model_name,
               "_Cladogenetic_Events.pdf")
pdf(file=pdffn, width=8.5, height=11)

plot(tr, sub = plot_subtitle)
title(paste0(model_name, ": Prob of types of cladogenetic events"))
axisPhylo()
# Add pie diagrams with proportion of type of cladogenetic event
nodelabels(pie = myTable, piecol = myColors, cex = 0.9)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

###  END  ###