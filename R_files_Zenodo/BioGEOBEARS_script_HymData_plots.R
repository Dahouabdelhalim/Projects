
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


#######################################################
# Plot per-area probabilities -- with COLOR
# 
# (This assumes you have run a standard BioGeoBEARS analysis
#  from the main example script first)
#######################################################

scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

areas = getareas_from_tipranges_object(tipranges)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

barwidth_proportion = 0.004
barheight_proportion = 0.0025
split_reduce_x = 0.85    # fraction of sizes for the split probabilities
split_reduce_y = 0.85    # fraction of sizes for the split probabilities

# Manual offsets, if desired (commented out below)
offset_nodenums = NULL
offset_xvals = NULL
offset_yvals = NULL

root.edge = TRUE

#numstates_from_numareas(numareas=8, maxareas=5, include_null_range=TRUE)

colors_list_for_states  = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666", rep("grey50", times=211))

# # # if plotted as PDF
# pdffn = "BioGeoBEARS_Hym_per-area_probabilities_DECjM1.pdf"
# pdf(pdffn, width=60, height=80)

# # #   DEC + J M1   # # #
# titletxt = "Probability of occupancy per area\\nunder DEC + J M1 analysis"
probs_each_area = plot_per_area_probs(tr, res=resDECjM1, areas=areas, states_list_0based=states_list_0based, titletxt=titletxt, cols_each_area=TRUE, barwidth_proportion=barwidth_proportion, barheight_proportion=barheight_proportion, offset_nodenums=offset_nodenums, offset_xvals=offset_xvals, offset_yvals=offset_yvals, root.edge=root.edge)

write.csv(probs_each_area, file = "BioGeoBEARS_DECjM1_areaprobs.csv")

# # # if plotted as PDF
# print("Waiting...")
# Sys.sleep(2)    # wait for this to finish
# dev.off()
# cmdstr = paste("open ", pdffn, sep="")
# system(cmdstr)

