### updated on 11Mar18


# ###description of program and list of what software and files are needed
# this script is written to perform a post-hoc simulation of data that matches empirical data
# the simulated data are then used as input for species delimitation with spedeSTEM
# in Jacobs etal 2018, we perform the simulation-spedeSTEM analysis 100 independent times
# after which, 10 simulated datasets are randomly drawn to perform species delimitation with BPP
 
# ###software that the script needs access to:
# R
# perl
# grep
# ms (to simulate geneologies: Hudson RR (2002) Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinformatics, 18, 337–338)
# seqgen (to simulate sequences: Rambaut A, Grass NC (1997) Seq-Gen: an application for the Monte Carlo simulation of DNA sequence evolution along phylogenetic trees. Bioinformatics, 13, 235–238.)
# phyutility (to concatenate simulated data: Smith SA, Dunn CW (2008) Phyutility: a phyloinformatics tool for trees, alignments and molecular data. Bioinformatics, 24, 715–716.)
# subsample.pl (to subsample simulated alignments: an altered version of Step1Paup.pl (from STEM distribution - Hird S, Kubatko L, Carstens B (2010) Rapid and accurate species tree estimation for phylogeographic investigations using replicated subsampling. Molecular Phylogenetics and Evolution, 57, 888-898.)
#   we utilize subsampling in the spedeSTEM analysis and our code incorporates this at step 3 (below)
#   if you do not wish to subsample, you will need to modify the code to omit step 3 (in possibly step 2, part c, where we sort sequences in order to subsample appropriately)
#	and make sure that spedeSTEM has access to your alignments for step 4
# spedeSTEM (performs species delimitation: Ence DD, Carstens BC (2010) SpedeSTEM: a rapid and accurate method for species delimitation. Molecular Ecology Resources, 11, 473–480.)


# ###input files and empirical information needed: refer to Jacobs etal 2018 supplemental data for more detail
# geneology simulations: number of populations and how many individuals per population, estimates of modern and ancestral population sizes, estimates of the timing of population merging
# geneology scaling: mean tree height of empirical trees
# 	you will need to make one file for each gene tree that you wish to scale as lined out in spedeSTEM documentation (these are used in step4 - ours are called scale_chl.txt and scale_nuc.txt)
# sequence simulation: number and length of partitions in empirical data, empirical models of nucleotide evolution

# ###general note:
# this script is not formatted for general use, or with variables that you can quickly plug your own information into
# we currently have all generated files named for our purposes
# all parameter settings correspond with our study system and WILL NEED TO BE CHANGED to suit your needs
# this script is only intended to provide a scaffold for you to build your own post-hoc simulations around


#
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
make_option(c("-s", "--sim"), type='integer', default=0, help="The number of simulations you want to do, yo [default %default]", dest="sim_number"),
make_option(c("-n", "--nreps"), type='integer', default=0, help="The number of replicates, within each simulation, you would like to perform [default %default]", dest="nreps_number")
)
opt <- parse_args(OptionParser(option_list=option_list))



#Setup the environment
# set up and designate the paths to the following
mywd <- getwd()
setwd(mywd)
perl <- "/opt/modules/devel/perl/5.18.0/bin/perl"
grep <- "/bin/grep"
ms_app <- "/mnt/lfs2/sjjacobs/apps/msdir/ms"
seqgen_app <- "/mnt/lfs2/sjjacobs/apps/Seq-Gen.v1.3.3/source/seq-gen"
phyutility_app <- "java -jar /mnt/lfs2/sjjacobs/apps/phyutility_old/phyutility.jar"



power_analysis <- function(sim,nreps) {
##sim=1
##nreps=2
  
  #simFolder is the term used to move to the primary Simulation folder.
  #Wherever the code is inititated from, a folder called 'Simulation' (plus a number indicating which simulation), will be added to the working directory path 
  #within the Simulation.'sim' folder, the folder 'filesUsed' will contain all files used during the analysis
  simFolder = file.path(getwd(),paste("Simulation",sim,sep="."))
  dir.create(simFolder)
  filesUsed = file.path(simFolder,"filesUsed")
  dir.create(filesUsed)
  spedeSTEM <- file.path(mywd,"/apps/SpedeSTEM-master/SpedeSTEM_2.py")
  
  #step 1, simulate a gene tree for both loci
  # currently set to estimate a chloroplast geneology and a nuclear geneology; add additional lines as needed
  setwd(simFolder)
  call.gen.chl <- paste(ms_app,"20 1 -T -t 1232.5 -I 4 13 3 3 1 -n 1 0.5 -n 2 0.001 -n 3 0.001 -n 4 0.5 -ej 0.2 3 1 -ej 0.55 2 1 -ej 3.7 4 1 | tail -n +4 | grep -v // | grep \\\\) >",file.path(simFolder,"cast.sim.gen.chl.tre"))
  call.gen.nuc <- paste(ms_app,"20 1 -T -t 5786.12 -I 4 13 3 3 1 -n 1 0.5 -n 2 0.001 -n 3 0.001 -n 4 0.5 -ej 0.1 3 1 -ej 0.275 2 1 -ej 1.85 4 1 | tail -n +4 | grep -v // | grep \\\\) >",file.path(simFolder,"cast.sim.gen.nuc.tre"))
  system(call.gen.chl)
  system(call.gen.nuc)
  setwd(mywd)
  
  
  #step 2, [a] simulate sequences for each partition based on the gene trees, [b] concatenate the partitions, and then [c] sort the simulated sequences sequentially (necessary for appropriate subsampling)
  # change part [a] according to how many partitions you need
  # change part [c] according to 
  ##[a]##
  call.seq.chl.par1 <- paste(seqgen_app,"-mHKY -l4386 -f0.3033,0.1286,0.158,0.4101 -i0.939 -t0.5 -s0.003 -on <",file.path(simFolder,"cast.sim.gen.chl.tre"), ">",file.path(simFolder,"cast.sim.seq.temp.chl.par1.nex"))
  call.seq.chl.par2 <- paste(seqgen_app,"-mHKY -l11921 -f0.3496,0.177,0.1676,0.3059 -i0.92 -t0.5 -s0.003 -on <",file.path(simFolder,"cast.sim.gen.chl.tre"), ">",file.path(simFolder,"cast.sim.seq.temp.chl.par2.nex"))
  call.seq.chl.par3 <- paste(seqgen_app,"-mHKY -l1643 -s0.003 -a0.768 -g4 -i0.947 -f0.3408,0.1523,0.13,0.3769 -t0.5 -on <",file.path(simFolder,"cast.sim.gen.chl.tre"), ">",file.path(simFolder,"cast.sim.seq.temp.chl.par3.nex"))
  call.seq.chl.par4 <- paste(seqgen_app,"-mHKY -l6339 -f0.2738,0.1983,0.1837,0.3441 -i0.95 -s0.003 -on <",file.path(simFolder,"cast.sim.gen.chl.tre"), ">",file.path(simFolder,"cast.sim.seq.temp.chl.par4.nex"))
  call.seq.chl.par5 <- paste(seqgen_app,"-mHKY -l508 -a0.272 -g4 -i0.252 -f0.3409,0.1491,0.1503,0.3597 -t0.5 -s0.003 -on <",file.path(simFolder,"cast.sim.gen.chl.tre"), ">",file.path(simFolder,"cast.sim.seq.temp.chl.par5.nex"))
  call.seq.chl.par6 <- paste(seqgen_app,"-mHKY -l554 -f0.3779,0.1557,0.1888,0.2776 -t0.5 -s0.003 -on <",file.path(simFolder,"cast.sim.gen.chl.tre"), ">",file.path(simFolder,"cast.sim.seq.temp.chl.par6.nex"))
  
  call.seq.nuc.par1 <- paste(seqgen_app,"-mGTR -l450 -s0.015 -a15.988 -g4 -f0.1579,0.2779,0.2859,0.2783 -r0.91238,1.63328,2.78612,0.32558,1.63328,1 -on <",file.path(simFolder,"cast.sim.gen.nuc.tre"),">",file.path(simFolder,"cast.sim.seq.temp.nuc.par1.nex"))
  call.seq.nuc.par2 <- paste(seqgen_app,"-mGTR -l689 -s0.015 -i0.763 -f0.1888,0.3222,0.3051,0.1839 -r0.6289,1.0201,0.9926,0.0134,2.2502,1 -on <",file.path(simFolder,"cast.sim.gen.nuc.tre"),">",file.path(simFolder,"cast.sim.seq.temp.nuc.par2.nex"))
  
  system(call.seq.chl.par1)
  system(call.seq.chl.par2)
  system(call.seq.chl.par3)
  system(call.seq.chl.par4)
  system(call.seq.chl.par5)
  system(call.seq.chl.par6)
  system(call.seq.nuc.par1)
  system(call.seq.nuc.par2)
  
  ##[b]##
  setwd(simFolder)
  call.phyutility.chl <- paste(phyutility_app, "-concat -in cast.sim.seq.temp.chl.par1.nex cast.sim.seq.temp.chl.par2.nex cast.sim.seq.temp.chl.par3.nex cast.sim.seq.temp.chl.par4.nex cast.sim.seq.temp.chl.par5.nex cast.sim.seq.temp.chl.par6.nex -out amb15.chl.sim.concat.nex", sep=" ")
  call.phyutility.nuc <- paste(phyutility_app, "-concat -in cast.sim.seq.temp.nuc.par1.nex cast.sim.seq.temp.nuc.par2.nex -out amb15.nuc.sim.concat.nex", sep=" ")
  system(call.phyutility.chl)
  system(call.phyutility.nuc)
  
  ##[c]##
  chl.nex <- readLines(file.path(simFolder,"amb15.chl.sim.concat.nex"))
  chl.nex[7:26] <- chl.nex[7:26][order(as.numeric(sapply(strsplit(chl.nex[7:26],split="\\\\t"),"[[",2L)))]
  writeLines(chl.nex,file.path(simFolder,"amb15.chl.sim.seq.sorted.nex"))
  nuc.nex <- readLines(file.path(simFolder,"amb15.nuc.sim.concat.nex"))
  nuc.nex[7:26] <- nuc.nex[7:26][order(as.numeric(sapply(strsplit(nuc.nex[7:26],split="\\\\t"),"[[",2L)))]
  writeLines(nuc.nex,file.path(simFolder,"amb15.nuc.sim.seq.sorted.nex"))
  
  ## clean up files
  file.rename((file.path(simFolder,"cast.sim.gen.chl.tre")), file.path(filesUsed,"cast.sim.gen.chl.tre"))
  file.rename((file.path(simFolder,"cast.sim.gen.nuc.tre")), file.path(filesUsed,"cast.sim.gen.nuc.tre"))
  file.remove(file.path(simFolder,"cast.sim.seq.temp.chl.par1.nex"))
  file.remove(file.path(simFolder,"cast.sim.seq.temp.chl.par2.nex"))
  file.remove(file.path(simFolder,"cast.sim.seq.temp.chl.par3.nex"))
  file.remove(file.path(simFolder,"cast.sim.seq.temp.chl.par4.nex"))
  file.remove(file.path(simFolder,"cast.sim.seq.temp.chl.par5.nex"))
  file.remove(file.path(simFolder,"cast.sim.seq.temp.chl.par6.nex"))
  file.remove(file.path(simFolder,"cast.sim.seq.temp.nuc.par1.nex"))
  file.remove(file.path(simFolder,"cast.sim.seq.temp.nuc.par2.nex"))
  
  #step 3, subsample the simulated sequences
  call.subsample_chl <- paste(perl, file.path(mywd,'subsampling.pl'), sim, mywd, "chl", nreps, simFolder) ## STILL NEED TO FIX CHL TO BE LIKE NUC
  call.subsample_nuc <- paste(perl, file.path(mywd,'subsampling.pl'), sim, mywd, "nuc", nreps, simFolder)
  system(call.subsample_chl)
  system(call.subsample_nuc) 
  setwd(simFolder)
  
  #step 4, prepare the subsampled trees for spedeSTEM analysis
  #this includes cleaning the genetrees in spedeSTEM (removes zero-length branches and adds a scaling factor to each tree),
  #merging a single tree from each locus into a new tree file
  system(paste("cp -r",file.path(mywd,"SpedeSTEM-master"),file.path(simFolder,"SpedeSTEM-master"),sep=" "))
  writeLines(c(rep('[6.8]',nreps)),"scale_chl.txt")
  clean.spedeSTEM <- paste("python ",spedeSTEM, " setup -s", " scale_chl.txt"," -c amb15.chl.sim.seq.sorted.ALL3.*",sep="")
  system(clean.spedeSTEM)  
  writeLines(c(rep('[1.0]',nreps)),"scale_nuc.txt")
  clean.spedeSTEM <- paste("python ",spedeSTEM, " setup -s", " scale_nuc.txt"," -c amb15.nuc.sim.seq.sorted.ALL3.*",sep="")
  system(clean.spedeSTEM)
  setwd(mywd)                
  for (i in seq.int(1,nreps)) {
    nuc <- readLines(file.path(simFolder,paste("cleaned.amb15.nuc.sim.seq.sorted.ALL3.",i,".phy",sep="")))                   
    chl <- readLines(file.path(simFolder,paste("cleaned.amb15.chl.sim.seq.sorted.ALL3.",i,".phy",sep="")))                   
    writeLines(c(chl,nuc),file.path(simFolder,paste("ReplicateTree_",i,".tre",sep="")))
  }
  #here, move all files -- except ReplicateTree files -- into the 'filesUsed' directory
  #this is necessary for the spedeSTEM analysis - the only files that can be in the primary folder are the tree files used for the analysis; otherwise, spedeSTEM complains
  setwd(simFolder)
  system(paste("mv amb15.* ",file.path(filesUsed),sep=""))
  system(paste("mv cleaned.* ",file.path(filesUsed),sep=""))
  system(paste("mv paupfeed.* ",file.path(filesUsed),sep=""))
  system(paste("mv log.log ", file.path(filesUsed),sep=""))
  system(paste("mv seedms ", file.path(filesUsed),sep=""))
  system(paste("mv scale_* ", file.path(filesUsed),sep=""))
  system(paste("mv phyutility.log ", file.path(filesUsed),sep=""))
  system(paste("mv out.txt ", file.path(filesUsed),sep=""))
  
  
  #  system(paste("cp -r",file.path(mywd,"SpedeSTEM-master"),file.path(simFolder,"SpedeSTEM-master"),sep=" "))
  #step 7, perform spedestem analysis on simulated data
  setwd(file.path(simFolder,"SpedeSTEM-master"))
  #  system(paste("mkdir spedestem_out.rep.",sim,sep=""))
  validate.spedeSTEM <- paste("python ",spedeSTEM, " validation -d ",file.path(simFolder)," -a ",file.path(mywd,"associations.stem")," -s ",file.path(mywd,"settings.stem"),sep="")
  system(validate.spedeSTEM)
  system(paste("mv boot.tre ",file.path(filesUsed),sep=""))
  system(paste("mv genetrees.tre ",file.path(filesUsed),sep=""))
  system(paste("mv itTable.txt ",file.path(filesUsed),sep=""))
  system(paste("mv mle.tre ",file.path(filesUsed),sep=""))
  system(paste("mv results.txt ",file.path(simFolder),sep=""))
  system(paste("mv settings ",file.path(filesUsed),sep=""))
  system(paste("mv stemOut.txt ",file.path(filesUsed),sep=""))
  setwd(simFolder)
  unlink("SpedeSTEM-master/",recursive=TRUE)  
  
  #step 8, parse results file and calculate statistics
  results.temp <- read.delim("results.txt", sep=";", header=FALSE)
  results.parsed <- as.data.frame(matrix(results.temp$V3,ncol=nreps))
  rownames(results.parsed) <- results.temp$V1[1:nrow(results.parsed)]
  colnames(results.parsed) <- paste("Rep",seq.int(1,nreps),sep="")
  
  #  results.parsed[c("lnLAvg", "k", "AIC", "delta-i", "Model likelihood", "wi")] <- NA   #add new rows at end of data frame for statistics calculations
  results.parsed["lnL Avg"] <- rowMeans((results.parsed), na.rm=TRUE)   #calculate average log likelihood
  results.parsed["k"] <- results.temp$V2[1:nrow(results.parsed)] #enter number of parameters for each lineage-composition model; this is equal to the number of lineages in model (not including outgroup)
  results.parsed$AIC <- ((-2 * as.numeric(results.parsed$"lnL Avg")) + (2 * as.numeric(results.parsed$"k")))  #calculate AIC scores
  results.parsed <- results.parsed[order(results.parsed$"AIC"),]  #rank rows according to AIC scores; lowest score at the top
  
  results.parsed$"delta-i" <- c(0,diff(results.parsed$AIC))
  results.parsed$"Model likelihood" <- exp(-0.5 * abs(as.numeric(results.parsed$"delta-i")))  #calculate model likelihood for each lineage-composition model
  
  results.parsed$"wi" <- (as.numeric(results.parsed$"Model likelihood"))/sum(results.parsed$"Model likelihood")
  write.csv(results.parsed, file="results.parsed.csv")
  system(paste("mv ReplicateTree_* ",file.path(filesUsed),sep=""))
  system(paste("mv results.txt ",file.path(filesUsed),sep=""))
}


#for (i in 1:10){
#power_analysis(i,10)
#}


power_analysis(opt$sim_number,opt$nreps_number)
