#load the pophelper library 
library(pophelper)

#####plotting delta K Evanno method
setwd("all_results_files/") #this folder contains 400 structure results files, one file for each K value and replicate (K1-20, 20 replicate runs each) 
sfiles<-c("Anolis_k1_run1_f","Anolis_k1_run2_f","Anolis_k1_run3_f","Anolis_k1_run4_f","Anolis_k1_run5_f","Anolis_k1_run6_f","Anolis_k1_run7_f",
          "Anolis_k1_run8_f","Anolis_k1_run9_f","Anolis_k1_run10_f","Anolis_k1_run11_f","Anolis_k1_run12_f","Anolis_k1_run13_f",
          "Anolis_k1_run14_f","Anolis_k1_run15_f","Anolis_k1_run16_f","Anolis_k1_run17_f","Anolis_k1_run18_f","Anolis_k1_run19_f","Anolis_k1_run20_f",
          "Anolis_k2_run1_f","Anolis_k2_run2_f","Anolis_k2_run3_f","Anolis_k2_run4_f","Anolis_k2_run5_f","Anolis_k2_run6_f","Anolis_k2_run7_f",
          "Anolis_k2_run8_f","Anolis_k2_run9_f","Anolis_k2_run10_f","Anolis_k2_run11_f","Anolis_k2_run12_f","Anolis_k2_run13_f",
          "Anolis_k2_run14_f","Anolis_k2_run15_f","Anolis_k2_run16_f","Anolis_k2_run17_f","Anolis_k2_run18_f","Anolis_k2_run19_f","Anolis_k2_run20_f",
          "Anolis_k3_run1_f","Anolis_k3_run2_f","Anolis_k3_run3_f","Anolis_k3_run4_f","Anolis_k3_run5_f","Anolis_k3_run6_f","Anolis_k3_run7_f",
          "Anolis_k3_run8_f","Anolis_k3_run9_f","Anolis_k3_run10_f","Anolis_k3_run11_f","Anolis_k3_run12_f","Anolis_k3_run13_f",
          "Anolis_k3_run14_f","Anolis_k3_run15_f","Anolis_k3_run16_f","Anolis_k3_run17_f","Anolis_k3_run18_f","Anolis_k3_run19_f","Anolis_k3_run20_f",
          "Anolis_k4_run1_f","Anolis_k4_run2_f","Anolis_k4_run3_f","Anolis_k4_run4_f","Anolis_k4_run5_f","Anolis_k4_run6_f","Anolis_k4_run7_f",
          "Anolis_k4_run8_f","Anolis_k4_run9_f","Anolis_k4_run10_f","Anolis_k4_run11_f","Anolis_k4_run12_f","Anolis_k4_run13_f",
          "Anolis_k4_run14_f","Anolis_k4_run15_f","Anolis_k4_run16_f","Anolis_k4_run17_f","Anolis_k4_run18_f","Anolis_k4_run19_f","Anolis_k4_run20_f",
          "Anolis_k5_run1_f","Anolis_k5_run2_f","Anolis_k5_run3_f","Anolis_k5_run4_f","Anolis_k5_run5_f","Anolis_k5_run6_f","Anolis_k5_run7_f",
          "Anolis_k5_run8_f","Anolis_k5_run9_f","Anolis_k5_run10_f","Anolis_k5_run11_f","Anolis_k5_run12_f","Anolis_k5_run13_f",
          "Anolis_k5_run14_f","Anolis_k5_run15_f","Anolis_k5_run16_f","Anolis_k5_run17_f","Anolis_k5_run18_f","Anolis_k5_run19_f","Anolis_k5_run20_f",
          "Anolis_k6_run1_f","Anolis_k6_run2_f","Anolis_k6_run3_f","Anolis_k6_run4_f","Anolis_k6_run5_f","Anolis_k6_run6_f","Anolis_k6_run7_f",
          "Anolis_k6_run8_f","Anolis_k6_run9_f","Anolis_k6_run10_f","Anolis_k6_run11_f","Anolis_k6_run12_f","Anolis_k6_run13_f",
          "Anolis_k6_run14_f","Anolis_k6_run15_f","Anolis_k6_run16_f","Anolis_k6_run17_f","Anolis_k6_run18_f","Anolis_k6_run19_f","Anolis_k6_run20_f",
          "Anolis_k7_run1_f","Anolis_k7_run2_f","Anolis_k7_run3_f","Anolis_k7_run4_f","Anolis_k7_run5_f","Anolis_k7_run6_f","Anolis_k7_run7_f",
          "Anolis_k7_run8_f","Anolis_k7_run9_f","Anolis_k7_run10_f","Anolis_k7_run11_f","Anolis_k7_run12_f","Anolis_k7_run13_f",
          "Anolis_k7_run14_f","Anolis_k7_run15_f","Anolis_k7_run16_f","Anolis_k7_run17_f","Anolis_k7_run18_f","Anolis_k7_run19_f","Anolis_k7_run20_f",
          "Anolis_k8_run1_f","Anolis_k8_run2_f","Anolis_k8_run3_f","Anolis_k8_run4_f","Anolis_k8_run5_f","Anolis_k8_run6_f","Anolis_k8_run7_f",
          "Anolis_k8_run8_f","Anolis_k8_run9_f","Anolis_k8_run10_f","Anolis_k8_run11_f","Anolis_k8_run12_f","Anolis_k8_run13_f",
          "Anolis_k8_run14_f","Anolis_k8_run15_f","Anolis_k8_run16_f","Anolis_k8_run17_f","Anolis_k8_run18_f","Anolis_k8_run19_f","Anolis_k8_run20_f",
          "Anolis_k9_run1_f","Anolis_k9_run2_f","Anolis_k9_run3_f","Anolis_k9_run4_f","Anolis_k9_run5_f","Anolis_k9_run6_f","Anolis_k9_run7_f",
          "Anolis_k9_run8_f","Anolis_k9_run9_f","Anolis_k9_run10_f","Anolis_k9_run11_f","Anolis_k9_run12_f","Anolis_k9_run13_f",
          "Anolis_k9_run14_f","Anolis_k9_run15_f","Anolis_k9_run16_f","Anolis_k9_run17_f","Anolis_k9_run18_f","Anolis_k9_run19_f","Anolis_k9_run20_f",
          "Anolis_k10_run1_f","Anolis_k10_run2_f","Anolis_k10_run3_f","Anolis_k10_run4_f","Anolis_k10_run5_f","Anolis_k10_run6_f","Anolis_k10_run7_f",
          "Anolis_k10_run8_f","Anolis_k10_run9_f","Anolis_k10_run10_f","Anolis_k10_run11_f","Anolis_k10_run12_f","Anolis_k10_run13_f",
          "Anolis_k10_run14_f","Anolis_k10_run15_f","Anolis_k10_run16_f","Anolis_k10_run17_f","Anolis_k10_run18_f","Anolis_k10_run19_f","Anolis_k10_run20_f",
          "Anolis_k11_run1_f","Anolis_k11_run2_f","Anolis_k11_run3_f","Anolis_k11_run4_f","Anolis_k11_run5_f","Anolis_k11_run6_f","Anolis_k11_run7_f",
          "Anolis_k11_run8_f","Anolis_k11_run9_f","Anolis_k11_run10_f","Anolis_k11_run11_f","Anolis_k11_run12_f","Anolis_k11_run13_f",
          "Anolis_k11_run14_f","Anolis_k11_run15_f","Anolis_k11_run16_f","Anolis_k11_run17_f","Anolis_k11_run18_f","Anolis_k11_run19_f","Anolis_k11_run20_f",
          "Anolis_k12_run1_f","Anolis_k12_run2_f","Anolis_k12_run3_f","Anolis_k12_run4_f","Anolis_k12_run5_f","Anolis_k12_run6_f","Anolis_k12_run7_f",
          "Anolis_k12_run8_f","Anolis_k12_run9_f","Anolis_k12_run10_f","Anolis_k12_run11_f","Anolis_k12_run12_f","Anolis_k12_run13_f",
          "Anolis_k12_run14_f","Anolis_k12_run15_f","Anolis_k12_run16_f","Anolis_k12_run17_f","Anolis_k12_run18_f","Anolis_k12_run19_f","Anolis_k12_run20_f",
          "Anolis_k13_run1_f","Anolis_k13_run2_f","Anolis_k13_run3_f","Anolis_k13_run4_f","Anolis_k13_run5_f","Anolis_k13_run6_f","Anolis_k13_run7_f",
          "Anolis_k13_run8_f","Anolis_k13_run9_f","Anolis_k13_run10_f","Anolis_k13_run11_f","Anolis_k13_run12_f","Anolis_k13_run13_f",
          "Anolis_k13_run14_f","Anolis_k13_run15_f","Anolis_k13_run16_f","Anolis_k13_run17_f","Anolis_k13_run18_f","Anolis_k13_run19_f","Anolis_k13_run20_f",
          "Anolis_k14_run1_f","Anolis_k14_run2_f","Anolis_k14_run3_f","Anolis_k14_run4_f","Anolis_k14_run5_f","Anolis_k14_run6_f","Anolis_k14_run7_f",
          "Anolis_k14_run8_f","Anolis_k14_run9_f","Anolis_k14_run10_f","Anolis_k14_run11_f","Anolis_k14_run12_f","Anolis_k14_run13_f",
          "Anolis_k14_run14_f","Anolis_k14_run15_f","Anolis_k14_run16_f","Anolis_k14_run17_f","Anolis_k14_run18_f","Anolis_k14_run19_f","Anolis_k14_run20_f",
          "Anolis_k15_run1_f","Anolis_k15_run2_f","Anolis_k15_run3_f","Anolis_k15_run4_f","Anolis_k15_run5_f","Anolis_k15_run6_f","Anolis_k15_run7_f",
          "Anolis_k15_run8_f","Anolis_k15_run9_f","Anolis_k15_run10_f","Anolis_k15_run11_f","Anolis_k15_run12_f","Anolis_k15_run13_f",
          "Anolis_k15_run14_f","Anolis_k15_run15_f","Anolis_k15_run16_f","Anolis_k15_run17_f","Anolis_k15_run18_f","Anolis_k15_run19_f","Anolis_k15_run20_f",
          "Anolis_k16_run1_f","Anolis_k16_run2_f","Anolis_k16_run3_f","Anolis_k16_run4_f","Anolis_k16_run5_f","Anolis_k16_run6_f","Anolis_k16_run7_f",
          "Anolis_k16_run8_f","Anolis_k16_run9_f","Anolis_k16_run10_f","Anolis_k16_run11_f","Anolis_k16_run12_f","Anolis_k16_run13_f",
          "Anolis_k16_run14_f","Anolis_k16_run15_f","Anolis_k16_run16_f","Anolis_k16_run17_f","Anolis_k16_run18_f","Anolis_k16_run19_f","Anolis_k16_run20_f",
          "Anolis_k17_run1_f","Anolis_k17_run2_f","Anolis_k17_run3_f","Anolis_k17_run4_f","Anolis_k17_run5_f","Anolis_k17_run6_f","Anolis_k17_run7_f",
          "Anolis_k17_run8_f","Anolis_k17_run9_f","Anolis_k17_run10_f","Anolis_k17_run11_f","Anolis_k17_run12_f","Anolis_k17_run13_f",
          "Anolis_k17_run14_f","Anolis_k17_run15_f","Anolis_k17_run16_f","Anolis_k17_run17_f","Anolis_k17_run18_f","Anolis_k17_run19_f","Anolis_k17_run20_f",
          "Anolis_k18_run1_f","Anolis_k18_run2_f","Anolis_k18_run3_f","Anolis_k18_run4_f","Anolis_k18_run5_f","Anolis_k18_run6_f","Anolis_k18_run7_f",
          "Anolis_k18_run8_f","Anolis_k18_run9_f","Anolis_k18_run10_f","Anolis_k18_run11_f","Anolis_k18_run12_f","Anolis_k18_run13_f",
          "Anolis_k18_run14_f","Anolis_k18_run15_f","Anolis_k18_run16_f","Anolis_k18_run17_f","Anolis_k18_run18_f","Anolis_k18_run19_f","Anolis_k18_run20_f",
          "Anolis_k19_run1_f","Anolis_k19_run2_f","Anolis_k19_run3_f","Anolis_k19_run4_f","Anolis_k19_run5_f","Anolis_k19_run6_f","Anolis_k19_run7_f",
          "Anolis_k19_run8_f","Anolis_k19_run9_f","Anolis_k19_run10_f","Anolis_k19_run11_f","Anolis_k19_run12_f","Anolis_k19_run13_f",
          "Anolis_k19_run14_f","Anolis_k19_run15_f","Anolis_k19_run16_f","Anolis_k19_run17_f","Anolis_k19_run18_f","Anolis_k19_run19_f","Anolis_k19_run20_f",
          "Anolis_k20_run1_f","Anolis_k20_run2_f","Anolis_k20_run3_f","Anolis_k20_run4_f","Anolis_k20_run5_f","Anolis_k20_run6_f","Anolis_k20_run7_f",
          "Anolis_k20_run8_f","Anolis_k20_run9_f","Anolis_k20_run10_f","Anolis_k20_run11_f","Anolis_k20_run12_f","Anolis_k20_run13_f",
          "Anolis_k20_run14_f","Anolis_k20_run15_f","Anolis_k20_run16_f","Anolis_k20_run17_f","Anolis_k20_run18_f","Anolis_k20_run19_f","Anolis_k20_run20_f")

#convert STRUCTURE run files to qlist
slist<-readQ(sfiles)

#tabulate the qlist
tr1 <- tabulateQ(qlist=slist)

#make a summary table from the tabulated dataframe
sr1 <- summariseQ(tr1)

#apply the Evanno method for inferring K and export pdf
evannoMethodStructure(data=sr1,exportplot=T,returnplot=F,returndata=F,imgtype="pdf")


#####plotting STRUCTURE membership coefficients for K=2 and K=5
setwd("../membership_plots/") #this folder contains the structure results files to be used for plotting membership coefficients
#input files "Anolis_k2.ordered" and "Anolis_k5.ordered" were obtained using the "reorder_samples_STRUCTURE.sh" script

sfiles_K2<-c("Anolis_k2.ordered")
sfiles_K5<-c("Anolis_k5.ordered")

#convert STRUCTURE run files to qlist
slist_K2<-readQ(sfiles_K2)
slist_K5<-readQ(sfiles_K5)

#read in the metadata file, which contains two columns (individual ID, population ID)
labels<-read.csv(file="metadata.csv", header=FALSE)
labels$V2 <- as.character(labels$V2)

#keep only population IDs
onelabset <- labels[,2,drop=FALSE]

#plot results with ordered samples for K=2
plotQ(slist_K2[1],returnplot=F,exportplot=T,quiet=T,basesize=11,width=15,grplab = onelabset,grplabsize=1, imgtype="pdf")

#plot results with ordered samples for K=5
plotQ(slist_K5[1],returnplot=F,exportplot=T,quiet=T,basesize=11,width=15,grplab = onelabset,grplabsize=1, imgtype="pdf")
