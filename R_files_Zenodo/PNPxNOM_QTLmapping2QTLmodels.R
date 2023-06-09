# PNPxNOM QTL mapping: analyse 2 QTL scans run on cluster 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# load libraries
library(qtl)
library(here)

# load scantwo output downloaded from cluster
load(file=as.character(paste("R_on_cluster/PNPxNOMoutscantwo1", sep="/")))
load(file=as.character(paste("R_on_cluster/PNPxNOMoutscantwo2", sep="/")))
load(file=as.character(paste("R_on_cluster/PNPxNOMoutscantwo3", sep="/")))
load(file=as.character(paste("R_on_cluster/PNPxNOMoutscantwo4", sep="/")))
load(file=as.character(paste("R_on_cluster/PNPxNOMoutscantwo5", sep="/")))
load(file=as.character(paste("R_on_cluster/PNPxNOMoutscantwo6", sep="/")))
load(file=as.character(paste("R_on_cluster/PNPxNOMoutscantwo7", sep="/")))
load(file=as.character(paste("R_on_cluster/PNPxNOMoutscantwo8", sep="/")))
# (plots are not very helpful here, summaries are easier to understand)

# load scantwo permutations downloaded from cluster
# run in four 250 permutation batches for each trait

# traits with covars sex and fam
permsHL1<-readRDS("R_on_cluster/permsHL1")
permsHL2<-readRDS("R_on_cluster/permsHL2")
permsHL3<-readRDS("R_on_cluster/permsHL3")
permsHL4<-readRDS("R_on_cluster/permsHL4")
# combine the 4 permutation runs:
permsHL1all<-c(permsHL1, permsHL2, permsHL3, permsHL4)
# generate summary
HLsummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 1, perms=permsHL1all, alpha=0.1, pvalues=TRUE) 

permsHW1<-readRDS("R_on_cluster/permsHW1")
permsHW2<-readRDS("R_on_cluster/permsHW2")
permsHW3<-readRDS("R_on_cluster/permsHW3")
permsHW4<-readRDS("R_on_cluster/permsHW4")
permsHW1all<-c(permsHW1, permsHW2, permsHW3, permsHW4)
HWsummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 2, perms=permsHW1all, alpha=0.1, pvalues=TRUE) 

permsLJL1<-readRDS("R_on_cluster/permsLJL1")
permsLJL2<-readRDS("R_on_cluster/permsLJL2")
permsLJL3<-readRDS("R_on_cluster/permsLJL3")
permsLJL4<-readRDS("R_on_cluster/permsLJL4")
permsLJL1all<-c(permsLJL1, permsLJL2, permsLJL3, permsLJL4)
LJLsummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 3, perms=permsLJL1all, alpha=0.1, pvalues=TRUE) 

permsLJW1<-readRDS("R_on_cluster/permsLJW1")
permsLJW2<-readRDS("R_on_cluster/permsLJW2")
permsLJW3<-readRDS("R_on_cluster/permsLJW3")
permsLJW4<-readRDS("R_on_cluster/permsLJW4")
permsLJW1all<-c(permsLJW1, permsLJW2, permsLJW3, permsLJW4)
LJWsummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 4, perms=permsLJW1all, alpha=0.1, pvalues=TRUE) 

permsSnL1<-readRDS("R_on_cluster/permsSnL1")
permsSnL2<-readRDS("R_on_cluster/permsSnL2")
permsSnL3<-readRDS("R_on_cluster/permsSnL3")
permsSnL4<-readRDS("R_on_cluster/permsSnL4")
permsSnL1all<-c(permsSnL1, permsSnL2, permsSnL3, permsSnL4)
SnLsummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 5, perms=permsSnL1all, alpha=0.1, pvalues=TRUE) 

permsSnW1<-readRDS("R_on_cluster/permsSnW1")
permsSnW2<-readRDS("R_on_cluster/permsSnW2")
permsSnW3<-readRDS("R_on_cluster/permsSnW3")
permsSnW4<-readRDS("R_on_cluster/permsSnW4")
permsSnW1all<-c(permsSnW1, permsSnW2, permsSnW3, permsSnW4)
SnWsummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 6, perms=permsSnW1all, alpha=0.1, pvalues=TRUE) 

permsEyL1<-readRDS("R_on_cluster/permsEyL1")
permsEyL2<-readRDS("R_on_cluster/permsEyL2")
permsEyL3<-readRDS("R_on_cluster/permsEyL3")
permsEyL4<-readRDS("R_on_cluster/permsEyL4")
permsEyL1all<-c(permsEyL1, permsEyL2, permsEyL3, permsEyL4)
EyLsummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 7, perms=permsEyL1all, alpha=0.1, pvalues=TRUE) 

permsEyD1<-readRDS("R_on_cluster/permsEyD1")
permsEyD2<-readRDS("R_on_cluster/permsEyD2")
permsEyD3<-readRDS("R_on_cluster/permsEyD3")
permsEyD4<-readRDS("R_on_cluster/permsEyD4")
permsEyD1all<-c(permsEyD1, permsEyD2, permsEyD3, permsEyD4)
EyDsummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 8, perms=permsEyD1all, alpha=0.1, pvalues=TRUE) 

permsPOW1<-readRDS("R_on_cluster/permsPOW1")
permsPOW2<-readRDS("R_on_cluster/permsPOW2")
permsPOW3<-readRDS("R_on_cluster/permsPOW3")
permsPOW4<-readRDS("R_on_cluster/permsPOW4")
permsPOW1all<-c(permsPOW1, permsPOW2, permsPOW3, permsPOW4)
POWsummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 9, perms=permsPOW1all, alpha=0.1, pvalues=TRUE) 

permsInnerRows1<-readRDS("R_on_cluster/permsInnerRows1")
permsInnerRows2<-readRDS("R_on_cluster/permsInnerRows2")
permsInnerRows3<-readRDS("R_on_cluster/permsInnerRows3")
permsInnerRows4<-readRDS("R_on_cluster/permsInnerRows4")
permsInnerRows1all<-c(permsInnerRows1, permsInnerRows2, permsInnerRows3, permsInnerRows4)
InnerRowssummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 10, perms=permsInnerRows1all, alpha=0.1, pvalues=TRUE) 

permsGap1<-readRDS("R_on_cluster/permsGap1")
permsGap2<-readRDS("R_on_cluster/permsGap2")
permsGap3<-readRDS("R_on_cluster/permsGap3")
permsGap4<-readRDS("R_on_cluster/permsGap4")
permsGap1all<-c(permsGap1, permsGap2, permsGap3, permsGap4)
Gapsummary<-summary(PNPxNOMoutscantwo1, lodcolumn = 11, perms=permsGap1all, alpha=0.1, pvalues=TRUE) 

# traits with covar sex
permsBD1<-readRDS("R_on_cluster/permsBD1")
permsBD2<-readRDS("R_on_cluster/permsBD2")
permsBD3<-readRDS("R_on_cluster/permsBD3")
permsBD4<-readRDS("R_on_cluster/permsBD4")
permsBD1all<-c(permsBD1, permsBD2, permsBD3, permsBD4)
BDsummary<-summary(PNPxNOMoutscantwo2, lodcolumn = 1, perms=permsBD1all, alpha=0.1, pvalues=TRUE) 

permsChD1<-readRDS("R_on_cluster/permsChD1")
permsChD2<-readRDS("R_on_cluster/permsChD2")
permsChD3<-readRDS("R_on_cluster/permsChD3")
permsChD4<-readRDS("R_on_cluster/permsChD4")
permsChD1all<-c(permsChD1, permsChD2, permsChD3, permsChD4)
ChDsummary<-summary(PNPxNOMoutscantwo2, lodcolumn = 2, perms=permsChD1all, alpha=0.1, pvalues=TRUE) 

permsPOD1<-readRDS("R_on_cluster/permsPOD1")
permsPOD2<-readRDS("R_on_cluster/permsPOD2")
permsPOD3<-readRDS("R_on_cluster/permsPOD3")
permsPOD4<-readRDS("R_on_cluster/permsPOD4")
permsPOD1all<-c(permsPOD1, permsPOD2, permsPOD3, permsPOD4)
PODsummary<-summary(PNPxNOMoutscantwo2, lodcolumn = 3, perms=permsPOD1all, alpha=0.1, pvalues=TRUE) 

permsIOW1<-readRDS("R_on_cluster/permsIOW1")
permsIOW2<-readRDS("R_on_cluster/permsIOW2")
permsIOW3<-readRDS("R_on_cluster/permsIOW3")
permsIOW4<-readRDS("R_on_cluster/permsIOW4")
permsIOW1all<-c(permsIOW1, permsIOW2, permsIOW3, permsIOW4)
IOWsummary<-summary(PNPxNOMoutscantwo2, lodcolumn = 4, perms=permsIOW1all, alpha=0.1, pvalues=TRUE) 


# traits with covar fam
permsCPL1<-readRDS("R_on_cluster/permsCPL1")
permsCPL2<-readRDS("R_on_cluster/permsCPL2")
permsCPL3<-readRDS("R_on_cluster/permsCPL3")
permsCPL4<-readRDS("R_on_cluster/permsCPL4")
permsCPL1all<-c(permsCPL1, permsCPL2, permsCPL3, permsCPL4)
CPLsummary<-summary(PNPxNOMoutscantwo3, lodcolumn = 1, perms=permsCPL1all, alpha=0.1, pvalues=TRUE) 

permsCPD1<-readRDS("R_on_cluster/permsCPD1")
permsCPD2<-readRDS("R_on_cluster/permsCPD2")
permsCPD3<-readRDS("R_on_cluster/permsCPD3")
permsCPD4<-readRDS("R_on_cluster/permsCPD4")
permsCPD1all<-c(permsCPD1, permsCPD2, permsCPD3, permsCPD4)
CPDsummary<-summary(PNPxNOMoutscantwo3, lodcolumn = 2, perms=permsCPD1all, alpha=0.1, pvalues=TRUE) 

permsgutlength1<-readRDS("R_on_cluster/permsgutlength1")
permsgutlength2<-readRDS("R_on_cluster/permsgutlength2")
permsgutlength3<-readRDS("R_on_cluster/permsgutlength3")
permsgutlength4<-readRDS("R_on_cluster/permsgutlength4")
permsgutlength1all<-c(permsgutlength1, permsgutlength2, permsgutlength3, permsgutlength4)
gutlengthsummary<-summary(PNPxNOMoutscantwo3, lodcolumn = 3, perms=permsgutlength1all, alpha=0.1, pvalues=TRUE) 

permstoothDensityL1<-readRDS("R_on_cluster/permstoothDensityL1")
permstoothDensityL2<-readRDS("R_on_cluster/permstoothDensityL2")
permstoothDensityL3<-readRDS("R_on_cluster/permstoothDensityL3")
permstoothDensityL4<-readRDS("R_on_cluster/permstoothDensityL4")
permstoothDensityL1all<-c(permstoothDensityL1, permstoothDensityL2, permstoothDensityL3, permstoothDensityL4)
toothDensityLsummary<-summary(PNPxNOMoutscantwo3, lodcolumn = 4, perms=permstoothDensityL1all, alpha=0.1, pvalues=TRUE) 

# traits with no covars
permstoothshapeL1<-readRDS("R_on_cluster/permstoothshapeL1")
permstoothshapeL2<-readRDS("R_on_cluster/permstoothshapeL2")
permstoothshapeL3<-readRDS("R_on_cluster/permstoothshapeL3")
permstoothshapeL4<-readRDS("R_on_cluster/permstoothshapeL4")
permstoothshapeL1all<-c(permstoothshapeL1, permstoothshapeL2, permstoothshapeL3, permstoothshapeL4)
toothshapeLsummary<-summary(PNPxNOMoutscantwo4, lodcolumn = 1, perms=permstoothshapeL1all, alpha=0.1, pvalues=TRUE) 

permstoothShapeU1<-readRDS("R_on_cluster/permstoothShapeU1")
permstoothShapeU2<-readRDS("R_on_cluster/permstoothShapeU2")
permstoothShapeU3<-readRDS("R_on_cluster/permstoothShapeU3")
permstoothShapeU4<-readRDS("R_on_cluster/permstoothShapeU4")
permstoothShapeU1all<-c(permstoothShapeU1, permstoothShapeU2, permstoothShapeU3, permstoothShapeU4)
toothShapeUsummary<-summary(PNPxNOMoutscantwo4, lodcolumn = 2, perms=permstoothShapeU1all, alpha=0.1, pvalues=TRUE) 

permstoothDensityU1<-readRDS("R_on_cluster/permstoothDensityU1")
permstoothDensityU2<-readRDS("R_on_cluster/permstoothDensityU2")
permstoothDensityU3<-readRDS("R_on_cluster/permstoothDensityU3")
permstoothDensityU4<-readRDS("R_on_cluster/permstoothDensityU4")
permstoothDensityU1all<-c(permstoothDensityU1, permstoothDensityU2, permstoothDensityU3, permstoothDensityU4)
toothDensityUsummary<-summary(PNPxNOMoutscantwo4, lodcolumn = 3, perms=permstoothDensityU1all, alpha=0.1, pvalues=TRUE) 

# transgressive
permstricuspid1<-readRDS("R_on_cluster/permstricuspid1")
permstricuspid2<-readRDS("R_on_cluster/permstricuspid2")
permstricuspid3<-readRDS("R_on_cluster/permstricuspid3")
permstricuspid4<-readRDS("R_on_cluster/permstricuspid4")
permstricuspid1all<-c(permstricuspid1, permstricuspid2, permstricuspid3, permstricuspid4)
tricuspidsummary<-summary(PNPxNOMoutscantwo5, lodcolumn = 1, perms=permstricuspid1all, alpha=0.1, pvalues=TRUE) 

# colour males, covar fam
permsRdorsalfin11<-readRDS("R_on_cluster/permsRdorsalfin11")
permsRdorsalfin12<-readRDS("R_on_cluster/permsRdorsalfin12")
permsRdorsalfin13<-readRDS("R_on_cluster/permsRdorsalfin13")
permsRdorsalfin14<-readRDS("R_on_cluster/permsRdorsalfin14")
permsRdorsalfin11all<-c(permsRdorsalfin11, permsRdorsalfin12, permsRdorsalfin13, permsRdorsalfin14)
Rdorsalfin1summary<-summary(PNPxNOMoutscantwo6, lodcolumn = 1, perms=permsRdorsalfin11all, alpha=0.1, pvalues=TRUE) 

permsRdorsalfin21<-readRDS("R_on_cluster/permsRdorsalfin21")
permsRdorsalfin22<-readRDS("R_on_cluster/permsRdorsalfin22")
permsRdorsalfin23<-readRDS("R_on_cluster/permsRdorsalfin23")
permsRdorsalfin24<-readRDS("R_on_cluster/permsRdorsalfin24")
permsRdorsalfin21all<-c(permsRdorsalfin21, permsRdorsalfin22, permsRdorsalfin23, permsRdorsalfin24)
Rdorsalfin2summary<-summary(PNPxNOMoutscantwo6, lodcolumn = 2, perms=permsRdorsalfin21all, alpha=0.1, pvalues=TRUE) 

permsRhead1<-readRDS("R_on_cluster/permsRhead1")
permsRhead2<-readRDS("R_on_cluster/permsRhead2")
permsRhead3<-readRDS("R_on_cluster/permsRhead3")
permsRhead4<-readRDS("R_on_cluster/permsRhead4")
permsRhead1all<-c(permsRhead1, permsRhead2, permsRhead3, permsRhead4)
Rheadsummary<-summary(PNPxNOMoutscantwo6, lodcolumn = 3, perms=permsRhead1all, alpha=0.1, pvalues=TRUE) 

permsRdorsum11<-readRDS("R_on_cluster/permsRdorsum11")
permsRdorsum12<-readRDS("R_on_cluster/permsRdorsum12")
permsRdorsum13<-readRDS("R_on_cluster/permsRdorsum13")
permsRdorsum14<-readRDS("R_on_cluster/permsRdorsum14")
permsRdorsum11all<-c(permsRdorsum11, permsRdorsum12, permsRdorsum13, permsRdorsum14)
Rdorsum1summary<-summary(PNPxNOMoutscantwo6, lodcolumn = 4, perms=permsRdorsum11all, alpha=0.1, pvalues=TRUE) 

permsRdorsum21<-readRDS("R_on_cluster/permsRdorsum21")
permsRdorsum22<-readRDS("R_on_cluster/permsRdorsum22")
permsRdorsum23<-readRDS("R_on_cluster/permsRdorsum23")
permsRdorsum24<-readRDS("R_on_cluster/permsRdorsum24")
permsRdorsum21all<-c(permsRdorsum21, permsRdorsum22, permsRdorsum23, permsRdorsum24)
Rdorsum2summary<-summary(PNPxNOMoutscantwo6, lodcolumn = 5, perms=permsRdorsum21all, alpha=0.1, pvalues=TRUE) 

permsYflanks1<-readRDS("R_on_cluster/permsYflanks1")
permsYflanks2<-readRDS("R_on_cluster/permsYflanks2")
permsYflanks3<-readRDS("R_on_cluster/permsYflanks3")
permsYflanks4<-readRDS("R_on_cluster/permsYflanks4")
permsYflanks1all<-c(permsYflanks1, permsYflanks2, permsYflanks3, permsYflanks4)
Yflankssummary<-summary(PNPxNOMoutscantwo6, lodcolumn = 6, perms=permsYflanks1all, alpha=0.1, pvalues=TRUE) 

permsYgillcover1<-readRDS("R_on_cluster/permsYgillcover1")
permsYgillcover2<-readRDS("R_on_cluster/permsYgillcover2")
permsYgillcover3<-readRDS("R_on_cluster/permsYgillcover3")
permsYgillcover4<-readRDS("R_on_cluster/permsYgillcover4")
permsYgillcover1all<-c(permsYgillcover1, permsYgillcover2, permsYgillcover3, permsYgillcover4)
Ygillcoversummary<-summary(PNPxNOMoutscantwo6, lodcolumn = 7, perms=permsYgillcover1all, alpha=0.1, pvalues=TRUE) 

# colour males, no covar 
permsRflanks11<-readRDS("R_on_cluster/permsRflanks11")
permsRflanks12<-readRDS("R_on_cluster/permsRflanks12")
permsRflanks13<-readRDS("R_on_cluster/permsRflanks13")
permsRflanks14<-readRDS("R_on_cluster/permsRflanks14")
permsRflanks11all<-c(permsRflanks11, permsRflanks12, permsRflanks13, permsRflanks14)
Rflanks1summary<-summary(PNPxNOMoutscantwo7, lodcolumn = 1, perms=permsRflanks11all, alpha=0.1, pvalues=TRUE) 

permsRcheek1<-readRDS("R_on_cluster/permsRcheek1")
permsRcheek2<-readRDS("R_on_cluster/permsRcheek2")
permsRcheek3<-readRDS("R_on_cluster/permsRcheek3")
permsRcheek4<-readRDS("R_on_cluster/permsRcheek4")
permsRcheek1all<-c(permsRcheek1, permsRcheek2, permsRcheek3, permsRcheek4)
Rcheeksummary<-summary(PNPxNOMoutscantwo7, lodcolumn = 2, perms=permsRcheek1all, alpha=0.1, pvalues=TRUE) 

permsRgillcover1<-readRDS("R_on_cluster/permsRgillcover1")
permsRgillcover2<-readRDS("R_on_cluster/permsRgillcover2")
permsRgillcover3<-readRDS("R_on_cluster/permsRgillcover3")
permsRgillcover4<-readRDS("R_on_cluster/permsRgillcover4")
permsRgillcover1all<-c(permsRgillcover1, permsRgillcover2, permsRgillcover3, permsRgillcover4)
Rgillcoversummary<-summary(PNPxNOMoutscantwo7, lodcolumn = 3, perms=permsRgillcover1all, alpha=0.1, pvalues=TRUE) 

permsRnose1<-readRDS("R_on_cluster/permsRnose1")
permsRnose2<-readRDS("R_on_cluster/permsRnose2")
permsRnose3<-readRDS("R_on_cluster/permsRnose3")
permsRnose4<-readRDS("R_on_cluster/permsRnose4")
permsRnose1all<-c(permsRnose1, permsRnose2, permsRnose3, permsRnose4)
Rnosesummary<-summary(PNPxNOMoutscantwo7, lodcolumn = 4, perms=permsRnose1all, alpha=0.1, pvalues=TRUE) 

permsYupperlip1<-readRDS("R_on_cluster/permsYupperlip1")
permsYupperlip2<-readRDS("R_on_cluster/permsYupperlip2")
permsYupperlip3<-readRDS("R_on_cluster/permsYupperlip3")
permsYupperlip4<-readRDS("R_on_cluster/permsYupperlip4")
permsYupperlip1all<-c(permsYupperlip1, permsYupperlip2, permsYupperlip3, permsYupperlip4)
Yupperlipsummary<-summary(PNPxNOMoutscantwo7, lodcolumn = 5, perms=permsYupperlip1all, alpha=0.1, pvalues=TRUE) 

permsYcheek1<-readRDS("R_on_cluster/permsYcheek1")
permsYcheek2<-readRDS("R_on_cluster/permsYcheek2")
permsYcheek3<-readRDS("R_on_cluster/permsYcheek3")
permsYcheek4<-readRDS("R_on_cluster/permsYcheek4")
permsYcheek1all<-c(permsYcheek1, permsYcheek2, permsYcheek3, permsYcheek4)
Ycheeksummary<-summary(PNPxNOMoutscantwo7, lodcolumn = 6, perms=permsYcheek1all, alpha=0.1, pvalues=TRUE) 

permsYpelvicfin1<-readRDS("R_on_cluster/permsYpelvicfin1")
permsYpelvicfin2<-readRDS("R_on_cluster/permsYpelvicfin2")
permsYpelvicfin3<-readRDS("R_on_cluster/permsYpelvicfin3")
permsYpelvicfin4<-readRDS("R_on_cluster/permsYpelvicfin4")
permsYpelvicfin1all<-c(permsYpelvicfin1, permsYpelvicfin2, permsYpelvicfin3, permsYpelvicfin4)
Ypelvicfinsummary<-summary(PNPxNOMoutscantwo7, lodcolumn = 7, perms=permsYpelvicfin1all, alpha=0.1, pvalues=TRUE) 

permsYnose1<-readRDS("R_on_cluster/permsYnose1")
permsYnose2<-readRDS("R_on_cluster/permsYnose2")
permsYnose3<-readRDS("R_on_cluster/permsYnose3")
permsYnose4<-readRDS("R_on_cluster/permsYnose4")
permsYnose1all<-c(permsYnose1, permsYnose2, permsYnose3, permsYnose4)
Ynosesummary<-summary(PNPxNOMoutscantwo7, lodcolumn = 8, perms=permsYnose1all, alpha=0.1, pvalues=TRUE) 

# stripes with addcovar family (normal)
permsNo_stripes1<-readRDS("R_on_cluster/permsNo_stripes1")
permsNo_stripes2<-readRDS("R_on_cluster/permsNo_stripes2")
permsNo_stripes3<-readRDS("R_on_cluster/permsNo_stripes3")
permsNo_stripes4<-readRDS("R_on_cluster/permsNo_stripes4")
permsNo_stripes1all<-c(permsNo_stripes1, permsNo_stripes2, permsNo_stripes3, permsNo_stripes4)
No_stripessummary<-summary(PNPxNOMoutscantwo8, lodcolumn = 1, perms=permsNo_stripes1all, alpha=0.1, pvalues=TRUE) 


# combine all summaries
PNPxNOMsummary2QTL<-rbind(HLsummary, HWsummary, LJLsummary, LJWsummary, SnLsummary, SnWsummary, EyLsummary, EyDsummary, POWsummary, InnerRowssummary, Gapsummary, 
                          BDsummary, ChDsummary, PODsummary, IOWsummary, 
                          CPLsummary, CPDsummary, gutlengthsummary, toothDensityLsummary, 
                          toothshapeLsummary, toothShapeUsummary, toothDensityUsummary, 
                          tricuspidsummary, 
                          Rdorsalfin1summary, Rdorsalfin2summary, Rheadsummary, Rdorsum1summary, Rdorsum2summary, Yflankssummary, Ygillcoversummary, 
                          Rflanks1summary, Rcheeksummary, Rgillcoversummary, Rnosesummary, Yupperlipsummary, Ycheeksummary, Ypelvicfinsummary, Ynosesummary, 
                          No_stripessummary)

# add names in new column with trait numbers
PNPxNOMsummary2QTL["trait"]<-rep(c("HL","HW","LJL","LJW","SnL","SnW","EyL","EyD","POW","InnerRows","Gap",
                                    "BD","ChD","POD","IOW",
                                    "CPL","CPD","gutlength","toothDensityL",
                                    "toothshapeL","toothShapeU","toothDensityU",
                                    "tricuspid",
                                    "Rdorsalfin1","Rdorsalfin2","Rhead","Rdorsum1","Rdorsum2","Yflanks","Ygillcover",
                                    "Rflanks1","Rcheek","Rgillcover","Rnose","Yupperlip","Ycheek","Ypelvicfin","Ynose",
                                    "No_stripes"), 
                                 sapply(list(HLsummary, HWsummary, LJLsummary, LJWsummary, SnLsummary, SnWsummary, EyLsummary, EyDsummary, POWsummary, InnerRowssummary, Gapsummary, 
                                             BDsummary, ChDsummary, PODsummary, IOWsummary, 
                                             CPLsummary, CPDsummary, gutlengthsummary, toothDensityLsummary, 
                                             toothshapeLsummary, toothShapeUsummary, toothDensityUsummary, 
                                             tricuspidsummary, 
                                             Rdorsalfin1summary, Rdorsalfin2summary, Rheadsummary, Rdorsum1summary, Rdorsum2summary, Yflankssummary, Ygillcoversummary, 
                                             Rflanks1summary, Rcheeksummary, Rgillcoversummary, Rnosesummary, Yupperlipsummary, Ycheeksummary, Ypelvicfinsummary, Ynosesummary, 
                                             No_stripessummary), nrow))
# write out as csv
write.csv(x = PNPxNOMsummary2QTL, file = "Output/PNPxNOMsummary2QTL_summary.csv")
