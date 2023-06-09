library(readxl) # To run this R-script, the R-library readxl needs to be installed and functional

# Change this to wherever you decide to store ellipsoid_evaluation.xlsx
root_folder="/path/to/excel_file"

setwd(root_folder)

ellipsoids=as.data.frame(read_excel("ellipsoid_evaluation.xlsx"))


ellipsoids$Substrate=factor(ellipsoids$Substrate,levels=c("TCP","PDMS","imprinted"))



breaks=c(seq(from=1,to=3,by=1/3),50)


for_barplot_ecc = matrix(ncol=length(levels(ellipsoids$Substrate)),nrow=length(breaks)-1)


colnames(for_barplot_ecc)=levels(ellipsoids$Substrate)

for_chisq=for_barplot_ecc

for(theSubstrate in colnames(for_barplot_ecc))
{
    
    for_barplot_ecc[,theSubstrate]=hist(ellipsoids$Eccentricity[ellipsoids$Substrate==theSubstrate],breaks=breaks,plot=FALSE)$counts
    
    for_barplot_ecc[,theSubstrate]=for_barplot_ecc[,theSubstrate]/sum(for_barplot_ecc[,theSubstrate])
    
    for_chisq[,theSubstrate]=hist(ellipsoids$Eccentricity[ellipsoids$Substrate==theSubstrate],breaks=breaks,plot=FALSE)$counts
    
    for_barplot_ecc[,theSubstrate]=for_barplot_ecc[,theSubstrate]/sum(for_barplot_ecc[,theSubstrate])
    
    
    
    
    
}



names.arg= round((breaks[1:(length(breaks)-1)]+breaks[2:(length(breaks))])/2*10)/10


barplot(t(for_barplot_ecc),beside=TRUE,legend=TRUE, col=c("#eeeeeeff","#aaaaaaff","#555555ff"),names.arg=names.arg,xlab="Eccentricity",ylab="Frequency",las=2,ylim=c(0,0.3))

# Chisquare test of the distribution against the TCP
# Consider all categories except for the very large eccentricities where there are almost no cells

cat("Chi-square test for the distribution, PDMS vs TCP")

theSubstrate="PDMS"

ei_sub=sum(for_chisq[-dim(for_chisq)[1],theSubstrate])*(for_barplot_ecc[-dim(for_chisq)[1],"TCP"]+for_barplot_ecc[-dim(for_chisq)[1],theSubstrate])/2
ei_TCP=sum(for_chisq[-dim(for_chisq)[1],"TCP"])*(for_barplot_ecc[-dim(for_chisq)[1],"TCP"]+for_barplot_ecc[-dim(for_chisq)[1],theSubstrate])/2

ei=c(ei_sub,ei_TCP)

oi=c(for_chisq[-dim(for_chisq)[1],theSubstrate],for_chisq[-dim(for_chisq)[1],"TCP"])

chisq=sum((oi-ei)^2/ei)

# degrees of freedom: for each ei substrate (or TCP) one comme P-value is estimated, and we use the two totals

df=length(oi)-2- length(ei_sub)

3*pchisq(chisq,df,lower.tail=FALSE)

cat("Chi-square test for the distribution, imprinted vs TCP")

theSubstrate="imprinted"

ei_sub=sum(for_chisq[-dim(for_chisq)[1],theSubstrate])*(for_barplot_ecc[-dim(for_chisq)[1],"TCP"]+for_barplot_ecc[-dim(for_chisq)[1],theSubstrate])/2
ei_TCP=sum(for_chisq[-dim(for_chisq)[1],"TCP"])*(for_barplot_ecc[-dim(for_chisq)[1],"TCP"]+for_barplot_ecc[-dim(for_chisq)[1],theSubstrate])/2

ei=c(ei_sub,ei_TCP)

oi=c(for_chisq[-dim(for_chisq)[1],theSubstrate],for_chisq[-dim(for_chisq)[1],"TCP"])

chisq=sum((oi-ei)^2/ei)

# degrees of freedom: for each ei substrate (or TCP) one comme P-value is estimated, and we use the two totals

df=length(oi)-2- length(ei_sub)

3*pchisq(chisq,df,lower.tail=FALSE)

cat("Chi-square test for the distribution, imprinted vs PDMS")

theSubstrate="imprinted"

ei_sub=sum(for_chisq[-dim(for_chisq)[1],theSubstrate])*(for_barplot_ecc[-dim(for_chisq)[1],"PDMS"]+for_barplot_ecc[-dim(for_chisq)[1],theSubstrate])/2
ei_PDMS=sum(for_chisq[-dim(for_chisq)[1],"PDMS"])*(for_barplot_ecc[-dim(for_chisq)[1],"PDMS"]+for_barplot_ecc[-dim(for_chisq)[1],theSubstrate])/2

ei=c(ei_sub,ei_PDMS)

oi=c(for_chisq[-dim(for_chisq)[1],theSubstrate],for_chisq[-dim(for_chisq)[1],"PDMS"])

chisq=sum((oi-ei)^2/ei)

# degrees of freedom: for each ei substrate (or TCP) one comme P-value is estimated, and we use the two totals

df=length(oi)-2- length(ei_sub)

3*pchisq(chisq,df,lower.tail=FALSE)






