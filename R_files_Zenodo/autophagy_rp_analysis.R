## requirements for analysis of ribosome profiling

library(edgeR)
library(gplots)
source("functions.R"); # Load all the of the code and functions in Babel 
library(gridExtra)
library(ggplot2)
library(devtools)
library(easyGgplot2)
library(dplyr)

#open the counts file
cov=read.table("autophagy_rp.orf.cov", header=T, sep ="\\t")

#open the babel output
load("all_groups.between.babel", .GlobalEnv)
load("all_groups.within.babel", .GlobalEnv)

#match gene names to gene id
	#gene names are in the counts file cov cov$symbol
	#all of the gene ids in the list between.babel are the same
	all(between.babel[[1]][,1]==between.babel[[2]][,1])
	#match
	gene_names=cov[,1:2]
	gene_names=gene_names[gene_names$gene_id %in% between.babel[[1]][,1],]
	
	for (i in 1:6)
	{
		between.babel[[i]]=cbind(gene_names$symbol, between.babel[[i]])
		
	}
	
# NOTE: in the datasets
### 1 = Atg12f/f in Control media
### 2 = Atg12f/f in HBSS
### 3 = Atg12KO in Control media
### 4 = Atg12KO in HBSS


#all the data
RP2v4 = between.babel[[5]]
RP1v3 = between.babel[[2]]
RP1v2 = between.babel[[1]]

#rename the data to better sort
names=c("genesymbol", "gene", "mRNA_LogFC", "mRNA_FDR", "Change_type", "Pval", "FDR", "Direction")

colnames(RP1v2) = names
colnames(RP1v3) = names
colnames(RP2v4) = names



##############################################################
####### Figure S2C ##########
##############################################################
#Pearson correlation between replicates from cov file

cor1v2 = cor.test(cov$control_fed_1.rna, cov$control_fed_2.rna)
cor1v3 = cor.test(cov$control_fed_1.rna, cov$control_fed_3.rna)
cor1v4 = cor.test(cov$control_fed_1.rna, cov$control_fed_4.rna)
cor2v3 = cor.test(cov$control_fed_2.rna, cov$control_fed_3.rna)
cor2v4 = cor.test(cov$control_fed_2.rna, cov$control_fed_4.rna)
cor3v4 = cor.test(cov$control_fed_3.rna, cov$control_fed_4.rna)

cor_control_fed.rna = c(cor1v2$p.value, cor1v3$p.value, cor1v4$p.value, cor2v3$p.value, cor2v4$p.value, cor3v4$p.value)
Cor_cor_control_fed.rna = c(cor1v2$estimate, cor1v3$estimate, cor1v4$estimate, cor2v3$estimate, cor2v4$estimate, cor3v4$estimate)


a = cor.test(cov$control_starved_1.rna, cov$control_starved_2.rna)
b = cor.test(cov$control_starved_1.rna, cov$control_starved_3.rna)
c = cor.test(cov$control_starved_1.rna, cov$control_starved_4.rna)
d= cor.test(cov$control_starved_2.rna, cov$control_starved_3.rna)
e=cor.test(cov$control_starved_2.rna, cov$control_starved_4.rna)
f=cor.test(cov$control_starved_3.rna, cov$control_starved_4.rna)

cor_control_starved.rna = c(a$p.value, b$p.value, c$p.value, d$p.value, e$p.value, f$p.value)
Cor_cor_control_starved.rna = c(a$estimate, b$estimate, c$estimate, d$estimate, e$estimate, f$estimate)


g=cor.test(cov$atg_minus_fed_1.rna, cov$atg_minus_fed_2.rna)
h=cor.test(cov$atg_minus_fed_1.rna, cov$atg_minus_fed_3.rna)
i=cor.test(cov$atg_minus_fed_1.rna, cov$atg_minus_fed_4.rna)
j=cor.test(cov$atg_minus_fed_2.rna, cov$atg_minus_fed_3.rna)
k=cor.test(cov$atg_minus_fed_2.rna, cov$atg_minus_fed_4.rna)
l=cor.test(cov$atg_minus_fed_3.rna, cov$atg_minus_fed_4.rna)

cor_atg_minus_fed.rna = c(g$p.value, h$p.value, i$p.value, j$p.value, k$p.value, l$p.value)
Cor_cor_atg_minus_fed.rna = c(g$estimate, h$estimate, i$estimate, j$estimate, k$estimate, l$estimate)


m=cor.test(cov$atg_minus_starved_1.rna, cov$atg_minus_starved_2.rna)
n=cor.test(cov$atg_minus_starved_1.rna, cov$atg_minus_starved_3.rna)
o=cor.test(cov$atg_minus_starved_1.rna, cov$atg_minus_starved_4.rna)
p=cor.test(cov$atg_minus_starved_2.rna, cov$atg_minus_starved_3.rna)
q=cor.test(cov$atg_minus_starved_2.rna, cov$atg_minus_starved_4.rna)
r=cor.test(cov$atg_minus_starved_3.rna, cov$atg_minus_starved_4.rna)

cor_atg_minus_starved.rna = c(m$p.value, n$p.value, o$p.value, p$p.value, q$p.value, r$p.value)
Cor_cor_atg_minus_starved.rna = c(m$estimate, n$estimate, o$estimate, p$estimate, q$estimate, r$estimate)

cor1v2=cor.test(cov$control_fed_1.rp, cov$control_fed_2.rp)
cor1v3=cor.test(cov$control_fed_1.rp, cov$control_fed_3.rp)
cor1v4=cor.test(cov$control_fed_1.rp, cov$control_fed_4.rp)
cor2v3=cor.test(cov$control_fed_2.rp, cov$control_fed_3.rp)
cor2v4=cor.test(cov$control_fed_2.rp, cov$control_fed_4.rp)
cor3v4=cor.test(cov$control_fed_3.rp, cov$control_fed_4.rp)

cor_control_fed.rp = c(cor1v2$p.value, cor1v3$p.value, cor1v4$p.value, cor2v3$p.value, cor2v4$p.value, cor3v4$p.value)
Cor_cor_control_fed.rp = c(cor1v2$estimate, cor1v3$estimate, cor1v4$estimate, cor2v3$estimate, cor2v4$estimate, cor3v4$estimate)

a=cor.test(cov$control_starved_1.rp, cov$control_starved_2.rp)
b=cor.test(cov$control_starved_1.rp, cov$control_starved_3.rp)
c=cor.test(cov$control_starved_1.rp, cov$control_starved_4.rp)
d=cor.test(cov$control_starved_2.rp, cov$control_starved_3.rp)
e=cor.test(cov$control_starved_2.rp, cov$control_starved_4.rp)
f=cor.test(cov$control_starved_3.rp, cov$control_starved_4.rp)

cor_control_starved.rp = c(a$p.value, b$p.value, c$p.value, d$p.value, e$p.value, f$p.value)
Cor_cor_control_starved.rp = c(a$estimate, b$estimate, c$estimate, d$estimate, e$estimate, f$estimate)


g=cor.test(cov$atg_minus_fed_1.rp, cov$atg_minus_fed_2.rp)
h=cor.test(cov$atg_minus_fed_1.rp, cov$atg_minus_fed_3.rp)
i=cor.test(cov$atg_minus_fed_1.rp, cov$atg_minus_fed_4.rp)
j=cor.test(cov$atg_minus_fed_2.rp, cov$atg_minus_fed_3.rp)
k=cor.test(cov$atg_minus_fed_2.rp, cov$atg_minus_fed_4.rp)
l=cor.test(cov$atg_minus_fed_3.rp, cov$atg_minus_fed_4.rp)

cor_atg_minus_fed.rp = c(g$p.value, h$p.value, i$p.value, j$p.value, k$p.value, l$p.value)
Cor_cor_atg_minus_fed.rp = c(g$estimate, h$estimate, i$estimate, j$estimate, k$estimate, l$estimate)

m=cor.test(cov$atg_minus_starved_1.rp, cov$atg_minus_starved_2.rp)
n=cor.test(cov$atg_minus_starved_1.rp, cov$atg_minus_starved_3.rp)
o=cor.test(cov$atg_minus_starved_1.rp, cov$atg_minus_starved_4.rp)
p=cor.test(cov$atg_minus_starved_2.rp, cov$atg_minus_starved_3.rp)
q=cor.test(cov$atg_minus_starved_2.rp, cov$atg_minus_starved_4.rp)
r=cor.test(cov$atg_minus_starved_3.rp, cov$atg_minus_starved_4.rp)

cor_atg_minus_starved.rp = c(m$p.value, n$p.value, o$p.value, p$p.value, q$p.value, r$p.value)
Cor_cor_atg_minus_starved.rp = c(m$estimate, n$estimate, o$estimate, p$estimate, q$estimate, r$estimate)

rownamescor= c("rep1v2", "rep1v3", "rep1v4", "rep2v3", "rep2v4", "rep3v4")

cor_pval_rna = cbind(cor_control_fed.rna, cor_control_starved.rna, cor_atg_minus_fed.rna, cor_atg_minus_starved.rna)
rownames(cor_pval_rna) = rownamescor

cor_pval_rp = cbind(cor_control_fed.rp, cor_control_starved.rp, cor_atg_minus_fed.rp, cor_atg_minus_starved.rp)
rownames(cor_pval_rp) = rownamescor

cor_rna = cbind(Cor_cor_control_fed.rna, Cor_cor_control_starved.rna, Cor_cor_atg_minus_fed.rna, Cor_cor_atg_minus_starved.rna)
rownames(cor_rna) = rownamescor

cor_rp = cbind(Cor_cor_control_fed.rp, Cor_cor_control_starved.rp, Cor_cor_atg_minus_fed.rp, Cor_cor_atg_minus_starved.rp)
rownames(cor_rp) = rownamescor

##############################################################
#Pearson correlation of pval of RP occupancy between replicates from within file

cor1v2 = cor.test(within.babel[[1]][,3], within.babel[[2]][,3])
cor1v3 = cor.test(within.babel[[1]][,3], within.babel[[3]][,3])
cor1v4 = cor.test(within.babel[[1]][,3], within.babel[[4]][,3])
cor2v3 = cor.test(within.babel[[2]][,3], within.babel[[3]][,3])
cor2v4 = cor.test(within.babel[[2]][,3], within.babel[[4]][,3])
cor3v4 = cor.test(within.babel[[3]][,3], within.babel[[4]][,3])

pearson_RPoccupancy_control_fed_pval = c(cor1v2$p.value, cor1v3$p.value, cor1v4$p.value, cor2v3$p.value, cor2v4$p.value, cor3v4$p.value)
pearson_RPoccupancy_control_fed= c(cor1v2$estimate, cor1v3$estimate, cor1v4$estimate, cor2v3$estimate, cor2v4$estimate, cor3v4$estimate)


a = cor.test(within.babel[[5]][,3], within.babel[[6]][,3])
b = cor.test(within.babel[[5]][,3], within.babel[[7]][,3])
c = cor.test(within.babel[[5]][,3], within.babel[[8]][,3])
d= cor.test(within.babel[[6]][,3], within.babel[[7]][,3])
e=cor.test(within.babel[[6]][,3], within.babel[[8]][,3])
f=cor.test(within.babel[[7]][,3], within.babel[[8]][,3])

pearson_RPoccupancy_control_starved_pval = c(a$p.value, b$p.value, c$p.value, d$p.value, e$p.value, f$p.value)
pearson_RPoccupancy_control_starved = c(a$estimate, b$estimate, c$estimate, d$estimate, e$estimate, f$estimate)


g=cor.test(within.babel[[9]][,3], within.babel[[10]][,3])
h=cor.test(within.babel[[9]][,3], within.babel[[11]][,3])
i=cor.test(within.babel[[9]][,3], within.babel[[12]][,3])
j=cor.test(within.babel[[10]][,3], within.babel[[11]][,3])
k=cor.test(within.babel[[10]][,3], within.babel[[12]][,3])
l=cor.test(within.babel[[11]][,3], within.babel[[12]][,3])

pearson_RPoccupancy_atg_minus_fed_pval = c(g$p.value, h$p.value, i$p.value, j$p.value, k$p.value, l$p.value)
pearson_RPoccupancy_atg_minus_fed = c(g$estimate, h$estimate, i$estimate, j$estimate, k$estimate, l$estimate)


m=cor.test(within.babel[[13]][,3], within.babel[[14]][,3])
n=cor.test(within.babel[[13]][,3], within.babel[[15]][,3])
o=cor.test(within.babel[[13]][,3], within.babel[[16]][,3])
p=cor.test(within.babel[[14]][,3], within.babel[[15]][,3])
q=cor.test(within.babel[[14]][,3], within.babel[[16]][,3])
r=cor.test(within.babel[[15]][,3], within.babel[[16]][,3])

pearson_RPoccupancy_atg_minus_starved_pval = c(m$p.value, n$p.value, o$p.value, p$p.value, q$p.value, r$p.value)
pearson_RPoccupancy_atg_minus_starved = c(m$estimate, n$estimate, o$estimate, p$estimate, q$estimate, r$estimate)

rownamescor= c("rep1v2", "rep1v3", "rep1v4", "rep2v3", "rep2v4", "rep3v4")

pearson_RPoccupancy_pval = cbind(pearson_RPoccupancy_control_fed_pval, pearson_RPoccupancy_control_starved_pval, pearson_RPoccupancy_atg_minus_fed_pval, pearson_RPoccupancy_atg_minus_starved_pval)
rownames(pearson_RPoccupancy_pval) = rownamescor

pearson_RPoccupancy = cbind(pearson_RPoccupancy_control_fed, pearson_RPoccupancy_control_starved, pearson_RPoccupancy_atg_minus_fed, pearson_RPoccupancy_atg_minus_starved)
rownames(pearson_RPoccupancy) = rownamescor

pearsondata=cbind(cor_rna, cor_pval_rna, cor_rp, cor_pval_rp, pearson_RPoccupancy, pearson_RPoccupancy_pval)


##############################################################
##### Figure 2A-C #################
##############################################################

# plot violin plots and histograms: plot the change in counts in fed vs starved
#in order to use ggplots, need to reshape the data so that all the values are in one column and the variable (call it fed and starved) are in column 2

# plot a histogram of the means of the counts

#compare fed versus starved for atg12f/f cells
control_fed.rp= rowMeans(cov[,20:23], na.rm=TRUE)
control_starved.rp = rowMeans(cov[,24:27], na.rm=TRUE)
control_fed.rp = data.frame(xx =control_fed.rp, yy= rep("fed"))
control_starved.rp = data.frame(xx= control_starved.rp, yy= rep("starved"))
starvehistdat= rbind(control_fed.rp, control_starved.rp)

#plot the histogram of the means of the RPF counts from the bioreplicates fed versus starved for atg12f/f cells

ggplot(starvehistdat, aes(x=xx), show.legend=TRUE) + 
geom_histogram(data=subset(starvehistdat,yy == "fed"), fill = "red", alpha = 0.2, show.legend=TRUE) +
geom_histogram(data=subset(starvehistdat,yy == "starved"), fill = "blue", alpha = 0.2, show.legend=TRUE) + 
scale_y_log10() + scale_x_log10() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_line(colour = "grey")) +scale_fill_grey()

#compare the change in RPF counts in atg12f/f vs atg12KO cells

control_fed.rp= rowMeans(cov[,20:23], na.rm=TRUE)
atg_minus_fed.rp = rowMeans(cov[,28:31], na.rm=TRUE)
control_fed.rp = data.frame(xx =control_fed.rp, yy= rep("control"))
atg_minus_fed.rp = data.frame(xx= atg_minus_fed.rp, yy= rep("atg minus"))
starvehistdat= rbind(control_fed.rp, atg_minus_fed.rp)

#plot the histogram of the means of RPF counts from the bioreplicates atg12f/f vs atg12KO cells in fed conditions

ggplot(starvehistdat, aes(x=xx), show.legend=TRUE) + 
geom_histogram(data=subset(starvehistdat,yy == "control"), fill = "red", alpha = 0.2, show.legend=TRUE) +
geom_histogram(data=subset(starvehistdat,yy == "atg minus"), fill = "blue", alpha = 0.2, show.legend=TRUE) + 
scale_y_log10() + scale_x_log10() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_line(colour = "grey")) +scale_fill_grey()

#compare the change in RPF counts in atg12f/f vs atg12KO starved cells

control_starved.rp = rowMeans(cov[,24:27], na.rm=TRUE)
atg_minus_starved.rp = rowMeans(cov[,32:35], na.rm=TRUE)
control_starved.rp = data.frame(xx =control_starved.rp, yy= rep("control"))
atg_minus_starved.rp = data.frame(xx= atg_minus_starved.rp, yy= rep("atg minus"))
atg_starvehistdat= rbind(control_starved.rp, atg_minus_starved.rp)

#plot the histogram of the means of RPF counts from the bioreplicates atg12f/f vs atg12KO starved cells

ggplot(atg_starvehistdat, aes(x=xx), show.legend=TRUE) + 
geom_histogram(data=subset(atg_starvehistdat,yy == "control"), fill = "red", alpha = 0.2, show.legend=TRUE) +
geom_histogram(data=subset(atg_starvehistdat,yy == "atg minus"), fill = "blue", alpha = 0.2, show.legend=TRUE) + 
scale_y_log10() + scale_x_log10() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_line(colour = "grey")) +scale_fill_grey()


#plot individual bioreplicates by violin plot

cov.rp = cov[,-(3:19)]
cov.rp2 = melt(data=cov.rp)
cov.rp2$variable = gsub(".rp","", cov.rp2$variable)
cov.rp2$variable = gsub("atg_minus", "KO", cov.rp2$variable)
conditions = strsplit(cov.rp2$variable, "_")
atg12 = sapply(conditions, "[",1)
nutrient = sapply(conditions, "[",2)
biorep = sapply(conditions, "[",3)
cov.rp2 = cbind(cov.rp2, atg12, nutrient, biorep)

#compare the change in RPF counts in atg12f/f vs atg12KO fed cells
atg12_fed_comp = cov.rp2[grep("fed",cov.rp2$nutrient),]

ggplot2.violinplot(data=atg12_fed_comp, xName='biorep', yName='value', groupName='atg12') + scale_y_log10() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_line(colour = "grey"))

#compare the change in RPF counts in atg12f/f vs atg12KO starved cells
atg12_starved_comp = cov.rp2[grep("starved",cov.rp2$nutrient),]

ggplot2.violinplot(data=atg12_starved_comp, xName='biorep', yName='value', groupName='atg12') + scale_y_log10() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_line(colour = "grey"))

#compare the change in RPF counts in atg12f/f fed cells vs atg12f/f starved cells
starved_comp = cov.rp2[grep("control",cov.rp2$atg12),]

ggplot2.violinplot(data=starved_comp, xName='biorep', yName='value', groupName='nutrient') + scale_y_log10() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_line(colour = "grey"))

##############################################################
######## Figure 2D-F ###########
##############################################################

#plot fold change mRNA v RPF counts

load("all_groups.within.babel", .GlobalEnv)
#for some reason the P-value name is causing trouble so replace it
names=c("Gene", "Direction", "oneside", "twoside", "FDR")

	for (i in 1:16)
	{
		colnames(within.babel[[i]]) = names
		
	}
	
cov=cov[cov$gene_id %in% within.babel[[1]][,1],]

#atg dependent translation changes (Atg12KO/Atg12ff in fed conditions)

atg_fed_foldchange_rna_1 = cov$atg_minus_fed_1.rna /cov$control_fed_1.rna 
atg_fed_foldchange_rna_2 = cov$atg_minus_fed_2.rna /cov$control_fed_2.rna 
atg_fed_foldchange_rna_3 = cov$atg_minus_fed_3.rna /cov$control_fed_3.rna 
atg_fed_foldchange_rna_4 = cov$atg_minus_fed_4.rna /cov$control_fed_4.rna 

atg_fed_foldchange_rna_mean = (atg_fed_foldchange_rna_1 + atg_fed_foldchange_rna_2 + atg_fed_foldchange_rna_3 + atg_fed_foldchange_rna_4) /4

atg_fed_foldchange_rp_1 = cov$atg_minus_fed_1.rp / cov$control_fed_1.rp 
atg_fed_foldchange_rp_2 = cov$atg_minus_fed_2.rp / cov$control_fed_2.rp 
atg_fed_foldchange_rp_3 = cov$atg_minus_fed_3.rp / cov$control_fed_3.rp 
atg_fed_foldchange_rp_4 = cov$atg_minus_fed_4.rp / cov$control_fed_4.rp 

atg_fed_foldchange_rp_mean = (atg_fed_foldchange_rp_1 + atg_fed_foldchange_rp_2 + atg_fed_foldchange_rp_3 + atg_fed_foldchange_rp_4) /4

gene_id = as.character(cov$gene_id)
symbol = as.character(cov$symbol)
atg_fed_foldchange = cbind(gene_id, symbol, atg_fed_foldchange_rp_mean, atg_fed_foldchange_rna_mean)

atg_fed_foldchange = as.data.frame(atg_fed_foldchange)

atg_fed_foldchange$atg_fed_foldchange_rp_mean = as.numeric(as.character(atg_fed_foldchange$atg_fed_foldchange_rp_mean))

atg_fed_foldchange$atg_fed_foldchange_rna_mean = as.numeric(as.character(atg_fed_foldchange$atg_fed_foldchange_rna_mean))

atg_fed_foldchange$atg_fed_foldchange_rna_mean_log10 = log10(atg_fed_foldchange$atg_fed_foldchange_rna_mean)
atg_fed_foldchange$atg_fed_foldchange_rp_mean_log10 = log10(atg_fed_foldchange$atg_fed_foldchange_rp_mean)

validated = c("Brca2", "Atf4")
ggplot(atg_fed_foldchange, aes(x = atg_fed_foldchange_rna_mean_log10, y = atg_fed_foldchange_rp_mean_log10)) +geom_point(alpha=0.5, colour="#998ec3") + geom_point(data = atg_fed_foldchange[which(atg_fed_foldchange$symbol%in%validated),], aes(x = atg_fed_foldchange_rna_mean_log10, y = atg_fed_foldchange_rp_mean_log10), colour="#f1a340", size=5)+ geom_text(aes(label=ifelse(atg_fed_foldchange$symbol%in%validated, as.character(symbol),'')), size=7, hjust=.1, color="black") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_line(colour = "grey")) +scale_fill_grey()


#atg dependent translation in starvation changes (Atg12KO/Atg12ff in starve conditions)

atg_starved_foldchange_rna_1 = cov$atg_minus_starved_1.rna /cov$control_starved_1.rna 
atg_starved_foldchange_rna_2 = cov$atg_minus_starved_2.rna /cov$control_starved_2.rna 
atg_starved_foldchange_rna_3 = cov$atg_minus_starved_3.rna /cov$control_starved_3.rna 
atg_starved_foldchange_rna_4 = cov$atg_minus_starved_4.rna /cov$control_starved_4.rna 

atg_starved_foldchange_rna_mean = (atg_starved_foldchange_rna_1 + atg_starved_foldchange_rna_2 + atg_starved_foldchange_rna_3 + atg_starved_foldchange_rna_4) /4

atg_starved_foldchange_rp_1 = cov$atg_minus_starved_1.rp / cov$control_starved_1.rp 
atg_starved_foldchange_rp_2 = cov$atg_minus_starved_2.rp / cov$control_starved_2.rp 
atg_starved_foldchange_rp_3 = cov$atg_minus_starved_3.rp / cov$control_starved_3.rp 
atg_starved_foldchange_rp_4 = cov$atg_minus_starved_4.rp / cov$control_starved_4.rp 

atg_starved_foldchange_rp_mean = (atg_starved_foldchange_rp_1 + atg_starved_foldchange_rp_2 + atg_starved_foldchange_rp_3 + atg_starved_foldchange_rp_4) /4


gene_id = as.character(cov$gene_id)
symbol = as.character(cov$symbol)
atg_starved_foldchange = cbind(gene_id, symbol, atg_starved_foldchange_rp_mean, atg_starved_foldchange_rna_mean)

atg_starved_foldchange = as.data.frame(atg_starved_foldchange)

atg_starved_foldchange$atg_starved_foldchange_rp_mean = as.numeric(as.character(atg_starved_foldchange$atg_starved_foldchange_rp_mean))

atg_starved_foldchange$atg_starved_foldchange_rna_mean = as.numeric(as.character(atg_starved_foldchange$atg_starved_foldchange_rna_mean))

atg_starved_foldchange$atg_starved_foldchange_rna_mean_log10 = log10(atg_starved_foldchange$atg_starved_foldchange_rna_mean)
atg_starved_foldchange$atg_starved_foldchange_rp_mean_log10 = log10(atg_starved_foldchange$atg_starved_foldchange_rp_mean)

validated = c("Cdc25a", "Cpt2", "Pfkfb3")

ggplot(atg_starved_foldchange, aes(x = atg_starved_foldchange_rna_mean_log10, y = atg_starved_foldchange_rp_mean_log10)) +geom_point(alpha=0.5, colour="#998ec3") + geom_point(data = atg_starved_foldchange[which(atg_starved_foldchange$symbol%in%validated),], aes(x = atg_starved_foldchange_rna_mean_log10, y = atg_starved_foldchange_rp_mean_log10), colour="#f1a340", size=5) + geom_text(aes(label=ifelse(atg_starved_foldchange$symbol%in%validated, as.character(symbol),'')), size=7, position=position_jitter(width=0.1, height=0.2),color="black") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_line(colour = "grey")) +scale_fill_grey()


#starvation dependent translation changes (Atg12ff starved/ Atg12ff fed conditions)

starved_foldchange_rna_1 = cov$control_starved_1.rna /cov$control_fed_1.rna 
starved_foldchange_rna_2 = cov$control_starved_2.rna /cov$control_fed_2.rna 
starved_foldchange_rna_3 = cov$control_starved_3.rna /cov$control_fed_3.rna 
starved_foldchange_rna_4 = cov$control_starved_4.rna /cov$control_fed_4.rna 

starved_foldchange_rna_mean = (starved_foldchange_rna_1 + starved_foldchange_rna_2 + starved_foldchange_rna_3 + starved_foldchange_rna_4) /4

starved_foldchange_rp_1 = cov$control_starved_1.rp / cov$control_fed_1.rp 
starved_foldchange_rp_2 = cov$control_starved_2.rp / cov$control_fed_2.rp 
starved_foldchange_rp_3 = cov$control_starved_3.rp / cov$control_fed_3.rp 
starved_foldchange_rp_4 = cov$control_starved_4.rp / cov$control_fed_4.rp 

starved_foldchange_rp_mean = (starved_foldchange_rp_1 + starved_foldchange_rp_2 + starved_foldchange_rp_3 + starved_foldchange_rp_4) /4


gene_id = as.character(cov$gene_id)
symbol = as.character(cov$symbol)
starved_foldchange = cbind(gene_id, symbol, starved_foldchange_rp_mean, starved_foldchange_rna_mean)

starved_foldchange = as.data.frame(starved_foldchange)

starved_foldchange$starved_foldchange_rp_mean = as.numeric(as.character(starved_foldchange$starved_foldchange_rp_mean))

starved_foldchange$starved_foldchange_rna_mean = as.numeric(as.character(starved_foldchange$starved_foldchange_rna_mean))

starved_foldchange$starved_foldchange_rna_mean_log10 = log10(starved_foldchange$starved_foldchange_rna_mean)
starved_foldchange$starved_foldchange_rp_mean_log10 = log10(starved_foldchange$starved_foldchange_rp_mean)

validated = c("Eef2")

ggplot(starved_foldchange, aes(x = starved_foldchange_rna_mean_log10, y = starved_foldchange_rp_mean_log10)) +geom_point(alpha=0.5, colour="#998ec3") + geom_point(data = starved_foldchange[which(starved_foldchange$symbol%in%validated),], aes(x = starved_foldchange_rna_mean_log10, y = starved_foldchange_rp_mean_log10), colour="#f1a340", size=5) + geom_text(aes(label=ifelse(starved_foldchange$symbol%in%validated, as.character(symbol),'')), size=7, hjust=0.1, color="black") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_line(colour = "grey")) +scale_fill_grey()

##############################################################
##### Figure S2B #######
##############################################################

#plot mRNA v RPF counts

load("all_groups.within.babel", .GlobalEnv)
#for some reason the P-value name is causing trouble so replace it
names=c("Gene", "Direction", "oneside", "twoside", "FDR")

	for (i in 1:16)
	{
		colnames(within.babel[[i]]) = names
		
	}
	
cov=cov[cov$gene_id %in% within.babel[[1]][,1],]

#Atg12ff control media
ggplot(cov, aes(x = control_fed_1.rna, y = control_fed_1.rp, colour= ifelse(within.babel$control_fed_1$oneside < 0.025 | within.babel$control_fed_1$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = control_fed_2.rna, y = control_fed_2.rp, colour= ifelse(within.babel$control_fed_2$oneside < 0.025 | within.babel$control_fed_2$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = control_fed_3.rna, y = control_fed_3.rp, colour= ifelse(within.babel$control_fed_3$oneside < 0.025 | within.babel$control_fed_3$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = control_fed_4.rna, y = control_fed_4.rp, colour= ifelse(within.babel$control_fed_4$oneside < 0.025 | within.babel$control_fed_4$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")

#Atg12ff HBSS
ggplot(cov, aes(x = control_starved_1.rna, y = control_starved_1.rp, colour= ifelse(within.babel$control_starved_1$oneside < 0.025 | within.babel$control_starved_1$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = control_starved_2.rna, y = control_starved_2.rp, colour= ifelse(within.babel$control_starved_2$oneside < 0.025 | within.babel$control_starved_2$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = control_starved_3.rna, y = control_starved_3.rp, colour= ifelse(within.babel$control_starved_3$oneside < 0.025 | within.babel$control_starved_3$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = control_starved_4.rna, y = control_starved_4.rp, colour= ifelse(within.babel$control_starved_4$oneside < 0.025 | within.babel$control_starved_4$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")


#Atg12KO control
ggplot(cov, aes(x = atg_minus_fed_1.rna, y = atg_minus_fed_1.rp, colour= ifelse(within.babel$atg_minus_fed_1$oneside < 0.025 | within.babel$atg_minus_fed_1$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = atg_minus_fed_2.rna, y = atg_minus_fed_2.rp, colour= ifelse(within.babel$atg_minus_fed_2$oneside < 0.025 | within.babel$atg_minus_fed_2$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = atg_minus_fed_3.rna, y = atg_minus_fed_3.rp, colour= ifelse(within.babel$atg_minus_fed_3$oneside < 0.025 | within.babel$atg_minus_fed_3$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = atg_minus_fed_4.rna, y = atg_minus_fed_4.rp, colour= ifelse(within.babel$atg_minus_fed_4$oneside < 0.025 | within.babel$atg_minus_fed_4$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")

#Atg12KO HBSS
ggplot(cov, aes(x = atg_minus_starved_1.rna, y = atg_minus_starved_1.rp, colour= ifelse(within.babel$atg_minus_starved_1$oneside < 0.025 | within.babel$atg_minus_starved_1$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = atg_minus_starved_2.rna, y = atg_minus_starved_2.rp, colour= ifelse(within.babel$atg_minus_starved_2$oneside < 0.025 | within.babel$atg_minus_starved_2$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = atg_minus_starved_3.rna, y = atg_minus_starved_3.rp, colour= ifelse(within.babel$atg_minus_starved_3$oneside < 0.025 | within.babel$atg_minus_starved_3$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")
ggplot(cov, aes(x = atg_minus_starved_4.rna, y = atg_minus_starved_4.rp, colour= ifelse(within.babel$atg_minus_starved_4$oneside < 0.025 | within.babel$atg_minus_starved_4$oneside > 0.975, "One Side P-val", ""))) +geom_point() + scale_x_log10() + scale_y_log10() + theme(legend.position= "none")

