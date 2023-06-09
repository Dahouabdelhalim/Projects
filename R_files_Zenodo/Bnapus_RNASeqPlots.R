#install.packages("UpsetR")
#install.packages("sjPlot")
#install.packages("here)

library(UpSetR)
library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggExtra)
library(gtools)
library(gridExtra)
library(grid)
library(plyr)
library(zoo)
library(sjPlot)
library(here)
library(cowplot)
devtools::source_gist("524eade46135f6348140", filename = "ggplot_smooth_func.R")

# FUNCTIONS-----------------------------------------------------------------------------------------
GetSyntelogPairs<- function(R500File,TO1000File){
  TO1000<-read.delim(file=TO1000File, header=TRUE)
  R500<-read.delim(file=R500File, header=TRUE)
  
  colnames(TO1000)<-c("TO1000ID","TO1000Chr","TO1000Strand","TO1000Start","TO1000Stop","TO1000_tName","TO1000Exons","TO1000Length","TO1000GeneID","TO1000GeneName","TO1000Cov","TO1000FPKM","TO1000TPM")
  
  colnames(R500)<-c("R500ID","R500Chr","R500Strand","R500Start","R500Stop","R500_tName","R500Exons","R500Length","R500GeneID","R500GeneName","R500Cov","R500FPKM","R500TPM")
  
  R500 %>% mutate_if(is.factor,as.character)-> R500
  TO1000 %>% mutate_if(is.factor,as.character)-> TO1000
  TO1000.R500 %>% mutate_if(is.factor,as.character)-> TO1000.R500
  
  
  Ortho<-TO1000 %>% inner_join(TO1000.R500,by="TO1000_tName")%>% inner_join(R500,by="R500_tName")
  return(Ortho)
  
}


# READ FILES--------------------------------------------------------------------------------------

TO1000.R500 <- read.delim(here("TO1000.R500.anchors"), header=FALSE, comment.char="#")
colnames(TO1000.R500)<-c("TO1000_tName","R500_tName","BlockSize")

R500SyntelogFiles<-list.files(path=here("RNASeq"),pattern="R500",full.names = TRUE)
TO1000SyntelogFiles<-list.files(path=here("RNASeq"),pattern="TO1000",full.names = TRUE)

BNapusRNASeqArray<-vector("list", length = length(R500SyntelogFiles))


for (Indiv in 1:length(R500SyntelogFiles)){
  BNapusRNASeqArray[[Indiv]] <-GetSyntelogPairs(R500SyntelogFiles[Indiv],TO1000SyntelogFiles[Indiv])
}

names<- c("RS-100S1","RS-100S10","RS-100S5","RS-1100S1","RS-1100S10","RS-1100S5","RS-200S1","RS-200S10","RS-200S5","RS-300S1","RS-300S10","RS-400S1","RS-400S5","RS-600S1","RS-600S10","RS-600S5","RS-Parent")

#names<-c("RS-BNapus","RS-Parent")
names(BNapusRNASeqArray)<-names
for(Individ in 1:length(BNapusRNASeqArray)){ 
  BNapusRNASeqArray[[Individ]]<-mutate(BNapusRNASeqArray[[Individ]],SumTPM=TO1000TPM+R500TPM)%>%filter(SumTPM>10)%>%mutate(L2FC=foldchange2logratio(foldchange(num=TO1000TPM+1,denom=R500TPM+1),base=2)) %>% mutate(Bias= ifelse(L2FC > 3.5,"BnC Biased", ifelse(L2FC < -3.5,"Bna Biased","Nonbiased")))
}


#Get 2:2 Homoeolog list-------------------------------------------------------------------------------------
TwotoTwoHE<-list.files(path=here("WGS"),pattern="TwotoTwo",full.names = T)
TwotoTwoGenes<-c()
for(x in 1:length(TwotoTwoHE)){TwotoTwoGenes[[x]]<-read.delim(file=TwotoTwoHE[x], header=T,sep = ",",stringsAsFactors = F)}
names(TwotoTwoGenes)<-gsub(TwotoTwoHE,pattern = ".csv",replacement = "")
names(TwotoTwoGenes)<-gsub(names(TwotoTwoGenes),pattern = here("WGS/"),replacement = "")

# LIST TO OBJECTS--------------------------------------------------------------------------------------


list2env(TwotoTwoGenes,envir =.GlobalEnv)
list2env(BNapusRNASeqArray,envir =.GlobalEnv)

#1.2 Upsetr of Homoeologous bias_________________________________________________________________

#filter out lowly expressed genes, and get only 2:2 homoeologs identified by WGS
TwotoTwoHomoeo<-c()
TwotoTwoHomoeoTPM<-c()
names<- c("RS-100S1","RS-100S10","RS-100S5","RS-1100S1","RS-1100S10","RS-1100S5","RS-200S1","RS-200S10","RS-200S5","RS-300S1","RS-300S10","RS-400S1","RS-400S5","RS-600S1","RS-600S10","RS-600S5","RS-Parent")

IndivNames<-mixedsort(gsub(names,pattern = "RS-",replacement = ""))

for(Individ in IndivNames){ 
  if(grepl("Parent",Individ)){
    TwotoTwoHomoeo[[Individ]]<-get(sprintf("RS-%s",Individ))
    TwotoTwoHomoeo[[Individ]]<-mutate(TwotoTwoHomoeo[[Individ]],SumTPM=TO1000TPM+R500TPM)%>%filter(SumTPM>10)%>%mutate(L2FC=foldchange2logratio(foldchange(num=TO1000TPM+1,denom=R500TPM+1),base=2)) %>% mutate(Bias= ifelse(L2FC > 3.5,"BnC Biased", ifelse(L2FC < -3.5,"Bna Biased","Nonbiased")))
    TwotoTwoHomoeo[[Individ]]$BnC_Density<-rollapply(TwotoTwoHomoeo[[Individ]]$Bias=="BnC Biased", width = 10, by = 1, FUN = sum, na.rm = TRUE,partial=T)
    TwotoTwoHomoeo[[Individ]]$BnA_Density<-rollapply(TwotoTwoHomoeo[[Individ]]$Bias=="Bna Biased", width = 10, by = 1, FUN = sum, na.rm = TRUE,partial=T)
    
        }
  else{
    TwotoTwoHomoeo[[Individ]]<-subset(get(sprintf("RS-%s",Individ)),TO1000_tName %in% get(sprintf("Ds-%sTwotoTwo",Individ))$TO1000)
    TwotoTwoHomoeo[[Individ]]<-mutate(TwotoTwoHomoeo[[Individ]],SumTPM=TO1000TPM+R500TPM)%>%filter(SumTPM>10)%>%mutate(L2FC=foldchange2logratio(foldchange(num=TO1000TPM+1,denom=R500TPM+1),base=2)) %>% mutate(Bias= ifelse(L2FC > 3.5,"BnC Biased", ifelse(L2FC < -3.5,"Bna Biased","Nonbiased")))
    TwotoTwoHomoeo[[Individ]]$BnC_Density<-rollapply(TwotoTwoHomoeo[[Individ]]$Bias=="BnC Biased", width = 10, by = 1, FUN = sum, na.rm = TRUE,partial=T)
    TwotoTwoHomoeo[[Individ]]$BnA_Density<-rollapply(TwotoTwoHomoeo[[Individ]]$Bias=="Bna Biased", width = 10, by = 1, FUN = sum, na.rm = TRUE,partial=T)
    
        }
}


## @knitr Upset Common 2:2 Gene
CommonGenesS1<-Reduce(intersect, list(TwotoTwoHomoeo$`100S1`$TO1000_tName, TwotoTwoHomoeo$`200S1`$TO1000_tName,
                                      TwotoTwoHomoeo$`300S1`$TO1000_tName, TwotoTwoHomoeo$`400S1`$TO1000_tName,
                                      TwotoTwoHomoeo$`600S1`$TO1000_tName, TwotoTwoHomoeo$`1100S1`$TO1000_tName,TwotoTwoHomoeo$Parent$TO1000_tName))

CommonGenesS5<-Reduce(intersect, list(TwotoTwoHomoeo$`100S5`$TO1000_tName, TwotoTwoHomoeo$`200S5`$TO1000_tName,
                                      TwotoTwoHomoeo$`400S5`$TO1000_tName,TwotoTwoHomoeo$`600S5`$TO1000_tName,
                                      TwotoTwoHomoeo$`1100S5`$TO1000_tName,TwotoTwoHomoeo$Parent$TO1000_tName))
                      

CommonGenesS10<-Reduce(intersect, list(TwotoTwoHomoeo$`100S10`$TO1000_tName, TwotoTwoHomoeo$`200S10`$TO1000_tName,
                                      TwotoTwoHomoeo$`300S10`$TO1000_tName, TwotoTwoHomoeo$`600S10`$TO1000_tName, 
                                      TwotoTwoHomoeo$`1100S10`$TO1000_tName,TwotoTwoHomoeo$Parent$TO1000_tName))

#BnC biased S1
Bol_Dominant100S1<-filter(TwotoTwoHomoeo$`100S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC > 3.5)
Bol_Dominant200S1<-filter(TwotoTwoHomoeo$`200S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC > 3.5)
Bol_Dominant300S1<-filter(TwotoTwoHomoeo$`300S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC > 3.5)
Bol_Dominant400S1<-filter(TwotoTwoHomoeo$`400S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC > 3.5)
Bol_Dominant600S1<-filter(TwotoTwoHomoeo$`600S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC > 3.5)
Bol_Dominant1100S1<-filter(TwotoTwoHomoeo$`1100S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC > 3.5)
Bol_DominantParentS1<-filter(TwotoTwoHomoeo$`Parent`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC > 3.5)

#BnA biased S1
Bra_Dominant100S1<-filter(TwotoTwoHomoeo$`100S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC < -3.5)
Bra_Dominant200S1<-filter(TwotoTwoHomoeo$`200S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC < -3.5)
Bra_Dominant300S1<-filter(TwotoTwoHomoeo$`300S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC < -3.5)
Bra_Dominant400S1<-filter(TwotoTwoHomoeo$`400S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC < -3.5)
Bra_Dominant600S1<-filter(TwotoTwoHomoeo$`600S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC < -3.5)
Bra_Dominant1100S1<-filter(TwotoTwoHomoeo$`1100S1`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC < -3.5)
Bra_DominantParentS1<-filter(TwotoTwoHomoeo$`Parent`, TO1000_tName %in% CommonGenesS1)%>%filter(L2FC < -3.5)

#BnC biased S5
Bol_Dominant100S5<-filter(TwotoTwoHomoeo$`100S5`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC > 3.5)
Bol_Dominant200S5<-filter(TwotoTwoHomoeo$`200S5`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC > 3.5)
Bol_Dominant400S5<-filter(TwotoTwoHomoeo$`400S5`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC > 3.5)
Bol_Dominant600S5<-filter(TwotoTwoHomoeo$`600S5`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC > 3.5)
Bol_Dominant1100S5<-filter(TwotoTwoHomoeo$`1100S5`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC > 3.5)
Bol_DominantParentS5<-filter(TwotoTwoHomoeo$`Parent`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC > 3.5)

#BnA biased S5
Bra_Dominant100S5<-filter(TwotoTwoHomoeo$`100S5`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC < -3.5)
Bra_Dominant200S5<-filter(TwotoTwoHomoeo$`200S5`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC < -3.5)
Bra_Dominant400S5<-filter(TwotoTwoHomoeo$`400S5`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC < -3.5)
Bra_Dominant600S5<-filter(TwotoTwoHomoeo$`600S5`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC < -3.5)
Bra_Dominant1100S5<-filter(TwotoTwoHomoeo$`1100S5`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC < -3.5)
Bra_DominantParentS5<-filter(TwotoTwoHomoeo$`Parent`, TO1000_tName %in% CommonGenesS5)%>%filter(L2FC < -3.5)

Bol_Dominant100S10<-filter(TwotoTwoHomoeo$`100S10`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC > 3.5)
Bol_Dominant200S10<-filter(TwotoTwoHomoeo$`200S10`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC > 3.5)
Bol_Dominant300S10<-filter(TwotoTwoHomoeo$`300S10`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC > 3.5)
Bol_Dominant600S10<-filter(TwotoTwoHomoeo$`600S10`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC > 3.5)
Bol_Dominant1100S10<-filter(TwotoTwoHomoeo$`1100S10`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC > 3.5)
Bol_DominantParentS10<-filter(TwotoTwoHomoeo$`Parent`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC > 3.5)

Bra_Dominant100S10<-filter(TwotoTwoHomoeo$`100S10`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC < -3.5)
Bra_Dominant200S10<-filter(TwotoTwoHomoeo$`200S10`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC < -3.5)
Bra_Dominant300S10<-filter(TwotoTwoHomoeo$`300S10`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC < -3.5)
Bra_Dominant600S10<-filter(TwotoTwoHomoeo$`600S10`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC < -3.5)
Bra_Dominant1100S10<-filter(TwotoTwoHomoeo$`1100S10`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC < -3.5)
Bra_DominantParentS10<-filter(TwotoTwoHomoeo$`Parent`, TO1000_tName %in% CommonGenesS10)%>%filter(L2FC < -3.5)

S1BolBiased<-list(Rs_100=Bol_Dominant100S1$TO1000_tName, Rs_200=Bol_Dominant200S1$TO1000_tName,
                  Rs_300=Bol_Dominant300S1$TO1000_tName, Rs_400=Bol_Dominant400S1$TO1000_tName, 
                  Rs_600=Bol_Dominant600S1$TO1000_tName, Rs_1100=Bol_Dominant1100S1$TO1000_tName,
                  Rs_Parent=Bol_DominantParentS1$TO1000_tName)

S1BraBiased <-list(Rs_100=Bra_Dominant100S1$TO1000_tName, Rs_200=Bra_Dominant200S1$TO1000_tName,
                   Rs_300=Bra_Dominant300S1$TO1000_tName, Rs_400=Bra_Dominant400S1$TO1000_tName, 
                   Rs_600=Bra_Dominant600S1$TO1000_tName, Rs_1100=Bra_Dominant1100S1$TO1000_tName,
                   Rs_Parent=Bra_DominantParentS1$TO1000_tName)

S5BolBiased<-list(Rs_100=Bol_Dominant100S5$TO1000_tName, Rs_200=Bol_Dominant200S5$TO1000_tName,
                  Rs_400=Bol_Dominant400S5$TO1000_tName,Rs_600=Bol_Dominant600S5$TO1000_tName, 
                  Rs_1100=Bol_Dominant1100S5$TO1000_tName,Rs_Parent=Bol_DominantParentS5$TO1000_tName)

S5BraBiased<-list(Rs_100=Bra_Dominant100S5$TO1000_tName, Rs_200=Bra_Dominant200S5$TO1000_tName,
                  Rs_400=Bra_Dominant400S5$TO1000_tName,Rs_600=Bra_Dominant600S5$TO1000_tName, 
                  Rs_1100=Bra_Dominant1100S5$TO1000_tName,Rs_Parent=Bra_DominantParentS5$TO1000_tName)
               

S10BolBiased<-list(Rs_100=Bol_Dominant100S10$TO1000_tName, Rs_200=Bol_Dominant200S10$TO1000_tName,
                  Rs_300=Bol_Dominant300S10$TO1000_tName, Rs_600=Bol_Dominant600S10$TO1000_tName, 
                  Rs_1100=Bol_Dominant1100S10$TO1000_tName,Rs_Parent=Bol_DominantParentS10$TO1000_tName)

S10BraBiased<-list(Rs_100=Bra_Dominant100S10$TO1000_tName, Rs_200=Bra_Dominant200S10$TO1000_tName,
                   Rs_300=Bra_Dominant300S10$TO1000_tName, Rs_600=Bra_Dominant600S10$TO1000_tName, 
                   Rs_1100=Bra_Dominant1100S10$TO1000_tName,Rs_Parent=Bra_DominantParentS10$TO1000_tName)

upset(fromList(S1BolBiased), sets = c("Rs_100", "Rs_200", "Rs_300", "Rs_400", "Rs_600", 
                                     "Rs_1100","Rs_Parent"), order.by = "freq",nsets = 7,keep.order = TRUE, 
      mainbar.y.label = "BnC Biased Homoeolog pairs S1",
      main.bar.color = "#E41A1C",text.scale = c(2.2, 2, 1.8, 1.8, 1.8,0),
      point.size = 3.5, line.size = 2, show.numbers = F)

upset(fromList(S1BraBiased), sets = c("Rs_100", "Rs_200", "Rs_300", "Rs_400", "Rs_600", 
                                      "Rs_1100","Rs_Parent"), order.by = "freq",nsets = 7,keep.order = TRUE, 
      mainbar.y.label = "BnA Biased Homoeolog pairs S1",
      main.bar.color = "#E41A1C",text.scale = c(2.2, 2, 1.8, 1.8, 1.8,0),
      point.size = 3.5, line.size = 2, show.numbers = F)

upset(fromList(S5BolBiased), sets = c("Rs_100", "Rs_200", "Rs_400", "Rs_600", 
                                      "Rs_1100","Rs_Parent"), order.by = "freq",nsets = 6,keep.order = TRUE, 
      mainbar.y.label = "BnC Biased Homoeolog pairs S5",
      main.bar.color = "#4DAF4A",text.scale = c(2.2, 2, 1.8, 1.8, 1.8,0),
      point.size = 3.5, line.size = 2, show.numbers = F)

upset(fromList(S5BraBiased), sets = c("Rs_100", "Rs_200", "Rs_400", "Rs_600", 
                                      "Rs_1100","Rs_Parent"), order.by = "freq",nsets = 6,keep.order = TRUE, 
      mainbar.y.label = "BnA Biased Homoeolog pairs S5",
      main.bar.color = "#4DAF4A",text.scale = c(2.2, 2, 1.8, 1.8, 1.8,0),
      point.size = 3.5, line.size = 2, show.numbers = F)

upset(fromList(S10BolBiased), sets = c("Rs_100", "Rs_200", "Rs_300", "Rs_600", 
                                      "Rs_1100","Rs_Parent"), order.by = "freq",nsets = 6,keep.order = TRUE, 
      mainbar.y.label = "BnC Biased Homoeolog pairs S10",
      main.bar.color = "#377EB8",text.scale = c(2.2, 2, 1.8, 1.8, 1.8,0),
      point.size = 3.5, line.size = 2, show.numbers = F)

upset(fromList(S10BraBiased), sets = c("Rs_100", "Rs_200", "Rs_300", "Rs_600", 
                                       "Rs_1100","Rs_Parent"), order.by = "freq",nsets = 6,keep.order = TRUE, 
      mainbar.y.label = "BnA Biased Homoeolog pairs S10",
      main.bar.color = "#377EB8",text.scale = c(2.2, 2, 1.8, 1.8, 1.8,0),
      point.size = 3.5, line.size = 2, show.numbers = F)

 #2 Histogram of biased HE pairs ________________________________________________________________
 
 BiasHist<-function(In){
dat<-In

ggplot(data=dat,aes(x=L2FC,fill=Bias))+xlim(-20,20)+ ylim(0,2200)+labs(x="Log2(Expr_BnC/Expr_BnA)")+
   geom_histogram(bins=30)+
  theme_cowplot(12)+
   scale_fill_manual(values=c("blue","red","grey"))+
  annotate(geom="text", x=0, y=1000, size=5, label=paste0("All pairs","\\n","N=", nrow(dat),"\\n","Bias =", round(median(dat$L2FC),3)))+
  annotate(geom="text", x=9, y=300, size=5, label=paste0("BnC \\n Biased pairs","\\n","N=",nrow(dat[dat$Bias=="BnC Biased",]), "\\n","Bias =", round(median(dat[dat$Bias=="BnC Biased",]$L2FC),3)))+
  annotate(geom="text", x=-9, y=300, size=5, label=paste0("BnA \\n Biased pairs","\\n","N=",nrow(dat[dat$Bias=="Bna Biased",]), "\\n","Bias =", round(median(dat[dat$Bias=="Bna Biased",]$L2FC),3)))+
  NULL
 }
 
 BiasHistPlots<-lapply(X=TwotoTwoHomoeo,FUN=BiasHist)
 
save_plot(plot_grid(BiasHistPlots[[17]],BiasHistPlots[[1]],BiasHistPlots[[2]],BiasHistPlots[[3]],nrow=4,labels = c("Parent","EL100S1","EL100S5","EL100S10")),base_aspect_ratio = .8, base_height = 12, filename ="EL100_HEBFig.pdf",device = "pdf",path = here())
save_plot(plot_grid(BiasHistPlots[[17]],BiasHistPlots[[4]],BiasHistPlots[[5]],BiasHistPlots[[6]],nrow=4,labels = c("Parent","EL200S1","EL200S5","EL200S10")),base_aspect_ratio = .8, base_height = 12, filename ="EL200_HEBFig.pdf",device = "pdf",path = here())
save_plot(plot_grid(BiasHistPlots[[17]],BiasHistPlots[[7]],BiasHistPlots[[8]],nrow=3,labels = c("Parent","EL300S1","EL300S10")),base_aspect_ratio = .8, base_height = 12, filename ="EL300_HEBFig.pdf",device = "pdf",path = here())
save_plot(plot_grid(BiasHistPlots[[17]],BiasHistPlots[[9]],BiasHistPlots[[10]],nrow=3,labels = c("Parent","EL400S1","EL400S5")),base_aspect_ratio = .8, base_height = 12, filename ="EL400_HEBFig.pdf",device = "pdf",path = here())
save_plot(plot_grid(BiasHistPlots[[17]],BiasHistPlots[[11]],BiasHistPlots[[12]],BiasHistPlots[[13]],nrow=4,labels = c("Parent","EL600S1","EL600S5","EL600S10")),base_aspect_ratio = .8,base_height = 12, filename ="EL600_HEBFig.pdf",device = "pdf",path =here())
save_plot(plot_grid(BiasHistPlots[[17]],BiasHistPlots[[14]],BiasHistPlots[[15]],BiasHistPlots[[16]],nrow=4,labels = c("Parent","EL1100S1","EL1100S5","EL1100S10")),base_aspect_ratio = .8, base_height = 12, filename ="EL1100_HEBFig.pdf",device = "pdf",path = here())

Chisq<-c()
lapply(TwotoTwoHomoeo,function(x){ testres<-wilcox.test(x$L2FC,mu=0,alternative = "g"); return(testres$p.value)})
lapply(TwotoTwoHomoeo,function(x){ testres<-wilcox.test(x$L2FC,TwotoTwoHomoeo$Parent$L2FC); return(testres$p.value)})
Chisq<-lapply(TwotoTwoHomoeo,function(x){ Chires<-chisq.test(c(sum(x$Bias=="BnC Biased"),sum(x$Bias=="Bna Biased")),p=c(0.5,0.5)); return(Chires)})
lapply(TwotoTwoHomoeo,function(x){ Chires<-chisq.test(c(sum(x$Bias=="BnC Biased"),sum(x$Bias=="Bna Biased")),p=c(3688/(3688+1796),1796/(3688+1796))); return(Chires)})

Sample<-c()
BnC.Obs<-c()
BnA.Obs<-c()
BnC.Exp<-c()
BnA.Exp<-c()
ChiSq<-c()
P.value<-c()
for(n in 1:length(Chisq)){
  Sample[n]<-names(Chisq[n])
  BnC.Obs[n]<-Chisq[[n]]$observed[1]
  BnA.Obs[n]<-Chisq[[n]]$observed[2]
  BnC.Exp[n]<-Chisq[[n]]$expected[1]
  BnA.Exp[n]<-Chisq[[n]]$expected[2]
  ChiSq[n]<-Chisq[[n]]$statistic
  P.value[n]<-Chisq[[n]]$p.value
}
ChiRes<-data.frame(Sample=Sample, `BnC Observed`=BnC.Obs,`BnC Expected`=BnC.Exp,`BnA Observed`=BnA.Obs,`BnA Expected`=BnA.Exp,`Chi Squared`=round(ChiSq,2),`P value`=P.value)

sjPlot::tab_df(x = as.data.frame(ChiRes), alternate.rows = T,title = "Homeolog Expression Bias Chi Squared table",file="HEBChiRes.doc")


