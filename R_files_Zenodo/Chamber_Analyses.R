#Chamber Number and Size Analyses

library(lme4)
library(nlme)
library(ape)
library(caper)
library(phytools)
library(phylolm)
library(MASS)
library(car)
library(scales)


#Create Nest_Chambers.csv" file from "ChamberWidths_perChamber.csv" file
WPC<-read.csv("....Chamber Widths_perChamber.csv")
NameKey<-read.csv("....filename to species name key.csv")


#MAKING Nest_Chambers.csv----
#Remove all Entrance Nodes
WPC_C<-WPC[!WPC$Structure.Type=="Entrance",]
WPC_C<-WPC_C[,c(1:4)]

WPC_MeanWidth<-aggregate(WPC_C$Maximum.Width,list(WPC_C$NestID),mean)
colnames(WPC_MeanWidth)<-c("Nest.ID","ch.width.avg")
WPC_ChNum<-aggregate(WPC_C$Structure.Type,list(WPC_C$NestID),length)
colnames(WPC_ChNum)<-c("Nest.ID","ch.num")

Nest.Chambers1<-merge(WPC_MeanWidth,WPC_ChNum,by="Nest.ID",all=TRUE)
filename.list<-paste(Nest.Chambers1$Nest.ID,".csv",sep="") #Make nest ID be csv so it will match with Name Key
Nest.Chambers2<-cbind(filename.list,Nest.Chambers1)
Nest.Chambers3<-merge(Nest.Chambers2,NameKey,by="filename.list",all.x = TRUE)

GenusSp<-data.frame(do.call("rbind", strsplit(as.character(Nest.Chambers3$Name), " ", fixed = TRUE)))
colnames(GenusSp)<-c("Genus","sp")
Nest.Chambers4<-cbind(GenusSp,Nest.Chambers3)

write.csv(Nest.Chambers4,"...Nest Chambers.csv")
#Creates a spreadsheet with chamber number and avg chamber width per nest

#GENERATE DataC----
#DataC is used in the analyses below
#Prepare Sean's Ant Length Measurements:
sean<-read.csv("....Sean - Morphometric measures.csv")
sean
sean$Weber.s.length<-as.numeric(sean$Weber.s.length)
mean_webers<-aggregate(sean$Weber.s.length,list(sean$Species.Name),mean)
mean_webers
colnames(mean_webers)<-c("Name","length.webers")

#Calculate Pheidole values and replace in spreadsheet
#Pheidole has discrete polymorphism, so calculating weighted average: 90% minors, 10% majors

pmorrisii<-sean$Weber.s.length[sean$Species.Name=="Pheidole morrisii"]
minor<-pmorrisii[1:3]
major<-pmorrisii[4:8]
weighted.avg.pmorrisii<-0.9*(mean(minor))+0.1*(mean(major))

pdentata<-sean$Weber.s.length[sean$Species.Name=="Pheidole dentata"]
minor<-pdentata[1:4]
major<-pdentata[5:7]
weighted.avg.pdentata<-0.9*(mean(minor))+0.1*(mean(major))

#Add weighted values to mean_webers spreadsheet summarizing weber's length value for each species
mean_webers$length.webers[mean_webers$Name=="Pheidole morrisii"]<-weighted.avg.pmorrisii
mean_webers$length.webers[mean_webers$Name=="Pheidole dentata"]<-weighted.avg.pdentata

#Add a genus and species column
mean_webers$Genus<-sapply(strsplit(as.character(mean_webers$Name), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
mean_webers$sp<-sapply(strsplit(as.character(mean_webers$Name), "\\\\ "), "[[", 2) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
write.csv(mean_webers,"....Ant Lengths - sean.csv")

#Open other datasheets
WidthpCh<-read.csv("....\\Chamber Widths_perChamber.csv")
NameKey<-read.csv("....filename to species name key.csv")
Col.Size<-read.csv("....Colony size.csv")
Num<-read.csv("....Nest Chambers.csv")
AntLength<-read.csv("....Ant Lengths - sean.csv")
names(Num)[names(Num) == "Nest.ID"] <- "NestID"

Data1<-merge(WidthpCh, Num,by="NestID",all=TRUE)
Data2<-merge(Data1,AntLength,by="Name",all.x=TRUE)

#Remove E1s
Data3<-Data2[!(Data2$Structure.ID=="E1")| is.na(Data2$Structure.ID),]

#Divide Chamber width by ant length (webers) - CWdithS per chamber
Data3$CWidth_S<-Data3$Maximum.Width/Data3$length.webers #S for standardized

#Calculate CV, max size, mean, SE for each nest ID using aggregate function
#Maximum Width
  maxChW<-aggregate(Data3[,"CWidth_S"],list(Data3$NestID),max)
  colnames(maxChW)<-c("NestID","MaxW_S")
  
#Mean
  meanChW<-aggregate(Data3[,"CWidth_S"],list(Data3$NestID),mean)
  colnames(meanChW)<-c("NestID","MeanW_S")
  
#Standard Deviation
  SDChW<-aggregate(Data3[,"CWidth_S"],list(Data3$NestID),sd)
  colnames(SDChW)<-c("NestID","SDW_S")
  
#Coefficient of Variation
  CV <- function(x){
    sd(x)/mean(x)
  }
  
  CVChW<-aggregate(Data3[,"CWidth_S"],list(Data3$NestID),CV)
  colnames(CVChW)<-c("NestID","CVW_S")

#Standard Error of Chamber Width 
  StE <- function(x) {
    sd(x)/sqrt(length(x))
  }
  
  SEChW<-aggregate(Data3[,"CWidth_S"],list(Data3$NestID),StE)
  colnames(SEChW)<-c("NestID","SEW_S")  
  

#MERGE DATA
  #Ref<-merge(NullTraitData,ColSize,by="Name",all.x=TRUE)
  ChData1<-merge(maxChW, meanChW,by="NestID",all.x=TRUE)
  ChData2<-merge(SDChW, CVChW,by="NestID",all.x=TRUE)
  ChData3<-merge(SEChW, ChData1,by="NestID",all.x=TRUE)
  ChData<-merge(ChData3, ChData2,by="NestID",all.x=TRUE)
  

#Combine with Ant Size and Colony size
  
  #Assemble all data into one data-frame
  
  #Add Chamber number count
  #names(Num)[names(Num)=="NestID"] <- "NestID"
  Chambers<-merge(ChData,Num,by="NestID",all.x=TRUE)
  # ChData$Nest.ID<-paste(ChData$Nest.ID,".csv",sep="")
  # names(NameKey)[names(NameKey)=="filename.list"] <- "Nest.ID"
  # Chambers<-merge(ChData,NameKey,by="Nest.ID",all.x=TRUE)
  DataC<-merge(Chambers,Col.Size,by="Name",all.x=TRUE)
  #Data<-merge(data.size,AntLength,by="Name",all.x=TRUE)
  str(DataC)
  write.csv(DataC,"....DataC.csv")
  
  
#Max W v Col Size
  
  plot(MaxW_S~log(AvgColSize),data=DataC)
  summary(lm(MaxW_S~log(AvgColSize),data=DataC))
 
  
#CV v Col Size
  
  plot(CVW_S~log(AvgColSize),data=DataC)
  summary(lm(CVW_S~log(AvgColSize),data=DataC))
  
#MEAN WIDTH v COLONY SIZE - linear and nonlinear  - calculate the PIC values----
  #Chamber Size~Col Size----
  #___________________________
  #*Preparing Tree----
  #____________________________
  setwd("....")
  ant.tree.moreau<-read.nexus("MoreauTree2016.nex")
  plot(ant.tree.moreau)
  str(ant.tree.moreau)
  ant.tree.moreau$tip.label
  
  #Prune tree to the genera in my dataset - Should be 19 genera
  tips<-c("Acromyrmex_versicolor","Aphaenogaster_occidentalis_NW","Camponotus_maritimus","Diacamma_rugosum","Dinoponera_australis","Dorymyrmex_bicolor","Ectatomma_opaciventre","Formica_moki" ,"Monomorium_pharaonis","Mycetagroicus_triangularis","Mycetarotes_acutus","Mycetophylax_conformis","Mycocepurus_goeldii","Odontomachus_coquereli","Pachycondyla_harpax","Pheidole_longispinosa","Pogonomyrmex_angustus", "Prenolepis_imparis" ,"Sericomyrmex_Sp","Trachymyrmex_arizonensis","Veromessor_andrei")
  pruned.tree.size<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])
  plot(pruned.tree.size)
  pruned.tree.size$tip.label
  
  #Rename tip labels to create genus level tree
  for (i in 1:length(pruned.tree.size$tip.label)) {
    split.tips<-strsplit(pruned.tree.size$tip.label[i],"_")
    genus.tip<-split.tips[[1]]
    pruned.tree.size$tip.label[i]<-genus.tip[1]
  }
  plot(pruned.tree.size) 
  
  #*Prepare Trait Data----
  #__________________________________
  DataC<-read.csv("....DataC.csv")
  
  mean_trait_size2<-aggregate(DataC[],list(DataC$Name),mean,na.rm=TRUE)
  
  # Rename Genus species column name to "species"
  colnames(mean_trait_size2)[colnames(mean_trait_size2)=="Group.1"] <- "species"
  
  
  #*Create data subsets----
  #_________________________________________
  #Must be removed
  mean_trait_size_sub0<-subset(mean_trait_size2,species !="Aphaenogaster" & species != "Pheidole oxyops" & species != "Mycocepurus goeldii" & species !="Myrmecocystus navajo" & species !="Formica japonica" &  species != "Forelius" & species != "Myrmecia dispar" &  species !="Dinoponera australis" & species !="Acromyrmex subterraneus"  & species!="Odontomachus chelifer")
  
  #Optional removal for diffeerent data subsets: 13 subsets for 20 genera
  mean_trait_size_sub1<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")  #archboldi
  mean_trait_size_sub2<-subset(mean_trait_size_sub0,species!="Acromyrmex balzani" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
  mean_trait_size_sub3<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster floridana" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
  mean_trait_size_sub4<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster floridana"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
  mean_trait_size_sub5<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes parallelus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
  mean_trait_size_sub6<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
  mean_trait_size_sub7<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
  mean_trait_size_sub8<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex mayri" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
  mean_trait_size_sub9<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
  mean_trait_size_sub10<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
  mean_trait_size_sub11<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex mayri" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
  mean_trait_size_sub12<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens") 
  mean_trait_size_sub13<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica archboldi" & species !="Formica subaenescens")       #pallidefulva
  mean_trait_size_sub14<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica archboldi" & species !="Formica pallidefulva" & species !="Formica subaenescens") #dolosa
  mean_trait_size_sub15<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica archboldi")      #subaenescens
  
  
  mean_trait_size_sub1$subset.no=rep(1,nrow(mean_trait_size_sub1))
  mean_trait_size_sub2$subset.no=rep(2,nrow(mean_trait_size_sub2))
  mean_trait_size_sub3$subset.no=rep(3,nrow(mean_trait_size_sub3))
  mean_trait_size_sub4$subset.no=rep(4,nrow(mean_trait_size_sub4))
  mean_trait_size_sub5$subset.no=rep(5,nrow(mean_trait_size_sub5))
  mean_trait_size_sub6$subset.no=rep(6,nrow(mean_trait_size_sub6))
  mean_trait_size_sub7$subset.no=rep(7,nrow(mean_trait_size_sub7))
  mean_trait_size_sub8$subset.no=rep(8,nrow(mean_trait_size_sub8))
  mean_trait_size_sub9$subset.no=rep(9,nrow(mean_trait_size_sub9))
  mean_trait_size_sub10$subset.no=rep(10,nrow(mean_trait_size_sub10))
  mean_trait_size_sub11$subset.no=rep(11,nrow(mean_trait_size_sub11))
  mean_trait_size_sub12$subset.no=rep(12,nrow(mean_trait_size_sub12))
  mean_trait_size_sub13$subset.no=rep(13,nrow(mean_trait_size_sub13))
  mean_trait_size_sub14$subset.no=rep(14,nrow(mean_trait_size_sub14))
  mean_trait_size_sub15$subset.no=rep(15,nrow(mean_trait_size_sub15))

  
  Data.Subset<-list(mean_trait_size_sub1,
                    mean_trait_size_sub2,
                    mean_trait_size_sub3,
                    mean_trait_size_sub4,
                    mean_trait_size_sub5,
                    mean_trait_size_sub6, 
                    mean_trait_size_sub7,
                    mean_trait_size_sub8,
                    mean_trait_size_sub9,
                    mean_trait_size_sub10, 
                    mean_trait_size_sub11,
                    mean_trait_size_sub12,
                    mean_trait_size_sub13,
                    mean_trait_size_sub14,
                    mean_trait_size_sub15)
  
      
 
      #Initialize output lists
      #_______________________
      chsize.pic.coef<-c()
      chsize.pic.lci<-c()
      chsize.pic.uci<-c()
      chsize.pic.r2<-c()
      chsize.pic.rlci<-c()
      chsize.pic.ruci<-c()
      chsize.pic.p<-c()
      chsize.pic.ss<-c()
      chsize.pic.aic<-c()
      
      tree.type<-c()
      
      #Run analysis for each data subset
      #___________________________________
      
      for (i in 1:length(Data.Subset)){
        
        data<-Data.Subset[[i]]
        tree.type<-c(tree.type,data$subset.no[i])
        
        #Edit "species" so that it's just genus name - take first element in the string
        data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
        
        rownames(data)<-data$genus 
        
        names(data)[names(data) == "genus"] <- "Species"
        
        pgls_size_data <- data[match(pruned.tree.size$tip.label,rownames(data)),]
        
        
        #PIC for chamber size  
        #________________________________
        #Make PIC data
        pic.chsize<-pic(pgls_size_data$MeanW_S,pruned.tree.size)
        pic.size<-pic(log(pgls_size_data$AvgColSize),pruned.tree.size)
        
        data_all=cbind(pic.chsize,pic.size)
        data_pic_log=as.data.frame(data_all)
        
        # Run linear Model
        fit_L=lm(pic.chsize~pic.size-1,data=data_pic_log) #model has no intercept
        chsize.pic.aic<-c(chsize.pic.aic,AIC(fit_L))
        
        ## get CIs based on lm
        ci_lm=confint(fit_L, level = 0.95)
        #confidence interval for beta
        #ci_lm_beta_low=ci_lm[1,1] 
        #ci_lm_beta_high=ci_lm[1,2]
        
        # calculate CI for the lm estimates using a bootstrap:
        iter=1000  # set numer of iterations
        
        # set up empty vectors to fill with bootstrap data:
        #boot_interc=rep(NA, iter)
        boot_beta_L=rep(NA, iter)
        boot_R2_L=rep(NA,iter)
        
        for (i in 1:iter){
          # sample data with replacement: (note: need to shuffle with replacement
          #the relationship, not each variable)
          boot_ix=sample(dim(data_pic_log)[1], size=dim(data_pic_log)[1], replace=TRUE) # create a shuffled index with replacement
          # compute lm of the resample:
          boot_fit_L=lm(data_pic_log$pic.chsize[boot_ix]~data_pic_log$pic.size[boot_ix]-1) #Runs the model with newly sampled data-sets
          # save the bootstrapped estimates:
          boot_beta_L[i]=boot_fit_L$coefficients[1]
          #save the bootstrapped R2:
          boot_R2_L[i]=summary(boot_fit_L)$r.squared
        }
        
        ## compute 95% CI for the bootstrap:
        ci_boot_beta_high=quantile(boot_beta_L, 0.975)
        ci_boot_beta_low=quantile(boot_beta_L, 0.025)
        
        ## compute 95% CI for r2 from the bootstrap:
        ci_boot_R2_high=quantile(boot_R2_L, 0.975)
        ci_boot_R2_low=quantile(boot_R2_L, 0.025)
        
        #Save results - NEEDS TO BE ADDED ABOVE AND BELOW
        chsize.pic.coef<-c(chsize.pic.coef,fit_L$coefficients[1])
        chsize.pic.lci<-c(chsize.pic.lci,ci_boot_beta_low)
        chsize.pic.uci<-c(chsize.pic.uci,ci_boot_beta_high)
        chsize.pic.r2<-c(chsize.pic.r2,summary(fit_L)$r.squared)
        chsize.pic.rlci<-c(chsize.pic.rlci,ci_boot_R2_low)
        chsize.pic.ruci<-c(chsize.pic.ruci,ci_boot_R2_high)
        chsize.pic.p<-c(chsize.pic.p,summary(fit_L)$coefficients[1,4])
        chsize.pic.ss<-c(chsize.pic.ss, length(resid(fit_L)))
        
      }
      
      chsize.pic.all.nonlinear<-data.frame(tree.type, chsize.pic.coef,
                                        chsize.pic.lci,
                                        chsize.pic.uci,
                                        chsize.pic.r2,
                                        chsize.pic.rlci,
                                        chsize.pic.ruci,
                                        chsize.pic.p,
                                        chsize.pic.ss,
                                        chsize.pic.aic)
      
      setwd(....)
      write.csv(chsize.pic.all.nonlinear,"LOG-NL_ChamberSize_PIC.csv")
      
   
    
  library(ggplot2)
      
  #*PLOTTING - chamber size----
  #cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  cbp1<-c("darkgoldenrod2","ivory4","magenta4","springgreen3","dodgerblue", "chocolate1")
      
  x=log(mean_trait_size_sub1$AvgColSize)
  mean_trait_size_sub1$SEW_S[is.na(mean_trait_size_sub1$SEW_S)]<-0  #set NA values to 0
 
  
  #Add Family
  NameKey<-read.csv("...filename to species name key.csv")
  names(NameKey)[2]<-"species"
  key<-unique(NameKey[,c(2,3)])
  mean_trait_size_sub1<-merge(x=mean_trait_size_sub1,y=key,by="species",all.x=TRUE,no.dups=TRUE)
  str(mean_trait_size_sub1)
  
  chambsize1 <- gls(MeanW_S~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE),na.action=na.omit)
  #Model coeff of pgls lambda = 1, not sig
  # (Intercept) log(AvgColSize) 
  #  1.5373073       0.2699538 
  
  par(mfrow=c(1,1))
  #Data needs to be fixed and lines to fit need to be added to num
  ggplot(mean_trait_size_sub1, aes(x=x, y=MeanW_S, group=Family.y, color=Family.y)) + 
  geom_pointrange(aes(ymin=MeanW_S-SEW_S, ymax=MeanW_S+SEW_S))+
    #scale_color_brewer(palette=pal)+
    scale_color_manual(values = cbp1)+
    theme_bw()+
    ggtitle("Mean Chamber Size") +
    xlab("Log[ Colony Size ]") + ylab("Mean Standardized Chamber Width +/- SE")+
    #geom_line(data=predicted_df, aes(x=x, y=eff_pred),linetype="dotted")
    geom_abline(aes(intercept =  1.5373073  , slope = 0.26995382  ,lty="Lambda = 1"),col="black",linetype="dashed")+
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),text = element_text(size=20))

  
 
#CHAMBER NUMBER V COLONY SIZE----
  
  #* Chamber Number~Col Size
  #_____________________________________
  #
  #*Preparing Tree----
  #_________________
  setwd("...")
  ant.tree.moreau<-read.nexus("MoreauTree2016.nex")
  plot(ant.tree.moreau)
  str(ant.tree.moreau)
  ant.tree.moreau$tip.label
  
  #Prune tree to the genera in my dataset - Should be 24 genera
  tips<-c("Acromyrmex_versicolor","Aphaenogaster_occidentalis_NW","Camponotus_maritimus","Diacamma_rugosum","Dinoponera_australis","Dorymyrmex_bicolor","Ectatomma_opaciventre","Forelius_pruinosus","Formica_moki","Monomorium_pharaonis","Mycetagroicus_triangularis","Mycetarotes_acutus","Mycetophylax_conformis","Mycocepurus_goeldii","Myrmecocystus_flaviceps","Myrmecia_pyriformis","Odontomachus_coquereli","Pachycondyla_harpax","Pheidole_longispinosa","Pogonomyrmex_angustus", "Prenolepis_imparis" ,"Sericomyrmex_Sp","Trachymyrmex_arizonensis","Veromessor_andrei")
  pruned.tree.size<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])
  plot(pruned.tree.size)
  pruned.tree.size$tip.label
  
  #Rename tip labels to create genus level tree
  for (i in 1:length(pruned.tree.size$tip.label)) {
    split.tips<-strsplit(pruned.tree.size$tip.label[i],"_")
    genus.tip<-split.tips[[1]]
    pruned.tree.size$tip.label[i]<-genus.tip[1]
  }
  plot(pruned.tree.size) #Has 24 tips
  

  #*Prepare Trait Data----
  #______________________
  mean_trait_size2<-aggregate(DataC[],list(DataC$Name),mean)
  #****NOTE!! T. holmgreni  chamber number data needs to be added manually??**** 
  #*Avg 4.2 chmabers per nest for n=50 nests
  mean_trait_size2$ch.num[mean_trait_size2$Group.1=="Trachymyrmex holmgreni"]<-4.2
  
  # Rename Genus species column name to "species"
  colnames(mean_trait_size2)[colnames(mean_trait_size2)=="Group.1"] <- "species"
  
  
  #*Create data subsets ----
  #_______________________________________
  #Must be removed
  mean_trait_chnum_sub0<-subset(mean_trait_size2,species !="Aphaenogaster"& species != "Pheidole oxyops" & species !="Acromyrmex subterraneus")
  
  #Optional removal for diffeerent data subsets: 17 subsets
  mean_trait_chnum_sub1<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"  & species !="Formica subaenescens" )
  mean_trait_chnum_sub2<-subset(mean_trait_chnum_sub0,species!="Acromyrmex balzani" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii"  & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub3<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster floridana" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii"  & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"  & species !="Formica subaenescens"  )
  mean_trait_chnum_sub4<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster floridana" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub5<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera quadriceps" & species!="Mycocepurus smithii"  & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub6<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus goeldii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens")
  mean_trait_chnum_sub7<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus brunneus" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub8<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes parallelus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub9<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub10<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub11<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex mayri" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub12<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub13<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub14<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex mayri" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub15<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" ) #japonica
  mean_trait_chnum_sub16<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica japonica" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" ) #archboldi
  mean_trait_chnum_sub17<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica japonica" & species!="Formica pallidefulva"   & species !="Formica subaenescens" ) #dolosa
  mean_trait_chnum_sub18<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica japonica"   & species !="Formica subaenescens" )       #pallidefulva
  mean_trait_chnum_sub19<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"  & species !="Formica japonica" ) #subaenescens
  
  
  
  mean_trait_chnum_sub1$subset.no=rep(1,nrow(mean_trait_chnum_sub1))
  mean_trait_chnum_sub2$subset.no=rep(2,nrow(mean_trait_chnum_sub2))
  mean_trait_chnum_sub3$subset.no=rep(3,nrow(mean_trait_chnum_sub3))
  mean_trait_chnum_sub4$subset.no=rep(4,nrow(mean_trait_chnum_sub4))
  mean_trait_chnum_sub5$subset.no=rep(5,nrow(mean_trait_chnum_sub5))
  mean_trait_chnum_sub6$subset.no=rep(6,nrow(mean_trait_chnum_sub6))
  mean_trait_chnum_sub7$subset.no=rep(7,nrow(mean_trait_chnum_sub7))
  mean_trait_chnum_sub8$subset.no=rep(8,nrow(mean_trait_chnum_sub8))
  mean_trait_chnum_sub9$subset.no=rep(9,nrow(mean_trait_chnum_sub9))
  mean_trait_chnum_sub10$subset.no=rep(10,nrow(mean_trait_chnum_sub10))
  mean_trait_chnum_sub11$subset.no=rep(11,nrow(mean_trait_chnum_sub11))
  mean_trait_chnum_sub12$subset.no=rep(12,nrow(mean_trait_chnum_sub12))
  mean_trait_chnum_sub13$subset.no=rep(13,nrow(mean_trait_chnum_sub13))
  mean_trait_chnum_sub14$subset.no=rep(14,nrow(mean_trait_chnum_sub14))
  mean_trait_chnum_sub15$subset.no=rep(15,nrow(mean_trait_chnum_sub15))
  mean_trait_chnum_sub16$subset.no=rep(16,nrow(mean_trait_chnum_sub16))
  mean_trait_chnum_sub17$subset.no=rep(17,nrow(mean_trait_chnum_sub17))
  mean_trait_chnum_sub18$subset.no=rep(18,nrow(mean_trait_chnum_sub18))
  mean_trait_chnum_sub19$subset.no=rep(19,nrow(mean_trait_chnum_sub19))
  
  
  
  Data.Subset<-list(mean_trait_chnum_sub1,
                    mean_trait_chnum_sub2,
                    mean_trait_chnum_sub3,
                    mean_trait_chnum_sub4,
                    mean_trait_chnum_sub5,
                    mean_trait_chnum_sub6, 
                    mean_trait_chnum_sub7,
                    mean_trait_chnum_sub8,
                    mean_trait_chnum_sub9,
                    mean_trait_chnum_sub10,
                    mean_trait_chnum_sub11,
                    mean_trait_chnum_sub12,
                    mean_trait_chnum_sub13,
                    mean_trait_chnum_sub14,
                    mean_trait_chnum_sub15,
                    mean_trait_chnum_sub16,
                    mean_trait_chnum_sub17,
                    mean_trait_chnum_sub18,
                    mean_trait_chnum_sub19)
  
  
  
  #Initialize output lists
  #_________________________
  Num1.int<-c()
  Num1.coef<-c()
  Num1.se<-c()
  Num1.p<-c()
  Num1.aic<-c()
  Num1.ll<-c()
  Num1.st<-c()
  Num1.ss<-c()
  
  Num0.int<-c()
  Num0.coef<-c()
  Num0.se<-c()
  Num0.p<-c()
  Num0.aic<-c()
  Num0.ll<-c()
  Num0.st<-c()
  Num0.ss<-c()
  
  num.int<-c()
  num.coef<-c()
  num.se<-c()
  num.p<-c()
  num.aic<-c()
  num.ll<-c()
  num.r2<-c()
  num.lam<-c()
  num.lci.u<-c()
  num.lci.l<-c()
  num.normal<-c()
  num.ss<-c()
  
  chnum.pic.coef<-c()
  chnum.pic.lci<-c()
  chnum.pic.uci<-c()
  chnum.pic.r2<-c()
  chnum.pic.rlci<-c()
  chnum.pic.ruci<-c()
  chnum.pic.p<-c()
  chnum.pic.ss<-c()
  chnum.aic<-c()
  
  tree.type<-c()
  
  #Run analysis for each data subset
  #___________________________________
  #RUN ALL CHAMBER NUMBER ANALYSES
  #_____________________________________
  for (i in 1:length(Data.Subset)){
    
    data<-Data.Subset[[i]]
    tree.type<-c(tree.type,data$subset.no[1])
    
    #Edit "species" so that it's just genus name - take first element in the string
    data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
    
    rownames(data)<-data$genus 
    
    names(data)[names(data) == "genus"] <- "Species"
    
    pgls_size_data <- data[match(pruned.tree.size$tip.label,rownames(data)),]
    
    #PGLS, Lambda=1
    #_______________________
    pgls.Num.1 <- gls(ch.num~ log(AvgColSize), data=pgls_size_data, correlation=corPagel(1,pruned.tree.size, fixed = TRUE,form=~Species),na.action=na.omit)
    
    Num1.int<-c(Num1.int,summary(pgls.Num.1)$tTable[1,1])
    Num1.coef<-c(Num1.coef,summary(pgls.Num.1)$tTable[2,1])
    Num1.se<-c(Num1.se,summary(pgls.Num.1)$tTable[2,2])
    Num1.p<-c(Num1.p,summary(pgls.Num.1)$tTable[2,4])
    Num1.aic<-c(Num1.aic,summary(pgls.Num.1)$AIC)
    Num1.ll<-c(Num1.ll,summary(pgls.Num.1)$logLik)
    Num1.st<-c(Num1.st,shapiro.test(residuals(pgls.Num.1))$p.value)
    Num1.ss<-c(Num1.ss,length(resid(pgls.Num.1)))
    
    #PGLS, Lambda=0
    #________________________
    pgls.Num.0 <- gls(ch.num~ log(AvgColSize), data=pgls_size_data, correlation=corPagel(0,pruned.tree.size, fixed = TRUE,form=~Species),na.action=na.omit)
    
    Num0.int<-c(Num0.int,summary(pgls.Num.0)$tTable[1,1])
    Num0.coef<-c(Num0.coef,summary(pgls.Num.0)$tTable[2,1])
    Num0.se<-c(Num0.se,summary(pgls.Num.0)$tTable[2,2])
    Num0.p<-c(Num0.p,summary(pgls.Num.0)$tTable[2,4])
    Num0.aic<-c(Num0.aic,summary(pgls.Num.0)$AIC)
    Num0.ll<-c(Num0.ll,summary(pgls.Num.0)$logLik)
    Num0.st<-c(Num0.st,shapiro.test(residuals(pgls.Num.0))$p.value)
    Num0.ss<-c(Num0.ss,length(resid(pgls.Num.0)))
    
    
    #PIC for chamber number  
    #________________________
    #Make PIC data
    pic.chnum<-pic(pgls_size_data$ch.num,pruned.tree.size)
    pic.size<-pic(log(pgls_size_data$AvgColSize),pruned.tree.size)
    
    data_all=cbind(pic.chnum,pic.size)
    data.pic=as.data.frame(data_all)
    
    # Run linear Model
    fit_L=lm(pic.chnum~pic.size-1,data=data.pic) #model has no intercept
    chnum.aic<-c(chnum.aic,AIC(fit_L))
    
    ## get CIs based on lm
    ci_lm=confint(fit_L, level = 0.95)
    #confidence interval for beta
    #ci_lm_beta_low=ci_lm[1,1] 
    #ci_lm_beta_high=ci_lm[1,2]
    
    # calculate CI for the lm estimates using a bootstrap:
    iter=1000  # set numer of iterations
    
    # set up empty vectors to fill with bootstrap data:
    #boot_interc=rep(NA, iter)
    boot_beta_L=rep(NA, iter)
    boot_R2_L=rep(NA,iter)
    
    for (k in 1:iter){
      # sample data with replacement: (note: need to shuffle with replacement
      #the relationship, not each variable)
      boot_ix=sample(dim(data.pic)[1], size=dim(data.pic)[1], replace=TRUE) # create a shuffled index with replacement
      # compute lm of the resample:
      boot_fit_L=lm(data.pic$pic.chnum[boot_ix]~data.pic$pic.size[boot_ix]-1) #Runs the model with newly sampled data-sets
      # save the bootstrapped estimates:
      boot_beta_L[k]=boot_fit_L$coefficients[1]
      #save the bootstrapped R2:
      boot_R2_L[k]=summary(boot_fit_L)$r.squared
    }
    
    ## compute 95% CI for the bootstrap:
    ci_boot_beta_high=quantile(boot_beta_L, 0.975)
    ci_boot_beta_low=quantile(boot_beta_L, 0.025)
    
    ## compute 95% CI for r2 from the bootstrap:
    ci_boot_R2_high=quantile(boot_R2_L, 0.975)
    ci_boot_R2_low=quantile(boot_R2_L, 0.025)
    
    #Save results - NEEDS TO BE ADDED ABOVE AND BELOW
    chnum.pic.coef<-c(chnum.pic.coef,fit_L$coefficients[1])
    chnum.pic.lci<-c(chnum.pic.lci,ci_boot_beta_low)
    chnum.pic.uci<-c(chnum.pic.uci,ci_boot_beta_high)
    chnum.pic.r2<-c(chnum.pic.r2,summary(fit_L)$r.squared)
    chnum.pic.rlci<-c(chnum.pic.rlci,ci_boot_R2_low)
    chnum.pic.ruci<-c(chnum.pic.ruci,ci_boot_R2_high)
    chnum.pic.p<-c(chnum.pic.p,summary(fit_L)$coefficients[1,4])
    chnum.pic.ss<-c(chnum.pic.ss, length(resid(fit_L)))
    
  
    
  }
  
  Num.small.ss<-data.frame(tree.type,Num1.int,Num1.coef, Num1.se, Num1.p, Num1.aic, Num1.ll,Num1.st,Num1.ss,Num0.int,Num0.coef,Num0.se, Num0.p, Num0.aic, Num0.ll, Num0.st,Num0.ss)
  Num.Results.PGLS<-data.frame(tree.type,num.int,num.coef,num.se, num.p,num.aic, num.ll,num.r2, num.lam,num.lci.u,num.lci.l,num.normal,num.ss)
  chnum.pic.all.nonlinear<-data.frame(tree.type, chnum.pic.coef,
                                   chnum.pic.lci,
                                   chnum.pic.uci,
                                   chnum.pic.r2,
                                   chnum.pic.rlci,
                                   chnum.pic.ruci,
                                   chnum.pic.p,
                                   chnum.pic.ss,
                                   chnum.aic)
  
  setwd("....")
  write.csv(chnum.pic.all.nonlinear,"LOG_ChamberNum_PIC.csv")
  write.csv(Num.small.ss,"LOG_ChamberNum_lambda_0or1.csv")
  write.csv(Num.Results.PGLS,"LOG_ChamberNum_PGLSMaxLL.csv")
  

  
  
 #*PLOTTING - chamber number-----
 #*Note: save pdf as 800x500, image as 600x400
  
  #LOG
  DataC_SE<-aggregate(DataC[,"ch.num"],list(DataC$Name),StE)
  colnames(DataC_SE)<-c("species","SE_chnum")
  DataC_SE$SE_chnum[is.na(DataC_SE$SE_chnum)]<-0  #set NA values to 0
  
  Data.plot<-merge(mean_trait_size2,DataC_SE,by="species",all.x=TRUE)
  
  mean_trait_chnum_sub0<-subset(Data.plot,species !="Aphaenogaster"& species != "Pheidole oxyops" & species !="Acromyrmex subterraneus")
  
  
  #Add Family
  NameKey<-read.csv("filename to species name key.csv")
  names(NameKey)[2]<-"species"
  key<-unique(NameKey[,c(2,3)])
  mean_trait_chnum_sub0<-merge(x=mean_trait_chnum_sub0,y=key,by="species",all.x=TRUE,no.dups=TRUE)
  str(mean_trait_chnum_sub0)
  
  #Optional removal for diffeerent data subsets:
  mean_trait_chnum_sub1<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"  & species !="Formica subaenescens" )
  
  #cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#D55E00" )
  
  cbp1<-c("darkgoldenrod2","ivory4","magenta4","chocolate1","springgreen3","dodgerblue")
  
  
  x=log(mean_trait_chnum_sub1$AvgColSize)
  ggplot(mean_trait_chnum_sub1, aes(x=x, y=ch.num,group=Family.y, color=Family.y)) + 
    geom_point()+
    geom_pointrange(aes(ymin=ch.num-SE_chnum, ymax=ch.num+SE_chnum))+
    scale_color_manual(values = cbp1)+
    theme_bw()+
    ggtitle("Chamber Number") +
    xlab("Log[ Colony Size ]") + ylab("Number of Chambers per Nest")+
    geom_abline(aes(intercept = -11.54222586, slope = 4.024218372,lty="Lambda = 1"),col="black") +
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),text = element_text(size=20))
  

  
    