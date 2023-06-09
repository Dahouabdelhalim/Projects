##########################################################################################################################
#R code used in The hidden legacy of megafaunal extinction:  Loss of functional diversity and ecosystem resilience over 
#the Late Quaternary at Hall's Cave by Carson P. Hedberg, S. Kathleen Lyons, and Felisa A. Smith
#All R code written by Carson P. Hedberg, chedberg@unm.edu

##########################################################################################################################


#Install required packages
install.packages("pacman", repos="https://cloud.r-project.org")
pacman::p_load(dplyr,
               ggplot2,
               ggtree,
               FD,
               picante,
               labdsv,
               hypervolume,
               BAT,
               ade4,
               foreach,
               doParallel,
               readr,
               reshape2,
               plyr,
               cluster,
               forcats,
               funrar,
               update= F)




###########################################################################################################################
###Data cleaning and Set-up###

#upload Input data from Dryad : LINK
All_Taxa <- read_csv("Input data/All_Taxa.csv")
PA_Chart <- read_csv("Input data/PA_Chart.csv")

#rename raw data 
Allcom<- as.data.frame(All_Taxa)
#set species names as rownames 
row.names(Allcom)<- (Allcom$Binomial)
Allcom<- Allcom[order(Allcom$Binomial),]

#Set rownames in presence absence chart as time bin names
PA_Chart<- as.data.frame(PA_Chart)
row.names(PA_Chart)<- PA_Chart$Time_bin
PA_Chart$Time_bin<- NULL
#Add dummy row that includes all taxa
PA_Chart[7,]<- rep(1, dim(PA_Chart)[2])
rownames(PA_Chart)[7]<- "All"

#Isolate traits and assign variable types
traits.all<- Allcom[, c(9:17, 19:21)]
traits.all$Running<- as.ordered(traits.all$Running)
traits.all$Arboreality<- as.ordered(traits.all$Arboreality)
traits.all$Soil_Disturbance<- as.ordered(traits.all$Soil_Disturbance)
traits.all$Activity_Period<- as.factor(traits.all$Activity_Period)
traits.all$Migratory<- as.numeric(traits.all$Migratory)

#create weights vector for traits
#body size = 1
#collective diet = 1
#other collective ecological traits = 1
weights<- c(1, 1/5,1/5,1/5,1/5,1/5, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6)


###########################################################################################################################



###########################################################################################################################
###Functional Diversity Metrics###


#Run function for assessing the quality of functional space from 2 to N number of dimensions fora given trait matrix. 
#Developed by Eva Maire and Sebastien Villager (https://doi.org/10.1111/geb.12299), R code for function found in supplementary 
#information
#Testing with a maximum of 12 dimensions because there are 12 separate trait values for each species
QF<-quality_funct_space(traits.all, nbdim=12, metric="Gower", plot="quality_funct_space")
QF$meanSD
QF$details_funct_space$alg_best_tree


#create Gower's distance matrix from functional traits matrix
distALL <- gowdis(traits.all, weights, ord="podani")


#Run principle coordinate analysis on Gower's distance matrix
pcoALL <- pco(distALL, k= 12)
pcopointsALL <- pcoALL$points

#caclulate variation explained by first 4 eigenvectors 
toteigALL <- sum(pcoALL$eig[pcoALL$eig>0])
eig1ALL <-round((pcoALL$eig[1]/toteigALL)*100,digits=0) #42
eig2ALL <-round((pcoALL$eig[2]/toteigALL)*100,digits=0) #14
eig3ALL<- round((pcoALL$eig[3]/toteigALL)*100,digits=0) #8
eig4ALL<- round((pcoALL$eig[4]/toteigALL)*100,digits=0) #7

#calculate how traits load onto each PCoA axis
efG.all <- envfit(pcopointsALL, traits.all, permu=999)
efG.all


#Make empty matrix to house all FD metrics
FD.metrics<- as.data.frame(matrix(data= NA, nrow= 7, ncol= 6, dimnames= list(c("TP", "EH", "MH", 
              "LH", "Hist", "Mod", "All"),  c("FRic", "FDis", "FDiv", "FEve", "FVol", "SpeciesRichness"))))


#check to make sure taxa are in the same order 
rownames(traits.all)== colnames(PA_Chart)
#if not all TRUE, run code below to put both in alphebetical order
#traits.all<- traits.all[order(rownames(traits.all)),]
#PA_Chart<- PA_Chart[,order(colnames(PA_Chart))]


#calculate FRic, FDis, FDiv, FEve, and community-wide means (CWM's) using dbFD function from FD package
#calculates Gower's distance among taxa based on trait matix and specified weights, then performs PCoA on distance matrix 
#to construct functional space
#ord= "podani" and corr = "calliez" specify methods for treatment of ordinal varibles and negative eigenvectors, respectivley
#stand.FRic determines whether FRic values should be standardized by FRic value when all species are included
#m = x designates number of PCoA axes to use for construction of functional space, using first 4 which represent 71% of total
#variation. Higher dimensionality marginally increases quality of functional space (see QF) up to 7 dimensions, however
#higher dimensionality can lead to overfitting and increased computation time, not recommended for hypervolume methods below
#Calc CWM = TRUE and CWM.type = "all" instructs the function to calcualte community wide means for each time bin and designates
#how categorical and ordinal variabels should be treated

FD.indices<- dbFD(x=traits.all, a=PA_Chart, w=weights, ord="podani",corr="cailliez", stand.FRic = F, 
                  m= 4, print.pco= F, calc.CWM= TRUE, CWM.type = "all") 


#pull out and save community wide means (CWM)
CWM<- as.data.frame(FD.indices$CWM)
CWM$Bin<- rownames(CWM)



#Calculate functinal richness via hypervolume
#isolate points from first 4 PCoA axes 
pcoALL4<- pcoALL$points[,1:4]
#explains 71% of variance

##Separate PCO coordinates for each time bin 
PA_Chart_trans<- as.data.frame(t(PA_Chart))
pcoTP<- pcoALL4[PA_Chart_trans$TP == "1",]
pcoEH<- pcoALL4[PA_Chart_trans$EH == "1",]
pcoMH<- pcoALL4[PA_Chart_trans$MH == "1",]
pcoLH<-pcoALL4[PA_Chart_trans$LH == "1",]
pcoHist<-pcoALL4[PA_Chart_trans$Hist == "1",]
pcoMod<- pcoALL4[PA_Chart_trans$Mod == "1",]


#Calculate hypervolumes for each temporal community 
#Using gaussian method with default quantile threshod cutoff (95%). Value of hypervolume will be slightly different each time
#function is run because obtained quantile threshold is not exacly 0.95. Run code to calculate each hypervolume until obtained
#threshold is >0.9450
#Bandwidth around points is fixed based on trait axes containing full Hall's Cave assemblage. 
#All hypervoumes constructed using same set of axes and thus all have same units, making volumes comparable. 
#Samples per point allowed to fluctuate based on default gaussian methods to maintain relatively equal sampling regardless 
#of number of observations

set.seed(4)

hyp_all<-hypervolume_gaussian(pcoALL4, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
get_volume(hyp_all)

hyp_TP<-hypervolume_gaussian(pcoTP, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
get_volume(hyp_TP)

hyp_EH<-hypervolume_gaussian(pcoEH, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
get_volume(hyp_EH)

hyp_MH<-hypervolume_gaussian(pcoMH, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
get_volume(hyp_MH)

hyp_LH<-hypervolume_gaussian(pcoLH, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
get_volume(hyp_LH)

hyp_Hist<-hypervolume_gaussian(pcoHist, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
get_volume(hyp_Hist)

hyp_Mod<-hypervolume_gaussian(pcoMod, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
get_volume(hyp_Mod)


#join hypervolumes in a list
hv_list<- hypervolume_join(hyp_TP, hyp_EH, hyp_MH, hyp_LH, hyp_Hist, hyp_Mod, hyp_all)


#Add all indices to functional metrics matrix 
FD.metrics$FRic<- FD.indices$FRic
FD.metrics$FDis <- FD.indices$FDis
FD.metrics$FDiv<- FD.indices$FDiv
FD.metrics$FEve<- FD.indices$FEve
FD.metrics$FVol<- get_volume(hv_list)
FD.metrics$SpeciesRichness<- colSums(PA_Chart_trans)

#save
write.csv(FD.metrics, "FD.metrics.csv")


###########################################################################################################################



###########################################################################################################################
###Null Model###

#Set up parallel cores to run null model faster 
registerDoParallel(cores = 2) #or however many your computer will allow
getDoParWorkers() #Check how many cores are called 

##create trait matrix connected to ID numbers instead of names
traits.all.num<- traits.all
rownames(traits.all.num)<- c(1:79)

#run null model 1000 times  
n<- 1000

#set up null model to run on n cores 
results<- foreach (e = 1:n, .combine = rbind, .packages= c('picante', 'FD', 'BAT', 'hypervolume'))  %dopar%  {
  
  
  #randomize PA chart with numbers instead of species names, randomizes what species end up in each bin but keeps 
  #species richness and turnover between bins the same
  PA_Chart_num<- PA_Chart[1:6,] #don't need row with all species present
  vector<- sample(1:79, 79, replace= FALSE)
  colnames(PA_Chart_num)<- vector
  trans<- as.data.frame(t(PA_Chart_num))
  trans$vector<- as.numeric(row.names(trans))
  trans<- trans[order(trans$vector),]
  trans$vector<- NULL
  PA_Chart_rand<-t(trans)
  rownames(PA_Chart_rand)<- c("TP","EH", "MH", "LH", "Hist", "Mod")
  
 
  RFD.indicies<- dbFD(x=traits.all.num, a=PA_Chart_rand  , w=weights, ord="podani",corr="cailliez", stand.FRic = F, 
                      m= 4, print.pco= F, calc.CWM= TRUE, CWM.type = "all") 
  
  #random hypervolumes for FVol
  RpcoTP<- pcoALL4[trans$TP == "1",]
  RpcoEH<- pcoALL4[trans$EH == "1",]
  RpcoMH<- pcoALL4[trans$MH == "1",]
  RpcoLH<- pcoALL4[trans$LH == "1",]
  RpcoHist<- pcoALL4[trans$Hist == "1",]
  RpcoMod<- pcoALL4[trans$Mod == "1",]
  
  Rhyp_TP<-hypervolume_gaussian(RpcoTP, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
  Rhyp_EH<-hypervolume_gaussian(RpcoEH, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
  Rhyp_MH<-hypervolume_gaussian(RpcoMH, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
  Rhyp_LH<-hypervolume_gaussian(RpcoLH, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
  Rhyp_Hist<-hypervolume_gaussian(RpcoHist, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
  Rhyp_Mod<-hypervolume_gaussian(RpcoMod, kde.bandwidth= estimate_bandwidth(pcoALL4, method = "silverman"))
  
  
  Rhv_list<- hypervolume_join(Rhyp_TP, Rhyp_EH, Rhyp_MH, Rhyp_LH, Rhyp_Hist, Rhyp_Mod)
  
  
  #output
  data.frame(
    Bin= as.factor(c("TP","EH", "MH", "LH", "Hist", "Mod")),
    FRic = RFD.indicies$FRic,
    FDis = RFD.indicies$FDis,
    FDiv =  RFD.indicies$FDiv,
    FEve = RFD.indicies$FEve, 
    FVol= get_volume(Rhv_list),
    Log10_Mass = RFD.indicies$CWM[,1],Vertebrates= RFD.indicies$CWM[,2], Invertebrates= RFD.indicies$CWM[,3], 
    Browse= RFD.indicies$CWM[,4],Graze= RFD.indicies$CWM[,5],Fruit_Grain= RFD.indicies$CWM[,6],
    R1= RFD.indicies$CWM[,7],R2= RFD.indicies$CWM[,8],R3= RFD.indicies$CWM[,9],R4= RFD.indicies$CWM[,10],
    A1= RFD.indicies$CWM[,11],A2= RFD.indicies$CWM[,12],A3= RFD.indicies$CWM[,13],A4= RFD.indicies$CWM[,14],
    SD1= RFD.indicies$CWM[,15],SD2= RFD.indicies$CWM[,16],SD3= RFD.indicies$CWM[,17], SD4= RFD.indicies$CWM[,18],
    Log10_GroupSize= RFD.indicies$CWM[,19],Migratory0= RFD.indicies$CWM[,20], Migratory1= RFD.indicies$CWM[,21],
    Activity_Period_1= RFD.indicies$CWM[,22], Activity_Period_2= RFD.indicies$CWM[,23], 
    Activity_Period_3= RFD.indicies$CWM[,24], Activity_Period_4= RFD.indicies$CWM[,25]) 
  
}



#Format results
RFD_metrics<- results[,1:6]
R_CWM<- results[,c(1, 7:31)]

#save data 
write.csv(RFD_metrics, "RFD.metrics.csv")
write.csv(R_CWM,"R_CWM.csv")


#turn off parallel cores 
stopCluster(cl)



#melt data
RFDmetrics_long<- melt(RFD_metrics, id.vars= "Bin")
#calculate median, 95% and 50% confidence intervals for each metric and time bin
indexsum95 <- ddply(RFDmetrics_long,.(Bin, variable), summarise,Median=quantile(value, 0.5,na.rm=T),
                    CIl95=quantile(value,0.025,na.rm=T),CIu95=quantile(value, 0.975,na.rm=T))
indexsum50<- ddply(RFDmetrics_long,.(Bin, variable), summarise,Median=quantile(value, 0.5,na.rm=T),
                   CIl50=quantile(value,0.25,na.rm=T),CIu50=quantile(value, 0.75,na.rm=T))


###Create Figure 2: Summary throught time ###
###FVol
indexsum95_Fvol<- indexsum95[indexsum95$variable == "FVol",]
indexsum50_Fvol<- indexsum50[indexsum50$variable == "FVol",]

FVol_plot<- ggplot(data=indexsum95_Fvol) + geom_line(aes(x= Bin, y= Median, group = 1), color= "darkorchid3", size = 1.5)+
  scale_x_discrete(limits= c("TP", "EH", "MH","LH", "Hist", "Mod"))+
  geom_ribbon(aes(x=Bin,ymin=CIl95,ymax=CIu95,group = 1), fill = "darkorchid3", alpha=.15)+
  geom_ribbon(data= indexsum50_Fvol, aes(x=Bin,ymin=CIl50,ymax=CIu50,group = 1), fill = "darkorchid3", alpha=.20)+
  geom_line(data = FD.metrics, aes(x= rownames(FD.metrics), y= FVol, group = 1), color= "black", size = 1.5)+
  labs(x = "Time Bin", y = "FVol")+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.grid.major=element_line(color= "grey"), panel.border = element_rect(color= "grey", fill= NA))

FVol_plot

###FRic
indexsum95_FRic<- indexsum95[indexsum95$variable == "FRic",]
indexsum50_FRic<- indexsum50[indexsum50$variable == "FRic",]

FRic.plot<- ggplot(data=indexsum95_FRic) + geom_line(aes(x= Bin, y= Median, group = 1), color= "orchid3", size = 1.5)+
  scale_x_discrete(limits= c("TP", "EH", "MH","LH", "Hist", "Mod"))+
  geom_ribbon(aes(x=Bin,ymin=CIl95,ymax=CIu95,group = 1), fill = "orchid3", alpha=.20)+
  geom_ribbon(data= indexsum50_FRic, aes(x=Bin,ymin=CIl50,ymax=CIu50,group = 1), fill = "orchid3", alpha=.25)+
  geom_line(data = FD.metrics, aes(x= rownames(FD.metrics), y= FRic, group = 1), color= "black", size = 1.5)+
  labs(x = "Time Bin", y = "FRic")+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.grid.major=element_line(color= "grey"), panel.border = element_rect(color= "grey", fill= NA))

FRic.plot


###FDis
indexsum95_FDis<- indexsum95[indexsum95$variable == "FDis",]
indexsum50_FDis<- indexsum50[indexsum50$variable == "FDis",]

FDis.plot<- ggplot(data=indexsum95_FDis) + geom_line(aes(x= Bin, y= Median, group = 1), color= "deepskyblue2", size = 2)+
  scale_x_discrete(limits= c("TP", "EH", "MH","LH", "Hist", "Mod"))+
  geom_ribbon(aes(x=Bin,ymin=CIl95,ymax=CIu95,group = 1), fill = "deepskyblue2", alpha=.20)+
  geom_ribbon(data= indexsum50_FDis, aes(x=Bin,ymin=CIl50,ymax=CIu50,group = 1), fill = "deepskyblue2", alpha=.25)+
  geom_line(data = FD.metrics, aes(x= rownames(FD.metrics), y= FDis, group = 1), color= "black", size = 2)+
  labs(x = "Time Bin", y = "FDis")+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.grid.major=element_line(color= "grey"), panel.border = element_rect(color= "grey", fill= NA))
FDis.plot


###FDiv
indexsum95_FDiv<- indexsum95[indexsum95$variable == "FDiv",]
indexsum50_FDiv<- indexsum50[indexsum50$variable == "FDiv",]

FDiv.plot<- ggplot(data=indexsum95_FDiv) + geom_line(aes(x= Bin, y= Median, group = 1), color= "turquoise3", size = 2)+
  scale_x_discrete(limits= c("TP", "EH", "MH","LH", "Hist", "Mod"))+
  geom_ribbon(aes(x=Bin,ymin=CIl95,ymax=CIu95,group = 1), fill = "turquoise3", alpha=.20)+
  geom_ribbon(data= indexsum50_FDiv, aes(x=Bin,ymin=CIl50,ymax=CIu50,group = 1), fill = "turquoise3", alpha=.25)+
  geom_line(data = FD.metrics, aes(x= rownames(FD.metrics), y= FDiv, group = 1), color= "black", size = 2)+
  labs(x = "Time Bin", y = "FDiv")+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.grid.major=element_line(color= "grey"), panel.border = element_rect(color= "grey", fill= NA))

FDiv.plot


###FEve
indexsum95_FEve<- indexsum95[indexsum95$variable == "FEve",]
indexsum50_FEve<- indexsum50[indexsum50$variable == "FEve",]

FEve.plot<- ggplot(data=indexsum95_FEve) + geom_line(aes(x= Bin, y= Median, group = 1), color= "seagreen3", size = 2)+
  scale_x_discrete(limits= c("TP", "EH", "MH","LH", "Hist", "Mod"))+
  geom_ribbon(aes(x=Bin,ymin=CIl95,ymax=CIu95,group = 1), fill = "seagreen3", alpha=.20)+
  geom_ribbon(data= indexsum50_FEve, aes(x=Bin,ymin=CIl50,ymax=CIu50,group = 1), fill = "seagreen3", alpha=.25)+
  geom_line(data = FD.metrics, aes(x= rownames(FD.metrics), y= FEve, group = 1), color= "black", size = 2)+
  labs(x = "Time Bin", y = "FEve")+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.grid.major=element_line(color= "grey"), panel.border = element_rect(color= "grey", fill= NA))

FEve.plot




###Calculate change in functional diveristy metrics (slopes) between time bins
#Create empty matrix for slopes
Slopes<- as.data.frame(matrix(data= NA,nrow= 5, ncol= 10, dimnames= list(c("TP-EH", "EH-MH", "MH-LH", "LH-Hist", "Hist-Mod"),
         c("FRic", "FDis", "FDiv", "FEve", "FVol", "FRic_Z", "FDis_Z", "FDiv_Z","FEve_Z","FVol_Z"))))

RFD_metrics_TP<- RFD_metrics[grep("TP", RFD_metrics$Bin), (2:6)]
RFD_metrics_EH<- RFD_metrics[grep("EH", RFD_metrics$Bin), (2:6)]
RFD_metrics_MH<- RFD_metrics[grep("MH", RFD_metrics$Bin), (2:6)]
RFD_metrics_LH<- RFD_metrics[grep("LH", RFD_metrics$Bin), (2:6)]
RFD_metrics_Hist<- RFD_metrics[grep("Hist", RFD_metrics$Bin), (2:6)]
RFD_metrics_Mod<- RFD_metrics[grep("Mod", RFD_metrics$Bin), (2:6)]


#Calculate empirical slopes
for(e in 1:5) { 
  for (g in 1:5){
    Slopes[g,e]<- FD.metrics[g+1,e]-FD.metrics[g,e]
    
  }
}


#write function to calculate slopes for null model values
slopes_fun<- function(data1, data2){
  output<- matrix(data= NA, nrow= dim(data1)[1], ncol= dim(data1)[2])
   for (k in 1:5) {
   output[,k]<- data2[,k]- data1[,k]
  }
  return(output)
}


R_TP_EH_Slopes<- slopes_fun(RFD_metrics_TP, RFD_metrics_EH)
R_EH_MH_Slopes<- slopes_fun(RFD_metrics_EH, RFD_metrics_MH)
R_MH_LH_Slopes<- slopes_fun(RFD_metrics_MH, RFD_metrics_LH)
R_LH_Hist_Slopes<- slopes_fun(RFD_metrics_LH, RFD_metrics_Hist)
R_Hist_Mod_Slopes<- slopes_fun(RFD_metrics_Hist, RFD_metrics_Mod)


#Calculate Z scores and add to Slopes dataframe
for (i in 1:5){
  Slopes[1,i+5]<- (Slopes[1,i]- mean(R_TP_EH_Slopes[,i], na.rm = TRUE))/sd(R_TP_EH_Slopes[,i], na.rm= TRUE)}
  
for (i in 1:5){
  Slopes[2,i+5]<- (Slopes[2,i]- mean(R_EH_MH_Slopes[,i],  na.rm = TRUE))/sd(R_EH_MH_Slopes[,i],  na.rm = TRUE)}

for (i in 1:5){
  Slopes[3,i+5]<- (Slopes[3,i]- mean(R_MH_LH_Slopes[,i], na.rm = TRUE))/sd(R_MH_LH_Slopes[,i],  na.rm = TRUE)}

for (i in 1:5){
  Slopes[4,i+5]<- (Slopes[4,i]- mean(R_LH_Hist_Slopes[,i],  na.rm = TRUE))/sd(R_LH_Hist_Slopes[,i],  na.rm = TRUE)}

for (i in 1:5){
  Slopes[5,i+5]<- (Slopes[5,i]- mean(R_Hist_Mod_Slopes[,i],  na.rm = TRUE))/sd(R_Hist_Mod_Slopes[,i],  na.rm = TRUE)}


#save
write.csv(Slopes, "Slopes_Z.csv")


###########################################################################################################################



###########################################################################################################################
###Functional Composition Metrics###


###Community Wide Means###
#For ordinal or categorical traits, CWM values are the proportion of taxa in each level or category. We calculated a 
#weighted average from the proportions in each level of the ordinal traits (cursoriality, arboreality, soil disturbance) 
#to represent each by a single value. Thus, a value closer to 1 indicates a greater prevalence of fauna with high scores 
#in a given ordinal trait. Because activity period is categorical, we could not condense this trait to a single value 
#representing the distribution across all categories, and each category is represented by a separate proportion. 

#create weighted averages for all ordinal traits
CWM.red<- CWM
CWM.red$Running_2<- CWM.red$Running_2*2
CWM.red$Running_3<- CWM.red$Running_3*3
CWM.red$Running_4<- CWM.red$Running_4*4
CWM.red$Running<- (CWM.red$Running_1 + CWM.red$Running_2+ CWM.red$Running_3 + CWM.red$Running_4)/4

CWM.red$Arboreality_2<- CWM.red$Arboreality_2*2
CWM.red$Arboreality_3<- CWM.red$Arboreality_3*3
CWM.red$Arboreality_4<- CWM.red$Arboreality_4*4
CWM.red$Arboreality<- (CWM.red$Arboreality_1+ CWM.red$Arboreality_2+ CWM.red$Arboreality_3+ CWM.red$Arboreality_4)/4

CWM.red$Soil_Disturbance_2<- CWM.red$Soil_Disturbance_2*2
CWM.red$Soil_Disturbance_3<- CWM.red$Soil_Disturbance_3*3
CWM.red$Soil_Disturbance_4<- CWM.red$Soil_Disturbance_4*4
CWM.red$Soil_Disturbance<- (CWM.red$Soil_Disturbance_1+ CWM.red$Soil_Disturbance_2+ CWM.red$Soil_Disturbance_3 + CWM.red$Soil_Disturbance_4)/4

#null out ordinal variables, leave weighted averages in their place
colnames(CWM.red)
CWM.red[7:18]<- NULL

#remove migratory=0, should be opposite proportion of migratory = 1, redundant
CWM.red$Migratory_0<- NULL
#remove column of time bin names
CWM.red$Bin<- NULL


###Compare empirical CWMs to null CWMs using Z scores

#Format random CWM matrix the same as CWM matrix
R_CWM_red<- R_CWM

R_CWM_red$R2<- R_CWM_red$R2*2
R_CWM_red$R3<- R_CWM_red$R3*3
R_CWM_red$R4<- R_CWM_red$R4*4
R_CWM_red$Running<- (R_CWM_red$R1+ R_CWM_red$R2+ R_CWM_red$R3+ R_CWM_red$R4)/4

R_CWM_red$A2<- R_CWM_red$A2*2
R_CWM_red$A3<- R_CWM_red$A3*3
R_CWM_red$A4<- R_CWM_red$A4*4
R_CWM_red$Arboreality<- (R_CWM_red$A1+ R_CWM_red$A2+ R_CWM_red$A3+ R_CWM_red$A4)/4

R_CWM_red$SD2<- R_CWM_red$SD2*2
R_CWM_red$SD3<- R_CWM_red$SD3*3
R_CWM_red$SD4<- R_CWM_red$SD4*4
R_CWM_red$Soil_Disturbance<- (R_CWM_red$SD1+ R_CWM_red$SD2+ R_CWM_red$SD3+ R_CWM_red$SD4)/4

colnames(R_CWM_red)
R_CWM_red[c(1, 8:19, 21)]<- NULL

#Check to make sure traits are in the same order(column names are not all the same)
colnames(R_CWM_red)
colnames(CWM.red)


#Separate CWM values between time bins
RCWM_TP<- as.data.frame(R_CWM_red[grep("TP", rownames(R_CWM_red)),])
RCWM_EH<- as.data.frame(R_CWM_red[grep("EH", rownames(R_CWM_red)),])
RCWM_MH<- as.data.frame(R_CWM_red[grep("MH", rownames(R_CWM_red)),])
RCWM_LH<- as.data.frame(R_CWM_red[grep("LH", rownames(R_CWM_red)),])
RCWM_His<- as.data.frame(R_CWM_red[grep("His", rownames(R_CWM_red)),])
RCWM_Mod<- as.data.frame(R_CWM_red[grep("Mod", rownames(R_CWM_red)),])


#Set up empty matrix
CWM_Z<- as.data.frame(matrix(data= NA,nrow= 6, ncol= 15, dimnames= list(c("TP", "EH", "MH", "LH", "Hist", "Mod"),
              c("Mass", "Vert", "Invert", "Browse", "Graze", "Fruit_Grain", "Group_Size","Mig1", "AP1", "AP2", "AP3", 
                "AP4", "Cursoriality","Arboreality", "Soil_Dist"))))


#Calculate Z scores for each trait in each temporal community 
for (i in 1:15){
  CWM_Z[1,i]<- ((CWM.red[1,i]- mean(RCWM_TP[,i]))/sd(RCWM_TP[,i]))}

for (i in 1:15){
  CWM_Z[2,i]<- ((CWM.red[2,i]- mean(RCWM_EH[,i]))/sd(RCWM_EH[,i]))}

for (i in 1:15){
  CWM_Z[3,i]<- ((CWM.red[3,i]- mean(RCWM_MH[,i]))/sd(RCWM_MH[,i]))}

for (i in 1:15){
  CWM_Z[4,i]<- ((CWM.red[4,i]- mean(RCWM_LH[,i]))/sd(RCWM_LH[,i]))}

for (i in 1:15){
  CWM_Z[5,i]<- ((CWM.red[5,i]- mean(RCWM_His[,i]))/sd(RCWM_His[,i]))}

for (i in 1:15){
  CWM_Z[6,i]<- ((CWM.red[6,i]- mean(RCWM_Mod[,i]))/sd(RCWM_Mod[,i]))}


#save data 
write.csv(CWM_Z, "CWM_Z.csv")


###########################################################################################################################



###########################################################################################################################
###Functional Redundancy Metrics###

###Functional Uniqueness
#Calculated via nearest taxon distance (NTD) in functional space
#use first 4 axes, matches number of dimensions used to calculate other functional metrics

#create distance matrix from pco points
dist.euc <- as.matrix(dist(pcoALL4, method="euclidean",diag=T,upper=T))

#create species lists for each temporal community 
sp_names<- Allcom[,5]
TPtaxa<- Allcom[PA_Chart_trans$TP==1,5]
EHtaxa<- Allcom[PA_Chart_trans$EH==1,5]
MHtaxa<- Allcom[PA_Chart_trans$MH==1,5]
LHtaxa<- Allcom[PA_Chart_trans$LH==1,5]
Histtaxa<- Allcom[PA_Chart_trans$Hist==1,5]
Modtaxa<- Allcom[PA_Chart_trans$Mod==1,5]


#function to calculate nearest neighbor pairs for each time bin
nn_pairs_fun<- function (species_list, dist1){
  bin_dist<- dist1[rownames(dist1) %in% species_list ,colnames(dist1) %in% species_list]
  pairs_list<- as.data.frame(matrix(data = NA, nrow= dim(bin_dist)[1], ncol=3, 
                                    dimnames = list(1:dim(bin_dist)[1], c("Pair1", "Dist1","Species"))))
  
  for(k in 1:dim(bin_dist)[2]) {
    taxadist<- sort(bin_dist[,k]) 
    pairs_list$Pair1[k]<- paste(names(taxadist[1]), names(taxadist[2]), sep= "-")
    pairs_list$Dist1[k]<- taxadist[2]
    pairs_list$Species[k]<- paste(names(taxadist[1]))}
  
  return(pairs_list)
  
}


#calculate for each time bin
Pairs.TP<- nn_pairs_fun(TPtaxa, dist.euc)
Pairs.EH<- nn_pairs_fun(EHtaxa, dist.euc)
Pairs.MH<- nn_pairs_fun(MHtaxa, dist.euc)
Pairs.LH<- nn_pairs_fun(LHtaxa, dist.euc)
Pairs.Hist<- nn_pairs_fun(Histtaxa, dist.euc)
Pairs.Mod<- nn_pairs_fun(Modtaxa, dist.euc)


#Figure 5: distribution of nearest neighbor distances in Terminal Pleistocene vs Early Holocene  and Terminal Pleistocene vs Modern

ndistTP_EH <- ggplot(data= Pairs.TP, aes(x = Dist1)) + geom_histogram(aes(y = ..density..), binwidth = 0.02, 
                                                                  color= "black", fill = "turquoise4", alpha= 0.5)+
  geom_density(color= "turquoise4", size = 2 )+
  geom_histogram(data= Pairs.EH, aes(x= Dist1, y = ..density..), binwidth = 0.02, color= "black",
                 fill = "purple3", alpha= 0.4)+
  geom_density(data= Pairs.EH, color= "purple3", size = 2 )+
  geom_vline(data= Pairs.TP, aes(xintercept = median(Dist1)), color = "turquoise4", linetype = 2, size = 2)+
  geom_vline(data= Pairs.EH, aes(xintercept = median(Dist1)), color = "purple3", linetype = 2, size = 2)+
  labs(x= "Nearest Neighbor Distance", y= "Frequency")+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.border = element_rect(color= "grey", fill= NA))+
  xlim(-0.01, 0.25)+
  ylim(0, 16)

ndistTP_EH 


ndistTP_Mod <- ggplot(data= Pairs.TP, aes(x = Dist1)) + geom_histogram(aes(y = ..density..), binwidth = 0.02, 
                                                                       color= "black", fill = "turquoise4", alpha= 0.5)+
  geom_density(color= "turquoise4", size = 2 )+
  geom_histogram(data= Pairs.Mod, aes(x= Dist1, y = ..density..), binwidth = 0.02, color= "black",
                 fill = "deeppink3", alpha= 0.4)+
  geom_density(data= Pairs.Mod, color= "deeppink3", size = 2 )+
  geom_vline(data= Pairs.TP, aes(xintercept = median(Dist1)), color = "turquoise4", linetype = 2, size = 2)+
  geom_vline(data= Pairs.Mod, aes(xintercept = median(Dist1)), color = "deeppink3", linetype = 2, size = 2)+
  labs(x= "Nearest Neighbor Distance", y= "Frequency")+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.border = element_rect(color= "grey", fill= NA))+
  xlim(-0.01, 0.25)+
  ylim(0, 16)

ndistTP_Mod 




###Functional Distinctiveness
#using distinctiveness function from funrar package, calcualtes mean functional distance for each species
#in each community
#using euclidena distance matrix based on first 4 PCoA axes to remain consistent with other metrics
PA_chart_m<- as.matrix(PA_Chart)
fdistinct<- distinctiveness(PA_chart_m, dist.euc)
fdistinct<- as.data.frame(t(fdistinct))
fdistinct$All<- NULL

#save
write.csv(fdistinct, "distinctiveness.csv")


#Plot distinctiveness values through time (Figure 4)
#make long
fdistinct_long<- melt(fdistinct)
distinct.plot<- ggplot(data= fdistinct_long, aes(x= variable, y= value, fill= variable))+ 
  geom_boxplot(alpha= 0.6)+
  geom_jitter(aes(color = variable), width = .1)+
  labs(x= "Time Bin" ,y = "Functional Distinctiveness")+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.border = element_rect(color= "grey", fill= NA))
distinct.plot


#Plot TP community separated by status (extinct/extant) (Figure 4)
distinctTP<- data.frame(distinct= fdistinct$TP, Binomial= rownames(fdistinct))
distinctTP.status<- merge(distinctTP, status, by= "Binomial")

Distinct_StatusTP<- ggplot(data= distinctTP.status, aes(x= Status, y = distinct, fill = Status))+ geom_boxplot(alpha = 0.6)+
  ylab("Functional Distinctiveness") + geom_jitter(aes(color = Status), width = .1)+ 
  scale_fill_manual(values=c("deepskyblue3","darkmagenta")) + scale_color_manual(values=c("deepskyblue3","darkmagenta"))+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.border = element_rect(color= "grey", fill= NA))
  
Distinct_StatusTP

#test for significant difference 
t.test(distinct~ Status, data= distinctTP.status)
#p= 0.008468



###########################################################################################################################



###########################################################################################################################
###Functional Diversity including Introduced Speices###

#upload data sets that includes introduced species
All_Taxa_Int <- read.csv("Input data/All_Taxa_Int.csv")
PA_Chart_Int <- read.csv("Input data/PA_Chart_Int.csv")

#Format data
Allcom.Int<- as.data.frame(All_Taxa_Int)
row.names(Allcom.Int)<- (Allcom.Int$Binomial)

PA_Chart_Int<- as.data.frame(PA_Chart_Int)
rownames(PA_Chart_Int)<- PA_Chart_Int$ï..Time_bin
PA_Chart_Int$ï..Time_bin<- NULL ####only run first time after loading spreadsheet in!
PA_Chart_Int<- PA_Chart_Int[,order(colnames(PA_Chart_Int))]#make sure columns are in alphbetical order

#Separate traits and assign variable types 
traits.all.int<- Allcom.Int[, c(9:17, 19:21)]
traits.all.int$Running<- as.ordered(traits.all.int$Running)
traits.all.int$Arboreality<- as.ordered(traits.all.int$Arboreality)
traits.all.int$Soil_Disturbance<- as.ordered(traits.all.int$Soil_Disturbance)
traits.all.int$Activity_Period<- as.factor(traits.all.int$Activity_Period)
traits.all.int$Migratory<- as.numeric(traits.all.int$Migratory)


###Compute Gower's distance matrix
distInt <- gowdis(traits.all.int, weights, ord="podani")

###Create PCoA
pcoInt <- pco(distInt, k= 12)
pcopoints_Int <- pcoInt$points
toteigInt <- sum(pcoInt$eig[pcoInt$eig>0])

eig1Int <- round((pcoInt$eig[1]/toteigInt)*100,digits=0) #42
eig2Int <- round((pcoInt$eig[2]/toteigInt)*100,digits=0) #13
eig3Int <- round((pcoInt$eig[3]/toteigInt)*100,digits=0) #8
eig4Int <- round((pcoInt$eig[4]/toteigInt)*100,digits=0) #7


#Checking loadings
efG.Int <- envfit(pcopoints_Int, traits.all.int, permu=999, choices= c(1,3))
efG.Int


#transpose PA chart
PA_Chart_Int_trans<- as.data.frame(t(PA_Chart_Int))
#pull out PCO points in first 4 dimensions for each time bin
pcoTP.Int<- as.data.frame(pcopoints_Int[PA_Chart_Int_trans$TP== "1", 1:4])
pcoTP.Int$Bin<- rep("TP", dim(pcoTP.Int)[1])
pcoEH.Int<- as.data.frame(pcopoints_Int[PA_Chart_Int_trans$EH== "1", 1:4])
pcoEH.Int$Bin<- rep("EH", dim(pcoEH.Int)[1])
pcoMH.Int<- as.data.frame(pcopoints_Int[PA_Chart_Int_trans$MH== "1", 1:4])
pcoMH.Int$Bin<- rep("MH", dim(pcoMH.Int)[1])
pcoLH.Int<- as.data.frame(pcopoints_Int[PA_Chart_Int_trans$LH== "1", 1:4])
pcoLH.Int$Bin<- rep("LH", dim(pcoLH.Int)[1])
pcoHist.Int<- as.data.frame(pcopoints_Int[PA_Chart_Int_trans$Hist== "1", 1:4])
pcoHist.Int$Bin<- rep("Hist", dim(pcoHist.Int)[1])
pcoMod.Int<- as.data.frame(pcopoints_Int[PA_Chart_Int_trans$Mod== "1", 1:4])
pcoMod.Int$Bin<- rep("Mod", dim(pcoMod.Int)[1])
pcoModI.Int<- as.data.frame(pcopoints_Int[PA_Chart_Int_trans$ModInt== "1", 1:4])
pcoModI.Int$Bin<- rep("ModI", dim(pcoModI.Int)[1])

#combine into one data frame
pcoInt.combined<- rbind(pcoTP.Int, pcoEH.Int, pcoMH.Int, pcoLH.Int, pcoHist.Int, pcoMod.Int, pcoModI.Int)
pcoInt.combined$Bin<- as.factor(pcoInt.combined$Bin)
pcoInt.combined$Bin<- fct_relevel(pcoInt.combined$Bin, "TP", "EH", "MH", "LH", "Hist", "Mod", "ModI")


#countour plot showing density taxa in first two dimensions of functional space
#iterations of this code used to generate Figure 1, Figures S2-S4
contour.all<- ggplot(data= pcoInt.combined, aes(x= V1, y= V2,))+
  geom_point(color= "black", size = 3)+
  stat_density2d_filled(aes(fill= ..ndensity..), geom= "tile", position= "identity", 
                        contour = FALSE, 
                        contour_var = "ndensity", alpha= 0.8, show.legend = FALSE )+
  stat_density2d(color= "gray80", size =0.4, bins= 10)+
  scale_fill_gradient(low= "white", high= "darkorchid4")+
  xlim(-0.45, 0.45)+
  ylim(-0.25, 0.3)+
  labs(x= "PC1 (42% of variance)", y= "PC2 (13% of variance)")+
  theme(axis.text = element_text(size= 16, color= "black"),
        axis.title= element_text(size= 18), panel.background= element_rect(fill= "white"), 
        panel.border = element_rect(color= "grey", fill= NA))+
  facet_wrap(vars(Bin), nrow=2)

contour.all



###calcualte functional diversity metrics
#Make empty matrix to house FD metrics
FD_metrics_Int<- as.data.frame(matrix(data= NA, nrow= 7, ncol= 6, dimnames= list(c("TP", "EH", "MH", 
                 "LH", "Hist", "Mod", "ModInt"),  c("FVol" ,"FRic","FDis", "FDiv", "FEve", "SpeciesRichness"))))


#check
rownames(traits.all.int)== colnames(PA_Chart_Int)
#if not, order traits matrix alphabetically 
#traits.all.int<- traits.all.int[order(rownames(traits.all.int)),]

#Calculate FRic, FDiv, FDis, FEve, RaoQ, and CWM's
FD.indices.Int<- dbFD(x=traits.all.int, a=PA_Chart_Int, w=weights, ord="podani",corr="cailliez", stand.FRic = F, 
                      m= 4, print.pco= F, calc.CWM= TRUE, CWM.type = "all") 

#save CWMs
CWM_Int<- FD.indices.Int$CWM


#calculate hypervolumes
pcoInt4<- pcoInt$points[,1:4]

hyp_TP_Int<-hypervolume_gaussian(pcoTP.Int[,1:4], kde.bandwidth= estimate_bandwidth(pcoInt4, method = "silverman"))
get_volume(hyp_TP_Int)

hyp_EH_Int<-hypervolume_gaussian(pcoEH.Int[,1:4], kde.bandwidth= estimate_bandwidth(pcoInt4, method = "silverman"))
get_volume(hyp_EH_Int) 

hyp_MH_Int<-hypervolume_gaussian(pcoMH.Int[,1:4], kde.bandwidth= estimate_bandwidth(pcoInt4, method = "silverman"))
get_volume(hyp_MH_Int)

hyp_LH_Int<-hypervolume_gaussian(pcoLH.Int[,1:4], kde.bandwidth= estimate_bandwidth(pcoInt4, method = "silverman"))
get_volume(hyp_LH_Int)

hyp_Hist_Int<-hypervolume_gaussian(pcoHist.Int[,1:4], kde.bandwidth= estimate_bandwidth(pcoInt4, method = "silverman"))
get_volume(hyp_Hist_Int)

hyp_Mod_Int<-hypervolume_gaussian(pcoMod.Int[,1:4], kde.bandwidth= estimate_bandwidth(pcoInt4, method = "silverman"))
get_volume(hyp_Mod_Int) 

hyp_ModI_Int<-hypervolume_gaussian(pcoModI.Int[,1:4],kde.bandwidth= estimate_bandwidth(pcoInt4, method = "silverman"))
get_volume(hyp_ModI_Int)


hv_list_Int<- hypervolume_join(hyp_TP_Int, hyp_EH_Int, hyp_MH_Int, hyp_LH_Int, hyp_Hist_Int,
                               hyp_Mod_Int, hyp_ModI_Int)


#Add indices to matrix 
FD_metrics_Int$FVol<- get_volume(hv_list_Int)
FD_metrics_Int$FRic<- FD.indices.Int$FRic
FD_metrics_Int$FDis <- FD.indices.Int$FDis
FD_metrics_Int$FDiv<- FD.indices.Int$FDiv
FD_metrics_Int$FEve<- FD.indices.Int$FEve
FD_metrics_Int$SpeciesRichness<- colSums(PA_Chart_Int_trans)



#calcualte beta diveristy between bins, decomposed into Brepl (shift in position) and Brich (change in volume)
beta_FD_Int <- kernel.beta(hv_list_Int)
beta_FD_Int$Brepl
beta_FD_Int$Brich
beta_FD_Int$Btotal
percent.contract.Int<- beta_FD_Int$Brich/beta_FD_Int$Btotal#percent of total beta diversity between bins attributed to 
#contraction (or expansion) of volume


#Plot showing Beta diveristy between Terminal Pleistocene community and subsequent communities, 
#including Mod + Introduced (Figure S5)
TP_Beta_Int<- data_frame(Diff= beta_FD_Int$Brich [1:6], Repl= beta_FD_Int$Brepl[1:6], 
                         Bin= as.factor(c("TP-EH", "TP-MH", "TP-LH", "TP-Hist", "TP-Mod", "TP-ModInt")))
TP_Beta_Int<- melt(TP_Beta_Int, id.var= "Bin")

TP_Beta_Int_plot<- ggplot(data= TP_Beta_Int, aes(x= Bin, y= value, fill= variable)) + geom_col(alpha= 1)+
  scale_fill_manual(values = c("turquoise3", "darkorchid4"))+
  scale_x_discrete(limits= c("TP-EH", "TP-MH","TP-LH", "TP-Hist", "TP-Mod", "TP-ModInt"))+
  scale_y_continuous(limits= c(0,0.50))+
  labs(x= "Time Bin" ,y = "Total Beta Diversity")+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.border = element_rect(color= "grey", fill= NA))
TP_Beta_Int_plot




###############Functional Compostion via Community Wide Means###################

#create weighted averages for all ordinal traits
CWM_Int_red<- CWM_Int
CWM_Int_red$Running_2<- CWM_Int_red$Running_2*2
CWM_Int_red$Running_3<- CWM_Int_red$Running_3*3
CWM_Int_red$Running_4<- CWM_Int_red$Running_4*4
CWM_Int_red$Running<- (CWM_Int_red$Running_1 + CWM_Int_red$Running_2+ CWM_Int_red$Running_3 + CWM_Int_red$Running_4)/4

CWM_Int_red$Arboreality_2<- CWM_Int_red$Arboreality_2*2
CWM_Int_red$Arboreality_3<- CWM_Int_red$Arboreality_3*3
CWM_Int_red$Arboreality_4<- CWM_Int_red$Arboreality_4*4
CWM_Int_red$Arboreality<- (CWM_Int_red$Arboreality_1+ CWM_Int_red$Arboreality_2+ CWM_Int_red$Arboreality_3+ CWM_Int_red$Arboreality_4)/4

CWM_Int_red$Soil_Disturbance_2<- CWM_Int_red$Soil_Disturbance_2*2
CWM_Int_red$Soil_Disturbance_3<- CWM_Int_red$Soil_Disturbance_3*3
CWM_Int_red$Soil_Disturbance_4<- CWM_Int_red$Soil_Disturbance_4*4
CWM_Int_red$Soil_Disturbance<- (CWM_Int_red$Soil_Disturbance_1+ CWM_Int_red$Soil_Disturbance_2+ CWM_Int_red$Soil_Disturbance_3 + CWM_Int_red$Soil_Disturbance_4)/4


#null out ordinal variables, leave weighted averages in their place
colnames(CWM_Int_red)
CWM_Int_red[7:18]<- NULL
#null out migratory=0, should be opposite proportion to migratory =1 so unnecessary 
CWM_Int_red$Migratory_0<- NULL

###Ordination using PCO to assess multivariate shifts in community composition
#using Euclidean distance because all CWM values are numeric
#each activity period category is weighted to 1/4 to give it collectively equal weight to the other CWMs 
distCWM.Int<-daisy(CWM_Int_red, metric = "euclidean", weights= c(1,1,1,1,1,1,1,1,1/4,1/4,1/4,1/4,1,1,1), stand= TRUE)
pcoCWM.Int <- pco(distCWM.Int, k= 6)
origscoreCWM.Int <- pcoCWM.Int$points
toteigCWM.Int <- sum(pcoCWM.Int$eig[pcoCWM.Int$eig>0])

eig1CWM.Int <-round((pcoCWM.Int$eig[1]/toteigCWM.Int)*100,digits=1) 
eig2CWM.Int <-round((pcoCWM.Int$eig[2]/toteigCWM.Int)*100,digits=1)
eig3CWM.Int<- round((pcoCWM.Int$eig[3]/toteigCWM.Int)*100,digits=1)

#Checking loadings
efG.cwm.int <- envfit(origscoreCWM.Int, CWM_Int_red , permu=999)
efG.cwm.int
#pull out loadings
cwm_vectors_pco<- as.data.frame(efG.cwm.int$vectors$arrows)

#pull out points from first two PCO axes
pcoCWM.Int2<- as.data.frame(origscoreCWM.Int[,1:2])
pcoCWM.Int2$Bin<- rownames(pcoCWM.Int2)


#plot communities (Figure 3)
CWM.plot.pco<- ggplot(data= pcoCWM.Int2, aes(x= V1, y= V2, label= Bin ))+ geom_point(size = 2)+
  geom_text(aes(label= Bin),hjust=1.2, vjust=1.2)+
  geom_segment(data= cwm_vectors_pco, aes(x= 0, y= 0, xend= Dim1*3.5, yend= Dim2*3.5), arrow= arrow(), 
               inherit.aes = FALSE, color = "turquoise3")+
  annotate("text", x= cwm_vectors_pco$Dim1*3.6, y= cwm_vectors_pco$Dim2*3.6, label= rownames(cwm_vectors_pco))+
  theme(axis.text = element_text(size= 23, color= "black"),
        axis.title= element_text(size= 26), panel.background= element_rect(fill= "white"), 
        panel.border = element_rect(color= "grey", fill= NA))+
  labs(x= "PC1 (63% of variance)", y= "PC2 (23% of variance)")+
  xlim(-6, 8) + ylim(-6, 4)


CWM.plot.pco


#Functional distinctiveness of Introduced species
#create distance matrix
dist.euc.int <- as.matrix(dist(pcoInt4, method="euclidean",diag=T,upper=T))

PA_chart_Int_m<- as.matrix(PA_Chart_Int)
fdistinct_int<- distinctiveness(PA_chart_Int_m, dist.euc.int)
fdistinct_int<- as.data.frame(t(fdistinct_int))

#pull out ModInt values and merge with status
Status_int<- All_Taxa_Int[,5:6]
distinctModI<- data.frame(distinct= fdistinct_int$ModInt, Binomial= rownames(fdistinct_int))
distinctModI.status<- merge(distinctModI, Status_int, by= "Binomial")

#test for significant difference 
t.test(distinct~ Status, data= distinctModI.status)
#p-value = 0.07351






##########################################################################################################################
##End of script###
#########################################################################################################################

