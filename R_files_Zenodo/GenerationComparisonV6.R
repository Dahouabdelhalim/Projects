
################# Full Model
library(MASS)
library(car)
library(emmeans)
library(rsq)
library(data.table)
library(tidyr)

setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/Phenovern3")

# This function runs the analyses for each trait. It also outputs marginal means. Both are written to files.
# The marginal means are corrected for false discovery at this step.
MakeAOVs_emmeans<-function(){
  Comb<-fread("Phenovern.combined.new.csv",data.table = FALSE)
  Comb$Tray<-paste(Comb$Generation,Comb$Tray,sep="_")
  #Comb<-Comb[-329,c(10,1:9)]
  #Comb$TID<-paste(Comb$Genotype,Comb$GrandmaternalTemp,Comb$GrandmaternalVern,sep=":")
  #table(Comb$TID,Comb$Tray,Comb$Generation)
  #Comb$TID<-NULL
  options(contrasts = c("contr.sum","contr.poly"))
  
  # Max leaf number
  Comb$MaxLeafNumber<-as.numeric(Comb$MaxLeafNumber)
  hist(log(Comb$MaxLeafNumber))
  
  fit1<-glm.nb(MaxLeafNumber~Genotype*
                        GrandmaternalTemp*
                        GrandmaternalVern*Generation,
                        data=Comb)
  #fit2<-lm(log(MaxLeafNumber)~Genotype*GrandmaternalTemp*GrandmaternalVern*Generation,data=Comb)
  #anova(fit1,fit2)
  x<-data.frame(Sum.Sq=NA,Anova(fit1,type="III"))
  x<-x[,c(1,3,2,4)]
  x<-data.frame(effect=rownames(x),x)
  #x$Pr..Chisq.<-p.adjust(x$Pr..Chisq.,method="holm")
  x$Pr..Chisq.<-round(x$Pr..Chisq.,4)
  colnames(x)<-c("effect","SumSq","Df","Fvalue/Lr.ChiSq","Pvalue")
  x$Dependent<-colnames(fit1$model)[1]
  fwrite(x,"Anovas.III.csv",quote=FALSE,col.names =TRUE)
  
  GetEMs<-function(){
    outcomes<-emmeans(fit1, pairwise ~ GrandmaternalTemp | Genotype+GrandmaternalVern+Generation)
    #outcomes$contrasts
    #x<-data.frame(outcomes$contrasts)
    #x$p.value<-round(x$p.value,3)
    #emtrends(fit1, pairwise ~ GrandmaternalTemp | Genotype+GrandmaternalVern+Generation)
    #emmip(fit1, Genotype ~ GrandmaternalTemp | GrandmaternalVern+Generation)
    #z<-data.frame(outcomes$contrasts)
    #joint_tests(fit1)
    #coef(outcomes)
    #z$contrasts
    z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
    z$contrast<-"Grandparental Temperature"
    z<-z[,-c(7:8)]
    z<-data.frame(colnames(fit1$model)[1],z)
    colnames(z)[1]<-"Phenotype"
    z$p.value<-p.adjust(z$p.value,method="fdr")
    z$p.value<-round(z$p.value,3)
    fwrite(z,"comparison.emMeans.temp.csv",quote=FALSE,row.names = FALSE,append=TRUE)
    
    outcomes<-emmeans(fit1, pairwise ~ GrandmaternalVern | Genotype+GrandmaternalTemp+Generation)
    z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
    z$contrast<-"Grandparental Vernalization"
    z<-z[,-c(7:8)]
    z<-data.frame(colnames(fit1$model)[1],z)
    colnames(z)[1]<-"Phenotype"
    z$p.value<-p.adjust(z$p.value,method="fdr")
    z$p.value<-round(z$p.value,3)
    fwrite(z,"comparison.emMeans.vern.csv",quote=FALSE,row.names = FALSE,append=TRUE)
    
    outcomes<-emmeans(fit1, pairwise ~ GrandmaternalTemp | Generation)
    z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
    z$contrast<-"Grandparental Temperature"
    #z<-z[,-c(7:8)]
    z<-data.frame(colnames(fit1$model)[1],z)
    colnames(z)[1]<-"Phenotype"
    z$p.value<-p.adjust(z$p.value,method="fdr")
    z$p.value<-round(z$p.value,3)
    fwrite(z,"comparison.emMeans.gen.csv",quote=FALSE,row.names = FALSE,append=TRUE)
    
    outcomes<-emmeans(fit1, pairwise ~ GrandmaternalVern | Generation)
    z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
    z$contrast<-"Grandparental Vernalization"
    #z<-z[,-c(7:8)]
    z<-data.frame(colnames(fit1$model)[1],z)
    colnames(z)[1]<-"Phenotype"
    z$p.value<-p.adjust(z$p.value,method="fdr")
    z$p.value<-round(z$p.value,3)
    fwrite(z,"comparison.emMeans.gen.csv",quote=FALSE,row.names = FALSE,append=TRUE)
  }
  GetEMs()
  
  # Leaf length
  
  Comb$MaxLength<-as.numeric(Comb$MaxLength)
  hist(Comb$MaxLength)
  
  #fit1<-lm(MaxLength~Genotype*GrandmaternalTemp*GrandmaternalVern*Generation,data=Comb)
  fit1<-lmer(MaxLength~Genotype*GrandmaternalTemp*GrandmaternalVern*Generation+(1|Tray),data=Comb)
  #anova(fit2,fit1)
  
  x<-data.frame(Sum.Sq=NA,Anova(fit1,type="III"))
  x<-x[,c(1,3,2,4)]
  x<-data.frame(effect=rownames(x),x)
  #x$Pr..Chisq.<-p.adjust(x$Pr..Chisq.,method="holm")
  x$Pr..Chisq.<-round(x$Pr..Chisq.,4)
  colnames(x)<-c("effect","SumSq","Df","Fvalue/Lr.ChiSq","Pvalue")
  x$Dependent<-"MaxLength"
  fwrite(x,"Anovas.III.csv",quote=FALSE,append = TRUE,col.names =FALSE)
  
  GetEMs.rnd<-function(effect){
    outcomes<-emmeans(fit1, pairwise ~ GrandmaternalTemp | Genotype+GrandmaternalVern+Generation)
    #outcomes$contrasts
    #x<-data.frame(outcomes$contrasts)
    #x$p.value<-round(x$p.value,3)
    #emtrends(fit1, pairwise ~ GrandmaternalTemp | Genotype+GrandmaternalVern+Generation)
    #emmip(fit1, Genotype ~ GrandmaternalTemp | GrandmaternalVern+Generation)
    #z<-data.frame(outcomes$contrasts)
    #joint_tests(fit1)
    #coef(outcomes)
    #z$contrasts
    z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
    z$contrast<-"Grandparental Temperature"
    z<-z[,-c(7:8)]
    z<-data.frame(effect,z)
    colnames(z)[1]<-"Phenotype"
    z$p.value<-p.adjust(z$p.value,method="fdr")
    z$p.value<-round(z$p.value,3)
    fwrite(z,"comparison.emMeans.temp.csv",quote=FALSE,row.names = FALSE,append=TRUE)
    
    outcomes<-emmeans(fit1, pairwise ~ GrandmaternalVern | Genotype+GrandmaternalTemp+Generation)
    z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
    z$contrast<-"Grandparental Vernalization"
    z<-z[,-c(7:8)]
    z<-data.frame(effect,z)
    colnames(z)[1]<-"Phenotype"
    z$p.value<-p.adjust(z$p.value,method="fdr")
    z$p.value<-round(z$p.value,3)
    fwrite(z,"comparison.emMeans.vern.csv",quote=FALSE,row.names = FALSE,append=TRUE)
    
    outcomes<-emmeans(fit1, pairwise ~ GrandmaternalTemp | Generation)
    z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
    z$contrast<-"Grandparental Temperature"
    #z<-z[,-c(7:8)]
    z<-data.frame(effect,z)
    colnames(z)[1]<-"Phenotype"
    z$p.value<-p.adjust(z$p.value,method="fdr")
    z$p.value<-round(z$p.value,3)
    fwrite(z,"comparison.emMeans.gen.csv",quote=FALSE,row.names = FALSE,append=TRUE)
    
    outcomes<-emmeans(fit1, pairwise ~ GrandmaternalVern | Generation)
    z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
    z$contrast<-"Grandparental Vernalization"
    #z<-z[,-c(7:8)]
    z<-data.frame(effect,z)
    colnames(z)[1]<-"Phenotype"
    z$p.value<-p.adjust(z$p.value,method="fdr")
    z$p.value<-round(z$p.value,3)
    fwrite(z,"comparison.emMeans.gen.csv",quote=FALSE,row.names = FALSE,append=TRUE)
  }
  GetEMs.rnd("MaxLength")
  
  # Bolting
  
  Comb$Bolting<-as.numeric(Comb$Bolting)
  hist(Comb$Bolting)
  fit1<-lmer(Bolting~Genotype*GrandmaternalTemp*GrandmaternalVern*Generation+(1|Tray),data=Comb)
  
  x<-data.frame(Sum.Sq=NA,Anova(fit1,type="III"))
  x<-x[,c(1,3,2,4)]
  x<-data.frame(effect=rownames(x),x)
  #x$Pr..Chisq.<-p.adjust(x$Pr..Chisq.,method="holm")
  x$Pr..Chisq.<-round(x$Pr..Chisq.,4)
  colnames(x)<-c("effect","SumSq","Df","Fvalue/Lr.ChiSq","Pvalue")
  x$Dependent<-"Bolting"
  fwrite(x,"Anovas.III.csv",quote=FALSE,append = TRUE,col.names =FALSE)
  
  GetEMs.rnd("Bolting")
  
  # Flower
  Comb$Flower<-as.numeric(Comb$Flower)
  hist(Comb$Flower)
  fit1<-lmer(Flower~Genotype*GrandmaternalTemp*GrandmaternalVern*Generation+(1|Tray),data=Comb)
  x<-data.frame(Sum.Sq=NA,Anova(fit1,type="III"))
  x<-x[,c(1,3,2,4)]
  x<-data.frame(effect=rownames(x),x)
  #x$Pr..Chisq.<-p.adjust(x$Pr..Chisq.,method="holm")
  x$Pr..Chisq.<-round(x$Pr..Chisq.,4)
  colnames(x)<-c("effect","SumSq","Df","Fvalue/Lr.ChiSq","Pvalue")
  x$Dependent<-"Flower"
  fwrite(x,"Anovas.III.csv",quote=FALSE,append = TRUE,col.names =FALSE)
  GetEMs.rnd("Flower")
  
  # First Fruit
  #Comb$FirstFruit<-as.numeric(Comb$FirstFruit)
  #hist(Comb$FirstFruit)
  #fit1<-lm(FirstFruit~Genotype*GrandmaternalTemp*GrandmaternalVern*Generation,data=Comb)
  #x<-data.frame(Anova(fit1,type="III"))
  #x$Pr..F.<-round(x$Pr..F.,4)
  #x$Dependent<-colnames(fit1$model)[1]
  #fwrite(x,"Anovas.III.csv",quote=FALSE,append = TRUE)
  #GetEMs()
  
  # First Yellow Fruit
  Comb$FirstYellowFruit<-as.numeric(Comb$FirstYellowFruit)
  hist(Comb$FirstYellowFruit)
  fit1<-lmer(FirstYellowFruit~Genotype*GrandmaternalTemp*GrandmaternalVern*Generation+(1|Tray),data=Comb)
  x<-data.frame(Sum.Sq=NA,Anova(fit1,type="III"))
  x<-x[,c(1,3,2,4)]
  x<-data.frame(effect=rownames(x),x)
  #x$Pr..Chisq.<-p.adjust(x$Pr..Chisq.,method="holm")
  x$Pr..Chisq.<-round(x$Pr..Chisq.,4)
  colnames(x)<-c("effect","SumSq","Df","Fvalue/Lr.ChiSq","Pvalue")
  x$Dependent<-"FirstYellowFruit"
  fwrite(x,"Anovas.III.csv",quote=FALSE,append = TRUE,col.names =FALSE)
  GetEMs.rnd("FirstYellowFruit")
  
  ### Germination
  
  library(pscl)
  library(lmerTest)
  library(lme4)
  library(car)
  library(MASS)
  library(rsq)
  library(data.table)
  setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/PhenoGerm3")
  
  Combined<-fread("GerminationData_Phenovern.csv",data.table = FALSE)
  
  library(ggplot2)
  #ggplot(Combined,aes(x=Generation,y=Max,fill=GrandmaternalTemp)) + geom_boxplot() +
  #  facet_grid(GerminationTemp~GrandmaternalVern)
  
  hist(Combined$Max)
  
  options(contrasts = c("contr.sum","contr.poly"))
  Combined$counts<-Combined$Max*20
  #hist(Combined$counts)
  fit1<-zeroinfl(counts ~ Genotype*GrandmaternalTemp*GrandmaternalVern*Generation*GerminationTemp | 1, data = Combined, dist="negbin")
  #fit1<-lm(asin(Max) ~ Genotype*GrandmaternalTemp*GrandmaternalVern*Generation*GerminationTemp, data = Combined)
  #fit1<-glm(Max ~ Genotype*GrandmaternalTemp*GrandmaternalVern*Generation*GerminationTemp, data = Combined,family="binomial")

  
  x<-data.frame(Sum.Sq=NA,Anova(fit1,type="III"))
  #x<-x[,c(1,3,2,4)]
  x<-data.frame(effect=rownames(x),x)
  #x$Pr..Chisq.<-p.adjust(x$Pr..Chisq.,method="holm")
  x$Pr..Chisq.<-round(x$Pr..Chisq.,4)
  colnames(x)<-c("effect","SumSq","Df","Fvalue/Lr.ChiSq","Pvalue")
  x$Dependent<-colnames(fit1$model)[1]
  fwrite(x,"Anovas.III.germ.csv",quote=FALSE,col.names =TRUE)
  
  library(emmeans)
  outcomes<-emmeans(fit1, pairwise ~ GrandmaternalTemp | Genotype+GrandmaternalVern+Generation+GerminationTemp)
  z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
  z$contrast<-"Grandparental Temperature"
  #z<-z[,-c(7:8)]
  z<-data.frame(colnames(fit1$model)[1],z)
  colnames(z)[1]<-"Phenotype"
  z$p.value<-p.adjust(z$p.value,method="fdr")
  z$p.value<-round(z$p.value,3)
  
  ## VERY IMPORTANT: effect sizes are reversed. fixed below:
  z$estimate<- -z$estimate
  
  fwrite(z,"comparison.gen12.emMeans.temp.csv",quote=FALSE,row.names = FALSE)
  
  outcomes<-emmeans(fit1, pairwise ~ GrandmaternalVern | Genotype+GrandmaternalTemp+Generation+GerminationTemp)
  z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
  z$contrast<-"Grandparental Vernalization"
  #z<-z[,-c(7:8)]
  z<-data.frame(colnames(fit1$model)[1],z)
  colnames(z)[1]<-"Phenotype"
  z$p.value<-p.adjust(z$p.value,method="fdr")
  z$p.value<-round(z$p.value,3)
  
  ## VERY IMPORTANT: effect sizes are reversed. fixed below:
  z$estimate<- -z$estimate
  
  fwrite(z,"comparison.gen12.emMeans.vern.csv",quote=FALSE,row.names = FALSE)
  
  outcomes<-emmeans(fit1, pairwise ~ GrandmaternalTemp | Generation)
  z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
  z$contrast<-"Grandparental Temperature"
  #z<-z[,-c(7:8)]
  z<-data.frame(colnames(fit1$model)[1],z)
  colnames(z)[1]<-"Phenotype"
  z$p.value<-p.adjust(z$p.value,method="fdr")
  z$p.value<-round(z$p.value,3)
  fwrite(z,"comparison.gen12.emMeans.gen.csv",quote=FALSE,row.names = FALSE)
  
  outcomes<-emmeans(fit1, pairwise ~ GrandmaternalVern | Generation)
  z<-data.frame(test(outcomes, side = "=", adjust = "tukey")$contrasts)
  z$contrast<-"Grandparental Vernalization"
  #z<-z[,-c(7:8)]
  z<-data.frame(colnames(fit1$model)[1],z)
  colnames(z)[1]<-"Phenotype"
  z$p.value<-p.adjust(z$p.value,method="fdr")
  z$p.value<-round(z$p.value,3)
  fwrite(z,"comparison.gen12.emMeans.gen.csv",quote=FALSE,row.names = FALSE,append=TRUE)
  
  ############## organize
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  
  
  ems<-fread("comparison.gen12.emmeans.temp.csv",data.table = FALSE)
  ems$GrandmaternalVern<-gsub("No Grandparental Vernalization","Unvernalized",ems$GrandmaternalVern)
  ems$GrandmaternalVern<-gsub("Grandparental Vernalization","Vernalized",ems$GrandmaternalVern)
  ems$Generation<-gsub("Two","Progeny",ems$Generation)
  ems$Generation<-gsub("Three","Grandprogeny",ems$Generation)
  ems$GerminationTemp<-gsub("Hot","22C",ems$GerminationTemp)
  ems$GerminationTemp<-gsub("Cold","10C",ems$GerminationTemp)
  ems$Phenotype<-"Maximum germination proportion"
  
  ems$Generation<-factor(ems$Generation,levels=c("Progeny","Grandprogeny"))
  
  ems<-ems[,-c(8:10)]
  ems$estimate<-as.numeric(ems$estimate)
  ems$estimate<-round(ems$estimate,2)
  ems<-ems[complete.cases(ems),]
  for(i in 1:nrow(ems)){
    print(i)
    if(ems[i,8] <= 0.05 & ems[i,8] > 0.01){
      ems[i,7]<-paste(ems[i,7],"*",sep="")
    } else if(ems[i,8] <= 0.01 & ems[i,8] > 0.001){
      ems[i,7]<-paste(ems[i,7],"**",sep="")
    } else if(ems[i,8] <= 0.001){
      ems[i,7]<-paste(ems[i,7],"***",sep="")
    }
  }
  ems$p.value<-NULL
  ems.temp<-ems %>% spread(key=1,7)
  
  ems<-fread("comparison.gen12.emmeans.vern.csv",data.table = FALSE)
  ems$GrandmaternalTemp<-gsub("Cold","Low",ems$GrandmaternalTemp)
  ems$GrandmaternalTemp<-gsub("Warm","High",ems$GrandmaternalTemp)
  ems$Generation<-gsub("Two","Progeny",ems$Generation)
  ems$Generation<-gsub("Three","Grandprogeny",ems$Generation)
  ems$GerminationTemp<-gsub("Hot","22C",ems$GerminationTemp)
  ems$GerminationTemp<-gsub("Cold","10C",ems$GerminationTemp)
  ems$Phenotype<-"Maximum germination proportion"
  
  #ems$Generation<-factor(ems$Generation,levels=c("Progeny","Grandprogeny"))
  
  ems<-ems[,-c(8:10)]
  ems$estimate<-as.numeric(ems$estimate)
  ems$estimate<-round(ems$estimate,2)
  ems<-ems[complete.cases(ems),]
  for(i in 1:nrow(ems)){
    print(i)
    if(ems[i,8] <= 0.05 & ems[i,8] > 0.01){
      ems[i,7]<-paste(ems[i,7],"*",sep="")
    } else if(ems[i,8] <= 0.01 & ems[i,8] > 0.001){
      ems[i,7]<-paste(ems[i,7],"**",sep="")
    } else if(ems[i,8] <= 0.001){
      ems[i,7]<-paste(ems[i,7],"***",sep="")
    }
  }
  ems$p.value<-NULL
  ems.vern<-ems %>% spread(key=1,7)
  colnames(ems.vern)<-colnames(ems.temp)
  ems<-rbind(ems.temp,ems.vern)
  colnames(ems)[3]<-"Alternate environmental cue"
  fwrite(ems,"germ.ems.org.csv",quote=FALSE,row.names = FALSE)
  
  
  
}
MakeAOVs_emmeans()

# These functions prep the type III anova tables, with stars for significance.
# These tests are also corrected for multiple testing here.


GetAOVs<-function(correction="holm"){
  library(data.table)
  library(tidyr)
  setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/Phenovern3")
  ANOVAs<-fread("Anovas.III.csv",data.table = FALSE)
  germ<-fread("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/Phenogerm3/Anovas.III.germ.csv",data.table = FALSE)
  colnames(germ)<-colnames(ANOVAs)
  ANOVAs<-rbind(ANOVAs,germ)
  rm(germ)
  
  
  # change here
  ANOVAs<-dplyr::filter(ANOVAs,!effect %in% c("(Intercept)","Residuals"))
  ANOVAs$Pvalue<-p.adjust(ANOVAs$Pvalue,method = correction)
  ANOVAs$`Fvalue/Lr.ChiSq`<-round(ANOVAs$`Fvalue/Lr.ChiSq`,2)
  
  for(i in 1:nrow(ANOVAs)){
    print(i)
    if(ANOVAs[i,5] <= 0.05 & ANOVAs[i,5] > 0.01){
      ANOVAs[i,4]<-paste(ANOVAs[i,4],"*",sep="")
    } else if(ANOVAs[i,5] <= 0.01 & ANOVAs[i,5] > 0.001){
      ANOVAs[i,4]<-paste(ANOVAs[i,4],"**",sep="")
    } else if(ANOVAs[i,5] <= 0.001){
      ANOVAs[i,4]<-paste(ANOVAs[i,4],"***",sep="")
    }
  }
  
  ANOVAs<-ANOVAs[,-c(2:3,5)]
  
  A2<-ANOVAs %>% spread(key=3,2)
  
  Order<-c(19,20,28,1,3,24,21,29,4,7,15,25,8,11,16,12)
  
  A3<-A2[c(Order),]
  
  A4<-A2[-c(Order),]
  A4$effect<-gsub(":GerminationTemp","",A4$effect)
  A4<-A4[,c(1,3)]
  colnames(A4)[2]<-"GermInt"
  
  out<-merge(A4,A3,all=TRUE,by="effect")
  out<-out[match(A3$effect, out$effect),]
  out<-out[,c(1,4,2,3,6,5,7,8)]
  
  fwrite(out,paste("Anovas.combined.",correction,".csv",sep=""),quote=FALSE,row.names = FALSE)
}
GetAOVs(correction="holm")
GetAOVs(correction="fdr")

# This function organizes the marginal means into a table for output.
EMS_organize<-function(){
  ## Organize ems
  
  setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/Phenovern3")
  ems<-fread("comparison.emmeans.temp.csv",data.table = FALSE)
  #ems$Phenotype<-gsub("FirstFruit","First fruit",ems$Phenotype)
  ems$Phenotype<-gsub("FirstYellowFruit","First mature fruit",ems$Phenotype)
  ems$Phenotype<-gsub("Flower","Flowering",ems$Phenotype)
  ems$Phenotype<-gsub("MaxLeafNumber","Leaf number at bolting",ems$Phenotype)
  ems$Phenotype<-gsub("MaxLength","Length of largest leaf at bolting",ems$Phenotype)
  #ems$Phenotype<-gsub("FruitInterval","Flowering - first yellow fruit  interval",ems$Phenotype)
  #ems$Phenotype<-gsub("FloweringInterval","Bolting - flowering interval",ems$Phenotype)
  #ems$Phenotype<-gsub("YFInterval","First fruit - first yellow fruit interval",ems$Phenotype)
  
  #ems$p.value<-p.adjust(ems$p.value,method="holm")
  
  class(ems$Phenotype)
  unique(ems$Phenotype)
  ems$Phenotype<-factor(ems$Phenotype,levels=unique(ems$Phenotype))
  ems$Generation<-factor(ems$Generation,levels=c("Progeny","Grandprogeny"))
  
  ems<-ems[,-c(7)]
  ems$estimate<-as.numeric(ems$estimate)
  ems$estimate<-round(ems$estimate,2)
  ems<-ems[complete.cases(ems),]
  for(i in 1:nrow(ems)){
    print(i)
    if(ems[i,7] <= 0.05 & ems[i,7] > 0.01){
      ems[i,6]<-paste(ems[i,6],"*",sep="")
    } else if(ems[i,7] <= 0.01 & ems[i,7] > 0.001){
      ems[i,6]<-paste(ems[i,6],"**",sep="")
    } else if(ems[i,7] <= 0.001){
      ems[i,6]<-paste(ems[i,6],"***",sep="")
    }
  }
  ems$p.value<-NULL
  library(dplyr)
  library(tidyr)
  
  ems<-ems %>% spread(key=1,6)
  
  ems2<-fread("comparison.emmeans.vern.csv",data.table = FALSE)
  #ems2$Phenotype<-gsub("FirstFruit","First fruit",ems2$Phenotype)
  ems2$Phenotype<-gsub("FirstYellowFruit","First mature fruit",ems2$Phenotype)
  ems2$Phenotype<-gsub("Flower","Flowering",ems2$Phenotype)
  ems2$Phenotype<-gsub("MaxLeafNumber","Leaf number at bolting",ems2$Phenotype)
  ems2$Phenotype<-gsub("MaxLength","Length of largest leaf at bolting",ems2$Phenotype)
  #ems2$Phenotype<-gsub("FruitInterval","Flowering - first yellow fruit  interval",ems2$Phenotype)
  #ems2$Phenotype<-gsub("FloweringInterval","Bolting - flowering interval",ems2$Phenotype)
  #ems2$Phenotype<-gsub("YFInterval","First fruit - first yellow fruit interval",ems2$Phenotype)
  
  class(ems2$Phenotype)
  unique(ems2$Phenotype)
  ems2$Phenotype<-factor(ems2$Phenotype,levels=unique(ems2$Phenotype))
  ems2$Generation<-factor(ems2$Generation,levels=c("Progeny","Grandprogeny"))
  
  ems2<-ems2[,-c(7)]
  ems2$estimate<-as.numeric(ems2$estimate)
  ems2$estimate<-round(ems2$estimate,2)
  ems2<-ems2[complete.cases(ems2),]
  for(i in 1:nrow(ems2)){
    print(i)
    if(ems2[i,7] <= 0.05 & ems2[i,7] > 0.01){
      ems2[i,6]<-paste(ems2[i,6],"*",sep="")
    } else if(ems2[i,7] <= 0.01 & ems2[i,7] > 0.001){
      ems2[i,6]<-paste(ems2[i,6],"**",sep="")
    } else if(ems2[i,7] <= 0.001){
      ems2[i,6]<-paste(ems2[i,6],"***",sep="")
    }
  }
  ems2$p.value<-NULL
  ems2<-ems2 %>% spread(key=1,6)
  
  fwrite(ems,"ems.temp.spread.csv",col.names = TRUE,row.names = FALSE,quote = FALSE)
  fwrite(ems2,"ems.vern.spread.csv",col.names = TRUE,row.names = FALSE,quote = FALSE)
  
  colnames(ems)[3]<-"Alternate environmental cue"
  colnames(ems2)<-colnames(ems)
  
  ems3<-rbind(ems,ems2)
  
  ems3$`Alternate environmental cue`<-gsub("High","Warm thermal regime",ems3$`Alternate environmental cue`)
  ems3$`Alternate environmental cue`<-gsub("Low","Cold thermal regime",ems3$`Alternate environmental cue`)
  
  
  #spread(ems3,1,6)
  fwrite(ems3,"ems.all.temp.csv",quote = FALSE,row.names = FALSE)
  
  
  # unite with germ
  
  setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/Phenovern3")
  Germ<-fread("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/Phenogerm3/germ.ems.org.csv",
              data.table = FALSE)
  
  ems<-fread("ems.all.temp.csv",data.table = FALSE)
  colnames(ems)[1]<-"Effect"
  colnames(Germ)[1:4]<-colnames(ems)[1:4]
  
  #Germ<-spread(Germ,5,6)
  Germ$Effect<-gsub("Grandparental","First-generation",Germ$Effect)
  Germ$Effect<-gsub("Temperature","thermal regime",Germ$Effect)
  Germ$Effect<-gsub("Vernalization","vernalization",Germ$Effect)
  
  ems$Effect<-gsub("Grandparental","First-generation",ems$Effect)
  ems$Effect<-gsub("Temperature","thermal regime",ems$Effect)
  ems$Effect<-gsub("Vernalization","vernalization",ems$Effect)
  
  ems$Generation<-gsub("Grandprogeny","Third",ems$Generation)
  ems$Generation<-gsub("Progeny","Second",ems$Generation)
  
  
  ems$ID<-paste(ems$Effect,ems$Genotype,ems$`Alternate environmental cue`,ems$Generation,sep="_")
  ems$ID<-gsub(" ","",ems$ID,fixed = TRUE)
  
  Germ$Generation<-gsub("Grandprogeny","Third",Germ$Generation)
  Germ$Generation<-gsub("Progeny","Second",Germ$Generation)
  Germ$`Alternate environmental cue`<-gsub("High","Warm thermal regime",Germ$`Alternate environmental cue`)
  Germ$`Alternate environmental cue`<-gsub("Low","Cold thermal regime",Germ$`Alternate environmental cue`)
  Germ$ID<-paste(Germ$Effect,Germ$Genotype,Germ$`Alternate environmental cue`,Germ$Generation,sep="_")
  Germ$ID<-gsub(" ","",Germ$ID,fixed = TRUE)
  
  Germ<-Germ[,c(7,5,6)]
  #ems<-ems[1:64,-c(5,6,14:20)]
  
  ems2<-merge(ems,Germ,by="ID",all=TRUE)
  ems2$ID<-NULL
  
  ems2<-ems2[,c(1:4,10:11,7:9,5:6)]
  
  ems2<-spread(ems2,5,6)
  ems2<-ems2[,c(1:4,10:11,5:9)]
  
  fwrite(ems2,"ems.complete.csv",quote = FALSE,row.names = FALSE)
}
EMS_organize()



### Testing for correlation between progeny and grandprogeny plasticity (persistance)
setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/Phenovern3")
library(data.table)
library(tidyr)

DataPrep<-function(){
  
  setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/Phenovern3")
  GetAndPrep<-function(file="comparison.emmeans.temp.csv"){
    ems<-fread(file,data.table = FALSE)
    #ems$Phenotype<-gsub("FirstFruit","First fruit",ems$Phenotype)
    ems$Phenotype<-gsub("FirstYellowFruit","First mature fruit",ems$Phenotype)
    ems$Phenotype<-gsub("Flower","Flowering",ems$Phenotype)
    ems$Phenotype<-gsub("MaxLeafNumber","Leaf number at bolting",ems$Phenotype)
    ems$Phenotype<-gsub("MaxLength","Length of largest leaf at bolting",ems$Phenotype)
    #ems$Phenotype<-gsub("FruitInterval","Flowering - first yellow fruit  interval",ems$Phenotype)
    #ems$Phenotype<-gsub("FloweringInterval","Bolting - flowering interval",ems$Phenotype)
    #ems$Phenotype<-gsub("YFInterval","First fruit - first yellow fruit interval",ems$Phenotype)
    
    ems$Phenotype<-factor(ems$Phenotype,levels=unique(ems$Phenotype))
    ems$Generation<-factor(ems$Generation,levels=c("Progeny","Grandprogeny"))
    
    #ems<-ems[,-c(7)]
    ems$estimate<-as.numeric(ems$estimate)
    #ems$estimate<-round(ems$estimate,2)
    ems<-ems[complete.cases(ems),]
    ems$p.value<-NULL
    
    ems.est<-ems[,-7]
    ems.est<-ems.est %>% spread(key=1,6)
    
    ems.se<-ems[,-6]
    ems.se<-ems.se %>% spread(key=1,6)
    
    return(list(ems.est,ems.se))
  }
  
  list.temp<-GetAndPrep()
  list.vern<-GetAndPrep("comparison.emmeans.vern.csv")
  
  ems.est.temp<-list.temp[[1]]
  ems.se.temp<-list.temp[[2]]
  
  ems.est.vern<-list.vern[[1]]
  ems.se.vern<-list.vern[[2]]
  
  # germ
  
  setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/PhenoGerm3")
  ems<-fread("comparison.gen12.emmeans.temp.csv",data.table = FALSE)
  ems$GrandmaternalVern<-gsub("No Grandparental Vernalization","Unvernalized",ems$GrandmaternalVern)
  ems$GrandmaternalVern<-gsub("Grandparental Vernalization","Vernalized",ems$GrandmaternalVern)
  ems$Generation<-gsub("Two","Progeny",ems$Generation)
  ems$Generation<-gsub("Three","Grandprogeny",ems$Generation)
  ems$GerminationTemp<-gsub("Hot","22C",ems$GerminationTemp)
  ems$GerminationTemp<-gsub("Cold","10C",ems$GerminationTemp)
  ems$Phenotype<-"Maximum germination proportion"
  
  ems$Generation<-factor(ems$Generation,levels=c("Progeny","Grandprogeny"))
  ems$p.value<-NULL
  ems<-ems[,-c(9:10)]
  ems$estimate<-as.numeric(ems$estimate)
  #ems$estimate<-round(ems$estimate,2)
  ems<-ems[complete.cases(ems),]
  
  germ.est<-ems[,-8]
  germ.est<-germ.est %>% spread(key=6,7)
  
  ems.est.temp$key<-paste(ems.est.temp[,1],ems.est.temp[,2],ems.est.temp[,3],ems.est.temp[,4],sep="_")
  germ.est$key<-paste(germ.est[,2],germ.est[,3],germ.est[,4],germ.est[,5],sep="_")
  germ.est<-germ.est[,c(8,6,7)]
  ems.est.temp<-merge(ems.est.temp,germ.est,by="key",all=TRUE)
  
  germ.se<-ems[,-7]
  germ.se<-germ.se %>% spread(key=6,7)
  ems.se.temp$key<-paste(ems.se.temp[,1],ems.se.temp[,2],ems.se.temp[,3],ems.se.temp[,4],sep="_")
  germ.se$key<-paste(germ.se[,2],germ.se[,3],germ.se[,4],germ.se[,5],sep="_")
  germ.se<-germ.se[,c(8,6,7)]
  ems.se.temp<-merge(ems.se.temp,germ.se,by="key",all=TRUE)
  rm(ems,germ.se,germ.est)
  
  
  ems<-fread("comparison.gen12.emmeans.vern.csv",data.table = FALSE)
  ems$GrandmaternalTemp<-gsub("Cold","Low",ems$GrandmaternalTemp)
  ems$GrandmaternalTemp<-gsub("Warm","High",ems$GrandmaternalTemp)
  ems$Generation<-gsub("Two","Progeny",ems$Generation)
  ems$Generation<-gsub("Three","Grandprogeny",ems$Generation)
  ems$GerminationTemp<-gsub("Hot","22C",ems$GerminationTemp)
  ems$GerminationTemp<-gsub("Cold","10C",ems$GerminationTemp)
  ems$Phenotype<-"Maximum germination proportion"
  
  ems$Generation<-factor(ems$Generation,levels=c("Progeny","Grandprogeny"))
  
  ems<-ems[,-c(9:10)]
  ems$estimate<-as.numeric(ems$estimate)
  #ems$estimate<-round(ems$estimate,2)
  ems<-ems[complete.cases(ems),]
  ems$p.value<-NULL
  germ.est<-ems[,-8]
  germ.est<-germ.est %>% spread(key=6,7)
  
  ems.est.vern$key<-paste(ems.est.vern[,1],ems.est.vern[,2],ems.est.vern[,3],ems.est.vern[,4],sep="_")
  germ.est$key<-paste(germ.est[,2],germ.est[,3],germ.est[,4],germ.est[,5],sep="_")
  germ.est<-germ.est[,c(8,6,7)]
  ems.est.vern<-merge(ems.est.vern,germ.est,by="key",all=TRUE)
  
  germ.se<-ems[,-7]
  germ.se<-germ.se %>% spread(key=6,7)
  ems.se.vern$key<-paste(ems.se.vern[,1],ems.se.vern[,2],ems.se.vern[,3],ems.se.vern[,4],sep="_")
  germ.se$key<-paste(germ.se[,2],germ.se[,3],germ.se[,4],germ.se[,5],sep="_")
  germ.se<-germ.se[,c(8,6,7)]
  ems.se.vern<-merge(ems.se.vern,germ.se,by="key",all=TRUE)
  rm(ems,germ.se,germ.est)
  
  
  colnames(ems.est.temp)[4]<-"Alternate environmental cue"
  colnames(ems.se.temp)[4]<-"Alternate environmental cue"
  colnames(ems.est.vern)[4]<-"Alternate environmental cue"
  colnames(ems.se.vern)[4]<-"Alternate environmental cue"
  
  ems.se<-rbind(ems.se.temp,ems.se.vern)
  ems.est<-rbind(ems.est.temp,ems.est.vern)
  
  ems.se$key<-NULL
  ems.est$key<-NULL
  
  ems.se$Generation<-gsub("Progeny","Second",ems.se$Generation)
  ems.se$Generation<-gsub("Grandprogeny","Third",ems.se$Generation)
  ems.est$Generation<-gsub("Progeny","Second",ems.est$Generation)
  ems.est$Generation<-gsub("Grandprogeny","Third",ems.est$Generation)
  
  ems.se$`Alternate environmental cue`<-gsub("Unvernalized","Not vernalized",ems.se$`Alternate environmental cue`)
  ems.est$`Alternate environmental cue`<-gsub("Unvernalized","Not vernalized",ems.est$`Alternate environmental cue`)
  
  ems.se$contrast<-gsub("Grandparental","First-generation",ems.se$contrast)
  ems.est$contrast<-gsub("Grandparental","First-generation",ems.est$contrast)
  ems.se$contrast<-gsub("Temperature","thermal regime",ems.se$contrast)
  ems.est$contrast<-gsub("Temperature","thermal regime",ems.est$contrast)
  ems.se$contrast<-gsub("Vernalization","vernalization",ems.se$contrast)
  ems.est$contrast<-gsub("Vernalization","vernalization",ems.est$contrast)
  
  ems.est$`Alternate environmental cue`<-gsub("Low","Cold thermal regime",ems.est$`Alternate environmental cue`)
  ems.est$`Alternate environmental cue`<-gsub("High","Warm thermal regime",ems.est$`Alternate environmental cue`)
  ems.se$`Alternate environmental cue`<-gsub("Low","Cold thermal regime",ems.se$`Alternate environmental cue`)
  ems.se$`Alternate environmental cue`<-gsub("High","Warm thermal regime",ems.se$`Alternate environmental cue`)
  
  
  setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/Phenovern3")
  
  return(list(ems.est,ems.se))
  
}

ems<-DataPrep()
ems.est<-ems[[1]]
ems.se<-ems[[2]]

Gen2vs3<-function(){
  Gen2<-dplyr::filter(ems.est,Generation == "Second")
  Gen3<-dplyr::filter(ems.est,!Generation == "Second")
  #Gen2se<-dplyr::filter(ems.se,Generation == "Second")
  #Gen3se<-dplyr::filter(ems.se,!Generation == "Second")
  
  out<-data.frame(Pheno=character(),Stat=numeric(),Estimate=numeric(),Pval=numeric())
  
  for(i in 5:11){
    #x<-summary(lm(Gen2[,i]~Gen3[,i]))
    #x2<-data.frame(colnames(Gen2)[i],t(x$coefficients[-1,-c(1:2)]),x$adj.r.squared)
    x<-cor.test(Gen2[,i],Gen3[,i],method="spearman")
    x2<-data.frame(colnames(Gen2)[i],x$statistic,x$estimate,x$p.value)
    colnames(x2)<-colnames(out)
    out<-rbind(x2,out)
  }
  return(out)
  
}
Gen2vs3minus<-function(){
  Gen2<-dplyr::filter(ems.est,Generation == "Second")
  Gen3<-dplyr::filter(ems.est,!Generation == "Second")
  
  out<-data.frame(Pheno=NA,Tstat=NA,Pval=NA)
  for(i in 5:9){
    x<-(cor.test(Gen2[,i],(Gen2[,i]-Gen3[,i])))
    out<-rbind(out,data.frame(Pheno=colnames(Gen2)[i],Tstat=x$statistic,Pval=x$p.value))
  }
  out<-out[-1,]
  return(out)
}



Generational<-Gen2vs3()


Generational$Pval<-p.adjust(Generational$Pval,method="holm")
max(abs(Generational$Estimate))

# developmental time test
library(car)
cor.df<-data.frame((ems.est[,c(1:4,7:9)]))
library(reshape2)
cor.df<-reshape2::melt(cor.df,id.vars=colnames(cor.df)[1:4])
cor.df$variable<- as.factor(cor.df$variable)
cor.df$variable<-as.numeric(cor.df$variable)

cor.df$value<-gsub("*","",cor.df$value,fixed = TRUE)
cor.df$value<-as.numeric(cor.df$value)
colnames(cor.df)[5:6]<-c("LifeHistoryStage","ModelEstimate")

ver<-dplyr::filter(cor.df,grepl("vernaliz",contrast))
ther<-dplyr::filter(cor.df,!grepl("vernaliz",contrast))

fit.cor<-lm(ModelEstimate~LifeHistoryStage,data=ver)
summary(fit.cor)

fit.cor<-lm(ModelEstimate~LifeHistoryStage,data=ther)
summary(fit.cor)

# instability test

setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/PhenoGerm3")
ems<-fread("comparison.gen12.emmeans.temp.csv",data.table = FALSE)
ems$GrandmaternalVern<-gsub("No Grandparental Vernalization","Unvernalized",ems$GrandmaternalVern)
ems$GrandmaternalVern<-gsub("Grandparental Vernalization","Vernalized",ems$GrandmaternalVern)
ems$Generation<-gsub("Two","Progeny",ems$Generation)
ems$Generation<-gsub("Three","Grandprogeny",ems$Generation)
ems$GerminationTemp<-gsub("Hot","22C",ems$GerminationTemp)
ems$GerminationTemp<-gsub("Cold","10C",ems$GerminationTemp)
ems$Phenotype<-"Maximum germination proportion"

ems$Generation<-factor(ems$Generation,levels=c("Progeny","Grandprogeny"))
ems$p.value<-NULL
ems<-ems[,-c(9:10)]
ems$estimate<-as.numeric(ems$estimate)
#ems$estimate<-round(ems$estimate,2)
ems<-ems[complete.cases(ems),]
germ.est<-ems[,-8]
germ.est<-germ.est %>% spread(key=6,7)

gen2<-dplyr::filter(germ.est, Generation == "Progeny")
gen3<-dplyr::filter(germ.est, !Generation == "Progeny")

cor.test(gen2$`10C`,gen2$`22C`,method="spearman")
