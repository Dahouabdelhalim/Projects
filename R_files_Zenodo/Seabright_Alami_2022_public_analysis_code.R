# library(devtools)
# install_github("edseab/Counterfact")
library(Counterfact)
library(sjPlot)
library(lme4)
library(brms)
library(dplyr)
library(ggplot2)
library(glmmTMB)


## Load in the datasets ##
dmom <- read.csv("final_analysis_table_mothers.csv")
d <- read.csv("final_analysis_table_children.csv")
social_group_plot <- read.csv("social_group_plot.csv",row.names=1)
activity_group_plot <- read.csv("activity_group_plot.csv",row.names=1)
sexbias <- read.csv("activity_sexbias.csv",row.names=1)

####################################
##                                ##
##            MODELS              ##
##                                ##
####################################

## First, the 4 models of social group size, poisson, negative binomial, and zero-inflated poisson
mm1 <- glmmTMB(out_sgroupsize ~ Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid)+ (1|Comunidad), data=dmom,family=poisson)

nb1 <- glmmTMB(out_sgroupsize ~ Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid)+ (1|Comunidad), data=dmom,family=nbinom2)

zip1nosquare <- glmmTMB(out_sgroupsize~Residence_Pattern+ age_centered + Morning   + (1|pid) + (1|Comunidad),data=dmom,ziformula= ~Residence_Pattern + age_centered + Morning + (1|pid) + (1|Comunidad),family=poisson)

zip1 <- glmmTMB(out_sgroupsize~Residence_Pattern+ age_centered +age_centered2 + Morning   + (1|pid) + (1|Comunidad),data=dmom,ziformula= ~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),family=poisson)



## Now active care and number of helpers
{
  mmcare <- glmmTMB(active_care ~ Residence_Pattern + age_centered +age_centered2 + nkids_under7 + (1|pid)+ Morning +  (1|Comunidad), data=dmom,family=binomial)
  
  mmhelpfood <- glmmTMB(help_with_food ~ Residence_Pattern + age_centered +age_centered2 + Morning +  (1|pid)+ (1|Comunidad), data=dmom,family=poisson)
  
  
  ziphelpfood <- glmmTMB(help_with_food~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),data=dmom,ziformula=~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),family=poisson)
  
  mmhelpwork <- glmmTMB(help_with_labor~ Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid)+ (1|Comunidad), data=dmom,family=poisson)
  
  ziphelpwork <- glmmTMB(help_with_labor~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),data=dmom,ziformula=~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),family=poisson)
  
  mmhelpresources <- glmmTMB(resource_acquisition_helpers~ Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid)+ (1|Comunidad), data=dmom,family=poisson)
  
  
  ziphelpresources <- glmmTMB(resource_acquisition_helpers~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),data=dmom,ziformula=~Residence_Pattern + Morning+ age_centered +age_centered2 + (1|pid) + (1|Comunidad),family=poisson)
  
  mmhelpmanufacture <- glmmTMB(help_with_manufacturing~ Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid)+ (1|Comunidad), data=dmom,family=poisson)
  
  ziphelpmanufacture <- glmmTMB(help_with_manufacturing~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),data=dmom,ziformula=~Residence_Pattern + age_centered +age_centered2 + Morning +  (1|pid) + (1|Comunidad),family=poisson)
}


## Now Activity group models
{
  mmpart <- glmmTMB(out_agroupsize~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),data=dmom,family=poisson)
  nbpart <- glmmTMB(out_agroupsize~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),data=dmom,family=nbinom2)
  
  zippart <- glmmTMB(out_agroupsize~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),data=dmom,ziformula=~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),family=poisson)
  
  mmpartfood <- glmmTMB(food_partners ~ Residence_Pattern + Morning + age_centered +age_centered2 + (1|pid)+ (1|Comunidad), data=dmom,family=poisson)
  
  zippartfood <- glmmTMB(food_partners~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),data=dmom,ziformula=~Residence_Pattern + Morning +  age_centered + (1|pid) + (1|Comunidad),family=poisson)
  
  mmpartwork <- glmmTMB(labor_partners~ Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid)+ (1|Comunidad), data=dmom,family=poisson)
  
  zippartwork <- glmmTMB(labor_partners~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),data=dmom,ziformula=~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid) + (1|Comunidad),family=poisson)
  
  mmpartresources <- glmmTMB(resource_partners~ Residence_Pattern + age_centered +age_centered2 + (1|Comunidad), data=dmom,family=poisson)
  
  zippartresources <- glmmTMB(resource_partners~ Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid)+ (1|Comunidad), data=dmom, ziformula = ~Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid)+ (1|Comunidad),family=poisson)
  
  
  mmpartmanu <- glmmTMB(manufacturing_partners~ Residence_Pattern + age_centered +age_centered2 + Morning + (1|pid)+ (1|Comunidad), data=dmom,family=poisson)
  
  zippartmanufacture <- glmmTMB(manufacturing_partners~Residence_Pattern + age_centered +age_centered2 + (1|pid) + (1|Comunidad),data=dmom,ziformula=~Residence_Pattern + age_centered +age_centered2+ (1|pid) + (1|Comunidad),family=poisson)
}

## Now childcare

allor <- glmmTMB(non_sibling_allocare ~ Residence_Pattern + age_centered + Morning + (1|pid) + (1|m_pid) +(1|Comunidad),family=binomial,data=d7)

allod <- glmmTMB(non_sibling_allocare ~ dist_to_mat_gramps + dist_to_pat_gramps + age_centered + Morning +(1|pid) + (1|m_pid) +(1|Comunidad),family=binomial,data=d7)

allorfull <- glmmTMB(all_direct_nonsib_allocare ~ Residence_Pattern + age_centered + Morning + (1|pid) + (1|m_pid) +(1|Comunidad),family=binomial,data=d7)

allodfull<- glmmTMB(all_direct_nonsib_allocare ~ log(dist_to_mat_gramps +1) + log(dist_to_pat_gramps + 1) + age_centered + Morning +(1|pid) + (1|m_pid) +(1|Comunidad),family=binomial,data=d7)

unsupr <- glmmTMB(unsupervized ~ Residence_Pattern + age_centered + Morning + (1|pid) + (1|m_pid) +(1|Comunidad),family=binomial,data=d7)

unsupd <- glmmTMB(unsupervized ~log(dist_to_mat_gramps +1) + log(dist_to_pat_gramps + 1) + age_centered + Morning + (1|pid) + (1|m_pid) +(1|Comunidad),family=binomial,data=d7)

# grpsize as function of distance to parents

dmom$log_dist_to_fam <- log(dmom$dist_to_fam + 1)
dmom$log_dist_to_inlaws <- log(dmom$dist_to_inlaws + 1)
# dmod1 <- brm(out_sgroupsize ~ log_dist_to_fam + log_dist_to_inlaws + age + (1|pid) + (1|Comunidad), data=dmom, family="poisson")
dnbhelp <- glmmTMB(out_sgroupsize ~ log_dist_to_fam + log_dist_to_inlaws + Morning + (1|pid) + (1|Comunidad), data=dmom, family=nbinom2)

dnbpart <- glmmTMB(out_agroupsize ~ log_dist_to_fam + log_dist_to_inlaws + age_centered + + Morning(1|pid) + (1|Comunidad), data=dmom, family=nbinom2)


dziphelp <- glmmTMB(out_sgroupsize ~ log_dist_to_fam + log_dist_to_inlaws + age_centered + age_centered2 +Morning + (1|pid) + (1|Comunidad),ziformula = ~ log_dist_to_fam + log_dist_to_inlaws + age_centered + age_centered2 + Morning + (1|pid) + (1|Comunidad), data=dmom, family=poisson)

dzippart <- glmmTMB(out_agroupsize ~ log_dist_to_fam + log_dist_to_inlaws + age_centered + age_centered2 +Morning + (1|pid) + (1|Comunidad),ziformula = ~ log_dist_to_fam + log_dist_to_inlaws + age_centered + age_centered2 +Morning + (1|pid) + (1|Comunidad), data=dmom, family=poisson)

####################################
##                                ##
##       graphs and figures       ##
##                                ##
####################################


# plot people in women's social groups by relatedness and sex #
  
  rel_palette <- c("#67a0d9","#9972C0","#CB6262","#8fd5c6")
  dev.new(width=16, height=14,unit="in",noRStudioGD = TRUE)
  
  barplot(as.matrix(social_group_plot), ylab="Proportion of social group partners",col=rel_palette,xlim=c(0,5),cex.names =0.9,width=(0.46),space = 0.6,cex.lab=1.2,xlab="Relation",family="serif")
  legend(3,1,rownames(social_group_plot),box.col = "white",fill=rel_palette)
  

  dev.new(width=16, height=14,unit="in",noRStudioGD = TRUE)
  
  b<-barplot(unlist(sexbias[1,]), ylab="Proportion women",col='#917cac',xlim=c(0,5),cex.names =1,ylim=c(0,1),cex.lab=1.3,space=1.0,width=(0.5),xlab="Activity type",family="serif")
  arrows(b,unlist(sexbias[1,])+qnorm(0.975)*unlist(sexbias[2,]), b,unlist(sexbias[1,])+qnorm(0.025)*unlist(sexbias[2,]) , angle=90, code=3,length=0.10)
  text(b,unlist(sexbias[1,])+0.06,labels=paste0('  ',round(unlist(sexbias[1,])*100,0),'%'))
  
  

# titles: "Sex bias in activity types across all scan observations (N= 18671)","Relatedness of socail group members by activity"

# Now activity group plots
{
  rel_palette <- c("#176BA0","#7D3AC1","#CD2323","#1de4bd")

  dev.new(width=15, height=9,unit="in",noRStudioGD = TRUE)
  barplot(partners_plot, ylab="Relatedness",col=rel_palette,xlim=c(0,5),cex.names =0.9,width=(0.6),cex.lab=1.2,xlab="Residence",space=0.6,legend=rownames(helpers_plot),family="serif",main="Partners")
}



## Code for the distance plot

# First generate counterfactual predictions along with confidence intervals (predict doesn't do this properly so this is coded into the function counterfact from package Counterfact at github.col/edseab)
s <- counterfact(dziphelp,x="log_dist_to_fam",other=list(log_dist_to_inlaws=0,Morning=0),CI=T,mixture=T)

s2 <-counterfact(dziphelp,x="log_dist_to_inlaws",other=list(log_dist_to_fam=0,Morning=0),CI=T,mixture=T)

# In order to plot the data points they need to be scaled to account for the control vartiables (essentially we want to plot residuals)

# For this we need to get both the effect on the lambda of the conditional model and the probability of the zeroinflation model
cond_scale <- exp(fixef(dziphelp)$cond["log_dist_to_inlaws"]*(dmom$log_dist_to_inlaws)+fixef(dziphelp)$cond["age_centered"]*(dmom$age_centered-mean(dmom$age_centered)) + fixef(dziphelp)$cond["age_centered2"]*(dmom$age_centered2-mean(dmom$age_centered2)) + fixef(dziphelp)$cond["Morning"]*dmom$Morning) 

odds_data <- exp(fixef(dziphelp)$zi["(Intercept)"] + fixef(dziphelp)$zi["log_dist_to_inlaws"]*(dmom$log_dist_to_inlaws)+fixef(dziphelp)$zi["age_centered"]*(dmom$age_centered) +fixef(dziphelp)$zi["age_centered2"]*(dmom$age_centered2) + fixef(dziphelp)$zi["Morning"]*dmom$Morning)



odds_baseline <- exp(fixef(dziphelp)$zi["(Intercept)"]+fixef(dziphelp)$zi["age_centered"]*mean(dmom$age_centered) + fixef(dziphelp)$zi["age_centered2"]*mean(dmom$age_centered2))

# we scale pby the probability ratio, where specifically we want (1-p) because the zi model is the probability of a structural 0 and we want the probability of not a structural 0
zi_scale <- (1-(odds_data/(1+odds_data)))/(1-(odds_baseline/(1+odds_baseline)))
scale <- cond_scale*zi_scale


resis <- dmom$out_sgroupsize/scale

# Now calculate the residuals for the second half of the plot


cond_scale2 <- exp(fixef(dziphelp)$cond["log_dist_to_fam"]*(dmom$log_dist_to_fam)+fixef(dziphelp)$cond["age_centered"]*(dmom$age_centered-mean(dmom$age_centered)) +fixef(dziphelp)$cond["age_centered2"]*(dmom$age_centered2-mean(dmom$age_centered2)) + fixef(dziphelp)$cond["Morning"]*dmom$Morning)

odds_data2 <- exp(fixef(dziphelp)$zi["(Intercept)"]+fixef(dziphelp)$zi["log_dist_to_fam"]*(dmom$log_dist_to_fam)+fixef(dziphelp)$zi["age_centered"]*(dmom$age_centered) +fixef(dziphelp)$zi["age_centered2"]*(dmom$age_centered2)+ fixef(dziphelp)$zi["Morning"]*dmom$Morning)

zi_scale2 <- (1-(odds_data2/(1+odds_data2)))/(1-(odds_baseline/(1+odds_baseline)))
scale2 <- cond_scale2*zi_scale2


resis2 <- dmom$out_sgroupsize/scale2



# Ok now we have residuals we can plot
dplot <- dev.new(noRStudioGD = T, width = 20, height = 14)
par(oma=c(0,1,0,0)+0.1,mfrow=c(2,1),mar=c(5,3,4,2)+0.1,family='serif')

plot(exp(s$xvalues),s$predicted.mean,
     xlab = "Distance from parents (m)",
     ylab = NA,
     type="l", xaxt="n",ylim=c(0,4),log="x")
axis(1,at=c(0,10,100,1000,20000,80000))

points(col=adjustcolor( "black",alpha.f = 0.3),pch=16,aggregate(dmom$dist_to_fam+1,by=list(dmom$pid),mean)[[2]],aggregate(resis,by=list(dmom$pid),mean)[[2]])


polygon(x=c(exp(s$xvalues),rev(exp(s$xvalues))),border=NA,y=c(s$LowerCI,rev(s$UpperCI)), col=adjustcolor("blue",alpha.f=0.2))

# Now the second graph doing the same for the distance to in-laws

plot(exp(s2$xvalues),s2$predicted.mean,
     ylab=NA,
     xlab = "Distance from in-laws (m)",
     type="l", xaxt="n",ylim=c(0,4),log="x")
axis(1,at=c(0,10,100,1000,20000,80000))

points(col=adjustcolor("black",alpha.f=0.3),pch=16,aggregate(dmom$dist_to_inlaws+1,by=list(dmom$pid),mean)[[2]],aggregate(resis2,by=list(dmom$pid),mean)[[2]])

polygon(x=c(exp(s2$xvalues),rev(exp(s2$xvalues))),border=NA,y=c(s2$LowerCI,rev(s2$UpperCI)), col=adjustcolor("blue",alpha.f=0.2))

mtext("Predicted social group size",side=2, outer=T, cex=1.2)


### Some more plots ###
par(family = "serif")
color_palette <- c("#2D7DD2","#97CC04","#EEB902","#F45D01","#C287E8","#7CF0BD","#E36397","#759EB8","#FAC05E","#563F1B","#058C42","#7C77B9")
activities <- prop.table(table(dmom$cat1,dmom$Residence_Pattern),margin=2)
par(family = "serif")
barplot(activities, col=color_palette,xlim=c(0,9),space=0.6,xlab="Residence",cex.lab=1.3,legend=rownames(activities),main="Relative time engaged in activities - Mothers", ylab="Proportion time observed")

par(family = "serif")

activities_kids <- prop.table(table(d7$cat1,d7$Residence_Pattern),margin=2)

barplot(activities_kids, col=color_palette,xlim=c(0,9),space=0.6,xlab="Residence",cex.lab=1.3,legend=rownames(activities_kids),main="Relative time engaged in activities - Children", ylab="Proportion time observed")

communities_residence_patterns <- prop.table(table(dmom_unique$Residence_Pattern,dmom_unique$Comunidad),margin=2)
plot.new()
barplot(communities_residence_patterns, col=color_palette[1:4],xlim=c(0,12),legend=rownames(communities_residence_patterns))


carers <- c(d7$carer1rel,d7$carer2rel,d7$carer3rel,d7$carer4rel)
carers[!is.na(carers) & !carers %in% c("mother","father","None","full sibling","aunt","uncle")] <- "other relative"
carers_t <- table(carers)[order(table(carers),decreasing=T)]
names(carers_t)[names(carers_t)=="None"] <- "unrelated"
par(family = "serif")
barplot(carers_t,col=color_palette,cex.lab=1.2,cex.names=1,main="Observed instances of direct childcare",ylab="Frequency",xlab="Carer")



barplot(table(dmom_unique$Residence_Pattern),ylim=c(0,50), col=colour_palette[1:4],xlab="Residence pattern")

dmom_unique <- dmom[!duplicated(dmom$pid),]
relatedness <- cbind(aggregate(dmom_unique$community_r,list(dmom_unique$Residence_Pattern), mean,na.rm=T),
                     aggregate(dmom_unique$community_affinal_r,list(dmom_unique$Residence_Pattern), mean, na.rm=T)[,2])
colnames(relatedness) <- c("Residence_Pattern","Mother_R","Father_R")
barplot(t(as.matrix(relatedness[,c(2,3)])),beside=T,ylab="Average relatedness to community members", col=color_palette[c(1,4)])
axis(1,at=c(2,5,8,11),labels=relatedness$Residence_Pattern)
legend(9,0.05,fill=color_palette[c(1,4)], legend=c("Mother r","Father r"))

boxplot(dmom_unique$age ~ dmom_unique$Residence_Pattern,ylab="Age",xlab="Residence")

### Average community wide relationships ###

calc_avg_pairwise_r <- function(x) {
  pids <- census_ta_data$pid[census_ta_data$origin_comid==x & census_ta_data$pid %in% colnames(r_matrix)]
  return(mean(r_matrix[pids,pids][lower.tri(r_matrix[pids,pids])],na.rm=T))}

com$average_pairwise_r <- sapply(com$comID, calc_avg_pairwise_r)
com$n_in_ta_censo <- sapply(com$comID,function(x)sum(census_ta_data$origin_comid==x))
com$n_in_r_matrix <- sapply(com$comID,function(x)sum(colnames(r_matrix) %in% census_ta_data$pid[which(census_ta_data$origin_comid==x)]))



####################################
##                                ##
##            #TABLES#            ##
##                                ##
####################################

# Table 1
t1_vars <- c("age","community_r","community_affinal_r","dist_to_fam","dist_to_inlaws")
dmom_unique <- aggregate(dmom[,t1_vars],by=list(dmom$pid),FUN=mean,na.rm=T)
dmom_unique <-  cbind(dmom_unique,dmom$Residence_Pattern[match(dmom_unique$Group.1,dmom$pid)])
colnames(dmom_unique)[c(1,ncol(dmom_unique))] <- c("pid","Residence")

aggregate(dmom$pid,list(dmom$Residence_Pattern),function(x)length(unique(x))) |>
  apply(2,unlist) |>
  cbind(aggregate(d7$pid,list(d7$Residence_Pattern),function(x)length(unique(x)))["x"])|> 
  cbind(aggregate(dmom_unique[,t1_vars],by=list(dmom_unique$Residence),FUN=mean,na.rm=T)[,-1]) |>
  cbind(sapply(levels(dmom$Residence_Pattern), function(x)sum(!is.na(dmom_unique$dist_to_fam[which(dmom_unique$Residence==x)]) & !is.na(dmom_unique$dist_to_fam[which(dmom_unique$Residence==x)])))) |>
  rbind(c("Total",length(unique(dmom$pid)),length(unique(d7$pid)),apply(dmom_unique[,c(t1_vars)],2,mean,na.rm=T),sum(!is.na(dmom_unique$dist_to_fam)&!is.na(dmom_unique$dist_to_inlaws)))) -> table1 

colnames(table1) <- c("Residence","N Women","N children <7", "Average woman's age","Mean %R in community (sd)","Mean husband's %R (sd)","Mean distance to parents in km (sd)","Mean distance to in-laws in km (sd)","N w/ complete parental GPS data")
round(as.numeric(table1$`Mean distance to in-laws in km (sd)`)/1000,2) |> paste0(' (',round(aggregate(dmom_unique[,'dist_to_inlaws'],by=list(dmom_unique$Residence),FUN=function(x)sd(x/1000,na.rm=T))$x,1), ')') -> table1$`Mean distance to in-laws in km (sd)`

round(as.numeric(table1$`Mean distance to parents in km (sd)`)/1000,1) |> paste0(' (',round(aggregate(dmom_unique[,'dist_to_fam'],by=list(dmom_unique$Residence),FUN=function(x)sd(x/1000,na.rm=T))$x,1), ')') -> table1$`Mean distance to parents in km (sd)`

round(as.numeric(table1$`Mean %R in community (sd)`)*100,1) |> paste0(' (',round(aggregate(dmom_unique[,'community_r'],by=list(dmom_unique$Residence),FUN=function(x)sd(x*100,na.rm=T))$x,1), ')') -> table1$`Mean %R in community (sd)`

round(as.numeric(table1$`Mean husband's %R (sd)`)*100,1) |> paste0(' (',round(aggregate(dmom_unique[,'community_affinal_r'],by=list(dmom_unique$Residence),FUN=function(x)sd(x*100,na.rm=T))$x,2), ')') -> table1$`Mean husband's %R (sd)`

table1$`Average woman's age` <- round(as.numeric(table1$`Average woman's age`),1)


write.csv(t(table1),"output/table_1.csv",row.names=T)


# Table S1
Model <- c('Poisson', 'Negative Binomial', 'Zero-Inflated Poisson')

compare1 <- cbind(Model,rbind(as.data.frame(anova(mm1,nb1,zip1))))
compare1$'Pr(>Chisq)'[2:3] <- '< E-20'
compare1[is.na(compare1)] <- '-'

compare2 <- cbind(Model,as.data.frame(anova(mmpart,nbpart,zippart)))
compare2$'Pr(>Chisq)'[2:3] <- '< E-20'
compare2[is.na(compare2)] <- '-'

s1 <- tab_dfs(list(compare1,compare2),titles = c('Social group','Activity group'),file='output/tableS1.html')
# Table S2
s2 <- tab_model(zip1,ziphelpresources,ziphelpmanufacture,ziphelpwork, ziphelpfood,file='output/tableS2.html')

# Table S3
s3 <- tab_model(zippart,zippartresources,zippartmanufacture,zippartwork, zippartfood,digits=1,file='output/tableS3.html')

# Table S4
s4 <- tab_model (dziphelp,dzippart,digits=1,file='output/tableS4.html')

# Table 2
t2 <- tab_model(allodfull,digits=1,file='output/table2.html')

# Table S5

s5 <- tab_model(mmcare,allorfull,unsupr,digits=1,file='output/tableS5.html')
