#Orzechowski, E.A., and S. Finnegan, 2020, Controls on range shifts of coastal Californian bivalves during the peak of the last interglacial and baseline predictions for today, Paleobiology. 

#setwd()
packages<- c("ggplot2", "reghelper", "dismo", "DescTools", "gbm", "plyr", "gridExtra", "ggrepel")
lapply(packages, require, character.only=T) 



############################
#Paleo-embayment v. paleo-open coast extralimital occurrences
############################



coasts<- read.csv("SI_coasts.csv", header = T, stringsAsFactors = F)
coasts<- coasts[,c("coast", "extra")]
coasts$extra<- as.numeric(as.character(coasts$extra))
co<- cbind(10, 3)
c<- cbind(0, 9)
coa<- rbind(co, c)
colnames(coa)<- c("non", "extra.occs")
rownames(coa)<- c("paleo-open", "paleo-bay")

##############
#G-test statistics compiled in Table 1
summary.coast<- GTest(coa)
summary.coast$statistic #chi-square test statistic
summary.coast$parameter #df
summary.coast$p.value
summary.coast$method #Log likelihood ratio (G-test) test of independence without correction
summary.coast$observed
summary.coast$expected #expected count under null hypothesis; the null hypothesis is: "G-test is performed on the null hypothesis that the joint distribution of the cell counts in a 2-dimensional contingency table is the product of the row and column marginals." Documentation for Gtest, Pete Hurd.



############################
#Range-through environmental analyses 
############################



#Stepwise regressions (range-through environmental conditions); statistics complied in Table 3
fauna_range_through <- read.csv("SI_fauna_range_through.csv",header=TRUE)
panamic.RTvar<- fauna_range_through[,c("extra", "rt_ann", "rt_sum", "rt_win", "lat.range", "n.islands", "crt_max", "crt_min")]
panamic.RTvar<- panamic.RTvar[complete.cases(panamic.RTvar), ]
fullmodel<- glm(extra~ ., data=panamic.RTvar, family=binomial(link="logit"))
nothing<- glm(extra~ 1, data=panamic.RTvar, family=binomial(link="logit"))
both<- step(nothing,scope=list(lower=formula(nothing),upper=formula(fullmodel)), direction="both")
best.model<- glm(extra~rt_ann + n.islands + crt_max + crt_min, data =panamic.RTvar, family=binomial(link="logit"))
formula(both)
summary(both) #forwards and backwards stepwise regressions
beta(best.model) #beta coefficients; listed in Table 3

##############
#Predictor (range-through environmental conditions) correlation statistics compiled in Table 2
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .5, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
corr_matrix <- abs(cor(panamic.RTvar[, -5])) <= .5
corr_matrix



############################
#Range-limiting environmental analyses 
############################



#Stepwise regressions (range-limiting environmental conditions); statistics complied in Table 3
fauna_range_limits <- read.csv("SI_fauna_range_limits.csv",header=TRUE)
panamic.RLvar<- fauna_range_limits[,c("extra", "rl_ann", "rl_sum", "rl_win", "lat.range", "n.islands", "crl_max", "crl_min")]
fullmodel<- glm(extra~ ., data=panamic.RLvar, family=binomial(link="logit"))
nothing<- glm(extra~ 1, data=panamic.RLvar, family=binomial(link="logit"))
both<- step(nothing,scope=list(lower=formula(nothing),upper=formula(fullmodel)), direction="both")
best.model<- glm(extra~rl_ann + n.islands + crl_max + crl_min, data =panamic.RLvar, family=binomial(link="logit"))
formula(both)
summary(both) #forwards and backwards stepwise regressions
beta(best.model) #beta coefficients; listed in Table 3

##############
#Predictor (range-limiting environmental conditions) correlation statistics compiled in Table 2
is.correlated <- function(i, j, data, conf.level = .95, cutoff = .5, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
vCorrelated <- Vectorize(is.correlated, c("i", "j"))
corr_matrix <- abs(cor(panamic.RLvar[, -5])) <= .5
corr_matrix

##############GBM models (range-limiting environmental conditions)
newdf <- panamic.RLvar
newdf <- newdf[c("extra","rl_ann","n.islands","crl_max","crl_min")]
mod <- gbm.step(data=newdf, gbm.x = c(2:5), gbm.y = 1,family = "bernoulli", tree.complexity = 3,learning.rate = 0.001, bag.fraction =.67,n.trees=100,max.trees=20000,n.folds=10,step.size=50,silent = FALSE,plot.main = TRUE)

partials <- list()
for(i in seq(1,4)){
 # i <- 5
pardep <- plot.gbm(mod, i.var = i, n.trees = mod$n.trees, continuous.resolution = 20, return.grid = TRUE, type="response") 
pardep$variable <- colnames(pardep)[1]
colnames(pardep) <- c("value","prediction","variable")
partials[[i]] <- pardep
}
partials <- do.call("rbind",partials)
partials$variable <- revalue(partials$variable,c("crl_max"="Chlor. Max","crl_min"="Chlor. Min","n.islands"="Num. Islands","rl_ann" = "Mean Annual Temp."))

#Figure 3: Partial dependence plots for all four predictors selected for inclusion in the GBM models. 
ggplot(partials,aes(value,prediction))+
  geom_line() +
  facet_wrap(~variable,scales="free_x",ncol=2)  +
  theme_bw() +
  xlab("Variable value") +
  ylab("Marginal effect on prob. extralimital") +
  theme(strip.background = element_rect(colour = 'white', fill = 'white'))

highest <- newdf[(newdf$rl_ann < quantile(newdf$rl_ann,.25) & newdf$crl_min > quantile(newdf$crl_min,.75) & newdf$crl_max > quantile(newdf$crl_max,.75) & newdf$n.islands > quantile(newdf$n.islands,.75)),]
lowest <- newdf[(newdf$rl_ann > quantile(newdf$rl_ann,.25) & newdf$crl_min < quantile(newdf$crl_min,.75) & newdf$crl_max < quantile(newdf$crl_max,.75) & newdf$n.islands < quantile(newdf$n.islands,.75)),]

##############
#code to generate Figure 4
predicted<- predict(mod, newdata=fauna_range_limits, type="response",n.trees = mod$n.trees)
pred.data<- cbind(predicted, fauna_range_limits)
pred.data<- pred.data[order(pred.data$predicted),]
pred.data<- pred.data[order(pred.data$predicted, decreasing=T),]
pred.data$seq1<- 1:nrow(pred.data)
pred.data$ScientificName<- pred.data$ScientificName
pred.data$ScientificName<- as.character(pred.data$ScientificName)
library(stringr)
pred.data$extra<- as.factor(pred.data$extra)
pred.data$extra<- pred.data$extra
pred.data$extra<- ifelse(pred.data$extra == 0, "N", "Y") 
pred.data<- pred.data[complete.cases(pred.data$extra),]
pred.data$predicted <-pred.data$predicted/max(pred.data$predicted)
migrations <- pred.data[(pred.data$seq1 <= 20),]
rownames(migrations) <- migrations$ScientificName
migrations <- migrations[c("predicted","extra")]
migrations$predicted <- round(migrations$predicted,2)
colnames(migrations) <- c("Predicted\\npotential","MIS 5e\\nextralimital?")

##############
#Figure 4: Predicted relative colonization potential for all species in the present-day northern Panamic species pool, based on the best-supported MIS 5e extralimital model. 

ggplot(pred.data,aes(predicted,seq1,colour=extra,label=ScientificName))+ geom_segment(aes(xend = 0,yend=seq1),size=.25) + scale_colour_manual(values=c("darkgrey", "black"),name = "MIS 5e\\nExtralimital?")+  theme_bw()  + theme(legend.title =element_text(colour="black", size=12,face="plain"),legend.text = element_text(colour="black", size = 10)) + xlab("Relative colonization potential") + theme(axis.title=element_text(size=12)) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())  + theme(legend.position=c(.2,.85)) + coord_cartesian(xlim=c(0,1.02),ylim=c(1,345),expand = FALSE)+ theme(axis.text=element_text(size=10),axis.title=element_text(size=12,face="plain")) + annotation_custom(tableGrob(migrations,theme =ttheme_minimal(base_size=6,hjust=1,padding = unit(c(1.2, 1.2), "mm"))), xmin=.45, xmax=.95, ymin=170, ymax=170) + annotate("text",x=.7,y=254,label="Relative colonization potential") + annotate("rect",xmin=0,xmax=1,ymin=1,ymax=20,fill=NA,colour="black") + annotate("rect",xmin=.42,xmax=.98,ymin=90,ymax=265,fill=NA,colour="black")+ylab("Rank order")+ annotate("segment",x=.7,xend=.7,y=90,yend=22,arrow = arrow(length = unit(0.2,"cm")))

geom_label_repel(data=migrations,aes(x=predicted,y=seq1,label=ScientificName),alpha=1,fill="white",fontface = 'bold',size=2,box.padding = unit(.2, "lines"),point.padding = unit(0, "lines"),segment.color = 'gray',force=.5,nudge_x = 1, hjust=0,direction="y",vjust=4)











