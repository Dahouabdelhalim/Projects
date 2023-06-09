

#analysis code for survival and predator community data for manuscript:
#Myers AT, Bahlai CA, Landis DA. Habitat type influences Danaus plexippus (Lepidoptera: Nymphalidae) oviposition and egg 
#survival on Asclepias syriaca (Gentianales: Apocynaceae). Environmental Entomology.





# Load packages
library(bbmle)
library(arm)
library(lme4)
library(piecewiseSEM)
library(ggplot2)
library(reshape2)
library(plyr)
library(ggthemes)
library(multcomp)
library(emmeans)


#analysis involves doing the same for three different experimental periods-- Aug 2016, July 2017, and Aug 2018

###### First working with Aug 2016 data ############

#bring data in
data_aug2016<-read.csv(file="Aug_2016_survival_manuscript_csv.csv", header=TRUE)

#make a plant patch column called 'patch'that is defined by block, treatment, and exclosure type
data_aug2016<-within(data_aug2016, patch <- paste(block,treatment,exclosure_treatment, sep ='.'))
# then make a column called "surviving" which is the the total divided by initial number of eggs
data_aug2016<-within(data_aug2016, surviving <- (total/initial_count))

#drop all but 72 h
data_aug2016_3day <- data_aug2016[ which(data_aug2016$hours_since_deployment == 71), ]

#make one with only Open treatment for use later
data_aug_2016_3day_Open <- data_aug2016_3day[ which(data_aug2016_3day$exclosure_treatment == 'Open'),]

###first make null model, and then make model with exclosure treatment
m.0.aug2016 = glmer(cbind(total, initial_count-total) ~
                          1 + (1|block/patch), family=binomial(), data=data_aug2016_3day)

m.3.aug2016 = glmer(cbind(total, initial_count-total) ~
                      exclosure_treatment + (1|block/patch) + (1|treatment), family=binomial(), data=data_aug2016_3day)

###then do LRT to see if exclosure trt significant
anova(m.0.aug2016, m.3.aug2016)

###exclosure treatment significant, so next use emmeans to do pairwise comparisons among excl treatments
library(emmeans)
emmeans.aug2016.excl = emmeans(m.3.aug2016, ~exclosure_treatment)
contrast(emmeans.aug2016.excl, method = "pairwise", adjust = "holm")


###next we model survival as a function of treatment in Open only
###only need block as random effect, because only one patch per block
###also have to make another null model for this

m.00.aug2016= glmer(cbind(total, initial_count-total) ~
                      1 + (1|block), family=binomial(), data=data_aug_2016_3day_Open)
m.2.aug2016 = glmer(cbind(total, initial_count-total) ~
                      treatment + (1|block), family=binomial(), data=data_aug_2016_3day_Open)

###same as above, LRT to see if treatment significant
anova(m.00.aug2016, m.2.aug2016)
###treatment significant, so now can do pairwise comparisons using emmeans
emmeans.aug2016.trt = emmeans(m.2.aug2016, ~treatment, transform="response")
contrast(emmeans.aug2016.trt, alpha=0.05, method="pairwise", adjust="holm")






###plotting 

##first do the ddply 

#only op to 72 hours. then average across block and treatment and calculate SEM.
data_aug2016_72_exclosure <- data_aug2016[ which(data_aug2016$hours_since_deployment <  80 ), ]
data_aug2016_72_exclosure <- ddply(data_aug2016_72_exclosure, .(hours_since_deployment, exclosure_treatment), summarize,
                                   N=length(surviving),
                                   mean=mean(surviving),
                                   sd   = sd(surviving),
                                   se   = sd / sqrt(N) )

#only up to 72 hours. this one is for treatment comparisons
data_aug2016_72_open <- data_aug2016[ which(data_aug2016$hours_since_deployment <  80 
                                       & data_aug2016$exclosure_treatment == "Open"), ]

data_aug2016_72_open<- ddply(data_aug2016_72_open, .(hours_since_deployment, treatment), summarize,
         N=length(surviving),
         mean=mean(surviving),
         sd   = sd(surviving),
         se   = sd / sqrt(N) )

#then make ggplots

cols2 <- c("Open" = "springgreen1", "Closed" = "lightcoral", "Sham" = "royalblue1")

ggplot(data_aug2016_72_exclosure, aes(x=hours_since_deployment, shape=exclosure_treatment, y=mean, colour=exclosure_treatment, fill=exclosure_treatment))+
  #ggtitle("Aug 2016")+
  geom_point(size=2, show.legend = FALSE)+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, linetype=NA), show.legend = FALSE, alpha=0.3)+
  scale_alpha(guide = "none")+
  geom_line(aes(color=exclosure_treatment), show.legend = FALSE, size=1)+
  scale_colour_manual(values=cols2)+ 
  scale_fill_manual(values=cols2)+
  #xlab("\\nHours Since Deployment")+
  #ylab("Proportion Surviving\\n")+
  labs(x=NULL, y=NULL) +
  theme_few()+
  guides(alpha=FALSE)+
  theme(text = element_text(size=14))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 75), breaks=c(0, 10, 20, 30, 40, 50, 60, 70))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
ggsave('ggplotsurvivalbyexclosureAug2016_binomial.png', dpi = 1200, width=3, height=2.5)


cols <- c("Corn" = "gold2", "Prairie" = "lightgreen", "Soybean" = "lightslateblue", "Bare" = "indianred3")

Aug2016_72_figure <- ggplot(data_aug2016_72_open,
  aes(x=hours_since_deployment, y=mean, shape=treatment, colour=treatment, fill=treatment))+
  geom_point(size=2, show.legend = FALSE)+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, linetype=NA), show.legend = FALSE, alpha=0.3)+
  scale_alpha(guide = "none")+
  geom_line(aes(color=treatment), show.legend = FALSE, size=1)+
  scale_colour_manual(values=cols)+ 
  scale_fill_manual(values=cols)+
  #xlab("Hours Since Deployment")+
  #ylab("Proportion Surviving")+
  labs(x=NULL, y=NULL) +
  theme_few()+
  theme(text = element_text(size=14, colour = "black"))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 75), breaks=c(0, 10, 20, 30, 40, 50, 60, 70))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
ggsave('ggplotsurvivalbyhabitatAug2016_binomial.png', dpi = 1200, width=3, height=2.5)











#############july 2017##################
## exact same thing as for Aug 2016
#bring data in

data_july2017<-read.csv(file="July_2017_survival_manuscript_csv.csv", header=TRUE)

#make a plant patch column called 'patch'that is defined by block, treatment, and exclosure type
data_july2017<-within(data_july2017, patch <- paste(block,treatment,exclosure_treatment, sep ='.'))
# then make a column called "surviving" which is the the total divided by initial number of eggs
data_july2017<-within(data_july2017, surviving <- (total/initial_count))

#all but 72 hr obs 
data_july2017_3day <- data_july2017[ which(data_july2017$hours_since_deployment == 71), ]

#make one with only Open trt
data_july_2017_3day_Open <- data_july2017_3day[ which(data_july2017_3day$exclosure_treatment == 'Open'),]

###first make null model, and then make model with exclosure treatment
m.0.july2017 = glmer(cbind(total, initial_count-total) ~
                      1 + (1|block/patch), family=binomial(), data=data_july2017_3day)

m.3.july2017 = glmer(cbind(total, initial_count-total) ~
                      exclosure_treatment + (1|block/patch) + (1|treatment), family=binomial(), data=data_july2017_3day)

###then do LRT to see if exclosure trt significant
anova(m.0.july2017, m.3.july2017)

###exclosure treatment significant, so next use emmeans to do pairwise comparisons among excl treatments
library(emmeans)
emmeans.july2017.excl = emmeans(m.3.july2017, ~exclosure_treatment)
contrast(emmeans.july2017.excl, method = "pairwise", adjust = "holm")

###next we model survival as a function of treatment in Open only
###only need block as random effect, because only one patch per block
###also have to make another null model for this

m.00.july2017= glmer(cbind(total, initial_count-total) ~
                      1 + (1|block), family=binomial(), data=data_july_2017_3day_Open)
m.2.july2017 = glmer(cbind(total, initial_count-total) ~
                      treatment + (1|block), family=binomial(), data=data_july_2017_3day_Open)

###same as above, LRT to see if treatment significant
anova(m.00.july2017, m.2.july2017)
###treatment significant, so now can do pairwise comparisons using emmeans
emmeans.july2017.trt = emmeans(m.2.july2017, ~treatment, transform="response")
contrast(emmeans.july2017.trt, alpha=0.05, method="pairwise", adjust="holm")











###plotting 

##first do the ddply

#average across block and treatment and calculate SEM.
data_july2017_72_exclosure <- data_july2017[ which(data_july2017$hours_since_deployment <  80), ]
data_july2017_72_exclosure <- ddply(data_july2017_72_exclosure, .(hours_since_deployment, exclosure_treatment), summarize,
                                   N=length(surviving),
                                   mean=mean(surviving),
                                   sd   = sd(surviving),
                                   se   = sd / sqrt(N) )

#only Open, and only up to 72 hours. this one is for treatment comparisons
data_july2017_72_open <- data_july2017[ which(data_july2017$hours_since_deployment <  80 
                                            & data_july2017$exclosure_treatment == "Open"), ]

data_july2017_72_open<- ddply(data_july2017_72_open, .(hours_since_deployment, treatment), summarize,
                             N=length(surviving),
                             mean=mean(surviving),
                             sd   = sd(surviving),
                             se   = sd / sqrt(N) )

#then make ggplots

ggplot(data_july2017_72_exclosure, aes(x=hours_since_deployment, shape=exclosure_treatment, y=mean, colour=exclosure_treatment, fill=exclosure_treatment))+
  #ggtitle("July 2017")+
  geom_point(size=2, show.legend = FALSE)+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, linetype=NA), show.legend = FALSE, alpha=0.3)+
  scale_alpha(guide = "none")+
  geom_line(aes(color=exclosure_treatment), show.legend = FALSE, size=1)+
  scale_colour_manual(values=cols2)+ 
  scale_fill_manual(values=cols2)+
  #xlab("\\nHours Since Deployment")+
  #ylab("Proportion Surviving\\n")+
  labs(x=NULL, y=NULL) +
  theme_few()+
  guides(alpha=FALSE)+
  theme(text = element_text(size=14))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 75), breaks=c(0, 10, 20, 30, 40, 50, 60, 70))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
ggsave('ggplotsurvivalbyexclosurejuly2017_binomial.png', dpi = 1200, width=3, height=2.5)


cols <- c("Corn" = "gold2", "Prairie" = "lightgreen", "Soybean" = "lightslateblue", "Bare" = "indianred3")

july2017_72_figure <- ggplot(data_july2017_72_open,
                            aes(x=hours_since_deployment, y=mean, shape=treatment, colour=treatment, fill=treatment))+
  geom_point(size=2, show.legend = FALSE)+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, linetype=NA), show.legend = FALSE, alpha=0.3)+
  scale_alpha(guide = "none")+
  geom_line(aes(color=treatment), show.legend = FALSE, size=1)+
  scale_colour_manual(values=cols)+ 
  scale_fill_manual(values=cols)+
  #xlab("Hours Since Deployment")+
  #ylab("Proportion Surviving")+
  labs(x=NULL, y=NULL) +
  theme_few()+
  theme(text = element_text(size=14, colour = "black"))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 75), breaks=c(0, 10, 20, 30, 40, 50, 60, 70))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
  ggsave('ggplotsurvivalbyhabitatJuly2017_binomial.png', dpi = 1200, width=3, height=2.5)
  july2017_72_figure


  
  
  
  
  
  
  
  

#############AUG 2017##################
#bring data in

data_aug2017<-read.csv(file="Aug_2017_survival_manuscript_csv.csv", header=TRUE)

#make a plant patch column called 'patch'that is defined by block, treatment, and exclosure type
data_aug2017<-within(data_aug2017, patch <- paste(block,treatment,exclosure_treatment, sep ='.'))
# then make a column called "surviving" which is the the total divided by initial number of eggs
data_aug2017<-within(data_aug2017, surviving <- (total/initial_count))

#all but 72 hr obs 
data_aug2017_3day <- data_aug2017[ which(data_aug2017$hours_since_deployment == 72), ]

#make one with only Open trt
data_aug_2017_3day_Open <- data_aug2017_3day[ which(data_aug2017_3day$exclosure_treatment == 'Open'),]

###first make null model, and then make model with exclosure treatment
m.0.aug2017 = glmer(cbind(total, initial_count-total) ~
                      1 + (1|block/patch), family=binomial(), data=data_aug2017_3day)

m.3.aug2017 = glmer(cbind(total, initial_count-total) ~
                      exclosure_treatment + (1|block/patch) + (1|treatment), family=binomial(), data=data_aug2017_3day)

###then do LRT to see if exclosure trt significant
anova(m.0.aug2017, m.3.aug2017)

###exclosure treatment significant, so next use emmeans to do pairwise comparisons among excl treatments
library(emmeans)
emmeans.aug2017.excl = emmeans(m.3.aug2017, ~exclosure_treatment)
contrast(emmeans.aug2017.excl, method = "pairwise", adjust = "holm")  

###next we model survival as a function of treatment in Open only
###only need block as random effect, because only one patch per block
###also have to make another null model for this

m.00.aug2017= glmer(cbind(total, initial_count-total) ~
                      1 + (1|block), family=binomial(), data=data_aug_2017_3day_Open)
m.2.aug2017 = glmer(cbind(total, initial_count-total) ~
                      treatment + (1|block), family=binomial(), data=data_aug_2017_3day_Open)

###same as above, LRT to see if treatment significant
anova(m.00.aug2017, m.2.aug2017)
###treatment significant, so now can do pairwise comparisons using emmeans
emmeans.aug2017.trt = emmeans(m.2.aug2017, ~treatment, transform="response")
contrast(emmeans.aug2017.trt, alpha=0.05, method="pairwise", adjust="holm")





###plotting 

##first do the ddply

#average across block and treatment and calculate SEM.
data_aug2017_72_exclosure <- data_aug2017[ which(data_aug2017$hours_since_deployment <  80), ]
data_aug2017_72_exclosure <- ddply(data_aug2017_72_exclosure, .(hours_since_deployment, exclosure_treatment), summarize,
                                   N=length(surviving),
                                   mean=mean(surviving),
                                   sd   = sd(surviving),
                                   se   = sd / sqrt(N) )

#do only Open only up to 72 hours. this one is for treatment comparisons
data_aug2017_72_open <- data_aug2017[ which(data_aug2017$hours_since_deployment <  80
                                            & data_aug2017$exclosure_treatment == "Open"), ]

data_aug2017_72_open<- ddply(data_aug2017_72_open, .(hours_since_deployment, treatment), summarize,
                             N=length(surviving),
                             mean=mean(surviving),
                             sd   = sd(surviving),
                             se   = sd / sqrt(N) )

#then make ggplots

ggplot(data_aug2017_72_exclosure, aes(x=hours_since_deployment, shape=exclosure_treatment, y=mean, colour=exclosure_treatment, fill=exclosure_treatment))+
  #ggtitle("aug 2017")+
  geom_point(size=2, show.legend = FALSE)+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, linetype=NA), show.legend = FALSE, alpha=0.3)+
  scale_alpha(guide = "none")+
  geom_line(aes(color=exclosure_treatment), show.legend = FALSE, size=1)+
  scale_colour_manual(values=cols2)+ 
  scale_fill_manual(values=cols2)+
  #xlab("\\nHours Since Deployment")+
  #ylab("Proportion Surviving\\n")+
  labs(x=NULL, y=NULL) +
  theme_few()+
  guides(alpha=FALSE)+
  theme(text = element_text(size=14))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 75), breaks=c(0, 10, 20, 30, 40, 50, 60, 70))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
ggsave('ggplotsurvivalbyexclosureaug2017_binomial.png', dpi = 1200, width=3, height=2.5)



cols <- c("Corn" = "gold2", "Prairie" = "lightgreen", "Soybean" = "lightslateblue", "Bare" = "indianred3")

Aug2017_72_figure <- ggplot(data_aug2017_72_open,
                            aes(x=hours_since_deployment, y=mean, shape=treatment, colour=treatment, fill=treatment))+
  geom_point(size=2, show.legend = TRUE)+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, linetype=NA), show.legend = TRUE, alpha=0.3)+
  scale_alpha(guide = "none")+
  geom_line(aes(color=treatment), show.legend = TRUE, size=1)+
  scale_colour_manual(values=cols)+ 
  scale_fill_manual(values=cols)+
  #xlab("Hours Since Deployment")+
  #ylab("Proportion Surviving")+
  labs(x=NULL, y=NULL) +
  theme_few()+
  theme(text = element_text(size=14, colour = "black"))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 75), breaks=c(0, 10, 20, 30, 40, 50, 60, 70))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
  ggsave('ggplotsurvivalbyhabitatAug2017_binomial.png', dpi = 1200, width=3, height=2.5)
  Aug2017_72_figure




  
  
  
  
  
  
  ####PREDATOR GGPLOTTING############  
  
  #read in predator survey files for 3 periods (excluding July 206, because it was only Prairie and soy)
  predators_Aug2016<-read.csv(file= "predator_surveys_Aug_2016.csv", header=TRUE) 
  predators_2017<-read.csv(file="predator_surveys_july_&_august_2017.csv", header=TRUE) 
  
  
  # first, only looking at data for first 72 hours and divide 2017 data into july and august deployments
  
  predators_Aug2016 <- predators_Aug2016[ which(predators_Aug2016$hours_since_deployment < 73), ]
  predators_July2017<- predators_2017[ which(predators_2017$hours_since_deployment < 73 & predators_2017$deployment == 1),]
  predators_Aug2017 <- predators_2017[ which(predators_2017$hours_since_deployment < 73 & predators_2017$deployment == 2),] 
  
  
  # for the anlaysis I used sums of predators, because glm.nb requires integers
  # but for my figure, it is better to have mean numbers of preds per stem....
  # so first here are separate data files for preds per stem
  predators_July2017_perstem=predators_July2017
  predators_July2017_perstem[, 7:28] <- predators_July2017[, 7:28]/3
  
  predators_July2017_perstem=predators_July2017
  predators_July2017_perstem[, 8:31] <- predators_July2017[, 8:31]/3
  
  predators_Aug2017_perstem=predators_Aug2017
  predators_Aug2017_perstem[, 8:31] <- predators_Aug2017[, 8:31]/3
  
  ###have to use melt function to make a column that shows predator id as a column
  predators_July2017_melt <- melt(predators_July2017_perstem, measure.vars = c("total_ants", "total_no_ants"), variable.name = "Predator_ID")
  predators_Aug2017_melt <- melt(predators_Aug2017_perstem, measure.vars = c("total_ants", "total_no_ants"), variable.name = "Predator_ID")
  
  predators_July2017_total_melt <- melt(predators_July2017_perstem, measure.vars = c("total"), variable.name = "Predator_ID")
  predators_Aug2017_total_melt <- melt(predators_Aug2017_perstem, measure.vars = c("total"), variable.name = "Predator_ID")
  
  # next using ddply to makemeans averaging across blocks for each treatment and exclosure
  
  predators_2016_means <- ddply(predators_Aug2016, .(treatment, exclosure_treatment), summarize,
                                N=4,
                                mean=mean(total_no_ants),
                                sd   = sd(total_no_ants),
                                se   = sd / sqrt(N) )
  
  predators_July2017_means <- ddply(predators_July2017_melt, .(treatment, exclosure_treatment, Predator_ID), summarize,
                                    N=4,
                                    mean=mean(value),
                                    sd   = sd(value),
                                    se   = sd / sqrt(N) )
  
  predators_Aug2017_means <- ddply(predators_Aug2017_melt, .(treatment, exclosure_treatment, Predator_ID), summarize,
                                   N=4,
                                   mean=mean(value),
                                   sd   = sd(value),
                                   se   = sd / sqrt(N) )
  
  #for total

  
  predators_July2017_means_total  <- ddply(predators_July2017_total_melt, .(treatment, exclosure_treatment, Predator_ID), summarize,
                                    N=4,
                                    mean=mean(value),
                                    sd   = sd(value),
                                    se   = sd / sqrt(N) )
  
  predators_Aug2017_means_total  <- ddply(predators_Aug2017_total_melt, .(treatment, exclosure_treatment, Predator_ID), summarize,
                                   N=4,
                                   mean=mean(value),
                                   sd   = sd(value),
                                   se   = sd / sqrt(N) )
  
  ggplot(data=predators_2016_means[which(predators_2016_means$exclosure_treatment=="open"),], aes(x=treatment, y=mean)) + 
    geom_bar(position="stack", stat="identity", colour = "black", width = 0.8, fill = "lightgray", show.legend = T) +
    geom_errorbar(aes(ymin=mean-0, ymax=mean+se), colour="black", width=.2, position=position_dodge(.9)) +
    #theme(panel.background = element_blank(), axis.text.x = element_blank(),  axis.ticks = element_blank())+
    #ggtitle("August 2017")+
    xlab("")+
    ylab("Mean Predators/stem/survey")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5))+
    theme_few()+
    theme(axis.text.x  = element_text(angle=30, vjust=0.5, size = 11))
  ggsave('totalpredatorsAug2016.png', dpi = 1200, width = 2.7, height = 3)
  
  
  ggplot(data= NULL, 
         aes(x=treatment, y=mean, fill=Predator_ID)) + 
    geom_bar(data=predators_July2017_means[which(predators_July2017_means$exclosure_treatment=="open"),],
             position="stack", stat="identity", colour = "black", show.legend = F, width = 0.8,)+ 
    scale_fill_manual(values=c("white", "#4c4c4c", "lightgray"), 
                      name = "",
                      breaks=c("total_ants", "total_no_ants", ""),
                      labels=c("ants", "non-ant predators", "")) +
    #geom_errorbar(aes(ymin=mean-0, ymax=mean+se), colour="black", width=.2, position=position_dodge(.9)) +
    #theme(panel.background = element_blank(), axis.text.x = element_blank(),  axis.ticks = element_blank())+
    #ggtitle("July 2017")+
    xlab("")+
    ylab("Mean Predators/stem/survey")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5))+
    theme_few()+
    theme(axis.text.x  = element_text(angle=30, vjust=0.5, size = 11))+
    geom_bar(data=predators_July2017_means_total[which(predators_July2017_means_total$exclosure_treatment=="open"),],
             stat="identity", colour = "black", show.legend = F, alpha = 0, width = 0.8,)+
    geom_errorbar(data=predators_July2017_means_total[which(predators_July2017_means_total$exclosure_treatment=="open"),],
                  aes(ymin=mean-0, ymax=mean+se), colour="black", width=.2, position=position_dodge(.9)) 
  
  ggsave('totalpredatorsJuly2017.png', dpi = 1200, width = 2.7, height = 3)
  
  
  
  ggplot(data= NULL, 
         aes(x=treatment, y=mean, fill=Predator_ID)) + 
    geom_bar(data=predators_Aug2017_means[which(predators_Aug2017_means$exclosure_treatment=="open"),],
             position="stack", stat="identity", colour = "black", show.legend = F, width = 0.8,)+ 
    scale_fill_manual(values=c("white", "#4c4c4c", "lightgray"), 
                      name = "",
                      breaks=c("total_ants", "total_no_ants", ""),
                      labels=c("ants", "non-ant predators", "")) +
    #geom_errorbar(aes(ymin=mean-0, ymax=mean+se), colour="black", width=.2, position=position_dodge(.9)) +
    #theme(panel.background = element_blank(), axis.text.x = element_blank(),  axis.ticks = element_blank())+
    #ggtitle("Aug 2017")+
    xlab("")+
    ylab("Mean Predators/stem/survey")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5))+
    theme_few()+
    theme(axis.text.x  = element_text(angle=30, vjust=0.5, size = 11))+
    geom_bar(data=predators_Aug2017_means_total[which(predators_Aug2017_means_total$exclosure_treatment=="open"),],
             stat="identity", colour = "black", show.legend = F, alpha = 0, width = 0.8,)+
    geom_errorbar(data=predators_Aug2017_means_total[which(predators_Aug2017_means_total$exclosure_treatment=="open"),],
                  aes(ymin=mean-0, ymax=mean+se), colour="black", width=.2, position=position_dodge(.9)) 
  
  ggsave('totalpredatorsAug2017.png', dpi = 1200, width = 2.7, height = 3)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


  #####PREDATOR ANALYSIS#######
  
  #read in predator survey files for 3 periods (excluding July 206, because it was only Prairie and soy)
  predators_Aug2016<-read.csv(file="predator_surveys_Aug_2016.csv", header=TRUE) 
  predators_2017<-read.csv(file="predator_surveys_july_&_august_2017.csv", header=TRUE) 
  
  
  # first, only looking at data for first 72 hours and divide 2017 data into july and august deployments
  
  predators_Aug2016 <- predators_Aug2016[ which(predators_Aug2016$hours_since_deployment < 73), ]
  predators_July2017<- predators_2017[ which(predators_2017$hours_since_deployment < 73 & predators_2017$deployment == 1),]
  predators_Aug2017 <- predators_2017[ which(predators_2017$hours_since_deployment < 73 & predators_2017$deployment == 2),] 
  
  ### want to compare the effect of treatment on total predator numbers
  ### doing open only
  ### doing separate analysis for each of the three experiments
  ### have to do negative binomial rather than simple anova, because data are not normally distributed
  ###data not normally distributed
  ###hist(predators_July2017_means_by_date[which(predators_July2017_means_by_date$exclosure_treatment == "open"),5])
  ###hist(predators_July2017_means_by_date[which(predators_July2017_means_by_date$exclosure_treatment == "open"),5])
  ###hist(predators_Aug2017_means_by_date[which(predators_Aug2017_means_by_date$exclosure_treatment == "open"),5])
  
  ####first have to sum across dates
  
  predators_Aug2016_sum_by_date <- ddply(predators_Aug2016, .(treatment, exclosure_treatment, block), summarize,
                                          sum.preds=sum(total_no_ants),
                                          observations=length(total_no_ants))
  
  predators_July2017_sum_by_date <- ddply(predators_July2017, .(treatment, exclosure_treatment, block), summarize,
                                          sum.preds=sum(total),
                                          sum.ants=sum(total_ants),
                                          observations=length(total))
  
  predators_Aug2017_sum_by_date <- ddply(predators_Aug2017, .(treatment, exclosure_treatment, block), summarize,
                                         sum.preds=sum(total),
                                         sum.ants=sum(total_ants),
                                         observations=length(total))
  
  
  #nb.GLM
  
  m.preds.Aug2016 <- glm.nb(sum.preds ~ treatment, data = predators_Aug2016_sum_by_date[which(predators_Aug2016_sum_by_date$exclosure_treatment == "open"),])
  m.preds.July2017 <- glm.nb(sum.preds  ~ treatment, data = predators_July2017_sum_by_date[which(predators_July2017_sum_by_date$exclosure_treatment == "open"),])
  m.preds.Aug2017 <- glm.nb(sum.preds  ~ treatment, data = predators_Aug2017_sum_by_date[which(predators_Aug2017_sum_by_date$exclosure_treatment == "open"),])
  
  null.preds.Aug2016 <- glm.nb(sum.preds ~ 1, data = predators_Aug2016_sum_by_date[which(predators_Aug2016_sum_by_date$exclosure_treatment == "open"),])
  null.preds.July2017 <- glm.nb(sum.preds  ~ 1, data = predators_July2017_sum_by_date[which(predators_July2017_sum_by_date$exclosure_treatment == "open"),])
  null.preds.Aug2017 <- glm.nb(sum.preds  ~ 1, data = predators_Aug2017_sum_by_date[which(predators_Aug2017_sum_by_date$exclosure_treatment == "open"),])
  
  ##LRT comparing models
  anova(m.preds.Aug2016, null.preds.Aug2016, test="Chisq")
  anova(m.preds.July2017, null.preds.July2017, test="Chisq")
  anova(m.preds.Aug2017, null.preds.Aug2017, test="Chisq")
  
  #pairwise comparisons for just Aug2017, which was the only significant one
  emmeans.preds.Aug2017<-emmeans(m.preds.Aug2017, ~treatment, type = "response")
  contrast.preds.Aug2017<-contrast(emmeans.preds.Aug2017, alpha=0.05, method="pairwise", adjust="holm")
  
  
  c##HERE DOING THE EXACT SAME THING AGAIN BUT WITH ANTS ONLY DATA AND FOR JUST 2017 WHEN ANTS WERE COUNTED
  m.ants.July2017 <- glm.nb(sum.ants  ~ treatment, data = predators_July2017_sum_by_date[which(predators_July2017_sum_by_date$exclosure_treatment == "open"),])
  m.ants.Aug2017 <- glm.nb(sum.ants  ~ treatment, data = predators_Aug2017_sum_by_date[which(predators_Aug2017_sum_by_date$exclosure_treatment == "open"),])
  
  null.ants.July2017 <- glm.nb(sum.ants  ~ 1, data = predators_July2017_sum_by_date[which(predators_July2017_sum_by_date$exclosure_treatment == "open"),])
  null.ants.Aug2017 <- glm.nb(sum.ants  ~ 1, data = predators_Aug2017_sum_by_date[which(predators_Aug2017_sum_by_date$exclosure_treatment == "open"),])
  
  ##LRT comparing models
  anova(m.ants.July2017, null.ants.July2017, test="Chisq")
  anova(m.ants.Aug2017, null.ants.Aug2017, test="Chisq")
  
  #pairwise comparisons for just Aug2017, which was the only significant one
  emmeans.ants.Aug2017<-emmeans(m.ants.Aug2017, ~treatment, type = "response")
  contrast.ants.Aug2017<-contrast(emmeans.ants.Aug2017, alpha=0.05, method="pairwise", adjust="holm")
  

  
  
  
  
  
  ### Creating a figure combining oviposition and survival
  
  #first generate the proper oviposition data table (same as from ovip code file)
  #read in 2016 data
  oviposition2016<-read.csv(file="oviposition2016.csv", header=TRUE) #read in oviposition2016 file
  oviposition2016<-na.omit(oviposition2016) #get rid of na's. There were several incidents when we were unable to count eggs (broken plants, plants were covered by exclosures, etc)
  #include only August data
  oviposition_aug2016<-oviposition2016[ which(oviposition2016$deployment == '3'), ]
  #####for days with more than one egg check, doug wants me to add up all the eggs and divide by the 
  ######number of plants present in the plot on that day (or average if it changed).
  #install and use the ddply function to find the sum of the number of eggs/patch/check and the number of plants present in each check
  library(plyr)
  oviposition_aug2016.avg <-ddply(oviposition_aug2016, .(treatment, date, time, block, deployment), summarize, 
                                  monarch_eggs.sum=sum(monarch_eggs),
                                  nplants=length(monarch_eggs))
  #average plants checked per day and sum all the eggs found per day
  oviposition_aug2016.avg.2 <-ddply(oviposition_aug2016.avg, .(treatment, date, block, deployment), summarize, 
                                    nplants.mean=mean(nplants),
                                    monarch_eggs.sum=sum(monarch_eggs.sum))
  #divide number of eggs seen in a day by average number of plants present that day
  oviposition_aug2016.avg.2 <-ddply(oviposition_aug2016.avg.2, .(treatment, date, block, deployment, monarch_eggs.sum, nplants.mean), summarize,
                                    monarch_eggs.per.plant=monarch_eggs.sum/nplants.mean)
  ##average across all dates, treating each date like a subsample, calculate the number of plant checks, sum all eggs found#####
  oviposition_aug2016.avg.3<-ddply(oviposition_aug2016.avg.2, .(treatment, block), summarize,
                                   monarch_eggs.mean=mean(monarch_eggs.per.plant),
                                   nplants.checks = (sum(nplants.mean)),
                                   monarch_eggs.sum = sum(monarch_eggs.sum),
                                   days.checked = length(treatment))
  
  #read in 2017 oviposition data
  oviposition2017<-read.csv(file="oviposition2017.csv", header=TRUE) #read in oviposition2017 file
  oviposition2017<-na.omit(oviposition2017) #get rid of na's. There were several incidents when we were unable to count eggs (broken plants, plants were covered by exclosures, etc)
  #include only August data
  oviposition_aug2017<-oviposition2017[ which(oviposition2017$deployment == '3'), ]
  #####for days with more than one egg check, doug wants me to add up all the eggs and divide by the 
  ######number of plants present in the plot on that day (or average if it changed).
  #install and use the ddply function to find the sum of the number of eggs/patch/check and the number of plants present in each check
  library(plyr)
  oviposition_aug2017.avg <-ddply(oviposition_aug2017, .(treatment, date, time, block, deployment), summarize, 
                                  monarch_eggs.sum=sum(monarch_eggs),
                                  nplants=length(monarch_eggs))
  #average plants checked per day and sum all the eggs found per day
  oviposition_aug2017.avg.2 <-ddply(oviposition_aug2017.avg, .(treatment, date, block, deployment), summarize, 
                                    nplants.mean=mean(nplants),
                                    monarch_eggs.sum=sum(monarch_eggs.sum))
  #divide number of eggs seen in a day by average number of plants present that day
  oviposition_aug2017.avg.2 <-ddply(oviposition_aug2017.avg.2, .(treatment, date, block, deployment, monarch_eggs.sum, nplants.mean), summarize,
                                    monarch_eggs.per.plant=monarch_eggs.sum/nplants.mean)
  ##average across all dates, treating each date like a subsample, calculate the number of plant checks, sum all eggs found#####
  oviposition_aug2017.avg.3<-ddply(oviposition_aug2017.avg.2, .(treatment, block), summarize,
                                   monarch_eggs.mean=mean(monarch_eggs.per.plant),
                                   nplants.checks = (sum(nplants.mean)),
                                   monarch_eggs.sum = sum(monarch_eggs.sum),
                                   days.checked = length(treatment))
  
  
  
  #merge ovipostion and survival by block and treatment horizontally
  ovip_survival_aug2016 <- merge(oviposition_aug2016.avg.3,data_aug_2016_3day_Open, by=c("block","treatment"))
  ovip_survival_aug2017 <- merge(oviposition_aug2017.avg.3,data_aug_2017_3day_Open, by=c("block","treatment"))
  
  #make a year column
  ovip_survival_aug2016<-within(ovip_survival_aug2016, Year <- paste("August 2016"))
  ovip_survival_aug2017<-within(ovip_survival_aug2017, Year <- paste("August 2017"))
  
  #vertically merge
  ovip_survival <- rbind(ovip_survival_aug2016, ovip_survival_aug2017)
  colnames(ovip_survival)[which(names(ovip_survival) == "treatment")] <- "Treatment"
  
  #plot
  cols <- c("Corn" = "gold2", "Prairie" = "lightgreen", "Soybean" = "lightslateblue", "Bare" = "indianred3")
  ggplot(ovip_survival, aes(y=monarch_eggs.mean,x=surviving))+
    geom_point(size=4,  aes(shape=Year, colour=Treatment))+
    geom_smooth(aes(group=Treatment, colour=Treatment), method ="lm", se=FALSE, fullrange = T)+  
    scale_colour_manual(values=cols)+
    theme_few()+
    xlim(0,1)+
    ylim(0,0.5)+
    ylab("Mean eggs/stem/day")+
    xlab("Proportion surviving\\n")
  ggsave('ovip_survival_figure.png', width = 6, height = 4)
  
  
  
  
  ### removing corn outlier
  
  ovip_survival_no_outlier <- ovip_survival[ which(ovip_survival$monarch_eggs.mean < 0.4), ]
  
  cols <- c("Corn" = "gold2", "Prairie" = "lightgreen", "Soybean" = "lightslateblue", "Bare" = "indianred3")
  ggplot(ovip_survival_no_outlier, aes(y=monarch_eggs.mean,x=surviving))+
    geom_point(size=4,  aes(shape=Year, colour=Treatment))+
    geom_smooth(aes(group=Treatment, colour=Treatment), method ="lm", se=FALSE, fullrange = T)+  
    scale_colour_manual(values=cols)+
    theme_few()+
    xlim(0,1)+
    ylim(0,0.3)+
    ylab("Mean eggs/stem/day")+
    xlab("Proportion surviving\\n")
  ggsave('ovip_survival_figure_no_outlier.png', width = 6, height = 4, dpi=1200)
  
  
