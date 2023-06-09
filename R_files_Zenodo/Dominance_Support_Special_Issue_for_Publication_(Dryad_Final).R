##############################################################
#                                                            #
#     Title: Developmental Sex Differences in Hyenas         #
#     Author: S. Kevin McCormick                             #
#     Date: 12/05/2022                                        #
#     Purpose: Comparisons of Winning Agnostic Encounters    #
#                Against the Opposite Sex                    #
#                                                            #
##############################################################

#set working directory
setwd("C:/Users/skmcc/Documents/R/R_wd")

#specify packages that you'll need in this script
library(dplyr)
library(ggplot2)
library(ggbreak)
library(GGally)
library(bbmle)
library(lme4)
library(ggeffects)
library(sjmisc)
library(sjlabelled)
library(sjPlot)
library(sjstats)
library(DHARMa)
library(emmeans)
library(scales)
citation(DHARMa)
citation("DHARMa")

#Import data
Data_File1 <- read.csv("tblAggression_AllClan_Age_Sex_Raw_Dryad.csv")

#Lets modify the data set to what we need
str(Data_File1)

#Remove counter attack aggressions
Data_File1 <- Data_File1[Data_File1$Context != "counterattack",]

#For now we will also remove mating a natal den location contexts
Data_File1 <- Data_File1[Data_File1$location != "n",]
Data_File1 <- Data_File1[Data_File1$location != "g",]

#There are females that where adults already when the study began and were never aged
Data_File1$Actor_Age_Class[Data_File1$Actor_Age_Class == "Unknown" &
                             Data_File1$Actor_Sex == "Female"] <- "Adult"
Data_File1$Recipient_Age_Class[Data_File1$Recipient_Age_Class == "Unknown" &
                                 Data_File1$Recipient_Sex == "Female"] <- "Adult"

######Define a winning aggression based on the first response

#First remove unknown responses (i.e. the entry is termed "unknown")
Data_File1 <- Data_File1[Data_File1$response1 != "unknown",]

#Then remove blank responses
Data_File1 <- Data_File1[Data_File1$response1 != "",]

#Define the wins based off of are submissive response
Data_File1$Win <- (Data_File1$response1 == "bo" |
                     Data_File1$response1 == "cc" |
                     Data_File1$response1 == "dp" |
                     Data_File1$response1 == "eb" |
                     Data_File1$response1 == "grin" |
                     Data_File1$response1 == "hb" |
                     Data_File1$response1 == "run" |
                     Data_File1$response1 == "sp" |
                     Data_File1$response1 == "squeal"|
                     Data_File1$response1 == "s1" |
                     Data_File1$response1 == "s2" |
                     Data_File1$response1 == "s3") == TRUE
Data_File1$Win[Data_File1$Win == TRUE] <- 1
Data_File1$Win[Data_File1$Win == FALSE] <- 0



#Remove strange sequence values/typos
Data_File1 <- Data_File1[Data_File1$Sequence < 2,]

#Remove strange negative sequence values/typos
Data_File1 <- Data_File1[Data_File1$Sequence > -1,]

#Pair up all the sex vs sex interactions
Data_File_Female_Male <- Data_File1[(Data_File1$Actor_Sex =="Female") & (Data_File1$Recipient_Sex =="Male"),]
Data_File_Male_Female <- Data_File1[(Data_File1$Actor_Sex =="Male") & (Data_File1$Recipient_Sex =="Female"),]
Data_File_Female_Female <- Data_File1[(Data_File1$Actor_Sex =="Female") & (Data_File1$Recipient_Sex =="Female"),]
Data_File_Male_Male <- Data_File1[(Data_File1$Actor_Sex =="Male") & (Data_File1$Recipient_Sex =="Male"),]

#Bind everything back together
Data_File <- rbind(Data_File_Female_Male,
                           Data_File_Male_Female,
                           Data_File_Female_Female,
                           Data_File_Male_Male)

#Lets calculate the number of hyenas present (both known IDs and unIDs)
known_hyenas <- Data_File$hyenas
coal.size <- function(known_hyenas_fuction){
  known_hyenas_fuction <- strsplit(known_hyenas_fuction, ',')
  return(unlist(lapply(known_hyenas_fuction, function(x)(return(sum(unique(x) != ''))))))
}
coal.size(known_hyenas)
Data_File$known_groupsize <- coal.size(Data_File$hyenas)
unknown_unidhyenas <- Data_File$unidhyenas
coal.size <- function(unknown_unidhyenas_fuction){
  unknown_unidhyenas_fuction <- strsplit(unknown_unidhyenas_fuction, ',')
  return(unlist(lapply(unknown_unidhyenas_fuction, function(x)(return(sum(unique(x) != ''))))))
}
coal.size(unknown_unidhyenas)
Data_File$unknown_groupsize <- coal.size(Data_File$unidhyenas)
Data_File <- Data_File[(Data_File$known_groupsize + Data_File$unknown_groupsize) > 2,]

#This function also works on the groupcomp column for the number of supporters
Support_Number <- Data_File$groupcomp
coal.size <- function(Support_Number_fuction){
  Support_Number_fuction <- strsplit(Support_Number_fuction, ',')
  return(unlist(lapply(Support_Number_fuction, function(x)(return(sum(unique(x) != ''))))))
}
coal.size(Support_Number)
Data_File$Support_Number <- coal.size(Data_File$groupcomp)

#Lets prep the data for analyses

Data_File$Supported <- Data_File$Sequence
Data_File$Supported[Data_File$Supported == 1] <- "Supported"
Data_File$Supported[Data_File$Supported == 0] <- "Unsupported"


#Create separate files for each ACTOR age class
Data_File_Cub <- Data_File[Data_File$Actor_Age_Class =="Cub",]
Data_File_Sub <- Data_File[Data_File$Actor_Age_Class =="Sub" &
                             (Data_File$Recipient_Age_Class =="Adult" |
                                Data_File$Recipient_Age_Class =="Sub" |
                                Data_File$Recipient_Age_Class =="Unknown"),]
Data_File_Adult <- Data_File[(Data_File$Actor_Age_Class =="Adult" | Data_File$Actor_Age_Class =="Unknown") &
                               (Data_File$Recipient_Age_Class =="Adult" | Data_File$Recipient_Age_Class =="Unknown"),]
Data_File_Adult$Actor_Sex[Data_File_Adult$Actor_Age_Class == "Unknown" &
                            Data_File_Adult$Actor_Sex == "Male"] <- "Natal_Male"
Data_File_Adult$Recipient_Sex[Data_File_Adult$Recipient_Age_Class == "Unknown" &
                                Data_File_Adult$Recipient_Sex == "Male"] <- "Natal_Male"
Data_File_Adult$Recipient_Sex[Data_File_Adult$Recipient_Sex == "Male"] <- "Immigrant_Male"
Data_File_Adult$Actor_Sex[Data_File_Adult$Actor_Sex == "Male"] <- "Immigrant_Male"

Data_File_Cub$Actor_Sex <- as.factor(Data_File_Cub$Actor_Sex)
Data_File_Cub$Recipient_Sex <- as.factor(Data_File_Cub$Recipient_Sex)
Data_File_Cub$Actor_Age_Class <- as.factor(Data_File_Cub$Actor_Age_Class)
Data_File_Sub$Actor_Sex <- as.factor(Data_File_Sub$Actor_Sex)
Data_File_Sub$Recipient_Sex <- as.factor(Data_File_Sub$Recipient_Sex)
Data_File_Sub$Actor_Age_Class <- as.factor(Data_File_Sub$Actor_Age_Class)
Data_File_Adult$Actor_Sex <- as.factor(Data_File_Adult$Actor_Sex)
Data_File_Adult$Recipient_Sex <- as.factor(Data_File_Adult$Recipient_Sex)
Data_File_Adult$Actor_Age_Class <- as.factor(Data_File_Adult$Actor_Age_Class)

write.csv(Data_File_Adult,"C:/Users/skmcc/Documents/R/R_wd/Data_File_Adult.csv", row.names = FALSE)
write.csv(Data_File_Cub,"C:/Users/skmcc/Documents/R/R_wd/Data_File_Cub.csv", row.names = FALSE)
write.csv(Data_File_Sub,"C:/Users/skmcc/Documents/R/R_wd/Data_File_Sub.csv", row.names = FALSE)



#######Basic Sex Difference in Winning Models


Model_AllClans_Mating_Cub <- glmer(Win ~ Actor_Sex*Recipient_Sex
                            + (1|ID:Recipient) + (1|Session),
                            family = binomial(link = "logit"),
                            data = Data_File_Cub,
                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

Model_AllClans_Mating_Sub <- glmer(Win ~ Actor_Sex*Recipient_Sex
                             + (1|Session),
                            family = binomial(link = "logit"),
                            data = Data_File_Sub,
                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Model overfit with pairs, so pairs removed

Model_AllClans_Mating_Adult <- glmer(Win ~ Actor_Sex*Recipient_Sex
                              + (1|ID:Recipient) + (1|Session),
                              family = binomial(link = "logit"),
                              data = Data_File_Adult,
                              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))



#Summarize the results of the best models
summary(Model_AllClans_Mating_Cub)
Model_AllClans_Mating_Cub_DF <- summary(Model_AllClans_Mating_Cub)
write.csv(Model_AllClans_Mating_Cub_DF$coefficients, "C:/Users/skmcc/Documents/R/R_wd/Model_AllClans_Mating_Cub_DF.csv")
summary(Model_AllClans_Mating_Sub)
Model_AllClans_Mating_Sub_DF <- summary(Model_AllClans_Mating_Sub)
write.csv(Model_AllClans_Mating_Sub_DF$coefficients, "C:/Users/skmcc/Documents/R/R_wd/Model_AllClans_Mating_Sub_DF.csv")
summary(Model_AllClans_Mating_Adult)
Model_AllClans_Mating_Adult_DF <- summary(Model_AllClans_Mating_Adult)
write.csv(Model_AllClans_Mating_Adult_DF$coefficients, "C:/Users/skmcc/Documents/R/R_wd/Model_AllClans_Mating_Adult_DF.csv")


#Look at the pairwise comparisons of the full model
emmeans(Model_AllClans_Mating_Cub, specs = pairwise ~ Recipient_Sex|Actor_Sex, type="response", adjust="tukey")
emmeans(Model_AllClans_Mating_Sub, specs = pairwise ~ Recipient_Sex|Actor_Sex, type="response", adjust="tukey")
emmeans(Model_AllClans_Mating_Adult, specs = pairwise ~ Recipient_Sex|Actor_Sex, type="response", adjust="tukey")



Cub <- ggeffects::ggemmeans(Model_AllClans_Mating_Cub, c("Recipient_Sex", "Actor_Sex"))
Cub_Plot <- plot(Cub)+
  ggtitle("B") +
  xlab("Recipient's Sex") +
  ylab("Probability of Recieving a Submissive Response")+
  theme_classic()
Plot_B <- Cub_Plot +labs(colour="Actor's Sex")
Sub <- ggeffects::ggemmeans(Model_AllClans_Mating_Sub, c("Recipient_Sex", "Actor_Sex"))
Sub_Plot <- plot(Sub)+
  ggtitle("C") +
  xlab("Recipient's Sex") +
  ylab("Probability of Recieving a Submissive Response")+
  theme_classic()
Plot_C <- Sub_Plot +labs(colour="Actor's Sex")
Adult <- ggeffects::ggemmeans(Model_AllClans_Mating_Adult, c("Recipient_Sex","Actor_Sex"))
Adult_Plot <- plot(Adult)+
  ggtitle("A") +
  xlab("Recipent's Sex") +
  ylab("Probability of Recieving a Submissive Response")+
  theme_classic()+
  #scale_y_continuous(labels = percent,
                     #limits = c(.4,1))+
  scale_y_break(c(.675,.875))

library(ggbreak)
Plot_A <- Adult_Plot +labs(colour="Actor's Sex")


library(ggpubr)
figure <- ggarrange(
  Plot_A + rremove("ylab") + rremove("xlab") +
    scale_y_break(c(.675,.875)),
  # Second row with box and dot plots
  ggarrange(Plot_B + rremove("ylab") + rremove("xlab"),
            Plot_C+ rremove("ylab") + rremove("xlab"),
            ncol = 2, labels = c("", ""),
            common.legend = TRUE, legend = "right"),
  nrow = 2,
  labels = ""       # Label of the line plot
)

library(grid)
annotate_figure(figure, left = textGrob("Probability of Recieving a Submissive Response", rot = 90, vjust = 0.75, gp = gpar(cex = 1.125)),
                bottom = textGrob("Recipient's Sex", hjust = 0.75, gp = gpar(cex = 1.125)))




#Now We Look at Supported Interactions to see if this is why females win more


#Models
Model_AllClans_Mating_Supported_Cub <- glmer(Win ~ Actor_Sex*Recipient_Sex*Supported
                                      + (1|ID:Recipient) + (1|Session),
                                      family = binomial(link = "logit"),
                                      data = Data_File_Cub,
                                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

Model_AllClans_Mating_Supported_Sub <- glmer(Win ~ Actor_Sex*Recipient_Sex*Supported
                                       + (1|Session),
                                      family = binomial(link = "logit"),
                                      data = Data_File_Sub,
                                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


### Overfitting, so random effect of ID:Recipient was removed

Model_AllClans_Mating_Supported_Adult <- glmer(Win ~ Actor_Sex*Recipient_Sex*Supported
                                        + (1|ID:Recipient) + (1|Session),
                                        family = binomial(link = "logit"),
                                        data = Data_File_Adult,
                                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))



#Summarize the results of the best models
summary(Model_AllClans_Mating_Supported_Cub)
Model_AllClans_Mating_Supported_Cub_DF <- summary(Model_AllClans_Mating_Supported_Cub)
write.csv(Model_AllClans_Mating_Supported_Cub_DF$coefficients, "C:/Users/skmcc/Documents/R/R_wd/Model_AllClans_Mating_Supported_Cub_DF.csv")
summary(Model_AllClans_Mating_Supported_Sub)
Model_AllClans_Mating_Supported_Sub_DF <- summary(Model_AllClans_Mating_Supported_Sub)
write.csv(Model_AllClans_Mating_Supported_Sub_DF$coefficients, "C:/Users/skmcc/Documents/R/R_wd/Model_AllClans_Mating_Supported_Sub_DF.csv")
summary(Model_AllClans_Mating_Supported_Adult)
Model_AllClans_Mating_Supported_Adult_DF <- summary(Model_AllClans_Mating_Supported_Adult)
write.csv(Model_AllClans_Mating_Supported_Adult_DF$coefficients, "C:/Users/skmcc/Documents/R/R_wd/Model_AllClans_Mating_Supported_Adult_DF.csv")


#Pairwise comparisons from full model
emmeans(Model_AllClans_Mating_Supported_Cub, specs = pairwise ~ Supported|Actor_Sex|Recipient_Sex, type="response", adjust="tukey")
emmeans(Model_AllClans_Mating_Supported_Sub, specs = pairwise ~ Supported|Actor_Sex|Recipient_Sex, type="response", adjust="tukey")
emmeans(Model_AllClans_Mating_Supported_Adult, specs = pairwise ~ Supported|Actor_Sex|Recipient_Sex, type="response", adjust="tukey")


##Plot the full model
Cub_Supported <- ggeffects::ggemmeans(Model_AllClans_Mating_Supported_Cub, c("Recipient_Sex", "Supported", "Actor_Sex"))
Cub_Plot_Supported <- plot(Cub_Supported)+
  ggtitle("B") +
  xlab("Recipient's Sex") +
  ylab("Probability of Recieving a Submissive Response")+
  theme_classic()
Plot_E <- Cub_Plot_Supported +labs(colour="Support")
Sub_Supported <- ggeffects::ggemmeans(Model_AllClans_Mating_Supported_Sub, c("Recipient_Sex", "Supported", "Actor_Sex"))
Sub_Plot_Supported <- plot(Sub_Supported)+
  ggtitle("C") +
  xlab("Recipient's Sex") +
  ylab("Probability of Recieving a Submissive Response")+
  theme_classic()
Plot_F <- Sub_Plot_Supported +labs(colour="Support")
Adult_Supported <- ggeffects::ggemmeans(Model_AllClans_Mating_Supported_Adult, c("Recipient_Sex", "Supported", "Actor_Sex"))
Adult_Plot_Supported <- plot(Adult_Supported)+
  ggtitle("A") +
  xlab("Recipient's Sex") +
  ylab("Probability of Recieving a Submissive Response")+
  theme_classic()
Plot_D <- Adult_Plot_Supported + theme(legend.position="none")

library(ggpubr)
figure.support <- ggarrange(
  Plot_D + rremove("ylab") + rremove("xlab") +
    scale_y_break(c(.675,.875)),
  # Second row with box and dot plots
  ggarrange(Plot_E + rremove("ylab") + rremove("xlab"),
            Plot_F+ rremove("ylab") + rremove("xlab"),
            ncol = 2, labels = c("", ""),
            common.legend = TRUE, legend = "right"),
  nrow = 2,
  labels = ""       # Label of the line plot
)

library(grid)
annotate_figure(figure.support, left = textGrob("Probability of Recieving a Submissive Response", rot = 90, vjust = 0.75, gp = gpar(cex = 1.125)),
                bottom = textGrob("Recipient's Sex", hjust = 0.75, gp = gpar(cex = 1.125)))




### Stratified models

library(broom)
library(broom.mixed)
library(patchwork)

#Stratify on Recipient

adult_female_recip_mating <- glmer(Win ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                              family = binomial(link = "logit"),
                              data = Data_File_Adult[Data_File_Adult$Recipient_Sex == 'Female',],
                              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

adult_nmale_recip_mating <- glmer(Win ~ Actor_Sex + (1|Session),
                            family = binomial(link = "logit"),
                            data = Data_File_Adult[Data_File_Adult$Recipient_Sex == 'Natal_Male',],
                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Pairs lead to overfitting, removed

adult_imale_recip_mating <- glmer(Win ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                           family = binomial(link = "logit"),
                           data = Data_File_Adult[Data_File_Adult$Recipient_Sex == 'Immigrant_Male',],
                           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

cub_female_recip_mating <- glmer(Win ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                           family = binomial(link = "logit"),
                           data = Data_File_Cub[Data_File_Cub$Recipient_Sex == 'Female',],
                           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

cub_male_recip_mating <- glmer(Win ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                          family = binomial(link = "logit"),
                          data = Data_File_Cub[Data_File_Cub$Recipient_Sex == 'Male',],
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

sub_female_recip_mating <- glmer(Win ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                          family = binomial(link = "logit"),
                          data = Data_File_Sub[Data_File_Sub$Recipient_Sex == 'Female',],
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

sub_male_recip_mating <- glmer(Win ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                        family = binomial(link = "logit"),
                        data = Data_File_Sub[Data_File_Sub$Recipient_Sex == 'Male',],
                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


female_recip_mating_est_dominance <- tidy(adult_female_recip_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = c('Adult\\nImmigrant\\nMales', 'Adult\\nNatal\\nMales'), recipient = 'Female', age.class = 'adult')

nmale_recip_mating_est_dominance <- tidy(adult_nmale_recip_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = c('Adult\\nImmigrant\\nMales', 'Adult\\nNatal\\nMales'), recipient = 'Natal Male', age.class = 'adult')

imale_recip_mating_est_dominance <- tidy(adult_imale_recip_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = c('Adult\\nImmigrant\\nMales', 'Adult\\nNatal\\nMales'), recipient = 'Immigrant Male', age.class = 'adult')

cub_female_recip_mating_est_dominance <- tidy(cub_female_recip_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = 'Male\\nCubs', recipient = 'Female', age.class = 'cub')

cub_male_recip_mating_est_dominance <- tidy(cub_male_recip_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = 'Male\\nCubs', recipient = 'Male', age.class = 'cub')

sub_female_recip_mating_est_dominance <- tidy(sub_female_recip_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = 'Male\\nSub-Adults', recipient = 'Female', age.class = 'sub')

sub_male_recip_mating_est_dominance <- tidy(sub_male_recip_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = 'Male\\nSub-Adults', recipient = 'Male', age.class = 'sub')


recip_strat_est_dominance <- rbind(cub_female_recip_mating_est_dominance, cub_male_recip_mating_est_dominance,
                                   sub_female_recip_mating_est_dominance, sub_male_recip_mating_est_dominance,
                                   female_recip_mating_est_dominance, imale_recip_mating_est_dominance, nmale_recip_mating_est_dominance)
recip_strat_est_dominance$plot_order <- letters[1:nrow(recip_strat_est_dominance)]
recip_strat_est_dominance$age.class <- factor(recip_strat_est_dominance$age.class, levels = c('cub', 'sub', 'adult'))


recip_strat_dominance_plot <- ggplot(recip_strat_est_dominance, aes(x = plot_order, y = estimate)) +
  geom_hline(yintercept = 1, color = 'red',
             linetype = 2) + # line at null behind coefs
  annotate('rect', xmin = 0, ymin =  0.0005, xmax = 10.5, ymax = 0.0025)+
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 8.5), linetype = 1, size = 1)+
  geom_point(aes(color = age.class), size = 2, alpha = 1,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(color = age.class, ymin=conf.low, ymax=conf.high), width=.1,
                position=position_dodge(.5)) +
  scale_color_manual(values = c('darkgreen', 'darkblue', 'brown')) +
  theme(text = element_text(size=10)) +
  theme(legend.position = 'none') +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = 'black',
                                 size = 1, linetype = 'solid')) +
  scale_x_discrete(labels = recip_strat_est_dominance$Actor, expand = expansion(add = c(1, 0)))+
  ylab(expression("Odds ratio of eliciting submission"))+
  xlab('Actor')+
  annotate('text', x = 0.05, y = 1.5, label = expression(italic('Female actors')), color = 'red', hjust = 0, size = 8/.pt)+
  annotate('text', x = 0.2, y = exp((log(0.0005) + log(0.0025))/2), label = expression(underline('Model')), angle = 90, color = 'white')+
  annotate('text', x = c(1,2,3,4,5.5, 7.5, 9.5), y = rep(exp((log(0.0005) + log(0.0025))/2), 7),
           label = c('Female\\nRecipients', 'Male\\nRecipients', 'Female\\nRecipients', 'Male\\nRecipients', 'Adult\\nFemale\\nRecipients', 'Adult\\nImmigrant Male\\nRecipients', 'Adult\\nNatal Male\\nRecipients'), color ='white', size = 8/.pt)+
  scale_y_continuous(trans='log10', limits = c(0.0005,20), expand = c(0,0), breaks = c(0.1,0.5, 1, 2, 8))
  

recip_strat_dominance_plot

## Stratify by actor

adult_female_actor_mating <- glmer(Win ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                            family = binomial(link = "logit"),
                            data = Data_File_Adult[Data_File_Adult$Actor_Sex == 'Female',],
                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

adult_nmale_actor_mating <- glmer(Win ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                           family = binomial(link = "logit"),
                           data = Data_File_Adult[Data_File_Adult$Actor_Sex == 'Natal_Male',],
                           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

adult_imale_actor_mating <- glmer(Win ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                           family = binomial(link = "logit"),
                           data = Data_File_Adult[Data_File_Adult$Actor_Sex == 'Immigrant_Male',],
                           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

cub_female_actor_mating <- glmer(Win ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                          family = binomial(link = "logit"),
                          data = Data_File_Cub[Data_File_Cub$Actor_Sex == 'Female',],
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

cub_male_actor_mating <- glmer(Win ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                        family = binomial(link = "logit"),
                        data = Data_File_Cub[Data_File_Cub$Actor_Sex == 'Male',],
                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

sub_female_actor_mating <- glmer(Win ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                          family = binomial(link = "logit"),
                          data = Data_File_Sub[Data_File_Sub$Actor_Sex == 'Female',],
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

sub_male_actor_mating <- glmer(Win ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                        family = binomial(link = "logit"),
                        data = Data_File_Sub[Data_File_Sub$Actor_Sex == 'Male',],
                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))




female_actor_mating_est_dominance <- tidy(adult_female_actor_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = c('Adult\\nImmigrant\\nMales', 'Adult\\nNatal\\nMales'), actor = 'Female', age.class = 'adult')

nmale_actor_mating_est_dominance <- tidy(adult_nmale_actor_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = c('Adult\\nImmigrant\\nMales', 'Adult\\nNatal\\nMales'), actor = 'Natal Male', age.class = 'adult')

imale_actor_mating_est_dominance <- tidy(adult_imale_actor_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = c('Adult\\nImmigrant\\nMales', 'Adult\\nNatal\\nMales'), actor = 'Immigrant Male', age.class = 'adult')

cub_female_actor_mating_est_dominance <- tidy(cub_female_actor_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = 'Male\\nCubs', actor = 'Female', age.class = 'cub')

cub_male_actor_mating_est_dominance <- tidy(cub_male_actor_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = 'Male\\nCubs', actor = 'Male', age.class = 'cub')

sub_female_actor_mating_est_dominance <- tidy(sub_female_actor_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = 'Male\\nSub-Adults', actor = 'Female', age.class = 'sub')

sub_male_actor_mating_est_dominance <- tidy(sub_male_actor_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = 'Male\\nSub-Adults', actor = 'Male', age.class = 'sub')



actor_strat_est_dominance <- rbind(cub_female_actor_mating_est_dominance, cub_male_actor_mating_est_dominance,
                                   sub_female_actor_mating_est_dominance, sub_male_actor_mating_est_dominance,
                                   female_actor_mating_est_dominance, imale_actor_mating_est_dominance, nmale_actor_mating_est_dominance)
str(actor_strat_est_dominance)
actor_strat_est_dominance$plot_order <- letters[1:nrow(actor_strat_est_dominance)]
actor_strat_est_dominance$Recipient[actor_strat_est_dominance$Recipient == "Male\\nCubs"] <- "Males"
actor_strat_est_dominance$Recipient[actor_strat_est_dominance$Recipient == "Male\\nSub-Adults"] <- "Males"
actor_strat_est_dominance$age.class <- factor(actor_strat_est_dominance$age.class, levels = c('cub', 'sub', 'adult'))

actor_strat_dominance_plot <- ggplot(actor_strat_est_dominance, aes(x = plot_order, y = estimate)) +
  geom_hline(yintercept = 1, color = 'red',
             linetype = 2) + # line at null behind coefs
  annotate('rect', xmin = 0, ymin =  0.05, xmax = 10.5, ymax = 0.17)+
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 8.5), linetype = 1, size = 1)+
  geom_point(aes(color = age.class), size = 2, alpha = 1,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(color = age.class, ymin=conf.low, ymax=conf.high), width=.1,
                position=position_dodge(.5)) +
  scale_color_manual(values = c('darkgreen', 'darkblue', 'brown')) +
  theme(text = element_text(size=10)) +
  theme(legend.position = 'none') +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = 'black',
                                 size = 1, linetype = 'solid')) +
  scale_x_discrete(labels = actor_strat_est_dominance$Recipient, expand = expansion(add = c(1, 0)))+
  ylab(expression("Odds ratio of offering submission"))+
  xlab('Recipient')+
  annotate('text', x = 0.05, y = 0.8, label = expression(italic('Female recipients')), color = 'red', hjust = 0, size = 8/.pt)+
  annotate('text', x = 0.2, y = exp((log(0.05) + log(0.17))/2), label = expression(underline('Model')), angle = 90, color = 'white')+
  annotate('text', x = c(1,2,3,4,5.5, 7.5, 9.5), y = rep(exp((log(0.05) + log(0.17))/2), 7),
           label = c('Cub\\nFemale\\nActors', 'Cub\\nMale\\nActors', 'Sub\\nFemale\\nActors', 'Sub\\nMale\\nActors', 'Adult\\nFemale\\nActors', 'Adult\\nImmigrant Male\\nActors', 'Adult\\nNatal Male\\nActors'), color ='white', size = 8/.pt)+
  scale_y_continuous(trans='log10', limits = c(0.05,150), expand = c(0,0), breaks = c(1, 2, 8, 32))
  

pdf(file = 'Figure 3 OR strat.pdf', width = 7.50, height = 9.00)
recip_strat_dominance_plot / actor_strat_dominance_plot
dev.off()
save(actor_strat_est_dominance, recip_strat_est_dominance, file = 'data_Figure_3.Rdata')

write.csv(actor_strat_est_dominance, "C:/Users/skmcc/Documents/R/R_wd/actor_strat_est_dominance.csv")
write.csv(recip_strat_est_dominance, "C:/Users/skmcc/Documents/R/R_wd/recip_strat_est_dominance.csv")


### Social support and dominance models


summary(Model_AllClans_Mating_Supported_Adult)


support_est_mating <- list()

## Cubs
for(actor in c('Female' ,'Male')){
  for(recipient in c('Female', 'Male')){

    df <- Data_File_Cub[Data_File_Cub$Recipient_Sex == recipient &
                            Data_File_Cub$Actor_Sex == actor,]

    df$Supported <- factor(df$Supported, levels = c('Unsupported', 'Supported'))
    support_model_mating <- glmer(Win ~ Supported +
                             (1|Session) + (1|ID:Recipient),
                           family = binomial(link = "logit"),
                           data = df,
                           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

    support_est_mating[[length(support_est_mating) + 1]] <- tidy(support_model_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
      filter(grepl('Support', term)) %>%
      select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
      mutate(Recipient = recipient, Actor = actor, age.class = 'cub')
  }
}

## Subs
for(actor in c('Female' ,'Male')){
  for(recipient in c('Female', 'Male')){

    df <- Data_File_Sub[Data_File_Sub$Recipient_Sex == recipient &
                          Data_File_Sub$Actor_Sex == actor,]

    df$Supported <- factor(df$Supported, levels = c('Unsupported', 'Supported'))
    support_model_mating <- glmer(Win ~ Supported +
                             (1|Session)+(1|ID:Recipient),
                           family = binomial(link = "logit"),
                           data = df,
                           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

    support_est_mating[[length(support_est_mating) + 1]] <- tidy(support_model_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
      filter(grepl('Support', term)) %>%
      select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
      mutate(Recipient = recipient, Actor = actor, age.class = 'sub')
  }
}

## Adults
for(actor in c('Female' ,'Immigrant_Male', 'Natal_Male')){
  for(recipient in c('Female', 'Immigrant_Male', 'Natal_Male')){

    df <- Data_File_Adult[Data_File_Adult$Recipient_Sex == recipient &
                            Data_File_Adult$Actor_Sex == actor,]

    df$Supported <- factor(df$Supported, levels = c('Unsupported', 'Supported'))
    support_model_mating <- glmer(Win ~ Supported +
                            (1|Session),
                          family = binomial(link = "logit"),
                          data = df,
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

    support_est_mating[[length(support_est_mating) + 1]] <- tidy(support_model_mating, effects = 'fixed', conf.int = T, exponentiate = T) %>%
      filter(grepl('Support', term)) %>%
      select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
      mutate(Recipient = recipient, Actor = actor, age.class = 'adult')
  }
}

support_est_mating <- bind_rows(support_est_mating)

support_est_mating$age.class <- factor(support_est_mating$age.class, levels = c('cub', 'sub', 'adult'))

support_est_mating <- arrange(support_est_mating, age.class, Recipient)
support_est_mating$plot_order <- letters[1:nrow(support_est_mating)]

## Limit to intersexual interactions
support_est_mating$interaction.class <- 'intra'
support_est_mating$interaction.class[support_est_mating$Recipient == 'Female' & grepl(pattern = 'Male', ignore.case = F,
                                                                        support_est_mating$Actor)] <- 'inter'
support_est_mating$interaction.class[grepl('Male', support_est_mating$Recipient, ignore.case = F) &
                                support_est_mating$Actor == 'Female'] <- 'inter'

support_est_mating_inter <- filter(support_est_mating, interaction.class == 'inter')


#
support_recip_mating_strat <- ggplot(support_est_mating_inter, aes(x = plot_order, y = estimate)) +
  geom_hline(yintercept = 1, color = 'grey20',
             linetype = 2) + # line at null behind coefs
  annotate('rect', xmin = 0, ymin =  0.001, xmax = 8.5, ymax = 0.005)+
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5), linetype = 1, size = 1)+
  geom_point(aes(color = age.class), size = 2, alpha = 1,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(color = age.class, ymin=conf.low, ymax=conf.high), width=.1,
                position=position_dodge(.5)) +
  scale_color_manual(values = c('darkgreen', 'darkblue','brown')) +
  #coord_flip() + # flip x and y axes
  theme(text = element_text(size=10)) +
  theme(legend.position = 'none') +
  #theme(axis.ticks = element_blank()) + # remove axis ticks
  # remove background color
  theme(panel.background = element_rect(fill = 'white')) +
  # add major axes
  theme(axis.line = element_line(colour = 'black',
                                 size = 1, linetype = 'solid')) +

  scale_x_discrete(labels = c('Male\\nCubs', 'Female\\nCubs', 'Male\\nSub-Adults', 'Female\\nSub-Adults',
                              'Adult\\nImmigrant\\nMales', 'Adult\\nNatal\\nMales', 'Adult\\nFemales', 'Adult\\nFemales'),
                   expand = expansion(add = c(1, 0)))+
  ylab(expression("Odds ratio of eliciting submission"))+
  xlab('Actor')+
  annotate('text', x = 0.1, y = 2.5, label = expression(italic('Unsupported actors')), color = 'grey20', hjust = 0, size = 8/.pt)+
  annotate('text', x = 0.2, y = exp((log(0.001) + log(0.005))/2), label = expression(underline('Model')), angle = 90, color = 'white')+
  annotate('text', x = c(1,2,3,4,5.5, 7, 8), y = rep(exp((log(0.001) + log(0.005))/2), 7),
           label = c('Female\\nRecipients', 'Male\\nRecipients', 'Female\\nRecipients', 'Male\\nRecipients', 'Adult\\nFemale\\nRecipients', 'Adult\\nImmigrant Male\\nRecipients', 'Adult\\nNatal Male\\nRecipients'), color ='white', size = 8/.pt)+
  scale_y_continuous(trans='log10', limits = c(0.001,12), expand = c(0,0), breaks = c(0.1, 0.2, 0.5, 1, 2, 8))




pdf(file = 'Figure 5 intersex intx.pdf', width = 8, height = 4)
support_recip_mating_strat
dev.off()
save(support_est_mating_inter, file = 'data_Figure_5.Rdata')
write.csv(support_est_mating_inter, "C:/Users/skmcc/Documents/R/R_wd/support_est_mating_inter.csv")
write.csv(support_est_mating, "C:/Users/skmcc/Documents/R/R_wd/support_est_mating.csv")



#######Sex Difference in Being Supported Models

Model_Supported_Cub2 <- glmer(Sequence ~ Actor_Sex+Recipient_Sex + (1|ID:Recipient) + (1|Session),
                              family = binomial(link = "logit"),
                              data = Data_File_Cub,
                              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

Model_Supported_Sub2 <- glmer(Sequence ~ Actor_Sex+Recipient_Sex + (1|ID:Recipient) + (1|Session),
                              family = binomial(link = "logit"),
                              data = Data_File_Sub,
                              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

Model_Supported_Adult2 <- glmer(Sequence ~ Actor_Sex+Recipient_Sex + (1|ID:Recipient) + (1|Session),
                                family = binomial(link = "logit"),
                                data = Data_File_Adult,
                                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


#Summarize the results of the best models
summary(Model_Supported_Cub2)
Model_Supported_Cub2_DF <- summary(Model_Supported_Cub2)
write.csv(Model_Supported_Cub2_DF$coefficients, "C:/Users/skmcc/Documents/R/R_wd/Model_Supported_Cub2_DF.csv")
summary(Model_Supported_Sub2)
Model_Supported_Sub2_DF <- summary(Model_Supported_Sub2)
write.csv(Model_Supported_Sub2_DF$coefficients, "C:/Users/skmcc/Documents/R/R_wd/Model_Supported_Sub2_DF.csv")
summary(Model_Supported_Adult2)
Model_Supported_Adult2_DF <- summary(Model_Supported_Adult2)
write.csv(Model_Supported_Adult2_DF$coefficients, "C:/Users/skmcc/Documents/R/R_wd/Model_Supported_Adult2_DF.csv")


#Look at the pairwise comparisons from the full model
emmeans(Model_Supported_Cub2, specs = pairwise ~ Recipient_Sex, type="response", adjust="tukey")
emmeans(Model_Supported_Sub2, specs = pairwise ~ Recipient_Sex, type="response", adjust="tukey")
emmeans(Model_Supported_Adult2, specs = pairwise ~ Recipient_Sex, type="response", adjust="tukey")
emmeans(Model_Supported_Cub2, specs = pairwise ~ Actor_Sex, type="response", adjust="tukey")
emmeans(Model_Supported_Sub2, specs = pairwise ~ Actor_Sex, type="response", adjust="tukey")
emmeans(Model_Supported_Adult2, specs = pairwise ~ Actor_Sex, type="response", adjust="tukey")

##Plot the full model
Cub <- ggeffects::ggemmeans(Model_Supported_Cub2, c("Recipient_Sex","Actor_Sex"))
Cub
plot(Cub)
Cub_Plot <- plot(Cub)+
  ggtitle("B") +
  xlab("Recipient's Sex") +
  ylab("Probability of Being Supported")+
  theme_classic()
Plot_B <- Cub_Plot +labs(colour="Actor's Sex")
Sub <- ggeffects::ggemmeans(Model_Supported_Sub2, c("Recipient_Sex","Actor_Sex"))
Sub
plot(Sub)
Sub_Plot <- plot(Sub)+
  ggtitle("C") +
  xlab("Recipient's Sex") +
  ylab("Probability of Being Supported")+
  theme_classic()
Plot_C <- Sub_Plot +labs(colour="Actor's Sex")
Adult <- ggeffects::ggemmeans(Model_Supported_Adult2, c("Recipient_Sex","Actor_Sex"))
Adult
plot(Adult)
Adult_Plot <- plot(Adult)+
  ggtitle("A") +
  xlab("Recipient's Sex") +
  ylab("Probability of Being Supported")+
  theme_classic()
Plot_A <- Adult_Plot +labs(colour="Actor's Sex")

library(ggpubr)
figure <- ggarrange(
  Plot_A + rremove("ylab") + rremove("xlab"),
  # Second row with box and dot plots
  ggarrange(Plot_B + rremove("ylab") + rremove("xlab"),
            Plot_C+ rremove("ylab") + rremove("xlab"),
            ncol = 2, labels = c("", ""),
            common.legend = TRUE, legend = "right"),
  nrow = 2,
  labels = ""       # Label of the line plot
)

library(grid)
annotate_figure(figure, left = textGrob("Probability of Being Supported", rot = 90, vjust = 0.75, gp = gpar(cex = 1.125)),
                bottom = textGrob("Recipient's Sex", hjust = 0.75, gp = gpar(cex = 1.125)))


### Stratified models by age class and recipient
## Adults
model_support_adult_female_recipient <- glmer(Sequence ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                                              family = binomial(link = "logit"),
                                              data = Data_File_Adult[Data_File_Adult$Recipient_Sex == 'Female',],
                                              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


model_support_adult_nmale_recipient <- glmer(Sequence ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                                             family = binomial(link = "logit"),
                                             data = Data_File_Adult[Data_File_Adult$Recipient_Sex == 'Natal_Male',],
                                             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_support_adult_imale_recipient <- glmer(Sequence ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                                             family = binomial(link = "logit"),
                                             data = Data_File_Adult[Data_File_Adult$Recipient_Sex == 'Immigrant_Male',],
                                             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

## Cubs

model_support_cub_female_recipient <- glmer(Sequence ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                                            family = binomial(link = "logit"),
                                            data = Data_File_Cub[Data_File_Cub$Recipient_Sex == 'Female',],
                                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_support_cub_male_recipient <- glmer(Sequence ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                                          family = binomial(link = "logit"),
                                          data = Data_File_Cub[Data_File_Cub$Recipient_Sex == 'Male',],
                                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

## Subadults

model_support_sub_female_recipient <- glmer(Sequence ~ Actor_Sex + (1|ID:Recipient) + (1|Session),
                                            family = binomial(link = "logit"),
                                            data = Data_File_Sub[Data_File_Sub$Recipient_Sex == 'Female',],
                                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_support_sub_male_recipient <- glmer(Sequence ~ Actor_Sex + (1|ID:Recipient),
                                          family = binomial(link = "logit"),
                                          data = Data_File_Sub[Data_File_Sub$Recipient_Sex == 'Male',],
                                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Model overfit with session, session removed



female_recip_est <- tidy(model_support_adult_female_recipient, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = c('Adult Immigrant Males', 'Adult Natal Males'), recipient = 'Female', age.class = 'adult')

nmale_recip_est <- tidy(model_support_adult_nmale_recipient, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = c('Adult Immigrant Males', 'Adult Natal Males'), recipient = 'Natal Male', age.class = 'adult')

imale_recip_est <- tidy(model_support_adult_imale_recipient, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = c('Adult Immigrant Males', 'Adult Natal Males'), recipient = 'Immigrant Male', age.class = 'adult')

female_cub_recip_est <- tidy(model_support_cub_female_recipient, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = 'Male Cubs', recipient = 'Female', age.class = 'cub')

male_cub_recip_est <- tidy(model_support_cub_male_recipient, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = 'Male Cubs', recipient = 'Male', age.class = 'cub')

female_sub_recip_est <- tidy(model_support_sub_female_recipient, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = 'Male Sub-Adults', recipient = 'Female', age.class = 'sub')

male_sub_recip_est <- tidy(model_support_sub_male_recipient, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Actor = 'Male Sub-Adults', recipient = 'Male', age.class = 'sub')

recip_strat_est <- rbind(female_recip_est, nmale_recip_est, imale_recip_est,
                         female_cub_recip_est, male_cub_recip_est, female_sub_recip_est,
                         male_sub_recip_est)


recip_strat_est$Actor <- factor(recip_strat_est$Actor, levels = c('Male Cubs','Male Sub-Adults','Adult Immigrant Males', 'Adult Natal Males'))

recip_strat_plot <- ggplot(recip_strat_est[], aes(x = Actor, y = estimate)) +
  geom_hline(yintercept = 1, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(aes(color = Actor), size = 6, alpha = 1,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(color = Actor, ymin=conf.low, ymax=conf.high), width=.1,
                position=position_dodge(.5)) +
  scale_color_manual(breaks = c('Male Cubs','Male Sub-Adults','Adult Immigrant Males', 'Adult Natal Males'),
                     values = c('darkgreen', 'darkblue', 'brown','brown')) +
  #coord_flip() + # flip x and y axes
  theme(plot.title = element_text(hjust = 0.5)) + # center title
  theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) +
  # bold and size title and axes labels
  theme(text = element_text(size=20, face = 'bold')) +
  theme(legend.position = 'none') +
  #theme(axis.ticks = element_blank()) + # remove axis ticks
  # remove background color
  theme(panel.background = element_rect(fill = 'white')) +
  # add major axes
  theme(axis.line = element_line(colour = 'black',
                                 size = 1, linetype = 'solid')) +
  # change axes font style, color, size, angle, and margin
  theme(axis.text.x = element_text(face='italic', color='black',
                                   size=15, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = -50, l = 0)),
        axis.text.y = element_text(face='bold', color='black',
                                   size=20, angle=0,
                                   margin = margin(t = 0, r = 0,
                                                   b = 0, l = 10))) +
  
  xlab(expression(italic("(Age stratified models)"))) +
  ylab(expression(atop(bold("Odds ratio (and 95% CI) of receiving social support"),
                       paste(italic("Actor sex relative to female actors"))))) +
  facet_wrap('recipient', scales = 'free')


### Stratified models by age class and actor
## Adults
model_support_adult_female_actor <- glmer(Sequence ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                                          family = binomial(link = "logit"),
                                          data = Data_File_Adult[Data_File_Adult$Actor_Sex == 'Female',],
                                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


model_support_adult_nmale_actor <- glmer(Sequence ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                                         family = binomial(link = "logit"),
                                         data = Data_File_Adult[Data_File_Adult$Actor_Sex == 'Natal_Male',],
                                         control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_support_adult_imale_actor <- glmer(Sequence ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                                         family = binomial(link = "logit"),
                                         data = Data_File_Adult[Data_File_Adult$Actor_Sex == 'Immigrant_Male',],
                                         control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

## Cubs

model_support_cub_female_actor <- glmer(Sequence ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                                        family = binomial(link = "logit"),
                                        data = Data_File_Cub[Data_File_Cub$Actor_Sex == 'Female',],
                                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_support_cub_male_actor <- glmer(Sequence ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                                      family = binomial(link = "logit"),
                                      data = Data_File_Cub[Data_File_Cub$Actor_Sex == 'Male',],
                                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

## Subadults

model_support_sub_female_actor <- glmer(Sequence ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                                        family = binomial(link = "logit"),
                                        data = Data_File_Sub[Data_File_Sub$Actor_Sex == 'Female',],
                                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

model_support_sub_male_actor <- glmer(Sequence ~ Recipient_Sex + (1|ID:Recipient) + (1|Session),
                                      family = binomial(link = "logit"),
                                      data = Data_File_Sub[Data_File_Sub$Actor_Sex == 'Male',],
                                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))




female_actor_est <- tidy(model_support_adult_female_actor, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = c('Adult Immigrant Males', 'Adult Natal Males'), actor = 'Female', age.class = 'adult')

nmale_actor_est <- tidy(model_support_adult_nmale_actor, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = c('Adult Immigrant Males', 'Adult Natal Males'), actor = 'Natal Male', age.class = 'adult')

imale_actor_est <- tidy(model_support_adult_imale_actor, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = c('Adult Immigrant Males', 'Adult Natal Males'), actor = 'Immigrant Male', age.class = 'adult')

female_cub_actor_est <- tidy(model_support_cub_female_actor, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = 'Male Cubs', actor = 'Female', age.class = 'cub')

male_cub_actor_est <- tidy(model_support_cub_male_actor, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = 'Male Cubs', actor = 'Male', age.class = 'cub')

female_sub_actor_est <- tidy(model_support_sub_female_actor, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = 'Male Sub-Adults', actor = 'Female', age.class = 'sub')

male_sub_actor_est <- tidy(model_support_sub_male_actor, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Recipient', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = 'Male Sub-Adults', actor = 'Male', age.class = 'sub')

actor_strat_est <- rbind(female_actor_est, nmale_actor_est, imale_actor_est,
                         female_cub_actor_est, male_cub_actor_est, female_sub_actor_est,
                         male_sub_actor_est)


actor_strat_est$Recipient <- factor(actor_strat_est$Recipient, levels = c('Male Cubs','Male Sub-Adults','Adult Immigrant Males', 'Adult Natal Males'))

actor_strat_plot <- ggplot(actor_strat_est[], aes(x = Recipient, y = estimate)) +
  geom_hline(yintercept = 1, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(aes(color = Recipient), size = 6, alpha = 1,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(color = Recipient, ymin=conf.low, ymax=conf.high), width=.1,
                position=position_dodge(.5)) +
  scale_color_manual(breaks = c('Male Cubs','Male Sub-Adults','Adult Immigrant Males', 'Adult Natal Males'),
                     values = c('darkgreen', 'darkblue', 'brown','brown')) +
  #coord_flip() + # flip x and y axes
  theme(plot.title = element_text(hjust = 0.5)) + # center title
  theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) +
  # bold and size title and axes labels
  theme(text = element_text(size=20, face = 'bold')) +
  theme(legend.position = 'none') +
  #theme(axis.ticks = element_blank()) + # remove axis ticks
  # remove background color
  theme(panel.background = element_rect(fill = 'white')) +
  # add major axes
  theme(axis.line = element_line(colour = 'black',
                                 size = 1, linetype = 'solid')) +
  # change axes font style, color, size, angle, and margin
  theme(axis.text.x = element_text(face='italic', color='black',
                                   size=15, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = -50, l = 0)),
        axis.text.y = element_text(face='bold', color='black',
                                   size=20, angle=0,
                                   margin = margin(t = 0, r = 0,
                                                   b = 0, l = 10))) +
  
  xlab(expression(italic("(Age stratified models)"))) +
  ylab(expression(atop(bold("Odds ratio (and 95% CI) of being targeted by social support"),
                       paste(italic("Recipient sex relative to female actors"))))) +
  facet_wrap('actor', scales = 'free')

### Unstratified models

full_adult_est <- tidy(Model_Supported_Adult2, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Sex', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = c(NA, NA, 'Adult Immigrant Males', 'Adult Natal Males'), Actor = c('Adult Immigrant Males', 'Adult Natal Males', NA, NA), age.class = 'adult')

full_cub_est <- tidy(Model_Supported_Cub2, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Sex', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = c(NA, 'Male Cubs'), Actor = c('Male Cubs', NA), age.class = 'cub')


full_sub_est <- tidy(Model_Supported_Sub2, effects = 'fixed', conf.int = T, exponentiate = T) %>%
  filter(grepl('Sex', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(Recipient = c(NA, 'Male Sub-Adults'), Actor = c('Male Sub-Adults', NA), age.class = 'sub')

unstratified_est <- rbind(full_adult_est, full_cub_est, full_sub_est)

unstratified_est$Recipient <- factor(unstratified_est$Recipient, levels = c(c('Male Cubs','Male Sub-Adults','Adult Immigrant Males', 'Adult Natal Males')))


recip_plot <- ggplot(unstratified_est[!is.na(unstratified_est$Recipient),], aes(x = Recipient, y = estimate)) +
  geom_hline(yintercept = 1, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_vline(xintercept = c(1.5, 2.5), size = 1)+
  geom_point(aes(color = Recipient), size = 2, alpha = 1,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(color = Recipient, ymin=conf.low, ymax=conf.high), width=.1,
                position=position_dodge(.5)) +
  scale_color_manual(breaks = c('Male Cubs','Male Sub-Adults','Adult Immigrant Males', 'Adult Natal Males'),
                     values = c('darkgreen', 'darkblue', 'brown','brown')) +
  #coord_flip() + # flip x and y axes
  theme(text = element_text(size=12)) +
  theme(legend.position = 'none') +
  #theme(axis.ticks = element_blank()) + # remove axis ticks
  # remove background color
  theme(panel.background = element_rect(fill = 'white')) +
  # add major axes
  theme(axis.line = element_line(colour = 'black',
                                 size = 1, linetype = 'solid')) +
  scale_x_discrete(labels = c('Male\\nCubs', 'Male\\nSub-Adults', 'Adult\\nImmigrant\\nMales',
                              'Adult\\nNatal\\nMales'), expand = expansion(add =c(1,0)))+
  ylab(expression("Odds ratio of being attacked with social support"))+
  annotate('text', x = 0.04, y = 1.02, label = expression(italic('Female recipients')), color = 'red', hjust = 0, size = 10/.pt)


unstratified_est$Actor <- factor(unstratified_est$Actor, levels = c(c('Male Cubs','Male Sub-Adults','Adult Immigrant Males', 'Adult Natal Males')))

write.csv(unstratified_est, "C:/Users/skmcc/Documents/R/R_wd/unstratified_est.csv")

actor_plot <- ggplot(unstratified_est[!is.na(unstratified_est$Actor),], aes(x = Actor, y = estimate)) +
  geom_hline(yintercept = 1, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_vline(xintercept = c(1.5, 2.5), size = 1)+
  geom_point(aes(color = Actor), size = 2, alpha = 1,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(color = Actor, ymin=conf.low, ymax=conf.high), width=.1,
                position=position_dodge(.5)) +
  scale_color_manual(breaks = c('Male Cubs','Male Sub-Adults','Adult Immigrant Males', 'Adult Natal Males'),
                     values = c('darkgreen', 'darkblue', 'brown','brown')) +
  #coord_flip() + # flip x and y axes
  theme(text = element_text(size=12)) +
  theme(legend.position = 'none') +
  #theme(axis.ticks = element_blank()) + # remove axis ticks
  # remove background color
  theme(panel.background = element_rect(fill = 'white')) +
  # add major axes
  theme(axis.line = element_line(colour = 'black',
                                 size = 1, linetype = 'solid')) +
  scale_x_discrete(labels = c('Male\\nCubs', 'Male\\nSub-Adults', 'Adult\\nImmigrant\\nMales',
                              'Adult\\nNatal\\nMales'), expand = expansion(add = c(1,0)))+
  ylab(expression("Odds ratio of receiving social support"))+
  annotate('text', x = 0.04, y = 1.1, label = expression(italic('Female actors')), color = 'red', hjust = 0, size = 10/.pt)


pdf(file = 'Figure_4.pdf', width = 10.00, height = 5.00)
actor_plot + recip_plot
dev.off()
