############SETUP############
#set working directory
setwd("~/Desktop/MAC19")

#load required packages
#data read in, wrangling, and cleanup
library(stringr);library(readxl);library(tidyverse);library(openxlsx);library(lubridate);library(janitor);library(ResourceSelection)
#calculating ROC
library(xts);library(TTR)
#linear models and post hoc tests
library(lme4);library(car);library(lmerTest);library(emmeans);library(lsmeans);library(rcompanion)
#data visualization
library(scales);library(patchwork);library(ggpmisc);library(ggplot2);library(magick);library(grid)
############

############TEMPERATURE DATA############
####read in temperature data ####
# acclimatization treatments
acclimatization.temps <-read_xlsx("Data_MASTER.xlsx", sheet = "Acclimatization_Profiles")%>%
  mutate_at("Time", as.character)%>%
  rename("Constant High" = "Constant")%>%
  separate(Time, into=c("garbage", "time"), sep= " ")%>%
  dplyr::select(-garbage)%>%
  unite(Timepoint, Date, time, sep = " ")%>%
  mutate_at("Timepoint", as_datetime)%>%
  gather(Treatment, Temperature, -Timepoint)

# 4 week recovery period 
recovery.temps <-read_xlsx("Data_MASTER.xlsx", sheet= "Recovery_Temperatures")%>%
  gather(Tank, Temperature, -Date, -Time)%>%
  select(everything(), -Time)%>%
  filter(Temperature < 28.5)%>%
  mutate_at("Date", as.Date)

# field temperatures
field.temps<-read_xlsx("Data_MASTER.xlsx", sheet="Field_Temperatures")%>%
  gather(Rack, Temperature, -Date, -Time)%>%
  select(everything(), -Time)%>%
  filter(Temperature >27)%>%
  mutate_at("Date", as.Date)%>%
  filter(Date < as.Date("2019-10-20"))%>%
  group_by(Date)%>%summarise(meantemps=mean(Temperature, na.rm=TRUE), sds=sd(Temperature, na.rm=TRUE))%>%
  mutate(Location = "Field")

#lab temperatures 
lab.temps <-read_xlsx("Data_MASTER.xlsx", sheet="Lab_Temperatures")%>%
  select(everything(), -Tank.22, -Tank.24)%>%
  gather(Tank, Temperature, -Date, -Time)%>%
  filter(Temperature >27)%>%
  select(everything(), -Time)%>%
  mutate_at("Date", as.Date)%>%
  filter(Date < as.Date("2019-10-20"))%>%
  full_join(.,recovery.temps, by=c("Date"="Date", "Temperature"="Temperature","Tank"="Tank"))%>%
  group_by(Date)%>%summarise(meantemps=mean(Temperature, na.rm=TRUE), sds=sd(Temperature, na.rm=TRUE))%>%
  mutate(Location = "Lab")%>%
  full_join(., field.temps, by=c("Date"="Date", "meantemps"="meantemps", "sds"="sds", "Location"="Location"))

####temperature plots ####
# recovery, lab and field temperatures during natural thermal stress
temperatures <- ggplot(lab.temps, aes(colour = Location))+
  scale_colour_manual(values = c("grey70", "black"))+
  scale_fill_manual(values = c("grey70", "black"))+
  geom_line(aes(x=Date, y=meantemps))+
  geom_ribbon(aes(x=Date, ymin=meantemps-sds, ymax=meantemps+sds, fill =Location, colour = NA), alpha = 0.3)+
  theme_classic()+
  ylim(27, 33)+
  scale_x_date(breaks = as.Date(c("2019-07-10","2019-07-15","2019-07-30","2019-08-14","2019-08-30","2019-09-14",
                                  "2019-09-28", "2019-10-10")), labels = date_format("%b %d"), limits =c(as.Date("2019-07-10"), as.Date("2019-10-15")))+
  annotate(geom = "segment", x = as.Date("2019-10-10"), xend = as.Date("2019-10-10"),
           y = 27, yend = 33, colour = "grey60", linetype = "dotted")+
  annotate(geom = "segment", x = as.Date("2019-08-14"), xend = as.Date("2019-08-14"),
           y = 27, yend = 33, colour = "grey60",linetype = "dotted")+
  annotate(geom = "segment", x = as.Date("2019-07-15"), xend = as.Date("2019-07-15"),
           y = 27, yend = 33, colour = "grey60",linetype = "dotted")+
  geom_hline(yintercept=27.5, color="grey30", linetype = "dashed") + 
  geom_hline(yintercept=28.5, color="grey30") +
  ggtitle("c")+
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.6,"line"),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 7),
    axis.title.x = element_text(size = 7, face = "bold")
  );temperatures

# average accilimatization temperatures every 2 hours 
hourseq = seq.POSIXt(min(acclimatization.temps$Timepoint), max(acclimatization.temps$Timepoint), by='2 hours')
temps.hourly = acclimatization.temps %>% group_by(Hourly = cut(Timepoint, breaks=hourseq), Treatment)%>%
  summarise(meantemps=mean(Temperature), sd = sd(Temperature))%>%ungroup(Hourly)%>%
  mutate_at("Hourly", as_datetime)
temps.hourly$Treatment <- factor(temps.hourly$Treatment,
                                 levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))

# expanded plot of acclimatization profiles
profiles<-ggplot(temps.hourly, aes(x = Hourly, y =meantemps, group =Treatment))+
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  scale_fill_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_ribbon(aes(x=Hourly, ymin=meantemps-sd, ymax=meantemps+sd, fill =Treatment), alpha = 0.3)+
  geom_line(aes(color = Treatment))+
  ylab("Temperature (°C)")+
  xlab("Date")+
  theme_classic()+
  scale_y_continuous(limits = c(27, 33),breaks = c(27,28,29,30,31,32,33))+
  scale_x_datetime(breaks = as_datetime(c("2019-07-10", "2019-07-15")), labels = date_format("%b %d"))+
  geom_hline(yintercept=27.5, color="grey30", linetype = "dashed") + 
  geom_hline(yintercept=28.5, color="grey30") +
  ggtitle("b")+
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.6,"line"),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 7, face = "bold"),
    axis.title.y = element_text(size = 7, face = "bold")
  );profiles

# compressed acclimatization profiles for overview plot
profiles.c<-ggplot(temps.hourly, aes(x = Hourly, y =meantemps, group =Treatment))+
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  scale_fill_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_ribbon(aes(x=Hourly, ymin=meantemps-sd, ymax=meantemps+sd, fill =Treatment), alpha = 0.3)+
  geom_line(aes(color = Treatment))+
  ylab("Temperature (°C)")+
  xlab("Date")+
  theme_classic()+
  scale_y_continuous(limits = c(27, 33),breaks = c(27,28,29,30,31,32,33))+
  scale_x_datetime(limits = as_datetime(c("2019-07-10", "2019-10-15")),
                   breaks = as_datetime(c("2019-07-10","2019-07-15","2019-07-30","2019-08-14","2019-08-30","2019-09-14",
                                          "2019-09-28", "2019-10-10")), labels = date_format("%b %d"))+
  geom_hline(yintercept=27.5, color="grey30", linetype = "dashed") + 
  geom_hline(yintercept=28.5, color="grey30") +
  annotate(geom = "segment", x = as_datetime("2019-07-10"), xend = as_datetime("2019-07-10"),
           y = 27, yend = 33, colour = "grey60", linetype = "dotted")+
  annotate(geom = "segment", x = as_datetime("2019-07-15"), xend = as_datetime("2019-07-15"),
           y = 27, yend = 33, colour = "grey60",linetype = "dotted")+
  ggtitle("c")+
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 7,angle = 90),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 7, face = "bold"),
    axis.title.y = element_text(size = 7, face = "bold")
  );profiles.c

# compile into one figure with space for map recovery/field temperatures to be added in photoshop
library(patchwork)
overview.fig <- (plot_spacer() | profiles) / profiles.c
overview.fig

####calculate average rate of change of acclimatization profiles ####
#get acclimatization temperatures and convert to .xts format
acclim.data <- acclimatization.temps%>%select(everything())
acclim.xts <- read.zoo(file = acclim.data,format = "%m/%d/%y %H:%M",FUN=as.POSIXct, split = "Treatment")

#calculate ROC of each treatmwnt per 4 15 minute increments (1 hour)
acclim.xts$constantROC <-ROC(acclim.xts$`Constant High`, n = 4)
acclim.xts$controlROC <-ROC(acclim.xts$`Control`, n = 4)
acclim.xts$pulseROC <-ROC(acclim.xts$`Pulse`, n = 4)
acclim.xts$pulsehighROC <-ROC(acclim.xts$`Pulse High`, n = 4)
acclim.xts$pulseincROC <-ROC(acclim.xts$`Pulse Increase`, n = 4)

#create dataframe with mean ROC per treatment
acclimROCs <- as.data.frame(acclim.xts)%>%select(everything(), -`Constant High`, -Control, -Pulse,
                                                 -`Pulse Increase`, -`Pulse High`)%>%
  mutate(`Constant High` = mean(abs(constantROC), na.rm = TRUE))%>%
  mutate(`Control` = mean(abs(controlROC), na.rm = TRUE))%>%
  mutate(`Pulse` = mean(abs(pulseROC), na.rm = TRUE))%>%
  mutate(`Pulse High` = mean(abs(pulsehighROC), na.rm = TRUE))%>%
  mutate(`Pulse Increase` = mean(abs(pulseincROC), na.rm = TRUE))%>%
  select(`Constant High`, Control, Pulse, `Pulse Increase`, `Pulse High`)%>%
  gather(Treatment, allROC)%>%group_by(Treatment)%>%summarise(ROC = mean(allROC))
############


############PAM DATA############
####read in PAM data ####
# acclimatization data, select relevant timepoints
acclimatization <- read.xlsx("Data_MASTER.xlsx", sheet = "Treatment_PAM")%>%
  dplyr::select(everything(), -Tank, -Row, -Column, -T0_record, -T1_record, -T2_record, -T3_record,
                -T4_record, -T5_record, -T6_record, -T7_record, -T8_record, -T9_record)%>%
  rename(Rack = Rack.Tag)%>%
  mutate_at("Rack", as.factor)%>%
  mutate_at("Treatment", as.factor)%>%
  mutate_at("Fragment", as.factor)%>%
  mutate_at("Colony", as.factor)%>%
  mutate_at("Phenotype", as.factor)%>%
  dplyr::select(Treatment, Fragment, Colony, Phenotype,T0_fvfm,T1_fvfm,T5_fvfm)%>%
  na.omit()

# lab coral data, select relevant timepoints
bleaching_lab <- read.xlsx("Data_MASTER.xlsx", sheet = "Lab_PAM")%>%
  dplyr::select(everything(), -Tank, -Row, -Column, -T10_record,-T11_record,-T12_record,-T13_record_1,-T13_record_2)%>%
  dplyr::rename(Rack=Rack.Tag)%>%
  mutate_at("Rack", as.factor)%>%
  mutate_at("Treatment", as.factor)%>%
  mutate_at("Fragment", as.factor)%>%
  mutate_at("Colony", as.factor)%>%
  mutate_at("Phenotype", as.factor)%>%
  dplyr::select(Treatment, Fragment, Colony, Phenotype, T9_fvfm,T12_fvfm)%>%
  mutate(Location = "lab")%>%
  na.omit()

# field coral data, select relevant timepoints and join with other PAM data, filter out junk data
bleaching <- read.xlsx("Data_MASTER.xlsx", sheet = "Field_PAM")%>%
  dplyr::select(everything(), -Rack,-Row, -Column, -T10_record,-T11_record,-T12_record,-T13_record_1,-T13_record_2)%>%
  dplyr::rename(Rack=Rack.Tag)%>%
  mutate_at("Rack", as.factor)%>%
  mutate_at("Treatment", as.factor)%>%
  mutate_at("Fragment", as.factor)%>%
  mutate_at("Colony", as.factor)%>%
  mutate_at("Phenotype", as.factor)%>%
  dplyr::select(Treatment, Fragment, Colony, Phenotype,T9_fvfm,T12_fvfm)%>%
  mutate(Location = "field")%>%
  na.omit()%>%
  full_join(.,bleaching_lab, by=c('Fragment'='Fragment','Colony'='Colony',
                                  'Phenotype'='Phenotype','Treatment'='Treatment',
                                  "T12_fvfm"="T12_fvfm", "Location"="Location", 
                                  "T9_fvfm"="T9_fvfm"))%>%
  full_join(.,acclimatization, by=c('Fragment'='Fragment','Colony'='Colony',
                                    'Phenotype'='Phenotype','Treatment'='Treatment'))%>%
  mutate_at("Fragment", as.factor)%>%
  mutate_at("T12_fvfm", as.numeric)%>%
  mutate_at("T1_fvfm", as.numeric)%>%
  mutate_at("T0_fvfm", as.numeric)%>%
  mutate_at("T9_fvfm", as.numeric)%>%
  mutate_at("T5_fvfm", as.numeric)%>%
  rename(T7 = T12_fvfm)%>%
  rename(T1 = T1_fvfm)%>%
  rename(T4 = T9_fvfm)%>%
  rename(T2 = T5_fvfm)%>%
  mutate_at("Location", as.factor)%>%
  na.omit()%>%
  gather(Date, fvfm, -Treatment, -Fragment, -Colony, -Phenotype, -Location)%>%
  filter(fvfm > 0.52)%>%
  filter(fvfm < 0.84)%>%
  spread(Date, fvfm)

####models for timeseries ####
#select relevant data and normalize to initial timepoint
bleachingmodel = bleaching%>%select(Fragment, Phenotype, Treatment, T1, T2, T4, T7, Colony, Location)%>%
  group_by(Fragment)%>%
  mutate(T1norm= T1/T1)%>%
  mutate(T2norm = T2/T1)%>%
  mutate(T4norm = T4/T1)%>%
  mutate(T7norm = T7/T1)%>%
  ungroup(Fragment)%>%
  select(everything(),-T1, -T2, -T4, -T7)%>%
  gather(Timepoint, fvfm, -Phenotype, -Fragment, -Treatment, -Colony, -Location)%>%
  filter(Timepoint != "T1norm")%>%
  na.omit()

#model with colony effects
overall.model = lmer(fvfm ~ Treatment*Timepoint*Colony +(1|Phenotype/Colony), data=bleachingmodel)

#check assumptions
qqPlot(residuals(overall.model))
leveneTest(residuals(overall.model)~bleachingmodel$Treatment*bleachingmodel$Timepoint*bleachingmodel$Colony)

summary(overall.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
rand(overall.model) #look at p values for random effects - whether or not they are significant doesn’t necessarily matter for your analysis b/c the model is accounting for them. But this can be informative depending on your question. 
anova(overall.model, type=3) #use this to view significant effects/interactions and get p-values

#post hoc testing
emm=emmeans(overall.model, ~Timepoint)
pairs(emm)

emm=emmeans(overall.model, ~Treatment)
pairs(emm)

emm=emmeans(overall.model, ~Treatment|Timepoint)
pairs(emm)

emm=emmeans(overall.model, ~Treatment|Colony)
pairs(emm)
#extract group letters from this interaction for data visualization
analysis.letters <- multcomp::cld(emm, Letters=c(LETTERS))

#model with phenotype effects
# need to transform data - boxcox transformation
likelihood<-MASS::boxcox(fvfm~Treatment*Timepoint*Phenotype, lambda=seq(-2,2,0.1), plotit=TRUE, data=bleachingmodel)
lambda<-likelihood$x
lik<-likelihood$y
bc<-cbind(lambda,lik)
sorted<-bc[order(-lik),]
head(sorted, n = 10)

bleachingmodel$tdata = log(bleachingmodel$fvfm)
phenotype.model = lmer(tdata ~ Treatment*Timepoint*Phenotype + (1|Colony), 
                       data = bleachingmodel, control = lmerControl(optimizer ="Nelder_Mead"))

#check assumptions - homogeneity of variances violated
qqPlot(residuals(phenotype.model)) #normality
leveneTest(residuals(phenotype.model)~bleachingmodel$Treatment*bleachingmodel$Timepoint*bleachingmodel$Phenotype)


####timeseries PAM plot and example colony plot####
#select relevant data and convert to date format
pam.plot.data <- bleaching%>%dplyr::select(T1, T2, T4, T7, Treatment, Phenotype, Fragment, Colony)%>%
  group_by(Fragment)%>%
  mutate(T1norm= T1/T1)%>%
  mutate(T2norm = T2/T1)%>%
  mutate(T4norm = T4/T1)%>%
  mutate(T7norm = T7/T1)%>%
  ungroup(Fragment)%>%
  select(everything(),-T1, -T2, -T4, -T7)%>%
  rename(`2019-07-09`=T1norm)%>%
  rename(`2019-07-14`=T2norm)%>%
  rename(`2019-08-14`=T4norm)%>%
  rename(`2019-10-10`=T7norm)%>%
  gather(Date, fvfm, -Phenotype, -Treatment, -Fragment, -Colony)%>%
  mutate_at("Date", as.Date)

#summarize by significant effects: timepoint and treatment
pam.plot<- plyr::ddply(pam.plot.data, c("Date","Treatment"), summarise,
                       N    = length(fvfm[!is.na(fvfm)]),
                       mean = mean(fvfm, na.rm=TRUE),
                       sd   = sd(fvfm, na.rm=TRUE),
                       se   = sd / sqrt(N)
);pam.plot
pam.plot$Treatment <- factor(pam.plot$Treatment,
                             levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))

#plot of normalized PAM data over time
pamPlot<-ggplot(data=pam.plot, aes(x=Date, y=mean, colour=Treatment)) +
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, group=Treatment), width=0.05, size=0.5, position=position_dodge(1.5), linetype=1)+
  theme_classic()+
  ylab("Normalized Fv/Fm")+
  scale_x_date(limits = as.Date(c("2019-07-09", "2019-10-10")),
               breaks = as.Date(c("2019-07-10","2019-07-15", "2019-08-14", "2019-10-10")), labels = date_format("%b %d"))+
  xlab("Date")+
  annotate("text", x = as.Date("2019-08-17"), y = 0.85, label = "Treatment p<0.001", size =2, fontface = 3)+
  annotate("text", x = as.Date("2019-08-17"), y = 0.84, label = "Timepoint p<0.001", size =2, fontface = 3)+
  annotate("text", x = as.Date("2019-08-17"), y = 0.83, label = "Treatment * Timepoint p<0.001", size =2, fontface = 3)+
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size =7, face = "bold"),
    axis.title.y = element_text(face = "bold",size=7),
    axis.title.x = element_text(face = "bold",size=7),
    axis.text.x = element_text(angle = 90, size=7),
    axis.text.y = element_text(angle = 90, size=7)
  );pamPlot

#select data for example colonies and extract letters information from post-hoc test
colony.plot.data <- bleaching%>%dplyr::select(T1, T2, T4, T7, Treatment, Phenotype, Fragment, Colony)%>%
  group_by(Fragment)%>%
  mutate(T1norm= T1/T1)%>%
  mutate(T2norm = T2/T1)%>%
  mutate(T4norm = T4/T1)%>%
  mutate(T7norm = T7/T1)%>%
  ungroup(Fragment)%>%
  select(everything(),-T1, -T2, -T4, -T7, -T1norm)%>%
  rename(`2019-07-14`=T2norm)%>%
  rename(`2019-08-14`=T4norm)%>%
  rename(`2019-10-10`=T7norm)%>%
  gather(Date, fvfm, -Phenotype, -Treatment, -Fragment, -Colony)%>%
  filter(Colony == "201" |Colony == "202")%>%
  mutate_at("Date", as.Date)

# select letters from post hoc testing and relabel colonies more clearly
plotlabels <- analysis.letters%>%select(Treatment,Colony, .group)%>%
  filter(Colony == "201"|Colony == "202")
plotlabels$Colony <- factor(plotlabels$Colony, levels = c("201", "202"),
                            labels = c("Colony 201", "Colony 202"))

#reorder treatment levels for correct color sequence and relabel colonies more clearly
colony.plot.data$Treatment <- factor(colony.plot.data$Treatment,
                                     levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))
colony.plot.data$Colony <- factor(colony.plot.data$Colony, levels = c("201", "202"),
                                  labels = c("Colony 201", "Colony 202"))
#boxplot of example colonies
colonyPlot<-ggplot(data=colony.plot.data, aes(x=Treatment, y=fvfm)) +
  facet_wrap(~Colony, nrow =2)+
  scale_fill_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_boxplot(aes(fill= Treatment,middle = mean(fvfm)), show.legend = FALSE)+
  theme_classic()+
  geom_text(data=plotlabels, aes(label=.group, y= 0.64), size = 2.5)+
  theme_classic()+
  ylab("Normalized Fv/Fm")+
  theme(
    axis.text.y = element_text(size=7),
    axis.title.y = element_text(size=7, face= "bold"),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  );colonyPlot

#compile into one figure
fvfm <- pamPlot|colonyPlot
fvfm + plot_annotation(tag_levels = 'a') + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

####model for peak temperature stress timepoint ####
# select relevant data
peakmodel = bleaching%>%select(Fragment, Phenotype, Treatment, T1, T2, T4, T7, Colony, Location)%>%
  group_by(Fragment)%>%
  mutate(T1norm= T1/T1)%>%
  mutate(T2norm = T2/T1)%>%
  mutate(T4norm = T4/T1)%>%
  mutate(T7norm = T7/T1)%>%
  ungroup(Fragment)%>%
  select(everything(),-T1, -T2, -T4, -T7)%>%
  gather(Timepoint, fvfm, -Phenotype, -Fragment, -Treatment, -Colony, -Location)%>%
  filter(Timepoint == "T7norm")%>%
  na.omit()

# model with phenotype as main effect
peak.model = lmer(fvfm ~ Treatment*Phenotype +(1|Colony) +(1|Location), data=peakmodel)

# check assumptions
qqPlot(residuals(peak.model)) #normality
leveneTest(residuals(peak.model)~peakmodel$Treatment*peakmodel$Phenotype)

summary(peak.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
rand(peak.model) #look at p values for random effects - whether or not they are significant doesn’t necessarily matter for your analysis b/c the model is accounting for them. But this can be informative depending on your question. 
anova(peak.model, type=3)#use this to view significant effects/interactions and get p-values

# post hoc testing
emm=emmeans(peak.model, ~Treatment|Phenotype)
pairs(emm)
# extract group letters for data visualization
analysis.letters2 <- multcomp::cld(emm, Letters=c(LETTERS))


####boxplot of treatment differences by phenotype ####
# select letters from post hoc testing for visualuzation
plotlabels2 <- analysis.letters2%>%select(Treatment,Phenotype, .group)
plotlabels2$Treatment <-factor(plotlabels2$Treatment, levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))

#rename phenotype labels to make them more clear
peakmodel$Phenotype <- factor(peakmodel$Phenotype, levels = c("B", "NB"),labels = c("Bleached", "Nonbleached"))
plotlabels2$Phenotype <- factor(plotlabels2$Phenotype, levels = c("B", "NB"),labels = c("Bleached", "Nonbleached"))

fvfm.plot <- ggplot(data = peakmodel, aes(x = Treatment, y = fvfm))+
  facet_wrap(~Phenotype)+
  scale_fill_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_boxplot(aes(fill= Treatment,middle = mean(fvfm)),show.legend = FALSE)+
  ylab("Normalized Fv/Fm")+
  theme_classic()+
  geom_text(data=plotlabels2, aes(label=.group, y= 0.67), size = 2.5)+
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=7),
    axis.title.y = element_text(face = "bold", size = 7),
  );fvfm.plot

####regression of ROC and DHW against PAM data ####
#calculate average fvfm at peak timepoint in each treatment
PAMregression <- peakmodel%>%
  select(Treatment, fvfm)%>%
  group_by(Treatment)%>%summarise(mean.fvfm = mean(fvfm, na.rm = TRUE))

# join with ROC data calculated earlier
ROCvPAM <- acclimROCs%>%
  full_join(.,PAMregression, by = c("Treatment"="Treatment"))

# run linear regression and extract relevant parameters
PAM.model <- lm(mean.fvfm~ ROC, data =ROCvPAM)
summary(PAM.model)  
coef(PAM.model)

ROCvPAM$Treatment <- factor(ROCvPAM$Treatment,
                            levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))

# plot regression with ROC
ROCreg <- ggplot(ROCvPAM, colour = Treatment)+
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_point(aes(x=ROC, y = mean.fvfm, color = Treatment))+
  theme_classic()+
  ylab("Normalized Fv/Fm")+
  xlab("Mean ROC")+ 
  annotate("text", x = 0.0048, y = 0.87, label = "p<0.05", size =2, fontface = 3)+
  annotate("text", x = 0.0048, y = 0.865, label = "R-squared=0.88", size =2, fontface = 3)+
  geom_abline(aes(intercept =0.9219403, slope=-4.3877931), colour = "grey60")+
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 7, face= "bold"),
    legend.text = element_text(size =7),
    axis.title.y = element_text(face = "bold", size =7),
    axis.title.x = element_text(face = "bold", size =7),
    axis.text.y = element_text(size =7),
    axis.text.x = element_text(size =7),
    );ROCreg

# dataframe with average fvfm data and approximated DHWs per treatment
DHWvPAM <- PAMregression%>%
  mutate(DHW = c(0.7, 3.2, 2.05, 2.45, 1.45))

DHWvPAM$Treatment <- factor(DHWvPAM$Treatment,
                            levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))

# run linear regression and extract relevant parameters
DHW.model <- lm(mean.fvfm~ DHW, data =DHWvPAM)
summary(DHW.model)  
coef(DHW.model)

# plot regression with DHW
DHWreg <- ggplot(DHWvPAM, colour = Treatment)+
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_point(aes(x=DHW, y = mean.fvfm, color = Treatment))+
  theme_classic()+
  ylab("Normalized Fv/Fm")+
  xlab("DHWs")+ 
  annotate("text", x = 1.1, y = 0.87, label = "p>0.05", size =2, fontface = 3)+
  annotate("text", x =1.1, y = 0.865, label = "R-squared=0.013", size =2, fontface = 3)+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 7, face= "bold"),
        legend.text = element_text(size =7),
        axis.title.y = element_text(face = "bold", size =7),
        axis.title.x = element_text(face = "bold", size =7),
        axis.text.y = element_text(size =7),
        axis.text.x = element_text(size =7))+
  geom_abline(aes(intercept =0.882958314, slope=0.002862967), colour = "grey60");DHWreg

#### make composite figure with peak boxplot and regressions
regressions = fvfm.plot/ (DHWreg | ROCreg)
regressions + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = 'bottom')


####timeseries PAM by colony ####
#extract relevant data
supp.data <-bleaching%>%select(T1, T4, T7, T2, Fragment, Phenotype, Treatment, Colony, Location)%>%
  mutate(T1norm= T1/T1)%>%
  mutate(T2norm = T2/T1)%>%
  mutate(T4norm = T4/T1)%>%
  mutate(T7norm = T7/T1)%>%
  select(everything(), -T1, -T2, -T4, -T7)%>%
  rename(`2019-07-09` = T1norm)%>%
  rename(`2019-07-14` = T2norm)%>%
  rename(`2019-08-14` = T4norm)%>%
  rename(`2019-10-10` = T7norm)%>%
  gather(Date, fvfm, -Fragment, -Phenotype, -Treatment, -Colony, -Location)%>%
  mutate_at("Date", as.Date)

#reorder factor levels
supp.data$Treatment <- factor(supp.data$Treatment,
                             levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))
supp.data$Colony <- factor(supp.data$Colony,
                          levels = c("11", "19", "201", "203", "211","12","20","202","214", "222"))
supp.data$Phenotype <- factor(supp.data$Phenotype, levels = c("B", "NB"),
                             labels = c("Bleached", "Nonbleached"))

#summarize by timepoint, treatment, phenotype, and colony
PAM.supp<- plyr::ddply(supp.data, c("Date", "Treatment", "Phenotype", "Colony"), summarise,
                  N    = length(fvfm[!is.na(fvfm)]),
                  mean = mean(fvfm, na.rm=TRUE),
                  sd   = sd(fvfm, na.rm=TRUE),
                  se   = sd / sqrt(N) 
);PAM.supp

#plot of normalized PAM data over time by colony
PAMPlot.supp<-ggplot(data=PAM.supp, aes(x=Date, y=mean, colour=Treatment)) +
  facet_wrap(~Colony, nrow =2)+
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  scale_linetype_manual(values=c("dotted", "solid"))+
  geom_line(aes(group=interaction(Treatment, Phenotype), linetype=Phenotype), position=position_dodge(1.5)) +
  geom_point(aes(color=Treatment, shape=Phenotype, group=interaction(Treatment, Phenotype)), position=position_dodge(1.5)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, group=interaction(Treatment, Phenotype)), width=0.05, size=0.5, position=position_dodge(1.5), linetype=1)+
  theme_classic()+
  ylab("Normalized Fv/Fm") +
  scale_x_date(breaks = as.Date(c( "2019-07-15","2019-08-14", "2019-10-10")), labels = date_format("%b %d"))+
  xlab("Date")+
  theme(
    legend.position = "bottom",
    text = element_text(size=7),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90)
  );PAMPlot.supp
############


############qPCR DATA############
# references:
#Cunning et al. 2012 - Excess algal symbionts increase the susceptibility of reef corals to bleaching
#Cunning et al. 2017 - Patterns of bleaching and recovery of Montipora capitata in Kaneohe Bay, Hawaii, USA

####upload and format raw qPCR data ####
#upload, cleanup 
rawdata<-read_excel("Data_MASTER.xlsx",sheet = "qPCR")                                         
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them
format_data$sample_name <- gsub('\\\\.',"-", format_data$sample_name);head(format_data)

#set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

#process data
#read in information about samples
frag_info<-read_excel("Data_MASTER.xlsx", sheet="Large_Frags")%>%
  dplyr::select(everything(), -Tag)%>%
  mutate_at("Treatment", as.factor)%>%
  mutate_at("Phenotype", as.factor)%>%
  mutate_at("Fragment", as.factor)%>%
  mutate_at("Tank", as.factor)%>%
  mutate_at("Colony", as.factor)

frag_info_field<-read_excel("Data_MASTER.xlsx", sheet="Small_Frags")%>%
  dplyr::select(everything(), -Tag)%>%
  rename(Fragment = Label)%>%
  mutate_at("Treatment", as.factor)%>%
  mutate_at("Phenotype", as.factor)%>%
  mutate_at("Fragment", as.factor)%>%
  mutate_at("Tank", as.factor)%>%
  mutate_at("Colony", as.factor)%>%
  full_join(., frag_info, by = c("Treatment"="Treatment", "Phenotype"="Phenotype","Fragment"="Fragment",
                                 "Tank"="Tank", "Colony"="Colony", "Fragment"="Fragment"))

data<-format_data%>%  
  dplyr::select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                       ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                   ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("con",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
Mcap<-filter(data,target_name=="Mcap")

#final formating of data, join with sample information, filter out junk data
final_data<-left_join(left_join(C,D,by="sample_name"),Mcap,by="sample_name")%>%
  dplyr::select(-target_name.x,-target_name.y,-target_name)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename(d_mean=ct_mean.y)%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  dplyr::rename(mcap_mean=ct_mean)%>%
  dplyr::rename(mcap_sd=ct_sd)%>%                                          
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(mcap_mean=mcap_mean-fluo.Mcap)%>%                                                 #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0~"CD",is.na(d_mean)~"C"))%>%                                                                             ###set presence/absence in new column based on ct values
  mutate(dm_ratio=(2^(mcap_mean-d_mean))/copy.n.ratioDM)%>%                                #from Cunning et al. 2012
  mutate(cm_ratio=(2^(mcap_mean-c_mean))/copy.n.ratioCM)%>%                                #from Cunning et al. 2012
  replace_na(list(dm_ratio = 0, cm_ratio = 0))%>%   
  mutate(cd_ratio=(cm_ratio/dm_ratio))%>%                                                  #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio                                                    ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(is.na(prop_c)~"1",TRUE ~ as.character(prop_c)))%>%               #from Cunning et al. 2012
  mutate(prop_c=as.numeric(prop_c))%>%
  mutate(prop_d=1-prop_c)%>%                                                                                                                 ###set proportion D off of proportion C
  mutate(abundance=case_when(prop_d>prop_c~"D>C",prop_c>prop_d~"C>D"))%>%                                                                    ###set abundance based on which proportion is dominant, create character factor
  mutate(sd_warning=case_when((c_sd>1&d_sd>1)~"cd*",c_sd>1~"c*",d_sd>1~"d*"))%>%                                                                                   ###set warnings where standard deviation of tech replicates is >1
  mutate(ct_warning=case_when((c_mean>38&d_mean>38)~"cd*",c_mean>38~"c*",d_mean>38~"d*"))%>%                                                                             ###set warnings where ct values is later than 34
  mutate(rep=case_when((is.na(d_sd)&d_mean>0)|(is.na(c_sd)&c_mean>0)~"1replicate"))%>%                                                       ###set warning if only one replicate processed (no standard deviation but mean present)
  separate(sample_name,into=c("Timepoint","Fragment"),sep="-")%>%
  mutate_at("Timepoint", as.factor)%>%
  mutate_at("Fragment", as.factor)%>%
  mutate(C=as.integer(prop_c*100))%>%
  mutate(D=100-C)%>%
  filter(is.na(rep))%>%
  left_join(.,frag_info_field, by="Fragment")%>%
  mutate_at("Fragment", as.factor)%>%
  filter(mcap_mean<34.08)%>%
  mutate_at("Fragment", as.character)%>%
  mutate_at("Fragment", as.numeric)%>%
  mutate(Location = case_when((Fragment<1000)~"Field", (Fragment>999)~"Lab"))%>%
  filter(Fragment !=1046)

####create log-log plots of symbiont community at each timepoint ####
#select relevant data for each timepoint
pre.data <-final_data%>%select(C,D,cm_ratio, dm_ratio, Colony, Phenotype, Timepoint, Fragment, Treatment, Location)%>%
  filter(Location == "Lab")%>%
  mutate_at("Colony", as.factor)%>%
  mutate_at("Phenotype", as.factor)%>%
  replace_na(list(dm_ratio = 0, cm_ratio = 0))%>%
  filter(Timepoint == "T1")
pre.data$Treatment <- factor(pre.data$Treatment,
                             levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))
pre.data$Phenotype <- factor(pre.data$Phenotype,
                            labels = c("Bleached", "Nonbleached"))

post.data <-final_data%>%select(C,D,cm_ratio, dm_ratio, Colony, Phenotype, Timepoint, Fragment, Treatment, Location)%>%
  filter(Location == "Lab")%>%
  mutate_at("Colony", as.factor)%>%
  mutate_at("Phenotype", as.factor)%>%
  replace_na(list(dm_ratio = 0, cm_ratio = 0))%>%
  filter(Timepoint == "T2")
post.data$Treatment <- factor(post.data$Treatment,
                              levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))
post.data$Phenotype <- factor(post.data$Phenotype,
                            labels = c("Bleached", "Nonbleached"))

recov.data <-final_data%>%select(C,D,cm_ratio, dm_ratio, Colony, Phenotype, Timepoint, Fragment, Treatment, Location)%>%
  filter(Location == "Lab")%>%
  mutate_at("Colony", as.factor)%>%
  mutate_at("Phenotype", as.factor)%>%
  replace_na(list(dm_ratio = 0, cm_ratio = 0))%>%
  filter(Timepoint == "T4")
recov.data$Treatment <- factor(recov.data$Treatment,
                               levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))
recov.data$Phenotype <- factor(recov.data$Phenotype,
                            labels = c("Bleached", "Nonbleached"))

bl.data <-final_data%>%select(cm_ratio, dm_ratio, Colony, Phenotype, Timepoint, Fragment, Treatment, Location)%>%
  filter(Location == "Lab")%>%
  mutate_at("Colony", as.factor)%>%
  mutate_at("Phenotype", as.factor)%>%
  replace_na(list(dm_ratio = 0, cm_ratio = 0))%>%
  filter(Timepoint == "T7")%>%
  mutate_at("Fragment", as.character)%>%
  mutate_at("Fragment", as.numeric)%>%
  mutate(diff = cm_ratio-dm_ratio)

bl.data$Treatment <- factor(bl.data$Treatment,
                            levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))
bl.data$Phenotype <- factor(bl.data$Phenotype,
                           labels = c("Bleached", "Nonbleached"))

#create log-log plot for each timepoint
pre.acclim<- ggplot(pre.data, aes(x = log10(dm_ratio), y = log10(cm_ratio), colour = Treatment, shape = Phenotype))+
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_point(aes(group = Phenotype),size = 2)+
  theme_classic()+
  ylim(-6, 1)+
  xlim(-6, 1)+
  coord_fixed()+
  geom_abline(linetype = "dashed")+
  ylab("log(C:H)")+
  xlab("log(D:H)")+
  theme(legend.position = "right",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)
  );pre.acclim

post.acclim<- ggplot(post.data, aes(x = log10(dm_ratio), y = log10(cm_ratio), colour = Treatment, shape = Phenotype))+
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_point(aes(group = Phenotype), size = 2)+
  theme_classic()+
  ylim(-6, 1)+
  xlim(-6, 1)+
  coord_fixed()+
  geom_abline(linetype = "dashed")+
  ylab("log(C:H)")+
  xlab("log(D:H)")+
  theme(legend.position = "right",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)
  );post.acclim

recovery<- ggplot(recov.data, aes(x = log10(dm_ratio), y = log10(cm_ratio), colour = Treatment, shape = Phenotype))+
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_point(aes(group = Phenotype), size = 2)+
  theme_classic()+
  ylim(-6, 1)+
  xlim(-6, 1)+
  coord_fixed()+
  geom_abline(linetype = "dashed")+
  ylab("log(C:H)")+
  xlab("log(D:H)")+
  theme(legend.position = "right",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)
  );recovery

bleach.symb<- ggplot(bl.data, aes(x = log10(dm_ratio), y = log10(cm_ratio),colour = Treatment, shape = Phenotype))+
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_point(aes(group = Phenotype), size = 2)+
  theme_classic()+
  ylim(-6, 1)+
  xlim(-6, 1)+
  coord_fixed()+
  geom_abline(linetype = "dashed")+
  ylab("log(C:H)")+
  xlab("log(D:H)")+
  theme(legend.position = "right",
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)
  );bleach.symb

# compile into a single figure
loglog <- (pre.acclim | post.acclim )/ (recovery | bleach.symb)
loglog + plot_annotation(tag_levels = 'a') + plot_layout(guides = "collect") & theme(legend.spacing = unit(0.1, 'cm'))

####model effect of treatment, timepoint on Cladocopium in mixed and C-only communities ####
#select relevant data
model.data <- final_data%>%dplyr::select(Phenotype, Colony, Treatment, Timepoint, cm_ratio, dm_ratio, Fragment, Location)%>%
  rename(Durusdinium = dm_ratio)%>%
  rename(Cladocopium = cm_ratio)%>%
  gather(Symbiont, Ratio, -Phenotype, -Colony, -Treatment, -Timepoint, -Fragment, -Location)%>%
  mutate_at("Location", as.factor)%>%
  mutate_at("Fragment", as.factor)%>%
  filter(Location == "Lab")%>%
  mutate_at("Treatment", as.factor)%>%
  mutate_at("Symbiont", as.factor)%>%
  filter(Symbiont == "Cladocopium")
summary(model.data)

#transform data
model.data$tdata = transformTukey(model.data$Ratio)

#model
symb.model = lmer(tdata ~ Timepoint*Treatment*Phenotype + (1|Colony/Fragment), data=model.data)

#check assumptions
qqPlot(residuals(symb.model)) #normality
leveneTest(residuals(symb.model)~model.data$Timepoint*model.data$Treatment*model.data$Phenotype)

summary(symb.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
rand(symb.model)
anova(symb.model, type = 3) # use this view significant effects and get p-values

#post-hoc testing
emm=emmeans(symb.model, ~Timepoint)
multcomp::cld(emm, Letters=c(LETTERS))

emm=emmeans(symb.model, ~Timepoint*Phenotype)
#extract letters for data visualization
analysis.letters3 <- multcomp::cld(emm, Letters=c(LETTERS))

####visualize C:H in the two phenotypes over time ####
# extract relevant letters from post-hoc tests
plotlabels1 <- analysis.letters3%>%select(Timepoint,Phenotype, .group)%>%
  spread(Timepoint, .group)%>%
  rename(`2019-07-09`=T1)%>%
  rename(`2019-07-14`=T2)%>%
  rename(`2019-08-14`=T4)%>%
  rename(`2019-10-10`=T7)%>%
  gather(Date, .group, -Phenotype)%>%
  mutate_at("Date", as.Date)%>%
  mutate_at(".group", as.factor)%>%
  filter(Phenotype == "B")%>%
  mutate(.group = str_trim(.group))
plotlabels1$Phenotype <- factor(plotlabels1$Phenotype, levels = c("B", "NB"),
                                labels = c("Bleached", "Nonbleached"))

plotlabels2 <- analysis.letters3%>%select(Timepoint,Phenotype, .group)%>%
  spread(Timepoint, .group)%>%
  rename(`2019-07-09`=T1)%>%
  rename(`2019-07-14`=T2)%>%
  rename(`2019-08-14`=T4)%>%
  rename(`2019-10-10`=T7)%>%
  gather(Date, .group, -Phenotype)%>%
  mutate_at("Date", as.Date)%>%
  filter(Phenotype == "NB")%>%
  mutate(.group = str_trim(.group))
plotlabels2$Phenotype <- factor(plotlabels2$Phenotype, levels = c("B", "NB"),
                                labels = c("Bleached", "Nonbleached"))

#extract relevant data, change format to date
symb.plot.data <- model.data%>%dplyr::select(Timepoint, Ratio, Treatment, Phenotype, Fragment)%>%
  spread(Timepoint, Ratio)%>%
  rename(`2019-07-09`=T1)%>%
  rename(`2019-07-14`=T2)%>%
  rename(`2019-08-14`=T4)%>%
  rename(`2019-10-10`=T7)%>%
  gather(Date, ratio, -Phenotype, -Treatment, -Fragment)%>%
  mutate_at("Date", as.Date)

#average by significant effects: timepoint and phenotype
symb.plot<- plyr::ddply(symb.plot.data, c("Date","Phenotype"), summarise,
                        N    = length(ratio[!is.na(ratio)]),
                        mean = mean(ratio, na.rm=TRUE),
                        sd   = sd(ratio, na.rm=TRUE),
                        se   = sd / sqrt(N)
);symb.plot

#rename phenotypes so they are easier to understand
symb.plot$Phenotype <- factor(symb.plot$Phenotype, levels = c("B", "NB"),
                              labels = c("Bleached", "Nonbleached"))

#plot timeseries of C:H ratio
symbPlot<-ggplot(data=symb.plot, aes(x=Date, y=log10(mean), shape = Phenotype)) +
  scale_linetype_manual(values=c("dotted", "solid"))+
  geom_line(aes(group=Phenotype, linetype=Phenotype), colour = "slateblue", position=position_dodge(1.5)) +
  geom_point(aes( shape=Phenotype, group=Phenotype), colour = "slateblue", position=position_dodge(1.5)) +
  geom_errorbar(aes(ymin=log10(mean-se), ymax=log10(mean+se), group=Phenotype), colour = "slateblue", width=0.05, size=0.5, position=position_dodge(1.5), linetype=1)+
  theme_classic()+
  ylab("log(C:H)") +
  annotate("text", x = as.Date("2019-07-30"), y = -3.5, label = "Phenotype p<0.001", size =2, fontface = 3)+
  annotate("text", x = as.Date("2019-07-30"), y = -3.3, label = "Timepoint p<0.001", size =2, fontface = 3)+
  geom_text(data=plotlabels1, aes(label=.group, y= -0.1), show.legend = FALSE, size = 2.5)+
  geom_text(data=plotlabels2, aes(label=.group, y= -3.1), show.legend = FALSE, size = 2.5)+
  scale_x_date(breaks = as.Date(c("2019-07-10", "2019-07-15", "2019-07-30", "2019-08-14",
                                  "2019-08-31", "2019-09-15", "2019-09-29", "2019-10-10")), labels = date_format("%b %d"))+
  xlab("Date")+
  guides(color = guide_legend(override.aes = list(color = "black"))) +
  theme(
    legend.position = "right",
    legend.title = element_text(size =7, face = "bold"),
    legend.text = element_text(size = 7),
    axis.title.y = element_text(face = "bold", size = 7),
    axis.title.x = element_text(face = "bold", size = 7),
    axis.text.x = element_text(angle = 90, size = 7),
    axis.text.y = element_text(size = 7)
  );symbPlot

####model of C and D in mixed community corals over time ####
#select relevant data
NBmodel.data <- final_data%>%dplyr::select(Phenotype, Colony, Treatment, Timepoint, cm_ratio, dm_ratio, Fragment, Location)%>%
  rename(Durusdinium = dm_ratio)%>%
  rename(Cladocopium = cm_ratio)%>%
  gather(Symbiont, Ratio, -Phenotype, -Colony, -Treatment, -Timepoint, -Fragment, -Location)%>%
  mutate_at("Location", as.factor)%>%
  mutate_at("Fragment", as.factor)%>%
  filter(Location == "Lab")%>%
  mutate_at("Treatment", as.factor)%>%
  mutate_at("Symbiont", as.factor)%>%
  filter(Phenotype == "NB")

#transform data
NBmodel.data$tdata = transformTukey(NBmodel.data$Ratio)

#model
NB.model = lmer(tdata ~ Timepoint*Treatment*Symbiont + (1|Colony/Fragment), data=NBmodel.data)

#check assumptions
qqPlot(residuals(NB.model)) #normality
leveneTest(residuals(NB.model)~NBmodel.data$Timepoint*NBmodel.data$Treatment*NBmodel.data$Symbiont)

summary(NB.model) #view coefficients and se, do not use the p-values from this output. In this output look for the amount of deviance accounted for by the random effects
rand(NB.model)
anova(NB.model, type = 3) # use this view significant effects and get p-values

#post-hoc testing
emm=emmeans(NB.model, ~Timepoint)
multcomp::cld(emm, Letters=c(LETTERS))

emm=emmeans(NB.model, ~Symbiont*Timepoint)
#extract letters for data visualization
analysis.letters4 <- multcomp::cld(emm, Letters=c(LETTERS))

####visualize NB symbiont data ####
#extract grouping letters from post-hoc tests for visualization
plotlabelsNB1 <- analysis.letters4%>%select(Timepoint,Symbiont, .group)%>%
  spread(Timepoint, .group)%>%
  rename(`2019-07-09`=T1)%>%
  rename(`2019-07-14`=T2)%>%
  rename(`2019-08-14`=T4)%>%
  rename(`2019-10-10`=T7)%>%
  gather(Date, .group, -Symbiont)%>%
  mutate_at("Date", as.Date)%>%
  mutate(.group = str_trim(.group))%>%
  filter(Symbiont == "Durusdinium")

plotlabelsNB2 <- analysis.letters4%>%select(Timepoint,Symbiont, .group)%>%
  spread(Timepoint, .group)%>%
  rename(`2019-07-09`=T1)%>%
  rename(`2019-07-14`=T2)%>%
  rename(`2019-08-14`=T4)%>%
  rename(`2019-10-10`=T7)%>%
  gather(Date, .group, -Symbiont)%>%
  mutate_at("Date", as.Date)%>%
  mutate(.group = str_trim(.group))%>%
  filter(Symbiont == "Cladocopium")

#extract relevant data, change format to date
NB.plot.data <- NBmodel.data%>%dplyr::select(Timepoint, Ratio, Treatment, Fragment, Symbiont)%>%
  spread(Timepoint, Ratio)%>%
  rename(`2019-07-09`=T1)%>%
  rename(`2019-07-14`=T2)%>%
  rename(`2019-08-14`=T4)%>%
  rename(`2019-10-10`=T7)%>%
  gather(Date, ratio, -Treatment, -Fragment, -Symbiont)%>%
  mutate_at("Date", as.Date)

#average by significant effects: timepoint and symbiont
NB.plot<- plyr::ddply(NB.plot.data, c("Date","Symbiont"), summarise,
                        N    = length(ratio[!is.na(ratio)]),
                        mean = mean(ratio, na.rm=TRUE),
                        sd   = sd(ratio, na.rm=TRUE),
                        se   = sd / sqrt(N)
);NB.plot

#plot timeseries of C:H ratio
NBPlot<-ggplot(data=NB.plot, aes(x=Date, y=log10(mean), colour = Symbiont)) +
  scale_color_manual(values=c("slateblue", "tomato3"))+
  geom_line(aes(group=Symbiont, colour=Symbiont), position=position_dodge(1.5)) +
  geom_point(aes(colour=Symbiont, group=Symbiont), shape = 17, position=position_dodge(1.5),show.legend = FALSE) +
  geom_errorbar(aes(ymin=log10(mean-se), ymax=log10(mean+se), group=Symbiont), width=0.05, size=0.5, position=position_dodge(1.5), linetype=1)+
  theme_classic()+
  ylab("log(S:H)") +
  annotate("text", x = as.Date("2019-07-30"), y = -3.5, label = "Symbiont p<0.001", size =2, fontface = 3)+
  annotate("text", x = as.Date("2019-07-30"), y = -3.3, label = "Timepoint p<0.001", size =2, fontface = 3)+
  annotate("text", x = as.Date("2019-07-30"), y = -3.7, label = "Timepoint * Symbiont p<0.001", size =2, fontface = 3)+
  geom_text(data=plotlabelsNB1, aes(label=.group, y= -0.2), show.legend = FALSE, size = 2.5)+
  geom_text(data=plotlabelsNB2, aes(label=.group, y= -3.1), show.legend = FALSE, size = 2.5)+
  scale_x_date(breaks = as.Date(c("2019-07-10", "2019-07-15", "2019-07-30", "2019-08-14",
                                  "2019-08-31", "2019-09-15", "2019-09-29", "2019-10-10")), labels = date_format("%b %d"))+
  xlab("Date")+
  theme(
    legend.position = "right",
    legend.title = element_text(size =7, face = "bold"),
    legend.text = element_text(size = 7),
    axis.title.y = element_text(face = "bold", size = 7),
    axis.title.x = element_text(face = "bold", size = 7),
    axis.text.x = element_text(angle = 90, size = 7),
    axis.text.y = element_text(size = 7)
  );NBPlot

#compile into figure with log-log plots and Cladocopium timeseries
symbionts <-symbPlot / NBPlot
symbionts + plot_annotation(tag_levels = 'a') + plot_layout(guides = "collect")
############


############PHOTOSYNTHETIC EFFICIENCY############
####get relevant PAM and qPCR data ####
#PAM data, find mean fvfm for each colony/treatment at each timepoint
fvfm_info <- bleaching%>%select(Fragment, Colony, Location, Phenotype, Treatment,T1, T2, T4, T7)%>%
  gather(Timepoint, fvfm, -Fragment, -Colony, -Treatment, -Location, -Phenotype)%>%
  group_by(Colony, Timepoint, Treatment, Phenotype)%>%mutate(col_fvfm = mean(fvfm, na.rm = TRUE))%>%
  select(Colony, Timepoint, col_fvfm, Treatment, Phenotype)%>%unique()

#proportion D at the first timepoint, join with relevant PAM data
d_ambient <- final_data%>%select(Fragment, Colony,  Phenotype, Treatment, Timepoint, prop_d)%>%
  mutate_at("Fragment", as.character)%>%
  mutate_at("Fragment", as.numeric)%>%
  mutate(Location = case_when((Fragment<1000)~"Field", (Fragment>999)~"Lab"))%>%
  mutate_at("Location", as.factor)%>%
  mutate_at("Fragment", as.factor)%>%
  select(Colony, Timepoint, prop_d, Treatment, Phenotype)%>%
  full_join(.,fvfm_info, by = c("Colony" = "Colony", "Timepoint"="Timepoint", "Treatment"="Treatment",
                                "Phenotype"= "Phenotype"))%>%
  filter(Timepoint == "T1")%>%
  filter(Phenotype == "NB")

#proprotion D at the final timepoint, including field samples, join with relevant PAM data
d_stressed <- final_data%>%select(Fragment, Colony,  Phenotype, Treatment, Timepoint, prop_d)%>%
  mutate_at("Fragment", as.character)%>%
  mutate_at("Fragment", as.numeric)%>%
  mutate(Location = case_when((Fragment<1000)~"Field", (Fragment>999)~"Lab"))%>%
  mutate_at("Location", as.factor)%>%
  mutate_at("Fragment", as.factor)%>%
  select(Colony, Timepoint, prop_d, Treatment, Phenotype)%>%
  full_join(.,fvfm_info, by = c("Colony" = "Colony", "Timepoint"="Timepoint", "Treatment"="Treatment",
                                "Phenotype"= "Phenotype"))%>%
  filter(Timepoint == "T7")%>%
  filter(Phenotype == "NB")

####run linear regressions at the two timepoints ####
# thermal stress timepoint
regression_stress <- lm(col_fvfm~prop_d, data = d_stressed) 

summary(regression_stress)
coef(regression_stress)

#reorder and relabel factors
d_stressed$Treatment <- factor(d_stressed$Treatment,
                               levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))
d_stressed$Phenotype <- factor(d_stressed$Phenotype, levels = c("B", "NB"), 
                                      labels = c("Bleached", "Nonbleached"))

# ambient timepoint
regression_ambient <- lm(col_fvfm~prop_d, data = d_ambient) 

summary(regression_ambient)
coef(regression_ambient)

d_ambient$Treatment <- factor(d_ambient$Treatment,
                              levels = c("Control", "Constant High", "Pulse High", "Pulse", "Pulse Increase"))
d_ambient$Phenotype <- factor(d_ambient$Phenotype, levels = c("B", "NB"), 
                               labels = c("Bleached", "Nonbleached"))

#### visualize the two timepoints ####
stress.plot<-ggplot(d_stressed, group = Treatment)+
  scale_color_manual(values =c('#0d0887ff', '#7d03a8ff', '#cb4679ff', '#f89441ff', '#cccc00'))+
  geom_point(aes(x=prop_d, y = col_fvfm, colour = Treatment),shape = 17)+
  geom_abline(aes(intercept =0.640265284, slope = -0.008521945), colour = "grey60")+
  annotate("text", x = 0.25, y = 0.69, label = "slope = -0.0085", size =2, fontface = 3)+
  theme_classic()+
  ylab("Fv/Fm")+
  ylim(0.59,0.78)+
  xlab("Proportion Durusdinium")+
  theme(legend.position = "right",
        legend.title = element_text(face ="bold", size = 7),
        legend.text = element_text(size =7),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.text.y=element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(face = "bold", size = 7));stress.plot

ambient.plot<-ggplot(d_ambient, group = Treatment)+
  geom_point(aes(x=prop_d, y = col_fvfm),shape = 17)+
  geom_abline(aes(intercept = 0.71038205, slope = 0.01479874), colour = "grey60")+
  annotate("text", x = 0.25, y = 0.69, label = "slope = 0.015", size =2, fontface = 3)+
  theme_classic()+
  ylab("Fv/Fm")+
  ylim(0.59,0.78)+
  xlim(0,1)+
  xlab("Proportion Durusdinium")+
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold", size = 7),
        axis.text.y=element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(face = "bold", size = 7));ambient.plot

# compile into one figure
photo.efficiency <- ambient.plot|stress.plot
photo.efficiency + plot_annotation(tag_levels = 'a')

