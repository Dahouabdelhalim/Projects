################################################# Experiment 1 - PATHOGEN DRIVERS ################################################################

# import data
library(readxl)
Distance <- read_excel("DistanceInds.xlsx")

Distance$Status <- as.factor(Distance$Status)
Distance$Name <- as.factor(Distance$Name)

# divide flies per treatment
infected <- subset(Distance, Distance$Status == "Infected")
susceptible <- subset(Distance, Distance$Status == "Susceptible")
InfSusc <- subset(Distance, Distance$Status == "Inf-Susc")


################################# INFECTED ##################################
library(lme4)
library(carData)
library(car)

modeli <- lmer (dist_mean ~ Name*Dose*Time + (1|Date), data = infected)
summary(modeli)

# ANOVA
Anova(modeli, test.statistic="F")
# pathogen and pathogen:dose are significant
summary(Anova(modeli, type =2))

# post hoc
library(emmeans)      
lsmeans(modeli, pairwise~Name*Dose, adjust="Tukey")

library(DescTools)
D1i <- subset(infected, infected$Dose != "D2")
summary(D1i)
D2i <- subset(infected, infected$Dose != "D1")
summary(D2i)

DunnettTest(x= D1i$dist_mean, g= D1i$Name)
DunnettTest(x= D2i$dist_mean, g= D2i$Name)
# infected flies stay closer together when compared to controls


############################### SUSCEPTIBLE #################################

library(lme4)
library(lmerTest)
library(car)

models <- lmer(dist_mean ~ Name*Dose*Time + (1|Date), data = susceptible) 
summary(models)

# ANOVA
Anova(models, test.statistic="F")
## time and pathogen:dose are significant
summary(Anova(models, type =2))

# post hoc
library(emmeans)  
lsmeans(models, pairwise~Name*Dose, adjust="Tukey")

library(DescTools)
D1s <- subset(susceptible, susceptible$Dose != "D2")
summary(D1s)
D2s <- subset(susceptible, susceptible$Dose != "D1")
summary(D2s)

DunnettTest(x= D1s$dist_mean, g= D1s$Name)
DunnettTest(x= D2s$dist_mean, g= D2s$Name)


############################## INF-SUSC ####################################

library(lme4)
library(lmerTest)
library(car)

modelis <- lmer(dist_mean ~ Name*Dose*Time + (1|Date), data = InfSusc)
summary(modelis)

#ANOVA
Anova(modelis, test.statistic="F")
# no significant effects 
summary(Anova(modelis, type =2))


# post hoc   ### NOT NEEDED
library(DescTools)
D1is <- subset(InfSusc, InfSusc$Dose != "D2")
summary(D1is)
D2is <- subset(InfSusc, InfSusc$Dose != "D1")
summary(D2is)

DunnettTest(x= D1is$dist_mean, g= D1is$Name)
DunnettTest(x= D2is$dist_mean, g= D2is$Name)


############################## GRAPHS ####################################
library(dplyr)
library(ggplot2)

## D1 = 0.01 (low dose) and D2 = 0.1 (high dose)

########### infected

# Overall rate of aggregation - boxplot
summaryboxi <- infected %>%
  group_by(Dose, Name, Replica) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))

summaryboxi$Dose <- factor(summaryboxi$Dose, 
                           levels = c("D1", "D2", "none"))
summaryboxi$Dose <- relevel(summaryboxi$Dose, "D1")

pi <- ggplot(summaryboxi, aes(y=mean_distance, x=Name, fill = Dose)) + 
  geom_point(aes(fill = Dose), size = 1, shape = 21, alpha=0.2, position = position_jitterdodge()) + 
  geom_boxplot(alpha=0.6,key_glyph = "rect") +
  xlab("Pathogen") + ylab("Mean pairwise distance (mm)") + 
  theme_classic() + scale_x_discrete(labels=c("Control_color" = "Control", "Efae" = "E. faecalis",
                                              "Pento" = "P. entomophila", "Pret" = "P. rettgeri", "Smar" = "S. marcescens"))+
  scale_fill_manual(values = c("#FC766AFF", "#184A45FF", "#B0B8B4FF"), name = "Dose", labels = c("O.D. 0.01", "O.D. 0.1", "None"))+ 
  ggtitle("Infected Group")+ theme(plot.title = element_text(hjust = 0.5, size = 14))
plot(pi)


# Social aggregation across time
pathogens.labs <- c("E. faecalis", "P. entomophila", "P. rettgeri", "S. marcescens")
names(pathogens.labs) <- c("Efae", "Pento", "Pret", "Smar")

InfD1 <- Distance %>%
  filter(Status == "Infected", Dose == "D1") %>%
  group_by(Time, Name) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))

InfD2 <- Distance %>%
  filter(Status == "Infected", Dose == "D2") %>%
  group_by(Time, Name) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))

Control.inf <- Distance %>%
  filter(Status == "Infected", Dose == "none") %>%
  group_by(Time, Name) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))


pinf <- ggplot(data = InfD1, aes(x = Time, y = mean_distance)) +
  geom_pointrange(data=transform(Control.inf,Name=NULL), #this line allows "control" to be plotted across the graphs
                  aes(ymin = mean_distance -se_distance, 
                      ymax = mean_distance + se_distance), size= 0.3, color="#B0B8B4FF") +
  geom_pointrange(data = InfD1, aes(ymin = mean_distance - se_distance,
                                    ymax = mean_distance + se_distance), size=0.3, color="#FC766AFF",alpha=0.6) +
  geom_pointrange(data = InfD2, aes(ymin = mean_distance - se_distance,
                                    ymax = mean_distance + se_distance), size=0.3, color="#184A45FF",alpha=0.6) + 
  geom_line(data=transform(Control.inf,Name=NULL), size= 0.7, color="#B0B8B4FF") +
  geom_line(data=InfD1, size = 0.7,  color="#FC766AFF",alpha=0.6) +
  geom_line(data=InfD2, size= 0.7, color="#184A45FF",alpha=0.6) + 
  xlab("Time point") + ylab("Mean pairwise distance (mm)") +
  facet_wrap("Name", labeller = labeller(Name = pathogens.labs)) +
  theme_bw()+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=c(2,4,6,8))+
  scale_y_continuous(breaks=c(18,21,24))+
  ggtitle("Infected Group") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))
plot(pinf)

########### susceptible

# Overall rate of aggregation - boxplot
summaryboxs <- susceptible %>%
  group_by(Dose, Name, Replica) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))

summaryboxs$Dose <- factor(summaryboxs$Dose, 
                           levels = c("D1", "D2", "none"))
summaryboxs$Dose <- relevel(summaryboxs$Dose, "D1")

ps <- ggplot(summaryboxs, aes(y=mean_distance, x=Name, fill = Dose)) + 
  geom_point(aes(fill = Dose), size = 1, shape = 21, alpha=0.2, position = position_jitterdodge()) +
  geom_boxplot(alpha=0.6,key_glyph = "rect") +                                                                                 
  xlab("Pathogen") + ylab("Mean pairwise distance (mm)") + 
  theme_classic() + scale_x_discrete(labels=c("Control_color" = "Control", "Efae" = "E. faecalis",
                                              "Pento" = "P. entomophila", "Pret" = "P. rettgeri", "Smar" = "S. marcescens")) +
  scale_fill_manual(values = c("#FC766AFF", "#184A45FF", "#B0B8B4FF"), name = "Dose", labels = c("O.D. 0.01", "O.D. 0.1", "None")) +
  ggtitle("Susceptible Group")+theme(plot.title = element_text(hjust = 0.5, size = 14))
plot(ps)


# Social aggregation across time
#pathogens.labs <- c("E. faecalis", "P. entomophila", "P. rettgeri", "S. marcescens")
#names(pathogens.labs) <- c("Efae", "Pento", "Pret", "Smar")

SuscD1 <- Distance %>%
  filter(Status == "Susceptible", Dose == "D1") %>%
  group_by(Time, Name) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))

SuscD2 <- Distance %>%
  filter(Status == "Susceptible", Dose == "D2") %>%
  group_by(Time, Name) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))

Control.susc <- Distance %>%
  filter(Status == "Susceptible", Dose == "none") %>%
  group_by(Time, Name) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))


psusc <- ggplot(data = SuscD1, aes(x = Time, y = mean_distance)) +
  geom_pointrange(data=transform(Control.susc,Name=NULL), #this line allows "control" to be plotted across the graphs
                  aes(ymin = mean_distance -se_distance, 
                      ymax = mean_distance + se_distance), size= 0.3, color="#B0B8B4FF") +
  geom_pointrange(data = SuscD1, aes(ymin = mean_distance - se_distance,
                                     ymax = mean_distance + se_distance), size=0.3, color="#FC766AFF",alpha=0.6) +
  geom_pointrange(data = SuscD2, aes(ymin = mean_distance - se_distance,
                                     ymax = mean_distance + se_distance), size=0.3, color="#184A45FF",alpha=0.6) + 
  geom_line(data = transform(Control.susc,Name=NULL), size= 0.7, color="#B0B8B4FF") +
  geom_line(data = SuscD1, size = 0.7,  color="#FC766AFF",alpha=0.6) +
  geom_line(data = SuscD2, size= 0.7, color="#184A45FF",alpha=0.6) + 
  xlab("Time point") + ylab("Mean pairwise distance (mm)") +
  facet_wrap("Name", labeller = labeller(Name = pathogens.labs)) +
  theme_bw()+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=c(2,4,6,8))+
  scale_y_continuous(breaks=c(22,24,26,28))+
  ggtitle("Susceptible Group") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))
plot(psusc)


########### infected-susceptible

# Overall rate of aggregation - boxplot
summaryboxis <- InfSusc %>%
  group_by(Dose, Name, Replica) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))

summaryboxis$Dose <- factor(summaryboxis$Dose, 
                           levels = c("D1", "D2", "none"))
summaryboxis$Dose <- relevel(summaryboxis$Dose, "D1")

pis <- ggplot(summaryboxis, aes(y=mean_distance, x=Name, fill = Dose)) + 
  geom_point(aes(fill = Dose), size = 1, shape = 21, alpha=0.2, position = position_jitterdodge()) +
  geom_boxplot(alpha=0.6,key_glyph = "rect") +                                                                                    
  xlab("Pathogen") + ylab("Mean pairwise distance (mm)") + 
  theme_classic() + scale_x_discrete(labels=c("Control_color" = "Control", "Efae" = "E. faecalis",
                                              "Pento" = "P. entomophila", "Pret" = "P. rettgeri", "Smar" = "S. marcescens"))+
  scale_fill_manual(values = c("#FC766AFF", "#184A45FF", "#B0B8B4FF"), name = "Dose", labels = c("O.D. 0.01", "O.D. 0.1", "None"))+
  ggtitle("Infected-Susceptible Group")+theme(plot.title = element_text(hjust = 0.5, size = 14))

plot(pis)


# Social aggregation across time
#pathogens.labs <- c("E. faecalis", "P. entomophila", "P. rettgeri", "S. marcescens")
#names(pathogens.labs) <- c("Efae", "Pento", "Pret", "Smar")

ISD1 <- Distance %>%
  filter(Status == "Inf-Susc", Dose == "D1") %>%
  group_by(Time, Name) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))

ISD2 <- Distance %>%
  filter(Status == "Inf-Susc", Dose == "D2") %>%
  group_by(Time, Name) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))

Control.IS <- Distance %>%
  filter(Status == "Inf-Susc", Dose == "none") %>%
  group_by(Time, Name) %>%
  summarise(mean_distance=mean(dist_mean), se_distance = (sd(dist_mean)/sqrt(length(dist_mean))))


pinfsusc <- ggplot(data = ISD1, aes(x = Time, y = mean_distance)) +
  geom_pointrange(data=transform(Control.IS,Name=NULL), #this line allows "control" to be plotted across the graphs
                  aes(ymin = mean_distance -se_distance, 
                      ymax = mean_distance + se_distance), size= 0.3, color="#B0B8B4FF") +
  geom_pointrange(data = ISD1, aes(ymin = mean_distance - se_distance,
                                     ymax = mean_distance + se_distance), size=0.3, color="#FC766AFF",alpha=0.6) +
  geom_pointrange(data = ISD2, aes(ymin = mean_distance - se_distance,
                                     ymax = mean_distance + se_distance), size=0.3, color="#184A45FF",alpha=0.6) + 
  geom_line(data = transform(Control.IS,Name=NULL), size= 0.7, color="#B0B8B4FF") +
  geom_line(data = ISD1, size = 0.7,  color="#FC766AFF",alpha=0.6) +
  geom_line(data = ISD2, size= 0.7, color="#184A45FF",alpha=0.6) + 
  xlab("Time point") + ylab("Mean pairwise distance (mm)") +
  facet_wrap("Name", labeller = labeller(Name = pathogens.labs)) +
  theme_bw()+
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=c(2,4,6,8))+
  scale_y_continuous(breaks=c(22,24,26,28))+
  ggtitle("Infected-Susceptible Group") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))
plot(pinfsusc)


library(cowplot)
plot_grid(pinf, pi, psusc, ps, pinfsusc, pis, labels = c('A', 'B', 'C', 'D', 'E', 'F'), label_size = 12, ncol = 2)


######################## Checking model assumptions ##############################

# check model
install.packages("performance")
install.packages("remotes")
remotes::install_github("easystats/performance")
library(remotes)
library(performance)
install.packages("see")  
library(see)
modeli <- lmer (dist_mean ~ Name+Dose+Time + (1|Date), data = infected)
check_model(modeli)
models <- lmer(dist_mean ~ Name+Dose+Time + (1|Date), data = susceptible)
check_model(models)
modelis <- lmer(dist_mean ~ Name+Dose+Time + (1|Date), data = InfSusc)
check_model(modelis)


############################################################# HOST DRIVERS ############################################################

# import data and create logNND column
NND <- read.csv("NND-Boyle.csv")
NND <- cbind(NND, log10(NND$Median.NND))
names(NND)[names(NND)=="log10(NND$Median.NND)"] <- "logNND"
names(NND)[names(NND)=="Ã¯..line"] <- "line"
NND$line <- as.factor(NND$line)

library(lme4)
library(car)
mymodel <- lmer(logNND ~ line*ssex*Infection.status + (1|day.replicate), data = NND)
summary(mymodel)

#ANOVA
Anova(mymodel, test.statistic="F", type =2)
# line, sex and sex*status are significant

# post hoc
library(emmeans)
#lsmeans(mymodel, pairwise~ssex, adjust="Tukey")
### on average females stay closer together (-0.0282, p=0.0427) retransform: 1.06708742mm
lsmeans(mymodel, pairwise~ssex*Infection.status, adjust="Tukey")
### infected females aggregated more closely than infected males; female inf:male inf= -0.05936 (df 418, p=0.0135)
### re transform (10^0.05936 = 1.1464628859 mm)


######################## Checking model assumptions ##############################
library(remotes)
library(performance)
library(see)
mymodel <- lmer(logNND ~ line + ssex + Infection.status + (1|day.replicate), data = NND)
check_model(mymodel)


################################## GRAPHS #######################################

library(RColorBrewer)
library(ggplot2)
farb <- c("#CE2929", "#5675D6")

## OPTION 1
sex.labs <- c("Female", "Male")
names(sex.labs) <- c("female", "male")

summaryNND <- NND %>%
  group_by(ssex, Infection.status, line) %>%
  summarise(mean_distance=mean(Median.NND), se_distance = (sd(Median.NND)/sqrt(length(Median.NND))))
ggplot(summaryNND, aes(y=mean_distance, x=line,  ymin=mean_distance-se_distance, 
                       ymax=mean_distance+se_distance,fill = ssex, alpha = Infection.status)) + geom_bar(stat="identity",
                                                                                                         position=position_dodge()) +
  xlab("Line") + ylab("Median NND (mm)") + facet_wrap("ssex", labeller = labeller(ssex = sex.labs)) + 
  labs(fill = "Sex", alpha="Infection Status") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  scale_fill_manual(values= farb)+
  scale_alpha_manual(values=c(0.5, 1))+
  geom_errorbar(data = summaryNND, position=position_dodge(width=0.9)) +
  guides(fill = FALSE)

# not infected: no differences in sex; infected: females get closer, males further apart
#reorder
NND$line <- factor(NND$line, levels=c("427", "324", "358", "712", "208", "21", "304", "375", "852", "28"))
#interaction plot
NND$ssex<- as.factor(NND$ssex)


## OPTION 2
#status.labs <- c("Control", "Infected")
#names(status.labs) <- c("control", "infected")
library("viridis")  

Fig2 <- ggplot(NND, aes(x=Median.NND, y= line, fill= ssex)) + 
  geom_boxplot (alpha = 0.6) +
  xlab("Median NND (mm)") + ylab("Genetic background") +
  facet_wrap(~ Infection.status) +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
                  # element.text = element_text(size = 14))+
  scale_fill_viridis(discrete = TRUE, option = "D")+
  #scale_fill_manual(values= c("#FC766AFF", "#184A45FF"))+
  #scale_alpha_manual(values=c(0.5, 1))+
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(title="Sex")) +
  theme(
    strip.background = element_rect(
      color="black", fill="grey")) + 
  theme(strip.text.x = element_text(size = 12))

ggsave("Exp2_Figure2.png")




## 2nd graph

### violin

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


Fig3 <- ggplot(NND, aes(x = Infection.status,  y = Median.NND, colour = ssex, group = ssex)) + 
  stat_summary(fun = mean, geom = "point") +  stat_summary(fun = mean, geom = "line", size=1.2)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  scale_color_viridis(discrete = TRUE, option = "D")+
  xlab("Infection Status") + ylab("Median NND (mm)") + labs(color = "Sex")+
  geom_point(alpha=0.3, position = "jitter")+ylim(0,12)  +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 10)) +
  guides(fill = guide_legend(title="Sex")) +
  theme(strip.background = element_rect(color="black", fill="grey", size=0.9, linetype="solid")) + 
  theme(strip.text.x = element_text(size = 14)) 
ggsave("Exp2_Figure3.png")


library(patchwork)

Fig2 - Fig3
