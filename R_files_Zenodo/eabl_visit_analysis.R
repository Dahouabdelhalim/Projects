############################################################
############### RFID Data Analysis 2020 ################
############################################################

#Installing feedR package

install.packages("devtools") # if not already installed
devtools::install_github("animalnexus/feedr")

#Loading necessary libraries

library(lme4) #creates linear mixed models
library(car) #Anova test
library(ggplot2) #create good graphs
library(dplyr) #data management
library(lubridate) #time and date variable management
library(bbmle) #AICc tests
library(lmerTest) #Calculates p-values for LMMs

setwd() #Set working directory


rfid = read.csv("rfid_visits_clean.csv", header = TRUE)
# Main Models - Females and Double broods  ----------------------

m1d = lmer(log_visits ~ treatment*nestling_age + 
             treatment*treat_time +
             jhatch_date +
             brood_size + #no interaction because treatment wouldn't affect it
             (1|animal_id), 
           data = rfid, 
           # family = "poisson"
           REML = FALSE
)
summary(m1d)



# Regression Diagnostics --------------------------------------------------
#Assessing Outliers
library(car)
library(predictmeans)
library(tidyverse)
library(broom)

outlierTest(m3d) # Bonferonni p-value for most extreme obs
qqPlot(m3d, main="QQ Plot") #qq plot for studentized resid 
leveragePlots(m3d) # leverage plots

assump<-function(model) {       #graphs for checking assumptions
r<-residuals(m3d)
ft<-fitted(m3d)
par(mfrow=c(2,2))
library(MASS)
truehist(r,main="Histogram of Residuals",xlab="Residuals")
curve(dnorm(x,mean=mean(r),sd=sd(r)),add=TRUE)
qqnorm(r, ylim=range(r), main="QQNorm Plot",ylab="Quantiles of (ie, ordered) residuals", xlab="Quantiles of normal distribution")
qqline(r,lty=2)
plot(r~ft,main="Residuals vs Fitted",xlab="Fitted (predicted) values",ylab="Residuals");abline(h=0)
qqnorm(ft, ylim=range(ft), main="QQNorm Plot",ylab="Quantiles of fitted values", xlab="Quantiles of normal distribution")
qqline(ft,lty=2)
acf(resid(model))
  }
assump(m3d)

par(mfrow = c(2,2))
plot(m7d)



#####Plotting Visitation Rate within Days (Per Hour)######
#Plotting visitation rates across hours

library(ggplot2)
library(ggthemes)
library(ggsci)
library(sjPlot)
library(sjmisc)

rfidMeans = rfid %>% #Finds Means of data for plotting purposes
  group_by(treatment, nestling_age) %>% #grouping by treatment and days since hatch
  mutate(treatment = dplyr::recode(treatment, 'control' = "Control", #need to have single quote marks on existing values
                                   'noise' = "Noise"),
                            # recode(animal_id = as.factor(animal_id))
         ) %>%
  dplyr::summarize(n = n(), #need to use dplyr:: because other libraries have the summarize function
                   meanVisits = mean(log(visits_per_ind_hour)), 
                   seVisits = (sd(log(visits_per_ind_hour))/sqrt(n))
  )

#Scatterplot for mean visitation rate across nestling age (days since hatch) - grey
ggplot(rfidMeans, aes(x=nestling_age, 
                            y=meanVisits, 
                            color = treatment)) + 
  geom_errorbar(aes(ymin=meanVisits-seVisits, 
                    ymax=meanVisits+seVisits), 
                width=.5) + 
  geom_line(size = 3) +
  geom_point(size=5) +
  labs(colour = "Treatment",
       x = "Nestling Age (Days)",
       y = "Mean Visitation Rate \\n(Log(Visits/hr))") +
  scale_colour_grey() +
  theme_classic(base_size = 24)

library(extrafont)
#font_import() only do this one time - it takes a while
loadfonts(device="win")
windowsFonts(Times=windowsFont("TT Times New Roman"))

ggplot(rfidMeans, aes(x=nestling_age, 
                            y=meanVisits, 
                            color = treatment)) + 
  geom_errorbar(aes(ymin=meanVisits-seVisits, 
                    ymax=meanVisits+seVisits), 
                width=.5) + 
  geom_line(size = 1) +
  geom_point(size=2) +
  labs(colour = "Treatment",
       x = "Nestling Age (Days)",
       y = "Mean Visitation Rate \\n(Log(Visits/hr))") +
  scale_colour_manual(values = c("#E69F00", #orange
                        "#56B4E9"))+ #blue +
  # theme_bw(base_size=12, base_family='TT Times New Roman')+
  # theme(panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank())
  theme_classic(base_size=12, 
                base_family='TT Times New Roman')
setwd("C:/Users/meely/OneDrive - University of Oklahoma/University of Oklahoma/Bridge Lab/EABL Anthro Noise/results")
ggsave("visits_nestling_age.tiff", dpi=300, height=4.5, width=6.5, units="in")
