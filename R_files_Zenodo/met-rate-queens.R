# ---- GOAL: METABOLIC RATE DIFFERENCES OF GROUPS OF QUEENS 
# ---- instal and load packages ----
library("lme4") #LM, LMM, GLM, GLMM
library("dplyr")
library("car") #anova for LM and LMM
library("multcomp")
# Rmisc loads 'plyr' which masks the group_by function from 
library("Rmisc")
library("ggplot2") #for all graphs

# ---- working directories shortcut ----
raw.dir <- ("G:/My Drive/R/met-rate")
out.sheet <- ("G:/My Drive/R/met-rate/data-output")
out.graph <- ("G:/My Drive/R/met-rate/graphs/queens-only")

# ---- Loading raw per file data and assigning features----
setwd(raw.dir)
list.files() #visualize files in the folder

# loading compiled metabolic rate CSV file and adding variables
# whole colony measurement
df.mr <- read.csv("mr-queens.csv")
df.mr <- df.mr %>%
  mutate(co2 = (co2.mean*air.flow*60)/1000000, #1ppm = 10^-6 g/ml, airflow = ml/min, 60 minute = 1 hour. CO2 flow in grams per hour
         co2.emission.rate = co2/mass, # ml/g*hr. CO2 emission rate per gram of queens
         o2.consmp = co2.emission.rate/rq, # ml/g*hr. #Oxygen consumption rate per gm of queens Species specific respiratory quotient to convert CO2 emission rate to O2 consumption rate.
         oxy.joule = (16 + (5.164*rq)), #Energy rate. mathematical formula to calculate joule heat energy produced per unit of O2 consumed
         met.rate.j = (o2.consmp*oxy.joule)/3600, #Joules/sec at measured temperature per gram of queens
         met.rate.w = met.rate.j*1000000,#microwatts at measured temperature from whole queens
         met.rate.25 = met.rate.w/(2^((temp-25)/10))) #micro watts per gram of queens at 25'C

# data for publication - removing extra columns
df.mr <- data.frame(subset(df.mr, select = -c(Date,
                                              first.baseline:last.baseline)))

# ---- Adding queen age  ----
df.mr$age.bin <- NA #blank column 

#if condition meets, then add age group, if not then move to other condition
df.mr <- df.mr %>%
  mutate(age.bin = ifelse(queen.age >= 40 & queen.age <= 45, "40-45", 
                          ifelse(queen.age >= 46 & queen.age <= 50, "46-50",
                                 ifelse(queen.age >= 51 & queen.age <= 55, "51-55",
                                        ifelse(queen.age >= 56 & queen.age <= 60, "56-60", 
                                               ifelse(queen.age >= 61 & queen.age <= 65, "61-65", "too old")))))) #if change in male < 30  then assign 1, if not then assign 2

# saving this file in output directory
setwd(out.sheet)
write.csv(df.mr, "mr2019-queens-compiled.csv")

df.mr <- read.csv("mr2019-queens-compiled.csv")

# ---- Sample size -----
# sample size analyzed per infection group and queen age.
mr.sample <- df.mr %>% 
  group_by(wolbachia, queen.age) %>% 
  summarize(col.sample.size = n())

mr.sample1 <- df.mr %>% 
  group_by(wolbachia, age.bin) %>% 
  summarize(col.sample.size = n())

# sample size analyzed per infection group
mr.sample.wol <- df.mr %>% 
  group_by(wolbachia) %>% 
  summarize(col.sample.size = n())

# ---- Testing different aspects of data ####################

#check disrtibution of data
hist(df.mr$met.rate.25) #looks normal

#Shapiro Wilks Test for sampling from normal distribution, p < 0.05: not drawn from normal distribution (use GLM
shapiro.test(df.mr$met.rate.25) #p > 0.05. 

#Batlett's test for unequal variance. p-value < 0.05 : variance can not be assumed equal
bartlett.test(df.mr$met.rate.25 ~ wolbachia, data=df.mr) #p > 0.05

#F-test for equal variance, p < 0.05: variance not equal
var.test(df.mr$met.rate.25 ~ wolbachia, data =df.mr ,alternative = "two.sided") ## p > 0.01

# --- LMER to identify effects #########
# metabolic rate with individual queen age
lm.mr.qn <- lm (met.rate.25~ wolbachia+ as.factor(queen.age) + as.factor(queen.age)*wolbachia,
                     data = df.mr)

drop1(lm.mr.qn, test = "Chisq")
car::Anova(lm.mr.qn, ddf="Satterthwaite")
summary(lm.mr.qn)
res <- residuals(lm.mr.qn)
hist(res)
qqPlot(res)
shapiro.test(res) #p<0.05

# metabolic rate with queen age binned
lm.mr.qn.bn <- lm (met.rate.25~ wolbachia+ age.bin + age.bin:wolbachia,
                data = df.mr)

drop1(lm.mr.qn.bn, test = "Chisq")
car::Anova(lm.mr.qn.bn, ddf="Satterthwaite")
summary(lm.mr.qn.bn)
res <- residuals(lm.mr.qn.bn)
hist(res)
qqPlot(res)
shapiro.test(res) #p>0.05

# mass of queens with individual queen age
lm.ms.qn <- lm (mass~ wolbachia+ as.factor(queen.age) + as.factor(queen.age):wolbachia,
                data = df.mr)

drop1(lm.ms.qn, test = "Chisq")
car::Anova(lm.ms.qn, ddf="Satterthwaite")
summary(lm.ms.qn)
res <- residuals(lm.ms.qn)
hist(res)
qqPlot(res)
shapiro.test(res) #p<0.05

# mass of queens with queen age binned
lm.ms.qn.bn <- lm (mass~ wolbachia+ age.bin + age.bin:wolbachia,
                data = df.mr)

drop1(lm.ms.qn.bn, test = "Chisq")
car::Anova(lm.ms.qn.bn, ddf="Satterthwaite")
summary(lm.ms.qn.bn)
res <- residuals(lm.ms.qn.bn)
hist(res)
qqPlot(res)
shapiro.test(res) #p<0.05



# ---- PLOT : Sample value with mean ---------
my.color <-c ("chocolate2","darkorchid4")

# Step 1: universal variable and plot text
df.mr$x <- df.mr$wolbachia
df.mr$y <- df.mr$met.rate.25
title <- "Meatbolic Rate of Queens \\nat 25'C"
x.lab <- "Wolbachia infectio"
y.lab <- "metabolic rate \\n(microwatts/gram)"

# Step 2: statistical summary
mr.mean <- summarySE(df.mr, measurevar="y",
                      groupvars=c("wolbachia"), na.rm = T)

# Step 3: re-order variables
df.mr$wolbachia<- factor(df.mr$wolbachia, 
                         levels = c("uninfected", "infected"))


mr.mean$wolbachia<- factor(mr.mean$wolbachia, 
                         levels = c("uninfected", "infected"))


# Step 4: Plot = sample values overlaid with mean and 95% CI. Texts were later edited in Inkscape
ggplot(df.mr, aes(x= x,
                  y = y,
                  group = queens.id,
                  color = wolbachia)) + 
  scale_colour_manual(values = my.color)+
  scale_fill_manual(values = my.color)+
  geom_dotplot(stackdir='center',
               position=position_dodge(0.8)) +
  ggtitle(title)+
  scale_x_discrete(limits = seq(0, 20, 4))+
  labs(x = x.lab,
       y = y.lab) +
  theme(strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(family = "Arial", colour = "black", face = "italic",size = "12"),
        plot.title = element_text(hjust = 0.5,
                                  family = "Arial", colour = "black", face = "bold", size = "14"),
        panel.background=element_rect(fill = NA),
        panel.spacing.x=unit(0.1, "in"),
        plot.background=element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        #legend.position = c(0.08,0.80),
        legend.text = element_text(family = "Arial", colour = "black", size = "10"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
        axis.text = element_text(family = "Arial", colour = "black", size = "12"),
        axis.ticks = element_line(colour = "black", size = 1))

#saving image in desired format
ggsave(file = "total-queens-wrap.png",
       path = gr.line,
       dpi = 300, 
       width = 10,
       height = 6, 
       units = "in")


# ---- PLOT : Box plot with raw values ---------
my.color <-c ("chocolate2","darkorchid4")

df.mr$m <- df.mr$mass
df.mr$wolbachia<- factor(df.mr$wolbachia, 
                             levels = c("uninfected", "infected"))

ggplot(df.mr, aes(x=wolbachia,
                      y = m,
                      fill = wolbachia)) + 
  geom_boxplot(notch = F,
               position= position_dodge(0.8))+
  facet_wrap(~age.bin)+
  scale_fill_manual(values = my.color)+
  labs(title = "Mass \\nof Groups of 15 Queens",
       y = "microwatts per gram")+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial",
                                  colour = "black",
                                  face = "italic",
                                  size = "12"),
        panel.background=element_rect(fill = NA),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        #legend.position = c(0.08,0.80),
        legend.text = element_text(family = "Arial",
                                   colour = "black",
                                   size = "10"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_text(family = "Arial",
                                  colour ="black",
                                  face = "bold", size = "11"),
        axis.text = element_text(family = "Arial",
                                 colour = "black", size = "12"),
        axis.title.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.background=element_rect(fill = NA),
        plot.title = element_text(family = "Arial", colour = "black",
                                  face = "bold", size = "12", 
                                  margin = margin(t='1', r='0.5', b='0.5', l='1', unit = 'cm')))

#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(file = "mr-queens-mass-age.png",
       path = out.graph,
       dpi = 300, 
       width = 3.5,
       height = 3.5, 
       units = "in")

