
########################Reproductive strategies affect telomere dynamics across the life course########################

library(lme4)
library(lsmeans)
setwd()
data <- read.csv("Quail_telomeres_offspring.csv")

head(data)
tail(data)
names(data)
summary(data)
str(data)

########################
# Telomere length day 3
m1 <- lmer(telomere.d3 ~ line+sex+line:sex+ mass.hatching+ (1|family)+(1|gel), data = data)
summary(m1)

# interaction effect
m1a <- lmer(telomere.d3 ~ line+sex+ mass.hatching+ (1|family)+(1|gel), data = data)
anova(m1a, m1)

# Contrasts
lsmeans(m1, pairwise ~ line | sex)

# body mass effect
m1b <- lmer(telomere.d3 ~ line+sex+ line:sex+(1|family)+(1|gel), data = data)
anova(m1b, m1)



##########################
# Telomere length day 30 

m2 <- lmer(telomere.d30 ~ line+sex+  line:sex+ mass.d30 +(1|family)+(1|gel), data = data)
summary(m2)

# interaction effect
m2a <- lmer(telomere.d30 ~ line+sex+  mass.d30 +(1|family)+(1|gel), data = data)
anova(m2a, m2)

# sex effect
m2b <- lmer(telomere.d30 ~ line+  mass.d30 +(1|family)+(1|gel), data = data)
anova(m2b, m2a)

# line effect
m2c <- lmer(telomere.d30 ~ sex+  mass.d30 +(1|family)+(1|gel), data = data)
anova(m2c, m2a)

# body mass effect
m2d <- lmer(telomere.d30 ~ line+sex +(1|family)+(1|gel), data = data)
anova(m2d, m2a)

# family effect
m2e <- lmer(telomere.d30 ~ line+sex+mass.d30+(1|gel), data = data)
anova(m2e, m2a)

# gel effect
m2f <- lmer(telomere.d30 ~ line+sex+mass.d30+ (1|family), data = data)
anova(m2f, m2a)


################################
# Hatching mass


m3 <- lmer(mass.hatching ~ line+sex+line:sex +(1|family), data = data)
summary(m3)

# interaction effect
m3a <- lmer(mass.hatching ~ line+sex +(1|family), data = data)
anova(m3a, m3)

# line effect
m3b <- lmer(mass.hatching ~ sex +(1|family), data = data)
anova(m3b, m3a)

# sex effect
m3c <- lmer(mass.hatching ~ line +(1|family), data = data)
anova(m3c, m3a)

# family effect
m3d <- lm(mass.hatching ~ line+sex, data = data)
anova(m3a, m3d)


################################
# Body mass gain hatching - day 30

data$bodymasschange<-(data$mass.d30-data$mass.hatching)


m4 <- lmer(bodymasschange ~ line+sex+line:sex +(1|family), data = data)
summary(m4)

# interaction effect
m4a <- lmer(bodymasschange ~ line+sex +(1|family), data = data)
anova(m4a, m4)

# line effect
m4b <- lmer(bodymasschange ~ sex +(1|family), data = data)
anova(m4b, m4a)

# sex effect
m4c <- lmer(bodymasschange ~ line +(1|family), data = data)
anova(m4c, m4a)

# family effect
m4d <- lm(bodymasschange ~ line+sex, data = data)
anova(m4a, m4d)



################################
# Change in telomere length day 3 - day 30

data$telomerechange<-(data$telomere.d30-data$telomere.d3)

m5 <- lmer(telomerechange ~ line+sex+  line:sex+ bodymasschange +(1|family)+(1|gel), data = data)
summary(m5)

# interaction effect
m5a <- lmer(telomerechange ~ line+sex+  bodymasschange +(1|family)+(1|gel), data = data)
anova(m5a, m5)

# sex effect
m5b <- lmer(telomerechange ~ line+  bodymasschange +(1|family)+(1|gel), data = data)
anova(m5b, m5a)

# line effect
m5c <- lmer(telomerechange ~ sex+  bodymasschange +(1|family)+(1|gel), data = data)
anova(m5c, m5a)

# body mass change effect
m5d <- lmer(telomerechange ~ line+sex +(1|family)+(1|gel), data = data)
anova(m5d, m5a)

# family effect
m5e <- lmer(telomerechange ~ line+sex+bodymasschange+(1|gel), data = data)
anova(m5e, m5a)

# gel effect
m5f <- lmer(telomerechange ~ line+sex+bodymasschange+ (1|family), data = data)
anova(m5a, m5f)



################################
# Adult telomere length

data1 <- read.csv("Quail_telomeres_adults.csv")

m6 <- lmer(telomere.adult ~ line+sex+line:sex+ body.mass+ (1|gel), data = data1)
summary(m6)

# interaction effect
m6a <- lmer(telomere.adult ~ line+sex+ body.mass+ (1|gel), data = data1)
anova(m6a, m6)

# line effect
m6b <- lmer(telomere.adult ~ sex+ body.mass+ (1|gel), data = data1)
anova(m6b, m6a)

# sex effect
m6c <- lmer(telomere.adult ~ line+ body.mass+ (1|gel), data = data1)
anova(m6c, m6a)

# body mass effect
m6d <- lmer(telomere.adult ~ line+sex+(1|gel), data = data1)
anova(m6d, m6a)

# gel effect
m6f <- lm(telomere.adult ~ line+sex+body.mass, data = data1)
anova(m6a,m6f)



################################
# Age adults

m7 <- lm(age.breeding ~ line+sex+line:sex, data = data1)
summary(m7)

m7a <- lm(age.breeding ~ line+sex, data = data1)
summary(m7a)


################################
# Egg size

m8 <- lm(mean.egg.mass  ~ line, data = data1)
summary(m8)


################################
# Number of eggs

m9 <- lm(eggs.laid  ~ line, data = data1)
summary(m9)


################################
# Prenatal development

m10 <- lmer(development.days ~ line +sex+line:sex+(1|family), data = data)
summary(m10)

# interaction effect
m10a <- lmer(development.days ~ line+sex+ (1|family), data = data)
anova(m10a, m10)

# line effect
m10b <- lmer(development.days ~ sex+ (1|family), data = data)
anova(m10b, m10a)

# sex effect
m10c <- lmer(development.days ~ line+ (1|family), data = data)
anova(m10c, m10a)



################################
# Figure 1

data_new$age <- factor(data_new$age,      
                       levels = c("hatching", "juvenile", "adult"))

t<-data_new %>% 
  ggplot(aes(x=line,y=telomere, fill=factor(sex))) +
  geom_boxplot(aes(middle = mean(telomere)), outlier.size=0.6, outlier.colour= alpha(0.6), notch=FALSE) + 
  labs(fill = "") + 
  xlab("Line")+ 
  ylab("Telomere length (kb)")+ 
  geom_point(colour = "grey26", position=position_jitterdodge(jitter.width = 0.2,
                                                              jitter.height = 0.1,
                                                              dodge.width = 0.75,
                                                              seed = NA),alpha=0.6, size=0.6) +
  theme_bw(base_size = 20)+
  scale_fill_grey(start=0.9, end=0.6 )+
  facet_wrap(~age, ncol = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank())

tag_facet(t)
ggsave("Fig 1.png", width = 26, height = 13, units = "cm",dpi=300 )



# Figure 2
                           
data_new$age <- factor(data_new$age,      
                       levels = c("hatching", "juvenile", "adult"))

data_new2 <- data_new[ which(data_new$age=='hatching'),]

a<-data_new2 %>% 
  ggplot(aes(x=line,y=body.mass, fill=factor(sex))) +
  geom_boxplot(aes(middle = mean(telomere)), show.legend = FALSE, outlier.size=0.6, outlier.colour= alpha(0.6), notch=FALSE) + 
  labs(fill = "") + 
  xlab("Line")+ 
  ylab("Body mass (g)")+
  geom_point(colour = "grey26", position=position_jitterdodge(jitter.width = 0.2,
                                                              jitter.height = 0.1,
                                                              dodge.width = 0.75,
                                                              seed = NA),alpha=0.6, size=0.6) +
  theme_bw(base_size = 18)+
  scale_fill_grey(start=0.9, end=0.6 )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank())+
  theme(legend.position = "none")

a


data_new3 <- data_new[ which(data_new$age=='juvenile'),]

b<-data_new3 %>% 
  ggplot(aes(x=line,y=body.mass, fill=factor(sex))) +
  geom_boxplot(aes(middle = mean(telomere)), outlier.size=0.6, outlier.colour= alpha(0.6), notch=FALSE) + 
  labs(fill = "") + 
  xlab("Line")+ 
  ylab("Body mass change (g)")+
  geom_point(colour = "grey26", position=position_jitterdodge(jitter.width = 0.2,
                                                              jitter.height = 0.1,
                                                              dodge.width = 0.75,
                                                              seed = NA),alpha=0.6, size=0.6) +
  theme_bw(base_size = 18)+
  scale_fill_grey(start=0.9, end=0.6 )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank())

b

plot_grid(a, b, labels = c('a)', 'b)'),label_fontface = "plain",align = "h", axis = "b", rel_widths = c(1, 1.5))


ggsave("Fig 2.png", width = 22, height = 13, units = "cm",dpi=300 )


# Figure 3

data_new$age <- factor(data_new$age,      
                       levels = c("hatching", "juvenile", "adult"))

data_new2 <- data_new[ which(data_new$age=='hatching'),]


c<-data_new2 %>% 
  ggplot(aes(body.mass,telomere, colour=line)) +
  geom_point(size=3)+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE,show_legend = FALSE, linetype="dashed", aes(fill="line"))+
  labs(fill = "") + 
  xlab("Hatching mass (g)")+ 
  ylab("Telomere length (kb)")+
  theme_bw(base_size = 18)+scale_colour_grey(start=0.1, end=0.7 )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank())+
  theme(legend.position = "none")

c

data_new$age <- factor(data_new$age,      
                       levels = c("hatching", "juvenile", "adult"))

data_new3 <- data_new[ which(data_new$age=='juvenile'),]


d<-data_new3 %>% 
  ggplot(aes(body.mass,telomere, colour=line)) +
  geom_point(size=3)+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE, linetype="dashed", aes(fill="line"), show_legend = FALSE)+
  guides(fill = FALSE) +
  labs(fill = "") + 
  xlab("Body mass change (g)")+ 
  ylab("Telomere length (kb)")+
  theme_bw(base_size = 18)+
  scale_colour_grey(start=0.1, end=0.7 )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank())

d


plot_grid(c, d, labels = c('a)', 'b)'),label_fontface = "plain",align = "h", axis = "b", rel_widths = c(1, 1.3))


ggsave("Fig 3.png", width = 22, height = 13, units = "cm",dpi=300 )


