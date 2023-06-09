rm(list=ls())

###change the working directory here
setwd("/Users/yinwu0407/Documents/PhD/Testosterone_donation_2019/data/")

raw_data <- read.csv("donation_raw_data.csv", header = T)


####session 1, public; session 0, private
raw_data$session[raw_data$session == 2] <- 0

###order, 1 public-private; 0 private-public

raw_data$order[raw_data$order == 2] <- 0

###, 1 testosterone; 0 placebo
raw_data$group[raw_data$group == 2] <- 0





library(lme4)
library(ggplot2)
library(languageR)
library(lmerTest)
library(lattice)  
library(pastecs)

options(scipen=999)


head(raw_data)

model1 <- glmer(donation ~ he + you + session*group + (1+ he + you |subject), family = binomial(), data = raw_data)

summary(model1)
coef(model1)


placebo <- raw_data[raw_data$group == 0, ]

model2 <- glmer(donation ~ he + you + session +  (1+ he + you|subject), family = binomial(), data = placebo  )
summary(model2)

testosterone <- raw_data[raw_data$group == 1, ]



model3 <- glmer(donation ~ he + you + session +  (1+ he + you |subject), family = binomial(), data = testosterone  )
summary(model3)


model4 <- glmer(donation ~   he + you + group + session + session:group + he:group + you:group + (1+ he + you|subject), family = binomial(), data = raw_data  )
summary(model4)


raw_data$rt_log <- log(raw_data$rt)
names(raw_data )


by(raw_data$rt, raw_data$session, mean)
by(raw_data$rt, raw_data$session, sd)

model5 <- lmer(rt_log ~   he + you + group + session + session:group + (1+ he + you|subject),  data = raw_data  )
summary(model5)

AIC(model5)
BIC(model5)




personality <- read.csv("personality_summary_osf.csv", header = T)
names(personality)


raw_data <- merge(raw_data, personality, by = c("subject"))
names(raw_data)


model6 <- glmer(donation ~ he + you + session*group.x + BIS_all + IRI_all + PPI_total +
                  AQ_total + MachIV + SIAS + 
                  (1+ he + you |subject), family = binomial(), data = raw_data)

summary(model6)

table(raw_data$group.x)

placebo <- raw_data[raw_data$group.x == 0, ]

model7 <- glmer(donation ~ he + you + session  + BIS_all + IRI_all + PPI_total +
                  AQ_total + MachIV + SIAS +  (1+ he + you|subject), family = binomial(), data = placebo  )
summary(model7)

testosterone <- raw_data[raw_data$group.x == 1, ]

model8 <- glmer(donation ~ he + you + session  + BIS_all + IRI_all + PPI_total +
                  AQ_total + MachIV + SIAS + (1+ he + you |subject), family = binomial(), data = testosterone  )
summary(model8)



raw_data$rt_log <- log(raw_data$rt)

model9 <- lmer(rt_log ~   he + you + group.x + session + session:group.x + BIS_all + IRI_all + PPI_total +
                 AQ_total + MachIV + SIAS + (1+ he + you|subject),  data = raw_data )
summary(model9)

AIC(model9)
BIC(model9)
