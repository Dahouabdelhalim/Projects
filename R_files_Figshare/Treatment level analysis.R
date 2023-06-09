#----------------------------------------------------------------------------
# This code is associated with  the paper 
# "Collaborative Dishonesty: A Meta-Analysis' by Margarita Leib, Nils C. Köbis, Ivan Soraperra, Ori Weisel, & Shaul Shalvi. 
# and contains the code for the treatment level analysis in the paper
# the code is associated with the data file 'treatment Level Data'
# the data set and read me files can be found in https://doi.org/10.21942/uva.14731593.v1  
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# load libraries 
library(readxl)
library(meta)
library(metafor)
library(dmetar)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# load data
alltasks <- read_excel("treatment Level Data.xlsx")
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
#  Ns for PRISMA CHART
#----------------------------------------------------------------------------

# n tasks

length(unique(alltasks$task.number[alltasks$decision.structure == "joint"]))
length(unique(alltasks$task.number[alltasks$decision.structure == "simultaneous"])) 
length(unique(alltasks$task.number[alltasks$decision.structure == "sequential"]))

# n papers

length(unique(alltasks$paper[alltasks$decision.structure == "joint"]))
length(unique(alltasks$paper[alltasks$decision.structure == "simultaneous"])) 
length(unique(alltasks$paper[alltasks$decision.structure == "sequential"]))

# n treatments

length(alltasks$treatment[alltasks$decision.structure == "joint"])
length(alltasks$treatment[alltasks$decision.structure == "simultaneous"]) 
length(alltasks$treatment[alltasks$decision.structure == "sequential"])


## N papers, participants, reports, treatments

length(unique(alltasks$paper)) # number of papers
length(unique(alltasks$paper[alltasks$published==0])) #number of unpublished papers
sum(alltasks$n.reports) #number of decisions 
length(alltasks$treatment) # number of treatments
length(unique(alltasks$task.number)) #number of tasks

# number of participants

#set the number of within to remove

alltasks$within.Subject.With[alltasks$within.Subject==1]


m1 <- alltasks$n.participants[alltasks$treatment== "PC, round 2, leader"] 
m2 <- alltasks$n.participants[alltasks$treatment== "PC, round 2, Ethical code"] 
m3 <- alltasks$n.participants[alltasks$treatment== "Out-group, aligned"] 
m4 <- alltasks$n.participants[alltasks$treatment== "Out-group, B fixed"] 


#number of participants (sum all minus those who repeat in within subject treatments)
sum(alltasks$n.participants) - sum(m1, m2, m3, m4)




#----------------------------------------------------------------------------
## This section calculates the frequency of country origin in the dataset  
#----------------------------------------------------------------------------

# re-name location 

alltasks$country[alltasks$country == "Mturk"] <- "Mturk/Prolific"
alltasks$country[alltasks$country == "Prolific"] <- "Mturk/Prolific"
alltasks$country[alltasks$country == "Mturk (USA)"] <- "USA"
alltasks$country[alltasks$country == "NL"] <- "The Netherlands"

# before calculating Number of participants per country, removing those who are within (to not have them counted twice)

alltasks$remove.for.country <- 0
alltasks$remove.for.country[alltasks$treatment== "Out-group, aligned"] <- 1 
alltasks$remove.for.country[alltasks$treatment== "Out-group, B fixed"] <- 1 
alltasks$remove.for.country[alltasks$treatment== "PC, round 2, leader"] <- 1 
alltasks$remove.for.country[alltasks$treatment== "PC, round 2, Ethical code"] <- 1 



# dataset for country frequency  

t <- aggregate(alltasks$n.participants[alltasks$remove.for.country ==0], by = list(alltasks$country[alltasks$remove.for.country ==0]), FUN= sum)

names(t) <- c("country", "N.participant")
t <- t[order(-t$N.participant),]
t$prop <- t$N.participant/sum(t$N.participant)
t$prop.round <- round(t$prop*100, digits = 2)

sum(t$N.participant)





#------------------------------------------------------------------------------------------------------------
#---------------                Note                       --------------------------------------------------
# The data set published online contains missing values for the standardized group report for 3 treatments.
# This is because we did not get approval from the authors to publish these values. 
# The values are for a sequential task (the dyadic die rolling task, task number 16) in a yet unpublished paper
# Thus, the results from running any section of the code where these values should be included 
# (e.g., all decision structures combined, only sequential tasks, all unpublished treatments)
# would lead to a result slightly different than the one reported in the paper
# (in the paper we do include the data from these 3 treatments)
#------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------
## This section calculates the descriptive information in Table 4 in the paper  
#----------------------------------------------------------------------------

# a function to generate mode + n of mode

Mode <- function(x) {
        ux <- table(x)
        cbind(names(ux[ux == max(ux)]), ux[ux == max(ux)])
}


#
# -------------------------------------------------------------------
### descriptive data for all decision structures combined 
#--------------------------------------------------------------------

# standardized group report 

#note that the data set published online contains 3 N/A values (in a sequential task, as we did not get approval from the authors to publish these values)
#thus the results from running this section of the code (for standardized group report) will not match the values reported in the table 4 in the paper

alltasks$standardized.group.report <- as.numeric(alltasks$standardized.group.report)
summary(alltasks$standardized.group.report)
sd(alltasks$standardized.group.report)
table(alltasks$standardized.group.report, useNA = "always")
Mode(alltasks$standardized.group.report)
length(alltasks$standardized.group.report)

# negative consequences for third parties 

# combine 1 (negative consequences for other participants) and 2 (negative consequences for charity) into one coding of negative consequences (0=no, 1=yes)

alltasks$negative.consequences.for.others[alltasks$negative.consequences.for.others == 2] <- 1

prop.table(table(alltasks$negative.consequences.for.others))
Mode(alltasks$negative.consequences.for.others)


# group size

summary(alltasks$group.size)
sd(alltasks$group.size)
Mode(alltasks$group.size)


# financial incentive to lie

# calculate financial incentive 

alltasks$max.minus.min.ppp <- alltasks$Max.group.can.earn.by.lying.USD.2015.PPP - alltasks$Min.group.can.earn.by.lying.USD.2015.PPP

summary(alltasks$max.minus.min.ppp)
sd(alltasks$max.minus.min.ppp)
Mode(alltasks$max.minus.min.ppp)

#payoff alignment 

prop.table(table(alltasks$aligned.payoff))
Mode(alltasks$aligned.payoff)

# study type 

# code study type lab in the field and field as one 

alltasks$study.type[alltasks$study.type == "lab in the field"] <- "field"

prop.table(table(alltasks$study.type))
Mode(alltasks$study.type)


# experimental deception 

prop.table(table(alltasks$experimetal.deception))
Mode(alltasks$experimetal.deception)


# repeated task 

# code for one shot vs. repeated

alltasks$repeated <- 0
alltasks$repeated[alltasks$N.rounds > 1] <- 1

prop.table(table(alltasks$repeated))
Mode(alltasks$repeated)


# among repeated tasks, number of rounds in the task

summary(alltasks$N.rounds[alltasks$repeated==1])
length(alltasks$N.rounds[alltasks$repeated==1])
sd(alltasks$N.rounds[alltasks$repeated==1])
Mode(alltasks$N.rounds[alltasks$repeated==1])

# year study run 

summary(alltasks$year.exp.run)
sd(alltasks$year.exp.run)
Mode(alltasks$year.exp.run)

# published 

prop.table(table(alltasks$published))
Mode(alltasks$published)


# Note that age and gender was calculated on the raw data (because we examine the average age/proportion of females within a group). Thus, the code associated with these analyses do not appear here

# -------------------------------------------------------------------
### descriptive data for joint tasks  
#--------------------------------------------------------------------

### restrict the data set for joint tasks only 

joint <- alltasks[alltasks$decision.structure == "joint", ]

# standardized group report 

summary(joint$standardized.group.report)
sd(joint$standardized.group.report)
Mode(joint$standardized.group.report)
length(joint$standardized.group.report)

# negative consequences for third parties 

prop.table(table(joint$negative.consequences.for.others))
Mode(joint$negative.consequences.for.others)


# group size

summary(joint$group.size)
sd(joint$group.size)
Mode(joint$group.size)


# financial incentive to lie

summary(joint$max.minus.min.ppp)
sd(joint$max.minus.min.ppp)
Mode(joint$max.minus.min.ppp)

#payoff alignment 

prop.table(table(joint$aligned.payoff))
Mode(joint$aligned.payoff)

# study type 

prop.table(table(joint$study.type))
Mode(joint$study.type)


# experimental deception 

prop.table(table(joint$experimetal.deception))
Mode(joint$experimetal.deception)


# repeated task 

prop.table(table(joint$repeated))
Mode(joint$repeated)


# year study run 

summary(joint$year.exp.run)
sd(joint$year.exp.run)
Mode(joint$year.exp.run)

# published 

prop.table(table(joint$published))
Mode(joint$published)


# Note that age and gender was calculated on the raw data (because we examine the average age/proportion of females within a group). Thus, the code associated with these analyses do not appear here

# -------------------------------------------------------------------
### descriptive data for simultaneous tasks  
#--------------------------------------------------------------------

### restrict the data set for simultaneous tasks only 

sim <- alltasks[alltasks$decision.structure == "simultaneous", ]

# standardized group report 

summary(sim$standardized.group.report)
sd(sim$standardized.group.report)
Mode(sim$standardized.group.report)
length(sim$standardized.group.report)

# negative consequences for third parties 

prop.table(table(sim$negative.consequences.for.others))
Mode(sim$negative.consequences.for.others)


# group size

summary(sim$group.size)
sd(sim$group.size)
Mode(sim$group.size)


# financial incentive to lie

summary(sim$max.minus.min.ppp)
sd(sim$max.minus.min.ppp)
Mode(sim$max.minus.min.ppp)

#payoff alignment 

prop.table(table(sim$aligned.payoff))
Mode(sim$aligned.payoff)

# stydy type 

prop.table(table(sim$study.type))
Mode(sim$study.type)


# experimental deception 

prop.table(table(sim$experimetal.deception))
Mode(sim$experimetal.deception)


# repeated task 

prop.table(table(sim$repeated))
Mode(sim$repeated)


# among repeated tasks, number of rounds in the task

summary(sim$N.rounds[sim$repeated==1])
length(sim$N.rounds[sim$repeated==1])
sd(sim$N.rounds[sim$repeated==1])
Mode(sim$N.rounds[sim$repeated==1])

# year study run 

summary(sim$year.exp.run)
sd(sim$year.exp.run)
Mode(sim$year.exp.run)

# published 

prop.table(table(sim$published))
Mode(sim$published)


# Note that age and gender was calculated on the raw data (because we examine the average age/proportion of females within a group). Thus, the code associated with these analyses do not appear here

# -------------------------------------------------------------------
### descriptive data for sequwntial tasks  
#--------------------------------------------------------------------

### restrict the data set for sequential tasks only 

seq <- alltasks[alltasks$decision.structure == "sequential", ]

# standardized group report 

#note that the dataset published online contains 3 N/A values (in a sequential task, as we did not get approval from the authors to publish these values)
#thus the results from running this section of the code (for standardized group report) will not match the values reported in the table 4 in the paper

summary(seq$standardized.group.report)
sd(seq$standardized.group.report)
Mode(seq$standardized.group.report)
length(seq$standardized.group.report)

# negative consequences for third parties 

prop.table(table(seq$negative.consequences.for.others))
Mode(seq$negative.consequences.for.others)

# group size

summary(seq$group.size)
sd(seq$group.size)
Mode(seq$group.size)


# financial incentive to lie

summary(seq$max.minus.min.ppp)
sd(seq$max.minus.min.ppp)
Mode(seq$max.minus.min.ppp)

#payoff alignment 

prop.table(table(seq$aligned.payoff))
Mode(seq$aligned.payoff)

# study type 

prop.table(table(seq$study.type))
Mode(seq$study.type)


# experimental deception 

prop.table(table(seq$experimetal.deception))
Mode(seq$experimetal.deception)


# repeated task 

prop.table(table(seq$repeated))
Mode(seq$repeated)


# among repeated tasks, number of rounds in the task

summary(seq$N.rounds[seq$repeated==1])
length(seq$N.rounds[seq$repeated==1])
sd(seq$N.rounds[seq$repeated==1])
Mode(seq$N.rounds[seq$repeated==1])

# year study run 

summary(seq$year.exp.run)
sd(seq$year.exp.run)
Mode(seq$year.exp.run)

# published 

prop.table(table(seq$published))
Mode(seq$published)


# Note that age and gender was calculated on the raw data (because we examine the average age/proportion of females within a group). Thus, the code associated with these analyses do not appear here


#----------------------------------------------------------------------------
## excluding Sutter, 2009; T2, groups of 3 because SD=0
alltasks <- alltasks[alltasks$sd.standardized.group.report !=0, ]
#----------------------------------------------------------------------------


## N papers, participants, reports after exclusion

length(unique(alltasks$paper))
length(alltasks$treatment)
sum(alltasks$n.participants) - sum(m1, m2, m3, m4)
sum(alltasks$n.reports)

#----------------------------------------------------------------------------
## analyses of Standardized group report -- overall effect and heterogeneity
#----------------------------------------------------------------------------

# make sure vars are numeric 

alltasks$standardized.group.report <- as.numeric(alltasks$standardized.group.report)
alltasks$sd.standardized.group.report <- as.numeric(alltasks$sd.standardized.group.report)


summary(alltasks$standardized.group.report)

# create a treatment per paper label 

alltasks$paper.treatment <- interaction(alltasks$paper, alltasks$treatment)
alltasks$paper.treatment <- droplevels(alltasks$paper.treatment)

## calculated SE = sd/sqrt(n groups)

alltasks$se <- alltasks$sd.standardized.group.report/sqrt(alltasks$n)  


### overall effect and heterogeneity

# joint decisions

m.dl.j <- metagen(standardized.group.report,
                  se,
                  data = alltasks[alltasks$decision.structure == "joint", ],
                  studlab = paste(paper.treatment),
                  comb.fixed = FALSE,
                  comb.random = TRUE,
                  method.tau = "DL",  # the default estimator
                  hakn = F,
                  prediction = TRUE,
                  sm = "SMD")

m.dl.j

# simultaneous decisions


m.dl.si <- metagen(standardized.group.report,
                   se,
                   data = alltasks[alltasks$decision.structure == "simultaneous", ],
                   studlab = paste(paper.treatment),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",  # the default estimator
                   hakn = F,
                   prediction = TRUE,
                   sm = "SMD")

m.dl.si

# sequential decisions

#note that the dataset published online contains 3 N/A values (in a sequential task, as we did not get approval from the authors to publish these values)
#thus the results from running this section of the code (for standardized group report) will not match the values reported in the table 4 in the paper

m.dl.se <- metagen(standardized.group.report,
                   se,
                   data = alltasks[alltasks$decision.structure == "sequential", ],
                   studlab = paste(paper.treatment),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",  # the default estimator
                   hakn = F,
                   prediction = TRUE,
                   sm = "SMD")

m.dl.se

# all decisions combined (for heterogeneity)
#note that the dataset published online contains 3 N/A values (in a sequential task, as we did not get approval from the authors to publish these values)
#thus the results from running this section of the code (for standardized group report) will not match the values reported in the table 4 in the paper


m.dl.all <- metagen(standardized.group.report,
                   se,
                   data = alltasks,
                   studlab = paste(paper.treatment),
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   method.tau = "DL",  # the default estimator
                   hakn = F,
                   prediction = TRUE,
                   sm = "SMD")


m.dl.all

## examples of  Standardized group report for tasks 1 and 2

# task 16 (dyadic die-rolling task)
# the value here will be different than the one in the paper as the 3 N/A vaules
#in the dataset used this task

metagen(standardized.group.report,
        se,
        data = alltasks[alltasks$task.number==16,],
        studlab = paste(paper.treatment),
        comb.fixed = FALSE,
        comb.random = TRUE,
        method.tau = "DL",  # the default estimator
        hakn = F,
        prediction = TRUE,
        sm = "SMD")

#task 7
metagen(standardized.group.report,
        se,
        data = alltasks[alltasks$task.number==7,],
        studlab = paste(paper.treatment),
        comb.fixed = FALSE,
        comb.random = TRUE,
        method.tau = "DL",  # the default estimator
        hakn = F,
        prediction = TRUE,
        sm = "SMD")

# examples from Gross and De Dreu (2021)
alltasks$standardized.group.report[alltasks$treatment == "HH (exp 2)"]
alltasks$standardized.group.report[alltasks$treatment == "LL (exp 2)"]


# -----------------------------------
# Outlier and influence diagnostics analyses (Viechtbauer & Cheung, 2010)
# https://www.metafor-project.org/doku.php/plots:plot_of_influence_diagnostics?s[]=outlier
# -----------------------------------


# note that the outcome will be slightly different because of the 3 N/A values in the data set published

m.dl.rma <- rma(yi = standardized.group.report,
                sei = se,
                data = alltasks, method = "DL", test = "knha")

inf.rma <- influence(m.dl.rma)
inf.rma
inf.rma$inf[100:120]


#----------------------------------------------------------------------------
# LEAVE ONE OUT ANALYSES 
#----------------------------------------------------------------------------


# all decisions structures 

# note that the outcome will be slightly different because of the 3 N/A values in the data set published


inf.analysis <- InfluenceAnalysis(x = m.dl.all,
                                  random = TRUE)
#I^2
min(inf.analysis$Data[5])
max(inf.analysis$Data[5])

#effect size
min(inf.analysis$Data[2])
max(inf.analysis$Data[2])

summary(inf.analysis)


## for joint decisions

inf.analysis <- InfluenceAnalysis(x = m.dl.j,
                                  random = TRUE)

#I^2
min(inf.analysis$Data[5])
max(inf.analysis$Data[5])

# EFFECT SIZE

min(inf.analysis$Data[2])
max(inf.analysis$Data[2])

summary(inf.analysis)

## for simultaneous decisions


inf.analysis <- InfluenceAnalysis(x = m.dl.si,
                                  random = TRUE)
#I^2

min(inf.analysis$Data[5])
max(inf.analysis$Data[5])

#effect size
min(inf.analysis$Data[2])
max(inf.analysis$Data[2])

summary(inf.analysis)


## for sequential decisions
# note that the outcome will be slightly different because of the 3 N/A values in the data set published

inf.analysis <- InfluenceAnalysis(x = m.dl.se,
                                  random = TRUE)
#I^2

min(inf.analysis$Data[5])
max(inf.analysis$Data[5])

#effect size
min(inf.analysis$Data[2])
max(inf.analysis$Data[2])

summary(inf.analysis)

#---------------------------------------------------------------------------------------
## publication bias (testing for difference between published and unpublished papers)
#---------------------------------------------------------------------------------------

## note that the outcome will be slightly different because of the 3 N/A values in the data set published

m <- rma.mv(standardized.group.report, 
            sd.standardized.group.report, 
            random = ~ 1 | paper, 
            tdist = TRUE, 
            data = alltasks,
            method = "REML", 
            mods = ~ published)


m

#unpublished (note that the outcome will be slightly different because of the 3 N/A values in the data set published)

metagen(standardized.group.report,
        se,
        data = alltasks[alltasks$published==0, ],
        studlab = paste(paper.treatment),
        comb.fixed = FALSE,
        comb.random = TRUE,
        method.tau = "DL",  # the default estimator
        hakn = F,
        prediction = TRUE,
        sm = "SMD")

#published

metagen(standardized.group.report,
        se,
        data = alltasks[alltasks$published==1,],
        studlab = paste(paper.treatment),
        comb.fixed = FALSE,
        comb.random = TRUE,
        method.tau = "DL",  # the default estimator
        hakn = F,
        prediction = TRUE,
        sm = "SMD")

# for joint decisions

rma.mv(standardized.group.report, 
        sd.standardized.group.report, 
        random = ~ 1 | paper, 
        tdist = TRUE, 
        data = alltasks[alltasks$decision.structure == "joint",],
        method = "REML", 
        mods = ~ published)

# for simultaneous decisions

rma.mv(standardized.group.report, 
       sd.standardized.group.report, 
       random = ~ 1 | paper, 
       tdist = TRUE, 
       data = alltasks[alltasks$decision.structure == "simultaneous",],
       method = "REML", 
       mods = ~ published)

# for sequential decisions
## note that the outcome will be slightly different because of the 3 N/A values in the data set published

rma.mv(standardized.group.report, 
       sd.standardized.group.report, 
       random = ~ 1 | paper, 
       tdist = TRUE, 
       data = alltasks[alltasks$decision.structure == "sequential",],
       method = "REML", 
       mods = ~ published)

