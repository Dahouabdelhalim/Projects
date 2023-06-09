# Packages ----------------------------------------------------------------
library(lme4)
library(plyr)
library(ggplot2)

# Body Length Data --------------------------------------------------------
bodylength <- read.csv("SA Body Lengths.csv")

bodylength$Line <- factor(bodylength$Line) # define fly line as a factor

BL.summarySE <- function(data=NULL, # function that calculates mean and SE
                         measurevar,
                         groupvars=NULL,
                         na.rm=FALSE,
                         conf.interval=.95,
                         .drop=TRUE){
  
  length2 <- function(x,
                      na.rm=FALSE){
    if (na.rm) sum(!is.na(x))
    else length(x)
  }
  
  BLmean <- ddply(bodylength,
                  groupvars,
                  .drop = .drop,
                  .fun = function(xx,col){
                    c(N = length2(xx[[col]],
                                  na.rm=na.rm),
                      mean = mean (xx[[col]],
                                   na.rm=na.rm),
                      sd = sd (xx[[col]],
                               na.rm=na.rm)
                    )
                  },
                  measurevar
  )
  BLmean <- rename(BLmean,
                   c("mean" = measurevar))
  
  BLmean$se <- BLmean$sd/sqrt(BLmean$N)
  
  ciMult <- qt(conf.interval/2 + .5,
               BLmean$N-1)
  BLmean$ci <- BLmean$se * ciMult
  
  return(BLmean)
}


Mean_BL.summarytable <- BL.summarySE(bodylength,  # calculate mean and SE for body length for each sex and genetic background combination
                                     measurevar = "Length.mm.",
                                     groupvars = c("Line",
                                                   "Sex")
)

Mean_BL.summarytable$Line <-  factor(Mean_BL.summarytable$Line) # define genetic background as a factor


# Graphing Body Length ----------------------------------------------------
ggplot(Mean_BL.summarytable,
       aes(x=reorder(Line,
                     Length.mm.),
           y=Length.mm.)
) +
  geom_bar(position = "dodge",
           stat = "identity",
           size=.5,
           color="black",
           aes(fill=Sex)) +
  geom_errorbar(aes(ymin = Length.mm.-se,
                    ymax = Length.mm.+se),
                position = position_dodge(.9),
                width=.2) +
  facet_grid(.~Sex) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60,
                                   hjust = 1),
        text = element_text(size=20)
  ) +
  labs(x=Genetic~Background,
       y=Body~Length~(mm))


# Model body length -------------------------------------------------------
model.BodyLength <- lm(Length.mm.~ Line * Sex,
                       data=bodylength)

summary(model.BodyLength)
anova(model.BodyLength)


# "--------------------" ----------------------------------------------------
# Social Aggregation Data -------------------------------------------------
# Upload Social Aggregation Experiment Dataframe/packages and Annotate
## Attach file "Social Aggregation Measurements.csv" as dataframe titled "SA"
SA <- read.csv("Social Aggregation Measurements.csv")

SA$Line <- factor(SA$Line) # Define genetic background as a factor


SA.summarySE <- function(data=NULL, # function to calculate mean and SE
                         measurevar,
                         groupvars=NULL,
                         na.rm=FALSE,
                         conf.interval=.95,
                         .drop=TRUE){
  
  length2 <- function(x,
                      na.rm=FALSE){
    if (na.rm) sum(!is.na(x))
    else length(x)
  }
  
  SA_mean <- ddply(SA,
                   groupvars,
                   .drop = .drop,
                   .fun = function(xx,col){
                     c(N = length2(xx[[col]],
                                   na.rm=na.rm),
                       mean = mean (xx[[col]],
                                    na.rm=na.rm),
                       sd = sd (xx[[col]],
                                na.rm=na.rm)
                     )
                   },
                   measurevar
  )
  SA_mean <- rename(SA_mean,
                    c("mean" = measurevar))
  
  SA_mean$se <- SA_mean$sd/sqrt(SA_mean$N)
  
  ciMult <- qt(conf.interval/2 + .5,
               SA_mean$N-1)
  SA_mean$ci <- SA_mean$se * ciMult
  
  return(SA_mean)
}



Median_SAmm.summarytable <- SA.summarySE(SA, # Calculate Median Nearest Neighbour Distance in millimetres
                                         measurevar = "MedianNND.mm.1",
                                         groupvars = c("Line",
                                                       "Sex",
                                                       "Infection")
)

Median_SAbl.summarytable <- SA.summarySE(SA, # Calculate Median Nearest Neighbour Distance in Body lengths
                                         measurevar = "MedianNND.BL.1",
                                         groupvars = c("Line",
                                                       "Sex",
                                                       "Infection")
)


Median_SAmm.summarytable$Line <-  factor(Median_SAmm.summarytable$Line) # define genetic background as factor for aggregation in millimetres
Median_SAbl.summarytable$Line <-  factor(Median_SAbl.summarytable$Line) # define genetic background as factor for aggregation in body lengths

# Graphing Social Aggregation (mm) -------------------------------------------------

ggplot(Median_SAmm.summarytable,
       aes(x=reorder(Line,
                     MedianNND.mm.1),
           y=MedianNND.mm.1,
           fill=Infection)
) +
  geom_bar(position = "dodge",
           stat = "identity",
           size=.5) +
  geom_errorbar(aes(ymin = MedianNND.mm.1-se,
                    ymax = MedianNND.mm.1+se),
                position = position_dodge(.9),
                width=.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60,
                                   hjust = 1),
        text = element_text(size=20)
  ) +
  facet_wrap(~Sex) +
  labs(x=Genetic~Background,
       y=Median~NND~(mm)
  )


# model Social Aggregation (mm) -------------------------------------------

model1.SocAggmm <- lm(MedianNND.mm.1 ~ Line * Sex * Infection,
                      data=SA)

summary(model1.SocAggmm)
anova(model1.SocAggmm)



# Graphing Social Aggregation Body Lengths --------------------------------

ggplot(Median_SAbl.summarytable,
       aes(x=reorder(Line,
                     MedianNND.BL.1),
           y=MedianNND.BL.1,
           fill=Infection)
) +
  geom_bar(position = "dodge",
           stat = "identity",
           size=.5) +
  geom_errorbar(aes(ymin = MedianNND.BL.1-se,
                    ymax = MedianNND.BL.1+se),
                position = position_dodge(.9),
                width=.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60,
                                   hjust = 1),
        text = element_text(size=20)
  ) +
  facet_wrap(~Sex) +
  labs(x=Genetic~Background,
       y=Median~NND~(Body~Lengths)
  )


# Model Social Aggregation ------------------------------------------------
model1.SocAggBL <- lm(MedianNND.BL.1 ~ Sex * Line *Infection,
                      data=SA)

summary(model1.SocAggBL)
anova(model1.SocAggBL)


# "--------------------" ----------------------------------------------------
# DAM total Activity data profile -----------------------------------------
## attach "DAM_total.csv" file and label dataframe: "df"
df <- read.csv("DAM_total.csv")

df2 <- df[!(df$Time_mins > 5770),] # exclude time steps greater than 5770
df2$Individual <- factor(df2$Individual) # define individual as a factor
# code to generate mean values --------------------------------------------------------------------
DAM.summarySE <- function(data=NULL, # function to calculate mean and SE
                          measurevar,
                          groupvars=NULL,
                          na.rm=FALSE,
                          conf.interval=.95,
                          .drop=TRUE){
  
  length2 <- function(x,
                      na.rm=FALSE){
    if (na.rm) sum(!is.na(x))
    else length(x)
  }
  
  DAM_mean <- ddply(df2,
                    groupvars,
                    .drop = .drop,
                    .fun = function(xx,col){
                      c(N = length2(xx[[col]],
                                    na.rm=na.rm),
                        mean = mean (xx[[col]],
                                     na.rm=na.rm),
                        sd = sd (xx[[col]],
                                 na.rm=na.rm)
                      )
                    },
                    measurevar
  )
  DAM_mean <- rename(DAM_mean,
                     c("mean" = measurevar))
  
  DAM_mean$se <- DAM_mean$sd/sqrt(DAM_mean$N)
  
  ciMult <- qt(conf.interval/2 + .5,
               DAM_mean$N-1)
  DAM_mean$ci <- DAM_mean$se * ciMult
  
  return(DAM_mean)
}



DAM.summarytable <- DAM.summarySE(df2, # calculate mean and SE values for activity counts at each time point for each combination of sex, genetic background and infection status
                                  measurevar = "Activity_Counts",
                                  groupvars = c("DGRP_Line",
                                                "Sex",
                                                "Infection",
                                                "Time_mins") 
                                  )


DAM.summarytable2 <- DAM.summarytable[!(DAM.summarytable$Sex=="Null"),] # remove empty vials from dataframe

# Plot of mean actuvuty counts -----------------------------------------------------

ggplot(DAM.summarytable2,
       aes(x=Time_mins,
           y=Activity_Counts)
) +
  geom_point(alpha=0.05) +
  geom_smooth(alpha=0.7,
              se=F,
              aes(color=Infection)
  ) +
  facet_grid(Sex~DGRP_Line) +
  theme_bw() +
  ylim(0,15)


# "--------------------" ----------------------------------------------------
# Dataframe for Summary Statistics ----------------------------------------
## Dataframe for summary statistics
# attach file titled 'DAM Summary Statistics.csv' as "DAM"
DAM <- read.csv("DAM Summary Statistics.csv")

DAM <- na.omit(DAM) # Remove empty cells from dtaaframe

# remove all values of Day except 'Sum' which combines all four days
DAM <- DAM[!(DAM$Day=="1"), ]
DAM <- DAM[!(DAM$Day=="2"), ]
DAM <- DAM[!(DAM$Day=="3"), ]
DAM <- DAM[!(DAM$Day=="4"), ]


DAM <- DAM[!(DAM$Alive==0), ] # Remove dead flies from experiment
DAM <- DAM[!(DAM$Line=="Null"), ] # Remove empty slots from dataframe
DAM <- DAM[!(DAM$Line=="Blank"), ] # Remove empty vials from dataframe
DAM <- DAM[!(DAM$Awake.Activity == "#DIV/0!"), ] # Remove flies with no activity from average activity when awake statistic


DAM <- DAM[!(DAM$Total.Activity==0), ] # Remove flies that were not active
DAM$Awake.Activity <- as.numeric(DAM$Awake.Activity) # Define average activity when awake as a numeric value

DAM$Sex <- factor(DAM$Sex,
                  levels=c("Male",
                           "Female"))

# DAM summary statistics --------------------------------------------------
# Graphing Total Activity ----------------------------------------------------------
## MEan±SE function
TA.summarySE <- function(data=NULL, # Function to calculate mean and SE
                         measurevar,
                         groupvars=NULL,
                         na.rm=FALSE,
                         conf.interval=.95,
                         .drop=TRUE){
  
  length2 <- function(x,
                      na.rm=FALSE){
    if (na.rm) sum(!is.na(x))
    else length(x)
  }
  
  DAM_mean <- ddply(DAM,
                    groupvars,
                    .drop = .drop,
                    .fun = function(xx,col){
                      c(N = length2(xx[[col]],
                                    na.rm=na.rm),
                        mean = mean (xx[[col]],
                                     na.rm=na.rm),
                        sd = sd (xx[[col]],
                                 na.rm=na.rm)
                      )
                    },
                    measurevar
  )
  DAM_mean <- rename(DAM_mean,
                     c("mean" = measurevar))
  
  DAM_mean$se <- DAM_mean$sd/sqrt(DAM_mean$N)
  
  ciMult <- qt(conf.interval/2 + .5,
               DAM_mean$N-1)
  DAM_mean$ci <- DAM_mean$se * ciMult
  
  return(DAM_mean)
}
Mean_TA.summarytable <- TA.summarySE(DAM, # calculate average total activity for each combination of sex, genetic background and ifnection treatment
                                     measurevar = "Total.Activity",
                                     groupvars = c("Line",
                                                   "Sex",
                                                   "Infection_treatment")
)

Mean_TA.summarytable$Line <-  factor(Mean_TA.summarytable$Line) # Define genetic background as a factor




# Graphing Total Activity Summary Stat ------------------------------------
ggplot(Mean_TA.summarytable,
                         aes(x=reorder(Line,
                                       Total.Activity),
                             y=Total.Activity,
                             fill=Infection_treatment)
) +
  geom_bar(position = "dodge",
           stat = "identity",
           size=.5,
           color="black") +
  geom_errorbar(aes(ymin = Total.Activity-se,
                    ymax = Total.Activity+se),
                position = position_dodge(.9)
                ) +
  facet_grid(.~Sex) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60,
                                   hjust = 1),
        text = element_text(size=20)
  ) +
  labs(x=Genetic~Background,
       y=Total~Activity)



# model of total activity -------------------------------------------------
model1.TotalActivity <- lm(log(Total.Activity+1) ~ Line*Sex*Infection_treatment,
            data=DAM)
summary(model1.TotalActivity)
anova(model1.TotalActivity)


# Graphing Proportion of time Awake  ----------------------------------------------------------------------
## MEan±SE
PA.summarySE <- function(data=NULL, # function to calculate mean and SE for proportion of time awake
                         measurevar,
                         groupvars=NULL,
                         na.rm=FALSE,
                         conf.interval=.95,
                         .drop=TRUE){
  
  length2 <- function(x,
                      na.rm=FALSE){
    if (na.rm) sum(!is.na(x))
    else length(x)
  }
  
  DAM_mean <- ddply(DAM,
                    groupvars,
                    .drop = .drop,
                    .fun = function(xx,col){
                      c(N = length2(xx[[col]],
                                    na.rm=na.rm),
                        mean = mean (xx[[col]],
                                     na.rm=na.rm),
                        sd = sd (xx[[col]],
                                 na.rm=na.rm)
                      )
                    },
                    measurevar
  )
  DAM_mean <- rename(DAM_mean,
                     c("mean" = measurevar))
  
  DAM_mean$se <- DAM_mean$sd/sqrt(DAM_mean$N)
  
  ciMult <- qt(conf.interval/2 + .5,
               DAM_mean$N-1)
  DAM_mean$ci <- DAM_mean$se * ciMult
  
  return(DAM_mean)
}
Mean_PA.summarytable <- PA.summarySE(DAM, # calculate mean and SE for proportion of time awake for each combination of line, sex and infection treatment
                                     measurevar = "Proportion.Activity",
                                     groupvars = c("Line",
                                                   "Sex",
                                                   "Infection_treatment")
)

Mean_PA.summarytable$Line <-  factor(Mean_PA.summarytable$Line) # Define genetic background as a factor


# Graphing Proportion of Time Spent Awake ---------------------------------
ggplot(Mean_PA.summarytable,
       aes(x=reorder(Line,
                     Proportion.Activity),
           y=Proportion.Activity,
           fill=Infection_treatment)
) +
  geom_bar(position = "dodge",
           stat = "identity",
           size=.5,
           color="black") +
  geom_errorbar(aes(ymin = Proportion.Activity-se,
                    ymax = Proportion.Activity+se),
                position = position_dodge(.9)
  ) +
  facet_grid(.~Sex) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60,
                                   hjust = 1),
        text = element_text(size=20)
  ) +
  labs(x=Genetic~Background,
       y=Proportion~Activity)



# model of proportion -----------------------------------------------------

model1.Prop.Activity <- lm(Proportion.Activity ~ Line*Sex*Infection_treatment,
            data=DAM)

summary(model1.Prop.Activity)

anova(model1.Prop.Activity)


# Graphing Activity when awake -----------------------------------------------------
## MEan±SE
AA.summarySE <- function(data=NULL, # function to calculate the mean and SE of average awake activitiy
                         measurevar,
                         groupvars=NULL,
                         na.rm=FALSE,
                         conf.interval=.95,
                         .drop=TRUE){
  
  length2 <- function(x,
                      na.rm=FALSE){
    if (na.rm) sum(!is.na(x))
    else length(x)
  }
  
  DAM_mean <- ddply(DAM,
                    groupvars,
                    .drop = .drop,
                    .fun = function(xx,col){
                      c(N = length2(xx[[col]],
                                    na.rm=na.rm),
                        mean = mean (xx[[col]],
                                     na.rm=na.rm),
                        sd = sd (xx[[col]],
                                 na.rm=na.rm)
                      )
                    },
                    measurevar
  )
  DAM_mean <- rename(DAM_mean,
                     c("mean" = measurevar))
  
  DAM_mean$se <- DAM_mean$sd/sqrt(DAM_mean$N)
  
  ciMult <- qt(conf.interval/2 + .5,
               DAM_mean$N-1)
  DAM_mean$ci <- DAM_mean$se * ciMult
  
  return(DAM_mean)
}
Mean_AA.summarytable <- AA.summarySE(DAM, # calculate mean and SE average awake activity for each treatment group combination of sex, genetic background and infection treatment
                                     measurevar = "Awake.Activity",
                                     groupvars = c("Line",
                                                   "Sex",
                                                   "Infection_treatment")
)

Mean_AA.summarytable$Line <-  factor(Mean_AA.summarytable$Line) # Define genetic background as a factor



# Graphing Average Awake Activity -----------------------------------------
ggplot(Mean_AA.summarytable,
       aes(x=reorder(Line,
                     Awake.Activity),
           y=Awake.Activity,
           fill=Infection_treatment)
) +
  geom_bar(position = "dodge",
           stat = "identity",
           size=.5,
           color="black") +
  geom_errorbar(aes(ymin = Awake.Activity-se,
                    ymax = Awake.Activity+se),
                position = position_dodge(.9)
  ) +
  facet_grid(.~Sex) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60,
                                   hjust = 1),
        text = element_text(size=20)
  ) +
  labs(x=Genetic~Background,
       y=Awake~Activity)



# Model Activity when awake  -----------------------------------------------
model1.AwakeAct <- lm(log(Awake.Activity+1) ~ Line*Sex*Infection_treatment,
            data=DAM)
summary(model1.AwakeAct)
anova(model1.AwakeAct)





