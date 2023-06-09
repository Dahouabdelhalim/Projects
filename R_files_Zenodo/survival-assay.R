# ------ GOAL 1 : SURVIAVL DIFFERENCES OF WOLBACHIA INFECTED AND UNINFECTED WORKERS
# ------ GOAL 2 : SURVIAVL DIFFERENCES OF WOLBACHIA INFECTED AND UNINFECTED QUEENS AND ARISING COLONIES
# install and load packages ----
library("lme4") #LM, LMM, GLM, GLMM
library("dplyr")
library("car") #anova for LM and LMM
library("multcomp")
library("MASS") #nelson-aalen estimator for cumulative hazard
library("ggfortify") #diagnostic plots of models

library("Rmisc") #correlations, summarize function, overrides DPLYR's group-by function
library("rstatix") #quick statistical comparisons such as Wilcox

library("Hmisc") #correlation significance
library("corrplot") #ggcorr plots
library("reshape2")
library("RColorBrewer")

library("survival") #cox-hazard and Kaplan-Meier survival tests
library("survminer") #better looking Kaplan-Meier plots
library("gtsummary") #extracting summaries from models
library("webshot") #saving gtsummary tables as png files

library("ggplot2") #for all graphs

# resources for consideration ----
# we are doing "time to event" analysis > time till the event (=death) happens
# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
# http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Survival/BS704_Survival_print.html
# https://www.theanalysisfactor.com/the-six-types-of-survival-analysis-and-challenges-in-learning-them/
# https://towardsdatascience.com/survival-analysis-part-a-70213df21c2e
# https://rpkgs.datanovia.com/survminer/index.html
# https://bookdown.org/sestelo/sa_financial/how-to-evaluate-the-ph-assumption.html
# http://www.sthda.com/english/wiki/survival-analysis-basics

# shortcuts for directories -----
#directory with raw unprocessed files
raw.dir <- ("G:/My Drive/R/2019-2020-survival/input")
#directory with processed files
out.dir <- ("G:/My Drive/R/2019-2020-survival/output")
gr.qn <- ("G:/My Drive/R/2019-2020-survival/graphs/queens")
gr.wrk <- ("G:/My Drive/R/2019-2020-survival/graphs/workers")

# statistics for GLMM/GLM ----
stats.glmer <- function(model){
  drop1(model,.~.,test="Chisq")
}

# quasi-models
stats.glm <- function(model){
  drop1(model,.~.,test="F")
}

  
# Q-Q plot of model ----
qq <- function(model){
  res <- residuals(model)
  qqnorm(res)
  qqline(res, probs=c(0.25,0.75))
}

# Histogam of residuals ---
res <- function(model) {
  res1 <- residuals(model)
  hist(res1)
}

# Normality of the residuals ---
shp.tst <- function(model) {
  res <- residuals(model)
  shapiro.test(res)
}

# Normality of the residuals ----
shp.tst <- function(model) {
  res <- residuals(model)
  shapiro.test(res)
}

# --------------------------- WORKER: SURVIVAL ANALYSIS -------
# loading processed file for analysis
setwd(out.dir)
df.wk <- read.csv("worker2019-survival.csv")
wk.time.death <- read.csv("worker2019-survival-time.csv")

# Raw file: variables ----
setwd(raw.dir)
list.files() #visualize files in the folder

df.wk <- read.csv("work-surv2019-long.csv")  #long-format for worker survival
df.wk$wolbachia <- as.factor(df.wk$wolbachia)
df.wk$colony.source <- as.factor(df.wk$colony.source)
df.wk$sub.colony <- as.factor(df.wk$sub.colony)

# change in queen number wrt to first census (overall) and wrt to previous time point
df.wk <- df.wk %>%
  dplyr::group_by(groupID) %>%
  dplyr::arrange(day, .by_group = TRUE) %>%
  dplyr::mutate(wk.dead = worker - dplyr::lag(worker, # lag = previous row
                              n = 1, # 1 row difference
                              default = first(worker)), # dead workers compared to previous time point
         wk.dead = abs(wk.dead), #making negative values as positive
         wk.dead.overall = first(worker)-worker, # number of dead workers since start
         total.workers = wk.dead.overall + worker, # confirming that the numbers match up# dead workers compared to starting
         wrk.surv.prp.overall = round((worker/first(worker)), digits = 3),# survival proportion wrt start
         wk.dead.prp = round((wk.dead/(worker+wk.dead)), digits = 3)) %>% # death proportion wrt previous time point
  dplyr::rename(wk.alive = worker)

# assigning death status
# if alive workers = 0 then death = 2 (group is dead), if number of alive workers at census > 0 then death = 1 (group is alive)
df.wk <- df.wk %>%
  dplyr::mutate(death = ifelse(wk.alive == 0, 2, 1))

#testing if it worked
df <- data.frame(subset(df.wk1, group == 12|group == 27|group == 46))

df <- data.frame(subset(df.wk1, death == 2))

# PLOT: trend of overall proportion of dead workers over time
data <- df.wk
data$x <- data$worker.age
data$y <- data$wrk.surv.prp.overall

ggplot(data, aes(x=x, y=y, group = groupID,
                 color = wolbachia, fill = wolbachia)) +
  geom_line()+
  facet_wrap(~groupID)

# Time till death of the group: data file -----
wk.time.death <- df.wk %>%
  dplyr::group_by(groupID) %>% # group by experimental group ID
  dplyr::arrange(day, .by_group = T ) %>% # arrange in increasing order of the day of census
  dplyr::filter(wk.dead.prp == 1) %>% # first time point at which the group was dead. Subsequent time points after death have 'NA' values
  dplyr::rename(worker.age.death = worker.age,  # age of workers at when the group died 
                time.death = day) # time in days passed till death of the group

#another way of doing this
df1 <- df.wk %>%
  dplyr::group_by(groupID) %>% # group by experimental group ID
  dplyr::arrange(day, .by_group = T ) %>% # arrange in increasing order of the day of census
  dplyr::filter(wk.alive == 0)

df1 <- df1[complete.cases(df1$wk.dead.prp), ]

identical(df1,wk.time.death)
all.equal(df1,wk.time.death) #are dataframes from both approaches same? (Names of 2 columns differ)

# Does the time till death of the group differ between infected and uninfected?
boxplot(worker.age.death~wolbachia, data = wk.time.death, notch = T)

# Worker : saving this processed file ----
# saving csv file
write.csv(df.wk, #data frame
          file.path(out.dir,# file path
                    "worker2019-survival-census.csv"), #file name
          row.names=T)

write.csv(wk.time.death, #data frame
          file.path(out.dir,# file path
                    "worker2019-survival-time.csv"), #file name
          row.names=T)

# renaming column names with more explanation
colnames(df.wk) <- c("group ID",
                     "day of census",
                     "estimated age of workers",
                     "source colony of workers",
                     "replicate created from source colony",
                     "wolbachia",
                     "workers alive at the time of census",
                     "workers dead since previous census",
                     "workers dead since start",
                     "number of workers in the colony (dead and alive)",
                     "proportion of workers surviving since start",
                     "proportion of workers dead since previous census",
                     "death (1=no, 2 = yes)")

write.csv(df.wk, #data frame
          file.path(out.dir,# file path
                    "worker2019-survival-longnames.csv"), #file name
          row.names=T)

# Kaplan-Meir (KP) method for fitted survival ----------
# first pass survival probability estimate
# non-parametric, univariate analysis of the effect of categorical variable
# since all colonies were not dead, using 'censoring': absence of event (death) due to any other factor such as experiment stopped at a time

# STEP 1: Creating survival object for worker group as individual ----
#(estimating the differences in proportion of colonies that are surviving)
# right censored data = event (death) is not observed over the course of experiment
# death = 1 if death observed (workers = 0 in colony). Death = 2 is censored when death not observed by the time
# avoid 0 since that means 'no observed value'

# create a survival object with time and event
data <- wk.time.death

srv.fit <- survfit(Surv(worker.age.death, death)~wolbachia, data = data)
srv.fit #extracting median surival time for Wolbachia infected and uninfected

# model summary
summary(srv.fit) #detailed
df.sum.srvfit <- surv_summary(srv.fit,
                              data = data) #data frame with summary survfit results (survival curve)

# P-value for difference in survival time: log-rank test
srv.diff <- survdiff(Surv(worker.age.death, death)~wolbachia, data = data) #p < 0.0001
summary(srv.diff)

pairwise_survdiff(Surv(worker.age.death, death) ~ wolbachia, data = data,
                  p.adjust.method = "bonferroni",
                  rho = 0) #0 = log-rank test

# saving survival curve output into csv file
write.csv(df.sum.srvfit, #data frame
          file.path(out.dir,# file path
                    "worker-survival-curve-summary.csv"), #file name
          row.names=T)

# PLOT: Survival Curves (KP) ---- 
ggsurv.wk <- ggsurvplot(srv.fit, # survfit object with calculated statistics
                     data = data,# data used to fit survival curves.
                     pval = TRUE,# show p-value of log-rank test.
                     conf.int = T,# show confidence intervals for point estimates of survival curves.
                     censor.shape="|",
                     censor.size = 4,
                     xlim = c(0,96),# customize X axis, but not affect survival estimates.
                     xlab = "Expected age of workers (days)",   # customize X axis label.
                     break.time.by = 6,     # break X axis in time intervals
                     risk.table = TRUE,# show risk table.
                     #risk.table.fontsize = 8, #font size to be used for the risk table.
                     risk.table.y.text.col = T,# colour risk table text annotations.
                     risk.table.height = 0.20, # the height of the risk table
                     risk.table.y.text = FALSE,# show bars instead of names in text annotations in legend of risk table.
                     ncensor.plot = TRUE,      # plot the number of censored subjects at time t
                     ncensor.plot.height = 0.25,
                     surv.median.line = "hv",  # add the median survival pointer.
                     legend.labs=c("Infected", "Uninfected"),
                     legend.title="Wolbachia", # change legend labels
                     palette=c("darkorchid4", "chocolate2"), # color for Wolbachia
                     ggtheme = theme(plot.background=element_rect(fill = NA),
                                     panel.background=element_rect(fill = NA),
                                     panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     legend.text = element_text(family = "Arial",colour = "black", size = "10"),
                                     legend.background = element_rect(fill = NA),
                                     legend.key = element_rect(color = NA, fill = NA))) # customize plot and risk table with a theme.

ggsurv.wk

# customizing survival plot
ggsurv.wk$plot <- ggpar(ggsurv.wk$plot,
                        title    = "Survival curves of Groups of Workers",
                        subtitle = "(Kaplan-Meier estimates)",
                        font.title    = c(14, "bold"),  #c(font size, text format, text color)
                        font.subtitle = c(12, "italic"),
                        ggtheme = theme (
                          axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "14"),
                          axis.text = element_text(family = "Arial", colour = "black", size = "11"),
                          axis.ticks = element_line(colour = "black", size = 1)) # customize plot and risk table with a theme.
                        )

# customizing risk table
ggsurv.wk$table <- ggpar(ggsurv.wk$table,
                        title    = "Risk Set",
                        font.title    = c(12, "bold"),  #c(font size, text format, text color)
                        ggtheme = theme (
                          axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
                          axis.text = element_text(family = "Arial", colour = "black", size = "10")))

# customizing survival plot - although no need for workers since all the groups had died by the end of the 
ggsurv.wk$ncensor.plot <- ggpar(ggsurv.wk$ncensor.plot,
                        title    = "Number of censoring",
                        font.title    = c(12, "bold"),  #c(font size, text format, text color)
                        ylab = "censoring \\nevents",
                        ggtheme = theme (
                          axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
                          axis.text = element_text(family = "Arial", colour = "black", size = "10"),
                          axis.ticks = element_line(colour = "black", size = 1)) # customize plot and risk table with a theme.
)

# looking at the edited plots
ggsurv.wk

#saving image in desired format (individual plots or all together)
ggsave("worker-survivalcurve-risktable.jpg", 
       plot = print(ggsurv.wk$table),
       path = gr.wrk,
       width = 5, 
       height = 1.5,
       dpi = 300)

# STEP 2: COX Model for survival of worker group ----
# step 1: build cox model
# step 2: diagnostics of model
# step 3: estimated survival curves


# when studying effect of more than 1 factor (include covariates). Same structure as GLM/LM
# model assumptions: proportional hazards, non-informative censoring
# non-parametric, multivariante analysis, catergorical (wolbachia) and quantitative predictor (counts) vairable

# STEP 1: HYPOTHESIS TESTING ---
#(estimating the differences in proportion of groups that are surviving)
data = wk.time.death
srv.fit <- survfit(Surv(worker.age.death, death) ~ wolbachia, data = data) #survival object

wk.srv.cox<- coxph(Surv(worker.age.death, death) ~ wolbachia, data = data) #cox model

# model summary
summary(wk.srv.cox) #hazard ratio and p-value
termplot(wk.srv.cox, terms = "wolbachia") #plot to show difference in hazard ratios

#tabular summary of hazard ratio and p-value
wk.hr <- wk.srv.cox %>%
  gtsummary::tbl_regression(exp = TRUE) 

# quoting log(HR) in paper
wk.hr.log <- wk.srv.cox %>%
  tbl_regression()

gt::gtsave(as_gt(wk.hr.log), 
           file = file.path(gr.wrk, "WkSurv-Cox-LogHazardRatio.png"))

#plot summary of hazard ratio with p-value
ggforest(wk.srv.cox, data = data)

# STEP 3: MODEL DIAGNOSTICS -----
# (https://bookdown.org/sestelo/sa_financial/how-to-evaluate-the-ph-assumption.html)

# LOG-LOG SURVIVAL CURVE: transformation of estimated survival curve that results from taking the natural log of an estimated survival probability twice
# log-log survial curve: 2 curves (grouping variable) should be parallel
# distance b/w 2 curve = linear expression of differences in predictor values which is doesn't involve time
wrk.cox <- ggsurvplot(srv.fit, data = data,
           fun = "cloglog",
           xlim=c(25, 100),
           title = "Diagnostic Plot for Worker Survival\\nLog-log survival by Wolbachia",
           subtitle = "should be parallel",
           ylab = "log-log survival" ,
           xlab = "Worker age (days) using log",
           legend.labs=c("Infected", "Uninfected"),
           legend.title="Wolbachia", # change legend labels
           palette=c("darkorchid4", "chocolate2"), # color for Wolbachia
           ggtheme = theme(plot.background=element_rect(fill = NA),
                           panel.background=element_rect(fill = NA),
                           panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           legend.text = element_text(family = "Arial",colour = "black", size = "10"),
                           legend.background = element_rect(fill = NA),
                           legend.key = element_rect(color = NA, fill = NA),
                           axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
                           axis.text = element_text(family = "Arial", colour = "black", size = "12"),
                           axis.ticks = element_line(colour = "black", size = 1)))

#saving image in desired format (individual plots or all together)
ggsave("worker-coxdiagnostics-loglogsurvival.jpg", 
       plot = print(wrk.cox),
       path = gr.wrk,
       width = 4.5, 
       height = 4,
       dpi = 300)

# GOODNESS-OF-FIT of Cox model (Testing for Proportional Hazards, p > 0.05)
# correlation (rho) b/w Schoenfeld residuals and survival time should be ~0 = good model (met model assumptions)
test.cox.ph <- cox.zph(wk.srv.cox)
test.cox.ph

# DEVIANCE RESIDUALS: check outliers.
# deviance residual is a normalized transform of the martingale residual.
# should be roughly symmetrically distributed about zero with a standard deviation of 1.
cox.dev<- ggcoxdiagnostics(wk.srv.cox,
                 type = "deviance", #change this to "deviance"
                 title = "Diagnostic Plot for Worker Survival: Deviance",
                 subtitle = "Cox-Propotional Hazard (Deviance): should be no trend",
                 linear.predictions = F,
                 font.subtitle = 9,
                 ggtheme = theme_bw())

# SCHOENFELD RESIDUALS : should be parallel to X-axis, flat and centered around 0
# SR represent the difference between the observed covariateS and the expected value given the risk set at that time
cox.sch <- ggcoxzph(test.cox.ph,
         title = " Worker Survival: Schoenfeld Residuals",
         subtitle = "Cox-Propotional Hazard (Schoenfeld Residuals): should be no trend",
         font.subtitle = 9,
         ggtheme = theme(plot.background=element_rect(fill = NA),
                         panel.background=element_rect(fill = NA),
                         panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         legend.text = element_text(family = "Arial",colour = "black", size = "10"),
                         legend.background = element_rect(fill = NA),
                         legend.key = element_rect(color = NA, fill = NA),
                         axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
                         axis.text = element_text(family = "Arial", colour = "black", size = "12"),
                         axis.ticks = element_line(colour = "black", size = 1)))

#saving image in desired format (individual plots or all together)
ggsave("WrkSurv-coxdiagnostics-deviance.jpg", 
       plot = print(cox.dev),
       path = gr.wrk,
       width = 4.5, 
       height = 4,
       dpi = 300)

# STEP 4: ESTIMATED SURVIVAL CURVES -----
wrk.cox2<- coxph(Surv(worker.age.death, death) ~ wolbachia, data = data)
newdf.wrk <- data.frame(wolbachia = levels(data$wolbachia), 2)

wrk.new.fit <- survfit(wrk.cox2, newdata = newdf.wrk)

#summary(fit) # to see the estimated values
ggsurv.wk <- ggsurvplot(wrk.new.fit, # survfit object with calculated statistics
                        data = newdf.wrk,# data used to fit survival curves
                        title    = "Estimated Survival curves for Worker Groups",
                        subtitle = "(Cox-Proportional Hazards Test)",
                        font.title    = c(12, "bold"),  #c(font size, text format, text color)
                        font.subtitle = c(12, "italic"),
                        pval = TRUE,# show p-value of log-rank test.
                        conf.int = T,# show confidence intervals for point estimates of survival curves.
                        censor.shape="|",
                        censor.size = 4,
                        xlim = c(0,96),# customize X axis, but not affect survival estimates.
                        xlab = "Time in days",   # customize X axis label.
                        break.time.by = 6,     # break X axis in time intervals
                        surv.median.line = "hv",  # add the median survival pointer.
                        legend.labs=c("Infected", "Uninfected"),
                        legend.title="Wolbachia", # change legend labels
                        palette=c("darkorchid4", "chocolate2"), # color for Wolbachia
                        ggtheme = theme(plot.background=element_rect(fill = NA),
                                        panel.background=element_rect(fill = NA),
                                        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        legend.text = element_text(family = "Arial",colour = "black", size = "10"),
                                        legend.background = element_rect(fill = NA),
                                        legend.key = element_rect(color = NA, fill = NA),
                                        axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
                                        axis.text = element_text(family = "Arial", colour = "black", size = "10"),
                                        axis.ticks = element_line(colour = "black", size = 1))) # customize plot and risk table with a theme.

ggsurv.wk


# GLMM for individual-level survival ----
data <- df.wk

data$group <- as.factor(data$group)

# proportion alive (wrt to start)
wk.srv.glm <- glmer(cbind(wk.alive, wk.dead.overall) ~ wolbachia + worker.age + (1|group), data = data,
                    family=binomial(link=logit), nAGQ = 0,
                    control=glmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=100000))) #using this model in paper

wk.srv.glm <- glm(cbind(wk.alive, wk.dead.overall) ~ wolbachia + worker.age + wolbachia:worker.age, data = data,
                    family=quasibinomial(link=logit))


model = wk.srv.glm
vif(model)
summary(model)
stats.glmer(model) 
stats.glm(model)
qq(model)
res(model)
shp.tst(model)
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# proportion alive (wrt to previous time point)
wk.srv.glm <- glmer(cbind(wk.alive, wk.dead) ~ wolbachia + worker.age + (1|colony), data = df.wk,
                    family=binomial(link=logit),
                    control=glmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=100000))) #using this model in paper

wk.srv.glm <- glm(cbind(wk.alive, wk.dead) ~ wolbachia + worker.age + wolbachia:worker.age, data = df.wk,
                  family=quasibinomial(link=logit))


model = wk.srv.glm
summary(model)
stats.glmer(model) 
stats.glm(model)
qq(model)
res(model)
shp.tst(model)
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# -------------------------- QUEENS: SURVIVAL ANALYSIS AND COLONY DYNAMICS -------
# loading processed file for analysis
setwd(out.dir)
df.qn <- read.csv("queen2020-survival.csv")
qn.time.event< - read.csv("queen2019-event-time.csv")

# loading raw files and generating variables
setwd(raw.dir)
list.files() #visualize files in the folder

df.qn <- read.csv("queen2020-survival.csv") #long-format for queen survival

# creating variables
df.qn <- df.qn %>%
  dplyr::mutate(young.brood = eggs + young.larvae, # number of eggs and young larvae
                older.brood = old.larvae + worker.pup, # number of late-instar larvae and worker pupae
                total.brood = young.brood + older.brood + worker.pre.pup, # total young and older brood
                adults = queens + worker, # number of queens and workers
                col.size = total.brood + adults,
                yb.prp = round(young.brood/col.size, digits = 3), # proportion of young.brood in the whole colony
                oldb.prp = round(older.brood/col.size, digits = 3) , # proportion of older.brood in the whole colony
                adlt.prp = round(adults/col.size, digits = 3)) # proportion of adults in the whole colony

# change in queen number wrt to first census (overall) and wrt to previous time point
df.qn <- df.qn %>%
  dplyr::group_by(colony) %>%
  dplyr::arrange(week, .by_group = TRUE) %>%
  dplyr::mutate(qn.dead = queens - dplyr::lag(queens, # lag = previous row
                                             n = 1, # 1 row difference
                                             default = first(queens)), # dead queens compared to previous time point
                qn.dead = abs(qn.dead), #absolute values, changing above -ve to +ve
                qn.dead.prp = round(qn.dead/(queens+qn.dead), digits = 3),# proportion of dead queens at each time point wrt previous time point
                qn.live.prp = round(queens/(queens+qn.dead), digits = 3),# proportion of live queens at each time point
                qn.dead.overall = first(queens)-queens, # dead queens wrt start (t=0)
                total.queen = qn.dead.overall + queens, # confirming that numbers match up
                qn.surv.prp.overall = round(queens/first(queens), digits = 3)) %>%# proportion of live queens wrt start (t=0)
  dplyr::rename(qn.alive = queens)

 

#testing if it worked
df <- data.frame(subset(df.qn, colony == 12|colony == 27|colony == 46))

df <- data.frame(subset(df.wk1, death == 2))

# PLOT: trend of overall proportion of dead workers over time
data <- df.wk
data$x <- data$worker.age
data$y <- data$wrk.surv.prp.overall

ggplot(data, aes(x=x, y=y, group = groupID,
                 color = wolbachia, fill = wolbachia)) +
  geom_line()+
  facet_wrap(~groupID)



# Raw file: variables ----
setwd(raw.dir)
list.files() #visualize files in the folder

df.qn <- read.csv("queen2020-survival.csv")  #long-format for queen survival
df.qn$wolbachia <- as.factor(df.qn$wolbachia)
df.qn$colony <- as.factor(df.qn$colony)

# change in queen number wrt to first census (overall) and wrt to previous time point
df.qn <- df.qn %>%
  dplyr::group_by(colony) %>%
  dplyr::arrange(week, .by_group = TRUE) %>%
  dplyr::mutate(qn.dead = dplyr::lag(queens, # lag = previous row
                                     n = 1, # 1 row difference
                                     default = first(queens))-queens, # dead queens compared to previous time point
                #qn.dead = abs(qn.dead), #making negative values as positive
                qn.dead.overall = first(queens)-queens, # number of dead queens since start
                total.queens = qn.dead.overall + queens, # confirming that the numbers match up# dead queens compared to starting
                qn.surv.prp.overall = round((queens/first(queens)), digits = 3),# survival proportion wrt start
                qn.dead.prp = round((qn.dead/(queens+qn.dead)), digits = 3)) %>% # death proportion wrt previous time point
  dplyr::rename(qn.alive = queens)

# assigning death status
# Some alive queens (death = 1, no, censored); No queens (death = 2, yes); New queen pupae (death = 2, reproduction) 
# censor 1, dead or reproduction = 2 
df.qn <- df.qn %>%
  dplyr::mutate(death = ifelse(qn.alive == 0, 2, ifelse(queen.pup > 0, 2, 1)))

#testing if it worked
df <- data.frame(subset(df.qn, colony == "33"|colony == "34"|colony == "38",
                 select = c(colony, wolbachia, week, qn.alive, qn.dead:qn.dead.prp)))

df <- data.frame(subset(df.qn, death != 1,
                 select = c(colony:week, qn.alive, queen.pup, qn.dead,death)))

# PLOT: trend of overall proportion of dead queens over time
data <- df.qn
data$x <- data$queen.age
data$y <- data$qn.alive

ggplot(data, aes(x=x, y=y, group = colony,
                 color = wolbachia, fill = wolbachia)) +
  geom_line()+
  facet_wrap(~colony)

# Time till death of the group: data file -----
# colonies which are dead or have produced new queens (death == 2) --
qn.time.death <- df.qn %>%
  dplyr::group_by(colony) %>% # group by experimental group ID
  dplyr::arrange(week, .by_group = T ) %>% # arrange in increasing order of the day of census
  dplyr::filter(qn.dead.prp == 1 | queen.pup > 0) %>% # first time point at which the group was dead. Subsequent time points after death have 'NA' values
  dplyr::slice(1) #first time the event of death happens

# Does the time till death of the group differ between infected and uninfected?
boxplot(queen.age~wolbachia, data = qn.time.death, notch = F)

# colonies that did not die till the end and are censored (death == 1) --
qn.censor <- df.qn %>%
  dplyr::group_by(colony) %>% # group by experimental group ID
  dplyr::arrange(week, .by_group = T ) %>% # arrange in increasing order of the day of census
  dplyr::slice(n()) %>% #last time point of census
  dplyr::filter(death == 1) # alive colonies at the last time point

# merging the two dataframes to know time to event or censor
qn.time.event <- rbind(qn.time.death, qn.censor)

# queen : saving this processed file ----
# saving csv file
write.csv(df.qn, #data frame
          file.path(out.dir,# file path
                    "queen2019-survival-census.csv"), #file name
          row.names=T)

write.csv(qn.censor, #data frame
          file.path(out.dir,# file path
                    "queen2019-censor-time.csv"), #file name
          row.names=T)

write.csv(qn.time.death, #data frame
          file.path(out.dir,# file path
                    "queen2019-survival-time.csv"), #file name
          row.names=T)

write.csv(qn.time.event, #data frame
          file.path(out.dir,# file path
                    "queen2019-event-time.csv"), #file name
          row.names=T)

# Kaplan-Meir (KP) method for fitted survival ----------
# first pass survival probability estimate
# non-parametric, univariate analysis of the effect of categorical variable
# since all colonies were not dead, using 'censoring': absence of event (death) due to any other factor such as experiment stopped at a time

# STEP 1: Creating survival object for queen group as individual ----
#(estimating the differences in proportion of colonies that are surviving)
# right censored data = event (death) is not observed over the course of experiment
# death = 1 if death observed (queens = 0 in colony). Death = 2 is censored when death not observed by the time
# avoid 0 since that means 'no observed value'

# create a survival object with time and event
data <- qn.time.event

srv.fit <- survfit(Surv(queen.age, death)~wolbachia, data = data)
srv.fit #extracting median surival time for Wolbachia infected and uninfected

# model summary
summary(srv.fit) #detailed
df.sum.srvfit <- surv_summary(srv.fit,
                              data = data) #data frame with summary survfit results (survival curve)

# P-value for difference in survival time: log-rank test
srv.diff <- survdiff(Surv(queen.age, death)~wolbachia, data = data) #p < 0.0001
summary(srv.diff)

pairwise_survdiff(Surv(queen.age, death) ~ wolbachia, data = data,
                  p.adjust.method = "bonferroni",
                  rho = 0) #0 = log-rank test

# saving survival curve output into csv file
write.csv(df.sum.srvfit, #data frame
          file.path(out.dir,# file path
                    "queen-survival-curve-summary.csv"), #file name
          row.names=T)

# PLOT: Survival Curves (KP) ---- 
ggsurv.qn <- ggsurvplot(srv.fit, # survfit object with calculated statistics
                        data = data,# data used to fit survival curves.
                        pval = TRUE,# show p-value of log-rank test.
                        conf.int = T,# show confidence intervals for point estimates of survival curves.
                        censor.shape="|",
                        censor.size = 4,
                        xlim = c(80,300),# customize X axis, but not affect survival estimates.
                        xlab = "Expected age of queens (days)",   # customize X axis label.
                        break.time.by = 20,     # break X axis in time intervals
                        risk.table = TRUE,# show risk table.
                        #risk.table.fontsize = 8, #font size to be used for the risk table.
                        risk.table.y.text.col = T,# colour risk table text annotations.
                        risk.table.height = 0.20, # the height of the risk table
                        risk.table.y.text = FALSE,# show bars instead of names in text annotations in legend of risk table.
                        ncensor.plot = TRUE,      # plot the number of censored subjects at time t
                        ncensor.plot.height = 0.25,
                        surv.median.line = "hv",  # add the median survival pointer.
                        legend.labs=c("Infected", "Uninfected"),
                        legend.title="Wolbachia", # change legend labels
                        palette=c("darkorchid4", "chocolate2"), # color for Wolbachia
                        ggtheme = theme(plot.background=element_rect(fill = NA),
                                        panel.background=element_rect(fill = NA),
                                        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        legend.text = element_text(family = "Arial",colour = "black", size = "10"),
                                        legend.background = element_rect(fill = NA),
                                        legend.key = element_rect(color = NA, fill = NA))) # customize plot and risk table with a theme.

ggsurv.qn

# customizing survival plot
ggsurv.qn$plot <- ggpar(ggsurv.qn$plot,
                        title    = "Survival curves of Groups of queens",
                        subtitle = "(Kaplan-Meier estimates)",
                        font.title    = c(14, "bold"),  #c(font size, text format, text color)
                        font.subtitle = c(12, "italic"),
                        pval = TRUE,# show p-value of log-rank test.
                        ggtheme = theme (
                          axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "14"),
                          axis.text = element_text(family = "Arial", colour = "black", size = "11"),
                          axis.ticks = element_line(colour = "black", size = 1)) # customize plot and risk table with a theme.
)

# customizing risk table
ggsurv.qn$table <- ggpar(ggsurv.qn$table,
                         title    = "Risk Set",
                         font.title    = c(12, "bold"),  #c(font size, text format, text color)
                         ggtheme = theme (
                           axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
                           axis.text = element_text(family = "Arial", colour = "black", size = "10")))

# customizing survival plot - although no need for queens since all the groups had died by the end of the 
ggsurv.qn$ncensor.plot <- ggpar(ggsurv.qn$ncensor.plot,
                                title    = "Number of censoring",
                                font.title    = c(12, "bold"),  #c(font size, text format, text color)
                                ylab = "censoring \\nevents",
                                ggtheme = theme (
                                  axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
                                  axis.text = element_text(family = "Arial", colour = "black", size = "10"),
                                  axis.ticks = element_line(colour = "black", size = 1)) # customize plot and risk table with a theme.
)

# looking at the edited plots
ggsurv.qn

#saving image in desired format (individual plots or all together)
ggsave("queen-survivalcurve-kp.jpg", 
       plot = print(ggsurv.qn$plot),
       path = gr.qn,
       width = 5, 
       height = 4,
       dpi = 300)

# STEP 2: COX Model for survival of queen group ----
# step 1: build cox model
# step 2: diagnostics of model
# step 3: estimated survival curves


# when studying effect of more than 1 factor (include covariates). Same structure as GLM/LM
# model assumptions: proportional hazards, non-informative censoring
# non-parametric, multivariante analysis, catergorical (wolbachia) and quantitative predictor (counts) vairable

# STEP 1: HYPOTHESIS TESTING ---
#(estimating the differences in proportion of groups that are surviving)
data = qn.time.event
srv.fit <- survfit(Surv(queen.age, death) ~ wolbachia , data = data) #survival object

qn.srv.cox<- coxph(Surv(queen.age, death) ~ wolbachia , data = data) #cox model

# model summary
summary(qn.srv.cox) #hazard ratio and p-value
termplot(qn.srv.cox, terms = "wolbachia") #plot to show difference in hazard ratios

#tabular summary of hazard ratio and p-value
qn.hr <- qn.srv.cox %>%
  gtsummary::tbl_regression(exp = TRUE) 

# quoting log(HR) in paper
qn.hr.log <- qn.srv.cox %>%
  tbl_regression()

gt::gtsave(as_gt(qn.hr), 
           file = file.path(gr.qn, "qnSurv-Cox-HazardRatio.png"))

#plot summary of hazard ratio with p-value
ggforest(qn.srv.cox, data = data)

# STEP 3: MODEL DIAGNOSTICS -----
# (https://bookdown.org/sestelo/sa_financial/how-to-evaluate-the-ph-assumption.html)

# LOG-LOG SURVIVAL CURVE: transformation of estimated survival curve that results from taking the natural log of an estimated survival probability twice
# log-log survial curve: 2 curves (grouping variable) should be parallel
# distance b/w 2 curve = linear expression of differences in predictor values which is doesn't involve time
qn.cox <- ggsurvplot(srv.fit, data = data,
                      fun = "cloglog",
                      xlim=c(120, 300),
                      title = "Diagnostic Plot for queen Survival\\nLog-log survival by Wolbachia",
                      subtitle = "should be parallel",
                      ylab = "log-log survival" ,
                      xlab = "queen age (days) using log",
                      legend.labs=c("Infected", "Uninfected"),
                      legend.title="Wolbachia", # change legend labels
                      palette=c("darkorchid4", "chocolate2"), # color for Wolbachia
                      ggtheme = theme(plot.background=element_rect(fill = NA),
                                      panel.background=element_rect(fill = NA),
                                      panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      legend.text = element_text(family = "Arial",colour = "black", size = "10"),
                                      legend.background = element_rect(fill = NA),
                                      legend.key = element_rect(color = NA, fill = NA),
                                      axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
                                      axis.text = element_text(family = "Arial", colour = "black", size = "12"),
                                      axis.ticks = element_line(colour = "black", size = 1)))

#saving image in desired format (individual plots or all together)
ggsave("queen-coxdiagnostics-loglogsurvival.jpg", 
       plot = print(qn.cox),
       path = gr.qn,
       width = 4.5, 
       height = 4,
       dpi = 300)

# GOODNESS-OF-FIT of Cox model (Testing for Proportional Hazards, p > 0.05)
# correlation (rho) b/w Schoenfeld residuals and survival time should be ~0 = good model (met model assumptions)
test.cox.ph <- cox.zph(qn.srv.cox)
test.cox.ph

# DEVIANCE RESIDUALS: check outliers.
# deviance residual is a normalized transform of the martingale residual.
# should be roughly symmetrically distributed about zero with a standard deviation of 1.
cox.dev<- ggcoxdiagnostics(qn.srv.cox,
                           type = "deviance", #change this to "deviance"
                           title = "Diagnostic Plot for queen Survival: Deviance",
                           subtitle = "Cox-Propotional Hazard (Deviance): should be no trend",
                           linear.predictions = F,
                           font.subtitle = 9,
                           ggtheme = theme_bw())

# SCHOENFELD RESIDUALS : should be parallel to X-axis, flat and centered around 0
# SR represent the difference between the observed covariateS and the expected value given the risk set at that time
cox.sch <- ggcoxzph(test.cox.ph,
                    title = " queen Survival: Schoenfeld Residuals",
                    subtitle = "Cox-Propotional Hazard (Schoenfeld Residuals): should be no trend",
                    font.subtitle = 9,
                    ggtheme = theme(plot.background=element_rect(fill = NA),
                                    panel.background=element_rect(fill = NA),
                                    panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    legend.text = element_text(family = "Arial",colour = "black", size = "10"),
                                    legend.background = element_rect(fill = NA),
                                    legend.key = element_rect(color = NA, fill = NA),
                                    axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
                                    axis.text = element_text(family = "Arial", colour = "black", size = "12"),
                                    axis.ticks = element_line(colour = "black", size = 1)))

#saving image in desired format (individual plots or all together)
ggsave("qnSurv-coxdiagnostics-schoenfeld.jpg", 
       plot = print(cox.sch),
       path = gr.qn,
       width = 4.5, 
       height = 4,
       dpi = 300)

# STEP 4: ESTIMATED SURVIVAL CURVES -----
qn.cox2<- coxph(Surv(queen.age, death) ~ wolbachia, data = data)
newdf.qn <- data.frame(wolbachia = levels(data$wolbachia), 2)

qn.new.fit <- survfit(qn.cox2, newdata = newdf.qn)

#summary(fit) # to see the estimated values
ggsurv.qn <- ggsurvplot(qn.new.fit, # survfit object with calculated statistics
                        data = newdf.qn,# data used to fit survival curves
                        title    = "Estimated Survival curves for queen Groups",
                        subtitle = "(Cox-Proportional Hazards Test)",
                        font.title    = c(12, "bold"),  #c(font size, text format, text color)
                        font.subtitle = c(12, "italic"),
                        pval = TRUE,# show p-value of log-rank test.
                        conf.int = T,# show confidence intervals for point estimates of survival curves.
                        censor.shape="|",
                        censor.size = 4,
                        xlim = c(80,300),# customize X axis, but not affect survival estimates.
                        xlab = "Time in days",   # customize X axis label.
                        break.time.by = 20,     # break X axis in time intervals
                        surv.median.line = "hv",  # add the median survival pointer.
                        legend.labs=c("Infected", "Uninfected"),
                        legend.title="Wolbachia", # change legend labels
                        palette=c("darkorchid4", "chocolate2"), # color for Wolbachia
                        ggtheme = theme(plot.background=element_rect(fill = NA),
                                        panel.background=element_rect(fill = NA),
                                        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        legend.text = element_text(family = "Arial",colour = "black", size = "10"),
                                        legend.background = element_rect(fill = NA),
                                        legend.key = element_rect(color = NA, fill = NA),
                                        axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "11"),
                                        axis.text = element_text(family = "Arial", colour = "black", size = "10"),
                                        axis.ticks = element_line(colour = "black", size = 1))) # customize plot and risk table with a theme.

ggsurv.qn

#saving image in desired format (individual plots or all together)
ggsave("queen-survivalcurve-cox.jpg", 
       plot = print(ggsurv.qn),
       path = gr.qn,
       width = 5, 
       height = 4,
       dpi = 300)


# GLMM for individual-level survival ----
data <- subset(df.qn, queen.pup <1)

data$colony <- as.factor(data$colony)

# proportion alive (wrt to start)
qn.srv.glm<- glmer(cbind(qn.alive, qn.dead.overall) ~ wolbachia + queen.age + (1|colony), data = data,
                    family=binomial(link=logit), nAGQ = 0,
                    control=glmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=100000))) #using this model in paper

qn.srv.glm <- glm(cbind(qn.alive, qn.dead.overall) ~ wolbachia + queen.age + wolbachia:queen.age, data = data,
                  family=quasibinomial(link=logit))


model = qn.srv.glm
vif(model)
summary(model)
stats.glmer(model) 
stats.glm(model)
qq(model)
res(model)
shp.tst(model)
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# ------------------------- PLOTS: RAW AND MEAN VALUES ---------------
# Group-level dynamics along with Mean --------------
my.color <-c ("darkorchid4","chocolate2" )

# Step 1: assign a variable and labels for the plot
data <- subset(df.qn, queen.pup <1)
data$colony <- as.factor(data$colony)

data$x <- data$queen.age #x-axis
data$y <- data$qn.surv.prp.overall #y-axis
title <- "Proportion of Alive Queens"
y.lab <- "proportion of alive queens \\n(since start)"
x.lab <- "expected age of queens (days)"

# Step 2: define common variable for colonyID-level data
data.mean <- summarySE(data, measurevar="y",
                       groupvars=c("wolbachia", "x"), na.rm = T)

# Step 3: reorder variablesa 
data$wolbachia<- factor(data$wolbachia, 
                        levels = c("infected", "uninfected"))

data.mean$wolbachia<- factor(data.mean$wolbachia, 
                             levels = c("infected", "uninfected"))

# Step 4: Plot = sample values overlaid with mean and 95% CI. Texts were later edited in Inkscape
ggplot(data, aes(x=x,
                 y = y,
                 group = colony)) + 
  scale_colour_manual(values = my.color)+
  scale_fill_manual(values = my.color)+
  geom_line(aes(color = wolbachia),
            size = 0.8,
            alpha = 0.2)+
  geom_errorbar(data = data.mean,
                aes(ymin= y-ci, ymax= y+ci,
                    colour = wolbachia,
                    group = wolbachia),
                #linetype = "longdash", 
                size = 0.8, width= 0.1)+
  geom_line(data = data.mean,
            aes(group = wolbachia,
                color = wolbachia),
            size = 2,
            position= position_dodge(0.2))+
  geom_point(data = data.mean,
             aes(fill = wolbachia,
                 group = wolbachia), shape = 21,
             size = 4)+
  #facet_wrap(~wolbachia)+
  ggtitle(title)+
  #scale_x_discrete(breaks = c("3", "12", "21", "30", "39", "48", "57", "66","75", "93"))+
  labs(x = x.lab,
       y = y.lab) +
  theme(strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(family = "Arial", colour = "black",
                                  face = "italic",size = "14"),
        plot.title = element_text(hjust = 0.5,
                                  family = "Arial", colour = "black",
                                  face = "bold", size = "18"),
        panel.background=element_rect(fill = NA),
        panel.spacing.x=unit(0.12, "in"),
        plot.background=element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        #legend.position = c(0.08,0.80),
        legend.text = element_text(family = "Arial", colour = "black", size = "12"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "16"),
        axis.text = element_text(family = "Arial", colour = "black", size = "14"),
        axis.ticks = element_line(colour = "black", size = 1))

#saving image in desired format
ggsave(file = "QnSurv-SurvProp-grouplevel.png",
       path = gr.wrk,
       dpi = 300, 
       width = 6.5,
       height = 4.5, 
       units = "in")

# PLOT : Box plot over time with raw values ---------
my.color <-c ("darkorchid4", "chocolate2")

#Step 1: create universal variable and 
data <- df.wk
data$y <- data$wrk.surv.prp.overall
data$x <- as.factor(data$worker.age)
y.lab <- "Proportion of Surviving Workers"
x.lab <- "Age of the workers (days) "
title <- "Proportion of Workers Surviving within Groups"

# Step 2: re-order
data$wolbachia <- factor(data$wolbachia, 
                         levels = c("infected", "uninfected"))

# Step 3: plot
ggplot(data,
       aes(x=x,
           y = y,
           fill = wolbachia)) + 
  geom_boxplot(notch = F,
               position= position_dodge(0.7),
               color = "black",
               alpha = 0.3,
               width = 0.5)+
  stat_boxplot(geom = "errorbar",
               position= position_dodge(0.7),
               width = 0.3)+ #width of the error bar
  geom_jitter(stat = "identity" ,
              shape = 21,
              size = 1.5,
              alpha = 0.6,
              position = position_dodge(0.7))+ 
  #facet_wrap(~q.age.month) +
  scale_fill_manual(values = my.color)+
  scale_color_manual(values = my.color)+
  labs(title =title,
       x = x.lab,
       y = y.lab)+
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
        legend.text = element_text(family = "Arial",
                                   colour = "black",
                                   size = "12"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_text(family = "Arial",
                                  colour ="black",
                                  face = "bold", size = "16"),
        axis.text = element_text(family = "Arial",
                                 colour = "black", size = "12"),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.background=element_rect(fill = NA),
        plot.title = element_text(family = "Arial", colour = "black",
                                  face = "bold", size = "18", 
                                  margin = margin(t='0.5', r='0.5', b='0.5', l='0.5', unit = 'cm')))


#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(file = "WkSurv-SurvProp-boxplot.png",
       path = gr.wrk,
       dpi = 300, 
       width = 7,
       height = 4, 
       units = "in")


# PLOT : Box plot per Wolbachia group with raw values ---------
my.color <-c ("darkorchid4", "chocolate2")

#Step 1: create universal variable and 
data <- subset(qn.time.death)

data <- wk.time.death
data$y <- data$queen.age
data$x <- data$wolbachia
y.lab <- "Queen age (days)"
title <- "Age of the Queens at the \\nTime of Death of the Group"

# Step 2: re-order
data$wolbachia <- factor(data$wolbachia, 
                         levels = c("infected", "uninfected"))

# Step 3: Mean
data.summary <- data %>%
  dplyr::group_by(x) %>%
  dplyr::summarise(mean = mean(y),
                   median = median(y),
                   sd = sd(y),
                   sample.size = n())
# Step 4: plot
ggplot(data,
       aes(x=x,
           y = y,
           fill = wolbachia)) + 
  geom_boxplot(notch = F,
               position= position_dodge(0.7),
               color = "black",
               alpha = 0.3,
               width = 0.5)+
  stat_boxplot(geom = "errorbar",
               position= position_dodge(0.7),
               width = 0.3)+ #width of the error bar
  geom_jitter(stat = "identity" ,
              shape = 21,
              size = 2.5,
              alpha = 0.8,
              width = 0.1)+ #width of the jitter 
  stat_summary(fun = mean, 
               geom = "point",
               shape = 24, 
               size = 3, 
               color = "black",
               fill = "black") +
  scale_fill_manual(values = my.color)+
  scale_color_manual(values = my.color)+
  scale_y_continuous(breaks = seq(110, 290, 20))+
  labs(title =title,
       y = y.lab)+
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
        legend.text = element_text(family = "Arial",
                                   colour = "black",
                                   size = "12"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_text(family = "Arial",
                                  colour ="black",
                                  face = "bold", size = "16"),
        axis.title.x = element_blank(),
        axis.text = element_text(family = "Arial",
                                 colour = "black", size = "12"),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.background=element_rect(fill = NA),
        plot.title = element_text(family = "Arial", colour = "black",
                                  face = "bold", size = "18", 
                                  margin = margin(t='0.5', r='0.5', b='0.5', l='0.5', unit = 'cm')))


#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(file = "QnSurv-time-death.png",
       path = gr.qn,
       dpi = 300, 
       width = 4,
       height = 3.5, 
       units = "in")


# PLOT: Individual colony data ----
my.color <-c ("darkorchid4","chocolate2" )

# Step 1: assign a variable and labels for the plot
data <- df.qn
data$x <- data$queen.age #x-axis
data$y <- data$qn.surv.prp.overall #y-axis
title <- "Proportion of Alive Queens"
y.lab <- "proportion of alive queens \\n(since start)"
x.lab <- "expected age of queens (days)"

# Step 2: reorder variablesa 
data$wolbachia<- factor(data$wolbachia, 
                         levels = c("infected",
                                    "uninfected"))

# Step 3: Plot = sample values overlaid with mean and 95% CI. Texts were later edited in Inkscape
ggplot(data, aes(x=x,
                  y = y)) + 
  scale_colour_manual(values = my.color)+
  geom_line(aes(color = wolbachia),
            size = 1)+
  facet_wrap(~colony)+ # producing each colony's graph
  ggtitle(title)+
  #scale_x_discrete(breaks = 4,limits = seq(1, 92))+ #y-axis ticks seq(lower limit, top limit, incremental steps
  labs(x = x.lab,
       y = y.lab) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", colour = "black",
                                  face = "italic",size = "12"),
        plot.title = element_text(hjust = 0.5, family = "Arial",
                                  colour = "black", face = "bold", size = "14"),
        plot.background=element_rect(fill = NA),
        panel.background=element_rect(fill = NA),
        panel.spacing.x=unit(0.1, "in"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        #legend.position = c(0.08,0.80),
        legend.text = element_text(family = "Arial", colour = "black", size = "10"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA))

#saving image in desired format
ggsave(file = "WkSurv-col-LiveWrk.png",
       path = gr.wk,
       dpi = 300, 
       width = 10,
       height = 8, 
       units = "in")



# PLOT: Mean value across coloneis --------------
my.color <-c ("darkorchid4","chocolate2" )

# Step 1: assign a variable and labels for the plot
data <- subset(df.qn, queen.pup <1)

data$colony <- as.factor(data$colony)
data$x <- data$queen.age #x-axis
data$y <- data$qn.surv.prp.overall #y-axis
title <- "Proportion of Alive Queens"
y.lab <- "proportion of alive queens \\n(since start)"
x.lab <- "expected age of queens (days)"

# Step 2: reorder variablesa 
data.mean <- summarySE(data, measurevar="y",
                       groupvars=c("wolbachia", "x"), na.rm = T)

# Step 3: reorder variablesa 
data.mean$wolbachia<- factor(data.mean$wolbachia, 
                             levels = c("infected",
                                        "uninfected"))

# Step 4: Plot = sample values overlaid with mean and 95% CI. Texts were later edited in Inkscape
ggplot(data.mean, aes(x=x,
                      y = y,
                      color = wolbachia,
                      group = wolbachia)) + 
  scale_colour_manual(values = my.color)+
  scale_fill_manual(values = my.color)+
  geom_errorbar(aes(ymin= y-ci, ymax= y+ci,
                    colour = wolbachia,
                    group = wolbachia),
                linetype = "longdash", 
                size = 0.8, width= 0.1,
                position= position_dodge(0.2))+
  geom_line(aes(group = wolbachia,
                color = wolbachia),
            size = 1,
            position= position_dodge(0.2))+
  geom_point(aes(fill = wolbachia,
                 group = wolbachia),shape = 21, 
             size = 4, position = position_dodge(0.2))+
  ggtitle(title)+
  #scale_x_discrete(limits = seq(70, 310, 30))+ #y-axis ticks seq(lower limit, top limit, incremental steps
  labs(x = x.lab,
       y = y.lab) +
  theme(strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(family = "Arial", colour = "black", face = "italic",size = "12"),
        plot.background=element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5,
                                  family = "Arial", colour = "black", face = "bold", size = "16"),
        panel.background=element_rect(fill = NA),
        panel.spacing.x=unit(0.1, "in"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(family = "Arial", colour ="black", face = "bold", size = "16"),
        axis.text = element_text(family = "Arial", colour = "black", size = "14"),
        axis.ticks = element_line(colour = "black", size = 1))

#saving image in desired format
ggsave(file = "QnSurv-SurvPrp-mean.png",
       path = gr.qn,
       dpi = 300, 
       width = 7,
       height = 4, 
       units = "in")
