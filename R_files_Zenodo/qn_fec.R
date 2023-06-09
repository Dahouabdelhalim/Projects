#GOAL: ARE THE EGG LAYING RATES DIFFERENT BETWEEN INFECTED AND UNINFECTED QUEENS?
# ---- instal and load packages ----
library("lme4") #LM, LMM, GLM, GLMM
library("car") #anova for LM and LMM
library("MASS")

library("multcomp")
library("dplyr")
library("rstatix") #outliers, annova
# Rmisc loads 'plyr' which masks the group_by function from dplyr
library("Rmisc") #stats (not correlation)
library("Hmisc") #correlation significance
library("permutes") #time series permutations

library("ggplot2") #for all graphs
library("ggfortify") #diagnostic plots of models

# shortcuts for directories ----
#directory with raw unprocessed files
raw.dir <- ("G:/My Drive/R/queen-fecundity/qn-fec-input")

#directory with processed files
out.dir <- ("G:/My Drive/R/queen-fecundity/qn-fec-output-files")
grph<- ("G:/My Drive/R/queen-fecundity/qn-fec-graphs")

# Loading processing files ----
setwd(out.dir)

df.fec <- read.csv("2019-fecundity-new-queen.csv") #newly eclosed and mated queens
df.ages <- read.csv("2017fecundity-queen-ages.csv") #egg laying by at different ages of the queen

# Raw Files: loading the file ----
setwd(raw.dir) #load the file path
list.files() #visualize files in the folder

df.fec <- read.csv("2019-fec-new-queens.csv") #newly eclosed and mated queens
df.ages <- read.csv("2017-fec-all-ages.csv") #egg laying by different queen ages

# Adding variables and cleaning data  ----

# Step 1: adding needed variables
df.fec <- df.fec %>%
  dplyr::mutate(date = as.Date(date,"%m/%d/%Y"),
                total.brood = eggs + yl + ol + worker.pp + worker.pup,
                total.adults = males + queens + workers,
                brd.prp = round(total.brood/(total.brood + total.adults), digits = 3),
                qn.fecundity = round(eggs/queens, digits = 0)) #eggs per queen in whole numbers

df.ages <- df.ages %>%
  dplyr::group_by(colonyID) %>%
  dplyr::rename(eggs48hr = eggs1,
                queens48hr = q1) %>%
  dplyr::mutate(total.brood.0 = yl0 + ol0 + worker.pp0 + worker.pup0,
         qn.fecundity.48hr = round(eggs48hr/queens48hr, digits = 0)) # eggs are laid in whole numbers

# Step 2: selecting only those colonies that have been counted all the time
df.fec <- df.fec %>% 
  dplyr::group_by(colonyID) %>% 
  dplyr::filter(n() >= 12)

# removing columns that are not needed
df.fec.start <- data.frame(subset(df.fec, select = -c(day),
                                  q.age.days == 4)) #census before to assess group composition

df.fec <- data.frame(subset(df.fec, select = -c(day),
                            q.age.days != 4)) #only census when queens are laying eggs

# Step 3: Time to an event
# assigning a value if a certain number of eggs were laid in the colony by a certain time.
# Event = 1, when more than 20 eggs in the colony, if not then '0'
df.fec <- df.fec %>%
  mutate(egg.laying = ifelse(eggs < 40, 0, 1))

# dropping further values when the eggs laid is > 40, i.e., the even occured
df.fec.egg <- df.fec %>%
  dplyr::group_by(colonyID) %>% 
  filter(cumsum(egg.laying == 1) == 1)

# Saving processed files ----
# Step 1: saving file for data analysis
write.csv(df.ages, #data frame
          file.path(out.dir,# file path
                  "2017fecundity-queen-ages.csv"))

write.csv(df.fec, #data frame
          file.path(out.dir,# file path
                    "2019-fecundity-new-queen.csv"))

write.csv(df.fec.egg, #data frame
          file.path(out.dir,# file path
                    "2019-fecundity-new-queen-first-egg.csv"))


# Step 2: saving file with detailed column names
colnames(df.ages) <- c("Colony ID","Observer",
                        "Age of the Queens (months)","Source Colony ID","Wolbachia infection",
                        "Young Larvae", "Late-instar larvae",
                       "Worker pre-pupae", "Sexual pre-pupae",
                       "Worker pupae", "Male pupae", "Queen pupae", "Workers",
                        "Eggs after 48 hrs", "Queens in Colony",
                        "Total amount of Brood", "Eggs per queen")
write.csv(df.ages,
          file.path(out.dir,# file path
                    "2017fecundity-queen-ages-detailednames.csv"))

colnames(df.fec) <- c("Date of census","Wolbachia infection", "Colony ID",
                      "Queen age (days)", "Males",
                     "Workers", "Queens", "Eggs", "Young larvae",
                      "Late-instar larvae", "Worker pre-pupae", "Worker pupae", "Total amount of brood", "Total adults",
                      "Proportion of Brood in the Colonies", "Eggs per queen", "More than 40 eggs")

write.csv(df.fec,
          file.path(out.dir,# file path
                    "2019fecundity-new-queen-detailednames.csv"))

colnames(df.fec.egg) <- c("Date of census","Wolbachia infection", "Colony ID",
                          "Queen age (days)", "Males",
                          "Workers", "Queens", "Eggs", "Young larvae",
                          "Late-instar larvae", "Worker pre-pupae", "Worker pupae", "Total amount of brood", "Total adults",
                          "Proportion of Brood in the Colonies", "Eggs per queen", "More than 40 eggs")
write.csv(df.fec,
  file.path(out.dir,# file path
            "2019-fecundity-new-queen-first-egg-detailednames.csv"))

# ---- Sample size -----
sample.size <- data.frame(aggregate(colonyID ~ wolbachia + q.age.month, df.ages,
                                   function(x) length(unique(x)))) 

sample.size <- data.frame(aggregate(colonyID ~ wolbachia + q.age.days, df.fec,
                                    function(x) length(unique(x)))) 

# --------------- RAW DISTRIBUTION OF DATA ----------
data <- df.fec
data$x = data$eggs

# since there are 0 values, creating a continuity corrected log value
# for visualization only
clog <- function(x) log(x + 0.5)

plot(clog(data$x))
hist(clog(data$x))

# density distribution overlapped with histogram of data
hist(data$x, freq = F)
lines(density(data$x))

# % of data that is 0
100*sum(data$x == 0)/nrow(data)

# Distribution for count data
# if mean is similar to variance = Poisson
# if mean <<< variance = Negative Binomial or Quasipoisson (overdispersed)
c(mean(data$x), var(data$x))


# predicting variance and mean to see which overlaps best with data
fit <- MASS::fitdistr(data$x, densfun = "negative binomial")

hist(data$x, prob=TRUE, main="")
curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)


# ---- Normality of raw data 
#Shapiro Wilks Test for sampling from normal distribution, p < 0.05: not drawn from normal distribution (use GLM)
shapiro.test(data$x) #test: if drawn from normal distribution
qqnorm(log(data$x))
qqline(log(data$x))

#Batlett's test for unequal variance. p-value < 0.05 : variance can not be assumed equal
bartlett.test(x ~ wolbachia, data=data)

# --------------- UNIVERSAL FUNCTIONS FOR MODEL ANALYSIS ------------------
# Goodness of fit for GLM ----
# Chi-square test on the residual deviance and degree of freedon
# ran GLM with Poisson, quasipoisson and negative binomial model
# null hypothesis = family (Poisson, quasipoisson or negative binomial) is an adquate fit for the data
# if P > 0.05, family specified in the model is a good fit, i.e., accept null hypothesis
p.gof <- function(model){
  dev <- deviance(model)
  rdf <- df.residual(model)
  pchisq(dev,
         df=rdf,
         lower.tail=FALSE)
}

# Overdispersion of GLM -----
overdisp_fun <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  } 
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model, type = "pearson") # computes pearson residuals
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}

# if overdispersed then use quasipoisson or negative binomial
# selected between the 2 by looking at 3 criteria in order of preference -
#(a) linear fit of residuals (q-q plot), (b) goodness of fit of modelm & (c) lack of over dispersion

# ------ Significance of model components ---
stats.glm <- function(model){
  stats.glm(model)
}

# Model diganostics ----
# Q-Q plot of model 
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

# Test statistics for model ----
# mixed-effect models
stats.glmer <- function(model){
  drop1(model,.~.,test="Chisq")
}

# quasi-models
stats.glm.q <- function(model){
  drop1(model,.~.,test="F")
}

# non-quasi models
stats.glm <- function(model){
  drop1(model,.~.,test="Chisq")
}

# linear models
stats.lm <- function(model) {
  car::Anova(model, ddf="Satterthwaite")
}
# ---------------- EGGS BY NEWLY MATED AND ECLOSED QUEENS  ----------------- 
# --- Eggs: overall effects ----
# repeated measures since colonies analyzed over a period of time.
df.fec$colID.factor <- as.factor(df.fec$colonyID)

egg.p <- glmer(eggs~wolbachia + queens + q.age.days  +(1|colID.factor),
               data=df.fec,family = poisson(link = log),
               nAGQ = 0, #prevent convergence issues by reducing accuracy of likelihood tests (https://stats.stackexchange.com/questions/77313/why-cant-i-match-glmer-family-binomial-output-with-manual-implementation-of-g)
               control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=3e5))) #fixed effects

eg.intr <- glm.nb (eggs~wolbachia + q.age.days + queens +
                   + wolbachia:q.age.days,
                   data = df.fec) #interaction


model <- eg.intr
vif(model) #testing for multicollinearity
p.gof(model)
overdisp_fun(model)
summary(model) 
stats.glmer(model)
qq(model)
res(model)
shp.tst(model)

# Permutation test to compare Wolbachia-specific change in eggs over time
eggs.perm <- permu.test(eggs~wolbachia + queens|q.age.days,
                        type = "regression",
                        data = df.fec)

eggs.perm[eggs.perm$p >= .05,'t'] <- NA
plot(eggs.perm)

# --- Eggs: Newly mated queen age specific ----
#step 1: creating datasets for each day
df8 <- data.frame(subset(df.fec, q.age.days == 8))
df11 <- data.frame(subset(df.fec, q.age.days == 11))
df14 <- data.frame(subset(df.fec, q.age.days == 14))
df17 <- data.frame(subset(df.fec, q.age.days == 17))
df20 <- data.frame(subset(df.fec, q.age.days == 20))
df23 <- data.frame(subset(df.fec, q.age.days == 23))
df35 <- data.frame(subset(df.fec, q.age.days == 35))
df39 <- data.frame(subset(df.fec, q.age.days == 39))
df43 <- data.frame(subset(df.fec, q.age.days == 43))
df50 <- data.frame(subset(df.fec, q.age.days == 50))

# step 2: age-specific model and model output
eg.time <- glm.nb (eggs~wolbachia + queens + workers,
                   data = df8)

model <- eg.time
summary(model) 
stats.glm(model) 
res(model)
shp.tst(model)
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

# ---------------- EGGS LAID BY AT DIFFERENT AGES ACROSS LIFE HISTORY  ------------
# (Fig 1b: 2017 dataset) Different queen ages : eggs in the colonies over time ----
# Step 1: age-specfic dataset
df1 <- data.frame((subset(df.ages, q.age.month ==1)))
df3 <- data.frame((subset(df.ages, q.age.month ==3)))
df4 <- data.frame((subset(df.ages, q.age.month ==4)))
df6 <- data.frame((subset(df.ages, q.age.month ==6)))
df9 <- data.frame((subset(df.ages, q.age.month ==9)))

egg.ages <- glm.nb (eggs48hr~wolbachia + ol0,
                   data = df9)

model <- egg.ages
vif(model)
summary(model) 
stats.glm(model) 
res(model)
shp.tst(model)
autoplot(model, which = 1:6, ncol = 3, label.size = 3) #diagnostic plots

#--------------- GRAPHS SHOWING TRENDS ------------
# PLOT : Box plot with raw values ---------
my.color <-c ("darkorchid4", "chocolate2")

#Step 1: create universal variable and 
data <- df.fec
data$y <- data$eggs
data$x <- as.factor(data$q.age.days)
y.lab <- "Number of eggs"
x.lab <- "Age of the queens (days) "
title <- "Total Number of Eggs in the Colonies"

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
                                 colour = "black", size = "14"),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.background=element_rect(fill = NA),
        plot.title = element_text(family = "Arial", colour = "black",
                                  face = "bold", size = "18", 
                                  margin = margin(t='0.5', r='0.5', b='0.5', l='0.5', unit = 'cm')))


#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(file = "2019fec-new-eggs.png",
       path = grph,
       dpi = 300, 
       width = 5.5,
       height = 4, 
       units = "in")


# PLOT : Correlation (Spearman Rank) plot ---------
my.color <-c ("chocolate2","darkorchid4")
yr.color <-c ("#52969A","#004549")

#Step 1: make universal variables
data <- df.col
data$y <- data$met.rate.25.col
data$x <- data$mass.before
y.lab <- "metabolic rates \\n(microwatts)"
x.lab <- "weight of the colony (gram)"
title <- "Effect of Colony Mass on \\nMetabolic Rates of \\nWhole Colonies (2018+2019)"

#Step 2: plot
ggplot(data, aes(x=x,y=y,
                 group =)) +
  geom_point(color = "#52969A")+
  geom_smooth(method=lm,
              formula = y~x,
              se = T,
              fullrange = F,
              alpha = 0.2,
              fill = "#52969A",
              color ="#52969A" )+
  labs(title=title,
       x = x.lab,
       y = y.lab) +
  stat_cor(label.y = 1700,
           method = "spearman")+ #position of R-square
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
        legend.position = "bottom",
        legend.text = element_text(family = "Arial",
                                   colour = "black",
                                   size = "12"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_text(family = "Arial",
                                  colour ="black",
                                  face = "bold", size = "16"),
        axis.text = element_text(family = "Arial",
                                 colour = "black", size = "14"),
        axis.ticks = element_line(colour = "black", size = 1),
        plot.background=element_rect(fill = NA),
        plot.title = element_text(family = "Arial", colour = "black",
                                  face = "bold", size = "16", 
                                  margin = margin(t='0.5', r='0.5', b='0.5', l='0.5', unit = 'cm')))

#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(file = "mr20182019-mW-humid.png",
       path = gr.col,
       dpi = 300, 
       width = 5,
       height =4.5, 
       units = "in")

# PLOT : Correlation Matrix ----
# reference:https://rstudio-pubs-static.s3.amazonaws.com/240657_5157ff98e8204c358b2118fa69162e18.html#visualization-of-a-correlation-matrix
# used the correlation matrix for this 

data <- corr18m

png(height= 8.5,
    width= 8,
    units = "in",
    res = 300,
    file="corr2018-sprm-colony.png")

corrplot(data$r, method = "color",
         type="lower", order="alphabet", #half of matrix and alphabetical ordering
         p.mat = data$P, sig.level = 0.05, insig = "blank", # p>0.05 are blank
         col = brewer.pal(n = 6, name = "BrBG"), #color scheme
         #outline = TRUE, #outline for each cell
         tl.col = "black", tl.srt = 45, #change label color and angle
         diag = FALSE, #hide correlation on the diaganol
         title = "Correlation Plot Colony with 2018 Dataset")
dev.off()

# good reference: https://www.guru99.com/r-pearson-spearman-correlation.html
#pick colors - https://htmlcolorcodes.com/

# plotting only correlations with Metabolic Rate
df.corr2018.mr <- data.frame(subset(df.corr2018, trait1 == "met.rate.25"))

# common variables
data = df.corr2018
x = data$trait1
y = data$trait2
value = data$spearman.corr
ggtitle = "Spearman Correlations \\n2019"

ggplot(data,aes(x,y))+
  geom_tile(aes(fill = value),
            color = "gray70")+
  ggtitle(ggtitle)+
  geom_text(aes(x, y, label = value),
            color = "black", size = 4) +
  scale_fill_gradient2(low = "#FF6347",# color scheme
                       mid = "white",
                       high = "cadetblue",
                       midpoint = 0,
                       limit = c(-1,1),
                       name="Correlations \\n(Pearson)",
                       guide = "colourbar") +
  theme_minimal()+ # minimal theme
  theme(plot.title = element_text(family = "Arial", colour = "black",
                                  face = "bold", size = "12", 
                                  margin = margin(t='0.5', r='0.5',
                                                  b='0.5', l='0.5', 
                                                  unit = 'cm')),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(family = "Arial",
                                   colour = "black",
                                   size = "10"),
        legend.background = element_rect(fill = NA,linetype = NULL),
        legend.key = element_rect(color = NA, fill = NA),
        axis.title = element_blank(),
        axis.text = element_text(family = "Arial",
                                 colour = "black",
                                 size = 11),
        axis.text.x = element_text(angle = 45,vjust = 1,size = 11,hjust = 1))        
#axis.ticks = element_line(colour = "black", size = 1)

#saving image in desired format to the file path set in the start of the script. Default saves last plot.
ggsave(filename = "corr2019-sprm-mr.png",
       path = gr.col,
       plot = last_plot(),
       dpi = 300, 
       width = 3.5,
       height = 4, 
       units = "in")

# PLOT: Individual colony data ----
my.color <-c ("darkorchid4","chocolate2" )

# Step 1: assign a variable and labels for the plot
data <- df.wk
data$x <- data$day
data$y <- data$wk.alive
title <- "Number of Alive  Workers since Last Time Point"
y.lab <- "number of alive workers"
x.lab <- "time (day) since starting the assay"

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
  scale_x_discrete(breaks = 4,
                   limits = seq(1, 92))+ #y-axis ticks seq(lower limit, top limit, incremental steps
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



# PLOT: Colony-level dynamics along with Mean --------------
my.color <-c ("darkorchid4","chocolate2" )

# Step 1: assign a variable and labels for the plot
data <- df.wk
data$x <- data$day #x-axis
data$y <- data$wrk.surv.prp.overall #y-axis
title <- "Proportion of Alive Workers"
y.lab <- "proportion of alive workers"
x.lab <- "time (days) since starting the assay"

# Step 2: reorder variablesa 
data.mean <- summarySE(data, measurevar="y",
                       groupvars=c("wolbachia", "day"), na.rm = T)

# Step 3: reorder variablesa 
data.mean$wolbachia<- factor(data.mean$wolbachia, 
                             levels = c("infected",
                                        "uninfected"))

# Step 4: Plot = sample values overlaid with mean and 95% CI. Texts were later edited in Inkscape
ggplot(data.mean, aes(x=day,
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
  scale_x_discrete(limits = seq(0, 100, 10))+ #y-axis ticks seq(lower limit, top limit, incremental steps
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
ggsave(file = "WkSurv-SurvPrp-mean.png",
       #path = gr.wk,
       dpi = 300, 
       width = 5,
       height = 4, 
       units = "in")
