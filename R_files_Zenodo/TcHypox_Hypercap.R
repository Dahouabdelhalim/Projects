`# Set WD to preferred location
setwd("~/Desktop/OneDrive - University of North Carolina at Chapel Hill/Copepod Research/Tc Hypox & Hypercap")

library(dplyr) #needed for %>%
library(car) #needed for vif
library(lme4) #needed for glmer
library(lmerTest) #needed for p values in glmer
library(numDeriv) #needed for checking scaling gradients in model assumptions
library(MuMIn) #needed for dredge
library(nloptr) #needed for optimizers
library(optimx) #needed for optimizers
library(ggplot2) #needed for plots
library(plotrix) #needed for std.error
library(nlme) #needed for gls
library(ggpubr) #needed for grid.arrange
library(grid) #needed for scales package
library(gridExtra) #needed for grid.arrange
library(scales) #needed for plot commands

# Low pH Analysis
#Check sample size, mean, and SD for treatment and control
hypercap.24 <- read.csv("Tc Hypercap indiv 24.csv")
group_by(hypercap.24, treatment) %>%
  summarise(
    count = n(),
    mean = mean(alive, na.rm = TRUE),
    sd = sd(alive, na.rm = TRUE)
  ) 
#Anova to see that treatment and control mean are different
hc.aov <- aov(alive ~ treatment, data = hypercap.24) 
summary(hc.aov)
#Running GLMER of pH treatment to measure effect sizes of fixed effects
hc.24 <- subset(hypercap.24, treatment == "pH")
#Length may be allometricly related to tolerance, taking the log-length
hc.24$log.length <- log(hc.24$length)
#Ensuring data are have correct collection information
pop.index <- c("AB", "BB", "FHL", "SCN", "SD", "SS")
lat.values <- c(33.44444, 38.18308, 48.54600, 36.94962, 32.74574, 34.00017)
coll.values <- c("5/24/11", "7/9/13", "4/3/18", "11/21/17", "11/26/17", "9/27/04")
hc.24$coll.date <- coll.values[match(hc.24$pop, pop.index)]
hc.24$lat <- lat.values[match(hc.24$pop, pop.index)]
hc.24$coll.date <- as.Date(hc.24$coll.date, "%m/%d/%y")
hc.24$year <- as.numeric(format(hc.24$coll.date,'%Y'))
#Centering and scaling continuous variables so better interpretation with categorical variables
cent.len <- function(x){
  y <- (x - mean(hc.24$log.length))/(2*sd(hc.24$log.length))
  return(y)
}
cent.lat <- function(x){
  y <- (x - mean(hc.24$lat))/(2*sd(hc.24$lat))
  return(y)
}
cent.year <- function(x){
  y <- (x - mean(hc.24$year))/(2*sd(hc.24$year))
  return(y)
}
hc.24$sc.len <- cent.len(hc.24$log.length)
hc.24$sc.lat <- cent.lat(hc.24$lat)
hc.24$sc.year <- cent.year(hc.24$year)
#Generating GLMER model
hc.complex.mod <- glmer(alive ~ sex + sc.len + sc.lat + sc.year + (1|position) + (1|batch/group), data = hc.24, family = binomial)
#Check for colinearity of fixed effects, singularity, residuals and overdispersion.
vif(hc.complex.mod)
isSingular(hc.complex.mod)
qqPlot(resid(hc.complex.mod))
hist(resid(hc.complex.mod))
overdisp_fun<- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(hc.complex.mod)
#Dredging the model to explore model space
options(na.action = na.fail)
hc.dredge <- dredge(hc.complex.mod, beta = "sd", rank = "AICc")
options(na.action = na.omit)
#Averaging the best fit models, calculating confidence intervals, and calculating importance of each predictor in averaged model
hc.avg <- model.avg(hc.dredge, subset = delta < 5, beta = "sd")
summary(hc.avg)
confint(hc.avg)
importance(hc.dredge)

# Low pH Plots
#Transforming raw data (0,1) into continuous response using model predictions
hc.24$pred <- predict(hc.avg, type = "response")
#Generating binomial plots
hc.24$lat.fact <- factor(hc.24$lat, levels = c("32.74574", "33.44444", "34.00017", "36.94962", "38.18308", "48.546"), labels = c("SD", "AB", "SS", "SCN", "BB", "FHL"))
hc.24$date.fact <- factor(hc.24$coll.date, levels = c("2004-09-27", "2011-05-24", "2013-07-09", "2017-11-21", "2017-11-26", "2018-04-03"), labels = c("SS", "AB", "BB", "SCN", "SD", "FHL"))
p3.hc <- ggplot(hc.24, aes(x = length, y = pred)) + geom_smooth(aes(y = as.numeric(alive), linetype = sex), colour="black", method = "glm",method.args = list(family = "binomial")) + facet_wrap(~lat.fact, labeller = labeller()) + labs(x = "Length (mm)", y = "Probability of Survival", color = "sex") + theme_classic() + geom_point(aes(x = length, y = as.numeric(hc.24$alive), shape = sex), color = "black")

# Hypoxia Analysis
anox.30 <- read.csv("Tc Anox indiv 30.csv")
#Check sample size, mean, and SD for treatment and control
group_by(anox.30, treatment) %>%
  summarise(
    count = n(),
    mean = mean(alive, na.rm = TRUE),
    sd = sd(alive, na.rm = TRUE)
  )
anox.aov <- aov(alive ~ treatment, data = anox.30)
summary(anox.aov)
#Organize data and ensure it has the correct collection information
ax.30 <- subset(anox.30, treatment == "anoxia")
ax.30 <- subset(ax.30, select = -notes)
pop.index <- c("AB", "BB", "FHL", "SCN", "SD", "SS")
lat.values <- c(33.44444, 38.18308, 48.54600, 36.94962, 32.74574, 34.00017)
coll.values <- c("5/24/11", "7/9/13", "4/3/18", "11/21/17", "11/26/17", "9/27/04")
ax.30$coll.date <- coll.values[match(ax.30$pop, pop.index)]
ax.30$lat <- lat.values[match(ax.30$pop, pop.index)]
ax.30$log.length <- log(ax.30$length)
ax.30$coll.date <- as.Date(ax.30$coll.date, "%m/%d/%y")
ax.30$year <- as.numeric(format(ax.30$coll.date,'%Y'))
ax.30$position <- as.factor(ax.30$position)
#Center and scale continuous variables
scl.len <- function(x){
  y <- (x - mean(ax.30$log.length))/(2*sd(ax.30$log.length))
  return(y)
}
scl.lat <- function(x){
  y <- (x - mean(ax.30$lat))/(2*sd(ax.30$lat))
  return(y)
}
scl.year <- function(x){
  y <- (x - mean(ax.30$year))/(2*sd(ax.30$year))
  return(y)
}
ax.30$sc.len <- scl.len(ax.30$log.length)
ax.30$sc.lat <- scl.lat(ax.30$lat)
ax.30$sc.year <- scl.year(ax.30$year)
#Generating non-interactive model to test colinearity of terms
ax.30.mod <- glmer(alive ~ sex + sc.len + sc.lat + sc.year + (1|position) + (1|batch/group), data = ax.30, family = binomial)
vif(ax.30.mod)
#Generating interactive model for analysis
ax.30.int <- glmer(alive ~ (sex + sc.len + sc.lat + sc.year)^2 + (1|position) + (1|batch/group), data = ax.30, family = binomial)
#Failure to converge warning, so checking singularity, absolute and scaled gradient, and optimizers
isSingular(ax.30.int)
#The model is not singular
dd <- update(ax.30.int, devFunOnly=TRUE)
pars <- unlist(getME(ax.30.int, c("theta", "fixef")))
grad2 <- grad(dd, pars)
hess2 <- hessian(dd, pars)
sc_grad2 <- solve(hess2, grad2)
max(pmin(abs(sc_grad2), abs(grad2)))
#Absolute and scaled gradient below tolerance threshold
m2 <- update(ax.30.int, control = glmerControl(optCtrl = list(maxfun = 2e4)))
m3 <- update(ax.30.int, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
m4 <- update(ax.30.int,control=glmerControl(optimizer="Nelder_Mead",
                                            optCtrl=list(maxfun=2e5)))
m5 <- update(ax.30.int, control = glmerControl(optimizer = 'optimx', optCtrl=list(method='L-BFGS-B')))
m6 <- update(ax.30.int,control=glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
m7 <- update(ax.30.int,control=glmerControl(optimizer='nloptwrap', 
                                            optCtrl=list(method = 'NLOPT_LN_NELDERMEAD')))
m8 <- update(ax.30.int,control=glmerControl(optimizer='nloptwrap', optCtrl=list(method = 'NLOPT_LN_BOBYQA')))
logLik(m3)
logLik(m5)
logLik(m6)
#All working optimizers equally functional so using 'bobyqa' because I wrote it first
ax.30.int <- glmer(alive ~ (sex + sc.len + sc.lat + sc.year)^2 + (1|position) + (1|batch/group), data = ax.30, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
#Checking residuals and overdispersion
qqPlot(resid(ax.30.int))
hist(resid(ax.30.int))
overdisp_fun(ax.30.int)
#Dredging the model to explore model space
options(na.action = na.fail)
ax.dredge <- dredge(ax.30.int, beta = "sd", rank = "AICc")
options(na.action = na.omit)
#Averaging the best fit models, calculating confidence intervals, and calculating importance of each predictor in averaged model
ax.avg <- model.avg(ax.dredge, subset = delta < 5, beta = "sd")
summary(ax.avg)
ax.avg.CI <- confint(ax.avg)
importance(ax.dredge)

# Hypoxia Plots
#Transforming raw data (0,1) into continuous response using model predictions
ax.30$pred <- predict(ax.avg, type = "response")
#Generating binomial plots
ax.30$lat.fact <- factor(ax.30$lat, levels = c("32.74574", "33.44444", "34.00017", "36.94962", "38.18308", "48.546"), labels = c("SD", "AB", "SS", "SCN", "BB", "FHL"))
ax.30$date.fact <- factor(ax.30$coll.date, levels = c("2004-09-27", "2011-05-24", "2013-07-09", "2017-11-21", "2017-11-26", "2018-04-03"), labels = c("SS", "AB", "BB", "SCN", "SD", "FHL"))
p3.ax <- ggplot(ax.30, aes(x = length, y = pred)) + geom_smooth(aes(y = as.numeric(alive), linetype = sex), colour="black", method = "glm",method.args = list(family = "binomial")) + facet_wrap(~lat.fact, labeller = labeller()) + labs(x = "Length (mm)", y = "Probability of Survival", color = "sex") + theme_classic() + geom_point(aes(x = length, y = as.numeric(ax.30$alive), shape = sex), color = "black")

# Physiochemical Analyses

field.2021 <- read.csv("~/Desktop/OneDrive - University of North Carolina at Chapel Hill/Copepod Research/Field Data/Summer 2021/2021 Summer Field Measurements - Sheet1.csv")
field.2021$date.time <- as.POSIXct(field.2021$date.time, format = "%m/%d/%Y %H:%M") 

pH.fit <- lm(pH ~ latitude:time.frame, field.2021)
pH.anova <- anova(pH.fit)
pH.CI <- confint(pH.fit)
DO.fit <- lm(DO ~ latitude:time.frame, field.2021)
DO.anova <- anova(DO.fit)
DO.CI <- confint(DO.fit)

# Latitude Plot
p1.ax <- ggplot(ax.30, aes(x = as.factor(lat), y = pred)) + geom_boxplot(aes(fill = as.factor(lat))) + theme_classic() + scale_x_discrete() + ylab(expression("Hypoxia Survival")) + scale_fill_brewer(palette="BrBG") + theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank()) + geom_text(aes(x = 0.6, y = 1, label = "A"))
p1.hc <- ggplot(hc.24, aes(x = as.factor(lat), y = pred)) + geom_boxplot(aes(fill = as.factor(lat))) + theme_classic() + scale_x_discrete(labels = c("SD\\n32.74574", "AB\\n33.44444", "SS\\n34.00017", "SCN\\n36.94962", "BB\\n38.18308", "FHL\\n48.546"))  + ylab(expression("Low pH Survival")) + xlab(expression("Population "(southern %->% northern))) + scale_fill_brewer(palette="BrBG") + theme(legend.position = "none") + geom_text(aes(x = 0.6, y = 1, label = "B"))
grid.arrange(p1.ax, p1.hc)

# Covariance Analysis
#Combining low pH and hypoxia data sets into single dataframe then find means by population and plot means with standard deviation
ax.corr <- group_by(ax.30, pop) %>% summarise(mean = mean(as.numeric(alive), na.rm=TRUE), sd = sd(as.numeric(alive), na.rm=TRUE))
pH.corr <- group_by(hc.24, pop) %>% summarise(mean = mean(as.numeric(alive), na.rm=TRUE), sd = sd(as.numeric(alive), na.rm=TRUE))
colnames(pH.corr) <- c("pop", "pH.mean", "pH.sd")
corr.df <- merge(ax.corr, pH.corr, by = "pop")
corr.df$pop.levels <- factor(corr.df$pop, levels = c("SD", "AB", "SS", "SCN", "BB", "FHL"))
cor.test(corr.df$mean, corr.df$pH.mean, method = "pearson")
ggplot(corr.df, aes(mean, pH.mean, color = pop.levels)) + geom_point() +  geom_errorbarh(aes(xmax = mean + sd, xmin = mean - sd)) + geom_errorbar(aes(ymax = pH.mean + pH.sd, ymin = pH.mean - pH.sd)) + xlab("Mean Hypoxia Tolerance") + ylab("Mean Acid Tolerance") + labs(color = "Populations") + theme_classic() + scale_color_brewer(palette = "BrBG")


# BB Rockpool pH and Temperature analysis
BB.meas <- read.csv("17-18 pH and Temp measurements.csv")
BB.meas <- BB.meas[-c(5),] #Removing row of pool 5 which evaporated after the first day
BB.meas$lab <- with(BB.meas, interaction(date, sun))
BB.meas$tp <- as.factor(BB.meas$tp)
BB.temp <- aggregate(temp ~ lab, data = BB.meas, FUN = function(x) c(mean = mean(x), se = std.error(x)))
BB.temp <- do.call(data.frame, BB.temp)
BB.temp$order <- as.factor(c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10))
BB.pH <- aggregate(pH ~ lab, data = BB.meas, FUN = function(x) c(mean = mean(x), se = std.error(x)))
BB.pH <- do.call(data.frame, BB.pH)
BB.pH$order <- as.factor(c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10))
BB.index <- BB.pH$lab
BB.ocean <- aggregate(ocean.pH ~ lab, data = BB.meas, FUN = mean)
BB.ocean$order <- c(2, 4, 6, 8, 10)
dates <- c("Dec 30", "", "Dec 31", "", "Jan 1", "", "Jan 2", "", "Jan 3", "")
ggplot(BB.pH, aes(x = order, y = pH.mean)) + geom_errorbar(aes(ymin = pH.mean - pH.se, ymax = pH.mean + pH.se), width = 0.5) + geom_line(group = 1) + geom_line(data = BB.ocean, aes(x = order, y = ocean.pH), color = "blue") + theme_classic() + theme(axis.title = element_text(size = 18), axis.text = element_text(size=16)) + labs(y = "pH", x = element_blank()) + scale_x_discrete(labels = dates)
#BB pH measurements
mean.pH <- aggregate(pH ~ sun, data = BB.meas, FUN = function(x) c(mean = mean(x), se = std.error(x)))
mean.pH
mean(BB.meas$ocean.pH, na.rm = TRUE)
std.error(BB.meas$ocean.pH, na.rm = TRUE)
BB.simp <- BB.meas[, c('tp', 'date', 'sun', 'pH')]
BB.simp %>%
  tidyr::spread(key = sun, value = pH) %>%
  dplyr::mutate(diff = set - rise)
mean.pH %>%
  tidyr::spread(key = sun, value = pH) %>%
  dplyr::mutate(diff = set - rise)

# Off-gas rate in low pH assays
hypercap.groups <- read.csv("Tc Hypercap groups.csv")
hc.groups <- subset(hypercap.groups, treatment == "pH")
hc.pH <- aggregate(pH ~ time.lapse, data = hc.groups, FUN = function(x) c(mean = mean(x), s.d. = sd(x)))
hc.pH <- do.call(data.frame, hc.pH)
hc.pH$init <- ifelse(hc.pH$time.lapse == 0:000, "sv", "tv")
hc.pH.plot <- ggplot(hc.pH, aes(x = time.lapse, y = pH.mean, color = init)) + geom_errorbar(aes(ymin = pH.mean - pH.s.d., ymax = pH.mean + pH.s.d.), width = 0.5) + geom_line(aes(group = init)) + theme_classic() + theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), legend.position = "none") + labs(y = "pH", x = "Time Lapse (hours)") + scale_color_manual(values=c("gray", "black")) + geom_text(x = 3, y = 8, label = "Pre-Injection pH", color = "gray")
hc.DO <- aggregate(DO ~ time.lapse, data = hc.groups, FUN = function(x) c(mean = mean(x), s.d. = sd(x)))
hc.DO <- do.call(data.frame, hc.DO)
hc.DO$init <- ifelse(hc.DO$time.lapse == 0:000, "sv", "tv")
hc.DO.plot <- ggplot(hc.DO, aes(x = time.lapse, y = DO.mean, color = init)) + geom_errorbar(aes(ymin = DO.mean - DO.s.d., ymax = DO.mean + DO.s.d.), width = 0.5) + geom_line(aes(group = init)) + theme_classic() + theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), legend.position = "none") + labs(y = "Dissolved Oxygen (mg/L)", x = "Time Lapse (hours)") + scale_color_manual(values=c("gray", "black")) + geom_text(x = 2.4, y = 6.2, label = "Pre-Injection\\nDO             ", color = "gray")
grid.arrange(hc.pH.plot, hc.DO.plot, ncol = 1)

# Find mean and SD of pH during control run
group_by(hypercap.groups, treatment) %>% 
  summarise(mean = mean(pH, na.rm = TRUE), 
            sd = sd(pH, na.rm = TRUE))

# Find mean and SD of DO during control run
anox.groups <- read.csv("Tc Anox groups.csv")
group_by(anox.groups, treatment) %>% 
  summarise(mean = mean(c(DO.pre, DO.20, DO.30), na.rm = TRUE), 
            sd = sd(c(DO.pre, DO.20, DO.30), na.rm = TRUE))

# 2021 Field Season Analysis and Plots
library(lme4)
library(lmerTest)
library(dplyr)
library(MuMIn)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(plotrix)
library(scales)

#full data set
field.2021 <- read.csv("~/Desktop/OneDrive - University of North Carolina at Chapel Hill/Copepod Research/Field Data/Summer 2021/2021 Summer Field Measurements - Sheet1.csv")
field.2021$date.time <- as.POSIXct(field.2021$date.time, format = "%m/%d/%Y %H:%M") 

# Determining descriptive statistics (mean and SD) by time and metric
field.morning <- subset(field.2021, time.frame == "morning")
field.morning.pH <- aggregate(pH ~ location, data = field.morning, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.morning.DO <- aggregate(DO ~ location, data = field.morning, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.morning.temp <- aggregate(temp ~ location, data = field.morning, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.morning.sal <- aggregate(salinity ~ location, data = field.morning, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.afternoon <- subset(field.2021, time.frame == "afternoon")
field.afternoon.pH <- aggregate(pH ~ location, data = field.afternoon, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.afternoon.DO <- aggregate(DO ~ location, data = field.afternoon, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.afternoon.temp <- aggregate(temp ~ location, data = field.afternoon, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.afternoon.sal <- aggregate(salinity ~ location, data = field.afternoon, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.evening <- subset(field.2021, time.frame == "evening")
field.evening.pH <- aggregate(pH ~ location, data = field.evening, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.evening.DO <- aggregate(DO ~ location, data = field.evening, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.evening.temp <- aggregate(temp ~ location, data = field.evening, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.evening.sal <- aggregate(salinity ~ location, data = field.evening, FUN = function(x) c(mean = mean(x), sd = sd(x)))

field.morning.pH
field.afternoon.pH
field.evening.pH
field.morning.DO
field.afternoon.DO
field.evening.DO
field.morning.temp
field.afternoon.temp
field.evening.temp
field.morning.sal
field.afternoon.sal
field.evening.sal

# Plots of metrics by location
field.2021.pH <- aggregate(pH ~ location + date.time, data = field.2021, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.2021.pH <- do.call(data.frame, field.2021.pH)

as.numeric(as.POSIXct("8/12/2021 00:00", format = "%m/%d/%Y %H:%M"))

pH.BB <- subset(field.2021.pH, location == "BB")
pH.BB.plot <- ggplot(pH.BB, aes(x = date.time, y = pH.mean)) + geom_errorbar(aes(ymin = pH.mean - pH.sd, ymax = pH.mean + pH.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = "pH", x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-12 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-15 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(7,9.5)) 

pH.SC <- subset(field.2021.pH, location == "SC")
pH.SC.plot <- ggplot(pH.SC, aes(x = date.time, y = pH.mean)) + geom_errorbar(aes(ymin = pH.mean - pH.sd, ymax = pH.mean + pH.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = "pH", x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-15 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-18 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(7,9.5))

pH.AB <- subset(field.2021.pH, location == "AB")
pH.AB.plot <- ggplot(pH.AB, aes(x = date.time, y = pH.mean)) + geom_errorbar(aes(ymin = pH.mean - pH.sd, ymax = pH.mean + pH.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = "pH", x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-18 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-22 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(7,9.5)) 

pH.SD <- subset(field.2021.pH, location == "SD")
pH.SD.plot <- ggplot(pH.SD, aes(x = date.time, y = pH.mean)) + geom_errorbar(aes(ymin = pH.mean - pH.sd, ymax = pH.mean + pH.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = "pH", x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-22 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-26 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(7,9.5)) 


field.2021.DO <- aggregate(DO ~ location + date.time, data = field.2021, FUN = function(x) c(mean = mean(x), sd = sd(x)))

field.2021.DO <- do.call(data.frame, field.2021.DO)

DO.BB <- subset(field.2021.DO, location == "BB")
DO.BB.plot <- ggplot(DO.BB, aes(x = date.time, y = DO.mean)) + geom_errorbar(aes(ymin = DO.mean - DO.sd, ymax = DO.mean + DO.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = "DO (mg/L)", x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-12 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-15 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(1,16)) + geom_text(x=1628744000, y=14, label="2")

DO.SC <- subset(field.2021.DO, location == "SC")
DO.SC.plot <- ggplot(DO.SC, aes(x = date.time, y = DO.mean)) + geom_errorbar(aes(ymin = DO.mean - DO.sd, ymax = DO.mean + DO.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = "DO (mg/L)", x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-15 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-18 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(1,16)) + geom_text(x=1629000300, y=14, label="3")

DO.AB <- subset(field.2021.DO, location == "AB")
DO.AB.plot <- ggplot(DO.AB, aes(x = date.time, y = DO.mean)) + geom_errorbar(aes(ymin = DO.mean - DO.sd, ymax = DO.mean + DO.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = "DO (mg/L)", x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-18 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-22 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(1,16)) + geom_text(x=1629259500, y=14, label="5")

DO.SD <- subset(field.2021.DO, location == "SD")
DO.SD.plot <- ggplot(DO.SD, aes(x = date.time, y = DO.mean)) + geom_errorbar(aes(ymin = DO.mean - DO.sd, ymax = DO.mean + DO.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = "DO (mg/L)", x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-22 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-26 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(1,16)) + geom_text(x=1629605300, y=14, label="6")

field.2021.temp <- aggregate(temp ~ location + date.time, data = field.2021, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.2021.temp <- do.call(data.frame, field.2021.temp)

temp.BB <- subset(field.2021.temp, location == "BB")
temp.BB.plot <- ggplot(temp.BB, aes(x = date.time, y = temp.mean)) + geom_errorbar(aes(ymin = temp.mean - temp.sd, ymax = temp.mean + temp.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = "Temperature (Â°C)", x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-12 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-15 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(14, 35)) + geom_text(x=1628744000, y=33, label="2")

temp.SC <- subset(field.2021.temp, location == "SC")
temp.SC.plot <- ggplot(temp.SC, aes(x = date.time, y = temp.mean)) + geom_errorbar(aes(ymin = temp.mean - temp.sd, ymax = temp.mean + temp.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = element_blank(), x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-15 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-18 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(14, 35)) + geom_text(x=1629000300, y=33, label="3")

temp.AB <- subset(field.2021.temp, location == "AB")
temp.AB.plot <- ggplot(temp.AB, aes(x = date.time, y = temp.mean)) + geom_errorbar(aes(ymin = temp.mean - temp.sd, ymax = temp.mean + temp.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = element_blank(), x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-18 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-22 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(14, 35)) + geom_text(x=1629259500, y=33, label="5")

temp.SD <- subset(field.2021.temp, location == "SD")
temp.SD.plot <- ggplot(temp.SD, aes(x = date.time, y = temp.mean)) + geom_errorbar(aes(ymin = temp.mean - temp.sd, ymax = temp.mean + temp.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = element_blank(), x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-22 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-26 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(14, 35)) + geom_text(x=1629605300, y=33, label="6")

field.2021.sal <- aggregate(salinity ~ location + date.time, data = field.2021, FUN = function(x) c(mean = mean(x), sd = sd(x)))
field.2021.sal <- do.call(data.frame, field.2021.sal)

sal.BB <- subset(field.2021.sal, location == "BB")
sal.BB.plot <- ggplot(sal.BB, aes(x = date.time, y = salinity.mean)) + geom_errorbar(aes(ymin = salinity.mean - salinity.sd, ymax = salinity.mean + salinity.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = "Salinity (ppt)", x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-12 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-15 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(30, 70)) 

sal.SC <- subset(field.2021.sal, location == "SC")
sal.SC.plot <- ggplot(sal.SC, aes(x = date.time, y = salinity.mean)) + geom_errorbar(aes(ymin = salinity.mean - salinity.sd, ymax = salinity.mean + salinity.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = element_blank(), x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-15 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-18 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(30, 70))

sal.AB <- subset(field.2021.sal, location == "AB")
sal.AB.plot <- ggplot(sal.AB, aes(x = date.time, y = salinity.mean)) + geom_errorbar(aes(ymin = salinity.mean - salinity.sd, ymax = salinity.mean + salinity.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = element_blank(), x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-18 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-22 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(30, 70)) 

sal.SD <- subset(field.2021.sal, location == "SD")
sal.SD.plot <- ggplot(sal.SD, aes(x = date.time, y = salinity.mean)) + geom_errorbar(aes(ymin = salinity.mean - salinity.sd, ymax = salinity.mean + salinity.sd), width = 10000) + geom_line() + theme_classic() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12), axis.text = element_text(size=10)) + labs(y = element_blank(), x = element_blank()) + scale_x_datetime(labels=date_format("%H:%m"), breaks = date_breaks(width = "12 hours"), expand=c(0,0)) + xlim(c(as.POSIXct('2021-08-22 00:00:00', format = "%Y-%m-%d %H:%M:%S"), as.POSIXct('2021-08-26 00:00:00', format = "%Y-%m-%d %H:%M:%S"))) + ylim(c(30, 70)) 

dev.new(width=10, height=18, unit="in")
x11(width=10, height=18)
grid.arrange(DO.BB.plot, pH.BB.plot, DO.SC.plot, pH.SC.plot, DO.AB.plot, pH.AB.plot, DO.SD.plot, pH.SD.plot, ncol = 2)
grid.arrange(temp.BB.plot, temp.SC.plot, temp.AB.plot, temp.SD.plot, sal.BB.plot, sal.SC.plot, sal.AB.plot, sal.SD.plot, ncol = 4)

# Linear models of DO and pH

DO.fit <- lm(DO ~ latitude:time.frame, data = field.2021)
summary(DO.fit)
anova(DO.fit)

pH.fit <- lm(pH ~ latitude:time.frame, data = field.2021)
summary(pH.fit)
anova(pH.fit)
