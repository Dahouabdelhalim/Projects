#### MFPS Prevalence 2015 field data ####

# Updated 3-17-20

data<-read.csv("MFPS_prev_3-17-20.csv")
head(data)
length(data$fan_num) # 1636 sea fans

# Normality/ distribution - binomial response

# Prevalence calculations
summary(data$mfps_fan) #
509/(509+1127)*100 # 31.1% prevalence 

#### REGRESSORS SET UP ########

# Numeric
data$temp<-as.numeric(data$temp) 
class(data$temp) 
class(data$temp_sd) # numeric
data$depth<-as.numeric(data$depth)
class(data$depth) # numeric
class(data$coralcov) #numeric
data$sf_num<-as.numeric(data$sf_num)
class(data$sf_num) # 
data$size<-as.numeric(data$size) # 
class(data$size) # numeric

# Factors
class(data$site) 
class(data$cu_all) 
class(data$mfps_fan) 

# CHECKING Data
plot(fan_num~site,data=data)
plot(transect~site,data=data) # 3 per
plot(temp~site,data=data) # 1
plot(temp_sd~site,data=data) # 1
plot(sf_num~site,data=data) # 1
plot(cu_all~site,data=data) # 1
plot(coralcov~transect,data=data) # 1 each
plot(depth~transect,data=data) # 1 each
plot(size~fan_num,data=data) # 1 each


## TEST FOR COLLINEARITY 
head(data)
round(cor(data[,c(7,8,9,10,12,13)]),2)  # only temp metrics, as expected

#### Rescale numeric predictors
library(arm)

data$temp<-rescale(data$temp, binary.inputs="center")
data$temp_sd<-rescale(data$temp_sd, binary.inputs="center")
data$coralcov<-rescale(data$coralcov, binary.inputs="center")
data$depth<-rescale(data$depth, binary.inputs="center")
data$size<-rescale(data$size, binary.inputs="center")
data$sf_num<-rescale(data$sf_num, binary.inputs="center")

### Contrasts - factors ## 
contrasts(data$site) # treatment relative to B 
contrasts(data$cu_all) 
contrasts(data$cu_all) <- cbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1)) # comparing all to level 1


################
# MODELS updated 5-28-19 
##############

# Predictors: size, coralcov, [temp, temp_sd,] cu_all, depth, sf_num_sum

null<-glmer(mfps_fan~1+(1|site),data=data,family='binomial')


# Temperature choice
temp1<-glmer(mfps_fan~temp+(1|site),data=data,family='binomial')
temp2<-glmer(mfps_fan~temp_sd+(1|site),data=data,family='binomial')

AIC(temp1,temp2) # temp1 is better: use mean instead of SD

## ADDITIVE ##
newmfps1<-glmer(mfps_fan~temp+cu_all+depth+coralcov+sf_num+size+(1|site),data=data,family='binomial')
drop1(newmfps1) # drop sf_num
newmfps1b<-glmer(mfps_fan~temp+cu_all+depth+coralcov+size+(1|site),data=data,family='binomial')
drop1(newmfps1b) # drop copper
newmfps1c<-glmer(mfps_fan~temp+depth+coralcov+size+(1|site),data=data,family='binomial')
drop1(newmfps1c) # drop none

### Test size quadratic fr. previous papers
newmfps2<-glmer(mfps_fan~size+I(size^2)+temp+cu_all+depth+coralcov+sf_num_sum+(1|site),data=data,family='binomial')
drop1(newmfps2) # drop sf_num
newmfps2b<-glmer(mfps_fan~size+I(size^2)+temp+cu_all+depth+coralcov+(1|site),data=data,family='binomial')
drop1(newmfps2b) # drop copper
newmfps2c<-glmer(mfps_fan~size+I(size^2)+temp+depth+coralcov+(1|site),data=data,family='binomial')
drop1(newmfps2c) # drop depth
newmfps2d<-glmer(mfps_fan~size+I(size^2)+temp+coralcov+(1|site),data=data,family='binomial')
drop1(newmfps2d) # drop none

### AIC Table ###
library(AICcmodavg) #
cand.set_lb1<-list(null,newmfps1c,newmfps2d)
aiclist<-aictab(cand.set_lb1,modnames=c("null","newmfps1c","newmfps2d"),second.ord=FALSE)
aiclist

### Best model: newmfps2d

summary(newmfps2d) 
plot(newmfps2d)

# Overdispersion
overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model)) - model.df
  # extracts the Pearson residuals
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}

overdisp_fun(newmfps2d) # not overdispersed

# Conf. intervals
library(MuMIn)
confint(newmfps2d,level=0.95, method="Wald") 


### VIF for lmer: Collinearity; legit for glm models### 
# HELPFUL: https://stats.stackexchange.com/questions/70679/which-variance-inflation-factor-should-i-be-using-textgvif-or-textgvif 

# VIF for lmer
vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

vif.mer(newmfps2d) # good


##############
### PLOTS ####
##############

### sjplots: Figure 4 (a) and (b)
mfps_num_mod<-glmer(mfps_num~size+I(size^2)+(1|site),data=data,family='binomial')
summary(data$mfps_num) # 31.1 reflects that percent of 1s
class(data$mfps_num)

library(sjPlot) 
# QQ plot 
sjp.glmer(mfps_num_mod, type = "re.qq") # random effects QQ plot

par(mfrow=c(2,2))
# Probability curve of fixed effects
sjp.glmer(mfps_num_mod, type = "fe.pc",show.ci=TRUE) 

## Showing in two plots
# Two plots - first for size, then size2; all sites on each
par(mfrow=c(2,1))
sjp.glmer(mfps_num_mod,
          type = "ri.pc",
          facet.grid = FALSE)

### Figure 4(c) - plot trend without site as random effect
mfps_glm<-glm(mfps_fan~size+I(size^2),data=data,family='binomial')
par(mfrow=c(1,1))
range(data$size) # -0.47 to 5
xsize <- seq(0, 6, 0.01)
ysize <- predict(mfps_glm, list(size=xsize),type="response")
plot(data$size, data$mfps_num, pch = 16, xlab = "Size (arbitrary unit)", ylab = "MFPS Prevalence")
lines(xsize, ysize)

