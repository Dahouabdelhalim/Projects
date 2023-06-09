# Load required packages
require(car)
require(lme4)
require(ggplot2)
require(psych)
require(AICcmodavg)
require(boot)

# This helper script imports and formats the data, then creates the distance
# matricies and calculates the local sex ratios for all individuals across
# all scales. It is not well vectorized and so it takes a long time to run.

source('formatData.R')


# Test whether mean sex ratio changes among levels using a general linear mixed
# model

varModel <- glmer(sexRatio ~ sex*level + (1|name), weights=density, 
                  data=varAnalysis, family=binomial(logit), 
                  control=glmerControl(optimizer="bobyqa"))
varAnalysis$resid <- resid(varModel)
Anova(varModel, type=3)

# Visually insepct the histogram of residuals and a plot of the 
# predicted values against the residuals to test for deviations from
# regression model assumptions
par(mfrow=c(1,2))
hist(resid(varModel))
plot(resid(varModel) ~ predict(varModel))
dev.off()

### calculate Cij, Pii, Pjj to compare the components of the selection
### differential

## Cij is the covariance between an individual's sex and social partner sex
## ratios calculated as the sum of products of the deviation of individual sex
## from mean sex of population and the deviation of local sex ratio from the 
## population sex ratio

cij <- as.character(levels(varAnalysis$level))
cij <- as.numeric(levels(varAnalysis$level))
cij <- rbind(cij, rep(NA, length(levels(varAnalysis$level))), 
             rep(NA, length(levels(varAnalysis$level))))
            
for (i in 1:length(levels(varAnalysis$level))) {
  cij[2,i] <- cov(varAnalysis$sexEncode[varAnalysis$level==levels(varAnalysis$level)[i]], varAnalysis$sexRatio[varAnalysis$level==levels(varAnalysis$level)[i]])
  cij[3,i] <- cor.test(varAnalysis$sexEncode[varAnalysis$level==levels(varAnalysis$level)[i]], varAnalysis$sexRatio[varAnalysis$level==levels(varAnalysis$level)[i]])$p.value
}

## Pii is just the variance in sex at each level. It should ideall be the same...

pii <- as.character(levels(varAnalysis$level))
pii <- as.numeric(levels(varAnalysis$level))
pii <- rbind(pii, rep(NA, length(levels(varAnalysis$level))))

for (i in 1:length(levels(varAnalysis$level))) {
  pii[2,i] <- var(varAnalysis$sexEncode[varAnalysis$level==levels(varAnalysis$level)[i]])
}

## Pjj is the variance in local sex ratios at each level

pjj <- as.character(levels(varAnalysis$level))
pjj <- as.numeric(levels(varAnalysis$level))
pjj <- rbind(pjj, rep(NA, length(levels(varAnalysis$level))))
pjj <- rbind(pjj, rep(NA, length(levels(varAnalysis$level))))
pjj <- rbind(pjj, rep(NA, length(levels(varAnalysis$level))))

for (i in 1:length(levels(varAnalysis$level))) {
  pjj[2,i] <- var(varAnalysis$sexRatio[varAnalysis$level==levels(varAnalysis$level)[i] & varAnalysis$sex=="hermaphrodite"])
  pjj[3,i] <- var(varAnalysis$sexRatio[varAnalysis$level==levels(varAnalysis$level)[i] & varAnalysis$sex=="female"])
  pjj[4,i] <- var(varAnalysis$sexRatio[varAnalysis$level==levels(varAnalysis$level)[i]])
}


## Cjj is the correlation of social phenotypes

cjj <- corr.test(localSex[,8:19])


######################
# Supplemental methods
######################

#Calculate what a hypothetical beta-N would be if females have
#a two-fold fitness advantage over hermaphrodites

# Dummy variable for fitness
localSex$fakeFit <- NA
localSex$sexCode <- NA

localSex$sex <- as.factor(localSex$sex)
localSex$sexCode[localSex$sex=="hermaphrodite"] = 1
localSex$sexCode[localSex$sex=="female"] = 0
localSex$fakeFit[localSex$sex=="hermaphrodite"] = 1
localSex$fakeFit[localSex$sex=="female"] = 2

# Calculate a variance-standardized value of the 0-1 index of sex
localSex$sexStd <- standardize(localSex$sexCode)

# Calculate mean-relativized value of fitness
localSex$relFit <- relativize(localSex$fakeFit)

# A simple linear regression of relative fitness on standardized sex
summary(lm(relFit ~ sexStd, data=localSex))


############################################## 
#Plotting functions
##############################################

### FIGURE 1
### Can use to plot the local sex ratio variance thing (like from arc)
sub <- varAnalysis[varAnalysis$level==0.5 | varAnalysis$level==1.5 | 
                     varAnalysis$level==3.5 | varAnalysis$level==4.0,]
sub$level <- factor(sub$level)
g1 <- (ggplot(data=sub, 
             aes(x=Easting, y=Northing))
      + geom_point(aes(fill = sexRatio), pch=21, color="#484848", size=1.5) 
      + scale_fill_gradientn(colours = gray.colors(n=2, start = 1, end = 0))
      + scale_y_continuous(breaks = c(4137080:4137129), labels=rep('', length(c(4137080:4137129))))
      + scale_x_continuous(breaks = c(557111:557177), labels=rep('', length(c(557111:557177))))
      + theme (panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank())
      + facet_wrap(~level, scales="free")
      )

g1 


### FIGURE 2
# Create a new data frame containing the factors that you want to plot
grid <- with(varAnalysis, 
             expand.grid(level = levels(level),
                         sex=levels(sex)))

# Use predictSE to estimate the LS means for all combinations of factor levels
preds <- predictSE(varModel, newdata = grid, type = 'link', se.fit=T)

# Append these LS mean estimates to the "grid" dataset
grid$means <- preds$fit

# Add upper and lower confidence intervals
grid = within(grid, {ucl <- grid$means + 1.96 * preds$se.fit
                     lcl <- grid$means - 1.96 * preds$se.fit})

# Back transform the LS means and confidence intervals from logistic
# to the original units of the data
grid$means <- inv.logit(grid$means)
grid$ucl <- inv.logit(grid$ucl)
grid$lcl <- inv.logit(grid$lcl)
colnames(grid)[2] <- "Sex"

# Dot plot with CIs from the above
pd <- position_dodge(0.4)
(ggplot(grid, aes(x=level, y=means, color=Sex, group=Sex)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1, position=pd) +
  geom_line(aes(group=Sex), position=pd) +
  geom_point(position=pd, size=3) +
  scale_color_manual(values=c("#909090", "black")) + 
  scale_fill_manual(values=c("#909090", "black")) +
  scale_x_discrete(name="Radius of social neighborhood (meters)", breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)) +
  scale_y_continuous(name="Hermaphrodite frequency", breaks = c(0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.8)) +
  theme(axis.title.x=element_text(size=18), 
           axis.title.y=element_text(size=18, angle=90), 
           axis.text.x=element_text(size=14, color='black'), 
           axis.text.y=element_text(size=14, color='black'),
           axis.line = element_line(color = 'black'),
           axis.ticks = element_line(color = 'black'),
           legend.title=element_text(size=18, face='bold', hjust=0),
           legend.text=element_text(size=18),
           strip.text.x = element_text(size=18, face="bold"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank()))

### FIGURE 3: Plot of Pjj
variancePlot <- data.frame(cbind(pjj[1,], pjj[2,], pjj[3,], pjj[4,]))

g3 = (ggplot(variancePlot, aes(x=X1, y=X4))
     + geom_line(size=1)
     + geom_point(size=3.5)
     + scale_x_continuous(name='Radius of social neighborhood (meters)', breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0))
     + scale_y_continuous(name='Variance in hermaphrodite frequency')
     + theme(axis.title.x=element_text(size=18), 
             axis.title.y=element_text(size=18, angle=90), 
             axis.text.x=element_text(size=14, color='black'), 
             axis.text.y=element_text(size=14, color='black'),
             axis.line = element_line(color = 'black'),
             legend.title=element_text(size=18, face='bold', hjust=0),
             legend.text=element_text(size=18),
             strip.text.x = element_text(size=18, face="bold"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             panel.background = element_blank()))

g3

### FIGURE 4: Plot of Cij
covariancePlot <- data.frame(cbind(cij[1,], cij[2,]))

g4 = (ggplot(covariancePlot, aes(x=X1, y=X2)) 
 + geom_point(size=3.5)
 + scale_x_continuous(name='Radius of social neighborhood (meters)', breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0))
 + scale_y_continuous(name='Phenotypic covarinace')
 + geom_line(size=1)
 + theme(axis.title.x=element_text(size=18), 
         axis.title.y=element_text(size=18, angle=90), 
         axis.text.x=element_text(size=14, color='black'), 
         axis.text.y=element_text(size=14, color='black'),
         axis.line = element_line(color = 'black'),
         axis.ticks = element_line(color = 'black'),
         legend.title=element_text(size=18, face='bold', hjust=0),
         legend.text=element_text(size=18),
         strip.text.x = element_text(size=18, face="bold"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank()))

g4

