##
# Analysis script by Lewis J Bartlett
# Coauthored - Elisa Visher and Mike Boots
# lewis.bartlett@uga.edu
# For associated study:
#'The target of selection matters: An established resistance - development-time negative genetic trade-off is not found when selecting on development time.'
##

#Load packages needed downstream for mixed modelling
library(boot)
library(lme4)
library(afex)
library(emmeans)
library(MASS)
library(npIntFactRep)

## Read in cleaned data
LHD <- read.csv(file = '~LHDMain.csv' , header=T)
IFD <- read.csv(file = '~IFDMain.csv' , header=T)

## Pull out number and names of selected lines,  doses used at assaying etc
Lines <- unique(LHD$Line)

Doses <- unique(IFD$Dose)

NumE <- length(unique(LHD$Line[which(LHD$Treatment == 'Early')]))

NumL <- length(unique(LHD$Line[which(LHD$Treatment == 'Late')]))

## Time for some analysis!

# Clean out damaged larvae (Mass = NA) from life history assats
LHD <- na.exclude(LHD)

ColE <- 'plum4'
ColL <- 'darkgoldenrod3'

# Quick function to flexibly make colours translucent
Transpa <- function(color, percent) {
  
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)
  
  return(t.col)
  
}

# Did selection alter development?

# Examine growth rate, calculated from mass and development time
LHD$GRate <- LHD$Mass/LHD$Days

mixed(GRate ~ Treatment + (1|Treatment/Line) + (1|Block), data = LHD, method = 'LRT')

# Treatment affected growth rate. What direction?
summary(glmer(GRate ~ Treatment + (1|Treatment/Line) + (1|Block), data = LHD))

#Arrange plots
par(mfrow = c(2,2), cex.axis = 1.25, cex.lab = 1.4, mar = c(5,8,2,2))
#Plotting parameter (colour transparency)
Opac <- 65

#Plot growth rate.
boxplot(LHD$GRate ~ LHD$Line, 
        ylab = expression(paste('Growth Rate (mg  ',day^-1,')', sep='')), 
        xlab = 'Line', notch = F, varwidth = F, lwd = 2, pch = 19,
        col = c(rep(Transpa(ColE, Opac), times = NumE), rep(Transpa(ColL, Opac), time = NumL)), border = c(rep(ColE, times = NumE), rep(ColL, time = NumL)))


# Plot and test for differences in development times
boxplot(LHD$Days ~ LHD$Line, ylab = 'Development Time (Days)', xlab = 'Line', notch = F, varwidth = F, lwd = 2, pch = 19,
        col = c(rep(Transpa(ColE, Opac), times = NumE), rep(Transpa(ColL, Opac), time = NumL)), border = c(rep(ColE, times = NumE), rep(ColL, time = NumL)))

# Poisson fit to 'days until pupation' (see data, effectively counts...) - means LRT rather than KR (see ?mixed)
# Use LRT throughout for consistency (will also be necessary for binomial fits for infection data)
mixed(Days ~ Treatment + (1|Treatment/Line) + (1|Block), data = LHD, family = poisson, method = 'LRT')

# Significant differences in development time for selection - check direction

summary(glmer(Days ~ Treatment + (1|Treatment/Line) - 1 + (1|Block), data = LHD, family = poisson))
# Selection happened in intended direction - significantly shorter development time in early / longer in late selected

# get indentity space difference
exp(3.09325) - exp(2.71338)

# Examine for differences in mass at pupation
boxplot(LHD$Mass ~ LHD$Line,  ylab = 'Pupal Mass (mg)', xlab = 'Line', notch = F, varwidth = F, lwd = 2, pch = 19,
        col = c(rep(Transpa(ColE, Opac), times = NumE), rep(Transpa(ColL, Opac), time = NumL)), border = c(rep(ColE, times = NumE), rep(ColL, time = NumL)))

# Same approach as above but normal fit (see data, mass is a typical normal example) and so can use preferred 'KR' method (see ?mixed)

mixed(Mass ~ Treatment + (1|Treatment/Line) + (1|Block), data = LHD, method = 'LRT')

# No real evidence for differences in pupal mass according to treatment

# Decent evidence indeed that late selected are slower growing (take longer, not bigger)

# Early-selected lines are significantly faster growing

# Add legend to plots

plot.new()

legend( x = 0, y = 0.7, c('Early Selected', 'Late Selected'), 
        col = c(ColE, ColL), pch = 19, border = NA, box.col= 'white', cex = 2)

par(mfrow=c(1,1))

### Selected has acted on development
### Do the populations differ in their susceptibility to the pathogen?
### Do we see differences between treatments?

## Test for changes in resistance, in line with the trade-off in this species previously documented

#Check for heterogeneity in dose response

mixed(Infected ~ Dose + Line + Dose:Line + (1|Treatment/Line) + (1|Block), data = IFD, family = binomial, method = 'LRT')

#No between-lines heterogeneity! Great.

mixed(Infected ~ Dose + Line + (1|Treatment/Line) + (1|Block), data = IFD, family = binomial, method = 'LRT')

# Did selection alter resistance to the pathogen?

mixed(Infected ~ Dose + Treatment + (1|Treatment/Line) + (1|Block), data = IFD, family = binomial, method = 'LRT')

# Inconclusive evidence that selection altered resistance
# Check direction.

summary(glmer(Infected ~ Dose + Treatment + (1|Treatment/Line) + (1|Block), data = IFD, family = binomial))

# Late selected (slower developing) more likely to be infected.
# Counter to hypothesis / expectation.

#Get some sense of effect size and CIs and create models for plotting

# Create model from above with treatment fitted

TreatMod <- glmer(Infected ~ Dose + Treatment + (1|Treatment/Line) + (1|Block), data = IFD, family = binomial)

# Create model with each line fitted

LinesMod <- glmer(Infected ~ Dose + Line + (1|Treatment/Line) + (1|Block), data = IFD, family = binomial)

# Get 'effects sizes' for treatment at a hypothetical average dose.

emmeans(object = TreatMod, specs = 'Treatment', by = 'Dose')

inv.logit(c(-2.26, -3.08, -1.446, -1.23, -1.99, -0.467))

# Plot infection data to get a sense of this graphically

# Plot each line individually

par(mfrow=c(2,max(c(NumE, NumL))), mar = c(6.0,6,2.5,2))

for(Line in sort(Lines)){
  
  
  Seb <- IFD[which(IFD$Line==Line),]
  
  if(unique(Seb$Treatment) == 'Early'){
    
    CT <- 'Early'
    
    PlotCol <- ColE
    
  }
  
  if(unique(Seb$Treatment) == 'Late'){
    
    CT <- 'Late'
    
    PlotCol <- ColL
    
  }
  
  CB <- unique(IFD$Block[which(IFD$Line == Line)])
  
  plot(x=seq(from = min(IFD$Dose), to = max(IFD$Dose), by = 0.1), 
       y = seq(from = 0, to = 1, length.out = length(seq(from = seq(from = min(IFD$Dose), to = max(IFD$Dose), by = 0.1)))), 
       type='n', xlab = expression(Log[10](Dose) + 1), ylab = 'Proportion Infected', main = Line,
       cex.main = 1.75, cex.lab = 1.635, cex.axis = 1.635)
  
  PlotData <- matrix(NA, ncol = 2, nrow = length(Doses))
  
  PlotData[,1] <- sort(unique(Seb$Dose))
  
  for(Dose in 1:length(unique(Seb$Dose))){
    
    PlotData[which(PlotData[,1] == (unique(Seb$Dose)[Dose])),2] <- sum(Seb$Infected[which(Seb$Dose == (unique(Seb$Dose)[Dose]))]) / NROW(Seb$Infected[which(Seb$Dose == (unique(Seb$Dose)[Dose]))])

  }
  
  points(x = PlotData[,1], y = PlotData[,2], type = 'p', col = PlotCol, lwd = 2 , pch = 4, cex=1.75)
  
  GV1 <- seq(from = min (IFD$Dose), to = max (IFD$Dose), length.out = 3000)
  GV2 <- rep(CT, time = NROW(GV1))
  GV3 <- rep(Line, time = NROW(GV1))
  GV4 <- rep('A', time = NROW(GV1))
  
  FakeData <- data.frame(Treatment = GV2, Dose = GV1, Line = GV3, Block = GV4)
  
  GV0 <- inv.logit(as.numeric(predict(TreatMod, FakeData, re.form = ~ (1|Treatment/Line)+ (1|Block) )))
  points(x = GV1, y = GV0 , type = 'l', lwd = 3 ,col = Transpa(PlotCol, 45), pch = 20, cex=1.75)
  
  
}


# Plot the 'overall effect of treatment'
par(mfrow=c(1,1), mar = c(6.5,6,2,2))

plot(x=seq(from = min(IFD$Dose), to = max(IFD$Dose), by = 0.1), 
     y = seq(from = 0, to = 1, length.out = length(seq(from = seq(from = min(IFD$Dose), to = max(IFD$Dose), by = 0.1)))), 
     type='n', xlab = expression(Log[10](Dose + 1)), ylab = ' Predicted Proportion Infected',
     cex.main = 1.75, cex.lab = 1.635, cex.axis = 1.635)

for( Treatment in sort(unique(LHD$Treatment))){
  
  if(Treatment == 'Early'){
    
    PlotCol <- ColE
    
  }
  
  if(Treatment == 'Late'){
    
    PlotCol <- ColL
    
  }
  
  GV1 <- seq(from = min (IFD$Dose), to = max (IFD$Dose), length.out = 3000)
  GV0 <- rep(Treatment, time = NROW(GV1))
  
  FakeData <- data.frame(Treatment = GV0, Dose = GV1)
  
  GV2 <- inv.logit(as.numeric(predict(TreatMod, FakeData, re.form = ~ 0)))
  
  points(x = GV1, y = GV2 , type = 'l', lwd = 4 ,col = PlotCol, pch = 20, cex=1.75)
  
}

legend( x = 0.3, y = 1, c('Early Selected', 'Late Selected'), col = c(ColE, ColL), pch = 19, border = NA, box.col= 'white', cex = 1.65)

## Combine these graphs into one composite figure for better publication clarity.



par(mfrow=c(1,1), mar = c(6.5,6,2,2))

plot(x=seq(from = min(IFD$Dose), to = max(IFD$Dose), by = 0.1), 
     y = seq(from = 0, to = 1, length.out = length(seq(from = seq(from = min(IFD$Dose), to = max(IFD$Dose), by = 0.1)))), 
     type='n', xlab = expression(Log[10](Dose) + 1), ylab = 'Proportion Infected', main = NULL,
     cex.main = 1.75, cex.lab = 1.635, cex.axis = 1.635)

for(Line in sort(Lines)){
  
  
  Seb <- IFD[which(IFD$Line==Line),]
  
  if(unique(Seb$Treatment) == 'Early'){
    
    CT <- 'Early'
    
    PlotCol <- ColE
    
  }
  
  if(unique(Seb$Treatment) == 'Late'){
    
    CT <- 'Late'
    
    PlotCol <- ColL
    
  }
  
  CB <- unique(IFD$Block[which(IFD$Line == Line)])
  
  PlotData <- matrix(NA, ncol = 2, nrow = length(Doses))
  
  PlotData[,1] <- sort(unique(Seb$Dose))
  
  for(Dose in 1:length(unique(Seb$Dose))){
    
    PlotData[which(PlotData[,1] == (unique(Seb$Dose)[Dose])),2] <- sum(Seb$Infected[which(Seb$Dose == (unique(Seb$Dose)[Dose]))]) / NROW(Seb$Infected[which(Seb$Dose == (unique(Seb$Dose)[Dose]))])
    
  }
  
  points(x = jitter(PlotData[,1], factor = 0.3), y = PlotData[,2], type = 'p', col = PlotCol, lwd = 2 , pch = 4, cex=1.75)
  
  GV1 <- seq(from = min (IFD$Dose), to = max (IFD$Dose), length.out = 3000)
  GV2 <- rep(CT, time = NROW(GV1))
  GV3 <- rep(Line, time = NROW(GV1))
  GV4 <- rep('A', time = NROW(GV1))
  
  FakeData <- data.frame(Treatment = GV2, Dose = GV1, Line = GV3, Block = GV4)
  
  GV0 <- inv.logit(as.numeric(predict(TreatMod, FakeData, re.form = ~ (1|Treatment/Line)+ (1|Block) )))
  points(x = GV1, y = GV0 , type = 'l', lwd = 3 ,col = Transpa(PlotCol, 45), pch = 20, cex=1.75)
  
  
}

for( Treatment in sort(unique(LHD$Treatment))){
  
  if(Treatment == 'Early'){
    
    PlotCol <- ColE
    
  }
  
  if(Treatment == 'Late'){
    
    PlotCol <- ColL
    
  }
  
  GV1 <- seq(from = min (IFD$Dose), to = max (IFD$Dose), length.out = 3000)
  GV0 <- rep(Treatment, time = NROW(GV1))
  
  FakeData <- data.frame(Treatment = GV0, Dose = GV1)
  
  GV2 <- inv.logit(as.numeric(predict(TreatMod, FakeData, re.form = ~ 0)))
  
  points(x = GV1, y = GV2 , type = 'l', lty = 'dashed' , lwd = 5 ,col = PlotCol, pch = 20, cex=1.75)
  
}

legend( x = 0.3, y = 1, c('Early Selected', 'Late Selected'), col = c(ColE, ColL), pch = 19, border = NA, box.col= 'white', cex = 1.65)

##############################

# Now undertake more detailed comparison to see if this finding (counter to expectations) exists using a 'trade-off-like' plot
# Let's do this with development time & growth rate


# Plot development time against infection susceptibility

# Need the dev. time model.

DTMM<- glmer(Days ~ Treatment + (1|Treatment/Line) + (1|Block), data = LHD, family = poisson)


LineData <- data.frame(Line = as.character(Lines))
LineData$DT <- NA
LineData$IC <- NA

for(J in 1:NROW(LineData)){
  
  FakeData <- data.frame(
    Treatment = unique(IFD$Treatment[which(IFD$Line == LineData$Line[J])]),
    Dose = max(IFD$Dose), 
    Line = LineData$Line[J], 
    #Block = unique(IFD$Block[which(IFD$Line == LineData$Line[J])])
    Block = 'A'
  )
  
  LineData$DT[J] <- exp(as.numeric(predict(DTMM, FakeData, re.form = ~ (1|Treatment/Line)+ (1|Block))))
  
  LineData$IC[J] <- inv.logit(as.numeric(predict(TreatMod, FakeData, re.form = ~ (1|Treatment/Line)+ (1|Block) )))
  
}


ColV <- rep(NA, length.out = NROW(LineData))

TempV <- unlist(sapply(X = LineData$Line, FUN = regexpr, pattern = 'E')) > 0

ColV[which(TempV)] <- ColE
ColV[which(!TempV)] <- ColL

rm(TempV)

LinMod <- lm(LineData$IC ~ LineData$DT)

PV1 <- seq(from = min(LineData$DT)*0.95, to = max(LineData$DT) *1.05, length.out = 100 )
PV2 <- seq(from = min(LineData$IC)*0.95, to = max(LineData$IC)* 1.05, length.out = 100 )

plot(PV2 ~ PV1, ylab = 'Susceptibility', xlab = expression(paste('Development Time (Days)', sep='')), 
     type = 'n', cex = 2, cex.lab = 1.5, cex.axis = 1.5)

points(LineData$IC ~ LineData$DT,
       pch = 4, col = ColV, cex = 1.5, lwd = 4)

V1 <- seq(from = par('usr')[1]*1.1, to = par('usr')[2]*0.95, length.out = 100)
V2 <- coef(LinMod)[[1]] + (V1*coef(LinMod)[[2]])

points(V2 ~ V1, type = 'l', lty = 2, lwd = 2, col = 'gray1')

legend( x = 12.5, y = 0.775, c('Early Selected', 'Late Selected'), col = c(ColE, ColL), 
        pch = 19, border = NA, box.col= 'white', cex = 1.2)


cor.test(LineData$DT, LineData$IC, method = "pearson")

# No good evidence that we can resolve this finding of faster growth = higher resistance at
# the level of the individual lines.


####
# Compare to populations when (approximately) assayed at start of experiment. See main text.
####

LHDS <- read.csv(file = '~/LHDStart.csv' , header=T)
IFDS <- read.csv(file = '~/IFDStart.csv' , header=T)

# Transform dose into log+1 space
IFDS$Dose <- log(IFDS$Dose+1)

#Life history analysis, see above, except treatment doesn't now matter as selection hasn't begun

# Clean out damaged larvae (Mass = NA) from life history assats
LHDS <- na.exclude(LHDS)

# Examine growth rate, calculated from mass and development time

LHDS$GRate <- LHDS$Mass/LHDS$Days

summary(aov(GRate ~ Line, data = LHDS))

# Lines differ in growth rate (but caveat: all assayed on different days)

# See if there is any evidence that it introduced a treatment bias:
mixed(GRate ~ Treatment + (1|Treatment/Line) , data = LHDS, method = 'LRT')

#No evidence for a treatment bias at set-up for larval growth rate

#Plot.

par(mfrow = c(1,1), cex.axis = 1.25, cex.lab = 1.4, mar = c(5,5,2,2))

boxplot(LHDS$GRate ~ LHDS$Line, 
        ylab = expression(paste('Growth Rate (mg  ',day^-1,')', sep='')), 
        xlab = 'Line', notch = F, varwidth = F, lwd = 2, pch = 19,
        col = c(rep(Transpa(ColE, Opac), times = NumE), rep(Transpa(ColL, Opac), time = NumL)), border = c(rep(ColE, times = NumE), rep(ColL, time = NumL)))

legend( x = 0.5, y = 0.875, c('Early Selected', 'Late Selected'), 
        col = c(ColE, ColL), pch = 19, border = NA, box.col= 'white', cex = 1)


## Test for changes in resistance in starting population

IFDS$Infected <- cbind(IFDS$I, IFDS$NI)

#Check for susceptibility of starting populations at generation = 1

mixed(Infected ~ Dose + Line + Dose:Line + (1|Treatment/Line) , data = IFDS, family = binomial, method = 'LRT')

summary(glmer(Infected ~ Dose + Line + Dose:Line + (1|Treatment/Line) , data = IFDS, family = binomial))

# Some differences in starting population apparent susceptibilities, but again, caveat of different assay days
# However for the data we have, the populations destined to be late-selected are all more resistant

# We can show that as an overall effect here

mixed(Infected ~ Dose + Treatment + (1|Treatment/Line) , data = IFDS, family = binomial, method = 'LRT')

summary(glmer(Infected ~ Dose + Treatment - 1 + (1|Treatment/Line) , data = IFDS, family = binomial))

# Late populations started off less susceptible (more resistant)

# Compare change in rank-order of susceptibility

# Get model
StartI <- glmer(Infected ~ Dose + Line + Dose:Line + (1|Treatment/Line) , data = IFDS, family = binomial)

IFRanks <- data.frame(Line = unique(IFDS$Line))
IFRanks$StartMS <- NA
IFRanks$EndMS <- NA

for(J in 1:NROW(unique(IFDS$Line))){
  
  FakeData <- data.frame(
    Treatment = unique(IFDS$Treatment[which(IFDS$Line == IFRanks$Line[J])]),
    Dose = median(IFDS$Dose), 
    Line = IFRanks$Line[J] 
  )
  
  IFRanks$StartMS[J] <- inv.logit(as.numeric(predict(StartI, FakeData, re.form = ~ (1|Treatment/Line))))
  
  IFRanks$EndMS[J] <- LineData$IC[which(LineData$Line == IFRanks$Line[J])]
  
  
}

IFRanks$StartR <- rank(IFRanks$StartMS)
IFRanks$EndR <- rank(IFRanks$EndMS)

IFRanks$Treatment <- sort(rep(c('Early', 'Late'), 4))

IFRanks
# Above table used to make Fig. S2 (image made 'by hand')

# Test whether selection regime alters rank change

IFRT <- data.frame(subj = IFRanks$Line, group = IFRanks$Treatment, TP1 = IFRanks$StartMS, TP2 = IFRanks$EndMS)
(npIntFactRep(IFRT))$RegularRanks

####
# Compare main experiment to populations re-assayed from Bartlett et al 2018 assayed in same lab with same virus stock
####

LHDIB <- read.csv(file = '~InbredLHD.csv' , header=T)
IFDIB <- read.csv(file = '~InbredIFD.csv' , header=T)

IFDIB$Infected <- cbind(IFDIB$I, IFDIB$NI)

IFDIB$Dose <- log10(IFDIB$Dose+1)

LHDIB <- na.exclude(LHDIB)

LHDIB$GRate <- LHDIB$Mass/LHDIB$Days

IIMod <- glm(Infected ~ Line + Dose + Line:Dose, data = IFDIB , family = binomial)
anova(IIMod, test = 'Chisq')
IIMod <- glm(Infected ~ Line + Dose, data = IFDIB , family = binomial)

head(LineData)

IBLD <- data.frame(Line = unique(LHDIB$Line))
IBLD$DT <- NA
IBLD$IC <- NA

for(J in 1:NROW(IBLD)){
  
  IBLD$DT[J] <- mean(LHDIB$Days[which(LHDIB$Line == IBLD$Line[J])])

  FakeData <- data.frame(
    Dose = max(IFDIB$Dose), 
    Line = IBLD$Line[J]
  )
  
  IBLD$IC[J] <- inv.logit(as.numeric(predict(IIMod, FakeData)))
  
}

ComPlot <- rbind(LineData, IBLD)

CompCols <- c(ColV, rep('green4', 3))

PV1 <- seq(from = min(ComPlot$DT)*0.95, to = max(ComPlot$DT) *1.05, length.out = 100 )
PV2 <- seq(from = min(ComPlot$IC)*1.05, to = max(ComPlot$IC)* 1.05, length.out = 100 )

plot(PV2 ~ PV1, ylab = 'Susceptibility', xlab = expression(paste('Development Time (Days)', sep='')), 
     type = 'n', cex = 2, cex.lab = 1.5, cex.axis = 1.5)

points(ComPlot$IC ~ ComPlot$DT,
       pch = 4, col = CompCols, cex = 1.5, lwd = 4)



legend( x = 25, y = 0.95, c('Early Selected', 'Late Selected', 'Inbred'), col = c(ColE, ColL, 'green4'), 
        pch = 19, border = NA, box.col= 'white', cex = 1.2)






