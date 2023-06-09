#### Working with the data ####
# Setup --------------------------------------------------------------------
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(polynom)

setwd("") # ENTER PATHNAME FOR LOCATION OF CSV 
full.dataframe <- read.csv("Salmon_ForceGape_Data.csv", header=T, skip=1)

# Create dataframes --------------------------------------------------------------------
#sp.G.F.data
sp.G.F.data <- full.dataframe[1:63, c("Fish.ID","S1.Gape.Normalized", "S1.Force.Normalized")] #includes all specimens
sp.G.F.data$Fish.ID <- sp.G.F.data$Fish.ID %>% as.factor()
specimens <- levels(sp.G.F.data$Fish.ID) 

#sp.G.F.data.woP3
sp.G.F.data.woP3 <-sp.G.F.data[1:58,] #excludes Pink 3
sp.G.F.data.woP3$Fish.ID <- sp.G.F.data.woP3$Fish.ID %>% as.character() %>% as.factor() #converting from character to factor refreshes the levels of Fish.ID. if we hadn't done this, it would have included Pink 3 even though it is trimmed from sp.G.F.data.woP3
specimens.woP3 <- levels(sp.G.F.data.woP3$Fish.ID) #check specimens.woP3 to make sure Pink 3 isn't included

#Be aware that sp.G.F.curve, sp.GHL.FRaw.curve, sp.fatigue.data, etc. are dataframes that are created later in the code and extract their data from the full.dataframe 

# Re-Normalizing regression curves to their max force --------------------------------------------------------------------
#When we normalize the raw force values to the maximum force for each individual and then calculate the polynomial regression for the normalized force-gape curves, the peaks of the polynomial regressions are not perfectly at 1 on the y axis (they range from .98 to 1.04) because the curves bent slightly above or below the data points.
#We needed the polynomial regression curves to peak exactly at 1 in order to do the analysis comparing curve widths (because those analyses are calculated at relative positions along the y axis). 
#To correct for this, we divided the normalized forces by the maximum of the polynomial regression (and called it S1.Force.ReNormalized), and then recomputed the polynomial regressions.

#Uses and alters sp.G.F.data and alters sp.G.F.data.woP3

predict.pts <- data.frame(S1.Gape.Normalized = seq(-0.5, 1.1, by = 0.005))

curve.Fmax.M <- matrix(nrow = length(specimens), ncol=2)

sp.G.F.data$S1.Force.ReNormalized <- 0 #I just need to create a new column

for(i in 1:length(specimens)) {
  #Generates predicted force values for the given seq of gape points.
  ithModel <- lm(S1.Force.Normalized ~ poly(S1.Gape.Normalized, 3, raw = T), data = sp.G.F.data[grep(specimens[i], sp.G.F.data$Fish.ID),]) #took this from modelsM
  predict.pts[,i+1] <- predict(ithModel, newdata = predict.pts)
  temp <- specimens[i]
  colnames(predict.pts)[i+1] <- paste0('S1.Force.Predicted.', temp) #[,1] is the seq; [,i+1] is the specimen's column
  
  #fills the curve.Fmax.M matrix with the specimen names and the value of their maximum force (which is close to 1, but not exact)
  curve.Fmax.M[i,1] <- specimens[i]
  predict.maxF.i <- predict.pts[which.max(predict.pts[,i+1]), i+1] #finds max force value. which.max give the index, i.e. row, and i+1 is the specimen's column of S1.Force.Predicted values   
  curve.Fmax.M[i,2] <- as.numeric(predict.maxF.i)
  
  #goes fish by fish and fills in values for S1.Force.ReNormalized (by dividing S1.Force.Normalized by the maximum of its S1.Force.Normalized polynomial fit curve)
  sp.G.F.data$S1.Force.ReNormalized[grep(specimens[i], sp.G.F.data$Fish.ID)] <- sp.G.F.data$S1.Force.Normalized[grep(specimens[i], sp.G.F.data$Fish.ID)] / as.numeric(curve.Fmax.M[i,2])
}

sp.G.F.data.woP3$S1.Force.ReNormalized <- sp.G.F.data$S1.Force.ReNormalized[grep("^Pink 3", sp.G.F.data$Fish.ID, invert=T)] 

#Now sp.G.F.data and sp.G.F.data.woP3 have a column S1.Force.ReNormalized, and curves that will be fit to that data will peak at 1


# creating a matrix of the F-G polynomial fit equations (and more) --------------------------------------------------------------------
#Creates a matrix that stores the fish ID, the coefficients of its polynomial regression, a variable (func#) that calls the complete equation, the optimal gape of the regression, the r-squared value, and the p-value for the polynomial regression 
#Uses sp.G.F.data

modelsM <- matrix(nrow = length(specimens), ncol=9)
for(i in 1:length(specimens)) { 
  modelsM[i,1] <- specimens[i]
  ithModel <- lm(S1.Force.ReNormalized ~ poly(S1.Gape.Normalized, 3, raw = T), data = sp.G.F.data[grep(specimens[i], sp.G.F.data$Fish.ID),]) 
  modelsM[i,2:5] <- coefficients(ithModel) 
  coE1 <- as.numeric(modelsM[i,2])
  coE2 <- as.numeric(modelsM[i,3])
  coE3 <- as.numeric(modelsM[i,4])
  coE4 <- as.numeric(modelsM[i,5])
  temp <- paste0('function(x) {',coE1,'+ x*',coE2,' + (x^2)*',coE3,' + (x^3)*',coE4,'}') #constructs the equation
  assign(paste0('func', i), eval(parse(text=temp))) #creates variables (func[1:10]) that store each of the eqn.s
  modelsM[i,6] <- paste0('func', i) #[,6] contains variables func[1:10], not the eqn.s themselves
  modelsM[i,7] <- solve(deriv(polynomial(c(coE1, coE2,coE3,coE4))))[predict(deriv(deriv(polynomial(c(coE1, coE2, coE3, coE4)))), solve(deriv(polynomial(c(coE1, coE2,coE3,coE4))))) < 0] #[,7] contains optimal gape (i.e. x-value of local maximum)
  modelsM[i,8] <- summary(ithModel)$r.squared 
  modelsM[i,9] <- pf(summary(ithModel)$fstatistic[1], summary(ithModel)$fstatistic[2], summary(ithModel)$fstatistic[3], lower.tail = F) #[,9] stores the calculated p-value.
}

#### Rate of muscle weakening ####

#NB: There are two columns in the Salmon_ForceGapeData.csv (S2.Fatigue.relative.to.S1 and S2.Time.after.end.of.initial.series)
#that were calculated by hand. We recreate these calculations in the code below, in the spirit of transparency and reproducibility. 

#The dataframe we originally read in (full.dataframe) contains, among other columns, tidy observations for the first force-gape curve 
#and tidy observations for the second force-gape curve, essentially bound together column-wise.
#To easily calculate muscle weakening, we want to convert this to a wide dataframe describing these two series, where each row corresponds to a unique pair of fish ID and gape.

#We begin by copying two separate dataframes from the original full dataframe, so they can be properly merged together.
sp.fatigue.S1 <- full.dataframe[1:63, c("Fish.ID", "S1.Gape..cm.", "S1.Force..N.", "S1.Time..min.")]
#We need to join by the Fish.ID and gape variables, so we need those variables to have the same names in both tables.
colnames(sp.fatigue.S1)[colnames(sp.fatigue.S1)=="S1.Gape..cm."] <- "Gape..cm." 

sp.fatigue.S2 <- full.dataframe[1:63, c("Fish.ID", "S2.Gape..cm.", "S2.Force..N.", "S2.Time..min.")]
colnames(sp.fatigue.S2)[colnames(sp.fatigue.S2)=="S2.Gape..cm."] <- "Gape..cm."
#We also need gape and force variables to be numeric.
sp.fatigue.S2$Gape..cm. <- sp.fatigue.S2$Gape..cm. %>% as.character() %>% as.numeric()
sp.fatigue.S2$S2.Force..N. <- sp.fatigue.S2$S2.Force..N. %>% as.character() %>% as.numeric()

#Now we can join them, left joining to keep rows where the first force-gape series has observations but the second series does not.
sp.fatigue.wide <- left_join(sp.fatigue.S1,sp.fatigue.S2, by=c("Fish.ID","Gape..cm."))

#We add a new variable indicating the time when the first force-gape series was completed (a property of each unique fish).
sp.fatigue.wide <- sp.fatigue.wide %>% group_by(Fish.ID) %>% mutate(S1.EndTime..min. = max(S1.Time..min.))

#Due to human error, one fish was measured at the same gape twice in succession. Here, we average over the forces and times of this duplicate measurement.
sp.fatigue.wide.remDup <- sp.fatigue.wide %>% group_by(Fish.ID, Gape..cm.) %>% summarise_all(funs(mean))

#For each gape size we calculate the difference in force production (in Newtons). Negative values correspond to decreases in force.
#Each difference in force production corresponds to a duration, which we calculated as the minutes between the end of the first force-gape curve to the time of the second force measurement.
sp.fatigue.wide.remDup <- sp.fatigue.wide.remDup %>% mutate(S1.S2.ForceDiff..N. = S2.Force..N.- S1.Force..N., Time.After.S1..min. = S2.Time..min. - S1.EndTime..min.)

#Before analyzing the muscle weakening we just calculated, we need to remove observations where the second series force (S2.Force..N.) is zero or missing,
#and we need to remove observations where the force difference is 0 (weakening has not yet begun). 
sp.fatigue.wide.nonzero <- sp.fatigue.wide.remDup %>% filter(S1.S2.ForceDiff..N. != 0, !is.na(S1.S2.ForceDiff..N.), S2.Force..N. != 0)


#We can now calculate the x-intercept of the linear regressions (decrease in force over time) for each fish. 
xinterceptsM <- matrix(nrow = 2, ncol = length(specimens)) 
for(i in 1:length(specimens)) { 
  xinterceptsM[1,i] <- specimens[i] #first row is fish IDs
  #stores the first coefficient (the y-intercept) of the linear regression between muscle weakening and time for each fish
  yintercept.i <- coefficients(lm (S1.S2.ForceDiff..N. ~ Time.After.S1..min., data = sp.fatigue.wide.nonzero %>% filter(Fish.ID == specimens[i])))[1]
  #stores the second coefficient (the slope) of the linear regression between muscle weakening and time for each fish
  slope.i <- coefficients(lm (S1.S2.ForceDiff..N. ~ Time.After.S1..min., data = sp.fatigue.wide.nonzero %>% filter(Fish.ID == specimens[i])))[2]
  #calculates the x-intercept of the linear regression between muscle weakening and time for each fish
  xinterceptsM[2,i] <- (yintercept.i*-1)/slope.i 
}
xinterceptsM

#### Stats ####
# MW-U test: KP scaled optimal gapes; p = 0.02381 --------------------------------------------------------------------
sp.g0.data <- data.frame(mM.Fish.ID = modelsM[,1], mM.g0 = (modelsM[,7])) #modelsM [,1] contains the fish ID, modelsM [,7] contains the optimal gape values #P3 is excluded in the code for the wilcox test below

#two-tailed, unpaired Mann-Whitney test
#compares Kings (n=3) to Pinks (n=6), excluding P3
wilcox.test(as.numeric(as.vector(sp.g0.data$mM.g0[grep("^King", sp.g0.data$mM.Fish.ID)])), as.numeric(as.vector(sp.g0.data$mM.g0[grep("^King|Pink 3", sp.g0.data$mM.Fish.ID, invert=T)])), alternative="two.sided", paired=F) 


# MW-U test: KP head lengths; p = 0.3682 --------------------------------------------------------------------
sp.hl.data <- full.dataframe[!is.na(full.dataframe$Head.Length..cm.), c("Fish.ID","Head.Length..cm.")] #filters the rows (aren't NA for head length)
sp.hl.data <- sp.hl.data[1:8,] #includes P1, P2, K3, K4, P5, P6, P8, K2 (excluding Pink 3). I didn't record the head length of P7, so that is missing here. 

avg.hl.pink <- mean (sp.hl.data$Head.Length..cm.[grep("^Pink",sp.hl.data$Fish.ID)]) #Mean head length of the Pink salmon is 11.14
avg.hl.king <- mean (sp.hl.data$Head.Length..cm.[grep("^King",sp.hl.data$Fish.ID)]) #Mean head length of the King salmon is 10.63

#two-tailed, unpaired Wilcoxon-Mann-Whitney test
#compares Kings (n=3) to Pinks (n=6), excluding P3
wilcox.test(sp.hl.data$Head.Length..cm.[grep("^King", sp.hl.data$Fish.ID)], sp.hl.data$Head.Length..cm.[grep("^Pink",sp.hl.data$Fish.ID)], alternative="two.sided", paired=F, exact=F) 
#wilcox.test assumes that the data is continuous, and exact=F corrects that here and allows for comparison of values with ties (i.e. two of the values for Pink are the same)

# MW-U test: KP absolute max gapes; p = 0.3621 --------------------------------------------------------------------
sp.gm.data <- full.dataframe[1:63, c("Fish.ID","S1.Gape..cm.")] #includes all specimens, but P3 is excluded in the code for the wilcox test

sp.gm.data$Fish.ID <- sp.gm.data$Fish.ID %>% as.factor() #convert to factor so that we can use levels()
specimens <- levels(sp.gm.data$Fish.ID)

gmM <- matrix(nrow = length(specimens), ncol=2)
for(i in 1:length(specimens)) {
  gmM[i,1] <- specimens[i]
  gmM[i,2] <- max(sp.gm.data$S1.Gape..cm.[grep(as.character(specimens[i]), sp.gm.data$Fish.ID)])
}
sp.gmM.data <- data.frame(gmM.Fish.ID = gmM[,1], gmM.gm = (gmM[,2])) #converts gmM to a dataframe that I can work with more easily

#two-tailed, unpaired Mann-Whitney test
#compares Kings (n=3) to Pinks (n=6), excluding P3
wilcox.test(as.numeric(as.vector(sp.gmM.data$gmM.gm[grep("^King", sp.gmM.data$gmM.Fish.ID)])), as.numeric(as.vector(sp.gmM.data$gmM.gm[grep("^King|Pink 3", sp.gmM.data$gmM.Fish.ID, invert=T)])), alternative="two.sided", paired=F, exact=F) 
#wilcox.test assumes that the data is continuous, and exact=F corrects that here and allows for comparison of values with ties (i.e. 2+ of the values are the same).

# MW-U test: KP absolute peak force; p = 0.2453 --------------------------------------------------------------------
sp.G.FRaw.data <- full.dataframe[1:58, c("Fish.ID","S1.Gape.Normalized", "S1.Force..N.")] #rows are trimmed to exclude Pink 3

sp.G.FRaw.data$Fish.ID <- sp.G.FRaw.data$Fish.ID %>% as.character() %>% as.factor() #converts to factor so that we can use levels()
specimens <- levels(sp.G.FRaw.data$Fish.ID) #Pink 3, Pink -, and Pink FT should be excluded
levels(sp.G.FRaw.data$Fish.ID)

sp.peakFRaw <- data.frame(Fish.ID = NA, G0 = NA, F0Raw = NA)

for(i in 1:length(specimens)) { 
  sp.peakFRaw[i,1] <- specimens[i]
  ithModel <- lm(S1.Force..N. ~ poly(S1.Gape.Normalized, 3, raw = T), data = sp.G.FRaw.data[grep(specimens[i], sp.G.FRaw.data$Fish.ID),]) #3rd order polynomial equation fit to normalized gape and non-normalized force
  coE1 <- as.numeric(coefficients(ithModel)[1])
  coE2 <- as.numeric(coefficients(ithModel)[2])
  coE3 <- as.numeric(coefficients(ithModel)[3])
  coE4 <- as.numeric(coefficients(ithModel)[4])
  eqn <- eval(parse(text=paste0('function(x) {',coE1,'+ x*',coE2,' + (x^2)*',coE3,' + (x^3)*',coE4,'}'))) #constructs the equation
  sp.peakFRaw[i,2] <- modelsM[i,7]
  sp.peakFRaw[i,3] <- eqn(as.numeric(sp.peakFRaw[i,2])) # evaluates function(x) for  x = optimal gape
}

sp.peakFRaw$G0 <- as.numeric(sp.peakFRaw$G0)
sp.peakFRaw$F0Raw <- as.numeric(sp.peakFRaw$F0Raw)

avg.fm.pink <- mean (sp.peakFRaw$F0Raw[grep("^P",sp.peakFRaw$Fish.ID)]) #4.097601
avg.fm.king <- mean (sp.peakFRaw$F0Raw[grep("^K",sp.peakFRaw$Fish.ID)]) #2.320697

#two-tailed, unpaired Mann-Whitney test
#compares Kings (n=3) to Pinks (n=6), excluding P3
wilcox.test(sp.peakFRaw$F0Raw[grep("^K",sp.peakFRaw$Fish.ID)], sp.peakFRaw$F0Raw[grep("^P",sp.peakFRaw$Fish.ID)], alternative="two.sided", paired=F, exact=F) 
#wilcox.test assumes that the data is continuous, and exact=F corrects that here and allows for comparison of values with ties (i.e. 2+ values are identical).

# MW-U test: KP lever ratio; p = 0.551 --------------------------------------------------------------------
sp.lr.data <- full.dataframe[!is.na(full.dataframe$Lever.Ratio..In.Out.), c("Fish.ID","Lever.Ratio..In.Out.")] #filters the rows (aren't NA for head length) and the columns (by name)
sp.lr.data <- sp.lr.data[1:8,] #includes P1, P2, K3, K4, P5, P6, P8, K2 (excluding Pink 3). I didn't record the lever ratio of P7, so that is missing here.

avg.lr.pink <- mean (sp.lr.data$Lever.Ratio..In.Out.[grep("^Pink",sp.lr.data$Fish.ID)]) #The mean lever ratio for Pink salmon is 0.228
avg.lr.king <- mean (sp.lr.data$Lever.Ratio..In.Out.[grep("^King",sp.lr.data$Fish.ID)]) #The mean lever ratio for King salmon is 0.237

#two-tailed, unpaired Wilcoxon-Mann-Whitney test
#compares Kings (n=3) to Pinks (n=6), excluding P3
wilcox.test(sp.lr.data$Lever.Ratio..In.Out.[grep("^King",sp.lr.data$Fish.ID)], sp.lr.data$Lever.Ratio..In.Out.[grep("^Pink",sp.lr.data$Fish.ID)], alternative="two.sided", paired=F, exact=F) 
#wilcox.test assumes that the data is continuous, and exact=F corrects that here and allows for comparison of values with ties (i.e. 2+ values are identical)

# MW-U test: KP muscle mass; p = 0.3711 --------------------------------------------------------------------
sp.mm.data <- full.dataframe[!is.na(full.dataframe$Muscle.mass..g.), c("Fish.ID","Muscle.mass..g.")] #filter the rows (aren't NA for head length) and the columns (by name) 
sp.mm.data <- sp.mm.data[1:8,] #includes P1, P2, K3, K4, P5, P6, P8, K2 (excluding Pink 3). I didn't record the muscle mass of P7, so that is missing here.

avg.mm.pink <- mean (sp.mm.data$Muscle.mass..g.[grep("^Pink",sp.mm.data$Fish.ID)]) #The mean muscle mass for Pink salmon is 2.94
avg.mm.king <- mean (sp.mm.data$Muscle.mass..g.[grep("^King",sp.mm.data$Fish.ID)]) #The mean muscle mass for King salmon is 2.66

#two-tailed, unpaired Wilcoxon-Mann-Whitney test
#compares Kings (n=3) to Pinks (n=6), excluding P3
wilcox.test(sp.mm.data$Muscle.mass..g.[grep("^King",sp.mm.data$Fish.ID)], sp.mm.data$Muscle.mass..g.[grep("^Pink",sp.mm.data$Fish.ID)], alternative="two.sided", paired=F, exact=F) 
#wilcox.test assumes that the data is continuous, and exact=F corrects that here and allows for comparison of values with ties (i.e. 2+ values are identical).

# comparing curve widths at 90% and 75% of P0 --------------------------------------------------------------------
#We calculated the widths of the polynomial regressions at 90% and 75% of peak force. 
#To do this, we generated a set of points on the polynomial regression and then located the x-

#Uses sp.G.F.data.woP3

#Creates new dataframes with 500 points (gape values range from -0.5 to 1.1)
predict.w.pts <- data.frame(S1.Gape.Normalized = seq(-0.5, 1.1, by = 0.005)) 

#Uses the polynomial regression of each fish to generates predicted force values for that individual (for the sequence of gape points in predict.w.pts)
for(i in 1:length(specimens.woP3)) { 
  ithModel <- lm(S1.Force.Normalized ~ poly(S1.Gape.Normalized, 3, raw = T), data = sp.G.F.data.woP3[grep(specimens.woP3[i], sp.G.F.data.woP3$Fish.ID),]) 
  predict.w.pts[,i+1] <- predict(ithModel, newdata = predict.w.pts) #predict.w.pts[,1] is the seq; predict.w.pts[,i+1] is the specimen's column
  temp <- specimens.woP3[i]
  colnames(predict.w.pts)[i+1] <- paste0('S1.Force.Re.Predicted.', temp)
}

#Creates a list (widthMatrices.list) with two elememts (named predict.width.90.M and predict.width.75.M)
pct.s <- c(0.90, 0.75)
widthMatrices.list <- vector("list", length = length(pct.s))
names(widthMatrices.list) <- c("predict.width.90.M", "predict.width.75.M")

#Makes each element in widthMatrices.list into a matrix that contains (for each individual fish) the fish ID, the maximum force, 90% or 75% of that maximum force, the gape and force values of the coordinate points (on the left and right sides of the curve) that are have force values closest to the calculated value of 90% or 75% of maximum force, and lastly, the difference in the gape values of those two points.  
for (j in 1:length(pct.s)) {
  widthMatrices.list[[j]] <- matrix(nrow = length(specimens.woP3), ncol=8) #widthMatrices.list[[j]] represents each of the predict.width.##.M matrices.
  for (i in 1:length(specimens.woP3)) { 
    widthMatrices.list[[j]][i,1] <- specimens.woP3[i]
    predict.maxF.i <- predict.w.pts[which.max(predict.w.pts[,i+1]), i+1] #finds max force value. which.max give the index, i.e. row, and i+1 is the specimen's column of S1.Force.Re.Predicted values.   
    widthMatrices.list[[j]][i,2] <- as.numeric(predict.maxF.i)
    widthMatrices.list[[j]][i,3] <- pct.s[j]*predict.maxF.i #max force times 0.90 or 0.75
    Lpt.idx.i <- which.min( abs( filter(predict.w.pts, predict.w.pts[,1] < .5, predict.w.pts[,i+1] > 0)[,i+1] - (pct.s[j]*predict.maxF.i))) #looks at left side of curve (i.e. where gape is less than 0.5) for where force is 90% or 75% of max and returns the index for the gape of that point
    Rpt.idx.i <- which.min( abs( filter(predict.w.pts, predict.w.pts[,1] > .5, predict.w.pts[,i+1] > 0)[,i+1] - (pct.s[j]*predict.maxF.i))) #looks at right side of curve (i.e. where gape is greater than 0.5) for where force is 90% or 75% of max and returns the index for the gape of that point
    widthMatrices.list[[j]][i,4] <- Lpt.i <- (filter(predict.w.pts, predict.w.pts[,1] < .5, predict.w.pts[,i+1] > 0)[,1])[Lpt.idx.i] #Lpt.i #uses the gape index to return the gape value
    widthMatrices.list[[j]][i,5] <- (filter(predict.w.pts, predict.w.pts[,1] < .5, predict.w.pts[,i+1] > 0)[,i+1])[Lpt.idx.i] #uses the gape index to return the force value
    widthMatrices.list[[j]][i,6] <- Rpt.i <- (filter(predict.w.pts, predict.w.pts[,1] > .5, predict.w.pts[,i+1] > 0)[,1])[Rpt.idx.i] #Rpt.i #uses the gape index to return the gape value
    widthMatrices.list[[j]][i,7] <- (filter(predict.w.pts, predict.w.pts[,1] > .5, predict.w.pts[,i+1] > 0)[,i+1])[Rpt.idx.i] #uses the gape index to return the force value
    widthMatrices.list[[j]][i,8] <- Rpt.i - Lpt.i 
  }
}

#widthMatrices.list
#widthMatrices.list[["predict.width.90.M"]]

# MW-U test: KP curve widths at 90% and 75%; p = 0.8955, p = 0.4367 --------------------------------------------------------------------
#creats a dataframe with the fish IDs, and the curve widths at 90% and 75% of maximum force
sp.w90.data <- data.frame(Fish.ID = widthMatrices.list$predict.width.90.M[,1] %>% as.factor(), Dif.90 = widthMatrices.list$predict.width.90.M[,8] %>% as.numeric()) 
sp.w75.data <- data.frame(Fish.ID = widthMatrices.list$predict.width.75.M[,1] %>% as.factor(), Dif.75 = widthMatrices.list$predict.width.75.M[,8] %>% as.numeric()) 

#two-tailed, unpaired Wilcoxon-Mann-Whitney test
#compares Kings (n=3) to Pinks (n=6), excluding P3
wilcox.test(sp.w90.data$Dif.90[grep("^King", sp.w90.data$Fish.ID)], sp.w90.data$Dif.90[grep("^Pink",sp.w90.data$Fish.ID)], alternative="two.sided", paired=F, exact=F) #p = 0.8955

wilcox.test(sp.w75.data$Dif.75[grep("^King", sp.w75.data$Fish.ID)], sp.w75.data$Dif.75[grep("^Pink",sp.w75.data$Fish.ID)], alternative="two.sided", paired=F, exact=F) #p = 0.4367
#wilcox.test assumes that the data is continuous, and exact=F corrects that here and allows for comparison of values with ties (i.e. two of the values for Pink are the same).

avg.w90.pink <- mean (sp.w90.data$Dif.90[grep("^Pink", sp.w90.data$Fish.ID)]) #The mean curve width at 90% of peak force for Pink salmon is 0.4408 (width is in terms of normalized gape)
avg.w90.king <- mean (sp.w90.data$Dif.90[grep("^King", sp.w90.data$Fish.ID)]) #The mean curve width at 90% of peak force for King salmon is 0.4083
avg.w90.pink
avg.w90.king
avg.w75.pink <- mean (sp.w75.data$Dif.75[grep("^Pink", sp.w75.data$Fish.ID)]) #The mean curve width at 75% of peak force for Pink salmon is 0.725
avg.w75.king <- mean (sp.w75.data$Dif.75[grep("^King", sp.w75.data$Fish.ID)]) #The mean curve width at 75% of peak force for King salmon is 0.65
avg.w75.pink
avg.w75.king
