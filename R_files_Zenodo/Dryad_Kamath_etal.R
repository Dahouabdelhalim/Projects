#### Table of contents ####
# This is the R Code that should re-create all the analyses and figures found in Kamath et al, 2020, Character displacement in the midst of background evolution in island populations of Anolis lizards: a spatiotemporal perspective

#If something is broken, or provides an answer inconsistent with what is reported in the paper, please email Yoel Stuart (ystuart@luc.edu).

# 0) Read in data for analysis 
# 1) Perch Height Analyses
# 2) Toepad analyses
# 3) New analyses added during review
# BESPOKE FUNCTIONS (run before size correcting)
# FIGURES
# end) Table of contents


#### 0) Read in data for analysis ####
#NB: Place the data files in a folder and generate a path to it. In the code below, replace "path.raw.data.southTwin" with the name of the path to your data folder.

#### 0) Read in data for analysis ####
#### 0a) Perch height data 2019 & 2010 ####
#2019
ph <- read.csv(paste(path.raw.data.southTwin, "perch.height.southTwin.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
ph.final <- ph[ph$dstrb == "no", c(2, 3, 4, 5, 7, 8, 9, 10, 11)] #dstrb == "no" subsets the data to behavioral measures that were not disturbed
ph.final <- ph.final[ph.final$sex == "m" | ph.final$sex == "f" , ] #just confirmed males and females
#2010
ph.2010 <- read.csv(paste(path.data.cleaned.southTwin, "101207_Mosquito Lagoon_Ecology Morphology_Master.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
ph.2010.final.tt <- ph.2010[ph.2010$sp == "carolinensis" & ph.2010$undist == 1 , ] #just undistrubed carolinensis
ph.2010.final.t <- ph.2010.final.tt[ph.2010.final.tt$sex == "m" | ph.2010.final.tt$sex == "f" , ]
islands.2010 <- c("Yang", "Channel", "Crescent", "Osprey", "South.Twin", "North.Twin", "Hornet") #subset of islands studied in 2019
ph.2010.final <- ph.2010.final.t[ph.2010.final.t$island %in% islands.2010, ]
# end 0a)

#### 0b) Toepad data 2019 & 2010 ####
#toepad data 2019
tp.2019.t <- read.csv(paste(path.data.cleaned.southTwin, "2019_toepad.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
tp.2019.t$liz.ID.univ <- substr(tp.2019.t$photoID, 1, 4)
tp.2019.t$right.left <- substr(tp.2019.t$photoID, 10, 10)
lizIDs <- unique(tp.2019.t$liz.ID.univ)
#average left and right
tp.2019.lr.av <- data.frame(matrix(NA, nrow = length(lizIDs), ncol = 3))
colnames(tp.2019.lr.av) <- c("liz.ID.univ", "toe.area.mm2.av", "lam.num.av")
for(i in 1:length(lizIDs)){
  focal.rows <- which(tp.2019.t$liz.ID.univ %in% lizIDs[i] )
  average.area <- mean(tp.2019.t$toe.area.mm2[focal.rows], na.rm = TRUE)
  average.lam <- mean(tp.2019.t$lamella.number[focal.rows], na.rm = TRUE)
  tp.2019.lr.av[i,1] <- lizIDs[i]
  tp.2019.lr.av[i,2] <- average.area
  tp.2019.lr.av[i,3] <- average.lam
}
#svl, mass, sex
sms.2019.t <- read.csv(paste(path.data.cleaned.southTwin, "2019_svl_mass.id.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# morph file
morph.2019 <- merge(sms.2019.t, tp.2019.lr.av, by = "liz.ID.univ")

#toepad data 2010
tp.2010.t <- read.csv(paste(path.raw.data.southTwin, "toes_male_female.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
tp.2010.tt <- tp.2010.t[tp.2010.t$morph.good.av == 1, c(2, 3, 6, 7,8, 17, 18)]
tp.2010.ttt <- tp.2010.tt[tp.2010.tt$island == "channel" | tp.2010.tt$island == "crescent" | tp.2010.tt$island == "hornet" | tp.2010.tt$island == "n.twin" | tp.2010.tt$island == "osprey" | tp.2010.tt$island == "s.twin" | tp.2010.tt$island == "yang", ]
tp.2010.ttt <- tp.2010.ttt[tp.2010.ttt$sex == "m" | tp.2010.ttt$sex == "f", ]
# end 0b)

##### 0c) toepad data size correction #####
#size correction. 
#Pool 2010 and 2019 for size correction
tp.2010.tttt <- tp.2010.ttt
tp.2010.tttt$year <- rep("2010", length(tp.2010.tttt$field.id))

morph.2019.t <- morph.2019
morph.2019.t$year <- rep("2019", length(morph.2019.t$liz.ID.univ))
morph.2019.t$island[which(morph.2019.t$island == "north.twin")] <- "n.twin"; morph.2019.t$island[which(morph.2019.t$island == "south.twin")] <- "s.twin"

pooled.for.size.correction <- data.frame(matrix(NA, nrow = (length(tp.2010.tttt$field.id) + length(morph.2019.t$liz.ID.univ)), ncol = 7))
names(pooled.for.size.correction) <- c("ID", "year", "svl", "area", "lam", "sex", "island")  
pooled.for.size.correction$ID <- c(tp.2010.tttt$field.id, morph.2019.t$liz.ID.univ)
pooled.for.size.correction$year <- c(tp.2010.tttt$year, morph.2019.t$year)
pooled.for.size.correction$svl <- c(tp.2010.tttt$svl, morph.2019.t$svl)
pooled.for.size.correction$sex <- c(tp.2010.tttt$sex, morph.2019.t$sex)
pooled.for.size.correction$lam <- c(tp.2010.tttt$count.av, morph.2019.t$lam.num.av)
pooled.for.size.correction$area <- c(tp.2010.tttt$area.mm2.av, morph.2019.t$toe.area.mm2.av)
pooled.for.size.correction$island <- c(tp.2010.tttt$island, morph.2019.t$island)

# Separate the sexes while pooling
#males
pooled.for.size.correction.m <- pooled.for.size.correction[pooled.for.size.correction$sex=="m" , ]
pooled.for.size.correction.m.sc <- f.sizecorrect.withoutsex.withyear(standard.length = pooled.for.size.correction.m$svl, vector.of.columns = c(4,5), dataframe.to.use = pooled.for.size.correction.m, watershed = pooled.for.size.correction.m$island, year = pooled.for.size.correction.m$year)
#females
pooled.for.size.correction.f <- pooled.for.size.correction[pooled.for.size.correction$sex=="f" , ]
pooled.for.size.correction.f.sc <- f.sizecorrect.withoutsex.withyear(standard.length = pooled.for.size.correction.f$svl, vector.of.columns = c(4,5), dataframe.to.use = pooled.for.size.correction.f, watershed = pooled.for.size.correction.f$island, year = pooled.for.size.correction.f$year)

#pool data
pooled.for.size.correction.sc <- merge(x = pooled.for.size.correction.m.sc, y = pooled.for.size.correction.f.sc, all = TRUE)
pooled.for.size.correction.sc.final <- pooled.for.size.correction.sc[ , 1:8]

morph.2019.sc <- pooled.for.size.correction.sc.final[pooled.for.size.correction.sc.final$year == 2019, ]
morph.2010.sc <- pooled.for.size.correction.sc.final[pooled.for.size.correction.sc.final$year == 2010, ]
# end 0c)

#### 0d) A. sagrei catch per unit effort 2019 ####
cpue <- read.csv(paste(path.raw.data.southTwin, "catch.effort.southTwin.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
cpue.t <- cpue[-which(cpue$date == "10-Jun-19"), ] #get rid of the bad weather day on 10 June. Yang
cpue.final <- cpue.t[, c(1, 2, 5, 6, 8)]
# end 0d)

# end 0)

##### 1) Perch Height Analyses #####
#### 1a) S. Twin 2010 vs. 2019. Did perch height change?####
stw.2010 <- ph.2010.final[ph.2010.final$island == "South.Twin",]
stw.2019 <- ph.final[ph.final$isle.code == "stw", ]
stw.both.years <- data.frame(c(stw.2010$perch.lizard, stw.2019$ph.liz.cm), c(rep("2010", length(stw.2010$perch.lizard)), rep("2019", length(stw.2019$ph.liz.cm))), c(stw.2010$sex, stw.2019$sex)) ; names(stw.both.years) <- c("phl", "year", "sex")

#### Question 2 models perch height ####
summary(model.1a.1 <- lm(stw.both.years$phl ~ stw.both.years$year + stw.both.years$sex))
#end 1a)

##### 1b) All islands, both years, linear model #####
pooled.ph.t <- ph.2010.final
pooled.ph.t$year <- rep(2010, length(pooled.ph.t$obs))
pooled.ph.t$island <- tolower(pooled.ph.t$island)
pooled.ph.tt <- ph.final
pooled.ph.tt$year <- rep(2019, length(pooled.ph.tt$liz.ID.univ))
pooled.ph.final <- data.frame(c(pooled.ph.t$island, pooled.ph.tt$island), c(pooled.ph.t$perch.lizard, pooled.ph.tt$ph.liz.cm), c(pooled.ph.t$year, pooled.ph.tt$year), c(pooled.ph.t$sex, pooled.ph.tt$sex), stringsAsFactors = FALSE); colnames(pooled.ph.final) <- c("island", "phl", "year", "sex")
#rename S.Twin so that it is alphabetically first and therefore is the baseline.
stws <- which(pooled.ph.final$island == "south.twin")
pooled.ph.final$island[stws] <- "a.south.twin"
##### Model for Table 6 , Table S1 #####
##### Model for Figure 1 #####
summary(model.1b.1 <- lm(pooled.ph.final$phl ~ as.character(pooled.ph.final$year) * pooled.ph.final$island + pooled.ph.final$sex))
anova(model.1b.1)
# end 1b)

##### 1c) Control islands, both years, linear mixed effects model #####
pooled.ph.final.lme <- pooled.ph.final[-(which(pooled.ph.final$island == "a.south.twin")), ]
pooled.ph.final.lme$category <- ifelse(pooled.ph.final.lme$island %in% c("hornet", "osprey", "crescent"), 0, 1)
##### Model for Table 3 #####
summary(model.1c.1 <- lme(phl ~ as.character(year) * as.character(category) + sex, data = pooled.ph.final.lme, random = ~1|island, na.action = "na.omit"))
# end 1c)

##### 1d) Effect of A. sagrei presence #####
##sagrei density##
isle.codes <- unique(cpue.final$isle.code)
mean.densities <- data.frame(matrix(NA, ncol = 1, nrow = length(unique(ph.final$isle.code)))); colnames(mean.densities) <- "density"
for(i in 1:length(isle.codes)){
  mean.densities[i, 1] <- mean(tapply(X = cpue.final$perch.height[cpue.final$isle.code == isle.codes[i]], INDEX = cpue.final$round[cpue.final$isle.code == isle.codes[i]], FUN = length))
  rownames(mean.densities)[i] <- isle.codes[i]
}
mean.densities[c(5:7), 1] <- 0; rownames(mean.densities)[5:7] <- c("osp", "cre", "hor")
ph.means.2019 <- tapply(X = ph.final$ph.liz.cm, INDEX = ph.final$island, FUN = mean, na.rm = TRUE)[c(7, 1, 6, 4, 5, 2, 3)]
density.vs.ph <- cbind(data.frame(ph.means.2019), mean.densities); colnames(density.vs.ph)[2] <- "mean.density"
##### Model for Figure 3 #####
summary(model.1d.1 <- lm(density.vs.ph$ph.means.2019 ~ density.vs.ph$mean.density))
# without uninvaded islands
summary(model.1d.2 <- lm(density.vs.ph$ph.means[c(1:4)] ~ density.vs.ph$mean.density[c(1:4)]))

##potential interactions##
##### Models for Table 9 #####
##### Models for Figure 4 ####
#2010
ph.Asonly.2010 <- ph.2010.final[is.na(ph.2010.final$as.int) == FALSE , ]
ph.Asonly.2010$as.binary <- ifelse(test=ph.Asonly.2010$as.int == "no", yes = 0, no = 1)
#summary(model.1d.3 <- lm(ph.Asonly.2010$perch.lizard ~ as.character(ph.Asonly.2010$as.binary) + ph.Asonly.2010$sex + ph.Asonly.2010$island))

model.1d.5 <- lme(perch.lizard ~ as.character(as.binary) + sex, data = ph.Asonly.2010, random = ~1|island, na.action = "na.omit")
summary(model.1d.5)

#2019
ph.Asonly.2019 <- ph.final[is.na(ph.final$as.int.sex) == FALSE, ]
ph.Asonly.2019$as.binary <- ifelse(test = ph.Asonly.2019$as.int.sex == "no", yes = 0, no = 1)
st <- which(ph.Asonly.2019$isle.code == "stw")
ph.Asonly.2019$isle.code[st] <- "ast"
#summary(model.1d.4 <- lm(ph.Asonly.2019$ph.liz.cm ~ as.character(ph.Asonly.2019$as.binary) + ph.Asonly.2019$sex + ph.Asonly.2019$isle.code))
#EtaSq((aov(model.1d.4)))

model.1d.6 <- lme(ph.liz.cm ~ as.character(as.binary) + sex, data = ph.Asonly.2019, random = ~1|island, na.action = "na.omit")
summary(model.1d.6)

# end 1d)
#end 1)


##### 2) Toepad analyses #####
#Toepad analysis uses size corrected toepad area. Not size corrected lamella number.
##### 2a) S. Twin 2010 vs. 2019. Toepad area and lamella number#####
stwin.2019.toe <- morph.2019.sc[morph.2019.sc$island == "s.twin", ] 
stwin.2010.toe <- morph.2010.sc[morph.2010.sc$island == "s.twin", ]
## Summary stats
#area
stwin.bothyears.area <- data.frame(c(stwin.2010.toe$area.sc, stwin.2019.toe$area.sc), c(rep("2010", length(stwin.2010.toe$area.sc)), rep("2019", length(stwin.2019.toe$area.sc))), c(stwin.2010.toe$sex, stwin.2019.toe$sex), c(stwin.2010.toe$svl, stwin.2019.toe$svl)) ; names(stwin.bothyears.area) <- c("area.sc", "year", "sex", "svl")
tapply(X = stwin.bothyears.area$area.sc, INDEX = as.character(stwin.bothyears.area$year), FUN = mean, na.rm = TRUE) - (tapply(X = stwin.bothyears.area$area.sc, INDEX = as.character(stwin.bothyears.area$year), FUN = sd, na.rm = TRUE)/sqrt(tapply(X = stwin.bothyears.area$area.sc, INDEX = as.character(stwin.bothyears.area$year), FUN = length)))
#lam num
stwin.bothyears.lam <- data.frame(c(stwin.2010.toe$lam, stwin.2019.toe$lam), c(rep("2010", length(stwin.2010.toe$lam)), rep("2019", length(stwin.2019.toe$lam))), c(stwin.2010.toe$sex, stwin.2019.toe$sex), c(stwin.2010.toe$svl, stwin.2019.toe$svl)) ; names(stwin.bothyears.lam) <- c("lam", "year", "sex", "svl")
tapply(X = stwin.bothyears.lam$lam, INDEX = as.character(stwin.bothyears.lam$year), FUN = mean, na.rm = TRUE) - (tapply(X = stwin.bothyears.lam$lam, INDEX = as.character(stwin.bothyears.lam$year), FUN = sd, na.rm = TRUE)/sqrt(tapply(X = stwin.bothyears.lam$lam, INDEX = as.character(stwin.bothyears.lam$year), FUN = length)))

#### Prediction 2 models toepads ####
summary(model.2a.1 <- lm(stwin.bothyears.area$area.sc ~ as.character(stwin.bothyears.area$year) + stwin.bothyears.area$sex))

summary(model.2a.2 <- lm(stwin.bothyears.lam$lam ~ as.character(stwin.bothyears.lam$year) + stwin.bothyears.lam$sex))
#end 2a)

##### 2b) Comparative analyses. Toepad area and lamella number #####
#build data frame that has all islands and combines both years.
pooled.for.lme <- pooled.for.size.correction.sc.final
pooled.for.lme$island[which(pooled.for.lme$island == "s.twin")] = "a.s.twin"
pooled.for.lme$category <- ifelse(pooled.for.lme$island %in% c("hornet", "osprey", "crescent"), yes = 0, no = 1)
pooled.for.lme$category[which(pooled.for.lme$island == "a.s.twin")] <- 2 

##### 2b.1 Linear Models with all islands#####
##### Model for Table 7, S2 #####
summary(model.2b.1 <- lm(area.sc ~ as.factor(year) * island + sex, data = pooled.for.lme))
anova(model.2b.1)

##### Model for Table 8, Table S3 #####
summary(model.2b.2 <- lm(lam ~ as.factor(year) * island + sex, data = pooled.for.lme))
anova(model.2b.2)

##### 2b.2 Linear Mixed Models with only control islands#####
##### Model for Table 4 #####
summary(model.2b.3 <- lme(area.sc ~ as.character(year) * as.character(category) + sex, data = pooled.for.lme[-which(pooled.for.lme$island == "a.s.twin"), ], random = ~1|island, na.action = "na.omit"))

##### Model for Table 5 #####
summary(model.2b.3 <- lme(lam ~ as.character(year) * as.character(category) + sex, data = pooled.for.lme[-which(pooled.for.lme$island == "a.s.twin"), ], random = ~1|island, na.action = "na.omit"))

# end 2b)
#end 2)

##### 3) New analyses added during review #####
##### 3a) Wilcoxon signed rank test. #####
#Does mismatch between toes and perch height improve over time
mean.2010.phl.for.SR <- tapply(X = pooled.ph.final$phl[pooled.ph.final$year==2010 & pooled.ph.final$island != "a.south.twin"], INDEX = pooled.ph.final$island[pooled.ph.final$year==2010 & pooled.ph.final$island != "a.south.twin"], FUN = mean, na.rm = TRUE)
mean.2010.area.for.SR <- tapply(X = pooled.for.lme$area.sc[pooled.for.lme$year==2010 & pooled.for.lme$island != "a.s.twin"], INDEX = pooled.for.lme$island[pooled.for.lme$year==2010 & pooled.for.lme$island != "a.s.twin"], FUN = mean, na.rm = TRUE)
mean.2010.lam.for.SR <- tapply(X = pooled.for.lme$lam[pooled.for.lme$year==2010 & pooled.for.lme$island != "a.s.twin"], INDEX = pooled.for.lme$island[pooled.for.lme$year==2010 & pooled.for.lme$island != "a.s.twin"], FUN = mean, na.rm = TRUE)

mean.2019.phl.for.SR <- tapply(X = pooled.ph.final$phl[pooled.ph.final$year==2019 & pooled.ph.final$island != "a.south.twin"], INDEX = pooled.ph.final$island[pooled.ph.final$year==2019 & pooled.ph.final$island != "a.south.twin"], FUN = mean, na.rm = TRUE)
mean.2019.area.for.SR <- tapply(X = pooled.for.lme$area.sc[pooled.for.lme$year==2019 & pooled.for.lme$island != "a.s.twin"], INDEX = pooled.for.lme$island[pooled.for.lme$year==2019 & pooled.for.lme$island != "a.s.twin"], FUN = mean, na.rm = TRUE)
mean.2019.lam.for.SR <- tapply(X = pooled.for.lme$lam[pooled.for.lme$year==2019 & pooled.for.lme$island != "a.s.twin"], INDEX = pooled.for.lme$island[pooled.for.lme$year==2019 & pooled.for.lme$island != "a.s.twin"], FUN = mean, na.rm = TRUE)

#phl vs. lam
#2010
wilcox.test(x = mean.2010.phl.for.SR, y = mean.2010.area.for.SR, paired = TRUE)
wilcox.test(x = mean.2010.phl.for.SR, y = mean.2010.lam.for.SR, paired = TRUE)

#2019
wilcox.test(x = mean.2019.phl.for.SR, y = mean.2019.area.for.SR, paired = TRUE)
wilcox.test(x = mean.2019.phl.for.SR, y = mean.2019.lam.for.SR, paired = TRUE)


#traits against each other
wilcox.test(x = mean.2010.phl.for.SR, y = mean.2019.phl.for.SR, paired = TRUE)
wilcox.test(x = mean.2010.area.for.SR, y = mean.2019.area.for.SR, paired = TRUE)
wilcox.test(x = mean.2010.lam.for.SR, y = mean.2019.lam.for.SR, paired = TRUE)
#end 5b)


##### 3b) Resampling Power analyses for supporting information. #####
##### 3b-1) Proprtion of observations with sagrei present, Control islands only, generalized linear mixed effects model #####
#Did the proportion of A. carolinensis observations with an A. sagrei nearby change from 2010 to 2019?
#COMBINE ph.Asonly.2010 AND 2019, REMOVE S TWIN, RUN GLM
ph.Asonly.2010.t<-subset(ph.Asonly.2010, select = c("island","sex","perch.lizard","as.binary"))
ph.Asonly.2010.t$year<-rep(2010, length(ph.Asonly.2010.t$island))

ph.Asonly.2019.t<-subset(ph.Asonly.2019, ph.Asonly.2019$island!="south.twin",select = c("island","sex","ph.liz.cm","as.binary"))
colnames(ph.Asonly.2019.t)=c("island","sex","perch.lizard","as.binary")
ph.Asonly.2019.t$year<-rep(2019, length(ph.Asonly.2019.t$island))
for(i in 1:length(ph.Asonly.2019.t$island)){
  if (ph.Asonly.2019.t[i,which(colnames(ph.Asonly.2019.t)=="island")]=="north.twin") {
    ph.Asonly.2019.t[i,which(colnames(ph.Asonly.2019.t)=="island")]="North.Twin"}
  if (ph.Asonly.2019.t[i,which(colnames(ph.Asonly.2019.t)=="island")]=="yang") {
    ph.Asonly.2019.t[i,which(colnames(ph.Asonly.2019.t)=="island")]="Yang"}
  if (ph.Asonly.2019.t[i,which(colnames(ph.Asonly.2019.t)=="island")]=="channel") {
    ph.Asonly.2019.t[i,which(colnames(ph.Asonly.2019.t)=="island")]="Channel"}
}
ph.Asonly.pooled <- as.data.frame(rbind(ph.Asonly.2010.t,ph.Asonly.2019.t))

#Model
ph.Asonly.pooled.lme <- ph.Asonly.pooled
summary(glmer(as.binary ~ as.character(year) + (1|island), family = binomial, data = ph.Asonly.pooled.lme,
              na.action = "na.omit"))
tapply(X = ph.Asonly.pooled$as.binary, INDEX = ph.Asonly.pooled$year, FUN = mean)
# end NH1)

#### 3b-2) - analyzing 2010 data with reduced sample size (to match our 2019 sampling effort)
###get all possible 3 island combos from 2010 (one and two species islands separately)###

#### 3b-2) TOEPADS #####
tp.2010.NH <- tp.2010.tt[tp.2010.tt$sex == "m" | tp.2010.tt$sex == "f", ]
tp.2010.NH$category <- ifelse(tp.2010.NH$island %in% c("hornet", "osprey", "crescent", "s.twin", "pine"), 0, 1)

#one species islands
onespec.2010 <- c("hornet", "osprey", "crescent", "s.twin", "pine")
onespec.matrix <- matrix(nrow=5000,ncol=3)
for (i in 1:nrow(onespec.matrix)){
  onespec.matrix[i,1:3]=sample(onespec.2010,3,replace = FALSE)
}
#sort all rows into alphabetical order and paste as columns in new matrix
onespec.combos <- apply(onespec.matrix[, 1:3], 1, function(x) paste0(sort(x))) 
onespec.combos.final <- unique.array(onespec.combos,MARGIN = 2) #MARGIN=1 for rows and 2 for columns; will identify only unique columns


#two species islands
twospec.2010=c("yang","n.twin","yin","lizard","hook","channel")
twospec.matrix <- matrix(nrow=5000,ncol=3)
for (i in 1:nrow(twospec.matrix)){
  twospec.matrix[i,1:3]=sample(twospec.2010,3,replace = FALSE)
}
#sort all rows into alphabetical order and paste as columns in new matrix
twospec.combos <- apply(twospec.matrix[, 1:3], 1, function(x) paste0(sort(x))) 
twospec.combos.final <-unique.array(twospec.combos,MARGIN = 2) #MARGIN=1 for rows and 2 for columns; will identify only unique columns


##THE SUPER LOOP:
#1) create dataframe for storing results
#2) choose a combo of 3-island sets (200 total combos)
#3) subset 2010 data based on chosen combo
#4) separate by sex and size correct
#5) pool sexes back together into single data drame
#6) run linear model to compare toepad traits - store effect size and p-value

into.the.toeverse<-matrix(nrow=200,ncol=10)
into.the.toeverse<-as.data.frame(into.the.toeverse)
colnames(into.the.toeverse)=c("isl1","isl2","isl3","isl4","isl5","isl6","area.ES","area.p","lam.ES","lam.p")
counter=as.numeric()

for(i in 1:ncol(onespec.combos.final)){
  for(j in 1:ncol(twospec.combos.final)){
    
    #Subset
    tp.2010.temp <- tp.2010.NH[tp.2010.NH$island == as.character(onespec.combos.final[1,i]) | 
                                 tp.2010.NH$island == as.character(onespec.combos.final[2,i]) |
                                 tp.2010.NH$island == as.character(onespec.combos.final[3,i]) |
                                 tp.2010.NH$island == as.character(twospec.combos.final[1,j]) |
                                 tp.2010.NH$island == as.character(twospec.combos.final[2,j]) |
                                 tp.2010.NH$island == as.character(twospec.combos.final[3,j]) ,]
    
    # Size Correction - AREA ONLY 
    tp.2010.f.temp <- tp.2010.temp [tp.2010.temp$sex=="f" , ]
    tp.2010.f.sc.temp <- f.sizecorrect(standard.length = tp.2010.f.temp$svl, vector.of.columns = 7, 
                                       dataframe.to.use = tp.2010.f.temp, watershed = tp.2010.f.temp$island)
    
    tp.2010.m.temp <- tp.2010.temp [tp.2010.temp$sex=="m" , ]
    tp.2010.m.sc.temp <- f.sizecorrect(standard.length = tp.2010.m.temp$svl, vector.of.columns = 7, 
                                       dataframe.to.use = tp.2010.m.temp, watershed = tp.2010.m.temp$island)
    
    tp.2010.pooled.sc.temp <- merge(x = tp.2010.f.sc.temp, y = tp.2010.m.sc.temp, all = TRUE)
    
    # Models 
    area.lme.temp <- summary(lme(area.mm2.av.sc ~ as.character(category) + sex, 
                                 data = tp.2010.pooled.sc.temp, random = ~1|island, na.action = "na.omit"))
    lam.lme.temp <-summary(lme(count.av ~ as.character(category) + sex,
                               data = tp.2010.pooled.sc.temp, random = ~1|island, na.action = "na.omit"))
    
    # Store everything 
    counter <- c(counter,1)
    userow <-as.numeric(length(counter))
    into.the.toeverse[userow,1]<-as.character(onespec.combos.final[1,i])
    into.the.toeverse[userow,2]<-as.character(onespec.combos.final[2,i])
    into.the.toeverse[userow,3]<-as.character(onespec.combos.final[3,i])
    into.the.toeverse[userow,4]<-as.character(twospec.combos.final[1,j])
    into.the.toeverse[userow,5]<-as.character(twospec.combos.final[2,j])
    into.the.toeverse[userow,6]<-as.character(twospec.combos.final[3,j])
    into.the.toeverse[userow,7]<-area.lme.temp$tTable[2,1] #effect size
    into.the.toeverse[userow,8]<-area.lme.temp$tTable[2,5] #p-value
    into.the.toeverse[userow,9]<-lam.lme.temp$tTable[2,1] #effect size
    into.the.toeverse[userow,10]<-lam.lme.temp$tTable[2,5] #p-value
    
  }
}

# Toepad change for the 6 control islands in 2019
morph.2019.NH <- morph.2019[-(which(morph.2019$island == "south.twin")), ]
morph.2019.NH$category <- ifelse(morph.2019.NH$island %in% c("hornet", "osprey", "crescent"), 0, 1)

morph.2019.NH.temp <- morph.2019.NH[morph.2019.NH$sex=="f" , ]
morph.2019.NH.f.sc <- f.sizecorrect(standard.length = morph.2019.NH.temp$svl, vector.of.columns = 7, 
                                    dataframe.to.use = morph.2019.NH.temp, watershed = morph.2019.NH.temp$island)

morph.2019.NH.temp <- morph.2019.NH[morph.2019.NH$sex=="m" , ]
morph.2019.NH.m.sc <- f.sizecorrect(standard.length = morph.2019.NH.temp$svl, vector.of.columns = 7, 
                                    dataframe.to.use = morph.2019.NH.temp, watershed = morph.2019.NH.temp$island)

morph.2019.NH.sc <- merge(x = morph.2019.NH.f.sc, y = morph.2019.NH.m.sc, all = TRUE)

###Models
tpAREA.2019.params<-summary(lme(toe.area.mm2.av.sc ~ as.character(category) + sex, 
                                data = morph.2019.NH.sc, 
                                random = ~1|island, na.action = "na.omit"))
tpAREA.2019.params
tpAREA.2019.params$tTable[2,1] #effect size for one- vs two-species islands
tpAREA.2019.params$tTable[2,5] #p-value

tapply(X = morph.2019.NH.sc$toe.area.mm2.av.sc, INDEX = morph.2019.NH.sc$category, FUN = mean)


tpLAM.2019.params<-summary(lme(lam.num.av ~ as.character(category) + sex, 
                               data = morph.2019.NH.sc, 
                               random = ~1|island, na.action = "na.omit"))
tpLAM.2019.params
tpLAM.2019.params$tTable[2,1] #effect size size for one- vs two-species islands
tpLAM.2019.params$tTable[2,5] #p-value

tapply(X = morph.2019.NH.sc$lam.num.av, INDEX = morph.2019.NH.sc$category, FUN = mean)

#How do model parameters for the 6 control islands in 2019 compare to analyses from 2010 subsetted data?

#this next line is needed to identify the row in into.the.toeverse that contains the exact 6 control islands from 2019
findem1=as.numeric(which(apply(into.the.toeverse, 1, function(r) {any(r %in% c("yang")) & any(r %in% c("n.twin")) & 
    any(r %in% c("channel")) & any(r %in% c("hornet")) & any(r %in% c("osprey")) & any(r %in% c("crescent"))})))

plot(into.the.toeverse$area.ES,into.the.toeverse$area.p,
     main="2010 Toepad Area Comparisons - Island Subsets",
     xlab="Effect Size",ylab="p-value",col="gray76")
abline(h=0.05,col="red")
points(x=tpAREA.2019.params$tTable[2,1],y=tpAREA.2019.params$tTable[2,5],col="blue",pch=15) #2019 control islands
points(x=into.the.toeverse$area.ES[findem1],y=into.the.toeverse$area.p[findem1],col="blue",pch=17) #exact 2019 control islands but in 2010

sum(into.the.toeverse$area.p<.05) #26 combos from 2010 had significant toepad area differences


plot(into.the.toeverse$lam.ES,into.the.toeverse$lam.p,
     main="2010 Lamellae Comparisons - Island Subsets",
     xlab="Effect Size",ylab="p-value",col="gray76")
abline(h=0.05,col="red")
points(x=tpLAM.2019.params$tTable[2,1],y=tpLAM.2019.params$tTable[2,5],col="blue",pch=15) #2019 control islands
points(x=into.the.toeverse$lam.ES[findem1],y=into.the.toeverse$lam.p[findem1],col="blue",pch=17) #exact 2019 control islands but in 2010

sum(into.the.toeverse$lam.p<.05) #45 combos from 2010 had significant toepad area differences
# end 3b-2)

#### 3b-3) PERCH HEIGHT #####

ph.2010.NH <- ph.2010.final.tt[ph.2010.final.tt$sex == "m" | ph.2010.final.tt$sex == "f" , ]
ph.2010.NH$category <- ifelse(ph.2010.NH$island %in% c("Hornet", "Osprey", "Crescent", "S.twin", "Pine"), 0, 1)


#one species islands - NEEDS TO BE RERUN BECAUSE OF DIFFERENCE IN CAPITALIZATION
unique(ph.2010.NH$island)
onespec.2010 <- c("Hornet", "Osprey", "Crescent", "South.Twin", "Pine")
onespec.matrix <- matrix(nrow=5000,ncol=3)
for (i in 1:nrow(onespec.matrix)){
  onespec.matrix[i,1:3]=sample(onespec.2010,3,replace = FALSE)
}
#sort all rows into alphabetical order and paste as columns in new matrix
onespec.combos <- apply(onespec.matrix[, 1:3], 1, function(x) paste0(sort(x))) 
onespec.combos.final <- unique.array(onespec.combos,MARGIN = 2) #MARGIN=1 for rows and 2 for columns; will identify only unique columns


#two species islands
twospec.2010=c("Yang","North.Twin","Yin","Lizard","Hook","Channel")
twospec.matrix <- matrix(nrow=5000,ncol=3)
for (i in 1:nrow(twospec.matrix)){
  twospec.matrix[i,1:3]=sample(twospec.2010,3,replace = FALSE)
}
#sort all rows into alphabetical order and paste as columns in new matrix
twospec.combos <- apply(twospec.matrix[, 1:3], 1, function(x) paste0(sort(x))) 
twospec.combos.final <-unique.array(twospec.combos,MARGIN = 2) #MARGIN=1 for rows and 2 for columns; will identify only unique columns


##THE SUPER LOOP:
#1) create dataframe for storing results
#2) choose a combo of 3-island sets (200 total combos)
#3) subset 2010 data based on chosen combo
#4) run linear model to compare perch heights - store effect size and p-value

into.the.perchverse<-matrix(nrow=200,ncol=8)
into.the.perchverse<-as.data.frame(into.the.perchverse)
colnames(into.the.perchverse)=c("isl1","isl2","isl3","isl4","isl5","isl6","ph.ES","ph.p")
counter=as.numeric()

for(i in 1:ncol(onespec.combos.final)){
  for(j in 1:ncol(twospec.combos.final)){
    
    ###Subset
    ph.2010.temp <- ph.2010.NH[ph.2010.NH$island == as.character(onespec.combos.final[1,i]) | 
                                 ph.2010.NH$island == as.character(onespec.combos.final[2,i]) |
                                 ph.2010.NH$island == as.character(onespec.combos.final[3,i]) |
                                 ph.2010.NH$island == as.character(twospec.combos.final[1,j]) |
                                 ph.2010.NH$island == as.character(twospec.combos.final[2,j]) |
                                 ph.2010.NH$island == as.character(twospec.combos.final[3,j]) ,]
    
    
    ##### Model 
    ph.lme.temp <- summary(lme(perch.lizard ~ as.character(category) + sex, 
                               data = ph.2010.temp, random = ~1|island, na.action = "na.omit"))
    
    
    ##### Store everything 
    counter <- c(counter,1)
    userow <-as.numeric(length(counter))
    into.the.perchverse[userow,1]<-as.character(onespec.combos.final[1,i])
    into.the.perchverse[userow,2]<-as.character(onespec.combos.final[2,i])
    into.the.perchverse[userow,3]<-as.character(onespec.combos.final[3,i])
    into.the.perchverse[userow,4]<-as.character(twospec.combos.final[1,j])
    into.the.perchverse[userow,5]<-as.character(twospec.combos.final[2,j])
    into.the.perchverse[userow,6]<-as.character(twospec.combos.final[3,j])
    into.the.perchverse[userow,7]<-ph.lme.temp$tTable[2,1] #effect size
    into.the.perchverse[userow,8]<-ph.lme.temp$tTable[2,5] #p-value
    
  }
}



# Perch height change for the 6 control islands in 2019 #

ph.final.NH <- ph.final[-(which(ph.final$island == "south.twin")), ]
ph.final.NH$category <- ifelse(ph.final.NH$island %in% c("hornet", "osprey", "crescent"), 0, 1)

###Model
ph.2019.params<-summary(lme(ph.liz.cm ~ as.character(category) + sex, 
                            data = ph.final.NH, 
                            random = ~1|island, na.action = "na.omit"))
ph.2019.params$tTable[2,1] #effect size for one- vs two-species islands
ph.2019.params$tTable[2,5] #p-value

tapply(X = ph.final.NH$ph.liz.cm, INDEX = ph.final.NH$category, FUN = mean)



#How do model parameters for the 6 control islands in 2019 compare to analyses from 2010 subsetted data?

#this next line is needed to identify the row in into.the.toeverse that contains the exact 6 control islands from 2019
findem2=as.numeric(which(apply(into.the.perchverse, 1, function(r) {any(r %in% c("Yang")) & any(r %in% c("North.Twin")) & 
    any(r %in% c("Channel")) & any(r %in% c("Hornet")) & any(r %in% c("Osprey")) & any(r %in% c("Crescent"))})))

plot(into.the.perchverse$ph.ES,into.the.perchverse$ph.p,
     main="2010 Perch Height Comparisons - Island Subsets",
     xlab="Effect Size",ylab="p-value",col="gray76",
     xlim = c(20,130), ylim=c(0,0.3))
abline(h=0.05,col="red")
points(x=ph.2019.params$tTable[2,1],y=ph.2019.params$tTable[2,5],col="blue",pch=15) #2019 control islands
points(x=into.the.perchverse$ph.ES[findem2],y=into.the.perchverse$ph.p[findem2],col="blue",pch=17) #exact 2019 control islands but in 2010

sum(into.the.perchverse$ph.p<.05) #97 combos from 2010 had significant pergh height differences

#end 3b-3)
#end 3b)
#end 3)

##### BESPOKE FUNCTIONS #####
#### 1 f.size.correct ####
# OKE ET AL COMMON GARDEN MANUSCRIPT METHOD
#"All relative warps and univariate shape traits were allometrically standardized to a common body size (Reist 1985; Lleonart et al. 2000) based on Ms = M0 (Ls / L0)^b, where Ms is the standardized trait value (mm), M0 is the non-standardized trait value (mm), Ls is the overall mean centroid size (for RWs) or standard length (mm, for univariate shape traits), and L0 is the centroid size or standard length (mm) of the individual. The common within-group slope, b, was calculated from a linear mixed model of log10(M0) regressed on log10(L0), with group included as a random factor (Reist 1985; Lleonart et al. 2000)." - Oke et al.    The size correction formula described above: Ms = M0(Ls/L0)^b

#variables:
#standard.length = the trait to be used for size correction. snout vent length
#watershed = the random factor in the linear mixed model. island
#dataframe.to.use = the data
# vector.of.columns = the columns in dataframe.to.use that are to be size corrected
#The output will have all the original columns from dataframe.to.use with the size-corrected columns added on at the end.

#include year but not sex in the mixed model. watershed is island
f.sizecorrect.withoutsex.withyear <- function(standard.length, vector.of.columns, dataframe.to.use, watershed, year) {
  #Calculate overall mean standard length (Ls)
  Ls <- mean(standard.length, na.rm = TRUE) 
  #Call individual standard length
  L0 <- standard.length
  #Calculate common, within-group slope, b, for each trait. 
  b.vector.lmm <- vector()
  for (i in vector.of.columns) {
    abcd <- (dataframe.to.use[i])
    #b.model <- lmer((abcd[,])~((standard.length))*sex + (1|watershed)) #originally didn't have log10 in this calculation of b. Added 16 Nov. 2019
    b.model <- lmer(log10(abcd[,])~(log10(standard.length)) + year + (1|watershed))
    b <- coef(summary(b.model))[2,1]
    b.vector.lmm <- c(b.vector.lmm, b)
  }
  # size correct
  xx <- dataframe.to.use  
  columnnames <- colnames(xx)
  j=1
  for (i in vector.of.columns) {
    M0 <- xx[,i] #grab the appropriate column of data
    Ms = M0 * ((Ls/L0)^b.vector.lmm[j]) #size correction formula
    j=j+1
    columnnames <- c(columnnames, paste(colnames(xx[i]), "sc", sep = "."))
    xx <- cbind(xx, Ms)
  }
  colnames(xx) <- columnnames # Rename thh columns in the temporary dataframe xx
  return(xx) #Output a new dataframe with the name provided in "outputfilename"
  return(print(b.vector.lmm))
}

##### FIGURES #####
##### Figure 1. Perch heights by island, 2010 to 2019 #####
#parameter table for 2010
#uses ph.2010.final which was created in analysis.southTwin.forPaper.v2.R, 0a)
ph.means.2010 <- round(tapply(X = ph.2010.final$perch.lizard, INDEX = ph.2010.final$island, FUN = mean, na.rm = TRUE), 1)
ph.sds.2010 <- round(tapply(X = ph.2010.final$perch.lizard, INDEX = ph.2010.final$island, FUN = sd, na.rm = TRUE), 1)
ph.n.2010 <- tapply(X = ph.2010.final$perch.lizard, INDEX = ph.2010.final$island, FUN = length)
ph.ses.2010 <- round(ph.sds.2010/sqrt(ph.n.2010), 1)
ph.table.for.plot.2010 <- data.frame(ph.means.2010, ph.sds.2010, ph.n.2010, ph.ses.2010)
ph.table.for.plot.sort.2010 <- ph.table.for.plot.2010[order(ph.table.for.plot.2010$ph.means), ]

#parameter table for 2019
#uses ph.final which was created in analysis.southTwin.forPaper.v2.R, 0a)
ph.for.plot <- ph.final
ph.means <- round(tapply(X = ph.for.plot$ph.liz.cm, INDEX = ph.for.plot$isle.code, FUN = mean, na.rm = TRUE), 1)
ph.sds <- round(tapply(X = ph.for.plot$ph.liz.cm, INDEX = ph.for.plot$isle.code, FUN = sd, na.rm = TRUE), 1)
ph.n <- tapply(X = ph.for.plot$ph.liz.cm, INDEX = ph.for.plot$isle.code, FUN = length)
ph.ses <- round(ph.sds/sqrt(ph.n), 1)
ph.table.for.plot <- data.frame(ph.means, ph.sds, ph.n, ph.ses)
ph.table.for.plot.sort <- ph.table.for.plot[order(ph.table.for.plot$ph.means), ]

#Figure
par(mfrow = c(1,1), mar = c(0, 3.5, 1, 1), oma = c(0,0,0,0))
plot(ph.table.for.plot.sort.2010$ph.means.2010 ~ c(seq(0.94, 1.06, by = 0.02)), xlim = c(0.9, 1.6), ylim = c(75, 275), pch = c(1, 1, 15, 1, 16, 16, 16), col = c(rep("grey50", 2), "black", rep("grey50", 4)), axes = FALSE, ylab = "", cex = 1.2)
segments(x0 = seq(0.94, 1.06, by = 0.02), x1 = seq(0.94, 1.06, by = 0.02), y0 = ph.table.for.plot.sort.2010$ph.means.2010, y1 = ph.table.for.plot.sort.2010$ph.means.2010 + ph.table.for.plot.sort.2010$ph.ses.2010, col = c(rep("grey50", 2), "black", rep("grey50", 4)))
segments(x0 = seq(0.94, 1.06, by = 0.02), x1 = seq(0.94, 1.06, by = 0.02), y0 = ph.table.for.plot.sort.2010$ph.means.2010, y1 = ph.table.for.plot.sort.2010$ph.means.2010 - ph.table.for.plot.sort.2010$ph.ses.2010, col = c(rep("grey50", 2), "black", rep("grey50", 4)))
axis(side = 2, at = c(75, 125, 175, 225, 275), las = 1)

points(x = seq(1.44, 1.56, by = 0.02), y = ph.table.for.plot.sort$ph.means, col = c(rep("grey50", 3), "black", rep("grey50", 3)), pch = c(1, 1, 1, 15, 16, 16, 16), cex = 1.2)
segments(x0 = seq(1.44, 1.56, by = 0.02), x1 = seq(1.44, 1.56, by = 0.02), y0 = ph.table.for.plot.sort$ph.means, y1 = ph.table.for.plot.sort$ph.means + ph.table.for.plot.sort$ph.ses, col = c(rep("grey50", 3), "black", rep("grey50", 3)))
segments(x0 = seq(1.44, 1.56, by = 0.02), x1 = seq(1.44, 1.56, by = 0.02), y0 = ph.table.for.plot.sort$ph.means, y1 = ph.table.for.plot.sort$ph.means - ph.table.for.plot.sort$ph.ses, col = c(rep("grey50", 3), "black", rep("grey50", 3)))
segments(x0 = seq(0.94, 1.06, by = 0.02), x1 = seq(1.44, 1.56, by = 0.02)[c(3, 1, 4, 2, 7, 5, 6)], y0 = ph.table.for.plot.sort.2010$ph.means.2010, y1 = ph.table.for.plot.sort$ph.means[c(3, 1, 4, 2, 7, 5, 6)], lty = (c(1,1,1,1,1,1,1)), col = c(rep("grey50", 2), "black", rep("grey50", 4)))
text(x = 1, y = 75, labels = "2010"); text(x = 1.5, y = 75, labels = "2019")
mtext(text = "Perch Height (cm)", side = 2, line = 2.5, las = 0)

text(x = seq(0.94, 1.06, by = 0.02)-0.03, y = ph.table.for.plot.sort.2010$ph.means.2010, labels = c("hor", "osp", "stw", "cre", "yan", "ntw", "cha"), cex = 1)
# end Figure 1


##### Figure 3. Effects of sagrei density on Perch Height #####
#uses "density.vs.ph" generated in analysis.southTwin.forPaper.v2.R, 1d)
par(mfrow = c(1,1), mar = c(2, 4.75, 1, 1), las = 0, oma = c(0,0,0,0))
density.vs.ph.with.error <- density.vs.ph
density.vs.ph.with.error$se <- ph.table.for.plot.sort$ph.ses[c(7, 6, 4, 5, 1, 2, 3)]

plot(density.vs.ph.with.error$ph.means.2019[density.vs.ph.with.error$mean.density>0] ~ density.vs.ph.with.error$mean.density[density.vs.ph.with.error$mean.density>0], pch = 16, ylim = c(75, 275), xlim = c(0, 45), xlab = "", ylab = "", axes = FALSE, cex = 1.2)
jittered.0 <- c(0, 0.3, 0.6)
points(x = jittered.0, y = density.vs.ph.with.error$ph.means.2019[which(density.vs.ph.with.error$mean.density==0)], pch = 16, cex = 1.2)
ex.0 <- c(density.vs.ph.with.error$mean.density[density.vs.ph.with.error$mean.density>0], jittered.0)
ex.1 <-c(density.vs.ph.with.error$mean.density[density.vs.ph.with.error$mean.density>0], jittered.0)
wy.0 <- density.vs.ph.with.error$ph.means.2019
wy.1 <- density.vs.ph.with.error$ph.means.2019 + density.vs.ph.with.error$se
segments(x0 = ex.0, x1 = ex.1, y0 = wy.0, y1 = wy.1)

segments(x0 = c(density.vs.ph.with.error$mean.density[density.vs.ph.with.error$mean.density>0], jittered.0), x1 = c(density.vs.ph.with.error$mean.density[density.vs.ph.with.error$mean.density>0], jittered.0), y0 = density.vs.ph.with.error$ph.means.2019, y1 = density.vs.ph.with.error$ph.means.2019 - density.vs.ph.with.error$se)

segments(x0 = 0, x1 = 43.0, y0 = model.1d.1$coefficients[1], y1 = model.1d.1$coefficients[1] + 43.0 * model.1d.1$coefficients[2])  #from analysis.southTwin.forPaper.v2.R 1d)
axis(side = 2, at = c(75, 125, 175, 225, 275), las = 1)
mtext(text = expression(atop(italic(A.~carolinensis), "mean perch height (cm)")), side = 2, line = 2.2)
axis(side = 1, at = c(0.0, 18.5, 32.5, 41.5, 43.0), labels = FALSE)
text(x = c("0.0", 18.5, 32.5), y = 72, labels = c(0.0, 18.5, 32.5))
text(x = c(18.5, 32.5), y = 78, labels = c("stw", "ntw"))
text(x = c(41.5,41.5), y = c(72, 78), labels = c(41.5, "cha"), adj = 1)
text(x = c(43.0, 43.0), y = c(72, 78), labels = c("43.0", "yan"), adj = 0)
mtext(text = expression("Average # of "*italic(A.~sagrei)*" seen during 1 hour searches"), side = 1, line = 1)
#end Figure 3

##### Figure 4. Effects of potential sagrei interactions on perch height #####
#A = 2019, B = 2010
par(mfrow = c(2,1), mar = c(1,4,0,0), oma = c(0,0,0,0))

#A
#Uses ph.final which was created in analysis.southTwin.forPaper.v2.R, 0a) 
ph.final.invaded <- ph.final[-which(is.na(ph.final$as.int.sex)), ]
ph.final.invaded$as.int.pool <- ifelse(ph.final.invaded$as.int.sex == "no", 0, 1)
stws.2 <- which(ph.final.invaded$island == "south.twin")
ph.final.invaded[stws.2, 1] <- "a.south.twin"


plot(y = ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "a.south.twin"], x = jitter(x = rep(1, length(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "a.south.twin"]))), xlim = c(0.75, 4.25), xaxt = "n", type = "n", ylim = c(0, 500), ylab = "", axes = FALSE)
axis(2, las = 1)
#mtext(text = c("South Twin", "North Twin", "Channel", "Yang"), side = 1, line = -1.5, at = c(1, 2, 3, 4), cex = 1.2)
#mtext(text = "Perch height (cm)", side = 2, line = 2.5, cex = 1.2)
mtext(text = expression(""*italic(A.~carolinensis)*" perch height (cm)"), side = 2, line = 2.5, cex = 1.2)

#text(x = 3, y = 600, labels = "green = sagrei not nearby, brown = sagrei nearby")
#South Twin
points(y = ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "a.south.twin" & ph.final.invaded$as.int.pool == 0], x = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "a.south.twin" & ph.final.invaded$as.int.pool == 0] == "f", 0.88, 0.92), pch = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "a.south.twin" & ph.final.invaded$as.int.pool == 0] == "f", 1, 2), cex = 1.2)
points(y = ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "a.south.twin" & ph.final.invaded$as.int.pool == 1], x = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "a.south.twin" & ph.final.invaded$as.int.pool == 1] == "f", 1.08, 1.12), pch = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "a.south.twin" & ph.final.invaded$as.int.pool == 1] == "f", 16, 17), cex = 1.2)
segments(x0 = 0.88, x1 = 1.08, y0 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "a.south.twin" & ph.final.invaded$as.int.pool == 0 & ph.final.invaded$sex == "f"]), y1 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "a.south.twin" & ph.final.invaded$as.int.pool == 1 & ph.final.invaded$sex == "f"]), lwd = 4)
segments(x0 = 0.92, x1 = 1.12, y0 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "a.south.twin" & ph.final.invaded$as.int.pool == 0 & ph.final.invaded$sex == "m"]), y1 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "a.south.twin" & ph.final.invaded$as.int.pool == 1 & ph.final.invaded$sex == "m"]), lwd = 4, lty = 2)


#North Twin
points(y = ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "north.twin" & ph.final.invaded$as.int.pool == 0], x = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "north.twin" & ph.final.invaded$as.int.pool == 0] == "f", 1.88, 1.92), pch = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "north.twin" & ph.final.invaded$as.int.pool == 0] == "f", 1, 2), cex = 1.2) #females are circles, males are triangles
points(y = ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "north.twin" & ph.final.invaded$as.int.pool == 1], x = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "north.twin" & ph.final.invaded$as.int.pool == 1] == "f", 2.08, 2.12), pch = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "north.twin" & ph.final.invaded$as.int.pool == 1] == "f", 16, 17), cex = 1.2)
segments(x0 = 1.88, x1 = 2.08, y0 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "north.twin" & ph.final.invaded$as.int.pool == 0 & ph.final.invaded$sex == "f"]), y1 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "north.twin" & ph.final.invaded$as.int.pool == 1 & ph.final.invaded$sex == "f"]), lwd = 4)
segments(x0 = 1.92, x1 = 2.12, y0 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "north.twin" & ph.final.invaded$as.int.pool == 0 & ph.final.invaded$sex == "m"]), y1 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "north.twin" & ph.final.invaded$as.int.pool == 1 & ph.final.invaded$sex == "m"]), lwd = 4, lty = 2)

#Channel
points(y = ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "channel" & ph.final.invaded$as.int.pool == 0], x = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "channel" & ph.final.invaded$as.int.pool == 0] == "f", 2.88, 2.92), pch = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "channel" & ph.final.invaded$as.int.pool == 0] == "f", 1, 2), cex = 1.2)
points(y = ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "channel" & ph.final.invaded$as.int.pool == 1], x = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "channel" & ph.final.invaded$as.int.pool == 1] == "f", 3.08, 3.12), pch = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "channel" & ph.final.invaded$as.int.pool == 1] == "f", 16, 17), cex = 1.2)
segments(x0 = 2.88, x1 = 3.08, y0 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "channel" & ph.final.invaded$as.int.pool == 0 & ph.final.invaded$sex == "f"]), y1 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "channel" & ph.final.invaded$as.int.pool == 1 & ph.final.invaded$sex == "f"]), lwd = 4)
segments(x0 = 2.92, x1 = 3.12, y0 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "channel" & ph.final.invaded$as.int.pool == 0 & ph.final.invaded$sex == "m"]), y1 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "channel" & ph.final.invaded$as.int.pool == 1 & ph.final.invaded$sex == "m"]), lwd = 4, lty = 2)

#Yang
points(y = ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "yang" & ph.final.invaded$as.int.pool == 0], x = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "yang" & ph.final.invaded$as.int.pool == 0] == "f", 3.88, 3.92), pch = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "yang" & ph.final.invaded$as.int.pool == 0] == "f", 1, 2), cex = 1.2)
points(y = ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "yang" & ph.final.invaded$as.int.pool == 1], x = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "yang" & ph.final.invaded$as.int.pool == 1] == "f", 4.08, 4.12), pch = ifelse(ph.final.invaded$sex[ph.final.invaded$island == "yang" & ph.final.invaded$as.int.pool == 1] == "f", 16, 17), cex = 1.2)
segments(x0 = 3.88, x1 = 4.08, y0 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "yang" & ph.final.invaded$as.int.pool == 0 & ph.final.invaded$sex == "f"]), y1 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "yang" & ph.final.invaded$as.int.pool == 1 & ph.final.invaded$sex == "f"]), lwd = 4)
segments(x0 = 3.92, x1 = 4.12, y0 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "yang" & ph.final.invaded$as.int.pool == 0 & ph.final.invaded$sex == "m"]), y1 = mean(ph.final.invaded$ph.liz.cm[ph.final.invaded$island == "yang" & ph.final.invaded$as.int.pool == 1 & ph.final.invaded$sex == "m"]), lwd = 4, lty = 2)

text(labels= "A. 2019", x =  2.5, y = 500, cex = 1.2)


#legend
points(x = c(3.7, 3.7, 3.9, 3.9), y = c(5, -10, 5, -10), pch = c(1, 16, 2, 17))
text(x = c(3.7, 3.9), y = c(20, 20), labels = c("f", "m"))
#text(x = c(3.5, 3.5), y = c(5, -10), labels = c("sagrei not nearby", "sagrei nearby"), adj = 0.9)
text(x = c(3.5, 3.5), y = c(5, -10), labels = c(expression(""*italic(A.~sagrei)*" not nearby"), expression(""*italic(A.~sagrei)*" nearby")), adj = 0.9)
segments(x0 = 3.6, x1 = 3.6, y0 = -15, y1 = 12.5)
segments(x0 = 3.6, x1 = 4, y0 = 12.5, y1 = 12.5)

#B
plot(y = ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "North Twin"], x = jitter(x = rep(2, length(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "North Twin"]))), xlim = c(0.75, 4.25), xaxt = "n", type = "n", ylim = c(0,600), ylab = "", xlab = "", axes = FALSE)
axis(2, las = 1)
mtext(text = c("South Twin", "North Twin", "Channel", "Yang"), side = 1, line = 0, at = c(1, 2, 3, 4), cex = 1.2)
mtext(text = expression(""*italic(A.~carolinensis)*" perch height (cm)"), side = 2, line = 2.5, cex = 1.2)

#text(x = 3, y = 600, labels = "green = sagrei not nearby, brown = sagrei nearby")

#North Twin
points(y = ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "North.Twin" & ph.Asonly.2010$as.binary == 0], x = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "North.Twin" & ph.Asonly.2010$as.binary == 0] == "f", 1.88, 1.92), pch = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "North.Twin" & ph.Asonly.2010$as.binary == 0] == "f", 1, 2), cex = 1.2)
points(y = ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "North.Twin" & ph.Asonly.2010$as.binary == 1], x = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "North.Twin" & ph.Asonly.2010$as.binary == 1] == "f", 2.08, 2.12), pch = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "North.Twin" & ph.Asonly.2010$as.binary == 1] == "f", 16, 17), cex = 1.2)
segments(x0 = 1.88, x1 = 2.08, y0 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "North.Twin" & ph.Asonly.2010$as.binary == 0 & ph.Asonly.2010$sex == "f"], na.rm = TRUE), y1 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "North.Twin" & ph.Asonly.2010$as.binary == 1 & ph.Asonly.2010$sex == "f"], na.rm = TRUE), lwd = 4)
segments(x0 = 1.92, x1 = 2.12, y0 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "North.Twin" & ph.Asonly.2010$as.binary == 0 & ph.Asonly.2010$sex == "m"], na.rm = TRUE), y1 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "North.Twin" & ph.Asonly.2010$as.binary == 1 & ph.Asonly.2010$sex == "m"], na.rm = TRUE), lwd = 4, lty = 2)

#Channel
points(y = ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Channel" & ph.Asonly.2010$as.binary == 0], x = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "Channel" & ph.Asonly.2010$as.binary == 0] == "f", 2.88, 2.92), pch = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "Channel" & ph.Asonly.2010$as.binary == 0] == "f", 1, 2), cex = 1.2)
points(y = ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Channel" & ph.Asonly.2010$as.binary == 1], x = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "Channel" & ph.Asonly.2010$as.binary == 1] == "f", 3.08, 3.12), pch = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "Channel" & ph.Asonly.2010$as.binary == 1] == "f", 16, 17), cex = 1.2)

segments(x0 = 2.88, x1 = 3.08, y0 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Channel" & ph.Asonly.2010$as.binary == 0 & ph.Asonly.2010$sex == "f"], na.rm = TRUE), y1 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Channel" & ph.Asonly.2010$as.binary == 1 & ph.Asonly.2010$sex == "f"], na.rm = TRUE), lwd = 4)
segments(x0 = 2.92, x1 = 3.12, y0 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Channel" & ph.Asonly.2010$as.binary == 0 & ph.Asonly.2010$sex == "m"], na.rm = TRUE), y1 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Channel" & ph.Asonly.2010$as.binary == 1 & ph.Asonly.2010$sex == "m"], na.rm = TRUE), lwd = 4, lty = 2)


#Yang
points(y = ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Yang" & ph.Asonly.2010$as.binary == 0], x = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "Yang" & ph.Asonly.2010$as.binary == 0] == "f", 3.88, 3.92), pch = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "Yang" & ph.Asonly.2010$as.binary == 0] == "f", 1, 2), cex = 1.2)
points(y = ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Yang" & ph.Asonly.2010$as.binary == 1], x = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "Yang" & ph.Asonly.2010$as.binary == 1] == "f", 4.08, 4.12), pch = ifelse(ph.Asonly.2010$sex[ph.Asonly.2010$island == "Yang" & ph.Asonly.2010$as.binary == 1] == "f", 16, 17), cex = 1.2)

segments(x0 = 3.88, x1 = 4.08, y0 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Yang" & ph.Asonly.2010$as.binary == 0 & ph.Asonly.2010$sex == "f"], na.rm = TRUE), y1 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Yang" & ph.Asonly.2010$as.binary == 1 & ph.Asonly.2010$sex == "f"], na.rm = TRUE), lwd = 4)
segments(x0 = 3.92, x1 = 4.12, y0 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Yang" & ph.Asonly.2010$as.binary == 0 & ph.Asonly.2010$sex == "m"], na.rm = TRUE), y1 = mean(ph.Asonly.2010$perch.lizard[ph.Asonly.2010$island == "Yang" & ph.Asonly.2010$as.binary == 1 & ph.Asonly.2010$sex == "m"], na.rm = TRUE), lwd = 4, lty = 2)

#legend continued
segments(x0 = 3.6, x1 = 4, y0 = 610, y1 = 610, lty = 1, lwd = 4); text(x = 3.5, y = 610, "female means", adj = .9, oma = TRUE)
segments(x0 = 3.6, x1 = 4, y0 = 590, y1 = 590, lty = 2, lwd = 4); text(x = 3.5, y = 590, "male means", adj = .9)

text(labels= "B. 2010", x =  2.5, y = 600, cex = 1.2)
#end Figure 4

##### Figure 2 Size corrected toepad area and lamella number change 2010 to 2019 #####
#make means and se tables
#morph.2010.sc and morph.2019.sc come from analysis.southTwin.forPaper.v2.R, 0c)
##2010
#area
tpa.means.2010 <- round(tapply(X = morph.2010.sc$area.sc, INDEX = morph.2010.sc$island, FUN = mean, na.rm = TRUE), 2)
tpa.sds.2010 <- round(tapply(X = morph.2010.sc$area.sc, INDEX = morph.2010.sc$island, FUN = sd, na.rm = TRUE), 2)
tpa.n.2010 <- tapply(X = morph.2010.sc$area.sc, INDEX = morph.2010.sc$island, FUN = length)
tpa.ses.2010 <- round(tpa.sds.2010/sqrt(tpa.n.2010), 3)
tpa.table.for.plot.2010 <- data.frame(tpa.means.2010, tpa.sds.2010, tpa.n.2010, tpa.ses.2010)
tpa.table.for.plot.sort.2010 <- tpa.table.for.plot.2010[order(tpa.table.for.plot.2010$tpa.means.2010), ]
#lamella number
tpl.means.2010 <- round(tapply(X = morph.2010.sc$lam, INDEX = morph.2010.sc$island, FUN = mean, na.rm = TRUE), 2)
tpl.sds.2010 <- round(tapply(X = morph.2010.sc$lam, INDEX = morph.2010.sc$island, FUN = sd, na.rm = TRUE), 2)
tpl.n.2010 <- tapply(X = morph.2010.sc$lam, INDEX = morph.2010.sc$island, FUN = length)
tpl.ses.2010 <- round(tpl.sds.2010/sqrt(tpl.n.2010), 3)
tpl.table.for.plot.2010 <- data.frame(tpl.means.2010, tpl.sds.2010, tpl.n.2010, tpl.ses.2010)
tpl.table.for.plot.sort.2010 <- tpl.table.for.plot.2010[order(tpl.table.for.plot.2010$tpl.means.2010), ]

##2019
#area
tpa.means.2019 <- round(tapply(X = morph.2019.sc$area.sc, INDEX = morph.2019.sc$island, FUN = mean, na.rm = TRUE), 2)
tpa.sds.2019 <- round(tapply(X = morph.2019.sc$area.sc, INDEX = morph.2019.sc$island, FUN = sd, na.rm = TRUE), 2)
tpa.n.2019 <- tapply(X = morph.2019.sc$area.sc, INDEX = morph.2019.sc$island, FUN = length)
tpa.ses.2019 <- round(tpa.sds.2019/sqrt(tpa.n.2019), 3)
tpa.table.for.plot.2019 <- data.frame(tpa.means.2019, tpa.sds.2019, tpa.n.2019, tpa.ses.2019)
tpa.table.for.plot.sort.2019 <- tpa.table.for.plot.2019[order(tpa.table.for.plot.2019$tpa.means.2019), ]
#lamella
tpl.means.2019 <- round(tapply(X = morph.2019.sc$lam, INDEX = morph.2019.sc$island, FUN = mean, na.rm = TRUE), 2)
tpl.sds.2019 <- round(tapply(X = morph.2019.sc$lam, INDEX = morph.2019.sc$island, FUN = sd, na.rm = TRUE), 2)
tpl.n.2019 <- tapply(X = morph.2019.sc$lam, INDEX = morph.2019.sc$island, FUN = length)
tpl.ses.2019 <- round(tpl.sds.2019/sqrt(tpl.n.2019), 3)
tpl.table.for.plot.2019 <- data.frame(tpl.means.2019, tpl.sds.2019, tpl.n.2019, tpl.ses.2019)
tpl.table.for.plot.sort.2019 <- tpl.table.for.plot.2019[order(tpl.table.for.plot.2019$tpl.means.2019), ]

#figure
par(mar = c(0, 4, 0, 0), mfrow = c(2,1), oma = c(0,0,0,0))
#toepad area
#oddly, rank order of size-corrected areas here is different from rank order of the size-corrected areas for 2014 plot, mostly because crescent jumps from 2 to 6. Regardless, analysis was done on non-size corrected so 2014 paper is fine. It's just that plot for Science Fig 3 2014 included a naive size-correction that wasn't the lme approach used here.
plot(tpa.table.for.plot.sort.2010$tpa.means.2010 ~ c(seq(0.94, 1.06, by = 0.02)), xlim = c(0.9, 1.6), ylim = c(3.0, 4.2), pch = c(15, 1, 16, 1, 1, 16, 16), col = c("black", rep("grey50", 6)), axes = FALSE, ylab = "", cex = 1.2)
segments(x0 = seq(0.94, 1.06, by = 0.02), x1 = seq(0.94, 1.06, by = 0.02), y0 = tpa.table.for.plot.sort.2010$tpa.means.2010, y1 = tpa.table.for.plot.sort.2010$tpa.means.2010 + tpa.table.for.plot.sort.2010$tpa.ses.2010, col = c("black", rep("grey50", 6)))
segments(x0 = seq(0.94, 1.06, by = 0.02), x1 = seq(0.94, 1.06, by = 0.02), y0 = tpa.table.for.plot.sort.2010$tpa.means.2010, y1 = tpa.table.for.plot.sort.2010$tpa.means.2010 - tpa.table.for.plot.sort.2010$tpa.ses.2010, col = c("black", rep("grey50", 6)))
axis(side = 2, at = c(3.0, 3.4, 3.8, 4.2), las = 1)
mtext(text = "Mean size-corrected Toepad Area (mm2)", side = 2, line = 2.5, at = 3.6)
text(x = 1, y = 3.0, labels = "2010"); text(x = 1.5, y = 3.0, labels = "2019", cex = 1)

points(x = seq(1.44, 1.56, by = 0.02), y = tpa.table.for.plot.sort.2019$tpa.means.2019, col = c(rep("grey50", 4), "black", rep("grey50", 2)), pch = c(1, 1, 1, 16, 15, 16, 16), cex = 1.2)
segments(x0 = seq(1.44, 1.56, by = 0.02), x1 = seq(1.44, 1.56, by = 0.02), y0 = tpa.table.for.plot.sort.2019$tpa.means.2019, y1 = tpa.table.for.plot.sort.2019$tpa.means.2019 + tpa.table.for.plot.sort.2019$tpa.ses.2019, col = c(rep("grey50", 4), "black", rep("grey50", 2)))
segments(x0 = seq(1.44, 1.56, by = 0.02), x1 = seq(1.44, 1.56, by = 0.02), y0 = tpa.table.for.plot.sort.2019$tpa.means.2019, y1 = tpa.table.for.plot.sort.2019$tpa.means.2019 - tpa.table.for.plot.sort.2019$tpa.ses.2019, col = c(rep("grey50", 4), "black", rep("grey50", 2)))
segments(x0 = seq(0.94, 1.06, by = 0.02), x1 = seq(1.44, 1.56, by = 0.02)[c(5, 1, 4, 2, 3, 7, 6)], y0 = tpa.table.for.plot.sort.2010$tpa.means.2010, y1 = tpa.table.for.plot.sort.2019$tpa.means.2019[c(5, 1, 4, 2, 3, 7, 6)], lty = rep(1, 7), lwd = c(1.0, rep(0.5, 6)))
#lty = c(1,1,1,1,1,1,1)
text(x = c(0.91, 0.98, 0.95, 1.02, 1, 1.06, 1.04), y = c(tpa.table.for.plot.sort.2010$tpa.means.2010[1], tpa.table.for.plot.sort.2010$tpa.means.2010[2] + 0.02, tpa.table.for.plot.sort.2010$tpa.means.2010[3], tpa.table.for.plot.sort.2010$tpa.means.2010[4] -.05, tpa.table.for.plot.sort.2010$tpa.means.2010[5] + 0.08, tpa.table.for.plot.sort.2010$tpa.means.2010[6] + 0.04, tpa.table.for.plot.sort.2010$tpa.means.2010[7]+0.04), labels = c("stw", "osp", "cha", "cre", "hor", "yan", "ntw"), cex = 1)

#lamella number
plot(tpl.table.for.plot.sort.2010$tpl.means.2010 ~ c(seq(0.94, 1.06, by = 0.02)), xlim = c(0.9, 1.6), ylim = c(22.9, 24.7), pch = c(1, 15, 16, 1, 1, 16, 16), col = c("grey50", "black", rep("grey50", 5)), axes = FALSE, ylab = "", cex = 1.2)
segments(x0 = seq(0.94, 1.06, by = 0.02), x1 = seq(0.94, 1.06, by = 0.02), y0 = tpl.table.for.plot.sort.2010$tpl.means.2010, y1 = tpl.table.for.plot.sort.2010$tpl.means.2010 + tpl.table.for.plot.sort.2010$tpl.ses.2010, col = c("grey50", "black", rep("grey50", 5)))
segments(x0 = seq(0.94, 1.06, by = 0.02), x1 = seq(0.94, 1.06, by = 0.02), y0 = tpl.table.for.plot.sort.2010$tpl.means.2010, y1 = tpl.table.for.plot.sort.2010$tpl.means.2010 - tpl.table.for.plot.sort.2010$tpl.ses.2010, col = c("grey50", "black", rep("grey50", 5)))
axis(side = 2, at = c(22.9, 23.4, 23.9, 24.4), las = 1)
mtext(text = "Mean Lamella Number", side = 2, line = 2.5, at = 23.65)

points(x = seq(1.44, 1.56, by = 0.02), y = tpl.table.for.plot.sort.2019$tpl.means.2019, col = c(rep("grey50", 4), "black", rep("grey50", 2)), pch = c(1, 1, 16, 1, 15, 16, 16), cex = 1.2)
segments(x0 = seq(1.44, 1.56, by = 0.02), x1 = seq(1.44, 1.56, by = 0.02), y0 = tpl.table.for.plot.sort.2019$tpl.means.2019, y1 = tpl.table.for.plot.sort.2019$tpl.means.2019 + tpl.table.for.plot.sort.2019$tpl.ses.2019, col = c(rep("grey50", 4), "black", rep("grey50", 2)))
segments(x0 = seq(1.44, 1.56, by = 0.02), x1 = seq(1.44, 1.56, by = 0.02), y0 = tpl.table.for.plot.sort.2019$tpl.means.2019, y1 = tpl.table.for.plot.sort.2019$tpl.means.2019 - tpl.table.for.plot.sort.2019$tpl.ses.2019, col = c(rep("grey50", 4), "black", rep("grey50", 2)))
segments(x0 = seq(0.94, 1.06, by = 0.02), x1 = seq(1.44, 1.56, by = 0.02)[c(1, 5, 3, 2, 4, 6, 7)], y0 = tpl.table.for.plot.sort.2010$tpl.means.2010, y1 = tpl.table.for.plot.sort.2019$tpl.means.2019[c(1, 5, 3, 2, 4, 6, 7)], lty = rep(1, 7), lwd = c(0.5, 1.0, rep(0.5, 5)))
#lty = c(1,1,1,1,1,1,3))
text(x = 1, y = 22.9, labels = "2010"); text(x = 1.5, y = 22.9, labels = "2019")
text(x = c(0.90, 0.99, 0.94, 1.50, 0.98, 1.01, 1.60), y = c(tpl.table.for.plot.sort.2010$tpl.means.2010[1], tpl.table.for.plot.sort.2010$tpl.means.2010[2]+.05, tpl.table.for.plot.sort.2010$tpl.means.2010[3], tpl.table.for.plot.sort.2019$tpl.means.2019[2], tpl.table.for.plot.sort.2010$tpl.means.2010[5], tpl.table.for.plot.sort.2010$tpl.means.2010[6]+0.1, tpl.table.for.plot.sort.2019$tpl.means.2019[7]), labels = c("hor", "stw", "cha", "osp", "cre", "ntw", "yan"), cex = 1)

# end Figure 2


##### Figure S1 2010 perch height, toepad, and lamella data, by island ##### 
#for all islands, to show where significance came from back then.

#perch height
ph.2010all <- ph.2010.final.t #from 0a) 
ph.2010all.means <- round(tapply(X = ph.2010all$perch.lizard, INDEX = ph.2010all$island, FUN = mean, na.rm = TRUE), 1)
ph.2010all.sds <- round(tapply(X = ph.2010all$perch.lizard, INDEX = ph.2010all$island, FUN = sd, na.rm = TRUE), 1)
ph.2010all.n <- tapply(X = ph.2010all$perch.lizard, INDEX = ph.2010all$island, FUN = length)
ph.2010all.ses <- round(ph.2010all.sds/sqrt(ph.2010all.n), 1)
ph.table.for.plot.2010all <- data.frame(ph.2010all.means, ph.2010all.sds, ph.2010all.n, ph.2010all.ses)
ph.table.for.plot.2010all.sort <- ph.table.for.plot.2010all[order(ph.table.for.plot.2010all$ph.2010all.means), ]
ph.table.for.plot.2010all.sort <- ph.table.for.plot.2010all.sort[-which(rownames(ph.table.for.plot.2010all.sort)=="Line.of.Cedars"), ]

#toepad data, all islands. From 0b)
tp.2010.tt.noLoc <- tp.2010.tt[-which(tp.2010.tt$island == "line.cedars") , ]

tp.2010.tt.noLoc.m <- tp.2010.tt.noLoc[which(tp.2010.tt.noLoc$sex == "m") , ]
tp.2010.tt.noLoc.f <- tp.2010.tt.noLoc[which(tp.2010.tt.noLoc$sex == "f") , ]

morph.2010all.sc.m <- f.sizecorrect(standard.length = tp.2010.tt.noLoc.m$svl, vector.of.columns = c(6,7), dataframe.to.use = tp.2010.tt.noLoc.m, watershed = tp.2010.tt.noLoc.m$island)
morph.2010all.sc.f <- f.sizecorrect(standard.length = tp.2010.tt.noLoc.f$svl, vector.of.columns = c(6,7), dataframe.to.use = tp.2010.tt.noLoc.f, watershed = tp.2010.tt.noLoc.f$island)

morph.2010all.sc <- merge(morph.2010all.sc.m, morph.2010all.sc.f, all = TRUE)

#toepad area
ta.2010all.means <- round(tapply(X = morph.2010all.sc$area.mm2.av.sc, INDEX = morph.2010all.sc$island, FUN = mean, na.rm = TRUE), 2)
ta.2010all.sds <- round(tapply(X = morph.2010all.sc$area.mm2.av.sc, INDEX = morph.2010all.sc$island, FUN = sd, na.rm = TRUE), 1)
ta.2010all.n <- tapply(X = morph.2010all.sc$area.mm2.av.sc, INDEX = morph.2010all.sc$island, FUN = length)
ta.2010all.ses <- round(ta.2010all.sds/sqrt(ta.2010all.n), 1)
ta.table.for.plot.2010all <- data.frame(ta.2010all.means, ta.2010all.sds, ta.2010all.n, ta.2010all.ses)
ta.table.for.plot.2010all.sort <- ta.table.for.plot.2010all[order(ta.table.for.plot.2010all$ta.2010all.means), ]

#lamella number
ln.2010all.means <- round(tapply(X = morph.2010all.sc$count.av, INDEX = morph.2010all.sc$island, FUN = mean, na.rm = TRUE), 1)
ln.2010all.sds <- round(tapply(X = morph.2010all.sc$count.av, INDEX = morph.2010all.sc$island, FUN = sd, na.rm = TRUE), 1)
ln.2010all.n <- tapply(X = morph.2010all.sc$count.av, INDEX = morph.2010all.sc$island, FUN = length)
ln.2010all.ses <- round(ln.2010all.sds/sqrt(ln.2010all.n), 1)
ln.table.for.plot.2010all <- data.frame(ln.2010all.means, ln.2010all.sds, ln.2010all.n, ln.2010all.ses)
ln.table.for.plot.2010all.sort <- ln.table.for.plot.2010all[order(ln.table.for.plot.2010all$ln.2010all.means), ]


#Plot
par(mar = c(1,4,0,0), mfrow = c(3,1), oma = c(1,1,1,1))

#Perch Height
plot(x = seq(1, nrow(ph.table.for.plot.2010all.sort), by = 1), y = ph.table.for.plot.2010all.sort$ph.2010all.means, axes = FALSE, ylab = "", xlab = "", pch = c(rep(1, 5), rep(16, 6)), ylim = c(min(ph.table.for.plot.2010all.sort$ph.2010all.means - ph.table.for.plot.2010all.sort$ph.2010all.ses), max(ph.table.for.plot.2010all.sort$ph.2010all.means + ph.table.for.plot.2010all.sort$ph.2010all.ses)), cex = 1.3)
segments(x0 = seq(1, nrow(ph.table.for.plot.2010all.sort)), x1 = seq(1, nrow(ph.table.for.plot.2010all.sort)), y0 = ph.table.for.plot.2010all.sort$ph.2010all.means, y1 = ph.table.for.plot.2010all.sort$ph.2010all.means + ph.table.for.plot.2010all.sort$ph.2010all.ses)
segments(x0 = seq(1, nrow(ph.table.for.plot.2010all.sort)), x1 = seq(1, nrow(ph.table.for.plot.2010all.sort)), y0 = ph.table.for.plot.2010all.sort$ph.2010all.means, y1 = ph.table.for.plot.2010all.sort$ph.2010all.means - ph.table.for.plot.2010all.sort$ph.2010all.ses)
axis(side = 2, at = c(75, 100, 125, 150, 175, 200, 225))
text(x = c(seq(1, nrow(ph.table.for.plot.2010all.sort)-1, by = 1) + 0.4, nrow(ph.table.for.plot.2010all.sort)-0.4), y = ph.table.for.plot.2010all.sort$ph.2010all.means, labels = c("hor", "osp", "stw", "", "cre", "yan", "", "ntw", "", "", "cha"), cex = 1.1)
mtext(side = 2, line = 2, text = "Perch Height")

#toepad area
plot(x = seq(1, nrow(ta.table.for.plot.2010all.sort), by = 1), y = ta.table.for.plot.2010all.sort$ta.2010all.means, axes = FALSE, ylab = "", xlab = "", pch = c(1,1,1,16,16,1,16,16,1,16,16), ylim = c(min(ta.table.for.plot.2010all.sort$ta.2010all.means - ta.table.for.plot.2010all.sort$ta.2010all.ses), max(ta.table.for.plot.2010all.sort$ta.2010all.means + ta.table.for.plot.2010all.sort$ta.2010all.ses)), cex = 1.3)
segments(x0 = seq(1, nrow(ta.table.for.plot.2010all.sort)), x1 = seq(1, nrow(ta.table.for.plot.2010all.sort)), y0 = ta.table.for.plot.2010all.sort$ta.2010all.means, y1 = ta.table.for.plot.2010all.sort$ta.2010all.means + ta.table.for.plot.2010all.sort$ta.2010all.ses)
segments(x0 = seq(1, nrow(ta.table.for.plot.2010all.sort)), x1 = seq(1, nrow(ta.table.for.plot.2010all.sort)), y0 = ta.table.for.plot.2010all.sort$ta.2010all.means, y1 = ta.table.for.plot.2010all.sort$ta.2010all.means - ta.table.for.plot.2010all.sort$ta.2010all.ses)
axis(side = 2, at = c(3.0, 3.3, 3.6, 3.9, 4.2, 4.5))
text(x = c(seq(1, nrow(ta.table.for.plot.2010all.sort)-1, by = 1) + 0.4, nrow(ta.table.for.plot.2010all.sort)-0.4), y = ta.table.for.plot.2010all.sort$ta.2010all.means, labels = c("", "stw", "osp", "", "cha", "cre", "", "", "hor", "yan", "ntw"), cex = 1.1)
mtext(side = 2, line = 2, text = "Toepad Area (mm2)")

#lamella number
plot(x = seq(1, nrow(ln.table.for.plot.2010all.sort), by = 1), y = ln.table.for.plot.2010all.sort$ln.2010all.means, axes = FALSE, ylab = "", xlab = "", pch = c(1,1,16,1,1,16,16,1,16,16,16), ylim = c(min(ln.table.for.plot.2010all.sort$ln.2010all.means - ln.table.for.plot.2010all.sort$ln.2010all.ses), max(ln.table.for.plot.2010all.sort$ln.2010all.means + ln.table.for.plot.2010all.sort$ln.2010all.ses)), cex = 1.3)
segments(x0 = seq(1, nrow(ln.table.for.plot.2010all.sort)), x1 = seq(1, nrow(ln.table.for.plot.2010all.sort)), y0 = ln.table.for.plot.2010all.sort$ln.2010all.means, y1 = ln.table.for.plot.2010all.sort$ln.2010all.means + ln.table.for.plot.2010all.sort$ln.2010all.ses)
segments(x0 = seq(1, nrow(ln.table.for.plot.2010all.sort)), x1 = seq(1, nrow(ln.table.for.plot.2010all.sort)), y0 = ln.table.for.plot.2010all.sort$ln.2010all.means, y1 = ln.table.for.plot.2010all.sort$ln.2010all.means - ln.table.for.plot.2010all.sort$ln.2010all.ses)
axis(side = 2, at = c(22.5, 23.0, 23.5, 24, 24.5))
text(x = c(seq(1, nrow(ln.table.for.plot.2010all.sort)-1, by = 1) + 0.4, nrow(ln.table.for.plot.2010all.sort)-0.4), y = ln.table.for.plot.2010all.sort$ln.2010all.means, labels = c("hor", "stw", "cha", "osp", "", "", "", "cre", "", "ntw", "yan"), cex = 1.1)
mtext(side = 2, line = 2, text = "Lamella Number")
legend(x = 8, y = 23.5, legend = c("one species", "two species"), pch = c(1, 16), cex = 1)

mtext(text = "2010", side = 3, outer = TRUE, line = -1)
# end Figure S1

