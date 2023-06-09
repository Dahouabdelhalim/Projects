### Analysis used in Backhouse et al 2021 "Differential geographic patterns in song components of male Albert's lyrebirds" ###

#Measurements were taken from individual elements in songs using Raven
#Measurements for individual elements are in file "ALB.elements.csv"

ALB.elements <- read.csv("ALB.elements.csv")
#Dataset includes principal components 1-3 from PCA in Minitab

## Create datasets for each song component and the full songs ##

#Introductory elements
intro.elements <- ALB.elements[(ALB.elements$Type=="intro"),]
intro.elements$PC1 <- NULL
intro.elements$PC2 <- NULL
intro.elements$PC3 <- NULL
intro.elements$log.dur <- log(intro.elements$Dur.90)
#write.csv(intro.elements, "intro.elements.csv")
#NB: PCA was conducted in Minitab to find principal components

#Final elements
final.elements <- ALB.elements[(ALB.elements$Order.pct==100),]
final.elements$PC1 <- NULL
final.elements$PC2 <- NULL
final.elements$PC3 <- NULL
final.elements$log.endfr <- log(4000-final.elements$End.freq)
#write.csv(final.elements, "final.elements.csv")
#NB: PCA was conducted in Minitab to find principal components

#Remove "buzz" elements for use in analysis
final.elements.nobuzz <- final.elements[!(final.elements$Type.b=="buzz"),]
#write.csv(final.elements.nobuzz, "final.elements.nobuzz.csv")
#NB: PCA was conducted in Minitab to find principal components and remove outliers

#Song body
library(plyr)

body.elements <- ALB.elements[(ALB.elements$Type.cat=="body"),]

#summarise element-level measurements to get song-level measurements
song.body <- ddply(body.elements, .(Bird.ID, Begin.File, Song, Location), summarise,
                   Length = (max(End.Time)-min(Begin.Time)), #song length including breaks between elements
                   Freq.range.90 = (max(Freq.95)-min(Freq.5)), #90% frequency range
                   Max.freq = max(Peak.Freq), #maximum peak frequency
                   Min.freq = min(Peak.Freq), #minimum peak frequency
                   CV.pf = sd(Peak.Freq)/mean(Peak.Freq),#coefficient of variation of peak frequency
                   CV.dur = sd(Dur.90)/mean(Dur.90)) #coefficient of variation of element duration
View(song.body)

#Create song slope
#first set the function
lin_fit <- function(dat){
  the_fit <-lm(Center.Freq #frequency measurement
               ~ Order.pct #time measurement
               , dat)
  setNames(coef(the_fit), c("intercept", "slope"))
}

#calculate slope and intercept for each song
body.coefs <- ddply(body.elements, .(Bird.ID, Begin.File, Song, Location), lin_fit)
View(body.coefs)

#merge slopes with rest of song vars
song.body <- merge(song.body, body.coefs, by=c("Location", "Bird.ID", "Begin.File", "Song"))

#write.csv(song.body, "song.body.csv")
#NB: PCA was conducted in Minitab to find principal components and remove outliers

#Full songs
#summarise element-level measurements to get song-level measurements
ALB.songs <- ddply(ALB.elements, .(Bird.ID, Begin.File, Song, Location), summarise,
                   Length = (max(End.Time)-min(Begin.Time)), #song length including breaks between elements
                   Freq.range.90 = (max(Freq.95)-min(Freq.5)), #90% frequency range
                   Max.freq = max(Peak.Freq), #maximum peak frequency
                   Min.freq = min(Peak.Freq), #minimum peak frequency
                   CV.pf = sd(Peak.Freq)/mean(Peak.Freq),#coefficient of variation of peak frequency
                   CV.dur = sd(Dur.90)/mean(Dur.90)) #coefficient of variation of element duration
View(ALB.songs)

#Create song slope
#first set the function
lin_fit <- function(dat){
  the_fit <-lm(Center.Freq #frequency measurement
               ~ Order.pct #time measurement
               , dat)
  setNames(coef(the_fit), c("intercept", "slope"))
}

#calculate slope and intercept for each song
song.coefs <- ddply(ALB.elements, .(Bird.ID, Begin.File, Song, Location), lin_fit)
View(song.coefs)

#merge slopes with rest of song vars
ALB.songs <- merge(ALB.songs, song.coefs, by=c("Location", "Bird.ID", "Begin.File", "Song"))

#write.csv(ALB.songs, "ALB.songs.csv")
#NB: PCA was conducted in Minitab to find principal components and remove outliers


## Run pDFA ##

#Function for permutated discriminant function analysis (pDFA) received from R Mundry (Mundry & Sommers 2007 Anim. Behav.)
#Modified function to save DFA output from original (unpermuted) data in "s.dfa" (full DFA result) and "dfas" (linear discriminants)

#Will need to request function from R Mundry to run this code

library(MASS)

#intro elements
intro.elements <- read.csv("intro.elements.csv")
View(intro.elements)
intro.scaled <- as.data.frame(scale(intro.elements[c(5:6, 8:10)]))
View(intro.scaled) #check variables are correct
#should include BW.90, Peak.Freq, Begin.Freq, End.freq, log.dur
intro.scaled$Location <- intro.elements$Location
intro.scaled$Bird.ID <- intro.elements$Bird.ID

#run pDFA
pdfa.intro <- pDFA.nested(test.fac="Location", contr.fac = "Bird.ID", 
                          variables=c("Peak.Freq", "Begin.freq", "End.freq", "BW.90", "log.dur"),
                          restrict.by=NULL, n.contr.fac.levels.to.sel=NULL,
                          n.to.sel.per.contr.fac.level=NULL, n.sel=100, n.perm=1000, 
                          pdfa.data=intro.scaled)
pdfa.intro$result #permutation result
pdfa.intro$conf.mat.cval #confusion matrix

pdfa.intro$dfas #linear discriminants
pdfa.intro$s.dfa #DFA result

#get discriminant functions for each song to use in Fig. 5 and Fig. 7
intro.df <- predict(pdfa.intro$s.dfa, newdata=intro.scaled); head(intro.df$x, 10)
zz1.i <- data.frame(intro.df$x)
intro.scaled$fn1 <- (zz1.i$LD1); intro.scaled$fn2 <- (zz1.i$LD2); intro.scaled$fn3 <- (zz1.i$LD3);
intro.scaled$fn4 <- (zz1.i$LD4); intro.scaled$fn5 <- (zz1.i$LD5);
View(intro.scaled)

#final notes, without buzz notes
final.elements.nobuzz <- read.csv("final.elements.nobuzz.csv")
View(final.elements.nobuzz)
final.scaled <- as.data.frame(scale(final.elements.nobuzz[c(5:8, 10)]))
#should include BW.90, Peak.Freq, Dur.90, Begin.Freq, log.endfr
View(final.scaled)
final.scaled$Location <- final.elements.nobuzz$Location
final.scaled$Bird.ID <- final.elements.nobuzz$Bird.ID

#run pDFA
pdfa.final <- pDFA.nested(test.fac="Location", contr.fac = "Bird.ID", 
                          variables=c("Peak.Freq", "Begin.freq", "log.endfr", "BW.90", "Dur.90"),
                          restrict.by=NULL, n.contr.fac.levels.to.sel=NULL,
                          n.to.sel.per.contr.fac.level=NULL, n.sel=100, n.perm=1000, 
                          pdfa.data=final.scaled)
pdfa.final$result #permutation result
pdfa.final$conf.mat.cval #confusion matrix

pdfa.final$dfas #linear discriminants
pdfa.final$s.dfa #DFA result

#get discriminant functions for each song to use in Fig. 5 and Fig. 7
final.df <- predict(pdfa.final$s.dfa, newdata=final.scaled); head(final.df$x, 10)
zz1.f <- data.frame(final.df$x)
final.scaled$fn1 <- (zz1.f$LD1); final.scaled$fn2 <- (zz1.f$LD2); final.scaled$fn3 <- (zz1.f$LD3);
final.scaled$fn4 <- (zz1.f$LD4); final.scaled$fn5 <- (zz1.f$LD5);
View(final.scaled)


#Song body
song.body <- read.csv("song.body.csv")
View(song.body)
body.scaled <- as.data.frame(scale(song.body[5:11]))
View(body.scaled)
body.scaled$Location <- song.body$Location
body.scaled$Bird.ID <- song.body$Bird.ID

#run pDFA
pdfa.body <- pDFA.nested(test.fac="Location", contr.fac = "Bird.ID", 
                         variables=c("Length", "Freq.range.90", "Max.freq", "Min.freq", 
                                     "CV.pf", "CV.dur", "slope"),
                         restrict.by=NULL, n.contr.fac.levels.to.sel=NULL,
                         n.to.sel.per.contr.fac.level=NULL, n.sel=100, n.perm=1000, 
                         pdfa.data=body.scaled)
pdfa.body$result
pdfa.body$conf.mat.cval

pdfa.body$dfas
pdfa.body$s.dfa

#get discriminant functions for each song
body.df <- predict(pdfa.body$s.dfa, newdata=body.scaled); head(body.df$x, 10)
zz1.b <- data.frame(body.df$x)
View(zz1.b)
body.scaled$fn1 <- (zz1.b$LD1); body.scaled$fn2 <- (zz1.b$LD2); body.scaled$fn3 <- (zz1.b$LD3);
body.scaled$fn4 <- (zz1.b$LD4); body.scaled$fn5 <- (zz1.b$LD5);
View(body.scaled)


#Full songs
ALB.songs <- read.csv("ALB.songs.csv")
View(ALB.songs)
songs.scaled <- as.data.frame(scale(ALB.songs[c(5:11)]))
View(songs.scaled)
songs.scaled$Location <- ALB.songs$Location
songs.scaled$Bird.ID <- ALB.songs$Bird.ID

#run pDFA
pdfa.songs <- pDFA.nested(test.fac="Location", contr.fac = "Bird.ID", 
                          variables=c("Length", "Freq.range.90", "Max.freq", "Min.freq", 
                                      "CV.pf", "CV.dur", "slope"),
                          restrict.by=NULL, n.contr.fac.levels.to.sel=NULL,
                          n.to.sel.per.contr.fac.level=NULL, n.sel=100, n.perm=1000, 
                          pdfa.data=songs.scaled)
pdfa.songs$result
pdfa.songs$conf.mat.cval

pdfa.songs$dfas
pdfa.songs$s.dfa

#get discriminant functions for each song
songs.df <- predict(pdfa.songs$s.dfa, newdata=songs.scaled); head(songs.df$x, 10)
zz1 <- data.frame(songs.df$x)
songs.scaled$fn1 <- (zz1$LD1); songs.scaled$fn2 <- (zz1$LD2); songs.scaled$fn3 <- (zz1$LD3);
songs.scaled$fn4 <- (zz1$LD4); songs.scaled$fn5 <- (zz1$LD5);
View(songs.scaled)

## Mahalanobis distance ##
#This analysis uses principal components calculated using Minitab
#The provided files in Dryad include these principal components

library(HDMD)
library(cluster)
library(ecodist)

#Intro elements
intro.elements <- read.csv("intro.elements.csv")

group = matrix(intro.elements$Location) #what is being compared
group = t(group[,1]) #prepare for pairwise.mahalanobis function

variables = c("PC1","PC2", "PC3") #variables (what is being used for comparison)
variables = as.matrix(intro.elements[,variables]) #prepare for pairwise.mahalanobis function

mahala_sq = pairwise.mahalanobis(x=variables, grouping=group) #get squared mahalanobis distances (see mahala_sq$distance).
names = rownames(mahala_sq$means) #capture labels

intro.mahala = sqrt(mahala_sq$distance) #mahalanobis distance
rownames(intro.mahala) = names #set rownames in the dissimilarity matrix
colnames(intro.mahala) = names #set colnames in the dissimilarity matrix

intro.mahala
intro.dist <- as.dist(intro.mahala); intro.dist

#Final elements
final.elements <- read.csv("final.elements.nobuzz.csv")
group = matrix(final.elements$Location) #what is being compared
group = t(group[,1]) #prepare for pairwise.mahalanobis function

variables = c("PC1","PC2", "PC3") #variables (what is being used for comparison)
variables = as.matrix(final.elements[,variables]) #prepare for pairwise.mahalanobis function

mahala_sq = pairwise.mahalanobis(x=variables, grouping=group) #get squared mahalanobis distances (see mahala_sq$distance).
names = rownames(mahala_sq$means) #capture labels

final.mahala = sqrt(mahala_sq$distance) #mahalanobis distance
rownames(final.mahala) = names #set rownames in the dissimilarity matrix
colnames(final.mahala) = names #set colnames in the dissimilarity matrix

final.mahala
final.dist <- as.dist(final.mahala); final.dist

#song body
song.body <- read.csv("song.body.csv")

group = matrix(song.body$Location) #what is being compared
group = t(group[,1]) #prepare for pairwise.mahalanobis function

variables = c("PC1","PC2", "PC3") #variables (what is being used for comparison)
variables = as.matrix(song.body[,variables]) #prepare for pairwise.mahalanobis function

mahala_sq = pairwise.mahalanobis(x=variables, grouping=group) #get squared mahalanobis distances (see mahala_sq$distance).
names = rownames(mahala_sq$means) #capture labels

body.mahala = sqrt(mahala_sq$distance) #mahalanobis distance
rownames(body.mahala) = names #set rownames in the dissimilarity matrix
colnames(body.mahala) = names #set colnames in the dissimilarity matrix

body.mahala
body.dist <- as.dist(body.mahala); body.dist

#Full song
ALB.songs <- read.csv("ALB.songs.csv")

group = matrix(ALB.songs$Location) #what is being compared
group = t(group[,1]) #prepare for pairwise.mahalanobis function

variables = c("PC1","PC2", "PC3") #variables (what is being used for comparison)
variables = as.matrix(ALB.songs[,variables]) #prepare for pairwise.mahalanobis function

mahala_sq = pairwise.mahalanobis(x=variables, grouping=group) #get squared mahalanobis distances (see mahala_sq$distance).
names = rownames(mahala_sq$means) #capture labels

songs.mahala = sqrt(mahala_sq$distance) #mahalanobis distance
rownames(songs.mahala) = names #set rownames in the dissimilarity matrix
colnames(songs.mahala) = names #set colnames in the dissimilarity matrix

songs.mahala
songs.dist <- as.dist(songs.mahala); songs.dist

## Mantel tests ##
library(ade4)
#Load geographic distances
distances <- read.csv("distances.csv")
#turn geographic distances in dist matrices
nams <- with(distances, unique(c(as.character(Site.x), as.character(Site.y))))
#straight line distance
SL.dist <- with(distances, structure(SL.dist,
                           Size=length(nams),
                           Labels = nams,
                           Diag = FALSE,
                           Upper = FALSE,
                           method = "user",
                           class= "dist"))
SL.dist

#least cost path distance weighted by resistance (used as measure of geographic separation)
LCP.dist <-  with(distances, structure(resist.length,
                            Size=length(nams),
                            Labels = nams,
                            Diag = FALSE,
                            Upper = FALSE,
                            method = "user",
                            class= "dist"))
LCP.dist

#mantel tests
mantel.rtest(intro.dist, SL.dist, nrepet=999)
mantel.rtest(intro.dist, LCP.dist, nrepet=999)

mantel.rtest(final.dist, SL.dist, nrepet=999)
mantel.rtest(final.dist, LCP.dist, nrepet=999)

mantel.rtest(body.dist, SL.dist, nrepet=999)
mantel.rtest(body.dist, LCP.dist, nrepet=999)

mantel.rtest(songs.dist, SL.dist, nrepet=999)
mantel.rtest(songs.dist, LCP.dist, nrepet=999)
