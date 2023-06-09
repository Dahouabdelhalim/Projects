################################to calculate effect sizes, yi and vi, can used lnRR or SMDH (used SMDH here)##################################
#mydat <- read.table(#DATA Final dataset#)
str(mydat)
names(mydat)


#to export
hedgesg<- escalc(measure = "SMDH", n1i = mydat$SYMB.N, n2i = mydat$CONT.N, m1i = mydat$SYMBIONT, 
                 m2i = mydat$CONTROL, sd1i = mydat$SYMB.SD, sd2i = mydat$CONT.SD)
write.table (hedgesg, file = "effect sizes", append = FALSE, quote = TRUE, sep = " ",
             eol = "\\n", na = "NA", dec = ".", row.names = TRUE,
             col.names = TRUE, qmethod = c("escape", "double"))

#to append
hedgesg<-escalc(measure = "SMDH", n1i = SYMB.N, n2i = CONT.N, m1i = SYMBIONT, 
                m2i = CONTROL, sd1i = SYMB.SD, sd2i = CONT.SD, data = mydat,
                append = TRUE)

############################################################### POWER ANALYSIS ###############################################################################################
#all.data <- read.table(#DATA Final dataset#)
attach(all.data)
names(all.data)

#apriori effect size 0.5
es <- (tapply(yi.SMDH, list(Dataset,RESPONSE),  length) / tapply(yi.SMDH, list(Dataset,RESPONSE),  length))/2; es
ave.N<- (CONT.N+SYMB.N)/2
as <-tapply(ave.N, list(Dataset,RESPONSE),  mean);as
mk <- tapply(yi.SMDH, list(Dataset,RESPONSE),  length);mk
hg <- (tapply(yi.SMDH, list(Dataset,RESPONSE),  length) / tapply(yi.SMDH, list(Dataset,RESPONSE),  length))/3; hg

eq1 <- ((as+as)/((as)*(as))) + ((es^2)/(2*(as+as)))
eq2 <- hg*(eq1)
eq3 <- eq2+eq1
eq4 <- eq3/mk
eq5 <- (es/sqrt(eq4))
Power <- (1-pnorm(1.96-eq5))*100 # Two-tailed
Power


#BY NAT VS EXPT
es<- (tapply(yi.SMDH, list(Dataset,RESPONSE,Nat.vs.Exp),  length) / tapply(yi.SMDH, list(Dataset,RESPONSE,Nat.vs.Exp),  length))/2; es
ave.N<- (CONT.N+SYMB.N)/2
as <-tapply(ave.N, list(Dataset,RESPONSE,Nat.vs.Exp),  mean);as
mk <- tapply(yi.SMDH, list(Dataset,RESPONSE,Nat.vs.Exp),  length);mk
hg <- (tapply(yi.SMDH, list(Dataset,RESPONSE,Nat.vs.Exp),  length) / tapply(yi.SMDH, list(Dataset,RESPONSE,Nat.vs.Exp),  length))/3; hg

eq1 <- ((as+as)/((as)*(as))) + ((es^2)/(2*(as+as)))
eq2 <- hg*(eq1)
eq3 <- eq2+eq1
eq4 <- eq3/mk
eq5 <- (es/sqrt(eq4))
Power <- (1-pnorm(1.96-eq5))*100 # Two-tailed
Power


#BY SYMBIONT
es<- (tapply(yi.SMDH, list(Dataset,RESPONSE,Symbiont),  length) / tapply(yi.SMDH, list(Dataset,RESPONSE,Symbiont),  length))/2; es
ave.N<- (CONT.N+SYMB.N)/2
as <-tapply(ave.N, list(Dataset,RESPONSE,Symbiont),  mean);as
mk <- tapply(yi.SMDH, list(Dataset,RESPONSE,Symbiont),  length);mk
hg <- (tapply(yi.SMDH, list(Dataset,RESPONSE,Symbiont),  length) / tapply(yi.SMDH, list(Dataset,RESPONSE,Symbiont),  length))/3; hg

eq1 <- ((as+as)/((as)*(as))) + ((es^2)/(2*(as+as)))
eq2 <- hg*(eq1)
eq3 <- eq2+eq1
eq4 <- eq3/mk
eq5 <- (es/sqrt(eq4))
Power <- (1-pnorm(1.96-eq5))*100 # Two-tailed
Power


#BY INSECT SPECIES
es<- (tapply(yi.SMDH, list(Dataset,RESPONSE,Insect.Species),  length) / tapply(yi.SMDH, list(Dataset,RESPONSE,Insect.Species),  length))/2; es
ave.N<- (CONT.N+SYMB.N)/2
as <-tapply(ave.N, list(Dataset,RESPONSE,Insect.Species),  mean);as
mk <- tapply(yi.SMDH, list(Dataset,RESPONSE,Insect.Species),  length);mk
hg <- (tapply(yi.SMDH, list(Dataset,RESPONSE,Insect.Species),  length) / tapply(yi.SMDH, list(Dataset,RESPONSE,Insect.Species),  length))/3; hg

eq1 <- ((as+as)/((as)*(as))) + ((es^2)/(2*(as+as)))
eq2 <- hg*(eq1)
eq3 <- eq2+eq1
eq4 <- eq3/mk
eq5 <- (es/sqrt(eq4))
Power <- (1-pnorm(1.96-eq5))*100 # Two-tailed
Power

########################################################## final data set  ###################################################################
#all.data <- read.table(#DATA Final dataset#)
names(all.data)
library(metafor)
library(dmetar)
library(meta)

aphid.data<- subset(all.data, Dataset=="Aphid")
aphid.metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, random = ~ 1 |Study.ID, method = "REML", data = aphid.data)
summary(aphid.metareg)

funnel(rma(yi.SMDH, vi.SMDH, mods= ~ Study.ID+RESPONSE-1, method = "REML", data = aphid.data))
regtest(rma(yi.SMDH, vi.SMDH, mods= ~ Study.ID+RESPONSE-1, method = "REML", data = aphid.data)) 
#test for funnel plot asymmetry: z = 1.4344, p = 0.1514
aphid.metagen <- metagen(TE=yi.SMDH,seTE=vi.SMDH, data=aphid.data,studlab=paste(Study.ID),
                  comb.fixed = FALSE,comb.random = TRUE, method.tau = "REML",hakn = TRUE,prediction = TRUE,sm = "SMD")
pcurve(aphid.metagen)




whitefly.data<- subset(all.data, Dataset=="Whitefly")
whitefly.metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, random = ~ 1 |Study.ID, method = "REML", data = whitefly.data)
summary(whitefly.metareg)

funnel(rma(yi.SMDH, vi.SMDH, mods= ~ Study.ID+RESPONSE-1, method = "REML", data = whitefly.data))
regtest(rma(yi.SMDH, vi.SMDH, mods= ~ Study.ID+RESPONSE-1, method = "REML", data = whitefly.data)) 
#test for funnel plot asymmetry: z = -1.0788, p = 0.2807
whitefly.metagen <- metagen(TE=yi.SMDH,seTE=vi.SMDH, data=whitefly.data,studlab=paste(Study.ID),
                         comb.fixed = FALSE,comb.random = TRUE, method.tau = "REML",hakn = TRUE,prediction = TRUE,sm = "SMD")
pcurve(whitefly.metagen)




Heteroptera.data<- subset(all.data, Dataset=="Heteroptera")
Heteroptera.metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, random = ~ 1 |Study.ID, method = "REML", data = Heteroptera.data)
summary(Heteroptera.metareg)

funnel(rma(yi.SMDH, vi.SMDH, mods= ~ Study.ID+RESPONSE-1, method = "REML", data = Heteroptera.data))
regtest(rma(yi.SMDH, vi.SMDH, mods= ~ Study.ID+RESPONSE-1, method = "REML", data = Heteroptera.data)) 
#test for funnel plot asymmetry: z = -0.3991, p = 0.6899
Heteroptera.metagen <- metagen(TE=yi.SMDH,seTE=vi.SMDH, data=Heteroptera.data,studlab=paste(Study.ID),
                         comb.fixed = FALSE,comb.random = TRUE, method.tau = "REML",hakn = TRUE,prediction = TRUE,sm = "SMD")
pcurve(Heteroptera.metagen)

########################################################## APHIDS  ############################################################################
#all.data <- read.table(#DATA Final dataset#)
aphid.data<- subset(all.data, Dataset=="Aphid")
library(metafor)
aphid.metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, random = ~ 1 |Study.ID, method = "REML", data = aphid.data)
summary(aphid.metareg)

#Compare aphid experimental and natural lines
expt.data<- subset(aphid.data, Nat.vs.Exp=="Expt");tapply(expt.data$yi.SMDH,expt.data$RESPONSE,  length)
nat.data<- subset(aphid.data, Nat.vs.Exp=="Nat");tapply(nat.data$yi.SMDH,nat.data$RESPONSE,  length)

metareg.expt <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, random = ~ 1 |Study.ID, method = "REML", data = expt.data)
summary(metareg.expt)
funnel(rma(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, method = "REML", data = expt.data), main ="Experimental")
regtest(rma(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, method = "REML", data = expt.data))
#test for funnel plot asymmetry: z = 0.7386, p = 0.4601


metareg.nat <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, random = ~ 1 |Study.ID, method = "REML", data = nat.data)
summary(metareg.nat)
funnel(rma(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, method = "REML", data = nat.data), main ="Natural")
regtest(rma(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, method = "REML", data = nat.data)) 
#test for funnel plot asymmetry: z = 1.8014, p = 0.0716

tapply(expt.data$yi.SMDH,list(expt.data$RESPONSE,expt.data$Symbiont),  length)
tapply(nat.data$yi.SMDH,list(nat.data$RESPONSE,nat.data$Symbiont),  length)


#analysis of individual response variables
names(aphid.data)
Bodysize.data<- subset(aphid.data, RESPONSE=="Bodysize")
Dev.time.data<- subset(aphid.data, RESPONSE=="Dev.time")
Fecundity.data<- subset(aphid.data, RESPONSE=="Fecundity")
Lifespan.data<- subset(aphid.data, RESPONSE=="Lifespan")
Parasitism.data<- subset(aphid.data, RESPONSE=="Parasitism")


metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Insect.Species-1, random = ~ 1 |Study.ID, method = "REML", data = Bodysize.data)
summary(metareg)
tapply(Bodysize.data$yi.SMDH,Bodysize.data$Insect.Species,  length)
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symbiont-1, random = ~ 1 |Study.ID, method = "REML", data = Bodysize.data)
summary(metareg)
tapply(Bodysize.data$yi.SMDH,Bodysize.data$Symbiont,  length)


metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Insect.Species-1, random = ~ 1 |Study.ID, method = "REML", data = Dev.time.data)
summary(metareg)
tapply(Dev.time.data$yi.SMDH,Dev.time.data$Insect.Species,  length)
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symbiont-1, random = ~ 1 |Study.ID, method = "REML", data = Dev.time.data)
summary(metareg)
tapply(Dev.time.data$yi.SMDH,Dev.time.data$Symbiont,  length)


metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Insect.Species-1, random = ~ 1 |Study.ID, method = "REML", data = Fecundity.data)
summary(metareg)
tapply(Fecundity.data$yi.SMDH,Fecundity.data$Insect.Species,  length)
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symbiont-1, random = ~ 1 |Study.ID, method = "REML", data = Fecundity.data)
summary(metareg)
tapply(Fecundity.data$yi.SMDH,Fecundity.data$Symbiont,  length)

metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Nat.vs.Exp-1, random = ~ 1 |Study.ID, method = "REML", data = Lifespan.data)
summary(metareg)
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Insect.Species-1, random = ~ 1 |Study.ID, method = "REML", data = Lifespan.data)
summary(metareg)
tapply(Lifespan.data$yi.SMDH,Lifespan.data$Insect.Species,  length)
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symbiont-1, random = ~ 1 |Study.ID, method = "REML", data = Lifespan.data)
summary(metareg)
tapply(Lifespan.data$yi.SMDH,Lifespan.data$Symbiont,  length)


metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Insect.Species-1, random = ~ 1 |Study.ID, method = "REML", data = Parasitism.data)
summary(metareg)
tapply(Parasitism.data$yi.SMDH,Parasitism.data$Insect.Species,  length)
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symbiont-1, random = ~ 1 |Study.ID, method = "REML", data = Parasitism.data)
summary(metareg)
tapply(Parasitism.data$yi.SMDH,Parasitism.data$Symbiont,  length)



############################################## HOST PLANT EFFECTS ##############################################################
library(metafor)
#plant.data <- read.table(#DATA Host plant effects_aphid#)
names(plant.data)
tapply(plant.data$yi.SMDH,plant.data$RESPONSE,  length)

tapply(plant.data$yi.SMDH, list(plant.data$Expt.host.plant,plant.data$RESPONSE),  length)
tapply(plant.data$yi.SMDH, list(plant.data$Expt.host.plant,plant.data$RESPONSE,plant.data$Symbiont),  length)

#insufficient data for all plant hosts, use only M sativa and V faba comparison across Lifespan, Fecundity and Parasitism
expthost.summary<- subset(plant.data, Expt.host.plant %in% c("Medicago.sativa", "Vicia.faba")) #N=120
tapply(expthost.summary$yi.SMDH, list(expthost.summary$Plant.host.of.aphid.collection,expthost.summary$RESPONSE),  length)

expthost.data<- subset(expthost.summary, Plant.host.of.aphid.collection =="Medicago.sativa") #N=40
tapply(expthost.data$yi.SMDH, list(expthost.data$RESPONSE,expthost.data$Expt.host.plant),  length)

metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, random = ~ 1 |Study.ID, method = "REML", data = expthost.data)
summary(metareg)
tapply(expthost.data$yi.SMDH, expthost.data$RESPONSE,  length)

lifespan.EHP.data<- subset(expthost.data, RESPONSE=="Lifespan") #n=6
fecundity.EHP.data<- subset(expthost.data, RESPONSE=="Fecundity") #n=17
parasitism.EHP.data<- subset(expthost.data, RESPONSE=="Parasitism") #n=8


tapply(lifespan.EHP.data$yi.SMDH,list(lifespan.EHP.data$Expt.host.plant,lifespan.EHP.data$Symbiont),  length)
#insufficient data to compare within symbionts
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Expt.host.plant-1, random = ~ 1 |Study.ID, method = "REML", data = lifespan.EHP.data)
summary(metareg)
tapply(lifespan.EHP.data$yi.SMDH,lifespan.EHP.data$Expt.host.plant,  length)


tapply(fecundity.EHP.data$yi.SMDH,list(fecundity.EHP.data$Expt.host.plant,fecundity.EHP.data$Symbiont),  length)
#insufficient data to compare within symbionts
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Expt.host.plant-1, random = ~ 1 |Study.ID, method = "ML", data = fecundity.EHP.data)
summary(metareg)
tapply(fecundity.EHP.data$yi.SMDH,fecundity.EHP.data$Expt.host.plant,  length)


tapply(parasitism.EHP.data$yi.SMDH,list(parasitism.EHP.data$Expt.host.plant,parasitism.EHP.data$Symbiont),  length)
#insufficient data to compare within symbionts
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Expt.host.plant-1, random = ~ 1 |Study.ID, method = "REML", data = parasitism.EHP.data)
summary(metareg)
tapply(parasitism.EHP.data$yi.SMDH,parasitism.EHP.data$Expt.host.plant,  length)


######################################################### Plant.host.of.aphid.collection OHP #############################################################
#plant.data <- read.table(#DATA Host plant effects_aphid#)
names(plant.data)
tapply(plant.data$yi.SMDH,plant.data$RESPONSE,  length)

tapply(plant.data$yi.SMDH, list(plant.data$Plant.host.of.aphid.collection,plant.data$RESPONSE),  length)

OHP.data<- subset(plant.data, Plant.host.of.aphid.collection %in% c("Medicago.sativa","Ononis.spinosa")) #N=68

lifespan.OHP.data<- subset(OHP.data, RESPONSE=="Lifespan") #n=13
fecundity.OHP.data<- subset(OHP.data, RESPONSE=="Fecundity") #n=28
parasitism.OHP.data<- subset(OHP.data, RESPONSE=="Parasitism") #n=14


metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Plant.host.of.aphid.collection-1, random = ~ 1 |Study.ID, method = "REML", data = lifespan.OHP.data)
summary(metareg)

metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Plant.host.of.aphid.collection-1, random = ~ 1 |Study.ID, method = "REML", data = fecundity.OHP.data)
summary(metareg)

metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Plant.host.of.aphid.collection-1, random = ~ 1 |Study.ID, method = "REML", data = parasitism.OHP.data)
summary(metareg)


##################################################### ORIGIN HOST PLANT HDEF ONLY FOR A PISUM APHIDS ############################################################
OHP.hdef.all<- subset(plant.data, Symbiont=="H.def") #N=34
tapply(OHP.hdef.all$yi.SMDH,list(OHP.hdef.all$Plant.host.of.aphid.collection,OHP.hdef.all$RESPONSE),  length)
tapply(OHP.hdef.all$yi.SMDH,list(OHP.hdef.all$Plant.host.of.aphid.collection,OHP.hdef.all$Aphid.Species),  length) #only A pisum
tapply(OHP.hdef.all$yi.SMDH,list(OHP.hdef.all$Expt.host.plant,OHP.hdef.all$RESPONSE),  length)

OHP.data<- subset(plant.data, Plant.host.of.aphid.collection %in% c("Medicago.sativa","Ononis.spinosa")) #N=68
OHP.hdef<- subset(OHP.data, Symbiont=="H.def") #N=34
tapply(OHP.hdef$yi.SMDH,list(OHP.hdef$Plant.host.of.aphid.collection,OHP.hdef$RESPONSE),  length)
tapply(OHP.hdef$yi.SMDH,list(OHP.hdef$Plant.host.of.aphid.collection,OHP.hdef$Aphid.Species),  length) #only A pisum
tapply(OHP.hdef$yi.SMDH,list(OHP.hdef$Expt.host.plant,OHP.hdef$RESPONSE),  length)

OHP.Hdef.Vfaba<- subset(OHP.hdef, Expt.host.plant=="Vicia.faba")


fecundity.OHP.hdef<- subset(OHP.Hdef.Vfaba, RESPONSE=="Fecundity") #n=9
parasitism.OHP.hdef<- subset(OHP.Hdef.Vfaba, RESPONSE=="Parasitism") #n=10

metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Plant.host.of.aphid.collection-1, random = ~ 1 |Study.ID, method = "REML", data = fecundity.OHP.hdef)
summary(metareg)
tapply(fecundity.OHP.hdef$yi.SMDH,fecundity.OHP.hdef$Plant.host.of.aphid.collection,  length)

metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Plant.host.of.aphid.collection-1, random = ~ 1 |Study.ID, method = "REML", data = parasitism.OHP.hdef)
summary(metareg)
tapply(parasitism.OHP.hdef$yi.SMDH,list(parasitism.OHP.hdef$Plant.host.of.aphid.collection, parasitism.OHP.hdef$Symbiont),  length)


############################################## MULTIPLE HOSTING OF SYMBIONTS ##############################################################

library(metafor)
#multi.data <- read.table(#DATA Aphid multiple hosting#)
names(multi.data)
paired<- subset(multi.data, Dataset=="pair") #n=88

lifespan.multi<- subset(paired, RESPONSE=="Lifespan") #n=16 (4 NA)
fecundity.multi<- subset(paired, RESPONSE=="Fecundity") #n=40
parasitism.multi<- subset(paired, RESPONSE=="Parasitism") #n=20

metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symb.focal-1, random = ~ 1 |Study.ID, method = "REML", data = lifespan.multi)
summary(metareg)
tapply(lifespan.multi$yi.SMDH,lifespan.multi$Symb.focal,  length)
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symb.pair-1, random = ~ 1 |Study.ID, method = "REML", data = lifespan.multi)
summary(metareg)
tapply(lifespan.multi$yi.SMDH,lifespan.multi$Symb.pair,  length)

metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symb.focal-1, random = ~ 1 |Study.ID, method = "REML", data = fecundity.multi)
summary(metareg)
tapply(fecundity.multi$yi.SMDH,fecundity.multi$Symb.focal,  length)
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symb.pair-1, random = ~ 1 |Study.ID, method = "REML", data = fecundity.multi)
summary(metareg)
tapply(fecundity.multi$yi.SMDH,fecundity.multi$Symb.pair,  length)

metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symb.focal-1, random = ~ 1 |Study.ID, method = "REML", data = parasitism.multi)
summary(metareg)
tapply(parasitism.multi$yi.SMDH,parasitism.multi$Symb.focal,  length)
metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ Symb.pair-1, random = ~ 1 |Study.ID, method = "REML", data = parasitism.multi)
summary(metareg)
tapply(parasitism.multi$yi.SMDH,parasitism.multi$Symb.pair,  length)


Hdef.comparison<- subset(multi.data, Symb.focal=="H.def") #n=26

metareg <- rma.mv(yi.SMDH, vi.SMDH, mods= ~ RESPONSE-1, random = ~ 1 |Study.ID, method = "REML", data = Hdef.comparison)
summary(metareg)
tapply(Hdef.comparison$yi.SMDH, Hdef.comparison$RESPONSE,  length)
