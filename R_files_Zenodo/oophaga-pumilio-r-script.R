############################################################################
############################################################################
########## 	Evaluating the conspecific attraction hypothesis 	##########
##########    with Oophaga pumilio in cacao plantations at 		##########
##########	   La Selva, Costa Rica -- B. Folt, 2017			##########
############################################################################
############################################################################

### Research question:
### Do immigrants settle closer to residents than expected by chance,
### in ways that are consistent with conspecific attraction?



# Clear the workspace
rm(list=ls())

# Set working directory
setwd("C:/Users/brian.folt/Dropbox/Auburn Ph.D. Conspecific Attraction/R")

# Load the datafiles
plotmaps = read.csv("plot-maps.csv", header=TRUE)
datum = read.csv("pumilio-cacao-v3.csv", header=TRUE)
head(datum)

# Subset data to 12 x 9 area that Mo sampled ~exhaustively
head(datum)
datum = droplevels(subset(datum, 
	datum$X < 12 & datum$X > 0 		# Only observations 0-12 x coords
	& datum$Y < 9 & datum$Y > 0))		# Only obserations 0-9 y coords
summary(datum$X)
summary(datum$Y)


##########################################################
######## Part 0 --				 	##########
######## Are cacao trees distributed uniformly 	##########
######## in the study plots?			      ##########
##########################################################

# Subset to plots, subset to structure, and use 
# Clark and Evan's test to evaluate regularity

library(spatstat)

# Cacao trees first

cacao = matrix(NA,4,2, 
	dimnames=list(c("Plot1","Plot2","Plot3","Plot4"),c("R","P-value")))
for (i in 1:4){		
	plot = subset(plotmaps, Plot==i)
	plot = subset(plot, Structure=="cacao")
	res = clarkevans.test(ppp(plot$X,plot$Y,c(0,12),c(0,9)), 
		correction=c("Donnelly"), alternative=c("greater"))
	cacao[i,c(1:2)] = c(res$statistic,res$p.value)
}
cacao # Uniform distribution in all four plots

# Logs
logs = matrix(NA,4,2, 
	dimnames=list(c("Plot1","Plot2","Plot3","Plot4"),c("R","P-value")))
for (i in 1:4){		
	plot = subset(plotmaps, Plot==i)
	plot = subset(plot, Structure=="log")
	res = clarkevans.test(ppp(plot$X,plot$Y,c(0,12),c(0,9)), 
		correction=c("Donnelly"))
	logs[i,c(1:2)] = c(res$statistic,res$p.value)
}
logs # Random distribution in all four plots


# Pejibaye
pejibaye = matrix(NA,4,2, 
	dimnames=list(c("Plot1","Plot2","Plot3","Plot4"),c("R","P-value")))
for (i in 1:4){		
	plot = subset(plotmaps, Plot==i)
	plot = subset(plot, Structure=="pejibaye")
	res = clarkevans.test(ppp(plot$X,plot$Y,c(0,12),c(0,9)), 
		correction=c("Donnelly"), alternative=c("greater"))
	pejibaye[i,c(1:2)] = c(res$statistic,res$p.value)
}
pejibaye # Uniform in 2 plots, random in two plots


# Buttresses
butts = matrix(NA,4,2, 
	dimnames=list(c("Plot1","Plot2","Plot3","Plot4"),c("R","P-value")))
for (i in 1:4){		
	plot = subset(plotmaps, Plot==i)
	plot = subset(plot, Structure=="buttress")
	res = clarkevans.test(ppp(plot$X,plot$Y,c(0,12),c(0,9)), 
		correction=c("Donnelly"))
	butts[i,c(1:2)] = c(res$statistic,res$p.value)
}
butts # Random distribution in all four plots


##########################################################
######## Part 1 --				 	##########
######## Describe the population structure 	##########
######## of each plot in each season	      ##########
##########################################################

library(gdata)
library(reshape)

### Create a for-loop which tabulates the number of juveniles, females, and males
### present in each plot in each season. 

# These matrices will save counts from the for-loop

strJ = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))
strF = strJ
strM = strJ 

# Initiate the for-loop

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to a plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to a season

group = droplevels(subset(season, Sex == "J"))	### Juveniles

if (length(table(group$FrogNo)) > 0){
strJ[i,j] = length(table(group$FrogNo))
}
else {
strJ[i,j] = 0}

group = droplevels(subset(season, Sex == "F"))	### Females

if (length(table(group$FrogNo)) > 0){
strF[i,j] = length(table(group$FrogNo))
}
else {
strF[i,j] = 0}

group = droplevels(subset(season, Sex == "M"))	### Males

if (length(table(group$FrogNo)) > 0){
strM[i,j] = length(table(group$FrogNo))
}
else {
strM[i,j] = 0}
	}
}

### Reshape and combine the three matrices 
### into a single dataframe called 'structure' (~Population Structure)

strJ = melt(as.data.frame(strJ))
Group = rep("J",length(strJ[,1]))
Plot = as.factor(rep(c(1:4),length(strJ[,1])/4))
(strJ = cbind(Group,Plot,strJ))

strF = melt(as.data.frame(strF))
Group = rep("F",length(strF[,1]))
Plot = as.factor(rep(c(1:4),length(strF[,1])/4))
(strF = cbind(Group,Plot,strF))

strM = melt(as.data.frame(strM))
Group = rep("M",length(strM[,1]))
Plot = as.factor(rep(c(1:4),length(strM[,1])/4))
(strM = cbind(Group,Plot,strM))

structure = rbind(strJ,strF,strM)
colnames(structure) = c("Group","Plot","Season","Number")

head(structure)

detach(package:reshape)
detach(package:gdata)


juvies = subset(structure, Group=="J")
ladies = subset(structure, Group=="F")
guys = subset(structure, Group=="M")

plot(Number ~ Season, data=juvies, pch=as.numeric(Plot))
plot(Number ~ Season, data=ladies, pch=as.numeric(Plot))
plot(Number ~ Season, data=guys, pch=as.numeric(Plot))


### Save this data-frame to the clipboard, copy it into Excel ('pumilio-cacao-v3'), 
### and create a histogram to describe mean abundance in plots ('Figure 1' tab).

write.table(structure, "clipboard", sep="\\t", row.names=FALSE)

# In Excel, sort the data by season and then group.
# Calculate mean abundance (+/- std. dev.) for each age-sex group 
# across plots in each season


### Perform an analysis to test for seasonal variation of age-structure
### across the plots (e.g., maybe juvenile entry is pulsed seasonally)

library(nlme)
library(AICcmodavg)

mod1 = lme(Number ~ 1, random=(~1|Plot), data=structure)			# Null
mod2 = lme(Number ~ Group, random=(~1|Plot), data=structure)		# Group
mod3 = lme(Number ~ Season, random=(~1|Plot), data=structure)		# Season
mod4 = lme(Number ~ Group + Season, random=(~1|Plot), data=structure)	# Group+Seasn
mod5 = lme(Number ~ Group*Season, random=(~1|Plot), data=structure)	# Group*Seasn
summary(mod5)

Cand.models = list(mod1, mod2, mod3, mod4, mod5)
Modnames = c("Null", "Group", "Season", "Group+Season", "Group*Season")

# Compute table
print(aictab(cand.set = Cand.models, 
	modnames = Modnames, second.ord = TRUE), digits = 4)

tab = aictab(cand.set = Cand.models, modnames = Modnames, second.ord = TRUE)

write.table(tab, "clipboard", sep="\\t", row.names=FALSE)


### I copied the table into the Excel file "pumilio-cacao" into a tab "Appendix I"
### This analysis could be included in the manuscript to justify a seasonal effect
### in models describing spatial 

detach(package:nlme)
detach(package:AICcmodavg)



########################################################
############### 		Part 2 -- 	     ###############
############### Are immigrants clustered ###############
############### relative to residents?   ###############
########################################################

#########
######### A) -- 
######### Test this idea by subsetting to different groups and measuring
######### aggregation with the Clark & Evans R statistic, where R > 1
######### indicates uniform distribution, and R < 1 indicates clumping.
#########

library(spatstat) # Access for spatial stats stuff
library(gdata)	# Access library for 'drop.levels()'
library(reshape)  # Acess for 'melt()' function later on

#?clarkevans.test()

# Use for-loops to calculate the R-value and P-value
# for each resident-immigrant grouping
# in each plot, in each season


### Juveniles

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))
P = R

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, Sex %in% c("J"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
R[i,j] = clarkevans.test(ppp(means$X,means$Y,c(0,12),c(0,9)),
	correction=c("Donnelly"))$statistic
P[i,j] = clarkevans.test(ppp(means$X,means$Y,c(0,12),c(0,9)),
	correction=c("Donnelly")nsim=999)$p.value
} 

else {
R[i,j] = NA
P[i,j] = NA
}
	}
}

(jR = R)	# Matrix of Clark-Evans R values for juveniles
(jP = P)	# Matrix of Clark-Evans p-values for juveniles

# I removed the calculation of the P-values from future iterations,
# to decrease run time for this analysis

### Males

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, Sex %in% c("M"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
R[i,j] = clarkevans.test(ppp(means$X,means$Y,c(0,12),c(0,9)),
	correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

(rmR = R)	


### Females

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, Sex %in% c("F"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
R[i,j] = clarkevans.test(ppp(means$X,means$Y,c(0,12),c(0,9)),
	correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

(fR = R)


### Juveniles-Resident Females

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("J","RF"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
R[i,j] = clarkevans.test(ppp(means$X,means$Y,c(0,12),c(0,9)),
	correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R = melt(as.data.frame(R))
Group = rep("J-RF",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(jrfR = cbind(Group,Plot,R))
colnames(jrfR) = c("Group","Plot","Season","R")


### Juveniles-Resident Males

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("J","RM"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
R[i,j] = clarkevans.test(ppp(means$X,means$Y,c(0,12),c(0,9)),
	correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R = melt(as.data.frame(R))
Group = rep("J-RM",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(jrmR = cbind(Group,Plot,R))
colnames(jrmR) = c("Group","Plot","Season","R")


### Immigrant females-Resident females

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("IF","RF"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
R[i,j] = clarkevans.test(ppp(means$X,means$Y,c(0,12),c(0,9)),
	correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R = melt(as.data.frame(R))
Group = rep("IF-RF",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(ifrfR = cbind(Group,Plot,R))
colnames(ifrfR) = c("Group","Plot","Season","R")


### Immigrant females - resident males

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("IF","RM"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
R[i,j] = clarkevans.test(ppp(means$X,means$Y,c(0,12),c(0,9)),
	correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R = melt(as.data.frame(R))
Group = rep("IF-RM",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(ifrmR = cbind(Group,Plot,R))
colnames(ifrmR) = c("Group","Plot","Season","R")


### Immigrant males - resident females

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("IM","RF"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
R[i,j] = clarkevans.test(ppp(means$X,means$Y,c(0,12),c(0,9)),
	correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R = melt(as.data.frame(R))
Group = rep("IM-RF",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(imrfR = cbind(Group,Plot,R))
colnames(imrfR) = c("Group","Plot","Season","R")


### Immigrant males - resident males

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("IM","RM"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
R[i,j] = clarkevans.test(ppp(means$X,means$Y,c(0,12),c(0,9)),
	correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R = melt(as.data.frame(R))
Group = rep("IM-RM",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(imrmR = cbind(Group,Plot,R))
colnames(imrmR) = c("Group","Plot","Season","R")


detach(package:spatstat)
detach(package:gdata)
detach(package:reshape)

clusterdata = na.omit(rbind(jrfR,jrmR,ifrfR,ifrmR,imrfR,imrmR))
clusterdata = droplevels(subset(clusterdata, R != "Inf"))
Type = rep("Observed", length(clusterdata[,1]))
clusterdata = cbind(clusterdata, Type)


#########
######### B) -- 
######### Simulating random centroids and measuring the R statistic
######### At random, values should be ~1, w/ no differences between groups.
#########

library(spatstat)
library(gdata)
library(reshape)

NoSims = 5	# Number different simulations for each grouping

### Juveniles-Resident Females

R0 = matrix(NA, ncol = 5, nrow = 1, 
	dimnames=list(c("na"),c("A","B","C","D","E")))

for (h in 1:NoSims){					# No of simulations

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("J","RF"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
rand = rpoint(length(means[,1]), win=owin(c(0,12),c(0,9)))
R[i,j] = clarkevans.test(rand, correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R0 = rbind(R0,R)

}	# Close simulation replication loop

R = R0[-c(1),]


R = melt(as.data.frame(R))
Group = rep("J-RF",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(jrfR = cbind(Group,Plot,R))
colnames(jrfR) = c("Group","Plot","Season","R")


### Juveniles-Resident Males

R0 = matrix(NA, ncol = 5, nrow = 1, 
	dimnames=list(c("na"),c("A","B","C","D","E")))

for (h in 1:NoSims){					# No of simulations

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("J","RM"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
rand = rpoint(length(means[,1]), win=owin(c(0,12),c(0,9)))
R[i,j] = clarkevans.test(rand, correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R0 = rbind(R0,R)

}	# Close simulation replication loop

R = R0[-c(1),]


R = melt(as.data.frame(R))
Group = rep("J-RM",length(R[,1]))
Plot = rep(c(1:4),(length(R[,1])/4))
(jrmR = cbind(Group,Plot,R))
colnames(jrmR) = c("Group","Plot","Season","R")


### Immigrant females-Resident females

R0 = matrix(NA, ncol = 5, nrow = 1, 
	dimnames=list(c("na"),c("A","B","C","D","E")))

for (h in 1:NoSims){					# No of simulations

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("IF","RF"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
rand = rpoint(length(means[,1]), win=owin(c(0,12),c(0,9)))
R[i,j] = clarkevans.test(rand, correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R0 = rbind(R0,R)

}	# Close simulation replication loop

R = R0[-c(1),]


R = melt(as.data.frame(R))
Group = rep("IF-RF",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(ifrfR = cbind(Group,Plot,R))
colnames(ifrfR) = c("Group","Plot","Season","R")


### Immigrant females - resident males

R0 = matrix(NA, ncol = 5, nrow = 1, 
	dimnames=list(c("na"),c("A","B","C","D","E")))

for (h in 1:NoSims){					# No of simulations

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("IF","RM"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
rand = rpoint(length(means[,1]), win=owin(c(0,12),c(0,9)))
R[i,j] = clarkevans.test(rand, correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R0 = rbind(R0,R)

}	# Close simulation replication loop

R = R0[-c(1),]


R = melt(as.data.frame(R))
Group = rep("IF-RM",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(ifrmR = cbind(Group,Plot,R))
colnames(ifrmR) = c("Group","Plot","Season","R")


### Immigrant males - resident females

R0 = matrix(NA, ncol = 5, nrow = 1, 
	dimnames=list(c("na"),c("A","B","C","D","E")))

for (h in 1:NoSims){					# No of simulations

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("IM","RF"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
rand = rpoint(length(means[,1]), win=owin(c(0,12),c(0,9)))
R[i,j] = clarkevans.test(rand, correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R0 = rbind(R0,R)

}	# Close simulation replication loop

R = R0[-c(1),]


R = melt(as.data.frame(R))
Group = rep("IM-RF",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(imrfR = cbind(Group,Plot,R))
colnames(imrfR) = c("Group","Plot","Season","R")


### Immigrant males - resident males

R0 = matrix(NA, ncol = 5, nrow = 1, 
	dimnames=list(c("na"),c("A","B","C","D","E")))

for (h in 1:NoSims){					# No of simulations

R = matrix(0, nrow = 4, ncol = 5, 
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to each plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors
	
season = subset(plot, Season == j)		# Subset to each season
group = droplevels(subset(season, ImmRes %in% c("IM","RM"))) 
							# Subset to immigr-res groups
if (length(row.names(group)) > 0){
means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=group, mean))
rand = rpoint(length(means[,1]), win=owin(c(0,12),c(0,9)))
R[i,j] = clarkevans.test(rand, correction=c("Donnelly"))$statistic
} 

else {
R[i,j] = NA
}
	}
}

R0 = rbind(R0,R)

}	# Close simulation replication loop

R = R0[-c(1),]


R = melt(as.data.frame(R))
Group = rep("IM-RM",length(R[,1]))
Plot = rep(c(1:4),length(R[,1])/4)
(imrmR = cbind(Group,Plot,R))
colnames(imrmR) = c("Group","Plot","Season","R")


detach(package:spatstat)
detach(package:gdata)
detach(package:reshape)

randomdata = na.omit(rbind(jrfR,jrmR,ifrfR,ifrmR,imrfR,imrmR))
randomdata = droplevels(subset(randomdata, R != "Inf"))
Type = rep("Random", length(randomdata[,1]))
randomdata = cbind(randomdata, Type)



#####
##### C) --
##### Bind the R-statistics for random and observed datasets.
##### And compare the groups statistically.
#####

# Explore statistical tests for deviation from random
# First, bind the OBSERVED and RANDOM data together

dataset = rbind(clusterdata,randomdata) 
	# Now use this dataset for tests

library(nlme)

### Paired tests between observed and random values
### for six hypothesized conspecific attraction groupings:

# J-RF
testJRF = lme(R ~ Type, random=(~1|Plot), data=subset(dataset, Group == "J-RF"))
summary(testJRF)

# J-RM
testJRM = lme(R ~ Type, random=(~1|Plot), data=subset(dataset, Group == "J-RM"))
summary(testJRM)

# IF-RF
testIFRF = lme(R ~ Type, random=(~1|Plot), data=subset(dataset, Group == "IF-RF"))
summary(testIFRF)



# IF-RM
testIFRM = lme(R ~ Type, random=(~1|Plot), data=subset(dataset, Group == "IF-RM"))
summary(testIFRM)

# IM-RF
testIMRF = lme(R ~ Type, random=(~1|Plot), data=subset(dataset, Group == "IM-RF"))
summary(testIMRF)

# IM-RM
testIMRM = lme(R ~ Type, random=(~1|Plot), data=subset(dataset, Group == "IM-RM"))
summary(testIMRM)


### Overall model (i.e., ANOVA) to test for differences among observed states

model = lme(R ~ Group, random=(~1|Plot), data=subset(dataset, Type=="Observed"))
summary(model)
anova(model)


# Create matrices summarizing the results of the analysis
# to produce a neat figure in Excel

SSobs = matrix(nrow=6,ncol=2,c(
	testJRF$coefficients$fixed[1],sqrt(summary(testJRF)$varFix)[1,1],
	testJRM$coefficients$fixed[1],sqrt(summary(testJRM)$varFix)[1,1],
	testIFRF$coefficients$fixed[1],sqrt(summary(testIFRF)$varFix)[1,1],	
	testIFRM$coefficients$fixed[1],sqrt(summary(testIFRM)$varFix)[1,1],
	testIMRF$coefficients$fixed[1],sqrt(summary(testIMRF)$varFix)[1,1],
	testIMRM$coefficients$fixed[1],sqrt(summary(testIMRM)$varFix)[1,1]),
	dimnames=list(c("JRF","JRM","IFRF","IFRM","IMRF","IMRM"),
			c("Parameter","S.E.")), byrow=TRUE)

SSrand = matrix(nrow=6,ncol=2,c(
	testJRF$coefficients$fixed[1]+testJRF$coefficients$fixed[2],sqrt(summary(testJRF)$varFix)[2,2],
	testJRM$coefficients$fixed[1]+testJRM$coefficients$fixed[2],sqrt(summary(testJRM)$varFix)[2,2],
	testIFRF$coefficients$fixed[1]+testIFRF$coefficients$fixed[2],sqrt(summary(testIFRF)$varFix)[2,2],	
	testIFRM$coefficients$fixed[1]+testIFRM$coefficients$fixed[2],sqrt(summary(testIFRM)$varFix)[2,2],
	testIMRF$coefficients$fixed[1]+testIMRF$coefficients$fixed[2],sqrt(summary(testIMRF)$varFix)[2,2],
	testIMRM$coefficients$fixed[1]+testIMRM$coefficients$fixed[2],sqrt(summary(testIMRM)$varFix)[2,2]),
	dimnames=list(c("JRF","JRM","IFRF","IFRM","IMRF","IMRM"),
			c("Parameter","S.E.")), byrow=TRUE)

# Copy these matrices and paste them into Excel
write.table(SSobs, "clipboard", sep="\\t", row.names=FALSE)
	# Paste into Excel

write.table(SSrand, "clipboard", sep="\\t", row.names=FALSE)
	# Paste into Excel




################################################################
############### 		Part 3 -- 	     	       ###############
###############  Do immigrants share sites with  ###############
###############  residents more than chance?     ###############
################################################################

#########
######### A) -- 
######### Test this idea by subsetting to different groups and measuring
######### the proportion of sites overlapped between the two groups,
######### and compare that to a random distribution.
#########

library(spatstat) # Access for spatial stats stuff
library(gdata)	# Access library for 'drop.levels()'
library(reshape)  # Acess for 'melt()' function later on

### Juveniles - resident females

res = matrix(0, nrow = 4, ncol = 65, 			# max surveys = 65 (plot 3)
	dimnames = list(c("One","Two","Three","Four"),  # 4 row for each plot
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("J","RF"))) 
							# Subset to immigr-res groups
overlap = subset(as.data.frame(table(droplevels(group$Locality2), 
	droplevels(group$ImmRes),droplevels(group$FrogNo))), Freq>0)
overlap2 = table(overlap$Var1, overlap$Var2)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap2) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap2[ which(overlap2[,1] > 0 & overlap2[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap2))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap2)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}
}
}	# Close the for loops

res 

res = melt(as.data.frame(res))
Group = rep("J-RF",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(jrf = cbind(Group,Plot,res))
colnames(jrf) = c("Group","Plot","Survey","PropShared")

### Juveniles - resident males

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("J","RM"))) 
							# Subset to immigr-res groups
overlap = subset(as.data.frame(table(droplevels(group$Locality2), 
	droplevels(group$ImmRes),droplevels(group$FrogNo))), Freq>0)
overlap2 = table(overlap$Var1, overlap$Var2)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap2) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap2[ which(overlap2[,1] > 0 & overlap2[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap2))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap2)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}							# Close the if-else loops

}}							# Close the for-loops		

res 

res = melt(as.data.frame(res))
Group = rep("J-RM",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(jrm = cbind(Group,Plot,res))
colnames(jrm) = c("Group","Plot","Survey","PropShared")
summary(jrm)


### Immigrant females - resident females

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("IF","RF"))) 
							# Subset to immigr-res groups
overlap = subset(as.data.frame(table(droplevels(group$Locality2), 
	droplevels(group$ImmRes),droplevels(group$FrogNo))), Freq>0)
overlap2 = table(overlap$Var1, overlap$Var2)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap2) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap2[ which(overlap2[,1] > 0 & overlap2[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap2))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap2)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}							# Close the if-else loops

}}							# Close the for-loops		

res 

res = melt(as.data.frame(res))
Group = rep("IF-RF",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(ifrf = cbind(Group,Plot,res))
colnames(ifrf) = c("Group","Plot","Survey","PropShared")
summary(ifrf)


### Immigrant female - resident males

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("IF","RM"))) 
							# Subset to immigr-res groups
overlap = subset(as.data.frame(table(droplevels(group$Locality2), 
	droplevels(group$ImmRes),droplevels(group$FrogNo))), Freq>0)
overlap2 = table(overlap$Var1, overlap$Var2)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap2) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap2[ which(overlap2[,1] > 0 & overlap2[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap2))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap2)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}							# Close the if-else loops

}}							# Close the for-loops		

res 

res = melt(as.data.frame(res))
Group = rep("IF-RM",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(ifrm = cbind(Group,Plot,res))
colnames(ifrm) = c("Group","Plot","Survey","PropShared")
summary(ifrm)


### Immigrant males - resident females

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("IM","RF"))) 
							# Subset to immigr-res groups
overlap = subset(as.data.frame(table(droplevels(group$Locality2), 
	droplevels(group$ImmRes),droplevels(group$FrogNo))), Freq>0)
overlap2 = table(overlap$Var1, overlap$Var2)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap2) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap2[ which(overlap2[,1] > 0 & overlap2[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap2))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap2)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}							# Close the if-else loops

}}							# Close the for-loops		

res 

res = melt(as.data.frame(res))
Group = rep("IM-RF",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(imrf = cbind(Group,Plot,res))
colnames(imrf) = c("Group","Plot","Survey","PropShared")
summary(imrf)


### Immigrant males - resident males

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("IM","RM"))) 
							# Subset to immigr-res groups
overlap = subset(as.data.frame(table(droplevels(group$Locality2), 
	droplevels(group$ImmRes),droplevels(group$FrogNo))), Freq>0)
overlap2 = table(overlap$Var1, overlap$Var2)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap2) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap2[ which(overlap2[,1] > 0 & overlap2[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap2))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap2)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}							# Close the if-else loops

}}							# Close the for-loops		

res 

res = melt(as.data.frame(res))
Group = rep("IM-RM",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(imrm = cbind(Group,Plot,res))
colnames(imrm) = c("Group","Plot","Survey","PropShared")
summary(imrm)



### Bind all the proportion of shared sites dataframes into a single file


ObsPropShared = rbind(jrf,jrm,ifrf,ifrm,imrf,imrm)

plot(PropShared ~ Group, data=na.omit(ObsPropShared))

res = lme(PropShared ~ Group, data=na.omit(ObsPropShared), random=(~1|Plot))
summary(res)


#########
######### B) -- 
######### Estimate the proportion of shared sites between groups  
######### that is expected by random chance
#########


### Juveniles - resident females

NoSims = 5	# Five different simulations for each grouping

res0 = matrix(NA, ncol = 65, nrow = 1, dimnames=list(c("na"),c(1:65)))

for (h in 1:NoSims){						# For each simulation

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("J","RF"))) 
							# Subset to immigr-res groups

RandSites = sample(c("A1","A12","A2","A23","A3","A34",	# Column A
			"AB1","AB12","AB2","AB23","AB3","AB34",	# Column AB
			"B1","B12","B2","B23","B3","B34",		# Column B
			"BC1","BC12","BC2","BC23","BC3","BC34",	# Column BC
			"C1","C12","C2","C23","C3","C34",		# Column C
			"CD1","CD12","CD2","CD23","CD3","CD34",	# Column CD
			"D1","D12","D2","D23","D3","D34",		# Column D
			"DE1","DE12","DE2","DE23","DE3","DE34"),	# Column DE
			nrow(group), replace=TRUE)			# ImmRes groups
										# Sample w/ replacement
RandSites = data.frame(cbind(group$ImmRes,RandSites))
colnames(RandSites) = c("Group","Site")

overlap = table(RandSites$Site, RandSites$Group)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap[ which(overlap[,1] > 0 & overlap[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}

}	# Close survey loop
}	# Close plot loop

res0 = rbind(res0,res)

}	# Close simulation replication loop

res0
res0 = res0[-c(1),]

res = melt(as.data.frame(res0))
Group = rep("J-RF",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(jrfRAND = cbind(Group,Plot,res))
colnames(jrfRAND) = c("Group","Plot","Survey","PropShared")



### Juveniles - resident males

res0 = matrix(NA, ncol = 65, nrow = 1, dimnames=list(c("na"),c(1:65)))

for (h in 1:NoSims){						# 5 different simulations

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("J","RM"))) 
							# Subset to immigr-res groups

RandSites = sample(c("A1","A12","A2","A23","A3","A34",	# Column A
			"AB1","AB12","AB2","AB23","AB3","AB34",	# Column AB
			"B1","B12","B2","B23","B3","B34",		# Column B
			"BC1","BC12","BC2","BC23","BC3","BC34",	# Column BC
			"C1","C12","C2","C23","C3","C34",		# Column C
			"CD1","CD12","CD2","CD23","CD3","CD34",	# Column CD
			"D1","D12","D2","D23","D3","D34",		# Column D
			"DE1","DE12","DE2","DE23","DE3","DE34"),	# Column DE
			nrow(group), replace=TRUE)			# ImmRes groups
										# Sample w/ replacement
RandSites = data.frame(cbind(group$ImmRes,RandSites))
colnames(RandSites) = c("Group","Site")

overlap = table(RandSites$Site, RandSites$Group)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap[ which(overlap[,1] > 0 & overlap[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}

}	# Close survey loop
}	# Close plot loop

res0 = rbind(res0,res)

}	# Close simulation replication loop

res0
res0 = res0[-c(1),]

res = melt(as.data.frame(res0))
Group = rep("J-RM",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(jrmRAND = cbind(Group,Plot,res))
colnames(jrmRAND) = c("Group","Plot","Survey","PropShared")



### Immigrant females - resident females

res0 = matrix(NA, ncol = 65, nrow = 1, dimnames=list(c("na"),c(1:65)))

for (h in 1:NoSims){	

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("IF","RF"))) 
							# Subset to immigr-res groups

RandSites = sample(c("A1","A12","A2","A23","A3","A34",	# Column A
			"AB1","AB12","AB2","AB23","AB3","AB34",	# Column AB
			"B1","B12","B2","B23","B3","B34",		# Column B
			"BC1","BC12","BC2","BC23","BC3","BC34",	# Column BC
			"C1","C12","C2","C23","C3","C34",		# Column C
			"CD1","CD12","CD2","CD23","CD3","CD34",	# Column CD
			"D1","D12","D2","D23","D3","D34",		# Column D
			"DE1","DE12","DE2","DE23","DE3","DE34"),	# Column DE
			nrow(group), replace=TRUE)			# ImmRes groups
										# Sample w/ replacement
RandSites = data.frame(cbind(group$ImmRes,RandSites))
colnames(RandSites) = c("Group","Site")

overlap = table(RandSites$Site, RandSites$Group)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap[ which(overlap[,1] > 0 & overlap[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}

}	# Close survey loop
}	# Close plot loop

res0 = rbind(res0,res)

}	# Close simulation replication loop

res0
res0 = res0[-c(1),]

res = melt(as.data.frame(res0))
Group = rep("IF-RF",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(ifrfRAND = cbind(Group,Plot,res))
colnames(ifrfRAND) = c("Group","Plot","Survey","PropShared")


### Immigrant females - resident males

res0 = matrix(NA, ncol = 65, nrow = 1, dimnames=list(c("na"),c(1:65)))

for (h in 1:NoSims){	

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("IF","RM"))) 
							# Subset to immigr-res groups

RandSites = sample(c("A1","A12","A2","A23","A3","A34",	# Column A
			"AB1","AB12","AB2","AB23","AB3","AB34",	# Column AB
			"B1","B12","B2","B23","B3","B34",		# Column B
			"BC1","BC12","BC2","BC23","BC3","BC34",	# Column BC
			"C1","C12","C2","C23","C3","C34",		# Column C
			"CD1","CD12","CD2","CD23","CD3","CD34",	# Column CD
			"D1","D12","D2","D23","D3","D34",		# Column D
			"DE1","DE12","DE2","DE23","DE3","DE34"),	# Column DE
			nrow(group), replace=TRUE)			# ImmRes groups
										# Sample w/ replacement
RandSites = data.frame(cbind(group$ImmRes,RandSites))
colnames(RandSites) = c("Group","Site")

overlap = table(RandSites$Site, RandSites$Group)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap[ which(overlap[,1] > 0 & overlap[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}

}	# Close survey loop
}	# Close plot loop

res0 = rbind(res0,res)

}	# Close simulation replication loop

res0
res0 = res0[-c(1),]

res = melt(as.data.frame(res0))
Group = rep("IF-RM",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(ifrmRAND = cbind(Group,Plot,res))
colnames(ifrmRAND) = c("Group","Plot","Survey","PropShared")




### Immigrant males - resident females

res0 = matrix(NA, ncol = 65, nrow = 1, dimnames=list(c("na"),c(1:65)))

for (h in 1:NoSims){	

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("IM","RF"))) 
							# Subset to immigr-res groups

RandSites = sample(c("A1","A12","A2","A23","A3","A34",	# Column A
			"AB1","AB12","AB2","AB23","AB3","AB34",	# Column AB
			"B1","B12","B2","B23","B3","B34",		# Column B
			"BC1","BC12","BC2","BC23","BC3","BC34",	# Column BC
			"C1","C12","C2","C23","C3","C34",		# Column C
			"CD1","CD12","CD2","CD23","CD3","CD34",	# Column CD
			"D1","D12","D2","D23","D3","D34",		# Column D
			"DE1","DE12","DE2","DE23","DE3","DE34"),	# Column DE
			nrow(group), replace=TRUE)			# ImmRes groups
										# Sample w/ replacement
RandSites = data.frame(cbind(group$ImmRes,RandSites))
colnames(RandSites) = c("Group","Site")

overlap = table(RandSites$Site, RandSites$Group)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap[ which(overlap[,1] > 0 & overlap[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}

}	# Close survey loop
}	# Close plot loop

res0 = rbind(res0,res)

}	# Close simulation replication loop

res0
res0 = res0[-c(1),]

res = melt(as.data.frame(res0))
Group = rep("IM-RF",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(imrfRAND = cbind(Group,Plot,res))
colnames(imrfRAND) = c("Group","Plot","Survey","PropShared")



### Immigrant males - resident males

res0 = matrix(NA, ncol = 65, nrow = 1, dimnames=list(c("na"),c(1:65)))

for (h in 1:NoSims){	

res = matrix(0, nrow = 4, ncol = 65, 			# 65 surveys in plot 3
	dimnames = list(c("One","Two","Three","Four"),
	c(1:65)))

for (i in 1:4){					# For each plot, 1:4
for (j in 1:65){					# For each survey 1:65	

plot = subset(datum, Plot == i)		# Subset to each plot

plot = drop.levels(plot)			# Drop the extra empty factors
Date4 = as.numeric(droplevels(factor(plot$Date2))) # Generate date as 1:65
plot = cbind(plot,Date4)			# Bind new date variable
	
survey = subset(plot, Date4 == j)		# Subset to each survey
group = droplevels(subset(survey, ImmRes %in% c("IM","RM"))) 
							# Subset to immigr-res groups

RandSites = sample(c("A1","A12","A2","A23","A3","A34",	# Column A
			"AB1","AB12","AB2","AB23","AB3","AB34",	# Column AB
			"B1","B12","B2","B23","B3","B34",		# Column B
			"BC1","BC12","BC2","BC23","BC3","BC34",	# Column BC
			"C1","C12","C2","C23","C3","C34",		# Column C
			"CD1","CD12","CD2","CD23","CD3","CD34",	# Column CD
			"D1","D12","D2","D23","D3","D34",		# Column D
			"DE1","DE12","DE2","DE23","DE3","DE34"),	# Column DE
			nrow(group), replace=TRUE)			# ImmRes groups
										# Sample w/ replacement
RandSites = data.frame(cbind(group$ImmRes,RandSites))
colnames(RandSites) = c("Group","Site")

overlap = table(RandSites$Site, RandSites$Group)

if (j > length(table(as.factor(plot$Date4)))){	# If-else statements to fill matrix
	res[i,j] = NA					# Didn't take this many samples
} else { if (length(table(group$ImmRes)) < 1){
	res[i,j] = NA					# No groups present
} else { if (ncol(overlap) < 2){			
	res[i,j] = 0					# Only one group present
} else 
shared = overlap[ which(overlap[,1] > 0 & overlap[,2] > 0),]
{ if (is.vector(shared)){
	res[i,j] = (1/nrow(overlap))			# Single shared site
} else { if (ncol(shared)){			
	res[i,j] = nrow(shared)/nrow(overlap)	# Multiple shared sites
} else {
	res[i,j] = 0					# No shared sites
}}}}}

}	# Close survey loop
}	# Close plot loop

res0 = rbind(res0,res)

}	# Close simulation replication loop

res0
res0 = res0[-c(1),]

res = melt(as.data.frame(res0))
Group = rep("IM-RM",length(res[,1]))
Plot = rep(c(1:4),length(res[,1])/4)
(imrmRAND = cbind(Group,Plot,res))
colnames(imrmRAND) = c("Group","Plot","Survey","PropShared")


###### C ---
###### Bind the random simulated shared sites into a single file, and then
###### Bind observed and random results into another file for analysis

RandPropShared = rbind(jrfRAND,jrmRAND,ifrfRAND,ifrmRAND,imrfRAND,imrmRAND)

SharedSites = rbind(ObsPropShared,RandPropShared)


## Add dummy-coded variables for each group

library(dummies)

x = dummy(SharedSites$Group)
ObsRand = c(rep(c("Obs"), times=nrow(ObsPropShared)),
	rep(c("Rand"), times=nrow(RandPropShared)))
x = cbind(x,ObsRand)
colnames(x) = c("JRF","JRM","IFRF","IFRM","IMRF","IMRM","ObsRand")
SharedSites = cbind(SharedSites,x)
SharedSites = na.omit(SharedSites)
head(SharedSites)

detach(package:dummies)


##### Analyses

library(nlme)

## Pairwise differences between observed and random proportion of shared sites

resJRF = lme(PropShared ~ ObsRand, random=(~1|Plot/Survey),
	 data=subset(SharedSites, JRF == 1))
summary(resJRF)

resJRM = lme(PropShared ~ ObsRand, random=(~1|Plot/Survey),
	 data=subset(SharedSites, JRM == 1))
summary(resJRM)

resIFRF = lme(PropShared ~ ObsRand, random=(~1|Plot/Survey),
	 data=subset(SharedSites, IFRF == 1))
summary(resIFRF)


resIFRM = lme(PropShared ~ ObsRand, random=(~1|Plot/Survey),
	 data=subset(SharedSites, IFRM == 1))
summary(resIFRM)

resIMRF = lme(PropShared ~ ObsRand, random=(~1|Plot/Survey),
	 data=subset(SharedSites, IMRF == 1))
summary(resIMRF)

resIMRM = lme(PropShared ~ ObsRand, random=(~1|Plot/Survey),
	 data=subset(SharedSites, IMRM == 1))
summary(resIMRM)


# Create matrices summarizing the results of the analysis
# to produce a neat figure in Excel

SSobs = matrix(nrow=6,ncol=2,c(
	resJRF$coefficients$fixed[1],sqrt(summary(resJRF)$varFix)[1,1],
	resJRM$coefficients$fixed[1],sqrt(summary(resJRM)$varFix)[1,1],
	resIFRF$coefficients$fixed[1],sqrt(summary(resIFRF)$varFix)[1,1],	
	resIFRM$coefficients$fixed[1],sqrt(summary(resIFRM)$varFix)[1,1],
	resIMRF$coefficients$fixed[1],sqrt(summary(resIMRF)$varFix)[1,1],
	resIMRM$coefficients$fixed[1],sqrt(summary(resIMRM)$varFix)[1,1]),
	dimnames=list(c("JRF","JRM","IFRF","IFRM","IMRF","IMRM"),
			c("Parameter","S.E.")), byrow=TRUE)

SSrand = matrix(nrow=6,ncol=2,c(
	resJRF$coefficients$fixed[1]+resJRF$coefficients$fixed[2],sqrt(summary(resJRF)$varFix)[2,2],
	resJRM$coefficients$fixed[1]+resJRM$coefficients$fixed[2],sqrt(summary(resJRM)$varFix)[2,2],
	resIFRF$coefficients$fixed[1]+resIFRF$coefficients$fixed[2],sqrt(summary(resIFRF)$varFix)[2,2],	
	resIFRM$coefficients$fixed[1]+resIFRM$coefficients$fixed[2],sqrt(summary(resIFRM)$varFix)[2,2],
	resIMRF$coefficients$fixed[1]+resIMRF$coefficients$fixed[2],sqrt(summary(resIMRF)$varFix)[2,2],
	resIMRM$coefficients$fixed[1]+resIMRM$coefficients$fixed[2],sqrt(summary(resIMRM)$varFix)[2,2]),
	dimnames=list(c("JRF","JRM","IFRF","IFRM","IMRF","IMRM"),
			c("Parameter","S.E.")), byrow=TRUE)

# Copy these matrices and paste them into Excel
write.table(SSobs, "clipboard", sep="\\t", row.names=FALSE)
	# Paste into Excel

write.table(SSrand, "clipboard", sep="\\t", row.names=FALSE)
	# Paste into Excel



## Overall model of differences in proportion of shared sites

res = lme(PropShared ~ Group, data=na.omit(ObsPropShared), random=(~1|Plot/Survey))
summary(res)
anova(res)

# dummy-coded model to test for pairwise differences
res2 = lme(PropShared ~ JRF + JRM + IFRF + IMRF + IMRM,
	data=subset(na.omit(SharedSites), ObsRand == "Obs"),
	random=(~1|Plot/Survey))
summary(res2)






##########################################################
######## Part 4 -- 					##########
######## Do immigrants settle in space closer	##########
######## to residents than expected by chance?  ##########
##########################################################

# I first explored using classic nearest-neighbor distances (NND) between
# different immigrant-resident groups. However, a limitation of this approach
# is that NND are confounded by the geometry of study window;
# e.g., true nearest neighbors may be off the plot!

# An alternative to this is to estimate the nearest-neighbor
# distribution function, G(r). This approach corrects for border issues
# and irregular geometry. It estimates NND across a continuous radii
# around individuals.

# Here, let's use the G(r) function to test how frequently
# the G(r) relationship is consistent with CLUSTERING.


#########
######### A) --
######### Use a for-loop to generate G(r) functions for 
######### each immigrant-resident group in each plot in each season.
#########

library(spatstat) 
library(gdata)
library(reshape)

results = matrix(0, nrow = 4, ncol = 5, 			
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))


ImmResGroups = list(c("J","RF"),c("J","RM"),
				c("IF","RF"),c("IF","RM"),
				c("IM","RF"),c("IM","RM"))

# Initiate for loops:

for (h in ImmResGroups){			# For each immigrant resident group
					
for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5	

plot = subset(datum, Plot == i)		# Subset to a plot

tab = table(plot$FrogNo)			# Tabulate no. of indiv. obs.
noRecaps = names(tab)[tab > 2]		# Find ind's seen >2 times
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])   # Subset to the above inds 
plot = drop.levels(plot)			# Drop the extra empty factors

group = droplevels(subset(plot, ImmRes %in% h))
							# Subset to imm-res group
season = subset(group, Season == j)		# Subset to a season

means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes, data=season, mean))
grid = ppp(means$X,means$Y,c(0,12),c(0,9))
marks(grid) = as.factor(means$ImmRes)	
	# Take average location of individuals and categorize them by groups

g = mad.test(grid, Gest, correction="km", alternative=c("greater"))
						# Nearest-neighbor distance function

results[i,j] = g$statistic$mad	# Store G(r) MAD statistic in dataset

	}
}				

res = melt(results)
group = data.frame(c(rep(paste0(h[1],h[2]), length(res[,1]))))
ObsRand = data.frame(c(rep("Observed", length(res[,1]))))
res = cbind(group, res, ObsRand)
colnames(res) = c("Group","Plot","Season","Gfxn","Random")

assign(paste0(h[1],h[2],"res"), res) 	
	# Creates unique name for each immigrant-resident group

	}	# End the immigrant-resident subset loop

# This loop created a series of new objects (e.g., JRMres, IMRMres, etc.)
# which describe the proportion of distances in nearest-neighbor distance
# analysis which were consistent with clustering.  
	# Larger MAD values indicate greater degrees of clustering
	# relative to smaller ones

nndObs = rbind(JRFres, JRMres, IFRFres, IFRMres, IMRFres, IMRMres)
	# Bind these together to facilitate comparison to random later on



detach(package:spatstat) 
detach(package:gdata)
detach(package:reshape)



######### 
######### B) -- 
######### Use a for-loop to calculate the same G-statistic function
######### estimating nearest-neighbor distances expected by RANDOM
######### 

library(spatstat) 
library(gdata)
library(reshape)

NoSims = 10	# Five simulations to better approximate true random means

rand0 = as.data.frame(matrix(NA, ncol = 5, nrow = 1,
	dimnames=list(c("na"),c("Group","Plot","Season","Gfxn","Random"))))
			# Dataframe to save simulated datasets

random = matrix(0, nrow = 4, ncol = 5,
	dimnames = list(c("One","Two","Three","Four"),
	c("A","B","C","D","E")))
			# Data matrix to store initial simulated estimates

ImmResGroups = list(c("J","RF"),c("J","RM"),c("IF","RF"),
			c("IF","RM"),c("IM","RF"),c("IM","RM"))
			# Immigr-resident groups for subsettting

for (h in ImmResGroups){			# For each immigrant resident group
for (g in 1:NoSims){				# Five simulation replications		
for (i in 1:4){					# For each plot, 1:4
for (j in 1:5){					# For each season 1:5						

plot = subset(datum, Plot == i)	

tab = table(plot$FrogNo)			
noRecaps = names(tab)[tab > 2]		
plot = data.frame(
	plot[plot$FrogNo %in% noRecaps,])    
plot = drop.levels(plot)	

group = droplevels(subset(plot, ImmRes %in% h))
		
season = subset(group, Season == j)		

means = data.frame(aggregate(cbind(X,Y) ~ FrogNo + ImmRes,
	 data=season, mean))	

rand = rpoint(length(means[,1]), win=owin(c(0,12),c(0,9)))
grid = ppp(rand$x,rand$y,c(0,12),c(0,9))
marks(grid) = as.factor(means$ImmRes)

g = mad.test(grid, Gest, alternative=c("less"))
						# Nearest-neighbor distance function

random[i,j] = g$statistic$mad		# Store G(r) p-value in dataset	

	}	# Close the season loop
}		# Close the plot loops

rand = melt(random)
group = data.frame(c(rep(paste0(h[1],h[2]), length(random[,1]))))
ObsRand = data.frame(c(rep("Random", length(random[,1]))))
rand = cbind(group, rand, ObsRand)
colnames(rand) = c("Group","Plot","Season","Gfxn","Random")

rand0 = rbind(rand0,rand)

}		# Close the simulation loop

rand = rand0[-c(1),]

assign(paste0(h[1],h[2],"rand"), rand) 	
	# Creates unique name for each immigrant-resident group

	}	# End the immigrant-resident subset loop


# Bind all the random simulated results together into a single file
nndRand = rbind(JRFrand,JRMrand,IFRFrand,IFRMrand,IMRFrand,IMRMrand)



detach(package:spatstat) 
detach(package:gdata)
detach(package:reshape)



#########
######### C) -- 
######### Compare the number of NND functions observed that were clustered
######### relative to those generated by RANDOM.
#########

### Combine the observed and random datasets here
analysis = rbind(nndObs,nndRand)


### Paired tests between observed and random values, using 'analysis' object
### for six hypothesized conspecific attraction groupings:

library(nlme)

# J-RF
testJRF = lme(Gfxn ~ Random, random=(~1|Plot/Season),
	data=subset(analysis, Group == "JRF"))
summary(testJRF)

# J-RM
testJRM = lme(Gfxn ~ Random, random=(~1|Plot/Season),
	data=subset(analysis, Group == "JRM"))
summary(testJRM)

# IF-RF
testIFRF = lme(Gfxn ~ Random, random=(~1|Plot/Season),
	data=subset(analysis, Group == "IFRF"))
summary(testIFRF)


# IF-RM
testIFRM = lme(Gfxn ~ Random, random=(~1|Plot/Season),
	data=subset(analysis, Group == "IFRM"))
summary(testIFRM)

# IM-RF
testIMRF = lme(Gfxn ~ Random, random=(~1|Plot/Season),
	data=subset(analysis, Group == "IMRF"))
summary(testIMRF)

# IM-RM
testIMRM = lme(Gfxn ~ Random, random=(~1|Plot/Season),
	data=subset(analysis, Group == "IMRM"))
summary(testIMRM)


### Extract effect sizes, save them into matrices.
### Then copy these into Excel for to generate a figure

NNDobs = matrix(nrow=6,ncol=2,c(
	testJRF$coefficients$fixed[1],sqrt(summary(testJRF)$varFix)[1,1],
	testJRM$coefficients$fixed[1],sqrt(summary(testJRM)$varFix)[1,1],
	testIFRF$coefficients$fixed[1],sqrt(summary(testIFRF)$varFix)[1,1],	
	testIFRM$coefficients$fixed[1],sqrt(summary(testIFRM)$varFix)[1,1],
	testIMRF$coefficients$fixed[1],sqrt(summary(testIMRF)$varFix)[1,1],
	testIMRM$coefficients$fixed[1],sqrt(summary(testIMRM)$varFix)[1,1]),
	dimnames=list(c("JRF","JRM","IFRF","IFRM","IMRF","IMRM"),
			c("Parameter","S.E.")), byrow=TRUE)

NNDrand = matrix(nrow=6,ncol=2,c(
	testJRF$coefficients$fixed[1]+testJRF$coefficients$fixed[2],sqrt(summary(testJRF)$varFix)[2,2],
	testJRM$coefficients$fixed[1]+testJRM$coefficients$fixed[2],sqrt(summary(testJRM)$varFix)[2,2],
	testIFRF$coefficients$fixed[1]+testIFRF$coefficients$fixed[2],sqrt(summary(testIFRF)$varFix)[2,2],	
	testIFRM$coefficients$fixed[1]+testIFRM$coefficients$fixed[2],sqrt(summary(testIFRM)$varFix)[2,2],
	testIMRF$coefficients$fixed[1]+testIMRF$coefficients$fixed[2],sqrt(summary(testIMRF)$varFix)[2,2],
	testIMRM$coefficients$fixed[1]+testIMRM$coefficients$fixed[2],sqrt(summary(testIMRM)$varFix)[2,2]),
	dimnames=list(c("JRF","JRM","IFRF","IFRM","IMRF","IMRM"),
			c("Parameter","S.E.")), byrow=TRUE)

write.table(NNDobs, "clipboard", sep="\\t", row.names=FALSE)
	# Paste into 'Figure 4' tab, Observed data

write.table(NNDrand, "clipboard", sep="\\t", row.names=FALSE)
	# Paste into 'Figure 4' tab, RANDOM data


#############################################################################
##### The End ###############################################################
#############################################################################
