###California Tiger Salamander (Ambystoma californiense) integral projection model and population viability analysis
###Arianne Messerman
###October 2021

#setwd("")

#Load source code for functions
source("CTS_Source.R")

##Load .csv of 500 samples from  posterior districbutions of survival, growth, maturity and fertility models
surv<-read.csv("survival-posterior-samples-CENTERED.csv", header = T)
grow<-read.csv("growth-posterior-samples-CENTERED.csv", header = T)
mat<-read.csv("maturity-posterior-samples-CENTERED.csv", header = T)
fert<-read.csv("fertility-posterior-samples-CENTERED.csv", header = T)

## Number of draws to take from the joint posterior distribution of the parameters. 
## Cannot be greater than the number of draws provided in the .csv file, which is 500.
Ndraws<-min(500,nrow(grow))

####################################################################################################
## Set up IPM parameter vectors
## 'cts' is a matrix where rows are vital rate coefficients and columns are posterior draws
## below, we will loop over columns, sending each set of coefficients into the stochastic IPM
####################################################################################################
cts<-matrix(NA,nrow=22,ncol=Ndraws)

## Params 01-03: Metamorph growth
cts[1,]<-grow$meta.inter1      	## growth intercept
cts[2,]<-grow$meta.slope1				## growth linear slope
cts[3,]<-grow$meta.sigma1				## SD of metamorph growth function

## Params 04-06: Adult growth
cts[4,]<-grow$ad.inter1      	  ## growth intercept
cts[5,]<-grow$ad.slope1				  ## growth linear slope
cts[6,]<-grow$ad.sigma1					## SD of adult growth function

## Params 07-11: Survival params
cts[7,]<-surv$mean.meta.eta.phi ## mean metamorph survival intercept for years 1,2,6, and 7
cts[8,]<-surv$mass.alpha				## survival on mass slope
cts[9,]<-surv$ad.eta.phi1			  ## adult survival intercept
cts[10,]<-surv$fidel.m	        ## metamorph site fidelity
cts[11,]<-surv$fidel.a	        ## adult site fidelity

## Params 12-13: Maturity params
cts[12,]<-mat$inter   		      ## maturity intercept
cts[13,]<-mat$slope			        ## maturity slope

## Params 14-15: Fertility params (clutch size)
cts[14,]<-fert$inter1    	      ## fertility intercept 
cts[15,]<-fert$slope1		        ## fertility linear slope

## Params 16-27: Misc params (bounds of continuous size domain in units of mean-centered ln(g))
cts[16,]<- -0.85                ## min size overall
cts[17,]<- 5.2                  ## max size overall
cts[18,]<- 2.24                 ## mean size of new metamorphs at low density
cts[19,]<- 0.345                ## SD size of new metamorphs at low density
cts[20,]<- 0.50                 ## Proportion of females in breeding population (Smith and Voss 2009)
cts[21,]<- 0.341                ## Proportion of mature females that return to breed in a given year (Trenham and Shaffer 2005)
cts[22,]<- 0.092                ## Maximum larval survival probability at low density (Trenham and Shaffer 2005)


###################################################################################################
## The full density-independent IPM kernel and meshpoints in the units of size
## params get passed to the vital rate functions in the SOURCE file
## matsize is the dimension of the approximating matrix
###################################################################################################
##Size of the approximating matrix to the IPM kernel
matsize<-122
lower<-cts[16]
upper<-cts[17]

n<-matsize
L<-lower
U<-upper
h<-(U-L)/n          #Bin size =0.05
b<-L+c(0:n)*h       #Lower boundaries of bins 
y<-0.5*(b[1:n]+b[2:(n+1)])  #Bin midpoints

#Save vital rate function outputs
g.m<-g.a<-s.m<-s.a<-p.m<-p.a<-ma<-fe<-f<-array(0, dim=c(1, length(y), Ndraws))
rec<-array(0, dim=c(1,length(y)))
for(i in 1:Ndraws){
 g.m[,,i]<-gxy.m(y,y,cts[,i])
 g.a[,,i]<-gxy.a(y,y,cts[,i])
 s.m[,,i]<-sx.m(y,cts[,i])
 s.a[,,i]<-sx.a(y,cts[,i])
 p.m[,,i]<-pxy.m(y,y,cts[,i])
 p.a[,,i]<-pxy.a(y,y,cts[,i])
 ma[,,i]<-mx(y,cts[,i])
 fe[,,i]<-fer(y,cts[,i])
 f[,,i]<-fx(y,cts[,i])
}

rec<-recruits(y,cts[,1])

f<-ifelse(f<0,0,f)#Fertility estimates <0 equal 0
for(i in 1:Ndraws){
  fe[,,i]<-ifelse(fe[,,i]<0,0,fe[,,i])#Clutch size estimates <0 equal 0
}

#Examine vital rate function outputs
plot(y,g.m[,,1], type='l', ylim=c(0,2), xlim=c(0.69, 4.2), ylab="Probability density", xlab="Ln(body mass (g))", bty="l", 
     cex.axis=1.4, cex.lab=1.4, col="skyblue3")
mtext("a", 3, .1, at=-0.8, cex=1.2)
for(i in 2:length(g.m[,,1])){
lines(y,g.m[,,i], col="Skyblue3")
}
for(i in 1:length(g.a[,,1])){
  lines(y,g.a[,,i], col="plum1")
}
g.m.med<-g.a.med<-c(1:length(y))
for (i in 1:length(y)){
  g.m.med[i]<-median(g.m[,i,])
  g.a.med[i]<-median(g.a[,i,])
}
lines(y, g.m.med, col="midnightblue", lwd=3)
lines(y, g.a.med, lty=2, col="plum4", lwd=3)

#pdf("Figure3.pdf", width = 12, height = 11)
#par(mfrow=c(3,3))
plot(y,(cts[1,1] + cts[2,1]*(y-2.39)), type='l', ylab="Body mass (g) in year t+1", xlim=c(0.69, 4.2), ylim=c(0.69, 4.2),
     xlab="Body mass (g) in year t", bty="l", cex.axis=1.4, cex.lab=1.4, col="skyblue3", xaxt = "n", yaxt="n")
axis(1, at=c(-1:5), labels=c("0.4", "1.0", "2.7", "7.4", "20.1", "54.6", "148.4"), cex.axis=1.4)
axis(2, at=c(-1:5), labels=c("0.4", "1.0", "2.7", "7.4", "20.1", "54.6", "148.4"), cex.axis=1.4)
mtext("a", 3, .1, at=.7, cex=1.2)
for(i in 2:length(g.m[,,1])){
  lines(y,(cts[1,i] + cts[2,i]*(y-2.39)), col="Skyblue3")
}
for(i in 1:length(g.a[,,1])){
  lines(y,(cts[4,i] + cts[5,i]*(y-2.95)), col="plum1")
}
g.m.med<-g.a.med<-c(1:length(y))
for (i in 1:length(y)){
  g.m.med[i]<-median((cts[1,] + cts[2,]*(y[i]-2.39)))
  g.a.med[i]<-median((cts[4,] + cts[5,]*(y[i]-2.95)))
}
lines(y, g.m.med, col="midnightblue", lwd=3)
lines(y, g.a.med, lty=2, col="plum4", lwd=3)
#points(M$ln.mass, rep(0,length(M$ln.mass)), col="grey",type='l', lwd=5)
#points(A$ln.mass, rep(0,length(A$ln.mass)), col="pink",type='l', lwd=5)
plot(y,s.m[,,1], type='l', ylim=c(0,1), ylab="True survival probability", xlab="Body mass (g)", bty="l", 
     cex.axis=1.4, cex.lab=1.4, col="skyblue3", xaxt = "n", xlim=c(0.69, 4.2))
axis(1, at=c(-1:5), labels=c("0.4", "1.0", "2.7", "7.4", "20.1", "54.6", "148.4"), cex.axis=1.4)
mtext("b", 3, .1, at=.7, cex=1.2)
for(i in 2:length(s.m[,,1])){
  lines(y,s.m[,,i], col="skyblue3")
}
for(i in 1:length(s.a[,,1])){
  lines(y,s.a[,,i], col="plum1")
}
s.m.med<-s.a.med<-c(1:length(y))
for (i in 1:length(y)){
  s.m.med[i]<-median(s.m[,i,])
  s.a.med[i]<-median(s.a[,i,])
}
lines(y, s.m.med, col="midnightblue", lwd=3)
lines(y, s.a.med, col="plum4", lwd=3, lty=2)

#plot(y,p.m[,,1], type='l', ylim=c(0,1.5), ylab="Survival*growth", xlab="Ln(body mass (g))")
#for(i in 2:length(p.m[,,1])){
#  lines(y,p.m[,,i], col=1)
#}
#for(i in 1:length(p.a[,,1])){
#  lines(y,p.a[,,i], col=2)
#}

plot(y,ma[,,1], type='l', ylim=c(0,1), ylab="Maturation probability", xlab="Body mass (g)", col="grey70", bty="l", 
cex.axis=1.4, cex.lab=1.4, xaxt = "n", xlim=c(0.69, 4.2))
axis(1, at=c(-1:5), labels=c("0.4", "1.0", "2.7", "7.4", "20.1", "54.6", "148.4"), cex.axis=1.4)
mtext("c", 3, .1, at=.7, cex=1.2)
for(i in 2:length(ma[,,1])){
  lines(y,ma[,,i], col="grey70")
}
ma.med<-c(1:length(y))
for (i in 1:length(y)){
  ma.med[i]<-median(ma[,i,])
}
lines(y, ma.med, col=1, lwd=1)

plot(y,fe[,,1], type='l', ylab="Clutch size", xlab="Body mass (g)", col="grey70", bty="l", 
     cex.axis=1.4, cex.lab=1.4, xaxt = "n", xlim=c(0.69, 4.2), ylim=c(0, 2000))
axis(1, at=c(-1:5), labels=c("0.4", "1.0", "2.7", "7.4", "20.1", "54.6", "148.4"), cex.axis=1.4)
mtext("d", 3, .1, at=.7, cex=1.2)
for(i in 2:length(fe[,,1])){
  lines(y,fe[,,i], col="grey70")
}
fe.med<-c(1:length(y))
for (i in 1:length(y)){
  fe.med[i]<-median(fe[,i,])
}
lines(y, fe.med, col=1, lwd=3)

plot(y,f[,,1], type='l', ylim=c(0,35), ylab="Fecundity (metamorphs produced)", xlab="Body mass (g)", col="grey70", bty="l", 
     cex.axis=1.4, cex.lab=1.4,xaxt = "n", xlim=c(0.69, 4.2))
axis(1, at=c(-1:5), labels=c("0.4", "1.0", "2.7", "7.4", "20.1", "54.6", "148.4"), cex.axis=1.4)
mtext("e", 3, .1, at=.7, cex=1.2)
for(i in 2:length(f[,,1])){
  lines(y,f[,,i], col="grey70")
}
f.med<-c(1:length(y))
for (i in 1:length(y)){
  f.med[i]<-median(f[,i,])
}
lines(y, f.med, col=1, lwd=3)

#plot(y,rec, type='l', ylab="Probability density", xlab="Ln(body mass (g))", bty="l", 
#     cex.axis=1.4, cex.lab=1.4, col=1)
#mtext("f", 3, .1, at=-0.8, cex=1.2)


#Create empty matrices that will be the quadrants of the overall transition matrix
Tmat<-array(NA,dim=c(2*length(y), 2*length(y), Ndraws))
Mpmat<-Apmat<-Mmat<-Amat<-Mrmat<-Armat<-Mfmat<-Afmat<-array(0,dim=c(length(y), length(y), Ndraws))

rec.prob<-rec/sum(rec)#recruitment probability at size y
rec.prob<-as.matrix(rec.prob)#one column, 122 rows

#Lower right and left quadrants of the transition matrix
for(i in 1:Ndraws){
  for(s in 1:n){
    #Metamorph to one-year-old growth and survival (lower left)
    Mmat[,,i]<-t(outer(y,y,pxy.m,params=cts[,i]))*h
    #Juvenile/adult growth and survival (lower right)
    Amat[,,i]<-t(outer(y,y,pxy.a,params=cts[,i]))*h
  }
}

#Intermediate step
for(i in 1:Ndraws){
  for(s in 1:n){
    #Metamorph to one-year-old growth and survival (without bin size correction)
    Mpmat[,,i]<-t(outer(y,y,pxy.m,params=cts[,i]))
    #Juvenile/adult growth and survival (without bin size correction)
    Apmat[,,i]<-t(outer(y,y,pxy.a,params=cts[,i]))
  }
}

#Intermediate step
for(i in 1:Ndraws){
  for(s in 1:n){
    #Metamorphs produced per year by one-year-old females of a size class (without bin size correction)
    Mfmat[s,,i]<-Mpmat[s,,i]*f[,s,i]
    #Metamorphs produced per year by adult females of a size class (without bin size correction)
    Afmat[s,,i]<-Apmat[s,,i]*f[,s,i]
  }
}

#Upper right and left quadrants of the transition matrix
for(i in 1:Ndraws){
  for(s in 1:n){
    #Assign distribution of metamorphs produced of a size class per year by one-year-olds of a size class (upper left)
    Mrmat[,s,i]<-(sum(Mfmat[,s,i]))*rec.prob*h
    #Assign distribution of metamorphs produced of a size class per year by adults of a size class (upper right)
    Armat[,s,i]<-(sum(Afmat[,s,i]))*rec.prob*h
  }
}

#Create the transition matrix
for(i in 1:Ndraws){
  for(s in 1:n){
    for(j in 1:n){
      Tmat[j,s,i]<-Mrmat[j,s,i]
      Tmat[j,s+122,i]<-Armat[j,s,i]
      Tmat[j+122,s,i]<-Mmat[j,s,i]
      Tmat[j+122,s+122,i]<-Amat[j,s,i]
    }#j dimension 2
  }#s dimension 1
}#i iteration of posterior draws


####################################################################################################
## Simulate population growth at low densities over 100 years
## Calculate deterministic growth rate at low densities
####################################################################################################

#Calculate deterministic growth rate
iter<-100  # how many years to sample

#Store lambda values (principle eigenvectors)
eig<-vector("numeric",length=Ndraws)

#Loop through for each of Ndraws transition matrices to calculate deterministic lambda distribution
for(i in 1:Ndraws) {
eig[i]<- Re(eigen(Tmat[,,i])$value[1])#Re() returns real part of value (i.e., non-imaginary)
}
lambdas<-if(all(Im(lambdas<-zapsmall(eig))==0)) as.numeric(lambdas) else eig #Removes"+0i"
hist(lambdas)
mean(lambdas)#1.48
median(lambdas)#1.47
sum(ifelse(lambdas <1,1,0))# 0/500 lambdas indicate decline

#Find 95% confidence intervals (12th and 488th values) around the population growth rate (lambda) when the population is at low density
lamb.s<-sort(lambdas)
left <- lamb.s[12]#1.22
right <- lamb.s[488]#1.73

## Simulate population growth over 100 years
N<-array(NA,dim=c(iter,2*length(y), Ndraws)) #Store population size distribution for each year of simulation
N[1,,]<-1 #Simulation begins with one individual in each size class

for(i in 1:Ndraws){
  for(t in 1:(iter-1)){
    N[t+1,,i]<- as.matrix(Tmat[,,i]) %*% N[t,,i]#Apply transition matrix to current population size distribution
  }#t
}#i

Nt<-matrix(NA, iter, Ndraws)
lambt<-rt<-matrix(NA, iter-1, Ndraws)
for(i in 1:Ndraws){
  for(t in 1:iter){
    Nt[t,i]<-sum(N[t,,i]) #Total population size per year
  }
}

hist(log10(Nt[100,]))

par(mfrow=c(1,1))
plot(c(1:100), log10(Nt[1:100,1]), ylim=c(0, 30), type='l', ylab="Log(population size)", xlab="Year")
for (i in 2:Ndraws){
  lines(c(1:100), log10(Nt[1:100,i]), col=i)
}

####################################################################################################
## Conduct a sensitivity and elasticity analysis of all matrix elements
## Code adapted from Merow et al.: 
## https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.12146&file=mee312146-sup-0001-AppendixS1.pdf
####################################################################################################
w.eigen<-stable.dist<-v.eigen<-repro.val<-array(NA,dim=c(2*length(y), Ndraws))
v.dot.w<-vector("numeric",length=Ndraws)
sens<-elas<-array(NA,dim=c(2*length(y), 2*length(y), Ndraws))

for(i in 1:Ndraws) {
  w.eigen[,i] <- Re(eigen(Tmat[,,i])$vectors[,1])
  stable.dist[,i] <- w.eigen[,i]/sum(w.eigen[,i])
  v.eigen[,i] <- Re(eigen(t(Tmat[,,i]))$vectors[,1])
  repro.val[,i] <- v.eigen[,i]/v.eigen[1,i]
  v.dot.w[i] <- sum(stable.dist[,i]*repro.val[,i])*h
  sens[,,i] <- outer(repro.val[,i],stable.dist[,i])/v.dot.w[i]
  elas[,,i] <- matrix(as.vector(sens[,,i])*as.vector(Tmat[,,i])/eig[i],nrow=(2*length(y)))
}

#Take means across 500 iterations
x.w.eig <- rowMeans(w.eigen)
x.dist <- rowMeans(stable.dist)
x.v.eig <- rowMeans(v.eigen)
x.repro <- rowMeans(repro.val)
x.v.w <- mean(v.dot.w)
x.sens<-rowMeans(sens,dims=2) 
x.elas<-rowMeans(elas,dims=2)

#Plot the stable size distribution, reproductive values, elasticity, and sensitivity
library(fields)
par(mfrow=c(2,3),mar=c(4,5,2,2))
image.plot(c(1:244),c(1:244),t(Tmat[,,1]), xlab="Size (t)",ylab="Size (t+1)",
           col=topo.colors(100), main="IPM matrix")
plot(c(1:244),x.dist,xlab="Size",type="l",main="Stable size distribution")
plot(c(1:244),x.repro,xlab="Size",type="l",main="Reproductive values")
image.plot(c(1:244), c(1:244),t(x.elas),xlab="Size (t)",ylab="Size (t+1)",main="Elasticity")
image.plot(c(1:244), c(1:244),t(x.sens),xlab="Size (t)",ylab="Size (t+1)", main="Sensitivity")

## Sensitivity
par(mfrow=c(2,2),mar=c(4,5,2,2))
image.plot(c(1:122), c(1:122), t(x.sens[c(1:122), c(1:122)]),xlab="Size (t)",ylab="Metamorph size (t+1)", col=topo.colors(200), main="Metamorph productivity", zlim=c(min(x.sens), max(x.sens)))
image.plot(c(1:122), c(1:122), t(x.sens[c(1:122), c(123:244)]),xlab="Size (t)",ylab="Metamorph size (t+1)", col=topo.colors(200), main="Adult productivity", zlim=c(min(x.sens), max(x.sens)))
image.plot(c(1:122), c(1:122), t(x.sens[c(123:244), c(1:122)]),xlab="Size (t)",ylab="Size (t+1)", col=topo.colors(200), main="Metamorph growth*survival", zlim=c(min(x.sens), max(x.sens)))
image.plot(c(1:122), c(1:122), t(x.sens[c(123:244), c(123:244)]),xlab="Size (t)",ylab="Size (t+1)", col=topo.colors(200), main="Adult growth*survival", zlim=c(min(x.sens), max(x.sens)))

## Elasticity
par(mfrow=c(2,2),mar=c(4,5,2,2))
image.plot(c(1:122), c(1:122), t(x.elas[c(1:122), c(1:122)]),xlab="Size (t)",ylab="Metamorph size (t+1)", col=topo.colors(200), main="Metamorph productivity", zlim=c(min(x.elas), max(x.elas)))
image.plot(c(1:122), c(1:122), t(x.elas[c(1:122), c(123:244)]),xlab="Size (t)",ylab="Metamorph size (t+1)", col=topo.colors(200), main="Adult productivity", zlim=c(min(x.elas), max(x.elas)))
image.plot(c(1:122), c(1:122), t(x.elas[c(123:244), c(1:122)]),xlab="Size (t)",ylab="Size (t+1)", col=topo.colors(200), main="Metamorph growth*survival", zlim=c(min(x.elas), max(x.elas)))
image.plot(c(1:122), c(1:122), t(x.elas[c(123:244), c(123:244)]),xlab="Size (t)",ylab="Size (t+1)", col=topo.colors(200), main="Adult growth*survival", zlim=c(min(x.elas), max(x.elas)))

library(tidyr)
library(ggplot2)
library(tibble)
library(hrbrthemes)
library(dplyr)
library(cowplot)

x.sens<-data.frame(x.sens)
colnames(x.sens)<-c(1:244)
x.elas<-data.frame(x.elas)
colnames(x.elas)<-c(1:244)

x.sens2<-x.sens %>% mutate(rowid=row_number())
sen.long <- x.sens2 %>% pivot_longer(cols = c(1:244),
                         names_to = "Size", 
                         values_to = "Sensitivity")
sen.long$Size<-as.integer(sen.long$Size)


x.elas2<-x.elas %>% mutate(rowid=row_number())
ela.long <- x.elas2 %>% pivot_longer(cols = c(1:244),
                                     names_to = "Size", 
                                     values_to = "Elasticity")
ela.long$Size<-as.integer(ela.long$Size)

classes<-as.character(c(1:122, 1:122))
# Sensitivity heatmap 
ggplot(data = sen.long, mapping = aes(y=rowid, x=Size, fill = Sensitivity)) +
  geom_tile() +
  xlab(label = "Ln(body mass (g))")+
  ylab(label = "Ln(body mass (g))")+
  scale_x_continuous(breaks=c(1, 61, 123, 183, 244), labels=c("-0.83", "2.15", "-0.83", "2.15", "5.18"))+
  scale_y_reverse(breaks=c(1, 61, 123, 183, 244), labels=c("-0.83", "2.15", "-0.83", "2.15", "5.18"))+
  geom_hline(yintercept=123)+
  geom_vline(xintercept=123)+
  theme_cowplot()+
  annotate("text", label="Offspring produced by metamorphs", y=5, x=61, size=4.5)+
  annotate("text", label="Offspring produced by juveniles/adults", y=5, x=183, size=4.5)+
  annotate("text", label="Metamorph growth*survival", y=127, x=61, size=4.5)+
  annotate("text", label="Juvenile/adult growth*survival", y=127, x=183, size=4.5)+
  scale_fill_distiller(palette="Spectral")

# Elasticity heatmap 
ggplot(data = ela.long, mapping = aes(y=rowid, x=Size, fill = Elasticity)) +
  geom_tile() +
  xlab(label = "Ln(body mass (g))")+
  ylab(label = "Ln(body mass (g))")+
  scale_x_continuous(breaks=c(1, 61, 123, 183, 244), labels=c("-0.83", "2.15", "-0.83", "2.15", "5.18"))+
  scale_y_reverse(breaks=c(1, 61, 123, 183, 244), labels=c("-0.83", "2.15", "-0.83", "2.15", "5.18"))+
  geom_hline(yintercept=123)+
  geom_vline(xintercept=123)+
  theme_cowplot()+
  annotate("text", label="Offspring produced by metamorphs", y=5, x=61, size=4.5)+
  annotate("text", label="Offspring produced by juveniles/adults", y=5, x=183, size=4.5)+
  annotate("text", label="Metamorph growth*survival", y=127, x=61, size=4.5)+
  annotate("text", label="Juvenile/adult growth*survival", y=127, x=183, size=4.5)+
  scale_fill_distiller(palette="Spectral")


####################################################################################################
## Conduct an elasticity analysis of all parameters
## Perturbing each parameter coefficient in the demographic functions by 1%
## Record the percent change in the principal eigenvalue
####################################################################################################
cts1<-array(cts, dim=c(length(cts[,1]), length(cts[1,]), length(cts[,1])))
f1<-array(NA, dim=c(1, length(y), Ndraws, length(cts[,1])))
rec1<-rec.prob1<-array(NA, dim=c(1,length(y),length(cts[,1])))

for (j in 1:length(cts[,1])){
  for (m in 1:Ndraws){
    #Iteratively perturb each parameter
    cts1[j,,j] <- cts[j,] - (cts[j,]*0.01) #perturb by 1%
    #Re-calculate vital rate functions
    f1[,,m,j]<-fx(y,cts1[,m,j])
    rec1[,,j]<-recruits(y,cts1[,m,j])
    rec.prob1[,,j]<-rec1[,,j]/sum(rec1[,,j])#recruitment probability at size y
  }
}

f1<-ifelse(f1<0,0,f1)#Fertility estimates <0 equal 0

Tmat1<-array(NA, dim=c(2*length(y), 2*length(y), Ndraws, length(cts[,1])))
Mpmat1<-Apmat1<-Mmat1<-Amat1<-Mrmat1<-Armat1<-Mfmat1<-Afmat1<-array(0,dim=c(length(y), length(y), Ndraws, length(cts[,1])))
eig1<-array(NA,dim=c(Ndraws, length(cts[,1])))

#Construct transition matrices and store eigenvalues **Requires long run time
for (j in 1:length(cts[,1])){
  for(i in 1:Ndraws){
    for(s in 1:n){
      #Metamorph to one-year-old growth and survival for lower left quadrant of transition matrix
      Mmat1[,,i,j]<-t(outer(y,y,pxy.m,params=cts1[,i,j]))*h
      #juvenile/adult growth and survival for lower right quadrant of transition matrix
      Amat1[,,i,j]<-t(outer(y,y,pxy.a,params=cts1[,i,j]))*h
    }
  }
}

for (j in 1:length(cts[,1])){
  for(i in 1:Ndraws){
    for(s in 1:n){    
      #Metamorph to one-year-old growth and survival (without bin size correction)
      Mpmat1[,,i,j]<-t(outer(y,y,pxy.m,params=cts1[,i,j]))
      #juvenile/adult growth and survival (without bin size correction)
      Apmat1[,,i,j]<-t(outer(y,y,pxy.a,params=cts1[,i,j]))
    }
  }
}

for (j in 1:length(cts[,1])){
  for(i in 1:Ndraws){
    for(s in 1:n){   
      #Metamorphs produced per year by one-year-old females of a size class (without bin size correction)
      Mfmat1[s,,i,j]<-Mpmat1[s,,i,j]*f1[,s,i,j]
      #Metamorphs produced per year by adult females of a size class (without bin size correction)
      Afmat1[s,,i,j]<-Apmat1[s,,i,j]*f1[,s,i,j]
    }
  }
}

for (j in 1:length(cts[,1])){
  for(i in 1:Ndraws){
    for(s in 1:n){  
      #Assign distribution of metamorphs produced of a size class per year by one-year-olds of a size class, upper right quadrant of transition matrix
      Mrmat1[,s,i,j]<-(sum(Mfmat1[,s,i,j]))*rec.prob1[,,j]*h
      #Assign distribution of metamorphs produced of a size class per year by adults of a size class, upper left quadrant of transition matrix
      Armat1[,s,i,j]<-(sum(Afmat1[,s,i,j]))*rec.prob1[,,j]*h
    }
  }
}

for (j in 1:length(cts[,1])){
  for(i in 1:Ndraws){
    for(s in 1:n){ 
      for(k in 1:n){
        #Create the transition matrix
        Tmat1[k,s,i,j]<-Mrmat1[k,s,i,j]
        Tmat1[k,s+122,i,j]<-Armat1[k,s,i,j]
        Tmat1[k+122,s,i,j]<-Mmat1[k,s,i,j]
        Tmat1[k+122,s+122,i,j]<-Amat1[k,s,i,j]
      }
    }#s dimension 1
  }#i draws from posteriors
}#j number of parameters

for (j in 1:length(cts[,1])){
  for(i in 1:Ndraws){
    #Loop through for each transition matrix to calculate deterministic lambda distribution
    eig1[i,j]<- Re(eigen(Tmat1[,,i,j])$value[1])#Re() returns real part of value (i.e., non-imaginary)
  }
}

lambdas1<-if(all(Im(lambdas1<-zapsmall(eig1))==0)) as.numeric(lambdas1) else eig1 #Removes"+0i"
lambdas1<-array(lambdas1, dim=c(Ndraws, length(cts[,1])))

lam.dif<-c(1:length(cts[,1]))

#Find percent difference in mean first eigenvector from original Tmat and matrices using each perturbed parameter
for (j in 1:length(cts[,1])){
  lam.dif[j]<-(abs(mean(lambdas)-mean(lambdas1[,j]))*100)/mean(lambdas)
}
max(abs(lam.dif)) #change of 4.44% at parameter 4, adult growth intercept, second at 1.08% for parameter 1, metamorph growth intercept
plot(c(1:length(cts[,1])), lam.dif)

#Find means across 500 matrices for each coefficient
#x.sens1<-rowMeans(sens1,dims=2) 
#x.elas1<-rowMeans(elas1,dims=2)
param<-c("growth.meta.inter", "growth.meta.slope", "growth.meta.sigma", "growth.adult.inter", "growth.adult.slope", "growth.adult.sigma", 
         "surv.meta.inter", "surv.mass.slope", "surv.adult.inter", "fidelity.meta", "fidelity.adult", 
         "maturity.inter", "maturity.slope", "fertility.inter", "fertility.slope", 
         "min.size", "max.size", "mean.meta.size", "sd.meta.size", "proportion.female", "proportion.breeding", "surv.larvae")
lam.pert<-data.frame(param, lam.dif)
lam.pert$param<-as.factor(lam.pert$param)
plot(lam.pert$param, lam.pert$lam.dif)

####################################################################################################
## Add density dependence on embryonic/larval survival and 
## body size at metamorphosis to 100 year simulation
####################################################################################################
#Larval survival by egg density data and function
lar.egg <- read.csv("Larval-Survival-Density.csv")
larv <- read.csv("larval-survival-posterior-samples.csv", header=T)

#Metamorph body size and egg density data
dens<-read.csv("meta-size-by-egg-density.csv", header=T)

cts2<-matrix(NA,nrow=26,ncol=Ndraws)
cts2[c(1:22),] <- cts

set.seed(23)
cts2[23,] <- rnorm(Ndraws, -0.272, 0.08)  ##Slope value and SD of metamorph size on egg density from Searcy et al. (2015) 
cts2[24,] <- rnorm(Ndraws, 2.28, 0.11)    ##Intercept value set as mean daily log-metamorph masses observed at Olcott over the study period, and 0.11=SE

cts2[25,] <- larv$lar.inter1   		       ## larval survival by egg density intercept
cts2[26,] <- larv$lar.slope1			       ## larval survival by egg density slope

#Test density-dependent functions
#Ln(egg density) calculated using Olcott volume of 101000 m^3
eggl<-sort(lar.egg$ln.dens, decreasing=FALSE)
egg<-sort(dens$ln.egg.dens, decreasing=FALSE)
sim.egg<-seq(0, max(egg), length.out=122)

#Save density dependent growth and survival function outputs
m.d<-array(0, dim=c(length(egg), Ndraws))
f2<-array(0, dim=c(length(y), 1, Ndraws))
for(i in 1:Ndraws){
  m.d[,i]<-met(egg,cts2[,i])
  f2[,,i]<-fx1(y,cts2[,i])
}
f2<-ifelse(f2<0,0,f2)#Fertility estimates <0 equal 0

l.d<-array(0, dim=c(length(egg), Ndraws))
for(i in 1:Ndraws){
  l.d[,i]<-lar(egg,cts2[,i])
}

rec4<-rec.prob4<-array(0, dim=c(1,length(y), length(sim.egg)))
for(j in 1:length(sim.egg)){
  rec4[,,j]<-recruits2(y, sim.egg[j], cts2[,1])
    #dnorm(y, mean = cts2[24,1] + cts2[23,1]*(sim.egg[j]-1.23), sd=0.227)
  rec.prob4[,,j]<-rec4[,,j]/sum(rec4[,,j])#recruitment probability at size y
}


#Examine function outputs
#par(mfrow=c(1,1))
plot(egg, m.d[,1], ylim=c(1.5,3), xlab="Egg density (eggs/m^3)", ylab="Body mass (g)", type='l', bty="l", 
     cex.axis=1.4, cex.lab=1.4, col="grey70", xaxt = "n", yaxt = "n")
axis(1, at=c(0,.5,1,1.5,2), labels=c("1.0", "1.7", "2.7", "4.5", "7.4"), cex.axis=1.4)
axis(2, at=c(1.5,2,2.5,3), labels=c("4.5", "7.4", "12.2", "20.1"), cex.axis=1.4)
mtext("f", 3, .1, at=0.0, cex=1.2)
for(i in 2:Ndraws){
  lines(egg, m.d[,i], col="grey70")
}
md.med<-c(1:8)
for (i in 1:8){
  md.med[i]<-median(m.d[i,])
}
lines(egg, md.med, col=1, lwd=3)

plot(egg, l.d[,1], ylim=c(-7,0), type='l', ylab="Larval survival probability", xlab="Egg density (eggs/m^3)", bty="l", 
     cex.axis=1.4, cex.lab=1.4, col="grey70", xaxt="n", yaxt="n")
axis(1, at=c(0,.5,1,1.5,2), labels=c("1.0", "1.7", "2.7", "4.5", "7.4"), cex.axis=1.4)
axis(2, at=c(-7,-5,-3,-1), labels=c("0.001","0.01", "0.05", "0.37"), cex.axis=1.4)
mtext("g", 3, .1, at=0.0, cex=1.2)
for(i in 2:Ndraws){
  lines(egg,l.d[,i], col="grey70")
}
ld.med<-c(1:8)
for (i in 1:8){
  ld.med[i]<-median(l.d[i,])
}
lines(egg, ld.med, col=1, lwd=3)


plot(y,f2[,,1], type='l', ylim=c(0,550), ylab="Fecundity (metamorphs produced excluding larval survival)", xlab="Ln(body mass (g))", col=2)
for(i in 2:length(f2[,,1])){
  lines(y,f2[,,i], col=2)
}
plot(y,rec4[,,1], type='l', ylab="Probability density", xlab="Ln(metamorph body mass (g))")
for(i in 2:length(rec4[,,1])){
  lines(y,rec4[,,i], col=1)
}
plot(y,rec.prob4[,,1], type='l', ylim=c(0,0.1), ylab="Probability density", xlab="Ln(metamorph body mass (g))")
for(i in 2:length(rec.prob4[,,1])){
  lines(y,rec.prob4[,,i], col=1)
}


# Create matrices to store the body size distribution of juveniles/adults (A1), 
# metamorphs (M1), and the total popuplation (N1) in each year
N1<-A1<-M1<-M1f<-array(NA, dim=c(iter, length(y), Ndraws))#row for every year (1:100), for 122 body size bins, for 500 matrices
rec2<-rec.prob2<-array(NA, dim=c(1, length(y), iter, Ndraws))#1 row, 122 predicted metamorph body size bins, for each annual egg density, for 500 matrices
E<-array(NA, dim=c(1,iter,Ndraws)) #summed egg density at t+1 for each matrix

# The simulation begins with one juvenile/adult in each size class
A1[1,,]<-1
M1[1,,]<-0
N1[1,,]<-0
E[1,,]<-0

#Run density-dependent model
for(i in 1:Ndraws){
  for(t in 1:(iter-1)){
      #Adult population size distribution
      A1[t+1,,i] <- (as.matrix(Mmat[,,i]) %*% M1[t,,i]) + (as.matrix(Amat[,,i]) %*% A1[t,,i]) #Apply growth*survival matrices to current population size distributions
      #Calculate egg density
      E[,t+1,i] <- log((sum(f2[,,i]*A1[t+1,,i])/101000)+1e-20) #ln(egg density) in each year given the number of adult females at t+1 and Olcott volume = 101000 m^3
      #Calculate recruiment probability
      rec2[,,t+1,i] <- dnorm(y, mean = cts2[24,i] + cts2[23,i]*(E[,t+1,i]-1.23), sd=0.227)
      rec.prob2[,,t+1,i] <- rec2[,,t+1,i]/(ifelse(sum(rec2[,,t+1,i])==0, 1e-20, sum(rec2[,,t+1,i])))
      #Density-dependent metamorph population size distribution
      M1f[t+1,,i] <- f2[,,i]*A1[t+1,,i] #Eggs produced per year by adult females of a size class in year t+1
      M1[t+1,,i] <- rec.prob2[,,t+1,i]*sum(f2[,,i]*A1[t+1,,i])*min(exp(lar(d=E[,t+1,i], cts2[,i])), 0.0917) #Assign estimated larval survival from 'lar' function given egg density, or as maximum observed larval survival=0.0917
      #Total population size distribution
      N1[t+1,,i] <- A1[t+1,,i] + M1[t+1,,i]
  }#t
}#i

mean(rec2[,,,1], na.rm=TRUE)
mean(rec.prob2[,,,1], na.rm=TRUE)
exp(mean(E[,100,])) #Back-transformed mean egg density at K = 5.59 eggs/m^3

#Calculate mean metamorph body size at K
M.size <- c(1:Ndraws)
for(i in 1:Ndraws){
  M.size[i] <- cts2[24,i] + cts2[23,i]*(E[,100,i]-1.23)
}
exp(mean(M.size)) #Back-transformed mean metamorph body size at K = 8.58 g

#Calculate mean larval survival at K
larv.surv<-array(NA, dim=c(1,iter,Ndraws)) #larval survival given egg density
for(i in 1:Ndraws){
  for(t in 1:(iter-1)){
    larv.surv[,t+1,i] <- min(exp(lar(d=E[,t+1,i], cts2[,i])), 0.0917)
  }
}
mean(larv.surv[,100,])*100 #1.97% mean larval survival at K

Nt1<-At1<-Mt1<-matrix(NA, iter, Ndraws)
for(i in 1:Ndraws){
  for(t in 1:iter){
    Nt1[t,i]<-sum(N1[t,,i]) #Total population sizes per year
    At1[t,i]<-sum(A1[t,,i])
    Mt1[t,i]<-sum(M1[t,,i])
  }
}
sum(ifelse(Nt1[100,]<122,1,0), na.rm=TRUE)# 13/500 populations declined compared to initial population size

hist(Nt1[100,])
hist(Mt1[100,])
hist(At1[100,])

hist(log(Nt1[100,]))
hist(log(Mt1[100,]))
hist(log(At1[100,]))

#Mean and median carrying capacity (K)
mean(Nt1[100,], na.rm=TRUE) #88875.61 for total population size
mean(At1[100,], na.rm=TRUE) #59605.18 for juvenile/adult population size
mean(Mt1[100,], na.rm=TRUE) #29270.43 for metamorph population size
sd(Nt1[100,], na.rm=TRUE) #364985.6
sd(At1[100,], na.rm=TRUE) #275423
sd(Mt1[100,], na.rm=TRUE) #89998.65
median(Nt1[100,], na.rm=TRUE) #60330.32 for total population size
median(At1[100,], na.rm=TRUE) #37506.07 for juvenile/adult population size
median(Mt1[100,], na.rm=TRUE) #21960.54 for metamorph population size

plot(c(1:100), rowMeans(Nt1[1:100,]), type='l', col="blue", ylab="Mean population size", xlab="Year")
lines(c(1:100), rowMeans(Mt1[1:100,]), type='l')
lines(c(1:100), rowMeans(At1[1:100,]), type='l', col=2)

#Try excluding populations that have crashed
mean(ifelse(Nt1[100,]>3, Nt1[100,], NA), na.rm=TRUE) #90874.85 for total population size
mean(ifelse(At1[100,]>3, At1[100,], NA), na.rm=TRUE) #60945.99 for juvenile/adult population size
mean(ifelse(Mt1[100,]>3, Mt1[100,], NA), na.rm=TRUE) #29990.19 for metamorph population size
median(ifelse(Nt1[100,]>3, Nt1[100,], NA), na.rm=TRUE) #61005.16 for total population size
median(ifelse(At1[100,]>3, At1[100,], NA), na.rm=TRUE) #38888.85 for juvenile/adult population size
median(ifelse(Mt1[100,]>3, Mt1[100,], NA), na.rm=TRUE) #22114.39 for metamorph population size

plot(c(1:100), Nt1[1:100,1], ylim=c(0,300000), type='l', ylab="Ln(population size)", xlab="Year")
for (i in 2:Ndraws){
  lines(c(1:100), Nt1[1:100,i], col=i)
}
abline(75560.3, 0, lwd=5)#Mean K
abline(58915.25, 0, lwd=5, col="blue")# Median K
plot(c(1:100), Mt1[1:100,1], ylim=c(0,130000), type='l', ylab="Ln(metamorph population size)", xlab="Year")
for (i in 2:Ndraws){
  lines(c(1:100), Mt1[1:100,i], col=i)
}
abline(25795.51, 0, lwd=5)#Mean K
abline(22235.32, 0, lwd=5, col="blue")# Median K
plot(c(1:100), At1[1:100,1], ylim=c(0,150000), type='l', ylab="Ln(adult population size)", xlab="Year")
for (i in 2:Ndraws){
  lines(c(1:100), At1[1:100,i], col=i)
}
abline(49764.79, 0, lwd=5)#Mean K
abline(36591.74, 0, lwd=5, col="blue")# Median K

#Examine log-transformed plots
plot(c(1:100), log(Nt1[1:100,1]), ylim=c(-60,20), type='l', ylab="Ln(population size)", xlab="Year")
for (i in 2:Ndraws){
  lines(c(1:100), log(Nt1[1:100,i]), col=i)
}
plot(c(1:100), log(Mt1[1:100,1]), ylim=c(-60,20), type='l', ylab="Ln(metamorph population size)", xlab="Year")
for (i in 2:Ndraws){
  lines(c(1:100), log(Mt1[1:100,i]), col=i)#incomplete lines result of log(0)
}
plot(c(1:100), log(At1[1:100,1]), ylim=c(-60,20), type='l', ylab="Ln(adult population size)", xlab="Year")
for (i in 2:Ndraws){
  lines(c(1:100), log(At1[1:100,i]), col=i)
}

####################################################################################################
## Conduct an elasticity analysis of all parameters under density-dependence
## Perturbing each parameter coefficient in the demographic functions by 1%
## Record the percent change in the carrying capacity
####################################################################################################
cts3<-array(cts2, dim=c(length(cts2[,1]), length(cts2[1,]), length(cts2[,1])))
Mmat2<-Amat2<-array(0,dim=c(length(y), length(y), Ndraws, length(cts2[,1])))
f3<-array(NA, dim=c(1, length(y), Ndraws, length(cts2[,1])))
rec3<-rec.prob3<-array(NA, dim=c(1, length(y), iter, Ndraws, length(cts2[,1])))
N2<-A2<-M2<-M2f<-array(NA, dim=c(iter, length(y), Ndraws, length(cts2[,1])))
E1<-array(NA, dim=c(1,iter,Ndraws, length(cts2[,1])))

for (j in 1:length(cts2[,1])){
  for (m in 1:Ndraws){
    #Iteratively perturb each parameter
    cts3[j,,j] <- cts2[j,] - (cts2[j,]*0.01) #perturb by 1% #perturb by 1%
    #Re-calculate fecundity function
    f3[,,m,j]<-fx1(y,cts3[,m,j])
  }
}

f3<-ifelse(f3<0,0,f3)#Fertility estimates <0 equal 0

#Each simulation begins with one juvenile/adult in each size class
A2[1,,,]<-1
M2[1,,,]<-0
N2[1,,,]<-0
E1[1,,,]<-0

#Construct transition matrices and store eigenvalues **Requires long run time
for (j in 1:length(cts2[,1])){
  for(i in 1:Ndraws){
    for(s in 1:n){
      #Metamorph to one-year-old growth and survival for lower left quadrant of transition matrix
      Mmat2[,,i,j]<-t(outer(y,y,pxy.m,params=cts3[,i,j]))*h
      #kuvenile/adult growth and survival for lower right quadrant of transition matrix
      Amat2[,,i,j]<-t(outer(y,y,pxy.a,params=cts3[,i,j]))*h
    }
  }
}

#Run density-dependent model with 1% perturbation
for (j in 1:length(cts2[,1])){
  for(i in 1:Ndraws){
    for(t in 1:(iter-1)){
      #Adult population size distribution
      A2[t+1,,i,j] <- (as.matrix(Mmat2[,,i,j]) %*% M2[t,,i,j]) + (as.matrix(Amat2[,,i,j]) %*% A2[t,,i,j]) #Apply growth*survival matrices to current population size distributions
      #Calculate egg density
      E1[,t+1,i,j] <- log((sum(f3[,,i,j]*A2[t+1,,i,j])/101000)+1e-20) #ln(egg density) in each year given the number of adult females at t+1 and Olcott volume = 101000 m^3
      #Calculate recruiment probability
      rec3[,,t+1,i,j] <- dnorm(y, mean = cts3[24,i,j] + cts3[23,i,j]*(E1[,t+1,i,j]-1.23), sd=0.227)
      rec.prob3[,,t+1,i,j] <- rec3[,,t+1,i,j]/(ifelse(sum(rec3[,,t+1,i,j])==0, 1e-20, sum(rec3[,,t+1,i,j])))
      #Density-dependent metamorph population size distribution
      M2f[t+1,,i,j] <- f3[,,i,j]*A2[t+1,,i,j] #Eggs produced per year by adult females of a size class in year t+1
      M2[t+1,,i,j] <- rec.prob3[,,t+1,i,j]*sum(f3[,,i,j]*A2[t+1,,i,j])*min(exp(lar(d=E1[,t+1,i,j], cts3[,i,j])), 0.0917) #Assign estimated larval survival from 'lar' function given egg density, or as maximum observed larval survival=0.0917
      #Total population size distribution
      N2[t+1,,i,j] <- A2[t+1,,i,j] + M2[t+1,,i,j]
    }#t
  }#i
}#j

Nt2<-At2<-Mt2<-array(NA, dim=c(iter, Ndraws, length(cts2[,1])))
for (j in 1:length(cts2[,1])){
  for(i in 1:Ndraws){
    for(t in 1:iter){
      Nt2[t,i,j]<-sum(N2[t,,i,j]) #Total population sizes per year
      At2[t,i,j]<-sum(A2[t,,i,j])
      Mt2[t,i,j]<-sum(M2[t,,i,j])
    }
  }
}
  
hist(Nt2[100,,])
hist(Mt2[100,,])
hist(At2[100,,])

Ndif<-Adif<-Mdif<-array(NA, dim=c(1, length(cts2[,1])))
#Percent change in median carrying capacity
for (j in 1:length(cts2[,1])){
  Ndif[,j] <- (abs(median(Nt1[100,])-median(Nt2[100,,j])))*100/median(Nt1[100,])
  Adif[,j] <- (abs(median(At1[100,])-median(At2[100,,j])))*100/median(At1[100,])  
  Mdif[,j] <- (abs(median(Mt1[100,])-median(Mt2[100,,j])))*100/median(Mt1[100,])
}

plot(c(1:length(cts2[,1])), Ndif)#parameter 4=juvenile/adult growth intercept, then parameter 25 = larval survival intercept
plot(c(1:length(cts2[,1])), Adif)
plot(c(1:length(cts2[,1])), Mdif)
Ndif

####################################################################################################
## Build in environmental stochasticity given
## the number of females breeding given December-January precipitation and
## the probability of reproductive/replacement success given December-January precipitation.
## Test for goodness of fit based on the predicted and observed numbers of metamorphs.
####################################################################################################
#Create matrices to store the body size distribution of juveniles/adults (Ak) starting from carrying capacity (A1[100,,]), 
#and the change in number of juveniles/adults between the two time steps (Adelta)
Ak<-array(NA, dim=c(2, length(y), Ndraws))#row for one year change from K, for 122 body size bins, for 500 matrices
AJk<-array(NA, dim=c(1, length(y), Ndraws))# Matrices to store the number of mature juveniles/adults (i.e., adults only) given body size (AJk)
Adelta<-array(NA, dim=c(1, length(y), Ndraws))
  
#Number of juveniles/adults in each size class at carrying capacity
AJk[1,,]<-A1[100,,]

for(i in 1:Ndraws){
  for(s in 1:n){
    #Estimate the number of adults only at carrying capacity (i.e., mature individuals)
    Ak[1,s,i] <- as.matrix(AJk[,s,i]*ma[,s,i]) #ma is maturity function output by body size bin created early in this script
  }
}

Aksum<-array(NA, dim=c(1, Ndraws))
for(i in 1:Ndraws){
    Aksum[,i]<-sum(Ak[1,,i]) #Total adult-only population size at carrying capacity
}

median(Aksum)

for(i in 1:Ndraws){
  #Adult-only population size distribution at carrying capacity (K)
  Ak[2,,i] <- as.matrix(Ak[1,,i]*s.a[,,i]) #Apply juvenile/adult survival arrays (s.a created above) to K population size distributions
  #Calculate adults lost between two time points in each body size class
  Adelta[,,i]<-Ak[1,,i]-Ak[2,,i]
}

#Calculate matrix-specific totals
Ak1<-matrix(NA, 2, Ndraws)
Adelta1<-matrix(NA, 1, Ndraws)

for(i in 1:Ndraws){
  for(t in 1:2){
    Ak1[t,i]<-sum(Ak[t,,i]) #Total adult population sizes at K and following time step
  }
}

for(i in 1:Ndraws){
  Adelta1[,i]<-sum(Adelta[,,i]) #Total change in adult population size
}

median(Adelta1) #3327.254 median adults lost annually from K

plot(density(Ak1[1,]), ylim=c(0, 7e-05)) #At K
lines(density(Ak1[2,]), col=2) #at next time step

hist(Adelta1) #Change in adult population size

Adelta2<-matrix(NA, 1, Ndraws)
for(i in 1:Ndraws){
  Adelta2[,i]<-Adelta1[,i]*0.159 #Total change in adult population size sampled by area of Olcott drift fence 15.9%
}

median(Adelta2) #529.03 median adults lost annually from K sampled from drift fence

##Load posterior samples of proportion femlaes breeding function, reproductive success inflection points, and replacement success inflection points
fem<-read.csv("females-precip-posterior-samples-CENTERED.csv", header=T)
replace<-read.csv("replace-success-infection-samples-CENTERED.csv", header=TRUE)
repro<-read.csv("repro-success-infection-samples-CENTERED.csv", header=TRUE)

#Calculate regression between each pair of inflection points
infl<-cbind.data.frame(repro$x, replace$x)*10 #Convert from cm to mm
infl$inter<-rep(NA, length(infl[,1]))
infl$slope<-rep(NA, length(infl[,1]))
probs<-c(0,1)

mean(infl[,1])#43.1 cm mean inflection point for reproductive success
mean(infl[,2])#56.5 cm mean inflection point for replacement

for(i in 1:length(infl[,1])){
  mod<-lm(probs~c(repro$x[i]*10, replace$x[i]*10))#Convert from cm to mm
  infl[i,c(3:4)]<-coef(mod)
}

#Specify stepwise function for inflection pairs where reproduction cutoff>replacement cutoff
#Such a situation is not ecologically possible
sum(ifelse(infl[,1]>=infl[,2],1,0))
infl[,3]<-ifelse(infl[,1]>=infl[,2], infl[,3]==NA, infl[,3])
infl[,4]<-ifelse(infl[,1]>=infl[,2], infl[,4]==NA, infl[,4])
infl[,2]<-ifelse(infl[,1]>=infl[,2], infl[,1]+.1, infl[,2])

plot(c(infl[1,1:2]), probs, xlim=c(300,700), bty="l", cex.axis=1.4, cex.lab=1.4, cex=1, col="white", 
     xlab="October-June precipitation (mm)", ylab="Probability of reproductive success")
abline(a=infl[1,3], b=infl[1,4], col="grey70")
mtext("i", 3, .1, at=300, cex=1.2)
for(i in 2:length(infl[,1])){
  points(infl[i,1:2], probs, col="grey70", cex=1)
}
infl.omit<-na.omit(infl)
for(i in 2:length(infl.omit[,1])){
  abline(a=infl.omit[i,3], b=infl.omit[i,4], col="grey70")
}
#Add pairwise connections
segments(x0=582.5, y0=-1, x1=582.5, y1=1.5, col="grey70")
segments(x0=596, y0=-1, x1=596, y1=1.5, col="grey70")
segments(x0=592, y0=-1, x1=592, y1=1.5, col="grey70")
segments(x0=620.5, y0=-1, x1=620.5, y1=1.5, col="grey70")
segments(x0=611, y0=-1, x1=612, y1=1.5, col="grey70")
#Add median lines
segments(x0=0, y0=0, x1=max(infl[,1]), y1=0, col="grey70", lwd=1)
segments(x0=min(infl[,2]), y0=1, x1=720, y1=1, col="grey70", lwd=1)
segments(x0=median(infl[,1]), y0=0, x1=median(infl[,2]), y1=1, col=1, lwd=3)
segments(x0=0, y0=0, x1=median(infl[,1]), y1=0, col=1, lwd=3)
segments(x0=median(infl[,2]), y0=1, x1=720, y1=1, col=1, lwd=3)
#dev.off()

##Add proportion female breeding function to IPM data frame
cts4<-matrix(NA,nrow=28,ncol=Ndraws)
cts4[c(1:26),] <- cts2
cts4[27,] <- fem$fem.inter1  		       ## proportion females breeding by Dec-Jan rainfall intercept
cts4[28,] <- fem$fem.slope1			       ## proportion females breeding by Dec-Jan rainfall slope


##Add environmental dependency to the density-dependent IPM
#Read in 96-year climate record data
clim<-read.table("Stochastic_Climate_Pool.txt",header=T)

#Set the number of years of simulation plus one (iter2)
iter2<-97

#Calculate precipitation- and density-dependent fecundity
f4<-array(NA, dim=c(length(y), 1, iter2-1, Ndraws))
fems<-array(NA, dim=c(1,iter2-1,Ndraws))

for(i in 1:Ndraws){
  for(t in 1:(iter2-1)){
    fems[,t,i]<-nfem.b1(clim$Dec.Jan[t], cts4[,i])#proportion females breeding given Dec-Jan rainfall
    f4[,,t,i]<-fx2(y, clim$Dec.Jan[t], cts4[,i])#fecundity given precipitation-dependent proportion females breeding
  }
}
mean(fems[])#Long-term mean females breeding = 37%
f4<-ifelse(f4<0,0,f4)#Fecundity estimates <0 equal 0

#Inspect outputs
plot(sort(clim$Dec.Jan), sort(fems[,,1]), ylim=c(0,1), xlim=c(0,650), ylab="Proportion of females breeding", xlab="December-January precipitation (mm)", 
     col="grey70", cex.axis=1.4, cex.lab=1.4, bty="l", type="l")
mtext("h", 3, .1, at=2, cex=1.2)
for(i in 1:Ndraws){
  lines(sort(clim$Dec.Jan),sort(fems[,,i]), col="grey70")
}
fems.med<-c(1:(iter2-1))
for(i in 1:(iter2-1)){
  fems.med[i] <-median(fems[,i,])
}
lines(sort(clim$Dec.Jan), sort(fems.med), col=1, lwd=3)



plot(y,f4[,,1,1], type='l', ylim=c(0,550), ylab="Fecundity (metamorphs produced excluding larval survival)", xlab="Ln(body mass (g))", col=2)
for(i in 1:Ndraws){
  for(t in 1:(iter2-1)){
    lines(y,f4[,,t,i], col=2)
  }
}

#Create matrices to hold the simulated size distributions of metamorphs (M3), adults (A3), and the total population (N3)
N3<-A3<-M3<-array(NA, dim=c(iter2, length(y), Ndraws))
rec4<-rec.prob4<-array(NA, dim=c(1, length(y), iter2, Ndraws))#1 row, 122 predicted metamorph body size bins, for each annual egg density, for 500 matrices
E2<-array(NA, dim=c(1,iter2,Ndraws)) #summed egg density at t+1 for each matrix

#Simulation begins with juveniles/adults at carrying capacity
A3[1,,]<-A1[100,,]
M3[1,,]<-0
N3[1,,]<-0
E2[1,,]<-0

#Run density- and environment-dependent model
for(i in 1:Ndraws){
  for(t in 1:(iter2-1)){
    #If Oct-Jun rainfall is less than reproduction inflection point [i], then the number of metamorphs recruited to the population is zero
    if(clim$Oct.Jun[t]<infl[i,1]){
      #Metamorph recruitment and population size distribution is zero (total breeding failure due to drying)
      M3[t+1,,i] <- rec4[,,t+1,i] <- rec.prob4[,,t+1,i] <- 0
      #Adult population size distribution
      A3[t+1,,i] <- (as.matrix(Mmat[,,i]) %*% M3[t,,i]) + (as.matrix(Amat[,,i]) %*% A3[t,,i]) #Apply growth*survival matrices to current population size distributions
      #Calculate egg density
      E2[,t+1,i] <- log((sum(f4[,,t,i]*A3[t+1,,i])/101000)+1e-20) #ln(egg density) in each year given the number of adult females at t+1 and Olcott volume = 101000 m^3 and proportion of females breeding
    }
    
    #If Oct-Jun rainfall is between reproduction inflection point [i] and replacement inflection point [i], then the 
    #number of metamorphs recruited to the popualtion is a fraction of the individuals that survived to late-stage larvae ((clim$Oct.Jun[t]-infl[i,1])/(infl[i,2]-infl[i,1]))
    if(clim$Oct.Jun[t]>=infl[i,1] && clim$Oct.Jun[t]<=infl[i,2]){
      #Adult population size distribution
      A3[t+1,,i] <- (as.matrix(Mmat[,,i]) %*% M3[t,,i]) + (as.matrix(Amat[,,i]) %*% A3[t,,i]) #Apply growth*survival matrices to current population size distributions
      #Calculate egg density given proportion females breeding
      E2[,t+1,i] <- log((sum(f4[,,t,i]*A3[t+1,,i])/101000)+1e-20) #ln(egg density) in each year given the number of adult females at t+1, Olcott volume = 101000 m^3, and proportion females breeding
      #Calculate recruiment probability
      rec4[,,t+1,i] <- dnorm(y, mean = cts4[24,i] + cts4[23,i]*(E2[,t+1,i]-1.23), sd=0.227)
      rec.prob4[,,t+1,i] <- rec4[,,t+1,i]/(ifelse(sum(rec4[,,t+1,i])==0, 1e-20, sum(rec4[,,t+1,i])))
      #Metamorph population size distribution given proportion females breeding
      M3[t+1,,i] <- ((clim$Oct.Jun[t]-infl[i,1])/(infl[i,2]-infl[i,1]))*rec.prob4[,,t+1,i]*sum(f4[,,t,i]*A3[t+1,,i])*min(exp(lar(d=E2[,t+1,i], cts4[,i])), 0.0917) #Assign estimated larval survival from 'lar' function given egg density, or as maximum observed larval survival=0.0917
    }
    
    #If Oct-Jun rainfall is greater than replacement inflection point [i], then the number of metamorphs recruited to
    #the population equals the number of individuals that survived to late-stage larvae
    if(clim$Oct.Jun[t]>infl[i,2]){
      A3[t+1,,i] <- (as.matrix(Mmat[,,i]) %*% M3[t,,i]) + (as.matrix(Amat[,,i]) %*% A3[t,,i]) #Apply growth*survival matrices to current population size distributions
      #Calculate egg density
      E2[,t+1,i] <- log((sum(f4[,,t,i]*A3[t+1,,i])/101000)+1e-20) #ln(egg density) in each year given the number of adult females at t+1 and Olcott volume = 101000 m^3
      #Calculate recruiment probability
      rec4[,,t+1,i] <- dnorm(y, mean = cts4[24,i] + cts4[23,i]*(E2[,t+1,i]-1.23), sd=0.227)
      rec.prob4[,,t+1,i] <- rec4[,,t+1,i]/(ifelse(sum(rec4[,,t+1,i])==0, 1e-20, sum(rec4[,,t+1,i])))
      #Metamorph population size distribution
      M3[t+1,,i] <- rec.prob4[,,t+1,i]*sum(f4[,,t,i]*A3[t+1,,i])*min(exp(lar(d=E2[,t+1,i], cts4[,i])), 0.0917) #Assign estimated larval survival from 'lar' function given egg density, or as maximum observed larval survival=0.0917
    }
    
    #Total population size distribution
    N3[t+1,,i] <- A3[t+1,,i] + M3[t+1,,i]
  }#t
}#i

fails <- matrix(NA, iter2, Ndraws)#Number of reproductive failures given matrix iterations across 96 years
fails.sum <- c(1:Ndraws)
for(i in 1:Ndraws){
  for(t in 1:iter2){
    fails[t,i] <- ifelse(clim$Oct.Jun[t]<infl[i,1], 1, 0)
    fails.sum[i] <-sum(fails[,i], na.rm = TRUE)
  }
}
mean(fails.sum/96) #On average across matrices, 24% of years feature complete reproductive failure

rec.fail <- matrix(NA, iter2, Ndraws)#Number of replacement failures given matrix iterations across 96 years
rec.fail.sum <- c(1:Ndraws)
for(i in 1:Ndraws){
  for(t in 1:iter2){
    rec.fail[t,i] <- ifelse(clim$Oct.Jun[t]<infl[i,2], 1, 0)
    rec.fail.sum[i] <-sum(rec.fail[,i], na.rm = TRUE)
  }
}
mean(rec.fail.sum/96) #On average across matrices, 50% of years feature failure to recruit enough metamorphs to achieve replacement

Nt3<-At3<-Mt3<-matrix(NA, iter2, Ndraws)
for(i in 1:Ndraws){
  for(t in 1:iter2){
    Nt3[t,i]<-sum(N3[t,,i]) #Total population sizes per year
    At3[t,i]<-sum(A3[t,,i])
    Mt3[t,i]<-sum(M3[t,,i])
  }
}

hist(Nt3[97,])
hist(Mt3[97,])
hist(At3[97,])

#Mean final abundance
mean(Nt3[97,], na.rm=TRUE) #33145.42 for total population size
mean(At3[97,], na.rm=TRUE) #31266.29 for juvenile/adult population size
mean(Mt3[97,], na.rm=TRUE) #1879.13 for metamorph population size
sd(Nt3[97,], na.rm=TRUE) #179069.7
sd(At3[97,], na.rm=TRUE) #177348.5
sd(Mt3[97,], na.rm=TRUE) #5677.984

#Median final abundance
median(Nt3[97,], na.rm=TRUE) #18384.8 for total population size
median(At3[97,], na.rm=TRUE) #17909.36 for juvenile/adult population size
median(Mt3[97,], na.rm=TRUE) #0 for metamorph population size

par(mfrow=c(2,1))
plot(c(1892, clim$Year), rowMeans(Nt3[1:97,]), type='l', col="blue", ylab="Mean population size", xlab="Year")
lines(c(1892, clim$Year), rowMeans(Mt3[1:97,]), type='l')
lines(c(1892, clim$Year), rowMeans(At3[1:97,]), type='l', col=2)
plot(clim$Year, clim$Oct.Jun, type='l', col=3, lty=1, ylab="October-June precipitation (mm)", xlab="Year")
abline(median(infl[,1]), 0, lty=2, col="black")#Median reproduction threshold
abline(median(infl[,2]), 0, lty=3, col="black")#Median replacement threshold

#Compare predicted to observed metamorphs for 2005-2013
pred.m<-c(1:9)
for (i in 1:9){
  pred.m[i] <- median(Mt3[i+88,])*0.159 #Account for percent of Olcott fenced
}

#Load in observed metamorph data (2007 and 2008 had total reproductive failure)
m05<-as.matrix(read.table("Olcott2005.txt", header=T))
m06<-as.matrix(read.table("Olcott2006.txt", header=T))
m09<-as.matrix(read.table("Olcott2009.txt", header=T))
m10<-as.matrix(read.table("Olcott2010.txt", header=T))
m11<-as.matrix(read.table("Olcott2011.txt", header=T))
m12<-as.matrix(read.table("Olcott2012.txt", header=T))
m13<-as.matrix(read.table("Olcott2013.txt", header=T))

#Create vector of annual total metamorphs observed
obs.m<-c(sum(m05), sum(m06), 0, 0, sum(m09), sum(m10), sum(m11), sum(m12), sum(m13))

pdf("Figure4.pdf", width = 9, height = 7)
par(mfrow=c(1,1))
plot(pred.m, obs.m, ylab="Observed number of metamorphs", xlab="Predicted number of metamorphs", bty="l",
     pch=19, cex=1.5, cex.axis=1.4, cex.lab=1.4)
abline(lm(obs.m~pred.m), lwd=3)
abline(0, 1, col=2, lty=2, lwd=3)
legend(2700,600, legend=c("Linear model", "1:1 line"), col=c(1,2), lty=c(1,2), lwd=3, cex=1.4, bty="n")
mtext(expression(R^2~"="~0.94), 1, -7, at=3080, cex=1.4)
mtext(expression(P~"<"~"0.0001"), 1, -5.5, at=3120, cex=1.4)
mtext(expression(P~"="~0.98), 1, -1.5, at=3050, cex=1.4, col=2)
dev.off()

gof.m<-lm(obs.m~pred.m)
summary(gof.m)#R^2=0.94

#Does the relationship between observed and predicted metamorph values significantly differ from the one-to-one line?
#Simulate one-to-one line
perfect.x <- perfect.y <- seq(min(obs.m), max(obs.m), 426.5)
plot(perfect.x, perfect.y)

#Create data frame with true and simulated data
y.gof<-c(obs.m, perfect.y)
x.gof<-c(pred.m, perfect.x)
group<-c(rep("true", 9), rep("perfect", 9))
test.gof<-cbind.data.frame(x.gof, y.gof, group)
test.gof$group<-as.factor(test.gof$group)

#Test whether true and simulated regression lines significantly differ
fit<-lm(y.gof~x.gof*group, test.gof)
summary(fit)
library(car)
Anova(fit, type="III")#No significant difference between slopes of true values and one-to-one perfect fit

#Examine all transition matrices
pred.m.all<-matrix(NA, 9, Ndraws)
for (j in 1:Ndraws){
  for (i in 1:9){
    pred.m.all[i,j] <- Mt3[i+88,j]*0.159 #Account for percent of Olcott fenced
  }
}

x.preds <- c(1:Ndraws)
for (i in 1:Ndraws){
  x.preds[i] <- mean(pred.m.all[,i])
}

sum(ifelse(x.preds>mean(obs.m), 1, 0)) #299 (60%) matrices have greater mean annual metamorph predictions than mean observed

inters<-slopes<-c(1:Ndraws)
for (i in 1:Ndraws){
  fit1<-lm(obs.m~pred.m.all[,i])
  coeffs<-c(coefficients(fit1))
  inters[i]<-as.numeric(coeffs[1])
  slopes[i]<-as.numeric(coeffs[2])
  slopes[i]<-ifelse(is.na(slopes[i])==TRUE, 0, slopes[i])#9 slopes NA due to total metamorph production failure across years
}

plot(pred.m.all[,1], obs.m, ylab="Observed number of metamorphs", xlab="Predicted number of metamorphs", col=0, ylim=c(0,30000))
for (i in 1:Ndraws){
  abline(inters[i], slopes[i])
}
abline(0, 1, col=2, lty=2, lwd=4)
legend(3200,30000, legend=c("Linear model", "1:1 line"), col=c(1,2), lty=c(1,2), lwd=c(1,4))

####################################################################################################
## Determine population size at which all PVA model runs should begin
####################################################################################################
#200 year simulation
iter3<-201

#Read in Nut Tree Airport-only climate record data
clim1<-read.table("Stochastic_Climate_Pool_Rev.txt",header=T)

#Sample across historic record for 200 year simulation, and repeat for 100 random iterations
historic<-array(NA, dim=c(iter3, 3, 100))
set.seed(9237)
for(h in 1:100){
  samples <- sample(nrow(clim1), iter3, replace=TRUE) #Sample historic climate record for 200-yr simulation
  historic[,,h] <- as.matrix(clim1[samples,])
}

#Create matrices to hold the simulated size distributions of metamorphs (M4), adults (A4), and the total population (N4)
f5<-array(NA, dim=c(length(y), 1))
N4<-array(NA, dim=c(iter3, length(y), Ndraws))
A4<-M4<-array(NA, dim=c(iter3, length(y), Ndraws))
rec5<-rec.prob5<-array(NA, dim=c(1, length(y)))#1 row, 122 predicted metamorph body size bins
E3<-array(NA, dim=c(1,iter3,Ndraws)) #summed egg density at t+1 for each matrix
Nt4<-array(NA, dim=c(iter3-1, Ndraws, 100)) #Store total population sizes per year given transition matrix and climate simulation
M4all<-A4all<-array(NA, dim=c(iter3, length(y), Ndraws, 100)) #Store metamorph and juvenile/adult population size distributions for all iterations for use in age at maturity calculation below
  
#Simulation begins with juveniles/adults at carrying capacity
A4[1,,]<-A1[100,,]
A4all[1,,,]<-A1[100,,]
M4[1,,]<-0
M4all[1,,,]<-0
N4[1,,]<-A1[100,,]
E3[1,,]<-0

#Run density- and environment-dependent model (Requires ~2 hr run time)
for (h in 1:100){
  for(i in 1:Ndraws){
    for(t in 1:(iter3-1)){
      #Calculate precipitation- and density-dependent fecundity across 100 randomly selected climate data sets
      f5<-fx2(y, historic[t,2,h], cts4[,i])
      f5<-ifelse(f5<0,0,f5)#Fecundity estimates <0 equal 0
      #If Oct-Jun rainfall is less than reproduction inflection point [i], then the number of metamorphs recruited to the population is zero
      if(historic[t,3,h]<infl[i,1]){
        #Metamorph recruitment and population size distribution is zero (total breeding failure due to drying)
        M4[t+1,,i] <- rec5 <- rec.prob5 <- 0
        #Adult population size distribution
        A4[t+1,,i] <- (as.matrix(Mmat[,,i]) %*% M4[t,,i]) + (as.matrix(Amat[,,i]) %*% A4[t,,i]) #Apply growth*survival matrices to current population size distributions
        #Calculate egg density
        E3[,t+1,i] <- log((sum(f5*A4[t+1,,i])/101000)+1e-20) #ln(egg density) in each year given the number of adult females at t+1 and Olcott volume = 101000 m^3 and proportion of females breeding
      }
      
      #If Oct-Jun rainfall is between reproduction inflection point [i] and replacement inflection point [i], then the 
      #number of metamorphs recruited to the popualtion is a fraction of the individuals that survived to late-stage larvae ((historic$Oct.Jun[t]-infl[i,1])/(infl[i,2]-infl[i,1]))
      if(historic[t,3,h]>=infl[i,1] && historic[t,3,h]<=infl[i,2]){
        #Adult population size distribution
        A4[t+1,,i] <- (as.matrix(Mmat[,,i]) %*% M4[t,,i]) + (as.matrix(Amat[,,i]) %*% A4[t,,i]) #Apply growth*survival matrices to current population size distributions
        #Calculate egg density given proportion females breeding
        E3[,t+1,i] <- log((sum(f5*A4[t+1,,i])/101000)+1e-20) #ln(egg density) in each year given the number of adult females at t+1, Olcott volume = 101000 m^3, and proportion females breeding
        #Calculate recruiment probability
        rec5 <- dnorm(y, mean = cts4[24,i] + cts4[23,i]*(E3[,t+1,i]-1.23), sd=0.227)
        rec.prob5 <- rec5/(ifelse(sum(rec5)==0, 1e-20, sum(rec5)))
        #Metamorph population size distribution given proportion females breeding
        M4[t+1,,i] <- ((historic[t,3,h]-infl[i,1])/(infl[i,2]-infl[i,1]))*rec.prob5*sum(f5*A4[t+1,,i])*min(exp(lar(d=E3[,t+1,i], cts4[,i])), 0.0917) #Assign estimated larval survival from 'lar' function given egg density, or as maximum observed larval survival=0.0917
        
      }
      
      #If Oct-Jun rainfall is greater than replacement inflection point [i], then the number of metamorphs recruited to
      #the population equals the number of individuals that survived to late-stage larvae
      if(historic[t,3,h]>infl[i,2]){
        A4[t+1,,i] <- (as.matrix(Mmat[,,i]) %*% M4[t,,i]) + (as.matrix(Amat[,,i]) %*% A4[t,,i]) #Apply growth*survival matrices to current population size distributions
        #Calculate egg density
        E3[,t+1,i] <- log((sum(f5*A4[t+1,,i])/101000)+1e-20) #ln(egg density) in each year given the number of adult females at t+1 and Olcott volume = 101000 m^3
        #Calculate recruiment probability
        rec5 <- dnorm(y, mean = cts4[24,i] + cts4[23,i]*(E3[,t+1,i]-1.23), sd=0.227)
        rec.prob5 <- rec5/(ifelse(sum(rec5)==0, 1e-20, sum(rec5)))
        #Metamorph population size distribution
        M4[t+1,,i] <- rec.prob5*sum(f5*A4[t+1,,i])*min(exp(lar(d=E3[,t+1,i], cts4[,i])), 0.0917) #Assign estimated larval survival from 'lar' function given egg density, or as maximum observed larval survival=0.0917
      }
      
      #Total population size distribution
      N4[t+1,,i] <- A4[t+1,,i] + M4[t+1,,i]
      #Total population sizes per year given transition matrix and climate simulation
      Nt4[t,i,h] <- sum(N4[t+1,,i])
      #Stage-specific body size distributions per year given transition matrix and climate simulation
      M4all[t+1,,i,h] <- M4[t+1,,i]
      A4all[t+1,,i,h] <- A4[t+1,,i]
    }#t
  }#i
  print("running") #allows you to evaluate progress of model
}#h

#Mean abundance of second 100 years of simulation
hist(Nt4[c(101:200),,])
mean(Nt4[c(101:200),,], na.rm=TRUE) #44421.76 for total population size
sd(Nt4[c(101:200),,], na.rm=TRUE) #118354.1

#Median abundance of second 100 years of simulation
median(Nt4[c(101:200),,], na.rm=TRUE) #53678.9 for total population size
median(ifelse(Nt4[c(101:200),,]>3, Nt4[c(101:200),,], NA), na.rm=TRUE) #31485.69 with crashed populations removed

NtpropK <- median(Nt4[c(101:200),,], na.rm=TRUE) / median(Nt1[100,], na.rm=TRUE) * 100 #53.88% overall median percent of K
median(ifelse(Nt4[c(101:200),,]>3, Nt4[c(101:200),,], NA), na.rm=TRUE) / median(ifelse(Nt1[100,]>3, Nt1[100,], NA), na.rm=TRUE) * 100 #54.94 with crashed populations removed

#Visualize some simulations
plot(c(101:200), Nt4[101:200,1,1], ylim=c(0,300000), type='l', ylab="Ln(population size)", xlab="Year")
for (i in 2:Ndraws){
  lines(c(101:200), Nt4[101:200,i,1], col=i)
}
abline(median(Nt1[100,], na.rm=TRUE), 0, lwd=5, col="blue")# Median K
abline(median(Nt4[c(101:200),,1], na.rm=TRUE), 0, lwd=5, col="red")# Median total population size


plot(c(101:200), Nt4[101:200,1,100], ylim=c(0,300000), type='l', ylab="Ln(population size)", xlab="Year")
for (i in 2:Ndraws){
  lines(c(101:200), Nt4[101:200,i,100], col=i)
}
abline(median(Nt1[100,], na.rm=TRUE), 0, lwd=5, col="blue")# Median K
abline(median(Nt4[c(101:200),,100], na.rm=TRUE), 0, lwd=5, col="red")# Median total population size


#Calculate climate simulation-specific median total population sizes and percent of K to examine variance
Ntmed<-Ntprob<-rep(NA, 100)
for (h in 1:100){
  Ntmed[h]<-median(Nt4[c(101:200),,h], na.rm = TRUE)#Medians by climate simulations
  Ntprob[h]<-median(Nt4[c(101:200),,h], na.rm = TRUE)/median(Nt1[100,], na.rm=TRUE)*100#Percents of K
}

hist(Ntmed, xlim=c(4000, 70000))
abline(v=median(Nt1[100,], na.rm=TRUE), lwd=5, col="blue")# Median K

hist(Ntprob, xlab="Percent of K")

mean(Ntprob)#Mean of climate simulation-specific medians = 50.94%

####################################################################################################
## Calculate mean age at maturity and post-metamorphic stage-specific survival
####################################################################################################
#M4all[c(101:200),,,] is the above-generated population size distribution of one-year old individuals given stochastic climate scenario and after burn-in
#Save summed individuals in body size bins (dim 2) for years 101-200 for each of 500 matrices (dim 3) across the 100 climate scenarios
Mall1 <- apply(M4all[c(101:200),,,], c(2,3), sum)

#M1ma is the population size distribution of mature one-year-olds after burn-in
#M1nma is the population size distribution of non-mature one-year-olds after burn-in
M1s<-M1ma<-M1nma<-Mall2<-array(NA, dim=c(10, length(y), Ndraws))#number mature one-year-olds (summed across 100 years and 100 climate simulations) in each of 10 consecutive years in 122 body size bins, for 500 matrices
Mall2[1,,]<-Mall1

#For loops for year one, when individuals start in the metamorph stage
for(i in 1:Ndraws){
  M1s[1,,i] <- as.matrix(Mmat[,,i]) %*% Mall2[1,,i] #Apply metamorph growth and survival function to all metamorphs
}

for(i in 1:Ndraws){
  for(s in 1:n){
    M1ma[1,s,i] <- as.matrix(M1s[1,s,i]*ma[,s,i]) #ma is maturity function output by body size bin created early in this script, applied to survivors
    M1nma[1,s,i] <- as.matrix(M1s[1,s,i] - M1ma[1,s,i]) #Immature individuals from first year move on as juveniles
  }
}

#For loop for years 2-10, when individuals from simulated cohort are in the juvenile/adult stage
for(t in 1:9){  
  for(i in 1:Ndraws){
    
    M1s[t+1,,i] <- as.matrix(Amat[,,i]) %*% M1nma[t,,i] #Apply juvenile/adult growth*survival matrix to current immature population size distributions
    
    for(s in 1:n){
      #Estimate the number of adults only at carrying capacity (i.e., mature individuals)
      M1ma[t+1,s,i] <- as.matrix(M1s[t+1,s,i]*ma[,s,i]) #Apply maturity function to juveniles that grew and survived to next year
      M1nma[t+1,s,i] <- M1s[t+1,s,i] - M1ma[t+1,s,i] #Immature individuals
    }
  }
}

#Sanity checks
plot(Mall1[,5], typ="l")
points(M1s[1,,5])
points(M1s[2,,5], col=2)
points(M1s[3,,5], col=3)
points(M1s[4,,5], col=4)
points(M1s[5,,5], col=5)
points(M1s[6,,5], col=6)

plot(Mall1[,5], typ="l")
points(M1nma[1,,5])
points(M1nma[2,,5], col=2)
points(M1nma[3,,5], col=3)
points(M1nma[4,,5], col=4)
points(M1nma[5,,5], col=5)
points(M1nma[6,,5], col=6)

plot(M1ma[1,,5], ylim=c(0,2100000))
points(M1ma[2,,5], col=2)
points(M1ma[3,,5], col=3)
points(M1ma[4,,5], col=4)
points(M1ma[5,,5], col=5)
points(M1ma[6,,5], col=6)
points(M1ma[7,,5], col=7)
points(M1ma[8,,5], col=8)
points(M1ma[9,,5], col=9)
points(M1ma[10,,5], col=10)


mat<-nmat<-matrix(NA, 10, Ndraws)
for(t in 1:10){
  mat[t,] <- apply(M1ma[t,,], 2, sum, na.rm=TRUE) #Number of mature individuals per year per matrix
  nmat[t,] <- apply(M1nma[t,,], 2, sum, na.rm=TRUE) #Number of immature individuals per year per matrix
}

age.mat <- data.frame(apply(mat, 1, mean, na.rm=TRUE)) #Calculate mean annual Number of mature individuals across 500 matrix iterations
age.mat$n <- c(1:10)
names(age.mat) <- c("mature", "year")
age.mat$mature <- as.integer(age.mat$mature)
age.mat1 <- rep(age.mat[,2], age.mat[,1])#Expand dataframe across frequency of mature individuals per year

#Calculate age at maturity across years
hist(age.mat1)
x.agemat <- mean(age.mat1) #mean age at maturity
median(age.mat1) #median age at maturity



### Calculate means for annual survival of metamorphs, juveniles, and adults
#Total number of individuals that enter the 10-year simulation
indT <- apply(Mall2[1,,], 2, sum, na.rm=TRUE)

#Total number of individuals that reach maturity across the 10-year simulation
maT<-c(1:Ndraws)
for (i in 1: Ndraws){
  maT[i] <- sum(mat[,i], na.rm=TRUE)
}

#Calculate pre-maturity survival (metamorphs+juveniles)
#Solve for mean(indT) * (mean annual survival[x]^x.agemat) = mean(maT)
  #x.agemat * log((mean annual survival[x])) = log(mean(maT)/mean(indT))
  #mean annual survival[x] = exp(log(mean(maT)/mean(indT))/x.agemat)
x.surv <- exp(log(mean(maT)/mean(indT))/x.agemat) #Prematurity-terrestrial survival


#Calculate adult and juvenile mean survivals
#Take the size distribution summed across years 101-200 and all climate simulations
Aall1 <- apply(A4all[c(101:200),,,], c(2,3), sum)

J1s<-A1s<-Aall2<-array(NA, dim=c(2, length(y), Ndraws))
Aall2[1,,]<-Aall1

#First year, determine life stages
for(i in 1:Ndraws){
  for(s in 1:n){
    A1s[1,s,i] <- as.matrix(Aall2[1,s,i] * ma[,s,i]) #Apply maturity function to determine adults
    J1s[1,s,i] <- as.matrix(Aall2[1,s,i] - A1s[1,s,i]) #Determine juveniles
  }
}

#Apply survival function to each life stage
for(i in 1:Ndraws){
  for(s in 1:n){
    #Estimate the number of adults only at carrying capacity (i.e., mature individuals)
    A1s[2,s,i] <- A1s[1,s,i]*s.a[,s,i] #Adult survival
    J1s[2,s,i] <- J1s[1,s,i]*s.a[,s,i] #Juvenile survival
  }
}

#Calculate ratio of stage specific survival vs. the total number of individuals that entered the simulation
#Total number of stage-specific individuals that enter thesimulation
aT1 <- apply(A1s[1,,], 2, sum, na.rm=TRUE)
jT1 <- apply(J1s[1,,], 2, sum, na.rm=TRUE)

#Total number of stage-specific individuals surviving the simulation
aT2 <- apply(A1s[2,,], 2, sum, na.rm=TRUE)
jT2 <- apply(J1s[2,,], 2, sum, na.rm=TRUE)

#Calculate ratios
a.surv<-mean(aT2/aT1, na.rm=TRUE) #Adult-only survival
j.surv<-mean(jT2/jT1, na.rm=TRUE) #Juvenile-only survival

js<-as<-wj<-wa<-xj<-xa<-matrix(NA, 8, Ndraws)
ms<-wm<-xm<-c(1:500)


#Repeat this approach to calculate mean metamorph survival in the stochastic model
M2s<-array(NA, dim=c(2, length(y), Ndraws))
M2s[1,,] <- Mall1

#Apply survival function to metamorph size distribution
for(i in 1:Ndraws){
  for(s in 1:n){
    M2s[2,s,i] <- M2s[1,s,i]*s.m[,s,i]
  }
}

#Total number of stage-specific individuals that enter thesimulation
mT1 <- apply(M2s[1,,], 2, sum, na.rm=TRUE)

#Total number of stage-specific individuals surviving the simulation
mT2 <- apply(M2s[2,,], 2, sum, na.rm=TRUE)

#Calculate ratios
m.surv<-mean(mT2/mT1, na.rm=TRUE) #Metamorph survival

####################################################################################################
## Population Viability Analysis
## Using the density- and environment-dependent IPM, calculate extinction probability of CTS populations
## under 80 different scenarios (20 preserve sizes x 2 pond sizes x 2 climatic regimes)
####################################################################################################

#20 preserve sizes from 100 to 2000 m in 100-m increments
#Load data for life stage density given distance from pond (calculated using methods in Searcy & Shaffer 2008; 2011)
dens.dist <- read.csv("density-distance.csv", header=TRUE)

#Two pond sizes of 101,000 m^3 (area of Olcott Lake, i.e., biggest) and 933 m^3 (area of Blomquist Pond, i.e., average)
pond <- c(101000, 933)

#Two climatic regimes:
# --Random draws of historic (1893-2013) weather data (see "clim1" data and "historic" simulated data object created above)
# --Future regime under clim1ate change based on trends in both the mean and variance of our two 
#   (December-January and October-June) precipitation variables. If such trends are present, 
#   we extrapolated their values in 2100 (see "future" simulated data object created below).

##Test for Dec-Jan trends
plot(clim1$Year, clim1$Dec.Jan)
abline(lm(clim1$Dec.Jan~clim1$Year))

mod.dec.jan <- lm(clim1$Dec.Jan~clim1$Year)
summary(mod.dec.jan)#No significant trend in Dec-Jan rainfall
resid.dec <- resid(mod.dec.jan) #Extract residuals
sigma(mod.dec.jan) #148.65

plot(clim1$Year, abs(resid.dec), xlim=c(1893,2100))
abline(lm(abs(resid.dec)~clim1$Year))

mod.dec.resid <- lm(abs(resid.dec)~clim1$Year)
summary(mod.dec.resid) #Significant positive trend in variance over time (p = 0.003; r^2 = 0.09)

#Extrapolate absolute residual variance value in 2100
var.dec <- (0.8445*2100)-1532.6465 #240.80

##Test for Oct-Jun trends
plot(clim1$Year, clim1$Oct.Jun)
abline(lm(Oct.Jun~Year, clim1))

mod.oct.jun<-lm(Oct.Jun~Year, clim1)
summary(mod.oct.jun)#No significant trend in Oct-Jun rainfall
resid.oct <- resid(mod.oct.jun) #Extract residuals
sigma(mod.oct.jun) #231.47

plot(clim1$Year, abs(resid.oct), xlim=c(1893,2100))
abline(lm(abs(resid.oct)~clim1$Year))

mod.oct.resid <- lm(abs(resid.oct)~clim1$Year)
summary(mod.oct.resid) #Significant positive trend in variance over time (p = 0.015; r^2 = 0.03)

var.oct <- (1.0323*2100)-1828.3297 #339.50

#Build Figure A3
pdf("FigureA3.pdf", width = 10, height = 10)
par(mfrow=c(2,2), mar = c(5, 5, 2.5, 1.8))
plot(clim1$Year, clim1$Dec.Jan, ylab="Cumulative precipitation (mm)", xlab="", bty="l",
     pch=19, cex=1.5, cex.axis=1.4, cex.lab=1.4, main="December-January")
abline(lm(clim1$Dec.Jan~clim1$Year), lwd=3)
mtext("A", 3, .1, at=1900, cex=1.4)
mtext(expression(R^2~"="~"0.006"), 3, -3, at=1910, cex=1.3)
mtext(expression(P~"="~"0.48"), 3, -4.5, at=1910, cex=1.3)
plot(clim1$Year, clim1$Oct.Jun, ylab="", xlab="", bty="l",
     pch=19, cex=1.5, cex.axis=1.4, cex.lab=1.4, main="October-June")
abline(lm(Oct.Jun~Year, clim1), lwd=3)
mtext("B", 3, .1, at=1900, cex=1.4)
mtext(expression(R^2~"="~"0.01"), 3, -3, at=1910, cex=1.3)
mtext(expression(P~"="~"0.30"), 3, -4.5, at=1910, cex=1.3)
plot(clim1$Year, abs(resid.dec), ylab="Absolute(residuals (mm))", xlab="Year", bty="l",
     pch=19, cex=1.5, cex.axis=1.4, cex.lab=1.4)
abline(lm(abs(resid.dec)~clim1$Year), lwd=3)
mtext("C", 3, .1, at=1900, cex=1.4)
mtext(expression(R^2~"="~"0.09"), 3, -3, at=1910, cex=1.3)
mtext(expression(P~"="~"0.003"), 3, -4.5, at=1910, cex=1.3)
plot(clim1$Year, abs(resid.oct), ylab="", xlab="Year", bty="l",
     pch=19, cex=1.5, cex.axis=1.4, cex.lab=1.4)
abline(lm(abs(resid.oct)~clim1$Year), lwd=3)
mtext("D", 3, .1, at=1900, cex=1.4)
mtext(expression(R^2~"="~"0.06"), 3, -3, at=1910, cex=1.3)
mtext(expression(P~"="~"0.015"), 3, -4.5, at=1910, cex=1.3)
dev.off()

#Relationship between two time frames of precipitation
mod.comp <- lm(Oct.Jun~Dec.Jan, clim1)
summary(mod.comp)

mean(clim1$Dec.Jan)
mean(clim1$Oct.Jun)
sd(clim1$Dec.Jan)
sd(clim1$Oct.Jun)

##Simulate 2100 future climate variation 100 times
set.seed(3463)
iter4<-101
future<-array(NA, dim=c(iter4, 2, 100))
for (f in 1:100){
  for (t in 1:iter4){
    #Dec-Jan rainfall in 2100 simulated as a random draw from a parametric distribution
    future[t,1,f] <- 254 + rnorm(1,0,370) #Where SD is value that produces desired residual in 2100
    #Correct any Dec-Jan rainfall values below zero to zero
    future[t,1,f] <- ifelse(future[t,1,f]<0, 0, future[t,1,f])
    #Oct-Jun rainfall in 2100 simulated as a random draw from a parametric distribution around the regression of Oct-Jun rainfall on Dec-Jan rainfall
    future[t,2,f] <- 1.15 * future[t,1,f] + rnorm(1,0,284) + 310 #Where SD is additional needed to produced desired residual in 2100, 310 is intercept and 1.15 is slope of the regression
    #If Oct-Jun rainfall is less than Dec-Jan rainfall, set equal to Dec-Jan rainfall
    future[t,2,f] <- ifelse(future[t,2,f] < future[t,1,f], future[t,1,f], future[t,2,f])
  }
}

plot(clim1$Year, clim1$Dec.Jan, xlim=c(1893, 2100), ylim=c(0,1000))
abline(lm(clim1$Dec.Jan~clim1$Year))
points(rep(2100,100), future[1,1,], col=2) #just examine first year

plot(clim1$Year, clim1$Oct.Jun, xlim=c(1893, 2100), ylim=c(0,2000))
abline(lm(clim1$Oct.Jun~clim1$Year))
points(rep(2100,100), future[1,2,], col=2) #just examine first year

#Do the simulated residuals approach predicted mean residual values?
mean(clim1$Dec.Jan) #253.95 mean Dec-Jan precip value
(1.15*253.95)+309.73 #601.77 mean Oct-Jun precip value as calculated from mean Dec-Jan precip
f.dec.resid<-future[,1,]-253.95
f.oct.resid<-future[,2,]-601.77
mean(abs(f.dec.resid))#240.03, target is 240.03
mean(abs(f.oct.resid))#339.46, target is 339.46


##########################################################################################################
##Run PVA for 100 historic climate simulations over 100 years
f5<-array(NA, dim=c(length(y), 1))
N5<-A5<-M5<-array(NA, dim=c(iter4, length(y), Ndraws))
rec5<-rec.prob5<-array(NA, dim=c(1, length(y)))#1 row, 122 predicted metamorph body size bins
E3<-array(NA, dim=c(1,iter4,Ndraws)) #summed egg density at t+1 for each matrix
At5<-Mt5<-Nt5<-array(NA, dim=c(iter4-1, Ndraws, 100, 20, 2)) #Store total population sizes per year given transition matrix, climate simulation, preserve size, and pond volume
A.min<-N.min<-array(NA, dim=c(1, Ndraws, 100, 20, 2)) #Store minimum population sizes given transition matrix, climate simulation, preserve size, and pond volume

#Simulation begins with juveniles/adults at 50.94% of carrying capacity (the historic median population size)
A5[1,,]<-0.5094*A1[100,,]
M5[1,,]<-0
N5[1,,]<-0.5094*A1[100,,]
E3[1,,]<-0

# Specify function for the fraction of juvenile (i.e., not mature; nm) individuals to incorporate density of juveniles given distance from pond
nm <- 1-ma #where ma is the maturity function output by body size bin created early in this script

#Run density- and environment-dependent model for historic climate simulation, both pond volumes, and all 20 preserve sizes
#Resquires long run time (~26 hrs)
for(p in 1:length(pond)){
  for(j in 1:length(dens.dist[,1])){
    for(h in 1:length(historic[1,1,])){
      for(i in 1:Ndraws){
        for(t in 1:(iter4-1)){
          #Calculate precipitation- and density-dependent fecundity across 100 randomly selected climate data sets
          f5<-fx2(y, historic[t,2,h], cts4[,i])
          f5<-ifelse(f5<0,0,f5)#Fecundity estimates <0 equal 0
          #If Oct-Jun rainfall is less than reproduction inflection point [i], then the number of metamorphs recruited to the population is zero
          if(historic[t,3,h]<infl[i,1]){
            #Metamorph recruitment and population size distribution is zero (total breeding failure due to drying)
            M5[t+1,,i] <- rec5 <- rec.prob5 <- 0
            #Adult/juvenile population size distribution
            A5[t+1,,i] <- (as.matrix(Mmat[,,i] %*% (M5[t,,i]*dens.dist[j,4]))) + (as.matrix(Amat[,,i] %*% (A5[t,,i]*ma[,,i]*dens.dist[j,2]))) + (as.matrix(Amat[,,i] %*% (A5[t,,i]*nm[,,i]*dens.dist[j,3]))) #Apply growth*survival matrices to current population size distributions of metamorphs, juveniles, and adults, accounting for stage-specific density given distance from pond
            #Calculate egg density
            E3[,t+1,i] <- log((sum(f5*A5[t+1,,i])/pond[p])+1e-20) #ln(egg density) in each year given the number of adult females at t+1, pond volume and proportion of females breeding
          }
          
          #If Oct-Jun rainfall is between reproduction inflection point [i] and replacement inflection point [i], then the 
          #number of metamorphs recruited to the popualtion is a fraction of the individuals that survived to late-stage larvae ((historic$Oct.Jun[t]-infl[i,1])/(infl[i,2]-infl[i,1]))
          if(historic[t,3,h]>=infl[i,1] && historic[t,3,h]<=infl[i,2]){
            #Adult/juvenile population size distribution
            A5[t+1,,i] <- (as.matrix(Mmat[,,i] %*% (M5[t,,i]*dens.dist[j,4]))) + (as.matrix(Amat[,,i] %*% (A5[t,,i]*ma[,,i]*dens.dist[j,2]))) + (as.matrix(Amat[,,i] %*% (A5[t,,i]*nm[,,i]*dens.dist[j,3]))) #Apply growth*survival matrices to current population size distributions of metamorphs, juveniles, and adults, accounting for stage-specific density given distance from pond
            #Calculate egg density given proportion females breeding
            E3[,t+1,i] <- log((sum(f5*A5[t+1,,i])/pond[p])+1e-20) #ln(egg density) in each year given the number of adult females at t+1, pond volume, and proportion females breeding
            #Calculate recruiment probability
            rec5 <- dnorm(y, mean = cts4[24,i] + cts4[23,i]*(E3[,t+1,i]-1.23), sd=0.227)
            rec.prob5 <- rec5/(ifelse(sum(rec5)==0, 1e-20, sum(rec5)))
            #Metamorph population size distribution given proportion females breeding
            M5[t+1,,i] <- ((historic[t,3,h]-infl[i,1])/(infl[i,2]-infl[i,1]))*rec.prob5*sum(f5*A5[t+1,,i])*min(exp(lar(d=E3[,t+1,i], cts4[,i])), 0.0917) #Assign estimated larval survival from 'lar' function given egg density, or as maximum observed larval survival=0.0917
          }
          
          #If Oct-Jun rainfall is greater than replacement inflection point [i], then the number of metamorphs recruited to
          #the population equals the number of individuals that survived to late-stage larvae
          if(historic[t,3,h]>infl[i,2]){
            #Adult/juvenile population size distribution
            A5[t+1,,i] <- (as.matrix(Mmat[,,i] %*% (M5[t,,i]*dens.dist[j,4]))) + (as.matrix(Amat[,,i] %*% (A5[t,,i]*ma[,,i]*dens.dist[j,2]))) + (as.matrix(Amat[,,i] %*% (A5[t,,i]*nm[,,i]*dens.dist[j,3]))) #Apply growth*survival matrices to current population size distributions of metamorphs, juveniles, and adults, accounting for stage-specific density given distance from pond
            #Calculate egg density
            E3[,t+1,i] <- log((sum(f5*A5[t+1,,i])/pond[p])+1e-20) #ln(egg density) in each year given the number of adult females at t+1 and pond volume
            #Calculate recruiment probability
            rec5 <- dnorm(y, mean = cts4[24,i] + cts4[23,i]*(E3[,t+1,i]-1.23), sd=0.227)
            rec.prob5 <- rec5/(ifelse(sum(rec5)==0, 1e-20, sum(rec5)))
            #Metamorph population size distribution
            M5[t+1,,i] <- rec.prob5*sum(f5*A5[t+1,,i])*min(exp(lar(d=E3[,t+1,i], cts4[,i])), 0.0917) #Assign estimated larval survival from 'lar' function given egg density, or as maximum observed larval survival=0.0917
          }
          
          #Total population size distribution
          N5[t+1,,i] <- A5[t+1,,i] + M5[t+1,,i]
          #Total population sizes per year given transition matrix, climate simulation, preserve size, and pond volume
          Nt5[t,i,h,j,p] <- sum(N5[t+1,,i])
          At5[t,i,h,j,p] <- sum(A5[t+1,,i])
          Mt5[t,i,h,j,p] <- sum(M5[t+1,,i])
        }#t (simulation years)
        
        #Minimum annual population size per matrix iteration and simulation
        N.min[,i,h,j,p] <- min(Nt5[,i,h,j,p])
        A.min[,i,h,j,p] <- min(At5[,i,h,j,p])#Juveniles/adults only
        
      }#i (transition matrix iterations)
      print(h) #allows you to evaluate progress of the model as it runs
    }#h (historic climate simulations)
    print("preserve")
  }#j (preserve sizes)
  print("pond")
}#p (pond volumes)

A.ext<-N.ext<-prob.ext<-array(NA, dim=c(1, Ndraws, 20, 2))
for(p in 1:length(pond)){
  for(j in 1:length(dens.dist[,1])){
    for(i in 1:Ndraws){
      #Times quasi-extinction threshold (3 individuals) reached for each matrix under each climate simulation
      N.ext[,i,j,p] <- sum(ifelse(N.min[,i,,j,p]<3, 1, NA), na.rm=TRUE)
      #Times quasi-extinction threshold (3 juveniles/adults) reached for each matrix under each climate simulation
      A.ext[,i,j,p] <- sum(ifelse(A.min[,i,,j,p]<3, 1, NA), na.rm=TRUE)
      #Extinction probability for each matrix, preserve size, and pond volume given 100 climate simulations
      prob.ext[,i,j,p]<-sum(A.ext[,i,j,p])/100 
    }
  }
}

pdf("Figure5.pdf", width = 10, height = 8)
par(mfrow=c(2,2), mar = c(5, 6, 2, 1.8))
plot(dens.dist$dist.m, prob.ext[,1,c(1:20), 2], col="grey70", typ='l', ylim=c(0, 1), cex.axis=1.4, cex.lab=1.4,
     ylab="Extinction probability (100 years)", xlab="", bty="l", main = "Average carrying capacity")#933 m^3 pond
mtext("Historic climate", side = 2, line = 4.7, cex=1.2, font=2)
for(i in 2:Ndraws){
  lines(dens.dist$dist.m, prob.ext[,i,c(1:20), 2], col="grey70")
}
med.prob.933 <- c(1:20)
for(i in 1:20){
  med.prob.933[i] <- median(prob.ext[,,i,2])
}
lines(dens.dist$dist.m, med.prob.933, col=1, pch=20, cex=3, lwd=4)
mtext("a", 3, .1, at=200, cex=1.5)

plot(dens.dist$dist.m, prob.ext[,1,c(1:20), 1], col="grey70", typ='l', ylim=c(0, 1), cex.axis=1.4, cex.lab=1.4,
     ylab="", xlab="", bty="l", main = "High carrying capacity")#101000 m^3 pond
for(i in 2:Ndraws){
  lines(dens.dist$dist.m, prob.ext[,i,c(1:20), 1], col="grey70")
}
med.prob.101 <- c(1:20)
for(i in 1:20){
  med.prob.101[i] <- median(prob.ext[,,i,1])
}
lines(dens.dist$dist.m, med.prob.101, col=1, pch=20, cex=3, lwd=4)
mtext("b", 3, .1, at=200, cex=1.5)


#Median extinction probabilities at each preserve radii
olcott.ext<-blomquist.ext<-matrix(NA, 20, 1)
for(i in 1:20){
  olcott.ext[i,] <- median(prob.ext[,,i,1])
  blomquist.ext[i,] <- median(prob.ext[,,i,2])
}

#Number of iterations very low extinction probability (<=2%) at each preserve radii
blom.min<-olc.min<-c(1:20)
for(i in 1:20){
  blom.min[i]<-sum(ifelse(prob.ext[,,i,2]<=.02,1,0))#97.4% of iterations very likely go extinct if 700 m radii, 100% if less
  olc.min[i]<-sum(ifelse(prob.ext[,,i,1]<=.02,1,0))#92.8% of iterations very likely go extinct if 300 m radii or less, 100% if less
}

##########################################################################################################
##Run PVA for 100 climate simulations given predicted 2100 precipitation variance over 100 years
f5<-array(NA, dim=c(length(y), 1))
N6<-A6<-M6<-array(NA, dim=c(iter4, length(y), Ndraws))
rec5<-rec.prob5<-array(NA, dim=c(1, length(y)))#1 row, 122 predicted metamorph body size bins
E4<-array(NA, dim=c(1,iter4,Ndraws)) #summed egg density at t+1 for each matrix
At6<-Mt6<-Nt6<-array(NA, dim=c(iter4-1, Ndraws, 100, 20, 2)) #Store total population sizes per year given transition matrix, climate simulation, preserve size, and pond volume
A.min1<-N.min1<-array(NA, dim=c(1, Ndraws, 100, 20, 2)) #Store minimum population sizes given transition matrix, climate simulation, preserve size, and pond volume

#Simulation begins with juveniles/adults at 53.88% of carrying capacity (the future median population size)
A6[1,,]<-0.5388*A1[100,,]
M6[1,,]<-0
N6[1,,]<-0.5388*A1[100,,]
E4[1,,]<-0

# Specify function for the fraction of juvenile (i.e., not mature; nm) individuals to incorporate density of juveniles given distance from pond
nm <- 1-ma #where ma is the maturity function output by body size bin created early in this script

#Run density- and environment-dependent model for future climate simulation, both pond volumes, and all 20 preserve sizes
#Resquires long run time (~26 hrs)
for(p in 1:length(pond)){
  for(j in 1:length(dens.dist[,1])){
    for(h in 1:length(future[1,1,])){
      for(i in 1:Ndraws){
        for(t in 1:(iter4-1)){
          #Calculate precipitation- and density-dependent fecundity across 100 randomly selected climate data sets
          f5<-fx2(y, future[t,1,h], cts4[,i])
          f5<-ifelse(f5<0,0,f5)#Fecundity estimates <0 equal 0
          #If Oct-Jun rainfall is less than reproduction inflection point [i], then the number of metamorphs recruited to the population is zero
          if(future[t,2,h]<infl[i,1]){
            #Metamorph recruitment and population size distribution is zero (total breeding failure due to drying)
            M6[t+1,,i] <- rec5 <- rec.prob5 <- 0
            #Adult/juvenile population size distribution
            A6[t+1,,i] <- (as.matrix(Mmat[,,i] %*% (M6[t,,i]*dens.dist[j,4]))) + (as.matrix(Amat[,,i] %*% (A6[t,,i]*ma[,,i]*dens.dist[j,2]))) + (as.matrix(Amat[,,i] %*% (A6[t,,i]*nm[,,i]*dens.dist[j,3]))) #Apply growth*survival matrices to current population size distributions of metamorphs, juveniles, and adults, accounting for stage-specific density given distance from pond
            #Calculate egg density
            E4[,t+1,i] <- log((sum(f5*A6[t+1,,i])/pond[p])+1e-20) #ln(egg density) in each year given the number of adult females at t+1, pond volume and proportion of females breeding
          }
          
          #If Oct-Jun rainfall is between reproduction inflection point [i] and replacement inflection point [i], then the 
          #number of metamorphs recruited to the popualtion is a fraction of the individuals that survived to late-stage larvae ((future[t,2,h]-infl[i,1])/(infl[i,2]-infl[i,1]))
          if(future[t,2,h]>=infl[i,1] && future[t,2,h]<=infl[i,2]){
            #Adult/juvenile population size distribution
            A6[t+1,,i] <- (as.matrix(Mmat[,,i] %*% (M6[t,,i]*dens.dist[j,4]))) + (as.matrix(Amat[,,i] %*% (A6[t,,i]*ma[,,i]*dens.dist[j,2]))) + (as.matrix(Amat[,,i] %*% (A6[t,,i]*nm[,,i]*dens.dist[j,3]))) #Apply growth*survival matrices to current population size distributions of metamorphs, juveniles, and adults, accounting for stage-specific density given distance from pond
            #Calculate egg density given proportion females breeding
            E4[,t+1,i] <- log((sum(f5*A6[t+1,,i])/pond[p])+1e-20) #ln(egg density) in each year given the number of adult females at t+1, pond volume, and proportion females breeding
            #Calculate recruiment probability
            rec5 <- dnorm(y, mean = cts4[24,i] + cts4[23,i]*(E4[,t+1,i]-1.23), sd=0.227)
            rec.prob5 <- rec5/(ifelse(sum(rec5)==0, 1e-20, sum(rec5)))
            #Metamorph population size distribution given proportion females breeding
            M6[t+1,,i] <- ((future[t,2,h]-infl[i,1])/(infl[i,2]-infl[i,1]))*rec.prob5*sum(f5*A6[t+1,,i])*min(exp(lar(d=E4[,t+1,i], cts4[,i])), 0.0917) #Assign estimated larval survival from 'lar' function given egg density, or as maximum observed larval survival=0.0917
          }
          
          #If Oct-Jun rainfall is greater than replacement inflection point [i], then the number of metamorphs recruited to
          #the population equals the number of individuals that survived to late-stage larvae
          if(future[t,2,h]>infl[i,2]){
            #Adult/juvenile population size distribution
            A6[t+1,,i] <- (as.matrix(Mmat[,,i] %*% (M6[t,,i]*dens.dist[j,4]))) + (as.matrix(Amat[,,i] %*% (A6[t,,i]*ma[,,i]*dens.dist[j,2]))) + (as.matrix(Amat[,,i] %*% (A6[t,,i]*nm[,,i]*dens.dist[j,3]))) #Apply growth*survival matrices to current population size distributions of metamorphs, juveniles, and adults, accounting for stage-specific density given distance from pond
            #Calculate egg density
            E4[,t+1,i] <- log((sum(f5*A6[t+1,,i])/pond[p])+1e-20) #ln(egg density) in each year given the number of adult females at t+1 and pond volume
            #Calculate recruiment probability
            rec5 <- dnorm(y, mean = cts4[24,i] + cts4[23,i]*(E4[,t+1,i]-1.23), sd=0.227)
            rec.prob5 <- rec5/(ifelse(sum(rec5)==0, 1e-20, sum(rec5)))
            #Metamorph population size distribution
            M6[t+1,,i] <- rec.prob5*sum(f5*A6[t+1,,i])*min(exp(lar(d=E4[,t+1,i], cts4[,i])), 0.0917) #Assign estimated larval survival from 'lar' function given egg density, or as maximum observed larval survival=0.0917
          }
          
          #Total population size distribution
          N6[t+1,,i] <- A6[t+1,,i] + M6[t+1,,i]
          #Total population sizes per year given transition matrix, climate simulation, preserve size, and pond volume
          Nt6[t,i,h,j,p] <- sum(N6[t+1,,i])
          At6[t,i,h,j,p] <- sum(A6[t+1,,i])
          Mt6[t,i,h,j,p] <- sum(M6[t+1,,i])
        }#t (simulation years)
        
        #Minimum annual population size per matrix iteration and simulation
        N.min1[,i,h,j,p] <- min(Nt6[,i,h,j,p])
        A.min1[,i,h,j,p] <- min(At6[,i,h,j,p])#Juveniles/adults only
        
      }#i (transition matrix iterations)
      print(h) #allows you to evaluate progress of the model as it runs
    }#h (future climate simulations)
    print("preserve")
  }#j (preserve sizes)
  print("pond")
}#p (pond volumes)

A.ext1<-N.ext1<-prob.ext1<-array(NA, dim=c(1, Ndraws, 20, 2))
for(p in 1:length(pond)){
  for(j in 1:length(dens.dist[,1])){
    for(i in 1:Ndraws){
      #Times quasi-extinction threshold (3 individuals) reached for each matrix under each climate simulation
      N.ext1[,i,j,p] <- sum(ifelse(N.min1[,i,,j,p]<3, 1, NA), na.rm=TRUE)
      #Times quasi-extinction threshold (3 juveniles/adults) reached for each matrix under each climate simulation
      A.ext1[,i,j,p] <- sum(ifelse(A.min1[,i,,j,p]<3, 1, NA), na.rm=TRUE)
      #Extinction probability for each matrix, preserve size, and pond volume given 100 climate simulations
      prob.ext1[,i,j,p]<-sum(A.ext1[,i,j,p])/100 
    }
  }
}

plot(dens.dist$dist.m, prob.ext1[,1,c(1:20), 2], col="grey70", typ='l', ylim=c(0, 1), cex.axis=1.4, cex.lab=1.4,
     ylab="Extinction probability (100 years)", xlab="Preserve radius (m)", bty="l")#933 m^3 pond
mtext("Future climate", side = 2, line = 4.7, cex=1.2, font=2)
for(i in 2:Ndraws){
  lines(dens.dist$dist.m, prob.ext1[,i,c(1:20), 2], col="grey70")
}
med.prob.933 <- c(1:20)
for(i in 1:20){
  med.prob.933[i] <- median(prob.ext1[,,i,2])
}
lines(dens.dist$dist.m, med.prob.933, col=1, pch=20, cex=3, lwd=4)
mtext("c", 3, .1, at=200, cex=1.5)

plot(dens.dist$dist.m, prob.ext1[,1,c(1:20), 1], col="grey70", typ='l', ylim=c(0, 1), cex.axis=1.4, cex.lab=1.4,
     ylab="", xlab="Preserve radius (m)", bty="l")#101000 m^3 pond
for(i in 2:Ndraws){
  lines(dens.dist$dist.m, prob.ext1[,i,c(1:20), 1], col="grey70")
}
med.prob.101 <- c(1:20)
for(i in 1:20){
  med.prob.101[i] <- median(prob.ext1[,,i,1])
}
lines(dens.dist$dist.m, med.prob.101, col=1, pch=20, cex=3, lwd=4)
mtext("d", 3, .1, at=200, cex=1.5)

dev.off()


#Median extinction probabilities at each preserve radii
olcott.ext1<-blomquist.ext1<-matrix(NA, 20, 1)
for(i in 1:20){
  olcott.ext1[i,] <- median(prob.ext1[,,i,1])
  blomquist.ext1[i,] <- median(prob.ext1[,,i,2])
}

#Number of iterations very low extinction probability (<=2%) at each preserve radii
blom.min1<-olc.min1<-c(1:20)
for(i in 1:20){
  blom.min1[i]<-sum(ifelse(prob.ext1[,,i,2]<=.02,1,0))#98.4% of iterations very likely go extinct if 800 m radii, >99% if less
  olc.min1[i]<-sum(ifelse(prob.ext1[,,i,1]<=.02,1,0))#98.8% of iterations very likely go extinct if 500 m radii or less, 100% if less
}
