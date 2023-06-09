# Analyses for Caruso et al. What are the environmental determinants of phenotypic selection? A meta-analysis of experimental studies. The American Naturalist 

### Main text analyses
#libraries needed: MCMCglmm; dclone; mvtnorm; plyr

# For analysis of selection gradients reporting standard errors described in the main text read in the combined dataset of experimental studies and spatiotemporally replicated studies (the later from Siepielski et al 2009 and taken from Kingsolver and Diamond 2011).

d<-read.csv("determinants of selection dataset_main analyses.csv", header=TRUE)

#subset to manipulative experiments dataset for objectives 1 & 2
data<-subset(d,d$dataset=="experimental")

##fitting of univariate and multivariate mixed models: equations 1 and 2

# assign unique experimental unit identity (i.e. selection gradients for each combination of trait, fitness, treatment within experiments, etc.)
data$id<-paste(data$study.id, data$generic.spp,data$generic.season,data$generic.trait, data$generic.sex,data$generic.age,data$generic.fitness.measure)


# reduce included traits to those with at least 15 estimates of selection. 
data<-subset(data,data$trait.reduced %in% c("LH","MO","PY","S"))

# checking data and factor levels
table(is.na(data$grad.linear.value)) # should all be FALSE
table(is.na(data$grad.linear.se)) # should all be FALSE
table(data$seln.agent.class)
table(data$experimental.context)
table(data$trait.reduced)
table(data$taxon.group)
table(data$fitness.component)

# data cloning 
library(dclone)

# this controls how long all the individual models run
# (useful for testing code, and then setting up for final # runs)
y<-20 # for final runs use 20
clones<-50 # for final runs use 50

data$id<-as.numeric(as.factor(as.character(data$id)))
data$row.id<-as.numeric(as.character(data$row.id))
data$trait.reduced<-as.factor(as.character(data$trait.reduced))

cloned.data<-dclone(data,n.clones=clones)
cloned.data$id <- dclone(dciid(data$id,iid="id"),n.clones=clones)
cloned.data$row.id <- dclone(dciid(data$row.id,iid="row.id"),n.clones=clones)

cloned.data$fitness.component <- as.factor(as.character(cloned.data$fitness.component))
cloned.data$seln.agent.class <- as.factor(as.character(cloned.data$seln.agent.class))
cloned.data$experimental.context <- as.factor(as.character(cloned.data$experimental.context))
cloned.data$trait.reduced <- as.factor(as.character(cloned.data$trait.reduced))
cloned.data$taxon.group <- as.factor(as.character(cloned.data$taxon.group))
cloned.data$id <- as.factor(as.character(cloned.data$id))
cloned.data$row.id <- as.factor(as.character(cloned.data$row.id))

cloned.mev<-cloned.data$grad.linear.se^2

# post cloning checks 
# these should be increased by 50X

length(table(data$id))
length(table(data$row.id))
length(table(cloned.data$id))
length(table(cloned.data$row.id))

# these should have the same number of levels
length(table(data$fitness.component))
length(table(cloned.data$fitness.component))
length(table(data$seln.agent.class))
length(table(cloned.data$seln.agent.class))
length(table(data$experimental.context))
length(table(cloned.data$experimental.context))
length(table(data$trait.reduced))
length(table(cloned.data$trait.reduced))
length(table(data$taxon.group))
length(table(cloned.data$taxon.group))

## 1. Full model with different modifier variables simultaneously: equation 2
library(MCMCglmm)

p1<-list(G=list( G1=list(V=0.001,nu=1),
                G2=list(V=diag(5)*0.001,nu=1),
                G3=list(V=diag(4)*0.001,nu=1),
                G4=list(V=diag(3)*0.001,nu=1),
                G5=list(V=diag(3)*0.001,nu=1) ),
        R=list(V=diag(4)*0.001,nu=1))
        
system.time(
  m1<-MCMCglmm(grad.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,    random=~id
     + idh(fitness.component):row.id
     + idh(seln.agent.class):row.id
     + idh(experimental.context):row.id
     + idh(taxon.group):row.id,
  rcov=~idh(trait.reduced):units,
  mev=cloned.mev,
  prior=p1,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)


## 2 Individual models for the modifier variables: equation 1

# fitness component
p2.1<-list(G=list( G1=list(V=0.001,nu=1)),
        R=list(V=diag(5)*0.001,nu=1))
system.time(
m.fit.comp.50clones<-MCMCglmm(grad.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,
  random=~id,
  rcov=~idh(fitness.component):units,
  mev=cloned.mev,
  prior=p2.1,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

# selective agent
p2.2<-list(G=list( G1=list(V=0.001,nu=1)),
        R=list(V=diag(4)*0.001,nu=1))
system.time(
m.sel.agent.50clones<-MCMCglmm(grad.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,
  random=~id,
  rcov=~idh(seln.agent.class):units,
  mev=cloned.mev,
  prior=p2.2,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

# experimental context
p2.3<-list(G=list( G1=list(V=0.001,nu=1)),
        R=list(V=diag(3)*0.001,nu=1))
system.time(
m.exp.context.50clones<-MCMCglmm(grad.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,
  random=~id,
  rcov=~idh(experimental.context):units,
  mev=cloned.mev,
  prior=p2.3,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

# trait type
p2.4<-list(G=list( G1=list(V=0.001,nu=1)),
        R=list(V=diag(4)*0.001,nu=1))
system.time(
m.trait.50clones<-MCMCglmm(grad.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,
  random=~id,
  rcov=~idh(trait.reduced):units,
  mev=cloned.mev,
  prior=p2.4,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

# taxa type
p2.5<-list(G=list( G1=list(V=0.001,nu=1)),
        R=list(V=diag(3)*0.001,nu=1))
system.time(
m.taxon.50clones<-MCMCglmm(grad.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,
  random=~id,
  rcov=~idh(taxon.group):units,
  mev=cloned.mev,
  prior=p2.5,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

##3. testing whether the magnitude of selection varies with mean fitness: equation 3

#needed packages
library(plyr)
library(MCMCglmm)
library(dclone)

# read in the combined dataset of experimental studies and studies from Siepielski 2009
d<-read.csv("determinants of selection dataset_main analyses.csv", header=TRUE)

# subset to manipulative experiments dataset
data<-subset(d,d$dataset=="experimental")

# assign unique experimental unit identity (i.e. selection gradients for each combination of trait, fitness, treatment within experiments, etc.)
data$id<-paste(data$study.id, data$generic.spp,data$generic.season,data$generic.trait, data$generic.sex,data$generic.age,data$generic.fitness.measure)

data$id<-as.numeric(as.factor(as.character(data$id)))
data$row.id<-as.numeric(as.character(data$row.id))
data$trait.reduced<-as.factor(as.character(data$trait.reduced))

# reduce included traits to those with at least 15 estimates of selection. 
data<-subset(data,data$trait.reduced %in% c("LH","MO","PY","S"))

# remove estimates with missing fitness values
data2<-data[!is.na(data$fitness.mean.num),]

# checking data
table(is.na(data2$grad.linear.value)) # should all be FALSE
table(is.na(data2$grad.linear.se)) # should all be FALSE
table(is.na(data2$fitness.mean.num)) # should all be FALSE

# calculate mean fitness, standardized among treatments within experiments
mean.fit<-tapply(data2$fitness.mean.num, data2$id, mean)
d3<-adply(mean.fit, c(1))
colnames(d3)<-c("id","mean.fit")
data2<-merge(data2,d3, by="id")
data2$standardFitness<-data2$fitness.mean.num/data2$mean.fit
data2$standardFitness<-data2$standardFitness-1
data2$id<-as.numeric(as.factor(as.character(data2$id)))

# this controls how long all the individual models run
# (useful for testing code, and then setting up for final # runs)
y<-20 # for final runs use 20
clones<-50 # for final runs use 50

# clone data
cloned.data2<-dclone(data2,n.clones=clones)
cloned.data2$id <- dclone(dciid(data2$id,iid="id"),n.clones=clones)
cloned.data2$id <- as.factor(cloned.data2$id)
cloned.mev2<-cloned.data2$grad.linear.se^2

#run model
p3<-list(R=list(V=0.001,nu=0.001),G=list(G1=list(V=diag(2)*0.001,nu=0.002)))
system.time(
  m.randReg<-MCMCglmm(grad.linear.value ~1+standardFitness,
      random=~us(1+standardFitness):id,
      data=cloned.data2,mev=cloned.mev2,
      prior=p3,
      verbose=FALSE,
      nitt=13000*y,thin=10*y,burnin=3000*y)
)

# 4 Comparing selection in manipulative experiments to selection in the wild: equation 4

#read the original data. Note that we are using the same data file here as above. For analyses up to now, we just subsetted out the spatiotemporal studies, now we want to retain them. 
d<-read.csv("determinants of selection dataset_main analyses.csv",header=TRUE,sep=',')

# assign unique experimental unit identity (i.e. selection grdients for each combination of trait, fitness, treatment within experiments as nescent.id

ds<-d[d$dataset=="spatiotemporal", ]
ds$id<-paste(ds$study.id, ds$generic.spp, 
 ds$generic.season, ds$generic.trait, ds$generic.sex, ds$generic.age, ds$generic.fitness.measure)
 
de<-d[d$dataset=="experimental", ] 
de$id<-paste(de$study.id, de$generic.year, de$generic.spp, 
 de$generic.season, de$generic.population, de$generic.trait, de$generic.sex, de$generic.age, de$generic.fitness.measure)

d3<-merge(de, ds, all=TRUE)

# reduce included traits to those with at least 15 rows in either dataset
d3<-subset(d3, d3$trait.reduced %in% c("LH","MO","PY","S"))

# format row and experimental unit unique id's and trait factor levels
d3$id<-as.numeric(as.factor(as.character(d3$id)))
d3$row.id<-as.numeric(as.character(d3$row.id))
d3$trait.reduced<-as.factor(as.character(d3$trait.reduced))

# this controls how long all the individual models run
# (useful for testing code, and then setting up for final # runs)
y<-20 # for final runs use 20
clones<-50 # for final runs use 50

# clone data
cloned.d3<-dclone(d3,n.clones=clones)
cloned.d3$id <- dclone(dciid(d3$id,iid="id"),n.clones=clones)
cloned.d3$row.id <- dclone(dciid(d3$row.id,iid="row.id"),n.clones=clones)
cloned.d3$fitness.component <- as.factor(as.character(cloned.d3$fitness.component))
cloned.d3$seln.agent.class <- as.factor(as.character(cloned.d3$seln.agent.class))
cloned.d3$experimental.context <- as.factor(as.character(cloned.d3$experimental.context))
cloned.d3$trait.reduced <- as.factor(as.character(cloned.d3$trait.reduced))
cloned.d3$taxon.group <- as.factor(as.character(cloned.d3$taxon.group))
cloned.d3$id <- as.factor(as.character(cloned.d3$id))
cloned.d3$row.id <- as.factor(as.character(cloned.d3$row.id))
cloned.d3$dataset <- as.factor(as.character(cloned.d3$dataset))
cloned.mev3<-cloned.d3$grad.linear.se^2

#run model
p4<-list(R=list(V=diag(2)*0.001,nu=0.002),G=list(G1=list(V=diag(2)*0.001,nu=0.002)))
system.time(
  mainDatasetComparisonModel<-MCMCglmm(grad.linear.value ~ dataset-1,
   random=~idh(dataset):id,
   rcov=~idh(dataset):units,
   mev=cloned.mev3,data= cloned.d3,
      verbose=FALSE,
   burnin=3000*y,nitt=13000*y,thin=10*y,prior=p4)
)
###################
###################

###for analysis of all selection gradients including those lacking estimates of standard error, described in Appendix D read in the dataset directly below. Use the script above and omit the "mev" argument from mcmcglmm models.
d<-read.csv("determinants of selection dataset.all.estimates_App D.csv",header=TRUE,sep=',')

###################
###################

### For analysis of selection differentials reporting standard error described in the Appendix A

#libraries needed 
library(MCMCglmm)
library(dclone)
library(mvtnorm)
library(plyr)
library(dclone)

#read in the dataset of selection differentials reporting standard error from experimental studies 
d<-read.csv("determinants of selection dataset.differentials_App A.csv")

#remove differentials lacking standard error estimates
table(is.na(d$diff.linear.sterr))
d<-d[!is.na(d$diff.linear.sterr), ]
table(is.na(d$diff.linear.sterr))

# data cloning 
library(dclone)

# this controls how long all the individual models run
# (useful for testing code, and then setting up for final # runs)
y<-20 # for final runs use 20
clones<-50 # for final runs use 50

# assign unique experimental unit identity (i.e. selection gradients for each combination of trait, fitness, treatment within experiments as nescent.id)
d$id<-paste(d$study.id, d$generic.spp,d$generic.season,d$generic.trait, d$generic.sex,d$generic.age,d$generic.fitness.measure)

d$id<-as.numeric(as.factor(as.character(d$id)))
d$row.id<-as.numeric(as.character(d$row.id))
d$trait.reduced<-as.factor(as.character(d$trait.reduced))

cloned.data<-dclone(d,n.clones=clones)
cloned.data$id <- dclone(dciid(d$id,iid="id"),n.clones=clones)
cloned.data$row.id <- dclone(dciid(d$row.id,iid="row.id"),n.clones=clones)

cloned.data$fitness.component <- as.factor(as.character(cloned.data$fitness.component))
cloned.data$seln.agent.class <- as.factor(as.character(cloned.data$seln.agent.class))
cloned.data$experimental.context <- as.factor(as.character(cloned.data$experimental.context))
cloned.data$trait.reduced <- as.factor(as.character(cloned.data$trait.reduced))
cloned.data$taxon.group <- as.factor(as.character(cloned.data$taxon.group))
cloned.data$id <- as.factor(as.character(cloned.data$id))
cloned.data$row.id <- as.factor(as.character(cloned.data$row.id))

cloned.mev<-cloned.data$diff.linear.sterr^2

# post cloning checks 
# these should be increased by 50X

length(table(d$id))
length(table(d$row.id))
length(table(cloned.data$id))
length(table(cloned.data$row.id))

# these should have the same number of levels
length(table(d$fitness.component))
length(table(cloned.data$fitness.component))
length(table(d$seln.agent.class))
length(table(cloned.data$seln.agent.class))
length(table(d$experimental.context))
length(table(cloned.data$experimental.context))
length(table(d$trait.class.reduced))
length(table(cloned.data$trait.reduced))
length(table(d$taxon.group))
length(table(cloned.data$taxon.group))

# fitting of univariate and multivariate mixed models: eqs 1 and eq 2
#needed packages
library(MCMCglmm)

# 1. Full model with different modifier variables simultaneously

p1<-list(G=list( G1=list(V=0.001,nu=1),
                G2=list(V=diag(4)*0.001,nu=1),
                G3=list(V=diag(4)*0.001,nu=1),
                G4=list(V=diag(3)*0.001,nu=1),
                G5=list(V=diag(3)*0.001,nu=1) ),
        R=list(V=diag(4)*0.001,nu=1))
        
system.time(
  m1<-MCMCglmm(diff.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,    random=~id
     + idh(fitness.component):row.id
     + idh(seln.agent.class):row.id
     + idh(experimental.context):row.id
     + idh(taxon.group):row.id,
  rcov=~idh(trait.reduced):units,
  mev=cloned.mev,
  prior=p1,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

# 2. Individual models for the modifier variables

# fitness component
p2.1<-list(G=list( G1=list(V=0.001,nu=1)),
        R=list(V=diag(4)*0.001,nu=1))
system.time(
m.fit.comp.50clones<-MCMCglmm(diff.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,
  random=~id,
  rcov=~idh(fitness.component):units,
  mev=cloned.mev,
  prior=p2.1,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

# selective agent
p2.2<-list(G=list( G1=list(V=0.001,nu=1)),
        R=list(V=diag(4)*0.001,nu=1))
system.time(
m.sel.agent.50clones<-MCMCglmm(diff.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,
  random=~id,
  rcov=~idh(seln.agent.class):units,
  mev=cloned.mev,
  prior=p2.2,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

# experimental context
p2.3<-list(G=list( G1=list(V=0.001,nu=1)),
        R=list(V=diag(3)*0.001,nu=1))
system.time(
m.exp.context.50clones<-MCMCglmm(diff.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,
  random=~id,
  rcov=~idh(experimental.context):units,
  mev=cloned.mev,
  prior=p2.3,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

# trait type
p2.4<-list(G=list( G1=list(V=0.001,nu=1)),
        R=list(V=diag(4)*0.001,nu=1))
system.time(
m.trait.50clones<-MCMCglmm(diff.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,
  random=~id,
  rcov=~idh(trait.class.reduced):units,
  mev=cloned.mev,
  prior=p2.4,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

# taxa type
p2.5<-list(G=list( G1=list(V=0.001,nu=1)),
        R=list(V=diag(3)*0.001,nu=1))
system.time(
m.taxon.50clones<-MCMCglmm(diff.linear.value~fitness.component+
     seln.agent.class + experimental.context+trait.reduced+
     taxon.group,
  random=~id,
  rcov=~idh(taxon.group):units,
  mev=cloned.mev,
  prior=p2.5,
  verbose=FALSE,
  nitt=13000*y,thin=10*y,burnin=3000*y,data= cloned.data)
)

###################
###################
###################
# 3. testing whether the magnitude of selection varies with mean fitness

#needed packages
library(plyr)
library(MCMCglmm)
library(dclone)

#read in the dataset of selection differentials reporting standard error from experimental studies 

d<-read.csv("determinants of selection dataset.differentials_App A.csv")

#remove differentials lacking standard error estimates
table(is.na(d$diff.linear.sterr))
d<-d[!is.na(d$diff.linear.sterr), ]
table(is.na(d$diff.linear.sterr))

# remove estimates with missing fitness values
d2<-d[!is.na(d$fitness.mean.num),]

# calculate mean fitness, standardized among treatments within experiments
mean.fit<-tapply(d2$fitness.mean.num, d2$id, mean)
d3<-adply(mean.fit, c(1))
colnames(d3)<-c("id","mean.fit")
d2<-merge(d2,d3, by="id")
d2$standardFitness<-d2$fitness.mean.num/d2$mean.fit
d2$standardFitness<-d2$standardFitness-1
d2$id<-as.numeric(as.factor(as.character(d2$id)))

# this controls how long all the individual models run
# (useful for testing code, and then setting up for final # runs)
y<-20 # for final runs use 20
clones<-50 # for final runs use 50

# clone data
cloned.d2<-dclone(d2,n.clones=clones)
cloned.d2$id <- dclone(dciid(d2$id,iid="id"),n.clones=clones)
cloned.d2$id <- as.factor(cloned.d2$id)
cloned.mev2<-cloned.d2$diff.linear.sterr^2

#run model
p3<-list(R=list(V=0.001,nu=0.001),G=list(G1=list(V=diag(2)*0.001,nu=0.002)))
system.time(
  mdiff.randReg<-MCMCglmm(diff.linear.value ~1+standardFitness,
      random=~us(1+standardFitness):id,
      data=cloned.d2,mev=cloned.mev2,
      prior=p3,
      verbose=FALSE,
      nitt=13000*y,thin=10*y,burnin=3000*y)
)