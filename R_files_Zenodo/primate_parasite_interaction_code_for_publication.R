library(dplyr);library(igraph);library(bipartite);
library(sna);library(network);library(ergm);library(ggplot2);library(bipartite)


#########################
##########  USE THIS FILE TO CALCULATE BIPARTITE NETS
###########################
byhostpar<-read.csv(file="Supplementary Appendix GMPD data feb 25 2019 Host-Parasite interactions after removing duplicated parasite sp.csv",sep=",",header=T)


##first just use the acast function with the aggregating argument being the mean
# to do that, we want to replace NA prevalence, and for the purposes of creating a pres / abs matrix, let's just insert a value of 0.07, the median, because we'll convert to 1 or 0 anyway, but if we want to use the prevalence info, would have to exclude data with no prevalence


##looks like some of the prevalence values have -99 for NA
which(colnames(byhostpar2)=="Prevalence")
byhostpar2[which(byhostpar2$Prevalence<0),9]<-NA
#View(byhostpar2)

###replace missing with mean, just for the purposes of keeping parasites for which prevalence is missing
colnames(byhostpar2)
byhostpar2[is.na(byhostpar2$Prevalence),9]<-0.07
#View(byhostpar2)
str(byhostpar2)
summary(byhostpar2$Prevalence)

##make the matrix using acast, fill in missing values with -99, use the prevalence column, and use the mean function
library(reshape2)
bipart2<-as.data.frame(acast(byhostpar2, 
                             Host.Corrected.Name.W.R05~Parasite.Corrected.Name, 
                           value.var='Prevalence',fill=-99,
                          fun.aggregate=mean, margins=F))
str(bipart2)
#View(bipart2)
dim(bipart2)
write.table(bipart2,file="bipart2.txt",sep="\\t")


pardattab<-as.matrix(bipart2)
pardattab[1:5,1:5]
dim(pardattab)
#View(pardattab)

write.table(bipart2,file="bipartitematrixwnasandmeanprev.txt",sep="\\t")

bipart2[bipart2==-99]<-NA

#convert nonzero prevalences to 0/1
bipart3[bipart3>0]<-1

##how many hosts per parasite?
hostsperparasite<-colSums(bipart3,na.rm=T)
#how many parasites per host?
parasitespperhost<-rowSums(bipart3,na.rm=T)

str(hostsperparasite)
summary(hostsperparasite)
summary(parasitespperhost)
hist(parasitespperhost)

dim(bipart3)

##remove any parasites/hosts with all 0s

bipart4<-bipart3
bipart4[,which(colSums(bipart4,na.rm=T)==0)]<-NULL

##omitted the species with all 0s for parasites - mostly NWM tested negative for plasmodium or OWM neg for STLV; also Hapa sp.

bipart4<-bipart4[-which(rowSums(bipart4,na.rm=T)==0),]

dim(bipart4)

#View(bipart3)

#double check, how many hosts per parasite?
hostsperparasite<-colSums(bipart4,na.rm=T)
hostsperparasite$parasitesp<-rownames(hostsperparasite)
#how many parasites per host?
parasitespperhost<-rowSums(bipart4,na.rm=T)
str(hostsperparasite)
str(parasitespperhost)
summary(hostsperparasite)
summary(parasitespperhost)



##convert NAs to 0s
bipart_bin<-bipart4
bipart_bin[is.na(bipart_bin)]<-0

write.table(bipart_bin,file="bipartmat_bin.txt",sep="\\t")

################################
################ FINAL BIPARTITE NETWORK BINARY
########### THIS FILE ALSO HAS REMOVED ALL ERRONEOUS PARASITES AFTER DOWNSTREAM CLEANING

bipart_bin<-read.csv(file="Supplementary Appendix Host-Parasite Matrix.csv",sep=",",header=T,row.names = 1)
dim(bipart_bin)

bipart_bin[1:5,1:5]
rownames(bipart_bin)<-bipart_bin$host
bipart_bin$host<-NULL
bipart_bin<-as.matrix(bipart_bin)


#################################################################
################ network analysis
##make bipartite graph

net<-network(bipart_bin,bipartite=T,matrix.type="bipartite")
net


net$val[1:213]
net$val[1:780]
plot(net)
##############################
#set the vertex attributes
#################################

#hosts first
hostnames<-as.data.frame(network.vertex.names(net)[1:213])
colnames(hostnames)<-"host"
parnames<-as.data.frame(network.vertex.names(net)[214:976])
colnames(parnames)<-"parasite"
allvertnames<-as.data.frame(network.vertex.names(net))
colnames(allvertnames)<-"species"


#################
###  assign the network vertex attributes
##basics
vert_dat<-read.csv(file="Supplementary Table all vertex traits.csv",sep=",",header=T)
str(vert_dat)
vertnames<-as.data.frame(net %v% 'vertex.names')
str(vertnames)

#verify network and traits match
results1 <- setdiff(vert_dat$species, vertnames$`net %v% "vertex.names"`)
results1 <- setdiff(vertnames$`net %v% "vertex.names"`,vert_dat$species)
d<-vert_dat[which(duplicated(vert_dat$species)),1]


type<-as.character(vert_dat$type)
names(type)<-vert_dat$species
length(type)
type<-type[ net %v% 'vertex.names' ]
network::set.vertex.attribute(net,"type",type)
network::get.vertex.attribute(net,attrname="type")


##primate tax group
taxgroup<-as.character(vert_dat$taxonomic.group)
names(taxgroup)<-vert_dat$species
length(taxgroup)
taxgroup<-taxgroup[ net %v% 'vertex.names' ]
network::set.vertex.attribute(net,"taxgroup",taxgroup)
network::get.vertex.attribute(net,attrname="taxgroup")


##mass
mass<-scale(log(vert_dat$mass_sex_mean+1))
summary(mass)
names(mass)<-vert_dat$species
length(mass)
mass<-mass[ net %v% 'vertex.names' ]
network::set.vertex.attribute(net,"mass",mass)
network::get.vertex.attribute(net,attrname="mass")

#sample effort, N
N<-scale(vert_dat$N)

summary(N)
names(N)<-vert_dat$species
length(N)
N<-N[ net %v% 'vertex.names' ]

network::set.vertex.attribute(net,"N",N)
network::get.vertex.attribute(net,attrname="N")

head(network::get.vertex.attribute(net,attrname="N"))

##geographic

#area
area<-scale(log(vert_dat$area+1))
names(area)<-vert_dat$species
length(area)
area<-area[ net %v% 'vertex.names' ]

network::set.vertex.attribute(net,"area",area)
head(network::get.vertex.attribute(net,"area"))

hosttraits<-as.data.frame(cbind(N,mass,area))
str(hosttraits)
write.table(hosttraits,file="hosttraits.txt",sep="\\t")

#climate pcs

pc1<-vert_dat$PC1
names(pc1)<-vert_dat$species
length(pc1)
pc1<-pc1[ net %v% 'vertex.names' ]
network::set.vertex.attribute(net,"clim.pc1",pc1)
network::get.vertex.attribute(net,"clim.pc1")

pc2<-vert_dat$PC2
names(pc2)<-vert_dat$species
length(pc2)
pc2<-pc2[ net %v% 'vertex.names' ]
network::set.vertex.attribute(net,"clim.pc2",pc2)
network::get.vertex.attribute(net,"clim.pc2")

pc3<-scale(vert_dat$PC3)
names(pc3)<-vert_dat$species
length(pc3)
pc3<-pc3[ net %v% 'vertex.names' ]
network::set.vertex.attribute(net,"clim.pc3",pc3)

network::get.vertex.attribute(net,"clim.pc3")

clim3<-as.numeric(scale(network::get.vertex.attribute(net,"clim.pc3")))

network::set.vertex.attribute(net,"clim.pc3",clim3)

network::get.vertex.attribute(net,"clim.pc3")


sprich<-vert_dat$sprich
names(sprich)<-vert_dat$species
length(sprich)
sprich<-sprich[ net %v% 'vertex.names' ]
network::set.vertex.attribute(net,"sprich",sprich)
network::get.vertex.attribute(net,"sprich")


##threat status
Threat<-vert_dat$threatened
summary(Threat)
Threat[1:213]
names(Threat)<-vert_dat$species
length(Threat)
Threat<-Threat[ net %v% 'vertex.names' ]
head(Threat)
network::set.vertex.attribute(net,"Threat",Threat)
network::get.vertex.attribute(net,"Threat")

##biogeo data, created below
biogeo<-vert_dat2$Continent
biogeo<-biogeo[ net %v% 'vertex.names' ]
network::set.vertex.attribute(net,"biogeo",as.character(biogeo))
head(network::get.vertex.attribute(net,"biogeo"))


##################################################
#parasite traits
##parasite types

##make sure parasite names match order of network
vert_dat<-vert_dat[order(match(vert_dat$species,net %v% 'vertex.names')),]


network::set.vertex.attribute(net,"par.N",as.numeric(vert_dat$N))
summary(network::get.vertex.attribute(net,"par.N"))


####parasite types
Parasite.Type<-as.character(vert_dat$Parasite.Type)
names(Parasite.Type)<-vert_dat$species
Parasite.Type<-Parasite.Type[ net %v% 'vertex.names' ]
str(Parasite.Type)

network::set.vertex.attribute(net,"Parasite.Type",Parasite.Type)
network::get.vertex.attribute(net,"Parasite.Type")


###add trans modes
Close<-as.character(vert_dat$Close)
Close
names(Close)<-vert_dat$species
length(Close)
Close<-Close[ net %v% 'vertex.names' ]

network::set.vertex.attribute(net,"Close",Close)
network::get.vertex.attribute(net,"Close")

Nonclose<-as.character(vert_dat$Nonclose)
Nonclose
names(Nonclose)<-vert_dat$species
length(Nonclose)
Nonclose<-Nonclose[ net %v% 'vertex.names' ]

network::set.vertex.attribute(net,"Nonclose",Nonclose)
network::get.vertex.attribute(net,"Nonclose")

Sexual<-as.character(vert_dat$Sexual)
Sexual
names(Sexual)<-vert_dat$species
length(Sexual)
Sexual<-Sexual[ net %v% 'vertex.names' ]

network::set.vertex.attribute(net,"Sexual",Sexual)
network::get.vertex.attribute(net,"Sexual")


Vector<-as.character(vert_dat$Vector)
Vector
names(Vector)<-vert_dat$species
length(Vector)
Vector<-Vector[ net %v% 'vertex.names' ]

network::set.vertex.attribute(net,"Vector",Vector)
network::get.vertex.attribute(net,"Vector")

Intermediate<-as.character(vert_dat$Intermediate)
Intermediate
names(Intermediate)<-vert_dat$species
length(Intermediate)
Intermediate<-Intermediate[ net %v% 'vertex.names' ]

network::set.vertex.attribute(net,"Intermediate",Intermediate)
network::get.vertex.attribute(net,"Intermediate")


#########ADD IN the parasite trait about whether they are found in other mammals or not
other<-vert_dat$other.mammals
names(other)<-net %v% 'vertex.names' 
length(other)
View(other)
network::set.vertex.attribute(net,"other.mammals",other)
network::get.vertex.attribute(net,"other.mammals")



############ADD PHYLOGENY
library(ape);library(geiger)
tree<-read.tree("springer.tre")

phylodist<-cophenetic.phylo(tree)
dim(phylodist)
summary(phylodist)
class(phylodist)
plot(tree4)
axisPhylo()

############################################
###########   ergms
###############################################3

##basic model with just edges, ie. intercept
#test if observed network has more or less edges than expected by random graph
Fm1<-ergm ( net ~ edges)
summary(Fm1)##neg effect of edges, meaning there are fewer edges than expected by chance - accounts for sparse matrix
fm1gof<-gof(Fm1)


##for other terms you can use with bipartite graphs, see:
search.ergmTerms(categories=c("bipartite"))
search.ergmTerms()

######testing host traits
#test if there are more edges for same tax group of hosts
Fm4<-ergm ( net ~ edges+b1nodematch("taxgroup"))
sumfm4<-summary(Fm4)##positive effect, nodes in same tax group have more edges to same parasites than nodes in diff tax groups

mcmc.diagnostics(Fm4)
fm4gof<-gof(Fm4) 

#test if there are more edges for same tax group of hosts and if some tax groups have more edges than others
Fm5<-ergm ( net ~ edges+b1nodematch("taxgroup")+b1factor("taxgroup"))
sumfm5<-summary(Fm5)#significant differences among tax groups - lemurs and lorises have fewer edges than expected comapred to catarhines, platys and tars not different


#test if mass covaries
Fm6<-ergm ( net ~ edges+b1cov("mass"))
summaryfm6<-summary(Fm6)#significant positive effect of mass - more edges with higher body mass

fm6gof<-gof(Fm6)
plot(fm6gof)

#test if mass and citations covaries,
Fm7<-ergm ( net ~ edges+b1cov("mass")+b1cov("N"))
summaryfm7<-summary(Fm7)#significant effects of mass and N - more edges for higher body mass and more citations

fm7gof<-gof(Fm7)
plot(fm7gof)


#test if mass and citations covaries, with taxgroup too
##error bcause of low edge number
Fm8<-ergm ( net ~ edges+b1cov("mass")+b1cov("N")+b1nodematch("taxgroup"))
summaryfm8<-summary(Fm8)#
fm8gof<-gof(Fm8)
plot(fm8gof)
mcmc.diagnostics(Fm8)

#test if area covaries,
Fm9<-ergm ( net ~ edges+b1cov("area"))
summary(Fm9)#significant effects of area - more edges for larger geo areas

#test if area, mass, and citations covaries,
Fm10<-ergm ( net ~ edges+b1cov("area")+b1cov("mass")+b1cov("N"))
fm10summary<-summary(Fm10)#significant effects of area, mass, and N - more edges for larger geo areas, larger body mass, and more samples
fm10gof<-gof(Fm10)
str(fm10gof)
plot(fm10gof)
fm10sims<-simulate(Fm10,nsim=100)
summary(fm10sims)


#test if area, mass, diet and citations covaries, and include tax group as node factor

Fm11<-ergm ( net ~ edges+b1cov("area")+b1cov("sprich")+
               b1cov("mass")+b1cov("N")+b1nodematch("diet")+
               b1nodematch("taxgroup"))
summaryfm11<-summary(Fm11)#significant effects of area, mass, and N - more edges for larger geo areas, larger body mass, and more samples
write.table(summaryfm11$coefs,file="hosttraitsres.txt",sep="\\t")

fm11gof<-gof(Fm11)
str(fm11gof)
plot(fm11gof)

##with nodefactor taxgroup instead of nodematch
Fm11.1<-ergm ( net ~ edges+b1cov("area")+b1cov("sprich")+
               b1cov("mass")+b1cov("N")+b1cov("clim.pc3")+
               b1factor("taxgroup"))
summaryfm11.1<-summary(Fm11.1)#significant effects of area, mass, and N - more edges for larger geo areas, larger body mass, and more samples
write.table(summaryfm11.1$coefs,file="hosttraitsres.txt",sep="\\t")

fm11.1gof<-gof(Fm11.1)
str(fm11gof)
plot(fm11gof)


#test if climate+area, mass, and citations covaries, and include tax group as node factor
net
Fmclim<-ergm ( net ~ edges+
                 b1cov("clim.pc3")+
                 b1cov("area")+b1cov("sprich")+
               b1cov("mass")+b1cov("N")+
               b1factor("taxgroup")+b1factor("Threat"))
summaryfmclim<-summary(Fmclim)#significant effects of area, mass, and N - more edges for larger geo areas, larger body mass, and more samples
fmclimgof<-gof(Fmclim)
str(fm11gof)
plot(fmclimgof)

View(summaryfmclim$coefs)

###w nodematch tax group
Fmclim2<-ergm ( net ~ edges+
                 b1cov("clim.pc3")+
                 b1cov("area")+b1cov("sprich")+
                 b1cov("mass")+b1cov("N")+
                 b1nodematch("taxgroup")+b1factor("Threat"))
summaryfmclim2<-summary(Fmclim2)#significant effects of area, mass, and N - more edges for larger geo areas, larger body mass, and more samples
###aic = 21738
fmclimgof2<-gof(Fmclim2)

View(summaryfmclim2$coefficients)

p.value <- 1 - pchisq(225890-21283, 162945-162937)


###w nodematch tax group
Fmclim3<-ergm ( net ~ edges+
                  b1cov("clim.pc3")+
                  b1cov("area")+b1cov("sprich")+
                  b1cov("mass")+b1cov("N")+
                  b1factor("taxgroup")+b1factor("Threat"))
summaryfmclim3<-summary(Fmclim3)#significant effects of area, mass, and N - more edges for larger geo areas, larger body mass, and more samples
###aic = 22847
fmclimgof3<-gof(Fmclim3)


############################################################################
########ERGMS testing parasite traits
##################################################################3
#test if there are more edges for same parasite types
Fm3<-ergm ( net ~ edges+
              b2nodematch("Parasite.Type"))
#aic = 24044
summaryfm3<-summary(Fm3)##
mcmcres3<-mcmc.diagnostics(Fm3)
fm3gof<-gof(Fm3)


#test if there are more edges for generalist/specialist
Fm3<-ergm ( net ~ edges+
              b2factor("other.mammals"))
#aic = 24044
summaryfm3<-summary(Fm3)##
mcmcres3<-mcmc.diagnostics(Fm3)
fm3gof<-gof(Fm3)

##################################################
#### test host and parasite traits
#####################################################

##combined host parasite traits
#test if area, mass, and citations covaries, and include tax group as node match
Fm12<-ergm ( net ~ edges+b1cov("clim.pc3")+
               b1cov("area")+
               b1nodematch("taxgroup")+
               b1factor("Threat")+
               b1cov("mass")+b1cov("N")+
               b2factor("Parasite.Type")+b2cov("par.N"))

summaryfm12<-summary(Fm12)#significant effects of area, mass, and N - more edges for larger geo areas, larger body mass, and more samples
write.table(summaryfm12$coefs,file="model results_23_6_standardized.txt",sep="\\t")

fm12gof<-gof(Fm12)
plot(fm12gof)
##convergence issues
mcmcfm12<-mcmc.diagnostics(Fm12)


#test if area, mass, and citations covaries, and include tax group as node match
Fm12.1<-ergm ( net ~ edges+b1cov("clim.pc3")+b1cov("area")+
                 b1nodematch("taxgroup")+b1factor("Threat")+
               b1cov("mass")+b1cov("N")+b2factor("Parasite.Type"))
summaryfm12.1<-summary(Fm12.1)#significant effects of area, mass, and N - more edges for larger geo areas, larger body mass, and more samples
fm12.1gof<-gof(Fm12.1)
plot(fm12.1gof)
write.table(summaryfm12.1$coefs,file="model results.txt",sep="\\t")

#test using nodematches and diet if area, mass, and citations covaries, and include tax group as node match
Fm12.2<-ergm ( net ~ edges+b1cov("clim.pc1")+b1cov("area")+
                 b1nodematch("taxgroup")+b1factor("Threat")+b1nodematch("diet")+
                 b1cov("mass")+b1cov("N")+b2factor("Parasite.Type"))
summaryfm12.2<-summary(Fm12.2)#significant effects of area, mass, and N - more edges for larger geo areas, larger body mass, and more samples
fm12.2gof<-gof(Fm12.2)
plot(fm12.1gof)
write.table(summaryfm12.2$coefs,file="model results.txt",sep="\\t")

###
Fm12.2$newnetwork
newnetmat<-as.matrix.network.adjacency(Fm12.2$newnetworks[[1]])
write.table(newnetmat,file="newnetmat.txt",sep="\\t")

###


#test if area, mass, and citations covaries, and include tax group as node match
Fm12.3<-ergm ( net ~ edges+b1cov("clim.pc3")+
               b1cov("area")+b1mindegree(2)+
               b1nodematch("taxgroup")+
                 b1nodematch("biogeo")+
                 b1factor("Threat")+
               b1cov("mass")+b1cov("N")+
               b2factor("Parasite.Type")+b2factor("other.mammals")+
                 b2mindegree(2)+b2cov("par.N"))
summaryfm12.3<-summary(Fm12.3)


summaryfm12.3

fm12.3gof<-gof(Fm12.3)
pdf("gof.pdf",height=6,width=22)
plot(fm12.3gof)
dev.off()

ci.res<-confint(Fm12.3)

write.table(summaryfm12.3$coefs,file="model results w primate parasite traits_23_2_22.txt",sep="\\t")
write.table(ci.res,file="CIs_23_2_22.txt",sep="\\t")

View(summaryfm12.3$coefs)
View(summaryfm12$coefs)
#compare with host tax group as factor

#test if area, mass, and citations covaries, and include tax group as node match
Fm12.4<-ergm ( net ~ edges+b1cov("clim.pc3")+
                 b1cov("area")+b1mindegree(2)+
                 b1nodematch("taxgroup")+
                 b1nodematch("biogeo")+
                 b1factor("Threat")+
                 b1cov("mass")+b1cov("N")+
                 b2factor("Parasite.Type")+b2factor("other.mammals")+
                 b2mindegree(2)+b2cov("par.N"))
summaryfm12.4<-summary(Fm12.4)



write.table(summaryfm12.4$coefficients,file="res w taxgroupfactor.txt",sep="\\t")
write.table(ci.res,file="confidenceintervals.txt",sep="\\t")

#####compare threat and no threat here

Fm12.4.nothreat<-ergm ( net ~ edges+b1cov("clim.pc3")+
                 b1cov("area")+b1mindegree(2)+
                 b1factor("taxgroup")+
                 b1nodematch("biogeo")+
                 b1cov("mass")+b1cov("N")+
                 b2factor("Parasite.Type")+b2factor("other.mammals")+
                 b2mindegree(2)+b2cov("par.N"))
summaryfm12.4.nothreat<-summary(Fm12.4.nothreat)

#no area
Fm12.4.noarea<-ergm ( net ~ edges+b1cov("clim.pc3")+
                 b1mindegree(2)+
                 b1factor("taxgroup")+
                 b1nodematch("biogeo")+
                 b1factor("Threat")+
                 b1cov("mass")+b1cov("N")+
                 b2factor("Parasite.Type")+b2factor("other.mammals")+
                 b2mindegree(2)+b2cov("par.N"),iterations=50)
summaryfm12.4.noarea<-summary(Fm12.4.noarea)

#compare without mindegree terms
Fm12.4<-ergm ( net ~ edges+b1cov("clim.pc3")+
                 b1cov("area")+
                 b1nodematch("taxgroup")+
                 b1factor("Threat")+
                 b1cov("mass")+b1cov("N")+
                 b2factor("Parasite.Type")+
                 b2cov("par.N"))
summaryfm12.4<-summary(Fm12.4)#significant effects of area, mass, and N - more edges for larger geo areas, larger body mass, and more samples

summaryfm12.4

Fm12.4$coef

mcmcfm12<-mcmc.diagnostics(Fm12.3)



###############
## unipartite

#####################################################
###################       test if phylogeny covaries
phynet<-network(phylodist)
#plot(phynet)

### do a unipartite projection for hosts and test if host phylo predicts unipartite project
##convert to igraph network
g <- igraph::graph.incidence(bipart_bin)
##do projection
projection<-bipartite_projection(g)
str(projection)
print(projection[[1]])
plot(projection[[1]])
hostmat<-as.matrix(as_adjacency_matrix(projection[[1]]))
hostmat[1:10,1:20]
nethosts<-network(hostmat,directed=F)


dim(hostmat)
hostmat[1:10,1:10]
plot(nethosts)
hostnames<-data.frame(host=network.vertex.names(nethosts))
network.vertex.names(net)[1:213]


 biogeo2<-vert.dat$Continent
  str(biogeo2)
 hostnames<-as.data.frame(hostnames)
 colnames(hostnames)[1]<-"host"
 subset(hostnames, !(host %in% biogeo2$taxon))
# 
# ##create a matrix of whether hosts are on same continent or not
 biogeo2<-biogeo2[order(match(biogeo2[,2],hostnames[,1])),]
 str(biogeo2)
 rownames(biogeo2)<-biogeo2[,2]
 View(biogeo2)
 biogeo2[,1:2]<-NULL
 library(cluster)
 str(biogeo2)
 biogeomat<-data.matrix(daisy(biogeo2,metric="Gower"))
 biogeomat[1:10,1:10]
# 
# ##daisy will give dissimilarity value of 1 for species on diff continents, so let's convert to 0 and 1
# ##kind of silly but it works
 biogeomat[biogeomat==1]<-2
 biogeomat[biogeomat==0]<-1
 biogeomat[biogeomat==2]<-0
 diag(biogeomat)<-1
View(biogeomat)

###proportion of geographic range overlap
rangeoverlap<-read.csv(file="Supplementary Appendix primate range overlap matrix based on IUCN ranges.csv",sep=",",header=T)

dim(rangeoverlap)
View(rangeoverlap)

str(rangeoverlap)
x <- as.matrix(rangeoverlap[,colnames(rangeoverlap) %in% rownames(rangeoverlap)])
rownames(x)
colnames(x)
str(x)
x[1:5,1:5]

rangeoverlapnet<-network(x)


###taking attributes from bipartite network, dropping nas
nethostsmass<-as.data.frame(net%v%"mass")
nethostsmass<-nethostsmass[!is.na(nethostsmass),]
nethostsarea<-as.data.frame(net%v%"area")
nethostsarea<-nethostsarea[!is.na(nethostsarea),]
nethostsn<-as.data.frame(net%v%"N")
nethostsn<-nethostsn[!is.na(nethostsn),]
nethoststax<-as.data.frame(net%v%"taxgroup")
nethoststax<-nethoststax[!is.na(nethoststax),]
length(nethoststax)
nethostssprich<-as.data.frame(net%v%"sprich")
nethostssprich<-nethostssprich[!is.na(nethostssprich),]
nethoststhreat<-as.data.frame(net%v%"Threat")
nethoststhreat<-nethoststhreat[!is.na(nethoststhreat),]

nethostsclim1<-as.data.frame(net%v%"clim.pc1")
nethostsclim1<-nethostsclim1[!is.na(nethostsclim1),]
nethostsclim2<-as.data.frame(net%v%"clim.pc2")
nethostsclim2<-nethostsclim2[!is.na(nethostsclim2),]
nethostsclim3<-as.data.frame(net%v%"clim.pc3")
nethostsclim3<-nethostsclim3[!is.na(nethostsclim3),]

##add on covariates
nethosts%v%"mass"<-nethostsmass
nethosts%v%"area"<-nethostsarea
nethosts%v%"N"<-nethostsn
nethosts%v%"taxgroup"<-as.character(nethoststax)
nethosts%v%"diet"<-as.character(nethostsdiet)
nethosts%v%"sprich"<-nethostssprich
nethosts%v%"diet"<-as.character(nethostsdiet)
nethosts%v%"Threat"<-as.character(nethoststhreat)
nethosts%v%"clim1"<-nethostsclim1
nethosts%v%"clim2"<-nethostsclim2
nethosts%v%"clim3"<-nethostsclim3

taxgr<-nethosts%v%"taxgroup"
unique(taxgr)

#take the inverse of phylodist, so that the matrix is phylo proximity
phylodist<-1/phylodist
diag(phylodist)<-0
#make a network out of the phylogeny
hostphynet<-network(phylodist,directed=F)

##order the phylo matrix to match the host network
# do this: 
library(xergm.common) 
# construct the matrix of t1 values for vertices 
# that appear in t2 and fill in NA for missing vertices 
phylodist2<-adjust(as.matrix(phylodist,attrname = 'value'),nethosts) 
phylodist2 # peek at matrix

phylodist2[1:10,1:10]

biogeomat<-adjust(as.matrix(biogeomat,attrname = 'value'),nethosts) 
biogeomat # peek at matrix

overlapmat<-adjust(as.matrix(x,attrname = 'value'),nethosts) 
overlapmat # peek at matrix

#########################
### ergm

#test just phylogeny
Fm14<-ergm ( nethosts ~ edges+dyadcov(phylodist2))
summary(Fm14)#significant effect
fm14gof<-gof(Fm14)
str(fm14gof)
plot(fm14gof)#really bad fit

##phylogeny with other network properties
Fm14.1<-ergm ( nethosts ~ edges+dyadcov(phylodist2)+kstar(2))
summaryfm14.1<-summary(Fm14.1)#significant effect
fm14.1gof<-gof(Fm14.1)
plot(fm14.1gof)#really bad fit

# phylo + mass, N, area
Fm15<-ergm ( nethosts ~ edges+nodecov("mass")+nodecov("N")+nodecov("area")+dyadcov(phylodist2))
summary(Fm15)#all significant effects and positive
fm15gof<-gof(Fm15)
str(fm15gof)
plot(fm15gof)##pretty aweful fit

##w biogeo
# phylo + mass, N, area, same continent or not
Fm16<-ergm ( nethosts ~ edges+nodecov("mass")+nodecov("N")+nodecov("area")+dyadcov(phylodist2)+dyadcov(biogeomat))
sumfm16<-summary(Fm16)#all significant effects and positive
fm16gof<-gof(Fm16)
str(fm16gof)
plot(fm16gof)##pretty aweful fit


##w range overlap + biogeo
# phylo + mass, N, area, same continent or not
Fm17<-ergm ( nethosts ~ edges+
               nodecov("mass")+
               nodecov("N")+
               nodecov("area")+
               nodecov("clim3")+
               nodematch("taxgroup")+
               dyadcov(phylodist2)+
               dyadcov(biogeomat)+
               dyadcov(overlapmat))
sumfm17<-summary(Fm17)#all significant effects and positive
sumfm17
fm17gof<-gof(Fm17)
str(fm16gof)
plot(fm17gof)##
write.table(sumfm17$coefs,file="unipartite res all full vars.txt",sep="\\t")

# phylo + mass, N, area, same continent or not
Fm18<-ergm ( nethosts ~ edges+
               nodecov("mass")+
               nodecov("N")+
               nodecov("area")+
               nodecov("clim3")+
               nodefactor("Threat")+
               dyadcov(phylodist2)+
               dyadcov(biogeomat)+
               dyadcov(overlapmat))
sumfm18<-summary(Fm18)#all significant effects except body mass, and positive
sumfm18
write.csv(sumfm18$coefs,file="unipartite all vars res.csv")
cisunipart<-confint(Fm18)
write.csv(cisunipart,file="unipartite CIS.csv")

fm18gof<-gof(Fm18)
str(fm16gof)

####plot coefficients


####plot of ergm coefficients from best unipartite model
ergm.coef<-read.csv(file="unipartite all vars res 3_6.csv",sep=",",header=T)
summary(ergm.coef)
View(ergm.coef)
colnames(ergm.coef)

ergm.coef$variable <- factor(ergm.coef$Variable,
                             levels = ergm.coef$Variable)
ergm.coef$variable  # notice the changed order of factor levels
ergm.coef1<-ergm.coef[-1,]
p<-ggplot(ergm.coef1, aes(x=variable, y=Estimate)) + 
  geom_errorbar(aes(ymin=LCI, 
                    ymax=HCI), 
                colour="black", width=0.5) +geom_hline(yintercept=0, linetype="dotted",colour="white")+
  geom_point(size=2,colour="black")+geom_hline(yintercept=0.01, linetype="dotted")+
  xlab("Node- and Dyad-level predictors")+ylab("Coefficient +/- 95% CI")

p+coord_flip()+
  theme(axis.text=element_text(size=60),
        axis.title=element_text(size=65,face="bold"))+theme_classic()

#p+theme_black()

ergm.coef$logodds<-plogis(abs(ergm.coef$coefficient))
p<-ggplot(ergm.coef, aes(x=variable, y=logodds)) + 
  geom_bar(aes(logodds))+
  theme_black() 
p  

