#########################################################################################
# This script provides a sample of code used to perform phylogenetic generalized linear 
# regression analyses testing the relationships among our phenotypic (PC axes 1-3 of body
# shape variation and body size) and ecological (buoyancy, depth, and PC axes 1-3 of diet
# variation) variables within icefishes and notoperches.
#########################################################################################

# load necessary libraries
library(ape)
library(nlme)
library(geiger)
library(phytools)

#read in trait data
dat<-read.csv("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/ice-trem_datafile.csv")

#read in phylogeny
noto.tree<-read.tree("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/notothenioid_timetree.tre")

############################################
# First complete analyses for icefishes
############################################

##isolate traits for icefishes
c.dat<-dat[1:14,]

#create unique variables for each of the individual traits
c.buoy<-as.matrix(c.dat$Buoyancy)
row.names(c.buoy)<-c.dat$Tree_label

c.depth<-as.matrix(c.dat$MeanDepth)
row.names(c.depth)<-c.dat$Tree_label

c.size<-as.matrix(c.dat$MaxSize)
row.names(c.size)<-c.dat$Tree_label

c.shapePC1<-as.matrix(c.dat$BodyShapePC1)
row.names(c.shapePC1)<-c.dat$Tree_label

c.shapePC2<-as.matrix(c.dat$BodyShapePC2)
row.names(c.shapePC2)<-c.dat$Tree_label

c.shapePC3<-as.matrix(c.dat$BodyShapePC3)
row.names(c.shapePC3)<-c.dat$Tree_label

c.dietPC1<-as.matrix(c.dat$DietPC1)
row.names(c.dietPC1)<-c.dat$Tree_label

c.dietPC2<-as.matrix(c.dat$DietPC2)
row.names(c.dietPC2)<-c.dat$Tree_label

c.dietPC3<-as.matrix(c.dat$DietPC3)
row.names(c.dietPC3)<-c.dat$Tree_label

###starting the phylogenetic GLS

##what is the relationship between buoyancy and each of the first three PC axes of body shape variation?

#first construct data frame including data on buoyancy and first three PC axes of body shape variation 
c.buoy.shape<-cbind(c.buoy,c.shapePC1,c.shapePC2,c.shapePC3)
#remove species with missing data
na.omit(c.buoy.shape)->c.buoy.shape
colnames(c.buoy.shape) = c("Buoyancy", "ShapePC1", "ShapePC2", "ShapePC3")
as.data.frame(c.buoy.shape)->c.buoy.shape

#prune tree such that tips match species for which we have trait data
row.names(c.buoy.shape)->c.buoy.shape.tax
c.buoy.shape.tree<-keep.tip(noto.tree, c.buoy.shape.tax)

#then define covariance structure - we'll allow covariance structure to match that expected under Brownian Motion (BM) process
c.bm<-corBrownian(1, c.buoy.shape.tree, form=~c.buoy.shape.tax)

#define model 1 - relationship between buoyancy and body shape PC axis 1 under Brownian Motion model of covariance
model1<-gls(Buoyancy~ShapePC1, data=c.buoy.shape, correlation=c.bm, method = "ML")
summary(model1)

# define model 2 - relationship between buoyancy and body shape PC axis 1 under Ornstein-Uhlenbeck model of covariance
c.ou<-corMartins(1, c.buoy.shape.tree, fixed=FALSE, form=~c.buoy.shape.tax)
model2<-gls(Buoyancy~ShapePC1, data=c.buoy.shape, correlation=c.ou, method = "ML")
summary(model2)

##what is the relationship between depth and maximum body size?

#construct data frame for relevant trait data
c.depth.size<-cbind(c.depth,c.size)
na.omit(c.depth.size)->c.depth.size
colnames(c.depth.size) = c("MeanDepth", "MaxSize")
as.data.frame(c.depth.size)->c.depth.size

#prune tree such that tips match species for which we have trait data
row.names(c.depth.size)->c.depth.size.tax
c.depth.size.tree<-keep.tip(noto.tree, c.depth.size.tax)

#define model 3 - relationship between depth and body size under Brownian Motion model of covariance
c.bm<-corBrownian(1, c.depth.size.tree, form=~c.depth.size.tax)
model3<-gls(MaxSize~MeanDepth, data=c.depth.size, correlation=c.bm, method = "ML")
summary(model3)

# define model 4 - relationship between depth and body size under Ornstein-Uhlenbeck model of covariance
c.ou<-corMartins(1, c.depth.size.tree, fixed=FALSE, form=~c.depth.size.tax)
model4<-gls(MaxSize~MeanDepth, data=c.depth.size, correlation=c.ou, method = "ML")
summary(model4)

##what is the relationship between body size and the first PC axes of variation in diet?

#construct data frame for relevant trait data
c.diet.size<-cbind(c.dietPC1,c.dietPC2,c.dietPC3,c.size)
na.omit(c.diet.size)->c.diet.size
colnames(c.diet.size) = c("DietPC1", "DietPC2", "DietPC3", "MaxSize")
as.data.frame(c.diet.size)->c.diet.size

#prune tree such that tips match species for which we have trait data
row.names(c.diet.size)->c.diet.size.tax
c.diet.size.tree<-keep.tip(noto.tree, c.diet.size.tax)

#define model 5 - relationship between depth and body size under Brownian Motion model of covariance
c.bm<-corBrownian(1, c.diet.size.tree, form=~c.diet.size.tax)
model5<-gls(MaxSize~DietPC1, data=c.diet.size, correlation=c.bm, method = "ML")
summary(model5)

# define model 6 - relationship between depth and body size under Ornstein-Uhlenbeck model of covariance
c.ou<-corMartins(1, c.diet.size.tree, fixed=FALSE, form=~c.diet.size.tax)
model6<-gls(MaxSize~DietPC1, data=c.diet.size, correlation=c.ou, method = "ML")
summary(model6)


############################################
# Next complete analyses for notoperches
############################################

##isolate traits for notoperches
t.dat<-dat[15:31,]

#create unique variables for each of the individual traits
t.buoy<-as.matrix(t.dat$Buoyancy)
row.names(t.buoy)<-t.dat$Tree_label

t.depth<-as.matrix(t.dat$MeanDepth)
row.names(t.depth)<-t.dat$Tree_label

t.size<-as.matrix(t.dat$MaxSize)
row.names(t.size)<-t.dat$Tree_label

t.shapePC1<-as.matrix(t.dat$BodyShapePC1)
row.names(t.shapePC1)<-t.dat$Tree_label

t.shapePC2<-as.matrix(t.dat$BodyShapePC2)
row.names(t.shapePC2)<-t.dat$Tree_label

t.shapePC3<-as.matrix(t.dat$BodyShapePC3)
row.names(t.shapePC3)<-t.dat$Tree_label

t.dietPC1<-as.matrix(t.dat$DietPC1)
row.names(t.dietPC1)<-t.dat$Tree_label

t.dietPC2<-as.matrix(t.dat$DietPC2)
row.names(t.dietPC2)<-t.dat$Tree_label

t.dietPC3<-as.matrix(t.dat$DietPC3)
row.names(t.dietPC3)<-t.dat$Tree_label

###starting the phylogenetic GLS

##what is the relationship between buoyancy and each of the first three PC axes of body shape variation?

#first construct data frame including data on buoyancy and first three PC axes of body shape variation 
t.buoy.shape<-cbind(t.buoy,t.shapePC1,t.shapePC2,t.shapePC3)
#remove species with missing data
na.omit(t.buoy.shape)->t.buoy.shape
colnames(t.buoy.shape) = c("Buoyancy", "ShapePC1", "ShapePC2", "ShapePC3")
as.data.frame(t.buoy.shape)->t.buoy.shape

#prune tree such that tips match species for which we have trait data
row.names(t.buoy.shape)->t.buoy.shape.tax
t.buoy.shape.tree<-keep.tip(noto.tree, t.buoy.shape.tax)

#then define covariance structure - we'll allow covariance structure to match that expected under Brownian Motion (BM) process
t.bm<-corBrownian(1, t.buoy.shape.tree, form=~t.buoy.shape.tax)

#define model 7 - relationship between buoyancy and body shape PC axis 1 under Brownian Motion model of covariance
model7<-gls(Buoyancy~ShapePC1, data=t.buoy.shape, correlation=t.bm, method = "ML")
summary(model7)

# define model 8 - relationship between buoyancy and body shape PC axis 1 under Ornstein-Uhlenbeck model of covariance
t.ou<-corMartins(1, t.buoy.shape.tree, fixed=FALSE, form=~t.buoy.shape.tax)
model8<-gls(Buoyancy~ShapePC1, data=t.buoy.shape, correlation=t.ou, method = "ML")
summary(model8)

##what is the relationship between depth and maximum body size?

#construct data frame for relevant trait data
t.depth.size<-cbind(t.depth,t.size)
na.omit(t.depth.size)->t.depth.size
colnames(t.depth.size) = c("MeanDepth", "MaxSize")
as.data.frame(t.depth.size)->t.depth.size

#prune tree such that tips match species for which we have trait data
row.names(t.depth.size)->t.depth.size.tax
t.depth.size.tree<-keep.tip(noto.tree, t.depth.size.tax)

#define model 9 - relationship between depth and body size under Brownian Motion model of covariance
t.bm<-corBrownian(1, t.depth.size.tree, form=~t.depth.size.tax)
model9<-gls(MaxSize~MeanDepth, data=t.depth.size, correlation=t.bm, method = "ML")
summary(model9)

# define model 10 - relationship between depth and body size under Ornstein-Uhlenbeck model of covariance
t.ou<-corMartins(1, t.depth.size.tree, fixed=FALSE, form=~t.depth.size.tax)
model10<-gls(MaxSize~MeanDepth, data=t.depth.size, correlation=t.ou, method = "ML")
summary(model10)

##what is the relationship between body size and the first PC axes of diet variation?

#construct data frame for relevant trait data
t.diet.size<-cbind(t.dietPC1,t.dietPC2,t.dietPC3,t.size)
na.omit(t.diet.size)->t.diet.size
colnames(t.diet.size) = c("DietPC1", "DietPC2", "DietPC3", "MaxSize")
as.data.frame(t.diet.size)->t.diet.size

#prune tree such that tips match species for which we have trait data
row.names(t.diet.size)->t.diet.size.tax
t.diet.size.tree<-keep.tip(noto.tree, t.diet.size.tax)

#define model 11 - relationship between depth and body size under Brownian Motion model of covariance
t.bm<-corBrownian(1, t.diet.size.tree, form=~t.diet.size.tax)
model11<-gls(MaxSize~DietPC1, data=t.diet.size, correlation=t.bm, method = "ML")
summary(model11)

# define model 12 - relationship between depth and body size under Ornstein-Uhlenbeck model of covariance
t.ou<-corMartins(1, t.diet.size.tree, fixed=FALSE, form=~t.diet.size.tax)
model12<-gls(MaxSize~DietPC1, data=t.diet.size, correlation=t.ou, method = "ML")
summary(model12)
