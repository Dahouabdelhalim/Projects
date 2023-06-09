## clean environment & plots
rm(list=ls())

#load packages 
library(tidyverse) 
library(ape)
library(geiger)
library(phytools)
library(Rmisc)
library(picante)
library(nlme)
library(MuMIn)
library(scales)
library(cowplot)

#set working directory

#load data
beast <- readxl::read_excel("wd/Kimmitt_etal_beast_summary.xlsx", na = "NA")

other.dat <- read_csv("wd/migration_distance_Winger&Pegan.csv")

other.dat <- other.dat %>%
  select(Species, Migratory_status, boreal_migdist)


beast$Species <- gsub(" ", "_", beast$species)
beast$Species
beast$Species[beast$Species == "Dryobates_villosus"] <- "Picoides_villosus"

dat <- beast %>%
  left_join(other.dat, by = "Species")

## set theme
th=theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(legend.position = "none")+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(strip.text=element_text(size=11),
        axis.text=element_text(size=11),
        axis.title=element_text(size=12)) 


#plot the variation in migratory distance as a histogram
ggplot(dat,aes(x=boreal_migdist))+geom_histogram(,color="black",fill="gray",bins=8)+th + xlab("Migratory distance (km)") +
  ylab("Number of species") 



### PGLS 
#load tree 
tree <- read.tree("/Users/abbykimmitt/Dropbox (University of Michigan)/winger_lab/demographic_history/beast/jetz_matched_names.tre")

# making sure species names match between dataset and tree 
dat$Species[!dat$Species %in% tree$tip.label]
dat$Species[dat$Species == "Regulus_calendula"] <- "Corthylio_calendula"
dat$Species[!dat$Species %in% tree$tip.label]

# prune tree down to species we need
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% dat$Species])
plot(tree, cex = 0.5)

# assigning row names so that branch lengths are recognized 
dat<- as.data.frame(dat)
row.names(dat) <- dat$Species


### Historic Ne (Ne1)

#removing missing data
dat.Ne1_ac2 <- dat[complete.cases(dat$Ne1_ac2, dat$boreal_migdist, dat$mass),]
tree.Ne1_ac2 <-  drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat.Ne1_ac2$Species])

dat.Ne1_ac4 <- dat[complete.cases(dat$Ne1_ac4, dat$boreal_migdist),]
tree.Ne1_ac4 <-  drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat.Ne1_ac4$Species])

dat.Ne1_3c2 <- dat[complete.cases(dat$Ne1_3c2, dat$boreal_migdist),]
tree.Ne1_3c2 <-  drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat.Ne1_3c2$Species])

dat.Ne1_3c4 <- dat[complete.cases(dat$Ne1_3c4, dat$boreal_migdist),]
tree.Ne1_3c4 <-  drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat.Ne1_3c4$Species])

dat.Ne1_cytb <- dat[complete.cases(dat$Ne1_cytb, dat$boreal_migdist),]
tree.Ne1_cytb <-  drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat.Ne1_cytb$Species])


#Migration distance 
#models without considering phylogeny 
Ne1_ac2 <-  lm(Ne1_ac2~ boreal_migdist, data = dat.Ne1_ac2)
Ne1_ac4 <-  lm(Ne1_ac4~ boreal_migdist, data = dat.Ne1_ac4)
Ne1_3c2 <-  lm(Ne1_3c2~ boreal_migdist, data = dat.Ne1_3c2)
Ne1_3c4 <-  lm(Ne1_3c4~ boreal_migdist, data = dat.Ne1_3c4)
Ne1_cytb <-  lm(Ne1_cytb~ boreal_migdist, data = dat.Ne1_cytb)

#checking phylogenetic signal
res_ac2 = residuals(Ne1_ac2)
psigNe1_ac2 <- phylosig(tree.Ne1_ac2, res_ac2, method="lambda", test=T)
psigNe1_ac2$P # not significant

res_ac4 = residuals(Ne1_ac4)
psigNe1_ac4 <- phylosig(tree.Ne1_ac4, res_ac4, method="lambda", test=T)
psigNe1_ac4$P #not significant

res_3c2 = residuals(Ne1_3c2)
psigNe1_3c2 <- phylosig(tree.Ne1_3c2, res_3c2, method="lambda", test=T)
psigNe1_3c2$P #not significant

res_3c4 = residuals(Ne1_3c4)
psigNe1_3c4 <- phylosig(tree.Ne1_3c4, res_3c4, method="lambda", test=T)
psigNe1_3c4$P # not significant

res_cytb = residuals(Ne1_cytb)
psigNe1_cytb <- phylosig(tree.Ne1_cytb, res_cytb, method="lambda", test=T)
psigNe1_cytb$P # not significant

#summarizing OLS without phylogenetic signal 
summary(Ne1_ac2) #not sig
summary(Ne1_ac4) #not sig
summary(Ne1_3c2) #not sig
summary(Ne1_3c4) #not sig 
summary(Ne1_cytb) #not sig 

##################################

### Pop expansion initiation (time1)

#removing missing data
dat.time1_ac2 <- dat[complete.cases(dat$time1_ac2, dat$boreal_migdist),]
tree.time1_ac2 <-  drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat.time1_ac2$Species])

dat.time1_ac4 <- dat[complete.cases(dat$time1_ac4, dat$boreal_migdist),]
tree.time1_ac4 <-  drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat.time1_ac4$Species])

dat.time1_3c2 <- dat[complete.cases(dat$time1_3c2, dat$boreal_migdist),]
tree.time1_3c2 <-  drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat.time1_3c2$Species])

dat.time1_3c4 <- dat[complete.cases(dat$time1_3c4, dat$boreal_migdist),]
tree.time1_3c4 <-  drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat.time1_3c4$Species])

dat.time1_cytb <- dat[complete.cases(dat$time1_cytb, dat$boreal_migdist),]
tree.time1_cytb <-  drop.tip(tree, tree$tip.label[!tree$tip.label%in%dat.time1_cytb$Species])

#Migration distance 
#models without considering phylogeny 

time1_ac2 <-  lm(time1_ac2~ boreal_migdist, data = dat.time1_ac2)
time1_ac4 <-  lm(time1_ac4~ boreal_migdist, data = dat.time1_ac4)
time1_3c2 <-  lm(time1_3c2~ boreal_migdist, data = dat.time1_3c2)
time1_3c4 <-  lm(time1_3c4~ boreal_migdist, data = dat.time1_3c4)
time1_cytb <-  lm(time1_cytb~ boreal_migdist, data = dat.time1_cytb)

#checking phylogenetic signal
res_ac2 = residuals(time1_ac2)
psigtime1_ac2 <- phylosig(tree.time1_ac2, res_ac2, method="lambda", test=T)
psigtime1_ac2$P #significant

res_ac4 = residuals(time1_ac4)
psigtime1_ac4 <- phylosig(tree.time1_ac4, res_ac4, method="lambda", test=T)
psigtime1_ac4$P #significant

res_3c2 = residuals(time1_3c2)
psigtime1_3c2 <- phylosig(tree.time1_3c2, res_3c2, method="lambda", test=T)
psigtime1_3c2$P #significant

res_3c4 = residuals(time1_3c4)
psigtime1_3c4 <- phylosig(tree.time1_3c4, res_3c4, method="lambda", test=T)
psigtime1_3c4$P #significant

res_cytb = residuals(time1_cytb)
psigtime1_cytb <- phylosig(tree.time1_cytb, res_cytb, method="lambda", test=T)
psigtime1_cytb$P # not significant


# full model pgls  
bm<-corBrownian(1, tree, form = ~Species)

gls_time1_ac2<- gls(time1_ac2 ~ boreal_migdist, data = dat.time1_ac2, correlation=bm, method="ML")
summary(gls_time1_ac2) #not sig

gls_time1_ac4<- gls(time1_ac4 ~ boreal_migdist, data = dat.time1_ac4, correlation=bm, method="ML")
summary(gls_time1_ac4) #not sig

gls_time1_3c2<- gls(time1_3c2 ~ boreal_migdist, data = dat.time1_3c2, correlation=bm, method="ML")
summary(gls_time1_3c2) #not sig

gls_time1_3c4<- gls(time1_3c4 ~ boreal_migdist, data = dat.time1_3c4, correlation=bm, method="ML")
summary(gls_time1_3c4) #not sig

summary(time1_cytb) #not sig

######################################

## figures for time of population expansion 


dat$title1<-"All codons, slow calibration"
dat$title2<-"All codons, fast calibration"
dat$title3<-"Third codon, slow calibration"
dat$title4<-"Third codon, fast calibration"
dat$title5<-"cytb, 2% rule"


p1<-ggplot(data=dat, aes(x=boreal_migdist,y=time1_ac2))+ geom_point()+scale_y_continuous(labels = label_number(scale = 1e-3), limits=c(0,250000)) + 
  ylab("Initiation of population size expansion (1000 ybp)") + xlab("Migration distance (km)") + facet_grid(. ~ title1)+ theme(strip.text.x = element_text(size = 12))
    
p2<-ggplot(data=dat, aes(x=boreal_migdist,y=time1_ac4))+ geom_point()+scale_y_continuous(labels = label_number(scale = 1e-3), limits=c(0,250000))+ 
  ylab("Initiation of population size expansion (1000 ybp)") + xlab("Migration distance (km)")+ facet_grid(. ~ title2)+ theme(strip.text.x = element_text(size = 12))

p3<-ggplot(data=dat, aes(x=boreal_migdist,y=time1_3c2))+ geom_point()+scale_y_continuous(labels = label_number(scale = 1e-3), limits=c(0,250000))+ 
  ylab("Initiation of population size expansion (1000 ybp)") + xlab("Migration distance (km)")+ facet_grid(. ~ title3)+ theme(strip.text.x = element_text(size = 12))

p4<-ggplot(data=dat, aes(x=boreal_migdist,y=time1_3c4))+ geom_point()+scale_y_continuous(labels = label_number( scale = 1e-3), limits=c(0,250000))+ 
  ylab("Initiation of population size expansion (1000 ybp)") + xlab("Migration distance (km)")+ facet_grid(. ~ title4)+ theme(strip.text.x = element_text(size = 12))

p5<-ggplot(data=dat, aes(x=boreal_migdist,y=time1_cytb))+ geom_point()+scale_y_continuous(labels = label_number(scale = 1e-3), limits=c(0,250000))+
ylab("Initiation of population size expansion (1000 ybp)") +xlab("Migration distance (km)")+ facet_grid(. ~ title5)+ theme(strip.text.x = element_text(size = 12))




plot1<-plot_grid(p1 + theme(axis.text.x = element_blank(),
                                      axis.title.x=element_blank()),
                   p2 + theme(axis.text.x = element_blank(),
                              axis.title.x=element_blank(),
                              axis.text.y = element_blank(),
                              axis.title.y = element_blank()),
                 p3, p4 + theme(axis.text.y = element_blank(),
                                axis.title.y = element_blank()),
                   nrow = 2,
                   labels = c('A', 'B', 'C', 'D'), label_size = 12, align="v")

plot2<-plot_grid(p5, labels=c('E'))

plot_grid(plot1, NULL, plot2, rel_widths = c(1,-0.2,1), nrow=1, scale =c(1,1,0.5))




#figures for effective population size 


p1<-ggplot(data=dat, aes(x=boreal_migdist,y=Ne1_ac2))+ geom_point()+scale_y_continuous(labels = label_number(scale = 1e-4), limits=c(0,4000000)) + 
  ylab(expression(median~italic(N[e])~~(~x10000))) + xlab("Migration distance (km)")+ facet_grid(. ~ title1)+ theme(strip.text.x = element_text(size = 12))

p2<-ggplot(data=dat, aes(x=boreal_migdist,y=Ne1_ac4))+ geom_point()+scale_y_continuous(labels = label_number(scale = 1e-4), limits=c(0,4000000)) + 
  ylab(expression(median~italic(N[e])~~(~x10000))) + xlab("Migration distance (km)")+ facet_grid(. ~ title2)+ theme(strip.text.x = element_text(size = 12))

p3<-ggplot(data=dat, aes(x=boreal_migdist,y=Ne1_3c2))+ geom_point()+scale_y_continuous(labels = label_number(scale = 1e-4), limits=c(0,4000000)) + 
  ylab(expression(median~italic(N[e])~~(~x10000))) + xlab("Migration distance (km)")+ facet_grid(. ~ title3)+ theme(strip.text.x = element_text(size = 12))

p4<-ggplot(data=dat, aes(x=boreal_migdist,y=Ne1_3c4))+ geom_point()+scale_y_continuous(labels = label_number(scale = 1e-4), limits=c(0,4000000)) + 
  ylab(expression(median~italic(N[e])~~(~x10000))) + xlab("Migration distance (km)")+ facet_grid(. ~ title4)+ theme(strip.text.x = element_text(size = 12))

p5<-ggplot(data=dat, aes(x=boreal_migdist,y=Ne1_cytb))+geom_point()+scale_y_continuous(labels = label_number(scale = 1e-4), limits=c(0,4000000)) + 
  ylab(expression(median~italic(N[e])~~(~x10000))) + xlab("Migration distance (km)")+ facet_grid(. ~ title5)+ theme(strip.text.x = element_text(size = 12))




plot1<-plot_grid(p1 + theme(axis.text.x = element_blank(),
                            axis.title.x=element_blank()),
                 p2 + theme(axis.text.x = element_blank(),
                            axis.title.x=element_blank(),
                            axis.text.y = element_blank(),
                            axis.title.y = element_blank()),
                 p3, p4 + theme(axis.text.y = element_blank(),
                                axis.title.y = element_blank()),
                 nrow = 2,
                 labels = c('A', 'B', 'C', 'D'), label_size = 12, align="v")

plot2<-plot_grid(p5, labels=c('E'))

plot_grid(plot1, NULL, plot2, rel_widths = c(1,-0.2,1), nrow=1, scale =c(1,1,0.5))
