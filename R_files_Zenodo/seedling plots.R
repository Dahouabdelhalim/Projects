##########################################################
#   Gawel et al. RSOS Contrasting Ecological Roles      ##
#   of non-native ungulates in a novel ecosystem.       ##
#   Scripts for analysing and plotting seedling         ##
#   plot data.                                          ##
##########################################################

########Analysis using generalized linear mixed effects models#######

##load packages##
library(lme4)
library(AICcmodavg)

##import dataset 

ungulate <- read.csv("seedlingplots_guam.csv", stringsAsFactors=F)

#assign class of each column
factor_cols <- c("island","site", "species", "trt", "deer.abundance", "pig.abundance")
numeric_cols <- c("numplant", "alive", "dead", "length.exposure")
ungulate[factor_cols] <- lapply(ungulate[factor_cols], as.factor)
ungulate[numeric_cols] <- lapply(ungulate[numeric_cols], as.numeric)

#standardize time
ungulate$length.exposure <-scale(ungulate$length.exposure)

#look at data
str(ungulate)
with(ungulate, table(species, length.exposure))

#create variable for survival
ungulate$survival<-cbind(ungulate$alive,ungulate$dead)

###Model selection between trt*time and trt*species###
timetrt <- glmer(survival~length.exposure*trt + (1|site), family = binomial, data=ungulate)
spptrt <- glmer(survival~species*trt + (1|site), family = binomial, data=ungulate)
aictab(list(timetrt, spptrt), modnames=c("timetrt", "spptrt"))

###Model selection for individual species###

#######AGLAIA#########
aglaia<-subset(ungulate,species=="aglaia")
survivalag<-cbind(aglaia$alive,aglaia$dead)

#candidate models
trtag<-glmer(survivalag~trt+(1|site),family=binomial,data=aglaia)
nullag<-glmer(survivalag~1+(1|site),family=binomial,data=aglaia)

#find best model
aictab(list(trtag, nullag), modnames=c("trtag", "nullag")) 

#######NEISO########
##***Neisosperma is now called Ochrosia, so will appear as "Ochrosia oppositifolia" in figures and text
neiso<-subset(ungulate,species=="neiso")
survivalne<-cbind(neiso$alive,neiso$dead)

#candidate models
trtne<-glmer(survivalne~trt+(1|site),family=binomial,data=neiso)
nullne<-glmer(survivalne~1+(1|site),family=binomial,data=neiso)

#find best model
aictab(list(trtne, nullne), modnames=c("trtne", "nullne"))

#####CARICA########
carica<-subset(ungulate,species=="carica")
survivalca<-cbind(carica$alive,carica$dead)

#candidate models
trtca<-glmer(survivalca~trt+(1|site),family=binomial,data=carica)
nullca<-glmer(survivalca~1+(1|site),family=binomial,data=carica)

#find best model
aictab(list(trtca, nullca), modnames=c("trtca", "nullca"))

########MORINDA###############

morinda<-subset(ungulate,species=="morinda")
survivalmo<-cbind(morinda$alive,morinda$dead)

#candidate models
trtmo<-glmer(survivalmo~trt+(1|site),family=binomial,data=morinda)
nullmo<-glmer(survivalmo~1+(1|site),family=binomial,data=morinda)

#find best model
aictab(list(trtmo,nullmo), modnames=c("trtmo", "nullmo"))

######PREMNA##########
premna<-subset(ungulate,species=="premna")
survivalpr<-cbind(premna$alive,premna$dead)

#candidate models
trtpr<-glmer(survivalpr~trt+(1|site),family=binomial,data=premna)
nullpr<-glmer(survivalpr~1+(1|site),family=binomial,data=premna)

#find best model
aictab(list(trtpr,nullpr), modnames=c("trtpr", "nullpr"))

########PSYCHOTRIA############

psychotria<-subset(ungulate,species=="psychotria")
survivalps<-cbind(psychotria$alive,psychotria$dead)

#candidate models
trtps<-glmer(survivalps~trt+(1|site),family=binomial,data=psychotria)
nullps<-glmer(survivalps~1+(1|site),family=binomial,data=psychotria)

#find best model
aictab(list(trtps,nullps), modnames=c("trtps", "nullps"))


########################################
#               GRAPHS                 #
########################################


###load packages for graphs/figures###
library(ggplot2)
library(binom)
library(plyr)
library(ggsignif)

#############
#create dataframe with binomial confidence intervals
newdf<-ddply(ungulate,.(trt, species),summarise,
             prop=sum(alive)/sum(numplant),
             low=prop.test(sum(alive),sum(numplant))$conf.int[1],
             upper=prop.test(sum(alive),sum(numplant))$conf.int[2])

##order x-axis by biggest difference between treatments to least##
##Order from most to least significant
newdf$species <- factor(newdf$species, levels = c("carica","psychotria","morinda","premna","aglaia","neiso"))

###Make graph###
p <- ggplot(newdf,aes(species,prop, ymin=low, ymax=upper, shape=trt))+
  geom_point(stat="identity", position=position_dodge(width=0.3), size=3)+
  geom_errorbar(width=0.2, position=position_dodge(width=0.3))+
  ylab("Proportion seedling survival")+
  ylim(0,1.1)+
  scale_x_discrete("Species", labels=c("Carica", "Psychotria","Morinda","Premna","Aglaia","Ochrosia"))+
  scale_shape_manual(values=c(17,19), breaks=c("fenced","ungulate"), labels=c("No ungulates", "Ungulates"))+
  theme_minimal()+
  theme(axis.title.y=element_text(size=10, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        strip.text=element_text(size=10, face="italic"),
        strip.background = element_rect(colour = "white", fill="white"))

###annotate to add asterisks and NS above bars###
p +annotate("text",x=1, y=0.9,label="*") +annotate("text",x=2,y=0.93,label="*") +annotate("text",x=3,y=0.94,label="*") +annotate("text",x=4,y=0.5,label="*")+annotate("text",x=5,y=0.93,label="NS")+annotate("text",x=6,y=0.97,label="NS")




