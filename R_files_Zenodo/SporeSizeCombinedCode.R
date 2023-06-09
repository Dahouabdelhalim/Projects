
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(cowplot)


#########################################################
#Analyzing field patterns

#Load data
MI <- read.table("MichiganField2015.txt", header=TRUE, sep="\\t", dec=".", strip.white=TRUE, stringsAsFactors = FALSE)
IN <- read.table("IndianaField2015.txt", header=TRUE, sep="\\t", dec=".", strip.white=TRUE, stringsAsFactors = FALSE)

#Filter to just lakes with Metschnikowia (in either host)
MLakes<-MI%>%group_by(Lake)%>%summarize(totalMetsch=sum(prevalence, na.rm=TRUE))
MLakes.1<-filter(MLakes, totalMetsch>0)
ML<-MLakes.1$Lake
MLakes.M<-filter(MI, Lake%in%ML)

ILakes<-IN%>%group_by(Lake)%>%summarize(totalMetsch=sum(prevalence, na.rm=TRUE))
ILakes.1<-filter(ILakes, totalMetsch>0)
IL<-ILakes.1$Lake
ILakes.M<-filter(IN, Lake%in%IL)

#Combine Michigan and Indiana data for lakes with some Metschnikowia
colnames(MLakes.M)<-c("Lake","Julian","Species","Total","prevalence")
MLakes.M$State<-rep("MI",length(MLakes.M$Lake))
ILakes.M$State<-rep("IN",length(ILakes.M$Lake))
In.Mi<-bind_rows(MLakes.M, ILakes.M)

In.Mi<-In.Mi[order(In.Mi$Julian),] #make sure julian days are in order


#Caluclate integrated area of infection
OUT<-NULL
STATE<-unique(In.Mi$State)
for (p in STATE){
  thisstate<-filter(In.Mi, State==p) #filter by state
  LAKE<-unique(thisstate$Lake)
  for(i in LAKE){
    thislake<-filter(thisstate, Lake==i) #filter by lake
    HOST<-unique(thislake$Species)
    for(j in HOST){
      thishost<-filter(thislake, Species==j) #filter by host species
      thishost1<-filter(thishost, !is.na(prevalence)) #take out NAs for prevalence (0s are still left in)
      TOTS<-NULL #make an empty matrix to include sequential prevalences. This will be used to calculate integrated areas.
      for (m in 1:length(thishost1$prevalence)){
        total<-thishost1[m, 5] #prevalence
        jul<-thishost1[m, 2] #julian day
        tots<-c(m, jul, total) #new row in dataframe with m, julian day, and total
        TOTS<-rbind(TOTS, tots) #bind all rows together (different one for each lake/year/parasite)
      }
      colnames(TOTS)<-c("m","date","total")
      rownames(TOTS)<-NULL
      TOTS<-as.data.frame(TOTS)
      CURVE<-NULL #make an empty matrix for calulating each chunk of the integrated area by trapezoid rule.
      for(n in 1:(length(TOTS$total)-1)){
        curve<-0.5*((TOTS[n,3]+TOTS[n+1,3])*(TOTS[n+1,2]-TOTS[n,2])) #trapezoid rule; this calculates area of each chunk
        CURVE<-c(CURVE, curve)} #write down each chunk
      infarea<-sum(CURVE) #take the sum
      output<-c(p, i, j, infarea)
      OUT<-rbind(OUT, output)
    }
  }
}

colnames(OUT)<-c("State","Lake","Species", "infected.area")
areas<-as.data.frame(OUT)
row.names(areas)<-NULL
areas$infected.area<-as.numeric(as.character(areas$infected.area))

#Make into a wide dataset
wide.prev<-spread(areas, key=Species, value=infected.area)
colnames(wide.prev)<-c("Region","Lake", "Cerio.infarea", "dent.infarea")

#Figure 1 - relationship between Daphnia integrated metsch prevalence and Ceriodaphnia integrated metsch prevalence
Fig1<-ggplot(wide.prev, aes(x=dent.infarea, y=Cerio.infarea, shape=Region))+
  geom_point(size=3, alpha=0.5)+
  labs(x=expression(atop("Outbreak size in"~italic("D. dentifera"), paste("(prevalence x days)"))),
       y=expression(atop("Outbreak size in "~italic("C. dubia"), paste("(prevalence x days)"))))+
  theme(axis.title = element_text(size=12))+
  theme(axis.text = element_text(size=10, color="black"))+
  theme(legend.position = c(0.70, 0.38))+
  theme(legend.background = element_rect(linetype=1, size = 0.5, colour = 1))+
  theme(legend.title = element_text(size=12))+
  theme(legend.text = element_text(size=10))

save_plot("Figure1.jpeg", Fig1, base_height = 3.14961, base_width = 3.14961)

#statistical models
#is the integrated area of daphnia prevalence related to the integrated area of Cerio prevalence?
M1<-lm(dent.infarea~Cerio.infarea, data=wide.prev) 
summary(M1)



####################################################
#Analyzing microsatellite data
library(adegenet)
library(ape)
library(magrittr)
library(treemap)
library(pegas)
library(mmod)
library(poppr)
library(dendextend)
library(ggtree)

#Load data
microsats.data <- read.table("microsats.data.txt",header=TRUE, sep="\\t",dec=".",strip.white=TRUE, na.strings="",fill=TRUE, stringsAsFactors = FALSE)

#how many differences between genotypes
id1<-c("WoodlandDaphnia10/15.1(MI)",
       "MillDaphnia9/8.1(MI)",
       "GoslingCerio10/26.1(MI)",
       "BenefielCerio10/15.2(IN)",
       "SycamoreCerio10/15.1(IN)",
       "BenefielCerio10/15.3(IN)")

fordistance <- filter(microsats.data, id%in% id1) #just look at these individuals, which represent different genotypes
fordistance <- select(fordistance,id,State_Host_Lake, p3,p7,p8,p9,primer10,p11,p12,p17,p19)

write.table(fordistance, "FD2.txt", col.names=TRUE, row.names=FALSE, sep="\\t", append=FALSE)

#Getting the data into the right format (genalex)
d <- read.loci("FD2.txt", col.loci=3:11, row.names=1,col.pop=2, header=TRUE, allele.sep="-")

D<-loci2genind(d, ploidy = 1)

genind2genalex(D, filename = "FD2.csv", quiet = FALSE, pop = NULL,
               allstrata = TRUE, geo = FALSE, geodf = "xy", sep = ",",
               sequence = FALSE, overwrite=TRUE)

D2<-read.genalex("FD2.csv")
A.distances<-provesti.dist(D2)
A.distances<-as.matrix(A.distances)
avgdistances<-mean(A.distances[ lower.tri( A.distances ) ])
9*avgdistances[1] #9 loci * average distance between genotypes.

locus_table(D2)


#All data - copepod
m.data.daphnia.cerio <- filter(microsats.data, Host!="copepod")
forpopgen <- select(m.data.daphnia.cerio,id,State_Host_Lake, p3,p7,p8,p9,primer10,p11,p12,p17,p19)
write.table(forpopgen, "microsats.data.popgen.txt", col.names=TRUE, row.names=FALSE, sep="\\t", append=FALSE)

z <- read.loci("microsats.data.popgen.txt", col.loci=3:11, row.names=1,col.pop=2, header=TRUE, allele.sep="-")
Z <- loci2genind(z, ploidy = 1) #Specify ploidy
#change format for genalex
genind2genalex(Z, filename = "microsats.data.popgen.csv", quiet = FALSE, pop = NULL,
               allstrata = TRUE, geo = FALSE, geodf = "xy", sep = ",",
               sequence = FALSE, overwrite=TRUE)
#Load data again
microsats.genalex<-read.genalex("microsats.data.popgen.csv")
locus_table(microsats.genalex)
info_table(microsats.genalex, type='missing',plot=TRUE)
mlg(microsats.genalex)

splitStrata(microsats.genalex) <- ~State/Host/Lake
setPop(microsats.genalex) <- ~Host
microsats.genalex%>% clonecorrect() %>% ia(sample = 999)

poppr(microsats.genalex)

#Get confidence intervals for heterozygosity; Pop 1 is Cerio, Pop 2 is Daphnia
OUT<-NULL
for (j in c(1,2)){
  Mpop <- microsats.genalex[pop = j]
  het <- function(i){x <- poppr(i);c(Hexp = mean(x$Hexp))}
  Msamples <- replicate(1000, shufflepop(Mpop, method = 2))
  res <- sapply(Msamples, het)
  res1<-as.matrix(res)
  CI<-apply(res1, 2, quantile, c(0.05, 0.95))
  output<-c(j, CI[1], CI[2])
  OUT<-rbind(OUT, output)
}
CI<-OUT

#To figure out the confidence intervals around total Hexp
het <- function(i){x <- poppr(i);tot<-filter(x, Pop=="Total"); c(Hexp = mean(tot$Hexp))}
Msamples <- replicate(1000, shufflepop(microsats.genalex, method = 2))
res <- sapply(Msamples, het)
res1<-as.matrix(res)
CI<-apply(res1, 2, quantile, c(0.05, 0.95))

setPop(microsats.genalex) <- ~Host/State/Lake

#AMOVA
AMOVA1 <- poppr.amova(microsats.genalex, ~Host/State/Lake,correction="lingoes",missing="ignore",within=FALSE, squared=FALSE)
set.seed(999)
AMOVA1.signif <- randtest(AMOVA1, nrepet=999)

AMOVA2 <- poppr.amova(microsats.genalex, ~State/Lake/Host,correction="lingoes",missing="ignore",within=FALSE)
set.seed(999)
AMOVA2.signif <- randtest(AMOVA2, nrepet=999)

#Build Tree (include the copepods here)
fortree <- select(microsats.data,id,State_Host_Lake, p3,p7,p8,p9,primer10,p11,p12,p17,p19)
write.table(fortree, "microsats.data.tree.txt", col.names=TRUE, row.names=FALSE, sep="\\t", append=FALSE)

z <- read.loci("microsats.data.tree.txt", col.loci=3:11, row.names=1,col.pop=2, header=TRUE, allele.sep="-")
Z <- loci2genind(z, ploidy = 1) #Specify ploidy
#change format for genalex
genind2genalex(Z, filename = "microsats.data.tree.csv", quiet = FALSE, pop = NULL,
               allstrata = TRUE, geo = FALSE, geodf = "xy", sep = ",",
               sequence = FALSE, overwrite=TRUE)
#Load data again
microsats.tree<-read.genalex("microsats.data.tree.csv")
microsats.genind <- genclone2genind(microsats.tree)

set.seed(999)
tree<-aboot(microsats.genind, dist = provesti.dist, sample = 1000, tree = "upgma", cutoff = 40, quiet = TRUE)

hostpalette <- c("firebrick","black", "dodgerblue2")
HOST <- factor(microsats.data$Host)
STATE <- factor(microsats.data$State)
statepalette <- c(19,17)
info<-as.data.frame(cbind(microsats.data$Host, microsats.data$State))

TREE<-ggtree(tree, size=1)+geom_tiplab(color=hostpalette[HOST], size=3, offset=0.02)+
  geom_treescale(width=0.05, x=0, y=43, fontsize=5)+ggplot2::xlim(-0.1, 1.2)+
  ggplot2::ylim(-5, 55)+
  geom_tippoint(color=hostpalette[HOST], shape=statepalette[STATE], size=3, position=position_nudge(x=.01))+
  geom_strip(4, 51, barsize=1, color='black', label=expression(italic("D. dentifera")*"-associated"), offset=0.5, fontsize = 4) +
  geom_strip(30, 3, barsize=1, color="black", label=expression(italic("C. dubia")*"-associated"), offset=0.5, fontsize = 4)+
  geom_strip(35, 3, barsize=1, color="black", label=expression("IN"), offset=0.4, fontsize = 4)+
  geom_strip(21, 30, barsize=1, color="black", label=expression("MI"), offset=0.4, fontsize = 4)+
  geom_nodelab(size = 3, col= "black", nudge_y = 0.8, nudge_x = -0.03)

#This is just so I can make a legend
stuff<-c(1,1,2,2,3,3,4,4,5)
colorst<-c("a","a","a","a","b","b","b","b","c")
shapes<-c("f","f","f","r","r","r","r","r","r")
Leg<-as.data.frame(cbind(stuff, colorst, shapes))
plotforLegend<-ggplot(Leg, aes(x=colorst,y=stuff, shape=shapes, color=colorst ))+
  geom_point(size=3)+
  scale_color_manual(values= c("firebrick","dodgerblue2","black"), name="Host species", 
                     labels=c(expression(italic("C. dubia")), expression(italic("D. dentifera")), "copepod"))+
  scale_shape_manual(values=c(19,17), name="State", labels=c("Indiana","Michigan"))+
  theme_bw()+theme(legend.title = element_text(size=16))+
  theme(legend.text = element_text(size=14))+theme(legend.text.align = 0)+
  theme(legend.key.size =  )
Leg<-get_legend(plotforLegend)

TREE3<- TREE + annotation_custom(grob = Leg, xmin=0, xmax=0.1, ymin=48, ymax=50)

save_plot("Figure2.jpg", TREE3, base_width = 6.53543, base_height = 8)


#########################################################
library(emmeans)

#Benefiel cross infection experiment analysis
Benefiel<-read.csv("BenefielInfections.csv", header=TRUE, sep=",")
Benefiel.sporecounts<-read.csv("BenefielSporeCounts.csv", header=TRUE, sep=",")
Goose<-read.csv("GooseInfections.csv", header=TRUE, sep=",")
Goose.sporecounts<-read.csv("GooseSporeCounts.csv", header=TRUE, sep=",")

spore.size <- read.table("Spore_Measurements.txt", sep="\\t", 
                         header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
colnames(spore.size)<-c("hostname","length","countdate")
spore.strain<-read.csv("FragAnalysisResults.Experiment.csv", sep=",", 
                       header=TRUE, stringsAsFactors = FALSE, fill=TRUE)


#Assign different clone numbers for shapes in figure
Benefiel$CloneNumber[Benefiel$Clone=="Ben Cerio 15"]<-1
Benefiel$CloneNumber[Benefiel$Clone=="Ben Cerio 6"]<-2
Benefiel$CloneNumber[Benefiel$Clone=="Ben Cerio 10"]<-3
Benefiel$CloneNumber[Benefiel$Clone=="Ben Cerio 1"]<-4
Benefiel$CloneNumber[Benefiel$Clone=="Ben Cerio 13"]<-5

Benefiel$CloneNumber[Benefiel$Clone=="Ben dent 6"]<-1
Benefiel$CloneNumber[Benefiel$Clone=="Ben dent 7"]<-2
Benefiel$CloneNumber[Benefiel$Clone=="Ben dent 16"]<-3
Benefiel$CloneNumber[Benefiel$Clone=="Ben dent 14"]<-4
Benefiel$CloneNumber[Benefiel$Clone=="Ben dent 4"]<-5

#Benefiel fig3a
fig3a<-ggplot(Benefiel, aes(x=HostType, y=FractionInf, color=MetschType))+
  ylab("Proportion infected")+
  geom_boxplot(aes(fill=MetschType), outlier.shape=NA,alpha=0.3)+
  geom_point(alpha=0.5, size=3, position=position_jitterdodge(), 
             aes(x=HostType, y=FractionInf, group=MetschType, 
                 colour=MetschType, fill=MetschType, shape=factor(CloneNumber)))+
  scale_colour_manual(values=c("firebrick","dodgerblue2"), labels=c("C. dubia", "D. dentifera"), name="Source host")+ 
  theme_classic()+ 
  scale_x_discrete(labels=c("C. dubia", "D. dentifera"))+ 
  theme(axis.text.x = element_text(face="italic"))+ 
  labs(x = "Exposed host")+
  theme(legend.text = element_text(face="italic"))+
  theme(axis.title = element_text(size=12))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(size=12))+
  theme(axis.text = element_text(size=10, color="black"))+
  theme(legend.position = "none")+
  scale_shape_manual(values=c(21:25))+
  scale_fill_manual(values=c("firebrick","dodgerblue2"), labels=c("C. dubia", "D. dentifera"), name="Source host")+
  annotate("text", x = 0.65, y = 0.53, label = "a", color="black", size=4)+
  annotate("text", x = 1.02, y = 0.325, label = "b", color="black", size=4)+
  annotate("text", x = 1.65, y = 0.41, label = "ab", color="black", size=4)+
  annotate("text", x = 2.02, y = 0.41, label = "ab", color="black", size=4)

model<-glmer(cbind(Benefiel$NumInf, Benefiel$uninf)~
               HostType*MetschType+(1|Clone), family=binomial, data=Benefiel)
drop1(model, test="Chisq")
emmeans(model, pairwise~HostType*MetschType, adjust="tukey")

#Goose Infection analysis
Goose$CloneNumber[Goose$Clone=="Goose Cerio B"]<-1
Goose$CloneNumber[Goose$Clone=="Goose Cerio A"]<-2
Goose$CloneNumber[Goose$Clone=="Goose Cerio C"]<-3
Goose$CloneNumber[Goose$Clone=="Goose Cerio I"]<-4
Goose$CloneNumber[Goose$Clone=="Goose Cerio J"]<-5

Goose$CloneNumber[Goose$Clone=="Goose dent A"]<-1
Goose$CloneNumber[Goose$Clone=="Goose dent H"]<-2
Goose$CloneNumber[Goose$Clone=="Goose dent E"]<-3
Goose$CloneNumber[Goose$Clone=="Goose dent C"]<-4
Goose$CloneNumber[Goose$Clone=="Goose dent D"]<-5

Fig4a<-ggplot(Goose, aes(x=HostType, y=FractionInf, color=MetschType))+
  ylab("Proportion infected")+
  geom_boxplot(aes(fill=MetschType), alpha=0.3, outlier.shape = NA)+
  geom_point(alpha=0.5, size=2,position=position_jitterdodge(jitter.width = 0.15), 
             aes(x=HostType, y=FractionInf, group=MetschType, 
                 colour=MetschType, fill=MetschType, shape=factor(CloneNumber)))+
  scale_colour_manual(values=c("firebrick","dodgerblue2"), labels=c("C. dubia", "D. dentifera"), name="Source Host")+ 
  scale_fill_manual(values=c("firebrick","dodgerblue2"), labels=c("C. dubia", "D. dentifera"), name="Source Host")+
  theme_classic()+ 
  scale_x_discrete(labels=c("C. dubia", "D. dentifera"))+ 
  theme(axis.text.x = element_text(face="italic"))+ 
  labs(x = "Exposed host")+
  theme(legend.text = element_text(face="italic"))+
  theme(legend.title = element_text(size=12)) +
  theme(legend.text=element_text(size=10))+
  theme(axis.title = element_text(size=12))+
  theme(axis.text = element_text(size=10, color="black"))+
  theme(legend.position = "none")+
  scale_shape_manual(values=c(21:25))

#within goose lake
#interaction not significant - nothing significant.
model<-glmer(cbind(Goose$NumInf, Goose$endNumber)~MetschType*HostType+(1|Clone)+(1|block), family=binomial, data=Goose)
drop1(model, test="Chisq")
summary(model)

#######################################
#spore counts
Benefiel.sporecounts$CloneNumber[Benefiel.sporecounts$Clone=="BeneCerio15"]<-1
Benefiel.sporecounts$CloneNumber[Benefiel.sporecounts$Clone=="BeneCerio6"]<-2
Benefiel.sporecounts$CloneNumber[Benefiel.sporecounts$Clone=="BeneCerio10"]<-3
Benefiel.sporecounts$CloneNumber[Benefiel.sporecounts$Clone=="BeneCerio1"]<-4
Benefiel.sporecounts$CloneNumber[Benefiel.sporecounts$Clone=="BeneCerio13"]<-5

Benefiel.sporecounts$CloneNumber[Benefiel.sporecounts$Clone=="BeneDent6"]<-1
Benefiel.sporecounts$CloneNumber[Benefiel.sporecounts$Clone=="BeneDent7"]<-2
Benefiel.sporecounts$CloneNumber[Benefiel.sporecounts$Clone=="BeneDent16"]<-3
Benefiel.sporecounts$CloneNumber[Benefiel.sporecounts$Clone=="BeneDent14"]<-4
Benefiel.sporecounts$CloneNumber[Benefiel.sporecounts$Clone=="BeneDent4"]<-5

fig3b<-ggplot(Benefiel.sporecounts, aes(x=Host.Species, y=real.count, color=Isolation.Species))+
  ylab((bquote('Spores/animal '~(x~10^4))))+
  geom_boxplot(aes(fill=Isolation.Species), outlier.shape=NA,alpha=0.3)+
  geom_point(alpha=0.5, size=2, position=position_jitterdodge(jitter.width = 0.2), 
             aes(x=Host.Species, y=real.count, group=Isolation.Species, 
                 colour=Isolation.Species, shape=factor(CloneNumber), fill=Isolation.Species))+
  scale_colour_manual(values=c("firebrick","dodgerblue2"), labels=c("C. dubia", "D. dentifera"), name="Source host")+ 
  scale_fill_manual(values=c("firebrick","dodgerblue2"), labels=c("C. dubia", "D. dentifera"), name="Source host")+
  theme_classic()+ 
  scale_x_discrete(labels=c("C. dubia", "D. dentifera"))+ 
  theme(axis.text.x = element_text(face="italic", size=10, color="black"))+ 
  labs(x = "Exposed host")+
  scale_y_continuous(limits=c(-200,85000), breaks=c(0,20000, 40000,60000,80000),labels=c(0,2,4,6,8))+
  theme(legend.text = element_text(face="italic", size=10))+
  theme(axis.text.y = element_text(size=10, color="black"))+
  theme(legend.title = element_text(size=12))+
  theme(axis.title = element_text(size=12))+
  scale_shape_manual(values=c(21,22,23,24,25), guide="none")+
  annotate("text", x = 0.65, y = 40800, label = "a", color="black", size=4)+
  annotate("text", x = 1.02, y = 22000, label = "b", color="black", size=4)+
  annotate("text", x = 1.65, y = 53500, label = "a", color="black", size=4)+
  annotate("text", x = 2.08, y = 49000, label = "a", color="black", size=4)

model<-lmer(real.count~Host.Species+Isolation.Species+(1|Clone)+(1|Beaker), data=Benefiel.sporecounts)
drop1(model, test="Chisq")
model.int<-lmer(real.count~Host.Species*Isolation.Species+(1|Clone)+(1|Beaker), data=Benefiel.sporecounts)
emmeans(model.int, pairwise~Host.Species*Isolation.Species, adjust="tukey")

#Goose
Goose.sporecounts$CloneNumber[Goose.sporecounts$Clone=="GooseCerioB"]<-1
Goose.sporecounts$CloneNumber[Goose.sporecounts$Clone=="GooseCerioA"]<-2
Goose.sporecounts$CloneNumber[Goose.sporecounts$Clone=="GooseCerioC"]<-3
Goose.sporecounts$CloneNumber[Goose.sporecounts$Clone=="GooseCerioI"]<-4
Goose.sporecounts$CloneNumber[Goose.sporecounts$Clone=="GooseCerioJ"]<-5

Goose.sporecounts$CloneNumber[Goose.sporecounts$Clone=="GooseDentA"]<-1
Goose.sporecounts$CloneNumber[Goose.sporecounts$Clone=="GooseDentH"]<-2
Goose.sporecounts$CloneNumber[Goose.sporecounts$Clone=="GooseDentE"]<-3
Goose.sporecounts$CloneNumber[Goose.sporecounts$Clone=="GooseDentC"]<-4
Goose.sporecounts$CloneNumber[Goose.sporecounts$Clone=="GooseDentD"]<-5

Fig4b<-ggplot(Goose.sporecounts, aes(x=Host.Species, y=real.count, color=Isolation.Species))+
  ylab(bquote("Spores/animal "~(x~10^4)))+
  geom_boxplot(aes(fill=Isolation.Species), alpha=0.3)+
  geom_point(alpha=0.5, size=2, position=position_jitterdodge(jitter.width = 0.15), 
             aes(x=Host.Species, y=real.count, group=Isolation.Species, 
                 colour=Isolation.Species, fill=Isolation.Species, shape=factor(CloneNumber)))+
  scale_colour_manual(values=c("firebrick","dodgerblue2"), labels=c("C. dubia", "D. dentifera"), name="Source host")+ 
  scale_fill_manual(values=c("firebrick","dodgerblue2"), labels=c("C. dubia", "D. dentifera"), name="Source host")+
  theme_classic()+ 
  scale_x_discrete(labels=c("C. dubia", "D. dentifera"))+ 
  theme(axis.text.x = element_text(face="italic"))+ 
  labs(x = "Exposed host")+
  theme(legend.text = element_text(face="italic")) +
  theme(legend.title=element_text(size=12))+
  theme(legend.text=element_text(size=10))+
  scale_shape_manual(values=c(21:25), guide="none")+
  theme(axis.title = element_text(size=12))+
  theme(axis.text=element_text(size=10, color="black"))+
  scale_y_continuous(limits=c(-100, 70000), breaks=c(0, 20000,40000,60000), labels=c(0,2,4,6))

#use this one
model<-lmer(real.count~Isolation.Species+Host.Species+(1|Clone)+(1|Beaker)+(1|Part), data=Goose.sporecounts)
drop1(model, test="Chisq")
#this is singular

#######################################
#Spore size
#This will take off the specific number of the measurement
spore.size$Animal.Name<-gsub("\\\\..*","", spore.size$hostname)
spore.size$Clone<-gsub( "_.*$", "", spore.size$Animal.Name )
spore.size.1<-select(spore.size, Animal.Name, length, Clone)

#Combine all the data
spore.counts<-rbind(Benefiel.sporecounts, Goose.sporecounts) #Previously loaded data on spore counts
Spores<-full_join(spore.counts, spore.size.1, by=c("Animal.Name","Clone"))
Spores<-filter(Spores, !is.na(Part))
Spores$Average<-as.numeric(as.character(Spores$Average))
SporesAndStrains<-full_join(Spores, spore.strain, by="Animal.Name")
SporesAndStrains<-filter(SporesAndStrains, !is.na(Part))
#calculate mean and ranges of number of spores measured per animal

#calculate mean spore length per individual
mean.spores<-SporesAndStrains%>%group_by(Animal.Name, Host.Lake, Parasite.Lake, Host.Species, Isolation.Species, strain, CloneNumber, Beaker, Part, Clone)%>%
  summarise(meanspore=mean(length), sd.spore=sd(length), number=length(length), meancount=mean(Average))
mean.spores$strain[is.na(mean.spores$strain)]<-"N"
mean.spores$se<-mean.spores$sd.spore/((mean.spores$number)^0.5)

Benefiel.size<-filter(mean.spores, Parasite.Lake=="Benefiel")
Goose.size<-filter(mean.spores, Parasite.Lake=="Goose")

#Fig 3c and d for Benefiel
fig3c<-ggplot(Benefiel.size, aes(x=Host.Species, y=meanspore, color=Isolation.Species))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size=2,position=position_jitterdodge(jitter.width = 0.2, seed=100), 
             aes(x=Host.Species, y=meanspore, group=Isolation.Species, 
                 colour=Isolation.Species, fill=strain, stroke=2, shape=factor(CloneNumber)))+
  ylab(expression("Mean spore length "~(mu*m)))+
  scale_colour_manual(values=c("firebrick","dodgerblue2"), 
                      labels=c("C. dubia", "D. dentifera"), 
                      name="Source host",
                      guide = guide_legend(override.aes = list(shape = 21)))+ 
  scale_fill_manual(values=c("black","dodgerblue2","firebrick","white"), labels=c("both","Daphnia-associated", "Ceriodaphnia-associated", "Not genotyped"), name="Fungal genotype", guide=FALSE)+
  theme_classic()+ 
  scale_x_discrete(labels=c("C. dubia", "D. dentifera"))+ 
  theme(axis.text.x = element_text(face="italic", size=10, color="black"))+ 
  labs(x = "Exposed host")+
  scale_shape_manual(values=c(21:25), guide="none")+
  theme(legend.text = element_text(face="italic", size=10)) +
  theme(legend.title = element_text(size=12))+
  theme(axis.text.y = element_text(size=10, color="black"))+
  theme(axis.title = element_text(size=12))+
  theme(legend.position = c(0.75, 0.2))+
  annotate("text", x = 0.65, y = 51, label = "a", color="black", size=4)+
  annotate("text", x = 1.02, y = 56.50, label = "b", color="black", size=4)+
  annotate("text", x = 1.65, y = 55.85, label = "b", color="black", size=4)+
  annotate("text", x = 2.02, y = 56.28, label = "b", color="black", size=4)

Benefiel.size$Burst<-Benefiel.size$meancount*10000*0.05
fig3d<-ggplot(Benefiel.size, aes(x=meanspore, y=Burst, group=Host.Species))+
  geom_point(size=2, aes(shape=factor(CloneNumber),fill=strain, color=Host.Species), stroke=1.3)+
  theme_classic()+
  scale_shape_manual(values=c(21:25), guide="none")+
  xlab("Mean spore length"~(mu*m))+
  ylab(bquote('Spores/animal '~(x~10^4)))+
  scale_colour_manual(values=c("firebrick","dodgerblue2"), labels=c("C. dubia", "D. dentifera"), 
                      name="Exposed host", 
                      guide = guide_legend(label.theme = element_text(angle=0, 
                                                                      face = "italic", size=10), 
                                           override.aes = list(shape = 21)))+
  theme(legend.text = element_text(face="italic", size=10))+
  theme(legend.title = element_text(size=12))+
  theme(axis.title = element_text(size=12))+
  theme(axis.text = element_text(size=10, color="black"))+
  scale_fill_manual(values=c("black","dodgerblue2","firebrick","white"), labels=c("both","D. dentifera-assoc.", "C. dubia-assoc.", "Not genotyped"), name="Fungal genotype",
                    guide = guide_legend(override.aes = list(shape = 21)))+
  scale_y_continuous(limits=c(-200, 85000), breaks=c(0,20000,40000,60000,80000), labels=c(0,2,4,6,8))+
  geom_smooth(method="lm", se=FALSE, aes(color=Host.Species))

LegTop<-get_legend(fig3b)
LegBot<-get_legend(fig3d)
Fig3<-plot_grid(fig3a,fig3b+theme(legend.position = "none"),LegTop, fig3c, fig3d+theme(legend.position = "none"),LegBot, labels = c('A', 'B','', 'C','D',''), 
                ncol=3, align="v", rel_widths = c(1,1,0.6, 1,1,0.6), axis="b")
save_plot("Figure3.jpg", Fig3, base_width = 6.53543, base_height = 6)

#used this one
model<-lmer(meanspore~Host.Species*Isolation.Species+(1|Clone)+(1|Beaker), data=Benefiel.size)
drop1(model, test="Chisq")
emmeans(model, pairwise~Host.Species*Isolation.Species, adjust="tukey")
#singular

burst1<-lmer(Burst~meanspore*Host.Species+(1|Clone)+(1|Beaker), data=Benefiel.size)
drop1(burst1, test="Chisq")

#Goose figure 4c & d
Fig4c<-ggplot(Goose.size, aes(x=Host.Species, y=meanspore, color=Isolation.Species))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(),size=2,stroke=2,
             aes(x=Host.Species, y=meanspore, group=Isolation.Species, 
                 colour=Isolation.Species, fill=strain, shape=factor(CloneNumber)))+
  ylab(expression("Mean spore length"~(mu*m)))+theme_bw()+
  theme_classic()+ 
  scale_x_discrete(labels=c("C. dubia", "D. dentifera"))+ 
  theme(axis.text.x = element_text(face="italic"))+ 
  labs(x = "Exposed host")+
  scale_shape_manual(values=c(21:25), guide="none")+
  theme(legend.text = element_text(face="italic", size=6))+
  theme(legend.title=element_text(size=8))+
  theme(axis.title = element_text(size=12))+
  theme(legend.background = element_rect(fill=alpha('white', 0)))+
  theme(axis.text = element_text(size=10, color="black"))+
  scale_colour_manual(values=c("firebrick","dodgerblue2"), 
                      labels=c("C. dubia", "D. dentifera"), 
                      name="Source host",
                      guide = guide_legend(override.aes = list(shape = 21)))+ 
  scale_fill_manual(values=c("dodgerblue2","white","black","firebrick"), labels=c("D. dentifera-assoc.","Not genotyped"), name="Fungal Genotype", guide=FALSE)+
  theme(legend.position=c(0.65,0.12))+
  guides(color = guide_legend(nrow = 1))


Goose.size$Burst<-Goose.size$meancount*10000*0.05
Fig4d<-ggplot(Goose.size, aes(x=meanspore, y=Burst))+
  geom_point(size=2,stroke=2, aes(shape=factor(CloneNumber), color=Host.Species, fill=strain))+
  theme_classic()+
  xlab(expression("Mean spore length"~(mu*m)))+
  ylab(bquote('Spores/animal '~(x~10^4)))+
  scale_colour_manual(values=c("firebrick","dodgerblue2"), 
                      labels=c("C. dubia", "D. dentifera"), 
                      name="Exposed host", 
                      guide = guide_legend(label.theme = element_text(angle=0, 
                                                                      face = "italic", size=10), 
                                                                      override.aes = list(shape = 21)))+
  theme(legend.text = element_text(face="italic"))+
  theme(legend.title = element_text(size=12))+
  theme(axis.title = element_text(size=12))+
  theme(legend.text=element_text(size=10))+
  theme(axis.text = element_text(size=10, color="black"))+
  theme(legend.text = element_text(size=10))+
  scale_fill_manual(values=c("dodgerblue2","white","black","firebrick"), labels=c("D. dentifera-assoc.","Not genotyped","",""), name="Fungal genotype",
                    guide = guide_legend(override.aes = list(shape = 21)))+
  scale_y_continuous(limits=c(-100, 70000), breaks=c(0,20000,40000,60000), labels=c(0,2,4,6))+
  scale_shape_manual(values=c(21:25),guide="none")


Legtop<-get_legend(Fig4b)
LegBottom<-get_legend(Fig4d)
Fig4<-plot_grid(Fig4a,Fig4b+theme(legend.position = "none"),Legtop, Fig4c, Fig4d+theme(legend.position = "none"),LegBottom, labels = c('A', 'B','', 'C','D',''), 
                ncol=3, align="v", rel_widths = c(1,1,0.6, 1,1,0.6), axis = "b")
save_plot("Figure4.jpg", Fig4, base_width = 6.53543, base_height = 6)

model2<-lmer(meanspore~Host.Species+Isolation.Species+(1|Clone)+(1|Part)+(1|Beaker), data=Goose.size)
drop1(model2, test="Chisq")

burst1<-lmer(Burst~meanspore+Host.Species+(1|Clone)+(1|Beaker)+(1|Part), data=Goose.size)
drop1(burst1, test="Chisq")

#Average difference in spore size
Dassoc<-filter(mean.spores, strain=="D")
Cassoc<-filter(mean.spores, strain=="IC")

mean(Dassoc$meanspore, na.rm=TRUE)
mean(Cassoc$meanspore, na.rm=TRUE)

#number of D-assoc spores in dent and cerio
DinCerio<-filter(Dassoc, Host.Species=="Cerio")
mean(DinCerio$meancount)
DinDent<-filter(Dassoc, Host.Species=="Dentifera")
mean(DinDent$meancount)
##


















