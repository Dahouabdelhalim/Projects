
#load packages (assuming they are downloaded already)
library(car) #to get pvalues from model
library(Rmisc) #to get basic stats
library(FSA) #to get basic stats
library(emmeans) #pairwise analysis
library(ggplot2)#build graphs
library(Hmisc)
library(splines)
library(sandwich)
library(effects)
library(RcmdrMisc)

#read data.xlsx from directory
data <- readXL("data.xlsx", rownames=FALSE, header=TRUE, na="", sheet="Compare ws and ds live.4grps", stringsAsFactors=TRUE)

#Setting up the data & basic calculations
attach(data)
head(data)
summary(data) 
sum<-summarySE(data=data, "total.dur", groupvars=NULL, na.rm=FALSE, conf.interval=0.95, .drop=TRUE)
sum
detach(data)

#### determine extreme outliers and plot by setting y limits to the most acceptable extreme value ####
#take note of up and low, to set y limits in plots, don't remove from actual data since outliers are true measurements, not due to instrument error or etc
q<-quantile(data$total.dur, probs=c(.25,.75), na.rm=FALSE)
q
iqr<-IQR(data$total.dur)
up<-q[2]+3*iqr
low<-q[1]-3*iqr

############# courtship duration 4 grps ################
names(data)
data$treatmentxseason1<-factor(data$treatmentxseason,levels=c("WT.ws", "yellow.ws", "WT.ds", "yellow.ds"),ordered=TRUE)
names(data)

sum<-summarySE(data=data, measurevar="total.dur", groupvars="treatmentxseason1", na.rm=FALSE, conf.interval=0.95, .drop=TRUE)
sum


plot<-ggplot(data,aes(x=treatmentxseason1,y=total.dur,color=treatmentxseason1))+
geom_jitter(position=position_jitter(0.1), size=2, shape=1, color="grey25")+
geom_pointrange(aes(ymin=total.dur-se, ymax=total.dur+se),data=sum, size=0.7)+
geom_errorbar(aes(ymin=total.dur-se, ymax=total.dur+se),data=sum, width=0.3, size = 1.1)+
theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)
        )+
scale_color_manual(values=c("#7E550E", "#E69F00", "#7E550E", "#E69F00"))+
scale_x_discrete(name="Treatment", labels=c("Wt (WS)", "Yellow (WS)", "Wt (DS)", "Yellow (DS)")) +
scale_y_continuous(name="Courtship Duration (s)", limits=c(0, 209))

plot
## remember to change name of output to decap or live accordingly
png("output.png",units="in",width=8, height=5, res=1280, bg="transparent")
plot
dev.off()

######### courtship frequency 4 grps #######

names(data)
###dont need the below if treatmentxseason1 is already done before/ established ###
data$treatmentxseason1<-factor(data$treatmentxseason,levels=c("WT.ws", "yellow.ws", "WT.ds", "yellow.ds"),ordered=TRUE)
names(data)

sum<-summarySE(data=data, measurevar="total.freq", groupvars="treatmentxseason1", na.rm=FALSE, conf.interval=0.95, .drop=TRUE)
sum

plot<-ggplot(data,aes(x=treatmentxseason1,y=total.freq,color=treatmentxseason1))+
geom_jitter(position=position_jitter(0.1),size=2, shape=1, color="grey25")+
geom_pointrange(aes(ymin=total.freq-se, ymax=total.freq+se),data=sum, size=0.7)+
geom_errorbar(aes(ymin=total.freq-se, ymax=total.freq+se),data=sum, width=0.3, size=1.1)+
theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)
        )+
scale_color_manual(values=c("#7E550E", "#E69F00", "#7E550E", "#E69F00"))+
scale_x_discrete(name="Treatment", labels=c("Wt (WS)", "Yellow (WS)", "Wt (DS)", "Yellow (DS)")) +
scale_y_continuous(name="Courtship Frequency", limits=c(0, 24))

plot

png("output.png",units="in",width=8, height=5, res=1280, bg="transparent")
plot
dev.off()
########## Mating duration 4 grps #############

names(data)
### needs to change treatmentxseason to mating specific column + duration.mated + level names #
data$treatment.mated1<-factor(data$treatment.mated,levels=c("WT.ws", "yellow.ws", "WT.ds", "yellow.ds"),ordered=TRUE)
names(data)

sum<-summarySE(data=data, measurevar="duration.mated", groupvars="treatment.mated1", na.rm=TRUE, conf.interval=0.95, .drop=TRUE)
sum1<-na.omit(sum)
sum1
xlim<-range(data$treatment.mated1, na.rm=TRUE)

plot<-ggplot(data,aes(x=treatment.mated1,y=duration.mated,color=treatment.mated1, ylim=120))+
scale_x_discrete(expand=c(0,0), limits=factor ("Wt (WS)", "Yellow (WS)", "Wt (DS", "Yellow (DS)"))+
scale_y_continuous(expand=c(0,0), limits=c(0, 80))+
expand_limits(x=0, y=0)+
expand_limits(x=c(0, NA), y=c(0,80))+
geom_jitter(position=position_jitter(0.1),size=2, shape=1, color="grey25", na.rm=TRUE)+
geom_pointrange(aes(ymin=duration.mated-se, ymax=duration.mated+se),data=sum1, size=0.7)+
geom_errorbar(aes(ymin=duration.mated-se, ymax=duration.mated+se),data=sum1, width=0.3, size=1.1)+
theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)
        )+
scale_color_manual(values=c("#7E550E", "#E69F00", "#7E550E", "#E69F00"))+
scale_x_discrete(name="Treatment", labels=c("Wt (WS)", "Yellow (WS)", "Wt (DS)", "Yellow (DS)")) +
scale_y_continuous(name="Mating Duration (min)", limits = c(0,80))

plot

png("output.png",units="in",width=8, height=5, res=1280, bg="transparent")
plot
dev.off()


########### Courtship latency 4 grps ############

names(data)
### needs to change treatmentxseason to mating specific column + latency.mated + level names ##########
data$treatment.mated1<-factor(data$treatment.mated,levels=c("WT.ws", "yellow.ws", "WT.ds", "yellow.ds"),ordered=TRUE)
names(data)

sum<-summarySE(data=data, measurevar="latency.mated", groupvars="treatment.mated1", na.rm=FALSE, conf.interval=0.95, .drop=TRUE)
sum1<-na.omit(sum)
sum1
xlim<-range(data$treatment.mated1, na.rm=TRUE)

plot<-ggplot(data,aes(x=treatment.mated1,y=latency.mated,color=treatment.mated1, ylim=64))+
expand_limits(x=0, y=0)+
expand_limits(x=c(0, NA), y=c(0,65))+
geom_jitter(position=position_jitter(0.1),size=2, shape=1, color="grey25", na.rm=TRUE)+
geom_pointrange(aes(ymin=latency.mated-se, ymax=latency.mated+se),data=sum1, size=0.7)+
geom_errorbar(aes(ymin=latency.mated-se, ymax=latency.mated+se),data=sum1, width=0.3)+
theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)
        )+
scale_color_manual(values=c("#7E550E", "#E69F00", "#7E550E", "#E69F00"))+
scale_x_discrete(name="Treatment", labels=c("Wt (WS)", "Yellow (WS)", "Wt (DS)", "Yellow (DS)")) +
scale_y_continuous(name="Courtship Latency (min)")

plot

png("output.png",units="in",width=8, height=5, res=1280, bg="transparent")
plot
dev.off()

############ Mating status 4 grps ###################

names(data)
data$treatmentxseason1<-factor(data$treatmentxseason,levels=c("WT.ws", "yellow.ws", "WT.ds", "yellow.ds"),ordered=TRUE)
names(data)
data$mating.status.chi1<-factor(data$mating.status.chi, levels=c("Mated", "Not mated"))

plot<-with(data, Barplot(mating.status.chi1, by=treatmentxseason1, style="parallel", xlab="Treatment", ylab="Percentage of Matings", 
scale="percent", col=c("#7E550E", "#E69F00", "#7E550E", "#E69F00"),ylim=c(0,60)))


