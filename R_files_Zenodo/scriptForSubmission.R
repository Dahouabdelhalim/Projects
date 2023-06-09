####################
## Load libraries ##
####################
library(plyr)
library(gmodels)
library(ggplot2)
library(scales)
library(reshape)
library(grid)
library(MASS)
library(lubridate)
library(multcomp)
library(Hmisc)
library(gtable)
library(xtable)
library(pander)
library(mgcv)
library(gridExtra)

##############################################################
## Define an inverse logit function for back-transformation ##
##############################################################
inverselogit <- function(x) {
 exp(x)/(exp(x)+1)
}


#######################
##  Read in the data ##
#######################

##,---------------------
##| Coral bleaching data
##`---------------------
data <- read.csv('data/predicting bleaching and mortality from zoox type v4_LB.csv', strip.white=TRUE)
head(data)
#visual <- read.csv('data/visual bleaching data_LB.csv', strip.white=TRUE)
visual <- read.csv('data/visual-bleaching-data_LB.csv', strip.white=TRUE)
head(visual)

##,-----------------
##| Temperature data
##`-----------------
temperature <- read.csv('data/Temperatures.csv', strip.white=TRUE)
head(temperature)
davies.thres <- read.csv('data/DaviesThresholds.csv', strip.white=TRUE)
davies.thres

keppels.thres <- read.csv('data/KeppelsThresholds.csv', strip.white=TRUE)
keppels.thres

maggie.thres <- read.csv('data/MaggieThresholds.csv', strip.white=TRUE)
maggie.thres


#####################
## Data processing ##
#####################

## 1. Recode `Population` as a factor (categorical predictor) and Time point as a factor and actual date
data$Population <- factor(data$Population, labels=c('Keppels','Davies','Magnetic Island'))
data$Time.point <- factor(data$Time.point)
data$Date <- as.Date(as.character(factor(data$Time.point, labels=c('2005-04-05','2005-09-14','2005-12-21','2006-01-16','2006-02-14','2006-04-30'))))
head(data)

visual$Population <- factor(visual$Population, levels=c('Keppels','Davies','Maggie'),labels=c('Keppels','Davies','Magnetic Island'))
visual$Time.point <- factor(visual$Time)
visual$Date <- as.Date(as.character(factor(visual$Time.point, labels=c('2005-04-05','2005-09-14','2005-12-21','2006-01-16','2006-02-14','2006-04-30'))))
visual$sample <- factor(visual$Sample)
head(visual)

## 2. Relabel  the S.H..C.D. column to D.Cratio for greater readibility
nms <- colnames(data)
nms[nms=='S.H..C.D.'] <- "D.Cratio"
colnames(data) <- nms

## 3. Add a column `S.H..C.D.` to represent the S:H ratio for both C and D clades
data$S.H..C.D. <- data$S.H..clade.C. + data$S.H..clade.D.

## 4. Some (five) of the transplanted Magnetic Island colonies were of sub-optimal starting health.
##    Five of the colonies had percentage mortality values greater than 0 at the first sampling period.  There inclusion (whilst not ideal)
##    is likely to have been to maintain sample sizes (which were already lower for Magnetic Island).
##
##    Nevertheless, as these were potentially compromised prior to the start of the experiment, there is good justification to exclude them as
##    atypical representatives.  I have therefore excluded all samples that had a starting mortality % greater than 0.
dim(data) 
sub <- data$Mortality>0 & data$Time.point==1  
##work out which sample it is and exclude this sample from all times
samples <- data[sub,'sample'] 
data <- subset(data,!(sample %in% samples))
dim(data)

## 5. Make alternative versions of the data to use for different analyses
a <- melt(subset(data, select=c('Population','Time.point', 'sample', 'S.H..C.D.',
                         'D.C.ratio..20bin.','Date')),id=c('Population','Time.point','sample','Date'))
data1 <- cast(a, Population+sample~variable+Date)
data1$SH_Dec <- data1$'S.H..C.D._2005-12-21'
data1$DC_Sep <- data1$'D.C.ratio..20bin._2005-09-14'
data1$DC_Apr <- data1$'D.C.ratio..20bin._2005-04-05'
head(data1)

a <- melt(subset(data, select=c('Population','Time.point', 'sample', 'Mortality',
                         'D.C.ratio..20bin.','Date')),id=c('Population','Time.point','sample','Date'))
data2 <- cast(a, Population+sample~variable+Date)
data2$Mortality_Apr <- data2$'Mortality_2006-04-30'
data2$DC_Sep <- data2$'D.C.ratio..20bin._2005-09-14'
head(data2)

a <- melt(subset(data, select=c('Population','Time.point', 'sample', 'Mortality','Upgraded.D','Unchanged','Downgraded.D','Date')),
          id=c('Population','Time.point','sample','Date'))
data3 <- cast(a, Population+sample~variable+Date)
data3$Mortality_Apr <- data3$'Mortality_2006-04-30'
data3$Upgraded_Sep <- data3$'Upgraded.D_2005-09-14'
data3$Unchanged_Sep <- data3$'Unchanged_2005-09-14'
data3$Downgraded_Sep <- data3$'Downgraded.D_2005-09-14'
head(data3)

a <- melt(subset(data, select=c('Population','Time.point', 'sample', 'S.H..clade.D.',
                         'D.C.ratio..20bin.','Date')),id=c('Population','Time.point','sample','Date'))
data4 <- cast(a, Population+sample~variable+Date)
data4$SH_Dec <- data4$'S.H..clade.D._2005-12-21'
data4$DC_Sep <- data4$'D.C.ratio..20bin._2005-09-14'
data4$DC_Apr <- data4$'D.C.ratio..20bin._2005-04-05'
head(data4)


## 6. Calculate the differences in temperatures experienced by each source from their typical long-term average.
##    - For the Magnetic Island sourced corals, this means subtracting the temperature of Geoffrey Bay from the long-term Geoffrey Bay average.
##    - For the Keppels sourced corals, this means subtracting the temperature of Geoffrey Bay from the long-term Keppels average.
##    - For the Davies Reef sourced corals, this means subtracting the temperature of Geoffrey Bay from the long-term Davies Reef average.
temperature$KeppelDiff <- temperature$Geoffrey.Bay.2005.06 - temperature$AvHalfFl
temperature$DavisDiff <- temperature$Geoffrey.Bay.2005.06 - temperature$AvDaviesFl
temperature$MaggieDiff <- temperature$Geoffrey.Bay.2005.06 - temperature$AvNelFl
temperature$Date <- as.Date(temperature$Experimental.date, format='%d/%m/%Y')


## 7. Calculate accumulated bleaching threshold exceedences over time for each population
A <- temperature$Geoffrey.Bay.2005.06
B <- keppels.thres$Av.daily.temp

#go with a 60 day window
AA <- A
Exp <- NULL
days <- 60
for (i in 1:length(A)) {
 d <- ifelse((i-60)<1,1,i-60)
 A<-AA[d:i]
 B <- keppels.thres$Av.daily.temp
 CC <- keppels.thres
 BB<-matrix(rep(B,length(A)), byrow=TRUE,nrow=length(A))
 a<-matrix(A>BB,nrow=length(A))
 freq<-colSums(a)
 B <- cbind(CC,f=freq)
 Exp <- c(Exp,(sum((B$f-B[,2])*(B$f>B[,2]), na.rm=TRUE)))
}
dt.k <- data.frame(temperature,Exp)

dt.exp <- dt.k$Date[which(dt.k$Exp==min(subset(dt.k, Exp>0)$Exp))[1]]

A <- temperature$Geoffrey.Bay.2005.06
B <- davies.thres$Av.daily.temp
CC <- davies.thres

#go with a 60 day window
AA <- A
Exp <- NULL
days <- 60
for (i in 1:length(A)) {
 d <- ifelse((i-60)<1,1,i-60)
 A<-AA[d:i]
 BB<-matrix(rep(B,length(A)), byrow=TRUE,nrow=length(A))
 a<-matrix(A>BB,nrow=length(A))
 freq<-colSums(a)
 temp <- cbind(CC,f=freq)
 Exp <- c(Exp,(sum((temp$f-temp[,2])*(temp$f>temp[,2]), na.rm=TRUE)))
}
dt.d <- data.frame(temperature,Exp)

dt.exp <- c(dt.exp,dt.d$Date[which(dt.d$Exp==min(subset(dt.d, Exp>0)$Exp))[1]])


A <- temperature$Geoffrey.Bay.2005.06
B <- maggie.thres$Av.daily.temp
CC <- maggie.thres

#go with a 60 day window
AA <- A
Exp <- NULL
days <- 60
for (i in 1:length(A)) {
 d <- ifelse((i-60)<1,1,i-60)
 A<-AA[d:i]
 BB<-matrix(rep(B,length(A)), byrow=TRUE,nrow=length(A))
 a<-matrix(A>BB,nrow=length(A))
 freq<-colSums(a)
 temp <- cbind(CC,f=freq)
 Exp <- c(Exp,(sum((temp$f-temp[,2])*(temp$f>temp[,2]), na.rm=TRUE)))
}
dt.m <- data.frame(temperature,Exp)

dt.exp <- c(dt.exp,dt.m$Date[which(dt.m$Exp==min(subset(dt.m, Exp>0)$Exp))[1]])
dt.exp <- data.frame(Population=c("Keppels","Davies","Magnetic Island"), Date=format(as.Date(dt.exp),'%d/%m/%Y'))

# Now for mortality
A <- temperature$Geoffrey.Bay.2005.06
B <- keppels.thres$Av.daily.temp

#go with a 60 day window
AA <- A
Exp <- NULL
days <- 60
for (i in 1:length(A)) {
 d <- ifelse((i-60)<1,1,i-60)
 A<-AA[d:i]
 B <- keppels.thres$Av.daily.temp
 CC <- keppels.thres
 BB<-matrix(rep(B,length(A)), byrow=TRUE,nrow=length(A))
 a<-matrix(A>BB,nrow=length(A))
 freq<-colSums(a)
 B <- cbind(CC,f=freq)
 Exp <- c(Exp,(sum((B$f-B[,3])*(B$f>B[,3]), na.rm=TRUE)))
}
mort.k <- data.frame(temperature,Exp)

mort.exp <- mort.k$Date[which(mort.k$Exp==min(subset(mort.k, Exp>0)$Exp))[1]]

A <- temperature$Geoffrey.Bay.2005.06
B <- davies.thres$Av.daily.temp
CC <- davies.thres

#go with a 60 day window
AA <- A
Exp <- NULL
days <- 60
for (i in 1:length(A)) {
 d <- ifelse((i-60)<1,1,i-60)
 A<-AA[d:i]
 BB<-matrix(rep(B,length(A)), byrow=TRUE,nrow=length(A))
 a<-matrix(A>BB,nrow=length(A))
 freq<-colSums(a)
 temp <- cbind(CC,f=freq)
 Exp <- c(Exp,(sum((temp$f-temp[,3])*(temp$f>temp[,3]), na.rm=TRUE)))
}
mort.d <- data.frame(temperature,Exp)

mort.exp <- c(mort.exp,mort.d$Date[which(mort.d$Exp==min(subset(mort.d, Exp>0)$Exp))[1]])

mort.exp <- c(mort.exp,NA)
mort.exp <- data.frame(Population=c("Keppels","Davies","Magnetic Island"), Date=format(as.Date(mort.exp),'%d/%m/%Y'))

dt.exp <- cbind(dt.exp,mort.exp[,2])
colnames(dt.exp) <- c('Population', 'Bleaching','Mortality')

########################################################################
## Summary table of when bleaching tolerance conditions were exceeded ##
## for each population                                                ##
########################################################################
pandoc.table(dt.exp)
print(xtable(dt.exp, caption = 'Dates at which bleaching tolerance conditions were exceeded for each population'),
      comment=FALSE,include.rownames=FALSE, booktabs=TRUE,caption.placement='top')


###############################
## Exploratory data analysis ##
###############################

##,--------------------------------------------------------------------------------------------------------------------------------------------
##| Temperature plot Figure 1. ================================================================================================================
##`--------------------------------------------------------------------------------------------------------------------------------------------
dt <- melt(temperature, id=c("Date"), measure.vars=c('KeppelDiff','DavisDiff','MaggieDiff'), variable_name='Population')
ymin <- min(dt$value)-0.5
ymax <- max(dt$value)
yticks <- -2:6
tmp1<-ggplot(dt, aes(y=value, x=Date, fill=Population)) + geom_ribbon(aes(ymin=0, ymax=value), alpha=.50)+
  geom_line(aes(col=Population), size=0.25)+
  coord_cartesian(ylim=c(ymin,9))+
  scale_y_continuous(expression(Temperature~difference~(degree~C)), limits=c(ymin,ymax), breaks=yticks)+
  scale_x_date("Date", expand=c(0,0))+
  scale_color_manual('',guide=FALSE,values=c('grey75','grey45','grey15'),breaks=c('KeppelDiff','DavisDiff','MaggieDiff'),lab=c("Keppels","Davies","Magnetic Island"))+
#  scale_fill_manual('',values=c('#ff000010','#00ff0010','#0000ff10'),breaks=c('KeppelDiff','DavisDiff','MaggieDiff'),lab=c("Keppels","Davies","Magnetic Island"))+
  scale_fill_manual('',values=c('grey80','grey50','grey20'),breaks=c('KeppelDiff','DavisDiff','MaggieDiff'),lab=c("Keppels","Davies","Magnetic Island"))+
  theme(panel.grid.major = element_blank(), # no major grid lines
  panel.grid.minor = element_blank(), # no minor grid lines
  panel.background = element_blank(), # no background
  panel.border = element_blank(), # no plot border
  axis.title.y=element_text(size=11, vjust=1,angle=90, hjust=0), # y-axis title
  axis.text.y=element_text(size=9), # y-axis labels
  axis.title.x=element_blank(),#element_text(size=11, vjust=-1), # x-axis title
  axis.text.x=element_text(size=9), # x-axis labels
  axis.line = element_line(),
  legend.position=c(1,0),
  legend.justification=c(0.95,0.2),
  legend.background=element_blank(),
  legend.direction="horizontal",
  plot.margin=unit(c(0.5,0.5,1,1),"lines"))#,text=element_text(family='Times',size=11))

g1 <- ggplot_gtable(ggplot_build(tmp1))

dtt <- temperature
dtt$Temp <- rescale(dtt$Geoffrey.Bay.2005.06, from=c(21,33),to=c(6,9))
yticks.2 <- rescale(seq(20,32,b=4), from=c(21,33),to=c(6,9))
yticksLab.2 <- seq(20,33,b=4)
ats <- rescale(seq(20,33,b=1), from=c(21,33),to=c(6,9))
tmp2<-ggplot(dtt, aes(y=Temp, x=Date)) +
    geom_hline(yintercept=ats, color='grey', size=0.25)+
    geom_line()+
    coord_cartesian(ylim=c(ymin,9))+
  scale_y_continuous(expression(Actual~temperature~(degree~C)), breaks=yticks.2, labels=yticksLab.2)+
  scale_x_date("Date", expand=c(0,0))+
  theme(panel.grid.major = element_blank(), # no major grid lines
  panel.grid.minor = element_blank(), # no minor grid lines
  panel.background = element_blank(), # no background
  panel.border = element_blank(), # no plot border
  axis.title.y=element_text(size=11,angle=90, debug=FALSE, hjust=1), # y-axis title
  axis.text.y=element_text(size=9), # y-axis labels
  axis.title.x=element_blank(),#element_text(size=11, vjust=-1), # x-axis title
  axis.text.x=element_text(size=9), # x-axis labels
  axis.line = element_line(),
  legend.position=c(0.5,0),
  legend.justification=c(0,0.2),
  legend.background=element_blank(),
  legend.direction="horizontal",
  plot.margin=unit(c(0.5,0.5,1,1),"lines"))#,text=element_text(family='Times',size=11))

tmp2
g2 <- ggplot_gtable(ggplot_build(tmp2))

combo_grob <- g1

panel_grob <- getGrob(g2$grobs[[4]], 'GRID.segments',grep = TRUE, global = TRUE)
combo_grob$grobs[[4]] <- addGrob(combo_grob$grobs[[4]], panel_grob)

panel_grob <- getGrob(g2$grobs[[4]], 'GRID.polyline',grep = TRUE, global = TRUE)
combo_grob$grobs[[4]] <- addGrob(combo_grob$grobs[[4]], panel_grob)

pos_a <- grep('axis-l', g2$layout$name)
axis <- g2$grobs[[pos_a]]
pp <- c(subset(g2$layout, name == 'panel', se = t:r))
## ax <- axis$children[[1]]
## #ax$widths <- rev(ax$widths)
## #ax$grobs <- rev(ax$grobs)
## ax$x <- ax$x - unit(1, "npc")# + unit(0.25, "cm")
## #ax$grobs[[2]]$x <- ax$grobs[[2]]$x - unit(1, "npc") + unit(0.75, "cm")
## #library(gtable)
## combo_grob <- gtable_add_grob(combo_grob, ax,  pp$t, length(combo_grob$widths) - 0, pp$b)

#axis line
ax <- axis$children[[1]]
ax$x <- ax$x - unit(1, "npc") + unit(0.00, "cm")
combo_grob <- gtable_add_grob(combo_grob, ax,  pp$t, length(combo_grob$widths) - 0, pp$b)

#ticks
ax <- axis$children[[2]]
ax$widths <- rev(ax$widths)
ax$grobs <- rev(ax$grobs)
ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") - unit(0.1, "cm")
#ax$grobs[[2]]$x <- ax$grobs[[2]]$x - unit(1, "npc") + unit(0.0, "cm") 
library(gtable)
combo_grob <- gtable_add_cols(combo_grob, g2$widths[g2$layout[pos_a,]$l], length(combo_grob$widths) - 0)
combo_grob <- gtable_add_grob(combo_grob, ax,  pp$t, length(combo_grob$widths) - 0, pp$b)

pp <- c(subset(g2$layout, name == 'ylab', se = t:r))
ia <- which(g2$layout$name == "ylab")
ga <- g2$grobs[[ia]]
ga$rot <- 270
#ga$y <- ga$y + unit(1.0, "npc") #+ unit(1.5, "cm")
#ga$x <- ga$x + unit(1.5, "lines")# + unit(1.5, "cm")
combo_grob <- gtable_add_cols(combo_grob, 1*g2$widths[g2$layout[ia,]$l])
combo_grob <- gtable_add_grob(combo_grob, ga, pp$t, length(combo_grob$widths), pp$b)

combo_grob$layout$clip <- "off"

grid.newpage()
grid.draw(combo_grob)
##External graphical formats
pdf(file='figures/Figure1.pdf', width=7, height=2.5)
grid.newpage()
grid.draw(combo_grob)
dev.off()

system('pdf2ps -dLanguageLevel=3 figures/Figure1.pdf figures/Figure1.eps')
system('pdftops -level3 -eps figures/Figure1.pdf figures/Figure1a.eps')
#system('convert figures/Figure1.pdf figures/Figure1.eps')

postscript(file='figures/Figure1p.eps', width=7, height=2.5, horizontal=FALSE, paper='special', fonts='Times')
grid.newpage()
grid.draw(combo_grob)
dev.off()

library(Cairo)
CairoFonts('Nimbus Roman No9 L')
CairoPS(file='figures/Figure1.eps', width=7, height=2.5, horizontal=FALSE, paper='special',pointsize=8)
grid.newpage()
grid.draw(combo_grob)
dev.off()

cairo_ps(file='figures/Figure1.eps', width=7, height=2.5)
grid.newpage()
grid.draw(combo_grob)
dev.off()
##,-----------------------------------------------------------------------------------------------------------------------------------
##| ==================================================================================================================================
##`-----------------------------------------------------------------------------------------------------------------------------------


##,-------------------
##| Visual scale model
##`-------------------
visual$Dt <- decimal_date(visual$Date)
visual$Vis <- visual$Color/5

aa.k <- gamm(Vis ~ s(Dt,k=4), random=list(sample=~1),data=subset(visual, Population=='Keppels'), family=binomial(logit), na.action=na.omit, correlation=corAR1())
#aa.k <- gam(Vis ~ s(Dt, k=4), random=list(sample=~1),data=subset(visual, Population=='Keppels'), family=binomial(logit), na.action=na.omit, correlation=corAR1())
newdata.k<-expand.grid(Population='Keppels', Dt=seq(min(visual$Dt), max(visual$Dt), l=100))
newdata.k$Date <- date_decimal(newdata.k$Dt)
fit=predict(aa.k$gam, newdata=newdata.k, type='link', se=TRUE)
newdata.k$fit=binomial()$linkinv(fit$fit)
newdata.k$lower=binomial()$linkinv(fit$fit-2*fit$se.fit)
newdata.k$upper=binomial()$linkinv(fit$fit+2*fit$se.fit)

aa.d <- gamm(Vis ~ s(Dt, k=4), random=list(sample=~1),data=subset(visual, Population=='Davies'), family=binomial(logit), na.action=na.omit, correlation=corAR1())
newdata.d<-expand.grid(Population='Davies', Dt=seq(min(visual$Dt), max(visual$Dt), l=100))
newdata.d$Date <- date_decimal(newdata.d$Dt)
fit=predict(aa.d$gam, newdata=newdata.d, type='link', se=TRUE)
newdata.d$fit=binomial()$linkinv(fit$fit)
newdata.d$lower=binomial()$linkinv(fit$fit-2*fit$se.fit)
newdata.d$upper=binomial()$linkinv(fit$fit+2*fit$se.fit)

aa.m <- gamm(Vis ~ s(Dt, k=4), random=list(sample=~1),data=subset(visual, Population=='Magnetic Island'), family=binomial(logit), na.action=na.omit, correlation=corAR1())
newdata.m<-expand.grid(Population='Magnetic Island', Dt=seq(min(visual$Dt), max(visual$Dt), l=100))
newdata.m$Date <- date_decimal(newdata.m$Dt)
fit=predict(aa.m$gam, newdata=newdata.m, type='link', se=TRUE)
newdata.m$fit=binomial()$linkinv(fit$fit)
newdata.m$lower=binomial()$linkinv(fit$fit-2*fit$se.fit)
newdata.m$upper=binomial()$linkinv(fit$fit+2*fit$se.fit)

newdata=rbind(newdata.k, newdata.d, newdata.m)

dts <- floor_date(unique(data$Date),'month')[-4]+15

p <- ggplot(visual, aes(Vis, x=Dt)) +
    geom_ribbon(data=newdata, aes(y=fit,ymin=lower, ymax=upper), fill='grey80', color=NA)+
        geom_line(data=newdata, aes(y=fit)) + 
    #stat_smooth(color='black',method='gam', formula=y~s(x, k=4), method.args=list(family="binomial")) +
    geom_point(position=position_jitter(width=0.02,height=0), colour="black", size=2)+facet_grid(~Population)+
    scale_y_continuous("Visual scale", breaks=c(1,2,3,4,5)/5, labels=c(1,2,3,4,5))+
    scale_x_continuous("Date", breaks=decimal_date(dts), labels=c(format(dts[1],format="%b\\n%Y"),strtrim(format(dts[c(-1,-5)],format="%b"),1),format(dts[5],format="%b\\n%Y")))+
      ggtitle('a)')+
  theme(axis.text.x=element_text(size = 9),
        panel.grid.major = element_blank(), # no major grid lines
        panel.grid.minor = element_blank(), # no minor grid lines
        panel.background = element_blank(), # no background
        panel.border = element_blank(), # no plot border
        axis.title.y=element_text(size=11, margin=margin(r=2, unit='line'),angle=90), # y-axis title
        axis.text.y=element_text(size=9), # y-axis labels
        axis.title.x=element_text(size=11, vjust=-2), # x-axis title
        axis.line = element_line(),
        legend.position=c(1,0),
        legend.justification=c(1,0),
        strip.background=element_blank(),
        strip.text=element_text(size=rel(1.0)),
        plot.margin=unit(c(0.5,0.5,2,2),"lines"),
        plot.title=element_text(hjust=0,vjust=0),
        panel.margin = unit(1, "lines"))#, text=element_text(family='Times'))

p

library(gtable)
vis.grob <- ggplotGrob(p)
grid.newpage()
grid.draw(rbind(vis.grob[1],vis.grob[4:5], vis.grob[1],vis.grob[3],size="first"))

##,----------------
##| Mortality model
##`----------------
data$Mort <- data$Mortality/100
data$Dt <- decimal_date(data$Date)
a<-glmmPQL(Mort ~ Dt*Population, random=~1|sample,data, family=binomial(logit), na.action=na.omit, correlation=corAR1())
newdata<-expand.grid(Population=levels(data$Population), Dt=seq(min(data$Dt), max(data$Dt), l=100))
newdata$Date <- date_decimal(newdata$Dt)
coefs <- fixef(a)
mm <- model.matrix(~ Dt*Population, newdata)
#newdata$fit <- inverselogit(as.vector(coefs %*% t(mm)))
newdata$fit <- binomial(logit)$linkinv(as.vector(coefs %*% t(mm)))*100
predvar <- diag(mm %*% vcov(a) %*% t(mm)) 
newdata$SE <- sqrt(predvar) 
newdata$SE2 <- sqrt(predvar+a$sigma^2)
newdata$lwr <- inverselogit(as.vector(coefs %*% t(mm))-newdata$SE2)*100
newdata$upr <- inverselogit(as.vector(coefs %*% t(mm))+newdata$SE2)*100
dts <- floor_date(unique(data$Date),'month')[-4]+15

p<-ggplot(newdata, aes(y=fit, x=Dt)) + 
  geom_point(data=data,aes(y=Mortality, x=Dt),position=position_jitter(0.02), colour="black", size=2)+
  geom_line(aes(x=Dt)) + facet_grid(~Population)+geom_ribbon(aes(ymin=lwr, ymax=upr, x=Dt), fill="grey", alpha=0.5)+
  scale_y_continuous("Mortality (%)")+
  #scale_x_continuous("Date", breaks=decimal_date(as.Date(pretty(newdata$Date))), labels=format(as.Date(pretty(newdata$Date)),format="%b\\n%Y"))+
  #scale_x_continuous("Date", breaks=decimal_date(dts), labels=c(format(dts[1],format="%b\\n%Y"),format(dts[c(-1,-5)],format="%b"),format(dts[5],format="%b\\n%Y")))+
  scale_x_continuous("Date", breaks=decimal_date(dts), labels=c(format(dts[1],format="%b\\n%Y"),strtrim(format(dts[c(-1,-5)],format="%b"),1),format(dts[5],format="%b\\n%Y")))+
  ggtitle('b)')+
  theme(axis.text.x=element_text(size = 9),
        panel.grid.major = element_blank(), # no major grid lines
        panel.grid.minor = element_blank(), # no minor grid lines
        panel.background = element_blank(), # no background
        panel.border = element_blank(), # no plot border
        axis.title.y=element_text(size=11, margin=margin(r=2, unit='line'),angle=90), # y-axis title
        axis.text.y=element_text(size=9), # y-axis labels
        axis.title.x=element_text(size=11), # x-axis title
        axis.line = element_line(),
        legend.position=c(1,0),
        legend.justification=c(1,0),
        strip.background=element_blank(),
        strip.text=element_text(size=rel(1.0)),
        plot.margin=unit(c(0.5,0.5,2,2),"lines"),
        plot.title=element_text(hjust=0,vjust=0),
        panel.margin = unit(1, "lines"))#, text=element_text(family='Times'))

library(gtable)
mort.grob <- ggplotGrob(p)
grid.newpage()
grid.draw(rbind(mort.grob[1],mort.grob[4:5], mort.grob[1],mort.grob[3],size="first"))

##,---------------------------------------------------------------------------------------------------------------------
##| Figure 2 ===========================================================================================================
##`---------------------------------------------------------------------------------------------------------------------
grid.newpage()
grid.draw(rbind(
vis.grob[1:2],vis.grob[4:5],    
mort.grob[1:2],mort.grob[4:5],mort.grob[3],size="first"))

## External versions
pdf(file='figures/Figure2.pdf', width=175/25.4, height=100/25.4)
grid.newpage()
grid.draw(rbind(
    vis.grob[1:2],vis.grob[4:5],    
    mort.grob[1:2],mort.grob[4:5],mort.grob[3],
    size="first"))
dev.off()

#system('pdf2ps -dLanguageLevel=3 figures/Figure2.pdf figures/Figure2.eps')
system('pdftops -level3 -eps figures/Figure2.pdf figures/Figure2a.eps')
#system('convert figures/Figure1.pdf figures/Figure1.eps')

postscript(file='figures/Figure2p.eps', width=175/25.4, height=100/25.4, horizontal=FALSE, paper='special', fonts='Times')
grid.newpage()
grid.draw(rbind(
    vis.grob[1:2],vis.grob[4:5],    
    mort.grob[1:2],mort.grob[4:5],mort.grob[3],
    size="first"))
dev.off()

library(Cairo)
CairoFonts('Nimbus Roman No9 L')
CairoPS(file='figures/Figure2.eps', width=175/25.4, height=100/25.4, horizontal=FALSE, paper='special',pointsize=8)
grid.newpage()
grid.draw(rbind(
    vis.grob[1:2],vis.grob[4:5],    
    mort.grob[1:2],mort.grob[4:5],mort.grob[3],
    size="first"))
dev.off()

cairo_ps(file='figures/Figure2.eps', width=175/25.4, height=100/25.4)
grid.newpage()
grid.draw(rbind(
    vis.grob[1:2],vis.grob[4:5],    
    mort.grob[1:2],mort.grob[4:5],mort.grob[3],
    size="first"))
dev.off()
##########################################################################################################################
## ==================================================================================================================== ##
##########################################################################################################################


##,----------------
##| D:C ratio model
##`----------------
data$DC <- data$D.C.ratio..20bin./20
#which samples have an initial DC<5
sa<-subset(data, Time.point=='1' & D.C.ratio..20bin.<5)$sample
data$DC5 <- as.factor(ifelse(data$sample %in% sa,1,0))

a<-glmmPQL(DC ~ Dt*Population, random=~1|sample,data, family=binomial(logit), na.action=na.omit, correlation=corAR1())
newdata<-expand.grid(Population=levels(data$Population), Dt=seq(min(data$Dt), max(data$Dt), l=100))
newdata$Date <- date_decimal(newdata$Dt)
coefs <- fixef(a)
mm <- model.matrix(~ Dt*Population, newdata)
newdata$fit <- binomial(logit)$linkinv(as.vector(coefs %*% t(mm)))*20
predvar <- diag(mm %*% vcov(a) %*% t(mm)) 
newdata$SE <- sqrt(predvar) 
newdata$SE2 <- sqrt(predvar+a$sigma^2)
newdata$lwr <- inverselogit(as.vector(coefs %*% t(mm))-newdata$SE2)*20
newdata$upr <- inverselogit(as.vector(coefs %*% t(mm))+newdata$SE2)*20

data$Dtj <- jitter(data$Dt, factor=2)
dts <- floor_date(unique(data$Date),'month')[-4]+15
p<-ggplot(newdata, aes(y=fit, x=Dt)) +
    geom_point(data=data, aes(y=D.C.ratio..20bin., x=Dt), position=position_jitter(width=0.02,height=0), size=2)+
        geom_line(aes(x=Dt)) +
            facet_grid(~Population)+
                geom_ribbon(aes(ymin=lwr, ymax=upr, x=Dt), fill="grey", alpha=0.5, show.legend=FALSE)+
                    scale_fill_manual(name="April 05 D:C ratio", values=c('red','blue'), breaks=c(1,0),labels=c('<5','>=5'))+
                        scale_shape_manual(name="April 05 D:C ratio", values=c(21,16), breaks=c(1,0), labels=c('<5','>=5'))+
                            scale_color_manual(name="April 05 D:C ratio", values=c('red','blue'), breaks=c(1,0), labels=c('<5','>=5'))+
                                scale_y_continuous("D:C ratio")+
                                    scale_x_continuous("Date", breaks=decimal_date(dts), labels=c(format(dts[1],format="%b\\n%Y"),strtrim(format(dts[c(-1,-5)],format="%b"),1),format(dts[5],format="%b\\n%Y")))+
                                        ggtitle('a)')+
                                            theme(axis.text.x=element_text(size = 9),
                                                  panel.grid.major = element_blank(), # no major grid lines
                                                  panel.grid.minor = element_blank(), # no minor grid lines
                                                  panel.background = element_blank(), # no background
                                                  panel.border = element_blank(), # no plot border
                                                  axis.title.y=element_text(size=11, margin=margin(r=2, unit='lines'),angle=90), # y-axis title
                                                  axis.text.y=element_text(size=9), # y-axis labels
                                                  axis.title.x=element_text(size=11, vjust=-2), # x-axis title
                                                  axis.line = element_line(),
                                                  legend.position=c(1,0),
                                                  legend.justification=c(1,0),
                                                  strip.background=element_blank(),
                                                  strip.text=element_text(size=rel(1.0)),
                                                  plot.margin=unit(c(0.5,0.5,2,2),"lines"),
                                                  plot.title=element_text(hjust=0,vjust=0),
                                                  panel.margin = unit(1, "lines"))#, text=element_text(family='Times'))

library(gtable)
dc.grob <- ggplotGrob(p)
grid.newpage()
grid.draw(rbind(dc.grob[1],dc.grob[4:5], dc.grob[1],dc.grob[3],size="first"))



##,----------------
##| S:H ratio model
##`----------------
data$Dt <- decimal_date(data$Date)
library(reshape)
data.melt <- melt(data[,c('Population','sample','S.H..clade.C.', 'S.H..clade.D.','Dt')],
                  id=c('Population','sample','Dt'),measure.var=c('S.H..clade.C.', 'S.H..clade.D.'))
dts <- floor_date(unique(data$Date),'month')[-4]+15

aa.k <- gamm((value)^0.25~s(Dt, k=6,by=variable)+variable, random=list(sample=~1),
             data=subset(data.melt, Population=='Keppels'),
             family=gaussian)
#gam.check(aa.k$gam)
#plot(aa.k$gam, pages=1)
newdata.k<-expand.grid(Population='Keppels',variable=levels(data.melt$variable), Dt=seq(min(data$Dt), max(data$Dt), l=100))
newdata.k$Date <- date_decimal(newdata.k$Dt)
fit=predict(aa.k$gam, newdata=newdata.k, se=TRUE, type='link')
newdata.k$fit=(fit$fit)^4 #-0.01
newdata.k$lower = (fit$fit-2*fit$se.fit)^4 #-0.01
newdata.k$upper = (fit$fit+2*fit$se.fit)^4 #- 0.01

aa.d <- gamm((value)^0.25~s(Dt, k=3,by=variable)+variable, random=list(sample=~1),
             data=subset(data.melt, Population=='Davies'),
             family=gaussian)
#gam.check(aa.d$gam)
#plot(aa.d$gam, pages=1)
newdata.d<-expand.grid(Population='Davies',variable=levels(data.melt$variable), Dt=seq(min(data$Dt), max(data$Dt), l=100))
newdata.d$Date <- date_decimal(newdata.d$Dt)
fit=predict(aa.d$gam, newdata=newdata.d, se=TRUE, type='link')
newdata.d$fit=(fit$fit)^4 #-0.01
newdata.d$lower = (fit$fit-2*fit$se.fit)^4#-0.01
newdata.d$upper = (fit$fit+2*fit$se.fit)^4#-0.01

aa.m <- gamm((value)^0.25~s(Dt, k=6,by=variable)+variable, random=list(sample=~1),
             data=subset(data.melt, Population=='Magnetic Island'),
             family=gaussian)
#gam.check(aa.m$gam)
#plot(aa.m$gam, pages=1)
newdata.m<-expand.grid(Population='Magnetic Island',variable=levels(data.melt$variable), Dt=seq(min(data$Dt), max(data$Dt), l=100))
newdata.m$Date <- date_decimal(newdata.m$Dt)
fit=predict(aa.m$gam, newdata=newdata.m, se=TRUE, type='link')
newdata.m$fit=(fit$fit)^4#-0.01
newdata.m$lower = (fit$fit-2*fit$se.fit)^4#-0.01
newdata.m$upper = (fit$fit+2*fit$se.fit)^4#-0.01

newdata <- rbind(newdata.k, newdata.d, newdata.m)
p <- ggplot(data.melt, aes(y=value, x=Dt, group=variable,color=variable,fill=variable, linetype=variable)) +
  #geom_point(data=data, aes(y=S.H..C.D., x=Dt), position=position_jitter(0.02)) +
  geom_point(position=position_jitter(0.02), show.legend=TRUE) +
                                        #geom_smooth(alpha=0.5,
                                        #color='black')+
      geom_ribbon(data=newdata, aes(y=fit, ymin=lower, ymax=upper), color=NA)+
      geom_line(data=newdata, aes(y=fit), color='black')+
      #stat_smooth(method='gam', formula='y~s(x)', method.args=list(family='binomial'))+
  facet_grid(~Population)+
  #scale_fill_manual('Symbiont Clade',values=c('#2b8cbe','#f03b20'),labels=c('C','D'))+
  #scale_color_manual('Symbiont Clade',values=c('#2b8cbe','#f03b20'),labels=c('C','D'))+
  scale_fill_manual('Symbiont Clade',values=c('grey80','grey30'),labels=c('C','D'))+
  scale_color_manual('Symbiont Clade',values=c('grey80','grey30'),labels=c('C','D'))+
  scale_linetype_manual('Symbiont Clade',values=c(1,2),labels=c('C','D'))+
  #scale_linetype_discrete('Symbiont Clade')+
  scale_y_continuous("S:H ratio",trans=sqrt_trans(), breaks=c(0,0.1,0.5,1,2,4,6,8))+
  scale_x_continuous("Date", breaks=decimal_date(dts), labels=c(format(dts[1],format="%b\\n%Y"),strtrim(format(dts[c(-1,-5)],format="%b"),1),format(dts[5],format="%b\\n%Y")))+
  ggtitle('b)')+
  theme(axis.text.x=element_text(size = 9),
        panel.grid.major = element_blank(), # no major grid lines
        panel.grid.minor = element_blank(), # no minor grid lines
        panel.background = element_blank(), # no background
        panel.border = element_blank(), # no plot border
        axis.title.y=element_text(size=11,margin=margin(r=2, unit='lines'),angle=90), # y-axis title
        axis.text.y=element_text(size=9), # y-axis labels
        axis.title.x=element_text(size=11, vjust=-2), # x-axis title
        axis.line = element_line(),
        legend.position=c(1,1.1),
        legend.justification=c(1,1),
        strip.background=element_blank(),
        strip.text=element_text(size=rel(1.0)),
        plot.margin=unit(c(0.5,0.5,2,2),"lines"),
        plot.title=element_text(hjust=0,vjust=0),
        panel.margin = unit(1, "lines"))#, text=element_text(family='Times'))



library(gtable)
sh.grob <- ggplotGrob(p)
grid.newpage()
grid.draw(rbind(sh.grob[1],sh.grob[4:5], sh.grob[1],sh.grob[3],size="first"))

##,----------
##| Dominance
##`----------
dat <- droplevels(subset(data, ITS2.DGGE.dominant.genotype!=""))
p<-ggplot(dat, aes(x=Dt, fill=ITS2.DGGE.dominant.genotype)) +
  geom_bar(aes(y = ..count..*100/sum(..count..)), position='fill', color='transparent',binwidth=0.1)+
      geom_bar(aes(y = ..count..*100/sum(..count..)), position='fill', color='black',binwidth=0.1, show_guide=FALSE)+
      facet_grid(~Population)+
  facet_grid(~Population)+
  #scale_fill_manual('Dominant genotype', values=c('#bdc9e1','#2b8cbe','#045a8d','#f03b20'))+
  scale_fill_manual('Dominant genotype', values=c('grey99','grey80','grey60','grey30'))+
  scale_y_continuous('Percentage of colonies', breaks=c(0,0.25,0.5,0.75,1), labels=c(0,25,50,75,100))+
  scale_x_continuous("Date", breaks=decimal_date(dts), labels=c(format(dts[1],format="%b\\n%Y"),strtrim(format(dts[c(-1,-5)],format="%b"),1),format(dts[5],format="%b\\n%Y")))+
  #scale_x_continuous("Date", breaks=decimal_date(dts), labels=c(format(dts[1],format="%b\\n%Y"),format(dts[c(-1,-5)],format="%b"),format(dts[5],format="%b\\n%Y")))+
  ggtitle('a)')+
theme(axis.text.x=element_text(size = 9),
        panel.grid.major = element_blank(), # no major grid lines
        panel.grid.minor = element_blank(), # no minor grid lines
        panel.background = element_blank(), # no background
        panel.border = element_blank(), # no plot border
        axis.title.y=element_text(size=11,margin=margin(r=2, unit='lines'),angle=90), # y-axis title
        axis.text.y=element_text(size=9), # y-axis labels
        axis.title.x=element_text(size=11, vjust=-2), # x-axis title
        axis.line = element_line(),
      legend.position=c(1,1.40),
      legend.key=element_rect(color='black'),
         #legend.position=c(1,-0.1),
        legend.justification=c(1,1),
      legend.direction='horizontal',
        strip.background=element_blank(),
        strip.text=element_text(size=rel(1.0)),
        plot.margin=unit(c(0.5,0.5,2,2),"lines"),
        plot.title=element_text(hjust=0,vjust=0),
        panel.margin = unit(1, "lines"))#, text=element_text(family='Times'))

library(gtable)
dom <- ggplotGrob(p)
gg <- dom

##,-------------------------------------------------------------------------------------------------------------------------------
##| Figure 3======================================================================================================================
##`-------------------------------------------------------------------------------------------------------------------------------
grid.newpage()
grid.draw(rbind(
dom[1:2],dom[4:5],
    sh.grob[1:2],sh.grob[4:5],sh.grob[3],
    size="first"))
## External graphics devices
pdf(file='figures/Figure3.pdf', width=175/25.4, height=100/25.4)
grid.newpage()
grid.draw(rbind(
    dom[1:2],dom[4:5],#dom[7],
    sh.grob[1:2],sh.grob[4:5],sh.grob[3],
    size="first"))
dev.off()

system('pdf2ps -dLanguageLevel=3 figures/Figure3.pdf figures/Figure3.eps')
system('pdftops -level3 -eps figures/Figure3.pdf figures/Figure3a.eps')
#system('convert figures/Figure1.pdf figures/Figure1.eps')

postscript(file='figures/Figure3p.eps', width=175/25.4, height=100/25.4, horizontal=FALSE, paper='special', fonts='Times')
grid.newpage()
grid.draw(rbind(
    dom[1:2],dom[4:5],#dom[7],
    sh.grob[1:2],sh.grob[4:5],sh.grob[3],
    size="first"))
dev.off()

library(Cairo)
CairoFonts('Nimbus Roman No9 L')
CairoPS(file='figures/Figure3.eps', width=175/25.4, height=100/25.4, horizontal=FALSE, paper='special',pointsize=8)
grid.newpage()
grid.draw(rbind(
    dom[1:2],dom[4:5],#dom[7],
    sh.grob[1:2],sh.grob[4:5],sh.grob[3],
    size="first"))
dev.off()

cairo_ps(file='figures/Figure3.eps', width=175/25.4, height=100/25.4)
grid.newpage()
grid.draw(rbind(
    dom[1:2],dom[4:5],#dom[7],
    sh.grob[1:2],sh.grob[4:5],sh.grob[3],
    size="first"))
dev.off()


## Table 2 compilations
t1 <- subset(data, Time.point==1,select=c('Population','sample','S.H..clade.C.','S.H..clade.D.','D.C.ratio..20bin.','S.H..C.D.','Mortality'))
a <- subset(data, Time.point==2,select=c('Population','sample','S.H..clade.C.','S.H..clade.D.','D.C.ratio..20bin.','S.H..C.D.'))
b <- subset(data, Time.point==6, select=c('Population','sample','Mortality'))
dat<-join(a,b)

library(plyr)
aa<-ddply(subset(dat,Mortality==100), ~Population, function(x) {
  data.frame(numcolwise(mean, na.rm=TRUE)(x), Lt5=length(subset(x,D.C.ratio..20bin.<5)$D.C.ratio..20bin.),Lt5P=100*length(subset(x,D.C.ratio..20bin.<5)$D.C.ratio..20bin.)/nrow(na.omit(x)),n=nrow(na.omit(x)))
})
bb<-ddply(subset(dat,Mortality<100), ~Population, function(x) {
  data.frame(numcolwise(mean, na.rm=TRUE)(x), Lt5=length(subset(x,D.C.ratio..20bin.<5)$D.C.ratio..20bin.),Lt5P=100*length(subset(x,D.C.ratio..20bin.<5)$D.C.ratio..20bin.)/nrow(na.omit(x)),n=nrow(na.omit(x)))
})

library(pander)
pandoc.table(rbind(aa,bb), split.tables=Inf, caption='Mean September 2005 S:H and D:C ratios as well as the number and proportion of colonies with D:C ratios less than 5 for colonies that displayed either 100% mortality or less than 100% mortality by April 2006.') 

tab <- rbind(aa,bb)
colnames(tab) <- c('Population','S:H ratio (Clade C)', 'S:H ratio (Clade D)', 'D:C ratio', 'S:H ratio (C D)','Mortality', 'lt5','Lt5P','n')
tab <- tab[,c(1,6,9,7,8,2:3,5)]
library(xtable)

library(Hmisc)
rownames(tab) <- c(paste(tab$Population[1:2]),paste(tab$Population[3:4],' ',sep=''))
tab<- tab[,-1:-2]


latex(object=tab,cdec=c(0,0,2,3,3,3),file="",align='cccccc', title='',
      cgroup=c('Apr 06','Sept 05 D:C \\\\[<\\\\]5','Mean Sept S:H ratio'), n.cgroup=c(1,2,3),
      colnamesTexCmd='bfseries',
      rowlabel='',
      dcolumn=TRUE,
      booktabs=TRUE,
      colheads=c('n','n','\\\\%','(Clade C)', '(Clade D)','(C and D)'),
      rgroup=c('100\\\\% Mortality','\\\\[<\\\\]100\\\\% Mortality'),
      caption='Mean September 2005 S:H and D:C ratios as well as the number and proportion of colonies with D:C ratios less than 5 for colonies that displayed either 100\\\\% mortality or less than 100\\\\% mortality by April 2006.',
      caption.loc='top'
      )
#April 05
t1 <- subset(data, Time.point==1,select=c('Population','sample','S.H..clade.C.','S.H..clade.D.','D.C.ratio..20bin.','S.H..C.D.'))
t6 <- subset(data, Time.point==6, select=c('Population','sample','Mortality'))
dat<-join(t1,t6)

library(plyr)
mort.dead<-ddply(subset(dat,Mortality==100), ~Population, function(x) {
  data.frame(numcolwise(mean, na.rm=TRUE)(x), Lt5=length(subset(x,D.C.ratio..20bin.<5)$D.C.ratio..20bin.),Lt5P=100*length(subset(x,D.C.ratio..20bin.<5)$D.C.ratio..20bin.)/nrow(na.omit(x)),n=nrow(na.omit(x)))
})
mort.alive<-ddply(subset(dat,Mortality<100), ~Population, function(x) {
  data.frame(numcolwise(mean, na.rm=TRUE)(x), Lt5=length(subset(x,D.C.ratio..20bin.<5)$D.C.ratio..20bin.),Lt5P=100*length(subset(x,D.C.ratio..20bin.<5)$D.C.ratio..20bin.)/nrow(na.omit(x)),n=nrow((x)))
})

pandoc.table(rbind(mort.dead,mort.alive), split.tables=Inf, caption='Mean September 2005 S:H and D:C ratios as well as the number and proportion of colonies with D:C ratios less than 5 for colonies that displayed either 100% mortality or less than 100% mortality by April 2006.') 

tab <- rbind(mort.dead,mort.alive)
colnames(tab) <- c('Population','S:H ratio (Clade C)', 'S:H ratio (Clade D)', 'D:C ratio', 'S:H ratio (C D)','Mortality', 'lt5','Lt5P','n')
tab <- tab[,c(1,6,9,7,8,2:3,5)]
library(Hmisc)
rownames(tab) <- c(paste(tab$Population[1:2]),paste(tab$Population[3:4],' ',sep=''))
tab<- tab[,-1:-2]


latex(object=tab,cdec=c(0,0,2,3,3,3),file="",align='cccccc', title='',
      cgroup=c('Apr 06','Apr 05 D:C \\\\[<\\\\]5','Mean Apr 05 S:H ratio'), n.cgroup=c(1,2,3),
      colnamesTexCmd='bfseries',
      rowlabel='',
      dcolumn=TRUE,
      booktabs=TRUE,
      colheads=c('n','n','\\\\%','(Clade C)', '(Clade D)','(C and D)'),
      rgroup=c('100\\\\% Mortality','\\\\[<\\\\]100\\\\% Mortality'),
      caption='Mean April 2005 S:H and D:C ratios as well as the number and proportion of colonies with D:C ratios less than 5 for colonies that displayed either 100\\\\% mortality or less than 100\\\\% mortality by April 2006.',
      caption.loc='top'
      )


##,----------------------
##| Suppluimentary figure
##`----------------------
data1$SH_April <- data1[,'S.H..C.D._2006-04-30']
data1$SH_April <- data1[,'S.H..C.D._2005-04-05'] #recent change

pN1 <- ggplot(data1, aes(y=SH_Dec, x=SH_April, color=Population, fill=Population)) +  
  geom_point(position=position_jitter(width=0.1,height=0.001), show_guide=FALSE) +
      geom_smooth(method='lm', show_guide=FALSE)+
  #facet_grid(~Population)+
  scale_y_sqrt("S:H ratio (Dec)")+
  scale_x_sqrt("S:H ratio (Apr)")+
  scale_fill_manual('Population', values=c('grey75','grey45','grey15'))+
  scale_color_manual('Population', values=c('grey75','grey45','grey15'))+
  theme_classic(9)+
  theme(panel.grid.major = element_blank(), # no major grid lines
  panel.grid.minor = element_blank(), # no minor grid lines
  panel.background = element_blank(), # no background
  panel.border = element_blank(), # no plot border
  axis.title.y=element_text(size=11, vjust=2,angle=90), # y-axis title
  axis.text.y=element_text(size=9), # y-axis labels
  axis.title.x=element_text(size=11, vjust=-0), # x-axis title
  axis.text.x=element_text(size=9), # x-axis labels
  axis.line = element_line(),
  legend.position=c(1,0),
  legend.justification=c(1,0),
  plot.margin=unit(c(0.5,0.5,0.5,1),"lines"))
pN1
summary(s2<-lm(sqrt(SH_Dec) ~ sqrt(SH_April)*Population, data=data1))

data1$SH_Sep <- data1[,'S.H..C.D._2005-09-14']
pN2 <- ggplot(data1, aes(y=SH_Dec, x=SH_Sep, color=Population,
fill=Population)) +
geom_point(position=position_jitter(width=0.1,height=0.001)) +
geom_smooth(method='lm')+
  #facet_grid(~Population)+
  scale_y_sqrt("S:H ratio (Dec)")+
  scale_x_sqrt("S:H ratio (Sept)")+
  scale_fill_manual('Population', values=c('grey75','grey45','grey15'))+
  scale_color_manual('Population', values=c('grey75','grey45','grey15'))+
  theme_classic(9)+
  theme(panel.grid.major = element_blank(), # no major grid lines
  panel.grid.minor = element_blank(), # no minor grid lines
  panel.background = element_blank(), # no background
  panel.border = element_blank(), # no plot border
  axis.title.y=element_text(size=11, vjust=2,angle=90), # y-axis title
  axis.text.y=element_text(size=9), # y-axis labels
  axis.title.x=element_text(size=11, vjust=-0), # x-axis title
  axis.text.x=element_text(size=9), # x-axis labels
  axis.line = element_line(),
  legend.position=c(0.5,1.1),
  legend.justification=c(0.5,1),
  plot.margin=unit(c(0.5,0.5,0.5,1),"lines"),legend.key.height=unit(0.75,'lines'))

summary(lm(sqrt(SH_Dec) ~ sqrt(SH_Sep)*Population, data=data1))


grid.arrange(pN1, pN2, ncol=2)

## External graphics devices
pdf(file='figures/Figure6New.pdf', width=175/25.4, height=100/25.4)
grid.arrange(pN1, pN2, ncol=2)
dev.off()

system('pdf2ps -dLanguageLevel=3 figures/Figure6New.pdf figures/Figure6New.eps')
system('pdftops -level3 -eps figures/Figure6New.pdf figures/Figure6Newa.eps')
#system('convert figures/Figure1.pdf figures/Figure1.eps')

postscript(file='figures/Figure6Newp.eps', width=175/25.4, height=100/25.4, horizontal=FALSE, paper='special', fonts='Times')
grid.arrange(pN1, pN2, ncol=2)
dev.off()

library(Cairo)
CairoFonts('Nimbus Roman No9 L')
CairoPS(file='figures/Figure6New.eps', width=175/25.4, height=100/25.4, horizontal=FALSE, paper='special',pointsize=8)
grid.arrange(pN1, pN2, ncol=2)
dev.off()

cairo_ps(file='figures/Figure6New.eps', width=175/25.4, height=100/25.4)
grid.arrange(pN1, pN2, ncol=2)
dev.off()



##,-----------------------------------------------------------------------------------------------------
##| Brocken stick models================================================================================
##`-----------------------------------------------------------------------------------------------------

## September---------------------------------------
## Keppels               
before <- function(x,bp) ifelse(x<bp, bp-x,0)
after <- function(x,bp) ifelse(x<bp, 0, x-bp)

data1.k <- subset(data1, Population=="Keppels")
foo <- function(bp) {
  mod <- glm(SH_Dec ~ before(DC_Sep,bp)+after(DC_Sep,bp),data1.k, family=quasipoisson, na.action=na.omit)
  mod$deviance
}
search.range <- c(min(data1.k$DC_Sep,na.rm=TRUE)+0.5,max(data1.k$DC_Sep, na.rm=TRUE)-0.5)
search.range <- c(2,10)
foo.opt <- optimize(foo, interval = search.range)
bp <- foo.opt$minimum
bp
mod.k <- glm(SH_Dec ~ before(DC_Sep,bp) + after(DC_Sep,bp),data1.k, family=quasipoisson, na.action=na.omit)
 
foo1 <- function(bp, data) {
  mod <- glm(SH_Dec ~ before(DC_Sep,bp)+after(DC_Sep,bp),data, family=quasipoisson, na.action=na.omit)
  mod$deviance
}
 
set.seed(4)
n <- 8
bps <- NULL
for (i in 1:nrow(data1.k)) {
  dt <-data1.k[sample(1:nrow(data1.k),size=nrow(data1.k)-n, replace=FALSE),]
  #  dt <- data1.k[-i,]
  search.range <- c(min(dt$DC_Sep,na.rm=TRUE)+0.5,max(dt$DC_Sep, na.rm=TRUE)-0.5)
  search.range <- c(2,10)
  foo.opt <- optimize(foo1, interval = search.range, data=dt)
  bp.m <- foo.opt$minimum
  bps <- c(bps,bp.m)
}
library(gmodels)
keppel.ci1 <- ci(bps)#[2:3]
keppel.m1<-quantile(bps)#[3]
keppel.thres1 <- c(keppel.m1,keppel.ci1)
  
names(mod.k$coefficients) <- c('(Intercept)','Before','After')
summary(mod.k)
library(multcomp)

summary(glht(mod.k, linfct=rbind(c(0,1,-1))))
newdata.k<-with(data1.k, expand.grid(DC_Sep=seq(min(DC_Sep,na.rm=TRUE),max(DC_Sep,na.rm=TRUE), l=1000)))

pred <- predict(mod.k, newdata=newdata.k, se=TRUE)
newdata.k$fit <- exp(pred$fit)
newdata.k$lwr <- exp(pred$fit-2*pred$se.fit)
newdata.k$upr <- exp(pred$fit+2*pred$se.fit)
newdata.k$BA <- factor(ifelse(newdata.k$DC_Sep>bp,1,0))
newdata.k$Population="Keppels"
bp

#Davies  
data1.d <- subset(data1, Population=="Davies")
foo <- function(bp) {
  mod <- glm(SH_Dec ~ before(DC_Sep,bp)+after(DC_Sep,bp),data1.d, family=quasipoisson, na.action=na.omit)
  mod$deviance
} 
search.range <- c(min(data1.d$DC_Sep,na.rm=TRUE),max(data1.d$DC_Sep, na.rm=TRUE))
search.range <- c(2,10)
foo.opt <- optimize(foo, interval = search.range)
bp <- foo.opt$minimum
bp
mod.d <- glm(SH_Dec ~ before(DC_Sep,bp) + after(DC_Sep,bp),data1.d, family=quasipoisson, na.action=na.omit)
names(mod.d$coefficients) <- c('(Intercept)','Before','After')
summary(mod.d)
summary(glht(mod.d, linfct=rbind(c(0,1,-1))))
mod.d <- glm(SH_Dec ~ DC_Sep,data1.d, family=quasipoisson, na.action=na.omit)
names(mod.d$coefficients) <- c('(Intercept)','Before')
newdata.d<-with(data1.d, expand.grid(DC_Sep=seq(min(DC_Sep,na.rm=TRUE),max(DC_Sep,na.rm=TRUE), l=1000)))
pred <- predict(mod.d, newdata=newdata.d, se=TRUE)
newdata.d$fit <- exp(pred$fit)
newdata.d$lwr <- exp(pred$fit-2*pred$se.fit)
newdata.d$upr <- exp(pred$fit+2*pred$se.fit)
newdata.d$BA <- factor(0)
newdata.d$Population="Davies"

#Magnetic Island  
data1.m <- subset(data1, Population=="Magnetic Island")
foo <- function(bp) {
  mod <- glm(SH_Dec ~ before(DC_Sep,bp)+after(DC_Sep,bp),data1.m, family=quasipoisson, na.action=na.omit)
  mod$deviance
}
search.range <- c(min(data1.m$DC_Sep,na.rm=TRUE),max(data1.m$DC_Sep, na.rm=TRUE)) 
#search.range <- c(2,10)
foo.opt <- optimize(foo, interval = search.range)
mod.m <- glm(SH_Dec ~ DC_Sep,data1.m, family=quasipoisson, na.action=na.omit)
names(mod.m$coefficients) <- c('(Intercept)','After')
summary(mod.m)
newdata.m<-with(data1.m, expand.grid(DC_Sep=seq(min(DC_Sep,na.rm=TRUE),max(DC_Sep,na.rm=TRUE), l=1000)))
pred <- predict(mod.m, newdata=newdata.m, se=TRUE)
newdata.m$fit <- exp(pred$fit)
newdata.m$lwr <- exp(pred$fit-2*pred$se.fit)
newdata.m$upr <- exp(pred$fit+2*pred$se.fit)
newdata.m$BA <- factor(1)
newdata.m$Population="Magnetic Island" 

#Full figure
newdata <- rbind(newdata.k, newdata.d, newdata.m) 
newdata$Population <- factor(newdata$Population, levels=c('Keppels', 'Davies', 'Magnetic Island'))
    
changepointFigSept<-ggplot(newdata, aes(y=fit, x=DC_Sep)) + 
  geom_point(data=data1, aes(y=SH_Dec, x=DC_Sep))+
  geom_line(aes(colour=BA)) + facet_grid(~Population)+
  geom_ribbon(aes(ymin=lwr, ymax=upr, x=DC_Sep, fill=BA), alpha=0.5)+
  scale_fill_manual(guide=FALSE, values=c('grey40','grey80'))+  scale_colour_manual(guide=FALSE, values=c('black','black'))+
#scale_fill_discrete(guide=FALSE)+  scale_colour_discrete(guide=FALSE)+
  scale_y_continuous("S:H ratio (Dec)")+coord_cartesian(ylim=c(0,0.4))+
  scale_x_continuous("D:C ratio (Sept)")+
theme(axis.text.x=element_text(size = rel(0.75)),
        panel.grid.major = element_blank(), # no major grid lines
        panel.grid.minor = element_blank(), # no minor grid lines
        panel.background = element_blank(), # no background
        panel.border = element_blank(), # no plot border
        axis.title.y=element_text(size=15, margin=margin(r=2, unit='lines'),angle=90), # y-axis title
        axis.text.y=element_text(size=10), # y-axis labels
        axis.title.x=element_text(size=15, vjust=-2), # x-axis title
        axis.line = element_line(),
        legend.position=c(1,0),
        legend.justification=c(1,0),
        strip.background=element_blank(),
        strip.text=element_text(size=rel(1.0)),
        plot.margin=unit(c(0.5,0.5,2,2),"lines"),
        panel.margin = unit(1, "lines"))
changepointFigSept



## April------------------------------------------------------------------------------------------------------------
#Keppels 
foo <- function(bp) {
  mod <- glm(SH_Dec ~ before(DC_Apr,bp)+after(DC_Apr,bp),data1.k, family=quasipoisson, na.action=na.omit)
  mod$deviance
}
search.range <- c(min(data1.k$DC_Apr,na.rm=TRUE)+0.5,max(data1.k$DC_Apr, na.rm=TRUE)-0.5)
search.range <- c(2,10)
foo.opt <- optimize(foo, interval = search.range)
bp <- foo.opt$minimum
bp
mod.k1 <- glm(SH_Dec ~ before(DC_Apr,bp) + after(DC_Apr,bp),data1.k, family=quasipoisson, na.action=na.omit)
names(mod.k1$coefficients) <- c('(Intercept)','Before','After')
summary(mod.k1)
library(multcomp)
summary(glht(mod.k1, linfct=rbind(c(0,1,-1))))
newdata.k<-with(data1.k, expand.grid(DC_Apr=seq(min(DC_Apr,na.rm=TRUE),max(DC_Apr,na.rm=TRUE), l=1000)))

pred <- predict(mod.k1, newdata=newdata.k, se=TRUE)
newdata.k$fit <- exp(pred$fit)
newdata.k$lwr <- exp(pred$fit-2*pred$se.fit)
newdata.k$upr <- exp(pred$fit+2*pred$se.fit)
newdata.k$BA <- factor(ifelse(newdata.k$DC_Apr>bp,1,0))
newdata.k$Population="Keppels"

foo1 <- function(bp, data) {
  mod <- glm(SH_Dec ~ before(DC_Apr,bp)+after(DC_Apr,bp),data, family=quasipoisson, na.action=na.omit)
  mod$deviance
}

set.seed(2)
bps <- NULL
for (i in 1:nrow(data1.k)) {
  dt <-data1.k[sample(1:nrow(data1.k),size=nrow(data1.k)-n, replace=FALSE),]
  #  dt <- data1.k[-i,]
  search.range <- c(min(dt$DC_Apr,na.rm=TRUE)+0.5,max(dt$DC_Apr, na.rm=TRUE)-0.5)
  search.range <- c(2,10)
  foo.opt <- optimize(foo1, interval = search.range, data=dt)
  bp.m <- foo.opt$minimum
  bps <- c(bps,bp.m)
}
#library(gmodels)
#ci(bps)

library(gmodels)
keppel.ci2 <- ci(bps)#[2:3]
keppel.m2<-quantile(bps)#[3]
keppel.thres2 <- c(keppel.m2,keppel.ci2)


#Davies 
data1.d <- subset(data1, Population=="Davies")
foo <- function(bp) {
  mod <- glm(SH_Dec ~ before(DC_Apr,bp)+after(DC_Apr,bp),data1.d, family=quasipoisson, na.action=na.omit)
  mod$deviance
}
search.range <- c(min(data1.d$DC_Apr,na.rm=TRUE),max(data1.d$DC_Apr, na.rm=TRUE))
#search.range <- c(2,10)
foo.opt <- optimize(foo, interval = search.range)
bp <- foo.opt$minimum
bp
mod.d1 <- glm(SH_Dec ~ before(DC_Apr,bp) + after(DC_Apr,bp),data1.d, family=quasipoisson, na.action=na.omit)
names(mod.d1$coefficients) <- c('(Intercept)','Before','After')
summary(mod.d1)
summary(glht(mod.d1, linfct=rbind(c(0,1,-1))))
mod.d1 <- glm(SH_Dec ~ DC_Apr,data1.d, family=quasipoisson, na.action=na.omit)
names(mod.d1$coefficients) <- c('(Intercept)','Before')
newdata.d<-with(data1.d, expand.grid(DC_Apr=seq(min(DC_Apr,na.rm=TRUE),max(DC_Apr,na.rm=TRUE), l=1000)))
pred <- predict(mod.d1, newdata=newdata.d, se=TRUE)
newdata.d$fit <- exp(pred$fit)
newdata.d$lwr <- exp(pred$fit-2*pred$se.fit)
newdata.d$upr <- exp(pred$fit+2*pred$se.fit)
newdata.d$BA <- factor(0)
newdata.d$Population="Davies"

#Magnetic Island 
data1.m <- subset(data1, Population=="Magnetic Island")
foo <- function(bp) {
  mod <- glm(SH_Dec ~ before(DC_Apr,bp)+after(DC_Apr,bp),data1.m, family=quasipoisson, na.action=na.omit)
  mod$deviance
}
search.range <- c(min(data1.m$DC_Apr,na.rm=TRUE),max(data1.m$DC_Apr, na.rm=TRUE)) 
#search.range <- c(2,10)
foo.opt <- optimize(foo, interval = search.range)
bp <- foo.opt$minimum
bp
mod.m1 <- glm(SH_Dec ~ DC_Apr,data1.m, family=quasipoisson, na.action=na.omit)
names(mod.m1$coefficients) <- c('(Intercept)','After')
summary(mod.m1)
newdata.m<-with(data1.m, expand.grid(DC_Apr=seq(min(DC_Apr,na.rm=TRUE),max(DC_Apr,na.rm=TRUE), l=1000)))
pred <- predict(mod.m1, newdata=newdata.m, se=TRUE)
newdata.m$fit <- exp(pred$fit)
newdata.m$lwr <- exp(pred$fit-2*pred$se.fit)
newdata.m$upr <- exp(pred$fit+2*pred$se.fit)
newdata.m$BA <- factor(1)
newdata.m$Population="Magnetic Island"


newdata <- rbind(newdata.k, newdata.d, newdata.m) 
newdata$Population <- factor(newdata$Population, levels=c('Keppels', 'Davies', 'Magnetic Island'))
  
changepointFig<-ggplot(newdata, aes(y=fit, x=DC_Apr)) + 
  geom_point(data=data1, aes(y=SH_Dec, x=DC_Apr))+
  geom_ribbon(aes(ymin=lwr, ymax=upr, x=DC_Apr, fill=BA), alpha=0.5)+
    geom_line(aes(colour=BA)) + facet_grid(~Population)+
    scale_fill_manual(guide=FALSE, values=c('grey40','grey80'))+  scale_colour_manual(guide=FALSE, values=c('black','black'))+
#  scale_fill_discrete(guide=FALSE)+  scale_colour_discrete(guide=FALSE)+
  scale_y_continuous("S:H ratio (Dec)")+coord_cartesian(ylim=c(0,0.4))+
  scale_x_continuous("D:C ratio (Apr)")+
  ggtitle('b)')+
theme(axis.text.x=element_text(size = 9),
        panel.grid.major = element_blank(), # no major grid lines
        panel.grid.minor = element_blank(), # no minor grid lines
        panel.background = element_blank(), # no background
        panel.border = element_blank(), # no plot border
        axis.title.y=element_text(size=11, margin=margin(r=2, unit='lines'),angle=90), # y-axis title
        axis.text.y=element_text(size=9), # y-axis labels
        axis.title.x=element_text(size=11, margin=margin(t=2, unit='lines')), # x-axis title
        axis.line = element_line(),
        legend.position=c(1,0),
        legend.justification=c(1,0),
        strip.background=element_blank(),
        strip.text=element_text(size=rel(1.0)),
        plot.margin=unit(c(0.5,0.5,2,2),"lines"),
        plot.title=element_text(hjust=0,vjust=0),
        panel.margin = unit(1, "lines"))#, text=element_text(family='Times'))
changepointFig

## the full figure
grid.newpage()
changepoint.grob <- ggplotGrob(changepointFig)

grid.draw(rbind(
dc.grob[3],dc.grob[1:2],dc.grob[4:5],    
changepoint.grob[1:2],changepoint.grob[4:7],
                size="first"))

#External devices
pdf(file='figures/Figure4.pdf', width=175/25.4, height=100/25.4)
grid.newpage()
grid.draw(rbind(
    dc.grob[3],dc.grob[1:2],dc.grob[4:5],    
    changepoint.grob[1:2],changepoint.grob[4:7],#changepoint.grob[3],
    size="first"))
dev.off()


#system('pdf2ps -dLanguageLevel=3 figures/Figure2.pdf figures/Figure2.eps')
system('pdftops -level3 -eps figures/Figure4.pdf figures/Figure4a.eps')
#system('convert figures/Figure1.pdf figures/Figure4.eps')

postscript(file='figures/Figure4p.eps', width=175/25.4, height=100/25.4, horizontal=FALSE, paper='special', fonts='Times')
grid.newpage()
grid.draw(rbind(
    dc.grob[3],dc.grob[1:2],dc.grob[4:5],    
    changepoint.grob[1:2],changepoint.grob[4:7],#changepoint.grob[3],
    size="first"))
dev.off()

library(Cairo)
CairoFonts('Nimbus Roman No9 L')
CairoPS(file='figures/Figure4.eps', width=175/25.4, height=100/25.4, horizontal=FALSE, paper='special',pointsize=8)
grid.newpage()
grid.draw(rbind(
    dc.grob[3],dc.grob[1:2],dc.grob[4:5],    
    changepoint.grob[1:2],changepoint.grob[4:7],#changepoint.grob[3],
    size="first"))
dev.off()

cairo_ps(file='figures/Figure4.eps', width=175/25.4, height=100/25.4)
grid.newpage()
grid.draw(rbind(
    dc.grob[3],dc.grob[1:2],dc.grob[4:5],    
    changepoint.grob[1:2],changepoint.grob[4:7],#changepoint.grob[3],
    size="first"))
dev.off()


## Change-point thresholds
## library(pander)
## thres <- c('D:C ratio (Sept)'= sprintf("%0.3f [%0.2f,%0.3f]",keppel.thres1[3],keppel.thres1[1],keppel.thres1[4]),
##            'D:C ratio (Apr)'= sprintf("%0.3f [%0.2f,%0.3f]",keppel.thres2[3],keppel.thres2[1],keppel.thres2[4]))

## thres.D <- c('D:C ratio (Sept)'= sprintf("%0.3f [%0.2f,%0.3f]",keppel.thres1.D[3],keppel.thres1.D[1],keppel.thres1.D[4]),
##            'D:C ratio (Apr)'= sprintf("%0.3f [%0.2f,%0.3f]",keppel.thres2.D[3],keppel.thres2.D[1],keppel.thres2.D[4]))
## thres <- rbind(thres,thres.D)
## rownames(thres) <- c('C \\\\& D symbionts', 'D symbionts')
## pandoc.table(thres)

library(texreg)
texreg(list(mod.k1, mod.k)
       ,custom.model.names=c('D/C ratio (Apr)','D/C ratio (Sep)')
       ,custom.coef.names=c('Intercept','Before','After')
       ,ci.force=TRUE, single.row=TRUE
       , caption="Piecewise regression coefficients (log odds scale) for associations between S:H ratio in December 2005 and D:C ratio in April 2005 and September 2005 for Keppels sourced corals."
       , caption.above=TRUE
       , include.bic=FALSE,include.aic=FALSE,include.loglik = FALSE,use.packages=FALSE,
       ,booktabs = TRUE, dcolumn = TRUE)
library(texreg)
texreg(list(mod.d1, mod.d)
       ,custom.model.names=c('D/C ratio (Apr)','D/C ratio (Sep)')
       ,custom.coef.names=c('Intercept','Before')
       ,ci.force=TRUE, single.row=TRUE
       , caption="Piecewise regression coefficients (log odds scale) for associations between S:H ratio in December 2005 and D:C ratio in April 2005 and September 2005 for Davies sourced corals."
       , caption.above=TRUE
       , include.bic=FALSE,include.aic=FALSE,include.loglik = FALSE,use.packages=FALSE,
       ,booktabs = TRUE, dcolumn = TRUE)
library(texreg)
texreg(list(mod.m1, mod.m)
       ,custom.model.names=c('D/C ratio (Apr)','D/C ratio (Sep)')
       ,custom.coef.names=c('Intercept','Before')
       ,ci.force=TRUE, single.row=TRUE
       , caption="Piecewise regression coefficients (log odds scale) for associations between S:H ratio in December 2005 and D:C ratio in April 2005 and September 2005 for Magnetic Island sourced corals."
       , caption.above=TRUE
       , include.bic=FALSE,include.aic=FALSE,include.loglik = FALSE,use.packages=FALSE,
       ,booktabs = TRUE, dcolumn = TRUE)




#####################
## Heat map figure ##
#####################

## sort by Mortality in Time.point 5
## cross for missing
data$Dt <- decimal_date(data$Date)
library(dplyr)
dd=data %>% group_by(sample) %>% droplevels() %>% do({
    x<-.
    t1_4=unique(x$ITS2.DGGE.dominant.genotype[x$Time.point %in% c(1,2,3,4)])
    t1_4=gsub('(.).*','\\\\1', t1_4)
    t1_4=unique(t1_4[t1_4!=""])
    t5=gsub('(.).*','\\\\1',x$ITS2.DGGE.dominant.genotype[x$Time.point==5])
    t1_5=unique(c(t1_4,t5))
    t1_5=t1_5[t1_5!=""]
    s='CD'
    ss=5
    #print(t1_5)
    if (length(t1_5)>1) {s='CD';ss=2;}
    if (length(t1_5)==1 & t1_5=='C') {s='CC';ss=3}
    if (length(t1_5)==1 & t1_5=='D') {s='DD';ss=1}
    xx = subset(x, Time.point==5)
    
    #data.frame(paste(t1_5,collapse=''),s)
    data.frame(Population=xx$Population, s,ss,Mortality=xx$Mortality)
}) %>% data.frame()

dd=dd %>% arrange(Population, Mortality, ss)

#data.T5 <- arrange(droplevels(subset(data, Time.point==5)), Population, !is.na(Mortality),100-Mortality)

#unique(data.T5$sample)
data.sort <- data
#data.sort$sample=factor(data.sort$sample, levels=unique(data.T5$sample))
data.sort$sample=factor(data.sort$sample, levels=rev(unique(dd$sample)))
ggplot(data.sort, aes(y=sample, x=Time.point)) +
    geom_tile(aes(color=D.C.ratio..20bin., fill=D.C.ratio..20bin.), size=0.5) +
        scale_y_discrete('')+
                        scale_x_discrete('', breaks=1:6,labels=format(unique(data$Date), format="%b\\n%Y"))+
        scale_fill_gradient('D:C ratio',low='grey',high='black', na.value=NA)+
        scale_color_gradient('D:C ratio', low='black',high='black', na.value=NA,guide=FALSE)+
            theme_classic()+
                theme(legend.position=c(1,0), legend.justification=c(1,0))
#                    facet_wrap(~Population)

#missing.samples <- subset(data, ITS2.DGGE.dominant.genotype=='' & !is.na(Mortality) & Mortality<100)
missing.samples <- subset(data.sort,is.na(data$D.C.ratio..20bin.) & (data$Mortality<100 | is.na(data$Mortality)))
#data.sort=data_old
data_old=data.sort
all <- data.sort

all$D.C.ratio..20bin. <- ifelse(is.na(data$D.C.ratio..20bin.) & (data$Mortality<100 | is.na(data$Mortality)),0,all$D.C.ratio..20bin.)
#subset(data, is.na(D.C.ratio..20bin.) & (Mortality<100 | is.na(Mortality)))

#& ITS2.DGGE.dominant.genotype=='' & !is.na(Mortality) & Mortality<100)

#data_old=data
#all <- data
#all$D.C.ratio..20bin. <- ifelse(is.na(all$D.C.ratio..20bin.) & all$ITS2.DGGE.dominant.genotype=='' & !is.na(all$Mortality) & all$Mortality<100, 0, all$D.C.ratio..20bin.)
data.sort=all
data.sort <- rbind(data.sort,
                   data.frame(Population='Davies', Time.point='3', sample='D29',S.H..clade.C.=NA, S.H..clade.D.=NA, D.Cratio=NA,
                              D.C.ratio..20bin.=0, D.C.ratio..4bin.=NA, Mortality=NA, ITS2.DGGE.dominant.genotype='',
                              Upgraded.D=NA, Downgraded.D=NA, Unchanged=NA,       Date=NA, S.H..C.D.=NA,       Dt=NA)
                   )
p1 <- ggplot(subset(data.sort,Population=='Keppels'), aes(y=sample, x=Time.point)) +
    geom_tile(aes(color=D.C.ratio..20bin., fill=D.C.ratio..20bin.), size=0.5) +
#        geom_text(data=subset(missing.samples,Population=='Keppels'), aes(label='////////////////'))+
            geom_text(data=subset(missing.samples,Population=='Keppels'), aes(label='XXXXXXX'), color='grey')+
                    geom_tile(aes(color=D.C.ratio..20bin.), fill=NA,size=0.5, show.legend = FALSE) +
            scale_y_discrete('Keppels')+
            scale_x_discrete('', breaks=1:6,labels=format(unique(data$Date), format="%b\\n%Y"))+
        scale_fill_gradient('D:C ratio',low='white',high='black', na.value=NA, guide=FALSE)+
        scale_color_gradient('D:C ratio', low='black',high='black', na.value=NA,guide=FALSE)+
            theme_classic()+
                theme(axis.title.x=element_blank(),
                      panel.margin.x=unit(0,'lines'),
                                            axis.ticks.x=element_blank(),
                      plot.margin=unit(c(0,0,0,0),'lines'),axis.text.x=element_blank(),
                      legend.position=c(1,0), legend.justification=c(1,0))#+facet_grid(Population~1, as.table=FALSE)


p2 <- ggplot(subset(data.sort,Population=='Davies'), aes(y=sample, x=Time.point)) +
    geom_tile(aes(color=D.C.ratio..20bin., fill=D.C.ratio..20bin.), size=0.5) +
#        annotate(geom='text', y='D29',x='3',
                                        #        label='////////////////')+
        annotate(geom='text', y='D29',x='3', label='XXXXXXX', color='grey')+
#            geom_text(data=subset(missing.samples,Population=='Davies'), aes(label='////////////////'))+
                            geom_text(data=subset(missing.samples,Population=='Davies'), aes(label='XXXXXXX'), color='grey')+
                    geom_tile(aes(color=D.C.ratio..20bin.), fill=NA,size=0.5, show.legend = FALSE) +

        scale_y_discrete('Davies')+
                    scale_x_discrete('', breaks=1:6,labels=format(unique(data$Date), format="%b\\n%Y"))+
        scale_fill_gradient('D:C ratio',low='white',high='black', na.value=NA)+
        scale_color_gradient('D:C ratio', low='black',high='black', na.value=NA,guide=FALSE)+
            theme_classic()+
                theme(axis.title.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      panel.margin.x=unit(0,'lines'),
                      plot.margin=unit(c(0,0,0,0),'lines'),axis.text.x=element_blank(),
                      legend.position=c(1,0), legend.justification=c(1,0))#+facet_grid(Population~1, as.table=FALSE)


p3 <- ggplot(subset(data.sort,Population=='Magnetic Island'), aes(y=sample, x=Time.point)) +
    geom_tile(aes(color=D.C.ratio..20bin., fill=D.C.ratio..20bin.), size=0.5) +
        scale_y_discrete('Magnetic Island')+
            geom_text(data=subset(missing.samples,Population=='Magnetic Island'), aes(label='XXXXXXX'), color='grey')+
                    geom_tile(aes(color=D.C.ratio..20bin.), fill=NA,size=0.5, show.legend = FALSE) +
                   scale_x_discrete('', breaks=1:6,labels=format(unique(data$Date), format="%b\\n%Y"))+
        scale_fill_gradient('D:C ratio',low='white',high='black', na.value=NA, guide=FALSE)+
        scale_color_gradient('D:C ratio', low='black',high='black', na.value=NA,guide=FALSE)+
            theme_classic()+
                theme(axis.title.x=element_blank(),panel.margin.x=unit(0,'lines'),
                      plot.margin=unit(c(0,0,1,0),'lines'),
                      legend.position=c(1,0), legend.justification=c(1,0))#+facet_grid(Population~1, as.table=FALSE)

#grid.arrange(p1,p2,p3,ncol=1,bottom='asdf')
cowplot:::plot_grid(p1,p2,p3, ncol=1, align='v', rel_heights = c(32,33,23))



pdf(file='figures/Figure7.pdf', width=175/25.4, height=300/25.4, onefile=FALSE)
cowplot:::plot_grid(p1,p2,p3, ncol=1, align='v', rel_heights = c(32,33,23))
dev.off()


#system('pdf2ps -dLanguageLevel=3 figures/Figure2.pdf figures/Figure2.eps')
system('pdftops -level3 -eps figures/Figure7.pdf figures/Figure7a.eps')
#system('convert figures/Figure1.pdf figures/Figure7.eps')

postscript(file='figures/Figure7p.eps', width=175/25.4, height=300/25.4, horizontal=FALSE, paper='special', fonts='Times')
cowplot:::plot_grid(p1,p2,p3, ncol=1, align='v', rel_heights = c(32,33,23))
dev.off()

library(Cairo)
CairoFonts('Nimbus Roman No9 L')
CairoPS(file='figures/Figure7.eps', width=175/25.4, height=300/25.4, horizontal=FALSE, paper='special',pointsize=8)
cowplot:::plot_grid(p1,p2,p3, ncol=1, align='v', rel_heights = c(32,33,23))
dev.off()

cairo_ps(file='figures/Figure7.eps', width=175/25.4, height=300/25.4)
cowplot:::plot_grid(p1,p2,p3, ncol=1, align='v', rel_heights = c(32,33,23))
dev.off()


## Now for the narrow version
p1 <- ggplot(subset(data.sort,Population=='Keppels'), aes(y=sample, x=Time.point)) +
    geom_tile(aes(color=D.C.ratio..20bin., fill=D.C.ratio..20bin.), size=0.5) +
#        geom_text(data=subset(missing.samples,Population=='Keppels'), aes(label='////////////////'))+
            geom_text(data=subset(missing.samples,Population=='Keppels'), aes(label='X'), color='grey')+
                    geom_tile(aes(color=D.C.ratio..20bin.), fill=NA,size=0.5, show.legend = FALSE) +
            scale_y_discrete('Keppels')+
            scale_x_discrete('', breaks=1:6,labels=format(unique(data$Date), format="%b\\n%Y"))+
        scale_fill_gradient('D:C ratio',low='white',high='black', na.value=NA, guide=FALSE)+
        scale_color_gradient('D:C ratio', low='black',high='black', na.value=NA,guide=FALSE)+
            theme_classic()+
                theme(axis.title.x=element_blank(),
                      panel.margin.x=unit(0,'lines'),
                                            axis.ticks.x=element_blank(),
                      plot.margin=unit(c(0,0,0,0),'lines'),axis.text.x=element_blank(),
                      legend.position=c(1,0), legend.justification=c(1,0))#+facet_grid(Population~1, as.table=FALSE)


p2 <- ggplot(subset(data.sort,Population=='Davies'), aes(y=sample, x=Time.point)) +
    geom_tile(aes(color=D.C.ratio..20bin., fill=D.C.ratio..20bin.), size=0.5, show.legend = FALSE) +
#        annotate(geom='text', y='D29',x='3',
                                        #        label='////////////////')+
        annotate(geom='text', y='D29',x='3', label='X', color='grey')+
#            geom_text(data=subset(missing.samples,Population=='Davies'), aes(label='////////////////'))+
                            geom_text(data=subset(missing.samples,Population=='Davies'), aes(label='X'), color='grey')+
                    geom_tile(aes(color=D.C.ratio..20bin.), fill=NA,size=0.5, show.legend = FALSE) +

        scale_y_discrete('Davies')+
                    scale_x_discrete('', breaks=1:6,labels=format(unique(data$Date), format="%b\\n%Y"))+
        scale_fill_gradient('D:C ratio',low='white',high='black', na.value=NA, guide=FALSE)+
        scale_color_gradient('D:C ratio', low='black',high='black', na.value=NA,guide=FALSE)+
            theme_classic()+
                theme(axis.title.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      panel.margin.x=unit(0,'lines'),
                      plot.margin=unit(c(0,0,0,0),'lines'),axis.text.x=element_blank(),
                      legend.position=c(1,0), legend.justification=c(1,0))#+facet_grid(Population~1, as.table=FALSE)

dts=unique(data.sort$Date)[1:6]
p3 <- ggplot(subset(data.sort,Population=='Magnetic Island'), aes(y=sample, x=Time.point)) +
    geom_tile(aes(color=D.C.ratio..20bin., fill=D.C.ratio..20bin.), size=0.5) +
        scale_y_discrete('Magnetic Island')+
            geom_text(data=subset(missing.samples,Population=='Magnetic Island'), aes(label='X'), color='grey')+
                    geom_tile(aes(color=D.C.ratio..20bin.), fill=NA,size=0.5, show.legend = FALSE) +
                        #scale_x_discrete('',
                        #breaks=1:6,labels=format(unique(data$Date),
                                        #format="%b\\n%Y"))+
                        scale_x_discrete('',breaks=1:6, labels=format(unique(data$Date), format='%b %y'))+
#                        scale_x_discrete('', breaks=1:6,labels=paste0(strtrim(format(unique(data$Date), format="%b"),1),'\\n',format(unique(data$Date), format="%y")))+
#                        scale_x_discrete('', breaks=1:6,labels=c(
#                                                            paste0(strtrim(format(dts[1],format="%b"),1),'\\n',format(dts[1],format="%y")),
#                                                            strtrim(format(dts[c(-1,-5)],format="%b"),1),
#                                                            paste0(strtrim(format(dts[5],format="%b"),1),'\\n',format(dts[5],format="%y"))))+
#                            scale_x_discrete('')+
        scale_fill_gradient('D:C ratio',low='white',high='black', na.value=NA, guide=FALSE)+
        scale_color_gradient('D:C ratio', low='black',high='black', na.value=NA,guide=FALSE)+
            theme_classic()+
                theme(axis.title.x=element_blank(),panel.margin.x=unit(0,'lines'),
                      axis.text.x=element_text(angle=90,vjust=0.5),
                      plot.margin=unit(c(0,0,1,0),'lines'),
                      legend.position=c(1,0), legend.justification=c(1,0))#+facet_grid(Population~1, as.table=FALSE)

#grid.arrange(p1,p2,p3,ncol=1,bottom='asdf')
cowplot:::plot_grid(p1,p2,p3, ncol=1, align='v', rel_heights = c(32,33,23))



pdf(file='figures/Figure7a.pdf', width=40/25.4, height=300/25.4, onefile=FALSE)
cowplot:::plot_grid(p1,p2,p3, ncol=1, align='v', rel_heights = c(32,33,23))
dev.off()


#system('pdf2ps -dLanguageLevel=3 figures/Figure2.pdf figures/Figure2.eps')
system('pdftops -level3 -eps figures/Figure7a.pdf figures/Figure7aa.eps')
#system('convert figures/Figure1.pdf figures/Figure7a.eps')

postscript(file='figures/Figure7ap.eps', width=40/25.4, height=300/25.4, horizontal=FALSE, paper='special', fonts='Times')
cowplot:::plot_grid(p1,p2,p3, ncol=1, align='v', rel_heights = c(32,33,23))
dev.off()

library(Cairo)
CairoFonts('Nimbus Roman No9 L')
CairoPS(file='figures/Figure7a.eps', width=40/25.4, height=300/25.4, horizontal=FALSE, paper='special',pointsize=8)
cowplot:::plot_grid(p1,p2,p3, ncol=1, align='v', rel_heights = c(32,33,23))
dev.off()

cairo_ps(file='figures/Figure7a.eps', width=40/25.4, height=300/25.4)
cowplot:::plot_grid(p1,p2,p3, ncol=1, align='v', rel_heights = c(32,33,23))
dev.off()


system('zip -ru FinalScriptAndFigures.zip *.R figures/*.pdf figures/*.eps')

