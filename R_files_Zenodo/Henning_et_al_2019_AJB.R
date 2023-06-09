library(car)
library(ggplot2)
library(multcomp)
library(e1071)
library(vegan)
library(compute.es)
library(ade4)
library(fBasics)
library(grid)
library(gridExtra)

########################################################################################################
#            Extra FUNCTIONS FOR GGPLOTS              #
########################################################################################################

stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="black", geom=geom, width=0.2, ...)
}  

# MULTIPLOT
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




################################################################################
################################################################################

################################################################################
# Single strain incoulations
################################################################################
################################################################################

################################################################################
singletons<-read.csv("Henning_et_al_2019_singleton_growth_data.csv")
names(singletons)


#############################
# Figure S1 - singleton biomass allocation patterns
#############################
library(lme4)
library(lmerTest)
mod1<-lmer(per_leaf~Treatment + (1|EP.run), data=singletons)
summary(mod1)
Anova(mod1)
difflsmeans(mod1)

mod2<-lmer(per_stem~Treatment + (1|EP.run), data=singletons)
summary(mod2)
Anova(mod2)
difflsmeans(mod2)

mod3<-lmer(per_root~Treatment + (1|EP.run), data=singletons)
summary(mod3)
Anova(mod3)
difflsmeans(mod3)

mod4<-lmer(SPAD~Treatment + (1|EP.run), data=singletons)
summary(mod4)
Anova(mod4)
difflsmeans(mod4)

names(singletons)
LMF <- ggplot(singletons, aes(x = Treatment, y = per_leaf))
SPADcontent <- ggplot(singletons, aes(x = Treatment, y = SPAD))
SMF <- ggplot(singletons, aes(x = Treatment, y = per_stem))
RMF <- ggplot(singletons, aes(x =Treatment, y = per_root))


singletons$Treatment
singletons$Treatment <- factor(singletons$Treatment, levels = c("C", "BT03", "GM17","GM30","GM41"))
SPADfig<- SPADcontent + geom_boxplot(fatten=2, lwd=1.3, color="black", fill="white")+theme_bw(base_size=13) + ylab(expression(paste("Chlorophyll Content (SPAD)")))  +  ylim(20,40) + 
  theme( panel.grid.major = element_blank(),
         plot.margin=unit(c(1,1,7,7), "pt"),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position='none',
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())  + geom_point(aes(size=1.4))


LMFfig<- LMF + geom_boxplot(fatten=2, lwd=1.3, color="black", fill="white")+theme_bw(base_size=13) + ylab(expression(paste("LMF (mg mg"^-1,")")))  +  ylim(30, 70) + 
  theme( panel.grid.major = element_blank(),
         plot.margin=unit(c(1,1,7,7), "pt"),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position='none',
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())  + geom_point(aes(size=1.4))

RMFfig <-RMF + geom_boxplot(fatten=2, lwd=1.3, color="black", fill="white")+theme_bw(base_size=13) + ylab(expression(paste("RMF (mg mg"^-1,")"))) + ylim(20, 60) + 
  theme( panel.grid.major = element_blank(),
         plot.margin=unit(c(1,1,10,7), "pt"),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position='none',
         axis.title.x=element_blank())  + geom_point(aes(size=1.4) ) 

SMFfig<-SMF + geom_boxplot(fatten=2, lwd=1.3, color="black", fill="white") +theme_bw(base_size=13)+ ylab(expression(paste("SMF (mg mg"^-1,")")))  + ylim(5,25)+ 
  theme( panel.grid.major = element_blank(),
         plot.margin=unit(c(1,1,7,7), "pt"),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position='none',
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())  + geom_point(aes(size=1.4) )



p1 <- ggplot_gtable(ggplot_build(SPADfig))
p2 <- ggplot_gtable(ggplot_build(LMFfig))
p3 <- ggplot_gtable(ggplot_build(SMFfig))
p4 <- ggplot_gtable(ggplot_build(RMFfig))

maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3], p3$widths[2:3], p4$widths[2:3])

p1$widths[2:3] <- maxWidth
p2$widths[2:3] <- maxWidth
p3$widths[2:3] <- maxWidth
p4$widths[2:3] <- maxWidth

tiff(filename="Henning_et_al_fig_S1.tiff", width=6.61, height=10.15, units='in', res=300)
grid.arrange(p1,p2,p3,p4,ncol = 1)
dev.off()


#############################################################################################
# nmds of looking at the phenotype differences of single-strain incoulations
# #figure 1 
#############################################################################################
singletons2<-singletons[,c("Treatment","per_stem", "per_leaf" ,"per_root",  "SPAD" )]
singletons2$Treatment
singletons3<-singletons[,c("per_stem", "per_leaf" ,"per_root",  "SPAD" )]
popcols<-c(rep("black",8),  rep("dodgerblue2", 8), rep("darkorange", 7), rep("gold", 8), rep("firebrick1", 10))

rdapca<-rda(singletons3, scale=TRUE)
scl=2
pca_scores<-scores(rdapca)
summary(rdapca)
names(rdapca)
rdapca$CA$u

str(scrs, max=1)
scrs<-scores(rdapca, display=c("sites", "species"), scaling=scl)
xlim <-with(scrs, range(species[,1]-0.25, sites[,1]))
ylim<-with(scrs, range(species[,2], sites[,2]))

tiff(filename="Henning_com_swi_fig1.tiff", width=10.66, height=8.35, units='in', res=300)
par(mar=c(6,6,1,1))
plot.new()
plot.window(xlim=xlim, ylim=ylim, asp=1)

abline(h=0, lty="dotted", lwd=2)
abline(v=0, lty="dotted", lwd=2)
with(arrows(0,0,0.9088733, -1.2913125,lwd=7, col='gray82'))
with(arrows(0,0,-1.7212708, -0.3455528,lwd=7, col='gray82'))
with(arrows(0,0,1.3285083,  1.1675822,lwd=7, col='gray82'))
with(arrows(0,0,0.9036858, -1.0759196,lwd=7, col='gray82'))
with(singletons3, points(scrs$sites, col=popcols, pch=19, cex=2))
with(scrs, text(species, labels=c("% Stem", "% Leaf", "% Root", "SPAD"), col="black", cex=2.2))
with(singletons3, legend("topleft", c("Control", "BT03", "GM30", "GM41", "GM17"), bty="n",col = c("black", "dodgerblue2","darkorange","gold", "firebrick1"), pch = 19, cex=2))
#ordihull(rdapca,popcols,draw=c("lines"), lwd=3, scaling=scl)
ordiellipse(rdapca, popcols, display="sites", kind = c("sd"), scaling=scl,draw = c("polygon"),
            col = c("black","darkorange", "dodgerblue2", "firebrick1","gold"),  label = FALSE, border = NULL, lty = NULL, lwd=2)
axis(side = 1, lwd=3, cex.axis=2)
axis(side = 2, lwd=3, cex.axis=2)
title(xlab = "PC 1 (50.3%)", ylab = "PC 2 (34.1%)", cex.lab=2)
box(lwd=4)
dev.off()


#####################################################################################
# Single strain inoculations as predictor of mixed communities
#####################################################################################

### simple predictive versus actual allocatiom model
allocation<-read.csv("Henning_et_al_2019_predictive_allocation_model_data.csv")
names(allocation)
allocation$Observation

t <- ggplot(allocation, aes(Observation,Stem_residual))
t +  geom_hline(yintercept=0)+ geom_point(aes(colour=Treatment)) +scale_color_brewer(palette="Dark2")+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background =element_blank(),axis.line.x =element_line(colour="black"), axis.line.y=element_line(colour="black"),axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold") )

p <- ggplot(allocation, aes(Observation,Leaf_residual))
p +  geom_hline(yintercept=0)+ geom_point(aes(colour=Treatment)) +scale_color_brewer(palette="Dark2")+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background =element_blank(),axis.line.x =element_line(colour="black"), axis.line.y=element_line(colour="black"),axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold") )

q <- ggplot(allocation, aes(Observation,Root_residual))
q +  geom_hline(yintercept=0)+ geom_point(aes(colour=Treatment)) +scale_color_brewer(palette="Dark2")+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background =element_blank(),axis.line.x =element_line(colour="black"), axis.line.y=element_line(colour="black"),axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold") )


mod3<-lm(Stem_residual~Treatment, data=allocation)
Anova(mod3)
summary(mod3)

t.test(subset(allocation$Stem_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$Stem_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$Stem_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$Stem_residual, allocation$Treatment=="30,41,bt"), mu=0)


mod1<-lm(Leaf_residual~Treatment, data=allocation)
Anova(mod1)
summary(mod1)

t.test(subset(allocation$Leaf_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$Leaf_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$Leaf_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$Leaf_residual, allocation$Treatment=="30,41,bt"), mu=0)


mod2<-lm(Root_residual~Treatment, data=allocation)
Anova(mod2)
summary(mod2)

t.test(subset(allocation$Root_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$Root_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$Root_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$Root_residual, allocation$Treatment=="30,41,bt"), mu=0)



################################################################################
################################################################################
################################################################################
# figures of  3 member communities
################################################################################
################################################################################
################################################################################



CSR<-read.csv("Henning_et_al_2019_community_data.csv")
CSR
names(CSR)

CSR$logSPAD<-log(CSR$SPAD+1)
CSR$logchange<-(log(CSR$change.in.wet.biomass..mg.+1))
CSR$logRL<-(log(CSR$root.leaf+1))
CSR$logBT<-(sqrt(CSR$BTcellswet.g.biomass))
CSR$log41<-(log(CSR$X41_cellswet.g.biomass+1))
CSR$log30<-(log(CSR$X30_cellswet.g.biomass+1))
CSR$log17<-(log(CSR$X17_cellswet.g.biomass+1))
CSR$logdryBT<-(sqrt(CSR$BT_cellsdrygbiomass))
CSR$logdry41<-(log(CSR$X41_cellsdrygbiomas+1))
CSR$logdry30<-(log(CSR$X30_cellsdrygbiomass+1))
CSR$logdry17<-(log(CSR$X17_cellsdrygbiomass+1))




####################################################################################################
####################################################################################################
####             Looking at effects on plant biomass allocation & chlorophyll content           #####
####################################################################################################
####################################################################################################

##############################################################################
#### FIGURE 2
##############################################################################
LMF <- ggplot(CSR, aes(x = Trt.code, y = Leaf))
SPADcontent <- ggplot(CSR, aes(x = Trt.code, y = SPAD))
SMF <- ggplot(CSR, aes(x = Trt.code, y = Stem))
RMF <- ggplot(CSR, aes(x =Trt.code, y = Root))
biomass<- ggplot(CSR, aes(x =Trt.code, y = total.biomass))

SPADfig<- SPADcontent + geom_boxplot(fatten=2, lwd=1.3, color="black", fill="white")+theme_bw(base_size=21) + ylab(expression(paste("Chlorophyll Content (SPAD)")))  +  ylim(20,32) + 
  theme( panel.grid.major = element_blank(),
         plot.margin=unit(c(1,1,7,7), "pt"),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position='none',
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())  + geom_point(aes(size=1.8))

totbio<- biomass + geom_boxplot(fatten=2, lwd=1.3, color="black", fill="white") + ylab(expression(paste("Total biomass (mg)")))  +theme_bw(base_size=21) +
  theme(panel.grid.major = element_blank(),
        plot.margin=unit(c(1,1,7,7), "pt"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  + geom_point(aes(size=1.8))

LMFfig<- LMF + geom_boxplot(fatten=2, lwd=1.3, color="black", fill="white")+theme_bw(base_size=21) + ylab(expression(paste("LMF (mg mg"^-1,")")))  +  ylim(45, 70) + 
  theme( panel.grid.major = element_blank(),
         plot.margin=unit(c(1,1,7,7), "pt"),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position='none',
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())  + geom_point(aes(size=1.8))

RMFfig <-RMF + geom_boxplot(fatten=2, lwd=1.3, color="black", fill="white")+theme_bw(base_size=21) + ylab(expression(paste("RMF (mg mg"^-1,")"))) + ylim(20, 60) + 
  theme( panel.grid.major = element_blank(),
         plot.margin=unit(c(1,1,10,7), "pt"),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position='none',
         axis.title.x=element_blank())  + geom_point(aes(size=1.8) ) + scale_x_discrete (breaks=c("A","B","C","D"), labels=c("GM17,GM30,GM41", "BT03,GM17,GM30", "BT03,GM17,GM41","BT03,GM30,GM41"))

SMFfig<-SMF + geom_boxplot(fatten=2, lwd=1.3, color="black", fill="white") +theme_bw(base_size=21)+ ylab(expression(paste("SMF (mg mg"^-1,")")))  + ylim(5,16)+ 
  theme( panel.grid.major = element_blank(),
         plot.margin=unit(c(1,1,7,7), "pt"),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position='none',
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())  + geom_point(aes(size=1.8) )

predpoints<- data.frame(allocation [,c("Trt.code","predictedstemalloc", "predictedleafalloc", "predictedrootalloc","per_root","per_leaf","per_stem" )])
coef(lm(predictedstemalloc ~ Trt.code, data = predpoints))
coef(lm(predictedleafalloc ~ Trt.code, data = predpoints))
coef(lm(predictedrootalloc ~ Trt.code, data = predpoints))

rootpredicfig<-RMFfig + geom_point(aes( x= predpoints$Trt.code, y= predpoints$predictedrootalloc), fill = "firebrick2", shape=23, size=8)
Leafpredicfig<-LMFfig + geom_point(aes( x= predpoints$Trt.code, y= predpoints$predictedleafalloc), fill = "firebrick2",shape=23,size=8)
Stempredicfig<-SMFfig + geom_point(aes( x= predpoints$Trt.code, y= predpoints$predictedstemalloc), fill = "firebrick2",shape=23,size=8)

rootpredicfig<-RMFfig
Leafpredicfig<-LMFfig 
Stempredicfig<-SMFfig 


p1 <- ggplot_gtable(ggplot_build(rootpredicfig))
p2 <- ggplot_gtable(ggplot_build(Leafpredicfig))
p3 <- ggplot_gtable(ggplot_build(Stempredicfig))
p4 <- ggplot_gtable(ggplot_build(totbio))
maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3], p3$widths[2:3])
maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3], p3$widths[2:3], p4$widths[2:3])

p1$widths[2:3] <- maxWidth
p2$widths[2:3] <- maxWidth
p3$widths[2:3] <- maxWidth
p4$widths[2:3] <- maxWidth

tiff(filename="Henning_allocation_fig_comm_switch_Fig2.tiff", width=10.67, height=9.24, units='in', res=300)
grid.arrange(p2,p3,p1,ncol = 1)
dev.off()


#####################################################################################################
#####################################################################################################
####                  Figure 3: Nutrient figure                                                            ####
#####################################################################################################
#####################################################################################################
shootn<- ggplot(CSR, aes(x =Trt.code, y = ShootN))

shootNfig<- shootn + geom_boxplot(fatten=2, lwd=1.2, color="black", fill="white") +theme_bw(base_size=22) + ylab(expression(paste("Shoot N (mg N g"^-1,")")))  +  ylim(17,28) + 
  theme(panel.grid.major = element_blank(),
        plot.margin=unit(c(1,1,7,7), "pt"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  + geom_point(aes(size=1.8))


shootp<- ggplot(CSR, aes(x =Trt.code, y = ShootP))
shootPfig<- shootp + geom_boxplot(fatten=2, lwd=1.2, color="black", fill="white") +theme_bw(base_size=22) + ylab(expression(paste("Shoot P (mg P g"^-1,")")))  +  ylim(1,4.3) + 
  theme(panel.grid.major = element_blank(),
        plot.margin=unit(c(1,1,7,7), "pt"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank())  + geom_point(aes(size=1.8))+ scale_x_discrete (breaks=c("A","B","C","D"), labels=c("GM17,GM30,GM41", "BT03,GM17,GM30", "BT03,GM17,GM41","BT03,GM30,GM41"))

rootn<- ggplot(CSR, aes(x =Trt.code, y = RootN))

rootnfig<- rootn + geom_boxplot(fatten=2, lwd=1.2, color="black", fill="white") +theme_bw(base_size=22) + ylab(expression(paste("Root N (mg N g"^-1,")")))  +  ylim(4.9,8.7) + 
  theme(panel.grid.major = element_blank(),
        plot.margin=unit(c(1,1,7,16), "pt"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position='none',
        axis.title.x=element_blank())  + geom_point(aes(size=1.8)) + scale_x_discrete (breaks=c("A","B","C","D"), labels=c("GM17,GM30,GM41", "BT03,GM17,GM30", "BT03,GM17,GM41","BT03,GM30,GM41"))


p1 <- ggplot_gtable(ggplot_build(shootNfig))
p2 <- ggplot_gtable(ggplot_build(shootPfig))
maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3])
p1$widths[2:3] <- maxWidth
p2$widths[2:3] <- maxWidth

tiff(filename="Henning_nutrient_comm_switch_fig3.tiff", width=10.47, height=6.90, units='in', res=300)
grid.arrange(p1,p2,ncol=1 )
dev.off()

####################################################################################################################
######    What if we assume that bacterial strains have the same effect regardless of biomass?
##               Appendix S8
############################################################################################################################################
names(allocation)


rootpredicfig<-RMFfig + geom_point(aes( x= allocation$Trt.code, y= allocation$UC_predictedrootalloc), fill = "red", shape=23, size=8)
Leafpredicfig<-LMFfig + geom_point(aes( x= allocation$Trt.code, y= allocation$UC_predictedleafalloc), fill = "red",shape=23,size=8)
Stempredicfig<-SMFfig + geom_point(aes( x= allocation$Trt.code, y= allocation$UC_predictedstemalloc), fill = "red",shape=23,size=8)
SPADpredicfig<-SPADfig + geom_point(aes( x= allocation$Trt.code, y= allocation$UC_predicted_SPAD), fill = "red",shape=23,size=8)


p1 <- ggplot_gtable(ggplot_build(rootpredicfig))
p2 <- ggplot_gtable(ggplot_build(Leafpredicfig))
p3 <- ggplot_gtable(ggplot_build(Stempredicfig))
#p4 <- ggplot_gtable(ggplot_build(SPADpredicfig))

maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3], p3$widths[2:3])

p1$widths[2:3] <- maxWidth
p2$widths[2:3] <- maxWidth
p3$widths[2:3] <- maxWidth
p4$widths[2:3] <- maxWidth

f
###
names(allocation)
mod3<-lm(NCresidSTEM~Treatment, data=allocation)
Anova(mod3)
summary(mod3)

t.test(subset(allocation$Stem_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$Stem_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$Stem_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$Stem_residual, allocation$Treatment=="30,41,bt"), mu=0)

t.test(subset(allocation$UC_Stem_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$UC_Stem_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$UC_Stem_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$UC_Stem_residual, allocation$Treatment=="30,41,bt"), mu=0)



mod1<-lm(NCresidLEAF~Treatment, data=uncorrect)
Anova(mod1)
summary(mod1)

t.test(subset(allocation$Root_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$Root_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$Root_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$Root_residual, allocation$Treatment=="30,41,bt"), mu=0)

t.test(subset(allocation$UC_Leaf_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$UC_Leaf_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$UC_Leaf_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$UC_Leaf_residual, allocation$Treatment=="30,41,bt"), mu=0)

mod2<-lm(NCresidROOT~Treatment, data=uncorrect)
Anova(mod2)
summary(mod2)

t.test(subset(allocation$Stem_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$Stem_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$Stem_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$Stem_residual, allocation$Treatment=="30,41,bt"), mu=0)

t.test(subset(allocation$UC_Root_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$UC_Root_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$UC_Root_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$UC_Root_residual, allocation$Treatment=="30,41,bt"), mu=0)

# SPAD

t.test(subset(allocation$spad_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$spad_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$spad_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$spad_residual, allocation$Treatment=="30,41,bt"), mu=0)

t.test(subset(allocation$UC_spad_residual, allocation$Treatment=="17,30,41"), mu=0)
t.test(subset(allocation$UC_spad_residual, allocation$Treatment=="17,30,bt"), mu=0)      
t.test(subset(allocation$UC_spad_residual, allocation$Treatment=="17,41,bt"), mu=0)
t.test(subset(allocation$UC_spad_residual, allocation$Treatment=="30,41,bt"), mu=0)

#############################################################################################
#   Is there a relationship between total plant biomass and # of bacterial cells?
#############################################################################################
cellfig <- ggplot(CSR, aes(x =log(total.cells), y = total.biomass))
SMFfig<-cellfig+ geom_point() +theme_bw(base_size=13)+ ylab(expression(paste("Plant biomass")))  + 
  theme( panel.grid.major = element_blank(),
         plot.margin=unit(c(1,1,7,7), "pt"),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position='none') + geom_point(aes(size=1.4) )

modbiomass<-lm(total.biomass~log(total.cells),data=CSR)
Anova(modbiomass)


CSR$Pseudcells<-(CSR$X41_cellsdrygbiomass+CSR$X30_cellsdrygbiomass+CSR$X17_cellsdrygbiomass)
Pcellfig <- ggplot(CSR, aes(x =log(Pseudcells), y = Leaf))
Pcellfig+ geom_point() +theme_bw(base_size=13)+ ylab(expression(paste("Plant biomass")))  + 
  theme( panel.grid.major = element_blank(),
         plot.margin=unit(c(1,1,7,7), "pt"),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.position='none') + geom_point(aes(size=1.4) )

modbiomass<-lm(Root~log(total.cells),data=CSR)
Anova(modbiomass)


#######################################################################################################
# Mixed community statistical models

############################################################################################################################
#################### SPAD MODEL * NS
############################################################################################################################

mod<-lm(logSPAD~startwetwt+Trt.code, data=CSR)
moda<-lm(logSPAD~Trt.code, data=CSR)
modb<-lm(logSPAD~startwetwt, data=CSR)
anova(mod,moda,modb)

mod<-lm(logSPAD~Trt.code, data=CSR)
Anova(mod)
summary(mod)
plot(moda)

############################################################################################################################
### total.biomass MODEL * not sig
############################################################################################################################

mod<-lm(total.biomass~startwetwt+Trt.code, data=CSR)
moda<-lm(total.biomass~Trt.code, data=CSR)
Anova(mod)
Anova(moda)
summary(moda)

mod<-aov(total.biomass~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)

modi<-lm(total.biomass~startwetwt*Trt.code, data=CSR)
Anova(modi) ### NO INTERACTION
############################################################################################################################
### dry.leaf MODEL * not sig
############################################################################################################################

mod<-lm(dry.leaves~startwetwt+Trt.code, data=CSR)
moda<-lm(dry.leaves~Trt.code, data=CSR)
Anova(mod)
Anova(moda)
summary(mod)
summary(moda)

mod<-aov(dry.leaves~startwetwt+Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)

modi<-lm(dry.leaves~startwetwt*Trt.code, data=CSR)
Anova(modi) ### NO INTERACTION
names(CSR)
############################################################################################################################
### dry.stem MODEL * SIG
############################################################################################################################

mod<-lm(dry.stem~startwetwt+Trt.code, data=CSR)
moda<-lm(dry.stem~Trt.code, data=CSR)
Anova(mod)
Anova(moda)
summary(moda)

mod<-aov(dry.stem~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)

modi<-lm(dry.stem~Trt.code, data=CSR)
Anova(modi) ### NO INTERACTION
names(CSR)
postHocs<-glht(modi, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)

boxplot(dry.stem~Trt.code, data=CSR)

############################################################################################################################
### dry.root MODEL * sig
############################################################################################################################

mod<-lm(dry.root~Trt.code, data=CSR)
moda<-lm(dry.root~startwetwt+Trt.code, data=CSR)
Anova(moda)
Anova(mod)
summary(mod)


mod<-aov(dry.root~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)

boxplot(dry.root~Trt.code, data=CSR)

modi<-lm(dry.root~startwetwt*Trt.code, data=CSR)
Anova(modi) ### NO INTERACTION
names(CSR)


shoots<-CSR$dry.stem + CSR$dry.leaves
mod<-lm(shoots~Trt.code, data=CSR)
Anova(mod)
summary(mod)

############################################################################################################################
### percent leaf  MODEL * sig
############################################################################################################################

mod<-lm(per_leaf~startwetwt+Trt.code, data=CSR)
moda<-lm(per_leaf~Trt.code, data=CSR)
Anova(moda)
Anova(mod)
summary(moda)

mod<-aov(per_leaf~startwetwt+Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)

modi<-lm(per_leaf~startwetwt*Trt.code, data=CSR)
Anova(modi) ### NO INTERACTION

############################################################################################################################
### percent stem MODEL * not sig
############################################################################################################################

mod<-lm(per_stem~startwetwt+Trt.code, data=CSR)
moda<-lm(per_stem~Trt.code, data=CSR)
Anova(moda)
Anova(mod)
summary(moda)

moda<-aov(per_stem~Trt.code, data=CSR)
postHocs<-glht(moda, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)

modi<-lm(per_stem~startwetwt*Trt.code, data=CSR)
Anova(modi) ### NO INTERACTION

############################################################################################################################
### percent roots MODEL * sig
############################################################################################################################

mod<-lm(per_root~Trt.code, data=CSR)
Anova(mod)
summary(mod)
mods<-lm(per_root~startwetwt+Trt.code, data=CSR)
Anova(mods)
summary(mods)

boxplot(per_root~Trt.code, data=CSR)
mod<-aov(per_root~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)

modi<-lm(per_root~startwetwt*Trt.code, data=CSR)
Anova(modi) ### NO INTERACTION

########### SHOOOT N   ##############
mod<-lm(ShootN~Trt.code, data=CSR)
boxplot(ShootN~Trt.code, data=CSR)
Anova(mod)
summary(mod)

mod<-aov(ShootN~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)


######## SHOOT  P   ##############
boxplot(ShootP~Trt.code, data=CSR)
mod<-lm(ShootP~Trt.code, data=CSR)
Anova(mod)
summary(mod)

mod<-aov(ShootP~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)

##### SHOOT N:P ########
mod<-aov(ShootNP~Trt.code, data=CSR)
Anova(mod)

########### ROOT N   ##############
mod<-lm(RootN~Trt.code, data=CSR)
boxplot(RootN~Trt.code, data=CSR)
Anova(mod)
summary(mod)

mod<-aov(RootN~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)

######## ROOOT  P   ##############
boxplot(RootP~Trt.code, data=CSR)
mod<-lm(RootP~Trt.code, data=CSR)
Anova(mod)
summary(mod)
mod<-aov(ShootP~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)
confint(postHocs)


############################################################################################################################
############ Populus phenotype - Trait PCA ###### ##
############################################################################################################################
# did not add in the root nutrient data because there is too much missing.
CSRP3<-CSR[,c("per_root", "per_leaf", "per_stem")]
TraitPCA<-CSRP3[complete.cases(CSRP2),]
popcols<-c(rep("dodgerblue2", 6), rep("darkorange", 7), rep("gold", 7), rep("firebrick1", 7))

rdapca<-rda(TraitPCA, scale=TRUE)
pca_scores<-scores(rdapca)
summary(rdapca)
names(rdapca)
rdapca$CA$u
scores(rdapca)
str(scrs, max=1)
scrs<-scores(rdapca, display=c("sites", "species"), scaling=scl)
xlim <-with(scrs, range(species[,1], sites[,1]))
ylim<-with(scrs, range(species[,2], sites[,2]))

tiff(filename="Henning_com_swi_fig3b.tiff", width=9.78, height=6.37, units='in', res=300)
par(mar=c(6,6,1,1))
plot.new()
plot.window(xlim=c(-2,2), ylim=c(-2,1), asp=1)


abline(h=0, lty="dotted", lwd=2)
abline(v=0, lty="dotted", lwd=2)
with(scrs, arrows(0,0,-1.689596, -0.2986387, lwd=7, col='gray82'))
with(scrs, arrows(0,0,1.606746 , 0.6019035, lwd=7, col='gray82'))
with(scrs, arrows(0,0,1.224921, -1.2014531, lwd=7, col='gray82'))

with(scrs, text(species, labels=c("% root", "% leaf", "% stem"), col="black", cex=2.2))
with(TraitPCA, legend("bottomright", bty="n",c("17,30,41", "17,30,bt", "17,41,bt", "30,41,bt"), col = c("dodgerblue2","darkorange","gold", "firebrick1"), pch = 19, cex=2))
ordihull(rdapca,popcols,draw=c("lines"), lwd=0.5, scaling=scl)
ordiellipse(rdapca, popcols, display="sites", kind = c("sd"), scaling=scl,draw = c("polygon"),
            col = c("darkorange", "dodgerblue2", "firebrick1","gold"),  label = FALSE, border = NULL, lty = NULL, lwd=2)
with(TraitPCA, points(scrs$sites, col=popcols, pch=19, cex=2))
axis(side = 1, lwd=3, cex.axis=2)
axis(side = 2, lwd=3, cex.axis=2)
title(xlab = "PC 1 (78.5%)", ylab = "PC 2 (21.5%)", cex.lab=2)
box(lwd=4)
dev.off()


































#####################################################################################################
#####################################################################################################
####                  Microbial community abundancne qPCR results                                ####
#####################################################################################################
#####################################################################################################

####################################################################################################
############      cells per wet root weight        #################################################
####################################################################################################
BTfig<-ggplot(CSR, aes(x=Trt.code, y=logBT))
BTfig2<- BTfig + geom_boxplot(fatten=2, lwd=1.8) +  theme(panel.grid.major = element_blank(),
                                                                        plot.margin=unit(c(1,8,7,7), "pt"),
                                                                        panel.grid.minor = element_blank(),
                                                                        panel.border = element_blank(),
                                                                        panel.background = element_blank(),
                                                                        axis.text = element_text(size=rel(1.4)),
                                                                        legend.position='none',
                                                                        axis.title.y = element_text(size = rel(1.5)),
                                                                        axis.title.x=element_blank(),
                                                                        axis.text.x=element_blank(),
                                                                        axis.line.x =element_line(colour="black"), axis.line.y=element_line(colour="black"),
                                                                        axis.ticks.x=element_blank())  




a41fig<-ggplot(CSR, aes(x=Trt.code, y=log41))
a41fig2<- a41fig + geom_boxplot(fatten=2, lwd=1.8) +  theme(panel.grid.major = element_blank(),
                                                                    plot.margin=unit(c(1,8,7,7), "pt"),
                                                                    panel.grid.minor = element_blank(),
                                                                    panel.border = element_blank(),
                                                                    panel.background = element_blank(),
                                                                    axis.text = element_text(size=rel(1.4)),
                                                                    legend.position='none',
                                                                    axis.title.y = element_text(size = rel(1.5)),
                                                                    axis.title.x=element_blank(),
                                                                    axis.text.x=element_blank(),
                                                                    axis.line.x =element_line(colour="black"), axis.line.y=element_line(colour="black"),
                                                                    axis.ticks.x=element_blank())  

a30fig<-ggplot(CSR, aes(x=Trt.code, y=log30))
a30fig2<- a30fig + geom_boxplot(fatten=2, lwd=1.8) +  theme(panel.grid.major = element_blank(),
                                                            plot.margin=unit(c(1,8,7,7), "pt"),
                                                            panel.grid.minor = element_blank(),
                                                            panel.border = element_blank(),
                                                            panel.background = element_blank(),
                                                            axis.text = element_text(size=rel(1.4)),
                                                            legend.position='none',
                                                            axis.title.y = element_text(size = rel(1.5)),
                                                            axis.title.x=element_blank(),
                                                            axis.text.x=element_blank(),
                                                            axis.line.x =element_line(colour="black"), axis.line.y=element_line(colour="black"),
                                                            axis.ticks.x=element_blank())  

a17fig<-ggplot(CSR, aes(x=Trt.code, y=log17))
a17fig2<- a17fig + geom_boxplot(fatten=2, lwd=1.8) +  theme(panel.grid.major = element_blank(),
                                                            plot.margin=unit(c(1,8,7,7), "pt"),
                                                            panel.grid.minor = element_blank(),
                                                            panel.border = element_blank(),
                                                            panel.background = element_blank(),
                                                            axis.text = element_text(size=rel(1.4)),
                                                            legend.position='none',
                                                            axis.title.y = element_text(size = rel(1.5)),
                                                            axis.title.x = element_text(size = rel(1.5)),
                                                            axis.line.x =element_line(colour="black"), axis.line.y=element_line(colour="black"),
                                                            axis.ticks.x=element_blank())  


p1 <- ggplot_gtable(ggplot_build(BTfig2))
p2 <- ggplot_gtable(ggplot_build(a41fig2))
p3 <- ggplot_gtable(ggplot_build(a30fig2))
p4 <- ggplot_gtable(ggplot_build(a17fig2))

maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3], p3$widths[2:3], p4$widths[2:3])

p1$widths[2:3] <- maxWidth
p2$widths[2:3] <- maxWidth
p3$widths[2:3] <- maxWidth
p4$widths[2:3] <- maxWidth


grid.arrange(p1,p2,p3,p4,ncol = 1)

# bt03
mod<-lm(logBT~Trt.code, data=CSR)
Anova(mod)
summary(mod)

mod<-aov(logBT~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)

# log41
mod<-lm(log41~Trt.code, data=CSR)
Anova(mod)
summary(mod)

mod<-aov(log41~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)

# log30
mod<-lm(log30~Trt.code, data=CSR)
Anova(mod)
summary(mod)

mod<-aov(log30~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)

# log17
mod<-lm(log17~Trt.code, data=CSR)
Anova(mod)
summary(mod)

mod<-aov(log17~Trt.code, data=CSR)
postHocs<-glht(mod, linfct=mcp(Trt.code="Tukey"))
summary(postHocs)



############################################################################################################################
###########       PCA of microbial community composition     ####################################
############################################################################################################################
CSRP2<-CSR[,c("logBT","log41","log30","log17")]
CommPCA<-CSRP2[complete.cases(CSRP2),]
commPCA

popcols<-c(rep("dodgerblue2", 6), rep("darkorange", 7), rep("gold", 7), rep("firebrick1", 7))

rdapca<-rda(commPCA, scale=TRUE)
scl=2
pca_scores<-scores(rdapca)
summary(rdapca)
names(rdapca)
rdapca$CA$u

str(scrs, max=1)
scrs<-scores(rdapca, display=c("sites", "species"), scaling=scl)
xlim <-with(scrs, range(species[,1], sites[,1]))
ylim<-with(scrs, range(species[,2], sites[,2]))

tiff(filename="Henning_com_swi_fig3a.tiff", width=9.78, height=6.37, units='in', res=300)
par(mar=c(6,6,1,1))
plot.new()
plot.window(xlim=xlim, ylim=ylim, asp=1)

abline(h=0, lty="dotted", lwd=2)
abline(v=0, lty="dotted", lwd=2)
with(scrs, arrows(0,0,-1.0124858,0.1728750, lwd=7, col='gray82'))
with(scrs, arrows(0,0,0.9975537, -0.2233056, lwd=7, col='gray82'))
with(scrs, arrows(0,0,0.2482928, -1.4098896, lwd=7, col='gray82'))
with(scrs, arrows(0,0,0.9776126,  0.7649843, lwd=7, col='gray82'))
with(commPCA, points(scrs$sites, col=popcols, pch=19, cex=2))
with(scrs, text(species, labels=c("BT03", "GM41", "GM30", "GM17"), col="black", cex=2.2))
with(commPCA, legend("bottomright", c("17,30,41", "17,30,bt", "17,41,bt", "30,41,bt"), bty="n",col = c("dodgerblue2","darkorange","gold", "firebrick1"), pch = 19, cex=2))
ordihull(rdapca,popcols,draw=c("lines"), lwd=0.5, scaling=scl)
ordiellipse(rdapca, popcols, display="sites", kind = c("sd"), scaling=scl,draw = c("polygon"),
            col = c("darkorange", "dodgerblue2", "firebrick1","gold"),  label = FALSE, border = NULL, lty = NULL, lwd=2)
axis(side = 1, lwd=3, cex.axis=2)
axis(side = 2, lwd=3, cex.axis=2)
title(xlab = "PC 1 (29.8%)", ylab = "PC 2 (26.0%)", cex.lab=2)
box(lwd=4)
dev.off()



############################################################################################################################
#                Co-inertia incorporating nutrients          #
#####                   FIGURE 4                        ######
############################################################################################################################
library(ade4)
PCABIO<- dudi.coa(commPCA, scannf=FALSE)
PCAENZ<- dudi.pca(TraitPCA,  row.w=PCABIO$lw, scannf=FALSE)

coin1<-coinertia(PCABIO, PCAENZ, scannf=FALSE)
plot(coin1)

summary(coin1)
rv1 <- RV.rtest(PCABIO$tab, PCAENZ$tab, 999)
plot(rv1)
rv1

### FIGURE
traitval<-(coin1$l1)
commval<-(coin1$c1)

tiff(filename="Henning_com_swi_fig3c.tiff", width=8.78, height=5.37, units='in', res=300)
par(mfrow = c(1,2),mar=c(1,1,1,0) )
plot(c(-1:1), xlim=c(-1,1), type='n',xaxt='n',yaxt='n',ann=FALSE)
abline(h=0, lty=3, lwd=1.5)
abline(v=0, lty=3, lwd=1.5)
arrows(0,0,0.5312060,  0.4238555, lwd=5) #per_root
text(0.5312060,  0.4738555, "Root", cex=1.5)
arrows(0,0,-0.4337560, -0.6077320, lwd=5) # per_leaf
text(-0.4337560, -0.6577320, "Leaf", cex=1.5)
arrows(0,0,-0.7277884,  0.6715715, lwd=5) #per_stem
text(-0.7277884,  0.7215715, "Stem", cex=1.5)
box(lwd=3)

##### adding in communinity members
par(mar=c(1,0,1,1))
plot(c(-25:25), xlim=c(-25,25), type='n',xaxt='n',yaxt='n',ann=FALSE)
abline(h=0, lty=3, lwd=1.5)
abline(v=0, lty=3, lwd=1.5)
arrows(0,0,  -0.06857851,   0.03736871, lwd=5) #BT03
text( -3.56857851,   1.03736871, "BT03", cex=1.5)
arrows(0,0, 1.28808211, -22.06321097, lwd=5) #GM41
text(1.28808211, -23.66321097, "GM41", cex=1.5)
arrows(0,0, 14.63940329,  -3.68999502, lwd=5) #GM30
text( 23.63940329 , -3.68999502, "GM30", cex=1.5, pos=2)
arrows(0,0,15.64375439,   6.68045857, lwd=5) #GM17
text(19.94375439,   6.68045857, "GM17", cex=1.5)
box(lwd=3)
dev.off()

########################################################################################################

