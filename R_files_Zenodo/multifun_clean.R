#cleaned version of the code for the trophic complexity manuscript

#calling the file
setwd("C:/Users/krish/Documents/Columbia/StatMod/cedar_creek")
bef<-read.csv("bef_manuscript.csv")

#initialising
col=c("chartreuse4", "goldenrod3", "chocolate4", "lemonchiffon4")
library("png")
library("corrplot")

#standardising functions----------------------------------------
bef$abvbiom<-bef$TOTBIOM01-bef$ROOT01 #calculating aboveground biomass
biomax<-mean(bef$abvbiom[order(bef$abvbiom, decreasing=T)][1:5])
#recmax<-mean(bef$rslnc[order(bef$rslnc, decreasing=T)][1:5])
recmax<-mean(bef$ABOVETOT02[order(bef$ABOVETOT02, decreasing=T)][1:5]) #should we use raw biomass
rtmax<-mean(bef$ROOT01[order(bef$ROOT01, decreasing=T)][1:5])
flomax<-mean(bef$FLOWTIME[order(bef$FLOWTIME, decreasing=T)][1:5])


bef$biomstd<-(bef$abvbiom-min(bef$abvbiom))/(biomax-min(bef$abvbiom))
bef$flostd<-(bef$FLOWTIME -min(bef$FLOWTIME))/(flomax-min(bef$FLOWTIME))
bef$recstd<-(bef$ABOVETOT02 -min(bef$ABOVETOT02))/(recmax-min(bef$ABOVETOT02)) #should we use raw biomass
#bef$recstd<-(bef$rslnc -min(bef$rslnc))/(recmax-min(bef$rslnc))
bef$rtstd<-(bef$ROOT01 -min(bef$ROOT01))/(rtmax-min(bef$ROOT01))

#correlations between functions--------------------------
fn.cor<-cor(bef[,c(20, 23, 21, 22)])
colnames(fn.cor)<-c("aboveground\\nbiomass", "root\\nbiomass", "water\\nretention", "biomass\\nrecovery")
rownames(fn.cor)<-c("aboveground\\nbiomass", "root\\nbiomass", "water\\nretention", "biomass\\nrecovery")
library(tiff)
tiff("manuscript/eco_evo_rev/fn_corrs_old.tiff", width=4, height=4, units="in", res=300)
corrplot(fn.cor, method="ellipse", type="upper", addCoef.col = "black", tl.cex=0.7, tl.col="black")
dev.off()


#Spider plots-----------------------------------------
#first make mult, the dataframe
library(plyr)
mult_sum<-ddply(bef, .(SPP,INTERACT), summarize, meanbiom=mean(biomstd), meanrt=mean(rtstd), meanflo=mean(flostd), meanrec=mean(recstd))

colnames(mult_sum)<-c("SPP", "INTERACT", "Aboveground biomass", "Root biomass", "Water retention", "Biomass recovery")

mult_sum_all<-ddply(bef, .(INTERACT), summarize, meanbiom=mean(biomstd), meanrt=mean(rtstd), meanflo=mean(flostd), meanrec=mean(recstd))

colnames(mult_sum_all)<-c("SPP", "Aboveground biomass", "Root biomass", "Water retention", "Biomass recovery")


# Library
library(fmsb)

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!

mult_sum_1<-mult_sum[mult_sum$SPP==1,3:6]
mult_sum_1<-rbind(rep(1, 4), rep(0, 4), mult_sum_1)

mult_sum_2<-mult_sum[mult_sum$SPP==2,3:6]
mult_sum_2<-rbind(rep(1, 4), rep(0, 4), mult_sum_2)

mult_sum_4<-mult_sum[mult_sum$SPP==4,3:6]
mult_sum_4<-rbind(rep(1, 4), rep(0, 4), mult_sum_4)

mult_sum_8<-mult_sum[mult_sum$SPP==8,3:6]
mult_sum_8<-rbind(rep(1, 4), rep(0, 4), mult_sum_8)

mult_sum_16<-mult_sum[mult_sum$SPP==16,3:6]
mult_sum_16<-rbind(rep(1, 4), rep(0, 4), mult_sum_16)

mult_sum_all<-rbind(rep(1, 5), rep(0, 5), mult_sum_all)

# Plot 2: Same plot with custom features
#colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
#colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}
## END

colors<-col
colors_border<-c(colors[4], colors[2], colors[3], colors[1])
#colors_in<-c( rgb(c(as.numeric(col2rgb(colors_border[1])[,1]), 0.4)), rgb(c(as.numeric(col2rgb(colors_border[2])[,1]), 0.4)) , rgb(c(as.numeric(col2rgb(colors_border[3])[,1]), 0.4)), rgb(c(as.numeric(col2rgb(colors_border[4])[,1]), 0.4)) )
colors_in<-c(t_col(colors_border[1], 40), t_col(colors_border[2], 40), t_col(colors_border[3], 40),
             t_col(colors_border[4], 40))


par(mar=c(1, 1, 1, 1))
radarchart( mult_sum_all[,2:5]  , axistype=1 ,
            #custom polygon
            pcol=colors_border , 
            #pfcol=colors_in , 
            plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            #custom labels
            #vlabels=c("Aboveground\\nbiomass", "Below-\\nground\\nbiomass", 
            #          "Water\\nretention", "Biomass\\nrecovery")
            vlabels=NA
            #vlcex=1.2
)

xpos<-c(0, -1.05, 0, 1.05)
ypos<-c(1.15, 0, -1.15, 0)
labs<-c("Aboveground\\nbiomass", "Belowground\\nbiomass", 
        "Water\\nretention", "Biomass\\nrecovery")
text(xpos, ypos, labs, cex=1.3)

?text
colnames(mult_sum_all)
?radarchart

#png("spider_plot_fill.png", width=8, height=12, units="in", res=300)
tiff("manuscript/eco_evo_rev/spider_plot_nofill_2.tiff", width=8, height=12, units="in", res=300)
par(mfrow=c(3,2), mar=c(1,1,1,1))
xpos<-c(0, -1, 0, 1.05)
ypos<-c(1.15, 0, -1.15, 0)
labs<-c("Aboveground\\nbiomass", "Root\\nbiomass", 
        "Water\\nretention", "Biomass\\nrecovery")

radarchart( mult_sum_all[,2:5]  , axistype=1 ,
            #custom polygon
            pcol=colors_border , 
            #pfcol=colors_in , 
            plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            #custom labels
            vlabels=NA
            #vlcex=1.2
)
legend("topleft", "a) All", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)

radarchart( mult_sum_1  , axistype=1 ,
            #custom polygon
            pcol=colors_border , 
            #pfcol=colors_in , 
            plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            #custom labels
            vlabels=NA
            #vlcex=1.2
)
legend("topleft", "b) 1 species", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)

radarchart( mult_sum_2  , axistype=1 ,
            #custom polygon
            pcol=colors_border , 
            #pfcol=colors_in , 
            plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            #custom labels
            vlabels=NA
            #vlcex=1.2
)
legend("topleft", "c) 2 species", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)

radarchart( mult_sum_4  , axistype=1 ,
            #custom polygon
            pcol=colors_border , 
            #pfcol=colors_in , 
            plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            #custom labels
            vlabels=NA
            #vlcex=1.2
)
legend("topleft", "d) 4 species", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)

radarchart( mult_sum_8  , axistype=1 ,
            #custom polygon
            pcol=colors_border , 
            #pfcol=colors_in , 
            plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            #custom labels
            vlabels=NA
            #vlcex=1.2
)
legend("topleft", "e) 8 species", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)

radarchart( mult_sum_16  , axistype=1 ,
            #custom polygon
            pcol=colors_border , 
            #pfcol=colors_in , 
            plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            #custom labels
            vlabels=NA
            #vlcex=1.2
)
legend("topleft", "f) 16 species", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)


#legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
par(xpd=NA)
legend(x=-1.7, y=5.2, bty="n", legend=c("NONE", "ABV", "LIT", "BOTH"), col=colors, lty=1, lwd=3, horiz = F, cex=1.3)
dev.off()


#BEF curves------------------------------------------
library(lme4)

abvbiom_glm0<-glm(abvbiom~log(SPP)+INTERACT+(log(SPP)*INTERACT), data=bef, family=Gamma(link=inverse))
#abvbiom_glm0<-glm(abvbiom~log(SPP)+INTERACT+(log(SPP)*INTERACT), data=bef)
root_glm0<-glm(ROOT01~log(SPP)+INTERACT+(log(SPP)*INTERACT), data=bef, family=Gamma(link="inverse"))
flo_glm0<-glm(FLOWTIME~log(SPP)+INTERACT+(log(SPP)*INTERACT), data=bef, family=Gamma(link="inverse"))
rslnc_glm0<-glm(ABOVETOT02~log(SPP)+INTERACT+(log(SPP)*INTERACT), data=bef, family=Gamma(link="inverse"))
rslnc_glm0<-glm(ABOVETOT02~log(SPP)+INTERACT+(log(SPP)*INTERACT), data=bef)


abvbiom_glm1<-glm(abvbiom~log(SPP)+INS+LIT+(INS*LIT)+(log(SPP)*INS*LIT), data=bef, family=Gamma(link=inverse))
root_glm1<-glm(ROOT01~log(SPP)+INS+LIT+(INS*LIT)+(log(SPP)*INS*LIT), data=bef, family=Gamma(link=inverse))
flo_glm1<-glm(FLOWTIME~log(SPP)+INS+LIT+(INS*LIT)+(log(SPP)*INS*LIT), data=bef, family=Gamma(link=inverse))
rslnc_glm1<-glm(ABOVETOT02~log(SPP)+INS+LIT+(INS*LIT)+(log(SPP)*INS*LIT), data=bef, family=Gamma(link=inverse))


abvbiom_glm2<-glm(abvbiom~log(SPP)+INTERACT, data=bef, family=Gamma(link="inverse"))
root_glm2<-glm(ROOT01~log(SPP)+INTERACT, data=bef, family=Gamma(link="inverse"))
flo_glm2<-glm(FLOWTIME~log(SPP)+INTERACT, data=bef, family=Gamma(link="inverse"))
rslnc_glm2<-glm(ABOVETOT02~log(SPP)+INTERACT, data=bef, family=Gamma(link="inverse"))

abvbiom_glm3<-glm(abvbiom~log(SPP), data=bef, family=Gamma(link="inverse"))
root_glm3<-glm(ROOT01~log(SPP), data=bef, family=Gamma(link="inverse"))
flo_glm3<-glm(FLOWTIME~log(SPP), data=bef, family=Gamma(link="inverse"))
rslnc_glm3<-glm(ABOVETOT02~log(SPP), data=bef, family=Gamma(link="inverse"))


summary(abvbiom_glm2)
summary(root_glm2)
summary(flo_glm2)
summary(rslnc_glm2)

library(bbmle)
AICctab(abvbiom_glm0, abvbiom_glm1, abvbiom_glm2, abvbiom_glm3)
AICctab(root_glm0, root_glm1, root_glm2, root_glm3)
AICctab(flo_glm0, flo_glm1, flo_glm2, flo_glm3)
AICctab(rslnc_glm0, rslnc_glm1, rslnc_glm2, rslnc_glm3)

library(plyr)
summary_allfn<-ddply(bef, .(SPP, INTERACT), summarise, 
                     abv.mean=mean(abvbiom, na.rm=T), abv.se=sd(abvbiom)/sqrt(length(abvbiom)),
                     root.mean=mean(ROOT01, na.rm=T), root.se=sd(ROOT01)/sqrt(length(ROOT01)),
                     flo.mean=mean(FLOWTIME, na.rm=T), flo.se=sd(FLOWTIME)/sqrt(length(FLOWTIME)),
                     rslnc.mean=mean(ABOVETOT02, na.rm=T), rslnc.se=sd(ABOVETOT02)/sqrt(length(ABOVETOT02)))
summary_allfn$col<-character(nrow(summary_allfn))
summary_allfn$col[which(summary_allfn$INTERACT=="NONE")]<-col[1]
summary_allfn$col[which(summary_allfn$INTERACT=="INS")]<-col[2]
summary_allfn$col[which(summary_allfn$INTERACT=="LIT")]<-col[3]
summary_allfn$col[which(summary_allfn$INTERACT=="BOTH")]<-col[4]

summary_allfn$abv.up<-summary_allfn$abv.mean+summary_allfn$abv.se
summary_allfn$abv.dn<-summary_allfn$abv.mean-summary_allfn$abv.se
summary_allfn$root.up<-summary_allfn$root.mean+summary_allfn$root.se
summary_allfn$root.dn<-summary_allfn$root.mean-summary_allfn$root.se
summary_allfn$flo.up<-summary_allfn$flo.mean+summary_allfn$flo.se
summary_allfn$flo.dn<-summary_allfn$flo.mean-summary_allfn$flo.se
summary_allfn$rslnc.up<-summary_allfn$rslnc.mean+summary_allfn$rslnc.se
summary_allfn$rslnc.dn<-summary_allfn$rslnc.mean-summary_allfn$rslnc.se


library(coefplot)

#p1<-coefplot(abvbiom_glm1, zeroColor="coral3", zeroType=1, 
#             newNames=c("log(SPP):INS:LIT"="log(SPP):BOTH", "INS:LIT"="BOTH", "INS"="ABV", "log(SPP):INS"="log(SPP):ABV"), title="Aboveground biomass")+ theme_bw()
#p2<-coefplot(root_glm1, zeroColor="coral3", zeroType=1, 
#             newNames=c("log(SPP):INS:LIT"="log(SPP):BOTH", "INS:LIT"="BOTH", "INS"="ABV", "log(SPP):INS"="log(SPP):ABV"), title="Root biomass")+theme_bw()
#p3<-coefplot(flo_glm1, zeroColor="coral3", zeroType=1, 
#             newNames=c("log(SPP):INS:LIT"="log(SPP):BOTH", "INS:LIT"="BOTH", "INS"="ABV", "log(SPP):INS"="log(SPP):ABV"), title="Water retention")+theme_bw()
#p4<-coefplot(rslnc_glm1, zeroColor="coral3", zeroType=1, 
#             newNames=c("log(SPP):INS:LIT"="log(SPP):BOTH", "INS:LIT"="BOTH", "INS"="ABV", "log(SPP):INS"="log(SPP):ABV"), title="Biomass recovery")+theme_bw()

library(sjPlot)
set_theme(base=theme_classic())
p1<-plot_model(abvbiom_glm1, show.values = T, value.offset = 0.3, 
             axis.labels=c("log(SPP):INS:LIT"="log(SPP):BOTH", "INS:LIT"="BOTH", "INS"="ABV", "log(SPP):INS"="log(SPP):ABV"), title="Aboveground biomass")
p2<-plot_model(root_glm1, show.values = T, value.offset = 0.3, 
               axis.labels=c("log(SPP):INS:LIT"="log(SPP):BOTH", "INS:LIT"="BOTH", "INS"="ABV", "log(SPP):INS"="log(SPP):ABV"), title="Root biomass")
p3<-plot_model(flo_glm1, show.values = T, value.offset = 0.3, 
             axis.labels=c("log(SPP):INS:LIT"="log(SPP):BOTH", "INS:LIT"="BOTH", "INS"="ABV", "log(SPP):INS"="log(SPP):ABV"), title="Water retention")
p4<-plot_model(rslnc_glm1, show.values = T, value.offset = 0.3, 
             axis.labels=c("log(SPP):INS:LIT"="log(SPP):BOTH", "INS:LIT"="BOTH", "INS"="ABV", "log(SPP):INS"="log(SPP):ABV"), title="Biomass recovery")

library(gridExtra)
library(gridBase)
library(ggplot2)
library(grid)

tiff("manuscript/eco_evo_rev/figure2_smoothers.tiff", width=8, height=16, units="in", res=300)

plot.new()
grid.newpage()
pushViewport(viewport(layout=grid.layout(4,2)))

pushViewport(viewport(layout.pos.row=1, layout.pos.col = 1))
print(p1, newpage=F)

popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col = 1))
print(p2, newpage=F)
popViewport()
pushViewport(viewport(layout.pos.row=3, layout.pos.col = 1))
print(p3, newpage=F)
popViewport()
pushViewport(viewport(layout.pos.row=4, layout.pos.col = 1))
print(p4, newpage=F)
popViewport()

x<-c(1, 2, 4, 8, 16)
pushViewport(viewport(layout.pos.row=1, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
par(mgp=c(2,1,0))
plot(abv.mean~SPP, summary_allfn, col=summary_allfn$col, pch=16, xlab="No. of species", ylab="",  
     cex=1.8, ylim=c(min(summary_allfn$abv.dn), max(summary_allfn$abv.up)), cex.lab=1.6, xaxt="n")
title(ylab=expression(paste("Aboveground biomass (g m"^"-1", ")", sep="")), cex.lab=1.6)
arrows(summary_allfn$SPP, summary_allfn$abv.mean, summary_allfn$SPP, summary_allfn$abv.up, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
arrows(summary_allfn$SPP, summary_allfn$abv.mean, summary_allfn$SPP, summary_allfn$abv.dn, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
#curve(predict(abvbiom_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="NONE")), col="chartreuse4", add=T, lwd=4, lty=1)
#curve(predict(abvbiom_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="INS")), col="goldenrod3", add=T, lwd=4, lty=1)
#curve(predict(abvbiom_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="LIT")), col="chocolate4", add=T, lwd=4, lty=1)
#curve(predict(abvbiom_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="BOTH")), col="lemonchiffon4", add=T, lwd=4, lty=1)
#curve(predict(abvbiom_glm3, type="response", newdata=data.frame(SPP=x)), col="black", add=T, lwd=4, lty=1)
#lines(smooth.spline(summary_allfn$SPP[which(summary_allfn$INTERACT=="NONE")], summary_allfn$abv.mean[which(summary_allfn$INTERACT=="NONE")]), lwd=1.8, col=col[1])
#lines(smooth.spline(summary_allfn$SPP[which(summary_allfn$INTERACT=="INS")], summary_allfn$abv.mean[which(summary_allfn$INTERACT=="INS")]), lwd=1.8, col=col[2])
#lines(smooth.spline(summary_allfn$SPP[which(summary_allfn$INTERACT=="LIT")], summary_allfn$abv.mean[which(summary_allfn$INTERACT=="LIT")]), lwd=1.8, col=col[3])
#lines(smooth.spline(summary_allfn$SPP[which(summary_allfn$INTERACT=="BOTH")], summary_allfn$abv.mean[which(summary_allfn$INTERACT=="BOTH")]), lwd=1.8, col=col[4])
#curve(predict(abvbiom_glm3, type="response", newdata=data.frame(SPP=x)), col="black", add=T, lwd=2, lty=1)
legend(-3,130, bty="n", legend=c("NONE", "ABV", "LIT", "BOTH"), col=col, pch=16, cex=1.2, xpd=NA, horiz = T)
axis(side=1, at=x, labels=x, cex.lab=1.2)
popViewport()

pushViewport(viewport(layout.pos.row=2, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
par(mgp=c(2,1,0))
plot(root.mean~SPP, summary_allfn, col=summary_allfn$col, pch=16, xlab="No. of species", ylab="", 
     cex=1.8, ylim=c(min(summary_allfn$root.dn), max(summary_allfn$root.up)), cex.lab=1.6, xaxt="n")
title(ylab=expression(paste("Root biomass (g m"^"-1", ")", sep="")), cex.lab=1.6)
arrows(summary_allfn$SPP, summary_allfn$root.mean, summary_allfn$SPP, summary_allfn$root.up, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
arrows(summary_allfn$SPP, summary_allfn$root.mean, summary_allfn$SPP, summary_allfn$root.dn, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
#curve(predict(root_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="NONE")), col="chartreuse4", add=T, lwd=4, lty=1)
#curve(predict(root_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="INS")), col="goldenrod3", add=T, lwd=4, lty=1)
#curve(predict(root_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="LIT")), col="chocolate4", add=T, lwd=4, lty=1)
#curve(predict(root_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="BOTH")), col="lemonchiffon4", add=T, lwd=4, lty=1)
#curve(predict(root_glm3, type="response", newdata=data.frame(SPP=x)), col="black", add=T, lwd=2, lty=1)
axis(side=1, at=x, labels=x, cex=1.2)
popViewport()

pushViewport(viewport(layout.pos.row=3, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
par(mgp=c(2,1,0))
plot(flo.mean~SPP, summary_allfn, col=summary_allfn$col, pch=16, xlab="No. of species", ylab="",
     cex=1.8, ylim=c(min(summary_allfn$flo.dn), max(summary_allfn$flo.up)), cex.lab=1.6, xaxt="n")
title(ylab="Water retention (s)", cex.lab=1.6)
arrows(summary_allfn$SPP, summary_allfn$flo.mean, summary_allfn$SPP, summary_allfn$flo.up, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
arrows(summary_allfn$SPP, summary_allfn$flo.mean, summary_allfn$SPP, summary_allfn$flo.dn, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
#curve(predict(flo_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="NONE")), col="chartreuse4", add=T, lwd=4, lty=1)
#curve(predict(flo_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="INS")), col="goldenrod3", add=T, lwd=4, lty=1)
#curve(predict(flo_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="LIT")), col="chocolate4", add=T, lwd=4, lty=1)
#curve(predict(flo_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="BOTH")), col="lemonchiffon4", add=T, lwd=4, lty=1)
#curve(predict(flo_glm3, type="response", newdata=data.frame(SPP=x)), col="black", add=T, lwd=2, lty=1)
axis(side=1, at=x, labels=x, cex=1)
popViewport()

pushViewport(viewport(layout.pos.row=4, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
par(mgp=c(2,1,0))
plot(rslnc.mean~SPP, summary_allfn, col=summary_allfn$col, pch=16, xlab="No. of species", ylab="", 
     cex=1.8, ylim=c(min(summary_allfn$rslnc.dn), max(summary_allfn$rslnc.up)), cex.lab=1.6, xaxt="n")
title(ylab=expression(paste("Biomass recovery (g m"^"-1", ")", sep="")), cex.lab=1.6)
arrows(summary_allfn$SPP, summary_allfn$rslnc.mean, summary_allfn$SPP, summary_allfn$rslnc.up, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
arrows(summary_allfn$SPP, summary_allfn$rslnc.mean, summary_allfn$SPP, summary_allfn$rslnc.dn, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
#curve(predict(rslnc_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="NONE")), col="chartreuse4", add=T, lwd=4, lty=1)
#curve(predict(rslnc_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="INS")), col="goldenrod3", add=T, lwd=4, lty=1)
#curve(predict(rslnc_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="LIT")), col="chocolate4", add=T, lwd=4, lty=1)
#curve(predict(rslnc_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="BOTH")), col="lemonchiffon4", add=T, lwd=4, lty=1)
#curve(predict(rslnc_glm3, type="response", newdata=data.frame(SPP=x)), col="black", add=T, lwd=2, lty=1)
axis(side=1, at=x, labels=x)
popViewport()

dev.off()

tiff("manuscript/eco_evo_rev/figure2_jitter_2.tiff", width=8, height=16, units="in", res=300)

jit<-jitter(summary_allfn$SPP)
x<-c(1, 2, 4, 8, 16)
plot.new()
grid.newpage()
pushViewport(viewport(layout=grid.layout(4,2)))

pushViewport(viewport(layout.pos.row=1, layout.pos.col = 1))
print(p1, newpage=F)
text("a.", x=-0.1, y=1.09, xpd=NA, cex=1.3)
popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col = 1))
print(p2, newpage=F)
text("c.", x=-0.1, y=0.785, xpd=NA, cex=1.3)
popViewport()
pushViewport(viewport(layout.pos.row=3, layout.pos.col = 1))
print(p3, newpage=F)
text("e.", x=-0.1, y=0.48, xpd=NA, cex=1.3)
popViewport()
pushViewport(viewport(layout.pos.row=4, layout.pos.col = 1))
print(p4, newpage=F)
text("g.", x=-0.1, y=0.175, xpd=NA, cex=1.3)
popViewport()

pushViewport(viewport(layout.pos.row=1, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
par(mgp=c(2,1,0))
plot(summary_allfn$abv.mean~jit, col=summary_allfn$col, pch=16, xlab="No. of species", ylab="",  
     cex=1.8, ylim=c(min(summary_allfn$abv.dn), max(summary_allfn$abv.up)), cex.lab=1.6, xaxt="n")
title(ylab=expression(paste("Aboveground biomass (g m"^"-1", ")", sep="")), cex.lab=1.6)
arrows(jit, summary_allfn$abv.mean, jit, summary_allfn$abv.up, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
arrows(jit, summary_allfn$abv.mean, jit, summary_allfn$abv.dn, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
legend(-3,130, bty="n", legend=c("NONE", "ABV", "LIT", "BOTH"), col=col, pch=16, cex=1.2, xpd=NA, horiz = T)
axis(side=1, at=x, labels=x, cex.lab=1.2)
text("b.", x=17, y=130, xpd=NA, cex=1.3)
popViewport()

pushViewport(viewport(layout.pos.row=2, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
par(mgp=c(2,1,0))
plot(summary_allfn$root.mean~jit, col=summary_allfn$col, pch=16, xlab="No. of species", ylab="", 
     cex=1.8, ylim=c(min(summary_allfn$root.dn), max(summary_allfn$root.up)), cex.lab=1.6, xaxt="n")
title(ylab=expression(paste("Root biomass (g m"^"-1", ")", sep="")), cex.lab=1.6)
arrows(jit, summary_allfn$root.mean, jit, summary_allfn$root.up, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
arrows(jit, summary_allfn$root.mean, jit, summary_allfn$root.dn, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
#curve(predict(root_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="NONE")), col="chartreuse4", add=T, lwd=4, lty=1)
#curve(predict(root_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="INS")), col="goldenrod3", add=T, lwd=4, lty=1)
#curve(predict(root_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="LIT")), col="chocolate4", add=T, lwd=4, lty=1)
#curve(predict(root_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="BOTH")), col="lemonchiffon4", add=T, lwd=4, lty=1)
#curve(predict(root_glm3, type="response", newdata=data.frame(SPP=x)), col="black", add=T, lwd=2, lty=1)
axis(side=1, at=x, labels=x, cex=1.2)
text("d.", x=17, y=6.5, xpd=NA, cex=1.3)
popViewport()

pushViewport(viewport(layout.pos.row=3, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
par(mgp=c(2,1,0))
plot(summary_allfn$flo.mean~jit, col=summary_allfn$col, pch=16, xlab="No. of species", ylab="",
     cex=1.8, ylim=c(min(summary_allfn$flo.dn), max(summary_allfn$flo.up)), cex.lab=1.6, xaxt="n")
title(ylab="Water retention (s)", cex.lab=1.6)
arrows(jit, summary_allfn$flo.mean, jit, summary_allfn$flo.up, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
arrows(jit, summary_allfn$flo.mean, jit, summary_allfn$flo.dn, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
#curve(predict(flo_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="NONE")), col="chartreuse4", add=T, lwd=4, lty=1)
#curve(predict(flo_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="INS")), col="goldenrod3", add=T, lwd=4, lty=1)
#curve(predict(flo_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="LIT")), col="chocolate4", add=T, lwd=4, lty=1)
#curve(predict(flo_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="BOTH")), col="lemonchiffon4", add=T, lwd=4, lty=1)
#curve(predict(flo_glm3, type="response", newdata=data.frame(SPP=x)), col="black", add=T, lwd=2, lty=1)
axis(side=1, at=x, labels=x, cex=1)
text("f.", x=17, y=9, xpd=NA, cex=1.3)
popViewport()

pushViewport(viewport(layout.pos.row=4, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
par(mgp=c(2,1,0))
plot(summary_allfn$rslnc.mean~jit, col=summary_allfn$col, pch=16, xlab="No. of species", ylab="", 
     cex=1.8, ylim=c(min(summary_allfn$rslnc.dn), max(summary_allfn$rslnc.up)), cex.lab=1.6, xaxt="n")
title(ylab=expression(paste("Biomass recovery (g m"^"-1", ")", sep="")), cex.lab=1.6)
arrows(jit, summary_allfn$rslnc.mean, jit, summary_allfn$rslnc.up, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
arrows(jit, summary_allfn$rslnc.mean, jit, summary_allfn$rslnc.dn, lwd=1.8, angle=90, col=summary_allfn$col, length=0)
#curve(predict(rslnc_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="NONE")), col="chartreuse4", add=T, lwd=4, lty=1)
#curve(predict(rslnc_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="INS")), col="goldenrod3", add=T, lwd=4, lty=1)
#curve(predict(rslnc_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="LIT")), col="chocolate4", add=T, lwd=4, lty=1)
#curve(predict(rslnc_glm0, type="response", newdata=data.frame(SPP=x, INTERACT="BOTH")), col="lemonchiffon4", add=T, lwd=4, lty=1)
#curve(predict(rslnc_glm3, type="response", newdata=data.frame(SPP=x)), col="black", add=T, lwd=2, lty=1)
axis(side=1, at=x, labels=x)
text("h.", x=17, y=60, xpd=NA, cex=1.3)
popViewport()

dev.off()


#multifunctionality at standard thresholds (Figure 3)-------------------------------
bef$l90<-rep(0, nrow(bef))
bef$l75<-rep(0, nrow(bef))
bef$l50<-rep(0, nrow(bef))
bef$l25<-rep(0, nrow(bef))

bef$l20<-rep(0, nrow(bef))
bef$l40<-rep(0, nrow(bef))
bef$l60<-rep(0, nrow(bef))
bef$l80<-rep(0, nrow(bef))


levs<-c(90, 75, 50, 25)
for(i in 1:nrow(bef)){
  for (j in 1:length(levs)){
    lev<-levs[j]
    count<-0
    for (k in 20:23){
      if (bef[i, k]>=0.01*lev) {count<-count+1}
    }
    bef[i, j+23]<-count
  }
}
levs<-c(20, 40, 60, 80)

for(i in 1:nrow(bef)){
  for (j in 1:length(levs)){
    lev<-levs[j]
    count<-0
    for (k in 20:23){
      if (bef[i, k]>=0.01*lev) {count<-count+1}
    }
    bef[i, j+27]<-count
  }
}
colnames(bef)
head(bef)

fit.90<-glm(l90~SPP+INDEX+SPP*INDEX, bef, family="poisson")
fit.75<-glm(l75~SPP+INDEX+SPP*INDEX, bef, family="poisson")
fit.50<-glm(l50~SPP+INDEX+SPP*INDEX, bef, family="poisson")
fit.25<-glm(l25~SPP+INDEX+SPP*INDEX, bef, family="poisson")

fit.20<-glm(l20~SPP+INDEX+SPP*INDEX, bef, family="poisson")
fit.40<-glm(l40~SPP+INDEX+SPP*INDEX, bef, family="poisson")
fit.60<-glm(l60~SPP+INDEX+SPP*INDEX, bef, family="poisson")
fit.80<-glm(l80~SPP+INDEX+SPP*INDEX, bef, family="poisson")

fit.20_0<-glm(l20~SPP, bef, family="poisson")
fit.40_0<-glm(l40~SPP, bef, family="poisson")
fit.60_0<-glm(l60~SPP, bef, family="poisson")
fit.80_0<-glm(l80~SPP, bef, family="poisson")


fit.20_1<-glm(l20~SPP+INTERACT+SPP*INTERACT, bef, family="poisson")
fit.40_1<-glm(l40~SPP+INTERACT+SPP*INTERACT, bef, family="poisson")
fit.60_1<-glm(l60~SPP+INTERACT+SPP*INTERACT, bef, family="poisson")
fit.80_1<-glm(l80~SPP+INTERACT+SPP*INTERACT, bef, family="poisson")

fit.20_2<-glm(l20~SPP+INS+LIT+INS*LIT, bef, family="poisson")
fit.40_2<-glm(l40~SPP+INS+LIT+INS*LIT, bef, family="poisson")
fit.60_2<-glm(l60~SPP+INS+LIT+INS*LIT, bef, family="poisson")
fit.80_2<-glm(l80~SPP+INS+LIT+INS*LIT, bef, family="poisson")

fit.20_3<-glm(l20~SPP+INS+LIT+INS*LIT+SPP*INS+SPP*LIT+SPP*INS*LIT, bef, family="poisson")
fit.40_2<-glm(l40~SPP+INS+LIT+INS*LIT, bef, family="poisson")
fit.60_2<-glm(l60~SPP+INS+LIT+INS*LIT, bef, family="poisson")
fit.80_2<-glm(l80~SPP+INS+LIT+INS*LIT, bef, family="poisson")


library(sjPlot)
plot_model(fit.90)
plot_model(fit.75)
plot_model(fit.50)
plot_model(fit.25)

plot_model(fit.20)
plot_model(fit.40)
plot_model(fit.60)
plot_model(fit.80)

plot_model(fit.20_2)
plot_model(fit.40)
plot_model(fit.60)
plot_model(fit.80)

library(tiff)
tiff("manuscript/eco_evo_rev/thresh_4_final2_glms.tiff", width=8, height=8, units="in", res=300)
par(xpd=F)

x<-c(1,2,4,8,16)
par(mfrow=c(2,2), oma=c(2,0,0,0), mai=c(1, 1, 0.1, 0.1), xpd=F)
plot(jitter(l90)~SPP, bef, pch=16, xlab="Number of species",
    ylab="Number of functions \\nabove threshold", xaxt="n", cex.lab=1.3,
     ylim=c(-0.5,4.5), 
     col="grey")
axis(1, x, labels=x)
curve(predict(fit.90, type="response", newdata=data.frame(SPP=x, INDEX=1)), col="chartreuse4", add=T, lwd=4, lty=1)
curve(predict(fit.90, type="response", newdata=data.frame(SPP=x, INDEX=2)), col="chocolate4", add=T, lwd=4, lty=1)
curve(predict(fit.90, type="response", newdata=data.frame(SPP=x, INDEX=3)), col="goldenrod3", add=T, lwd=4, lty=1)
curve(predict(fit.90, type="response", newdata=data.frame(SPP=x, INDEX=4)), col="lemonchiffon4", add=T, lwd=4, lty=1)
legend("topright", "90%", cex=1.5, bty="n")


plot(jitter(l75)~SPP, bef, pch=16, #main="75% threshold", 
     xlab="Number of species", ylab="Number of functions \\nabove threshold", cex.lab=1.3,
     ylim=c(-0.5, 4.5), xaxt="n", 
     col="grey")
curve(predict(fit.75, type="response", newdata=data.frame(SPP=x, INDEX=1)), col="chartreuse4", add=T, lwd=4, lty=1)
curve(predict(fit.75, type="response", newdata=data.frame(SPP=x, INDEX=2)), col="chocolate4", add=T, lwd=4, lty=1)
curve(predict(fit.75, type="response", newdata=data.frame(SPP=x, INDEX=3)), col="goldenrod3", add=T, lwd=4, lty=1)
curve(predict(fit.75, type="response", newdata=data.frame(SPP=x, INDEX=4)), col="lemonchiffon4", add=T, lwd=4, lty=1)
axis(1, x, labels=x)
legend("topright", "75%", cex=1.5, bty="n")

plot(jitter(l50)~SPP, bef, pch=16, #main="50% threshold", 
     xlab="Number of species", ylab="Number of functions \\nabove threshold", cex.lab=1.3,
     ylim=c(-0.5, 4.5), xaxt="n",
     col="grey")
legend("topright", "50%", cex=1.5, bty="n")
curve(predict(fit.50, type="response", newdata=data.frame(SPP=x, INDEX=1)), col="chartreuse4", add=T, lwd=4, lty=1)
curve(predict(fit.50, type="response", newdata=data.frame(SPP=x, INDEX=2)), col="chocolate4", add=T, lwd=4, lty=1)
curve(predict(fit.50, type="response", newdata=data.frame(SPP=x, INDEX=3)), col="goldenrod3", add=T, lwd=4, lty=1)
curve(predict(fit.50, type="response", newdata=data.frame(SPP=x, INDEX=4)), col="lemonchiffon4", add=T, lwd=4, lty=1)
axis(1, x, labels=x)

plot(jitter(l25)~SPP, bef, pch=16, #main="25% threshold", 
     xlab="Number of species", ylab="Number of functions \\nabove threshold", cex.lab=1.3, ylim=c(-0.5, 4.5), xaxt="n", col="grey")
legend("topright", "25%", cex=1.5, bty="n")
curve(predict(fit.25, type="response", newdata=data.frame(SPP=x, INDEX=1)), col="chartreuse4", add=T, lwd=4, lty=1)
curve(predict(fit.25, type="response", newdata=data.frame(SPP=x, INDEX=2)), col="chocolate4", add=T, lwd=4, lty=1)
curve(predict(fit.25, type="response", newdata=data.frame(SPP=x, INDEX=3)), col="goldenrod3", add=T, lwd=4, lty=1)
curve(predict(fit.25, type="response", newdata=data.frame(SPP=x, INDEX=4)), col="lemonchiffon4", add=T, lwd=4, lty=1)
axis(1, x, labels=x)

par(xpd=NA)
legend(x=-22, y=-2.4, bty="n", legend=c("NONE", "ABV", "LIT", "BOTH"), col=col, lty=1, lwd=3, horiz = T, cex=1.3)
dev.off()

#Figure 3 alt-------------------------------
set_theme(base=theme_classic())
p1<-plot_model(fit.20, title="20%", transform="exp", show.values=T, value.offset=0.3, axis.labels = c("SPP"="Plant richness", "INDEX"="Trophic treatment", "SPP:INDEX"="Interaction"))
p2<-plot_model(fit.40, title="40%", show.values=T, value.offset=0.3, axis.labels = c("SPP"="Plant richness", "INDEX"="Trophic treatment", "SPP:INDEX"="Interaction"))
p3<-plot_model(fit.60, title="60%", show.values=T, value.offset=0.3, axis.labels = c("SPP"="Plant richness", "INDEX"="Trophic treatment", "SPP:INDEX"="Interaction"))
p4<-plot_model(fit.80, title="80%", show.values=T, value.offset=0.3, axis.labels = c("SPP"="Plant richness", "INDEX"="Trophic treatment", "SPP:INDEX"="Interaction"))

p1<-plot_model(fit.20_2, title="20%", transform="exp", show.values=T, value.offset=0.3, axis.labels = c("SPP"="Plant richness", "INS"="ABV", "INS:LIT"="BOTH"))
p2<-plot_model(fit.40_2, title="40%", show.values=T, value.offset=0.3, axis.labels = c("SPP"="Plant richness", "INS"="ABV", "INS:LIT"="BOTH"))
p3<-plot_model(fit.60_2, title="60%", show.values=T, value.offset=0.3, axis.labels = c("SPP"="Plant richness", "INS"="ABV", "INS:LIT"="BOTH"))
p4<-plot_model(fit.80_2, title="80%", show.values=T, value.offset=0.3, axis.labels = c("SPP"="Plant richness", "INS"="ABV", "INS:LIT"="BOTH"))


tiff("manuscript/eco_evo_rev/Figure3_alt_2_labs.tiff", width=8, height=16, units="in", res=300)
#par(xpd=F)

plot.new()
grid.newpage()
pushViewport(viewport(layout=grid.layout(4,2)))

pushViewport(viewport(layout.pos.row=1, layout.pos.col = 1))
print(p1, newpage=F)
text("a.", x=-0.1, y=1.09, xpd=NA, cex=1.3)
popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col = 1))
print(p2, newpage=F)
text("c.", x=-0.1, y=0.785, xpd=NA, cex=1.3)
popViewport()
pushViewport(viewport(layout.pos.row=3, layout.pos.col = 1))
print(p3, newpage=F)
text("e.", x=-0.1, y=0.48, xpd=NA, cex=1.3)
popViewport()
pushViewport(viewport(layout.pos.row=4, layout.pos.col = 1))
print(p4, newpage=F)
text("g.", x=-0.1, y=0.175, xpd=NA, cex=1.3)
popViewport()

x<-c(1, 2, 4, 8, 16)
pushViewport(viewport(layout.pos.row=1, layout.pos.col = 2))
par(fig=gridFIG(), new=T)

#par(mfcol=c(4,2), oma=c(2,0,0,0), mai=c(1, 1, 0.1, 0.1), xpd=F)

plot(jitter(l20)~SPP, bef, pch=16,
     xlab="Number of species", ylab="Multifunctionality", cex.lab=1.3,
     ylim=c(-0.5,4.5), xaxt="n",
     col="grey")
curve(predict(fit.20_2, type="response", newdata=data.frame(SPP=x, INS=0, LIT=0)), col="chartreuse4", add=T, lwd=4, lty=1)
curve(predict(fit.20_2, type="response", newdata=data.frame(SPP=x, INS=0, LIT=1)), col="chocolate4", add=T, lwd=4, lty=1)
curve(predict(fit.20_2, type="response", newdata=data.frame(SPP=x, INS=1, LIT=0)), col="goldenrod3", add=T, lwd=4, lty=1)
curve(predict(fit.20_2, type="response", newdata=data.frame(SPP=x, INS=1, LIT=1)), col="lemonchiffon4", add=T, lwd=4, lty=1)
axis(1, x, labels=x)
legend("topright", "20%", cex=1.5, bty="n")
text("b.", x=17, y=6, xpd=NA, cex=1.3)
popViewport()

pushViewport(viewport(layout.pos.row=2, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
plot(jitter(l40)~SPP, bef, pch=16, #main="75% threshold", 
     xlab="Number of species", ylab="Multifunctionality", cex.lab=1.3,
     ylim=c(-0.5, 4.5), xaxt="n",
     col="grey")
curve(predict(fit.40_2, type="response", newdata=data.frame(SPP=x, INS=0, LIT=0)), col="chartreuse4", add=T, lwd=4, lty=1)
curve(predict(fit.40_2, type="response", newdata=data.frame(SPP=x, INS=0, LIT=1)), col="chocolate4", add=T, lwd=4, lty=1)
curve(predict(fit.40_2, type="response", newdata=data.frame(SPP=x, INS=1, LIT=0)), col="goldenrod3", add=T, lwd=4, lty=1)
curve(predict(fit.40_2, type="response", newdata=data.frame(SPP=x, INS=1, LIT=1)), col="lemonchiffon4", add=T, lwd=4, lty=1)
legend("topright", "40%", cex=1.5, bty="n")
axis(1, x, labels=x)
text("d.", x=17, y=6, xpd=NA, cex=1.3)
popViewport()

pushViewport(viewport(layout.pos.row=3, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
plot(jitter(l60)~SPP, bef, pch=16, #main="50% threshold", 
     xlab="Number of species", ylab="Multifunctionality", cex.lab=1.3,
     ylim=c(-0.5, 4.5), xaxt="n",
     col="grey")
legend("topright", "60%", cex=1.5, bty="n")
curve(predict(fit.60_2, type="response", newdata=data.frame(SPP=x, INS=0, LIT=0)), col="chartreuse4", add=T, lwd=4, lty=1)
curve(predict(fit.60_2, type="response", newdata=data.frame(SPP=x, INS=0, LIT=1)), col="chocolate4", add=T, lwd=4, lty=1)
curve(predict(fit.60_2, type="response", newdata=data.frame(SPP=x, INS=1, LIT=0)), col="goldenrod3", add=T, lwd=4, lty=1)
curve(predict(fit.60_2, type="response", newdata=data.frame(SPP=x, INS=1, LIT=1)), col="lemonchiffon4", add=T, lwd=4, lty=1)
axis(1, x, labels=x)
text("f.", x=17, y=6, xpd=NA, cex=1.3)
popViewport()

pushViewport(viewport(layout.pos.row=4, layout.pos.col = 2))
par(fig=gridFIG(), new=T)
plot(jitter(l80)~SPP, bef, pch=16, #main="25% threshold", 
     xlab="Number of species", ylab="Multifunctionality", cex.lab=1.3, ylim=c(-0.5, 4.5), xaxt="n", col="grey")
legend("topright", "80%", cex=1.5, bty="n")
curve(predict(fit.80_2, type="response", newdata=data.frame(SPP=x, INS=0, LIT=0)), col="chartreuse4", add=T, lwd=4, lty=1)
curve(predict(fit.80_2, type="response", newdata=data.frame(SPP=x, INS=0, LIT=1)), col="chocolate4", add=T, lwd=4, lty=1)
curve(predict(fit.80_2, type="response", newdata=data.frame(SPP=x, INS=1, LIT=0)), col="goldenrod3", add=T, lwd=4, lty=1)
curve(predict(fit.80_2, type="response", newdata=data.frame(SPP=x, INS=1, LIT=1)), col="lemonchiffon4", add=T, lwd=4, lty=1)
axis(1, x, labels=x)
text("h.", x=17, y=6, xpd=NA, cex=1.3)
popViewport()

par(xpd=NA)
legend(x=-22, y=-2.4, bty="n", legend=c("NONE", "ABV", "LIT", "BOTH"), col=col, lty=1, lwd=3, horiz = T, cex=1.3)
dev.off()

#correlations between functions and MF------------------
colnames(bef[28:31])<-c("l20", "l40", "l60", "l80")

mfcorrs<-matrix(0, 4, 4)
colnames(mfcorrs)<-c("20%", "40%", "60%", "80%")
rownames(mfcorrs)<-c("aboveground\\nbiomass", "root\\nbiomass", "water\\nretention", "biomass\\nrecovery")

mfcorrs[1,1]<-cor(bef$biomstd, bef[,28])
mfcorrs[1,2]<-cor(bef$biomstd, bef[,29])
mfcorrs[1,3]<-cor(bef$biomstd, bef[,30])
mfcorrs[1,4]<-cor(bef$biomstd, bef[,31])

mfcorrs[2,1]<-cor(bef$rtstd, bef[,28])
mfcorrs[2,2]<-cor(bef$rtstd, bef[,29])
mfcorrs[2,3]<-cor(bef$rtstd, bef[,30])
mfcorrs[2,4]<-cor(bef$rtstd, bef[,31])

mfcorrs[3,1]<-cor(bef$flostd, bef[,28])
mfcorrs[3,2]<-cor(bef$flostd, bef[,29])
mfcorrs[3,3]<-cor(bef$flostd, bef[,30])
mfcorrs[3,4]<-cor(bef$flostd, bef[,31])

mfcorrs[4,1]<-cor(bef$recstd, bef[,28])
mfcorrs[4,2]<-cor(bef$recstd, bef[,29])
mfcorrs[4,3]<-cor(bef$recstd, bef[,30])
mfcorrs[4,4]<-cor(bef$recstd, bef[,31])

tiff("manuscript/eco_evo_rev/mult_corrs.tiff", width=4, height=4, units="in", res=300)
#par(mfrow=c(1,2))
#corrplot(fn.cor, method="ellipse", type="upper", addCoef.col = "black", tl.cex=0.7, tl.col="black")
corrplot(mfcorrs, method="ellipse", addCoef.col = "black", tl.col="black")
dev.off()

#Figure 4-------------------------------

thresh<-bef

par(mfrow=c(1,1))
x<-seq(0, 100)
plot(1, type="n", xlim=c(0, 100), ylim=c(-0.1, 0.1), xlab="Threshold", ylab="Biodiversity-multifunctionality effect (BMF)", cex.lab=1.5)
slope.store<-numeric(0)
ind<-numeric(0)
slope.ind<-numeric(0)
thresh.ind<-numeric(0)
slope.sigma<-numeric(0)
ind.sigma<-numeric(0)

for(i in 1:length(x)){
  t<-x[i]
  thresh.1<-thresh
  thresh.1$val<-rep(t, nrow(thresh.1))
  thresh.1$mult<-rep(0, nrow(thresh.1))
  for(j in 1:nrow(thresh.1)){
    count<-0
    for (k in 20:23){
      if (thresh.1[j, k]>=0.01*t) {count<-count+1}
    }
    thresh.1$mult[j]<-count
  }
  #main.lm<-lm(mult~SPP, thresh.1)
  main.lm<-glm(mult~SPP, thresh.1, family="poisson")
  slope<-main.lm$coefficients[2]
  sig<-confint(main.lm)[2,]
  slope.store<-c(slope.store, slope)
  slope.sigma<-rbind(slope.sigma, sig)
  for (l in 1:4){
    sub<-subset(thresh.1, INDEX==l)
    #ind.lm<-lm(mult~SPP, sub)
    ind.lm<-glm(mult~SPP, sub, family="poisson")
    slp<-ind.lm$coefficients[2]
    sig.ind<-confint(ind.lm)[2,]
    ind<-c(ind, l)
    thresh.ind<-c(thresh.ind, t)
    slope.ind<-c(slope.ind, slp)
    ind.sigma<-rbind(ind.sigma, sig.ind)
    
  }
}

col<-c(col, "black")
tiff("manuscript/eco_evo_rev/diversity_multifun_poisson.tiff", width=8, height=8, units="in", res=300)

par(mfrow=c(1,1), xpd=F, mai=c(1,1,0.5,0.5))
plot(x, slope.store, ylim=c(-0.06, 0.1), xlab="Threshold", ylab="Biodiversity-multifunctionality effect (BMF)",
     cex=0.75, col="black", pch=16, cex.lab=2, type="n")
polygon(c(x, rev(x)), c(slope.sigma[,2], rev(slope.sigma[,1])),border=NA, col=paste(rgb(t(col2rgb("lemonchiffon2"))/255), "AA", sep=""))

lines(lowess(x, slope.store, f=0.15), col="black", lwd=3, lty=3)
ind.all<-cbind(thresh.ind, ind, slope.ind, ind.sigma)
for(i in 1:4){
  sub<-ind.all[ind.all[,2]==i,]
  points(x, sub[,3], cex=1.2, col=col[i], pch=16)
  lines(lowess(x, sub[,3], f=0.15), lwd=4, col=col[i])
}
abline(h=0)
legend("topright", c("PLANT ONLY", "ABV", "LIT", "BOTH"), col=c(col), pch=16, cex=1.3)
dev.off()  
#Wilcoxon tests-------------------------
slope.none<-slope.ind[1:404 %% 4 ==1]
slope.lit<-slope.ind[1:404 %% 4 ==2]
slope.ins<-slope.ind[1:404 %% 4 ==3]
slope.both<-slope.ind[1:404 %% 4 ==0]

plot(x, slope.none, ylim=c(-0.09, 0.09), pch=16, col="chartreuse4")
lines(x, ind.sigma[1:404 %% 4 ==1,1], col=col[1], lwd=2)
lines(x, ind.sigma[1:404 %% 4 ==1,2], col=col[1], lwd=2)
points(x, slope.ins, pch=16, col="goldenrod3")
points(x, slope.lit, pch=16, col="chocolate4")
points(x, slope.both, pch=16, col="lemonchiffon4")
abline(h=0)

ins.none<-wilcox.test(slope.ins, slope.none, conf.int = T)
lit.none<-wilcox.test(slope.lit, slope.none, conf.int=T)
both.none<-wilcox.test(slope.both, slope.none, conf.int=T)
