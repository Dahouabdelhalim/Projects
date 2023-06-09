require(Cairo)
#Loads other font including "Arial"
require("extrafont")

#Sets working directory
setwd("/Users/ajblake/Documents/Projects/Polarized_Light/Manuscripts/Behavioral Mechanisms Manuscript/")

#Import Spectra
trans<-read.csv("dryad/aquarium glass.csv")
qwp<-read.csv("dryad/QWP Edmund Optics 88-253.csv")
kraft<-read.csv("dryad/kraft paper.csv")
fluor_tubes<-read.csv("dryad/fluorescent tubes.csv")
monitors<-read.csv("dryad/monitors RGB 255 spectra.csv")



#Restricts spectral data to 350-700 nm
trans<-trans[trans[,1]>=350 & trans[,1]<=700,]
qwp<-qwp[qwp[,1]>=350 & qwp[,1]<=700,]
kraft<-kraft[kraft[,1]>=350 & kraft[,1]<=700,]
fluor_tubes<-fluor_tubes[fluor_tubes[,1]>=350 & fluor_tubes[,1]<=700,]
monitors<-monitors[monitors[,1]>=350 & monitors[,1]<=700,]

#Calculates spectral average of the two monitors
monitors<-with(monitors, data.frame(WL, 
								R.255=rowMeans(cbind(M1_R, M2_R)),
								G.255=rowMeans(cbind(M1_G, M2_G)),
								B.255=rowMeans(cbind(M1_B, M2_B)),
								Wh.255=rowMeans(cbind(M1_Wh, M2_Wh))))

#Multiplies the transmission of the QWP and aquarium glass
trans$AQGL<-100*(trans$AQGL/100)*(qwp$QWP/100) 

#--------Modeling illumination spectra of monitor
#Creates a function that will estimate the spectra of a monitor based on the RGB pixel values. It uses a gamma function to scale the spectra of the full brightness of the individual color channels and sums them for a final spectrum. It also takes into account the transmission of aquarium glass and the quater-wave retarder film.

monitor.spec<-function(R,G,B){
	#Corrects monitor monitors for transmission through QWP and aquarium glass
	monitors[,2:4]<-monitors[,2:4]*trans$AQGL/100
	
	#The encoding gamma numbers below were estimated from monitor monitors at 
	#different pixel values. This scales the full brightness values of the monitor 
	#monitors to levels appropriate to the RGB pixel values
	ɣ<-1.90
	spec<-cbind(monitors$WL,
		(monitors$R.255*((1/(1+0.005519431))*(R/255)^ɣ+0.005519431)
		+monitors$G.255*((1/(1+0.005519431))*(G/255)^ɣ+0.005519431)
		+monitors$B.255*((1/(1+0.005519431))*(B/255)^ɣ+0.005519431)))
}


#intensity-vs-color discrimination experiment
m_int_0.64<-monitor.spec(72,82,64)
m_int_0.93<-monitor.spec(105,120,93)
m_int_1.00<-monitor.spec(112,129,100)
m_int_1.15<-monitor.spec(129,148,115)
m_red<-monitor.spec(237,93,82)
#color-removal experiment
m_rgb<-monitor.spec(112,129,100)
m_rg<-monitor.spec(112,129,0)
m_gb<-monitor.spec(0,129,100)
m_g<-monitor.spec(0,129,0)

#intensity-vs-color discrimination experiment
m_int_0.64.col<-rgb(72,82,64, maxColorValue=255)
m_int_0.93.col<-rgb(105,120,93, maxColorValue=255)
m_int_1.00.col<-rgb(112,129,100, maxColorValue=255)
m_int_1.15.col<-rgb(129,148,115, maxColorValue=255)
m_red.col<-rgb(237,93,82, maxColorValue=255)
#color-removal experiment
m_rgb.col<-rgb(112,129,100, maxColorValue=255)
m_rg.col<-rgb(112,129,0, maxColorValue=255)
m_gb.col<-rgb(0,129,100, maxColorValue=255)
m_g.col<-rgb(0,129,0, maxColorValue=255)

#------------- Graph Code Begins ------------

#Set figure width and calculate plot sizes and figure height
#All inputs in inches
f_height<-7.5
m_bottom_in<-0.25
m_bottom_out<-0.10
m_left<-0.5
m_top<-0.1
m_right<-0.1
m_inner<-0.25
p_height<-(f_height-3*m_bottom_in-m_bottom_out-m_top)/3
f_width<-7.09
p_width<-f_width-2*m_left-m_right
xlim<-c(350,700)
ylim1<-c(0, 100)
ylim2<-c(0, 4.4e12)
ylim3<-c(0, 3.3e12)
ylim4<-c(0, 7e11)
ylim5<-c(0, 1.32e12)
ylim6<-ylim4
mtext_x_line<-1.5
mtext_y_line<-2.00
mtext_cex<-1
sub_lab_cex<-1.5
leg_cex<-1
axis_cex<-1.1765
axis1_mgp<-c(3,0.33,0)
axis2_mgp<-c(3,0.66,0)
tck<-(-0.02)
y.intersp<-1.25
subp_lab_pos_x<-285
subp_lab_pos_y<-0.95


#Export figure to a pdf with the width and height from above
cairo_pdf("Figures/Fig. S1 raw.pdf", width=f_width, height=f_height, 
	family="Arial", pointsize=8)

par(	mfcol = c(3, 2), 
		omi=c(m_bottom_out,0,m_top,m_right), 
		mai=c(m_bottom_in,m_left,0,0),
		cex=1,
		cex.axis=1,
		xpd=NA)

#--------------- Plot 1 --------------
plot.new()
plot.window(ylim=ylim1, xlim=xlim, yaxs="i")
lines(trans$WL, trans$AQGL, col=hsv(0,0,0.5))
lines(kraft$WL, kraft$kraft_paper, col=hsv(0.106,0.68,0.58))
axis(1, las=1, lwd=0, tck=tck, lwd.ticks=1, mgp=axis1_mgp, 
	at=seq(xlim[1],xlim[2], 50))
axis(2, las=1, lwd=0, tck=tck, lwd.ticks=1, mgp=axis2_mgp, 
	at=seq(ylim1[1],ylim1[2], 10))
box(bty="l")
#Subpanel Label
text(subp_lab_pos_x, subp_lab_pos_y*ylim1[2], substitute(bold(A)), adj = c(0,0), cex=sub_lab_cex)
mtext("reflectance or transmission (%)", side = 2, line=mtext_y_line, cex=mtext_cex)

legend(475, 50, bty="n", cex=leg_cex,
	c("transmission of stimulus window \\n+ λ/4 retarder film",
		"reflectance of brown kraft paper"), 
	lty=c(1,1), 
	col=c(hsv(1,0,0),col=hsv(0.106,0.68,0.58)), 
	seg.len=1.25, 
	xjust=0.5, 
	yjust=0.5, 
	y.intersp=y.intersp)

#--------------- Plot 2 --------------	
plot.new()
plot.window(ylim=ylim2, xlim=xlim, yaxs="i")
lines(fluor_tubes, col=hsv(0,0,0.5), lwd=1)
axis(1, las=1, lwd=0, tck=tck, lwd.ticks=1, mgp=axis1_mgp, 
	at=seq(xlim[1],xlim[2], 50))
axis(2, las=0, lwd=0, tck=tck, lwd.ticks=1, mgp=axis2_mgp, 
	at=seq(ylim2[1],ylim2[2], 1e12))
box(bty="l")
#Subpanel Label
text(subp_lab_pos_x, subp_lab_pos_y*ylim2[2], substitute(bold(B)), adj = c(0,0), cex=sub_lab_cex)
mtext("irradiance (photons/cm2/s/nm)", side = 2, line=mtext_y_line, cex=mtext_cex)

legend(xlim[1], ylim2[2]+3, bty="n", cex=leg_cex,
	c("fluorescent lamps"), 
	lty=c(1), 
	col=hsv(0,0,0.5), 
	seg.len=1.25, 
	xjust=0, 
	yjust=1,
	y.intersp=y.intersp)

#--------------- Plot 3 --------------	
plot.new()
plot.window(ylim=ylim3, xlim=xlim, yaxs="i")
polygon(rbind(c(400,0),monitors[,c(1,2)],c(700,0)), border=NA, col=rgb(1,0,0,0.4))
polygon(rbind(c(400,0),monitors[,c(1,4)],c(700,0)), border=NA, col=rgb(0,0,1,0.4))
polygon(rbind(c(400,0),monitors[,c(1,3)],c(700,0)), border=NA, col=rgb(0,1,0,0.4))
lines(monitors$WL, monitors$Wh.255, col=hsv(0,0,0), lwd=1)
axis(1, las=1, lwd=0, tck=tck, lwd.ticks=1, mgp=axis1_mgp, 
	at=seq(xlim[1],xlim[2], 50))
axis(2, las=0, lwd=0, tck=tck, lwd.ticks=1, mgp=axis2_mgp, 
	at=seq(ylim3[1],ylim3[2], 1e12))
box(bty="l")
#Subpanel Label
text(subp_lab_pos_x, subp_lab_pos_y*ylim3[2], substitute(bold(C)), adj = c(0,0), cex=sub_lab_cex)
mtext("relative irradiance (%)", side = 2, line=mtext_y_line, cex.lab=mtext_cex)
mtext("wavelength (nm)", side = 1, line=mtext_x_line, cex=mtext_cex)


legend(xlim[1]+11, ylim3[2], bty="n", cex=leg_cex,
	c("pure white pixels", "pure blue pixels", 
			"pure green pixels", "pure red pixels"), 
	lty=c(1,NA,NA,NA), lwd=c(1,0,0,0), pch=c(NA, 15, 15, 15),
	col=c(hsv(0,0,0), rgb(0,0,1,0.4), rgb(0,1,0,0.4), rgb(1,0,0,0.4)),
	pt.cex=2.5, seg.len=1.25, xjust=0, yjust=1, y.intersp=y.intersp)

#--------------- Plot 4 --------------	
plot.new()
plot.window(ylim=ylim4, xlim=xlim, yaxs="i")
lines(m_int_0.64, lwd=1, lty=2, col=m_int_0.64.col)
lines(m_int_0.93, lwd=1, lty=4, col=m_int_0.93.col)
lines(m_int_1.00, lwd=1, lty=1, col=m_int_1.00.col)
lines(m_int_1.15, lwd=1, lty=3, col=m_int_1.15.col)
axis(1, las=1, lwd=0, tck=tck, lwd.ticks=1, mgp=axis1_mgp, 
	at=seq(xlim[1],xlim[2], 50))
axis(2, las=0, lwd=0, tck=tck, lwd.ticks=1, mgp=axis2_mgp, 
	at=seq(ylim4[1],ylim4[2], 2e11))
box(bty="l")
#Subpanel Label
text(subp_lab_pos_x, subp_lab_pos_y*ylim4[2], substitute(bold(D)), adj = c(0,0), cex=sub_lab_cex)
mtext("irradiance (photons/cm2/s/nm)", side = 2, line=mtext_y_line, cex.lab=mtext_cex)

legend(xlim[1]+11, ylim4[2], bty="n", cex=leg_cex,
	c("treatment image", "130% intensity", "100% intensity", 
			"87% intensity", "44% intensity"), 
	col=c(NA, m_int_1.15.col, m_int_1.00.col, m_int_0.93.col, m_int_0.64.col), 
	lwd=1, lty=c(NA,3,1,4,2), seg.len=1.25, 
	xjust=0, yjust=1, y.intersp=y.intersp)

#--------------- Plot 5 --------------	
plot.new()
plot.window(ylim=ylim5, xlim=xlim, yaxs="i")
lines(m_red, 			lwd=1, col=m_red.col)
axis(1, las=1, lwd=0, tck=tck, lwd.ticks=1, mgp=axis1_mgp, 
	at=seq(xlim[1],xlim[2], 50))
axis(2, las=0, lwd=0, tck=tck, lwd.ticks=1, mgp=axis2_mgp, 
	at=seq(ylim5[1],ylim5[2], 4e11))
box(bty="l")
#Subpanel Label
text(subp_lab_pos_x, subp_lab_pos_y*ylim5[2], substitute(bold(E)), adj = c(0,0), cex=sub_lab_cex)
mtext("relative irradiance (%)", side = 2, line=mtext_y_line, cex.lab=mtext_cex)

legend(xlim[1]+11, ylim5[2], bty="n", cex=leg_cex,
	c("control image", "red image"), 
	col=c(NA, m_red.col), 
	lwd=1, lty=c(NA,1), seg.len=1.25, 
	xjust=0, yjust=1, y.intersp=y.intersp)
	
#--------------- Plot 6 --------------	
plot.new()
plot.window(ylim=ylim6, xlim=xlim, yaxs="i")
polygon(rbind(c(350,0),m_rg[,c(1,2)],c(700,0)), border=NA, col=m_rg.col)
polygon(rbind(c(350,0),m_gb[,c(1,2)],c(700,0)), border=NA, col=m_gb.col)
polygon(rbind(c(350,0),m_g[,c(1,2)],c(700,0)), border=NA, col=m_g.col)
lines(m_g, lwd=1, col="white")
lines(m_rgb, lwd=1, col="black")
axis(1, las=1, lwd=0, tck=tck, lwd.ticks=1, mgp=axis1_mgp, 
	at=seq(xlim[1],xlim[2], 50))
axis(2, las=0, lwd=0, tck=tck, lwd.ticks=1, mgp=axis2_mgp, 
	at=seq(ylim6[1],ylim6[2], 2e11))
box(bty="l")
#Subpanel Label
text(subp_lab_pos_x, subp_lab_pos_y*ylim6[2], substitute(bold(F)), adj = c(0,0), cex=sub_lab_cex)
mtext("irradiance (photons/cm2/s/nm)", side = 2, line=mtext_y_line, cex.lab=mtext_cex)
mtext("wavelength (nm)", side = 1, line=mtext_x_line, cex=mtext_cex)

legend(xlim[1]+11, ylim6[2], bty="n", cex=leg_cex,
	c("R + G + B", "R + G", "G + B", "G"),
	col=c(m_rgb.col, m_rg.col, m_gb.col, m_g.col),
	pt.cex=2.5, pch=15, seg.len=1.25, xjust=0, yjust=1, y.intersp=y.intersp)

#Close Graph dev
dev.off()
