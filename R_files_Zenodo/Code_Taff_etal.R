##########################################################################################################
# Code for: Taff, Zimmer, Vitousek. Efficacy of negative feedback predicts recovery from acute challenges.
# Written by Conor Taff: cct63@cornell.edu or cct663@gmail.com Run in R v 3.3.3; 
# Updated: 27 Feb 2018 (CCT), 5 Apr 2018 (CCT)
##########################################################################################################

## Explanation of column names in the main data file
	# Band: unique metal band number for the bird
	# SampleID: unique identifier from database; not relevant for these analyses
	# Unit: Identifier for local site where bird box was located
	# Box: Unique box number identifier nested within 'unit'
	# Sex: only females are included in this analysis
	# Treatment: identifies treatment group as control, dmso, long cort, or short cort (see paper)
	# FinalValue: not relevant for this analysis
	# PreBase: baseline cort measurement from pre-treatment capture
	# PreStress: stress induced cort measurement from pre-treatment capture
	# PreDex: cort measurement after dexamethasone injection from pre-treatment capture
	# PreStressResp: PreStress - PreBase
	# PreDexResp: PreStress - PreDex
	# PostBase: baseline cort measurement from second capture (immediately post treatment)
	# Post2Base: baseline cort measurement from third capture
	# Post2Stress: stress induced cort measurement from third capture
	# Post2Dex: post dexamethasone cort measurement from third capture
	# Post2StressResp: Post2Stress-Post2Base
	# Post2DexResp: Post2Stress-Post2Dex

### Load in libraries to be used
	library(MuMIn)			# Used in AIC table construction
	library(rethinking)		# Just to get the color transparency function
	library(plotrix)		# Just to get the corner plot label function

### Load in data. Set working directory to folder that contains data first.
	d<-read.csv("Taff_etal_Data.csv")

## Set a palette of colors that is color blind friendly
	cbbPalette <- c("#000000","#0072B2","#D55E00","#009E73")
	
## Set up data objects for each of the separate response variables
	d.pb<-na.omit(d[,c("PostBase","Treatment","PreDex","PreBase","PreStressResponse")])
	d.pb2<-na.omit(d[,c("Post2Base","Treatment","PreDex","PreBase","PreStressResponse")])
	d.ps2<-na.omit(d[,c("Post2Stress","Treatment","PreBase","PreDex","PreStress")])
	d.psr2<-na.omit(d[,c("Post2StressResp","Treatment","PreBase","PreDex","PreStressResponse")])
	d.pd2<-na.omit(d[,c("Post2Dex","Treatment","PreBase","PreDex","PreStressResponse")])
	d.pdr2<-na.omit(d[,c("Post2DexResp","Treatment","PreBase","PreDexResp","PreStressResponse")])
	
	data.list<-list(d.pb,d.pb2,d.ps2,d.psr2,d.pd2,d.pdr2)
	
## Rename response variable to 'resp' for each of these data frames
	for(i in 1:6){
		colnames(data.list[[i]])[1]<-"resp"
	}
	
## Define model list for AIC fitting
	modlist.pb<-c("resp~1","resp~Treatment","resp~PreBase","resp~Treatment+PreBase","resp~Treatment+PreDex",
		"resp~Treatment+PreStressResponse","resp~Treatment*PreBase","resp~Treatment*PreDex",
		"resp~Treatment*PreStressResponse","resp~Treatment*PreDex+PreBase","resp~Treatment*PreStressResponse+PreBase")
	modlist.pb2<-c("resp~1","resp~Treatment","resp~PreBase","resp~Treatment+PreBase","resp~Treatment+PreDex",
		"resp~Treatment+PreStressResponse","resp~Treatment*PreBase","resp~Treatment*PreDex",
		"resp~Treatment*PreStressResponse","resp~Treatment*PreDex+PreBase","resp~Treatment*PreStressResponse+PreBase")
	modlist.ps2<-c("resp~1","resp~Treatment","resp~PreStress","resp~Treatment+PreBase",
		"resp~Treatment+PreDex","resp~Treatment+PreStress","resp~Treatment*PreBase","resp~Treatment*PreDex",
		"resp~Treatment*PreStress","resp~Treatment*PreDex+PreStress","resp~Treatment*PreBase+PreStress")
	modlist.psr2<-c("resp~1","resp~Treatment","resp~PreStressResponse","resp~Treatment+PreBase",
		"resp~Treatment+PreDex","resp~Treatment+PreStressResponse","resp~Treatment*PreBase",
		"resp~Treatment*PreDex","resp~Treatment*PreStressResponse","resp~Treatment*PreDex+PreStressResponse",
		"resp~Treatment*PreBase+PreStressResponse")
	modlist.pd2<-c("resp~1","resp~Treatment","resp~PreDex","resp~Treatment+PreBase","resp~Treatment+PreDex",
		"resp~Treatment+PreStressResponse","resp~Treatment*PreBase","resp~Treatment*PreDex",
		"resp~Treatment*PreStressResponse","resp~Treatment*PreBase+PreDex","resp~Treatment*PreStressResponse+PreDex")
	modlist.pdr2<-c("resp~1","resp~Treatment","resp~PreDexResp","resp~Treatment+PreBase","resp~Treatment+PreDexResp",
		"resp~Treatment+PreStressResponse","resp~Treatment*PreBase","resp~Treatment*PreDexResp",
		"resp~Treatment*PreStressResponse","resp~Treatment*PreBase+PreDexResp","resp~Treatment*PreStressResponse+PreDexResp")
		
	modlists<-list(modlist.pb,modlist.pb2,modlist.ps2,modlist.psr2,modlist.pd2,modlist.pdr2)
		
## Define response variables
	resps<-c("PostBase","Post2Base","Post2Stress","Post2StressResp","Post2Dex","Post2DexResp")
	
## Fit candidate model set to each response variable
	for(k in 1:length(data.list)){
		modlist<-modlists[k]
		mod.output<-list()
		mod.output2<-list()
			for(i in 1:length(modlist[[1]])){
				dat<-data.list[[k]]
				model<-paste("fig.mod",i,sep="_")
				mod<-lm(as.formula(modlist[[1]][i]),data=dat)
				assign(model[1],mod)
				mod.output[[i]]<-get(model)
				mod.output2[[i]]<-summary(get(model))
			}
		AICtab<-as.data.frame(modlist)
		for(j in 1:nrow(AICtab)){
			AICtab$AICc[j]<-AICc(mod.output[[j]])
			AICtab$LL[j]<-logLik(mod.output[[j]])
			AICtab$K[j]<-attr(logLik(mod.output[[j]]),"df")
		}
		minim<-min(AICtab$AICc)
		AICtab$D.AICc<-AICtab$AICc-minim
		AICtab$rellik<-exp(-.5* AICtab$D.AICc)
		sumAICs<-sum(AICtab$rellik)
		for(u in 1:nrow(AICtab)){
			AICtab$Wi[u]<-round(AICtab$rellik[u]/sumAICs,2)
		}
		AICtab<-AICtab[order(-AICtab$rellik),]
		assign(paste("AIC.",resps[k],sep=""),AICtab)
		assign(paste("mods.",resps[k],sep=""),mod.output)
		assign(paste("mods.s.",resps[k],sep=""),mod.output2)
	}	
	
## Make a list of al the AIC tables for each response variable
	colnames(AIC.PostBase)[1]<-"Models"
	colnames(AIC.Post2Base)[1]<-"Models"
	colnames(AIC.Post2Stress)[1]<-"Models"
	colnames(AIC.Post2StressResp)[1]<-"Models"
	colnames(AIC.Post2Dex)[1]<-"Models"
	colnames(AIC.Post2DexResp)[1]<-"Models"
	AIC.tabs<-rbind(AIC.PostBase,AIC.Post2Base,AIC.Post2Stress,AIC.Post2StressResp,AIC.Post2Dex,AIC.Post2DexResp)
	AIC.tabs$Response<-c(rep("PostBase",11),rep("Post2Base",11),rep("Post2Stress",11),
		rep("Post2StressResp",11),rep("Post2Dex",11),rep("Post2DexResp",11))
	write.csv(AIC.tabs,"Output_All_AIC_Tables_Alt.csv")
	
## Relable fit model lists
	formulas<-c(gsub("resp","PostBase",modlist.pb),gsub("resp","Post2Base",modlist.pb2),
		gsub("resp","Post2Stress",modlist.ps2),gsub("resp","Post2StressResp",modlist.psr2),
		gsub("resp","Post2Dex",modlist.pd2),gsub("resp","Post2DexResp",modlist.pdr2))
	
## Write fit models to a file
	Mod.outputs<-c(mods.s.PostBase,mods.s.Post2Base,mods.s.Post2Stress,
		mods.s.Post2StressResp,mods.s.Post2Dex,mods.s.Post2DexResp)
	for(i in 1:length(formulas)){
		Mod.outputs[[i]][1]<-formulas[i]
	}
	sink("Output_All_Fit_Models.txt")
		Mod.outputs
	sink()
	
## Plot interaction
	for(i in 1:nrow(d)){
		ifelse(d$Treatment[i]=="Control",d$col[i]<-cbbPalette[1],
			ifelse(d$Treatment[i]=="DMSO",d$col[i]<-cbbPalette[2],
			ifelse(d$Treatment[i]=="Long Cort",d$col[i]<-cbbPalette[3],
			d$col[i]<-cbbPalette[4])))
	}
	
	tiff("Output_Fig2.tiff",width=6.2,height=5.8,units="in",res=300)
	palette<-cbbPalette
	plot(d$PreDex,d$Post2Base,col=d$col,pch=16,xlab="Pre-treatment dexamethasone-induced corticosterone (ng/mL)",
		ylab="2nd recapture baseline corticosterone (ng/mL)",cex=1.5)
	legend(x="topleft",legend=c("Control","Vehicle","Long-Cort","Short-Cort"),col=cbbPalette,pch=16,lwd=2,bty="n")
	abline(lm(subset(d$Post2Base,d$Treatment=="Long Cort")~subset(d$PreDex,d$Treatment=="Long Cort")),
		lwd=2,col=cbbPalette[3])
	abline(lm(subset(d$Post2Base,d$Treatment=="Short Cort")~subset(d$PreDex,d$Treatment=="Short Cort")),
		lwd=2,col=cbbPalette[4])
	abline(lm(subset(d$Post2Base,d$Treatment=="Control")~subset(d$PreDex,d$Treatment=="Control")),
		lwd=2,col=cbbPalette[1])
	abline(lm(subset(d$Post2Base,d$Treatment=="DMSO")~subset(d$PreDex,d$Treatment=="DMSO")),
		lwd=2,col=cbbPalette[2])
	dev.off()
	
## Plot main effects
	## Gathering averages and SEs to plot
		p1.c.b.m<-mean(na.omit(subset(d$PostBase,d$Treatment=="Control")))
		p1.c.b.se<-sd(na.omit(subset(d$PostBase,d$Treatment=="Control")))/
			sqrt(length(na.omit(subset(d$PostBase,d$Treatment=="Control"))))
		p1.d.b.m<-mean(na.omit(subset(d$PostBase,d$Treatment=="DMSO")))
		p1.d.b.se<-sd(na.omit(subset(d$PostBase,d$Treatment=="DMSO")))/
			sqrt(length(na.omit(subset(d$PostBase,d$Treatment=="DMSO"))))
		p1.6.b.m<-mean(na.omit(subset(d$PostBase,d$Treatment=="Long Cort")))
		p1.6.b.se<-sd(na.omit(subset(d$PostBase,d$Treatment=="Long Cort")))/
			sqrt(length(na.omit(subset(d$PostBase,d$Treatment=="Long Cort"))))	
		p1.3.b.m<-mean(na.omit(subset(d$PostBase,d$Treatment=="Short Cort")))
		p1.3.b.se<-sd(na.omit(subset(d$PostBase,d$Treatment=="Short Cort")))/
			sqrt(length(na.omit(subset(d$PostBase,d$Treatment=="Short Cort"))))
			
		p2.c.b.m<-mean(na.omit(subset(d$Post2Base,d$Treatment=="Control")))
		p2.c.b.se<-sd(na.omit(subset(d$Post2Base,d$Treatment=="Control")))/
			sqrt(length(na.omit(subset(d$Post2Base,d$Treatment=="Control"))))
		p2.d.b.m<-mean(na.omit(subset(d$Post2Base,d$Treatment=="DMSO")))
		p2.d.b.se<-sd(na.omit(subset(d$Post2Base,d$Treatment=="DMSO")))/
			sqrt(length(na.omit(subset(d$Post2Base,d$Treatment=="DMSO"))))
		p2.6.b.m<-mean(na.omit(subset(d$Post2Base,d$Treatment=="Long Cort")))
		p2.6.b.se<-sd(na.omit(subset(d$Post2Base,d$Treatment=="Long Cort")))/
			sqrt(length(na.omit(subset(d$Post2Base,d$Treatment=="Long Cort"))))	
		p2.3.b.m<-mean(na.omit(subset(d$Post2Base,d$Treatment=="Short Cort")))
		p2.3.b.se<-sd(na.omit(subset(d$Post2Base,d$Treatment=="Short Cort")))/
			sqrt(length(na.omit(subset(d$Post2Base,d$Treatment=="Short Cort"))))
			
		p2.c.d.m<-mean(na.omit(subset(d$Post2Dex,d$Treatment=="Control")))
		p2.c.d.se<-sd(na.omit(subset(d$Post2Dex,d$Treatment=="Control")))/
			sqrt(length(na.omit(subset(d$Post2Dex,d$Treatment=="Control"))))
		p2.d.d.m<-mean(na.omit(subset(d$Post2Dex,d$Treatment=="DMSO")))
		p2.d.d.se<-sd(na.omit(subset(d$Post2Dex,d$Treatment=="DMSO")))/
			sqrt(length(na.omit(subset(d$Post2Dex,d$Treatment=="DMSO"))))
		p2.6.d.m<-mean(na.omit(subset(d$Post2Dex,d$Treatment=="Long Cort")))
		p2.6.d.se<-sd(na.omit(subset(d$Post2Dex,d$Treatment=="Long Cort")))/
			sqrt(length(na.omit(subset(d$Post2Dex,d$Treatment=="Long Cort"))))	
		p2.3.d.m<-mean(na.omit(subset(d$Post2Dex,d$Treatment=="Short Cort")))
		p2.3.d.se<-sd(na.omit(subset(d$Post2Dex,d$Treatment=="Short Cort")))/
			sqrt(length(na.omit(subset(d$Post2Dex,d$Treatment=="Short Cort"))))
	
	tiff("Output_Fig1a.tiff",width=5,height=5.4,units="in",res=300)
		h1<-c(0,p1.c.b.m,p1.d.b.m,p1.6.b.m,p1.3.b.m,0,p2.c.b.m,p2.d.b.m,p2.6.b.m,p2.3.b.m,0)
		h1.u<-c(p1.c.b.m+p1.c.b.se,p1.d.b.m+p1.d.b.se,p1.6.b.m+p1.6.b.se,p1.3.b.m+p1.3.b.se,
			p2.c.b.m+p2.c.b.se,p2.d.b.m+p2.d.b.se,p2.6.b.m+p2.6.b.se,p2.3.b.m+p2.3.b.se)
		h1.l<-c(p1.c.b.m-p1.c.b.se,p1.d.b.m-p1.d.b.se,p1.6.b.m-p1.6.b.se,p1.3.b.m-p1.3.b.se,
			p2.c.b.m-p2.c.b.se,p2.d.b.m-p2.d.b.se,p2.6.b.m-p2.6.b.se,p2.3.b.m-p2.3.b.se)
		xer<-c(1.5,2.5,3.5,4.5,6.5,7.5,8.5,9.5)
		barplot(h1,space=0,ylim=c(0,60),col=rep(c("white",cbbPalette[1],cbbPalette[2],cbbPalette[3],cbbPalette[4]),2),
			ylab="Baseline corticosterone (ng/mL)")
		for(i in 1:length(xer)){
			lines(rep(xer[i],2),c(h1.u[i],h1.l[i]))
		}
		legend("topleft",legend=c("Control","Vehicle","Long-Cort","Short-Cort"),bty="n",
				fill=c(cbbPalette[1],cbbPalette[2],cbbPalette[3],cbbPalette[4]),cex=0.8)
		axis(1,at=c(-1,4,8,15),labels=c("","1st Recapture","2nd Recapture",""))	
		corner.label("(A)",x=1,y=1)
	dev.off()
		
	tiff("Output_Fig1b.tiff",width=3,height=5.4,units="in",res=300)
		h2<-c(0,p2.c.d.m,p2.d.d.m,p2.6.d.m,p2.3.d.m,0)
		h2.u<-c(p2.c.d.m+p2.c.d.se,p2.d.d.m+p2.d.d.se,p2.6.d.m+p2.6.d.se,p2.3.d.m+p2.3.d.se)
		h2.l<-c(p2.c.d.m-p2.c.d.se,p2.d.d.m-p2.d.d.se,p2.6.d.m-p2.6.d.se,p2.3.d.m-p2.3.d.se)
		xer2<-c(1.5,2.5,3.5,4.5)
		barplot(h2,space=0,col=c("white",cbbPalette[1],cbbPalette[2],cbbPalette[3],cbbPalette[4]),
			ylab="Dexamethasone-induced corticosterone (ng/mL)",ylim=c(0,60))
		for(i in 1:length(xer)){
			lines(rep(xer2[i],2),c(h2.u[i],h2.l[i]))
		}
		axis(1,at=c(-1,3,10),labels=c("","2nd Recapture",""))
		corner.label("(B)",x=1,y=1)
	dev.off()

## Plot cort value versus capture by treatment
	tiff("Output_Cort_Treatment.tiff",width=8.1,height=5,units="in",res=300)
		# Set jittering of points and transparency for all points	
			jit<-0.035
			alph<-0.3
		# Subset data by treatment
			sub1<-subset(d,d$Treatment=="Control")
			sub2<-subset(d,d$Treatment=="DMSO")
			sub3<-subset(d,d$Treatment=="Long Cort")
			sub4<-subset(d,d$Treatment=="Short Cort")
		# Begin plotting
			plot(1,1,type="n",ylim=c(-5,120),xlim=c(0.5,9.5),ylab="Corticosterone (ng/mL)",xlab="Capture number",xaxt="n")
			axis(1,at=c(0,2,5,8,10),labels=c("","First (pre-treatment)","Second","Third",""))
			legend("topleft",c("Control","Vehicle","Long-Cort","Short-Cort"),bty="n",
				pch=15,col=c(cbbPalette[1],cbbPalette[2],cbbPalette[3],cbbPalette[4]),lwd=2)
		# Add in capture 1 baseline values
			points(rep(.85,nrow(sub1))+runif(nrow(sub1),-jit,jit),sub1$PreBase,pch=16,col=col.alpha(cbbPalette[1],alph))
			points(rep(.95,nrow(sub2))+runif(nrow(sub2),-jit,jit),sub2$PreBase,pch=16,col=col.alpha(cbbPalette[2],alph))
			points(rep(1.05,nrow(sub3))+runif(nrow(sub3),-jit,jit),sub3$PreBase,pch=16,col=col.alpha(cbbPalette[3],alph))
			points(rep(1.15,nrow(sub4))+runif(nrow(sub4),-jit,jit),sub4$PreBase,pch=16,col=col.alpha(cbbPalette[4],alph))
		# Add in capture 1 mean values by group
			points(.85,mean(na.omit(sub1$PreBase)),pch=15,col=cbbPalette[1],cex=1.3)
			points(.95,mean(na.omit(sub2$PreBase)),pch=15,col=cbbPalette[2],cex=1.3)
			points(1.05,mean(na.omit(sub3$PreBase)),pch=15,col=cbbPalette[3],cex=1.3)	
			points(1.15,mean(na.omit(sub4$PreBase)),pch=15,col=cbbPalette[4],cex=1.3)
		# Add in capture 1 stress values by group
			points(rep(1.85,nrow(sub1))+runif(nrow(sub1),-jit,jit),sub1$PreStress,pch=16,col=col.alpha(cbbPalette[1],alph))
			points(rep(1.95,nrow(sub2))+runif(nrow(sub2),-jit,jit),sub2$PreStress,pch=16,col=col.alpha(cbbPalette[2],alph))
			points(rep(2.05,nrow(sub3))+runif(nrow(sub3),-jit,jit),sub3$PreStress,pch=16,col=col.alpha(cbbPalette[3],alph))
			points(rep(2.15,nrow(sub4))+runif(nrow(sub4),-jit,jit),sub4$PreStress,pch=16,col=col.alpha(cbbPalette[4],alph))
		# Add in capture 1 mean stress values by group
			points(1.85,mean(na.omit(sub1$PreStress)),pch=15,col=cbbPalette[1],cex=1.3)
			points(1.95,mean(na.omit(sub2$PreStress)),pch=15,col=cbbPalette[2],cex=1.3)
			points(2.05,mean(na.omit(sub3$PreStress)),pch=15,col=cbbPalette[3],cex=1.3)	
			points(2.15,mean(na.omit(sub4$PreStress)),pch=15,col=cbbPalette[4],cex=1.3)
		# Add in capture 1 dex values by group
			points(rep(2.85,nrow(sub1))+runif(nrow(sub1),-jit,jit),sub1$PreDex,pch=16,col=col.alpha(cbbPalette[1],alph))
			points(rep(2.95,nrow(sub2))+runif(nrow(sub2),-jit,jit),sub2$PreDex,pch=16,col=col.alpha(cbbPalette[2],alph))
			points(rep(3.05,nrow(sub3))+runif(nrow(sub3),-jit,jit),sub3$PreDex,pch=16,col=col.alpha(cbbPalette[3],alph))
			points(rep(3.15,nrow(sub4))+runif(nrow(sub4),-jit,jit),sub4$PreDex,pch=16,col=col.alpha(cbbPalette[4],alph))
		# Add in capture 1 mean dex values by group
			points(2.85,mean(na.omit(sub1$PreDex)),pch=15,col=cbbPalette[1],cex=1.3)
			points(2.95,mean(na.omit(sub2$PreDex)),pch=15,col=cbbPalette[2],cex=1.3)
			points(3.05,mean(na.omit(sub3$PreDex)),pch=15,col=cbbPalette[3],cex=1.3)	
			points(3.15,mean(na.omit(sub4$PreDex)),pch=15,col=cbbPalette[4],cex=1.3)
		# Add in capture 1 dex values by group
			points(rep(4.85,nrow(sub1))+runif(nrow(sub1),-jit,jit),sub1$PostBase,pch=16,col=col.alpha(cbbPalette[1],alph))
			points(rep(4.95,nrow(sub2))+runif(nrow(sub2),-jit,jit),sub2$PostBase,pch=16,col=col.alpha(cbbPalette[2],alph))
			points(rep(5.05,nrow(sub3))+runif(nrow(sub3),-jit,jit),sub3$PostBase,pch=16,col=col.alpha(cbbPalette[3],alph))
			points(rep(5.15,nrow(sub4))+runif(nrow(sub4),-jit,jit),sub4$PostBase,pch=16,col=col.alpha(cbbPalette[4],alph))
		# Add in capture 1 mean dex values by group
			points(4.85,mean(na.omit(sub1$PostBase)),pch=15,col=cbbPalette[1],cex=1.3)
			points(4.95,mean(na.omit(sub2$PostBase)),pch=15,col=cbbPalette[2],cex=1.3)
			points(5.05,mean(na.omit(sub3$PostBase)),pch=15,col=cbbPalette[3],cex=1.3)	
			points(5.15,mean(na.omit(sub4$PostBase)),pch=15,col=cbbPalette[4],cex=1.3)
		# Add in capture 1 baseline values
			points(rep(6.85,nrow(sub1))+runif(nrow(sub1),-jit,jit),sub1$Post2Base,pch=16,col=col.alpha(cbbPalette[1],alph))
			points(rep(6.95,nrow(sub2))+runif(nrow(sub2),-jit,jit),sub2$Post2Base,pch=16,col=col.alpha(cbbPalette[2],alph))
			points(rep(7.05,nrow(sub3))+runif(nrow(sub3),-jit,jit),sub3$Post2Base,pch=16,col=col.alpha(cbbPalette[3],alph))
			points(rep(7.15,nrow(sub4))+runif(nrow(sub4),-jit,jit),sub4$Post2Base,pch=16,col=col.alpha(cbbPalette[4],alph))
		# Add in capture 1 mean values by group
			points(6.85,mean(na.omit(sub1$Post2Base)),pch=15,col=cbbPalette[1],cex=1.3)
			points(6.95,mean(na.omit(sub2$Post2Base)),pch=15,col=cbbPalette[2],cex=1.3)
			points(7.05,mean(na.omit(sub3$Post2Base)),pch=15,col=cbbPalette[3],cex=1.3)	
			points(7.15,mean(na.omit(sub4$Post2Base)),pch=15,col=cbbPalette[4],cex=1.3)
		# Add in capture 1 stress values by group
			points(rep(7.85,nrow(sub1))+runif(nrow(sub1),-jit,jit),sub1$Post2Stress,pch=16,col=col.alpha(cbbPalette[1],alph))
			points(rep(7.95,nrow(sub2))+runif(nrow(sub2),-jit,jit),sub2$Post2Stress,pch=16,col=col.alpha(cbbPalette[2],alph))
			points(rep(8.05,nrow(sub3))+runif(nrow(sub3),-jit,jit),sub3$Post2Stress,pch=16,col=col.alpha(cbbPalette[3],alph))
			points(rep(8.15,nrow(sub4))+runif(nrow(sub4),-jit,jit),sub4$Post2Stress,pch=16,col=col.alpha(cbbPalette[4],alph))
		# Add in capture 1 mean stress values by group
			points(7.85,mean(na.omit(sub1$Post2Stress)),pch=15,col=cbbPalette[1],cex=1.3)
			points(7.95,mean(na.omit(sub2$Post2Stress)),pch=15,col=cbbPalette[2],cex=1.3)
			points(8.05,mean(na.omit(sub3$Post2Stress)),pch=15,col=cbbPalette[3],cex=1.3)	
			points(8.15,mean(na.omit(sub4$Post2Stress)),pch=15,col=cbbPalette[4],cex=1.3)
		# Add in capture 1 dex values by group
			points(rep(8.85,nrow(sub1))+runif(nrow(sub1),-jit,jit),sub1$Post2Dex,pch=16,col=col.alpha(cbbPalette[1],alph))
			points(rep(8.95,nrow(sub2))+runif(nrow(sub2),-jit,jit),sub2$Post2Dex,pch=16,col=col.alpha(cbbPalette[2],alph))
			points(rep(9.05,nrow(sub3))+runif(nrow(sub3),-jit,jit),sub3$Post2Dex,pch=16,col=col.alpha(cbbPalette[3],alph))
			points(rep(9.15,nrow(sub4))+runif(nrow(sub4),-jit,jit),sub4$Post2Dex,pch=16,col=col.alpha(cbbPalette[4],alph))
		# Add in capture 1 mean dex values by group
			points(8.85,mean(na.omit(sub1$Post2Dex)),pch=15,col=cbbPalette[1],cex=1.3)
			points(8.95,mean(na.omit(sub2$Post2Dex)),pch=15,col=cbbPalette[2],cex=1.3)
			points(9.05,mean(na.omit(sub3$Post2Dex)),pch=15,col=cbbPalette[3],cex=1.3)	
			points(9.15,mean(na.omit(sub4$Post2Dex)),pch=15,col=cbbPalette[4],cex=1.3)
		## Add text labels
			text(1,-5,"Base",cex=.8)
			text(2,-5,"Stress",cex=.8)
			text(3,-5,"Dex",cex=.8)
			text(5,-5,"Base",cex=.8)
			text(7,-5,"Base",cex=.8)
			text(8,-5,"Stress",cex=.8)
			text(9,-5,"Dex",cex=.8)
		## Add lines		
			lines(c(.85,1.85,2.85),c(mean(na.omit(sub1$PreBase)),mean(na.omit(sub1$PreStress)),
				mean(na.omit(sub1$PreDex))),col=cbbPalette[1],lwd=2)
			lines(c(.95,1.95,2.95),c(mean(na.omit(sub2$PreBase)),mean(na.omit(sub2$PreStress)),
				mean(na.omit(sub2$PreDex))),col=cbbPalette[2],lwd=2)
			lines(c(1.05,2.05,3.05),c(mean(na.omit(sub3$PreBase)),mean(na.omit(sub3$PreStress)),
				mean(na.omit(sub3$PreDex))),col=cbbPalette[3],lwd=2)
			lines(c(1.15,2.15,3.15),c(mean(na.omit(sub4$PreBase)),mean(na.omit(sub4$PreStress)),
				mean(na.omit(sub4$PreDex))),col=cbbPalette[4],lwd=2)
			lines(c(6.85,7.85,8.85),c(mean(na.omit(sub1$Post2Base)),mean(na.omit(sub1$Post2Stress)),
				mean(na.omit(sub1$Post2Dex))),col=cbbPalette[1],lwd=2)
			lines(c(6.95,7.95,8.95),c(mean(na.omit(sub2$Post2Base)),mean(na.omit(sub2$Post2Stress)),
				mean(na.omit(sub2$Post2Dex))),col=cbbPalette[2],lwd=2)
			lines(c(7.05,8.05,9.05),c(mean(na.omit(sub3$Post2Base)),mean(na.omit(sub3$Post2Stress)),
				mean(na.omit(sub3$Post2Dex))),col=cbbPalette[3],lwd=2)
			lines(c(7.15,8.15,9.15),c(mean(na.omit(sub4$Post2Base)),mean(na.omit(sub4$Post2Stress)),
				mean(na.omit(sub4$Post2Dex))),col=cbbPalette[4],lwd=2)
	dev.off()
	
## Test for differences between treatment groups in pre-treatment characteristics
	anova(lm(d$PreStress~d$Treatment))
	anova(lm(d$PreStressResponse~d$Treatment))
	anova(lm(d$PreDex~d$Treatment))	
	
	

