require(data.table)
require(ggplot2)
library(dplyr)
library(grid)
library(gtable)
library(multcomp)
library(lsmeans)
require(Rcpp)
require(scales)
library(boot)


dr = read.table("/Users/remi/Documents/biology/Research/PlasticityRobustness/Simulations/archive/DR_archive.txt", header=TRUE)
 
dr$ids = paste0(dr$Project, dr$Replicate)

drTO1 = dr[seq(1,nrow(dr), 2),]
drTO2 = dr[seq(2,nrow(dr), 2),]
if (any(dim(drTO1) != dim(drTO2)))
{
	stop("error")
}
drTO1$meanP_otherEnv = drTO2$meanP
C_pos = grep("C",drTO1$Project)
S_pos = grep("S",drTO1$Project)
T_pos = grep("T",drTO1$Project)
if (length(C_pos) + length(S_pos) + length(T_pos) != nrow(drTO1))
{
	stop("error")
}
drTO1$meanW[C_pos] = drTO1$meanW[C_pos]
drTO1$meanW[S_pos] = (drTO1$meanW[S_pos] + drTO2$meanW[S_pos])/2
drTO1$meanW[T_pos] = sqrt(drTO1$meanW[T_pos] * drTO2$meanW[T_pos])
dr = drTO1
drTO1=NULL
drTO2=NULL

nbGenerationsSampled = length(unique(dr$Generation))
dr$SimID = paste0(dr$Project, "_", dr$Replicate)
dr$RNslope = (dr$meanP_otherEnv - dr$meanP) / 2000
dr$Plastic = dr$RNslope > 0.2



########################
###### Prepare DR ######
########################

dr$PlasticCategory__End = ""
dr$PlasticCategory_50k_ = ""
dr$PlasticCategory_5k_ = ""
dr$PlasticCategory_4k_ = ""
dr$PlasticCategory_start_ = ""
dr$PlasticCategory_current_End = ""
dr$PlasticCategory_start_End = ""
dr$PlasticCategory_50k_End = ""
dr$PlasticCategory_5k_End = ""
dr$PlasticCategory_4k_End = ""

dr = subset(dr, Generation <= 200000)


for (i in 1:nrow(dr))
{
	SimIDMatchRows = dr$SimID==dr$SimID[i]
	if (length(which(SimIDMatchRows)) > nbGenerationsSampled)
	{
		print(paste("length(which(SimIDMatchRows)) = ", length(which(SimIDMatchRows))))
		stop()
	}
	generations 								= dr$Generation[SimIDMatchRows]
	areOtherGenerationsPlastic 	= dr$Plastic[SimIDMatchRows]
	

	if (all(generations <= dr$Generation[i])) # is last generation
	{
		if (dr$Plastic[i])
		{
			dr$PlasticCategory__End[SimIDMatchRows] = "P"	
		} else
		{
			dr$PlasticCategory__End[SimIDMatchRows] = "NP"
		}
	} else if (dr$Generation[i] == 2)
	{
		if (dr$Plastic[i])
		{
			dr$PlasticCategory_start_[SimIDMatchRows] = "P"
		} else
		{
			dr$PlasticCategory_start_[SimIDMatchRows] = "NP"
		}
	}	else if (dr$Generation[i] == 4000)
	{
		if (dr$Plastic[i])
		{
			dr$PlasticCategory_4k_[SimIDMatchRows] = "P"
		} else
		{
			dr$PlasticCategory_4k_[SimIDMatchRows] = "NP"
		}
	} else if (dr$Generation[i] == 5000)
	{
		if (dr$Plastic[i])
		{
			dr$PlasticCategory_5k_[SimIDMatchRows] = "P"
		} else
		{
			dr$PlasticCategory_5k_[SimIDMatchRows] = "NP"
		}
	} else if (dr$Generation[i] == 50000)
	{
		if (dr$Plastic[i])
		{
			dr$PlasticCategory_50k_[SimIDMatchRows] = "P"
		} else
		{
			dr$PlasticCategory_50k_[SimIDMatchRows] = "NP"
		}
	}
}

if (any(dr$PlasticCategory__End=="")){stop("error")}
if (any(dr$PlasticCategory_start_=="")){stop("error")}
if (any(dr$PlasticCategory_50k_=="")){stop("error")}
if (any(dr$PlasticCategory_5k_=="")){stop("error")}
if (any(dr$PlasticCategory_4k_=="")){stop("error")}

dr$PlasticCategory_current_End = paste0(ifelse(dr$Plastic, "P", "NP"), " - ", dr$PlasticCategory__End)
dr$PlasticCategory_start_End = paste0(dr$PlasticCategory_start_, " - ", dr$PlasticCategory__End)
dr$PlasticCategory_50k_End = paste0(dr$PlasticCategory_50k_, " - ", dr$PlasticCategory__End)
dr$PlasticCategory_5k_End = paste0(dr$PlasticCategory_5k_, " - ", dr$PlasticCategory__End)
dr$PlasticCategory_4k_End = paste0(dr$PlasticCategory_4k_, " - ", dr$PlasticCategory__End)



dr$environment = ""
dr$environment[grepl("S",dr$Project)] = "Spatial heterogeneity"
dr$environment[grepl("T",dr$Project)] = "Temporal heterogeneity"
dr$environment[grepl("C.*F",dr$Project)] = "Constant low environment"
dr$environment[grepl("C.*G",dr$Project)] = "Constant high environment"
dr$environment[grepl("C.*K",dr$Project)] = "Constant environment (inv)"
dr$environment = as.factor(dr$environment)
dr$environment = factor(dr$environment, levels = unique(dr$environment))

dr$signal = ""
dr$signal[grepl("no",dr$Project)] = "No signal"
dr$signal[grepl("so",dr$Project)] = "Env. Signal"
dr$signal[grepl("np",dr$Project)] = "Perf. feedback"
dr$signal[grepl("sp",dr$Project)] = "Both signals"
dr$signal = as.factor(dr$signal)
dr$signal = factor(dr$signal, levels = unique(dr$signal)[c(1,3,2,4)])


table(dr$signal, dr$environment, dr$Generation)



################################
###### Graphing functions ######
################################

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

rmXAxis = theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
rmYAxis = theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

dr = subset(dr, signal != "Both signals")
dr$isLastGeneration = dr$Generation == 100000
dr$isLastGeneration[dr$environment == "Spatial heterogeneity" & dr$Replicate <= 100] = dr$Generation[dr$environment == "Spatial heterogeneity" & dr$Replicate <= 100] == 200000

if (sum(dr$isLastGeneration) != length(unique(dr$ids)))
{
	stop("sum(dr$isLastGeneration) != length(unique(dr$ids))")
}


f.w = function(pheno, env) {exp(-(pheno - env)^2/(5e5))}
f.wLowEnv = function(pheno){f.w(pheno,1000)}
f.wHighEnv = function(pheno){f.w(pheno,3000)}

# function f.wLowEnv defined below
f.wLowEnv(qnorm(p=0.25,mean=1000,sd=250)) # 0.9447197
f.wLowEnv(3000)  # 0.0003354626
f.wHighEnv(1000) # 0.0003354626

# In Mathematica: Mean[TransformedDistribution[  Exp[-((z - 1000)^2/500000)], {z \\[Distributed] NormalDistribution[1000, 250]}]] // N



####################################
## Some means and standard errors ##
####################################


drSub = subset(dr, Generation == 100000 & !Plastic & environment!="Constant low environment")
as.data.table(drSub)[,list(mean=mean(sdP), SE=sd(sdP)/sqrt(.N)),list(environment)]

drSub = subset(dr, Generation == 100000 & Plastic & environment!="Constant low environment")
as.data.table(drSub)[,list(mean=mean(sdP), SE=sd(sdP)/sqrt(.N)),list(environment)]


################################
###### DR over treatments ######
################################

drAverage = as.data.frame(as.data.table(dr)[,
	{
		meanmeanW = mean(meanW)
		SEmeanW = sd(meanW) / sqrt(.N)
		meansdP = mean(sdP)
		SEsdP = sd(sdP) / sqrt(.N)
		list(
			meanW 	= meanmeanW,
			meanW_lowerCI = meanmeanW - 2 * SEmeanW,
			meanW_upperCI = meanmeanW + 2 * SEmeanW,
			SEmeanW = SEmeanW,

			sdP 	= meansdP,
			sdP_lowerCI = meansdP - 2 * SEsdP,
			sdP_upperCI = meansdP + 2 * SEsdP,
			sdP_lowerBar = meansdP - SEsdP,
			sdP_upperBar = meansdP + SEsdP,
			SEsdP = SEsdP,
			nbObservations = .N
		)
	},
	list(
		Generation,
		PlasticCategory__End,
		Project,
		environment,
		signal
	)
])

drAverage$Project = factor(drAverage$Project, levels=unique(drAverage$Project)[c(1,3,2,4,6,5,7,8,10,9,11)])

pos = position_dodge(0.1)

d = drAverage
#d = subset(d, grepl("NP - ",PlasticCategory_start_End))

dM = subset(d, signal != "Both signals" & ((environment == "Spatial heterogeneity" & Generation == 200000) | (environment == "Constant low environment"  & Generation == 100000)))
dM <- dM %>%
  group_by(Project) %>%
  mutate(
    width = 0.03 * n()
  )

ggplot(dM, aes(x=0, y=sdP, colour=PlasticCategory__End, width=width)) + geom_point(position = pos, size=3.3) + geom_errorbar(aes(ymax = sdP_lowerBar, ymin = sdP_upperBar),position = pos, size=0.7) + ggtitle("") + facet_grid(.~environment+signal) + theme_classic(12) + rmXAxis + ylab("")+ theme(legend.title=element_blank()) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) + coord_cartesian(ylim = c(210, 430)) + scale_y_continuous(breaks=c(250,300,350, 400)) + theme(axis.text.y = element_text(size=15))



##########################################################
###### Reaction norms in main treatment of interest ######
##########################################################


d = subset(dr, isLastGeneration & (environment == "Constant low environment" | environment == "Constant environment (inv)")) #| Generation==100000)
d = data.frame(
	pheno  = c(d$meanP, d$meanP_otherEnv),
	optima = c(rep(1000,nrow(d)), rep(3000,nrow(d))),
	Project = rep(d$Project,2),
	ID = rep(1:nrow(d), 2),
	trueID = rep(d$ids,2),
	Plastic = rep(abs(d$RNslope) > 0.2,2),
	Generation = rep(d$Generation, 2),
	environment = rep(d$environment, 2),
	signal = rep(d$signal, 2),
	RNslope = rep(d$RNslope, 2)
)
d$optima[d$optima==1000 & d$environment == "Constant environment (inv)"] = -1
d$optima[d$optima==3000 & d$environment == "Constant environment (inv)"] = 1000
d$optima[d$optima==-1 & d$environment == "Constant environment (inv)"] = 3000

d$Plastic[d$RNslope > 0.17 & d$environment == "Constant environment (inv)"] = TRUE


mainplot = ggplot(d, aes(x=optima, y=pheno, group=ID, color=Plastic))+ geom_point() + geom_line() + facet_grid(.~environment+signal, scales="free") + scale_y_continuous(name="Average phenotype",breaks=c(1000,3000), labels=c("zopt0","zopt1"), limits=c(500,3500)) + scale_x_continuous(name="Developmental Environment",breaks=c(1000,3000), labels=c("E0","E1")) + theme_classic(15) + scale_colour_manual(name="", values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) + theme(legend.position = c(0, 1.15), legend.justification = c(0, 1.1), legend.background=element_blank()) + theme(plot.margin = unit(c(1,0,1,1), "cm")) + geom_hline(yintercept = c(1000,3000), color="black", linetype="dashed", size=1)

sideplot = ggplot(data.frame(pheno = 0), aes(x=pheno)) + stat_function(fun = f.wLowEnv, color="black", linetype="dashed",) + stat_function(fun = f.wHighEnv, color="black", linetype="dashed",) + theme_classic(15) + scale_x_continuous(limits=c(500,3500)) + scale_y_continuous(name="Fitness", breaks=c(0,0.5,1), labels=c(0,0.5,1), limits=c(0,1)) + rmYAxis + coord_flip() + theme(plot.margin = unit(c(2.65,1,1.08,0), "cm")) + geom_vline(xintercept = c(1000,3000), color="black", linetype="dashed", size=1) + theme(axis.line.y=element_blank())

multiplot(mainplot, sideplot, layout=matrix(c(1,1,1,1,2),nrow=1))



###################################################
###### appendix B. Constant high environment ######
###################################################


d = subset(dr, Generation == 5e3 & ((environment == "Constant low environment" | environment == "Constant high environment") & signal == "Perf. feedback")) #| )
stopifnot(  abs((d$meanP_otherEnv - d$meanP) - (d$RNslope * 2000)) < 0.0001 )
plot(y = d$meanP_otherEnv - d$meanP, x = d$RNslope * 2000)
d = data.frame(
	pheno  = c(d$meanP, d$meanP_otherEnv),
	optima = c(rep(1000,nrow(d)), rep(3000,nrow(d))),
	Project = rep(d$Project,2),
	ID = rep(1:nrow(d), 2),
	trueID = rep(d$ids,2),
	Plastic = rep(abs(d$RNslope) > 0.2,2),
	Generation = rep(d$Generation, 2),
	environment = rep(d$environment, 2),
	signal = rep(d$signal, 2),
	RNslope = rep(d$RNslope, 2)
)
stopifnot(  abs((d$meanP_otherEnv - d$meanP) - (d$RNslope * 2000)) < 0.0001 )

stopifnot(abs(subset(d, Plastic)$RNslope) > 0.2)
stopifnot(abs(subset(d, !Plastic)$RNslope) < 0.2)




ggplot(d, aes(x=optima, y=pheno, group=ID, color=Plastic))+ geom_point() + geom_line() + facet_grid(.~environment+signal, scales="free") + scale_y_continuous(name="Average phenotype",breaks=c(1000,4000), labels=c(expression('Z'['opt,0']), expression('Z'['opt,1'])), limits=c(500,4500)) + scale_x_continuous(name="Developmental Environment",breaks=c(1000,3000), labels=c(expression('E'[0]),expression('E'[1]))) + theme_classic(15) + scale_colour_manual(name="", values=c("grey", "black")) + theme(plot.margin = unit(c(1,1,1,1), "cm")) + geom_hline(yintercept = c(1000,3000), color="black", linetype="dashed", size=1)


###########################################
###### Fraction of plastic genotypes ######
###########################################


mean((subset(dr, Generation == 5e3 & (environment == "Constant low environment" & signal == "Perf. feedback"))$RNslope) > 0.2)
mean((subset(dr, Generation == 5e3 & (environment == "Constant high environment" & signal == "Perf. feedback"))$RNslope) > 0.2)
tmp = subset(dr, Generation == 5e3 & (environment == "Constant high environment" & signal == "Perf. feedback"))
sum(tmp$RNslope > 0.2)
sum(tmp$RNslope < -0.2)


##########################################
###### Mean and SE of reaction norm ######
##########################################


as.data.table(subset(dr, isLastGeneration & signal == "Env. Signal" & (environment == "Constant low environment" | environment == "Constant environment (inv)")))[,
	list(
		fractionPP = mean(Plastic),
		fractionOfPPNegSlope=mean(RNslope[abs(RNslope) > 0.2] < 0),
		fractionOfPPPosSlope=mean(RNslope[abs(RNslope) > 0.2] > 0),
		fractionOfPPNegSlope=sum(RNslope[abs(RNslope) > 0.2] < 0),
		fractionOfPPPosSlope=sum(RNslope[abs(RNslope) > 0.2] > 0),
		N=.N,
		meanRNslope = mean(RNslope),
		SERNslope = sd(RNslope) / .N
	),
	list(environment, signal)
]

t.test(RNslope ~ environment, data=subset(dr, isLastGeneration & signal == "Env. Signal" & (environment == "Constant low environment" | environment == "Constant environment (inv)")))



############################################################################
###### Developmental robustness - every single simulation represented ######
############################################################################


drAverage = as.data.frame(as.data.table(subset(dr, isLastGeneration))[,
	{
		meanmeanW = mean(meanW)
		SEmeanW = sd(meanW) / sqrt(.N)
		meansdP = mean(sdP)
		SEsdP = sd(sdP) / sqrt(.N)
		list(
			meanW 	= meanmeanW,
			meanW_lowerCI = meanmeanW - 1.96 * SEmeanW,
			meanW_upperCI = meanmeanW + 1.96 * SEmeanW,
			SEmeanW = SEmeanW,

			sdP 	= meansdP,
			sdP_lowerCI = meansdP - 1.96 * SEsdP,
			sdP_upperCI = meansdP + 1.96 * SEsdP,
			sdP_lowerBar = meansdP - SEsdP,
			sdP_upperBar = meansdP + SEsdP,
			SEsdP = SEsdP,
			nbObservations = .N
		)
	},
	list(
		Generation,
		isLastGeneration,
		PlasticCategory__End,
		#PlasticCategory_start_End,
		Project,
		environment,
		signal
	)
])

drAverage$Project = factor(drAverage$Project, levels=unique(drAverage$Project)[c(1,3,2,4,6,5,7,8,10,9,11)])

d = drAverage
#d = subset(d, grepl("NP - ",PlasticCategory_start_End))

dM = subset(d, signal != "Both signals" & ((environment == "Constant low environment" & isLastGeneration) | (environment == "Spatial heterogeneity" & Generation == 200000)))
dM <- dM %>%
  group_by(Project) %>%
  mutate(
    width = 0.16 * n()
  )

ggplot(subset(dr, signal != "Both signals" & ((environment == "Constant low environment" & isLastGeneration) | (environment == "Spatial heterogeneity" & Generation == 200000))), aes(x=0, y=sdP, colour=PlasticCategory__End)) + geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width=0.5 ), size=0.5, alpha=0.5) + ggtitle("") + facet_grid(.~environment+signal) + theme_classic(12) + rmXAxis + ylab("")+ theme(legend.title=element_blank()) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) + geom_point(data=dM, position = position_dodge(0.5), size=3.3) + geom_errorbar(data=dM, aes(ymax = sdP_lowerBar, ymin = sdP_upperBar, width=width),position = position_dodge(0.5), size=1)


###################################
###### Same but with boxplot ######
###################################


dip = subset(dr, signal != "Both signals" & ((environment == "Constant low environment" & isLastGeneration) | (environment == "Spatial heterogeneity" & Generation == 200000)))

dip2 = subset(dip, )

anova(lm(dip$nbGenes ~ dip$PlasticCategory__End + dip$signal + dip$environment))

as.data.table(dip)[,list(mean = mean(nbGenes), var= var(nbGenes), max = max(nbGenes), min = min(nbGenes)),list(PlasticCategory__End, environment, signal)]

ggplot(dip, aes(x=0, y=nbGenes, colour=PlasticCategory__End)) + geom_boxplot() + ggtitle("") + facet_grid(.~environment+signal) + theme_classic(12) + rmXAxis + ylab("")+ theme(legend.title=element_blank()) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic"))





#######################################################################
###### Developmental robustness - main treatments and stat tests ######
#######################################################################


drAverage = as.data.frame(as.data.table(subset(dr, isLastGeneration))[,
	{
		meanmeanW = mean(meanW)
		SEmeanW = sd(meanW) / sqrt(.N)
		meansdP = mean(sdP)
		SEsdP = sd(sdP) / sqrt(.N)
		list(
			meanW 	= meanmeanW,
			meanW_lowerCI = meanmeanW - 1.96 * SEmeanW,
			meanW_upperCI = meanmeanW + 1.96 * SEmeanW,
			SEmeanW = SEmeanW,

			sdP 	= meansdP,
			sdP_lowerCI = meansdP - 1.96 * SEsdP,
			sdP_upperCI = meansdP + 1.96 * SEsdP,
			sdP_lowerBar = meansdP - SEsdP,
			sdP_upperBar = meansdP + SEsdP,
			SEsdP = SEsdP,
			nbObservations = .N
		)
	},
	list(
		Generation,
		PlasticCategory__End,
		#PlasticCategory_start_End,
		Project,
		environment,
		signal
	)
])

drAverage$Project = factor(drAverage$Project, levels=unique(drAverage$Project)[c(1,3,2,4,6,5,7,8,10,9,11)])

pos = position_dodge(0.5)

d = drAverage
#d = subset(d, grepl("NP - ",PlasticCategory_start_End))

dM = subset(d, signal != "Both signals" & signal != "No signal" & environment == "Spatial heterogeneity" & Generation == 200000)
dM <- dM %>%
  group_by(Project) %>%
  mutate(
    width = 0.08 * n()
  )

ggplot(dM, aes(x=0, y=sdP, colour=PlasticCategory__End, width=width)) + geom_point(position = pos, size=3.3) + geom_errorbar(aes(ymax = sdP_lowerBar, ymin = sdP_upperBar),position = pos, size=0.7) + ggtitle("") + facet_grid(.~environment+signal) + theme_classic(12) + rmXAxis + ylab("")+ theme(legend.title=element_blank()) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic"))



ggplot(subset(dr, isLastGeneration & signal != "Both signals" & signal != "No signal" & environment == "Spatial heterogeneity" & Generation == 200000), aes(x=0, y=sdP, colour=PlasticCategory__End)) + geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width=0.5 ), size=0.5) + ggtitle("") + facet_grid(.~environment+signal) + theme_classic(12) + rmXAxis + ylab("")+ theme(legend.title=element_blank()) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) + geom_point(data=dM, position = pos, size=3.3) + geom_errorbar(data=dM, aes(ymax = sdP_lowerBar, ymin = sdP_upperBar, width=width),position = pos, size=0.7)


ggplot(dM, aes(x=0, y=meanW, colour=PlasticCategory__End, width=width)) + geom_point(position = pos, size=3.3) + geom_errorbar(aes(ymax = meanW_upperCI, ymin = meanW_lowerCI),position = pos, size=0.7) + ggtitle("") + facet_grid(.~environment+signal) + theme_classic(12) + rmXAxis + ylab("")+ theme(legend.title=element_blank()) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic"))




dtest = subset(dr, grepl("NP - ",PlasticCategory_start_End) & signal != "Both signals" & signal != "No signal" & environment == "Spatial heterogeneity" & Generation == 200000)
dtestA= subset(dtest, signal == "Env. Signal" )
t.test(formula=sdP~PlasticCategory_start_End, data=dtestA)
wilcox.test(formula=sdP~PlasticCategory_start_End, data=dtestA)

dtestB= subset(dtest, signal == "Perf. feedback" )
t.test(formula=sdP~PlasticCategory_start_End, data=dtestB)
wilcox.test(formula=sdP~PlasticCategory_start_End, data=dtestB)

dtestC= subset(dtest, PlasticCategory_start_End == "NP - P")
t.test(formula=sdP~signal, data=dtestC)
wilcox.test(formula=sdP~signal, data=dtestC)



#################################################################
###### Developmental robustness - all treatments and stats ######
#################################################################


dS = subset(d, environment == "Spatial heterogeneity" & Generation == 200000 | environment != "Spatial heterogeneity" & Generation == 100000)
dS = subset(dS, environment != "Constant high environment")
dS <- dS %>%
  group_by(Project) %>%
  mutate(
    width = 0.08 * n()
  )  

ggplot(dS, aes(x=0, y=sdP, colour=PlasticCategory__End, width=width*1.5)) + geom_point(position = pos, size=3.8) + geom_errorbar(aes(ymax = sdP_lowerBar, ymin = sdP_upperBar),position = pos, size=1) + ggtitle("") + facet_grid(.~environment+signal) + theme_classic(12) + theme(axis.text.y = element_text(size=16)) + rmXAxis + ylab("")+ theme(legend.title=element_blank()) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) + scale_y_continuous(breaks=c(300,550,800))



dd = subset(dS, (PlasticCategory_start_End == "NP - P" | PlasticCategory_start_End == "NP - NP") & environment == "Constant low environment" & signal == "Perf. feedback")
dd = subset(dr, Generation == 2 & (PlasticCategory_start_End == "NP - P" | PlasticCategory_start_End == "NP - NP") & environment == "Constant low environment" & signal == "Perf. feedback")
t.test(sdP~PlasticCategory_start_End, data=dd)



dd = subset(dr, isLastGeneration & PlasticCategory__End=="P" & environment == "Constant low environment" & signal == "Perf. feedback")
anova(lm(sdP~RNslope,data=dd))
plot(y=dd$sdP, x=dd$RNslope)


###############################################################
###### Differences in developmental robustness over time ######
###############################################################


dd = subset(dr, isLastGeneration & (PlasticCategory_start_End == "NP - NP" | PlasticCategory_start_End == "NP - P") & environment == "Constant low environment" & signal != "No signal")


diffs = c()
dd = subset(dr, PlasticCategory_start_End == "NP - P" & environment == "Constant low environment" & signal != "No signal")
for (id in unique(dd$ids))
{
	ddd = subset(dd, ids==id)
	previous_sdP = NA
	for (gen in sort(unique(ddd$Generation)))
	{
		dddd = subset(ddd, Generation == gen)
		stopifnot(nrow(dddd) == 1)
		current_sdP = dddd$sdP[1]
		if (dddd$Plastic[1])
		{
			diffs = c(diffs, current_sdP - previous_sdP)
			break
		}
		previous_sdP = current_sdP
	}
}


stopifnot(all(subset(dr, isLastGeneration & PlasticCategory_start_End == "P - P")$Plastic))
stopifnot(all(subset(dr, isLastGeneration & PlasticCategory_start_End == "NP - P")$Plastic))
stopifnot(!any(subset(dr, isLastGeneration & PlasticCategory_start_End == "P - NP")$Plastic))
stopifnot(!any(subset(dr, isLastGeneration & PlasticCategory_start_End == "NP - NP")$Plastic))

stopifnot(all(subset(dr, Generation==2 & PlasticCategory_start_End == "P - P")$Plastic))
stopifnot(!any(subset(dr, Generation==2 & PlasticCategory_start_End == "NP - P")$Plastic))
stopifnot(all(subset(dr, Generation==2 & PlasticCategory_start_End == "P - NP")$Plastic))
stopifnot(!any(subset(dr, Generation==2 & PlasticCategory_start_End == "NP - NP")$Plastic))



dd = subset(dr, PlasticCategory_start_End == "P - P")
lost = rep(FALSE, length(unique(dd$ids)))
for (id_index in 1:length(unique(dd$ids)))
{
	id = unique(dd$ids)[id_index]
	ddd = subset(dd, ids==id)
	
	if (any(!ddd$Plastic))
	{
		lost[id_index] = TRUE
	}
}


table(dd$PlasticCategory_start_End)

t.test(dd$sdP~dd$PlasticCategory_start_End)


ggplot(dr, aes(x=0, y=meanW, colour=PlasticCategory__End)) + geom_boxplot() + ggtitle("") + facet_grid(.~environment+signal) + theme_classic(12) + theme(axis.text.y = element_text(size=16)) + rmXAxis + ylab("")+ theme(legend.title=element_blank()) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) 



##############################################
###### Table for description in results ######
##############################################


write.table(
	as.data.table(subset(dr, isLastGeneration))[,{nbPP = sum(grepl("end up P",PlasticCategory_start_End)); list(nbPP=nbPP, nbNP=.N-nbPP, fractionPP=nbPP / .N)},list(environment, signal)]
	,"Table2.txt",row.names=FALSE, quote=FALSE,sep='\\t')

as.data.table(subset(dr, isLastGeneration))[,
	{
		nbPP = sum(grepl("end up P",PlasticCategory_start_End))
		list(nbPP=nbPP, nbNP=.N-nbPP, fractionPP=nbPP / .N)
	},
	list(environment, signal)
]

fisher.test(matrix(c(0,200,76,124),  byrow=T, ncol=2)) # No signal vs perf. signal
fisher.test(matrix(c(0,200,15,185), byrow=T, ncol=2)) # No signal vs env. signal
fisher.test(matrix(c(15,185,76,124), byrow=T, ncol=2)) # env. signal vs perf. signal





###################################################
###### Evolution of developmental robustness ######
###################################################

drAverage = as.data.frame(as.data.table(dr)[,
	{
		meanmeanW = mean(meanW)
		SEmeanW = sd(meanW) / sqrt(.N)
		meansdP = mean(sdP)
		SEsdP = sd(sdP) / sqrt(.N)
		list(
			meanW 	= meanmeanW,
			meanW_lowerCI = meanmeanW - 1.96 * SEmeanW,
			meanW_upperCI = meanmeanW + 1.96 * SEmeanW,
			SEmeanW = SEmeanW,

			sdP 	= meansdP,
			sdP_lowerCI = meansdP - 1.96 * SEsdP,
			sdP_upperCI = meansdP + 1.96 * SEsdP,
			SEsdP = SEsdP,
			nbObservations = .N
		)
	},
	list(
		Generation,
		PlasticCategory_start_End,
		Project,
		environment,
		signal
	)
])



d = subset(drAverage, grepl("NP -",PlasticCategory_start_End))
d = subset(d, environment == "Spatial heterogeneity")
d = subset(d, signal != "No signal" & signal != "Both signals")
d <- d %>%
  group_by(Project) %>%
  mutate(
    width = 0.001 * n()
  )


ggplot(d, aes(x=Generation, y=sdP, colour=PlasticCategory_start_End, width=width)) + geom_point(size=2,position = position_dodge(0.5)) + ggtitle("") + scale_linetype_manual(values=c("solid", "solid", "solid", "solid")) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) + facet_grid(.~environment+signal) + theme_classic(22)  + ylab("")+ theme(legend.title=element_blank()) + geom_line() + theme(legend.position = c(0.85, 0.8),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) + geom_errorbar(aes(ymax = sdP_lowerCI, ymin = sdP_upperCI), size=0.5, position = pos) + scale_x_continuous(breaks=c(0,50e3,100e3, 150e3, 200e3), labels = c("0","50k", "100k", "150k", "200k")) + scale_x_log10()




ggplot(subset(d, Generation %in% c(2, 5e3, 5e4, 1e5, 2e5)), aes(x=Generation, y=sdP, colour=PlasticCategory_start_End, width=width)) + geom_point(size=2,position = position_dodge(0.5)) + ggtitle("") + scale_linetype_manual(values=c("solid", "solid", "solid", "solid")) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) + facet_grid(.~environment+signal) + theme_classic(22)  + ylab("")+ theme(legend.title=element_blank()) + geom_line() + theme(legend.position = c(0.85, 0.8),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) + geom_errorbar(aes(ymax = sdP_lowerCI, ymin = sdP_upperCI), size=0.5, position = pos) + scale_x_continuous(breaks=c(0,50e3,100e3, 150e3, 200e3), labels = c("0","50k", "100k", "150k", "200k")) 

ggplot(d, aes(x=Generation, y=sdP, colour=PlasticCategory_start_End, width=width)) + geom_point(size=2,position = position_dodge(0.5)) + ggtitle("") + scale_linetype_manual(values=c("solid", "solid", "solid", "solid")) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) + facet_grid(.~environment+signal) + theme_classic(22)  + ylab("")+ theme(legend.title=element_blank()) + geom_line() + theme(legend.position = c(0.85, 0.8),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) + geom_errorbar(aes(ymax = sdP_lowerCI, ymin = sdP_upperCI), size=0.5, position = pos) + scale_x_continuous(breaks=c(0,50e3,100e3, 150e3, 200e3), labels = c("0","50k", "100k", "150k", "200k"))

ggplot(subset(d, Generation >= 5e3), aes(x=Generation, y=sdP, colour=PlasticCategory_start_End, width=width)) + geom_point(size=2,position = position_dodge(0.5)) + ggtitle("") + scale_linetype_manual(values=c("solid", "solid", "solid", "solid")) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) + facet_grid(.~environment+signal) + theme_classic(22)  + ylab("")+ theme(legend.title=element_blank()) + geom_line() + theme(legend.position = c(0.85, 0.8),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) + geom_errorbar(aes(ymax = sdP_lowerCI, ymin = sdP_upperCI), size=0.5, position = pos) + scale_x_continuous(breaks=c(0,50e3,100e3, 150e3, 200e3), labels = c("0","50k", "100k", "150k", "200k"))




dtmp = subset(dr, signal == "Env. Signal" & environment == "Spatial heterogeneity" & (Generation == 100000 | Generation == 50000))
dtmpEarly = subset(dtmp, Generation == 50000)
dtmpLate = subset(dtmp, Generation == 100000)
dtmp2 = data.frame(
	sdP_at_50k = dtmpEarly$sdP,
	sdP_at_100k = dtmpLate$sdP
)

dtmp2$category = ""
for (row in 1:nrow(dtmp2))
{
	if (!dtmpEarly$Plastic[row] & !dtmpLate$Plastic[row])
	{
		dtmp2$category[row] = "NPP->NPP"
	} else if (dtmpEarly$Plastic[row] & dtmpLate$Plastic[row])
	{
		dtmp2$category[row] = "PP->PP"
	} else if (!dtmpEarly$Plastic[row] & dtmpLate$Plastic[row])
	{
		dtmp2$category[row] = "NPP->PP"
	} else
	{
		dtmp2$category[row] = "PP->NPP"
	}
}
	
ggplot(subset(dtmp2, category == "NPP->NPP" | category == "PP->PP"), aes(x=sdP_at_50k, y=sdP_at_100k, color=category)) + geom_point() + geom_abline(intercept=0, slope=1) + xlab("Developmental noise at 50k generations") + ylab("Developmental noise at 100k generations") + ggtitle("Environmental signal")





dtmp = subset(dr, signal == "Env. Signal" & environment == "Spatial heterogeneity" & (Generation == 100000 | Generation == 4000))
dtmpEarly = subset(dtmp, Generation == 4000)
dtmpLate = subset(dtmp, Generation == 100000)
dtmp2 = data.frame(
	sdP_at_50k = dtmpEarly$sdP,
	sdP_at_100k = dtmpLate$sdP
)

dtmp2$category = ""
for (row in 1:nrow(dtmp2))
{
	if (!dtmpEarly$Plastic[row] & !dtmpLate$Plastic[row])
	{
		dtmp2$category[row] = "NPP->NPP"
	} else if (dtmpEarly$Plastic[row] & dtmpLate$Plastic[row])
	{
		dtmp2$category[row] = "PP->PP"
	} else if (!dtmpEarly$Plastic[row] & dtmpLate$Plastic[row])
	{
		dtmp2$category[row] = "NPP->PP"
	} else
	{
		dtmp2$category[row] = "PP->NPP"
	}
}
	
ggplot(subset(dtmp2, category == "NPP->NPP" | category == "PP->PP"), aes(x=sdP_at_50k, y=sdP_at_100k, color=category)) + geom_point() + geom_abline(intercept=0, slope=1) + xlab("Developmental noise at 4k generations") + ylab("Developmental noise at 100k generations") + ggtitle("Environmental signal")




drAverage = as.data.frame(as.data.table(subset(dr, grepl("NP -",PlasticCategory_start_End)))[,
	{
		meanmeanW = mean(meanW)
		SEmeanW = sd(meanW) / sqrt(.N)
		meansdP = mean(sdP)
		SEsdP = sd(sdP) / sqrt(.N)
		list(
			meanW 	= meanmeanW,
			meanW_lowerCI = meanmeanW - 1.96 * SEmeanW,
			meanW_upperCI = meanmeanW + 1.96 * SEmeanW,
			SEmeanW = SEmeanW,

			sdP 	= meansdP,
			sdP_lowerCI = meansdP - 1.96 * SEsdP,
			sdP_upperCI = meansdP + 1.96 * SEsdP,
			SEsdP = SEsdP,
			nbObservations = .N
		)
	},
	list(
		Generation,
		Plastic,
		Project,
		environment,
		signal
	)
])



d = subset(drAverage, environment == "Spatial heterogeneity")
d = subset(d, signal != "No signal" & signal != "Both signals")
d <- d %>%
  group_by(Project) %>%
  mutate(
    width = 0.08 * n()
  )

ggplot(d, aes(x=Generation, y=sdP, colour=Plastic, width=width)) + geom_point(size=2,position = position_dodge(0.5)) + ggtitle("") + scale_linetype_manual(values=c("solid", "solid", "solid", "solid")) + scale_colour_manual(values=c("grey", "black"), labels=c("Not Plastic", "Plastic")) + facet_grid(.~environment+signal) + theme_classic(22)  + ylab("")+ theme(legend.title=element_blank()) + geom_line() + theme(legend.position = c(0.85, 0.8),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) + geom_errorbar(aes(ymax = sdP_lowerCI, ymin = sdP_upperCI), size=0.5, position = pos) + scale_x_continuous(breaks=c(0,50e3,100e3, 150e3, 200e3), labels = c("0","50k", "100k", "150k", "200k"))





#################################################
###### Fraction of plastic genotypes again ######
#################################################


d = as.data.frame(as.data.table(subset(dr, isLastGeneration & grepl("Constant", environment)))[,
	{
		nbCases = .N
		nbPlastic = sum(grepl("- P",PlasticCategory_start_End))
		p = nbPlastic / nbCases
		pp = (nbPlastic + 2) / (nbCases + 4)
		SE = sqrt((pp*(1-pp))/(nbCases+4))

		list(
			proportionPlastic = p,
			upSE = pp + SE,
			downSE = pp - SE
			)
	},
	list(signal)
	]
)


ggplot(d, aes(y=proportionPlastic, x=signal)) + geom_bar(stat="identity", fill="grey", colour='black') + geom_errorbar(aes(ymin = upSE, ymax=downSE), width=0.3) + theme_classic(20) + ylab("Proportion of plastic genotypes") + xlab(NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), title=element_text(size=15)) + scale_x_discrete(labels=c("No Signal","Env. Signal","Perf. Signal"))

