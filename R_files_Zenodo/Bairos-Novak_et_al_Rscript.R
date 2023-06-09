#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-     Analysis acompanying: Trust thy neighbour in times of trouble:     #-#-
#-#-                  background risk alters how tadpoles                   #-#-
#-#-                release and respond to disturbance cues                 #-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#
# Kevin Bairos-Novak, 2017
# Run on a Mac using R v.3.2.3
# Please do not re-distribute 
# Request permission to use by emailing: kevin.bairos-novak'at'usask.ca
#
### Initial set-up -------------------------------------------------------------
##
## House-keeping:
# 
setwd("~/Desktop") # Change to your directory with dataset + .R script
options(contrasts = c("contr.sum", "contr.poly")) # Set proper ANOVA Type III SS
require('gdata') # for loading data from an Excel spreadsheet
require('magrittr') # for data manipulation
require('dplyr') # for data manipulation
require('tidyr') # for data manipulation
require('car') # R-M ANOVA design
require('afex') # wrapper for easy R-M ANOVA design with car package
require('ggplot2') # for plotting data
require('psych') # for summarizing data
require('lattice') # for outlier detection plots
#
## Function to return test statistics + P-values from afex/car object:
resultz <- function (model, row=NULL, report=c("overall","within"), 
		     format="short") {
	a <- NULL
	b <- NULL
	format <- match.arg(format, c(1,2,"short", "long"))
	report <- match.arg(report, 
			    choices = c("all", "overall","within","none",
			    	    "n","N","k", "sample.size",
			    	    "r","R","r^2","r.squared","adj.r.squared"), 
			    several.ok=T)
	model.type <- match.arg(is(model),choices=c("lm", "aov", "anova", "Anova.mlm", "data.frame", "htest"), several.ok=T)
	## Within part for Anova.mlm type model or summary of Anova.mlm model
	if("htest" %in% model.type) {
		b <- model %$% 
			paste0("t(",(parameter %>% as.numeric),") = ",
			       (statistic %>% as.numeric %>% round(2)),
			       ", P = ", (p.value %>% as.numeric %>% signif(2)) )
		return(return(cat(b)))
	}
	if("Anova.mlm" %in% model.type | "anova" %in% model.type){
		if("Anova.mlm" %in% model.type){
			model <- summary(model, multivariate=T)$univariate	
		}
		x <- matrix(nrow=dim(model)[1], ncol=4)
		x[,1] <- model[,which(colnames(model)=="num Df")]
		x[,2] <- model[,which(colnames(model)=="den Df")]
		x[,3] <- model[,which(colnames(model)=="F")]
		x[,4] <- model[,which(colnames(model)=="Pr(>F)")]
		rownames(x) <- rownames(model)
		if(!is.null(row)){x <- x[row,,drop=F]}
		c. <- NULL
		d <- NULL
	}
	## Within part reading from ezANOVA and aov_car outputs
	if("data.frame" %in% model.type & !("anova" %in% model.type)){
		as.matrix(model)
		x <- matrix(nrow=dim(model)[1], ncol=4)
		x.df <- model[grep("df", colnames(model), ignore.case=T)]
		if(is.character(x.df[1,1])){
			x[,1:2] <- ( apply(x.df,1,strsplit, split=", ") %>% unlist %>% 
				     	as.numeric %>% matrix(nrow=dim(model)[1], byrow=TRUE) )
			x[,3] <- model[,which(colnames(model)=="F")] %>% matrix %>% apply(FUN=substr,1, start=1,stop=4) %>% as.numeric
			suppressWarnings(x[,4] <- model[,which(colnames(model)=="p.value")] %>% as.numeric)
			for(i in 1:length(which(is.na(x[,4])))) {
				j <- which(is.na(x[,4]))[i]
				if(model[j,which(colnames(model)=="p.value")]==c("<.0001")) x[j,4] <- 0.0001
			}
		} else {
			x[,1:2] <- x.df %>% as.matrix
			x[,3] <- model[,which(colnames(model)=="F")]
			x[,4] <- model[,which(colnames(model)=="p")]
		}
		rownames(x) <- model[,which(colnames(model)=="Effect")]
		if(!is.null(row)){x <- x[row,,drop=F]}
		c. <- NULL
		d <- NULL
	} 
	if("lm" %in% model.type) {## Model if from lm object
		## A: overall F-statistic
		fval <- summary(model)$fstatistic
		pval <- 1 - pf(fval[1],fval[2],fval[3])
		if("all" %in% report | "overall" %in% report){
			a <- c()
			# can also do: 1- do.call(pf, unname(as.list(fval)))
			if(format == "short" | format == 1){
				a <- paste0("overall model F(",fval[2],",",fval[3],") = ",
					    round(fval[1],2),", P = ", signif(pval,2), if(pval < 0.05) {
					    	c("  *\\n\\n")} else {"\\n\\n"})
			}
			if(format == "long" | format == 2){
				a <- paste0("overall model F = ",round(fval[1],2),"; df = ",fval[2],
					    ", ", fval[2], "; P = ", signif(pval,2), if(pval < 0.05) {
					    	c("  *\\n\\n")} else {"\\n\\n"})
			}
		} else {a <- NULL}	
		## B: within-model
		x <- subset(drop1(model, ~., test = "F"), # Type III SS for B
			    !is.na(Df), select = c(1,5,6))
		x <- cbind(x, Df2 = rep(summary(model)$df[2]))
		x <- x[,c(1,4,2,3), drop=F]
		if(!is.null(row)){x <- x[row,, drop=F]}
		## C: r^2 and adj-r^2 values (optional)
		if("all" %in% report | "r" %in% report | "R" %in% report | "r^2" %in% report | "r.squared" %in% report | "adj.r.squared" %in% report){
			c. <- paste0("\\nr^2 = ",round(summary(model)$r.squared,3),
				     "; adj. r^2 = ",round(summary(model)$adj.r.squared,3),"\\n")
		} else {c. <- NULL}
		## D: sample sizes
		if("all" %in% report | "n" %in% report | "N" %in% report| "k" %in% report | "sample.size" %in% report) {
			d <- paste0("\\nN total = ",fval[2] + fval[3] + 1,
				    "; k parameters = ",fval[2],"\\n")
		} else {d <- NULL}
	}
	
	## B: Print Type III SS results for all model types
	if("all" %in% report | "within" %in% report){
		b <- c()
		if(format == "short" | format == 1){
			if(!is.null(row)){ # If row is specified, don't use name
				b <- paste0("F(", x[,1],",",x[,2], ") = ", 
					    round(x[,3], 2),", P = ", 
					    signif(x[,4],2)
				)
			} else { # If row is not specified, display names
				for(i in 1:dim(x)[1]) {
					b[i] <- paste0(
						rownames(x)[i], " F(", x[i,1],",", 
						x[i,2], ") = ", round(x[i,3], 2), 
						", P = ", signif(x[i,4],2), 
						if(x[i,4] < 0.05) {c("  *\\n")} else {"\\n"} # Signif indicator
					)
				}	
			}
			
		}
		if(format == "long" | format == 2){
			if(!is.null(row)){ # If row is specified, don't use name
				b <- paste0("F = ", round(x[,3], 2), 
					    "; df = ", x[,1],", ", x[,2],
					    "; P = ", 
					    signif(x[,4],2)
				)
			} else { # If row is not specified, display names
				for(i in 1:dim(x)[1]) {
					b[i] <- paste0(
						rownames(x)[i], " F = ", round(x[i,3], 2), 
						"; df = ", x[i,1],", ", x[i,2],
						"; P = ", signif(x[i,4],2), 
						if(x[i,4] < 0.05) {c("  *\\n")} else {"\\n"} # Signif indicator
					)
				}
			}
		}
	} else {b <- NULL}
	cat(a,b,c.,d)
	
}
#
## Load data and pre-processing:
data <- read.xls("Bairos-Novak_et_al_rawdata.xlsx", sheet = 1, header = TRUE)
data <- data[!is.na(data$pre),,drop=F] # non-existant row removal
data %<>% within( { # factorize categorical variables appearing as integers
	trial   <- factor(trial); day <- factor(day); who <- factor(who)
	receiver.risk <- factor(receiver.risk);	donor.risk  <- factor(donor.risk)
	disturb <- factor(disturb)
})
data <- data[,c("trial", "receiver.risk", "donor.risk", "disturb", "pre", "post"),drop=F] # remove extraneous columns
data2 <- gather(data,"time","lines.crossed", pre, post) # Construct vertical dataset for R-M approach
data2$time %<>% factor(labels=c("pre","post"), ordered=TRUE)
#
#
### Outliers and assumptions ---------------------------------------------------
# ##
# ## Multi-panel barplots:
# data3 <- data2 %>% within( { # factorize categorical variables
# 	receiver.risk <- factor(receiver.risk, levels = c(0,1), labels = c("Low receiver risk", "High receiver risk"))
# 	donor.risk  <- factor(donor.risk, levels = c(0,1), labels = c("Low cue donor risk", "High cue donor risk"))
# 	disturb <- factor(disturb, levels = c(0,1), labels = c("Undisturbed", "Disturbance cue"))
# })
# bwplot(lines.crossed ~   receiver.risk | donor.risk, data = data3,
#        # strip = strip.custom(bg = 'white'),
#        cex = .5, layout = c(1, 2),
#        xlab = "", ylab = "Lines Crossed",
#        par.settings = list(
#        	box.rectangle = list(col = 1),
#        	box.umbrella  = list(col = 1),
#        	plot.symbol   = list(cex = .5, col = 1)),
#        scales = list(x = list(relation = "same"),
#        	      y = list(relation = "same")))
# bwplot(lines.crossed ~ disturb  | receiver.risk * donor.risk, data = data3,
#        # strip = strip.custom(bg = 'white'),
#        cex = .5, layout = c(2, 2),
#        xlab = "", ylab = "Lines Crossed",
#        par.settings = list(
#        	box.rectangle = list(col = 1),
#        	box.umbrella  = list(col = 1),
#        	plot.symbol   = list(cex = .5, col = 1)),
#        scales = list(x = list(relation = "same"),
#        	      y = list(relation = "same")))
# bwplot(lines.crossed ~ disturb | time * receiver.risk * donor.risk, data = data3,
#        # strip = strip.custom(bg = 'white'),
#        cex = .5, layout = c(2, 2),
#        xlab = "", ylab = "Lines Crossed",
#        par.settings = list(
#        	box.rectangle = list(col = 1),
#        	box.umbrella  = list(col = 1),
#        	plot.symbol   = list(cex = .5, col = 1)),
#        scales = list(x = list(relation = "same"),
#        	      y = list(relation = "same")))
# #
# ## Distributions:
# par(mfrow=c(2,2))
# data %>% filter(receiver.risk==0 & donor.risk==0 & disturb==0) %$% 
# 	hist(pre, col=rgb(0.1,0.1,0.1,0.5), xlim = c(0,60),
# 	     main = "Pre-exposure low-risk receivers")
# data %>% filter(receiver.risk==0 & donor.risk==0 & disturb==1) %$% 
# 	hist(pre, add=T, col=rgb(0.8,0.8,0,0.5))
# data %>% filter(receiver.risk==0 & donor.risk==1 & disturb==0) %$% 
# 	hist(pre, add=T, col=rgb(0,0,0.9,0.5))
# data %>% filter(receiver.risk==0 & donor.risk==1 & disturb==1) %$% 
# 	hist(pre, add=T, col=rgb(0.9,0,0,0.5))
# data %>% filter(receiver.risk==1 & donor.risk==0 & disturb==0) %$% 
# 	hist(pre, col=rgb(0.1,0.1,0.1,0.5), xlim = c(0,50),
# 	     main = "Pre-exposure high-risk receivers")
# data %>% filter(receiver.risk==1 & donor.risk==0 & disturb==1) %$%
# 	hist(pre, add=T, col=rgb(0.8,0.8,0,0.5))
# data %>% filter(receiver.risk==1 & donor.risk==1 & disturb==0) %$% 
# 	hist(pre, add=T, col=rgb(0,0,0.9,0.5))
# data %>% filter(receiver.risk==1 & donor.risk==1 & disturb==1) %$% 
# 	hist(pre, add=T, col=rgb(0.9,0,0,0.5))
# data %>% filter(receiver.risk==0 & donor.risk==0 & disturb==0) %$% 
# 	hist(post, col=rgb(0.1,0.1,0.1,0.5), xlim = c(0,60),
# 	     main = "Post-exposure low-risk receivers")
# data %>% filter(receiver.risk==0 & donor.risk==0 & disturb==1) %$% 
# 	hist(post, add=T, col=rgb(0.8,0.8,0,0.5))
# data %>% filter(receiver.risk==0 & donor.risk==1 & disturb==0) %$% 
# 	hist(post, add=T, col=rgb(0,0,0.9,0.5))
# data %>% filter(receiver.risk==0 & donor.risk==1 & disturb==1) %$% 
# 	hist(post, add=T, col=rgb(0.9,0,0,0.5))
# data %>% filter(receiver.risk==1 & donor.risk==0 & disturb==0) %$% 
# 	hist(post, col=rgb(0.1,0.1,0.1,0.5), xlim = c(0,40),
# 	     main = "Post-exposure high-risk receivers")
# data %>% filter(receiver.risk==1 & donor.risk==0 & disturb==1) %$%
# 	hist(post, add=T, col=rgb(0.8,0.8,0,0.5))
# data %>% filter(receiver.risk==1 & donor.risk==1 & disturb==0) %$% 
# 	hist(post, add=T, col=rgb(0,0,0.9,0.5))
# data %>% filter(receiver.risk==1 & donor.risk==1 & disturb==1) %$% 
# 	hist(post, add=T, col=rgb(0.9,0,0,0.5))
# par(mfrow=c(1,1))
# #
# ## Homogeneity of variances assumption:
# unite_(data2, "time_rec_cue_dist",c("time", "receiver.risk","donor.risk","disturb")) %$% 
# 	leveneTest(y = lines.crossed, group = time_rec_cue_dist)
# unite_(data2, "time_rec_cue_dist",c("time", "receiver.risk","donor.risk","disturb")) %$% 
# 	boxplot(lines.crossed ~ time_rec_cue_dist, col=c(rep(2,8),rep(3,8))) # Christmas box plot
# #
# ## Residual normality assumption:
# # Not the best, but ANOVA is robust to small violations of assumptions
# m1.resid <- aov_car(
# 	lines.crossed ~ receiver.risk * donor.risk * disturb + Error(trial/time)
# 	, data = data2, return="lm")$resid
# qqnorm(m1.resid); qqline(m1.resid)
# plot(m1.resid[,1],data$trial, col=4); points(m1.resid[,2],data$trial, col=2); abline(v=0, lwd=2)
# #
# #
### 3-way Pre-exposure data analysis --------------------------------------------------

## No difference in pre-exposure line crosses among treatments:
m0 <- lm(pre ~ receiver.risk*donor.risk*disturb, data=data)
# No significant overall model: F(7,258) = 0.54, P = 0.80
resultz(m0, report="overall")
#
# No difference in pre-exposure line crosses between receiver risk levels:
m0.rr <- lm(pre ~ receiver.risk, data=data)
resultz(m0.rr, row="receiver.risk","within")
# F(1,261) = 1.91, P = 0.17
# 
# 
### 4-way (1 within-subjects, 3 between-subjects) R-M ANOVA --------------------
##
## Test for 4-way interaction:
m1 <- aov_car(lines.crossed ~ receiver.risk * donor.risk * disturb + Error(trial/time)
	      , data = data2, return=afex_options(return_aov="nice"))
# Significant 4-way interaction: F(1,258) = 6.42, P = 0.012
resultz(m1, "receiver.risk:donor.risk:disturb:time") # ges = .007
#
#
## Test 3-way interaction for only undisturbed cues:
m1.uc <- aov_car(
	lines.crossed ~ receiver.risk * donor.risk + Error(trial/time)
	, data = subset(data2, disturb==0), return=afex_options(return_aov="nice"))
# Non-significant 3-way interaction (i.e. slopes are the same): F(1,129) = 0.34, P = 0.56
resultz(m1.uc, "receiver.risk:donor.risk:time")
#
# No other 2-way interactions significant, except time:
resultz(m1.uc)
# Time variable effect size calculation:
data %>% subset(disturb==0) %>% 
	with(., ((mean(pre) - mean(post))/mean(pre))*100) %>% 
	round(digits=2) %>%
	paste0("% decrease in lines crossed from pre to post") %>% cat
# 18.96% decrease in lines crossed from pre to post for undisturbed cues
# 
#
## Test 3-way interaction using only disturbance cues:
m1.dc <- aov_car(lines.crossed ~ receiver.risk * donor.risk + Error(trial/time)
		 , data=subset(data2, disturb==1), return=afex_options(return_aov="nice"))
# Significant 3-way interaction: F(1,129) = 8.59, P = 0.004
resultz(m1.dc, "receiver.risk:donor.risk:time") # ges = .02
#
#
## R-M ANOVA/t-tests to compare differences in time for disturbed vs. undisturbed:
# Note: t-value^2 = F-value
#
# High receiver risk, High DC donor risk through time:
data.HrrHcr <- subset(data, disturb==1 & receiver.risk==1 & donor.risk==1)
t.HrrHcr <- with(data.HrrHcr, t.test(pre, post, paired=T))
data2.HrrHcr <- subset(data2, disturb==1 & receiver.risk==1 & donor.risk==1)
m1.dc.HrrHcr <- aov_car(lines.crossed ~ Error(trial/time), 
			data = data2.HrrHcr, return=afex_options(return_aov="nice"))
# Significant effect of time for high rr, high cr: 
resultz(m1.dc.HrrHcr, "time") 	# F(1,33) = 19.84, P = 9.1e-05
resultz(t.HrrHcr) 		# t(33) = 4.45, P = 9.1e-05
# Effect size of time:
with(data.HrrHcr, ((mean(pre) - mean(post))/mean(pre))*100) %>% round(digits=2) %>%
	paste0("% decrease in lines crossed from pre to post") %>% cat
# 34.98% decrease in lines crossed from pre to post
# 
## Low receiver risk, High DC donor risk through time
data.LrrHcr <- subset(data, disturb==1 & receiver.risk==0 & donor.risk==1)
t.LrrHcr <- with(data.LrrHcr, t.test(pre, post, paired=T))
data2.LrrHcr <- subset(data2, disturb==1 & receiver.risk==0 & donor.risk==1)
m1.dc.LrrHcr <- aov_car(lines.crossed ~ Error(trial/time), 
			data = data2.LrrHcr, return=afex_options(return_aov="nice"))
# Significant effect of time for low rr, high cr: 
resultz(m1.dc.LrrHcr, "time") 	# F(1,32) = 44.37, P = 1.6e-07
resultz(t.LrrHcr) 		# t(32) = 6.66, P = 1.6e-07
# Effect size of time:
with(data.LrrHcr, ((mean(pre) - mean(post))/mean(pre))*100) %>% round(digits=2) %>%
	paste0("% decrease in lines crossed from pre to post") %>% cat
# 48.64% decrease in lines crossed from pre to post
# 
## High receiver risk, Low DC donor risk through time
data.HrrLcr <- subset(data, disturb==1 & receiver.risk==1 & donor.risk==0)
t.HrrLcr <- with(data.HrrLcr, t.test(pre, post, paired=T))
data2.HrrLcr <- subset(data2, disturb==1 & receiver.risk==1 & donor.risk==0)
m1.dc.HrrLcr <- aov_car(lines.crossed ~ Error(trial/time), 
			data = data2.HrrLcr, return=afex_options(return_aov="nice"))
# Significant effect of time for high rr, low cr: 
resultz(m1.dc.HrrLcr, "time") 	# F(1,32) = 21.5, P = 5.7e-05
resultz(t.HrrLcr) 		# t(32) = 4.64, P = 5.7e-05
# Effect size of time:
with(data.HrrLcr, ((mean(pre) - mean(post))/mean(pre))*100) %>% round(digits=2) %>%
	paste0("% decrease in lines crossed from pre to post") %>% cat
# 45.18% decrease in lines crossed from pre to post
# 
## Low receiver risk, Low DC donor risk through time
data.LrrLcr <- subset(data, disturb==0 & receiver.risk==0 & donor.risk==0)
t.LrrLcr <- with(data.LrrLcr, t.test(pre, post, paired=T))
data2.LrrLcr <- subset(data2, disturb==1 & receiver.risk==0 & donor.risk==0)
m1.dc.LrrLcr <- aov_car(lines.crossed ~ Error(trial/time), 
			data = data2.LrrLcr, return=afex_options(return_aov="nice"))
# No effect of time for low rr, low cr: 
resultz(m1.dc.LrrLcr, "time") 	# F(1,32) = 0.87, P = 0.36
resultz(t.LrrLcr) 		# t(32) = 0.93, P = 0.36
#
#
## Test whether positive responses differed in magnitude using 2-way RM-ANOVAs:
#
# Compare within high cue donor risk: 
m1.dc.Hcr <- aov_car(lines.crossed ~ receiver.risk + Error(trial/time), 
			data = subset(data2, disturb==1 & donor.risk==1), 
			return=afex_options(return_aov="nice"))
# No difference in receiver responses to DCs from high-risk cue donors 
# (remember, P-value doubled due to multiple comparisons)
resultz(m1.dc.Hcr, "receiver.risk:time") # F(1,65) = 4.28, P = 0.043
#
# Compare within high receiver risk:
m1.dc.Hrr <- aov_car(lines.crossed ~ donor.risk + Error(trial/time), 
		     data = subset(data2, disturb==1 & receiver.risk==1), 
		     return=afex_options(return_aov="nice"))
# No difference in high-risk receiver responses to DCs
resultz(m1.dc.Hrr, "donor.risk:time") # F(1,65) = 1.22, P = 0.27
#

# # Effect size of equally-significant DC treatments:
# data %>% subset(disturb==1 & !(receiver.risk==0 & donor.risk==0) ) %>% 
# 	with(., ((mean(pre) - mean(post))/mean(pre))*100) %>% 
# 	round(digits=2) %>%
# 	paste0("% decrease in lines crossed from pre to post") %>% cat
# # 43.31% decrease in lines crossed from pre to post
#
#
### Plot R-M ANOVA results -----------------------------------------------------
##
## Set-up of data:
#
# Create a data frame for mean and standard errors of each bar in graph
data2 <- gather(data,"time","lines.crossed", pre, post)
data2$time <- factor(data2$time, levels=c("pre","post"), labels=c(0,1))
df  <- data2 %$% 
	describeBy(lines.crossed, group=list(disturb, donor.risk, receiver.risk, time), digits=2,
		   mat=T)[,c("group1", "group2", "group3", "group4", "mean","se", "n")]
colnames(df) <- c("disturb", "donor.risk", "receiver.risk",  "time", "mean", "se", "n")
#
# Create a custom theme for plot:
dotplot.theme <- 
	theme_bw() +
	theme(
		plot.margin = unit(c(20,20,4,20), "pt"), # margin space around entire plot
		panel.grid.major = element_blank(), # major panel grid lines
		panel.grid.minor = element_blank(), # minor grid lines
		panel.border = element_rect(colour = "black", size = 1), # panel borders
		panel.background = element_blank(), # background colour of plot
		axis.title.y = element_text(size = 17, margin = unit(c(0,20,0,0), "pt")),
		axis.title.x = element_blank(), # remove x-axis title
		axis.text.y = element_text(size = 13, margin = unit(c(0,2,0,5), "pt")), # y-axis numbers
		axis.text.x = element_blank(), # x-axis categories
		axis.ticks.x = element_blank(),
		axis.line = element_line(), # weren't showing up after facetting
		strip.background = element_blank(), # facet panel background
		strip.text.x = element_text(size = 17, margin = unit(c(0,0,20,0), "pt")), # x-axis facet variable title
		strip.text.y = element_text(size = 17, margin = unit(c(0,0,0,20), "pt")), # y-axis facet variable title
		panel.spacing = unit(-0.09, "lines"), # space between facet panels
		legend.position = "bottom",
		legend.text = element_text(size = 15),
		legend.key = element_rect(colour = "white", fill = NA), # remove grey rectangle around legend symbols
		legend.key.size = unit(50, "pt") # increase legend symbol size
	)
#
## Begin creating plotting object:
#
# Base plot of data frame 
p <- ggplot(df, aes(x = time, y = mean, shape = disturb, group = disturb))
#
# Add themes, geoms, and scaling
dp <- p +
	guides(size=FALSE) + # removes extreneous legend elements (size in this case)
	dotplot.theme +
	geom_point(aes(shape = disturb), size = 3) +
	geom_line(aes(linetype = disturb)) +
	facet_grid(receiver.risk ~ donor.risk,
		   #switch = "y", # switch y or x labels to opposite side
		   labeller = labeller(
		   	receiver.risk = as_labeller(c('0' = "Low-Risk Receivers",
		   				      '1' = "High-Risk Recievers")),
		   	donor.risk = as_labeller(c('0' = "Low-Risk Cue Donor",
		   				 '1' = "High-Risk Cue Donor")))) +
	geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
		      colour = 'black', width = 0.1) +
	scale_y_continuous(limits = c(7,20.5), breaks = seq(0,20, by=5)) +
	# scale_x_discrete(breaks = c("0", "1"),labels = c("Pre", "Post")) +
	scale_shape_manual(values = c(1,19),
			   name = "",
			   labels = c("Undisturbed Cue    ", "Disturbance Cue")) +
	scale_linetype_manual(values = c(2,1),
			      name = "",
			      labels = c("Undisturbed Cue    ", "Disturbance Cue")) +
	labs(list(x=NULL, y = "Number of Lines Crossed")) +
	geom_text(data = data.frame(
		disturb = factor(c(0,0,0), levels = levels(df$disturb)),
		donor.risk = factor(c(1,0,1), levels = levels(df$donor.risk)),
		receiver.risk = factor(c(0,1,1), levels = levels(df$receiver.risk)),
		time = factor(c(1,1,1), levels = levels(df$time)),
		mean = 19, Label = c("*")), 
		aes(label = Label), size = 10, nudge_x=-0.47) +
	annotate("text", label = "Pre", size = 5, x = 1, y = 7.1) +
	annotate("text", label = "Post", size = 5, x = 2, y = 7.1)

df <- df[order(df$time,df$disturb,df$receiver.risk),,drop=F] # sorting important for sample sizes
df.uc <- df[9:12,]
df.dc <- df[13:16,]
#
# Add 'a' through 'd' to the corner of each plot:
dp <- dp +
	geom_text( 
		data = data.frame(mean=20,
				  disturb = factor(c(0,0,0,0), levels = levels(df.dc$disturb)),
				  donor.risk = factor(c(0,1,0,1), levels = levels(df.dc$donor.risk)),
				  receiver.risk = factor(c(0,0,1,1), levels = levels(df.dc$receiver.risk)),
				  time = factor(c(0,0,0,0), levels = levels(df.dc$time)),# will have to nudge
				  Label = c("a", "b", "c", "d")
		), aes(label = Label), size = 8, nudge_x=-0.45)
dp
#
# Add sample sizes
dp2 <- dp +	geom_text( # Add sample size labels
	data = data.frame(mean=df.dc$mean,
			  disturb = factor(c(0,0,0,0), levels = levels(df.dc$disturb)),
			  donor.risk = factor(c(0,1,0,1), levels = levels(df.dc$donor.risk)),
			  receiver.risk = factor(c(0,0,1,1), levels = levels(df.dc$receiver.risk)),
			  time = factor(c(1,1,1,1), levels = levels(df.dc$time)),# will have to nudge
			  Label = with(df.dc, c(paste("n=",n[1]),paste("n=",n[2]),paste("n=",n[3]),paste("n=",n[4])))
	), aes(label = Label), size = 5, nudge_x=0.35) +
	geom_text(
		data = data.frame(mean=df.uc$mean,
				  disturb = factor(c(0,0,0,0), levels = levels(df.uc$disturb)),
				  donor.risk = factor(c(0,1,0,1), levels = levels(df.uc$donor.risk)),
				  receiver.risk = factor(c(0,0,1,1), levels = levels(df.uc$receiver.risk)),
				  time = factor(c(1,1,1,1), levels = levels(df.uc$time)),# will have to nudge
				  Label = with(df.uc, c(paste("n=",n[1]),paste("n=",n[2]),paste("n=",n[3]),paste("n=",n[4])))
		), aes(label = Label), size = 5, nudge_x=0.35)
dp2

## Use the following
# setEPS()
# postscript("dotplot.eps")
# dp
# dev.off()
# setEPS()
# postscript("dotplot2.eps")
# dp2
# dev.off()
# 
