#########################################################################################################################################
### Script carrying out the analyses published in:																					  ###
### Deguines N, Brashares JS, Prugh LR - J. of Animal Ecology - Precipitation alters interactions in a grassland ecological community ###
#########################################################################################################################################


# Table of contents:
	#### A. Data loading
	#### B. Mixed-effects modelling as a preliminary step to inform and support our structural equation modelling approach
	#### C. Piecewise structural equation model
		# C.1 Evaluation of initial SEM
		# C.2 Path addition
		# C.3 Path deletion
		# C.4 Calculating the indirect, plant-mediated, effects of precipitation
	#### D. Linear regressions of annual standardized coefficients against precipitation
	
# Literature cited:
	# Bolker, B.M. (2015) Linear and generalized linear mixed models. Ecological Statistics: Contemporary Theory and Application, pp. 309–333. Oxford University Press.
	# Grace, J.B., Scheiner, S.M. & Schoolmaster, D.R. (2015) Structural equation modeling: building and evaluating causal models. Ecological statistics: contemporary theory and application, pp. 168–199. Oxford University Press, Oxford, UK.
	# Lefcheck, J.S. (2016) piecewiseSEM: Piecewise structural equation modelling in r for ecology, evolution, and systematics. Methods in Ecology and Evolution, 7, 573–579.
	# Zuur, A., Ieno, E.N. & Smith, G.M. (2007) Analyzing Ecological Data. Springer Science+Business Media, New York, USA.
	# Zuur, A., Ieno, E.N., Walker, N., Saveliev, A.A. & Smith, G.M. (2009) Mixed Effects Models and Extensions in Ecology with R. Springer Science+Business Media, New York, USA.

# Delete objects from previous R session:
rm(list=ls(all=TRUE))						

# Load required packages:
	# package piecewiseSEM: install latest version following steps at https://github.com/jslefche/piecewiseSEM
		# at time of publication, the steps were the following:
			#	install.packages("devtools",dep=TRUE)
			#	library(devtools)
			#	install.packages("stringi",dep=TRUE)
			#	library(stringi)													# may not be required
			#	install.packages("quadprog",repos="http://cran.r-project.org")
			#	library(quadprog)													# may not be required
			#	install.packages("mnormt",dep=TRUE)	# for r.squaredGLMM()
			#	library(mnormt)														# may not be required
			#	install.packages("pbivnorm",repos="http://cran.r-project.org")
			#	library(pbivnorm)													# may not be required
			#	install_github("jslefche/piecewiseSEM")
	library(piecewiseSEM)
	#install.packages("car",dep=TRUE)
	library(car)
	#install.packages("Hmisc",dep=TRUE)
	library(Hmisc)
	#install.packages("MuMIn",dep=TRUE)
	library(MuMIn)
	#install.packages("nlme",dep=TRUE)
	library(nlme)
						



####
#   A. Data loading
####

# Table Carrizo Data importation
carrizo=read.csv("Deguines_etal_JAE_carrizo.csv")
head(carrizo)

# Dataset description:
	# year
	unique(carrizo$year)			# Seven years of data (not used in models)
	# fyear
	unique(carrizo$fyear)			# year as a categorical variable ; used in lme() models
	# plot_id
	summary(carrizo$plot_id)		# identity of the 30 2ha plots ; categorical variable ; incl. as a random effect in the lme() models as each have 7 measures (one per year)
	# local_cond_id
	unique(carrizo$local_cond_id)	# identity of the two sites (not used in models)
	# local_cond_code
	unique(carrizo$local_cond_code)	# variables representing the identity of the two sites ; numerical variable (0 or 1, for CenterWell or Swain respectively) 
	# Precip
	summary(carrizo$Precip)			# Precipitation during the growing season (October - April) in cm
	# Plant_biomass
	summary(carrizo$Plant_biomass)	# biomass measured in April ; data were log() transformed
	tapply(carrizo$Plant_biomass,carrizo$fyear,mean)
	# Krats
	summary(carrizo$Krats)			# density of the giant kangaroo rat	(Dipodomys ingens)
	tapply(carrizo$Krats,carrizo$fyear,mean)
	# Ants
	summary(carrizo$Ants)			# abundance of Formicidae ; data were log(...+1) transformed
	tapply(carrizo$Ants,carrizo$fyear,mean)
	# Beetles
	summary(carrizo$Beetles)		# abundance of Coleoptera ; data were log(...+1) transformed
	tapply(carrizo$Beetles,carrizo$fyear,mean)
	# Orthopterans
	summary(carrizo$Orthopterans)	# abundance of Orthoptera ; data were log(...+1) transformed
	tapply(carrizo$Orthopterans,carrizo$fyear,mean)
	# Squirrels
	summary(carrizo$Squirrels)		# density of the San Joaquin antelope squirrel (Ammospermophilus nelsoni) ; data were log(...+1) transformed
	tapply(carrizo$Squirrels,carrizo$fyear,mean)
	# Lizards
	summary(carrizo$Lizards)		# abundance of the common side-blotched lizard (Uta stansburiana) ; data were log(...+1) transformed
	tapply(carrizo$Lizards,carrizo$fyear,mean)

# Multicollinearity among our variables ?
collinearity_tab=cbind(carrizo$Precip
			,carrizo$local_cond_code
			,carrizo$Plant_biomass,carrizo$Krats,carrizo$Ants,carrizo$Beetles,carrizo$Orthopterans,carrizo$Squirrels,carrizo$Lizards)
colnames(collinearity_tab)=c("Precip","local_cond_code","Plant_biomass","Krats","Ants","Beetles","Orthopterans","Squirrels","Lizards")
source("HighstatLibV6.R")	# "HighstatLibV6.R" is a script for the corvif function (Zuur et al. 2009). It is available online at http://www.highstat.com/book2.htm
corvif(collinearity_tab)
	#=> all GVIF values are similar and <5: no indication of collinearity among variables used in the models (Zuur et al. 2007, 2009).



	
####
#	B. Mixed-effects modelling as a preliminary step to inform and support our structural equation modelling approach
####

# Specifying control values for lme fit:
lmecontr_foc=lmeControl(maxIter=50,msMaxIter=50) #the default


	###
	# Plant biomass model
	###

m.prod.lme = lme( Plant_biomass ~ Precip*local_cond_code
	,random = ~ 1|plot_id
	,data=carrizo,method = "REML",control=lmecontr_foc)
mod=m.prod.lme

# Homosedasticity of the residuals? (Zuur et al. 2009)
plot(mod)

	# source of heteroscedasticity (if any)?
	par(mfrow=c(1,3))
	plot(carrizo$Precip,resid(mod),xlab="Precip",ylab="residuals")
	plot(carrizo$local_cond_code,resid(mod),xlab="local_cond_code",ylab="residuals")
	plot(carrizo$fyear,resid(mod),xlab="fyear",ylab="residuals")
	
		#Precip
		m.prod_hetf1=lme(Plant_biomass ~ Precip*local_cond_code
				,weights = varFixed(~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.prod_hetp1=lme(Plant_biomass ~ Precip*local_cond_code
				,weights = varPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.prod_hetcp1=lme(Plant_biomass ~ Precip*local_cond_code
				,weights = varConstPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.prod_hete1=lme(Plant_biomass ~ Precip*local_cond_code
				,weights = varExp(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.prod.lme,m.prod_hetf1,m.prod_hetp1,m.prod_hetcp1,m.prod_hete1)
		# => AIC improves best with hete1, varExp(form =~Precip).
		AIC(m.prod.lme)-AIC(m.prod_hete1) # dAIC=12.21 > 2, substantial improvement.

		#local_cond_code (weight as varIdent because a binary variable)
		m.prod_heti2=lme(Plant_biomass ~ Precip*local_cond_code
				,weights = varIdent(form = ~ 1 | local_cond_code)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.prod.lme,m.prod_heti2)
		# => AIC is not improved.

		#fyear
		m.prod_heti4=lme(Plant_biomass ~ Precip*local_cond_code
				,weights = varIdent(form = ~ 1 | fyear)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.prod.lme,m.prod_heti4)
		# => AIC improves with heti4, varIdent(form = ~ 1 | fyear)
		AIC(m.prod.lme)-AIC(m.prod_heti4) # dAIC=15.52 > 2, substantial improvement.
		
		#Combinations:
		hete1, varExp(form =~Precip)
		heti4, varIdent(form = ~ 1 | fyear)		
				
		m.prod_hete1i4=lme(Plant_biomass ~ Precip*local_cond_code
				,weights = varComb(varExp(form =~Precip),varIdent(form = ~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		tabAIC=anova(m.prod.lme,m.prod_hete1,m.prod_heti4,m.prod_hete1i4)
		tabAIC[order(tabAIC$AIC,decreasing=F),]

	# Best AIC2 model?
	 x model with no weight is:
		orimod=m.prod.lme
	 x lowest AIC is:
		low1AICmod=m.prod_heti4
		AIC(orimod)-AIC(low1AICmod) #=> dAIC = 15.52 > 2: low1AICmod model is better than orimod model
	 x 2nd lowest AIC model among the simpler ones remaining (i.e. lower df than lowest AIC model) is:
		lowsimpl2AICmod=m.prod_hete1
		AIC(orimod)-AIC(lowsimpl2AICmod) #=> dAIC = 12.21
		AIC(lowsimpl2AICmod)-AIC(low1AICmod) #=> dAIC = 3.31 > 2, low1AICmod model is better.
	#==> m.prod_heti4 is the best model.

mod=m.prod_heti4
plot(mod)	# no sign of heteroscedasticity (Zuur et al. 2009)

# Normality of the residuals?
qqnorm(mod,abline=c(0,1)) # note that "non-normality is generally less of a concern than heteroscedasticity [in linear and generalized linear mixed models]" (Bolker 2015)
	
# Model results	
anova(mod,type="marginal")
summary(mod)

# Removal of non-significant "Precip x Predictor variable" interactions ?
m.prod_heti4ML=lme(Plant_biomass ~ Precip*local_cond_code
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)	#ML
anova(m.prod_heti4ML,type="marginal") 
	# no interaction to drop.

# Model used in the following initial SEM:
m.prod_heti4REML=lme(Plant_biomass ~ Precip*local_cond_code
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	#REML
mod=m.prod_heti4REML

# Model assumptions:
plot(mod)					# Homoscedasticity of the residuals
qqnorm(mod,abline=c(0,1))	# Normality of the residuals








	###
	# Krats model
	###

m.gkr.lme = lme( Krats ~ Precip*Plant_biomass
	,random = ~ 1|plot_id
	,data=carrizo,method = "REML",control=lmecontr_foc)
mod=m.gkr.lme

# Homosedasticity of the residuals? (Zuur et al. 2009)
plot(mod)

	# source of heteroscedasticity (if any)?
	par(mfrow=c(1,3))
	plot(carrizo$Precip,resid(mod),xlab="Precip",ylab="residuals")
	plot(carrizo$Plant_biomass,resid(mod),xlab="Plant_biomass",ylab="residuals")
	plot(carrizo$fyear,resid(mod),xlab="fyear",ylab="residuals")
	
		#Precip
		m.gkr_hetf1=lme(Krats ~ Precip*Plant_biomass
				,weights = varFixed(~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.gkr_hetp1=lme(Krats ~ Precip*Plant_biomass
				,weights = varPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.gkr_hetcp1=lme(Krats ~ Precip*Plant_biomass
				,weights = varConstPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.gkr_hete1=lme(Krats ~ Precip*Plant_biomass
				,weights = varExp(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.gkr.lme,m.gkr_hetf1,m.gkr_hetp1,m.gkr_hetcp1,m.gkr_hete1)
		# => AIC is not improved.

		#Plant_biomass
		m.gkr_hetf2=lme(Krats ~ Precip*Plant_biomass
				,weights = varFixed(~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.gkr_hetp2=lme(Krats ~ Precip*Plant_biomass
				,weights = varPower(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.gkr_hetcp2=lme(Krats ~ Precip*Plant_biomass
				,weights = varConstPower(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.gkr_hete2=lme(Krats ~ Precip*Plant_biomass
				,weights = varExp(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.gkr.lme,m.gkr_hetf2,m.gkr_hetp2,m.gkr_hetcp2,m.gkr_hete2)
		# => AIC is not improved.

		#fyear
		m.gkr_heti3=lme(Krats ~ Precip*Plant_biomass
				,weights = varIdent(form = ~ 1 | fyear)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.gkr.lme,m.gkr_heti3)
		# => AIC improves with heti3, varIdent(form = ~ 1 | fyear):
		AIC(m.gkr.lme)-AIC(m.gkr_heti3) # dAIC=9.50 > 2, substantial improvement.

	#==> m.gkr_heti3 is the best model.	
		
mod=m.gkr_heti3
plot(mod)	# no sign of heteroscedasticity (Zuur et al. 2009)

# Normality of the residuals?
qqnorm(mod,abline=c(0,1)) # note that "non-normality is generally less of a concern than heteroscedasticity [in linear and generalized linear mixed models]" (Bolker 2015)

# Model results
anova(mod,type="marginal")
summary(mod)

# Removal of non-significant "Precip x Predictor variable" interactions ?
m.gkr_heti3ML=lme(Krats ~ Precip*Plant_biomass
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)	#ML		
anova(m.gkr_heti3ML,type="marginal")
	# no interaction to drop.
	
# Model used in the following initial SEM:
m.gkr_heti3REML=lme(Krats ~ Precip*Plant_biomass
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	#REML
mod=m.gkr_heti3REML

# Model assumptions:
plot(mod)					# Homoscedasticity of the residuals
qqnorm(mod,abline=c(0,1))	# Normality of the residuals








	###
	# Ant abundance model
	###

m.ant.lme = lme( Ants ~ Precip*Plant_biomass + Precip*Krats
	,random = ~ 1|plot_id
	,data=carrizo,method = "REML",control=lmecontr_foc)
mod=m.ant.lme

# Homosedasticity of the residuals? (Zuur et al. 2009)
plot(mod)

	# source of heteroscedasticity (if any)?
	par(mfrow=c(2,2))
	plot(carrizo$Precip,resid(mod),xlab="Precip",ylab="residuals")
	plot(carrizo$Plant_biomass,resid(mod),xlab="Plant_biomass",ylab="residuals")
	plot(carrizo$Krats,resid(mod),xlab="Krats",ylab="residuals")
	plot(carrizo$fyear,resid(mod),xlab="fyear",ylab="residuals")
			
		#Precip
		m.ant_hetf1=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varFixed(~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ant_hetp1=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ant_hetcp1=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varConstPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.ant_hete1=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varExp(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.ant.lme,m.ant_hetf1,m.ant_hetp1,m.ant_hetcp1,m.ant_hete1)
		# => AIC is not improved

		#Plant_biomass
		m.ant_hetf2=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varFixed(~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ant_hetp2=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varPower(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ant_hetcp2=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varConstPower(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.ant_hete2=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varExp(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.ant.lme,m.ant_hetf2,m.ant_hetp2,m.ant_hetcp2,m.ant_hete2)
		# => AIC is not improved

		#Krats
		m.ant_hetf3=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varFixed(~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ant_hetp3=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varPower(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ant_hetcp3=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varConstPower(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.ant_hete3=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varExp(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.ant.lme,m.ant_hetf3,m.ant_hetp3,m.ant_hetcp3,m.ant_hete3)
		# => AIC is not improved

		#fyear
		m.ant_heti4=lme(Ants ~ Precip*Plant_biomass + Precip*Krats
				,weights = varIdent(form =~ 1 | fyear)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.ant.lme,m.ant_heti4)
		# => AIC is not improved

	#==> m.ant.lme is the best model.
mod=m.ant.lme
plot(mod)	# no sign of heteroscedasticity (Zuur et al. 2009)

# Normality of the residuals?
qqnorm(mod,abline=c(0,1)) # note that "non-normality is generally less of a concern than heteroscedasticity [in linear and generalized linear mixed models]" (Bolker 2015)

# Model results	
anova(mod,type="marginal")
summary(mod)	

# Removal of non-significant "Precip x Predictor variable" interactions ?
m.ant.lmeML = lme( Ants ~ Precip*Plant_biomass + Precip*Krats
	,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)	# ML
anova(m.ant.lmeML,type="marginal")
	# removing Precip x Plant_biomass or Precip x Krats?
	m.ant_s1a=lme( Ants ~ Plant_biomass + Precip*Krats
		,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
	m.ant_s1b=lme( Ants ~ Precip*Plant_biomass + Krats
		,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
	anova(m.ant.lmeML,m.ant_s1a,m.ant_s1b)
	anova(m.ant_s1b,type="marginal")
		# removing Precip x Plant_biomass ?
		m.ant_s2=lme( Ants ~ Precip + Plant_biomass + Krats
			,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
		anova(m.ant_s1b,m.ant_s2)	# dAIC<2, simpler model is better
		anova(m.ant_s2,type="marginal")

# Model used in the following initial SEM:
m.ant_s2REML=lme( Ants ~ Precip + Plant_biomass + Krats
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	# REML
mod=m.ant_s2REML

# Model assumptions:
plot(mod)					# Homoscedasticity of the residuals
qqnorm(mod,abline=c(0,1))	# Normality of the residuals








	###
	# Beetle abundance model
	###

m.col.lme = lme( Beetles ~ Precip*Plant_biomass + Precip*Krats
	,random = ~ 1|plot_id
	,data=carrizo,method = "REML",control=lmecontr_foc)

mod=m.col.lme

# Homosedasticity of the residuals? (Zuur et al. 2009)
plot(mod)

	# source of heteroscedasticity (if any)?
	par(mfrow=c(2,2))
	plot(carrizo$Precip,resid(mod),xlab="Precip",ylab="residuals")
	plot(carrizo$Plant_biomass,resid(mod),xlab="Plant_biomass",ylab="residuals")
	plot(carrizo$Krats,resid(mod),xlab="Krats",ylab="residuals")
	plot(carrizo$fyear,resid(mod),xlab="fyear",ylab="residuals")
	
		#Precip
		m.col_hetf1=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varFixed(~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.col_hetp1=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.col_hetcp1=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varConstPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.col_hete1=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varExp(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.col.lme,m.col_hetf1,m.col_hetp1,m.col_hetcp1,m.col_hete1)
		# => AIC improves with hete1, varExp(form =~Precip):
		AIC(m.col.lme)-AIC(m.col_hete1) # dAIC=2.71 > 2, substantial improvement

		#Plant_biomass
		m.col_hetf2=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varFixed(~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.col_hetp2=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varPower(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.col_hetcp2=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varConstPower(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.col_hete2=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varExp(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.col.lme,m.col_hetf2,m.col_hetp2,m.col_hetcp2,m.col_hete2)
		# => AIC is not improved

		#Krats
		m.col_hetf3=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varFixed(~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.col_hetp3=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varPower(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.col_hetcp3=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varConstPower(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.col_hete3=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varExp(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.col.lme,m.col_hetf3,m.col_hetp3,m.col_hetcp3,m.col_hete3)
		# => AIC improves best with hete3, varExp(form =~Krats):
		AIC(m.col.lme)-AIC(m.col_hete3) # dAIC=1.9994 < 2, no substantial improvement

		#fyear
		m.col_heti4=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varIdent(form =~ 1 | fyear)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.col.lme,m.col_heti4)
		# => AIC improves with heti4, varIdent(form =~ 1 | fyear):
		AIC(m.col.lme)-AIC(m.col_heti4) # dAIC=10.65 > 2, substantial improvement
		
		# combinations of factor?
		m.col_hete1i4=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb(varExp(form =~Precip),varIdent(form =~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	tabAIC=anova(m.col.lme,m.col_hete1,m.col_heti4,m.col_hete1i4)
	tabAIC[order(tabAIC$AIC,decreasing=F),]
	
	# Best AIC2 model?
	 x model with no weight is:
		orimod=m.col.lme
	 x lowest AIC is:
		low1AICmod=m.col_heti4
		AIC(orimod)-AIC(low1AICmod) #=> dAIC = 10.65: low1AICmod model is better than orimod model
	 x 2nd lowest AIC model among the simpler ones remaining (i.e. lower df than lowest AIC model) is:
		lowsimpl2AICmod=m.col_hete1
		AIC(orimod)-AIC(lowsimpl2AICmod) #=> dAIC = 2.71
		AIC(lowsimpl2AICmod)-AIC(low1AICmod) #=> dAIC = 7.93 > 2, lower AIC model is better
	==> m.col_heti4 is the best model.

mod=m.col_heti4
plot(mod)	# no sign of heteroscedasticity (Zuur et al. 2009)

# Normality of the residuals?
qqnorm(mod,abline=c(0,1)) # note that "non-normality is generally less of a concern than heteroscedasticity [in linear and generalized linear mixed models]" (Bolker 2015)

# Model results	
anova(mod,type="marginal")
summary(mod)	

# Removal of non-significant "Precip x Predictor variable" interactions ?
m.col_heti4ML=lme(Beetles ~ Precip*Plant_biomass + Precip*Krats
	,weights = varIdent(form =~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)	#ML
anova(m.col_heti4ML,type="marginal")
	#=> removing Precip x Krats?
	m.col_s1a=lme(Beetles ~ Precip*Plant_biomass + Krats
		,weights = varIdent(form =~ 1 | fyear)
		,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
	anova(m.col_heti4ML,m.col_s1a)
	anova(m.col_s1a,type="marginal")
	# no interaction to drop.

# Model used in the following initial SEM:
m.col_s1aREML=lme(Beetles ~ Precip*Plant_biomass + Krats
		,weights = varIdent(form =~ 1 | fyear)
		,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	#REML
mod=m.col_s1aREML

# Model assumptions:
plot(mod)					# Homoscedasticity of the residuals
qqnorm(mod,abline=c(0,1))	# Normality of the residuals








	###
	# Orthopteran abundance model
	###

m.ort.lme = lme( Orthopterans ~ Precip*Plant_biomass + Precip*Krats
	,random = ~ 1|plot_id
	,data=carrizo,method = "REML",control=lmecontr_foc)
mod=m.ort.lme

# Homosedasticity of the residuals? (Zuur et al. 2009)
plot(mod)

	# source of heteroscedasticity (if any)?
	par(mfrow=c(2,2))
	plot(carrizo$Precip,resid(mod),xlab="Precip",ylab="residuals")
	plot(carrizo$Plant_biomass,resid(mod),xlab="Plant_biomass",ylab="residuals")
	plot(carrizo$Krats,resid(mod),xlab="Krats",ylab="residuals")
	plot(carrizo$fyear,resid(mod),xlab="fyear",ylab="residuals")

		#Precip
		m.ort_hetf1=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varFixed(~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hetp1=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hetcp1=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varConstPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.ort_hete1=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varExp(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.ort.lme,m.ort_hetf1,m.ort_hetp1,m.ort_hetcp1,m.ort_hete1)
		# => AIC improves with hete1, varExp(form =~Precip):
		AIC(m.ort.lme)-AIC(m.ort_hete1) # dAIC=5.19 > 2, substantial improvement.

		#Plant_biomass
		m.ort_hetf2=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varFixed(~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hetp2=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varPower(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hetcp2=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varConstPower(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.ort_hete2=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varExp(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.ort.lme,m.ort_hetf2,m.ort_hetp2,m.ort_hetcp2,m.ort_hete2)
		# => AIC improves best with hete2, varExp(form =~Plant_biomass):
		AIC(m.ort.lme)-AIC(m.ort_hete2) # dAIC=11.19 > 2, substantial improvement.

		#Krats
		m.ort_hetf3=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varFixed(~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hetp3=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varPower(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hetcp3=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varConstPower(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.ort_hete3=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varExp(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.ort.lme,m.ort_hetf3,m.ort_hetp3,m.ort_hetcp3,m.ort_hete3)
		# => AIC improves best with hete3, varExp(form =~Krats):
		AIC(m.ort.lme)-AIC(m.ort_hete3) # dAIC=9.27 > 2, substantial improvement.

		#Year as factor
		m.ort_heti4=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varIdent(form = ~ 1 | fyear)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.ort.lme,m.ort_heti4)
		# => AIC improves with heti4, varIdent(form = ~ 1 | fyear):
		AIC(m.ort.lme)-AIC(m.ort_heti4) # dAIC=132.29 > 2, substantial improvement.
		
		#Combinations:
		hete1, varExp(form =~Precip)
		hete2, varExp(form =~Plant_biomass)
		hete3, varExp(form =~Krats)
		heti4, varIdent(form =~ 1 | fyear)
	
		m.ort_hete1e2=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb(varExp(form =~Precip) , varExp(form =~Plant_biomass))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)		
		m.ort_hete1e3=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb(varExp(form =~Precip) ,varExp(form =~Krats) )
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hete1i4=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb(varExp(form =~Precip) , varIdent(form = ~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.ort_hete2e3=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb(varExp(form =~Plant_biomass) ,varExp(form =~Krats) )
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hete2i4=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb(varExp(form =~Plant_biomass) , varIdent(form = ~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.ort_hete3i4=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb(varExp(form =~Krats) , varIdent(form = ~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hete1e2e3=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb( varExp(form =~Precip), varExp(form =~Plant_biomass) ,varExp(form =~Krats) )
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hete1e2i4=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb( varExp(form =~Precip), varExp(form =~Plant_biomass) , varIdent(form =~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.ort_hete1e3i4=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb(varExp(form =~Precip) ,varExp(form =~Krats)  , varIdent(form =~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.ort_hete2e3i4=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb( varExp(form =~Plant_biomass),varExp(form =~Krats)  , varIdent(form =~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.ort_hete1e2e3i4=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
				,weights = varComb( varExp(form =~Precip),varExp(form =~Plant_biomass)  ,varExp(form =~Krats), varIdent(form =~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
				
	tabAIC=anova(m.ort.lme,m.ort_hete1,m.ort_hete2,m.ort_hete3,m.ort_heti4
			,m.ort_hete1e2,m.ort_hete1e3,m.ort_hete1i4,m.ort_hete2e3,m.ort_hete2i4,m.ort_hete3i4,m.ort_hete1e2e3,m.ort_hete1e2i4,m.ort_hete1e3i4,m.ort_hete2e3i4,m.ort_hete1e2e3i4
			)
	tabAIC[order(tabAIC$AIC,decreasing=F),]

	# Best AIC2 model?
	 x model with no weight is:
		orimod=m.ort.lme
	 x lowest AIC is:
		low1AICmod=m.ort_hete3i4
		AIC(orimod)-AIC(low1AICmod) #=> dAIC = 133.02 > 2: low1AICmod model is better than orimod model
	 x 2nd lowest AIC model among the simpler ones remaining (i.e. lower df than lowest AIC model) is:
		lowsimpl2AICmod=m.ort_heti4
		AIC(orimod)-AIC(lowsimpl2AICmod) #=> dAIC = 132.29
		AIC(lowsimpl2AICmod)-AIC(low1AICmod) #=> dAIC = 0.74 < 2, simpler model preferred
	 x 3rd lowest AIC model among the simpler ones remaining (i.e. lower df than lowest AIC model) is:
		lowsimpl3AICmod=m.ort_hete1e2e3
		AIC(orimod)-AIC(lowsimpl3AICmod) #=> dAIC = 21.82
		AIC(lowsimpl3AICmod)-AIC(lowsimpl2AICmod) #=> dAIC = 110.47 >2, lower AIC model preferred
	#==> m.ort_heti4 is the best model.

mod=m.ort_heti4
plot(mod)	# no sign of heteroscedasticity (Zuur et al. 2009)

# Normality of the residuals?
qqnorm(mod,abline=c(0,1)) # note that "non-normality is generally less of a concern than heteroscedasticity [in linear and generalized linear mixed models]" (Bolker 2015)

# Model results	
anova(mod,type="marginal")
summary(mod)	

# Removal of non-significant "Precip x Predictor variable" interactions ?
m.ort_heti4ML=lme(Orthopterans ~ Precip*Plant_biomass + Precip*Krats
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)	#ML		
Anova(m.ort_heti4ML,type="III")
anova(m.ort_heti4ML,type="marginal")
	# removing Precip x Plant_biomass ?
	m.ort_s1=lme(Orthopterans ~ Plant_biomass + Precip*Krats
		,weights = varIdent(form = ~ 1 | fyear)
		,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)	
	anova(m.ort_heti4ML,m.ort_s1)
	anova(m.ort_s1,type="marginal")
	# no interaction to drop.

# Model used in the following initial SEM:
m.ort_s1REML=lme(Orthopterans ~ Plant_biomass + Precip*Krats
		,weights = varIdent(form = ~ 1 | fyear)
		,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	#REML		
mod=m.ort_s1REML
# Model assumptions:
plot(mod)					# Homoscedasticity of the residuals
qqnorm(mod,abline=c(0,1))	# Normality of the residuals








	###
	# Squirrels model
	###			

m.sjas.lme = lme( Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
	,random = ~ 1|plot_id
	,data=carrizo,method = "REML",control=lmecontr_foc)
mod=m.sjas.lme

# Homosedasticity of the residuals? (Zuur et al. 2009)
plot(mod)

	# source of heteroscedasticity (if any)?
	par(mfrow=c(2,4))
	plot(carrizo$Precip,resid(mod),xlab="Precip",ylab="residuals")
	plot(carrizo$Plant_biomass,resid(mod),xlab="Plant_biomass",ylab="residuals")
	plot(carrizo$Krats,resid(mod),xlab="Krats",ylab="residuals")
	plot(carrizo$Ants,resid(mod),xlab="Ants",ylab="residuals")
	plot(carrizo$Beetles,resid(mod),xlab="Beetles",ylab="residuals")
	plot(carrizo$Orthopterans,resid(mod),xlab="Orthopterans",ylab="residuals")
	plot(carrizo$fyear,resid(mod),xlab="fyear",ylab="residuals")

		#Precip
		m.sjas_hetf1=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varFixed(~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetp1=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetcp1=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varConstPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.sjas_hete1=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varExp(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.sjas.lme,m.sjas_hetf1,m.sjas_hetp1,m.sjas_hetcp1,m.sjas_hete1)
		# => AIC improves with hetp1, varPower(form =~Precip):
		AIC(m.sjas.lme)-AIC(m.sjas_hetp1) # dAIC=0.49 < 2, no substantial improvement

		#Plant_biomass
		m.sjas_hetf2=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varFixed(~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetp2=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varPower(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetcp2=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varConstPower(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.sjas_hete2=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varExp(form =~Plant_biomass)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.sjas.lme,m.sjas_hetf2,m.sjas_hetp2,m.sjas_hetcp2,m.sjas_hete2)
		# => AIC is not improved.

		#Krats
		m.sjas_hetf3=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varFixed(~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetp3=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varPower(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetcp3=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varConstPower(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.sjas_hete3=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varExp(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.sjas.lme,m.sjas_hetf3,m.sjas_hetp3,m.sjas_hetcp3,m.sjas_hete3)
		# => AIC improves best with hete3, varExp(form =~Krats)
		AIC(m.sjas.lme)-AIC(m.sjas_hete3) # dAIC=1.51 < 2, no substantial improvement

		#Year as factor
		m.sjas_heti4=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varIdent(form = ~ 1 | fyear)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.sjas.lme,m.sjas_heti4)
		# => AIC improves with heti4, varIdent(form = ~ 1 | fyear):
		AIC(m.sjas.lme)-AIC(m.sjas_heti4) # dAIC=2.37 > 2, substantial improvement

		#Ants
		m.sjas_hetf5=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varFixed(~Ants)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetp5=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varPower(form =~Ants)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetcp5=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varConstPower(form =~Ants)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.sjas_hete5=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varExp(form =~Ants)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.sjas.lme,m.sjas_hete5)	# ,m.sjas_hetf5,m.sjas_hetp5,m.sjas_hetcp5 did not work
		# => AIC improves with hete5, varExp(form =~Ants):
		AIC(m.sjas.lme)-AIC(m.sjas_hete5) # dAIC=3.21 > 2, substantial improvement.

		#Beetles
		m.sjas_hetf6=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varFixed(~Beetles)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetp6=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varPower(form =~Beetles)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetcp6=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varConstPower(form =~Beetles)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.sjas_hete6=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varExp(form =~Beetles)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.sjas.lme,m.sjas_hetf6,m.sjas_hetp6,m.sjas_hetcp6,m.sjas_hete6)
		# => AIC improves best with hete6, varExp(form =~Beetles):
		AIC(m.sjas.lme)-AIC(m.sjas_hete6) # dAIC=3.79 > 2, substantial improvement.
		
		#Orthopterans
		m.sjas_hetf7=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varFixed(~Orthopterans)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetp7=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varPower(form =~Orthopterans)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_hetcp7=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varConstPower(form =~Orthopterans)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.sjas_hete7=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varExp(form =~Orthopterans)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.sjas.lme,m.sjas_hetf7,m.sjas_hetp7,m.sjas_hetcp7,m.sjas_hete7)
		# => AIC is not improved.
	
		#Combinations?
		heti4, varIdent(form = ~ 1 | fyear)
		hete5, varExp(form =~Ants)
		hete6, varExp(form =~Beetles)

		m.sjas_heti4e5=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varComb ( varIdent(form = ~ 1 | fyear) , varExp(form =~Ants) )
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.sjas_heti4e6=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varComb ( varIdent(form = ~ 1 | fyear) , varExp(form =~Beetles) )
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)										
		m.sjas_hete5e6=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varComb ( varExp(form =~Ants) , varExp(form =~Beetles) )
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.sjas_heti4e5e6=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
				,weights = varComb (  varIdent(form = ~ 1 | fyear) ,    varExp(form =~Ants) ,varExp(form =~Beetles) )
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		
		tabAIC=anova(m.sjas.lme,m.sjas_heti4,m.sjas_hete5,m.sjas_hete6
				,m.sjas_heti4e5
				,m.sjas_heti4e6
				,m.sjas_hete5e6
				,m.sjas_heti4e5e6
				)
		tabAIC[order(tabAIC$AIC,decreasing=F),]

	# Best AIC2 model?
	 x model with no weight is:
		orimod=m.sjas.lme
	 x lowest AIC is:
		low1AICmod=m.sjas_hete5e6
		AIC(orimod)-AIC(low1AICmod) #=> dAIC = 4.51 > 2: low1AICmod model is better than orimod model
	 x 2nd lowest AIC model among the simpler ones remaining (i.e. lower df than lowest AIC model) is:
		lowsimpl2AICmod=m.sjas_hete6
		AIC(orimod)-AIC(lowsimpl2AICmod) #=> dAIC = 3.79
		AIC(lowsimpl2AICmod)-AIC(low1AICmod) #=> dAIC = 0.71 <2, simpler model is preferred.
	 #==> m.sjas_hete6 is the best model.

mod=m.sjas_hete6
plot(mod)	# no sign of heteroscedasticity (Zuur et al. 2009)

# Normality of the residuals?
qqnorm(mod,abline=c(0,1)) # note that "non-normality is generally less of a concern than heteroscedasticity [in linear and generalized linear mixed models]" (Bolker 2015)

# Model results	
anova(mod,type="marginal")
summary(mod)	

# Removal of non-significant "Precip x Predictor variable" interactions ?
m.sjas_hete6ML=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
	,weights = varExp(form =~Beetles)
	,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)	#ML
anova(m.sjas_hete6ML,type="marginal")
	# removing Precip x Krats or Precip x Orthopterans ?
	m.sjas_s1a=lme(Squirrels ~ Precip*Plant_biomass + Krats + Precip*Ants + Precip*Beetles + Precip*Orthopterans
		,weights = varExp(form =~Beetles)
		,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
	m.sjas_s1b=lme(Squirrels ~ Precip*Plant_biomass + Precip*Krats + Precip*Ants + Precip*Beetles + Orthopterans
		,weights = varExp(form =~Beetles)
		,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
	anova(m.sjas_hete6ML,m.sjas_s1a,m.sjas_s1b)
	anova(m.sjas_s1b,type="marginal")
		# removing Precip x Krats ?
		m.sjas_s2=lme(Squirrels ~ Precip*Plant_biomass + Krats + Precip*Ants + Precip*Beetles + Orthopterans
			,weights = varExp(form =~Beetles)
			,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
		anova(m.sjas_s1b,m.sjas_s2)
		anova(m.sjas_s2,type="marginal")
	# no interaction to drop.

# Model used in the following initial SEM:
m.sjas_s2REML=lme(Squirrels ~ Precip*Plant_biomass + Krats + Precip*Ants + Precip*Beetles + Orthopterans
	,weights = varExp(form =~Beetles)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	#REML
mod=m.sjas_s2REML

# Model assumptions:
plot(mod)					# Homoscedasticity of the residuals
qqnorm(mod,abline=c(0,1))	# Normality of the residuals








	###
	# Lizards model
	###

m.uta.lme = lme( Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
	,random = ~ 1|plot_id
	,data=carrizo,method = "REML",control=lmecontr_foc)
mod=m.uta.lme

# Homosedasticity of the residuals? (Zuur et al. 2009)
plot(mod)

	# source of heteroscedasticity (if any)?
	par(mfrow=c(2,4))
	plot(carrizo$Precip,resid(mod),xlab="Precip",ylab="residuals")
	plot(carrizo$Krats,resid(mod),xlab="Krats",ylab="residuals")
	plot(carrizo$Ants,resid(mod),xlab="Ants",ylab="residuals")
	plot(carrizo$Beetles,resid(mod),xlab="Beetles",ylab="residuals")
	plot(carrizo$Orthopterans,resid(mod),xlab="Orthopterans",ylab="residuals")
	plot(carrizo$Squirrels,resid(mod),xlab="Squirrels",ylab="residuals")
	plot(carrizo$fyear,resid(mod),xlab="fyear",ylab="residuals")

		#Precip
		m.uta_hetf1=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varFixed(~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetp1=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetcp1=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varConstPower(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.uta_hete1=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varExp(form =~Precip)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.uta.lme,m.uta_hetf1,m.uta_hetp1,m.uta_hetcp1,m.uta_hete1)
		# => AIC improves with hete1, varExp(form =~Precip):
		AIC(m.uta.lme)-AIC(m.uta_hete1) # dAIC=0.23 < 2, no substantial improvement

		#Krats
		m.uta_hetf3=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varFixed(~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetp3=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varPower(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetcp3=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varConstPower(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.uta_hete3=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varExp(form =~Krats)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		anova(m.uta.lme,m.uta_hetf3,m.uta_hetp3,m.uta_hetcp3,m.uta_hete3)
		# => AIC and BIC model improve with hetp3, varPower(form =~Krats):
		AIC(m.uta.lme)-AIC(m.uta_hetp3) # dAIC=6.88 > 2, substantial improvements.

		#Year as factor
		m.uta_heti4=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varIdent(form = ~ 1 | fyear)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.uta.lme,m.uta_heti4)
		# => AIC improves with heti4, varIdent(form = ~ 1 | fyear):
		AIC(m.uta.lme)-AIC(m.uta_heti4) # dAIC=14.75 > 2, substantial improvements.
		
		#Ants
		m.uta_hetf5=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varFixed(~Ants)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetp5=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varPower(form =~Ants)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetcp5=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varConstPower(form =~Ants)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.uta_hete5=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varExp(form =~Ants)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.uta.lme,m.uta_hete5)	# ,m.uta_hetf5,m.uta_hetp5,m.uta_hetcp5 did not work
		# => AIC is not improved.		
		
		#Beetles
		m.uta_hetf6=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varFixed(~Beetles)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetp6=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varPower(form =~Beetles)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetcp6=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varConstPower(form =~Beetles)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.uta_hete6=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varExp(form =~Beetles)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.uta.lme,m.uta_hetf6,m.uta_hetp6,m.uta_hetcp6,m.uta_hete6)
		# => AIC improves with hete6, varExp(form =~Beetles):
		AIC(m.uta.lme)-AIC(m.uta_hete6) # dAIC=0.31 < 2, no substantial improvements.
		
		#Orthopterans
		m.uta_hetf7=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varFixed(~Orthopterans)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetp7=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varPower(form =~Orthopterans)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetcp7=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varConstPower(form =~Orthopterans)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.uta_hete7=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varExp(form =~Orthopterans)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.uta.lme,m.uta_hetf7,m.uta_hetp7,m.uta_hetcp7,m.uta_hete7)
		# => AIC improves best with hete7, varExp(form =~Orthopterans):
		AIC(m.uta.lme)-AIC(m.uta_hete7) # dAIC=3.69 > 2, substantial improvements.

		#Squirrels
		m.uta_hetf8=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varFixed(~Squirrels)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetp8=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varPower(form =~Squirrels)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetcp8=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varConstPower(form =~Squirrels)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)				
		m.uta_hete8=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varExp(form =~Squirrels)
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	anova(m.uta.lme,m.uta_hete8)	# ,m.uta_hetf8,m.uta_hetp8,m.uta_hetcp8 did not work
		# => AIC is not improved
		
		#Combinations:
		hetp3, varPower(form =~Krats)
		heti4, varIdent(form = ~ 1 | fyear)
		hete7, varExp(form =~Orthopterans)

		m.uta_hetp3i4=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetp3e7=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varComb( varPower(form =~Krats) ,  varExp(form =~Orthopterans))
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_heti4e7=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varComb( varIdent(form = ~ 1 | fyear) , varExp(form =~Orthopterans) )
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		m.uta_hetp3i4e7=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varComb( varPower(form =~Krats) ,varIdent(form = ~ 1 | fyear)  , varExp(form =~Orthopterans) )
				,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
				
		tabAIC=anova(m.uta.lme,m.uta_hetp3,m.uta_heti4,m.uta_hete7
				,m.uta_hetp3i4
				,m.uta_hetp3e7
				,m.uta_heti4e7
				,m.uta_hetp3i4e7
				)
	tabAIC[order(tabAIC$AIC,decreasing=F),]

	# Best AIC2 model?
	 x model with no weight is:
		orimod=m.uta.lme
	 x lowest AIC is:
		low1AICmod=m.uta_hetp3i4e7
		AIC(orimod)-AIC(low1AICmod) #=> dAIC = 23.01 > 2: low1AICmod model is better than orimod model
	 x 2nd lowest AIC model among the simpler ones remaining (i.e. lower df than lowest AIC model) is:
		lowsimpl2AICmod=m.uta_hetp3i4
		AIC(orimod)-AIC(lowsimpl2AICmod) #=> dAIC = 22.50 > 2
		AIC(lowsimpl2AICmod)-AIC(low1AICmod) #=> dAIC = 0.51 < 2, simpler model is better
	 x 3rd lowest AIC model among the simpler ones remaining (i.e. lower df than lowest AIC model) is:
		lowsimpl3AICmod=m.uta_hetp3e7
		AIC(orimod)-AIC(lowsimpl3AICmod) #=> dAIC = 15.06 > 2
		AIC(lowsimpl3AICmod)-AIC(lowsimpl2AICmod) #=> dAIC = 7.45 > 2, lower AIC model is better
	#==> m.uta_hetp3i4 is the best model.

mod=m.uta_hetp3i4
plot(mod)	# no sign of heteroscedasticity (Zuur et al. 2009)

# Normality of the residuals?
qqnorm(mod,abline=c(0,1)) # note that "non-normality is generally less of a concern than heteroscedasticity [in linear and generalized linear mixed models]" (Bolker 2015)

# Model results	
anova(mod,type="marginal")
summary(mod)	

# Removal of non-significant "Precip x Predictor variable" interactions ?
m.uta_hetp3i4ML=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
	,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
	,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc) #ML
anova(m.uta_hetp3i4ML,type="marginal")
	#=> removing Precip x Ants, Precip x Beetles, Precip x Orthopterans or Precip x Squirrels ?
	m.uta_s1a=lme(Lizards ~ Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
		,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
		,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
	m.uta_s1b=lme(Lizards ~ Precip*Ants + Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
		,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
		,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
	m.uta_s1c=lme(Lizards ~ Precip*Ants + Precip*Beetles + Orthopterans + Precip*Krats + Precip*Squirrels
		,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
		,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
	m.uta_s1d=lme(Lizards ~ Precip*Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Squirrels
		,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
		,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
	anova(m.uta_hetp3i4ML,m.uta_s1a,m.uta_s1b,m.uta_s1c,m.uta_s1d)
	anova(m.uta_s1a,type="marginal")
		#=> removing Precip x Beetles, Precip x Orthopterans or Precip x Squirrels ?
		m.uta_s2a=lme(Lizards ~ Ants + Beetles + Precip*Orthopterans + Precip*Krats + Precip*Squirrels
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
			,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
		m.uta_s2b=lme(Lizards ~ Ants + Precip*Beetles + Orthopterans + Precip*Krats + Precip*Squirrels
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
			,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
		m.uta_s2c=lme(Lizards ~ Ants + Precip*Beetles + Precip*Orthopterans + Precip*Krats + Squirrels
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
			,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
		anova(m.uta_s1a,m.uta_s2a,m.uta_s2b,m.uta_s2c)
		anova(m.uta_s2a,type="marginal")
			#=> removing Precip x Orthopterans or Precip x Squirrels?
			m.uta_s3a=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Precip*Squirrels
				,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
			m.uta_s3b=lme(Lizards ~ Ants + Beetles + Precip*Orthopterans + Precip*Krats + Squirrels
				,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
				,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
			anova(m.uta_s2a,m.uta_s3a,m.uta_s3b)
			anova(m.uta_s3b,type="marginal")
				#=> removing Precip x Orthopterans ?
				m.uta_s4=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels
					,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
					,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)
				anova(m.uta_s3b,m.uta_s4)	# dAIC<2, simpler model is better
				anova(m.uta_s4,type="marginal")
	# no interaction to drop.

# Model used in the following initial SEM:
m.uta_s4REML=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels
	,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
	,random = ~ 1|plot_id,data=carrizo,method = "ML",control=lmecontr_foc)	#REML
mod=m.uta_s4REML
# Model assumptions:
plot(mod)					# Homoscedasticity of the residuals
qqnorm(mod,abline=c(0,1))	# Normality of the residuals		








#####
#### C. Piecewise structural equation model
#####
	
# Set of sub-models, obtained from the preliminary analysis, to be used as the initial SEM:
	# Plant biomass model
m.prod_heti4REML=lme(Plant_biomass ~ Precip*local_cond_code
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Krats model
m.gkr_heti3REML=lme(Krats ~ Precip*Plant_biomass
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Ant abundance model
m.ant_s2REML=lme( Ants ~ Precip + Plant_biomass + Krats
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Beetle abundance model
m.col_s1aREML=lme(Beetles ~ Precip*Plant_biomass + Krats
	,weights = varIdent(form =~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Orthopteran abundance model
m.ort_s1REML=lme(Orthopterans ~ Plant_biomass + Precip*Krats
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Squirrels model
m.sjas_s2REML=lme(Squirrels ~ Precip*Plant_biomass + Krats + Precip*Ants + Precip*Beetles + Orthopterans
	,weights = varExp(form =~Beetles)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Lizards model
m.uta_s4REML=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels
	,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)

	
#
# C.1 Evaluation initial SEM
#

#list of models:
carrizo_initialSEM = list(
	# Plant biomass model
Plant_biomass = lme(Plant_biomass ~ Precip*local_cond_code
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	
	# Krats model
,Krats = lme(Krats ~ Precip*Plant_biomass
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Ant abundance model
,Ants = lme( Ants ~ Precip + Plant_biomass + Krats
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)		
	# Beetle abundance model
,Beetles = lme(Beetles ~ Precip*Plant_biomass + Krats
	,weights = varIdent(form =~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	
	# Orthopteran abundance model
,Orthopterans = lme(Orthopterans ~ Plant_biomass + Precip*Krats
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	
	# Squirrels model
,Squirrels = lme(Squirrels ~ Precip*Plant_biomass + Krats + Precip*Ants + Precip*Beetles + Orthopterans
	,weights = varExp(form =~Beetles)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)		
	# Lizards model
,Lizards = lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels
	,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
)	
	# model evaluation
modeval=carrizo_initialSEM
modfit=sem.fit(modeval, carrizo
		,corr.errors=c("Ants ~~ Beetles","Ants ~~ Orthopterans","Beetles ~~ Orthopterans")
		,model.control=list(rep(lmecontr_foc,7)))
modfit
		# Fish.C=241.21,pval<0.001, AIC=413.21
		# significant lack of fit between initial SEM and the data

	# Model complexity?
# sample size:
	modfit$AIC$n
# max nb of paths authorized (Grace et al. 2015):
	modfit$AIC$n/5
# nb of paths:
	unstdcoef=sem.coefs(modeval, carrizo ,corr.errors=c("Ants ~~ Beetles","Ants ~~ Orthopterans","Beetles ~~ Orthopterans")	,standardize = "none")
	nrow(unstdcoef[which(is.na(unstdcoef$std.error)==FALSE),])	# Nb of paths
# Ratio sample size/paths should be higher than 5:
	modfit$AIC$n/nrow(unstdcoef[which(is.na(unstdcoef$std.error)==FALSE),])>5


#
# C.2 Path addition
#

#list of models:
carrizo_addingpathsSEM = list(				# the models below are the results of 11 path additions ; see details below
	# Plant biomass model
Plant_biomass = lme(Plant_biomass ~ Precip*local_cond_code
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	
	# Krats model
,Krats = lme(Krats ~ Precip*Plant_biomass + Precip*local_cond_code
	,weights = varComb(varIdent(form = ~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	
	# Ant abundance model
,Ants = lme( Ants ~ Precip + Plant_biomass + Krats + local_cond_code
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)		
	# Beetle abundance model
,Beetles = lme(Beetles ~ Precip*Plant_biomass + Krats + Precip*local_cond_code
	,weights = varComb(varIdent(form =~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Orthopteran abundance model
,Orthopterans = lme(Orthopterans ~ Plant_biomass + Precip*Krats + local_cond_code
	,weights = varComb(varIdent(form = ~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Squirrels model
,Squirrels = lme(Squirrels ~ Precip*Plant_biomass + Krats + Precip*Ants + Precip*Beetles + Orthopterans + Precip*local_cond_code
	,weights = varExp(form =~Beetles)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)		
	# Lizards model
,Lizards = lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels + Precip*local_cond_code + Plant_biomass
	,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
)
	# model evaluation
modeval=carrizo_addingpathsSEM
modfit=sem.fit(modeval, carrizo
		,corr.errors=c("Ants ~~ Beetles","Ants ~~ Orthopterans","Beetles ~~ Orthopterans")
		,model.control=list(rep(lmecontr_foc,7)))
modfit
	# Missing paths?
modmis=sem.missing.paths(modeval, carrizo,corr.errors=c("Ants ~~ Beetles","Ants ~~ Orthopterans","Beetles ~~ Orthopterans"),model.control=list(rep(lmecontr_foc,7)))
modmis[order(modmis$p.value,decreasing=F),]

	# Addition of (lowest pvalue first; path creating cycle not added, path with stronger estimate chosen when equal pvalues):
		# none: 											Fish.C=241.21, pval=0.000, AIC=413.21
		# first path proposed is "Beetles ~ Precip*local_cond_code". This first requires adding "Beetles ~ local_cond_code" (also significantly suggested)
			# "Beetles ~ local_cond_code"					Fish.C=212.01, pval=0.000
			# "Beetles ~ Precip*local_cond_code":			Fish.C=182.28, pval=0.000
		# "Lizards ~ local_cond_code":						Fish.C=153.65, pval=0.000
		# "Ants ~ local_cond_code":							Fish.C=128.63 ,pval=0.000
		# path proposed is "Krats ~ Precip*local_cond_code". This first requires adding "Krats ~ local_cond_code" (also significantly suggested)
			# "Krats ~ local_cond_code"						Fish.C=127.36 ,pval=0.000
			# "Krats ~ Precip*local_cond_code":				Fish.C=95.46  ,pval=0.000
		# "Squirrels ~ local_cond_code":					Fish.C=75.84  ,pval=0.000
		# "Lizards ~ Plant_biomass":						Fish.C=59.02  ,pval=0.009
		# "Lizards ~ local_cond_code*Precip":				Fish.C=46.06  ,pval=0.081
		# "Squirrels ~ local_cond_code*Precip":				Fish.C=35.97  ,pval=0.288
		# "Orthopterans ~ local_cond_code":					Fish.C=25.39  ,pval=0.706
	#=> no more paths to add.
	#=> No significant lack of fit between the SEM and the data; 

	# is sample size enough for model complexity? Ratio sample size/paths should be higher than 5:
	unstdcoef=sem.coefs(modeval, carrizo, corr.errors=c("Ants ~~ Beetles","Ants ~~ Orthopterans","Beetles ~~ Orthopterans"), standardize = "none")
	modfit$AIC$n/nrow(unstdcoef[which(is.na(unstdcoef$std.error)==FALSE),])>=5
		nrow(unstdcoef[which(is.na(unstdcoef$std.error)==FALSE),])	# Two paths too much; max nb of paths = modfit$AIC$n/5
		
	# The 99 lines below were used to check for the need for updating structure of the variance, each time a new path was suggested:
		#"Beetles ~ local_cond_code":
		modsam.col=lme(Beetles ~ Precip*Plant_biomass + Krats + local_cond_code
			,weights = varIdent(form =~ 1 | fyear)
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		modsam.col.iX=lme(Beetles ~ Precip*Plant_biomass + Krats + local_cond_code
			,weights = varComb(varIdent(form =~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		anova(modsam.col,modsam.col.iX)
			#=> AIC improves with hetiX, varIdent(form =~ 1 | local_cond_code)
			AIC(modsam.col)-AIC(modsam.col.iX) # dAIC=36.90 >2, updating weights is required.

		#"Beetles ~ Precip*local_cond_code":
		modsam.col=lme(Beetles ~ Precip*Plant_biomass + Krats + Precip*local_cond_code
			,weights = varComb(varIdent(form =~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)			
			
		#"Lizards ~ local_cond_code":
		modsam.uta=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels + local_cond_code
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		modsam.uta.iX=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels + local_cond_code
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		anova(modsam.uta,modsam.uta.iX)
			#=> AIC is not improved.
				
		#"Ants ~ local_cond_code":
		modsam.ant=lme( Ants ~ Precip + Plant_biomass + Krats + local_cond_code
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		modsam.ant.iX=lme( Ants ~ Precip + Plant_biomass + Krats + local_cond_code
			,weights = varIdent(form =~ 1 | local_cond_code)
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		anova(modsam.ant,modsam.ant.iX)
			#=> AIC is not improved.

		#"Krats ~ local_cond_code":
		modsam.gkr=lme(Krats ~ Precip*Plant_biomass + local_cond_code
			,weights = varIdent(form = ~ 1 | fyear)
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		modsam.gkr.iX=lme(Krats ~ Precip*Plant_biomass + local_cond_code
			,weights = varComb(varIdent(form = ~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		anova(modsam.gkr,modsam.gkr.iX)
			#=> AIC improves with hetiX, varIdent(form =~ 1 | local_cond_code)
			AIC(modsam.gkr)-AIC(modsam.gkr.iX) # dAIC=2.46 > 2, updating weights is required.
			
		#"Krats ~ Precip*local_cond_code":
		modsam.gkr=lme(Krats ~ Precip*Plant_biomass + Precip*local_cond_code
			,weights = varComb(varIdent(form = ~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	
		#"Squirrels ~ local_cond_code":
		modsam.sjas=lme(Squirrels ~ Precip*Plant_biomass + Krats + Precip*Ants + Precip*Beetles + Orthopterans + local_cond_code
			,weights = varExp(form =~Beetles)
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		modsam.sjas.iX=lme(Squirrels ~ Precip*Plant_biomass + Krats + Precip*Ants + Precip*Beetles + Orthopterans + local_cond_code
			,weights = varComb(varExp(form =~Beetles),varIdent(form =~ 1 | local_cond_code))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		anova(modsam.sjas,modsam.sjas.iX)
			#=> AIC is not improved.

		#"Lizards ~ Plant_biomass":
		modsam.uta=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels + local_cond_code + Plant_biomass
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		modsam.uta.fX=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels + local_cond_code + Plant_biomass
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear), varFixed(~Plant_biomass))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		modsam.uta.pX=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels + local_cond_code + Plant_biomass
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear), varPower(form=~Plant_biomass))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		modsam.uta.cpX=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels + local_cond_code + Plant_biomass
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear), varConstPower(form=~Plant_biomass))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		modsam.uta.eX=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels + local_cond_code + Plant_biomass
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear),varExp(form=~Plant_biomass))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		anova(modsam.uta,modsam.uta.fX,modsam.uta.pX,modsam.uta.eX)	# ,modsam.uta.cpX did not work
			#=> AIC is not improved.
	
		#"Lizards ~ Precip*local_cond_code":
		modsam.uta=lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels + Precip*local_cond_code + Plant_biomass
			,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)

		#"Squirrels ~ Precip*local_cond_code":
		modsam.sjas=lme(Squirrels ~ Precip*Plant_biomass + Krats + Precip*Ants + Precip*Beetles + Orthopterans + Precip*local_cond_code
			,weights = varExp(form =~Beetles)
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
			
		#"Orthopterans ~ local_cond_code":
		modsam.ort=lme(Orthopterans ~ Plant_biomass + Precip*Krats + local_cond_code
			,weights = varIdent(form = ~ 1 | fyear)
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		modsam.ort.iX=lme(Orthopterans ~ Plant_biomass + Precip*Krats + local_cond_code
			,weights = varComb(varIdent(form = ~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
			,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
		anova(modsam.ort,modsam.ort.iX)
			#=> AIC improves with hetiX, varIdent(form =~ 1 | local_cond_code)
			AIC(modsam.ort)-AIC(modsam.ort.iX) # dAIC=7.72 > 2, updating weights is required.

	
#
# C.3 Path deletion
#

#list of models:
carrizo_delpathsSEM =  list(
	# Plant biomass model
Plant_biomass = lme(Plant_biomass ~ Precip*local_cond_code
	,weights = varIdent(form = ~ 1 | fyear)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	
	# Krats model
,Krats = lme(Krats ~ Plant_biomass + Precip*local_cond_code
		# removed paths: Precip x Plant_biomass
	,weights = varComb(varIdent(form = ~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)	
	# Ant abundance model
,Ants = lme( Ants ~ Plant_biomass + Krats + local_cond_code
		# removed paths: Precip
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)		
	# Beetle abundance model
,Beetles = lme(Beetles ~ Plant_biomass + Krats + Precip*local_cond_code
		# removed paths: Precip x Plant_biomass
	,weights = varComb(varIdent(form =~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Orthopteran abundance model
,Orthopterans = lme(Orthopterans ~ Plant_biomass + Precip*Krats + local_cond_code
	,weights = varComb(varIdent(form = ~ 1 | fyear),varIdent(form =~ 1 | local_cond_code))
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
	# Squirrels model
,Squirrels = lme(Squirrels ~ Plant_biomass + Krats + Precip*Ants + Precip*Beetles + Orthopterans + Precip*local_cond_code
		# removed paths: Precip x Plant_biomass
	,weights = varExp(form =~Beetles)
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)		
	# Lizards model
,Lizards = lme(Lizards ~ Ants + Beetles + Orthopterans + Precip*Krats + Squirrels + Precip*local_cond_code + Plant_biomass
	,weights = varComb( varPower(form =~Krats) ,  varIdent(form = ~ 1 | fyear))
	,random = ~ 1|plot_id,data=carrizo,method = "REML",control=lmecontr_foc)
)
	# model evaluation
modeval=carrizo_delpathsSEM
modfit=sem.fit(modeval, carrizo
		,corr.errors=c("Ants ~~ Beetles","Ants ~~ Orthopterans","Beetles ~~ Orthopterans")
		,model.control=list(rep(lmecontr_foc,7)))
modfit
	# is sample size enough for model complexity? Ratio sample size/paths should be higher than 5:
	unstdcoef=sem.coefs(modeval, carrizo,corr.errors=c("Ants ~~ Beetles","Ants ~~ Orthopterans","Beetles ~~ Orthopterans"),standardize = "none")
	nrow(unstdcoef[which(is.na(unstdcoef$std.error)==FALSE),])	# nb of paths
	modfit$AIC$n/5	# max nb of paths authorized

	# Deletion of non-significant paths:
		# highest pvalue first
		# main effects involved in two-way interactions are not removed
		# non-significant "Precip x Predictor variable" two-way interactions are removed to improve interpretability of the main effect.
	unstdcoef[order(unstdcoef$p.value,decreasing=T),]
		# none												Fish.C=25.39 ,pval=0.706, AIC=225.39
		# "Squirrels ~ Krats"								Fish.C=42.83 ,pval=0.096, AIC=240.83 => dAIC>2,	lowerAIC model is preferred. Path IS NOT removed
		# "Orthopterans ~ Plant_biomass"					Fish.C=33.81 ,pval=0.380, AIC=231.81 => dAIC>2,	lowerAIC model is preferred. Path IS NOT removed		
	# "Ants ~ Precip"									Fish.C=26.24 ,pval=0.753, AIC=224.24 => dAIC<2, simpler & lowerAIC model is preferred. Path IS removed
		# "Lizards ~ Squirrels"								Fish.C=40.34 ,pval=0.210, AIC=236.34 => dAIC>2,	lowerAIC model is preferred. Path IS NOT removed
	# "Krats ~ Precip*Plant_biomass"					Fish.C=25.84 ,pval=0.841, AIC=221.84 => dAIC<2, simpler & lowerAIC model is preferred. Path IS removed		
		# "Lizards ~ Plant_biomass"							Fish.C=44.53 ,pval=0.156, AIC=238.53 => dAIC>2,	lowerAIC model is preferred. Path IS NOT removed
	# "Squirrels ~ Precip*Plant_biomass"				Fish.C=29.91 ,pval=0.753, AIC=223.91 => dAIC=2.07; n.s. two-way interaction IS removed to improve interpretation of main effect.
		# "Squirrels ~ Plant_biomass"						Fish.C=46.18 ,pval=0.170, AIC=218.94 => dAIC>2,	lower AIC model preferred. Path IS NOT removed.
		# "Squirrels ~ Orthopterans"						Fish.C=44.66 ,pval=0.212, AIC=236.66 => dAIC>2,	lowerAIC model is preferred. Path IS NOT removed
		# "Lizards ~ Orthopterans"							Fish.C=43.23 ,pval=0.258, AIC=235.23 => dAIC>2, lowerAIC model is preferred. Path IS NOT removed
		# "Krats ~ Plant_biomass"							Fish.C=57.97 ,pval=0.020, AIC=249.97 => dAIC>2,	lowerAIC model is preferred. Path IS NOT removed
	# "Beetles ~ Precip*Plant_biomass"					Fish.C=19.55 ,pval=0.812, AIC=211.55 => dAIC<2,	simpler & lowerAIC model is preferred. Path IS removed
		# "Lizards ~ Ants"									Fish.C=38.83 ,pval=0.083, AIC=228.84 => dAIC>2, lowerAIC model is preferred. Path IS NOT removed
	#=> no other non-significant paths remain to be tested.
	#=> 40 paths
	#=> This is the final SEM: it fits adequately the data.


# A table of the path coefficients:
	unstdcoef=sem.coefs(modeval, carrizo,corr.errors=c("Ants ~~ Beetles","Ants ~~ Orthopterans","Beetles ~~ Orthopterans"),standardize = "none")
		# Unstandardized paths coefficients	
	scalcoef=sem.coefs(modeval, carrizo,corr.errors=c("Ants ~~ Beetles","Ants ~~ Orthopterans","Beetles ~~ Orthopterans"),standardize = "scale")
		# Standardized paths coefficients (i.e. after z-transformation of the data)
	ounstdcoef=unstdcoef[order(unstdcoef$response,unstdcoef$predictor),]
	oscalcoef=scalcoef[order(scalcoef$response,scalcoef$predictor),]
	colnames(ounstdcoef)[3:5]=paste("unstd",colnames(ounstdcoef)[3:5],sep="_")
	colnames(oscalcoef)[3:5]=paste("scal",colnames(oscalcoef)[3:5],sep="_")

	mod_coefs=cbind(ounstdcoef,oscalcoef[,3:5])
	mod_coefs

	# R² for individual models
sem.model.fits(modeval)
				#					Class	Family		Link		N  	Marginal	Conditional
				#	Plant_biomass   lme		gaussian	identity	210 0.2739314   0.3577692
				#	Krats			lme		gaussian	identity	210 0.3745781   0.6834798
				#	Ants			lme		gaussian	identity	210 0.2666920   0.3665653
				#	Beetles			lme		gaussian	identity	210 0.3973810   0.4492278
				#	Orthopterans	lme		gaussian	identity	210 0.1193517   0.1224732
				#	Squirrels		lme		gaussian	identity	210 0.2309355   0.3114217
				#	Lizards			lme		gaussian	identity	210 0.5799765   0.9185089

#
# C.4 Calculating the indirect, plant-mediated, effects of precipitation 
#

## Direct precipation effects on each variable
predictor.variable="Precip"
	mod_coefs[which(mod_coefs$predictor=="Precip"),c("response","predictor","scal_estimate")]	# direct effects
	
## Indirect, plant-mediated, precipitation effects on each variable
	# Plant_biomass: no indirect effects.

	# Indirect effect on Krats composed of:
	mod_coefs[which(mod_coefs$response=="Krats"),c("response","predictor","scal_estimate")]
		#=>Indirect effect on Krats composed of:
			#Precipitation to Plant_biomass to Krats
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]

	# Indirect effect on Ants composed of:
	mod_coefs[which(mod_coefs$response=="Ants"),c("response","predictor","scal_estimate")]
		#=>Indirect effect on Ants composed of:
			#Precipitation to Plant_biomass to Ants
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Ants"),"scal_estimate"]	+
			#Precipitation to Plant_biomass to Krats to Ants
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Ants"),"scal_estimate"]

	# Indirect effect on Beetles composed of:
	mod_coefs[which(mod_coefs$response=="Beetles"),c("response","predictor","scal_estimate")]
		#=>Indirect effect on Beetles composed of:
			#Precipitation to Plant_biomass to Beetles
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Beetles"),"scal_estimate"]	+
			#Precipitation to Plant_biomass to Krats to Beetles
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Beetles"),"scal_estimate"]
	
	# Indirect effect on Orthopterans composed of:
	mod_coefs[which(mod_coefs$response=="Orthopterans"),c("response","predictor","scal_estimate")]
		#=>Indirect effect on Orthopterans composed of:
			#Precipitation to Plant_biomass to Orthopterans
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Orthopterans"),"scal_estimate"]	+
			#Precipitation to Plant_biomass to Krats to Orthopterans
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Orthopterans"),"scal_estimate"]	
	
	# Indirect effect on Squirrels composed of:
	mod_coefs[which(mod_coefs$response=="Squirrels"),c("response","predictor","scal_estimate")]
		#=>Indirect effect on Krats composed of:
			#Precipitation to Plant_biomass to Squirrels
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Squirrels"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Krats to Squirrels
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Squirrels"),"scal_estimate"]	+			
			
			#Precipitation to Plant_biomass to Ants to Squirrels
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Ants"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Ants" & mod_coefs$response=="Squirrels"),"scal_estimate"]	+
	
			#Precipitation to Plant_biomass to Beetles to Squirrels
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Beetles"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Beetles" & mod_coefs$response=="Squirrels"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Orthopterans to Squirrels
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Orthopterans"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Orthopterans" & mod_coefs$response=="Squirrels"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Krats to Ants to Squirrels
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Ants"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Ants" & mod_coefs$response=="Squirrels"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Krats to Beetles to Squirrels
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Beetles"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Beetles" & mod_coefs$response=="Squirrels"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Krats to Orthopterans to Squirrels
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Orthopterans"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Orthopterans" & mod_coefs$response=="Squirrels"),"scal_estimate"]
			
	# Indirect effect on Lizards composed of:
	mod_coefs[which(mod_coefs$response=="Lizards"),c("response","predictor","scal_estimate")]
		#=>Indirect effect on Krats composed of:
			#Precipitation to Plant_biomass to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Lizards"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Krats to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Lizards"),"scal_estimate"]	+
			
			#Precipitation to Plant_biomass to Ants to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Ants"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Ants" & mod_coefs$response=="Lizards"),"scal_estimate"]	+
			
			#Precipitation to Plant_biomass to Beetles to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Beetles"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Beetles" & mod_coefs$response=="Lizards"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Orthopterans to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Orthopterans"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Orthopterans" & mod_coefs$response=="Lizards"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Squirrels to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Squirrels"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Squirrels" & mod_coefs$response=="Lizards"),"scal_estimate"]	+
			
			#Precipitation to Plant_biomass to Krats to Ants to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Ants"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Ants" & mod_coefs$response=="Lizards"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Krats to Beetles to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Beetles"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Beetles" & mod_coefs$response=="Lizards"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Krats to Orthopterans to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Krats"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Krats" & mod_coefs$response=="Orthopterans"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Orthopterans" & mod_coefs$response=="Lizards"),"scal_estimate"] +

			#Precipitation to Plant_biomass to Ants to Squirrels to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Ants"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Ants" & mod_coefs$response=="Squirrels"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Squirrels" & mod_coefs$response=="Lizards"),"scal_estimate"]	+
			
			#Precipitation to Plant_biomass to Beetles to Squirrels to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Beetles"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Beetles" & mod_coefs$response=="Squirrels"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Squirrels" & mod_coefs$response=="Lizards"),"scal_estimate"]	+

			#Precipitation to Plant_biomass to Orthopterans to Squirrels to Lizards
			mod_coefs[which(mod_coefs$predictor=="Precip" & mod_coefs$response=="Plant_biomass"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Plant_biomass" & mod_coefs$response=="Orthopterans"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Orthopterans" & mod_coefs$response=="Squirrels"),"scal_estimate"]*
			mod_coefs[which(mod_coefs$predictor=="Squirrels" & mod_coefs$response=="Lizards"),"scal_estimate"]
			
	

#####
#### D. Linear regressions of annual standardized coefficients against precipitation
#####

# Load outputs from AMOS multi-group analysis:
amos.outputs=read.csv("Deguines_etal_JAE_linearregressions.csv")
head(amos.outputs)
	# "stdcoeff": standardized path coefficient estimated with the multi-group analysis in Amos
	# "stdcoeff_low" and "stdcoeff_low": lower and upper bounds of the 95% bias-corrected confidence interval around stdcoeff.

unique(as.character(amos.outputs$path))	# these are the interactions the final SEM found to be significantly modified by precipitation
	
# Are standardized paths associated with precipitation levels?
allyears=sort(unique(amos.outputs$year))
allpaths=unique(amos.outputs$path)

xvar="Precip"
yvar="stdcoeff"
yvar_low="stdcoeff_low"
yvar_up="stdcoeff_up"

# The loop below makes the linear regressions for Figure 3 and Figure 4, and save each panel individually:
for(p in 1:length(allpaths))
	{
	P=allpaths[p]
	P
	amos.outputs_P=amos.outputs[which(amos.outputs$path==P),]
	amos.outputs_P

	# Figure
	plot(amos.outputs_P[,xvar],amos.outputs_P[,yvar],xlab=xvar,ylab=yvar,main=P,xlim=c(5,max(amos.outputs[,xvar])),ylim=c(-1.3,1.3),pch=16,axes=F)	
	axis(side=2,las=1,tck = 0.02)
	axis(side=1,las=1,tck = 0.02)	
	abline(h=0,lty=3,col="gray70")
	errbar(amos.outputs_P[,xvar],amos.outputs_P[,yvar],yminus=amos.outputs_P[,yvar_low],yplus=amos.outputs_P[,yvar_up]
	,xlab="",ylab="",add=TRUE)

	# Linear regression
	esti_range=amos.outputs_P[,yvar_up]-amos.outputs_P[,yvar_low]
	mod=lm(amos.outputs_P[,yvar] ~ amos.outputs_P[,xvar],weights=1/esti_range)
	anomod=Anova(mod,type="III")

	mtext(paste(
		"F val=",round(anomod["amos.outputs_P[, xvar]","F value"],3)
		,", pval=",round(anomod["amos.outputs_P[, xvar]","Pr(>F)"],3)
		,", R²=",round(r.squaredGLMM(mod)["R2c"],2),	sep="")		
		, side = 3, line = 0, outer = FALSE,cex = 1.0)
		
	mydata<- data.frame("amos.outputs_P[,xvar]"=amos.outputs_P[,xvar])
	prec_esti<- predict(mod, mydata,type="response")
	lines(amos.outputs_P[,xvar],prec_esti,lwd=2,type="l",lty=ifelse(anomod["amos.outputs_P[, xvar]","Pr(>F)"]<0.05,1,2),col="black")
	
	plotname=paste(P," (",yvar,") vs ",xvar,", years ",paste(allyears[1],"-",allyears[length(allyears)],sep=""),".bmp",sep="")
	plotname
	savePlot(filename=plotname,type="bmp")
	}
	
	
#########################################################################################################################################
### End of script for analyses published in:																						  ###
### Deguines N, Brashares JS, Prugh LR - J. of Animal Ecology - Precipitation alters interactions in a grassland ecological community ###
#########################################################################################################################################