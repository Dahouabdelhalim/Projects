### load packages
library(MCMCglmm)  #loads coda too
library(parallel)
library(lattice)
library(tictoc)

### load data and pedigree
load("cebus.RData")

### ENSO phases
cebus$ENSO_ <- 'Average/Neutral'
cebus$ENSO_[cebus$ENSO<=-0.5] <- 'Cool/La_Niña'
cebus$ENSO_[cebus$ENSO>=0.5] <- 'Warm/El_Niño'

### create sine and cosine of month to capture annual seasonality
require(plyr)
cebus$rad_12 <- mapvalues(
	cebus$Month, 
	from=c("01/Jan", "02/Feb", "03/Mar", "04/Apr", "05/May",
				 "06/Jun", "07/Jul", "08/Aug", "09/Sep", "10/Oct",
				 "11/Nov", "12/Dec"), 
	to=seq(12)*(2*pi/12))
cebus$rad_12 <- as.numeric(cebus$rad_12)

### format factor variables (data)
cols <- c("id", "animal", "Mother", "Group", "Sex", "Year", "Month", "ENSO_")
cebus[cols] <- lapply(cebus[cols], factor)

### format factor variables (pedigree)
cols <- c("animal", "Mother", "Father")
Ped[cols] <- lapply(Ped[cols], factor)

rm(cols)

###############################################################################
###############################################################################
###############################################################################

### NOTE: to set the prior
# one G structure needs to be specified for each random effect
# one R structure (for residuals) needs to be specified; 
### univariate models have only one

### Parameter expanded prior ###
### taken from Tutorial by Pierre de Villemereuil
### Estimation of a biological trait Heritability using the animal model 
#### and MCMCglmm (version 2), 2021-09-22

prior_m9a <- list(
	R=list(V=1, nu=0.002),  
	G=list(
		G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),  
		G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G8=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G9=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

## weakly informative inverse-Wishart prior (V=1, nu=0.002)
## performs poorly if variance is close to zero
## in case better to use a stronger prior (increase nu)
# lower values for nu create fatter tails 
## less conservative, will sample a broader parameter space
## danger of weird sampling issues at the tails
prior_m9b <- list(R=list(V=1, nu=0.002),  
									G=list(G1=list(V=1, nu=0.002),  
												 G2=list(V=1, nu=0.002),
												 G3=list(V=1, nu=0.002),
												 G4=list(V=1, nu=0.002),
												 G5=list(V=1, nu=0.002),
												 G6=list(V=1, nu=0.002),
												 G7=list(V=1, nu=0.002),
												 G8=list(V=1, nu=0.002),
												 G9=list(V=1, nu=0.002)))
prior_m9c <- list(R=list(V=1, nu=0.02),  
									G=list(G1=list(V=1, nu=0.02),  
												 G2=list(V=1, nu=0.02),
												 G3=list(V=1, nu=0.02),
												 G4=list(V=1, nu=0.02),
												 G5=list(V=1, nu=0.02),
												 G6=list(V=1, nu=0.02),
												 G7=list(V=1, nu=0.02),
												 G8=list(V=1, nu=0.02),
												 G9=list(V=1, nu=0.02)))

## Priors for poisson models ##
###############################

k <- 14 # number of fixed effects plus intercept

prior_p9a <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),  
		G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G8=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G9=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
# set strong prior for log of sampling effort
# Note: k-3 is because log(n) is 3rd to last
prior_p9a$B$mu[k-3] <- 1 
prior_p9a$B$V[k-3,k-3] <- 1e-7

prior_p9b <- list(B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
									R=list(V=1, nu=0.002),  
									G=list(
										G1=list(V=1, nu=0.002),  
										G2=list(V=1, nu=0.002),
										G3=list(V=1, nu=0.002),
										G4=list(V=1, nu=0.002),
										G5=list(V=1, nu=0.002),
										G6=list(V=1, nu=0.002),
										G7=list(V=1, nu=0.002),
										G8=list(V=1, nu=0.002),
										G9=list(V=1, nu=0.002)))
# set strong prior for log of sampling effort
prior_p9b$B$mu[k-3] <- 1 
prior_p9b$B$V[k-3,k-3] <- 1e-7

prior_p9c <- list(B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
									R=list(V=1, nu=0.02),  
									G=list(
										G1=list(V=1, nu=0.02),  
										G2=list(V=1, nu=0.02),
										G3=list(V=1, nu=0.02),
										G4=list(V=1, nu=0.02),
										G5=list(V=1, nu=0.02),
										G6=list(V=1, nu=0.02),
										G7=list(V=1, nu=0.02),
										G8=list(V=1, nu=0.02),
										G9=list(V=1, nu=0.02)))
# set strong prior for log of sampling effort
prior_p9c$B$mu[k-3] <- 1 
prior_p9c$B$V[k-3,k-3] <- 1e-7

###############################################################################
###############################################################################
###############################################################################

### multinomial2 models ###

set.seed(1977)  #Star Wars :) 
tic() 
m09a <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12) + 
			ENSO_, 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 1000, burnin = 100000, nitt = 4100000,
		prior = prior_m9a)
}, mc.cores=4)
toc()  #369617.741 sec elapsed

### multinomial2 model ###

set.seed(1977)  #Star Wars :) 
tic() 
m09b <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12) + 
			ENSO_, 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000, 
		prior = prior_m9b)
}, mc.cores=4)
toc()  #82709.315 sec elapsed

set.seed(1977)  #Star Wars :) 
tic() 
m09c <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12) + 
			ENSO_, 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000, 
		prior = prior_m9c)
}, mc.cores=4)
toc()  #84670.983 sec elapsed

### poisson models ###

set.seed(2001)  #: A Space Odyssey :) 
tic() 
p09a <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12) + 
			ENSO_ + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 500, burnin = 50000, nitt = 2050000, 
		prior = prior_p9a)
}, mc.cores=4)
toc()  #152544.88 sec elapsed

set.seed(2001)  #: A Space Odyssey :) 
tic() 
p09b <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12) + 
			ENSO_ + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9b)
}, mc.cores=4)
toc()  #45395.58 sec elapsed

set.seed(2001)  #: A Space Odyssey :) 
tic() 
p09c <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12) + 
			ENSO_ + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9c)
}, mc.cores=4)
toc()  #44983.867 sec elapsed

###############################################################################
###############################################################################
###############################################################################

### model diagnostics ###
## visual inspection ##

# you can load in the previously run models
load("Models1_fullModels.RData")

#####

bindChains <- function(model, pars){
	cbind(as.mcmc.list(model[[1]][[pars]]),
				as.mcmc.list(model[[2]][[pars]]),
				as.mcmc.list(model[[3]][[pars]]),
				as.mcmc.list(model[[4]][[pars]]))
}

chains.m.vcv <- bindChains(m09a, "VCV")
chains.m.sol <- bindChains(m09a, "Sol")
chains.p.vcv <- bindChains(p09a, "VCV")
chains.p.sol <- bindChains(p09a, "Sol")

## function for plotting chains together, trace plot next to density plot
color_scheme_set("viridis")
genPlt1 <- function(chain, component, ylab2=NULL){
	p_trace <- mcmc_trace(chain, pars = component) + 
		theme(legend.position = "none")
	if(is.null(ylab2) != TRUE) {
		p_trace <- p_trace + ylab(ylab2)
	}
	p_pdens <- mcmc_dens_overlay(chain, pars = component) + 
		theme(legend.position = "none") +
		theme(axis.title.x = element_blank())
	p <- plot_grid(p_trace, p_pdens, 
								 ncol = 2, nrow = 1) +
		theme(axis.text=element_text(size=12),
					axis.title=element_text(size=12))
	return(p)
}

rlist <- names(posterior.mode(m09a[[1]]$VCV))
tmp <- mcmc_dens_overlay(chains.m.vcv, pars = "id")
legend <- cowplot::get_legend(tmp)

## function to aggregate the plots into one graph
genPltVCV <- function(m){  #m=model_chain
	p <- plot_grid(
		genPlt1(m,"id"), genPlt1(m,"id:Year"), 
		genPlt1(m,"animal"), genPlt1(m,"Mother"),
		genPlt1(m,"Mother:Group:Year"), genPlt1(m,"Group"), 
		genPlt1(m,"Group:Year"), genPlt1(m,"Year"), 
		genPlt1(m,"Month:Year"), genPlt1(m,"units"), 
		ncol = 2, nrow = 5) +
		theme(axis.text=element_text(size=12),
					axis.title=element_text(size=12)
		)
	return(p)
}
## function to aggregate the plots into one graph
genPltSol <- function(m){  #m=model_chain
	p <- plot_grid(
		genPlt1(m,"(Intercept)"), 
		genPlt1(m,"SexMale"), 
		genPlt1(m,"poly(scale(age), degree = 3, raw = TRUE)1", "scale(age)ˆ1"), 
		genPlt1(m,"poly(scale(age), degree = 3, raw = TRUE)2", "scale(age)ˆ2"), 
		genPlt1(m,"poly(scale(age), degree = 3, raw = TRUE)3", "scale(age)ˆ3"), 
		genPlt1(m,"poly(scale(age), degree = 3, raw = TRUE)1:SexMale", "scale(age)ˆ1:Male"),
		genPlt1(m,"poly(scale(age), degree = 3, raw = TRUE)2:SexMale", "scale(age)ˆ2:Male"), 
		genPlt1(m,"poly(scale(age), degree = 3, raw = TRUE)3:SexMale", "scale(age)ˆ3:Male"), 
		genPlt1(m,"scale(grp_size)", "scale(group_size)"),
		genPlt1(m,"sin(rad_12)", "sin(month')"), 
		genPlt1(m,"cos(rad_12)", "cos(month')"), 
		genPlt1(m,"ENSO_Cool/La_Niña"), 
		genPlt1(m,"ENSO_Warm/El_Niño"), 
		legend, 
		ncol = 2, nrow = 7) +
		theme(axis.text=element_text(size=12),
					axis.title=element_text(size=12)
		)
	return(p)
}

library(Cairo)  #so that special characters print to PDF (i.e., ˆ)
library(gridGraphics)

grDevices::cairo_pdf("plots/supp/SI_File2a_Convergence_m09a.pdf", 
										 width = 9, height = 12,
										 onefile = T)
genPltVCV(chains.m.vcv)
genPltSol(chains.m.sol)
dev.off()

grDevices::cairo_pdf("plots/supp/SI_File2b_Convergence_p09a.pdf", 
										 width = 9, height = 12,
										 onefile = T)
genPltVCV(chains.p.vcv)
genPltSol(chains.p.sol)
dev.off()

#####

chains.m.vcv <- bindChains(m09b, "VCV")
chains.m.sol <- bindChains(m09b, "Sol")
chains.p.vcv <- bindChains(p09b, "VCV")
chains.p.sol <- bindChains(p09b, "Sol")

grDevices::cairo_pdf("plots/supp/SI_File2a_Convergence_m09b.pdf", 
										 width = 9, height = 12,
										 onefile = T)
genPltVCV(chains.m.vcv)
genPltSol(chains.m.sol)
dev.off()

grDevices::cairo_pdf("plots/supp/SI_File2b_Convergence_p09b.pdf", 
										 width = 9, height = 12,
										 onefile = T)
genPltVCV(chains.p.vcv)
genPltSol(chains.p.sol)
dev.off()

#####

chains.m.vcv <- bindChains(m09c, "VCV")
chains.m.sol <- bindChains(m09c, "Sol")
chains.p.vcv <- bindChains(p09c, "VCV")
chains.p.sol <- bindChains(p09c, "Sol")

grDevices::cairo_pdf("plots/supp/SI_File2a_Convergence_m09c.pdf", 
										 width = 9, height = 12,
										 onefile = T)
genPltVCV(chains.m.vcv)
genPltSol(chains.m.sol)
dev.off()

grDevices::cairo_pdf("plots/supp/SI_File2b_Convergence_p09c.pdf", 
										 width = 9, height = 12,
										 onefile = T)
genPltVCV(chains.p.vcv)
genPltSol(chains.p.sol)
dev.off()

###############################################################################
###############################################################################
###############################################################################

### model diagnostics ###

## more formal inspection ##
mods <- list(m09a, p09a)

### Gelman and Rubin's convergence diagnostic
for (i in 1:length(mods)){
	print("===== ===== ===== ===== =====")
	print(deparse(substitute(mods[[i]])))
	m <- mods[[i]]
	mclist.vcv <- mcmc.list(m[[1]]$VCV, m[[2]]$VCV, m[[3]]$VCV, m[[4]]$VCV)
	print("VCV")
	print(gelman.diag(mclist.vcv))
	mclist.sol <- mcmc.list(m[[1]]$Sol, m[[2]]$Sol, m[[3]]$Sol, m[[4]]$Sol)
	print("----- ----- ----- ----- ----- ----- ----- ----- -----")
	print("Sol")
	print(gelman.diag(mclist.sol))
}

### Heidelberger and Welch's convergence diagnostic
for (i in 1:length(mods)){
	print("===== ===== ===== ===== =====")
	print(deparse(substitute(mods[[i]])))
	m <- mods[[i]]
	for (chain in 1:4){
		print("----- ----- ----- ----- -----")
		print(paste0("Chain: ", chain))
		print(heidel.diag(m[[chain]]$VCV))
		print(heidel.diag(m[[chain]]$Sol))
	}
}

###############################################################################
###############################################################################
###############################################################################

### model diagnostics ###

## more formal inspection ##

mods <- list(m09b, p09b)

### Gelman and Rubin's convergence diagnostic
for (i in 1:length(mods)){
	print("===== ===== ===== ===== =====")
	print(deparse(substitute(mods[[i]])))
	m <- mods[[i]]
	mclist.vcv <- mcmc.list(m[[1]]$VCV, m[[2]]$VCV, m[[3]]$VCV, m[[4]]$VCV)
	print("VCV")
	print(gelman.diag(mclist.vcv))
	mclist.sol <- mcmc.list(m[[1]]$Sol, m[[2]]$Sol, m[[3]]$Sol, m[[4]]$Sol)
	print("----- ----- ----- ----- ----- ----- ----- ----- -----")
	print("Sol")
	print(gelman.diag(mclist.sol))
}

### Heidelberger and Welch's convergence diagnostic
for (i in 1:length(mods)){
	print("===== ===== ===== ===== =====")
	print(deparse(substitute(mods[[i]])))
	m <- mods[[i]]
	for (chain in 1:4){
		print("----- ----- ----- ----- -----")
		print(paste0("Chain: ", chain))
		print(heidel.diag(m[[chain]]$VCV))
		print(heidel.diag(m[[chain]]$Sol))
	}
}

###############################################################################
###############################################################################
################## SI Fig 3: Prior comparison ############################
###############################################################################
###############################################################################

library(ggplot2)
library(cowplot)
library(patchwork)

### prior comparison

### Define function for extracting MCMCglmm model information
### and placing into a data frame
genTabRE1 <- function(model){ #m=model
	mname <- deparse(substitute(model)) #store name of model
	m <- model[[1]]  #take first chain
	rlist <- names(posterior.mode(m$VCV)) #store list of random effects
	df_tmp <- data.frame( #create dataframe
		"Model" = mname, 
		"DIC" = m$DIC, 
		"Component" = rlist, #list
		"eff.samp" = effectiveSize(m$VCV), #list
		"posterior.mode" = posterior.mode(m$VCV), #list
		"HPDI95" = HPDinterval(m$VCV), #of posterior.mode
		"sum.VP.modes" = sum(posterior.mode(m$VCV)),
		"VP.mode" = posterior.mode(rowSums(m$VCV)), #posterior mode of sum of VCV components
		"latent.mode" = posterior.mode(m$VCV/rowSums(m$VCV)), #list, prop var
		"latent.HPDI95.lower" = HPDinterval(m$VCV/rowSums(m$VCV))[,1], #of prop var
		"latent.HPDI95.upper" = HPDinterval(m$VCV/rowSums(m$VCV))[,2], #of prop var
		row.names = 1:length(rlist)
	)
	df_tmp
}

genTabRE2 <- function(dre){
	dre$Model <- substr(dre$Model,1,4)
	dre$rank <- substr(dre$Model,4,4)
	#relabel
	dre$model <- "blank"
	dre$model[dre$rank=="a"] <- "parameter-expanded"
	dre$model[dre$rank=="b"] <- "inverse Wishart (V=1, nu=0.002)"
	dre$model[dre$rank=="c"] <- "inverse Wishart (V=1, nu=0.02)"
	#relabel
	dre$Component <- as.character(dre$Component)
	dre$Component[dre$Component=='animal'] <- 'Additive genetic'
	dre$Component[dre$Component=='Group'] <- 'GroupAlpha'
	dre$Component[dre$Component=='Group:Year'] <- 'GroupAlpha:Year'
	dre$Component[dre$Component=='Mother:Group:Year'] <- 'Mother:GroupAlpha:Year'
	dre$Component[dre$Component=='units'] <- 'residual'
	dre$Component[dre$Component=='id'] <- 'ID'
	dre$Component[dre$Component=='id:Year'] <- 'ID:Year'
	# reorder
	dre$sort[dre$Component=='ID'] <- 10
	dre$sort[dre$Component=='Mother'] <- 9
	dre$sort[dre$Component=='Additive genetic'] <- 8
	dre$sort[dre$Component=='ID:Year'] <- 7
	dre$sort[dre$Component=='Mother:GroupAlpha:Year'] <- 6
	dre$sort[dre$Component=='GroupAlpha'] <- 5
	dre$sort[dre$Component=='GroupAlpha:Year'] <- 4
	dre$sort[dre$Component=='Month:Year'] <- 3
	dre$sort[dre$Component=='Year'] <- 2
	dre$sort[dre$Component=='residual'] <- 1
	dre$Component <- factor(
		dre$Component, 
		levels = unique(dre$Component[order(dre$sort,dre$rank)]))
	dre$model <- factor(
		dre$model, 
		levels = unique(dre$model[order(dre$rank)]))
	dre <- dre[order(dre$Model, dre$sort),]
}

genPlotRE <- function(dre){
	plt <- ggplot(dre, 
								aes(y = Component, 
										x = latent.mode,  
										colour = model),
								poition = position_dodge(padding=1.0)) + 
		coord_cartesian(xlim = c(0, 1.0)) + 
		scale_color_manual(values = clst) + 
		geom_point(position = position_dodge(width = 1.0)) + 
		geom_tile(color = "grey") + 
		xlab("proportion variance explained (latent scale)") + 
		ylab("") + 
		labs(color = "prior for random effects") + 
		geom_linerange(aes(xmin = latent.HPDI95.lower, 
											 xmax = latent.HPDI95.upper, 
											 colour = model),
									 position = position_dodge(width = 1.0),
									 size = 1) + 
		theme_classic() + 
		theme(legend.position = c(0.7, 0.75)) + 
		theme(text = element_text(size=14)) + 
		guides(color = guide_legend(override.aes = list(size=3),
																reverse=TRUE)) + 
		theme(
			panel.grid.major.y = element_blank(),
			panel.grid.major.x = element_blank(),
			panel.grid.minor.x = element_blank())
}

clst <- c(
	"parameter-expanded" = "grey30",
	"inverse Wishart (V=1, nu=0.002)" = "grey50", 
	"inverse Wishart (V=1, nu=0.02)" = "grey80")

######

dreM <- rbind(genTabRE1(m09a),
							genTabRE1(m09b),
							genTabRE1(m09c))

table_RE_M <- genTabRE2(dreM)
pltM <- genPlotRE(table_RE_M) + 
	ggtitle("a) social vs alone")
(pltM)

######

dreP <- rbind(genTabRE1(p09a),
							genTabRE1(p09b),
							genTabRE1(p09c))

table_RE_P <- genTabRE2(dreP)
pltP <- genPlotRE(table_RE_P) + 
	ggtitle("b) number of partners")
(pltP)

plt <- pltM + pltP
plt[[1]] = plt[[1]] + theme(legend.position = "none")
plt[[2]] = plt[[2]] + theme(axis.text.y = element_blank(),
														axis.title.y = element_blank())
(plt)

######

pdf("plots/supp/SI_Fig3_PriorComparison.pdf", height = 7, width = 12)
plt
dev.off()

### save all full models in .RData file
save.image('Models1_fullModels.RData')