### load packages
library(MCMCglmm)
library(viridis)
library(bayesplot)
library(coda)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(patchwork)
if (!require("postMCMCglmm")) {
  devtools::install_github("JWiley/postMCMCglmm")
  library("postMCMCglmm", force=TRUE) # for extracting fixed effects
}

###############################################################################
###############################################################################
####### Figure 1: Distributions (raw data) #####################################
###############################################################################
###############################################################################

### load in original data
load("cebus.RData")

df <- cebus
df$p_soc <- df$social / df$n
df$p_par <- df$partners / df$n

###############################################################################

## plot distribution of age, group size, and sociality in population
theme_set(theme_classic(base_size = 7, base_family = 'Helvetica'))

### b) sampling ###
###################

ind1 <- as.data.frame(table(df$id, df$Year)) #sampling across diff years
ind1 <- ind1 %>% filter(Freq > 0)
ind1 <- as.data.frame(table(ind1$Var1))      #calendar years per id

ind1$Years <- "1-4"
ind1$Years[ind1$Freq>4] <- "5-9"
ind1$Years[ind1$Freq>9] <- "10+"
ind1$lev <- 1
ind1$lev[ind1$Freq>4] <- 2
ind1$lev[ind1$Freq>9] <- 3
ind1$Years <- factor(ind1$Years, levels = unique(ind1$Years[order(ind1$lev)]))
table(ind1$Freq)
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
# 1 32 52 38 18 33 36 19 22 36 20  9 10 17 10  2  6  2  9  4 

# plot distribution of calendar years per subjects
p3_yrs <- ggplot(data=ind1, aes(Freq, fill=Years)) + 
	geom_bar(col="black", alpha=.5, size=0.2) +
	scale_x_continuous(breaks=seq(0,20,5)) + 
	labs(title="b) Sampling", 
			 x="years per subject") + 
	scale_fill_manual(values = c("grey90","grey60","grey10")) + 
	theme(axis.text=element_text(size=5), 
				axis.title.y=element_text(size=5), 
				legend.key.size = unit(0.2, 'cm'), 
				legend.position = c(0.9, 0.9))
(p3_yrs)

ind2 <- as.data.frame(table(df$id)) #number of months per individual
mean(ind2$Freq)  #53.51862
range(ind2$Freq) #6 185

### c) group size ###
#####################

p5_grp <- ggplot(data=df, 
							aes(x = grp_size)) + 
	geom_bar(col="black", 
					 fill="#999999", 
					 size=0.2, 
					 alpha=.5) + 
	scale_x_continuous(breaks=seq(0,50,5)) + 
	labs(title="c) Group size", 
			 x="group size") +
	theme(axis.text=element_text(size=5), 
				axis.title.y=element_text(size=5))
(p5_grp)

# d) tenure length #
####################

alpha <- as.data.frame(table(df$Group,df$Year))
alpha <- alpha %>% filter(Freq > 0)
alpha <- as.data.frame(table(alpha$Var1)) #calendar years per alpha
mean(alpha$Freq)  #3.571429
range(alpha$Freq)  #1 14

p6_alp <- ggplot(
	data=alpha, 
	aes(x = Freq)) + 
	geom_bar(col="black", 
					 fill="#999999", 
					 size=0.2, 
					 alpha=.5) + 
	scale_x_continuous(breaks=seq(1,15,2)) + 
	scale_y_continuous(breaks=seq(0,15,3)) + 
  labs(title="d) Alpha tenure length", 
  		 x="alpha tenure length (years)") +
	theme(axis.text=element_text(size=5), 
				axis.title.y=element_text(size=5))
(p6_alp)

##### combine 3 figures together #####
#####################################

plt_3 <- plot_grid(p3_yrs, p5_grp, p6_alp, 
               ncol = 1, nrow = 3)
(plt_3)

##### ##### ##### ##### ##### #####

# generate capuchin image
library(sketcher)
img <- sketcher::im_load('MEE.2022.02.19.NSC.IMG_9201.JPG')
MEE <- sketcher::sketch(img, 
												lineweight = 4, 
												smooth = 0.05)
plot(MEE)

# plot distributions of sociality scores
p <- ggplot(df, aes(x=p_soc, y=p_par, size=n)) + 
	geom_point(alpha=0.1) + scale_size_area(max_size=1) + 
	theme_bw(base_size=7, base_family='Helvetica') + 
	xlab("proportion of time social") + 
	ylab("average number of partners") + 
	ggtitle("a) Sociality") +
	theme(axis.text=element_text(size=5),
				legend.position="none")
p1 <- ggExtra::ggMarginal(p, type = "histogram")

plt_soc <- cowplot::ggdraw() +
	draw_plot(p1) +
	draw_image(plot(MEE), 
						 x = 0.131, y = 0.53, 
						 width = 0.4/1.25, height = 0.3/1.25)

Fig1 <- plot_grid(plt_soc, plt_3, 
									ncol = 2, nrow = 1)
pdf("plots/Fig1_distributionsNew.pdf", width = 7.08661, height = 3.543305)
Fig1
dev.off()

###############################################################################
###############################################################################
###### Figure 2: Age by Sex interaction plot ######
###############################################################################
###############################################################################

### load in data
load("cebus.RData")               #original data
load("Models1_fullModels.RData")  #full models

logit2prob <- function(logit){
	odds <- exp(logit)
	prob <- odds / (1 + odds)
	return(prob)
}

xvals <- seq(min(scale(cebus$age)), 
						 max(scale(cebus$age)), 
						 length.out=length(df$age))

df <- cebus
df$p_soc <- df$social / df$n
df$p_par <- df$partners / df$n

########

m <- m09a[[1]]

pop.int <- posterior.mode(m$Sol[, "(Intercept)"])
pop.slope <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)1"])
pop.quad <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)2"])
pop.cub <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)3"])
pop.mal <- posterior.mode(m$Sol[, "SexMale"])
pop.slope.mal <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)1:SexMale"])
pop.quad.mal <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)2:SexMale"])
pop.cub.mal <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)3:SexMale"])

yvals_fem <- pop.int + xvals*pop.slope + I(xvals^2)*pop.quad + I(xvals^3)*pop.cub
yvals_mal <- pop.int + pop.mal + xvals*pop.slope + I(xvals^2)*pop.quad + I(xvals^3)*pop.cub + 
	xvals*pop.slope.mal + I(xvals^2)*pop.quad.mal + I(xvals^3)*pop.cub.mal

m_yvals_fem <- logit2prob(yvals_fem)
m_yvals_mal <- logit2prob(yvals_mal)

########

m <- p09a[[1]]

pop.int <- posterior.mode(m$Sol[, "(Intercept)"])
pop.slope <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)1"])
pop.quad <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)2"])
pop.cub <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)3"])
pop.mal <- posterior.mode(m$Sol[, "SexMale"])
pop.slope.mal <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)1:SexMale"])
pop.quad.mal <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)2:SexMale"])
pop.cub.mal <- posterior.mode(m$Sol[, "poly(scale(age), degree = 3, raw = TRUE)3:SexMale"])

yvals_fem <- pop.int + xvals*pop.slope + I(xvals^2)*pop.quad + I(xvals^3)*pop.cub
yvals_mal <- pop.int + pop.mal + xvals*pop.slope + I(xvals^2)*pop.quad + I(xvals^3)*pop.cub + 
	xvals*pop.slope.mal + I(xvals^2)*pop.quad.mal + I(xvals^3)*pop.cub.mal

p_yvals_fem <- exp(yvals_fem)
p_yvals_mal <- exp(yvals_mal)

########

library(ggplot2)
library(scales)  #for hue_pal()

p <- ggplot(df, aes(x=scale(age), y=p_soc, color=Sex, size=n)) + 
	geom_point(alpha=0.2) + scale_size_area(max_size=1) + 
	geom_line(aes(x=xvals, y=m_yvals_fem), 
						color = hue_pal()(2)[1], 
						size = 0.5) + 
	geom_line(aes(x=xvals, y=m_yvals_mal), 
						color = hue_pal()(2)[2], 
						size = 0.5) + 
	scale_x_continuous(
		sec.axis = dup_axis(breaks = c(-1,0,1,2,3), 
												labels = c(2.3, 9.3, 16.2, 23.1, 30.1),
												name = "age (years)")) + 
	theme_bw(base_size=7, base_family="Helvetica") + 
	xlab("age (scaled)") + 
	ylab("proportion of time social") + 
	ggtitle("a) social versus alone") + 
	theme(legend.position = 'none') + 
	guides(color=guide_legend(override.aes=list(fill=NA)))
p1 <- ggExtra::ggMarginal(p, 
													type = "histogram", 
													margins = 'x', 
													binwidth = 0.05, 
													size = 4,  #main plot 4x larger than hist
													groupColour = T, 
													groupFill = T, 
													alpha = 0.5)
(p1)

p <- ggplot(df, aes(x=scale(age), y=p_par, color=Sex, size=n)) + 
	geom_point(alpha=0.2) + scale_size_area(max_size=1) + 
	geom_line(aes(x=xvals, y=p_yvals_fem), 
						color = hue_pal()(2)[1], 
						size = 0.5) + 
	geom_line(aes(x=xvals, y=p_yvals_mal), 
						color = hue_pal()(2)[2], 
						size = 0.5) + 
	scale_x_continuous(
		sec.axis = dup_axis(breaks = c(-1,0,1,2,3), 
												labels = c(2.3, 9.3, 16.2, 23.1, 30.1),
												name = "age (years)")) + 
	theme_bw(base_size=7, base_family="Helvetica") + 
	xlab("age (scaled)") + 
	ylab("average number of partners") + 
	ggtitle("b) number of partners") + 
	theme(legend.title = element_blank(), 
				legend.justification = c(1, 1), 
				legend.position = c(0.95, 0.95), 
				legend.key = element_rect(colour = "transparent", fill = "white")) + 
	guides(size=FALSE)
p2 <- ggExtra::ggMarginal(p, 
													type = "histogram", 
													margins = 'x', 
													binwidth = 0.05, 
													size = 4,  #main plot 4x larger than hist
													groupColour = T, 
													groupFill = T, 
													alpha = 0.5)
(p2)

plt_age <- cowplot::plot_grid(p1, p2, 
															ncol = 2, nrow = 1)

width = 7.08661
pdf("plots/Fig2_Soc_by_AgeSex.pdf", width = width, height = width/2)
plt_age
dev.off()

###############################################################################
###############################################################################
####### SI Table 1: Random Effects table  #####################################
###############################################################################
###############################################################################

# load previously run models
load("Models1_fullModels.RData")     #full models
load("Models2_reducedModels.RData")  #reduced models

# functions taken from de Villemeruil 2021 tutorial
## for adding fixed effects variance to denominator
compute_varpred <- function(beta, design_matrix) {
	var(as.vector(design_matrix %*% beta))
}

### Define function for extracting MCMCglmm model information
### and placing into a data frame
genTabRE <- function(m){ #m=model
	mname <- deparse(substitute(m)) #store name of model
	rlist <- names(posterior.mode(m$VCV)) #store list of random effects
	vf <- apply(m[["Sol"]], 1, compute_varpred, design_matrix = m[["X"]])
	df_tmp <- data.frame( #create dataframe
		"Model" = mname, 
		"Component" = rlist, #list
		"eff.samp" = effectiveSize(m$VCV), #list
		"posterior.mode" = posterior.mode(m$VCV), #list
		"HPDI95.lower" = HPDinterval(m$VCV)[,1], #of posterior.mode
		"HPDI95.upper" = HPDinterval(m$VCV)[,2], #of posterior.mode
		"latent.mode" = posterior.mode(m$VCV/rowSums(m$VCV)), #list, prop var
		"latent.HPDI95.lower" = HPDinterval(m$VCV/rowSums(m$VCV))[,1], #of prop var
		"latent.HPDI95.upper" = HPDinterval(m$VCV/rowSums(m$VCV))[,2], #of prop var
		"FE.latent.mode" = posterior.mode(m$VCV/(rowSums(m$VCV) + vf)), #list, prop var
		"FE.latent.HPDI95.lower" = HPDinterval(m$VCV/(rowSums(m$VCV) + vf))[,1], #of prop var
		"FE.latent.HPDI95.upper" = HPDinterval(m$VCV/(rowSums(m$VCV) + vf))[,2], #of prop var
		row.names = 1:length(rlist)
	)
	df_tmp
}

## generate table with results of variance components from models ##
## use first chain only
latent <- rbind(genTabRE(m09a[[1]]), genTabRE(m09a.0[[1]]), 
						 genTabRE(m09a.1[[1]]), genTabRE(m09a.2[[1]]), 
						 genTabRE(m09a.3[[1]]), genTabRE(m09a.4[[1]]), 
						 genTabRE(p09a[[1]]), genTabRE(p09a.0[[1]]), 
						 genTabRE(p09a.1[[1]]), genTabRE(p09a.2[[1]]), 
						 genTabRE(p09a.3[[1]]), genTabRE(p09a.4[[1]]))

latent$Model <- substr(latent$Model, 1, nchar(latent$Model)-5)
latent$eff.samp <- round(latent$eff.samp, 0)

### 

# load qgglmm objects
load("QGglmmTransform.RData")

# generate table with posterior mode and 95HPDI
genTabICC <- function(df){
	df_tmp <- data.frame()
	mnames <- unique(df[["model"]])
	for (m in mnames){ #subset by model
		tmp1 <- subset(df,model==m)
		tmp1 <- droplevels(tmp1)
		rlst <- unique(tmp1$component)
		for (r in rlst){ #subset by variance component
			tmp2 <- subset(tmp1,component==r)
			tmp2 <- droplevels(tmp2)
			tmp3 <- data.frame( #create dataframe
				"Model" = m, 
				"Component" = r, 
				"data.posterior.mode" = posterior.mode(as.mcmc(tmp2[["icc.obs"]])), 
				"data.HPDI95.lower" = HPDinterval(as.mcmc(tmp2[["icc.obs"]]))[,1], 
				"data.HPDI95.upper" = HPDinterval(as.mcmc(tmp2[["icc.obs"]]))[,2], 
				row.names = 1:length(r))
			df_tmp <- rbind(df_tmp,tmp3)
		}
	}
	df_tmp
}

FE_no <- rbind(
	genTabICC(M_ICC1), genTabICC(M_ICC1.0), genTabICC(M_ICC1.1), 
	genTabICC(M_ICC1.2), genTabICC(M_ICC1.3), genTabICC(M_ICC1.4),
	genTabICC(P_ICC1), genTabICC(P_ICC1.0), genTabICC(P_ICC1.1), 
	genTabICC(P_ICC1.2), genTabICC(P_ICC1.3), genTabICC(P_ICC1.4)
)
FE_no$Model <- substr(FE_no$Model, 1, nchar(FE_no$Model)-5)
head(FE_no)

dre <- merge(latent, FE_no, by=c("Model","Component"), all.x = T)

###

M_ICC2_sub   <- rbind(M_ICC2[[1]],M_ICC2[[2]],M_ICC2[[3]],M_ICC2[[4]],M_ICC2[[5]],
						          M_ICC2[[6]],M_ICC2[[7]],M_ICC2[[8]],M_ICC2[[9]],M_ICC2[[10]])
M_ICC2.0_sub <- rbind(M_ICC2.0[[1]],M_ICC2.0[[2]],M_ICC2.0[[3]],M_ICC2.0[[4]],M_ICC2.0[[5]],
										  M_ICC2.0[[6]],M_ICC2.0[[7]],M_ICC2.0[[8]],M_ICC2.0[[9]],M_ICC2.0[[10]])
M_ICC2.1_sub <- rbind(M_ICC2.1[[1]],M_ICC2.1[[2]],M_ICC2.1[[3]],M_ICC2.1[[4]],M_ICC2.1[[5]],
										  M_ICC2.1[[6]],M_ICC2.1[[7]],M_ICC2.1[[8]],M_ICC2.1[[9]],M_ICC2.1[[10]])
M_ICC2.2_sub <- rbind(M_ICC2.2[[1]],M_ICC2.2[[2]],M_ICC2.2[[3]],M_ICC2.2[[4]],M_ICC2.2[[5]],
										  M_ICC2.2[[6]],M_ICC2.2[[7]],M_ICC2.2[[8]],M_ICC2.2[[9]],M_ICC2.2[[10]])
M_ICC2.3_sub <- rbind(M_ICC2.3[[1]],M_ICC2.3[[2]],M_ICC2.3[[3]],M_ICC2.3[[4]],M_ICC2.3[[5]],
										  M_ICC2.3[[6]],M_ICC2.3[[7]],M_ICC2.3[[8]],M_ICC2.3[[9]],M_ICC2.3[[10]])
M_ICC2.4_sub <- rbind(M_ICC2.4[[1]],M_ICC2.4[[2]],M_ICC2.4[[3]],M_ICC2.4[[4]],M_ICC2.4[[5]],
										  M_ICC2.4[[6]],M_ICC2.4[[7]],M_ICC2.4[[8]],M_ICC2.4[[9]],M_ICC2.4[[10]])

FE_Yes <- rbind(genTabICC(M_ICC2_sub), genTabICC(M_ICC2.0_sub), 
								genTabICC(M_ICC2.1_sub), genTabICC(M_ICC2.2_sub), 
								genTabICC(M_ICC2.3_sub), genTabICC(M_ICC2.4_sub), 
								genTabICC(P_ICC2), 
								genTabICC(P_ICC2.1),  genTabICC(P_ICC2.2),
								genTabICC(P_ICC2.3),  genTabICC(P_ICC2.4))
FE_Yes$Model <- substr(FE_Yes$Model, 1, nchar(FE_Yes$Model)-5)
colnames(FE_Yes) <- c("Model", "Component", "FE.data.posterior.mode", 
											"FE.data.HPDI95.lower", "FE.data.HPDI95.upper")

dre <- merge(dre, FE_Yes, by=c("Model","Component"), all.x = T)

###

cols <- c("posterior.mode", "HPDI95.lower", "HPDI95.upper", 
					"latent.mode", "latent.HPDI95.lower", "latent.HPDI95.upper",
					"FE.latent.mode", "FE.latent.HPDI95.lower", "FE.latent.HPDI95.upper",
					"data.posterior.mode", "data.HPDI95.lower", "data.HPDI95.upper", 
					"FE.data.posterior.mode", "FE.data.HPDI95.lower", "FE.data.HPDI95.upper")
dre[cols] <- lapply(dre[cols], round, 5)


olst <- data.frame(
	Component = dimnames(m09a[[1]]$VCV)[[2]],
	sort = c(1:10)
)
dre <- merge(dre, olst)

#relabel
dre$Component <- as.character(dre$Component)
dre$Component[dre$Component=='animal'] <- 'Additive genetic'
dre$Component[dre$Component=='Group'] <- 'GroupAlpha'
dre$Component[dre$Component=='Group:Year'] <- 'GroupAlpha:Year'
dre$Component[dre$Component=='Mother:Group:Year'] <- 'Mother:GroupAlpha:Year'
dre$Component[dre$Component=='units'] <- 'residual'
dre$Component[dre$Component=='id'] <- 'ID'
dre$Component[dre$Component=='id:Year'] <- 'ID:Year'

dre$Component <- factor(dre$Component, 
												levels = unique(dre$Component[order(dre$sort,dre$Model)]))
dre <- dre[order(dre$Model, dre$sort),]
head(dre,20)

write.table(dre, "tables/SI_Table1_RE.txt", row.names = F)

###############################################################################
###############################################################################
##### Figure 3: Fixed effects estimates from models ####
###############################################################################
###############################################################################

flabels <- c("(Intercept)" = "(Intercept)", 
						 "poly(scale(age), degree = 3, raw = TRUE)1" = expression(age^1), 
						 "poly(scale(age), degree = 3, raw = TRUE)2" = expression(age^2), 
						 "poly(scale(age), degree = 3, raw = TRUE)3" = expression(age^3), 
						 "poly(scale(age), degree = 3, raw = TRUE)1:SexMale" = expression(age^1~"*"~male),
						 "poly(scale(age), degree = 3, raw = TRUE)2:SexMale" = expression(age^2~"*"~male),
						 "poly(scale(age), degree = 3, raw = TRUE)3:SexMale" = expression(age^3~"*"~male),
						 "SexMale" = "male", 
						 "scale(grp_size)" = "group size", 
						 "sin(rad_12)" = "sin(month')",   
						 "cos(rad_12)" = "cos(month')", 
						 "ENSO_Cool/La_Niña" = "ENSO: La Niña (cool)", 
						 "ENSO_Warm/El_Niño" = "ENSO: El Niño (warm)"
)
forders <- c("(Intercept)", 
						 "poly(scale(age), degree = 3, raw = TRUE)1", 
						 "poly(scale(age), degree = 3, raw = TRUE)2", 
						 "poly(scale(age), degree = 3, raw = TRUE)3", 
						 "poly(scale(age), degree = 3, raw = TRUE)1:SexMale",
						 "poly(scale(age), degree = 3, raw = TRUE)2:SexMale",
						 "poly(scale(age), degree = 3, raw = TRUE)3:SexMale",
						 "SexMale", 
						 "scale(grp_size)", 
						 "sin(rad_12)",   
						 "cos(rad_12)", 
						 "ENSO_Cool/La_Niña", 
						 "ENSO_Warm/El_Niño"
)

library(bayesplot)
library(Cairo)  #so that special characters print to PDF (i.e., ˆ)

pltFE <- function(m, plt_title){
	bayesplot::mcmc_intervals(
		m$Sol,
		prob = 0.95,
		prob_outer = 0.99,
		inner_size = 2,
		point_size = 1
	) + 
		ggplot2::theme_bw(base_size=7, base_family="Helvetica") + 
		ggplot2::scale_y_discrete(labels = flabels, limits = rev(forders)) + 
		ggplot2::coord_cartesian(xlim=c(-0.6,0.6)) + 
		ggplot2::scale_x_continuous(breaks=seq(-0.6, 0.6, 0.2), 
																labels=c(-.6,-.4,-.2,0,.2,.4,.6)) + 
		ggplot2::ggtitle(plt_title) + 
		ggplot2::theme(axis.text=element_text(size=5),
									 panel.grid.major.x = element_line(colour = "gray80"), 
									 panel.grid.major.y = element_blank(),
									 panel.grid.minor.y = element_blank())
}
fixed.m09a <- pltFE(m09a[[1]], "a) social versus alone")
fixed.p09a <- pltFE(p09a[[1]], "b) number of partners")

p <- fixed.m09a + fixed.p09a
p[[2]] <- p[[2]] + 
	theme(axis.text.y = element_blank(),
				axis.title.y = element_blank())

(p)

width = 7.08661
grDevices::cairo_pdf("plots/Fig3_FE_mcmcIntervals.pdf", 
										 width = width, height = width/3)
p
dev.off()

###############################################################################
###############################################################################
####### SI Table 2: Fixed Effects table  #####################################
###############################################################################
###############################################################################

# load models
load("Models1_fullModels.RData")
load("Models2_reducedModels.RData")

### Define function for extracting MCMCglmm model information
### and placing into a data frame
genTabFE <- function(m){ #m=model
	mname <- deparse(substitute(m)) #store name of model
	flist <- names(posterior.mode(m$Sol)) #store list of fixed effects
	df_tmp <- data.frame( #create dataframe
		"Model" = mname, 
		"Component" = flist,  #list
		"eff.samp" = effectiveSize(m$Sol),  #list
		"posterior.mode" = posterior.mode(m$Sol),  #list
		"HPDI95.lower" = HPDinterval(m$Sol)[,1],  #of posterior.mode
		"HPDI95.upper" = HPDinterval(m$Sol)[,2],  #of posterior.mode
		row.names = 1:length(flist)
	)
	df_tmp
}

dfe <- rbind(genTabFE(m09a[[1]]), genTabFE(m09a.0[[1]]), genTabFE(m09a.1[[1]]), 
						 genTabFE(m09a.2[[1]]), genTabFE(m09a.3[[1]]), genTabFE(m09a.4[[1]]), 
						 genTabFE(p09a[[1]]), genTabFE(p09a.0[[1]]), genTabFE(p09a.1[[1]]), 
						 genTabFE(p09a.2[[1]]), genTabFE(p09a.3[[1]]), genTabFE(p09a.4[[1]]))

cols <- c("posterior.mode", "HPDI95.lower", "HPDI95.upper")
dfe[cols] <- lapply(dfe[cols], round, 4)
dfe$eff.samp <- round(dfe$eff.samp, 0)
dfe$Model <- substr(dfe$Model, 1, nchar(dfe$Model)-5)

#relabel
dfe$Component <- as.character(dfe$Component)
dfe$Component[dfe$Component=='scale(grp_size)'] <- 'group size'
dfe$Component[dfe$Component=='poly(scale(age), degree = 3, raw = TRUE)1'] <- "age^1"
dfe$Component[dfe$Component=='poly(scale(age), degree = 3, raw = TRUE)2'] <- "age^2"
dfe$Component[dfe$Component=='poly(scale(age), degree = 3, raw = TRUE)3'] <- "age^3"
dfe$Component[dfe$Component=='SexMale'] <- "male"
dfe$Component[dfe$Component=='poly(scale(age), degree = 3, raw = TRUE)1:SexMale'] <- "age^1:male"
dfe$Component[dfe$Component=='poly(scale(age), degree = 3, raw = TRUE)2:SexMale'] <- "age^2:male"
dfe$Component[dfe$Component=='poly(scale(age), degree = 3, raw = TRUE)3:SexMale'] <- "age^3:male"

dfe$orderR[dfe$Component=="(Intercept)"] <- 1
dfe$orderR[dfe$Component=='age^1'] <- 2
dfe$orderR[dfe$Component=='age^2'] <- 3
dfe$orderR[dfe$Component=='age^3'] <- 4
dfe$orderR[dfe$Component=='age^1:male'] <- 5
dfe$orderR[dfe$Component=='age^2:male'] <- 6
dfe$orderR[dfe$Component=='age^3:male'] <- 7 
dfe$orderR[dfe$Component=='male'] <- 8
dfe$orderR[dfe$Component=='group size'] <- 9
dfe$orderR[dfe$Component=="sin(rad_12)"] <- 10
dfe$orderR[dfe$Component=="cos(rad_12)"] <- 11
dfe$orderR[dfe$Component=="ENSO_Cool/La_Niña"] <- 12
dfe$orderR[dfe$Component=="ENSO_Warm/El_Niño"] <- 13
dfe$orderR[dfe$Component=="log(n)"] <- 14

dfe$Component <- factor(dfe$Component, levels = unique(dfe$Component[order(dfe$orderR,dfe$Model)]))
dfe <- dfe[order(dfe$Model, dfe$orderR),]
head(dfe,20)

plotFE <- subset(dfe, Component != "log(n)")
plotFE <- droplevels(plotFE)

dfe_txt <- dfe[order(dfe$Model, dfe$orderR),]
cols <- c("Model", "Component", "eff.samp", "posterior.mode", "HPDI95.lower", "HPDI95.upper")
dfe_txt <- dfe_txt[,cols]
head(dfe_txt, 15)

write.table(dfe_txt, "tables/SI_Table1_FE.txt", row.names = F)

###############################################################################
###############################################################################
###############################################################################
##### Figure 4: RE proportion of variance explained (latent scale) ####
###############################################################################
###############################################################################
###############################################################################

library(ggplot2)
library(patchwork)
library(bayesplot)

load("Models1_fullModels.RData")

rlst <- c("id",
					"Mother", 
					"animal",
					"id:Year", 
					"Mother:Group:Year",
					"Group",
					"Group:Year",
					"Year",
					"Month:Year",
					"units")
rlab <- c("ID",
					"Mother",
					"additive genetic",
					"ID:Year", 
					"Mother:GroupAlpha:Year",
					"GroupAlpha",
					"GroupAlpha:Year",
					"Year",
					"Month:Year",
					"residual")

library(bayesplot)
library(Cairo)  #so that special characters print to PDF (i.e., ˆ)

pltRE <- function(m, plt_title){
	bayesplot::mcmc_intervals(
		m$VCV/rowSums(m$VCV),
		prob = 0.95,
		prob_outer = 0.99,
		inner_size = 2,
		point_size = 1
	) + 
		ggplot2::theme_bw(base_size=7, base_family="Helvetica") + 
		ggplot2::scale_y_discrete(labels = rev(rlab), limits = rev(rlst)) + 
		ggplot2::coord_cartesian(xlim=c(0,1)) + 
		ggplot2::scale_x_continuous(breaks=seq(0, 1, 0.2)) + 
		ggplot2::ggtitle(plt_title) + 
		ggplot2::theme(
					axis.text = element_text(size = 5), 
					panel.grid.major.x = element_line(colour = "gray80"), 
					panel.grid.major.y = element_blank(),
					panel.grid.minor.y = element_blank()) + 
		ggplot2::xlab("proportion of variance explained (latent scale)")
}
random.m09a <- pltRE(m09a[[1]], "a) social versus alone")
random.p09a <- pltRE(p09a[[1]], "b) number of partners")

p <- random.m09a + random.p09a
p[[2]] <- p[[2]] + 
	theme(axis.text.y = element_blank(),
				axis.title.y = element_blank())
(p)

width = 7.08661
pdf("plots/Fig4_RE_mcmcIntervals.pdf", 
		width = width, height = width/3)
p
dev.off()

###############################################################################
###############################################################################
####################### Figure 5  ###########################################
###############################################################################
###############################################################################

## GENERATE REPEATABILITLY ESTIMATES ###

library(QGglmm)  #for estimates on observed data scale

load("cebus.RData")  #original data
load("Models1_fullModels.RData")

# sketch capuchin
library(sketcher)
img <- sketcher::im_load('UTT(left)_VYY(right).2022.03.14.NSC.IMG_2314.jpg')
UTT <- sketcher::sketch(img, 
												lineweight = 4, 
												smooth = 0.05)
plot(UTT)

###############################################################################

# generate posterior distributions for ICC estimates
## poisson model,
## not accounting for fixed effects variance in denominator
m <- p09a[[1]]  #use first chain
var.mu    <- m[["Sol"]][, "(Intercept)"]
var.tot   <- rowSums(m[["VCV"]])
## components
var.rep.l <- m[["VCV"]][,"animal"] + m[["VCV"]][,"Mother"] + m[["VCV"]][,"id"]
var.rep.s <- m[["VCV"]][,"animal"] + m[["VCV"]][,"Mother"] + m[["VCV"]][,"id"] + 
	m[["VCV"]][,"Mother:Group:Year"] + m[["VCV"]][,"id:Year"]
var.grp.l <- m[["VCV"]][,"Group" ]
var.grp.s <- m[["VCV"]][,"Group" ] + m[["VCV"]][, "Group:Year"]

# function for generating posterior distributions for ICC estimates
## poisson models,
## not accounting for fixed effects variance in denominator
library(purrr)  #for pmap_dfr function
genICC_P <- function(var_mu, var_comp, var_tot, lab1, lab2) {
	#m=model
	tmp <- pmap_dfr(
		list(
			mu       = var_mu,
			var.comp = var_comp,
			var.p    = var_tot
		),
		QGicc,
		model = "Poisson.log",
		verbose = FALSE
	)
	tmp$model <- "p09a[[1]]"
	tmp$component <- lab1
	tmp$behavior <- lab2
	tmp$FE.var <- "No"
	tmp$model.type <- "Poisson.log"
	tmp$scale <- "data scale"
	return(tmp)
}

repP_long_data <- genICC_P(var.mu, var.rep.l, var.tot,
													 "long-term", "Repeatability (number of partners)")
repP_shor_data <- genICC_P(var.mu, var.rep.s, var.tot,
													 "short-term", "Repeatability (number of partners)")
grpP_long_data <- genICC_P(var.mu, var.grp.l, var.tot,
													 "long-term", "Group of residence (number of partners)")
grpP_shor_data <- genICC_P(var.mu, var.grp.s, var.tot,
													 "short-term", "Group of residence (number of partners)")

dataP_scale_rep <- rbind(repP_long_data, repP_shor_data, 
												 grpP_long_data, grpP_shor_data)

# generate posterior distributions for ICC estimates
## poisson model,
## not accounting for fixed effects variance in denominator
m       <- m09a[[1]]  #use first chain
var.mu  <- m[["Sol"]][, "(Intercept)"]
var.tot <- rowSums(m[["VCV"]])
## components
var.rep.l <- m[["VCV"]][,"animal"] + m[["VCV"]][,"Mother"] + m[["VCV"]][,"id"]
var.rep.s <- m[["VCV"]][,"animal"] + m[["VCV"]][,"Mother"] + m[["VCV"]][,"id"] + 
	m[["VCV"]][,"Mother:Group:Year"] + m[["VCV"]][,"id:Year"]
var.grp.l <- m[["VCV"]][,"Group" ]
var.grp.s <- m[["VCV"]][,"Group" ] + m[["VCV"]][, "Group:Year"]

# function for generating posterior distributions for ICC estimates
## multinomial2 models,
## not accounting for fixed effects variance in denominator
obs <- round(mean(cebus$n)) ### mean #obs per subject per month (i.e., 32)
genICC_M <- function(var_mu, var_comp, var_tot, lab1, lab2) {
	#m=model, col=variance component
	tmp <- pmap_dfr(
		list(
			mu       = var_mu,
			var.comp = var_comp,
			var.p    = var_tot
		),
		QGicc,
		model = "binomN.logit",
		n.obs = obs,
		verbose = FALSE
	)
	tmp$model <- "m09a[[1]]"
	tmp$component <- lab1
	tmp$behavior <- lab2
	tmp$FE.var <- "No"
	tmp$model.type <- "binomN.logit"
	tmp$scale <- "data scale"
	return(tmp)
}

repM_long_data <- genICC_M(var.mu, var.rep.l, var.tot,
													 "long-term", "Repeatability (social vs alone)")
repM_shor_data <- genICC_M(var.mu, var.rep.s, var.tot,
													 "short-term", "Repeatability (social vs alone)")
grpM_long_data <- genICC_M(var.mu, var.grp.l, var.tot,
													 "long-term", "Group of residence (social vs alone)")
grpM_shor_data <- genICC_M(var.mu, var.grp.s, var.tot,
													 "short-term", "Group of residence (social vs alone)")

dataM_scale_rep <- rbind(repM_long_data, repM_shor_data, 
												 grpM_long_data, grpM_shor_data)

scale_data_rep <- rbind(dataM_scale_rep, dataP_scale_rep)
scale_data_rep <- scale_data_rep[,c('icc.obs','component','behavior','scale')]
names(scale_data_rep)[names(scale_data_rep) == "icc.obs"] <- "var1"

## latent scale

m <- m09a[[1]]
repM_long <- as.data.frame(
	(m[["VCV"]][,"id"]+m[["VCV"]][,"Mother"]+m[["VCV"]][,"animal"])/rowSums(m$VCV)
)
repM_long$component <- "long-term"
repM_short <- as.data.frame(
	(m[["VCV"]][,"id:Year"]+m[["VCV"]][,"Mother:Group:Year"]+m[["VCV"]][,"id"]+
	 	m[["VCV"]][,"Mother"]+m[["VCV"]][,"animal"])/rowSums(m$VCV)
)
repM_short$component <- "short-term"

repM <- rbind(repM_long, repM_short)
repM$behavior <- "Repeatability (social vs alone)"

m <- p09a[[1]]
repP_long <- as.data.frame(
	(m[["VCV"]][,"id"]+m[["VCV"]][,"Mother"]+m[["VCV"]][,"animal"])/rowSums(m$VCV)
)
repP_long$component <- "long-term"

repP_short <- as.data.frame(
	(m[["VCV"]][,"id:Year"]+m[["VCV"]][,"Mother:Group:Year"]+m[["VCV"]][,"id"]+
	 	m[["VCV"]][,"Mother"]+m[["VCV"]][,"animal"])/rowSums(m$VCV)
)
repP_short$component <- "short-term"

repP <- rbind(repP_long, repP_short)
repP$behavior <- "Repeatability (number of partners)"

scale_latent_rep <- rbind(repM, repP)
scale_latent_rep$scale <- "latent scale"

## 

fig <- rbind(scale_latent_rep, scale_data_rep)
fig <- fig[grep("Repeatability", fig$behavior), ]
fig <- droplevels(fig)

fig$order <- 1
fig$order[fig$behavior=="Repeatability (number of partners)"] <- 2
fig$behavior <- factor(fig$behavior, 
												levels = unique(fig$behavior[order(fig$order)]))

fig$order <- 1
fig$order[fig$scale=="data scale"] <- 2
fig$scale <- factor(fig$scale, 
											 levels = unique(fig$scale[order(fig$order)]))

# manually set colors
clst <- c("long-term" = "grey30", 
					"short-term" = "grey80")

# function for generating supplementary figure
pltRep <- function(df) {
	ggplot(
		df, 
		aes(
			x=var1, 
			fill=component)) + 
		theme_bw(base_size=7, base_family="Helvetica") + 
		geom_density(alpha=0.7)  + 
		facet_grid(scale~behavior) + 
		coord_cartesian(xlim = c(0, 0.5)) + 
		scale_fill_manual(values=clst) + 
		theme(
			panel.grid.major.x = element_line(
				colour='grey', size=0.25, linetype="solid"),
			panel.grid.major.y = element_blank(),
			panel.grid.minor.y = element_blank(),
			panel.background = element_blank(),
			panel.spacing.x=unit(1.0, "lines"),
			panel.spacing.y=unit(1.0, "lines")) + 
		theme(
			legend.title = element_blank(),
			legend.position = c(.99, .99),
			legend.justification = c("right", "top"),
			legend.box.just = "right",
			legend.key.height = unit(0.5, "cm"),
			legend.key.width = unit(0.5, "cm"),
			legend.text.align = 0, 
			legend.background = element_rect(fill=alpha('white', 0.4))) + 
		labs(
			title="",
			x    ="proportion variance explained")
}

# function for generating figure 5 (repeatability)
pltRep5 <- function(df) {
	ggplot(df, aes(x=var1, fill=component)) + 
		theme_bw(base_size=7, base_family="Helvetica") + 
		geom_density(alpha=0.7)  + 
		facet_wrap(vars(behavior), nrow=2, ncol=1) + 
		coord_cartesian(xlim = c(0, 0.5)) + 
		scale_fill_manual(values=clst) + 
		theme(
			panel.grid.major.x = element_line(
				colour='grey', size=0.25, linetype="solid"),
			panel.grid.major.y = element_blank(),
			panel.grid.minor.y = element_blank(),
			panel.background = element_blank(),
			panel.spacing.x=unit(1.0, "lines"),
			panel.spacing.y=unit(1.0, "lines")) + 
		theme(
			legend.title = element_blank(),
			legend.position = c(.99, .99),
			legend.justification = c("right", "top"),
			legend.box.just = "right",
			legend.key.height = unit(0.5, "cm"),
			legend.key.width = unit(0.5, "cm"),
			legend.text.align = 0, 
			legend.background = element_rect(fill=alpha('white', 0.4))) + 
		labs(
			title="",
			x    ="proportion variance explained (latent scale)")
}

# generate figure 5
fig5 <- subset(fig, scale == "latent scale")
fig5 <- droplevels(fig5)
p <- pltRep5(fig5)
(p)

# add capuchin sketch to figure 5
plt_rep <- cowplot::ggdraw() +
	draw_plot(p) +
	draw_image(plot(UTT), 
						 x = 0.622, y = 0.15, 
						 width = 0.4/1.25, height = 0.3/1.25)

width = 3.46457
pdf("plots/Fig5_Repeatability_latentNew.pdf", 
		width = width, height = width)
plt_rep
dev.off()

###############################################################################
###############################################################################
############ Figure 6: FE on RE ################################
###############################################################################
###############################################################################

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
		"latent.mode" = posterior.mode(m$VCV/rowSums(m$VCV)), #list, prop var
		"latent.HPDI95.lower" = HPDinterval(m$VCV/rowSums(m$VCV))[,1], #of prop var
		"latent.HPDI95.upper" = HPDinterval(m$VCV/rowSums(m$VCV))[,2], #of prop var
		row.names = 1:length(rlist)
	)
	df_tmp
}

genTabRE2 <- function(dre){
	dre$group <- substr(dre$Model,1,4)
	dre$rank <- substr(dre$Model,6,6)
	dre$rank[dre$rank==""] <- 5
	#relabel
	dre$model <- "blank"
	dre$model[dre$rank==0] <- "Intercept-only"
	dre$model[dre$rank==1] <- "Age"
	dre$model[dre$rank==2] <- "Age*Sex"
	dre$model[dre$rank==3] <- "Age*Sex, group size"
	dre$model[dre$rank==4] <- "Age*Sex, group size, seasonality"
	dre$model[dre$rank==5] <- "Age*Sex, group size, seasonality, ENSO"
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

###############################################################################

load("Models1_fullModels.RData")
load("Models2_reducedModels.RData")

dreM <- rbind(genTabRE1(m09a.0),
							genTabRE1(m09a.1),
							genTabRE1(m09a.2),
							genTabRE1(m09a.3),
							genTabRE1(m09a.4),
							genTabRE1(m09a),
							genTabRE1(m09b.0),
							genTabRE1(m09b.1),
							genTabRE1(m09b.2),
							genTabRE1(m09b.3),
							genTabRE1(m09b.4),
							genTabRE1(m09b),
							genTabRE1(m09c.0),
							genTabRE1(m09c.1),
							genTabRE1(m09c.2),
							genTabRE1(m09c.3),
							genTabRE1(m09c.4),
							genTabRE1(m09c)
							)
table_RE_M <- genTabRE2(dreM)

dreP <- rbind(genTabRE1(p09a.0),
							genTabRE1(p09a.1),
							genTabRE1(p09a.2),
							genTabRE1(p09a.3),
							genTabRE1(p09a.4),
							genTabRE1(p09a),
							genTabRE1(p09b.0),
							genTabRE1(p09b.1),
							genTabRE1(p09b.2),
							genTabRE1(p09b.3),
							genTabRE1(p09b.4),
							genTabRE1(p09b),
							genTabRE1(p09c.0),
							genTabRE1(p09c.1),
							genTabRE1(p09c.2),
							genTabRE1(p09c.3),
							genTabRE1(p09c.4),
							genTabRE1(p09c))
table_RE_P <- genTabRE2(dreP)

# plot of variance components
genPlotRE <- function(dre, title){
	plt <- ggplot(dre, aes(y = Component, 
												 x = latent.mode,  
												 colour = model),
								poition = position_dodge(padding=1.0)) + 
		coord_cartesian(xlim = c(0, 1.0)) + 
		scale_color_manual(values = clst) + 
		geom_point(position = position_dodge(width = 1.0), size = 0.75) + 
		geom_tile(color = "grey") + 
		xlab("proportion variance explained (latent scale)") + 
		ylab("") + 
		labs(color = "fixed effects") + 
		geom_linerange(aes(xmin = latent.HPDI95.lower, 
											 xmax = latent.HPDI95.upper, 
											 colour = model),
									 position = position_dodge(width = 1.0),
									 size = 0.5) + 
		theme_classic(base_size=7, base_family="Helvetica") + 
		theme(legend.position = c(0.7, 0.75), 
					legend.key.size = unit(0.3, 'cm'), 
					legend.text = element_text(size=5)) + 
		guides(color = guide_legend(override.aes = list(size=2),
																reverse=TRUE)) + 
		theme(
			panel.grid.major.y = element_blank(),
			panel.grid.major.x = element_blank(),
			panel.grid.minor.x = element_blank()) +
		ggtitle(title)
}

clst <- c(
	"Intercept-only" = "grey80",
	"Age" = "grey50", 
	"Age*Sex" = "grey30", 
	"Age*Sex, group size" = "orange1", 
	"Age*Sex, group size, seasonality" = "red",
	"Age*Sex, group size, seasonality, ENSO" = "red4")

pltM_a <- genPlotRE(subset(table_RE_M, 
													 group=='m09a'), 
										"a) social vs alone")
pltM_b <- genPlotRE(subset(table_RE_M, 
													 group=='m09b'), 
										"a) social vs alone")
pltM_c <- genPlotRE(subset(table_RE_M, 
													 group=='m09c'), 
										"a) social vs alone")

pltP_a <- genPlotRE(subset(table_RE_P, 
													 group=='p09a'), 
										"b) number of partners")
pltP_b <- genPlotRE(subset(table_RE_P, 
													 group=='p09b'), 
										"b) number of partners")
pltP_c <- genPlotRE(subset(table_RE_P, 
													 group=='p09c'),  
										"b) number of partners")

### add in DIC info

dre <- rbind(table_RE_M, table_RE_P)

dic <- dre %>% 
	distinct(Model, DIC, .keep_all = TRUE)
dic$type <- substr(dic$Model, 1, 1)
dic$prior <- substr(dic$Model, 4, 4)
dic <- dic[,c('Model','group','type','prior','rank','DIC','model')]

dic <- dic %>%
	group_by(group) %>%
	mutate(DIC_max = max(DIC))

dic$DIC_diff <- round(dic$DIC - dic$DIC_max, 0)

# function to generate plots of change in DIC scores
genPltDIC <- function(df, y1, y2){  #lolipop chart
	plt <- ggplot(df, aes(x=model, y=DIC_diff, color=model)) + 
		theme_classic(base_size = 7, base_family = "Helvetica") + 
		scale_color_manual(values = clst) + 
		geom_segment(aes(y = y2, 
										 x = model, 
										 yend = DIC_diff, 
										 xend = model), 
								 color = "black", alpha = 0.5, size = 0.25) + 
		geom_point(stat='identity', size=1)  +
		labs(title=expression(paste(Delta, " DIC")),
				 xlab="") + 
		ylim(y1, y2) +
		theme(legend.position = "none") + 
		theme(axis.title.y = element_blank(),
					axis.title.x = element_blank(),
					axis.text.y = element_blank(),
					axis.ticks.y = element_blank()) +  
		coord_flip()
	return(plt)
}

dic_sub <- subset(dic, type=='m' & prior=='a')
pltDIC_Ma <- genPltDIC(dic_sub, min(dic_sub$DIC_diff), max(dic_sub$DIC_diff))
(pltDIC_Ma)

dic_sub <- subset(dic, type=='m' & prior=='b')
pltDIC_Mb <- genPltDIC(dic_sub, min(dic_sub$DIC_diff), max(dic_sub$DIC_diff))
(pltDIC_Mb)

dic_sub <- subset(dic, type=='m' & prior=='c')
pltDIC_Mc <- genPltDIC(dic_sub, min(dic_sub$DIC_diff), max(dic_sub$DIC_diff))
(pltDIC_Mc)

dic_sub <- subset(dic, type=='p' & prior=='a')
pltDIC_Pa <- genPltDIC(dic_sub, min(dic_sub$DIC_diff), max(dic_sub$DIC_diff))
(pltDIC_Pa)

dic_sub <- subset(dic, type=='p' & prior=='b')
pltDIC_Pb <- genPltDIC(dic_sub, min(dic_sub$DIC_diff), max(dic_sub$DIC_diff))
(pltDIC_Pb)

dic_sub <- subset(dic, type=='p' & prior=='c')
pltDIC_Pc <- genPltDIC(dic_sub, min(dic_sub$DIC_diff), max(dic_sub$DIC_diff))
(pltDIC_Pc)

###

plt_a <- pltM_a + pltP_a
plt_a[[1]] = plt_a[[1]] + theme(legend.position = "none")
plt_a[[2]] = plt_a[[2]] + theme(axis.text.y = element_blank(),
														axis.title.y = element_blank())

plt_b <- pltM_b + pltP_b
plt_b[[1]] = plt_b[[1]] + theme(legend.position = "none")
plt_b[[2]] = plt_b[[2]] + theme(axis.text.y = element_blank(),
																axis.title.y = element_blank())

plt_c <- pltM_c + pltP_c
plt_c[[1]] = plt_c[[1]] + theme(legend.position = "none")
plt_c[[2]] = plt_c[[2]] + theme(axis.text.y = element_blank(),
																axis.title.y = element_blank())
###

plt2_a <- ggdraw() +
	draw_plot(plt_a) +
	draw_plot(pltDIC_Ma, 
						x = 0.4, y = 0.15, 
						width = 0.15, height = 0.3) +
	draw_plot(pltDIC_Pa, 
						x = 0.82, y = 0.15, 
						width = 0.15, height = 0.3)
(plt2_a)

width = 7.08661
pdf("plots/Fig6_FEonRE_a.pdf", 
		width = width, height = width/2)
plt2_a
dev.off()

plt2_b <- ggdraw() +
	draw_plot(plt_b) +
	draw_plot(pltDIC_Mb, 
						x = 0.4, y = 0.15, 
						width = 0.15, height = 0.3) +
	draw_plot(pltDIC_Pb, 
						x = 0.82, y = 0.15, 
						width = 0.15, height = 0.3)

pdf("supp/Fig6_FEonRE_b.pdf", 
		width = width, height = width/2)
plt2_b
dev.off()

plt2_c <- ggdraw() +
	draw_plot(plt_c) +
	draw_plot(pltDIC_Mc, 
						x = 0.4, y = 0.15, 
						width = 0.15, height = 0.3) +
	draw_plot(pltDIC_Pc, 
						x = 0.82, y = 0.15, 
						width = 0.15, height = 0.3)

pdf("supp/Fig6_FEonRE_c.pdf", 
		width = width, height = width/2)
plt2_c
dev.off()

##############################################################################
##############################################################################
##############################################################################

printEst <- function(m, col){
	print("====================================================")
	print(col)
	print(posterior.mode(m[["VCV"]][,col]/rowSums(m[["VCV"]])))
	print(HPDinterval(m[["VCV"]][,col]/rowSums(m[["VCV"]])))
}

m <- m09a[[1]]

for (col in colnames(m$VCV)){
	printEst(m, col)
}

###### Repeatability

vp <- rowSums(m[["VCV"]])
vi_within <- m[["VCV"]][,"id:Year"]+m[["VCV"]][,"Mother:Group:Year"]
posterior.mode(vi_within/vp)  #0.1202318 
HPDinterval(vi_within/vp)
#           lower    upper
# var1 0.09698419 0.143342
# attr(,"Probability")
# [1] 0.95

vi_long <- m[["VCV"]][,"id"]+m[["VCV"]][,"Mother"]+m[["VCV"]][,"animal"]
posterior.mode(vi_long/vp)  #0.2070075 
HPDinterval(vi_long/vp)
#          lower     upper
# var1 0.1686143 0.2654926
# attr(,"Probability")
# [1] 0.95

vi_short <- m[["VCV"]][,"id"]+m[["VCV"]][,"Mother"]+m[["VCV"]][,"animal"]+
	m[["VCV"]][,"id:Year"]+m[["VCV"]][,"Mother:Group:Year"]
posterior.mode(vi_short/vp)  #0.343289
HPDinterval(vi_short/vp)
#          lower     upper
# var1 0.2777952 0.3909005
# attr(,"Probability")
# [1] 0.95

m <- p09a[[1]]
vp <- rowSums(m[["VCV"]])

for (col in colnames(m$VCV)){
	printEst(m, col)
}

vi_within <- m[["VCV"]][,"id:Year"]+m[["VCV"]][,"Mother:Group:Year"]
posterior.mode(vi_within/vp)  #0.06384525 
HPDinterval(vi_within/vp)
#           lower      upper
# var1 0.05221002 0.07766866
# attr(,"Probability")
# [1] 0.95

vi_long <- m[["VCV"]][,"id"]+m[["VCV"]][,"Mother"]+m[["VCV"]][,"animal"]
posterior.mode(vi_long/vp)  #0.1436724 
HPDinterval(vi_long/vp)
#         lower     upper
# var1 0.112994 0.1810743
# attr(,"Probability")
# [1] 0.95

vi_short <- m[["VCV"]][,"id"]+m[["VCV"]][,"Mother"]+m[["VCV"]][,"animal"]+
	m[["VCV"]][,"id:Year"]+m[["VCV"]][,"Mother:Group:Year"]
posterior.mode(vi_short/vp)  #0.2126618 
HPDinterval(vi_short/vp)
#          lower     upper
# var1 0.1718745 0.2517314
# attr(,"Probability")
# [1] 0.95

###### group of residence

m <- m09a[[1]]
v_alpha <- m[["VCV"]][,"Group"]+m[["VCV"]][,"Group:Year"]
vp <- rowSums(m[["VCV"]])
posterior.mode(v_alpha/vp)  #0.1632801 
HPDinterval(v_alpha/vp)
#          lower     upper
# var1 0.1076107 0.2460715
# attr(,"Probability")
# [1] 0.95

m <- p09a[[1]]
v_alpha <- m[["VCV"]][,"Group"]+m[["VCV"]][,"Group:Year"]
vp <- rowSums(m[["VCV"]])
posterior.mode(v_alpha/vp)  #0.3317198 
HPDinterval(v_alpha/vp)
#          lower     upper
# var1 0.2448509 0.4363897
# attr(,"Probability")
# [1] 0.95

###### maternal

m <- m09a[[1]]
vm <- m[["VCV"]][,"Mother"]+m[["VCV"]][,"Mother:Group:Year"]
vp <- rowSums(m[["VCV"]])
posterior.mode(vm/vp)  #0.04546752 
HPDinterval(vm/vp)
#           lower     upper
# var1 0.02309197 0.0853172
# attr(,"Probability")
# [1] 0.95

m <- p09a[[1]]
vm <- m[["VCV"]][,"Mother"]+m[["VCV"]][,"Mother:Group:Year"]
vp <- rowSums(m[["VCV"]])
posterior.mode(vm/vp)  #0.04446521
HPDinterval(vm/vp)
#           lower      upper
# var1 0.02191045 0.07337441
# attr(,"Probability")
# [1] 0.95

######
