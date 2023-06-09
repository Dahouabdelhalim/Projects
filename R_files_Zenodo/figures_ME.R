########################################################################################

packages <- c("ggplot2", "dplyr", "lavaan", "plyr", "cowplot", "rmarkdown", 
              "readr", "caTools", "bitops")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library(gridExtra)
library(cowplot)
library(dplyr)
library(readr)
library(ggpubr)
setwd("W:\\\\Ph_D\\\\RainCloudPlots-master\\\\tutorial_R")
getwd()

source("R_rainclouds.R")
source("summarySE.R")
source("simulateData.R")

raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))


################################################
# phi_{ad}
################################################
n <- length(sim.results[1,1,,1])
cmr <- cbind(phi.results[2,1,,1], rep('c',n), rep('0.1',n))
ipm1 <- cbind(sim.results.omega[2,1,,1], rep('b',n), rep('0.1',n))
ipm2 <- cbind(sim.results[2,1,,1], rep('a',n), rep('0.1',n))
dat1 <- data.frame(rbind(cmr,ipm1,ipm2))

cmr <- cbind(phi.results[2,2,,1], rep('c',n), rep('0.0',n))
ipm1 <- cbind(sim.results.omega[2,2,,1], rep('b',n), rep('0.0',n))
ipm2 <- cbind(sim.results[2,2,,1], rep('a',n), rep('0.0',n))
dat2 <- data.frame(rbind(cmr,ipm1,ipm2))
dat <- rbind(dat1,dat2)
names(dat) <- c('est','Model','sim')

dat$est <- as.numeric(as.character(dat$est))
dat$est <- format(round(dat$est, 4), nsmall = 4)
dat$est <- as.numeric(as.character(dat$est))
str(dat)


################################################
# plot with simulated data as an example for marker loss
################################################
labs <- NULL

labs[1] <- expression(IPM['no-imm'])
labs[2] <- expression(IPM['imm'])
labs[3] <- expression('CMR')
phi.ad <- ggplot(dat, aes(x = sim, y = est)) +
  geom_flat_violin(aes(fill = Model),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(sim)-.15, y = est, colour = Model),position = position_jitter(width = .05), size = 3, shape = 20)+
  geom_boxplot(aes(x = sim, y = est, fill = Model),
               outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  # scale_colour_brewer(palette = "Dark2")+
  theme(text = element_text(size = 16, family = 'serif'), axis.title.x = element_blank(),
        legend.position = 'none')+
  # scale_fill_brewer(palette = "Dark2") + 
  ylab(expression(phi['ad'])) +
  geom_hline(yintercept = 0.5, linetype = 'dotted') +
  scale_colour_manual(values = c('a' = 'red', 'b' = 'green', 'c' = 'blue'),
                          labels = c(labs[1],labs[2],labs[3])) +
  guides(fill = FALSE)
phi.ad


















################################################
# phi_{juv}
################################################
n <- length(sim.results[1,1,,1])
cmr <- cbind(phi.results[1,1,,1], rep('c',n), rep('0.1',n))
ipm1 <- cbind(sim.results.omega[1,1,,1], rep('b',n), rep('0.1',n))
ipm2 <- cbind(sim.results[1,1,,1], rep('a',n), rep('0.1',n))
dat1 <- data.frame(rbind(cmr,ipm1,ipm2))

cmr <- cbind(phi.results[1,2,,1], rep('c',n), rep('0.0',n))
ipm1 <- cbind(sim.results.omega[1,2,,1], rep('b',n), rep('0.0',n))
ipm2 <- cbind(sim.results[1,2,,1], rep('a',n), rep('0.0',n))
dat2 <- data.frame(rbind(cmr,ipm1,ipm2))
dat <- rbind(dat1,dat2)
names(dat) <- c('est','Model','sim')

dat$est <- as.numeric(as.character(dat$est))
dat$est <- format(round(dat$est, 4), nsmall = 4)
dat$est <- as.numeric(as.character(dat$est))
str(dat)


################################################
# plot with simulated data as an example for marker loss
################################################
labs <- NULL
labs[1] <- expression(IPM['no-imm'])
labs[2] <- expression(IPM['imm'])
labs[3] <- expression('CMR')
phi.juv <- ggplot(dat, aes(x = sim, y = est)) +
  geom_flat_violin(aes(fill = Model),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(sim)-.15, y = est, colour = Model),position = position_jitter(width = .05), size = 3, shape = 20)+
  geom_boxplot(aes(x = sim, y = est, fill = Model),
               outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  # scale_colour_brewer(palette = "Dark2")+
  theme(text = element_text(size = 16, family = 'serif'), axis.title.x = element_blank(),
        legend.position = 'none')+
  # scale_fill_brewer(palette = "Dark2") + 
  ylab(expression(phi['juv'])) +
  geom_hline(yintercept = 0.3, linetype = 'dotted') +
  scale_colour_manual(values = c('a' = 'red', 'b' = 'green', 'c' = 'blue'),
                      labels = c(labs[1],labs[2],labs[3])) +
  guides(fill = FALSE)
phi.juv













################################################
# p
################################################
n <- length(sim.results[1,1,,1])
cmr <- cbind(phi.results[3,1,,1], rep('c',n), rep('0.1',n))
ipm1 <- cbind(sim.results.omega[4,1,,1], rep('b',n), rep('0.1',n))
ipm2 <- cbind(sim.results[4,1,,1], rep('a',n), rep('0.1',n))
dat1 <- data.frame(rbind(cmr,ipm1,ipm2))

cmr <- cbind(phi.results[3,2,,1], rep('c',n), rep('0.0',n))
ipm1 <- cbind(sim.results.omega[4,2,,1], rep('b',n), rep('0.0',n))
ipm2 <- cbind(sim.results[4,2,,1], rep('a',n), rep('0.0',n))
dat2 <- data.frame(rbind(cmr,ipm1,ipm2))
dat <- rbind(dat1,dat2)
names(dat) <- c('est','Model','sim')

dat$est <- as.numeric(as.character(dat$est))
dat$est <- format(round(dat$est, 4), nsmall = 4)
dat$est <- as.numeric(as.character(dat$est))
str(dat)


################################################
# plot with simulated data as an example for marker loss
################################################
labs <- NULL
labs[1] <- expression(IPM['no-imm'])
labs[2] <- expression(IPM['imm'])
labs[3] <- expression('CMR')
p <- ggplot(dat, aes(x = sim, y = est)) +
  geom_flat_violin(aes(fill = Model),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(sim)-.15, y = est, colour = Model),position = position_jitter(width = .05), size = 3, shape = 20)+
  geom_boxplot(aes(x = sim, y = est, fill = Model),
               outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  # scale_colour_brewer(palette = "Dark2")+
  theme(text = element_text(size = 16, family = 'serif'), axis.title.x = element_blank(),
        legend.position = 'none')+
  # scale_fill_brewer(palette = "Dark2") + 
  ylab(expression(p)) +
  geom_hline(yintercept = 0.6, linetype = 'dotted') +
  scale_colour_manual(values = c('a' = 'red', 'b' = 'green', 'c' = 'blue'),
                      labels = c(labs[1],labs[2],labs[3])) +
  guides(fill = FALSE)
p











################################################
# sigma.y
################################################
n <- length(sim.results[1,1,,1])
cmr <- cbind(ssm.results[1,1,,1], rep('c',n), rep('0.1',n))
ipm1 <- cbind(sim.results.omega[5,1,,1], rep('b',n), rep('0.1',n))
ipm2 <- cbind(sim.results[5,1,,1], rep('a',n), rep('0.1',n))
dat1 <- data.frame(rbind(cmr,ipm1,ipm2))

cmr <- cbind(ssm.results[1,2,,1], rep('c',n), rep('0.0',n))
ipm1 <- cbind(sim.results.omega[5,2,,1], rep('b',n), rep('0.0',n))
ipm2 <- cbind(sim.results[5,2,,1], rep('a',n), rep('0.0',n))
dat2 <- data.frame(rbind(cmr,ipm1,ipm2))
dat <- rbind(dat1,dat2)
names(dat) <- c('est','Model','sim')

dat$est <- as.numeric(as.character(dat$est))
dat$est <- format(round(dat$est, 4), nsmall = 4)
dat$est <- as.numeric(as.character(dat$est))
str(dat)


################################################
# plot with simulated data as an example for marker loss
################################################
labs <- NULL
labs[1] <- expression(IPM['no-imm'])
labs[2] <- expression(IPM['imm'])
labs[3] <- expression('SSM')
sig.y <- ggplot(dat, aes(x = sim, y = est)) +
  geom_flat_violin(aes(fill = Model),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(sim)-.15, y = est, colour = Model),position = position_jitter(width = .05), size = 3, shape = 20)+
  geom_boxplot(aes(x = sim, y = est, fill = Model),
               outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  # scale_colour_brewer(palette = "Dark2")+
  theme(text = element_text(size = 16, family = 'serif'), axis.title.x = element_blank(),
        legend.position = 'none')+
  # scale_fill_brewer(palette = "Dark2") + 
  ylab(expression(sigma['y'])) +
  geom_hline(yintercept = 100, linetype = 'dotted') +
  scale_colour_manual(values = c('a' = 'red', 'b' = 'green', 'c' = 'blue'),
                      labels = c(labs[1],labs[2],labs[3])) +
  guides(fill = FALSE) +
  ylim(0,175)
sig.y

# boxplot(sim.results.omega[5,1,,1], sim.results[5,1,,1], ssm.results[1,1,,1])














################################################
# fecundity
################################################
n <- length(sim.results[1,1,,1])
cmr <- cbind(fec.results[1,1,,1], rep('c',n), rep('0.1',n))
ipm1 <- cbind(sim.results.omega[3,1,,1], rep('b',n), rep('0.1',n))
ipm2 <- cbind(sim.results[3,1,,1], rep('a',n), rep('0.1',n))
dat1 <- data.frame(rbind(cmr,ipm1,ipm2))

cmr <- cbind(fec.results[1,2,,1], rep('c',n), rep('0.0',n))
ipm1 <- cbind(sim.results.omega[3,2,,1], rep('b',n), rep('0.0',n))
ipm2 <- cbind(sim.results[3,2,,1], rep('a',n), rep('0.0',n))
dat2 <- data.frame(rbind(cmr,ipm1,ipm2))
dat <- rbind(dat1,dat2)
names(dat) <- c('est','Model','sim')

dat$est <- as.numeric(as.character(dat$est))
dat$est <- format(round(dat$est, 2), nsmall = 4)
dat$est <- as.numeric(as.character(dat$est))
str(dat)


################################################
# plot with simulated data as an example for marker loss
################################################
labs <- NULL
labs[1] <- expression(IPM['no-imm'])
labs[2] <- expression(IPM['imm'])
labs[3] <- expression('PRM')
fec <- ggplot(dat, aes(x = sim, y = est)) +
  geom_flat_violin(aes(fill = Model),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(sim)-.15, y = est, colour = Model),position = position_jitter(width = .05), size = 3, shape = 20)+
  geom_boxplot(aes(x = sim, y = est, fill = Model),
               outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  # scale_colour_brewer(palette = "Dark2")+
  theme(text = element_text(size = 16, family = 'serif'), axis.title.x = element_blank(),
        legend.position = 'none')+
  # scale_fill_brewer(palette = "Dark2") + 
  ylab(expression(f)) +
  geom_hline(yintercept = 1.6, linetype = 'dotted') +
  scale_colour_manual(values = c('a' = 'red', 'b' = 'green', 'c' = 'blue'),
                      labels = c(labs[1],labs[2],labs[3])) +
  guides(fill = FALSE)
fec


















########################################################################
# immigration
########################################################################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 3
cols = gg_color_hue(n)

# plot(1:n, pch = 16, cex = 2, col = cols)

ipm <- cbind(sim.results.omega[6,1,,1], rep('a',n))
ssm <- cbind(sim.results.omega[6,2,,1], rep('b',n))
dat <- data.frame(rbind(ipm,ssm))
names(dat) <- c('est','Model')
dat$est <- as.numeric(as.character(dat$est))
labs <- NULL
labs[1] <- expression('0.1')
labs[2] <- expression('0.0')
imm <- ggplot(dat, aes(x = Model, y = est)) +
  geom_flat_violin(fill = cols[2],position = position_nudge(x = .1, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(Model)-.15, y = est), color = 'green',
             position = position_jitter(width = .05), size = 3, shape = 20)+
  geom_boxplot(aes(x = Model, y = est), fill = cols[2],
               outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  # scale_colour_brewer(palette = "Dark2")+
  theme(text = element_text(size = 16, family = 'serif'), axis.title.x = element_blank(),
        legend.position = 'none')+
  # scale_fill_brewer(palette = "Dark2") + 
  ylab(expression(omega)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  # scale_colour_manual(values = c('a' = 'red','b' = 'blue'),
  #                     labels = c(labs[1],labs[2])) +
  guides(fill = FALSE) +
  scale_x_discrete('', c('a','b'), c('0.1','0.0'))
imm








fig <- ggarrange(phi.ad,phi.juv,p,sig.y,fec,imm, nrow = 3, ncol = 2,
                 labels = c("A","B","C","D",'E','F'))
# fig

pdf("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\new_TEX\\\\Figure2.pdf", width=10, height=10)
annotate_figure(fig, bottom = text_grob('Marker Loss Rate', size = 16, family = 'serif'))
dev.off()











###################################################################
# estimate coverage
###################################################################
# \\phi_{ad}
length(which(sim.results[2,1,,2] < 0.5 & sim.results[2,1,,3] > 0.5))
length(which(sim.results.omega[2,1,,2] < 0.5 & sim.results.omega[2,1,,3] > 0.5))
length(which(phi.results[2,1,,2] < 0.5 & phi.results[2,1,,3] > 0.5))

length(which(sim.results[2,2,,2] < 0.5 & sim.results[2,2,,3] > 0.5))
length(which(sim.results.omega[2,2,,2] < 0.5 & sim.results.omega[2,2,,3] > 0.5))
length(which(phi.results[2,2,,2] < 0.5 & phi.results[2,2,,3] > 0.5))


# \\phi_{juv}
length(which(sim.results[1,1,,2] < 0.3 & sim.results[1,1,,3] > 0.3))
length(which(sim.results.omega[1,1,,2] < 0.3 & sim.results.omega[1,1,,3] > 0.3))
length(which(phi.results[1,1,,2] < 0.3 & phi.results[1,1,,3] > 0.3))

length(which(sim.results[1,2,,2] < 0.3 & sim.results[1,2,,3] > 0.3))
length(which(sim.results.omega[1,2,,2] < 0.3 & sim.results.omega[1,2,,3] > 0.3))
length(which(phi.results[1,2,,2] < 0.3 & phi.results[1,2,,3] > 0.3))


# p
length(which(sim.results[4,1,,2] < 0.6 & sim.results[4,1,,3] > 0.6))
length(which(sim.results.omega[4,1,,2] < 0.6 & sim.results.omega[4,1,,3] > 0.6))
length(which(phi.results[3,1,,2] < 0.6 & phi.results[3,1,,3] > 0.6))

length(which(sim.results[4,2,,2] < 0.6 & sim.results[4,2,,3] > 0.6))
length(which(sim.results.omega[4,2,,2] < 0.6 & sim.results.omega[4,2,,3] > 0.6))
length(which(phi.results[3,2,,2] < 0.6 & phi.results[3,2,,3] > 0.6))

# \\simga_{y}
length(which(sim.results[5,1,,2] < 100 & sim.results[5,1,,3] > 100))
length(which(sim.results.omega[5,1,,2] < 100 & sim.results.omega[5,1,,3] > 100))
length(which(ssm.results[1,1,,2] < 100 & ssm.results[1,1,,3] > 100))

length(which(sim.results[5,2,,2] < 100 & sim.results[5,2,,3] > 100))
length(which(sim.results.omega[5,2,,2] < 100 & sim.results.omega[5,2,,3] > 100))
length(which(ssm.results[1,2,,2] < 100 & ssm.results[1,2,,3] > 100))

# f
length(which(sim.results[3,1,,2] < 1.6 & sim.results[3,1,,3] > 1.6))
length(which(sim.results.omega[3,1,,2] < 1.6 & sim.results.omega[3,1,,3] > 1.6))
length(which(fec.results[1,1,,2] < 1.6 & fec.results[1,1,,3] > 1.6))

length(which(sim.results[3,2,,2] < 1.6 & sim.results[3,2,,3] > 1.6))
length(which(sim.results.omega[3,2,,2] < 1.6 & sim.results.omega[3,2,,3] > 1.6))
length(which(fec.results[1,2,,2] < 1.6 & fec.results[1,2,,3] > 1.6))


# \\omega
length(which(sim.results.omega[6,1,,2] < 0 & sim.results.omega[6,1,,3] > 0))
length(which(sim.results.omega[6,2,,2] < 0 & sim.results.omega[6,2,,3] > 0))








###################################################################
# determine relative/absolute difference
###################################################################
###################################################################
#  CJS
###################################################################
# juvenile survival: IPM_{no_imm}
mean((phi.results[1,1,,1] - 0.3)/0.3)
mean((phi.results[1,2,,1] - 0.3)/0.3)
mean(phi.results[1,1,,1] - 0.3)
mean(phi.results[1,2,,1] - 0.3)

# adult survival: IPM_{no_imm}
mean((phi.results[2,1,,1] - 0.5)/0.5)
mean((phi.results[2,2,,1] - 0.5)/0.5)
mean(phi.results[2,1,,1] - 0.5)
mean(phi.results[2,2,,1] - 0.5)

# p: IPM_{no_imm}
mean((phi.results[3,1,,1] - 0.6)/0.6)
mean((phi.results[3,2,,1] - 0.6)/0.6)
mean(phi.results[3,1,,1] - 0.6)
mean(phi.results[3,2,,1] - 0.6)


###################################################################
#  SSM
###################################################################
mean((ssm.results[1,1,,1] - 100)/100)
mean((ssm.results[1,2,,1] - 100)/100)
mean(ssm.results[1,1,,1] - 100)
mean(ssm.results[1,2,,1] - 100)



###################################################################
#  fec
###################################################################
mean((fec.results[1,1,,1] - 1.6)/1.6)
mean((fec.results[1,2,,1] - 1.6)/1.6)
mean(fec.results[1,1,,1] - 1.6)
mean(fec.results[1,2,,1] - 1.6)


###################################################################
#  IPM_{imm}
###################################################################
# juvenile survival: IPM_{no_imm}
mean((sim.results.omega[1,1,,1] - 0.3)/0.3)
mean((sim.results.omega[1,2,,1] - 0.3)/0.3)
mean(sim.results.omega[1,1,,1] - 0.3)
mean(sim.results.omega[1,2,,1] - 0.3)

# adult survival: IPM_{no_imm}
mean((sim.results.omega[2,1,,1] - 0.5)/0.5)
mean((sim.results.omega[2,2,,1] - 0.5)/0.5)
mean(sim.results.omega[2,1,,1] - 0.5)
mean(sim.results.omega[2,2,,1] - 0.5)

# p: IPM_{no_imm}
mean((sim.results.omega[4,1,,1] - 0.6)/0.6)
mean((sim.results.omega[4,2,,1] - 0.6)/0.6)
mean(sim.results.omega[4,1,,1] - 0.6)
mean(sim.results.omega[4,2,,1] - 0.6)

# sigma_y: IPM_{no_imm}
mean((sim.results.omega[5,1,,1] - 100)/100)
mean((sim.results.omega[5,2,,1] - 100)/100)
mean(sim.results.omega[5,1,,1] - 100)
mean(sim.results.omega[5,2,,1] - 100)

# fecundity: IPM_{no_imm}
mean((sim.results.omega[3,1,,1] - 1.6)/1.6)
mean((sim.results.omega[3,2,,1] - 1.6)/1.6)
mean(sim.results.omega[3,1,,1] - 1.6)
mean(sim.results.omega[3,2,,1] - 1.6)

# immigration
mean(sim.results.omega[6,1,,1])
mean(sim.results.omega[6,2,,1])








###################################################################
#  IPM_{no-imm}
###################################################################
# juvenile survival: IPM_{no_imm}
mean((sim.results[1,1,,1] - 0.3)/0.3)
mean((sim.results[1,2,,1] - 0.3)/0.3)
mean(sim.results[1,1,,1] - 0.3)
mean(sim.results[1,2,,1] - 0.3)

# adult survival: IPM_{no_imm}
mean((sim.results[2,1,,1] - 0.5)/0.5)
mean((sim.results[2,2,,1] - 0.5)/0.5)
mean(sim.results[2,1,,1] - 0.5)
mean(sim.results[2,2,,1] - 0.5)

# p: IPM_{no_imm}
mean((sim.results[4,1,,1] - 0.6)/0.6)
mean((sim.results[4,2,,1] - 0.6)/0.6)
mean(sim.results[4,1,,1] - 0.6)
mean(sim.results[4,2,,1] - 0.6)

# sigma_y: IPM_{no_imm}
mean((sim.results[5,1,,1] - 100)/100)
mean((sim.results[5,2,,1] - 100)/100)
mean(sim.results[5,1,,1] - 100)
mean(sim.results[5,2,,1] - 100)

# fecundity: IPM_{no_imm}
mean((sim.results[3,1,,1] - 1.6)/1.6)
mean((sim.results[3,2,,1] - 1.6)/1.6)
mean(sim.results[3,1,,1] - 1.6)
mean(sim.results[3,2,,1] - 1.6)








