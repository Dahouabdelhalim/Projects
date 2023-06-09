# Data Analysis for MS, A multivariate test of disease risk reveals conditions leading to disease amplification

# File prepared by FH on 15 June 2017

# load required libraries ----
require(dplyr)
require(tidyr)
require(ggplot2)
require(grid)
require(gridExtra)
require(lsmeans)
require(gtable)
require(nlme)
require(MuMIn)
require(stringr)
require(piecewiseSEM)


# figure theme ----
theme_figs <- theme_classic() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) 
theme_set(theme_figs)
pd <- position_dodge(0.6)
pd2 <- position_dodge(width = 1)
update_geom_defaults('point', c(size = 2))
update_geom_defaults('errorbar', c(size = 0.5))



# read data ----
drdre.dat <- read.csv("DRE_Disease_DF.csv")

# metadata ----
# Plot ID identifies the unique ID for each plot
# Blk identifies experimental blocks (replicates)
# Plot identifies the combination of treatments for each plot
    # First five letters identifies the species that is either missing from polyculture or was planted in monoculture
    # The next letter identifies whether this plot is a polyculture (P) or a monoculture (M)
    # The next letter identifies whether it is ambient (O) or fertilized (R)
    # The next letter identifies the enemy exclusion treatment -- all plots are labeled O for no enemy exclusion
# Div identifies diversity treatment (M = 1 host species, P = 5 host species)
# Resource identifies fertility treatment (0 = ambient, 1 = fertilized)
# Load is the community weighted damage caused by all parasites
# RarS is the rarefied richness of all parasites
# load.mic is the community weighted damage caused by microbial parasites
# load.ins is the community weighted damage caused by insect parasites
# RarS.Mic is the rarefied richness of microbial parasites
# RarS.Ins is the rarefied richness of insect parasites

# Note that rarefaction was performed separately for different groups of parasites
  # therefore RarS.Mic and RarS.Ins do not add up to RarS



# Univariate analyses (Results Fig 2, Table S1) ----
dre.df <- drdre.dat %>% 
  mutate(ihs.load = asinh(load), #hyperbolic arcsine transformation (inverse hyperbolic sine)
         fblock = factor(Blk),
         fdiv = factor(ifelse(Div == "M", "Monoculture", "Polyculture")),
         fResource = factor(ifelse(Resource == 0, "Ambient", "Fertilized")),
         #identify different composition treatments for each richness 
         plant.trt = substr(Plot,1,5), 
         PlotID = factor(PlotID)) %>% data.frame()

# Univariate models not accounting for composition
out.rich <- gls(RarS ~ fblock + fResource*fdiv, data = dre.df, 
                weights = varIdent(form = ~ 1 | fdiv),
                control = list(maxIter=5000000, msMaxIter=500000), method = "ML")

out.abund <- gls(ihs.load ~ fblock + fResource*fdiv, data = dre.df, 
                 weights = varIdent(form = ~ 1 | fdiv),
                 control = list(maxIter=5000000, msMaxIter=500000), method = "ML")

# Univariate models acounting for composition
out.rich.comp <- lme(RarS ~ fblock + fResource*fdiv, data = dre.df,
                     random = ~ 1 | plant.trt/PlotID,
                     weights = varIdent(form = ~ 1 | fdiv),
                     control = list(maxIter=5000000, msMaxIter=500000), method = "ML")

out.abund.comp <- lme(ihs.load ~ fblock + fResource*fdiv, data = dre.df,
                      random = ~ 1 | plant.trt/PlotID,
                      weights = varIdent(form = ~ 1 | fdiv),
                      control = list(maxIter=5000000, msMaxIter=500000), method = "ML")

# Results presented in Table S1 ----
anova(out.abund)
anova(out.rich)

anova(out.abund.comp)
anova(out.rich.comp)

# plot results in Fig 2 ----
result.plot <- function(mod, response, comp) {
  
  if(response == "Parasite richness") {
    f.m <- function(x) {x}
  }
  
  else if (response == "Parasite abundance"){
    f.m <- function(x) {sinh(x)}
  }
  
  tmp <- lsmeans(mod, ~ fResource * fdiv) %>% summary() %>% as.data.frame() %>% 
    dplyr::mutate(est = f.m(lsmean),
                  up95 = f.m(upper.CL), 
                  low95 = f.m(lower.CL)) %>% 
    ggplot(., aes(x = fdiv, y = est, color = fResource, shape = fResource)) +
    geom_errorbar(aes(ymax = up95, ymin = low95), position = pd, width = 0, show.legend = FALSE) +
    geom_point(position = pd, fill = 'white', size = 2) +
    labs(y = ifelse(response == "Parasite richness", paste(response, "\\n\\n"),
                    paste(response, "\\n(% leaf area damaged) \\n")),
         x = "") + 
    ggtitle(paste("Single-response regression \\n", ifelse(comp == "yes", "accounting for composition\\n", "\\n"))) +
    scale_colour_manual(values = c('black', 'red')) + 
    theme(legend.position = c(.5,.95))
  
}

grid.arrange(ggplotGrob(result.plot(out.abund, "Parasite abundance", "no") + 
                          scale_y_continuous(limits = c(0,4.5), breaks = c(0,2,4), expand = c(0,0)) +
                          scale_x_discrete(labels = c("", "")) +
                          annotate("text", x = 0.6, y = 0.5, label = "A)")),
             ggplotGrob(result.plot(out.abund.comp, "Parasite abundance", "yes") + 
                          scale_y_continuous(limits = c(0,4.5), breaks = c(0,2,4),
                                             labels = c("", "", ""), expand = c(0,0)) +
                          theme( legend.position = "none") + labs(y = "\\n") +
                          scale_x_discrete(labels = c("", "")) +
                          annotate("text", x = 0.6, y = 0.5, label = "B)")),
             ggplotGrob(result.plot(out.rich, "Parasite richness", "no") + 
                          scale_y_continuous(limits = c(0,7), breaks = c(0,3,6), expand = c(0,0)) +
                          labs(y = "\\nRarefied parasite richness\\n") + ggtitle("") +
                          theme(legend.position = "none") +
                          annotate("text", x = 0.6, y = 0.5, label = "C)")),
             ggplotGrob(result.plot(out.rich.comp, "Parasite richness", "yes") + 
                          scale_y_continuous(limits = c(0,7), breaks = c(0,3,6),
                                                          labels = c("", "", ""), expand = c(0,0)) +
                          theme( legend.position = "none") + labs(y = "\\n") + ggtitle("") +
                          annotate("text", x = 0.6, y = 0.5, label = "D)")),
             ncol=2, nrow = 2)



# Multivariate analyses (Results Fig 3, Table S2, Table S3, Fig S1)----

dre.long <- dre.df %>% dplyr::select(PlotID, fResource, fblock, fdiv, load.ins, load.mic, RarS.Ins, RarS.Mic, plant.trt) %>% 
  # Create scaled variables for insect and microbial abundance and richness
  mutate(I.abundance = asinh(load.ins) * 100 / max(asinh(load.ins)), #scaled insect parasite abundance variable
         I.richness = RarS.Ins * 100 / max(RarS.Ins), #scaled insect parasite richness variable
         M.abundance = asinh(load.mic) * 100 /max(asinh(load.mic)),
         M.richness = RarS.Mic * 100 / max(RarS.Mic)) %>%  
  # Format data frame for multivariate analysis
  gather(key = Var, value = Y, I.abundance, I.richness, M.abundance, M.richness) %>% 
  mutate(Var = factor(Var),
         PlotID = factor(PlotID),
         VarNum = ifelse(Var == "I.abundance", 1,
                         ifelse(Var == "I.Richness", 2,
                                ifelse(Var == "M.abundance", 3, 4))))

# Multivariate response models
out.m2 <- gls(Y ~ fblock + Var + fResource:Var + fdiv:Var + fResource:fdiv:Var -1, 
                         data = dre.long, 
                         na.action = na.omit,
                         correlation = corSymm(form = ~ -1| PlotID),
                         weights = varComb(varIdent(form = ~ 1 | fdiv), varIdent(form = ~ -1 | VarNum)),
                         method = "ML")

out.m2.comp <- lme(Y ~ fblock + Var + fResource:Var + fdiv:Var + fResource:fdiv:Var -1, 
                              data = dre.long, 
                              na.action = na.omit,
                              random = ~ 1 | plant.trt/PlotID,
                              correlation = corSymm(form = ~ -1| plant.trt/PlotID),
                              weights = varComb(varIdent(form = ~ 1 | fdiv), varIdent(form = ~ -1 | VarNum)),
                              control = list(maxIter=5000000, msMaxIter=500000), method = "ML")

# Multivariate response ANOVAs for Table S2 ----
anova(out.m2)
anova(out.m2.comp)

# Simplify model by removing non-significant interaction between resources and diversity
out.m2b <- update(out.m2, .~. -fResource:fdiv:Var)
out.m2b.comp <- update(out.m2.comp, .~. -fResource:fdiv:Var)

# Estimated model terms for Table S3 ---- 
# coefficients - not accounting for composition
out.lsm.m2 <- lsmeans(out.m2b, ~ fResource * fdiv * Var)
out.lsm.m2.comp <- lsmeans(out.m2b.comp, ~ fResource * fdiv * Var)

summary(out.lsm.m2)

# parwise tests - not controlling for composition
# note alpha for significance should be interpreted at 0.0125 following bonferroni correction
pairs(out.lsm.m2, adjust = "none")
pairs(out.lsm.m2.comp, adjust = "none")

# Compare model predictions within a response
cld(out.lsm.m2, adjust = "none", alpha = 0.0125, Letters = "12345678", by = "Var")

# Comparisons for each treatment, only across responses
cld(out.lsm.m2, adjust = "none", alpha = 0.0125, Letters = "abcdefghijklmn", by = c("fResource", "fdiv"))

# extract correlation structure among response variables
# 1 = ins abund, 2 = ins rich, 3 = microbe abund, 4 = microbe rich
summary(out.m2b)$modelStruct$corStruct


# extract variance covariance matrix
VClist <- list()
for(n in 1:length(unique(dre.long$PlotID))){
  tmp <- getVarCov(out.m2b, individual = dre.long$PlotID[n])[1:4,1:4]
  VClist[[n]] <- tmp
}
temp_array <- abind::abind(VClist, along = 3)
varcov <- apply(temp_array, 1:2, mean)
varcov

# Construct multivariate response figure from reduced model (Fig 3)----
# Extract pseudo R-squared
r2.3mi2.var <- piecewiseSEM::sem.model.fits(out.m2b)[5]
lb3mi2.var <- paste("Marginal~R^{2}  == ", round(r2.3mi2.var,4))

r2.3mi2.var.comp <- r.squaredGLMM(out.m2b.comp)[1]
lb3mi2.var.comp <- paste("Marginal~R^{2}  == ", round(r2.3mi2.var.comp,4))

r2.3c2.var.comp <- r.squaredGLMM(out.m2b.comp)[2]
lb3c2.var.comp <- paste("Conditional~R^{2}  == ", round(r2.3c2.var.comp,4))

# placement for pseudo R-squared on figures
ann_text_m4 <- data.frame(fdiv = "Monoculture", fResource = "Ambient",
                          var = factor("Richness: Microbes",levels = c("Abundance: Microbes","Abundance: Insects", "Richness: Microbes", "Richness: Insects")),
                          est = 1)
ann_text_c4 <- ann_text_m4 %>% mutate(est = -2)


# not accounting for composition
m3.plot <- lsmeans(out.m2b, ~ fResource * fdiv * Var) %>% 
  cld(Letters = "abcdefghijklmn", adjust = "none", alpha = 0.0125) %>%
  dplyr::mutate(.group = stringr::str_trim(.group), # removes white space
                est = lsmean, up95 = upper.CL, low95 = lower.CL,
                var = ifelse(Var == "M.abundance", "Abundance: Microbes",
                             ifelse(Var == "M.richness", "Richness: Microbes",
                                    ifelse(Var == "I.abundance", "Abundance: Insects", "Richness: Insects")))) %>% 
  ggplot(., aes(x = fdiv, y = est, color = fResource, shape = fResource)) +
  facet_wrap(~var, nrow = 1, strip.position = "bottom")+
  geom_errorbar(aes(ymax = up95, ymin = low95), position = pd, width = 0, show.legend = FALSE) +
  geom_point(position = pd, fill = "white", size = 3) +
  labs(y = 'Standardized abundance or richness (%) \\n', x = '\\n Treatment') + 
  ggtitle("Multivariate response regression \\n\\n")+
  geom_text(aes(label = .group, x = fdiv, group = fResource, y = upper.CL + 3), position = pd, show.legend = FALSE) +
  scale_colour_manual(values = c('black', 'red')) + 
  theme(legend.position = c(.15,.85),
        plot.title = element_text(hjust = 0.01,
                                  margin=margin(b = -20, unit = "pt")),
        strip.placement = "outside",
        strip.text.x = element_text(hjust = 0.5, face = "bold"),
        panel.spacing.x = unit(0.5, 'cm')) +
  geom_text(data = ann_text_m4, label = lb3mi2.var, parse = TRUE, nudge_x = .4, nudge_y = 4, show.legend = FALSE) + 
  geom_vline(xintercept = 0.4, color = "black", size = 0.5, alpha = 0.5) + 
  annotate(geom="segment", y=c(0,20,40,60,80), yend = c(0,20,40,60,80),
           x=0.4, xend= 0.45, color = "grey20", size = 0.5, alpha = 0.5) +
  scale_y_continuous(expand = c(0,0))

# accounting for composition
m4.plot <- lsmeans(out.m2b.comp, ~ fResource * fdiv * Var) %>% 
  cld(Letters = "abcdefghijklmnopqrst", adjust = "none", alpha = 0.0125) %>%
  dplyr::mutate(.group = stringr::str_trim(.group), # removes white space
                est = lsmean, up95 = upper.CL, low95 = lower.CL,
                var = ifelse(Var == "M.abundance", "Abundance: Microbes",
                             ifelse(Var == "M.richness", "Richness: Microbes",
                                    ifelse(Var == "I.abundance", "Abundance: Insects", "Richness: Insects"))),
                fdiv = ifelse(fdiv == "Monoculture", "Mono", "Poly")) %>% 
  ggplot(., aes(x = fdiv, y = est, color = fResource, shape = fResource)) +
  facet_wrap(~var, nrow = 1, strip.position = "bottom")+
  geom_errorbar(aes(ymax = up95, ymin = low95), position = pd, width = 0) +
  geom_point(position = pd, fill = "white", size = 3) +
  labs(y = 'Standardized abundance or richness (%) \\n', x = '') + 
  ggtitle("Multivariate response regression accounting for composition \\n\\n")+
  geom_text(aes(label = .group, x = fdiv, group = fResource, y = upper.CL + 3), position = pd, show.legend = FALSE) +
  scale_colour_manual(values = c('black', 'red')) + 
  theme(legend.position = c(.15,.95),
        plot.title = element_text(hjust = 0.01,
                                  margin=margin(b = -29, unit = "pt")),
        panel.spacing.x = unit(0.5, 'cm'),
        strip.placement = "outside",
        strip.text.x = element_text(hjust = 0.5, face = "bold")) + 
  geom_text(data = ann_text_m4 %>% 
              mutate(fdiv = "Mono"), label = lb3mi2.var.comp, parse = TRUE, nudge_x = .4, nudge_y = 6, show.legend = FALSE) + 
  geom_text(data = ann_text_c4 %>% 
              mutate(fdiv = "Mono"),
            label = lb3c2.var.comp, parse = TRUE, show.legend = FALSE, nudge_x = .5, nudge_y = 4) + 
  geom_vline(xintercept = 0.4, color = "black", size = 0.5, alpha = 0.5) + 
  annotate(geom="segment", y=c(0,20,40,60,80), yend = c(0,20,40,60,80),
           x=0.4, xend= 0.45, color = "grey20", size = 0.5, alpha = 0.5) +
  scale_y_continuous(expand = c(0,0))

grid.arrange(m3.plot + scale_x_discrete(labels = c("","")) +
               theme(strip.text.x = element_blank()) +
               labs(x = ""),
             m4.plot + theme(legend.position = "none"), nrow=2)


# Construct multivariate response figure from full model (Fig S1)----
# extract pseudo r2
r2.3mi2.var.full <- piecewiseSEM::sem.model.fits(out.m2)[5]
lb3mi2.var.full <- paste("Marginal~R^{2}  == ", round(r2.3mi2.var.full,4))

r2.3mi2.var.comp.full <- r.squaredGLMM(out.m2.comp)[1]
lb3mi2.var.comp.full <- paste("Marginal~R^{2}  == ", round(r2.3mi2.var.comp.full,4))

r2.3c2.var.comp.full <- r.squaredGLMM(out.m2.comp)[2]
lb3c2.var.comp.full <- paste("Conditional~R^{2}  == ", round(r2.3c2.var.comp.full,4))

# Not accounting for composition
m3.plot.full <- cld(lsmeans(out.m2, ~ fResource * fdiv * Var), Letters = "abcdefghijklmn", adjust = "none", alpha = 0.0125) %>%
  dplyr::mutate(.group = stringr::str_trim(.group), # removes white space
                est = lsmean, up95 = upper.CL, low95 = lower.CL,
                var = ifelse(Var == "M.abundance", "Abundance: Microbes",
                             ifelse(Var == "M.richness", "Richness: Microbes",
                                    ifelse(Var == "I.abundance", "Abundance: Insects", "Richness: Insects")))) %>% 
  ggplot(., aes(x = fdiv, y = est, color = fResource, shape = fResource)) +
  facet_wrap(~var, nrow = 1, strip.position = "bottom")+
  geom_errorbar(aes(ymax = up95, ymin = low95), position = pd2, width = 0, show.legend = FALSE) +
  geom_point(position = pd2, fill = "white", size = 3) +
  labs(y = 'Standardized abundance or richness (%) \\n', x = '\\n Treatment') + 
  ggtitle("Multivariate response regression \\n\\n")+
  geom_text(aes(label = .group, x = fdiv, group = fResource, y = upper.CL + 3), position = pd2, show.legend = FALSE) +
  scale_colour_manual(values = c('black', 'red')) + 
  theme(legend.position = c(.15,.85),
        plot.title = element_text(hjust = 0.01,
                                  margin=margin(b = -20, unit = "pt")),
        strip.placement = "outside",
        strip.text.x = element_text(hjust = 0.5),
        panel.spacing.x = unit(-0.5, 'cm')) +
  geom_text(data = ann_text_m4, label = lb3mi2.var.full, parse = TRUE, nudge_x = .4, nudge_y = 4, show.legend = FALSE) + 
  geom_vline(xintercept = 0.4, color = "black", size = 0.5, alpha = 0.5) + 
  annotate(geom="segment", y=c(0,20,40,60,80), yend = c(0,20,40,60,80),
           x=0.4, xend= 0.45, color = "grey20", size = 0.5, alpha = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(NA,89))

# Accounting for composition
m4.plot.full <- cld(lsmeans(out.m2.comp, ~ fResource * fdiv * Var), Letters = "abcdefghijklmnopqrst", adjust = "none", alpha = 0.0125) %>%
  dplyr::mutate(.group = stringr::str_trim(.group), # removes white space
                est = lsmean, up95 = upper.CL, low95 = lower.CL,
                var = ifelse(Var == "M.abundance", "Abundance: Microbes",
                             ifelse(Var == "M.richness", "Richness: Microbes",
                                    ifelse(Var == "I.abundance", "Abundance: Insects", "Richness: Insects"))),
                fdiv = ifelse(fdiv == "Monoculture", "Mono", "Poly")) %>% 
  ggplot(., aes(x = fdiv, y = est, color = fResource, shape = fResource)) +
  facet_wrap(~var, nrow = 1, strip.position = "bottom")+
  geom_errorbar(aes(ymax = up95, ymin = low95), position = pd2, width = 0) +
  geom_point(position = pd2, fill = "white", size = 3) +
  labs(y = 'Standardized abundance or richness (%) \\n', x = '') + 
  ggtitle("Multivariate response regression accounting for composition \\n\\n")+
  geom_text(aes(label = .group, x = fdiv, group = fResource, y = upper.CL + 3), position = pd2, show.legend = FALSE) +
  scale_colour_manual(values = c('black', 'red')) + 
  theme(legend.position = c(.15,.95),
        plot.title = element_text(hjust = 0.01,
                                  margin=margin(b = -29, unit = "pt")),
        panel.spacing.x = unit(-0.5, 'cm'),
        strip.placement = "outside",
        strip.text.x = element_text(hjust = 0.5)) + 
  geom_text(data = ann_text_m4 %>% 
              mutate(fdiv = "Mono"), label = lb3mi2.var.comp.full, parse = TRUE, nudge_x = .4, nudge_y = 6, show.legend = FALSE) + 
  geom_text(data = ann_text_c4 %>% 
              mutate(fdiv = "Mono"),
            label = lb3c2.var.comp.full, parse = TRUE, show.legend = FALSE, nudge_x = .5, nudge_y = 4) + 
  geom_vline(xintercept = 0.4, color = "black", size = 0.5, alpha = 0.5) + 
  annotate(geom="segment", y=c(0,20,40,60,80), yend = c(0,20,40,60,80),
           x=0.4, xend= 0.45, color = "grey20", size = 0.5, alpha = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(NA,89))

    grid.arrange(m3.plot.full + scale_x_discrete(labels = c("","")) +
                   theme(strip.text.x = element_blank()) +
                   labs(x = ""),
                 m4.plot.full + theme(legend.position = "none"), nrow=2)
