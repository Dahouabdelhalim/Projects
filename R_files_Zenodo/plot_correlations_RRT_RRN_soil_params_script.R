### plot correlation between RRT and RRN with soil parameters

# This script creates Figure 5. 
# Relationships were tested based on indicated Spearman's coefficients
# of correlation according to the correlogram of Figure 4. 

library(tidyverse)
library(RColorBrewer)
library(ggpmisc)
library(ggpubr)
library(cowplot)
#--------------------------

# load script, which loads needed data
source("correlograms_RRT_RRN_soil_param_script.R")  

# formatter is a function to define decimal places in facets
# https://stackoverflow.com/questions/37669064/facet-wrap-title-wrapping-decimal-places-on-free-y-axis-ggplot2
formatter <- function(...){
  function(x) format(round(x, 1), ...)
}

# based on indicated Spearman correlations according to correlograms 
my.formula <- y~x
my.formula.log <- y ~ log(x)
my.formula.exp <- y ~ exp(x)
my.formula.pol <- y ~ poly(x, 2, raw = TRUE)

#### RRN ####

#### RRN Cgrowth ~ F:B ratio ####

fit.linear <- lm(Cgrowth_ngC_gDW_h ~ F.B, data = gensoilRRTN.sub.N)
fit.log <- lm(Cgrowth_ngC_gDW_h ~ log(F.B), data = gensoilRRTN.sub.N)
fit.exp <- lm(Cgrowth_ngC_gDW_h ~ exp(F.B), data = gensoilRRTN.sub.N)
fit.poly <- lm(Cgrowth_ngC_gDW_h ~ poly(F.B, 2, raw = TRUE), data = gensoilRRTN.sub.N)
anova.modelfit <- anova(fit.linear, fit.log, fit.exp, fit.poly)
anova.modelfit
AIC(fit.linear, fit.log, fit.exp, fit.poly) 
model.list <- list(fit.linear, fit.log, fit.exp, fit.poly)
lapply(model.list, FUN = summary)

RRNCgrowth.FB.fig <- 
  ggplot(data=gensoilRRTN.sub.N, aes(x=F.B, y=Cgrowth_ngC_gDW_h, fill = RR)) +
  geom_hline(yintercept=1, linetype="dotted") +
  stat_smooth(method="lm", formula = my.formula.log, se=TRUE, linetype = "dotted", color = "grey35") +
  geom_point(aes(shape=landuse), size = 2.5, color="grey35") +
  stat_poly_eq(formula = my.formula.log, 
               aes(label =  paste(
                 stat(adj.rr.label), 
                 stat(p.value.label), 
                 sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  theme_bw(base_size = 13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  xlab(expression("F:B"~ratio)) +
  ylab(expression(RR[Cgrowth])) +
  scale_fill_manual(values = c("#0073C2FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  geom_text(aes(label=ifelse(Cgrowth_ngC_gDW_h>1.2,as.character(sample_ID),'')),hjust=0,vjust=0) +
  ylim(0.5,1.3)  # one outlier removed
RRNCgrowth.FB.fig


#### RRN cumulresp ~ C:N ####
fit.1 <- lm(cumul_resp_incub_ngC_g_total_time ~ CN_ratio, data = gensoilRRTN.sub.N[gensoilRRTN.sub.N$cumul_resp_incub_ngC_g_total_time<=1.2,])
fit.2 <- lm(cumul_resp_incub_ngC_g_total_time ~ log(CN_ratio), data = gensoilRRTN.sub.N[gensoilRRTN.sub.N$cumul_resp_incub_ngC_g_total_time<=1.2,])
fit.3 <- lm(cumul_resp_incub_ngC_g_total_time ~ exp(CN_ratio), data = gensoilRRTN.sub.N[gensoilRRTN.sub.N$cumul_resp_incub_ngC_g_total_time<=1.2,])
fit.4 <- lm(cumul_resp_incub_ngC_g_total_time ~ poly(CN_ratio, 2, raw = TRUE), data = gensoilRRTN.sub.N[gensoilRRTN.sub.N$cumul_resp_incub_ngC_g_total_time<=1.2,])
anova(fit.1, fit.2, fit.3, fit.4)
AIC(fit.1, fit.2, fit.3, fit.4)

RRNcumulresp.CN.fig <- 
  ggplot(data=gensoilRRTN.sub.N, aes(x=CN_ratio, y=cumul_resp_incub_ngC_g_total_time, fill = RR)) +
  geom_hline(yintercept=1, linetype="dotted") +
  stat_smooth(method="lm", formula = my.formula.log, se=TRUE, linetype = "dotted", color = "grey35") +
  geom_point(aes(shape=landuse), size = 2.5, color="grey35") +
  stat_poly_eq(formula = my.formula.log, 
               aes(label =  paste(
                 stat(adj.rr.label), 
                 stat(p.value.label), 
                 sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  theme_bw(base_size = 13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom") +
  guides(color="none") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  xlab(expression("C:N"~ratio)) +
  ylab(expression(RR[cumulative~respiration])) +
  scale_fill_manual(values = c("#0073C2FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  geom_text(aes(label=ifelse(cumul_resp_incub_ngC_g_total_time>1.2,as.character(sample_ID),'')),hjust=0,vjust=0) +
  ylim(0.55,1.2) # two outliers removed


#### Figure 5 - drivers of RRN ####
empty.plot <- ggplot()
legend <- get_legend(RRNCgrowth.FB.fig + theme(legend.position = "bottom") + guides(fill = "none"))
plot_N <- plot_grid(RRNcumulresp.CN.fig +
                         theme(legend.position = "none",
                               axis.title.y = element_text(size=13),
                               axis.title.x = element_text(size=13)),
                       empty.plot,
                       RRNCgrowth.FB.fig +
                         theme(legend.position = "none",
                               axis.title.y = element_text(size=13),
                               axis.title.x = element_text(size=13)),
                       labels = c("a","","b"), 
                       vjust = 1,
                       ncol = 1,
                       nrow = 3,
                       align = "v", 
                       axis = "rl",
                       rel_heights = c(1,0.05,1))

plot_grid(plot_N,
          legend,
          nrow = 2,
          rel_heights = c(1, 0.1))





#### RRT ####

#### RRT cumulresp ~ Fungi ####
fit.1 <- lm(cumul_resp_incub_ngC_g_total_time ~ F_cn_gsoil, data = gensoilRRTN.sub.T)
fit.2 <- lm(cumul_resp_incub_ngC_g_total_time ~ log(F_cn_gsoil), data = gensoilRRTN.sub.T)
#fit.3 <- lm(cumul_resp_incub_ngC_g_total_time ~ exp(F_cn_gsoil), data = gensoilRRTN.sub.T)
fit.4 <- lm(cumul_resp_incub_ngC_g_total_time ~ poly(F_cn_gsoil, 2, raw = TRUE), data = gensoilRRTN.sub.T)
anova(fit.1, fit.2, fit.4)
AIC(fit.1, fit.2, fit.4)

RRTcumulresp.F.fig <- ggplot(data=gensoilRRTN.sub.T, aes(x=F_cn_gsoil, y=cumul_resp_incub_ngC_g_total_time, fill = RR)) +
  geom_hline(yintercept=1, linetype="dotted") +
  stat_smooth(method="lm", formula = my.formula.log, se=TRUE, linetype = "dotted", color = "grey35") +
  geom_point(aes(shape=landuse), size = 2.5, color = "grey35") +
  stat_poly_eq(formula = my.formula.log, 
               aes(label =  paste(
                 stat(adj.rr.label), 
                 stat(p.value.label), 
                 sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  theme_bw(base_size = 13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom") +
  guides(fill="none") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  xlab(expression(Fungi~"["~number~of~gene~copies~g^-1~soil~DW~"]")) +
  ylab(expression(RR[cumulative~respiration])) +
  scale_fill_manual(values = c("#EFC000FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  geom_text(aes(label=ifelse(cumul_resp_incub_ngC_g_total_time<1.2,as.character(sample_ID),'')),hjust=1.1,vjust=0) #+
  #ylim(1.1,2.5) # to remove one outlier
RRTcumulresp.F.fig

# RRT CUE over CN soil   # does not reveal interesting correlation
# RRT CUE over Ntotal    # does not reveal interesting correlation



#### RRT ~ clay content ####

# Interesting corrleations as seen from correlograms:
# RRT CUE, Cgrowth ~ clay
# Correlations were driven by one site and it was concluded, that they are not trustworthy.

#### RRT CUE ~ clay ####
fit.1 <- lm(CUE_mean ~ clay_perc, data = compsoilmeanRRT)
fit.2 <- lm(CUE_mean ~ log(clay_perc), data = compsoilmeanRRT)
fit.3 <- lm(CUE_mean ~ exp(clay_perc), data = compsoilmeanRRT)
fit.4 <- lm(CUE_mean ~ poly(clay_perc, 2, raw = TRUE), data = compsoilmeanRRT)
anova(fit.1, fit.2, fit.3, fit.4)
AIC(fit.1, fit.2, fit.3, fit.4) 
model.list <- list(fit.1, fit.2, fit.3, fit.4)
lapply(model.list, FUN = summary)

RRCUEclay.fig <- ggplot(data=compsoilmeanRRT, aes(x=clay_perc, y=CUE_mean, fill = RR)) +
  geom_errorbar(aes(ymin=CUE_mean-CUE_sd, 
                    ymax=CUE_mean+CUE_sd),
                width=0, color = "grey35") +
  stat_smooth(method="lm", formula = my.formula, se=TRUE, linetype = "dotted", color = "grey35") +
  geom_point(aes(shape=landuse, fill=site), size = 2.5, color = "grey35") +
  stat_poly_eq(formula = my.formula, 
               aes(label =  paste(
                 stat(adj.rr.label), 
                 stat(p.value.label), 
                 sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(clay~"["~"%"~"]")) +
  ylab(expression(RR[CUE])) +
  scale_fill_manual(values = c("#8F7700FF","#999999","#EFC000FF","#EFC000FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  facet_wrap(.~RR, scales = "free")
RRCUEclay.fig 

# Exponential fit shows that the trend is mostly driven by site CD
RRCUEclay.fig.exp <- ggplot(data=compsoilmeanRRT, aes(x=clay_perc, y=CUE_mean, fill = RR)) + 
  geom_errorbar(aes(ymin=CUE_mean-CUE_sd, 
                    ymax=CUE_mean+CUE_sd),
                width=0, color = "grey35") +
  stat_smooth(method="lm", formula = my.formula.exp, se=TRUE, linetype = "dotted", color = "grey35") +
  geom_point(aes(shape=landuse, fill=site), size = 2.5, color = "grey35") +
  stat_poly_eq(formula = my.formula.exp, 
               aes(label =  paste(
                 stat(adj.rr.label), 
                 stat(p.value.label), 
                 sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(clay~"["~"%"~"]")) +
  ylab(expression(RR[CUE])) +
  scale_fill_manual(values = c("#8F7700FF","#999999","#EFC000FF","#EFC000FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  facet_wrap(.~RR, scales = "free")
RRCUEclay.fig.exp

#### RRT Cgrowth ~ clay ####
fit.1 <- lm(Cgrowth_ngC_gDW_h_mean ~ clay_perc, data = compsoilmeanRRT)
fit.2 <- lm(Cgrowth_ngC_gDW_h_mean ~ log(clay_perc), data = compsoilmeanRRT)
fit.3 <- lm(Cgrowth_ngC_gDW_h_mean ~ exp(clay_perc), data = compsoilmeanRRT)
fit.4 <- lm(Cgrowth_ngC_gDW_h_mean ~ poly(clay_perc, 2, raw = TRUE), data = compsoilmeanRRT)
anova(fit.1, fit.2, fit.3, fit.4)
AIC(fit.1, fit.2, fit.3, fit.4) 
model.list <- list(fit.1, fit.2, fit.3, fit.4)
lapply(model.list, FUN = summary)


RRCgrowthclay.fig <- ggplot(data=compsoilmeanRRT, aes(x=clay_perc, y=Cgrowth_ngC_gDW_h_mean, fill = RR)) + 
  geom_errorbar(aes(ymin=Cgrowth_ngC_gDW_h_mean-Cgrowth_ngC_gDW_h_sd, 
                    ymax=Cgrowth_ngC_gDW_h_mean+Cgrowth_ngC_gDW_h_sd),
                width=0, color = "grey35") +
  stat_smooth(method="lm", formula = my.formula.exp, se=TRUE, linetype = "dotted", color = "grey35") +
  geom_point(aes(shape=landuse, fill=site), size = 2.5, color = "grey35") +
  stat_poly_eq(formula = my.formula.exp, 
               aes(label =  paste(
                 stat(adj.rr.label), 
                 stat(p.value.label), 
                 sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  theme_bw(base_size = 13) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(clay~"["~"%"~"]")) +
  ylab(expression(RR[Cgrowth])) +
  scale_fill_manual(values = c("#8F7700FF","#999999","#EFC000FF","#EFC000FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  facet_wrap(.~RR, scales="free") +
  scale_y_continuous(labels = formatter(nsmall = 1)) 
RRCgrowthclay.fig


#### RRT cumulresp ~ clay content ####
fit.1 <- lm(cumul_resp_incub_ngC_g_total_time_mean ~ clay_perc, data = compsoilmeanRRT)
fit.2 <- lm(cumul_resp_incub_ngC_g_total_time_mean ~ log(clay_perc), data = compsoilmeanRRT)
fit.3 <- lm(cumul_resp_incub_ngC_g_total_time_mean ~ exp(clay_perc), data = compsoilmeanRRT)
fit.4 <- lm(cumul_resp_incub_ngC_g_total_time_mean ~ poly(clay_perc, 2, raw = TRUE), data = compsoilmeanRRT)
anova(fit.1, fit.2, fit.3, fit.4)
AIC(fit.1, fit.2, fit.3, fit.4) 
model.list <- list(fit.1, fit.2, fit.3, fit.4)
lapply(model.list, FUN = summary)

RRcumulclay.fig <- ggplot(data=compsoilmeanRRT, aes(x=clay_perc, y=cumul_resp_incub_ngC_g_total_time_mean, fill = RR)) + 
  geom_errorbar(aes(ymin=cumul_resp_incub_ngC_g_total_time_mean-cumul_resp_incub_ngC_g_total_time_sd, 
                    ymax=cumul_resp_incub_ngC_g_total_time_mean+cumul_resp_incub_ngC_g_total_time_sd),
                width=0, color = "grey35") +
  stat_smooth(method="lm", formula = my.formula.log, se=TRUE, linetype = "dotted", color = "grey35") +
  geom_point(aes(shape=landuse, fill=site), size = 2.5, color = "grey35") +
  stat_poly_eq(formula = my.formula.log, 
               aes(label =  paste(
                 stat(adj.rr.label), 
                 stat(p.value.label), 
                 sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  theme_bw(base_size = 13) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  ylim(1,2.5) +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(clay~"["~"%"~"]")) +
  ylab(expression(RR[cumulative~respiration])) +
  scale_fill_manual(values = c("#8F7700FF","#999999","#EFC000FF","#EFC000FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  facet_wrap(.~RR, scales ="free")
RRcumulclay.fig


#### Figure RRT ~ clay ####
empty.plot <- ggplot()
legend <- get_legend(RRCUEclay.fig + theme(legend.position = "bottom"))# + guides(fill = "none"))
plot_clay <- plot_grid(RRCUEclay.fig +
                         theme(legend.position = "none",
                               axis.title.y = element_text(size=13),
                               axis.text.x = element_blank(),
                               axis.title.x = element_blank(),
                               axis.ticks.x = element_blank()),
                       empty.plot,
                       RRCgrowthclay.fig +
                         theme(legend.position = "none",
                               axis.title.y = element_text(size=13),
                               axis.text.x = element_blank(),
                               axis.title.x = element_blank(),
                               axis.ticks.x = element_blank()),
                       empty.plot,
                       RRcumulclay.fig +
                         theme(legend.position = "none",
                               axis.title.y = element_text(size=13),
                               axis.text.x = element_text(size=13)),
                       labels = c("a","","b","","c"), 
                       vjust = 1,
                       ncol = 1,
                       nrow = 5,
                       align = "v", 
                       axis = "rl",
                       rel_heights = c(1,0.05,1,0.05,1.2))
# 400x1000 (for 3 facets without legend)

plot_grid(plot_clay,
          legend,
          nrow = 2,
          rel_heights = c(1, 0.1))
