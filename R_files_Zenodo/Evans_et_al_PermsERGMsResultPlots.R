require(dplyr)
require(ggplot2)
require(grid)
require(gridExtra)


allresults = read.csv("Evans_et_al_summarised_data.csv")

summary(allresults)

#type = whether network analysis method was an exponential random graph model "ERGM", or used controlled permutations "Perm"
#sex.degree = parameter controlling whether males had fewer (negative) or more (positive) or the same (zero) interaction strengths as females
#sex.assort = parameter controlling whether same-sex interactions are stronger (greater than 1), weaker (less than 1), or the same (equals 1) as opposite sex interactions
#i.dens = parameter controlling the strength of within-group interactions
#o.dens = parameter controlling the strength of between-group interactions
#d.eff = parameter controlling how greatly distance reduces between-group interactions
#obs.eff = parameter controlling the observer effort, with 1 being perfect and lower numbers indicating a greater number of social interactions missed
#dyad.fail.rate = the failure rate of models fitted to dyad based interactions. Rates are out of 100. If the effect was absent, higher values indicate more models returned a significant result. If the effect was present, higher values indicate more models returned non-significant results
#group.fail.rate = the failure rate of models fitted to grouping-event based interactions. Rates are out of 100. If the effect was absent, higher values indicate more models returned a significant result. If the effect was present, higher values indicate more models returned non-significant results

#False positives ####

false_pos = allresults %>% filter(sex.degree == 0, obs.eff==1)

#Dyad-based networks

#x = assort effect neg, 0, pos + type
#panels = effect and network
FP_sexeff_dyad = ggplot(false_pos, aes(x=interaction(sex.assort,type), 
                                         y=dyad.fail.rate, 
                                         fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  ylab("Failure %") + 
  theme_classic() + 
  annotation_custom(grob = textGrob(label="False positives in dyad-based networks", 
                                    just="left", gp=gpar(fontsize=12)), 
                    xmin = -0.5, xmax = -0.5,
                    ymin = 105, ymax = 105) +
  annotation_custom(grob = textGrob(label="a) i", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80),  clip = 'off') +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(2,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("-", "0\\nERGMs", "+",
                            "-", "0\\nPermutations", "+"))

FP_idens_dyad = ggplot(false_pos, aes(x=interaction(i.dens,type), 
                                        y=dyad.fail.rate, 
                                        fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  ylab("Failure %") +  
  theme_classic() + 
  annotation_custom(grob = textGrob(label="a) ii", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("0.4", "0.8\\nERGMs","1.2",
                            "0.4", "0.8\\nPermutations","1.2"))

FP_odens_dyad = ggplot(false_pos, aes(x=interaction(o.dens,type), 
                                        y=dyad.fail.rate, 
                                        fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  ylab("Failure %") + 
  theme_classic() + 
  annotation_custom(grob = textGrob(label="a) iii", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("0.1", "0.2\\nERGMs","0.4",
                            "0.1", "0.2\\nPermutations","0.4"))

FP_deff_dyad = ggplot(false_pos, aes(x=interaction(d.eff,type), 
                                       y=dyad.fail.rate, 
                                       fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  ylab("Failure %") + 
  theme_classic() + 
  annotation_custom(grob = textGrob(label="a) iv", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("0", "4\\nERGMs","8",
                            "0", "4\\nPermutations","8"))

#Grouping-event based networks

FP_sexeff_grp = ggplot(false_pos, aes(x=interaction(sex.assort,type), 
                                        y=group.fail.rate, 
                                        fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic()+ labs(fill= "Method") +
  annotation_custom(grob = textGrob(label="False positives in grouping-based networks", 
                                    just="left", gp=gpar(fontsize=12)), 
                    xmin = -0.5, xmax = -0.5,
                    ymin = 105, ymax = 105) +
  annotation_custom(grob = textGrob(label="b) i", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(2,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations"))+
  scale_x_discrete(labels=c("-", "0\\nERGMs", "+",
                            "-", "0\\nPermutations", "+"))

FP_idens_grp = ggplot(false_pos, aes(x=interaction(i.dens,type), 
                                       y=group.fail.rate, 
                                       fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic()+ labs(fill= "Method") +
  annotation_custom(grob = textGrob(label="b) ii", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations"))+
  scale_x_discrete(labels=c("0.4", "0.8\\nERGMs","1.2",
                            "0.4", "0.8\\nPermutations","1.2"))

FP_odens_grp = ggplot(false_pos, aes(x=interaction(o.dens,type), 
                                       y=group.fail.rate, 
                                       fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic()+ labs(fill= "Method") +
  annotation_custom(grob = textGrob(label="b) iii", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations"))+
  scale_x_discrete(labels=c("0.1", "0.2\\nERGMs","0.4",
                            "0.1", "0.2\\nPermutations","0.4"))

FP_deff_grp = ggplot(false_pos, aes(x=interaction(d.eff,type), 
                                      y=group.fail.rate, 
                                      fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic()+ labs(fill= "Method") +
  annotation_custom(grob = textGrob(label="b) iv", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations"))+
  scale_x_discrete(labels=c("0", "4\\nERGMs","8",
                            "0", "4\\nPermutations","8"))

#False negatives ####

false_neg = allresults %>% filter(sex.degree != 0, obs.eff == 1) 

FN_sexeff_dyad = ggplot(false_neg, aes(x=interaction(sex.assort,type), 
                                         y=dyad.fail.rate, 
                                         fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic() + labs(fill= "Method") +
  annotation_custom(grob = textGrob(label="False negatives in dyad-based networks", 
                                    just="left", gp=gpar(fontsize=12)), 
                    xmin = -0.5, xmax = -0.5,
                    ymin = 105, ymax = 105) +
  annotation_custom(grob = textGrob(label="c) i", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(2,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("-", "0\\nERGMs", "+",
                            "-", "0\\nPermutations", "+"))


FN_idens_dyad = ggplot(false_neg, aes(x=interaction(i.dens,type), 
                                        y=dyad.fail.rate, 
                                        fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic() + labs(fill= "Method") +
  annotation_custom(grob = textGrob(label="c) ii", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("0.4", "0.8\\nERGMs","1.2",
                            "0.4", "0.8\\nPermutations","1.2"))


FN_odens_dyad = ggplot(false_neg, aes(x=interaction(o.dens,type), 
                                        y=dyad.fail.rate, 
                                        fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic() + labs(fill= "Method") +
  annotation_custom(grob = textGrob(label="c) iii", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("0.1", "0.2\\nERGMs","0.4",
                            "0.1", "0.2\\nPermutations","0.4"))


FN_deff_dyad = ggplot(false_neg, aes(x=interaction(d.eff,type), 
                                       y=dyad.fail.rate, 
                                       fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic() + labs(fill= "Method") +
  annotation_custom(grob = textGrob(label="c) iv", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("0", "4\\nERGMs","8",
                            "0", "4\\nPermutations","8"))



FN_sexeff_grp = ggplot(false_neg, aes(x=interaction(sex.assort,type), 
                                        y=group.fail.rate, 
                                        fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  labs(tag = "Assortativity\\neffect") +
  theme_classic()+ 
  annotation_custom(grob = textGrob(label="False negatives in grouping-based networks", 
                                    just="left", gp=gpar(fontsize=12)), 
                    xmin = -0.5, xmax = -0.5,
                    ymin = 105, ymax = 105) +
  annotation_custom(grob = textGrob(label="d) i", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(2,1,0.5,1), "lines"),
        legend.position = "none",
        plot.tag.position = "right",
        plot.tag = element_text(vjust = 10, size = 12)) +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations"))+
  scale_x_discrete(labels=c("-", "0\\nERGMs", "+",
                            "-", "0\\nPermutations", "+"))


FN_idens_grp = ggplot(false_neg, aes(x=interaction(i.dens,type), 
                                       y=group.fail.rate, 
                                       fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  labs(tag = "Within\\ngroup\\ninteraction\\ndensity") +
  theme_classic()+ 
  annotation_custom(grob = textGrob(label="d) ii", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none",
        plot.tag.position = "right",
        plot.tag = element_text(vjust = 10, size = 12)) +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations"))+
  scale_x_discrete(labels=c("0.4", "0.8\\nERGMs","1.2",
                            "0.4", "0.8\\nPermutations","1.2"))

FN_odens_grp = ggplot(false_neg, aes(x=interaction(o.dens,type), 
                                       y=group.fail.rate, 
                                       fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  labs(tag = "Between\\ngroup\\ninteraction\\ndensity") +
  theme_classic()+ 
  annotation_custom(grob = textGrob(label="d) iii", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none",
        plot.tag.position = "right",
        plot.tag = element_text(vjust = 10, size = 12)) +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations"))+
  scale_x_discrete(labels=c("0.1", "0.2\\nERGMs","0.4",
                            "0.1", "0.2\\nPermutations","0.4"))


FN_deff_grp = ggplot(false_neg, aes(x=interaction(d.eff,type), 
                                      y=group.fail.rate, 
                                      fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  labs(tag = "Distance\\neffect") +  
  theme_classic()+ 
  annotation_custom(grob = textGrob(label="d) iv", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  coord_cartesian(xlim = c(1,6), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(1,1,0.5,1), "lines"),
        legend.position = "none",
        plot.tag.position = "right",
        plot.tag = element_text(vjust = 10, size = 12)) +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations"))+
  scale_x_discrete(labels=c("0", "4\\nERGMs","8",
                            "0", "4\\nPermutations","8"))



#pdf("NetSims1 Fig 4x4 plots.pdf")
grid.arrange(FP_sexeff_dyad, FP_sexeff_grp, FN_sexeff_dyad, FN_sexeff_grp,
             FP_idens_dyad, FP_idens_grp, FN_idens_dyad, FN_idens_grp,
             FP_odens_dyad, FP_odens_grp, FN_odens_dyad, FN_odens_grp,
             FP_deff_dyad, FP_deff_grp, FN_deff_dyad, FN_deff_grp, ncol=4,
             heights = c(1.2,1,1,1),
             widths = c(1.1,1,1,1.3))
#dev.off()
#export as 37.3 x 19.2 cm

#Observation effort ####

false.pos.ob = allresults %>% filter(sex.degree == 0)
false.neg.ob = allresults %>% filter(sex.degree != 0)


FP_obseff_dyad = ggplot(false.pos.ob, aes(x=interaction(obs.eff,type), 
                                           y=dyad.fail.rate, 
                                           fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic() + 
  ylab("Failure %") +
  annotation_custom(grob = textGrob(label="False positives in dyad-based networks", 
                                    just="left", gp=gpar(fontsize=12)), 
                    xmin = -0.5, xmax = -0.5,
                    ymin = 105, ymax = 105) +
  annotation_custom(grob = textGrob(label="a)", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  annotation_custom(grob = textGrob(label="ERGMs", just="left", gp=gpar(fontsize=10)), 
                    xmin = 2, xmax = 2,
                    ymin = -30, ymax = -30) +
  annotation_custom(grob = textGrob(label="Permutations", just="left", gp=gpar(fontsize=10)), 
                    xmin = 5.5, xmax = 5.5,
                    ymin = -30, ymax = -30) +
  coord_cartesian(xlim = c(1,8), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(2,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("0.3", "0.6\\n","0.9","1.0",
                            "0.3", "0.6","0.9","1.0"))


FN_obseff_dyad = ggplot(false.neg.ob, aes(x=interaction(obs.eff,type), 
                                           y=dyad.fail.rate, 
                                           fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic() + labs(fill= "Method") +
  annotation_custom(grob = textGrob(label="False negatives in dyad-based networks", 
                                    just="left", gp=gpar(fontsize=12)), 
                    xmin = -0.5, xmax = -0.5,
                    ymin = 105, ymax = 105) +
  annotation_custom(grob = textGrob(label="b)", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  annotation_custom(grob = textGrob(label="ERGMs", just="left", gp=gpar(fontsize=10)), 
                    xmin = 2, xmax = 2,
                    ymin = -30, ymax = -30) +
  annotation_custom(grob = textGrob(label="Permutations", just="left", gp=gpar(fontsize=10)), 
                    xmin = 5.5, xmax = 5.5,
                    ymin = -30, ymax = -30) +
  coord_cartesian(xlim = c(1,8), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(2,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("0.3", "0.6\\n","0.9","1.0",
                            "0.3", "0.6","0.9","1.0"))

FP_obseff_grp = ggplot(false.pos.ob, aes(x=interaction(obs.eff,type), 
                                          y=group.fail.rate, 
                                          fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic() + 
  ylab("Failure %") +
  annotation_custom(grob = textGrob(label="False positives in grouping-based networks", 
                                    just="left", gp=gpar(fontsize=12)), 
                    xmin = -0.5, xmax = -0.5,
                    ymin = 105, ymax = 105) +
  annotation_custom(grob = textGrob(label="c)", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  annotation_custom(grob = textGrob(label="ERGMs", just="left", gp=gpar(fontsize=10)), 
                    xmin = 2, xmax = 2,
                    ymin = -30, ymax = -30) +
  annotation_custom(grob = textGrob(label="Permutations", just="left", gp=gpar(fontsize=10)), 
                    xmin = 5.5, xmax = 5.5,
                    ymin = -30, ymax = -30) +
  coord_cartesian(xlim = c(1,8), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(2,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("0.3", "0.6\\n","0.9","1.0",
                            "0.3", "0.6","0.9","1.0"))


FN_obseff_grp = ggplot(false.neg.ob, aes(x=interaction(obs.eff,type), 
                                          y=group.fail.rate, 
                                          fill=type)) + 
  geom_violin(draw_quantiles=c(0.25,0.75), colour="white") +
  geom_violin(draw_quantiles=c(0.5),  alpha=0.1, colour="black") +
  theme_classic() + labs(fill= "Method") +
  annotation_custom(grob = textGrob(label="False negatives in grouping-based networks", 
                                    just="left", gp=gpar(fontsize=12)), 
                    xmin = -0.5, xmax = -0.5,
                    ymin = 105, ymax = 105) +
  annotation_custom(grob = textGrob(label="d)", just="left", gp=gpar(fontsize=12)), 
                    xmin = 0.5, xmax = 0.5,
                    ymin = 85, ymax = 85) +
  annotation_custom(grob = textGrob(label="ERGMs", just="left", gp=gpar(fontsize=10)), 
                    xmin = 2, xmax = 2,
                    ymin = -30, ymax = -30) +
  annotation_custom(grob = textGrob(label="Permutations", just="left", gp=gpar(fontsize=10)), 
                    xmin = 5.5, xmax = 5.5,
                    ymin = -30, ymax = -30) +
  coord_cartesian(xlim = c(1,8), ylim=c(0,80), clip = 'off') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        plot.margin = unit(c(2,1,0.5,1), "lines"),
        legend.position = "none") +
  scale_fill_manual(values=c("cornflowerblue","darkorange"),
                    labels = c("ERGMs", "Permutations")) +
  scale_x_discrete(labels=c("0.3", "0.6\\n","0.9","1.0",
                            "0.3", "0.6","0.9","1.0"))

grid.arrange(FP_obseff_dyad, FN_obseff_dyad, FP_obseff_grp, FN_obseff_grp, ncol=2,
             widths = c(1.2,1))
#export as 19.6 x 10.9 cm