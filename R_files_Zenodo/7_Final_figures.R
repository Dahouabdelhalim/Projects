######
###### Graph final figures
######

library(cowplot)

### Figure 1
load("Graphs/Final/Figure 1A.RData") # map_plot
load("Graphs/Final/Figure 1B.RData") # group_plot
load("Graphs/Final/Figure 1C.RData") # whit_plot
load("Graphs/Final/Figure 1D.RData") # rgdist_plot
rw_onsite_metadata <- read_csv("~/Data/FLUXNET2015/On-site RW data/Chronologies/rw_onsite_metadata.csv") 
rw_onsite_sites <- unique(rw_onsite_metadata$Site)
toprow <- cowplot::plot_grid(map_plot, 
                             group_plot+ylab("Paired GPPâ€“RW obs.\\n(site x yr)"), 
                             labels = c('A', "B"), rel_widths = c(3,1))
bottomrow <- cowplot::plot_grid(whit_plot, 
                                r_gdist_plot+ylab(expression(GPP-RW[region]~corr.~(r[region]))), 
                                labels = c("C", "D")) # Load beginning of Regional RW-FLUX correlations.R
plot_grid(toprow, bottomrow, ncol = 1, rel_heights = c(1,1.2)) # load Climatic space.R for whit_plot
ggsave("Graphs/Final/Figure 1.png", height = 6, width = 8)

### Figure 2
load("Graphs/Final/Figure 2A.RData") # r_movwin_plot
load("Graphs/Final/Figure 2B.RData") # r_biome_plot
load("Graphs/Final/Figure 2C.RData") # rpl_plot
load("Graphs/Final/Figure 2D.RData") # resid_plot
fig2_top <- cowplot::plot_grid(r_movwin_plot+theme(axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank()), 
                               r_biome_plot+scale_x_discrete(position = "top", name = "Biome")+theme(axis.text.x = element_text(angle = -30, hjust = 0.9)), 
                               rel_widths = c(3,1.2), align = "h", axis = "b", labels = c("A", "B"))
fig2_bottom <- cowplot::plot_grid(rpl_plot+theme(axis.text.x.top = element_blank(), strip.text = element_blank()), 
                                  resid_plot+scale_x_discrete(name = "Effect"), 
                                  rel_widths = c(3,1.2), align = "h", axis = "t", labels = c("C", "D"))
plot_grid(fig2_top, fig2_bottom, ncol = 1)
ggsave("Graphs/Final/Figure 2.png", width = 8, height = 6)


### Figure 3
gpp_pr_m <- readRDS("Graphs/Final/pr GPP-climate.RDS") %>% mutate(Dep = "GPP", Month = as.factor(Month))
rw_pr_m <- readRDS("Graphs/Final/pr RWI-climate.RDS") %>% mutate(Dep = "RWI", Month = as.factor(Month))
pr_m <- bind_rows(gpp_pr_m, rw_pr_m)
season <- c("Winter", "Spring", "Summer", "Fall")

pr_m %>% 
  ungroup %>% 
  mutate(Dep = replace(Dep, Dep == "RWI", "RW")) %>% 
  mutate(Lag = factor(Lag, levels = c("Previous Year", "Current Year")))%>%
  filter(name %in% c("(Intercept)")) %>%
  mutate(Month = factor(Month, labels = season)) %>%
  mutate(Var = factor(Var, labels = c("T[mean]", "PDSI", "Rad"))) %>%
  mutate(estimate = logitr_inv(estimate), se = logitr_inv(se)) %>% 
  rowwise %>% 
  mutate(signif = signif (pval)) %>% 
  ggplot(aes(x = Month, y = estimate, fill = Month, group = name))+
  geom_col(show.legend = F, position = 'dodge', aes(alpha = Lag))+
  geom_errorbar(aes(ymin = estimate-se, ymax = estimate+se, width = 0.5))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_label(aes(label = signif, y = pmax(0,estimate+se)), show.legend = F, label.size = 0, fill = "transparent", size = 5, position = position_dodge(width=0.9))+
  facet_grid(Var~Dep, labeller = label_parsed)+
  coord_cartesian(ylim = c(-0.12,0.5), expand = T)+
  scale_fill_brewer(type = "div", palette = 1, direction = -1)+
  scale_alpha_manual(values = c(1,1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9), panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype =2))+
  ylab("Partial correlation")+xlab(NULL)
ggsave("Graphs/Final/Figure 3.png", width = 4, height = 5)
 
