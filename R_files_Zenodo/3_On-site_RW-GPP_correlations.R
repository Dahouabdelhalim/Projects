#####
##### Calculate correlations between GPP and on-site RW data
#####

### Load data
# Metadata
flux_metadata <- read_csv("Data/flux_metadata.csv") 
rw_onsite_metadata <- read_csv("Data/rw_onsite_metadata.csv") 

# Observations
rw_onsite <- read_csv("Data/rw_onsite.csv")
flux_onsite <- read_csv("Data/flux_onsite.csv") 

# Merge RW and flux data
flux_rw_onsite <- rw_onsite %>% 
  inner_join(flux_onsite %>% 
               mutate(Year = year(Date), Month = month(Date)) %>% 
               group_by(Site, Year) %>% 
               filter(mean(QC)>0) %>% # Filter out low quality data
               filter(n()>=11)) %>% # Filter out site years with unsufficient data
  filter(!is.na(GPP))

# Preview of RW and GPP time series (undetrended)
flux_rw_onsite %>%
  group_by(Site, Year, IDcrn) %>%
  summarise(across(c(RWI, GPP), mean)) %>% # Yearly average
  group_by(Site, IDcrn) %>%
  mutate(across(c(RWI, GPP), std)) %>% # Standardization
  left_join(rw_onsite_metadata) %>% # MergÃ©
  group_by(Site, Year) %>% 
  summarise(across(c(RWI, GPP), ~weighted.mean(.,RW_proportion))) %>% # Average across species weighting by species composition 
  group_by(Site) %>%
  # filter(n()>=4) %>%
  mutate(r = cor.test(RWI, GPP)$estimate,
         pval = cor.test(RWI, GPP)$p.value %>% map(signif) %>% str_replace("ns", ""))  %>%
  pivot_longer(c(GPP, RWI), names_to = "Variable", values_to = "Z-score") %>% 
  ggplot(aes(x = Year, y = `Z-score`, color = Variable))+
  geom_line()+
  geom_label(data = . %>% group_by(Site, pval, r) %>% summarise, aes(label = paste0(round(r, 2), pval)), x = -Inf, y = Inf, color = "black", label.size = 0, fill = NA, vjust = 1, hjust = 0)+
  scale_color_manual(values = c("red", "black"))+
  scale_x_continuous(breaks = c(1995, 2005, 2015))+
  facet_wrap(~Site)+
  theme(legend.position = c(1,0), legend.justification = c(1,0), legend.direction = "horizontal")
ggsave("Graphs/Overlay GPP&RWI time series - by site.png", width = 8, height = 6)
             
### Analyses
## Preliminary analysis of the role of pluri-annual vs year-to-year GPP variations in explaining RW
trd_m <- flux_rw_onsite %>% 
  left_join(rw_onsite_metadata %>% select(Site, IDcrn, RW_proportion, RW_species)) %>%
  group_by(Site, Year) %>%
  summarise(across(c(RWI, GPP), ~weighted.mean(.,RW_proportion))) %>%
  group_by(Site) %>%
  filter(n()>=10) %>%
  mutate(GPP = std(GPP), RWI = std(RWI), 
         GPP_dtrd = dtrd(GPP), GPP_trd = trd(GPP)) %>%
  lme4::lmer(RWI~GPP_dtrd+GPP_trd-1+(GPP_dtrd+GPP_trd-1|Site), data = .)
coefs.lmer(trd_m)
MuMIn::r.squaredGLMM(trd_m)
anova(trd_m)

## Calculate correlations for each site and each time window
r <- rw_onsite %>%
  left_join(rw_onsite_metadata) %>%
  group_by(Site, Year) %>% 
  summarise(across(c(RWI, RWI_dtrd), ~weighted.mean(.,RW_proportion))) %>% 
  nest(data = -Site) %>% 
  mutate(r = map2(data,Site, function(x,y) {
    cat(".")
    res <- data.frame(Lag = -12:12, Period = rep(1:12, each = 25)) %>%
      filter(!Lag - Period < -12) %>%
      mutate(r = map2(Lag, Period, function(l, p) {
        x %>% 
          inner_join(
            flux_onsite %>%
              group_by(Site, Year) %>%
              filter(!is.na(GPP)) %>% 
              filter(Site %in% y) %>% 
              mutate(Date = Date + duration(l, "month") + duration(1, "day")) %>% 
              mutate(Year = year(Date), Month = month(Date)) %>% 
              group_by(Site, Year) %>%
              filter(Month %in% 1:p) %>% 
              filter(mean(QC)>0.5) %>% # More stringent filtering of high quality data in order to ensure less noise
              filter(n() >= p-1) %>% 
              summarise(across(c(GPP, GPP_dtrd), mean), .groups = "drop"),
            by = "Year"
          ) %>% 
          summarise(r = cor(RWI_dtrd, GPP_dtrd), n = n())
      }))
  })) %>% 
  select(-data) %>% 
  unnest(r) %>% unnest(r)

## Temporal plot
r_movwin_plot <- r %>% 
  mutate(Weight = n-2) %>% 
  group_by(Period, Lag) %>%
  summarise(pval = try(weights::wtd.t.test(r, 0, Weight)$coefficients["p.value"]),
            r = weighted.mean(r, Weight, na.rm = T)) %>%
  ungroup %>% 
  mutate(Year = c(Lag < 1) %>% factor(., labels = c("Previous year", "Current year"))) %>% 
  ggplot(aes(x = -Lag, y = Period, fill = r))+
  geom_tile()+
  geom_label(data = . %>% filter(pval<.05, pval>.01), aes(label = round(r, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pval<.01), aes(label = round(r, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  coord_cartesian(expand = F)+
  scale_fill_distiller(type ="div", palette = 4, limits = c(-.4,.4), oob = scales::squish, name = expression(r[on-site]))+
  scale_y_continuous(breaks = (1:6)*2, name = "Window length (month)", sec.axis = dup_axis(name = NULL))+
  scale_x_continuous(breaks = (-6:11)*2, labels = rep(c("Ja", "Mr", "My", "Jl", "Se", "No"), 3), name = 'Window onset (MOY)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", 
        panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = c(0.99,0.99), legend.justification = c(1,1))+
  facet_grid(.~Year, scales = "free_x")
r_movwin_plot
save(r_movwin_plot, file = "Graphs/Final/Figure 2A.RData")

## Spatial variations
r_biome <- r %>% 
  filter(Period == 10, Lag == 2) %>% # Consider the case with maximum average correlation
  group_by(Site) %>% 
  summarise(r = weighted.mean(r, n-2),
            n = weighted.mean(n, n-2)) %>% 
  filter(n>=4) %>% 
  ungroup %>% 
  left_join(flux_metadata %>%
              mutate(biome = KoppenGeiger_find(Lon, Lat) %>%
                       str_extract("[:alpha:]{2}") %>%
                       c("BS" = "Med./Sav.", "Cf" = "Temperate", "Cs" = "Med./Sav.", "Df" = "Bor./Cont.")[.]
                     )) %>%
  mutate(biome = factor(biome, levels = c("Bor./Cont.", "Temperate", "Med./Sav.")))
  
# Between site variations
Hmisc::wtd.quantile(r_biome$r, r_biome$n, probs = c(0.1, 0.9))

# Differences between biomes
r_biome_m <- r_biome %>% 
  lm(logitr(r)~biome, weights = n-2, data = .)
anova(r_biome_m) # No biome effect
emmeans::emmeans(r_biome_m, pairwise~biome) # No pairwise differences between 
emmeans::emmeans(r_biome_m, ~biome) # No pairwise differences between 

r_biome_plot <- r_biome %>% 
  ggplot(aes(x = biome, y = r, size = n, fill = biome))+
  geom_boxplot(show.legend = F, outlier.color = NA)+
  geom_point(position = position_dodge(width = 0.8), show.legend = F, shape = 21, color = "black", fill = "grey", alpha = 0.5)+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_fill_brewer(type = "qual", palette = 3, aesthetics = c("color", "fill"))+
  scale_alpha_manual(values =  c(0.4,1))+
  ylab(expression(r[on-site]))+
  scale_size(range = c(0.5,3))+
  theme(axis.text.x = element_text(angle = 30, hjust = 0.9))
r_biome_plot
save(r_biome_plot, file = "Graphs/Final/Figure 2B.RData")
