#####
##### ITRDB and FLUXNET correlations
#####

#### RW and flux data overlap 
flux_rw_overlap <- read_csv("Data/flux_rw_overlap.csv")

flux_years <- flux_rw_overlap %>% 
  group_by(Site = RW_site) %>% 
  summarise(Flux_start = min(Start), 
            Flux_end = max(End))

rw_years <- flux_rw_overlap %>% 
  group_by(Site) %>% 
  summarise(RW_start = min(RW_start), 
            RW_end = max(RW_end))

#### RW data
rw <- read_csv("Data/rw.csv") %>% 
  left_join(flux_years) %>% 
  group_by(Site) %>% 
  filter(Year>=Flux_start, Year<=Flux_end+1) %>% 
  select(-contains("Flux")) %>% 
  mutate(RWI_dtrd = dtrd(RWI))

rw_metadata <- read_csv("Data/rw_metadata.csv")
rw_onsite_metadata <- read_csv("Data/rw_onsite_metadata.csv") 
rw_onsite_sites <- unique(rw_onsite_metadata$Site)

#### Fluxnet data 
flux_metadata <- read_csv("Data/flux_metadata.csv") %>% 
  mutate(Group = if_else(str_detect(IGBP, "N"), "Gymnosperm", "Angiosperm"), 
         Group = if_else(str_detect(IGBP, "M"), "Mixed", Group)) 

flux <- read_csv("Data/flux.csv") %>% 
  left_join(flux_metadata %>% select(Site, Lat)) %>% 
  mutate(Date = if_else(Lat<0, Date-duration(181, "days"), Date), # Tree rings in the SH are dated such that the reported year corresponds to the year of when growth started
         Date = round_date(Date, "month"),                        # Thus GPP needs to be back-shifted by six months, such that July of year X becomes Jan of year X
         Year = year(Date), Month = month(Date)) %>%              # 
  select(-Lat) %>% 
  left_join(rw_years) %>% 
  group_by(Site) %>% 
  filter(Year>=RW_start, Year<=RW_end+1) %>% 
  select(-contains("RW")) %>% 
  mutate(GPP_dtrd = dtrd(GPP))

flux_sp <- read_csv("Data/flux_species.csv")
flux_richness <- flux_sp %>% 
  filter(!is.na(Group)) %>% 
  group_by(Site) %>% 
  summarise(GymnoProp = sum(Proportion[Group == "Gymnosperm"]),
            Richness = n(),
            Shannon = -sum(Proportion*log(Proportion)))

#### Merge data
flux_rw <- flux_rw_overlap %>%
  left_join(rw %>% rename(RW_site = Site)) %>% 
  left_join(flux %>%
              group_by(Site, Year) %>%
              filter(mean(QC)>0.5) %>% # Filter out bad quality data
              mutate(Date = floor_date(Date, "12 months")) %>%
              mutate(Year = year(Date), Month = month(Date)) %>% 
              group_by(Site, Month, Year) %>%
              summarise(across(c(NEE, GPP, R, GPP_dtrd), mean, na.rm = T))
            ) %>% 
  group_by(Site, RW_site) %>% 
  filter(!is.na(GPP), !is.na(RWI), !is.na(RW_climdist)) 

######
###### Calculate correlation between RWI and annual GPP with no lag

## Dataframe of correlation per site
r <- flux_rw %>%
  group_by(Month, Site, RW_site) %>%
  mutate(n = n()) %>% 
  filter(n >= 3) %>%
  nest %>%
  mutate(r = map(data, function(x) cor(x$RWI_dtrd, x$GPP_dtrd, method = "pearson")) %>% unlist,
         n = map(data, function(x) unique(x$n))) %>%
  select(-data) %>% 
  unnest(c(r,n)) %>% 
  left_join(flux_rw_overlap)
  
## Fit correlation extinction model and make model prediction for controlled distance variations
r_lmer <- r %>%
  filter(abs(r)<0.99) %>%
  group_by(Month) %>%
  mutate(across(c(RW_geodist, RW_phylodist, RW_climdist), function(x) x/quantile(x, 0.95))) %>%
  nest %>%
  mutate(m = map(data, function(x) lme4::lmer(logitr(r)~RW_phylodist*RW_climdist*RW_geodist+(1|Site)+(1|Site:RW_species), data = x, weights = n-2)),
         RW_geodist = list(rep(seq(0,1,length.out = 100), each = 1)),
         RW_climdist = list(rep(seq(0,0.65,length.out = 100), each = 1)),  # clim distance is ~0.65 for geodist = 1000
         RW_phylodist = list(rep(seq(0,0.278,length.out = 100), each = 1)), # phylo distance is 0.278 on average for onsite series (1.39 in absolute value)
         boot = map(m, function(x) lme4::bootMer(x,function(x){
           predict(x, re.form = ~0,
         newdata = data.frame(RW_geodist = unlist(RW_geodist),
                              RW_climdist = unlist(RW_climdist),
                              RW_phylodist = unlist(RW_phylodist))) %>% logitr_inv
         }, nsim=500)),
         r = map(boot, function(x) x$t %>% data.frame %>% summarise_all(mean) %>% unlist),
         se = map(boot, function(x) x$t %>% data.frame %>% summarise_all(sd) %>% unlist)
        )

## Model stats
# AIC(r_lmer$m[[1]])
# coefs.lmer(r_lmer$m[[1]])
# MuMIn::r.squaredGLMM(r_lmer$m[[1]])
# plot(r_lmer$m[[1]], resid(.,scaled = T)~fitted(.))
# lattice::qqmath(r_lm$m[[1]])

## Plot results
r_gdist_plot <- r %>% 
  ungroup %>% 
  select(RW_geodist, r) %>% 
  ggplot(aes(x = RW_geodist, y = r))+
  geom_ribbon(data = r_lmer %>% 
                unnest(c(r, se, RW_geodist)) %>% 
                select(c(r, se, RW_geodist)) %>% 
                mutate(RW_geodist = RW_geodist*1000),  
              aes(ymin = r-se, ymax = r+se), alpha = 0.3, fill = "forestgreen")+
  geom_line(data = r_lmer %>%
              ungroup %>% 
              unnest(c(r, RW_geodist)) %>%
              select(c(r, RW_geodist)) %>% 
              mutate(RW_geodist = RW_geodist*1000), 
            aes(color = "2", size = "2"))+
  stat_summary_bin(bins = 20, fun.data = mean_se, size = 0.4)+
  geom_smooth(method = "gam", formula = y~s(x, k = 3), aes(color = "1", size = "1"), level = 0.70)+
  geom_hline(yintercept = 0, linetype = 2)+
  annotate("segment", x = 0, xend = 250, y = 0.24, yend = 0.26, arrow = arrow(ends = "first", type = "closed", length = unit(6, "points")))+
  annotate("text", x = 250, y = 0.26, label = "r[D==0]", parse = T, hjust = 0)+
  coord_cartesian(expand = F)+
  ylab(expression(r[region]))+xlab("Geographic distance (km)")+
  scale_color_manual(values = c("black", "forestgreen"), name = "Model", labels = c(expression(D[geo]), "Full"))+
  scale_size_manual(values = c(1, 2), name = "Model", labels = c(expression(D[geo]), "Full"))+
  theme(axis.text.x = element_text(hjust = 0.9), legend.position = c(0.99,0.99), legend.justification = c(1,1))
r_gdist_plot
save(r_gdist_plot, file = "Graphs/Final/Figure 1D.RData")

#####
##### Moving window correlations
## Function to extract predictions from model object
pred <- function(x) predict(x, newdata = data.frame(RW_phylodist = 0, RW_climdist = 0, RW_geodist = 0), re.form = ~0)

## Function to calculate correlation and fit model for a given time window
rmov <- function(lag = 0, period = 12){
  cat(".")
  # Merge flux data of period of interest
  flux_rw <- flux_rw_overlap %>% 
    filter(RW_overlap>=5) %>%
    left_join(rw %>% 
                rename(RW_site = Site),
              by = "RW_site"
    ) %>%
    left_join(flux %>% 
                group_by(Site, Year) %>%
                mutate(Date = Date + duration(lag, "months") + duration(1, "day")) %>% # Apply lag
                mutate(Year = year(Date), Month = month(Date)) %>% 
                group_by(Site, Year) %>% 
                filter(Month %in% c(1:period)) %>% # Filter period of interest
                filter(mean(QC)>0.5) %>% # Filter out bad quality data
                filter(n() >= period-1) %>% # Filter out site-years with not enough obseravations
                summarise(across(c(GPP, GPP_dtrd), mean), .groups = "drop"),
              by = c("Site", "Year")
    ) %>% 
    filter(!is.na(GPP) & !is.na(RWI)) %>% 
    filter(!is.na(RW_climdist))
  
  # Calculate correlation
  r <- flux_rw %>% 
    group_by(Site, RW_site, IGBP, RW_group, RW_geodist, RW_phylodist, RW_climdist, RW_species) %>%
    mutate(n = n()) %>% 
    filter(n >= 3) %>% 
    summarise(r = cor(RWI_dtrd, GPP_dtrd, method = "pearson"),
              n = unique(n),
              .groups = "drop")

  # Fit model and calculate interecpt (rd=0)
  r0 <- r %>% 
    filter(abs(r)<0.99) %>%
    mutate(across(c(RW_geodist, RW_phylodist, RW_climdist), function(x) x/quantile(x, 0.95))) %>%
    group_by(group = 1) %>% nest %>% 
    mutate(m = map(data, function(x) lme4::lmer(logitr(r)~RW_phylodist*RW_climdist*RW_geodist+(1|Site)+(1|Site:RW_species), data = x, weights = n-2)),
           r0 = map(m, function(x) lme4::bootMer(x,function(x) pred(x) %>% logitr_inv, nsim=25)),
           se = map(r0, function(x) x$t %>% sd),
           r0 = map(m, function(m) pred(m) %>% logitr_inv),
           pval = map2(r0, se, function(x,y) pnorm(0, mean = x, sd = y))
    ) %>%
    select(-data) %>% 
    unnest(c(group, r0, se, pval)) %>%
    mutate(signif = map(pval,signif) %>% unlist)
  
  return(r0)
}

# Run correlation function and collect result in data frame
rpl <- tibble(period = rep((1:12),37),lag = rep((-24:12),each = 12)) %>% 
  filter(!lag - period < -12) %>%
  mutate(r0 = map2(lag, period, ~rmov(.x,.y)))

# Plot of the calculated intercept rd=0
rpl_plot <- rpl %>%
  unnest(r0) %>% 
  group_by(period, lag, group, r0, se, pval, signif) %>% 
  summarise %>% ungroup %>% 
  mutate(Year = c(lag < 1) %>% factor(., labels = c("Previous year", "Current year"))) %>% 
  ggplot(aes(x = -lag+1, y = period, fill = r0))+
  geom_tile(size = 1)+
  geom_label(data = . %>% filter(pval<=0.1, pval>0.01), aes(label = round(r0, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent")+
  geom_label(data = . %>% filter(pval<=0.01), aes(label = round(r0, 2) %>% str_replace("0","")), label.size = 0, size=3, fill = "transparent", fontface = "bold")+
  scale_fill_distiller(type = "div", palette = 5, limits = c(-0.4,0.4), name=expression(r[D==0]), oob = scales::squish)+
  scale_color_manual(values = c("transparent", "black"))+
  facet_grid(~Year, scales= "free_x")+
  coord_cartesian(expand = F)+
  scale_y_continuous(breaks = (1:6)*2, name = "Window length (month)", sec.axis = dup_axis(name = NULL))+
  scale_x_continuous(breaks = (-6:11)*2+1, labels = rep(c("Ja", "Mr", "My", "Jl", "Se", "No"), 3), name = 'Window onset (MOY)', sec.axis = dup_axis(name = NULL))+
  theme(panel.spacing.x = unit(0,"lines"), strip.placement = "outside", panel.grid = element_line(linetype = 3, color = "lightgrey"),
        legend.position = c(0.99,0.99), legend.justification = c(1,1))
rpl_plot
save(rpl_plot, file = "Graphs/Final/Figure 2C.RData")

#####
##### Model residuals as function of site characterisics

## Extract residuals
resid <- rpl %>% 
  # filter(period > 3) %>%
  unnest(r0) %>%
  mutate(frame = map(m, ~.@frame),
         residuals = map(m, ~residuals(.))) %>% 
  select(-m, -r0, -se, -pval, -signif, -group) %>% 
  unnest(cols = c(frame, residuals)) %>% 
  left_join(flux_rw_overlap %>% group_by(RW_species, RW_group) %>% summarise) %>%
  left_join(flux_metadata)

## Fit model
resid_m <- resid %>%
  filter(!is.na(LAI)) %>%
  left_join(flux_richness) %>% 
  group_by(Lag = lag < 1) %>% # Group by current vs previous year
  mutate(MACWD = sqrt(MACWD), MAP = sqrt(MAP)) %>% # Transform MACWD and MAP to respect normal distribution
  mutate(across(c(MAT, MAP, MACWD, MAR, MADF, LAI, Lat, Elev, GymnoProp, Richness, Shannon), std)) %>% # Standardize variables
  nest %>%
  mutate(m = map(data, function(x) lme4::lmer(`logitr(r)` ~ -1 + (MAT+MACWD+LAI+GymnoProp+Richness)+(MAT+MACWD+LAI+GymnoProp+Richness):(RW_geodist*RW_climdist*RW_phylodist)+(1|period:lag), data = x, weights = `(weights)`)),
         coefs = map(m, coefs.lmer))
## Calculate AIC
resid_m %>% mutate(stat = map(m, AIC) %>% unlist) %>% .$stat

## Plot result
resid_plot <- resid_m %>% 
  select(-data, -m) %>% 
  unnest(coefs) %>%
  filter(!str_detect(name,"dist")) %>% 
  mutate(name = str_replace(name, "RW_group", "")) %>% 
  mutate(name = str_replace(name, "sperm", ".")) %>%
  mutate(name = str_replace(name, "GymnoProp", "Gymno. prop.")) %>%
  mutate(name = str_replace(name, "Richness", "Sp. richness")) %>%
  mutate(name = factor(name, levels = c("Gymno. prop.", "Sp. richness", "LAI", "MACWD", "MAT", "GST"))) %>% 
  mutate(Lag = factor(Lag, levels = c(FALSE, TRUE), labels = c("Prv. yr", "Crt. yr"))) %>% 
  mutate(estimate = logitr_inv(estimate), se = logitr_inv(se)) %>%
  ggplot(aes(x = name, fill = name, color = name, y = estimate, alpha = Lag, group = Lag))+ 
  geom_col(show.legend = T, position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = estimate-se, ymax = estimate+se), width = 0.5, position = position_dodge(width = 0.9), color = "black", linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 2)+
  scale_fill_brewer(type = "qual", palette = 6, direction = -1, aesthetics = c("fill", "color"), guide = F)+
  scale_alpha_manual(values = c(0.5,1), name = NULL)+
  xlab(NULL)+ylab("Effect size")+
  theme(legend.position = c(0.25,1.05), legend.justification = c(0,1), legend.background = element_blank())+
  theme(axis.text.x = element_text(angle = 30, hjust = 0.9))
resid_plot
save(resid_plot, file = "Graphs/Final/Figure 2D.RData")

