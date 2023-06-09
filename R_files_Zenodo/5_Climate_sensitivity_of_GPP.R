######
###### Effect of diffuse light on GPP and tree ring width
######

library(ppcor) # Partial correlations
select <- dplyr::select
library(lme4) # Mixed models
library(MuMIn) 

#### Load data
flux_metadata <- read_csv("Data/flux_metadata.csv") 

flux <- read_csv("Data/flux.csv") %>% 
  group_by(Site, Year) %>% filter(sum(QC)>0) %>% 
  left_join(flux_metadata %>% select(Site, Lat)) %>%
  mutate(Date = if_else(Lat<0, Date - duration(181,"days"), Date), # 6 months lag for southern sites
         Date = round_date(Date, "month")) %>%
  select(-Lat) %>%
  mutate(Date = Date + duration(32, "days")) %>% # 1 month lag for season definition perpose
  mutate(Date = floor_date(Date, "3 months")) %>%           # Average across 3-months periods
  mutate(Year = year(Date), Month = month(Date)) %>%        # 
  group_by(Site, Year, Month) %>%                           #
  summarise(across(c("NEE", "GPP", "R"), mean, na.rm = T))  #

# Calculate the proportion of total GPP participated by each season for weighting of correlations later on
flux_season <- flux %>% 
  mutate(GPP = pmax(0,GPP)) %>% 
  group_by(Site, Year) %>% 
  filter(n()==4) %>% 
  mutate(GPPtot = sum(GPP)) %>% 
  group_by(Site, Month) %>% 
  summarise(GPP = mean(GPP/GPPtot)) %>% 
  group_by(Site) %>% 
  mutate(GPP = GPP/mean(GPP))

flux_rw_overlap <- read_csv("Data/flux_rw_overlap.csv")
flux_rw_overlap_sites <- unique(flux_rw_overlap$Site)

# Met data
met <- read_csv("Data/flux_met.csv")
# Aggregate met data at the three monthly scale
met_agg <- met %>% 
  left_join(flux_metadata %>% select(Site, Lat)) %>%
  mutate(Date = if_else(Lat<0, Date - duration(181,"days"), Date), # 6 months lag for southern sites
         Date = round_date(Date, "month")) %>%
  select(-Lat) %>%
  mutate(Date = Date + duration(32, "days")) %>% # 1 month lag to include previous year Dec
  mutate(Date = floor_date(Date, "3 months")) %>% 
  group_by(Site, Year = year(Date), Month = month(Date)) %>%        # Average within 3-months periods
  select(-Date) %>%                                                 #
  summarise(across(everything(), function(x) mean(x, na.rm = T)))   #

# Include 1 year lag (not used)
met_agg <- bind_rows(met_agg %>% mutate(Lag = "Current Year"), 
                     met_agg %>% mutate(Year = Year+1, Lag = "Previous Year"))

#### Merge flux data with rad and met data
flux_met <- flux %>% 
  filter(Site %in% flux_rw_overlap_sites) %>%
  left_join(met_agg) %>%
  filter(!is.na(tmin))

#### Prep data for partial correlations and linear models
m_data  <- flux_met %>%
  mutate(tvar = tmax-tmin) %>%
  mutate(DIFF = sqrt(DIFF)) %>%
  group_by(Site, Month, Lag) %>% 
  mutate(across(c(NEE:R,tmin:DIFF), std)) %>% 
  filter(!is.na(GPP), !is.na(PARTOT))

## Evaluate multicolinearity between predictors
# vif <- m_data %>% 
#   select(tmean, PDSI, PARTOT) %>%
#   group_by(Site, Month) %>% nest %>%
#   mutate(VIF = map(data, VIF)) %>%
#   unnest(cols = VIF)
# 
# vif %>% ggplot(aes(x = Var, y = VIF))+
#   geom_boxplot()+
#   scale_y_log10()+
#   geom_hline(yintercept = 5, linetype = 2)+
#   facet_wrap(~Month)

#### Run partial correlations
pr <- m_data %>%
  select(Site, Month, Lag, Dep = GPP, tmean, PDSI, PARTOT) %>%
  group_by(Site, Month, Lag) %>% 
  filter(n()>=6) %>% # Partial correlations need at least 2+n degrees of freedom, where n is the number of variables
  nest %>%  
  mutate(
    pr = map(data, function(x) pcor(x, method = "pearson")),
    pval = map(pr, function(x) x$p.value[-1,1,drop=F] %>% as.data.frame %>% rename(pval = Dep)),
    n = map(pr, function(x) x$n),
    pr = map(pr, function(x) x$estimate[-1,1,drop=F] %>% as.data.frame %>% rownames_to_column(var = "Var") %>% rename(pr = Dep))
  ) %>%
  select(-data) %>%
  unnest(cols = c(pval, n, pr)) %>%
  mutate(Var = factor(Var, levels = c("tmean", "PDSI", "PARTOT"))) %>% 
  left_join(flux_season) %>%
  mutate(pr = pr*GPP) %>%
  select(-GPP)

month <- unique(pr$Month)
season <- c("Winter", "Spring", "Summer", "Fall")
## Graph of the raw partial correlations results
# pr %>%
#   ungroup %>%
#   left_join(flux_metadata) %>%
#   filter(Lag == "Current Year") %>%
#   mutate(Fac = MAT) %>% 
#   mutate(Fac = cut(Fac, quantile(Fac, p=c(0,0.5,1)),include.lowest = T)) %>% 
#   mutate(Month = factor(Month, label = season)) %>%
#   mutate(Lag = factor(Lag, levels = c("Previous Year", "Current Year"))) %>% 
#   ggplot(aes(x = Month, y = pr, fill = factor(Month), alpha = Fac))+
#   geom_hline(yintercept = 0)+
#   stat_summary(fun = mean, show.legend = F, geom = "col", position = "dodge")+
#   stat_summary(fun.data = mean_se, col = "black", show.legend = F, geom = "errorbar", width = 0.5, position = position_dodge(width=0.9))+
#   facet_grid(Var~.)+
#   scale_fill_brewer(type = "div", palette = 1, direction = -1)+
#   scale_alpha_discrete(range = c(0.6,1))+
#   xlab("")+ylab("Partial correlation coefficient")+
#   # coord_cartesian(expand = T, ylim = c(-0.2, 0.5))+
#   theme(panel.grid.major.y = element_line(color = "grey", size = 0.2, linetype = 2))
# ggsave("Graphs/Partial correlation GPP~(tmean,PDSI,PARTOT,PARDF).pdf", width = 4, height = 4)

### Normalize partial correlations by accounting for the effect of site climate 
pr_lm <- pr%>%
  filter(Lag == "Current Year") %>%
  left_join(flux_metadata) %>%
  group_by(Var, Month, Lag) %>% 
  mutate(MACWD = sqrt(MACWD)) %>%
  mutate(across(c(MAT, MAR, MACWD), std)) %>% # Standardize predictors
  nest %>%
  mutate(m = map(data, function(x) lm(logitr(pr/2)*2~(MACWD+MAT), weights = n-4, data =x)), # division per two to have r vary between -1 and 1. Directly multiplying by 2 the result instead of estimated coefficients yields almost exactly the same result and simplifies the code
         coef = map(m, coefs.lmer),
         R2 = map(m, function(x) rsq::rsq(x)[1]) %>% unlist) %>%
  select(-data,-m) %>%
  unnest(cols = c(coef, R2))
# Save results for final graph
saveRDS(pr_lm, "Graphs/Final/pr GPP-climate.RDS")

# Plot result
pr_lm %>% 
  ungroup %>%
  filter(name %in% c("(Intercept)")) %>%
  mutate(estimate = logitr_inv(estimate)) %>%
  rowwise %>% 
  mutate(signif = signif(pval))%>%
  ggplot(aes(x = Month, y = estimate, fill = factor(Month), group = name))+
  geom_col(show.legend = F, position = 'dodge', aes(alpha = pval != ""))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_label(aes(label = signif, y = pmax(0,estimate)), show.legend = F, label.size = 0, fill = "transparent", size = 5, position = position_dodge(width=0.9))+
  facet_grid(Var~"GPP")+
  coord_cartesian(ylim = c(-0.15,0.65), expand = F)+
  scale_fill_brewer(type = "div", palette = 1, direction = -1)+
  scale_alpha_manual(values = c(1,1))+
  ylab("Partial correlation")+xlab(NULL)+
  theme(panel.grid.major.y = element_line(size = 0.2, linetype = 2, color = "grey"))
# ggsave("Graphs/lmer(pr GPPmayoct ~ MAT+MACWD+MAR+Group+1|IGBP) - Intercept.pdf", width = 4, height = 6)
