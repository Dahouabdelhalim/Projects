######
###### Effect of diffuse light on tree ring width
######

library(ppcor)
library(lme4)
library(MuMIn)
select <- dplyr::select

#### Load data
rw_metadata <- read_csv("Data/rw_metadata.csv") 
rw <- read_csv("Data/rw.csv")

flux_rw_overlap <- read_csv("Data/flux_rw_overlap.csv")
flux_rw_overlap_sites <- unique(flux_rw_overlap$RW_site)

# Met data
met <- read_csv("Data/rw_met.csv")
# Aggregate met data at the 3-monthly scale
met_agg <- met %>% 
  left_join(rw_metadata %>% select(Site, Lat)) %>%
  mutate(Date = if_else(Lat<0, Date - duration(181,"days"), Date), # 6 month lag for southern sites
         Date = round_date(Date, "month")) %>%                     #
  mutate(Date = Date + duration(32,"days")) %>% # 1 month shift for season definition
  select(-Lat) %>%
  mutate(Date = floor_date(Date, "3 months")) %>%                     # 3 months average
  group_by(Site, Year = year(Date), Month = factor(month(Date))) %>%  #
  select(-Date) %>%                                                   #
  summarise(across(everything(), function(x) mean(x, na.rm = T)))     #
# Add previous year (not used)
met_agg <- bind_rows(met_agg %>% mutate(Lag = "Current Year"), 
                     met_agg %>% mutate(Year = Year+1, Lag = "Previous Year"))


### Merge RW data with rad and met data
rw_met <- rw %>% 
  filter(Site %in% flux_rw_overlap_sites) %>% 
  left_join(met_agg) %>%
  filter(!is.na(tmin))

### Prep data for partial correlations and linear models
m_data  <- rw_met %>%
  group_by(Site, Month, Lag) %>% 
  mutate(across(!c(Year), std)) %>%
  filter(sum(!is.na(RWI) & !is.na(PARTOT)) > 11) 

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
  select(Site, Month, RWI, tmean, PDSI, PARTOT) %>%
  group_by(Site, Month, Lag) %>%
  filter(n()>=6) %>% # Partial correlations need at least 2+n degrees of freedom, where n is the number of variables
  nest %>%  
  mutate(pr = map(data, function(x) pcor(x, method = "pearson")), 
         pval = map(pr, function(x) x$p.value[-1,1,drop=F] %>% as.data.frame %>% rename(pval = RWI)),
         n = map(pr, function(x) x$n),
         pr = map(pr, function(x) x$estimate[-1,1,drop=F] %>% as.data.frame %>% rownames_to_column(var = "Var") %>% rename(pr = RWI))) %>%
  select(-data) %>%
  unnest(cols = c(pval, n, pr)) %>%
  mutate(Var = factor(Var, levels = c("tmean", "PDSI", "PARTOT")))

month <- unique(pr$Month)
season <- c("Winter", "Spring", "Summer", "Fall")
# Plot raw partial correlations
pr %>%
  ungroup %>%
  filter(Lag == "Current Year") %>% 
  mutate(Month = factor(Month, labels = season)) %>% 
  ggplot(aes(x = Month, y = pr, fill = Month))+
  geom_hline(yintercept = 0)+
  stat_summary(fun = mean, show.legend = F, geom = "col", position='dodge')+
  stat_summary(fun.data = mean_cl_normal, col = "black", show.legend = F, geom = "errorbar", width = 0.5, position = position_dodge(width=0.9))+
  facet_grid(Var~.)+
  scale_alpha_manual(values = c(0.7, 1))+
  scale_fill_brewer(type = "div", palette = 1, direction = -1)+
  xlab("")+ylab("Partial correlation coefficient")+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9), panel.grid.major.y = element_line(color = "grey", size = 0.2, linetype = 2))

## Normalize partial correlations by accounting for the effect of site climate 
pr_lm <- pr %>%
  filter(Lag == "Current Year") %>%
  left_join(rw_metadata) %>%
  group_by(Var, Month, Lag) %>% 
  mutate(MACWD = sqrt(MACWD)) %>%
  mutate(MAT = MAT-7.847) %>%             # Center on the mean climate of flux sites
  mutate(MACWD = MACWD-sqrt(240.71)) %>%  # to enable comparison of partial correlations
  mutate(across(c(MAT, MACWD), function(x) x/(sd(x, na.rm = T)))) %>% # divide by standard deviation
  nest %>%
  mutate(m = map(data, function(x) lm(logitr(pr)~(MACWD+MAT), weight = n-4, data =x)),
         coef = map(m, coefs.lmer),
         R2 = map(m, function(x) rsq::rsq(x)[1]) %>% unlist) %>%
  select(-data,-m) %>%
  unnest(cols = c(coef, R2)) %>% 
  mutate(name = factor(name, levels = c("(Intercept)", "MACWD", "MAT")))
# Save model result for final plot later on
saveRDS(pr_lm, "Graphs/Final/pr RWI-climate.RDS")

# Preliminary plot of results
pr_lm %>% 
  ungroup %>% 
  filter(name %in% c("(Intercept)")) %>%
  mutate(Month = factor(Month, labels = season)) %>%
  mutate(estimate = logitr_inv(estimate)) %>% 
  rowwise %>% 
  mutate(signif = signif(pval)) %>% 
  ggplot(aes(x = Month, y = estimate, fill = Month, group = name))+
  geom_col(show.legend = F, position = 'dodge', aes(alpha = Lag))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_label(aes(label = signif, y = pmax(0,estimate)), show.legend = F, label.size = 0, fill = "transparent", size = 5, position = position_dodge(width=0.9))+
  facet_grid(Var~"RWI")+
  coord_cartesian(ylim = c(-0.15,0.65), expand = F)+
  scale_fill_brewer(type = "div", palette = 1, direction = -1)+
  scale_alpha_manual(values = c(1,1))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype =2))+
  ylab("Partial correlation")+xlab(NULL)
