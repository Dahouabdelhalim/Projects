library(tidyverse)
library(sjPlot)
library(sjmisc)
library(ggeffects)
library(cowplot)

# load and combine data
spm <- list.files("input/records_through_time/specimen")
spm <- lapply(file.path("input/records_through_time/specimen", spm), "read_csv")
spm <- bind_rows(spm)

tot <- list.files("input/records_through_time/total")
tot <- lapply(file.path("input/records_through_time/total", tot), "read_csv")
tot <- bind_rows(tot)

inp <- bind_rows(tot, spm, .id = "data_type") %>% 
  mutate(data_type = recode(data_type, `1` = "total", `2` = "specimens"))

tmi <- read_csv("output/colonial_ties.csv")
tmi <- select(tmi, iso3, indep_year, Name)


# limit to time after independence
li <- unique(inp$country)
dat <- list()

for(i in 1:length(li)){
  tmi_sub <- tmi %>% filter(iso3 == li[i])
  
  sub <- inp %>% 
    filter(country == li[i]) %>% 
    filter(year > tmi_sub$indep_year) %>% 
    mutate(time_since_independence = year - tmi_sub$indep_year) %>% 
    mutate(country_name = tmi_sub$Name)
  
  dat[[i]] <- sub
}

dat <- bind_rows(dat)

# prepare relevant data
dat <- dat %>% 
  select(data_type, year, time_since_independence, country, country_name, records, records_domestic, records_coloniser) %>% 
  mutate(records_other = records - records_domestic - records_coloniser) %>% 
  mutate(fraction_domestic = records_domestic / records) %>% 
  mutate(fraction_domestic = ifelse(is.infinite(fraction_domestic), 0, fraction_domestic)) %>% 
  mutate(fraction_coloniser = records_coloniser/ records) %>% 
  mutate(fraction_coloniser = ifelse(is.infinite(fraction_coloniser), 0, fraction_coloniser))

# fit two mixed-effect models
dom <- dat %>% filter(data_type == "total")
dom_me <- glmer(as.matrix(dom[, c("records_domestic", "records")]) ~ 1 + time_since_independence + (1 + time_since_independence | country),
             family = binomial(link = "logit"),
             data = dom)

summary(dom_me)
cc <- confint(dom_me,parm="beta_")  ## slow (~ 11 seconds)
ctab <- cbind(est=fixef(dom_me),cc)
rtab <- exp(ctab)
print(rtab,digits=3)

coef(dom_me)
ranef(dom_me)

# The same for records from the former colonial suppressor
supp <- dat %>% filter(data_type == "specimens")
supp_me <- glmer(as.matrix(supp[, c("records_coloniser", "records")]) ~ 1 + time_since_independence + (1 + time_since_independence | country),
                family = binomial(link = "logit"),
                data = supp)
plot(supp_me)
summary(supp_me)
cc <- confint(supp_me, parm = "beta_")  ## slow (~ 11 seconds)
ctab <- cbind(est = fixef(supp_me), cc)
rtab <- exp(ctab)
print(rtab, digits = 3)

coef(supp_me)
ranef(supp_me)

# table
tab_model(dom_me)
tab_model(supp_me)

#Overall plot
# raw data
plo <- dat %>% 
  select(time_since_independence, fraction_coloniser, fraction_domestic, country) %>% 
  pivot_longer(contains("fraction"), names_to = "type", values_to = "fraction")  %>% 
  mutate(type = recode(type,
                       fraction_domestic = "A) Domestic records",
                       fraction_coloniser = "B) Preserved specimen from former colonial suppressor"))



ggplot(data = plo, aes(x = time_since_independence, y = fraction* 100))+
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black')+
  geom_point()+ 
  xlab("Time since independence")+
  ylab("Proportion of records")+
  facet_grid(type~.)+
  theme_bw()

ggsave(file = "output/supplementary_figure_raw_data_colonial_models.jpg", height = 10, width = 8)

# overall model
dat_dom <- ggpredict(dom_me, terms = c("time_since_independence"))
dat_supp <- ggpredict(supp_me, terms = c("time_since_independence"))

pred <- bind_rows(dat_dom, dat_supp, .id = "type") %>% 
  mutate(type = recode(type,
                       `1` = "A) Domestic records",
                       `2` = "B) Preserved specimen from former colonial suppressor"))

ggplot()+
  geom_point(data = plo, aes(x = time_since_independence, y = fraction ), alpha = 0.3, color = "grey")+
  geom_line(data = pred, aes(x = x, y = predicted ))+
  xlab("Years since independence")+
  ylab("Proportion of records")+
  facet_grid(type~.)+
  theme_bw()

ggsave(file = "output/supplementary_figure_overall_predictions.jpg", height = 10, width = 8)

# random effects
dom_rs <- plot_model(dom_me,
                     type="pred",
                     terms=c("time_since_independence","country"),
                     pred.type="re",
                     colors = "gs")+
  ggtitle("Domestic records")+
  xlab("Time since independence")+
  ylab("Proportion domestic")+
  theme_bw()

supp_rs <- plot_model(supp_me,
                      type="pred",
                     terms=c("time_since_independence","country"),
                     pred.type="re",
                     colors = "gs")+
  ggtitle("Preserved specimen from former colonial suppressor")+
  xlab("Years since independence")+
  ylab("Proportion of records")+
  theme_bw()

plot_grid(dom_rs, supp_rs,ncol = 1)

ggsave(file = "output/supplementary_figure_random_predictions.jpg", height = 10, width = 8)

# random slopes

# loop to fit one linear models and  produce output table
# li <- sort(unique(dat$country))
# out <- list()
# 
# pdf("output/regression_plots_countries.pdf")
# 
# for(i in 1:length(li)){
#   # get subset
#   sub <- filter(dat, country == li[i])
# 
#   # all records domestic
#   sub_dom <- filter(sub, data_type == "total")
#   mod_tot <- glm(as.matrix(sub_dom[, c("records_domestic", "records")])  ~ sub_dom$time_since_independence,
#                  family = binomial) %>%
#     tidy() %>%
#     mutate(country = li[i]) %>%
#     mutate(independence_year = tmi %>% filter(iso3 == li[i]) %>% pull(indep_year)) %>%
#     mutate(records = sum(filter(sub_dom, data_type == "total") %>% pull(records))) %>%
#     mutate(records_domestic = sum(filter(sub_dom, data_type == "total") %>% pull(records_domestic))) %>%
#     mutate(records_coloniser = sum(filter(sub_dom, data_type == "total") %>% pull(records_coloniser))) %>%
#     mutate(data_source = "all_records")%>%
#     mutate(dependent_variable = "domestic_occurrences")
# 
#   ## specimen only domestic
#   # mod_spm_dom <- lm(fraction_domestic ~ year, data = filter(sub, data_type == "specimens")) %>%
#   #   tidy() %>%
#   #   mutate(dependent_variable = "domestic_occurrences")
#   #
#   ## specimens colonial suppressor
#   sub_supp <- filter(sub, data_type == "specimens")
#   mod_spm_col <- glm(as.matrix(sub_supp[, c("records_coloniser", "records")])  ~ sub_supp$time_since_independence,
#                      family = binomial) %>%
#     tidy() %>%
#     mutate(dependent_variable = "coloniser_occurrences")
# 
#   # mod_spm <- bind_rows(mod_spm_dom, mod_spm_col)%>%
#   mod_spm <- mod_spm_col %>%
#     mutate(country = li[i]) %>%
#     mutate(independence_year  = tmi %>% filter(iso3 == li[i]) %>% pull(indep_year)) %>%
#     mutate(records = sum(filter(sub_supp, data_type == "specimens") %>% pull(records))) %>%
#     mutate(records_domestic = sum(filter(sub_supp, data_type == "specimens") %>% pull(records_domestic))) %>%
#     mutate(records_coloniser = sum(filter(sub_supp, data_type == "specimens") %>% pull(records_coloniser))) %>%
#     mutate(data_source = "specimens_records")
# 
#   #return output data
#   out_sub <- bind_rows(mod_tot,
#                        mod_spm)
# 
#   out[[i]] <- out_sub
# 
#  #plot
#   plo <- sub %>%
#     select(year, fraction_domestic, fraction_coloniser, data_type) %>%
#     pivot_longer(cols = contains("fraction_"), values_to = "fraction", names_to = "type") %>%
#     mutate(type = paste(data_type, type, sep = "_")) %>%
#     mutate(type = recode(type,
#                          total_fraction_domestic = "Any source, domestic",
#                          total_fraction_coloniser = "Any source, former colonial suppressor",
#                          specimens_fraction_domestic = "Specimens, domestic",
#                          specimens_fraction_coloniser = "Specimens, former colonial suppressor")) %>%
#     filter(type %in% c("Any source, domestic", "Specimens, former colonial suppressor"))
# 
#   plo <- ggplot()+
#     geom_point(data = plo,
#                aes(x = year, y = fraction))+
#     geom_smooth(data = plo,
#                 aes(x = year,
#                     y = fraction),
#                 method = "glm",
#                 method.args = list(family = "binomial"),
#                 se = FALSE)+
#     ylab("Fraction of records")+
#     xlab("Year")+
#     theme_bw()+
#     theme(legend.position = "none")+
#     ggtitle(paste(li[i], unique(sub$country_name), sep = " - "))+
#     facet_wrap(type~., scales = "free", ncol = 1)
#   print(plo)
# }
# dev.off()
# 
# 
# # add false positive detection corrections
# out <- bind_rows(out) %>%
#   filter(term != "(Intercept)") %>%
#   mutate(bf_sig = p.value < (0.05 / (nrow(.)/2))) %>%
#   mutate(bh_sig = p.adjust(p.value, method = "BH", n = 80) < 0.05)
# 
# 
# # generate output data.frame and save to disk
# write_csv(out, file = "output/supplement_results_country_regressions.csv")
# 
# # relevant stats for the manuscript
# ## number of regression with a significant effect
# ###Benjaminiâ€“Hochberg
# out %>%
#   mutate(effect_direction = ifelse(estimate <0, "negative", "positive")) %>%
#   group_by(dependent_variable,bh_sig, effect_direction) %>%
#   count()
# 
# ### Bonferroni
# out %>%
#   mutate(effect_direction = ifelse(estimate <0, "negative", "positive")) %>%
#   group_by(dependent_variable, bf_sig ,effect_direction) %>%
#   count()
# 
# # mean effect sizes of significant cases
# out %>%
#   group_by(dependent_variable) %>%
#   filter(bf_sig) %>%
#   summarize(mean(estimate))
# 