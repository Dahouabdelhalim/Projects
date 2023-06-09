################################################################################
#                           Conduct survival analysis                          #
#                                                                              #
#                                 Eli Strauss                                  #
#                                                                              #
#                                                                              #
#                               February 2019                                  #
################################################################################


######Load packages and set global options
rm(list = ls())
library(dplyr)
library(here)
library(survival)
library(survminer)
library(coxme)
library(viridis)
library(sjPlot)
library(lme4)
options(stringsAsFactors = FALSE)
options(scipen = 8)

load('04.cub_dev_vars.RData')
source('00.define_functions.R')

### Color scheme for plotting
colors <- viridis(5)[c(1,4)]
colors4 <- viridis(5)[c(1:4)]

####Calculate end date
all_cub_survival <- list()
for(i in 1:nrow(cub_dev_vars)){
  ID <- cub_dev_vars$id[i]
  if(is.na(cub_dev_vars$disappeared[i])){
    #times <- as.numeric(cub_dev_vars$last.seen.sessions[i] - cub_dev_vars$birthdate[i])
    times <- as.numeric(as.Date('2019-06-30') - cub_dev_vars$birthdate[i])
    events <- 0
    #if(is.na(times)) next
  }else{
    times <- as.numeric(cub_dev_vars$disappeared[i] - cub_dev_vars$birthdate[i])
    events <- 1
  }
  cub_dat <- cub_dev_vars[i,]
  #if(length(events)-1) cub_dat <- rbind(cub_dat,cub_dev_vars[i,])
  all_cub_survival[[i]] <- data.frame(time = times, etype = events, cub_dat)
}
all_cub_survival <- do.call(rbind, all_cub_survival)

### Remove individuals without time, etype, end_diff, and mom_rank
all_survival <- filter(all_cub_survival, !is.na(time), !is.na(etype), !is.na(end_diff), !is.na(mom_rank),
                           birthdate >= "1988-06-26")

##### Right censor data for males that survive to 2 years old at 2####
all_survival[!is.na(all_survival$time) & all_survival$time >= (365*2) & all_survival$sex == 'm','etype'] <- 0
all_survival[!is.na(all_survival$time) & all_survival$time >= (365*2) & all_survival$sex == 'm','time'] <- 365*2

#### Right censor hyenas who fission to Talek East clan (not well studied) at date of fission ####
east.membership <- read.csv('data/ClanMembership.csv')
east.membership <- filter(east.membership, Membership == 'e' |
                            Mom %in% filter(east.membership, Membership == 'e'))
all_survival[all_survival$id %in% east.membership$ID & (is.na(all_survival$disappeared) | all_survival$disappeared >= '2000-01-01'),'etype'] <- 0
all_survival[all_survival$id %in% east.membership$ID & (is.na(all_survival$disappeared) | all_survival$disappeared >= '2000-01-01'),'time'] <- as.Date('2000-01-01') - 
  all_survival[all_survival$id %in% east.membership$ID & (is.na(all_survival$disappeared) | all_survival$disappeared >= '2000-01-01'),'birthdate']


all_survival$mom_survive_to_2 <- factor(all_survival$mom_survive_to_2, levels = c(1,0),
                                        labels = c('survived', 'dead'))

###Restrict to hyenas that survive to graduation
all_grad <- filter(all_survival, survive_to_grad == 1)
all_grad$age <- all_grad$time/365


#### Some descriptives ####
mean(all_grad$end_diff)
sd(all_grad$end_diff)
sum(!is.na(all_grad$end_diff))

mean(all_grad$postgrad_diff, na.rm = TRUE)
sd(all_grad$postgrad_diff, na.rm = TRUE)
sum(!is.na(all_grad$postgrad_diff))

mean(all_grad$adult_diff, na.rm = TRUE)
sd(all_grad$adult_diff, na.rm = TRUE)
sum(!is.na(all_grad$adult_diff))

#### Center variables
all_grad$end_diff_centered <- scale(all_grad$end_diff)
all_grad$num_intx_centered <- scale(all_grad$num_intx)
all_grad$end_obs_centered <- scale(all_grad$end_obs)
all_grad$postgrad_diff_centered <- scale(all_grad$postgrad_diff)
all_grad$adult_diff_centered <- scale(all_grad$adult_diff)

### Dichotomize variables
all_grad$diff_class <- cut(all_grad$end_diff, breaks = c(-1000, 0, 1000), labels = c('Elo < expected', 'Elo ≥ expected'), right = FALSE)
all_grad$diff_class <- factor(all_grad$diff_class, levels = c('Elo ≥ expected', 'Elo < expected'))

all_grad$rank_class <- cut(all_grad$mom_rank, breaks = c(-1.1, 0, 1.1), labels = c('low rank', 'high rank'))
all_grad$rank_class <- factor(all_grad$rank_class, levels = c('high rank', 'low rank'))

all_grad$obs_class <- cut(all_grad$end_obs, breaks = c(-1000, mean(all_grad$end_obs), 1000), labels = c('below average', 'above average'))
all_grad$obs_class <- factor(all_grad$obs_class, levels = c('above average', 'below average'))

#### Explore different models relating elo deviance at den independence and survival####
## Models include number of interactions, rank class (high/low), mom survive to 2

##
elo.dev.mod <- coxme(Surv(age, etype) ~ end_diff_centered + num_intx_centered + mom_rank + mom_survive_to_2 + (1|clan), data = all_grad)
diff.class.mod <- coxme(Surv(age, etype) ~ diff_class + num_intx_centered + mom_rank + mom_survive_to_2 + (1|clan), data = all_grad)
end.obs.mod <- coxme(Surv(age, etype) ~ end_obs_centered + num_intx_centered + mom_rank + mom_survive_to_2 + (1|clan), data = all_grad)
obs.class.mod <- coxme(Surv(age, etype) ~ obs_class + num_intx_centered + mom_rank + mom_survive_to_2 + (1|clan), data = all_grad)
both.class.mod <- coxme(Surv(age, etype) ~ diff_class + obs_class + num_intx_centered + mom_rank + mom_survive_to_2 + (1|clan), data = all_grad)
null.mod <- coxme(Surv(age, etype) ~ num_intx_centered + mom_rank + mom_survive_to_2 + (1|clan), data = all_grad)
MuMIn::AICc(elo.dev.mod, end.obs.mod, obs.class.mod, both.class.mod, diff.class.mod, null.mod)

MuMIn::AICc(elo.dev.mod, diff.class.mod)[1,2] - MuMIn::AICc(elo.dev.mod, diff.class.mod)[2,2]
MuMIn::AICc(obs.class.mod, diff.class.mod)[1,2] - MuMIn::AICc(obs.class.mod, diff.class.mod)[2,2]
MuMIn::AICc(end.obs.mod, diff.class.mod)[1,2] - MuMIn::AICc(end.obs.mod, diff.class.mod)[2,2]
MuMIn::AICc(null.mod, diff.class.mod)[1,2] - MuMIn::AICc(null.mod, diff.class.mod)[2,2]


#### Build primary model ####
mod <- coxme(Surv(age, etype) ~ diff_class + num_intx_centered + mom_rank + mom_survive_to_2 + (1|clan), data = all_grad)
mod.without.numintx <- coxme(Surv(age, etype) ~ diff_class + mom_rank + mom_survive_to_2 + (1|clan), data = all_grad)
MuMIn::AICc(mod, mod.without.numintx)

primary.mod <- mod
summary(primary.mod)


table1 <- cox.tab(primary.mod, fixef.labs = c('Elo deviance (< expected)',
                                              'Number of interactions',
                                              'Maternal rank',
                                              'Maternal death before adulthood (dead)'))
table1.caption <- c('Table 1. Cox mixed-effects model of survival as a function of Elo-deviance at den independence and other covariates')


### Primary model with only females for supplemental material #### 
all_grad_females <- filter(all_grad, sex == 'f')
mod.female <- coxme(Surv(age, etype) ~ diff_class + num_intx_centered + mom_rank + mom_survive_to_2 + (1|clan), data = all_grad_females)
summary(mod.female)

table2 <- cox.tab(mod.female, fixef.labs = c('Elo deviance (< expected)',
                                             'Number of interactions',
                                             'Maternal rank',
                                             'Maternal death before adulthood (dead)'))
table2.caption <-  c('Table 2. Cox mixed-effects model of survival in females as a function of Elo-deviance at den independence and other covariates')
#########

p <- ggsurvplot(surv_fit(Surv(age, etype) ~ diff_class, data = all_grad), conf.int = F, pval = F, data = all_grad, risk.table = F,
                xlab = 'Age (Years)', legend.labs = c('Expected and above', 'Below expected'), break.time.by = 1,
                palette = colors, size = 1)

cairo_pdf('plots/5_deviance_at_den_indpendence.pdf', 4, 4)
ggplot(data = p$data.survplot, aes(x = time, y = surv, color = strata)) + 
  geom_step(size = 1) + 
  theme_survminer() +
  theme(legend.position = c(0.6, 0.8))+
  xlab('Age (Years)')+
  ylab('Survival probability')+
  scale_colour_manual(name = " ",
                      labels = c('Elo ≥ expected', 'Elo < expected'),
                      values = colors)
dev.off()


#### Are there interaction effects? ####   # no
primary.mod.intx <- coxme(Surv(age, etype) ~ diff_class * mom_rank + mom_rank * mom_survive_to_2 + diff_class*mom_survive_to_2 + num_intx_centered + (1|clan), data = all_grad)
MuMIn::AICc(primary.mod, primary.mod.intx)


### Model with first adult rank instead of maternal rank ###
mod.first.rank <- coxme(Surv(age, etype) ~ diff_class + num_intx_centered + first_adult_rank + mom_survive_to_2 + (1|clan), data = all_grad)
summary(mod.first.rank)

table3 <- cox.tab(mod.first.rank, fixef.labs = c('Elo deviance (< expected)',
                                              'Number of interactions',
                                              'First adult rank',
                                              'Maternal death before adulthood (dead)'))
table3.caption <- c('Table 3. Cox mixed-effects model of survival as a function of Elo-deviance at den independence using first adult rank rather than maternal rank')


#### Deviance at adulthood
all_grad_postgrad <- filter(all_grad, (clan %in% c('talek', 'happy.zebra') & birthdate <= (as.Date('2016-06-19') - 365*2)) |
         (clan == 'serena.n' & birthdate <= (as.Date('2016-12-31') - 365*2)) |
         (clan == 'serena.s' & birthdate <= (as.Date('2017-03-05') - 365*2)))


all_grad_postgrad$postgrad_diff_class <- cut(all_grad_postgrad$postgrad_diff, breaks = c(-1000, 0, 1000), labels = c('Elo < expected', 'Elo ≥ expected'), right = FALSE)
all_grad_postgrad$postgrad_diff_class <- factor(all_grad_postgrad$postgrad_diff_class, levels = c('Elo ≥ expected', 'Elo < expected'))
all_grad_postgrad$postgrad_intx_centered <- scale(all_grad_postgrad$postgrad_intx)

coxme(Surv(age, etype) ~ postgrad_diff_class + postgrad_intx_centered + mom_rank + (1|clan), data = all_grad_postgrad) %>%
  cox.tab()


#### Deviance at 3
all_grad_adult <- filter(all_grad, (clan %in% c('talek', 'happy.zebra') & birthdate <= (as.Date('2016-06-19') - 365*3)) |
                              (clan == 'serena.n' & birthdate <= (as.Date('2016-12-31') - 365*3)) |
                              (clan == 'serena.s' & birthdate <= (as.Date('2017-03-05') - 365*3)))

all_grad_adult$adult_diff_class <- cut(all_grad_adult$adult_diff, breaks = c(-1000, 0, 1000), labels = c('Elo < expected', 'Elo ≥ expected'), right = FALSE)
all_grad_adult$adult_diff_class <- factor(all_grad_adult$adult_diff_class, levels = c('Elo ≥ expected', 'Elo < expected'))
all_grad_adult$adult_intx_centered <- scale(all_grad_adult$adult_intx)

coxme(Surv(age, etype) ~ adult_diff_class + adult_intx_centered + mom_rank + (1|clan), data = all_grad_adult) %>%
  cox.tab()


##### Are deviance measures at different life-history stages correlated?
cor.test(all_grad_postgrad$postgrad_diff_centered, all_grad_postgrad$end_diff_centered, use = 'complete.obs')
cor.test(all_grad_adult$adult_diff_centered, all_grad_adult$end_diff_centered, use = 'complete.obs')

#### Stratify by rank and Elo deviance together
all_grad$both_class <- paste0(all_grad$diff_class, all_grad$rank_class)
all_grad$both_class <- factor(all_grad$both_class, levels = c('Elo ≥ expectedhigh rank', 'Elo < expectedhigh rank', 'Elo ≥ expectedlow rank', 'Elo < expectedlow rank'))
coxme(Surv(time = age, event = etype) ~ both_class + mom_survive_to_2 + (1|clan), data = all_grad)
p <- ggsurvplot(surv_fit(Surv(age, etype) ~ rank_class + diff_class, data = all_grad), conf.int = F, pval = F, data = all_grad, risk.table = F,
           xlab = 'Age (Years)', break.time.by = 1,
           color = 'strata',
           linetype = 'strata')

levs <- levels(p$data.survplot$strata)
levs
p$data.survplot$diff_class <- ifelse(as.numeric(p$data.survplot$strata) %in% c(1,3), 'Elo ≥ expected', 'Elo < expected')
p$data.survplot$rank_class <- ifelse(as.numeric(p$data.survplot$strata) %in% c(1,2), 'High rank', 'Low Rank')

cairo_pdf('plots/5_rank_learning_commparison.pdf', 5, 5)
ggplot(data = p$data.survplot, aes(x = time, y = surv, color = strata, linetype = strata)) +
  geom_step(size = 1) +
  theme_survminer() +
  scale_colour_manual(name = " ",
                      labels = c("\\nHigh rank\\nElo ≥ expected\\n", "\\nHigh rank\\nElo < expected\\n", "\\nLow rank\\nElo ≥ expected\\n", "\\nLow rank\\nElo < expected\\n"),
                      values = c(colors[1], colors[2], colors[1], colors[2])) +
  scale_linetype_manual(name = " ",
                     labels = c("\\nHigh rank\\nElo ≥ expected\\n", "\\nHigh rank\\nElo < expected\\n", "\\nLow rank\\nElo ≥ expected\\n", "\\nLow rank\\nElo < expected\\n"),
                     values = c(1, 1, 3, 3))+
  theme(legend.position = c(0.75, 0.65))+
  xlab('Age (Years)')+
  ylab('Survival probability')

dev.off()




### Are adverse effects cumulative?
all_grad$adversity_count <- 0
all_grad[all_grad$diff_class == 'Elo < expected',]$adversity_count <- all_grad[all_grad$diff_class == 'Elo < expected',]$adversity_count + 1
all_grad[all_grad$rank_class == 'low rank',]$adversity_count <- all_grad[all_grad$rank_class == 'low rank',]$adversity_count + 1
all_grad[all_grad$mom_survive_to_2 == 'dead',]$adversity_count <- all_grad[all_grad$mom_survive_to_2 == 'dead',]$adversity_count + 1



cumulative <- coxme(Surv(time = age, event = etype) ~ adversity_count + (1|clan), data = all_grad)

table4 <- cox.tab(cumulative, fixef.labs = c('Number of adverse conditions'))
table4.caption <- c('Table 4. Cox mixed-effects model of survival as a function of the number of adverse conditions')



indiv <- coxme(Surv(time = age, event = etype) ~ mom_survive_to_2 + rank_class + diff_class + (1|clan), data = all_grad)
indiv.intx <- coxme(Surv(time = age, event = etype) ~ mom_survive_to_2 * rank_class * diff_class + (1|clan), data = all_grad)

MuMIn::AICc(cumulative, indiv, indiv.intx)
MuMIn::AICc(cumulative, indiv)[1,2] - MuMIn::AICc(cumulative, indiv)[2,2]
MuMIn::AICc(indiv.intx, indiv)[1,2] - MuMIn::AICc(indiv.intx, indiv)[2,2]

p <- ggsurvplot(surv_fit(Surv(age, etype) ~ adversity_count, data = all_grad), conf.int = F, pval = F, data = all_grad, risk.table = F,
                xlab = 'Age (Years)', break.time.by = 1,
                palette = colors4, size = 1)

cairo_pdf(file = here('plots/5_adverse_conditions.pdf'), 4, 4)

ggplot(data = p$data.survplot, aes(x = time, y = surv, color = strata)) + 
  geom_step(size = 1) + 
  theme_survminer() +
  theme(legend.position = c(0.8, 0.7), legend.text.align = 1)+
  xlab('Age (Years)')+
  ylab('Survival probability')+
  scale_colour_manual(name = "# of adverse\\nconditions",
                      labels = 0:3,
                      values = colors4)

dev.off()


#### Lifetime reproductive success - add lifespan as requested by reviewers
lrs.dat <-filter(all_grad, etype == 1, sex == 'f')
lrs.dat$time.scaled <- scale(lrs.dat$time)
lrs.mod <- glmer(data = lrs.dat,
              lrs ~ diff_class + num_intx_centered + mom_rank + mom_survive_to_2 + (1|clan),
              family = 'poisson')
summary(lrs.mod)

table5 <- glm.tab(lrs.mod, fixef.labs = c('Intercept',
                                            'Elo deviance (< expected)',
                                            'Number of interactions',
                                            'Maternal rank',
                                            'Maternal death before adulthood (dead)'))
table5.caption <- c('Table 5. Poisson GLMM of lifetime reproductive success as a function of Elo-deviance at den independence and other covariates')

lrs.mod.lifespan <- glmer(data = lrs.dat,
                 lrs ~ diff_class + num_intx_centered + mom_rank + mom_survive_to_2 + time.scaled + (1|clan),
                 family = 'poisson')
summary(lrs.mod.lifespan)

table6 <- glm.tab(lrs.mod.lifespan, fixef.labs = c('Intercept',
                                                    'Elo deviance (< expected)',
                                                    'Number of interactions',
                                                    'Maternal rank',
                                                    'Maternal death before adulthood (dead)',
                                                    'Lifespan'))
table6.caption <- c('Table 6. Poisson GLMM of lifetime reproductive success as a function of Elo-deviance at den independence, lifespan, and other covariates')

lrs.mod.postgrad <- glmer(data = filter(all_grad_postgrad, etype == 1, sex == 'f'),
              lrs ~ postgrad_diff_class + mom_rank + postgrad_intx_centered + (1|clan),
              family = 'poisson')
summary(lrs.mod.postgrad)

lrs.mod.adult <- glmer(data = filter(all_grad_adult, etype == 1, sex == 'f'),
                       lrs ~ adult_diff_class + mom_rank + adult_intx_centered + (1|clan),
                     family = 'poisson')
summary(lrs.mod.adult)

cairo_pdf(file = here('plots/5_lrs.pdf'), 6, 4)
ggplot(data = filter(all_grad, etype == 1, sex == 'f'), aes(x = both_class, y = lrs, color = both_class, linetype = both_class, shape = both_class))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(height = 0, width = 0.3, alpha = 0.7)+
  theme_survminer()+
  theme(legend.position = 'bottom', axis.text.x = element_blank(), axis.title.x = element_blank())+
  scale_colour_manual(name = " ",
                      labels = c("\\nHigh rank\\nElo ≥ expected\\n", "\\nHigh rank\\nElo < expected\\n", "\\nLow rank\\nElo ≥ expected\\n", "\\nLow rank\\nElo < expected\\n"),
                      values = c(colors[1], colors[2], colors[1], colors[2])) +
  scale_linetype_manual(name = " ",
                        labels = c("\\nHigh rank\\nElo ≥ expected\\n", "\\nHigh rank\\nElo < expected\\n", "\\nLow rank\\nElo ≥ expected\\n", "\\nLow rank\\nElo < expected\\n"),
                        values = c(1, 1, 3, 3))+
scale_shape_manual(name = " ",
                      labels = c("\\nHigh rank\\nElo ≥ expected\\n", "\\nHigh rank\\nElo < expected\\n", "\\nLow rank\\nElo ≥ expected\\n", "\\nLow rank\\nElo < expected\\n"),
                      values = c(19, 17, 19, 17))+
  ylab('Lifetime reproductive success')
dev.off()



#### Save data
rank.acquisition <- all_grad
save(rank.acquisition, file = '06.rank.acquisition.RData')


#### Output tables ####
tab_dfs(x = list(table1[[1]], table2[[1]], table3[[1]], table4[[1]], table5[[1]], table6[[1]]),
        titles = list(table1.caption, table2.caption, table3.caption, table4.caption, table5.caption, table6.caption),
        footnotes = list(table1[[2]], table2[[2]], table3[[2]], table4[[2]], table5[[2]], table6[[2]]),
        show.footnote = TRUE, file = '~/Documents/Writing/Papers/2019 StraussShizukaHolekamp Juvenile rank acquisition influences fitness independent of adult rank/SecondRevision/SupplementalMaterials2_Tables.doc')

