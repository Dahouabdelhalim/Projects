################################################################################
#                              Rank acquisition                                #
#                                                                              #
#                                 Eli Strauss                                  #
#                                                                              #
#                                                                              #
#                               February 2019                                  #
################################################################################

######Load packages and set global options
rm(list = ls())
library(dplyr)
library(aniDom)
library(here)
library(ggplot2)
library(survminer)
library(viridis)
library(ggpubr)
options(stringsAsFactors = FALSE)

set.seed(1989)

load('02.cohortInfo.RData')
load('06.rank.acquisition.RData')
load('data/raw_data.RData')
source('00.define_functions.R')
colors <- viridis(5)[c(1,4)]


##### Rank learning - calculate Elo deviance every month
rank_learning_builder <- list()
counter <- 1
for(i in 1:nrow(cohortInfo)){
  id.intx <- filter(aggsFull, aggressor == cohortInfo$id[i] | recip == cohortInfo$id[i])
  if(!nrow(id.intx))next
  max.month <- ceiling(as.numeric(max(id.intx$date, na.rm = TRUE) - cohortInfo$birthdate[i])/30)
  if(max.month <= 0)next
  for(m in seq(from = 1, to = max.month, by = 1)){
    intx <- filter(id.intx, aggressor == cohortInfo$id[i] | recip == cohortInfo$id[i],
                   agg.rank != recip.rank, 
                   date > cohortInfo$birthdate[i]+(30*(m-1)),
                   date <= cohortInfo$birthdate[i] + (30*(m)))
    
    
    
    if(nrow(intx) < 1) next
    rank_learning_builder[[counter]] <- data.frame(deviance = elo_deviance(intx, cohortInfo$id[i]),
                                             age.months = m,
                                             id = cohortInfo$id[i],
                                             sex = cohortInfo$sex[i])
    counter <- counter + 1
  }
}
rank_learning <- do.call(rbind, rank_learning_builder)

rank_learning_summary <- rank_learning %>% 
  group_by(age.months) %>%
  summarize(sd.deviance = sd(deviance), num.ids = length(id)) %>%
  filter(num.ids >= 20)

segmented::davies.test(lm(data = rank_learning_summary, sd.deviance ~ age.months), seg.Z = ~age.months)
seg <- segmented::segmented(lm(data = rank_learning_summary, sd.deviance ~ age.months), seg.Z = ~age.months,
                             psi = c(18))

## Model of before rank acquisition
summary(lm(data = filter(rank_learning_summary, age.months <= seg$psi[2]), sd.deviance ~ age.months))

## Model of before rank acquisition
summary(lm(data = filter(rank_learning_summary, age.months > seg$psi[2]), sd.deviance ~ age.months))


pdf('plots/6_rank_acquisition_timing.pdf', 3.5,3.5)
ggplot(rank_learning_summary, aes(x = age.months, y = sd.deviance)) + 
  geom_jitter() +
  theme_survminer()+
  xlab('Age (months)')+
  ylab('Std. deviation of Elo deviance')+
  geom_vline(xintercept = 18, lty = 2)+
  geom_vline(xintercept = seg$psi[2], lty = 3)+
  geom_smooth(data = filter(rank_learning_summary, age.months <= round(seg$psi[2])), method = 'lm',
              col = 'red', fill = 'red')+
  geom_smooth(data = filter(rank_learning_summary, age.months >= round(seg$psi[2])), method = 'lm',
              col = 'red', fill = 'red')
dev.off()
  

#### Do juveniles acquire ranks as expected?

first.ranks <- tblFemaleRanks %>%
  group_by(id) %>%
  summarize(first.year = min(year), 
            first.rank = as.numeric(stan_rank[which.min(year)]))

first.ranks$mom <- left_join(first.ranks, tblHyenas, by = 'id')$mom
first.ranks$mom.rank <- as.numeric(left_join(first.ranks, tblFemaleRanks, 
                                             by = c('mom' = 'id', 'first.year' = 'year'))$stan_rank)
first.ranks$diff_class <- left_join(first.ranks, rank.acquisition, by = 'id')$diff_class
first.ranks <- filter(first.ranks, !is.na(diff_class))
first.ranks$mom.cub.diff <- first.ranks$mom.rank - first.ranks$first.rank
first.ranks$adult_expected_first_rank <- left_join(first.ranks, rank.acquisition, by = 'id')$adult_expected_first_rank
first.ranks$first_adult_rank <- left_join(first.ranks, rank.acquisition, by = 'id')$first_adult_rank
first.ranks$adult_first_rank_diff <- left_join(first.ranks, rank.acquisition, by = 'id')$adult_first_rank_diff

cor.test(first.ranks$first.rank, first.ranks$mom.rank)

rank.acquisition$mri <- ifelse(rank.acquisition$adult_first_rank_diff == 0, 'mri',
                               ifelse(rank.acquisition$adult_first_rank_diff < 0, 'below', 'above'))

mri.table <- table(rank.acquisition[,c('mri', 'diff_class')])
chisq.test(mri.table)


#cairo_pdf('plots/6_mri.pdf', 3.5, 3.5)
mri <- ggplot(first.ranks, aes(x = mom.rank, y = first.rank, col = diff_class)) +
  geom_point(alpha = 0.7) +
  theme_survminer()+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  scale_color_manual(values = colors,
                     name = '',
                     labels = c('Elo ≥ expected',
                                'Elo < expected')
                     )+
  ylab('Rank at onset of adulthood')+
  xlab('Mother\\'s rank')+
  theme(legend.position = c(0.35,0.9), legend.key.size = unit(4, 'pt'))
  # geom_text(data = data.frame(), aes(-0.9, 1, label = 'b)'), inherit.aes = FALSE,
  #           fontface = 2)
#dev.off()

#### Histogram of scores
#cairo_pdf('plots/5_hist_of_deviance_at_den_indpendence.pdf', 4, 4)
hist <- ggplot(data = rank.acquisition, aes(x = end_diff, fill = diff_class))+
  geom_histogram(bins = 30)+theme_survminer()+
  xlab('Deviance at den independence')+
  theme(legend.position = c(0.27, 0.8))+
  scale_fill_manual(name = " ",
                      labels = c('Elo ≥ expected', 'Elo < expected'),
                      values = colors)+
  ylab('')
  # geom_text(data = data.frame(), aes(-315, 170, label = 'a)'), inherit.aes = FALSE,
  #           fontface = 2)
#dev.off()

cairo_pdf('plots/figure_1.pdf', 8, 4)
ggarrange(hist, mri, labels = c('a)', 'b)'))
dev.off()

   