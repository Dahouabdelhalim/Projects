library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# read excel, grab only the count part of the matrix and sum every 3 rows together to combine
# each slide to a single count for plotting
df <- readxl::read_excel('Counts Fungi Probes 445_01 and 445_02V2.xlsx') %>% 
  dplyr::select(9:23) %>% 
  dplyr::mutate_all(as.integer) %>% 
  dplyr::mutate('slide' = rep(1:12, each=3)) %>%
  dplyr::group_by(slide) %>% 
  dplyr::summarise_all(sum) %>%
  dplyr::select(2:16)

plot_df <- df %>% dplyr::select(2:15) %>% 
  replace(is.na(.), 0) %>%
  tidyr::gather(key = 'Interaction', value = 'Observations') 


sum_df <- df %>% dplyr::mutate(N=rowSums(.)) %>% dplyr::select(2:16) %>% summarise_all(sum)
N <- sum_df %>% dplyr::select(N)
N <- N[[1]]

sum_df <- sum_df %>% 
  dplyr::select(1:14) %>% 
  tidyr::gather(key = 'Interaction', value='Sum')

# bonferroni corrected p-values
binomial_test <- function(X, N, prop){
  p_values <- rep(0, length(X))
  for (i in 1:length(X)){
    p_values[i] <- binom.test(X[i], N, prop, alternative = 'greater')$p.value * length(X)
  }
  return(p_values)
}

binom_df <- sum_df %>% 
  dplyr::mutate(p_value = binomial_test(Sum, N, 0.01)) %>%
  dplyr::mutate(signif = dplyr::if_else(p_value <= 0.05, 'p <= 0.05', 'p > 0.05'))

plot_df <- dplyr::full_join(plot_df, binom_df, by='Interaction')

p <- ggplot(plot_df, aes(x=reorder(Interaction, Observations, FUN=median), y=Observations, color=signif)) + 
  geom_boxplot() + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(color='Bonferroni Corrected\\n1-Sided Exact\\nBinomial Test (0.01)') +
  xlab('NCLC1:Diatom Interaction') + 
  ylab('Observations per Slide') +
  ggtitle('NCLC1:Diatom Interactions Observed using FISH Microscopy') 

p

ggsave('figure_2f.pdf', p)
ggsave('figure_2f.svg', p)
