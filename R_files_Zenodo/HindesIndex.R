## Juvenile join-leave analysis
## Ali Galezo

options(stringsAsFactors = FALSE)
library(tidyverse)

## Set seed for reproducibility.
set.seed(123)

## Load data.
events <- read.csv("juvenile_jlevents.csv")

## Set number of simulations.
num_sims <- 1000

## Set cutoff for number of interactions a dyad needs to be included in analysis.
dyad_cutoff <- 5

## Identify all dyads with at least 5 interactions.
all <-
  events %>%
  group_by(event_type) %>%
  count(actor, recipient) %>%
  ungroup() %>%
  complete(nesting(actor, recipient), event_type, fill = list(n = 0)) %>%
  mutate(dyad = ifelse(actor > recipient, paste(actor, recipient, sep = "_"), paste(recipient, actor, sep = "_"))) %>%
  group_by(dyad) %>%
  mutate(total = sum(n)) %>%
  filter(total >= dyad_cutoff) %>%
  spread(key = event_type, value = n) %>%
  group_by(dyad) %>%
  arrange(dyad) %>%
  mutate(joins = sum(JOIN),
         leaves = sum(LEAVE),
         join_prop = ifelse(joins != 0, JOIN / joins, 0),
         leave_prop = ifelse(leaves != 0, LEAVE / leaves, 0)) %>%
  ungroup()

## List of animals represented in dyads with at least 5 interactions, and their sexes.
animals <- unique(c(all$actor, all$recipient))
sex <- data.frame(id = animals, sex = NA)
for (i in sex$id){
  sex[sex$id == i, "sex"] <- unique(c(events[events$actor == i, "sex_actor"], events[events$recipient == i, "sex_recipient"]))
}

## Function to calculate Hinde's index for male-female dyads.
hindeCalculator <- function(sex){
  matrix <-
    all %>%
    merge(x = ., y = sex, by.x = "actor", by.y = "id", all.x = TRUE) %>%
    merge(x = ., y = sex, by.x = "recipient", by.y = "id", all.x = TRUE, suffixes = c("_actor","_recipient")) %>%
    filter(sex_actor != sex_recipient) %>%
    filter(sex_actor == "MALE") %>%
    mutate(hinde = join_prop - leave_prop)
  return(list(mean = mean(matrix$hinde),
              se = sd(matrix$hinde) / sqrt(nrow(matrix)),
              matrix = matrix))
}

## Observed Hinde's index: mean Hinde's index for all actual male-female pairs with at least 5 interactions.
obs <- hindeCalculator(sex)

## Expected Hinde's index: 1) Randomize sexes of animals. 2) Calculate mean Hinde's index for simulated male-female pairs with at least 5 interactions.
exp <- replicate(num_sims,
                 hindeCalculator(data.frame(id = sex$id,
                                            sex = sample(sex$sex))), 
                 simplify = FALSE)

## P-value.
extreme <- sum(obs$mean < unlist(lapply(exp, "[[", "mean")))
(extreme + 1) / (num_sims +1)

## Summary stats.
length(unique(all$dyad)) # number of dyads with 5+ interactions
all %>% distinct(dyad, joins, leaves) %>% summarize(joins = sum(joins), leaves = sum(leaves)) # number of joins and leaves between dyads with 5+ interactions
nrow(obs$matrix) # number of male-female dyads with 5+ interactions
obs$matrix %>% summarize(joins = sum(joins), leaves = sum(leaves)) # number of joins and leaves between male-female dyads with 5+ interactions
sex %>% count(sex) # counts of sexes of animals from dyads with 5+ interactions
events %>% count(event_type, sex_actor, sex_recipient) %>% filter(sex_actor != sex_recipient) # number of joins/leaves between males and females
# % of dyads in which male vs female is responsible for maintaining proximity
obs$matrix %>%
  mutate(responsible = ifelse(hinde > 0, "male",
                              ifelse(hinde == 0, "mutual", "female"))) %>%
  count(responsible) %>%
  mutate(prop = n/sum(n))

## Plot expected vs. observed.
ggplot() +
  geom_density(data = data.frame(x = c(obs$mean, unlist(lapply(exp, "[[", "mean"))),
                                 class = c("Observed", rep("Expected", num_sims))),
               aes(x = x, fill = class)) +
  scale_fill_manual(values = c("gray","black")) +
  geom_vline(xintercept = obs$mean, size = 0.7) +
  geom_vline(xintercept = obs$mean + obs$se, size = 0.7, linetype = "dashed") +
  geom_vline(xintercept = obs$mean - obs$se, size = 0.7, linetype = "dashed") +
  theme_classic() +
  xlab("Hinde's Index") +
  ylab("Density") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-0.25,0.25), breaks = c(-0.2,-0.1,0,0.1,0.2)) +
  theme(legend.position=c(0.18,0.9)) +
  guides(fill=guide_legend(title=NULL))
ggsave("JuvieHindes.tif", plot = last_plot(), device = "tiff", width = 85, height = 85, units = "mm", dpi = 300)

