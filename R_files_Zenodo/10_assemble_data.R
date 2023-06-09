# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Simonet & McNally 2020                     #
#        Assembling data frame to use in statistical analysis   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#local_project_dir='/path/to/where/repo/is/cloned'
setwd(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables'))


# Data to assemble are:
# - relatedness
# - relative abundance [species.host, within_host_relative_abundance]
# - sporulation scores
# - GO categories cooperation
# - secretome
# - gram profiles


relatedness<- read.table('relatedness.txt', header=TRUE, sep = '\\t') %>%
  mutate(species.host = paste0(species_id, '.', host))

abundance<- read.table('relative_abundance.txt', header=TRUE, sep = '\\t')

spo.scores<- read.table('sporulation_scores.txt', header=TRUE, sep = '\\t') %>%
  rename(species_id = species,
         sporulation_score = score_30_raw) %>%
  select(species_id, sporulation_score)

gos<- read.table('go_cooperation_categories.txt', header=TRUE, sep = '\\t', colClasses = c('character', rep('numeric', 7))) %>%
  rename(species_id = species)

ss<- read.table('secretome.txt', header=TRUE, sep = '\\t', colClasses = c('character', 'numeric')) %>%
  select(species_id, nb_extracellular)
# 97 species only because 4 don't have gram assigned so ss cannot be computed

grams<- read.table('../../data/species_info_files/gram_profiles_db.txt', colClasses = 'character', stringsAsFactors = FALSE, col.names = c('species_id', 'gram_profile')) %>%
  filter(species_id %in% spo.scores$species_id)


# ASSEMBLE ALL

sp.focus<- data.frame(species_id = unique(relatedness$species_id)) # my 101 species

dat<- sp.focus %>%
  left_join(gos, by = 'species_id') %>%
  left_join(grams, by = 'species_id') %>%
  left_join(ss, by = 'species_id') %>%
  left_join(spo.scores, by = 'species_id') %>%
  left_join(relatedness, by = 'species_id') %>%
  left_join(abundance, by = 'species.host')


write.table(dat, paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/ANALYSIS_DATA_ASSEMBLED.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\\t')

