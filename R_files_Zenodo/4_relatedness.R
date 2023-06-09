# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Simonet & McNally 2020                     #
#            Compute relatedness from diversity estimates       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Set path to directory containing the output of diversityanalysis (./output/midas/diversity)
# We use the within host and between host diversity measures to compute 1 - pi which is within and between hosts genomic similarity, and from this relatedness = (sim_within - sim_between)/(1 - sim_between)

# STEPS ARE:
# 1) Assembling all the tables together
# 2) Compute genomic similarity as  1 - diversity, do so for within and between
# 3) Compute relatedness and assemble table


# ~~~~~~~~~~~~~~~~~~~~~~ #
#  1) Assemble tables    #
# ~~~~~~~~~~~~~~~~~~~~~~ #


# Load within and between diversity tables
#local_project_dir='/path/to/where/repo/is/cloned'
setwd(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/midas/diversity/'))
source(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/scripts/analysis/sourced_packages.R'))

pi.within.files<- list.files(pattern = 'within')
pi.between.files<- list.files(pattern = 'between')
pi.within.ls<- vector('list', length = length(pi.within.files))
pi.between.ls<- vector('list', length = length(pi.between.files))

for(i in 1:length(pi.within.files)){
  pi.within.ls[[i]]<- read.table(pi.within.files[i], header=TRUE, stringsAsFactors = FALSE) %>%
    mutate(species = strsplit(pi.within.files[i], '.', fixed = TRUE)[[1]][1])
  
  pi.between.ls[[i]]<- read.table(pi.between.files[i], header=TRUE, stringsAsFactors = FALSE) %>%
    mutate(species = strsplit(pi.between.files[i], '.', fixed = TRUE)[[1]][1])
}

pi.within<- do.call('rbind', pi.within.ls)
pi.between<- do.call('rbind', pi.between.ls)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  2) Compute genomic similarity    #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Compute within and between genomic similarity
sim.within.MM<- pi.within %>%
  mutate(sim_within = 1-pi_bp) %>%
  select(species, sample_id, sim_within, sites) %>%
  rename(host = sample_id,
         nb_site_within = sites)

sim.between.MM<- pi.between %>%
  mutate(sim_between = 1-pi_bp) %>%
  select(species, sim_between, sites, samples) %>%
  rename(nb_host = samples,
         nb_site_between = sites)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  3) Compute relatedness    #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~ #

sims.HMPnew.MM<- left_join(sim.within.MM, sim.between.MM, by = 'species') %>%
  mutate(within_host_relatedness = (sim_within - sim_between)/(1 - sim_between)) %>% # point estimate of relatedness (i.e. within sample)
  group_by(species) %>%
  mutate(mean_relatedness = mean(within_host_relatedness)) %>% # mean relatedness
  select(species, host, sim_within, nb_site_within, sim_between, nb_site_between, within_host_relatedness, mean_relatedness, nb_host) %>%
  as.data.frame() %>%
  rename(species_id = species) %>%
  filter(
    nb_site_within > 1000, # species for which no core-genome site could be identified
    nb_host > 1) # need at least 2 hosts to estimate a between host diversity!)


length(unique(sims.HMPnew.MM$species_id)) # 141 species with no filtering at all 
length(unique(sims.HMPnew.MM$species_id)) # 140 after filtering for number of core genomic sites [nb_site_within > 1000,]
length(unique(sims.HMPnew.MM$species_id)) # 101 after further filtering for at least two hosts [nb_host > 1]

length(unique(sims.HMPnew.MM$host)) # 239 hosts remain after filtering for [nb_host > 1] (means that that host contained only species that were present only in that host)


write.table(sims.HMPnew.MM, paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/relatedness.txt'), col.names = TRUE, row.names = FALSE, sep = '\\t')


# Keep tables about pattern of diversity to add in supplementary dataset

pi.within.keep<- pi.within %>%
  rename(host = sample_id,
         species_id = species,
         nb_sites = sites) %>%
  mutate(species.host = paste0(species_id, '.', host)) %>%
  select(species_id, host, species.host, nb_sites, depth, snps, pi, snps_kb, pi_bp)


pi.between.keep<- pi.between %>%
  rename(species_id = species,
         nb_host = samples,
         nb_sites = sites) %>%
  select(species_id, nb_host, nb_sites, snps, pi, snps_kb, pi_bp)


colnames(pi.within.keep)[4:9]<- paste0(colnames(pi.within.keep)[4:9], '_within')
colnames(pi.between.keep)[3:7]<- paste0(colnames(pi.between.keep)[3:7], '_between')

pi.keep<- left_join(pi.within.keep, pi.between.keep, by = 'species_id')

kept.species.host<- paste0(sims.HMPnew.MM$species_id, '.', sims.HMPnew.MM$host)
pi.keep.trim<- pi.keep[pi.keep$species.host %in% kept.species.host,]


write.table(pi.keep.trim, paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/diversity_patterns.txt'), col.names = TRUE, row.names = FALSE, sep = '\\t')









