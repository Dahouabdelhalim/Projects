# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Simonet & McNally 2020                     #
#                     Relative abundance                        #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# The species relative abundance is included in the phylogenetic mixed model
# We extract these directly from MIDAS merge output
#local_project_dir='/path/to/where/repo/is/cloned'


setwd(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/'))


# relative abundance table from midas merge output
ra.from.merge<- read.table('./output/midas/merge/species_merged/relative_abundance.txt', header=TRUE, sep = '\\t')

# use table of relatedness to get the 101 species
rel<- read.table('./output/tables/relatedness.txt', header=TRUE, sep = '\\t')
sps<- unique(rel[,c('species_id', 'mean_relatedness')])

# reformat
ra.from.merge.long<- gather(ra.from.merge, 'host', 'within_host_relative_abundance', 2:241)

# remove hosts that got filtered out by nb_host > 1 when computing relatedness
ra.from.merge.long<- ra.from.merge.long[ra.from.merge.long$host %in% rel$host,]

# keep only our 101 focus species
rel$species.host<- paste0(rel$species_id, '.', rel$host)
ra.from.merge.long$species.host<- paste0(ra.from.merge.long$species_id, '.', ra.from.merge.long$host)

# intersect
rel2<- select(rel, species_id, host, species.host)
ra.trim<- left_join(rel2, ra.from.merge.long, by = 'species.host')
ra.trim<- select(ra.trim, species.host, within_host_relative_abundance)

write.table(ra.trim, file = paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/relative_abundance.txt'), sep = '\\t', col.names = TRUE, row.names = FALSE)

