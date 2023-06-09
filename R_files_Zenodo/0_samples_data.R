# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#            Simonet & McNally 2020                   #
#    Extract access links for unique HMP samples      #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# We used healthy subjects stool metagenomes from the Human Metagenome Project (HMP portal, accessed April 2020: under Project > HMP, Body Site > feces, Studies > WGS-PP1, File Type > WGS raw sequences set, File format > FASTQ. From this, download "manifest" file and associated "metadata" file listing the available fastq files.

# This script parses the manifest file to keep the earliest time point available per subject_id
# (+ retains a single links in cases several urls are available for the same sample from the same subject)


#local_project_dir='/path/to/where/repo/is/cloned'
setwd(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/'))


# Manifest file
hmp.manifest<- read.csv('./data/metagenomes/HMP_manifest_68fb5dcb39.tsv', sep = '\\t', colClasses = 'character') %>%
  arrange(sample_id)
# Metadata
hmp.md<- read.csv('./data/metagenomes/HMP_manifest_metadata_1955206147.tsv', sep = '\\t', colClasses = 'character') %>% 
  mutate(visit_number = as.numeric(visit_number)) %>% 
  arrange(sample_id)

sum(hmp.manifest$sample_id != hmp.md$sample_id) # sample ids correspond.

# retain first visit & link
hmp.first<- cbind(hmp.manifest, hmp.md[,-1]) %>%
  arrange(subject_id, visit_number) %>% 
  mutate(first_occurence = !duplicated(subject_id)) %>%
  filter(first_occurence == TRUE) # That's a total of 251 samples

# This retains the earliest visit number per sample id
# in case several links are available for a same metagenomic sample, that retains a single urls for it (the first occurence in the table)
# (sometimes two version of same samples, I imagine it's re-sequencing of same sample coming from same subject ...? Not explicited in HMP portal. Does not really matter though, according to "sample_id" and "subject_id", we have 251 unique sample and 251 unique subject which mean it's sure that each downloaded metagenome come from different individual)

length(unique(hmp.first$sample_id))
length(unique(hmp.first$subject_id))
length(unique(hmp.first$subject_uuid))


# Make internal host name
# = actual fastq file name as host name, taken from the link.
# to ensure tracking of the actual fastq file used for future  reference
hmp.first$host<- gsub('.tar.bz2', '', stri_reverse(gsub('/.*', '', stri_reverse(do.call('rbind', strsplit(hmp.first$urls, ','))[,1]))))

# write table
hmp.first.final<- hmp.first %>%
  select(sample_id, host, subject_id, subject_uuid, urls, file_id, md5, size, sample_body_site, subject_gender, visit_number)

write.table(hmp.first.final, paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/HMP_first_visit_samples.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\\t')


