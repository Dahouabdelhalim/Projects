# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Simonet & McNally 2020                     #
#                Measure GO cooperation categories              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Using the bacteria sociality GO terms identified previously, we now browse the GO annotations of our focal 101 species to quantify the number of genes in their genomes involved in cooperation, and further categorise it into the 5 cooperative behaviour categories (biofilm, cooperative antibiotic degradation, quorum-sensing, siderophores, secretion systems)

# STEPS ARE:
# 1) Define a function that will process the GO annotation of a given species, for a given focal cooperation category, from the list  of bacteria social GO, to output a table with the number of genes in that species falling into that focal cooperation cetagory
# 2) Apply function on each species for each cooperation cetagories (101 * 5)
# 3) Assemble & write table
# + added following review: same code to output per gene breakdown. See second section of the script.



#local_project_dir='/path/to/where/repo/is/cloned'


setwd(local_project_dir)
source(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/scripts/analysis/sourced_packages.R'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  1) Define function to quantify GO categories    #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# for a given species
# read GO annotaion (pannzer output), keep best hit per gene
# read the feature table
# assemble the two, verify the peg names match [note, pannzeer output need not contain all the pegs that were originally listed in feature files, if a peg had not GO match]
# intersect with a list of social GO terms, taken from social_go_list table, for focal_behaviour
# collapse BP/CC/MF ontologies (each gene counted a single time)
# sum up the number of genes that were annotated with one of those social GO terms
# compile a table to see which GO terms contributed to the hits (for records)

quant_go_sociality_2<- function(species, focal_behaviour, social_go_list){
  
  social_gos<- social_go_list
  
  # Reading panzzer output (i.e. GO annotation of CDS), keep best hit only, do some housekeeping
  # Each peg can be present in the table up to three times: once for each of the three GO ontologies
  sp<- read.csv(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/pannzer/', species, '.GO.out'), header=TRUE, sep = '\\t', colClasses = c(rep('character', 4), rep('numeric', 3)))
  sp<- sp[sp$ARGOT_rank == 1,]
  sp<- sp[,1:4]
  colnames(sp)<- c('peg', 'Ontology', 'GO_id', 'Description')
  sp$GO_id<- paste0('GO:', sp$GO_id)
  
  
  # open PATRIC features table of that species [theone that came along the fasta file feeding into Pannzer]
  sp_cds<- read.csv(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/data/patric/features/', species, '.features'), header = TRUE, sep = '\\t')
  sp_cds<- sp_cds[sp_cds$feature_type == 'CDS', c('patric_id', 'product', 'go')]
  sp_cds$patric_id<- as.character(sp_cds$patric_id)
  sp_cds$product<- as.character(sp_cds$product)
  sp_cds$go<- as.character(sp_cds$go)
  colnames(sp_cds)<- c('peg', 'product_patric', 'go_patric')
  
  
  # Intersect with panzzer table to get for each peg, the GO assigned by panzzer
  sp$peg<- paste0(do.call('rbind', strsplit(sp$peg, '\\\\|'))[,1], '|', do.call('rbind', strsplit(sp$peg, '\\\\|'))[,2])
  # Check that all peg names in the panzzer output table are in the peg names of the PATRIC table
  checkpoint<- ifelse(length(which(sp$peg %in% sp_cds$peg == FALSE)) == 0, 'ok', 'not_ok')
  print(paste0('checkpoint 1 : ', checkpoint))
  sp_cds_annot<- full_join(sp_cds, sp, 'peg')
  
  
  # record proportion of pegs now annotated
  tmp<- sp_cds_annot[,c('peg', 'GO_id')] # take just pegs and GO annotation by panzzer
  tmp<- tmp[is.na(tmp$GO_id) == FALSE,] # remove all those with annotation
  
  
  # intersect with social GO list = quantify sociality
  sp_cds_annot$is_focal_behaviour <- sp_cds_annot$GO_id %in% social_gos[social_gos$behaviour == focal_behaviour,1]
  
  # Also record for each term the number of time it was hit in that species --> to check the contribution of each term to a behaviour quantification
  social_gos_focus<- social_gos[social_gos$behaviour == focal_behaviour, ]
  go_contribution<- as.data.frame(table(sp_cds_annot$GO_id)) #%>% rename(GO_id = Var1) %>% right_join(social_gos_focus, 'GO_id')
  
  colnames(go_contribution)<- c('GO_id', 'Freq')
  
  go_contribution<- right_join(go_contribution, social_gos_focus, 'GO_id')
  
  
  go_contribution[is.na(go_contribution$Freq) == TRUE,'Freq']<- 0
  names(go_contribution)[2]<- species
  go_contribution<- go_contribution[,c(1,3,5,2)]
  
  
  sp_cds_annot2<- sp_cds_annot[,c('peg', 'Description', 'is_focal_behaviour')]
  sp_cds_annot2$is_annotated<- ifelse(is.na(sp_cds_annot2$Description) == TRUE, 0, 1)
  
  
  # Each peg can be present up to three times, because we retain the top hit GO match for all three ontologies
  # But when a peg is assigned to e.g. biofilm by both its BP and CC for example, we don't count it as twice biofilm
  # basically if either of the ontologies GO of a given peg falls in one social behaviour this peg is counted as being part of that social behaviour
  # The following thus converts those potential 'multiple hits' ACROSS THE 3 ONTOLOGIES into binary 0/1
  
  
  test<- sp_cds_annot2 %>% group_by(peg) %>% summarise(focal_behaviour_counts = sum(is_focal_behaviour),
                                                       annotated_counts = sum(is_annotated))
  
  
  test2<- data.frame(peg = test$peg, ifelse(test[,c('focal_behaviour_counts', 'annotated_counts')] > 0, 1, 0))
  
  
  quant_sociality<- data.frame(
    species = species,
    focal_behaviour = sum(test2$focal_behaviour_counts),
    total_cds = nrow(test2),
    annotated_cds = sum(test2$annotated_counts)
  )
  
  return(list(quant_go = quant_sociality, go_term_contribution = go_contribution))
  
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  2) Apply on each species/cooperation category    #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# LOAD DATA
# Final bacterial social GO list
social_go_list<- as.data.frame(read_excel(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/social_go_list_final.xls')))


# Get list of final species from the relatedness table
dat<- read.csv(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/relatedness.txt'), sep = '\\t') %>%
  select(species_id, mean_relatedness) %>%
  unique()

colnames(dat)<- c('species', 'mean_relatedness')


# Code wrapper to apply the previous function and runs plotting commands on each species/cooperation category
run_trait_quantification<- function(focal_behaviour){
  
  trait_quantification<- vector('list', length = nrow(dat))
  trait_go_terms_contribution<- vector('list', length = nrow(dat))
  for(i in 1:nrow(dat)){
    trait_quantification[[i]]<- quant_go_sociality_2(dat$species[i], focal_behaviour, social_go_list)[[1]]
    trait_go_terms_contribution[[i]]<- quant_go_sociality_2(dat$species[i], focal_behaviour, social_go_list)[[2]]
    print(i)
  }
  
  # Measure of the behaviour
  trait_quantification_df<- do.call('rbind', trait_quantification)
  
  
  # GO terms contribution
  trait_go_terms_contribution_df<- trait_go_terms_contribution[[1]]
  for(i in 2:length(trait_go_terms_contribution)){
    trait_go_terms_contribution_df<- cbind(trait_go_terms_contribution_df, trait_go_terms_contribution[[i]][,4])
    colnames(trait_go_terms_contribution_df)<- c(colnames(trait_go_terms_contribution_df)[c(1:ncol(trait_go_terms_contribution_df)-1)], colnames(trait_go_terms_contribution[[i]])[4])
  }
  trait_go_terms_contribution_df<- trait_go_terms_contribution_df[,c(2, 4:ncol(trait_go_terms_contribution_df))]
  rownames(trait_go_terms_contribution_df)<- trait_go_terms_contribution_df$description
  trait_go_terms_contribution_df<- trait_go_terms_contribution_df[,-1]
  trait_go_terms_contribution_df <- as.data.frame(t(trait_go_terms_contribution_df))
  trait_go_terms_contribution_df$species<- rownames(trait_go_terms_contribution_df)
  trait_go_terms_contribution_df2<- gather(trait_go_terms_contribution_df, 'GO_id', 'hits', 1:(ncol(trait_go_terms_contribution_df)-1))
  
  plot_ids<- read.table("./data/species_info_files/species_plot_names.txt", header=TRUE, sep = '\\t')
  
  trait_go_terms_contribution_df2<- left_join(trait_go_terms_contribution_df2, plot_ids, 'species')
  
  heatmap<- ggplot(trait_go_terms_contribution_df2, aes(x = GO_id, y = plot_names))+
    geom_tile(aes(fill = log(1+hits))) +
    scale_fill_gradient(low = "white", high = "darkred")+xlab('')+ylab('')+
    theme(#legend.position="none",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.7, "cm"),
      legend.key.width = unit(0.4,"cm") ,
      panel.border= element_blank(),
      axis.text.y = element_text(colour="black", size=5),
      axis.text.x = element_text(colour="black", face = "bold", size=5, angle = 45, vjust=1, hjust=1),
      axis.line.y = element_line(color="black", size = 0.3),
      axis.line.x = element_line(color="black", size = 0.3),
      axis.ticks.y = element_line(color="black", size = 0.3),
      axis.ticks.x = element_line(color="black", size = 0.3),
      plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))
  
  
  return(list('contribution_heatmap' = heatmap,
              'trait_quantification' = trait_quantification_df))
}

qt_biofilm<- run_trait_quantification('biofilm')
qt_ab_degradation<- run_trait_quantification('antibiotic_degradation')
qt_quorum_sensing<- run_trait_quantification('quorum_sensing')
qt_siderophores<- run_trait_quantification('siderophores')
qt_secretion_system_no4<- run_trait_quantification('secretion_system')


# (plots of GO contributions done moved to SI.R script)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  3) Assemble & write table    #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Assembling tables for dataset to use in analysis
qt_biofilm$trait_quantification<- qt_biofilm$trait_quantification %>% rename (biofilm = focal_behaviour) 
qt_ab_degradation$trait_quantification<- qt_ab_degradation$trait_quantification %>% rename (ab_degradation = focal_behaviour) 
qt_quorum_sensing$trait_quantification<- qt_quorum_sensing$trait_quantification %>% rename (quorum_sensing = focal_behaviour) 
qt_siderophores$trait_quantification<- qt_siderophores$trait_quantification %>% rename (siderophores = focal_behaviour) 
qt_secretion_system_no4$trait_quantification<- qt_secretion_system_no4$trait_quantification %>% rename (secretion_system_no4 = focal_behaviour) 

traits_all<- cbind(qt_biofilm$trait_quantification[,c('species', "total_cds", "annotated_cds")],
                   qt_biofilm$trait_quantification[,'biofilm'],
                   qt_ab_degradation$trait_quantification[,'ab_degradation'],
                   qt_quorum_sensing$trait_quantification[,'quorum_sensing'],
                   qt_siderophores$trait_quantification[,'siderophores'],
                   qt_secretion_system_no4$trait_quantification[,'secretion_system_no4'])


colnames(traits_all)<- c('species', 'total_cds', 'annotated_cds', 'biofilm', 'ab_degradation', 'quorum_sensing', 'siderophores', 'secretion_system_no4')


write.table(traits_all, paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/go_cooperation_categories.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\\t')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Simonet & McNally 2020                     #
#             Check per-gene break down of GO categories        #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Following reviews, adding the per-gene breakdown of the six measures of cooperation in the supplementary material tables and an associated supplementary figure.

# This simply re-runs some intermediate steps of the secretome and go_categories scripts
# to output a table with the per-gene details rather than collapsing at the species level


#local_project_dir='/path/to/where/repo/is/cloned'
setwd(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/'))
dat<- read.table(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/ANALYSIS_DATA_ASSEMBLED.txt'), header=TRUE)
dat<- unique(dat[,c(1:11, 18, 19)])


social_go_list<- as.data.frame(read_excel('./output/tables/social_go_list_final.xls'))
social_gos<- social_go_list
cooperative_behaviours<- c('biofilm', 'antibiotic_degradation', 'quorum_sensing', 'siderophores', 'secretion_system')


per_peg_annotations<- vector('list')

for(s in 1:length(dat$species_id)){
  
  
  species<- dat$species_id[s]
  pegs_assignation<- vector('list')
  
  for(k in 1:length(cooperative_behaviours)){
    
    focal_behaviour<- cooperative_behaviours[k]
    
    
    sp<- read.csv(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/pannzer/', species, '.GO.out'), header=TRUE, sep = '\\t', colClasses = c(rep('character', 4), rep('numeric', 3)))
    sp<- sp[sp$ARGOT_rank == 1,]
    sp<- sp[,1:4]
    colnames(sp)<- c('peg', 'Ontology', 'GO_id', 'Description')
    sp$GO_id<- paste0('GO:', sp$GO_id)
    
    
    # open PATRIC features table of that species [theone that came along the fasta file feeding into Pannzer]
    sp_cds<- read.csv(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/data/patric/features/', species, '.features'), header = TRUE, sep = '\\t')
    sp_cds<- sp_cds[sp_cds$feature_type == 'CDS', c('patric_id', 'product', 'go')]
    sp_cds$patric_id<- as.character(sp_cds$patric_id)
    sp_cds$product<- as.character(sp_cds$product)
    sp_cds$go<- as.character(sp_cds$go)
    colnames(sp_cds)<- c('peg', 'product_patric', 'go_patric')
    
    
    # Intersect with panzzer table to get for each peg, the GO assigned by panzzer
    sp$peg<- paste0(do.call('rbind', strsplit(sp$peg, '\\\\|'))[,1], '|', do.call('rbind', strsplit(sp$peg, '\\\\|'))[,2])
    
    sp_cds_annot<- full_join(sp_cds, sp, 'peg')
    
    
    sp_cds_annot$is_focal_behaviour <- sp_cds_annot$GO_id %in% social_gos[social_gos$behaviour == focal_behaviour,1]
    sp_cds_annot2<- sp_cds_annot[,c('peg', 'Description', 'is_focal_behaviour')]
    #test<- sp_cds_annot2 %>% group_by(peg) %>% summarise(focal_behaviour_counts = sum(is_focal_behaviour))
    test<- sp_cds_annot2 %>% group_by(peg) %>% mutate(focal_behaviour_counts = sum(is_focal_behaviour)) %>% select(peg, focal_behaviour_counts) %>% unique()
    test2<- data.frame(peg = test$peg, ifelse(test[,c('focal_behaviour_counts')] > 0, 1, 0))
    colnames(test2)<- c('peg', focal_behaviour)
    pegs_assignation[[k]]<- test2
    rm(test2)
    rm(test)
    rm(sp_cds_annot)
    rm(sp_cds_annot2)
    rm(sp)
    rm(sp_cds)
  }
  
  all<- left_join(pegs_assignation[[1]], pegs_assignation[[2]], by = 'peg') %>%
    left_join(pegs_assignation[[3]], 'peg') %>%
    left_join(pegs_assignation[[4]], 'peg') %>%
    left_join(pegs_assignation[[5]], 'peg')
  
  
  # ADD PER-PEG SECRETOME ANNOTATION
  
  no_gram<- c('Clostridiales_bacterium_61057', 'Clostridiales_bacterium_56470', 'Lachnospiraceae_bacterium_56833', 'Lachnospiraceae_bacterium_51870', 'Guyana_massiliensis_60772')
  
  if(species %in% no_gram){
    coop.keep<- data.frame(peg = all$peg,
                           secretome = NA)
  }else{
    
    
    coop<- read.csv(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/psortb/psortb_output/', species, '.psortb.out'), sep='\\t', row.names = NULL)
    colnames(coop)<- c(colnames(coop)[-1], 'foo')
    coop.keep<- coop[,c('SeqID', 'Final_Localization')]
    colnames(coop.keep)<- c('peg', 'secretome')
    
    
    coop.keep$peg<- paste0(do.call('rbind', strsplit(coop.keep$peg, '\\\\|'))[,1], '|', do.call('rbind', strsplit(coop.keep$peg, '\\\\|'))[,2])
    
    coop.keep$peg<- gsub(' ', '', coop.keep$peg)
    
    #head(coop.keep$peg)
    #head(all$peg)
    #sum(!coop.keep$peg %in% all$peg)
    
    
    coop.keep$secretome<- ifelse(coop.keep$secretome == 'Extracellular', 1, 0)
    
  }
  
  all<- left_join(all, coop.keep, by = 'peg')
  
  all$sum_GO<- rowSums(all[,2:6])
  all$sum_ANY<- rowSums(all[,2:7])
  
  
  if(max(range(all$sum_GO)) > 1){
    print(paste0(species, ' has peg(s) assigned to more than one GO'))
  }else{
    print(paste0(species, ' --- GO ok'))
  }
  
  
  if(is.na(max(range(all$sum_ANY))) == TRUE){
    print(paste0(species, ' --- Secretome NA'))
  }else{
    if(max(range(all$sum_ANY)) > 1){
      print(paste0(species, ' has peg(s) assigned to GO and Secretome'))
    }else{
      print(paste0(species, ' --- Secretome ok'))
    }
  }
  
  all<- cbind(species, all)
  
  per_peg_annotations[[s]]<- all
  
}

per_peg_annotations_flat<- do.call('rbind', per_peg_annotations)
length(unique(per_peg_annotations_flat$species))
length(unique(per_peg_annotations_flat$peg))


write.table(per_peg_annotations_flat, paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/per_gene_annotation.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\\t')



