# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Simonet & McNally 2020           #
#           SI file figures                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#local_project_dir='/path/to/where/repo/is/cloned'
setwd(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/'))
source('./scripts/analysis/sourced_ggthemes.R')
source('./scripts/analysis/sourced_packages.R')
library(kableExtra)



# GO TRAITS CONTRIBUTION HEAMAPS ----

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

social_go_list<- as.data.frame(read_excel(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/social_go_list_final.xls')))

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




pdf("./output/figures/SI_FIG_ContribBiofilm.pdf", width = 17.3/2.54, height = 25/2.54) 
print(qt_biofilm$contribution_heatmap)
dev.off()

pdf("./output/figures/SI_FIG_ContribQuorumSensing.pdf", width = 17.3/2.54, height = 25/2.54) 
print(qt_quorum_sensing$contribution_heatmap)
dev.off()

pdf("./output/figures/SI_FIG_ContribABdegradation.pdf", width = 17.3/2.54, height = 25/2.54) 
print(qt_ab_degradation$contribution_heatmap)
dev.off()

pdf("./output/figures/SI_FIG_ContribSiderophores.pdf", width = 17.3/2.54, height = 25/2.54) 
print(qt_siderophores$contribution_heatmap)
dev.off()

pdf("./output/figures/SI_FIG_ContribSecretionSyst.pdf", width = 17.3/2.54, height = 25/2.54) 
print(qt_secretion_system_no4$contribution_heatmap)
dev.off()










# PER GENES ANNOTATION VENN DIAGRAMS ----


per_gene<- read.table('./output/tables/per_gene_annotation.txt', header=TRUE)
test<- per_gene[!is.na(per_gene$secretome),] # remove the four species for which there is no secretome data
test<- test[test$sum_ANY != 0,] # focus on overlap among genes that are annotated

pegs.biofilm<- test[test$biofilm == 1, 'peg']
pegs.ab<- test[test$antibiotic_degradation == 1, 'peg']
pegs.qs<- test[test$quorum_sensing == 1, 'peg']
pegs.sid<- test[test$siderophores == 1, 'peg']
pegs.ss<- test[test$secretion_system == 1, 'peg']
pegs.secretome<- test[test$secretome == 1, 'peg']



# We know from SUM_GO that the GO do not overlap with each other at all
# So figure out the overlaps between secretome and the 5 GO categories
length(pegs.secretome)
length(pegs.biofilm)
length(pegs.ab)
length(pegs.qs)
length(pegs.sid)
length(pegs.ss)

length(intersect(pegs.secretome, pegs.biofilm))
length(intersect(pegs.secretome, pegs.ab))
length(intersect(pegs.secretome, pegs.qs))
length(intersect(pegs.secretome, pegs.sid))
length(intersect(pegs.secretome, pegs.ss))
# only Biofilm and AB overalp with secretome

# Use above to fill in a VennDiagram
# Do it this way and nto with venn.diagram() because with >2 sets, this is the only way to have a scale diagram.

library(eulerr)


# OVERALL DIAGRAM
pdf(file = paste0("./output/figures/SI_Fig_VennsDiagrams.pdf"), width=(11.3/2.54), height=(6.5/2.54))

VennDiag <- euler(c("Secretome" = length(pegs.secretome),
                    "Biofilm" = length(pegs.biofilm),
                    "Antibiotic" = length(pegs.ab),
                    "Quorum-sensing" = length(pegs.qs),
                    "Siderophores" = length(pegs.sid),
                    "Sec.Syst." = length(pegs.ss),
                    "Secretome&Biofilm" = length(intersect(pegs.secretome, pegs.biofilm)),
                    "Secretome&Antibiotic" = length(intersect(pegs.secretome, pegs.ab))))

plot(VennDiag,
     counts = TRUE,
     edges = TRUE,#quantities = TRUE, 
     font=2,
     labels = list(cex = .5),
     quantities = list(cex = .5),
     cex=2,
     alpha=0.5,
     fill=c("darkorchid", "forestgreen", "magenta", "yellow", "darkorange", "navy"))


dev.off()




# PER SPECIES DIAGRAM


test<- per_gene[!is.na(per_gene$secretome),] # remove the four species for which there is no secretome data
test<- test[test$sum_ANY != 0,] # focus on overlap among genes that are annotated
mysp = unique(test$species)

pdf(file = paste0("./output/figures/additional_figs/VennsDiagrams_per_species.pdf"), width=(11.3/2.54), height=(6.5/2.54))
for(s in 1:length(mysp)){
  
  
  print(s) 
  
  
  test.s = test[test$species == mysp[s],]
  
  
  pegs.biofilm<- test.s[test.s$biofilm == 1, 'peg']
  pegs.ab<- test.s[test.s$antibiotic_degradation == 1, 'peg']
  pegs.qs<- test.s[test.s$quorum_sensing == 1, 'peg']
  pegs.sid<- test.s[test.s$siderophores == 1, 'peg']
  pegs.ss<- test.s[test.s$secretion_system == 1, 'peg']
  pegs.secretome<- test.s[test.s$secretome == 1, 'peg']
  
  # we know from the overall dataframe that overlap always happen between secretome/biofilm or secretome/antibiotic, never between seceretome and anything else, or the GO categories between each other
  # so we can simply re-apply this code looping over each species
  
  VennDiag <- euler(c("Secretome" = length(pegs.secretome),
                      "Biofilm" = length(pegs.biofilm),
                      "Antibiotic" = length(pegs.ab),
                      "Quorum-sensing" = length(pegs.qs),
                      "Siderophores" = length(pegs.sid),
                      "Sec.Syst." = length(pegs.ss),
                      "Secretome&Biofilm" = length(intersect(pegs.secretome, pegs.biofilm)),
                      "Secretome&Antibiotic" = length(intersect(pegs.secretome, pegs.ab))))
  
  
  print(plot(VennDiag,
             counts = TRUE,
             edges = TRUE,#quantities = TRUE, 
             font=2,
             labels = list(cex = .5),
             quantities = list(cex = .5),
             cex=1,
             alpha=0.5,
             main = list(cex = 0.7, label = mysp[s], font = 2),
             fill=c("darkorchid", "forestgreen", "magenta", "yellow", "darkorange", "navy")))
  
  
}
dev.off()







# SUPPLEMENTARY TABLES ----

format_effects<- function(string){
  
  string_edited<- string %>%
    gsub(pattern = "(Intercept)", replacement = "Intercept", ., fixed = TRUE) %>%
    gsub(pattern = "mean_relatedness", replacement = "Mean relatedness", ., fixed = TRUE) %>%
    gsub(pattern = "log(nb_cds_not_involved_in_response)", replacement = "Log(genome size)", ., fixed = TRUE) %>%
    gsub(pattern = "gram_profilep", replacement = "Gram positive", ., fixed = TRUE) %>%
    gsub(pattern = "species", replacement = "Phylogenetic", ., fixed = TRUE) %>%
    gsub(pattern = "units", replacement = "Residual (non-phylogenetic)", ., fixed = TRUE) %>%
    gsub(pattern = "mean_relative_abundance", replacement = "Mean relative abundance", ., fixed = TRUE) %>%
    gsub(pattern = "sporulation_score", replacement = "Sporulation score", ., fixed = TRUE)
  return(string_edited)  
}

format_effects.model3<- function(string){
  
  string_edited<- string %>%
    gsub(pattern = "(Intercept)", replacement = "Intercept", ., fixed = TRUE) %>%
    gsub(pattern = "mean_relatedness", replacement = "Mean relatedness", ., fixed = TRUE) %>%
    gsub(pattern = "log(nb_cds_not_involved_in_response)", replacement = "Log(genome size)", ., fixed = TRUE) %>%
    gsub(pattern = "gram_profilep", replacement = "Gram positive", ., fixed = TRUE) %>%
    gsub(pattern = "species.ide", replacement = "Species, non-phylogenetic", ., fixed = TRUE) %>%
    gsub(pattern = "species", replacement = "Species, phylogenetic", ., fixed = TRUE,  ignore.case = FALSE) %>%
    gsub(pattern = "units", replacement = "Residual", ., fixed = TRUE) %>%
    gsub(pattern = "mean_relative_abundance", replacement = "Mean relative abundance", ., fixed = TRUE) %>%
    gsub(pattern = "sporulation_score", replacement = "Sporulation score", ., fixed = TRUE) %>%
    gsub(pattern = "host", replacement = "Host", ., fixed = TRUE) %>%
    gsub(pattern = "within_Host_relative_abundance", replacement = "Within host relative abundance", ., fixed = TRUE) %>%
    gsub(pattern = "biofilm", replacement = "Biofilm", ., fixed = TRUE) %>%
    gsub(pattern = "ab_degradation", replacement = "Antibiotic degradation", ., fixed = TRUE) %>%
    gsub(pattern = "quorum_sensing", replacement = "Quorum sensing", ., fixed = TRUE) %>%
    gsub(pattern = "siderophores", replacement = "Siderophores", ., fixed = TRUE) %>%
    gsub(pattern = "secretion_system_no4", replacement = "Secretion systems", ., fixed = TRUE) %>%
    gsub(pattern = "nb_extracellular", replacement = "Secretome", ., fixed = TRUE)
  
    
    return(string_edited)  
}

format_full_summary<- function(model, trait){
  
  fixed<- summary(model)$solutions %>%
    as.data.frame()%>%
    rownames_to_column('Effect') %>%
    mutate(Structure = 'Fixed effect',
           pMCMC = ifelse(pMCMC<0.01,
                          formatC(pMCMC, digit = 2, format = 'e'),
                          formatC(pMCMC, digit = 3, format = 'f')))
  
  random<- summary(model)$Gcovariances %>%
    as.data.frame()%>%
    rownames_to_column('Effect') %>%
    mutate(pMCMC = '') %>%
    mutate(Structure = '(Co)variance')
  
  
  unit<- summary(model)$Rcovariances %>%
    as.data.frame()%>%
    rownames_to_column('Effect') %>%
    mutate(pMCMC = '') %>%
    mutate(Structure = '(Co)variance')
  
  
  tab<- rbind(fixed, random, unit) %>%
    mutate(Effect = format_effects(Effect),
           Model = trait) %>%
    mutate(Structure = ifelse(duplicated(Structure) == TRUE, '', Structure),
           Model = ifelse(duplicated(Model) == TRUE, '', Model)) %>%
    select(Model, Structure, Effect, post.mean,`l-95% CI`, `u-95% CI`, eff.samp, pMCMC) %>%
    rename(`Posterior\\n mean` = post.mean,
           `CI95% lower` = `l-95% CI`,
           `CI95% upper` = `u-95% CI`,
           `Effective\\n sampling` = eff.samp) %>%
    as.data.frame() %>%
    mutate(`Effective\\n sampling` = formatC(`Effective\\n sampling`, digit = 0, format = 'f')) %>%
    mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f')))
  
  
  return(tab)
  
}


format_effects_model2<- function(string){
  
  string_edited<- string %>%
    gsub("at.level(trait, 2).species.ide:at.level(trait, 2).species.ide", "Relatedness_species.ide",  ., fixed = TRUE) %>%
    
    gsub("at.level(first, \\"TRUE\\"):at.level(trait, 1).species.ide:at.level(trait, 2).species.ide", "Relatedness,Trait_species.ide",  ., fixed = TRUE) %>%
    gsub("at.level(trait, 2).species.ide:at.level(first, \\"TRUE\\"):at.level(trait, 1).species.ide", "Trait,Relatedness_species.ide",  ., fixed = TRUE) %>%
    gsub("at.level(first, \\"TRUE\\"):at.level(trait, 1).species.ide:at.level(first, \\"TRUE\\"):at.level(trait, 1).species.ide", "Trait_species.ide",  ., fixed = TRUE) %>%
    
    gsub("at.level(trait, 2):at.level(trait, 2).species", "Relatedness_species",  ., fixed = TRUE) %>%
    gsub("at.level(first, \\"TRUE\\"):at.level(trait, 1):at.level(trait, 2).species", "Relatedness,Trait_species",  ., fixed = TRUE) %>%
    gsub("at.level(trait, 2):at.level(first, \\"TRUE\\"):at.level(trait, 1).species", "Trait,Relatedness_species",  ., fixed = TRUE) %>%
    gsub("at.level(first, \\"TRUE\\"):at.level(trait, 1):at.level(first, \\"TRUE\\"):at.level(trait, 1).species", "Trait_species",  ., fixed = TRUE) %>%
    gsub("at.level(trait, 2).units", "Relatedness_units",  ., fixed = TRUE) %>%
    gsub("at.level(first, \\"FALSE\\"):at.level(trait, 1).units", "construct",  ., fixed = TRUE) %>%
    gsub("at.level(trait, 1):at.level(first, \\"TRUE\\"):gram_profilen", "Gram profile (negative)_fixed",  ., fixed = TRUE) %>%
    gsub("at.level(trait, 1):at.level(first, \\"TRUE\\"):log(nb_cds_not_involved_in_response)", "Log(genome size)_fixed",  ., fixed = TRUE) %>%
    gsub("at.level(trait, 2)", "Intercept relatedness_fixed",  ., fixed = TRUE) %>%
    gsub("at.level(trait, 1):at.level(first, \\"TRUE\\")", "Intercept trait_fixed",  ., fixed = TRUE)
  
  return(string_edited)
}
format_full_summary_model2<- function(model, trait){
  
  fixed<- summary(model)$solutions %>%
    as.data.frame()%>%
    rownames_to_column('Effect') %>%
    mutate(Structure = 'Fixed effect',
           pMCMC = ifelse(pMCMC<0.01,
                          formatC(pMCMC, digit = 2, format = 'e'),
                          formatC(pMCMC, digit = 3, format = 'f')))
  
  random<- summary(model)$Gcovariances %>%
    as.data.frame()%>%
    rownames_to_column('Effect') %>%
    mutate(pMCMC = '') %>%
    mutate(Structure = 'G')
  # Add the pMCMC for the (phylogenetic) covariance term
  species.cov.post<- model$VCV[,which(colnames(model$VCV) == 'at.level(first, "TRUE"):at.level(trait, 1):at.level(trait, 2).species')]
  random$pMCMC[which(random$Effect == 'at.level(first, "TRUE"):at.level(trait, 1):at.level(trait, 2).species')]<- (2*sum(species.cov.post<0))/length(species.cov.post)
  
  
  
  unit<- summary(model)$Rcovariances %>%
    as.data.frame()%>%
    rownames_to_column('Effect') %>%
    mutate(pMCMC = '') %>%
    mutate(Structure = 'R')
  # Add the pMCMC for the (mom-phylogenetic) covariance term
  ide.cov.post<- model$VCV[,which(colnames(model$VCV) == 'at.level(first, "TRUE"):at.level(trait, 1).species.ide:at.level(trait, 2).species.ide')]
  unit$pMCMC[which(unit$Effect == 'at.level(first, "TRUE"):at.level(trait, 1).species.ide:at.level(trait, 2).species.ide')]<- (2*sum(ide.cov.post<0))/length(ide.cov.post)
  
  
  tab<- rbind(fixed, random, unit) %>%
    mutate(Effect = format_effects_model2(Effect),
           Model = trait) %>%
    mutate(Structure = do.call('rbind', strsplit(Effect, '_'))[,2]) %>%
    filter(!duplicated(post.mean),
           Structure != 'construct') %>%
    mutate(Structure = gsub('fixed', 'Fixed effects', Structure)) %>%
    mutate(Structure = gsub('species.ide', 'Species (co)-variances', Structure)) %>%
    mutate(Structure = gsub('species', 'Phylogenetic (co)-variances', Structure)) %>%
    mutate(Structure = gsub('units', 'Residual variance', Structure)) %>%
    mutate(Structure = ifelse(duplicated(Structure), '', Structure)) %>%
    mutate(Effect = do.call('rbind', strsplit(Effect, '_'))[,1]) %>%
    select(Model, Structure, Effect, post.mean,`l-95% CI`, `u-95% CI`, eff.samp, pMCMC) %>%
    rename(`Posterior\\n mean` = post.mean,
           `CI95% lower` = `l-95% CI`,
           `CI95% upper` = `u-95% CI`,
           `Effective\\n sampling` = eff.samp) %>%
    as.data.frame() %>%
    mutate(`Effective\\n sampling` = formatC(`Effective\\n sampling`, digit = 0, format = 'f')) %>%
    mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f')))
  
  
  return(tab)
  
}

library(kableExtra)

load("./output/analyses/MODEL1_CHAIN_1.RData")
load("./output/analyses/MODEL2_CHAIN_1.RData")
load("./output/analyses/MODEL3_CHAIN_1.RData")
load("./output/analyses/MODEL4_CHAIN_1.RData")
load("./output/analyses/MODEL5_CHAIN_1.RData")


# MODEL 1

tab<-rbind(
  format_full_summary(mods.R$siderophores, ''),
  format_full_summary(mods.R$biofilm, ''),
  format_full_summary(mods.R$ab_degradation, ''),
  format_full_summary(mods.R$secretome, ''),
  format_full_summary(mods.R$secretion_system_no4, ''),
  format_full_summary(mods.R$quorum_sensing, '')
)



tab.s1<- kable(tab, "latex", booktabs = T, caption = 'Model summaries for the phylogenetic mixed models of cooperation') %>%
  footnote(c("CI95%: 95% credible interval of the posterior distribution",
             "pMCMC: taken as twice the posterior probability that the estimate is negative"),
           fixed_small_size = TRUE, general_title = "") %>%
  kable_styling() %>%
  pack_rows("Siderophores", 1, 5) %>%
  pack_rows("Biofilm", 6, 10) %>%
  pack_rows("Antibiotic degradation", 11, 15) %>%
  pack_rows("Secretome", 16, 21) %>%
  pack_rows("Secretion systems", 22, 26) %>%
  pack_rows("Quorum sensing", 27, 31)


fileConn<-file("./output/figures/TABLE_S1.tex")
writeLines(tab.s1, fileConn)
close(fileConn)


# MODEL 2

tab.model2<-rbind(
  format_full_summary_model2(mods.R.UNCERTAINTY$siderophores, ''),
  format_full_summary_model2(mods.R.UNCERTAINTY$biofilm, ''),
  format_full_summary_model2(mods.R.UNCERTAINTY$ab_degradation, ''),
  format_full_summary_model2(mods.R.UNCERTAINTY$secretome, ''),
  format_full_summary_model2(mods.R.UNCERTAINTY$secretion_system_no4, ''),
  format_full_summary_model2(mods.R.UNCERTAINTY$quorum_sensing, '')
)


tab.model2 <- tab.model2 %>%
  rename(`Post. mean` = `Posterior\\n mean`,
         `Eff. samp.` = `Effective\\n sampling`) %>%
  mutate(Structure = gsub('Phylogenetic', 'Phyl.', Structure, fixed = TRUE))


tab.s2<- kable(tab.model2, "latex", longtable = T, booktabs = T, caption = 'Model summaries for phylogenetic mixed models of cooperation, for the six forms of cooperation, when accounting for uncertainty in relatedness estimates. The total regression coefficient of the response trait over relatedness is the sum of the phylogenetic and non-phylogenetic (residual) covariances divided by the sum of the phylogenetic and non-phylogenetic (residual) variances') %>%
  kable_styling(latex_options = c("HOLD_position", "repeat_header")) %>%
  footnote(c("CI95%: 95% credible interval of the posterior distribution",
             "pMCMC: taken as twice the posterior probability that the estimate is negative"),
           fixed_small_size = TRUE, general_title = "") %>%
  pack_rows("Siderophores", 1, 10) %>%
  pack_rows("Biofilm", 11, 20) %>%
  pack_rows("Antibiotic degradation", 21, 30) %>%
  pack_rows("Secretome", 31, 41) %>%
  pack_rows("Secretion systems", 42, 51) %>%
  pack_rows("Quorum sensing", 52, 61)


fileConn<-file("./output/figures/TABLE_S2.tex")
writeLines(tab.s2, fileConn)
close(fileConn)


# MODEL 3

fixed<- summary(m3)$solutions %>% as.data.frame()%>%
  rownames_to_column('Effect') %>%
  mutate(Structure = 'Fixed effects',
         pMCMC = ifelse(pMCMC<0.01,
                        formatC(pMCMC, digit = 2, format = 'e'),
                        formatC(pMCMC, digit = 3, format = 'f')))

random<- summary(m3)$Gcovariances %>% as.data.frame()%>%
  rownames_to_column('Effect') %>%
  mutate(pMCMC = '') %>%
  mutate(Structure = 'Variances')


unit<- summary(m3)$Rcovariances %>% as.data.frame()%>%
  rownames_to_column('Effect') %>%
  mutate(pMCMC = '') %>%
  mutate(Structure = 'Variances')


tab.model3<- rbind(fixed, random, unit) %>%
  mutate(Effect = format_effects.model3(Effect)) %>%
  mutate(Structure = ifelse(duplicated(Structure) == TRUE, '', Structure)) %>%
  select(Structure, Effect, post.mean,`l-95% CI`, `u-95% CI`, eff.samp, pMCMC) %>%
  rename(`Posterior\\n mean` = post.mean,
         `CI95% lower` = `l-95% CI`,
         `CI95% upper` = `u-95% CI`,
         `Effective\\n sampling` = eff.samp) %>%
  as.data.frame() %>%
  mutate(`Effective\\n sampling` = formatC(`Effective\\n sampling`, digit = 0, format = 'f')) %>%
  as.data.frame() %>%
  mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f')))
  

tab.model3<- tab.model3[c(1:4,6,5,7),]


tab.s3<- kable(tab.model3, "latex", booktabs = T, caption = 'Model summary for the phylogenetic mixed model of relatedness (drivers of relatedness)', row.names = FALSE, linesep = "") %>%
  footnote(c("CI95%: 95% credible interval of the posterior distribution",
             "pMCMC: taken as twice the posterior probability that the estimate is negative"),
           fixed_small_size = TRUE, general_title = "")


fileConn<-file("./output/figures/TABLE_S3.tex")
writeLines(tab.s3, fileConn)
close(fileConn)



# MODEL 4
tab.model4<-rbind(
  format_full_summary(mods.R.RA.SPO$siderophores, ''),
  format_full_summary(mods.R.RA.SPO$biofilm, ''),
  format_full_summary(mods.R.RA.SPO$ab_degradation, ''),
  format_full_summary(mods.R.RA.SPO$secretome, ''),
  format_full_summary(mods.R.RA.SPO$secretion_system_no4, ''),
  format_full_summary(mods.R.RA.SPO$quorum_sensing, '')
)

tab.s4<- kable(tab.model4, "latex", booktabs = T, caption = 'Model summaries for the phylogenetic mixed models of cooperation when including sporulation scores and relative abundance as predictors') %>%
  footnote(c("CI95%: 95% credible interval of the posterior distribution",
             "pMCMC: taken as twice the posterior probability that the estimate is negative"),
           fixed_small_size = TRUE, general_title = "")%>%
  kable_styling() %>%
  pack_rows("Siderophores", 1, 7) %>%
  pack_rows("Biofilm", 8, 14) %>%
  pack_rows("Antibiotic degradation", 15, 21) %>%
  pack_rows("Secretome", 22, 29) %>%
  pack_rows("Secretion systems", 30, 36) %>%
  pack_rows("Quorum sensing", 37, 43)


fileConn<-file("./output/figures/TABLE_S4.tex")
writeLines(tab.s4, fileConn)
close(fileConn)



# MODEL 5

fixed<- summary(m5)$solutions %>% as.data.frame()%>%
  rownames_to_column('Effect') %>%
  mutate(Structure = 'Fixed effects',
         pMCMC = ifelse(pMCMC<0.01,
                        formatC(pMCMC, digit = 2, format = 'e'),
                        formatC(pMCMC, digit = 3, format = 'f')))

random<- summary(m5)$Gcovariances %>% as.data.frame()%>%
  rownames_to_column('Effect') %>%
  mutate(pMCMC = '') %>%
  mutate(Structure = 'Variances')


unit<- summary(m5)$Rcovariances %>% as.data.frame()%>%
  rownames_to_column('Effect') %>%
  mutate(pMCMC = '') %>%
  mutate(Structure = 'Variances')


tab.model5<- rbind(fixed, random, unit) %>%
  mutate(Effect = format_effects.model3(Effect)) %>%
  mutate(Structure = ifelse(duplicated(Structure) == TRUE, '', Structure)) %>%
  select(Structure, Effect, post.mean,`l-95% CI`, `u-95% CI`, eff.samp, pMCMC) %>%
  rename(`Posterior\\n mean` = post.mean,
         `CI95% lower` = `l-95% CI`,
         `CI95% upper` = `u-95% CI`,
         `Effective\\n sampling` = eff.samp) %>%
  as.data.frame() %>%
  mutate(`Effective\\n sampling` = formatC(`Effective\\n sampling`, digit = 0, format = 'f')) %>%
  as.data.frame() %>%
  mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f')))


tab.model5<- tab.model5[c(1:9,10,12,11,13),]


tab.s5<- kable(tab.model5, "latex", booktabs = T, caption = 'Model summary for the phylogenetic mixed model of relatedness (drivers of relatedness) with cooperation (all six forms) included as fixed predictors', row.names = FALSE, linesep = "") %>%
  footnote(c("CI95%: 95% credible interval of the posterior distribution",
             "pMCMC: taken as twice the posterior probability that the estimate is negative"),
           fixed_small_size = TRUE, general_title = "")

fileConn<-file("./output/figures/TABLE_S5.tex")
writeLines(tab.s5, fileConn)
close(fileConn)


# META-ANALYSES

tab.meta<- rbind(MA.MODELS_1, MA.MODELS_2, MA.MODELS_4) %>%
  as.data.frame() %>%
  select(predictor, estimate, se, ci.lower, ci.upper, z.value, p.value) %>%
  rename(Predictor = predictor,
         Estimate = estimate,
         `Std. Err.` = se,
         `CI95% lower` =  ci.lower,
         `CI95% upper` = ci.upper,
         `z value` = z.value,
         `p value` = p.value) %>%
  mutate(Predictor = format_effects(Predictor)) %>%
  mutate(Predictor = gsub('Genome size', 'Log(genome size)', Predictor))



library(kableExtra)
tab.s6<- kable(tab.meta, "latex", booktabs = T, caption = 'Meta-analysis model summaries') %>%
  footnote(c("CI95%: 95% confidence intervals",
                 "pMCMC: taken as twice the posterior probability that the estimate is negative",
                 "Model 1: meta-analysis over the models of cooperation with mean relatedness as predictor",
                 "Model 2: meta-analysis over the models of cooperation accounting for uncertainty in relatedness estimates",
                 "Model 3: meta-analysis over the models of cooperation with mean relatedness and sporulation scores and relative abundance as predictors"), fixed_small_size = TRUE, general_title = "")%>%
  kable_styling() %>%
  pack_rows("Model 1", 1, 2) %>%
  pack_rows("Model 2", 3, 4) %>%
  pack_rows("Model 3", 5, 8)



fileConn<-file("./output/figures/TABLE_S6.tex")
writeLines(tab.s6, fileConn)
close(fileConn)






