# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Simonet & McNally 2020                     #
#                     Statistical analyses                      #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# We run 5 models:
# (1) mean relatedness
# (1.2) meta-analysis on models 1
# --> saved as MODEL1_CHAIN_X.RData

# (2) Within host relatedness (uncertainty)
# (2.2) meta-analysis on models 2
# --> saved as MODEL2_CHAIN_X.RData

# (3) Drivers of relatedness
# --> saved as MODEL3_CHAIN_X.RData

# (4) mean relatedness, with sporulation and relative abundance
# (4.2) meta-analysis on models 4
# --> saved as MODEL4_CHAIN_X.RData

# (5) Drivers of relatedness, with cooperative measures included
# --> saved as MODEL5_CHAIN_X.RData


#local_project_dir='/path/to/where/repo/is/cloned'
setwd(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/'))
source(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/scripts/analysis/sourced_packages.R'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# PREP DATA ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# d: full dataset, to use for models using within host values, i.e. (#2) uncertainty model and (#3) drivers of relatedness
# d_109: same as d but for models with secretome size as response: we remove species for which PSORTb courld not be run because gram profile is not be determined.
# d_mean: relatedness and abundance averaged by species, for model (#1) and (#4)
# d_mean_109: same as d_mean but for models with secretome size as response


d<- read.table(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/ANALYSIS_DATA_ASSEMBLED.txt'), header=TRUE, stringsAsFactors = FALSE) %>% mutate(first = !duplicated(species_id)) %>% rename(species = species_id)


toremove<- unique(d[d$gram_profile == 'gram0','species'])

d_109<- d %>% filter(!is.element(species, toremove)) %>% as.data.frame()
                        

d_mean<- d %>%
  group_by(species) %>%
  mutate(mean_relative_abundance = mean(within_host_relative_abundance)) %>%
  select(-host, -within_host_relatedness, -within_host_relative_abundance) %>%
  subset(first == TRUE) %>%
  ungroup() %>%
  as.data.frame()

d_mean_109<- d_mean %>% filter(!is.element(species, toremove)) %>% as.data.frame()


# PHYLOGENY
library("ape")
library('MCMCglmm')

midas.tree<- read.tree('./data/species_info_files/midas_tree_renamed.newick')

phylogeny<- drop.tip(midas.tree, midas.tree$tip.label[which(!is.element(midas.tree$tip.label, d_mean$species))])
phylogeny<- chronopl(phylogeny, lambda = 0)
phylogeny<-makeNodeLabel(phylogeny)

phylogeny_109<- drop.tip(phylogeny, phylogeny$tip.label[which(phylogeny$tip.label %in% toremove)])
phylogeny<-makeNodeLabel(phylogeny)
phylogeny_109<-makeNodeLabel(phylogeny_109)

Ainv<-inverseA(phylogeny, scale=FALSE)$Ainv
Ainv_109<-inverseA(phylogeny_109, scale=FALSE)$Ainv


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1) MEAN R ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# RESPONSE: number of genes involved in a cooperative trait, either secretome size, biofilm, siderophores, QS, antibiotic degradation, secretion systems
# PREDICTORS:  mean relatedness, log(genome_size), gram_profile. Genome_size is the number of CDS NOT involved in the cooperative trait. Gram_profile is included only when the response is the secretome size (because of PSORTb)


# MODEL CODE WRAPPER
# Model code wrapped in a function to re-run over each cooperative behaviour
# Function allows to run either model (#1) or (#4) by specifying the 'predictors' argument which can take two values: 'r_only' or 'r_ra_spo'.
# If 'r_only' is specified, the model function runs the model (#1) with only mean relatedness/genome_size/gram_profile (I.A.) as predictors
# if 'r_ra_spo' is specified, it runs the model (#4) by adding sporulation score and mean relative abundance in the predictors
# Gram_profile is included only for the model on secretome size. It's the user that specifies in the arguments to include it or not. Set it to TRUE when the response is secretome size. Set it to FALSE otherwise.

# The only thing that change between these different model is the inclusion or not of some main predictors. We use the default priors for the main effects, so the prior to specify is the same regardless of which model is run. Can use the same prior for all models.


pmm.new<- function(response, predictors, dat, phylogeny_inverse_matrix, prior, nitt, thin, burnin, include_gram){
  
  
  # ADD THE GENOME SIZE WHICH MUST BE INCLUDED AS PREDICTOR VARIABLE - Total n.o. CDS - n.o. CDS involved in response
  dat$nb_cds_not_involved_in_response<- dat$total_cds - dat[,which(colnames(dat) == response)]
  names(dat)[names(dat) == response] <- 'response_of_interest'
  
  
  if(include_gram == TRUE){
    if(predictors == 'r_only'){
      mod<- MCMCglmm(response_of_interest ~ mean_relatedness + log(nb_cds_not_involved_in_response) + gram_profile,
                     random = ~species,
                     ginverse = list(species=phylogeny_inverse_matrix),
                     data = dat,
                     prior=prior,
                     family=c("poisson"),
                     start=list(QUASI=FALSE),
                     pl = TRUE, pr = TRUE, nodes = 'ALL', DIC = TRUE, nitt=nitt, thin=thin, burnin=burnin)
      
    }else if(predictors == 'r_ra_spo'){
      
      mod<- MCMCglmm(response_of_interest ~ mean_relatedness + mean_relative_abundance + sporulation_score + log(nb_cds_not_involved_in_response) + gram_profile,
                     random = ~species,
                     ginverse = list(species=phylogeny_inverse_matrix),
                     data = dat,
                     prior=prior,
                     family=c("poisson"),
                     start=list(QUASI=FALSE),
                     pl = TRUE, pr = TRUE, nodes = 'ALL', DIC = TRUE, nitt=nitt, thin=thin, burnin=burnin)
    }else{
      print('predictor argument is invalid')
    }
  }else if(include_gram == FALSE){
    if(predictors == 'r_only'){
      mod<- MCMCglmm(response_of_interest ~ mean_relatedness + log(nb_cds_not_involved_in_response),
                     random = ~species,
                     ginverse = list(species=phylogeny_inverse_matrix),
                     data = dat,
                     prior=prior,
                     family=c("poisson"),
                     start=list(QUASI=FALSE),
                     pl = TRUE, pr = TRUE, nodes = 'ALL', DIC = TRUE, nitt=nitt, thin=thin, burnin=burnin)
      
    }else if(predictors == 'r_ra_spo'){
      mod<- MCMCglmm(response_of_interest ~ mean_relatedness + mean_relative_abundance + sporulation_score + log(nb_cds_not_involved_in_response),
                     random = ~species,
                     ginverse = list(species=phylogeny_inverse_matrix),
                     data = dat,
                     prior=prior,
                     family=c("poisson"),
                     start=list(QUASI=FALSE),
                     pl = TRUE, pr = TRUE, nodes = 'ALL', DIC = TRUE, nitt=nitt, thin=thin, burnin=burnin)
      
    }else{
      print('predictor argument is invalid')
    }
  }
  
  
}


# PRIOR
# default uninformative prior for all main effects (so not specified)

prior.1<- list(R=list(R1 = list(V=diag(1), nu = 0.002)),  # a single residual variance
               G=list(G1 = list(V=diag(1), nu = 0.002)))  # a single random effect


# MODEL RUN SPECS:
nitt<- 1e+06
burnin<- 100000
thin<- 100

# RUNNING MODELS
pmm.new.secretome.R<- pmm.new('nb_extracellular', 'r_only', d_mean_109, Ainv_109, prior.1, nitt, thin, burnin, TRUE)
pmm.new.biofilm.R<- pmm.new('biofilm', 'r_only', d_mean, Ainv, prior.1, nitt, thin, burnin, FALSE)
pmm.new.ab_degradation.R<- pmm.new('ab_degradation', 'r_only', d_mean, Ainv, prior.1, nitt, thin, burnin, FALSE)
pmm.new.quorum_sensing.R<- pmm.new('quorum_sensing', 'r_only', d_mean, Ainv, prior.1, nitt, thin, burnin, FALSE)
pmm.new.siderophores.R<- pmm.new('siderophores', 'r_only', d_mean, Ainv, prior.1, nitt, thin, burnin, FALSE)
pmm.new.secretion_system_no4.R<- pmm.new('secretion_system_no4', 'r_only', d_mean, Ainv, prior.1, nitt, thin, burnin, FALSE)

# SAVE MODEL OUTPUTS

mods.R<- list(ab_degradation = pmm.new.ab_degradation.R, 
              siderophores = pmm.new.siderophores.R,
              secretion_system_no4 = pmm.new.secretion_system_no4.R,
              quorum_sensing = pmm.new.quorum_sensing.R,
              biofilm = pmm.new.biofilm.R,
              secretome = pmm.new.secretome.R)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# META-ANALYSIS ON #1 ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#load('output/MODELS_1_chain1.RData')

# small code wrapper to exact mean, sd and HPD intervals out of the Sol posterior
extract_mev.2<-function(model, cooperative_trait, predictor_variable){
  post = model$Sol[,predictor_variable]
  effect = mean(post)
  sd_post = sd(post)
  hpd_lower = HPDinterval(post)[,1]
  hpd_higher = HPDinterval(post)[,2]
  return(data.frame(cooperative_trait = cooperative_trait, effect = effect, sd_post = sd_post, hpd_lower = hpd_lower, hpd_higher = hpd_higher, predictor_variable = predictor_variable))
}


# Extract those effects and sd from MCMC posteriors
dat.R<- rbind(extract_mev.2(mods.R$secretome, 'secretome', 'mean_relatedness'),
              extract_mev.2(mods.R$secretion_system_no4, 'secretion_system', 'mean_relatedness'),
              extract_mev.2(mods.R$biofilm, 'biofilm', 'mean_relatedness'),
              extract_mev.2(mods.R$quorum_sensing, 'quorum_sensing', 'mean_relatedness'),
              extract_mev.2(mods.R$siderophores, 'siderophores', 'mean_relatedness'),
              extract_mev.2(mods.R$ab_degradation, 'ab_degradation', 'mean_relatedness'),
              
              extract_mev.2(mods.R$secretome, 'secretome', 'log(nb_cds_not_involved_in_response)'),
              extract_mev.2(mods.R$secretion_system_no4, 'secretion_system', 'log(nb_cds_not_involved_in_response)'),
              extract_mev.2(mods.R$biofilm, 'biofilm', 'log(nb_cds_not_involved_in_response)'),
              extract_mev.2(mods.R$quorum_sensing, 'quorum_sensing', 'log(nb_cds_not_involved_in_response)'),
              extract_mev.2(mods.R$siderophores, 'siderophores', 'log(nb_cds_not_involved_in_response)'),
              extract_mev.2(mods.R$ab_degradation, 'ab_degradation', 'log(nb_cds_not_involved_in_response)')
)



# Small code wrapepr to run meta-analysis on a given predictor variable
# To be run on a dataframe that extracted the effect and sd of the posterior of the models, with a column called 'predictor_variable' that specifies which effect is reported
run.meta.analysis<- function(models.effects.df, predictor_variable){
  library(metafor)
  library(formattable)
  foo<- models.effects.df[models.effects.df$predictor_variable == predictor_variable,]
  randomMA<- rma(yi = foo$effect, sei = foo$sd_post)
  
  
  myp<- ifelse(randomMA$pval< 0.001, as.character(formattable(randomMA$pval, digits = 3, format = "e")), as.character(round(randomMA$pval, 3)))
  
  out<- c(predictor = predictor_variable,
          estimate = round(randomMA$b, 3),
          se = round(randomMA$se,3),
          z.value = round(randomMA$zval,3),
          p.value = myp,
          ci.lower = round(randomMA$ci.lb,3),
          ci.upper = round( randomMA$ci.ub,3))
  
  return(out)
}

MA.MODELS_1<- rbind(
  run.meta.analysis(dat.R, 'mean_relatedness'),
  run.meta.analysis(dat.R, 'log(nb_cds_not_involved_in_response)')
)
MA.MODELS_1


rm(pmm.new.secretome.R)
rm(pmm.new.secretion_system_no4.R)
rm(pmm.new.biofilm.R)
rm(pmm.new.quorum_sensing.R)
rm(pmm.new.siderophores.R)
rm(pmm.new.ab_degradation.R)


save.image('./output/analyses/MODEL1_CHAIN_2.RData')
print('Section 1 done!')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 2) WITHIN HOST R (UNCERTAINTY MODEL) ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# This asks exactly the same question as model (#1) but accounts for the uncertainty in the estimates of relatedness. 
# EXPLAINED IN SHORT: We use a bivariate model formulation to estimate (i) the covariance matrix between the species level variance for relatedness and the residual variance for secretome size and (ii) the covariance matrix between the inter-species phylogenetic variance for relatedness and the response cooperative trait. This gives two regression coefficient of the response cooperative trait and relatedness: the covariance from (i) divided by the residual variance in the response cooperative trait gives the non-phylogenetic regression coefficient and the covariance divided by the phylogenetic variance on the response cooperative trait gives the phylogenetic regression coefficient. The total regression coefficient is [cov(i) + cov(ii)]/[phylogenetic variance in Y + residual variance in Y]
# RESPONSE: number of genes involved in a cooperative trait, either secretome size, biofilm, siderophores, QS, antibiotic degradation, secretion systems
# PREDICTORS: within host relatedness, log(genome_size), gram_profile I.A.


# MODEL FUNCTION WRAPPED
# As aboved, wrapped the model code in a function to re-run for each cooperative trait
# This one has always the same predictors because we do not run it w/wo the ecological traits
# It is run w/wo gram_profile though so it does have the 'include_gram' argument to be specified as previously

hrpmm<- function(response, dataset, phylogeny_inverse_matrix, prior, nitt, thin, burnin, include_gram){
  
  # Set the genome size (total nb of CDS - number of CDS involved in the focal behaviour (response))
  dataset$nb_cds_not_involved_in_response<- dataset[,'total_cds'] - dataset[,which(colnames(dataset) == response)]
  
  # Some house-keeping to be able to run this bivariate model formulation.
  # We add a fake variable called "first" wich allow us to select the correct rows in the dataframe depending on which matrix is being estimated
  # it will estimate a residual variance for this fake variable, but it can be sinply disregarded
  dataset$first<-as.factor(!duplicated(dataset$species))
  dataset[which(dataset$first == 'FALSE'),which(colnames(dataset) == response)]<- NA
  
  # Add species.ide for non phylogenetic effects variable
  dataset$species.ide<-dataset$species
  
  # Rename the column for the focus response
  names(dataset)[names(dataset) == response] <- 'response_of_interest' # rename the column so that this one is used in the model. The column renamed will be the one set by user in the "response" argument
  
  
  if(include_gram == FALSE){
    
    model<-  MCMCglmm(cbind(response_of_interest, within_host_relatedness) ~ -1+
                        at.level(trait,1):at.level(first,"FALSE")+
                        at.level(trait,2)+                                                              # relatedness mean
                        at.level(trait,1):at.level(first,"TRUE")+                                       # behaviour of interest intercept
                        at.level(first,"TRUE"):at.level(trait,1):log(nb_cds_not_involved_in_response),  # effect of genome size
                      
                      random = ~us(at.level(first,"TRUE"):at.level(trait,1)+at.level(trait,2)):species+ # phylogenetic covariance between the two traits
                        idh(at.level(trait,2)):species.ide, # non-phylogenetic variance in relatedness
                      
                      
                      # rcov used to estimate a covariance between a phylogenetic effect and a residual effect
                      # This is not normally done
                      # here it is the case because there is a non-phylogenetic variance + a residual variance for relatedness
                      # While non-phylogenetic and residual variance are the same thing for the response trait
                      # the non-phylogenetic variance in relatedness is estimated above in the idh() formula
                      # the non-phylo/residual variance in response + residual variance in relatedness are estimated is the (co)-variance matrix below (rcov)
                      
                      rcov = ~idh(at.level(first,"TRUE"):at.level(trait,1)):species.ide+ # non-phylogenetic/residual variance in response
                        idh(at.level(trait,2)):units+                                     # residual variance in relatedness
                        idh(at.level(first,"FALSE"):at.level(trait,1)):units,             # residual variance of our fake variable
                      
                      
                      ginverse = list(species=phylogeny_inverse_matrix),
                      data = dataset,
                      prior = prior,
                      start=list(QUASI=FALSE),
                      family=c("poisson", "gaussian"), # poisson for the response trait (a number of genes), gaussiand for the relatedness
                      pl = TRUE, pr = TRUE, DIC = TRUE, nodes = 'ALL',
                      nitt=nitt, thin=thin, burnin=burnin)
    
  }else if(include_gram == TRUE){
    
    # The only thing that changes is the additional main predictor: at.level(first,"TRUE"):at.level(trait,1):gram_profile
    # Everything else (the variance/covariance structure) is exactly the same as above
    
    model<-  MCMCglmm(cbind(response_of_interest, within_host_relatedness) ~ -1+
                        at.level(trait,1):at.level(first,"FALSE")+
                        at.level(trait,2)+                                                      
                        at.level(trait,1):at.level(first,"TRUE")+                                
                        at.level(first,"TRUE"):at.level(trait,1):gram_profile +                  # effect of being gram negative on behaviour of interest
                        at.level(first,"TRUE"):at.level(trait,1):log(nb_cds_not_involved_in_response),
                      
                      random = ~us(at.level(first,"TRUE"):at.level(trait,1)+at.level(trait,2)):species+ 
                        idh(at.level(trait,2)):species.ide,
                      
                      rcov = ~idh(at.level(first,"TRUE"):at.level(trait,1)):species.ide+ 
                        idh(at.level(trait,2)):units+                                   
                        idh(at.level(first,"FALSE"):at.level(trait,1)):units,          
                      
                      
                      ginverse = list(species=phylogeny_inverse_matrix),
                      data = dataset,
                      prior = prior,
                      start=list(QUASI=FALSE),
                      family=c("poisson", "gaussian"),
                      pl = TRUE, pr = TRUE, DIC = TRUE, nodes = 'ALL',
                      nitt=nitt, thin=thin, burnin=burnin)
    
    
    
    
    
    
  }else{
    print('Include_gram argument badly set. Must be either TRUE (if model runs on secretome size) or FALSE (for other models)')
  }
  
  return(model)
}


# PRIOR
prior.2<-list(R=list(R1=list(V=diag(2), nu=0.002, covu=TRUE), # non-phylogenetic covariance, i.e. non-phylogenetic variance in relatedness and non-phylogenetic/residual variance in response trait. The covu is for the model to estimate the covariance between a random effect and a residual effect
                     R2=list(V=1, nu=0.002), # residual variance in relatedness
                     R3=list(V=1, fix=1)),  # residual variance in the fake variable
              G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*1000))) #  Phylogenetic covariance

# MODEL SPECS
nitt<- 13000*50
burnin<- 3000*50
thin<- 10*5


# RUNNING MODELS
pmm_secretome <- hrpmm('nb_extracellular', d_109, Ainv_109, prior.2, nitt, thin, burnin, include_gram = TRUE)
pmm_secretion_system <- hrpmm('secretion_system_no4', d, Ainv, prior.2, nitt, thin, burnin, include_gram = FALSE)
pmm_biofilm<- hrpmm('biofilm', d, Ainv, prior.2, nitt, thin, burnin, include_gram = FALSE)
pmm_QS <- hrpmm('quorum_sensing', d, Ainv, prior.2, nitt, thin, burnin, include_gram = FALSE)
pmm_siderophores <- hrpmm('siderophores', d, Ainv, prior.2, nitt, thin, burnin, include_gram = FALSE)
pmm_ab_degradation <- hrpmm('ab_degradation', d, Ainv, prior.2, nitt, thin, burnin, include_gram = FALSE)

# SAVE MODEL OUTPUTS
mods.R.UNCERTAINTY<- list(ab_degradation = pmm_ab_degradation,
                          siderophores = pmm_siderophores,
                          secretion_system_no4 = pmm_secretion_system,
                          quorum_sensing = pmm_QS,
                          biofilm = pmm_biofilm,
                          secretome = pmm_secretome)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# META-ANALYSIS ON #2 ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
extract_mev.genomeSize<-function(model, cooperative_trait){
  post = model$Sol[,grep('nb_cds_not_involved_in_response', colnames(model$Sol))]
  effect = mean(post)
  sd_post = sd(post)
  hpd_lower = HPDinterval(post)[,1]
  hpd_higher = HPDinterval(post)[,2]
  
  return(data.frame(cooperative_trait = cooperative_trait, effect = effect, sd_post = sd_post, hpd_lower = hpd_lower, hpd_higher = hpd_higher, predictor_variable = 'Genome size'))
}
extract_mev.relatedness<- function(model, cooperative_trait){
  
  model_output = model
  effect = mean(mcmc(rowSums(model_output$VCV[,c(2,6)])/rowSums(model_output$VCV[,c(4,8)])))
  sd_post = sd(mcmc(rowSums(model_output$VCV[,c(2,6)])/rowSums(model_output$VCV[,c(4,8)])))
  hpd_lower = HPDinterval(mcmc(rowSums(model_output$VCV[,c(2,6)])/rowSums(model_output$VCV[,c(4,8)])))[,1]
  hpd_higher = HPDinterval(mcmc(rowSums(model_output$VCV[,c(2,6)])/rowSums(model_output$VCV[,c(4,8)])))[,2]
  
  return(data.frame(cooperative_trait = cooperative_trait, effect = effect, sd_post = sd_post, hpd_lower = hpd_lower, hpd_higher = hpd_higher, predictor_variable = 'Within host relatedness'))
  
}

dat.R.UNCERTAINTY<- rbind(extract_mev.relatedness(mods.R.UNCERTAINTY$secretome, 'secretome'),
                          extract_mev.relatedness(mods.R.UNCERTAINTY$secretion_system_no4, 'secretion_system_no4'),
                          extract_mev.relatedness(mods.R.UNCERTAINTY$biofilm, 'biofilm'),
                          extract_mev.relatedness(mods.R.UNCERTAINTY$quorum_sensing, 'quorum_sensing'),
                          extract_mev.relatedness(mods.R.UNCERTAINTY$siderophores, 'siderophores'),
                          extract_mev.relatedness(mods.R.UNCERTAINTY$ab_degradation, 'ab_degradation'),
                          
                          extract_mev.genomeSize(mods.R.UNCERTAINTY$secretome, 'secretome'),
                          extract_mev.genomeSize(mods.R.UNCERTAINTY$secretion_system_no4, 'secretion_system_no4'),
                          extract_mev.genomeSize(mods.R.UNCERTAINTY$biofilm, 'biofilm'),
                          extract_mev.genomeSize(mods.R.UNCERTAINTY$quorum_sensing, 'quorum_sensing'),
                          extract_mev.genomeSize(mods.R.UNCERTAINTY$siderophores, 'siderophores'),
                          extract_mev.genomeSize(mods.R.UNCERTAINTY$ab_degradation, 'ab_degradation'))




MA.MODELS_2<- rbind(
  run.meta.analysis(dat.R.UNCERTAINTY, 'Within host relatedness'),
  run.meta.analysis(dat.R.UNCERTAINTY, 'Genome size')
)

MA.MODELS_2


rm(pmm_secretome)
rm(pmm_secretion_system)
rm(pmm_biofilm)
rm(pmm_QS)
rm(pmm_siderophores)
rm(pmm_ab_degradation)


save.image('./output/analyses/MODEL2_CHAIN_2.RData')
print('section 2 done!')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 3) DRIVERS OF RELATEDNESS ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Estimate:
# residual variance (R1)
# Species non-phylogenetics effects (species.ide) (G1)
# host effects variances  (G2)
# Species phylogenetics effects variance (G3)


nitt<-1005000
burnin<-5000
thin<-50

d$species.ide<-d$species

dm3<- d[,c('species', 'host', 'within_host_relatedness', 'species.ide', 'sporulation_score', 'within_host_relative_abundance')]


# We use a parameter expanded prior
prior.0.expanded<- list(R=list(V = diag(1), nu = 0.002),
                        G=list(
                          G1=list(V= 1, nu= 1, alpha.mu= 0, alpha.V= 1000),
                          G2=list(V= 1, nu= 1, alpha.mu= 0, alpha.V= 1000),
                          G3=list(V= 1, nu= 1, alpha.mu= 0, alpha.V= 1000)
                        ))


# We can use the same prior, nitt, burnin and thin. Simply has two main effect added in (for which default prior is used)
m3<- MCMCglmm(within_host_relatedness ~ 1 + sporulation_score + within_host_relative_abundance,
              random = ~species.ide+host+species,
              ginverse = list(species=Ainv),
              data = dm3, prior=prior.0.expanded,
              family=c("gaussian"),
              #start=list(QUASI=FALSE),
              pl = TRUE, nodes = 'ALL', nitt=nitt, thin=thin, burnin=burnin)


summary(m3)


save.image('./output/analyses/MODEL3_CHAIN_2.RData')

print('section 3 done!')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 4) MEAN R RA SPO ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# As in (1), uses mean relatedness as predictor BUT HERE ALSO INCLUDES ecological factors as predictors
# RESPONSE: number of genes involved in a cooperative trait, either secretome size, biofilm, siderophores, QS, antibiotic degradation, secretion systems
# PREDICTORS:  mean relatedness, log(genome_size), gram_profile AND sporulation score AND mean_relative_abundance


# PRIOR
# Is prior.1 defined in section (1)

# MODEL RUN SPECS:
nitt<-  1e+06
burnin<-  100000
thin<- 100

# RUNNING MODELS
pmm.new.secretome.R.RA.SPO<- pmm.new('nb_extracellular', 'r_ra_spo', d_mean_109, Ainv_109, prior.1, nitt, thin, burnin, TRUE)
pmm.new.biofilm.R.RA.SPO<- pmm.new('biofilm', 'r_ra_spo', d_mean, Ainv, prior.1, nitt, thin, burnin, FALSE)
pmm.new.ab_degradation.R.RA.SPO<- pmm.new('ab_degradation', 'r_ra_spo', d_mean, Ainv, prior.1, nitt, thin, burnin, FALSE)
pmm.new.quorum_sensing.R.RA.SPO<- pmm.new('quorum_sensing', 'r_ra_spo', d_mean, Ainv, prior.1, nitt, thin, burnin, FALSE)
pmm.new.siderophores.R.RA.SPO<- pmm.new('siderophores', 'r_ra_spo', d_mean, Ainv, prior.1, nitt, thin, burnin, FALSE)
pmm.new.secretion_system_no4.R.RA.SPO<- pmm.new('secretion_system_no4', 'r_ra_spo', d_mean, Ainv, prior.1, nitt, thin, burnin, FALSE)


# SAVE MODEL OUTPUTS
mods.R.RA.SPO<- list(ab_degradation = pmm.new.ab_degradation.R.RA.SPO,
                     siderophores = pmm.new.siderophores.R.RA.SPO,
                     secretion_system_no4 = pmm.new.secretion_system_no4.R.RA.SPO,
                     quorum_sensing = pmm.new.quorum_sensing.R.RA.SPO,
                     biofilm = pmm.new.biofilm.R.RA.SPO,
                     secretome = pmm.new.secretome.R.RA.SPO)



summary(pmm.new.secretome.R.RA.SPO)
summary(pmm.new.biofilm.R.RA.SPO)
summary(pmm.new.siderophores.R.RA.SPO)
summary(pmm.new.ab_degradation.R.RA.SPO)
summary(pmm.new.quorum_sensing.R.RA.SPO)
summary(pmm.new.secretion_system_no4.R.RA.SPO)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# META-ANALYSIS ON #4 ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Extract effects and sd from posterior of the models
dat.R.RA.SPO<- rbind(extract_mev.2(mods.R.RA.SPO$secretome, 'secretome', 'mean_relatedness'),
                     extract_mev.2(mods.R.RA.SPO$secretion_system_no4, 'secretion_system', 'mean_relatedness'),
                     extract_mev.2(mods.R.RA.SPO$biofilm, 'biofilm', 'mean_relatedness'),
                     extract_mev.2(mods.R.RA.SPO$quorum_sensing, 'quorum_sensing', 'mean_relatedness'),
                     extract_mev.2(mods.R.RA.SPO$siderophores, 'siderophores', 'mean_relatedness'),
                     extract_mev.2(mods.R.RA.SPO$ab_degradation, 'ab_degradation', 'mean_relatedness'),
                     
                     extract_mev.2(mods.R.RA.SPO$secretome, 'secretome', 'mean_relative_abundance'),
                     extract_mev.2(mods.R.RA.SPO$secretion_system_no4, 'secretion_system', 'mean_relative_abundance'),
                     extract_mev.2(mods.R.RA.SPO$biofilm, 'biofilm', 'mean_relative_abundance'),
                     extract_mev.2(mods.R.RA.SPO$quorum_sensing, 'quorum_sensing', 'mean_relative_abundance'),
                     extract_mev.2(mods.R.RA.SPO$siderophores, 'siderophores', 'mean_relative_abundance'),
                     extract_mev.2(mods.R.RA.SPO$ab_degradation, 'ab_degradation', 'mean_relative_abundance'),
                     
                     extract_mev.2(mods.R.RA.SPO$secretome, 'secretome', 'sporulation_score'),
                     extract_mev.2(mods.R.RA.SPO$secretion_system_no4, 'secretion_system', 'sporulation_score'),
                     extract_mev.2(mods.R.RA.SPO$biofilm, 'biofilm', 'sporulation_score'),
                     extract_mev.2(mods.R.RA.SPO$quorum_sensing, 'quorum_sensing', 'sporulation_score'),
                     extract_mev.2(mods.R.RA.SPO$siderophores, 'siderophores', 'sporulation_score'),
                     extract_mev.2(mods.R.RA.SPO$ab_degradation, 'ab_degradation', 'sporulation_score'),
                     
                     extract_mev.2(mods.R.RA.SPO$secretome, 'secretome', 'log(nb_cds_not_involved_in_response)'),
                     extract_mev.2(mods.R.RA.SPO$secretion_system_no4, 'secretion_system', 'log(nb_cds_not_involved_in_response)'),
                     extract_mev.2(mods.R.RA.SPO$biofilm, 'biofilm', 'log(nb_cds_not_involved_in_response)'),
                     extract_mev.2(mods.R.RA.SPO$quorum_sensing, 'quorum_sensing', 'log(nb_cds_not_involved_in_response)'),
                     extract_mev.2(mods.R.RA.SPO$siderophores, 'siderophores', 'log(nb_cds_not_involved_in_response)'),
                     extract_mev.2(mods.R.RA.SPO$ab_degradation, 'ab_degradation', 'log(nb_cds_not_involved_in_response)')
                     
)


MA.MODELS_4<- rbind(
  run.meta.analysis(dat.R.RA.SPO, 'mean_relatedness'),
  run.meta.analysis(dat.R.RA.SPO, 'mean_relative_abundance'),
  run.meta.analysis(dat.R.RA.SPO, 'sporulation_score'),
  run.meta.analysis(dat.R.RA.SPO, 'log(nb_cds_not_involved_in_response)')
)


MA.MODELS_4


rm(pmm.new.secretome.R.RA.SPO)
rm(pmm.new.secretion_system_no4.R.RA.SPO)
rm(pmm.new.biofilm.R.RA.SPO)
rm(pmm.new.quorum_sensing.R.RA.SPO)
rm(pmm.new.siderophores.R.RA.SPO)
rm(pmm.new.ab_degradation.R.RA.SPO)


save.image('./output/analyses/MODEL4_CHAIN_2.RData')

print('section 4 done!')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 5) DRIVERS AND COOPERATIVE TRAITS ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Same as model 3 but adding the cooperative traits along with the ecological factors
# We run it on d_109 thought because NA values in secretome size (97 species instead of 101)
# The prior remains the same
# residual variance (R1)
# Species non-phylogenetics effects (species.ide) (G1)
# host effects variances  (G2)
# Species phylogenetics effects variance (G3)

nitt<- 1005000
burnin<-5000
thin<-50

d_109$species.ide<-d_109$species

# We use a parameter expanded prior
prior.0.expanded<- list(R=list(V = diag(1), nu = 0.002),
                        G=list(
                          G1=list(V= 1, nu= 1, alpha.mu= 0, alpha.V= 1000),
                          G2=list(V= 1, nu= 1, alpha.mu= 0, alpha.V= 1000),
                          G3=list(V= 1, nu= 1, alpha.mu= 0, alpha.V= 1000)
                        ))


# We can use the same prior, nitt, burnin and thin. Simply has two main effect added in (for which default prior is used)
m5<- MCMCglmm(within_host_relatedness ~ 1 + sporulation_score + within_host_relative_abundance + biofilm + ab_degradation + quorum_sensing + siderophores + secretion_system_no4 + nb_extracellular,
              random = ~species.ide+host+species,
              ginverse = list(species=Ainv),
              data = d_109, prior=prior.0.expanded,
              family=c("gaussian"),
              #start=list(QUASI=FALSE),
              pl = TRUE, nodes = 'ALL', nitt=nitt, thin=thin, burnin=burnin)


summary(m5)


add.cor.title<- function(y, x, n){
  
  y2 = y[sample(x = seq(1:length(y)), size = n, replace = FALSE)]
  x2 = x[sample(x = seq(1:length(x)), size = n, replace = FALSE)]
  cor = cor.test(y = y2, x = x2)
  
  mtext(side = 3,
        paste0('cor = ', formatC(cor[[4]], digit = 3, format = 'f'), ' p = ', formatC(cor[[3]], digit = 3, format = 'f')), 
        col = ifelse(cor[[3]] < 0.05, 'red', 'black'))
  
  abline(h = 0, v = 0, lty = 2, col = "grey")
}


foo<- m5$Sol %>% as.data.frame()
par(mfrow = c(4,4))
colnames(foo)
plot(foo[,4] ~ foo[,5], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")
plot(foo[,4] ~ foo[,6], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")
plot(foo[,4] ~ foo[,7], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")
plot(foo[,4] ~ foo[,8], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")
plot(foo[,4] ~ foo[,9], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")

plot(foo[,5] ~ foo[,6], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")
plot(foo[,5] ~ foo[,7], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")
plot(foo[,5] ~ foo[,8], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")
plot(foo[,5] ~ foo[,9], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")

plot(foo[,6] ~ foo[,7], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")
plot(foo[,6] ~ foo[,8], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")
plot(foo[,6] ~ foo[,9], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")

plot(foo[,7] ~ foo[,8], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")
plot(foo[,7] ~ foo[,9], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")

plot(foo[,8] ~ foo[,9], pch = 16, col = adjustcolor('black', alpha = .1)); abline(h = 0, v = 0, lty = 2, col = "grey")


m5.joint.test<- aod::wald.test(cov(m5$Sol), colMeans(m5$Sol), Terms=4:9)$result$chi2


save.image('./output/analyses/MODEL5_CHAIN_2.RData')

print("section 5 done :)")
