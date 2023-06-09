#############################################################
#####               rubias assignments                #######
#####                Elacatinus lori                  #######
#####                   Jan 2021                      #######
####                   CC D'Aloia                     #######
#####   Acknowledgment: This code is based on the     #######
#####   rubias tutorial written by Dr. Eric C.        #######
####    Anderson, and accessed in Nov 2020 at:        #######
# https://cran.r-project.org/web/packages/rubias/vignettes/rubias-overview.html #
#################################################################################

###############################
###  SET WORKING DIRECTORY  ###
###############################

# If you set wd to folder containing both rubias genotype .txt files
# and this script, the code will run as written

########################
###   LOAD PACKAGES  ###
########################
library(rubias)
library(tidyverse)

########################
###   READ IN DATA   ###
########################

# Genotype files for 40 microsatellite loci in HWE and >15 reads per locus

# Adults. Reporting units are as defined in manuscript (3 major genetic clusters)
goby_micro <- read.table("rubias_adult_genotypes.txt", sep = "\\t", header=TRUE,
                          stringsAsFactors = FALSE)

# Settlers for mixture analyses
goby_sett <- read.table("rubias_settler_genotypes.txt", sep = "\\t", header=TRUE,
                        stringsAsFactors = FALSE)

##########################
###   Self-Assignments ###
##########################

### (1) Start with self-assignment among adults

sa_adults1 <- self_assign(reference=goby_micro, gen_start_col = 5)                 

#Summarize output by repunit (scaled_likelihood = posterior probability of assigning fish to the inferred collection
# given equal prior on every collection in reference; here, sum probabilities across reporting units)
sa_to_repu1 <- sa_adults1 %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))


#Return largest scaled_likelihood value
sa_filtered1 <- sa_adults1 %>%
  group_by(indiv) %>%
  filter(scaled_likelihood==max(scaled_likelihood)) %>%
  select(-missing_loci)
 sa_filtered1                        

 ###########################################
 ###   Explore accuracy of marker panels ###
 ###########################################

 ##Assess accuracy based on markers and groupings of pops into reporting units
 # Control simulated mixtures to reflect sampling
  rep_prop1 <- goby_micro  %>%
   group_by( repunit)  %>%
   summarize(repprop =n() ) %>%
   mutate(proportion = repprop/sum(repprop))
 
#Use these proportions to supply Dirichlet random variable
 arep1 <- rep_prop1 %>%
   ungroup() %>%
   mutate(dirichlet = 10 * proportion) %>%
   select(repunit, dirichlet)
 arep1
 
 #Output when return_indiv_posteriors=TRUE
 set.seed(500)
 goby_msat_sim_prop1b <- assess_reference_loo(reference = goby_micro, 
                                              gen_start_col = 5, 
                                              reps = 100,                       #How many simulated mixtures to create
                                              mixsize = 500,                    #Num individuals in each simulated mixture
                                              alpha_repunit = arep1,
                                              return_indiv_posteriors = TRUE)   #Returns posteriors (PofZs for simulated indivs)

#Summarize first tibble output of mixing proportions
prop <- goby_msat_sim_prop1b$mixing_proportions  %>%
  group_by(iter, repunit) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n)) %>%
  mutate(repu_n_prop = repu_n / sum(repu_n))


ggplot(prop, aes(x = true_repprop, y = reprop_posterior_mean, colour = repunit)) +
  geom_point() +
  scale_color_manual(values=c("#F1DC21","#4DAF4A","#377EB8")) +
  geom_abline(intercept = 0, slope = 1) +
  labs(y= "Simulated Mixture Proportion", x = "True Mixture Proportion") +
  facet_wrap(~ repunit) +
  theme_bw()


#Look at distribution of posteriors to correct reporting unit for fish from diff simulated collections
repu_pofzs <- goby_msat_sim_prop1b$indiv_posteriors %>%
  filter(repunit == simulated_repunit) %>%
  group_by(iter, indiv, simulated_collection, repunit) %>%  # first aggregate over reporting units
  summarise(repu_PofZ = sum(PofZ)) %>%
  ungroup() %>%
  arrange(repunit, simulated_collection) %>%
  mutate(simulated_collection = factor(simulated_collection, levels = unique(simulated_collection)))%>%
  mutate(simulated_repunit = factor(repunit, levels = unique(repunit)))

# Get the number of simulated individuals from each collection
num_simmed <-  goby_msat_sim_prop1b$indiv_posteriors %>%
  group_by(iter, indiv) %>%
  slice(1) %>%
  ungroup() %>%
  count(simulated_repunit) 

# Plot simulated individual assignments 
ggplot(repu_pofzs, aes(x = as.factor(simulated_repunit), y = repu_PofZ)) +
  geom_boxplot(aes(colour = repunit)) +
  scale_color_manual(values=c("#F1DC21","#4DAF4A","#377EB8")) +
  geom_text(data = num_simmed, mapping = aes(y = 1.025, label = n), angle = 90, hjust = 0, vjust = 0.5, size = 4) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, vjust = 0.5)) +
  ylim(c(NA, 1.05)) +
  theme_bw()+
  labs(x = "Simulated Reporting Unit", y = "Probability of assignment to reporting unit")

#Show  values below 0.9 for each simulated reporting unit
num_below09 <- goby_msat_sim_prop1b$indiv_posteriors %>%
  filter(repunit == simulated_repunit) %>%
  group_by(iter, indiv, simulated_collection, repunit) %>%  # first aggregate over reporting units
  summarise(repu_PofZ = sum(PofZ))  %>%
  filter(repu_PofZ < 0.9) %>%
  arrange(repunit, repu_PofZ) %>%
  group_by(repunit) %>%
  tally 

##############################
#Mixture analysis
###############################

#PB Method 
mix_est_pb <- infer_mixture(reference = goby_micro, 
                            mixture = goby_sett, 
                            gen_start_col = 5,
                            method = "PB",              # Parametric bootstrapping method
                            reps=20000,                 # Num interations in MCMC
                            burn_in = 1000,             # Burn-in reps to discard
                            pb_iter = 100)              # 100 boostrapped data sets for boostrap correction (default)             

#Mixture among settlers (including bs corrected)
mixprop <- mix_est_pb$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi)) %>%
  left_join(mix_est_pb$bootstrapped_proportions) %>%
  ungroup() %>%
  arrange(desc(repprop))
mixprop

#Individual assignment probabilities to each reporting unit
rep_indiv_ests_pb <- mix_est_pb$indiv_posteriors %>%
  group_by(mixture_collection, indiv, repunit) %>%
  summarise(rep_pofz = sum(PofZ))

#Return top assignment for each indiv
mix_filtered_pb <- rep_indiv_ests_pb  %>%
  group_by(indiv) %>%
  filter(rep_pofz==max(rep_pofz)) #%>%
mix_filtered_pb