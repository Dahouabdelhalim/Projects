#'---
#'author: Allan Baino
#'date: 30th June, 2019
#'title: Research Analysis 
#'---

#'
#'###Background
#'**We used C/N/S isotope ratios in vulture blood, feathers and ungulate muscle 
#'tissue to estimate diet composition and relative constribution of prey items 
#'over time within and across Rukwa, Selous Game Reserve and Serengeti National 
#'Park in Tanzania. 
#'---
#' 

#'###Research objectives
#'**To determine diet and prey source in vultures over time within Selous Game 
#'Reserve (SGR) and Serengeti National Park (SER) and Rukwa Game Reserve (RGR)**
#'**To estimate the relative contribution of prey items in Serengeti vultures 
#'over time**
#'---
#'

#'###Hypotheses
#'**There is a difference in TDFs,TEFs between vulture blood and 
#'feathers**
#'
#'**There is a difference in diet composition (d13C) in vultures across 
#'habitats**
#'
#'**There is no difference in diet composition (d13C) in vultures over 
#'time**
#'
#'**There is a difference in trophic level of prey items (d15N) within vulture 
#'diet across habitats**
#'
#'**There is no difference in trophic level of prey items (d15N) within vulture 
#'diet over time**
#'
#'**There is no difference in prey source (d34S) in vultures 
#'across habitats**
#'
#'**There is no difference in prey source (d34S) for vultures 
#'over time**
#'
#'**The relative contribution of grazing herbivores vulture diet is larger than 
#'that of browsing herbivores over space and time**
#'---
#'  

#'load libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(SIDER))
suppressPackageStartupMessages(library(MixSIAR))
suppressPackageStartupMessages(library(tidyverse))

#'First objective
#'we used SIDER a software developed for use in R (Healy et al, 2018),to 
#'impute Trophic Discrimination/Enrichment Factors (TDFs,TEFs) for d13C,d15N. 
#'The software contains published TDF estimates for bird and mammal phylogenies 
#'and performs a bayesian analysis on our species of interests' ecological 
#'information, to generate posterior probability distributions for our three 
#'biotracer TDFs,TEFs.
#'---
#'

#'import data set to impute TDF values
SIDER.data<-scrumpSider(iso.data = "all")
#'upload distribution of phylogenies 
combined.trees <- scrumpSider(tree = "all")
#'checking our provided unknown obs. phylogeny against two underlying phylogenies
#'in SIDER
#'African white-backed vultures (AWB)    
new.data.af<-recipeSider(species = "Gyps_africanus", 
                         habitat = "terrestrial", 
                         taxonomic.class = "aves", 
                         tissue = "feather", 
                         diet.type = "carnivore", 
                         tree = combined.trees)
#'we then add our checked obs. to the SIDER data set by pruning the corresponding
#'phylogeny to include only species relevant for d13C TDF estimation 
tdf.data.afc<- prepareSider(data.estimate = new.data.af, 
                            data.isotope = SIDER.data, 
                            tree = combined.trees, 
                            isotope = "carbon")
#'set formula for the fixed and random parts of the phylogenetic regression model
#'the response variable is set as d13C which varies as a function of habitat and
#'diet type (fixed effects). Phylogeny, tissue type and species (to account
#'for multiple observations reported in some species) are all set as random terms
formula.c<- delta13C ~ diet.type + habitat 
random.terms<-( ~ animal + species + tissue)
#'set non-informative priors
prior<-list(R = list(V = 1, nu=0.002), 
            G = list(G1=list(V = 1, nu=0.002),
                     G2=list(V = 1, nu=0.002), 
                     G3=list(V = 1, nu=0.002)))             
#'set default parameters for MCMC exploration
nitt<-c(1200000)#n# of iterations, for chains to converge
burnin<-c(200000)#burn-in rate 
thin<-c(500)#sample thinning
parameters <- c(nitt, thin, burnin)#collect vectors into obj
no.chains <- c(2)#set n# of Markovian chains
convergence =  c(1.1)#set convergence
ESS = c(1000)#set effective sample sizes
#'priors and default parameters have been optimised on a basic model that runs
#'on the full SIDER dataset with no species to be imputed so that all chains 
#'converge and have an effective sample size above 1000.
#'---
#'

#'call MCMCglmm run models that will provide posterior distribution of the focal 
#'TDF estimate and standard diagnostics of the MCMC chain convergence and 
#'effective sample size for each model parameter
TDF.est.afc <- imputeSider(mulTree.data = tdf.data.afc, 
                           formula = formula.c, 
                           random.terms = random.terms,
                           prior = prior, 
                           output = "AWB_cf_run",
                           parameters = parameters,
                           chains = no.chains, 
                           convergence =  convergence, 
                           ESS = ESS)
#'All MCMCglmm chains in models converged under 1.1 and effective sample sizes 
#'were greater than 1000
#'---
#'

#'get bayesian credible intervals of the posterior estimate
CI.afc<-hdrcde::hdr(TDF.est.afc$tdf_global,prob = c(50,95,99))
#'CI at 99% were -2.28 to 4.904, CI at 95% were -1.372 to 4.09, CI at 50% were 
#'0.436 to 2.270
#'---
#' 

#'plot the estimated posterior probability
coda::densplot(TDF.est.afc$tdf_global, main = "Estimated Posterior Probabilty")

#'summary statistics for posterior probability estimate
summary(TDF.est.afc$tdf_global)
#'mean TDF for d13C is 1.38‰in feathers for the African white-backed 
#'at 95% confidence
#'---
#' 

#'repeat procedure to estimate posterior probability for d15N for our obs. 
tdf.data.afn<- prepareSider(data.estimate = new.data.af, 
                            data.isotope = SIDER.data, 
                            tree = combined.trees, 
                            isotope = "nitrogen")
#
formula.n<- delta15N ~ diet.type + habitat
random.terms<-( ~ animal + species + tissue)
#
prior<-list(R = list(V = 1, nu=0.002), 
            G = list(G1=list(V = 1, nu=0.002),
                     G2=list(V = 1, nu=0.002), 
                     G3=list(V = 1, nu=0.002)))             
#
nitt<-c(1200000)
burnin<-c(200000) 
thin<-c(500)
parameters <- c(nitt, thin, burnin)
no.chains <- c(2)
convergence =  c(1.1)
ESS = c(1000)
#
TDF.est.afn<- imputeSider(mulTree.data = tdf.data.afn, 
                          formula = formula.n, 
                          random.terms = random.terms,
                          prior = prior, 
                          output = "AWB_nf_run",
                          parameters = parameters,
                          chains = no.chains, 
                          convergence =  convergence, 
                          ESS = ESS)
#'All MCMCglmm chains in models converged under 1.1 and effective sample sizes 
#'were greater than 1000
#'---
#'

#'get bayesian credible intervals of the posterior estimate
CI.afn<-hdrcde::hdr(TDF.est.afn$tdf_global,prob = c(50,95,99))
#'CI at 99% were -0.056 to 6.49, CI at 95% were 0.74 to 5.702, 
#'CI at 50% were 2.39 to 4.05
#'---
#' 

#'plot the estimated posterior prbability
coda::densplot(TDF.est.afn$tdf_global, main = "Estimated Posterior Probabilty")

#'summary statistics for posterior probability estimate
summary(TDF.est.afn$tdf_global)
#'mean TDF for d15N is 3.21‰ in feathers for the African white-backed 
#'at 95% confidence
#'---
#' 

#'repeat the procedures above for the same species with blood d13C,d15N 
new.data.ab<-recipeSider(species = "Gyps_africanus", 
                         habitat = "terrestrial", 
                         taxonomic.class = "aves", 
                         tissue = "blood", 
                         diet.type = "carnivore", 
                         tree = combined.trees)
#
tdf.data.abc<- prepareSider(data.estimate = new.data.ab, 
                            data.isotope = SIDER.data, 
                            tree = combined.trees, 
                            isotope = "carbon")
#
formula.c<- delta13C ~ diet.type + habitat
random.terms<-( ~ animal + species + tissue)
#
prior<-list(R = list(V = 1, nu=0.002), 
            G = list(G1=list(V = 1, nu=0.002),
                     G2=list(V = 1, nu=0.002), 
                     G3=list(V = 1, nu=0.002)))             
#
nitt<-c(1200000)
burnin<-c(200000) 
thin<-c(500)
parameters <- c(nitt, thin, burnin)
no.chains <- c(2)
convergence =  c(1.1)
ESS = c(1000)
#
TDF.est.abc <- imputeSider(mulTree.data = tdf.data.abc, 
                           formula = formula.c, 
                           random.terms = random.terms,
                           prior = prior, 
                           output = "AWB_cb_run",
                           parameters = parameters,
                           chains = no.chains, 
                           convergence =  convergence, 
                           ESS = ESS)
#'All MCMCglmm chains in models converged under 1.1 and effective sample sizes 
#'were greater than 1000
#'---
#'

#
CI.abc<-hdrcde::hdr(TDF.est.abc$tdf_global,prob = c(50,95,99))
#'CI at 99% were -3.12 to 3.85, CI at 95% were -2.27 to 2.90, CI at 50% were 
#'-0.54 to 1.18
#'---
#' 

#'plot of the estimated posterior prbability
coda::densplot(TDF.est.abc$tdf_global, main = "Estimated Posterior Probabilty")
#'

#'summary statistics
summary(TDF.est.abc$tdf_global)
#'mean TDF for d13C is 0.29‰ in blood for the African white-backed 
#'at 95% confidence
#'---
#' 

#
tdf.data.abn<- prepareSider(data.estimate = new.data.ab, 
                            data.isotope = SIDER.data, 
                            tree = combined.trees, 
                            isotope = "nitrogen")
#
formula.n<- delta15N ~ diet.type + habitat
random.terms<-( ~ animal + species + tissue)
#
prior<-list(R = list(V = 1, nu=0.002), 
            G = list(G1=list(V = 1, nu=0.002),
                     G2=list(V = 1, nu=0.002), 
                     G3=list(V = 1, nu=0.002)))             
#
nitt<-c(1200000)
burnin<-c(200000) 
thin<-c(500)
parameters <- c(nitt, thin, burnin)
no.chains <- c(2)
convergence =  c(1.1)
ESS = c(1000)
#
TDF.est.abn<- imputeSider(mulTree.data = tdf.data.abn, 
                          formula = formula.n, 
                          random.terms = random.terms,
                          prior = prior, 
                          output = "AWB_nb_run",
                          parameters = parameters,
                          chains = no.chains, 
                          convergence =  convergence, 
                          ESS = ESS)
#'All MCMCglmm chains in models converged under 1.1 and effective sample sizes 
#'were greater than 1000
#'---
#'

#
CI.abn<-hdrcde::hdr(TDF.est.abn$tdf_global,prob = c(50,95,99))
#'CI at 99% were -1.0 to 5.55, CI at 95% were -0.124 to 4.68, CI at 50% were
#'1.435 to 3.06 
#'---
#' 

#'plot of the estimated posterior prbability
coda::densplot(TDF.est.abn$tdf_global, main = "Estimated Posterior Probabilty")
#'

#'summary statistics
summary(TDF.est.abn$tdf_global)
#'mean TDF for d15N is 2.23‰in blood for the African white-backed 
#'at 95% confidence 

#'repeat all procedures above for carbon and nitrogen TDF imputation in Ruppell 
#'griffon vultures blood and feather
new.data.rf<-recipeSider(species = "Gyps_rueppellii", 
                         habitat = "terrestrial", 
                         taxonomic.class = "aves", 
                         tissue = "feather", 
                         diet.type = "carnivore", 
                         tree = combined.trees)
#
tdf.data.rfc<- prepareSider(data.estimate = new.data.rf, 
                            data.isotope = SIDER.data, 
                            tree = combined.trees, 
                            isotope = "carbon")
#
formula.c<- delta13C ~ diet.type + habitat 
random.terms<-( ~ animal + species + tissue)
#
prior<-list(R = list(V = 1, nu=0.002), 
            G = list(G1=list(V = 1, nu=0.002),
                     G2=list(V = 1, nu=0.002), 
                     G3=list(V = 1, nu=0.002)))             
#
nitt<-c(1200000)
burnin<-c(200000) 
thin<-c(500)
parameters <- c(nitt, thin, burnin)
no.chains <- c(2)
convergence =  c(1.1)
ESS = c(1000)
#
TDF.est.rfc <- imputeSider(mulTree.data = tdf.data.rfc, 
                           formula = formula.c, 
                           random.terms = random.terms,
                           prior = prior, 
                           output = "RPV_cf_run",
                           parameters = parameters,
                           chains = no.chains, 
                           convergence =  convergence, 
                           ESS = ESS)
#'get bayesian credible intervals
CI.rfc<-hdrcde::hdr(TDF.est.rfc$tdf_global,prob = c(50,95,99))
#'CI at 99% were -2.158 to 4.877, CI at 95% were -1.310 to 3.953, CI at 50% were 
#'0.458 to 2.215
#'---
#' 

#'plot of the estimated posterior prbability
coda::densplot(TDF.est.rfc$tdf_global, main = "Estimated Posterior Probabilty")

#'summary statistics
summary(TDF.est.rfc$tdf_global)
#'mean TDF for d13C is 1.23‰ in blood for the Ruppells griffon vultures 
#'at 95% confidence
#'---
#' 

#'TDF for nitrogen in RPV feathers
tdf.data.rfn<- prepareSider(data.estimate = new.data.rf, 
                            data.isotope = SIDER.data, 
                            tree = combined.trees, 
                            isotope = "nitrogen")
#
formula.c<- delta15N ~ diet.type + habitat 
random.terms<-( ~ animal + species + tissue)
#
prior<-list(R = list(V = 1, nu=0.002), 
            G = list(G1=list(V = 1, nu=0.002),
                     G2=list(V = 1, nu=0.002), 
                     G3=list(V = 1, nu=0.002)))             
#
nitt<-c(1200000)
burnin<-c(200000) 
thin<-c(500)
parameters <- c(nitt, thin, burnin)
no.chains <- c(2)
convergence =  c(1.1)
ESS = c(1000)
#
TDF.est.rfn <- imputeSider(mulTree.data = tdf.data.rfn, 
                           formula = formula.c, 
                           random.terms = random.terms,
                           prior = prior, 
                           output = "RPV_nf_run",
                           parameters = parameters,
                           chains = no.chains, 
                           convergence =  convergence, 
                           ESS = ESS)
#'get bayesian credible intervals
CI.rfn<-hdrcde::hdr(TDF.est.rfn$tdf_global,prob = c(50,95,99))
#'CI at 99% were -0.146 to 6.5, CI at 95% were 0.731 to 5.645, CI at 50% were 
#'2.433 to 4.091
#'---
#' 

#'plot of the estimated posterior probability
coda::densplot(TDF.est.rfn$tdf_global, main = "Estimated Posterior Probabilty")

#'summary statistics
summary(TDF.est.rfn$tdf_global)
#'mean TDF for d15N is 3.28‰ in blood for the Ruppells griffon vultures 
#'at 95% confidence
#'---
#' 

#'repeat procedure for RPV blood
new.data.rb<-recipeSider(species = "Gyps_rueppellii", 
                         habitat = "terrestrial", 
                         taxonomic.class = "aves", 
                         tissue = "blood", 
                         diet.type = "carnivore", 
                         tree = combined.trees)
#
tdf.data.rbc<- prepareSider(data.estimate = new.data.rb, 
                            data.isotope = SIDER.data, 
                            tree = combined.trees, 
                            isotope = "carbon")
#
formula.c<- delta13C ~ diet.type + habitat 
random.terms<-( ~ animal + species + tissue)
#
prior<-list(R = list(V = 1, nu=0.002), 
            G = list(G1=list(V = 1, nu=0.002),
                     G2=list(V = 1, nu=0.002), 
                     G3=list(V = 1, nu=0.002)))             
#
nitt<-c(1200000)
burnin<-c(200000) 
thin<-c(500)
parameters <- c(nitt, thin, burnin)
no.chains <- c(2)
convergence =  c(1.1)
ESS = c(1000)
#
TDF.est.rbc <- imputeSider(mulTree.data = tdf.data.rbc, 
                           formula = formula.c, 
                           random.terms = random.terms,
                           prior = prior, 
                           output = "RPV_cb_run",
                           parameters = parameters,
                           chains = no.chains, 
                           convergence =  convergence, 
                           ESS = ESS)
#'get bayesian credible intervals
CI.rbc<-hdrcde::hdr(TDF.est.rbc$tdf_global,prob = c(50,95,99))
#'CI at 99% were -3.120 to 3.684, CI at 95% were -2.308 to 2.88, CI at 50% were 
#'-0.535 to 1.188
#'---
#' 

#'plot of the estimated posterior prbability
coda::densplot(TDF.est.rbc$tdf_global, main = "Estimated Posterior Probabilty")

#'summary statistics
summary(TDF.est.rbc$tdf_global)
#'mean TDF for d13C is 0.30‰ in blood for the Ruppells griffon vultures 
#'at 95% confidence
#'---
#' 

#'TDF for nitrogen in RPV blood
tdf.data.rbn<- prepareSider(data.estimate = new.data.rb, 
                            data.isotope = SIDER.data, 
                            tree = combined.trees, 
                            isotope = "nitrogen")
#
formula.c<- delta15N ~ diet.type + habitat 
random.terms<-( ~ animal + species + tissue)
#
prior<-list(R = list(V = 1, nu=0.002), 
            G = list(G1=list(V = 1, nu=0.002),
                     G2=list(V = 1, nu=0.002), 
                     G3=list(V = 1, nu=0.002)))             
#
nitt<-c(1200000)
burnin<-c(200000) 
thin<-c(500)
parameters <- c(nitt, thin, burnin)
no.chains <- c(2)
convergence =  c(1.1)
ESS = c(1000)
#
TDF.est.rbn <- imputeSider(mulTree.data = tdf.data.rbn, 
                           formula = formula.c, 
                           random.terms = random.terms,
                           prior = prior, 
                           output = "RPV_nb_run",
                           parameters = parameters,
                           chains = no.chains, 
                           convergence =  convergence, 
                           ESS = ESS)
#'get bayesian credible intervals
CI.rbn<-hdrcde::hdr(TDF.est.rbn$tdf_global,prob = c(50,95,99))
#'CI at 99% were -1.098 to 5.648, CI at 95% were -0.164 to 4.76, CI at 50% were 
#'1.52 to 3.17
#'---
#' 

#'plot of the estimated posterior prbability
coda::densplot(TDF.est.rbn$tdf_global, main = "Estimated Posterior Probabilty")

#'summary statistics
summary(TDF.est.rbn$tdf_global)
#'mean TDF for d15N is 2.39‰ in blood for the Ruppells griffon vultures 
#'at 95% confidence
#'---
#' 

#'Second objective
#'we estimated diet and diet source as given our d13C,d15N,d34S compositions in
#'African white-backed (AWB) and Ruppells griffon (RPV) vulture blood and 
#'feathers, we also compared isotopic compositions in prey items and AWB vulture
#'tissues. 
#'---
#'

#'import data
bf.iso_data<-read.csv("bf_iso_data.csv", header = T)
#'add color and symbols to data frame
#'add empty column for colors
bf.iso_data$color<-NA
#'set colors for the different species
bf.iso_data$color[bf.iso_data$species=='AWB']<- "red"
bf.iso_data$color[bf.iso_data$species=='RPV']<- "blue"
#'add empty column for symbols
bf.iso_data$symbol<-NA
#'set symbols for the different locations sampled
bf.iso_data$symbol[bf.iso_data$location=='RGR']<-0
bf.iso_data$symbol[bf.iso_data$location=='SGR']<-1
bf.iso_data$symbol[bf.iso_data$location=='SER']<-2
#'check data
glimpse(bf.iso_data)
#'convert all d13C, d15N, d34S to numeric values
bf.iso_data$d13C<-as.numeric(as.character(bf.iso_data$d13C))
bf.iso_data$d15N<-as.numeric(as.character(bf.iso_data$d15N))
bf.iso_data$d34S<-as.numeric(as.character(bf.iso_data$d34S))
bf.iso_data$weight<-as.numeric(as.character(bf.iso_data$weight))
bf.iso_data$C.N.ratio<-as.numeric(as.character(bf.iso_data$C.N.ratio))
#'delete all NAs in data set
bf.iso_data<-na.omit(bf.iso_data)

#'basic plots to visualize data
#'visualize d13C signals in vulture tissues
ggplot(data = bf.iso_data[bf.iso_data$location != 'RGR',], aes(x=location, 
  y=d13C, fill=t.subset))+
  geom_boxplot()+
  ylab("d13C in vulture tissues (‰)")+
  xlab("Sampled sites")+
  theme_bw()+
  facet_wrap(~species)+
  scale_fill_manual(name='Tissue type', values=alpha(c("lightgreen","red",
                                                       "lightblue")))
#'d13C signals in blood and feathers indicated similarities in diet for different
#'species of vultures sampled within the same region over time, there appeared 
#'to be no differences in diet for vultures sampled across SER and SGR regions 
#'over time however, AWBs in SGR appeared to have a varied diet over longer 
#'periods of time as indicated by d13C in feathers. Diet in vultures sampled 
#'from RGR apeared to be different from those sampled in SER and SGR.
#'---
#'

#'visualize d15N signals in vulture tissues
ggplot(data = bf.iso_data[bf.iso_data$location != 'RGR',], aes(x=location, 
  y=d15N, fill=t.subset))+
  geom_boxplot()+
  ylab("d15N in vulture tissues (‰)")+
  xlab("Sampled sites")+
  theme_bw()+
  facet_wrap(~species)+
  scale_fill_manual(name='Tissue type', values=alpha(c("lightgreen","red",
                                                       "lightblue")))
#'d15N signals in blood and feathers for vultures species sampled within the 
#'same region indicate differences in trophic ecology over time, however, this 
#'did not appear to vary across SER and SGR regions. The trophic ecology for 
#'vultures in RGR appeared to be different from those in SER and SGR
#'---
#'  

#'visualize d34S signals in vulture tissues
ggplot(data = bf.iso_data[bf.iso_data$location != 'RGR',], aes(x=location, 
  y=d34S, fill=t.subset))+
  geom_boxplot()+
  ylab("d34S in vulture tissues (‰)")+
  xlab("Sampled sites")+
  theme_bw()+
  facet_wrap(~species)+
  scale_fill_manual(name='Tissue type', values=alpha(c("lightgreen","red",
                                                       "lightblue")))
#'d34S signals in blood and feathers for vultures sampled within the same region
#'did not indicate differences in diet source over time. Diet source for vultures
#'in SGR appeared to different from sampled in SER and RGR.
#'---
#'

#'We then calculated the absolute difference in d13C,d15N,d34S for feather 
#'subsets (basal & proximal ends) to identify individual shifts in diet and prey
#'source
#'view unique animal.IDs in our data
unique(bf.iso_data$animal.ID)
#'calculate and collect values for abs differences in d13C in a list
list_diff <- c()
for (id in unique(bf.iso_data$animal.ID)){
  #for rows where t.subset basal.barbs and proximal.barbs has a true value
  ifelse(nrow(bf.iso_data[bf.iso_data$animal.ID == id & bf.iso_data$t.subset == 
                            "basal.barbs",]) == 1 &
           nrow(bf.iso_data[bf.iso_data$animal.ID==id & bf.iso_data
                            $t.subset == "proximal.barbs",]) == 1,
         #calculate their difference
         list_diff <- c(list_diff,as.numeric(bf.iso_data[bf.iso_data$animal.ID==
                                                           id & bf.iso_data$
                                                           t.subset == "basal.barbs",
                                                         ]$d13C - bf.iso_data
                                             [bf.iso_data$animal.ID == id & 
                                                 bf.iso_data$t.subset
                                               == "proximal.barbs",]$d13C)),
         list_diff <- c(list_diff,NA))
}
#'d15N
list_diffN <- c()
for (id in unique(bf.iso_data$animal.ID)){
  #for rows where t.subset basal.barbs and proximal.barbs has a true value
  ifelse(nrow(bf.iso_data[bf.iso_data$animal.ID == id & bf.iso_data$t.subset == 
                            "basal.barbs",]) == 1 &
           nrow(bf.iso_data[bf.iso_data$animal.ID==id & bf.iso_data
                            $t.subset == "proximal.barbs",]) == 1,
         #calculate their difference
         list_diffN <- c(list_diffN,as.numeric(bf.iso_data[bf.iso_data$animal.ID==
                                                             id & bf.iso_data$
                                                             t.subset == "basal.barbs",
                                                           ]$d15N - bf.iso_data
                                               [bf.iso_data$animal.ID == id & 
                                                   bf.iso_data$t.subset
                                                 == "proximal.barbs",]$d15N)),
         list_diffN <- c(list_diffN,NA))
}
#'d34S
list_diffS <- c()
for (id in unique(bf.iso_data$animal.ID)){
  #for rows where t.subset basal.barbs and proximal.barbs has a true value
  ifelse(nrow(bf.iso_data[bf.iso_data$animal.ID == id & bf.iso_data$t.subset == 
                            "basal.barbs",]) == 1 &
           nrow(bf.iso_data[bf.iso_data$animal.ID==id & bf.iso_data
                            $t.subset == "proximal.barbs",]) == 1,
         #calculate their difference
         list_diffS <- c(list_diffS,as.numeric(bf.iso_data[bf.iso_data$animal.ID==
                                                             id & bf.iso_data$
                                                             t.subset == "basal.barbs",
                                                           ]$d34S - bf.iso_data
                                               [bf.iso_data$animal.ID == id & 
                                                   bf.iso_data$t.subset
                                                 == "proximal.barbs",]$d34S)),
         list_diffS <- c(list_diffS,NA))
}

#'put lists in data frame
d13C.abs<-as.data.frame(cbind(unique(as.character(bf.iso_data$animal.ID)),
                              as.numeric(list_diff)))
d15N.abs<-as.data.frame(cbind(unique(as.character(bf.iso_data$animal.ID)),
                              as.numeric(list_diffN)))
d34S.abs<-as.data.frame(cbind(unique(as.character(bf.iso_data$animal.ID)),
                              as.numeric(list_diffS)))
#'convert all variables to numeric
d13C.abs$V2<-as.numeric(as.character(d13C.abs$V2))
d15N.abs$V2<-as.numeric(as.character(d15N.abs$V2))
d34S.abs$V2<-as.numeric(as.character(d34S.abs$V2))
#'remove NAs
d13C.abs<-na.omit(d13C.abs)
d15N.abs<-na.omit(d15N.abs)
d34S.abs<-na.omit(d34S.abs)

#'plot absolute differences to look at distributions
#'d13C
hist(d13C.abs$V2,
  breaks=10,
  main="distribution of d13C in feathers of vultures sampled",
  xlab="d13C absolute values (‰)",
  ylab="Frequency",
  col="lightblue",
  las=1,
  prob=T)
lines(density(d13C.abs$V2))
#'d15N
hist(d15N.abs$V2,
  breaks=10,
  main="distribution of d15N in feathers of vultures sampled",
  xlab="d15N absolute values (‰)",
  ylab="Frequency",
  col="lightblue",
  las=1,
  prob=T)
lines(density(d15N.abs$V2))
#'d34S
hist(d34S.abs$V2,
  breaks=10,
  main="distribution of d34S in feathers of vultures sampled",
  xlab="d34S absolute values (‰)",
  ylab="Frequency",
  col="lightblue",
  las=1,
  prob=T)
lines(density(d34S.abs$V2))
#'d13C,d15N and d34S values for indvs. vultures sampled were centered around the
#'mean value, this meant that diet and diet source in vultures as given by abs
#'isotope values in feathers did not change individually, however sample sizes 
#'were not high enough to confirm this   
#'---
#'

#'forage category biomass estimation for Serengeti national park
#'import data
SER.biomass<-read.csv("Biomass_SER.csv", header = T)
SGR.biomass<-read.csv("Biomass_SGR.csv", header = T)

#'view data
glimpse(SER.biomass)
glimpse(SGR.biomass)

#'basic plots to visualize differences in biomass
#'SER 
ggplot(SER.biomass, aes(x=forage.category, y=log(count), fill=forage.category))+
  geom_boxplot()+
  ggtitle("Categorized absolute counts for common herbivores in Serengeti N.P")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Logged herbivore numbers")+
  xlab("Herbivore forage categories")+
  theme_classic()

#'Grazing herbivores had higher absolute numbers in Serengeti National Park
#'compared to browsing and mixed feeding herbivores
#'---
#' 

#'SGR
ggplot(SGR.biomass, aes(x=forage.category, y=log(count), fill=forage.category))+
  geom_boxplot()+
  ggtitle("Categorized absolute counts for common herbivores in Selous G.R")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Logged herbivore numbers")+
  xlab("Herbivore forage categories")+
  theme_classic()

#'Grazing herbivores had higher absolute numbers in Selous Game Reserve
#'compared to browsing and mixed feeding herbivores
#'---
#' 


#'Regression analyses
#'run regression models to look at d13C,d15N and d34S variation in vulture 
#'feathers and relationship between d13C and C/N ratios.
#'we suspect differences in isotopic compositions are driven by species and 
#'location
f.AWB<-filter(bf.iso_data, species=='AWB')
f.RPV.ser<-filter(bf.iso_data, species=='RPV')
f.AWB.rgr<-filter(f.AWB, location=='RGR')
f.AWB.sgr<-filter(f.AWB, location=='SGR')
f.AWB.ser<-filter(f.AWB, location=='SER')
#'---
#'

#'AWBs
#'d13C model
d13C.RGR.md<-lm(f.AWB.rgr[f.AWB.rgr$t.subset=="basal.barbs",
                            ]$d13C~f.AWB.rgr[f.AWB.rgr$t.subset=="proximal.barbs",]
                $d13C)
#'d15N model
d15N.RGR.md<-lm(f.AWB.rgr[f.AWB.rgr$t.subset=="basal.barbs",
                            ]$d15N~f.AWB.rgr[f.AWB.rgr$t.subset=="proximal.barbs",]
                $d15N)

#'d34S model
d34S.RGR.md<-lm(f.AWB.rgr[f.AWB.rgr$t.subset=="basal.barbs",
                            ]$d34S~f.AWB.rgr[f.AWB.rgr$t.subset=="proximal.barbs",]
                $d34S)
#'d13C model
d13C.SER.md<-lm(f.AWB.ser[f.AWB.ser$t.subset=="basal.barbs",
                          ]$d13C~f.AWB.ser[f.AWB.ser$t.subset=="proximal.barbs",]
                $d13C)
#'d15N model
d15N.SER.md<-lm(f.AWB.ser[f.AWB.ser$t.subset=="basal.barbs",
                          ]$d15N~f.AWB.ser[f.AWB.ser$t.subset=="proximal.barbs",]
                $d15N)

#'d34S model
d34S.SER.md<-lm(f.AWB.ser[f.AWB.ser$t.subset=="basal.barbs",
                          ]$d34S~f.AWB.ser[f.AWB.ser$t.subset=="proximal.barbs",]
                $d34S)
#'d13C model
d13C.SGR.md<-lm(f.AWB.sgr[f.AWB.sgr$t.subset=="basal.barbs",
                          ]$d13C~f.AWB.sgr[f.AWB.sgr$t.subset=="proximal.barbs",]
                $d13C)
#'d15N model
d15N.SGR.md<-lm(f.AWB.sgr[f.AWB.sgr$t.subset=="basal.barbs",
                          ]$d15N~f.AWB.sgr[f.AWB.sgr$t.subset=="proximal.barbs",]
                $d15N)

#'d34S model
d34S.SGR.md<-lm(f.AWB.sgr[f.AWB.sgr$t.subset=="basal.barbs",
                          ]$d34S~f.AWB.sgr[f.AWB.sgr$t.subset=="proximal.barbs",]
                $d34S)
#'RPV
#'d13C model2
d13C.SER.md2<-lm(f.RPV.ser[f.RPV.ser$t.subset=="basal.barbs",
                      ]$d13C~f.RPV.ser[f.RPV.ser$t.subset=="proximal.barbs",]
                $d13C)
#'d15N model
d15N.SER.md2<-lm(f.RPV.ser[f.RPV.ser$t.subset=="basal.barbs",
                      ]$d15N~f.RPV.ser[f.RPV.ser$t.subset=="proximal.barbs",]
                $d15N)

#'d34S model
d34S.SER.md2<-lm(f.RPV.ser[f.RPV.ser$t.subset=="basal.barbs",
                      ]$d34S~f.RPV.ser[f.RPV.ser$t.subset=="proximal.barbs",]
                $d34S)

#'d13C vs C:N model for AWBs
dAWB_CN<-lm(d13C~C.N.ratio, data = f.AWB[f.AWB$t.subset !='basal.barbs',])

#'d13C vs C:N model for RPVs
dRPV_CN<-lm(d13C~C.N.ratio, data = f.RPV.ser[f.RPV.ser$t.subset !='basal.barbs',])

#'extract coeff(s) from regression models
coeff1<-coefficients(d13C.RGR.md)
coeff2<-coefficients(d15N.RGR.md)
coeff3<-coefficients(d34S.RGR.md)
coeff4<-coefficients(d13C.SER.md)
coeff5<-coefficients(d15N.SER.md)
coeff6<-coefficients(d34S.SER.md)
coeff7<-coefficients(d13C.SGR.md)
coeff8<-coefficients(d15N.SGR.md)
coeff9<-coefficients(d34S.SGR.md)
coeff10<-coefficients(d13C.SER.md2)
coeff11<-coefficients(d15N.SER.md2)
coeff12<-coefficients(d34S.SER.md2)
coeff13<-coefficients(dAWB_CN)
coeff14<-coefficients(dRPV_CN)

#'put coeff(s) into equations
eq1<-paste0("y = ", round(coeff1[2],1), "*x+", round(coeff1[1],1))
eq2<-paste0("y = ", round(coeff2[2],1), "*x+", round(coeff2[1],1))
eq3<-paste0("y = ", round(coeff3[2],1), "*x+", round(coeff3[1],1))
eq4<-paste0("y = ", round(coeff4[2],1), "*x+", round(coeff4[1],1))
eq5<-paste0("y = ", round(coeff5[2],1), "*x+", round(coeff5[1],1))
eq6<-paste0("y = ", round(coeff6[2],1), "*x+", round(coeff6[1],1))
eq7<-paste0("y = ", round(coeff7[2],1), "*x+", round(coeff7[1],1))
eq8<-paste0("y = ", round(coeff8[2],1), "*x+", round(coeff8[1],1))
eq9<-paste0("y = ", round(coeff9[2],1), "*x+", round(coeff9[1],1))
eq10<-paste0("y = ", round(coeff10[2],1), "*x+", round(coeff10[1],1))
eq11<-paste0("y = ", round(coeff11[2],1), "*x+", round(coeff11[1],1))
eq12<-paste0("y = ", round(coeff12[2],1), "*x+", round(coeff12[1],1))
eq13<-paste0("y = ", round(coeff13[2],1), "*x+", round(coeff13[1],1))
eq14<-paste0("y = ", round(coeff14[2],1), "*x+", round(coeff14[1],1))

#'plot regressions
#'AWBs in RGR
#'d13C at feather proximal and basal barbs
plot(f.AWB.rgr[f.AWB.rgr$t.subset=="f.root",
                 ]$d13C~f.AWB.rgr[f.AWB.rgr$t.subset=="f.tip",]$d13C, 
     sub=eq1,
     xlab = "d13C in roots (‰)",
     ylab = "d13C in tips (‰)",
     col = f.AWB.rgr$color,pch=f.AWB.rgr$symbol)
legend("topleft",legend=c("AWB","RGR"),
       col=c("red","black"),
       pch=c(15,0))
abline(d13C.RGR.md, col='blue', lty=2)

#'d15N at feather vane tips and roots
plot(f.AWB.rgr[f.AWB.rgr$t.subset=="f.root",
                 ]$d15N~f.AWB.rgr[f.AWB.rgr$t.subset=="f.tip",]$d15N, 
     sub=eq2,
     xlab = "d15N in roots (‰)",
     ylab = "d15N in tips (‰)",
     col = f.AWB.rgr$color,pch=f.AWB.rgr$symbol)
legend("topleft",legend=c("AWB","RGR"),
       col=c("red","black"),
       pch=c(15,0))
abline(d15N.RGR.md, col='blue', lty=2)

#'d34S at feather vane tips and roots
plot(f.AWB.rgr[f.AWB.rgr$t.subset=="f.root",
                 ]$d34S~f.AWB.rgr[f.AWB.rgr$t.subset=="f.tip",]$d34S,
     sub=eq3,
     xlab = "d34S in root (‰)",
     ylab = "d34S in tip (‰)",
     col = f.AWB.rgr$color,pch=f.AWB.rgr$symbol)
legend("topleft",legend=c("AWB","RGR"),
       col=c("red","black"),
       pch=c(15,0))
abline(d34S.RGR.md, col='blue', lty=2)

#'AWBs in SER
#'d13C feather vane tips and roots
plot(f.AWB.ser[f.AWB.ser$t.subset=="f.root",
               ]$d13C~f.AWB.ser[f.AWB.ser$t.subset=="f.tip",]$d13C, 
     sub=eq4,
     xlab = "d13C in roots (‰)",
     ylab = "d13C in tips (‰)",
     col = f.AWB.ser$color,pch=f.AWB.ser$symbol)
legend("topleft",legend=c("AWB","SER"),
       col=c("red","black"),
       pch=c(15,2))
abline(d13C.SER.md, col='blue', lty=2)

#'d15N at feather vane tips and roots
plot(f.AWB.ser[f.AWB.ser$t.subset=="f.root",
               ]$d15N~f.AWB.ser[f.AWB.ser$t.subset=="f.tip",]$d15N, 
     sub=eq5,
     xlab = "d15N in roots (‰)",
     ylab = "d15N in tips (‰)",
     col = f.AWB.ser$color,pch=f.AWB.ser$symbol)
legend("topleft",legend=c("AWB","SER"),
       col=c("red","black"),
       pch=c(15,2))
abline(d15N.SER.md, col='blue', lty=2)

#'d34S at feather vane tips and roots
plot(f.AWB.ser[f.AWB.ser$t.subset=="f.root",
               ]$d34S~f.AWB.ser[f.AWB.ser$t.subset=="f.tip",]$d34S,
     sub=eq6,
     xlab = "d34S in root (‰)",
     ylab = "d34S in tip (‰)",
     col = f.AWB.ser$color,pch=f.AWB.ser$symbol)
legend("topleft",legend=c("AWB","SER"),
       col=c("red","black"),
       pch=c(15,2))
abline(d34S.SER.md, col='blue', lty=2)

#'AWBs in SGR
#'d13C feather vane tips and roots
plot(f.AWB.sgr[f.AWB.sgr$t.subset=="f.root",
               ]$d13C~f.AWB.sgr[f.AWB.sgr$t.subset=="f.tip",]$d13C, 
     sub=eq7,
     xlab = "d13C in roots (‰)",
     ylab = "d13C in tips (‰)",
     col = f.AWB.sgr$color,pch=f.AWB.sgr$symbol)
legend("topleft",legend=c("AWB","SGR"),
       col=c("red","black"),
       pch=c(15,1))
abline(d13C.SGR.md, col='blue', lty=2)

#'d15N at feather vane tips and roots
plot(f.AWB.sgr[f.AWB.sgr$t.subset=="f.root",
               ]$d15N~f.AWB.sgr[f.AWB.sgr$t.subset=="f.tip",]$d15N, 
     sub=eq8,
     xlab = "d15N in roots (‰)",
     ylab = "d15N in tips (‰)",
     col = f.AWB.sgr$color,pch=f.AWB.sgr$symbol)
legend("topleft",legend=c("AWB","SGR"),
       col=c("red","black"),
       pch=c(15,1))
abline(d15N.SGR.md, col='blue', lty=2)

#'d34S at feather vane tips and roots
plot(f.AWB.sgr[f.AWB.sgr$t.subset=="f.root",
               ]$d34S~f.AWB.sgr[f.AWB.sgr$t.subset=="f.tip",]$d34S,
     sub=eq9,
     xlab = "d34S in root (‰)",
     ylab = "d34S in tip (‰)",
     col = f.AWB.sgr$color,pch=f.AWB.sgr$symbol)
legend("topleft",legend=c("AWB","SGR"),
       col=c("red","black"),
       pch=c(15,1))
abline(d34S.SGR.md, col='blue', lty=2)

#'RPVs
#'d13C at feather vane tips and roots
plot(f.RPV.ser[f.RPV.ser$t.subset=="f.root",
           ]$d13C~f.RPV.ser[f.RPV.ser$t.subset=="f.tip",]$d13C, 
     sub=eq10,
     xlab = "d13C in roots (‰)",
     ylab = "d13C in tips (‰)",
     col = f.RPV.ser$color,pch=f.RPV.ser$symbol)
legend("topleft",legend=c("RPV","SER"),
       col=c("blue","black"),
       pch=c(15,2))
abline(d13C.SER.md2, col='red', lty=2)

#'d15N at feather vane tips and roots
plot(f.RPV.ser[f.RPV.ser$t.subset=="f.root",
           ]$d15N~f.RPV.ser[f.RPV.ser$t.subset=="f.tip",]$d15N, 
     sub=eq11,
     xlab = "d15N in roots (‰)",
     ylab = "d15N in tips (‰)",
     col = f.RPV.ser$color,pch=f.RPV.ser$symbol)
legend("topleft",legend=c("RPV","SER"),
       col=c("blue","black"),
       pch=c(15,2))
abline(d15N.SER.md2, col='red', lty=2)

#'d34S at feather vane tips and roots
plot(f.RPV.ser[f.RPV.ser$t.subset=="f.root",
           ]$d34S~f.RPV.ser[f.RPV.ser$t.subset=="f.tip",]$d34S,
     sub=eq12,
     xlab = "d34S in root (‰)",
     ylab = "d34S in tip (‰)",
     col = f.RPV.ser$color,pch=f.RPV.ser$symbol)
legend("topleft",legend=c("RPV","SER"),
       col=c("blue","black"),
       pch=c(15,2))
abline(d34S.SER.md2, col='red', lty=2)

#'d13C vs. CN plots for AWBs
attach(f.AWB[f.AWB$t.subset !='basal.barbs',])
plot(f.AWB$C.N.ratio, f.AWB$d13C, pch=1, col="blue", 
  xlab = "C:N ratios in proximal barbs",
  ylab = "d13C in proximal barbs")
abline(dAWB_CN, col='red', lty=2)

#'d13c vs. CN plots for RPVs
attach(f.RPV.ser[f.RPV.ser$t.subset !='basal.barbs',])
plot(f.RPV.ser$C.N.ratio, f.RPV.ser$d13C, pch=1, col="blue",
  xlab = "C:N ratios in proximal barbs",
  ylab = "d13C in proximal barbs")
abline(dRPV_CN, col='red', lty=2)

#'summary of regression models for AWBs in RGR
#'d13C
summary(d13C.RGR.md)
#'d13C values for AWBs in RGR were not significantly different btn feather vane 
#'tips and roots. F-stat 0.16 on 1 and 3 DF, P value = 0.72, R2 = 0.05
#'---
#'

#'d15N
summary(d15N.RGR.md)
#'d15N values for AWBs in RGR indicated a significant difference of 0.76‰ 
#'more d15N in feather vane tips than roots. F-stat 116.1 on 1 and 3 DF, 
#'P value = 0.0017, R2 = 0.97 
#'---
#'    

#'d34S
summary(d34S.RGR.md)
#'d34S values for AWBs in RGR indicated no difference in d34S in feather vane 
#'tips and roots. F-stat 5.418 on 1 and 3 DF, P value = 0.1, R2 = 0.64
#'---
#'

#'summary of regression models for AWBs in SER
#'d13C
summary(d13C.SER.md)
#'d13C values for AWBs in SER were not significantly different btn feather vane 
#'tips and roots. F-stat 2.409 on 1 and 10 DF, P value = 0.15, R2 = 0.19
#'---
#'

#'d15N
summary(d15N.SER.md)
#'d15N values for AWBs in SER indicated no difference in d15N btn feather vane 
#'tips than roots. F-stat 3.499 on 1 and 10 DF, P value = 0.09, R2 = 0.26 
#'---
#'    

#'d34S
summary(d34S.SER.md)
#'d34S values for AWBs in SER indicated a significant difference of 0.6‰ more
#'d34S in feather vane tips than roots. F-stat 8.866 on 1 and 10 DF, 
#'P value = 0.014, R2 = 0.47
#'---
#'

#'summary of regression models for AWBs in SGR
#'d13C
summary(d13C.SGR.md)
#'d13C values for AWBs in SGR were not significantly different btn feather vane 
#'tips and roots. F-stat 0.05 on 1 and 4 DF, P value = 0.83, R2 = 0.012
#'---
#'

#'d15N
summary(d15N.SGR.md)
#'d15N values for AWBs in SGR indicated no difference in d15N btn feather vane 
#'tips than roots. F-stat 2.45 on 1 and 4 DF, P value = 0.2, R2 = 0.38 
#'---
#'    

#'d34S
summary(d34S.SGR.md)
#'d34S values for AWBs in SGR indicated a significant difference of 0.7‰ more
#'d34S in feather vane tips than roots. F-stat 23.06 on 1 and 4 DF, 
#'P value = 0.0086, R2 = 0.85
#'---
#'

#'summary of regression models for RPVs in SER
#'d13C
summary(d13C.SER.md2)
#'d13C values for RPVs in SER were not significantly different btn feather vane 
#'tips and roots. F-stat 0.29 on 1 and 7 DF, P value = 0.6, R2 = 0.04
#'---
#'

#'d15N
summary(d15N.SER.md2)
#'d15N values for RPVs in SER indicated no difference in d15N btn feather vane 
#'tips than roots. F-stat 2.07 on 1 and 7 DF, P value = 0.2, R2 = 0.23 
#'---
#'    

#'d34S
summary(d34S.SER.md2)
#'d34S values for RPVs in SER indicated no significant difference in d34S in
#'feather vane tips than roots. F-stat 4.88 on 1 and 7 DF, P value = 0.06, 
#'R2 = 0.41
#'---
#'

#'General linear models
#'run linear models to estbalish diet composition and source in vultures within
#'sampled PAs
#'model with d13C signals from vulture tissues to establish whether diet varies
#'as a function of location, vulture species, tissue type and an interaction 
#'btn tissue type & vulture species

#'run model
d13C.model1<-glm(d13C~location+species+t.subset+(t.subset*species), 
  data = bf.iso_data[bf.iso_data$location !='RGR',])

#'subset the data to run model with feathers as only tissue type
bf.iso_data.f<-filter(bf.iso_data, tissue.type=="Feather")

#'model with feathers as the only tissue type #After review there was no need!
d13C.model11<-glm(d13C~location+species+t.subset+(t.subset*species),
  data = bf.iso_data.f[bf.iso_data.f$location !='RGR',])

#'summary of the model
summary(d13C.model1)
#'diet in vultures (-11.62‰) in SER differed significantly from those in SGR
#'by (-1.58‰) (P value = 8.42e-05). Diet did not vary in the different 
#'vulture species sampled (P value = 0.57). Diet composition in blood (-11.62‰) 
#'was significantly differently from that informed by feather vane roots by 
#'+1.22‰ (P value = 0.00995). RPVs had significantly more carbon at feather vane
#'tips compared to AWBs (P value = 0.04166) 
#'---
#' 

#'summary of model11
summary(d13C.model11)
#'this model was run to look at the robustness of feathers in identifying 
#'temporal diet differences.
#'There were consistencies in identifying differences in diet across space SER 
#'(-10.4712‰) and SGR (-1.3779‰) P = 0.005. There was significantly less 
#'carbon in feather vane tips (-0.9032‰) than roots P = 0.03, and RPVs vultures
#'had significantly more carbon in feather vane tips +2.30‰ compared to AWBs
#'P = 0.004.
#'---
#' 

#'model with d15N signals from vulture tissues to establish whether trophic 
#'ecology varies as a function of location, vulture species, tissue subset and 
#'an interaction btn tissue subset & vulture species. 
d15N.model2<-glm(d15N~location+species+t.subset+(t.subset*species),
                 data = bf.iso_data[bf.iso_data$location !='RGR',])

#'model with feathers as the only tissue type #After review there was no need!
d15N.model22<-glm(d15N~location+species+t.subset+(t.subset*species),
                  data = bf.iso_data.f[bf.iso_data.f$location !='RGR',])

#'summary of model
summary(d15N.model2)
#'Average trophic level of prey items fed on by vultures is SER (10.80‰) was not 
#'significantly different from those fed on by vultures in SGR (P = 0.85).
#'Overtime there was a significant difference between trophic level of prey items
#'in blood and feather vane roots (+1.41‰) P = 5.92e-05), it was also 
#'significantly different from feather vane tips (+1.06‰) P = 0.002. 
#'---
#'
#'summary of model22
summary(d15N.model22)
#'this model was run to look at the robustness of feathers in identifying 
#'temporal differences in the average trophic level of prey items fed on by 
#'vultures.
#'Consistent with the blood model, there were no differences on the average 
#'trophic level of prey items fed by vultures across SER and SGR. Overtime given 
#'by feather vane tipos and roots there was no difference in the trophic level
#'of prey items fed on by vultures.
#'---
#' 

#'model with d34S signals from vulture feathers to establish diet source for 
#'vultures sampled as a function of location, vultures species, tissue type and 
#'an interaction btn tissue type & vulture species  
d34S.model3<-glm(d34S~location+species+t.subset+(t.subset*species),
                 data = bf.iso_data[bf.iso_data$location !='RGR',])

#'model with feathers as the only tissue type #After review there was no need!
d34S.model33<-glm(d34S~location+species+t.subset+(t.subset*species),
                  data = bf.iso_data.f[bf.iso_data.f$location !='RGR',])
#'summary of model
summary(d34S.model3)
#'diet source (10.1‰) for vultures in SER differed significantly from those in 
#'SGR by (+3.21‰) P = 2.17e-14. Overtime given by blood and feathers there was 
#'difference in diet source for vultures in SER and SGR.
#'---
#'  
#'summary of model33
summary(d34S.model33)
#'this model was run to look at the robustness of feathers in identifying 
#'differences in the location of prey items fed on by vultures in SER and SGR.
#'Consistent with model with blood, there was a significant difference in diet 
#'source across locations SER (9.13‰) and SGR (+3.60‰) P = 2.89e-10. There were 
#'differences in diet source across species and over time as was seen in model 
#'with blood. 
#'---
#' 

#'model diagnostics
#'check how much variation in response variable was accounted for by explanatory 
#'variables in model
#'d13C
autoplot(d13C.model1, smooth.colour = NA)
#'Top right plot indicated that there was no difference between expectations  
#'under normal distribution and residuals i.e. our residuals met the assumptions
#'of normality. 
#'---
#'

#'d15N
autoplot(d15N.model2, smooth.colour = NA)
#'Top right plot indicated that there was no difference between expectations  
#'under normal distribution and residuals i.e. our residuals met the assumptions
#'of normality. 
#'---
#'

#'d34S
autoplot(d34S.model3, smooth.colour = NA)
#'Top right plot indicated that there was no difference between expectations  
#'under normal distribution and residuals i.e. our residuals met the assumptions
#'of normality.
#'The use of the top right plot in our model diagnostics was a useful comparison
#'for homoscedasticity because our data consisted of small sample sizes
#'---
#'

#'Third objective
#'calculate means and standard deviation for browsing and grazing sources
Sources<-read.csv("V.sources.csv", header = T)
#'d13C, 
Sources%>%
  group_by(source)%>%
  summarise(meand13C = mean(d13C),
            SDd13C = sd(d13C))
#'d15N
Sources%>%
  group_by(source)%>%
  summarise(meand15N = mean(d15N),
            SDd15N = sd(d15N))
#'d34S
Sources%>%
  group_by(source)%>%
  summarise(meand34S = mean(d34S),
            SDd34S = sd(d34S))
#'save generated means and SD in different files for mixSIAR use
#'---
#'

#'Background of data for third objective
#'we have 1 fixed independent categorical variable (species of vulture) for which
#'we wish to estimate the mixing proportion of prey items (browsers and grazers)
#'in their diet using 3 biotracers or isotopes (d13C, d15N, d34S)
#'---
#'

#'add consumer blood data
ARmix<-load_mix_data(filename = "AR.blood.csv", 
                      iso_names = c("d13C","d15N","d34S"),
                      factors = c("species"),
                      fac_random = c(F),
                      fac_nested = c(F),
                      cont_effects = NULL)
#'add sources data
ARsources<-load_source_data(filename = "sources.csv", 
                             source_factors = NULL,
                             conc_dep = F,
                             data_type = "means", ARmix)
#'add blood discrimination data
ARdiscr<-load_discr_data(filename = "AR.blood.discr.csv", ARmix)

#'make isospace plot
plot_data(filename = "isospace_plot_blood", plot_save_pdf = T,
          plot_save_png = T, ARmix,ARsources,ARdiscr)

#'plot prior
#'used the generalist uninformative prior/alpha.prior = 1
plot_prior(alpha.prior = 1, ARsources)
#'Plot indicate that all combination of the sources are equally likely 
#'i.e. the probability of having browser or grazer as a source is the same,
#'hence the plot being flat - all value proportions are equally likely.
#'---
#'

#'JAGS model
#'set parameters for JAGS model
model_filename<-"MixSIAR_model.txt"
resid_err<- T
process_err<- F
write_JAGS_model(model_filename, resid_err, process_err, ARmix, ARsources)

#'conduct a run test run to check if data has loaded correctly and model is 
#'specified correctly
jags.1<-run_model(run = "test", ARmix, ARsources, ARdiscr, model_filename,
                  alpha.prior = 1, resid_err, process_err)
#'the test model ran, all the data is in and specified correctly
#'---
#'

#'we now run the full model
#'specify MCMC parameters for model
#'we desired the 'long' setting, chainlength = 300,000, Burnin = 200,000,
#'sample thin = 100, chains = 3
jags.1<-run_model(run = "long",ARmix, ARsources, ARdiscr, model_filename,
                  alpha.prior = 1, resid_err, process_err)

#'specify outputs from jags model
output_options<-list(summary_save = T,
                     summary_name = "summary_statistics_blood",
                     sup_post = F,
                     plot_post_save_pdf = T,
                     plot_post_name = "posterior_density_blood",
                     sup_pairs = F,
                     plot_pairs_save_pdf = T,
                     plot_pairs_name = "pairs_plot",
                     sup_xy = T,
                     plot_xy_save_pdf = T,
                     plot_xy_name = "xy_plot",
                     gelman = T,
                     heidel = F,
                     geweke = T,
                     diag_save = T,
                     diag_name = "diagnostics_blood",
                     indiv_effect = F,
                     plot_post_save_png = T,
                     plot_pairs_save_png = T,
                     plot_xy_save_png = T)
#'retrieve outputs from jags model
output_JAGS(jags.1, ARmix, ARsources, output_options) 

#'add consumer feather data
ARf.mix<-load_mix_data(filename = "AR.feathers.csv", 
                     iso_names = c("d13C","d15N","d34S"),
                     factors = c("species"),
                     fac_random = c(F),
                     fac_nested = c(F),
                     cont_effects = NULL)
#'add sources data
ARf.sources<-load_source_data(filename = "sources.csv", 
                            source_factors = NULL,
                            conc_dep = F,
                            data_type = "means", ARf.mix)

#'add feathers discrimination data
ARf.discr<-load_discr_data(filename = "AR.feathers.discr.csv", ARf.mix)

#'make isospace plot
plot_data(filename = "isospace_plot_feathers", plot_save_pdf = T,
          plot_save_png = T, ARf.mix,ARf.sources,ARf.discr)

#'plot prior
#'used the generalist uninformative prior/alpha.prior = 1
plot_prior(alpha.prior = 1, ARf.sources)
#'Plot indicate that all combination of the sources are equally alike 
#'i.e. the probability of having browser or grazer as a source is the same,
#'hence the plot being flat - all value proportions are equally likely.
#'---
#'

#'JAGS model
#'set parameters for JAGS model
model_filename<-"MixSIAR_model2.txt"
resid_err<- T
process_err<- F
write_JAGS_model(model_filename, resid_err, process_err, ARf.mix, ARf.sources)

#'we now run the full model
#'specify MCMC parameters for model
#'we desired the 'long' setting, chainlength = 300,000, Burnin = 200,000,
#'sample thin = 100, chains = 3
jags.2<-run_model(run = "long", ARf.mix, ARf.sources, ARf.discr, model_filename,
                  alpha.prior = 1, resid_err, process_err)
#'specify outputs from jags model
output_options2<-list(summary_save = T,
                     summary_name = "summary_statistics_feathers",
                     sup_post = F,
                     plot_post_save_pdf = T,
                     plot_post_name = "posterior_density_feathers",
                     sup_pairs = F,
                     plot_pairs_save_pdf = T,
                     plot_pairs_name = "pairs_plot_feathers",
                     sup_xy = T,
                     plot_xy_save_pdf = T,
                     plot_xy_name = "xy_plot2",
                     gelman = T,
                     heidel = F,
                     geweke = T,
                     diag_save = T,
                     diag_name = "diagnostics_feathers",
                     indiv_effect = F,
                     plot_post_save_png = T,
                     plot_pairs_save_png = T,
                     plot_xy_save_png = T)
#'retrieve outputs from jags model
output_JAGS(jags.2, ARf.mix, ARf.sources, output_options2) 
#'---
#'

#'Get citations from all packages that have been used
citations <- function(includeURL = TRUE, includeRStudio = TRUE) {
  if(includeRStudio == TRUE) {
    ref.rstudio <- RStudio.Version()$citation
    if(includeURL == FALSE) {
      ref.rstudio$url = NULL;
    }
    print(ref.rstudio, style = 'text')
    cat('\\n')
  }
  
  cit.list <- c('base', names(sessionInfo()$otherPkgs))
  for(i in 1:length(cit.list)) {
    ref <- citation(cit.list[i])
    if(includeURL == FALSE) {
      ref$url = NULL;
    }
    print(ref, style = 'text')
    cat('\\n')
  }
}

####################
##  END OF SCRIPT ##
####################
