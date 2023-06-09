#####--------------------------#####
# Habitat and phylogenetic analyses
#####--------------------------#####

# Note: if downloading and working through code, I've set up a way to make it easier to move through sections. All subsections that start with 5 "#" signs and a description should have another 5 "#" signs in following row. If opening in R, there should be a downward arrow next to the script line number next to this row - clicking that will hide all of the code in between that line and the next line that only has 5 "#" signs. 

##### Libraries
#####
library(tidyr)
library(ggbiplot)
library(reshape2)
library(ggplot2)
library(ggfortify)
library("wesanderson")
library(Hmisc)
library(PerformanceAnalytics)
library(ggrepel)
library(rjags)
library(ape)
library(phylotate)
library(phytools)
library(nlme)
library(dplyr)
library(RColorBrewer)
library(geiger)
library("gridExtra")
#####

##### Formulas
#####
  # calculate standard error
  sem<-function(x){sd(x,na.rm=T)/sqrt(length(x))}

  # calculate 25-D Euclidean distances
  Euc_dist_25_PCs<-function(x){
    sqrt(
      abs(x[x$soil_type=="S","PC1"]-x[x$soil_type=="NS","PC1"])^2 +
        abs(x[x$soil_type=="S","PC2"]-x[x$soil_type=="NS","PC2"])^2 +
        abs(x[x$soil_type=="S","PC3"]-x[x$soil_type=="NS","PC3"])^2 +
        abs(x[x$soil_type=="S","PC4"]-x[x$soil_type=="NS","PC4"])^2 +
        abs(x[x$soil_type=="S","PC5"]-x[x$soil_type=="NS","PC5"])^2 +
        abs(x[x$soil_type=="S","PC6"]-x[x$soil_type=="NS","PC6"])^2 +
        abs(x[x$soil_type=="S","PC7"]-x[x$soil_type=="NS","PC7"])^2 +
        abs(x[x$soil_type=="S","PC8"]-x[x$soil_type=="NS","PC8"])^2 +
        abs(x[x$soil_type=="S","PC9"]-x[x$soil_type=="NS","PC9"])^2 +
        abs(x[x$soil_type=="S","PC10"]-x[x$soil_type=="NS","PC10"])^2 +
        abs(x[x$soil_type=="S","PC11"]-x[x$soil_type=="NS","PC11"])^2 +
        abs(x[x$soil_type=="S","PC12"]-x[x$soil_type=="NS","PC12"])^2 +
        abs(x[x$soil_type=="S","PC13"]-x[x$soil_type=="NS","PC13"])^2 +
        abs(x[x$soil_type=="S","PC14"]-x[x$soil_type=="NS","PC14"])^2 +
        abs(x[x$soil_type=="S","PC15"]-x[x$soil_type=="NS","PC15"])^2 +
        abs(x[x$soil_type=="S","PC16"]-x[x$soil_type=="NS","PC16"])^2 +
        abs(x[x$soil_type=="S","PC17"]-x[x$soil_type=="NS","PC17"])^2 +
        abs(x[x$soil_type=="S","PC18"]-x[x$soil_type=="NS","PC18"])^2 +
        abs(x[x$soil_type=="S","PC19"]-x[x$soil_type=="NS","PC19"])^2 +
        abs(x[x$soil_type=="S","PC20"]-x[x$soil_type=="NS","PC20"])^2 +
        abs(x[x$soil_type=="S","PC21"]-x[x$soil_type=="NS","PC21"])^2 +
        abs(x[x$soil_type=="S","PC22"]-x[x$soil_type=="NS","PC22"])^2 +
        abs(x[x$soil_type=="S","PC23"]-x[x$soil_type=="NS","PC23"])^2 +
        abs(x[x$soil_type=="S","PC24"]-x[x$soil_type=="NS","PC24"])^2 +
        abs(x[x$soil_type=="S","PC25"]-x[x$soil_type=="NS","PC25"])^2 
    )
  }  
#####
  
##### Upload data
#####
  setwd("your_directory")
  
  # soil data
  soil <- read.csv("soil_chemistry_texture.csv")
  
  # bare ground data
  comp <- read.csv("bare_ground.csv")
  
  # create geographic distance dataframe (from Table 1)
  geo.dist <- data.frame(pair_name = c("CACO","CABR","PLER","MGUT","COSP","COHT","TWILD","NAPB","NAHX","NAJP","NARS","CAGT","CLDV","LADI","MNUD","COGR","CABE"),
                         geo.dist = c(1.85,3.03,1.75,6.17,2.48,2.23,1.57,1.44,37.76,37.58,20.87,52.05,11.06,6.19,6.4,1.17,4.58))

  geo.dist

#####
  
#####-------------------------------######
#####-------------------------------######
#####  Phylogenetic analyses
#####-------------------------------######
  
  # read in Mr Bayes concensus tree
  tree.out<-read_annotated("tree.con.tre", format = "nexus")

  # Reroot tree such that Asterids and Rosids are in separate clades
  tree.root<-reroot(tree.out,20)
  
  # Ultrametricize tree with Grafen (1989) method
  tree9 <- compute.brlen(tree.root, method = "Granfen", power = 1) 
  plot(tree9)
  
  # Assign this tree as tree to use in soil analyses
  soil.tree <- tree9

  # Drop the two Navarretia pairs not used in the bare ground analyses
  tree18<-drop.tip(tree9, c("NARS","NAHX"))
  
  plot(tree18)
  
  # Assign this tree as tree to use in bare ground analyses
  bare.tree <- tree18

#####-------------------------------######  
#####-------------------------------######
#####  Habitat analyses: Do serpentine endemics occur in harsher serpentine soils than serpentine tolerators?
#####-------------------------------######
  
  ##### Dataframe ("aff") prep (subset out serpentine populations)
  #####
  aff <- soil %>% filter(soil_type=="S")
  aff$sp
  aff %>% group_by(E_T_NT) %>% dplyr::summarise(n()) # 17 taxa total, 8 E and 9 T
  
  aff <-aff %>%
    dplyr::select (sp, E_T_NT, CaMg,Ca_mg.kg,Mg_mg.kg,P_mg.kg,K_mg.kg,NH4_NP_mg.kg,Ni_mg.kg,Cr_mg.kg,Co_mg.kg,percent_sand,percent_silt,percent_clay)
  
  # put into long format
  affL <-aff %>%
    gather (.,"variable","value",3:14)
  
  # add spelled out variable names for figures
  affL$full_variable<-affL$variable
  affL$full_variable<-gsub("\\\\<CaMg\\\\>","Ca:Mg",affL$full_variable)
  affL$full_variable<-gsub("\\\\<Ca_mg.kg\\\\>","Ca (mg/kg) *",affL$full_variable)
  affL$full_variable<-gsub("\\\\<Mg_mg.kg\\\\>","Mg (mg/kg)",affL$full_variable)
  affL$full_variable<-gsub("\\\\<P_mg.kg\\\\>","P (mg/kg)",affL$full_variable)
  affL$full_variable<-gsub("\\\\<K_mg.kg\\\\>","K (mg/kg)",affL$full_variable)
  affL$full_variable<-gsub("\\\\<NH4_NP_mg.kg\\\\>","N (mg/kg)",affL$full_variable)
  affL$full_variable<-gsub("\\\\<Ni_mg.kg\\\\>","Ni (mg/kg)",affL$full_variable)
  affL$full_variable<-gsub("\\\\<Cr_mg.kg\\\\>","Cr (mg/kg)",affL$full_variable)
  affL$full_variable<-gsub("\\\\<Co_mg.kg\\\\>","Co (mg/kg)",affL$full_variable)
  affL$full_variable<-gsub("\\\\<percent_sand\\\\>","% Sand",affL$full_variable) 	
  affL$full_variable<-gsub("\\\\<percent_silt\\\\>","% Silt",affL$full_variable)
  affL$full_variable<-gsub("\\\\<percent_clay\\\\>","% Clay",affL$full_variable)
  
  # Relevel factors --
  affL$variable = factor(affL$variable, levels=c("CaMg","Ca_mg.kg","Mg_mg.kg","P_mg.kg","K_mg.kg","NH4_NP_mg.kg","Ni_mg.kg","Cr_mg.kg","Co_mg.kg","percent_sand","percent_silt","percent_clay"))
  affL$E_T_NT = factor(affL$E_T_NT,levels=c("T","E"))
  
  affL$full_variable = factor(affL$full_variable, levels=c("Ca (mg/kg) *","Mg (mg/kg)","Ca:Mg","P (mg/kg)","K (mg/kg)","N (mg/kg)","Ni (mg/kg)","Cr (mg/kg)","Co (mg/kg)","% Sand","% Silt","% Clay"))
  #####  
  
  ##### Figure 1 code
  #####
  # Facet plot of individual variables
  affL %>%
    ggplot (.,aes(E_T_NT,log(value),fill=E_T_NT)) +
    geom_boxplot (width = 0.6) +
    facet_wrap (~full_variable, scales="free_y",nrow=4)+
    xlab ("") +
    ylab ("") +
    theme_classic() +
    scale_fill_grey(start=.85,end=.4,name="Species type",labels=c("Tolerator","Endemic")) +
    theme(legend.key.size = unit(1,"line"), legend.position = "none", axis.text.x = element_text(angle = 45, vjust = .95, hjust = .95,size = rel(1.2)), axis.text.y = element_text(size = rel(1)), strip.text = element_text(size = rel(.8)), plot.margin = margin(t = .2, r = .4, b = .2, l = .2, unit = "cm"))+
    scale_x_discrete(labels = c("Tolerators","Endemics"))+
    ylab("Log-transformed variables")
  #####
  
##### PGLS 
  
  ##### Data frame prep
  #####
  # double check tree and data frame line up
  name.check(soil.tree, data.names = aff$sp) # looks good
  
  # Create log-transformed soil variables dataframe
  colnames(aff)
  aff.log<-log(aff[3:14])
  aff.log$E_T_NT <- aff$E_T_NT
  aff.log$sp <- aff$sp
  aff.log
  
  # Create data table to hold PGLS values
  pgls.table<-data.frame(Soil_variable=colnames(aff.log)[c(1:12)],
                         numDF=NA,
                         denDF=NA,
                         intercept_E_coeff=NA,
                         slope_T_coeff=NA,
                         F_value=NA,
                         p_value=NA)
  
  # Make sure dataframe rows match tree tips
  soil.tree$tip.label
  aff.log$sp
  aff.log$sp <- factor(aff.log$sp, levels=soil.tree$tip.label)
  aff.log<-aff.log[order(aff.log$sp),]
  row.names(aff.log)<-aff.log$sp
  #####
  
  ##### Run PGLS analyses
  #####
  # Mg 
  mg.pgls<-gls(Mg_mg.kg ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(mg.pgls)
  anova(mg.pgls)
  pgls.table[3,2:7]<-c(as.data.frame(anova(mg.pgls))[2,1], 15, as.data.frame(summary(mg.pgls)$tTable)[1,1], as.data.frame(summary(mg.pgls)$tTable)[2,1], as.data.frame(anova(mg.pgls))[2,2], as.data.frame(anova(mg.pgls))[2,3])
  
  
  # Ca
  ca.pgls<-gls(Ca_mg.kg ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(ca.pgls)
  anova(ca.pgls)
  as.data.frame(summary(ca.pgls)$tTable)[2,4] # to extract pvalue
  pgls.table[2,2:7]<-c(as.data.frame(anova(ca.pgls))[2,1], 15, as.data.frame(summary(ca.pgls)$tTable)[1,1], as.data.frame(summary(ca.pgls)$tTable)[2,1], as.data.frame(anova(ca.pgls))[2,2], as.data.frame(anova(ca.pgls))[2,3])
  
  
  # CaMg
  camg.pgls<-gls(CaMg ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(camg.pgls)
  pgls.table[1,2:7]<-c(as.data.frame(anova(camg.pgls))[2,1], 15, as.data.frame(summary(camg.pgls)$tTable)[1,1], as.data.frame(summary(camg.pgls)$tTable)[2,1], as.data.frame(anova(camg.pgls))[2,2], as.data.frame(anova(camg.pgls))[2,3])
  
  # P
  p.pgls<-gls(P_mg.kg ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(p.pgls)
  pgls.table[4,2:7]<-c(as.data.frame(anova(p.pgls))[2,1], 15, as.data.frame(summary(p.pgls)$tTable)[1,1], as.data.frame(summary(p.pgls)$tTable)[2,1], as.data.frame(anova(p.pgls))[2,2], as.data.frame(anova(p.pgls))[2,3])
  
  
  # K
  K.pgls<-gls(K_mg.kg ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(K.pgls)
  pgls.table[5,2:7]<-c(as.data.frame(anova(K.pgls))[2,1], 15, as.data.frame(summary(K.pgls)$tTable)[1,1], as.data.frame(summary(K.pgls)$tTable)[2,1], as.data.frame(anova(K.pgls))[2,2], as.data.frame(anova(K.pgls))[2,3])
  
  # NH4
  nh4.pgls<-gls(NH4_NP_mg.kg ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(nh4.pgls)
  pgls.table[6,2:7]<-c(as.data.frame(anova(nh4.pgls))[2,1], 15, as.data.frame(summary(nh4.pgls)$tTable)[1,1], as.data.frame(summary(nh4.pgls)$tTable)[2,1], as.data.frame(anova(nh4.pgls))[2,2], as.data.frame(anova(nh4.pgls))[2,3])
  
  
  # Ni
  ni.pgls<-gls(Ni_mg.kg ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(ni.pgls)
  pgls.table[7,2:7]<-c(as.data.frame(anova(ni.pgls))[2,1], 15, as.data.frame(summary(ni.pgls)$tTable)[1,1], as.data.frame(summary(ni.pgls)$tTable)[2,1], as.data.frame(anova(ni.pgls))[2,2], as.data.frame(anova(ni.pgls))[2,3])
  
  # Cr
  cr.pgls<-gls(Cr_mg.kg ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(cr.pgls)
  pgls.table[8,2:7]<-c(as.data.frame(anova(cr.pgls))[2,1], 15, as.data.frame(summary(cr.pgls)$tTable)[1,1], as.data.frame(summary(cr.pgls)$tTable)[2,1], as.data.frame(anova(cr.pgls))[2,2], as.data.frame(anova(cr.pgls))[2,3])
  
  
  # Co
  co.pgls<-gls(Co_mg.kg ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(co.pgls)
  pgls.table[9,2:7]<-c(as.data.frame(anova(co.pgls))[2,1], 15, as.data.frame(summary(co.pgls)$tTable)[1,1], as.data.frame(summary(co.pgls)$tTable)[2,1], as.data.frame(anova(co.pgls))[2,2], as.data.frame(anova(co.pgls))[2,3])
  
  # % sand
  sand.pgls<-gls(percent_sand ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(sand.pgls)
  sand.pgls
  pgls.table[10,2:7]<-c(as.data.frame(anova(sand.pgls))[2,1], 15, as.data.frame(summary(sand.pgls)$tTable)[1,1], as.data.frame(summary(sand.pgls)$tTable)[2,1], as.data.frame(anova(sand.pgls))[2,2], as.data.frame(anova(sand.pgls))[2,3])
  
  # % silt
  silt.pgls<-gls(percent_silt ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(silt.pgls)
  pgls.table[11,2:7]<-c(as.data.frame(anova(silt.pgls))[2,1], 15, as.data.frame(summary(silt.pgls)$tTable)[1,1], as.data.frame(summary(silt.pgls)$tTable)[2,1], as.data.frame(anova(silt.pgls))[2,2], as.data.frame(anova(silt.pgls))[2,3])
  
  # % clay
  clay.pgls<-gls(percent_clay ~ E_T_NT,correlation = corBrownian(phy = soil.tree), data = aff.log, method = "ML")
  summary(clay.pgls)
  pgls.table[12,2:7]<-c(as.data.frame(anova(clay.pgls))[2,1], 15, as.data.frame(summary(clay.pgls)$tTable)[1,1], as.data.frame(summary(clay.pgls)$tTable)[2,1], as.data.frame(anova(clay.pgls))[2,2], as.data.frame(anova(clay.pgls))[2,3])
  #####
  
  ##### PGLS summary table
  #####
  pgls.table
  pgls.table<-arrange(pgls.table, p_value)
  
  # sequential bonferonni
  bon <- rep(0,12)
  for ( i in 1:12){
    bon[i]<- round(0.05/(12 - i + 1),5)
  }
  pgls.table$bon.p <- bon
  pgls.table # Ca comes out as the only one that is significant after sequential bonferoni
  #####
  
#####-------------------------------######  
#####-------------------------------######
#####  Habitat analyses: Do serpentine endemics occur in less competitive serpentine habitats than serpentine tolerators?
#####-------------------------------###### 

  ##### Set up coancestry matrix
  #####
  G <- vcv(bare.tree, corr = T)
  G
  Ginv<-solve(G)
  Ginv
  rownames(Ginv) 
  #####

  ##### filter just the S populations and code data frame
  #####  
  comp.S <- comp %>% filter(., soil == "S")
  
  # Change the E/T to 1/0 -> 0 for T and 1 for E
  comp.S$E_T<-comp.S$pair_type
  comp.S$E_T<-gsub("\\\\<E\\\\>",1,comp.S$E_T)
  comp.S$E_T<-gsub("\\\\<T\\\\>",0,comp.S$E_T)
  comp.S$E_T <- as.numeric(comp.S$E_T)
  
  # tree built with no outgroup - depending on tree 
  bare.tree$tip.label
  comp.S$pair.no<-comp.S$pair_name
  comp.S$pair.no<-gsub("\\\\<CABE\\\\>","1",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<CABR\\\\>","2",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<CACO\\\\>","3",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<CAGT\\\\>","4",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<CLDV\\\\>","6",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<COGR\\\\>","10",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<COHT\\\\>","11",comp.S$pair.no) 
  comp.S$pair.no<-gsub("\\\\<COSP\\\\>","12",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<LADI\\\\>","9",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<MGUT\\\\>","14",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<MNUD\\\\>","13",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<NAJP\\\\>","7",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<NAPB\\\\>","8",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<PLER\\\\>","15",comp.S$pair.no)
  comp.S$pair.no<-gsub("\\\\<TWILD\\\\>","5",comp.S$pair.no)
  comp.S$pair.no <- as.numeric(comp.S$pair.no)
  
  # number of quadrats taken per population   
  comp.S %>%
    group_by(sp) %>%
    dplyr::summarise(n())
  #####  

  ##### JAGS dataframe  
  #####
  # need an E_T list for each species
  head(comp.S)
  et.species1<-comp.S[comp.S$quad==1,c("E_T","pair.no")] %>% arrange(.,pair.no)
  et.species1  
  
  data= list(
    #n.sp = length(as.double(unique(comp.S$species))),
    n = as.double(nrow(comp.S)),
    n.pts = as.double(comp.S$sum),
    y = as.double(comp.S$bare),
    E_T = as.double(et.species1$E_T),
    species = as.double(comp.S$pair.no),
    mu.0 = as.double(rep(0,15)),
    Ginv = Ginv
  )
  
  ## Figuring out the limits on sigma for the alpha and beta parameters in the beta distribution
  
  find.sigma<-function(mu) {
    (mu^2 - mu^3)/(mu)
  }
  
  mus<-seq(0.01,1,.01)
  mus
  find.sigma(mus)
  par(mfrow = c(1,1))
  plot(mus,find.sigma(mus)) # sigma has to be between 0 and .25 for alpha, and 0 and 1 for beta)
  #####
  
  ##### JAGS model
  #####
  sink("bare_phylo.R")
  cat("
      model{
      
      # Model actual data as binomial process
      
      for (i in 1:n){
      y[i] ~ dbin(phi[species[i]],n.pts[i])
      }
      
      # MODEL PROBABILITY OF BARE GROUND IN EACH POP (PHI) AS BETA DIST W/ MEAN FROM REGRESSION
      for (j in 1:15){
      mu[j] <- ilogit(beta.0[j] + beta.1 * E_T[j])
      
      a[j] <- (mu[j]^2 - mu[j]^3 - mu[j]*sigma2) / sigma2
      b[j] <- (mu[j] - 2 * mu[j]^2 + mu[j]^3 - sigma2 + mu[j]*sigma2) / sigma2
      
      phi[j] ~ dbeta(a[j],b[j])
      }
      
      # PRIORS
      beta.0 ~ dmnorm(mu.0, Ginv)
      
      beta.1 ~ dnorm(0, tau)
      tau <- 1 / sigma2b
      sigma2b ~ dunif(0,100)
      sigma2 ~ dunif(0,0.25) # MAX OF .25 SO THAT a[j] IS >0
      
      }
      ",fill = TRUE)
  sink()
  
  #####
  
  ##### Run JAGS
  #####
  n.adapt = 5000 # JAGS searches forright starting values, MCMC sampler methods, etc
  n.update = 10000 # burnin
  n.iter = 50000 # run length
  
  jm1 = jags.model("bare_phylo.R", data = data, n.chains = 3, n.adapt = n.adapt)
  update(jm1, n.iter = n.update)
  zm1 = coda.samples(jm1, variable.names = c("phi", "beta.0", "beta.1"), n.iter = n.iter, n.thin = 10)
  zj1 = jags.samples(jm1, variable.names = c("phi", "beta.0", "beta.1"), n.iter = n.iter, n.thin = 10)
  #####

  ##### Diagnostics
  #####
  summary(zm1)
  #plot(zm1) # trace and density plots of all parameters
  gelman.diag(zm1) 
  #####

  ##### Summary table for MCMC combined chains
  #####
  head(zm1[[1]])
  str(zm1)
  
  phi.sp1<-as.data.frame(rbind(zm1[[1]], zm1[[2]], zm1[[3]]) )
  head(phi.sp1)
  
  key<-comp.S[comp.S$quad == 1, c("pair.no","sp","pair_type")] %>% arrange(pair.no) # used for graphs
  #####
  

  ##### Graph of phi posterior for all species on same axis - just lines (Figure 2)
  #####
  comp.sum<-comp.S %>%
    mutate(perc.bare = bare / sum) %>%
    dplyr::group_by(sp, pair.no, E_T, pair_type) %>%
    summarise(avg = mean (perc.bare), low.sd = mean(perc.bare) - sd(perc.bare), high.sd = mean(perc.bare) + sd(perc.bare))
  
  comp.sum.E <- comp.sum %>% filter(., E_T ==1) %>% as.data.frame()
  comp.sum.T <- comp.sum %>% filter(., E_T ==0) %>% as.data.frame()
  
  par(mfrow = c(1,1))
  par(mar=c(3.5,3.5,1,1.5))
  plot(NA,  xlab = " ", xlim = c(0,1), ylim = c(0,20), ylab = " ",bty = "n",cex.axis=.8, cex.lab=2)
  mtext(side = 1, line = 2.3, expression(paste("Estimated proportion bare ground ","(",phi,")")), size = 2)
  mtext(side = 2, line = 2.3, "Probability density")
  lines(density(phi.sp1[,17+(comp.sum.E[1,"pair.no"]-1)]), col = "blue", lwd = 2)
  
  for(i in 2:8){
    lines(density(phi.sp1[,17+(comp.sum.E[i,"pair.no"]-1)]), col = "blue", lwd = 2)
  }
  for(i in 1:8){
    lines(density(phi.sp1[,17+(comp.sum.T[i,"pair.no"]-1)]), col = "orange", lwd = 2)
  }
  legend(.01, 20.5, legend = c("Tolerator","Endemic"), col = c("orange","blue"), lwd = c(2,2), cex = .8)
  #####
  
  ##### Facets of all the tolerator phis, w/ empirical mean overlayed (Appendix S1)
  #####
  par(mfrow=c(4,2))
  par(mar=c(4,4,3,2))
  for(i in 1:8){
    hist(phi.sp1[,17+(comp.sum.T[i,"pair.no"]-1)], main = paste(comp.sum.T[i,"sp"]), breaks = 15, freq = F, xlim = c(0,1), xlab = "",cex.lab = 1.5, cex.main = 1.5, ylab="", cex.axis = 1.5)
    mtext(side = 1, line = 2.7, expression(paste("Estimated proportion bare ground ","(",phi,")")))
    mtext(side = 2, line = 2.5, "Probability density")
    lines(density(phi.sp1[,17+(comp.sum.T[i,"pair.no"]-1)]), col = "orange", lwd = 2.5)
    abline(v = comp.sum.T[i, "avg"], col = "#E31A1C", lwd = 4)
  }
  #####

  ##### Facets of all endemic phis, w/ empirical mean overlayed (Appendix S2)
  #####
  par(mfrow=c(4,2))
  par(mar=c(4,4.2,3,2))
  for(i in 1:7){
    hist(phi.sp1[,17+(comp.sum.E[i,"pair.no"]-1)], main = paste(comp.sum.E[i,"sp"]), breaks = 12, freq = F, xlim = c(0,1),xlab = " ", cex.lab = 1.5, cex.main = 1.5, ylab=" ", cex.axis = 1.5)
    mtext(side = 1, line = 2.7, expression(paste("Estimated proportion bare ground ","(",phi,")")))
    mtext(side = 2, line = 2.5, "Probability density")
    lines(density(phi.sp1[,17+(comp.sum.E[i,"pair.no"]-1)]), col = "blue", lwd = 2.5)
    abline(v = comp.sum.E[i, "avg"], col = "#E31A1C", lwd = 2.5)
  }
  #####
  
  ##### Graph of beta.1 w/ credible intervals
  #####
  
  par(mfrow = c(1,1))
  par(mar = c(5,5,5,5))
  hist(phi.sp1[,"beta.1"], main = "beta.1 posterior", xlab = "beta.1", breaks = 100, freq = F,ylab="PDF")
  lines(density(phi.sp1[,"beta.1"]), col = "blue", lwd = "2")
  abline(v = quantile(zj1$beta.1,c(.025,.975))[[1]], col = "blue", lty = 2)
  abline(v = quantile(zj1$beta.1,c(.025,.975))[[2]], col = "blue", lty = 2)
  #####
  
  ##### % of beta.1 posterior greater than zero
  1-ecdf(zj1$beta.1)(0)

    
#####-------------------------------######  
#####-------------------------------######
#####  Habitat analyses: Do serpentine endemic sister taxa pairs have more divergence in soil chemistry and texture than serpentine tolerator sister taxa pairs?
#####-------------------------------###### 
  
  #####-----------------Multivariate soil distance-----------------#####
  
  ##### Get "pairs" dataframe ready 
  #####
  soil %>%
    group_by(pair_name, soil_type) %>%
    dplyr::summarise(n = n()) %>% as.data.frame() # double check all pairs (except shared ones) have 2 entries
  colnames(soil)
  
  pairs<-soil[c(1,3,4,5,9:34)]
  colnames(pairs)
  head(pairs)
  
  soil.tree$tip.label
  
  unique(pairs$pair_name)
  
  pairs$pair.no<-pairs$pair_name
  pairs$pair.no<-gsub("\\\\<CABE\\\\>","1",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<CABR\\\\>","2",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<CACO\\\\>","3",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<CAGT\\\\>","4",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<CLDV\\\\>","6",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<COGR\\\\>","10",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<COHT\\\\>","11",pairs$pair.no) 	
  pairs$pair.no<-gsub("\\\\<COSP\\\\>","12",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<LADI\\\\>","9",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<MGUT\\\\>","14",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<MNUD\\\\>","13",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<NAJP\\\\>","7",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<NAPB\\\\>","8",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<PLER\\\\>","15",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<NARS\\\\>","16",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<NAHX\\\\>","17",pairs$pair.no)
  pairs$pair.no<-gsub("\\\\<TWILD\\\\>","5",pairs$pair.no)
  
  pairs$pair_type<-pairs$pair_name
  pairs$pair_type<-gsub("\\\\<CABE\\\\>","E",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<CABR\\\\>","T",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<CACO\\\\>","T",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<CAGT\\\\>","E",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<CLDV\\\\>","E",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<COGR\\\\>","E",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<COHT\\\\>","T",pairs$pair_type) ;	
  pairs$pair_type<-gsub("\\\\<COSP\\\\>","T",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<LADI\\\\>","E",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<MGUT\\\\>","T",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<MNUD\\\\>","E",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<NAJP\\\\>","E",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<NAPB\\\\>","T",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<PLER\\\\>","T",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<TWILD\\\\>","T",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<NARS\\\\>","E",pairs$pair_type);
  pairs$pair_type<-gsub("\\\\<NAHX\\\\>","T",pairs$pair_type)
  
  head(pairs)
  #####
  
  ##### Generate PCA
  #####
  colnames(pairs)
  pca_soil_pairs<-prcomp(pairs[,c(5:29)],center=T,scale.=T)
  print(pca_soil_pairs)  
  pca_soil_pairs$rotation
  
  
  soil_loading_matrix<-as.data.frame(unlist(print(pca_soil_pairs$rotation)))  
  #####
  
  ##### Loading matrix
  #####
  soil_loading_matrix
  soil_loading_matrix$variable <- rownames(soil_loading_matrix)
  head(soil_loading_matrix)
  colnames(soil_loading_matrix)
  soil_loading_matrix<-soil_loading_matrix[,c(26,1:25)]
  
  paper.loadings<-arrange(soil_loading_matrix, -abs(PC1)) %>% select(1:6)
  
  #####
  
  ##### PCA plots
  #####
  pairs$E_T_NT = droplevels(pairs$E_T_NT)
  
  pca.1<-autoplot(pca_soil_pairs,data=pairs,shape="soil_type",colour="black", fill = "pair_type",size=5 ) +
    theme_classic()+
    theme(axis.text = element_text(size=11, color="black"),plot.title = element_text(hjust = 0.5),legend.key.size = unit(1.25,"line"),legend.spacing = unit(1,"line"),legend.box.spacing = unit(1,"line"), legend.text = element_text(size = 11), legend.title = element_text(size = 12), axis.title = element_text(size = 12), legend.position = c(.35, .20), legend.box = "horizontal") +
    scale_shape_manual(values=c(21,24), name="Soil",labels=c("Nonserpentine","Serpentine"))+
  scale_fill_grey(start=.4,end=.85,name = "Pair type", labels = c("Endemic","Tolerator"))+
    guides(fill = guide_legend( override.aes = list(shape = c(15,15), color = c("gray30","gray72"), size = 5)))+
    scale_linetype_discrete(guide = FALSE)+
    geom_line(aes(group = pair_name, lty = pair_type), lwd = .5) +
    xlab("PC1 (26.62% variance explained)")+
    ylab("PC2 (19.76% variance explained)")
  pca.1
  

  PCA_GH_full_loadings<-autoplot(pca_soil_pairs,data=pairs,shape="soil_type",colour="E_T_NT",label=T,label.label="pair_name", label.repel=T,loadings=T,loadings.colour="#FC4E2A", loadings.label.colour="#FC4E2A", loadings.label=T,loadings.label.repel=T,size=4) +
    theme_classic(base_size=14)+theme(plot.title = element_text(hjust = 0.5)) +
    scale_shape_manual(values=c(9,19),name = "Soil", labels=c("NS","S"))+
    scale_colour_grey(start=.25,end=.55,name="Pair type",labels=c("Tolerator","Endemic"))
  PCA_GH_full_loadings
  #####
  
  ##### Calculate Euclidean Distance
  #####
  # Make dataframe with ind PC values
  ind_values<-as.data.frame(pca_soil_pairs$x)
  ind_values
  colnames(pairs)
  id<-pairs[,c(1:4,30:32)]
  values<-bind_cols(id,ind_values)
  values  
  colnames(values)
  head(values)
  
  # Duplicate rows of NS taxon for COSP, MGUT and NAHX pairs 
  cosp<-filter(values,sp=="COSP"&soil_type=="NS")
  cosp$pair.no<-"10"
  cosp$pair_name <- "COGR"
  cosp$pair_type <- "E"
  
  mgut<-filter(values,sp=="MGUT"&soil_type=="NS")  
  mgut$pair.no<-"13"
  mgut$pair_name <- "MNUD"
  mgut$pair_type <- "E"
  
  nahx<-filter(values,sp=="NAHX"&soil_type=="NS")  
  nahx$pair.no<- "16"
  nahx$pair_name <- "NARS"
  nahx$pair_type <- "E"
  
  values<-bind_rows(values,cosp,mgut,nahx)
  values
  
  # Create empty distance dataframe to store values  
  distances_full<-data.frame (pair.no=c(1:17), distance=NA)
  summary(pca_soil_pairs)
  
  # Loop to calculate distances for each pair
  values$pair.no <- as.numeric(values$pair.no)
  values[values$pair.no==1,]
  for (i in 1:17){
    a<-values[values$pair.no==i,]
    distances_full[i,2]<-Euc_dist_25_PCs(a)
    i=i+1
  }
  
  distances_full$pair.no<-as.character(distances_full$pair.no)
  
  # beef up distance frame
  colnames(pairs)
  ab<-pairs[pairs$soil_type=="S",c(1,2,3,4,30:32)]
  distances_full<-full_join(ab,distances_full,by="pair.no")
  distances_full
  #####
  
  ##### Graph distances (figure 3)
  #####

  head(distances_full)
  
  # graph Average distances  
  distances_full$E_T_NT<-factor(distances_full$E_T_NT,levels=c("T","E"))
  distance.box<-distances_full %>%
    ggplot(.,aes(x=E_T_NT,y=distance,fill=E_T_NT)) +
    geom_boxplot()+
    xlab("Pair Type") + 
    ylab("Euclidean distance") + 
    #ggtitle("Average pairwise distance in all-variable soil PCA")+
    theme_classic()+
    scale_fill_grey(start=.85,end=.4,name="Pair type",labels=c("Tolerator","Endemic"))+
    scale_x_discrete(labels=c("T" = "Tolerator", "E" = "Endemic")) +
    theme(axis.text = element_text(size=11, color="black"), axis.title=element_text(size=12),legend.position = "none", axis.text.x = element_text(angle = 0, vjust = .5), plot.margin = margin(.5,.5,.5,.5, "cm"))+
    scale_y_continuous(limits = c(0,10),breaks =c(0,2,4,6,8,10), expand = c(0,0))
  distance.box
  
  # Combine the PCA and euclidian distances plots (Figure 3)
  grid.arrange(pca.1,distance.box, ncol = 2, widths = c(5,2.5))
  #####
  
  ##### PGLS with geographic distance as a covariate
  #####
  soil.geo <- full_join(distances_full,geo.dist, by = "pair_name")
  
  head(soil.geo)
  soil.geo$pair_name
  row.names(soil.geo)<-soil.geo$pair_name
  name.check(soil.geo, soil.tree)

   soil.geo.gls <- gls(distance ~ E_T_NT + geo.dist, correlation = corBrownian(phy = soil.tree), data = soil.geo, method = "ML")
  anova(soil.geo.gls)

  #####
  
  #####-----------------Divergence in individual soil variables-----------------#####
  
  ##### Pairwise differences in individual harshness variables
  #####
  pairs1<-pairs
  
  head(pairs1)
  table(pairs1$pair_name)
  pairs1 %>% select(sp, pair_name, pair.no, pair_type)
  colnames(pairs1)
  
  
  cosp<-filter(pairs1,sp=="COSP"&soil_type=="NS")
  cosp$pair.no<-"10"  
  cosp$pair_name<-"COGR"
  cosp$pair_type<-"E"
  
  mgut<-filter(pairs1,sp=="MGUT"&soil_type=="NS")  
  mgut$pair.no<-"13"
  mgut$pair_name <-"MNUD"
  mgut$pair_type<-"E"
  
  
  nahx<-filter(pairs1,sp=="NAHX"&soil_type=="NS")  
  nahx$pair.no<-"16"
  nahx$pair_name <-"NARS"
  nahx$pair_type<-"E"
  
  
  pairs1<-bind_rows(pairs1,cosp,mgut,nahx)
  
  table(pairs1$pair_name, pairs1$pair_type)
  
  
  colnames(pairs1)
  soil.difference<-pairs1 %>%
    select(c(1:4,8,9,10,13,14,12,20,21,22,27,28,29,30,31,32)) %>%
    group_by(pair_name) %>%
    gather(.,"soil_var","value",5:16) %>%
    select(-c(E_T_NT, county, sp)) %>%
    group_by(pair_name, pair.no, pair_type) %>%
    spread(soil_type,value) ## This is the dataframe to use to subtract values
  
  colnames(soil.difference)
  soil.difference$full_variable<-soil.difference$soil_var
  soil.difference$full_variable<-gsub("\\\\<CaMg\\\\>","Ca:Mg",soil.difference$full_variable)
  soil.difference$full_variable<-gsub("\\\\<Ca_mg.kg\\\\>","Ca",soil.difference$full_variable)
  soil.difference$full_variable<-gsub("\\\\<Mg_mg.kg\\\\>","Mg",soil.difference$full_variable)
  soil.difference$full_variable<-gsub("\\\\<P_mg.kg\\\\>","P",soil.difference$full_variable)
  soil.difference$full_variable<-gsub("\\\\<K_mg.kg\\\\>","K",soil.difference$full_variable)
  soil.difference$full_variable<-gsub("\\\\<NH4_NP_mg.kg\\\\>","N",soil.difference$full_variable)
  soil.difference$full_variable<-gsub("\\\\<Ni_mg.kg\\\\>","Ni",soil.difference$full_variable)
  soil.difference$full_variable<-gsub("\\\\<Cr_mg.kg\\\\>","Cr",soil.difference$full_variable)
  soil.difference$full_variable<-gsub("\\\\<Co_mg.kg\\\\>","Co",soil.difference$full_variable)
  soil.difference$full_variable<-gsub("\\\\<percent_sand\\\\>","% Sand",soil.difference$full_variable) 	
  soil.difference$full_variable<-gsub("\\\\<percent_silt\\\\>","% Silt",soil.difference$full_variable)
  soil.difference$full_variable<-gsub("\\\\<percent_clay\\\\>","% Clay",soil.difference$full_variable)
  
  #change long variable and E/T to factors with right pair levels
  soil.difference$soil_var = factor(soil.difference$soil_var, levels=c("CaMg","Ca_mg.kg","Mg_mg.kg","P_mg.kg","K_mg.kg","NH4_NP_mg.kg","Ni_mg.kg","Cr_mg.kg","Co_mg.kg","percent_sand","percent_silt","percent_clay"))
  #change long variable and E/T to factors with right pair levels
  soil.difference$full_variable = factor(soil.difference$full_variable, levels=c("Ca:Mg","Ca","Mg","P","K","N","Ni","Cr","Co","% Sand","% Silt","% Clay"))
  soil.difference$pair_type<-factor(soil.difference$pair_type,levels = c("T","E"))
  #####
  
  ##### Graph of proportional pairwise differences (Appendix S3)
  #####
  soil.difference %>%
    mutate(diff = S/NS) %>%
    ggplot(.,aes(pair_type,diff, fill = pair_type))+
    geom_boxplot(width = .6)+
    facet_wrap(~full_variable,scales="free_y",nrow = 4)+
    theme_classic() +
    geom_hline(yintercept = 1, lty = 2) +
    ylab("Proportional pairwise divergence (S taxon / NS taxon)")+
    xlab("")+
    scale_x_discrete(labels = c("Tolerators","Endemics"))+
    scale_fill_grey(start=.5,end=.95,name="Species type",labels=c("Tolerator","Endemic"))+
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 20, vjust = .5, size = rel(1.4)), 
          axis.text.y = element_text(size = rel(1.25)), 
          strip.text = element_text(size = rel(1.2)), 
          axis.title.y = element_text(size = rel(1.2)),
          plot.margin = margin(t = .6, r = .6, b = .6, l = .6, unit = "cm"))
  #####
  
  #####Pairwise analysis dataframe
  #####
  head(soil.difference)
  
  prop.diff<-soil.difference %>%
    ungroup() %>%
    mutate(prop.dif = S/NS) %>%
    select(pair_name,pair.no,pair_type,soil_var,prop.dif) %>%
    spread(soil_var,prop.dif)
  
  prop.diff$pair_type<-as.factor(prop.diff$pair_type)
  str(prop.diff$pair_type)
  prop.diff<-as.data.frame(prop.diff)
  head(prop.diff)
  
  prop.geo <- full_join(prop.diff, geo.dist, by = "pair_name")
  
  prop.geo$pair_name <- factor(prop.geo$pair_name, levels = c("MGUT","CACO","PLER","COHT","TWILD","COSP","NAPB","NAHX","CABR","CAGT","COGR","NAJP","CLDV","MNUD","NARS","LADI","CABE"))
  
  #####
  
  ##### PGLS on proportional difference with geographic distance as a covariate
  #####
  
  row.names(prop.geo)<-prop.geo$pair_name
  name.check(soil.tree, data.names = prop.geo$pair_name)
  
  head(prop.geo)
  str(prop.geo$pair_type)
  pgls.dif.geo.table<-data.frame(Soil_variable=colnames(prop.geo)[4:15],
                         numDF=NA,
                         denDF=NA,
                         intercept_E_coeff=NA,
                         slope_T_coeff=NA,
                         geo.dist.effect = NA,
                         F_value_pair=NA,
                         F_value_dist = NA,
                         p_value_pair=NA,
                         p_value_dist= NA)
  
  
   prop.geo.camg<-gls(CaMg ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.camg)
  anova(prop.geo.camg)
  pgls.dif.geo.table[1,2:10]<-c(as.data.frame(anova(prop.geo.camg))[2,1], 14, as.data.frame(summary(prop.geo.camg)$tTable)[1,1], as.data.frame(summary(prop.geo.camg)$tTable)[2,1], as.data.frame(summary(prop.geo.camg)$tTable)[3,1],as.data.frame(anova(prop.geo.camg))[2,2],as.data.frame(anova(prop.geo.camg))[3,2], as.data.frame(summary(prop.geo.camg)$tTable)[2,4],as.data.frame(summary(prop.geo.camg)$tTable)[3,4])

    prop.geo.ca<-gls(Ca_mg.kg ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.ca)
  anova(prop.geo.ca)
  pgls.dif.geo.table[2,2:10]<-c(as.data.frame(anova(prop.geo.ca))[2,1], 14, as.data.frame(summary(prop.geo.ca)$tTable)[1,1], as.data.frame(summary(prop.geo.ca)$tTable)[2,1], as.data.frame(summary(prop.geo.ca)$tTable)[3,1],as.data.frame(anova(prop.geo.ca))[2,2], as.data.frame(anova(prop.geo.ca))[3,2],as.data.frame(summary(prop.geo.ca)$tTable)[2,4],as.data.frame(summary(prop.geo.ca)$tTable)[3,4])
  
  
    prop.geo.mg<-gls(Mg_mg.kg ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.mg)
  anova(prop.geo.mg)
  pgls.dif.geo.table[3,2:10]<-c(as.data.frame(anova(prop.geo.mg))[2,1], 14, as.data.frame(summary(prop.geo.mg)$tTable)[1,1], as.data.frame(summary(prop.geo.mg)$tTable)[2,1], as.data.frame(summary(prop.geo.mg)$tTable)[3,1],as.data.frame(anova(prop.geo.mg))[2,2],as.data.frame(anova(prop.geo.mg))[3,2], as.data.frame(summary(prop.geo.mg)$tTable)[2,4],as.data.frame(summary(prop.geo.mg)$tTable)[3,4])
  
  prop.geo.p<-gls(P_mg.kg ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.p)
  anova(prop.geo.p)
  pgls.dif.geo.table[4,2:10]<-c(as.data.frame(anova(prop.geo.p))[2,1], 14, as.data.frame(summary(prop.geo.p)$tTable)[1,1], as.data.frame(summary(prop.geo.p)$tTable)[2,1], as.data.frame(summary(prop.geo.p)$tTable)[3,1],as.data.frame(anova(prop.geo.p))[2,2], as.data.frame(anova(prop.geo.p))[3,2],as.data.frame(summary(prop.geo.p)$tTable)[2,4],as.data.frame(summary(prop.geo.p)$tTable)[3,4])
  
  prop.geo.k<-gls(K_mg.kg ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.k)
  anova(prop.geo.k)
  pgls.dif.geo.table[5,2:10]<-c(as.data.frame(anova(prop.geo.k))[2,1], 14, as.data.frame(summary(prop.geo.k)$tTable)[1,1], as.data.frame(summary(prop.geo.k)$tTable)[2,1], as.data.frame(summary(prop.geo.k)$tTable)[3,1],as.data.frame(anova(prop.geo.k))[2,2],as.data.frame(anova(prop.geo.k))[3,2], as.data.frame(summary(prop.geo.k)$tTable)[2,4],as.data.frame(summary(prop.geo.k)$tTable)[3,4])
  
  prop.geo.n<-gls(NH4_NP_mg.kg ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.n)
  anova(prop.geo.n)
  pgls.dif.geo.table[6,2:10]<-c(as.data.frame(anova(prop.geo.n))[2,1], 14, as.data.frame(summary(prop.geo.n)$tTable)[1,1], as.data.frame(summary(prop.geo.n)$tTable)[2,1], as.data.frame(summary(prop.geo.n)$tTable)[3,1],as.data.frame(anova(prop.geo.n))[2,2],as.data.frame(anova(prop.geo.n))[3,2], as.data.frame(summary(prop.geo.n)$tTable)[2,4],as.data.frame(summary(prop.geo.n)$tTable)[3,4])
  
  prop.geo.ni<-gls(Ni_mg.kg ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.ni)
  anova(prop.geo.ni)
  pgls.dif.geo.table[7,2:10]<-c(as.data.frame(anova(prop.geo.ni))[2,1], 14, as.data.frame(summary(prop.geo.ni)$tTable)[1,1], as.data.frame(summary(prop.geo.ni)$tTable)[2,1], as.data.frame(summary(prop.geo.ni)$tTable)[3,1],as.data.frame(anova(prop.geo.ni))[2,2],as.data.frame(anova(prop.geo.ni))[3,2], as.data.frame(summary(prop.geo.ni)$tTable)[2,4],as.data.frame(summary(prop.geo.ni)$tTable)[3,4])
  
  
  prop.geo.cr<-gls(Cr_mg.kg ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.cr)
  anova(prop.geo.cr)
  pgls.dif.geo.table[8,2:10]<-c(as.data.frame(anova(prop.geo.cr))[2,1], 14, as.data.frame(summary(prop.geo.cr)$tTable)[1,1], as.data.frame(summary(prop.geo.cr)$tTable)[2,1], as.data.frame(summary(prop.geo.cr)$tTable)[3,1],as.data.frame(anova(prop.geo.cr))[2,2],as.data.frame(anova(prop.geo.cr))[3,2], as.data.frame(summary(prop.geo.cr)$tTable)[2,4],as.data.frame(summary(prop.geo.cr)$tTable)[3,4])
  
  
  prop.geo.co<-gls(Co_mg.kg ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.co)
  anova(prop.geo.co)
  pgls.dif.geo.table[9,2:10]<-c(as.data.frame(anova(prop.geo.co))[2,1], 14, as.data.frame(summary(prop.geo.co)$tTable)[1,1], as.data.frame(summary(prop.geo.co)$tTable)[2,1], as.data.frame(summary(prop.geo.co)$tTable)[3,1],as.data.frame(anova(prop.geo.co))[2,2],as.data.frame(anova(prop.geo.co))[3,2], as.data.frame(summary(prop.geo.co)$tTable)[2,4],as.data.frame(summary(prop.geo.co)$tTable)[3,4])
  
  
  prop.geo.sand<-gls(percent_sand ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.sand)
  anova(prop.geo.sand)
  pgls.dif.geo.table[10,2:10]<-c(as.data.frame(anova(prop.geo.sand))[2,1], 14, as.data.frame(summary(prop.geo.sand)$tTable)[1,1], as.data.frame(summary(prop.geo.sand)$tTable)[2,1], as.data.frame(summary(prop.geo.sand)$tTable)[3,1],as.data.frame(anova(prop.geo.sand))[2,2],as.data.frame(anova(prop.geo.sand))[3,2], as.data.frame(summary(prop.geo.sand)$tTable)[2,4],as.data.frame(summary(prop.geo.sand)$tTable)[3,4])
  
  
  prop.geo.silt<-gls(percent_silt ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.silt)
  anova(prop.geo.silt)
  pgls.dif.geo.table[11,2:10]<-c(as.data.frame(anova(prop.geo.silt))[2,1], 14, as.data.frame(summary(prop.geo.silt)$tTable)[1,1], as.data.frame(summary(prop.geo.silt)$tTable)[2,1], as.data.frame(summary(prop.geo.silt)$tTable)[3,1],as.data.frame(anova(prop.geo.silt))[2,2],as.data.frame(anova(prop.geo.silt))[3,2], as.data.frame(summary(prop.geo.silt)$tTable)[2,4],as.data.frame(summary(prop.geo.silt)$tTable)[3,4])
  
  
  prop.geo.clay<-gls(percent_clay ~ pair_type + geo.dist,correlation = corBrownian(phy = soil.tree), data = prop.geo, method = "ML")
  summary(prop.geo.clay)
  anova(prop.geo.clay)
  pgls.dif.geo.table[12,2:10]<-c(as.data.frame(anova(prop.geo.clay))[2,1], 14, as.data.frame(summary(prop.geo.clay)$tTable)[1,1], as.data.frame(summary(prop.geo.clay)$tTable)[2,1], as.data.frame(summary(prop.geo.clay)$tTable)[3,1],as.data.frame(anova(prop.geo.clay))[2,2],as.data.frame(anova(prop.geo.clay))[3,2], as.data.frame(summary(prop.geo.clay)$tTable)[2,4],as.data.frame(summary(prop.geo.clay)$tTable)[3,4])
  
  #####
  
  ##### PGLS proportional difference summary table
  #####
  pgls.dif.geo.table
  pgls.dif.geo.table<-arrange(pgls.dif.geo.table, p_value_pair)
  pgls.dif.geo.table$bon.p <- bon
  pgls.dif.geo.table
  #####
  
#####-------------------------------######  
#####-------------------------------######
#####  Habitat analyses: Do serpentine endemic sister taxa pairs have more divergence in competitive environment than serpentine tolerator sister taxa pairs? 2 - step model.
#####-------------------------------###### 
  
  ##### Coancestry matrix
  #####
  G
  Ginv
  #####
  
  ##### Set up the data frame
  #####
  head(comp)
  cosp <- comp %>% filter(sp == "COSP" & soil == "NS")
  cosp$pair_name <- "COGR"
  cosp$pair_type <- "E"

  mgut <- comp %>% filter(sp == "MGUT" & soil == "NS")
  mgut$pair_name <- "MNUD"
  mgut$pair_type <- "E"

  comp <- bind_rows(comp, cosp)
  comp <- bind_rows( comp, mgut)
  head(comp)
  
  comp$E_T<-comp$pair_type
  comp$E_T<-gsub("\\\\<E\\\\>",1,comp$E_T)
  comp$E_T<-gsub("\\\\<T\\\\>",0,comp$E_T)
  comp$E_T <- as.numeric(comp$E_T)
  
  comp$soil.type<-comp$soil
  comp$soil.type<-gsub("\\\\<S\\\\>",1,comp$soil.type)
  comp$soil.type<-gsub("\\\\<NS\\\\>",2,comp$soil.type)
  comp$soil.type <- as.numeric(comp$soil.type)
  
  rownames(Ginv)
  comp$pair.no<-comp$pair_name
  comp$pair.no<-gsub("\\\\<CABE\\\\>","1",comp$pair.no)
  comp$pair.no<-gsub("\\\\<CABR\\\\>","2",comp$pair.no)
  comp$pair.no<-gsub("\\\\<CACO\\\\>","3",comp$pair.no)
  comp$pair.no<-gsub("\\\\<CAGT\\\\>","4",comp$pair.no)
  comp$pair.no<-gsub("\\\\<CLDV\\\\>","6",comp$pair.no)
  comp$pair.no<-gsub("\\\\<COGR\\\\>","10",comp$pair.no)
  comp$pair.no<-gsub("\\\\<COHT\\\\>","11",comp$pair.no)
  comp$pair.no<-gsub("\\\\<COSP\\\\>","12",comp$pair.no)
  comp$pair.no<-gsub("\\\\<LADI\\\\>","9",comp$pair.no)
  comp$pair.no<-gsub("\\\\<MGUT\\\\>","14",comp$pair.no)
  comp$pair.no<-gsub("\\\\<MNUD\\\\>","13",comp$pair.no)
  comp$pair.no<-gsub("\\\\<NAJP\\\\>","7",comp$pair.no)
  comp$pair.no<-gsub("\\\\<NAPB\\\\>","8",comp$pair.no)
  comp$pair.no<-gsub("\\\\<PLER\\\\>","15",comp$pair.no)
  comp$pair.no<-gsub("\\\\<TWILD\\\\>","5",comp$pair.no)
  comp$pair.no<-as.numeric(comp$pair.no)
  
  
  comp %>%
    group_by(pair_name) %>%
    filter(quad==1) %>%
    dplyr::summarise(n())
  
  et.pair<-comp[comp$quad==1&comp$soil=="S",c("E_T","pair.no")] %>% arrange(.,pair.no)
  et.pair
  
  head(comp)
  
  #####
  
  ##### Make JAGS model 1 data
  #####
  head(comp)
  paired.m1 = list(
    pair.no = as.double(comp$pair.no),
    soil.type = as.double(comp$soil.type),
    n = nrow(comp),
    y = as.double(comp$bare),
    n.pts = as.double(comp$sum)
  )
  #####
  
  ##### JAGS model 1
  #####
  sink("paired_bare_m1.R")
  cat("
      model{
      
      # phi prior for all pairs
      for (k in 1:15){
      for (j in 1:2){
      phi[j,k] ~ dbeta(0.001, 0.001)
      }
      }
      
      # every taxon per pair gets its own phi
      for (i in 1:n){
      y[i] ~ dbin(phi[soil.type[i], pair.no[i]],n.pts[i])
      }
      
      # derived value - pair.dif (phi.S - phi.NS)
      for(k in 1:15){
      pair.dif[k] <- phi[1,k] - phi[2,k]
      }
      }
      ",fill = TRUE)
  sink()
  
  n.adapt = 5000
  n.update = 10000
  n.iter = 20000
  #####
  
  ##### Run JAGS model 1
  #####
  m1 = jags.model("paired_bare_m1.R", data = paired.m1, n.chains = 3, n.adapt = n.adapt)
  update(m1, n.iter = n.update)
  m1c = coda.samples(m1, variable.names = c("phi","pair.dif"), n.iter = n.iter, n.thin = 1)
  m1j = jags.samples(m1, variable.names = c("phi","pair.dif"), n.iter = n.iter, n.thin = 1)
  #####
  
    ##### diagnostics model 1
    #####
    gelman.diag(m1c, multivariate = F)
    #plot(m1c)
    #####
  
    ##### summary of all three chains model 1
    #####
    model1 <- as.data.frame(rbind(m1c[[1]], m1c[[2]], m1c[[3]]))
    colnames(model1)
    
    model1$E_mean <- ((model1[,1] + model1[,4] + model1[,6] + model1[,7] + model1[,9] + model1[,10] + model1[,13]) / 7)
   
    model1$T_mean <- ((model1[,2] + model1[,3] + model1[,5] + model1[,8] + model1[,11] + model1[,12] + model1[,14] + model1[,15]) / 8)

    bare_div_pairmeans<-model1 %>%
      select(E_mean, T_mean) %>%
      gather("pairtype", "mean") %>%
      group_by(pairtype) %>%
      summarise(avg = mean(mean), Q25 = quantile(mean, .025), Q50 = quantile(mean,.975))
    bare_div_pairmeans$y <- c(14,14)
    Emean <- bare_div_pairmeans %>% filter (pairtype == "E.mean")
    Tmean <- bare_div_pairmeans %>% filter (pairtype == "T.mean")  
    #####
  
    ##### All graphs on same axis (Figure 4)
    #####
    par(mfrow = c(1,1))
    par(mar=c(3.5,3.5,1,1.5))
    plot(NA, xlab = "", xlim = c(-1,1), ylim = c(0,15),bty="n", ylab = "", cex.axis = .8, cex.lab = 2)
    abline(v = 0, lty = 4, lwd = 1.75, col = "gray47")
    mtext(side = 1, line = 2.3, expression(paste("Divergence in bare ground (",phi[S] - phi[NS],")")))
    mtext(side = 2, line = 2, "Probability density")
    for (i in c(2,3,5,8,11,12,14,15)){
      lines(density(model1[,1 + (i-1)]), col = "orange", lwd = 1.75)
    }
    for (i in c(1,4,6,7,9,10,13)){
      lines(density(model1[,1 + (i-1)]), col = "blue", lwd = 1.75)
    }
    legend(-1.05, 14, legend = c("Tolerator","Endemic"), col = c("orange","blue"), lwd = c(2,2), cex = .8) 
    # arrows( Emean$Q25, Emean$y, Emean$Q50,Emean$y,  length = 0.06, code = 3, angle = 90, col = "blue", lwd = 2)
    # arrows( Tmean$Q25, Tmean$y, Tmean$Q50,Tmean$y,  length = 0.06, code = 3, angle = 90, col = "orange", lwd = 2)
    points(data = filter(bare_div_pairmeans, pairtype == "E_mean"), y~avg,  type = "p", col = "black", pch = 23, cex = 1, bg = "blue") 
    points(data = filter(bare_div_pairmeans, pairtype == "T_mean"), y~avg,  type = "p", col = "black", pch = 23, cex = 1, bg = "orange")
    
    #####
    
    ##### Facet graph of divergence in the tolerator pairs (Appendix S6)
    #####
    par(mfrow = c(4,2))
    hist(model1[,2], freq = F, xlab = " ", main = paste(comp[comp$pair.no==2 & comp$quad ==1 & comp$soil=="S","pair_name"]), xlim = c(-1,1), cex.lab = 1.5, cex.main = 1.5, ylab = " ")
    lines(density(model1[,2]), col = "orange", lwd = 2)
    abline(v = 0)
    mtext(side = 1, line = 2.5, expression(paste("Divergence in bare ground (",phi[S] - phi[NS],")")))
    mtext(side = 2, line = 2.5, "Probability density")
    for (i in c(3,5,8,11,12,14,15)){
      hist(model1[,1 +(i-1)], freq = F, xlab = " ", main = paste(comp[comp$pair.no==i & comp$quad ==1 & comp$soil=="S","pair_name"]),xlim = c(-1,1),cex.lab = 1.5, cex.main = 1.5,ylab = " ")
      lines(density(model1[,1 + (i-1)]), col = "orange", lwd = 2)
      abline(v = 0)
      mtext(side = 1, line = 2.5, expression(paste("Divergence in bare ground (",phi[S] - phi[NS],")")))
      mtext(side = 2, line = 2.5, "Probability density")
    }
    #####
    
    ##### Facet graph of divergence in the endemic pairs (Appendix S7)
    #####
    par(mfrow = c(4,2))
    par(mar=c(4,4,3,2))
    colnames(model1)
    hist(model1[,1], freq = F, xlab = " ", main = paste(comp[comp$pair.no==1 & comp$quad ==1 & comp$soil=="S","pair_name"]), xlim = c(-1,1), cex.lab = 1.5, cex.main = 1.5, ylab = " ")
    lines(density(model1[,1]), col = "blue", lwd = 2)
    mtext(side = 1, line = 2.5, expression(paste("Divergence in bare ground (",phi[S] - phi[NS],")")))
    mtext(side = 2, line = 2.5, "Probability density")
    abline(v = 0)
    for (i in c(4,6,7,9,10,13)){
      hist(model1[,1 +(i-1)], freq = F, xlab = " ", main = paste(comp[comp$pair.no==i & comp$quad ==1 & comp$soil=="S","pair_name"]),xlim = c(-1,1),cex.lab = 1.5, cex.main = 1.5,ylab = " ")
      lines(density(model1[,1 + (i-1)]), col = "blue", lwd = 2)
      abline(v = 0)
      mtext(side = 1, line = 2.5, expression(paste("Divergence in bare ground (",phi[S] - phi[NS],")")))
      mtext(side = 2, line = 2.5, "Probability density")
    }
    
    
    #####
    
##### data for JAGS model 2
#####

head(model1)
mu.var <- model1 %>%
  select(1:15) %>%
  gather("pair","mcmc") %>%
  group_by(pair) %>%
  summarise(mu = mean(mcmc), var = var(mcmc)) %>% as.data.frame()

mu.var$pair.no = c(1,10,11,12,13,14,15,2,3,4,5,6,7,8,9)
mu.var <- arrange(mu.var,pair.no)    
mu.var
#calculate tau
mu.var<-mu.var %>%
  mutate(tau = 1/var)

mu.var <- comp %>% dplyr::select(pair_name,pair.no) %>% group_by(pair_name, pair.no) %>% summarise(n = n()) %>% dplyr::select(-n) %>% right_join(.,mu.var)
mu.var <- left_join(mu.var, geo.dist, by = "pair_name") %>% as.data.frame()
mu.var


data.paired.m2 = list(
  mu.dif = as.double(mu.var$mu),
  tau.dif = as.double (mu.var$tau),
  mu.0 = as.double(rep(0,15)),
  Ginv = Ginv,
  E_T = as.double(et.pair$E_T),
  geo.dist = as.double(mu.var$geo.dist)
)
#####
  
##### JAGS model 2 with geographic distance as covariate
#####
sink("paired_bare_m2.R")
cat("
    model{
    
    # data structure - every k is a row - columns mean of pair.dif [k] (=mu.dif[k]) and sd of pair.dif[k] (=sd.dif[k]), and tau for each pair (=1/sd)
  
    for (k in 1:15){
    mu.dif[k] ~ dnorm(mu[k], tau.dif[k])
    
    mu[k] ~ dnorm(alpha[k], tau.mu)

    alpha[k] <- b0[k] + b1 * E_T[k] + b2 * geo.dist[k]
    #alpha[k] <- b1 * E_T[k] + b2 * geo.dist[k]
    }

    # Priors on tau.mu (process error)
    tau.mu <- 1/sigma.error
    sigma.error ~ dunif(0,100)
    
    # Priors on regression coefficients
    b0 ~ dmnorm(mu.0,Ginv)
    b1 ~ dnorm(0,tau.b1)
    tau.b1 <- 1 / sigma2b
    sigma2b ~ dunif(0,100)
    b2 ~ dnorm(0,tau.b2)
    tau.b2 <- 1 / sigma2b2
    sigma2b2 ~ dunif(0,100)
    }
    ",fill = TRUE)
sink()

n.adapt = 50000
n.update = 100000
n.iter = 100000

#####

##### Run JAGS model 2 
#####
m2 = jags.model("paired_bare_m2.R", data = data.paired.m2, n.chains = 3, n.adapt = n.adapt)
update(m2, n.iter = n.update)
m2c = coda.samples(m2, variable.names = c("b1","mu","alpha","b2","b0"), n.iter = n.iter, n.thin = 1)
#m2j = jags.samples(m2, variable.names = c("b1","mu","alpha","b2","b0"), n.iter = n.iter, n.thin = 1)
#####
    
    ##### diagnostics
    #####
    gelman.diag(m2c, multivariate = F)
    #plot(m2c)
    #####
  
    ##### summary of all three chains
    #####
    model2 <- as.data.frame(rbind(m2c[[1]], m2c[[2]], m2c[[3]]))
    head(model2)
    #####
    
    ##### b1 posterior
    #####
    colnames(model2)
    
    par(mfrow = c(1,1))
    hist(model2[,31], freq = FALSE, main = "B1 coefficient", xlab = "b1", ylim = c(0,2))
    lines(density(model2[,31]), col = "blue", lwd = 2)
    abline(v = quantile(model2[,31],c(0.0275, .975)), col = "red", lty = 3, lwd = 2)
    #####
    
    ##### % of the b1 distribution over 0
    1-ecdf(model2$b1)(0)
    #####
    
    ##### b2 posterior - effect of geographic distance
    #####
    colnames(model2)
    
    par(mfrow = c(1,1))
    hist(model2[,32], freq = FALSE, main = "B2 coefficient", xlab = "b2")
    lines(density(model2[,32]), col = "blue", lwd = 2)
    abline(v = quantile(model2[,32],c(0.0275, .975)), col = "red", lty = 3, lwd = 2)
    
    quantile(model2[,32],c(0.0275, .5, .975))
    1-ecdf(model2$b2)(0)
    median(model2[,32])
    mean(model2[,32])
    #####
  
  
#####-------------------------------######  
#####-------------------------------######
#####  Habitat analyses: Is soil chemistry/texture correlated with microhabitat bareness?
#####-------------------------------###### 
    
  
  ## JUST SERPENTINE TAXA

  ##### make S only df
  #####
  # data frame of median phi values for each serpentine population   
  med.phi<-data.frame(sp.no = key$pair.no,
                      sp = key$sp,
                      phi.med = NA)
  for (i in 1:15){
    med.phi[i,3] <-quantile(phi.sp1[,17 + (i-1)],.5)[[1]]
  }
  med.phi
  
  # Soil data frame for all serpentine populations
  
  head(pairs) # soil data for all populations
  
  S.soil <- pairs %>% filter(soil_type=="S")
  colnames(S.soil)
  S.soil$pair.no <- as.numeric(S.soil$pair.no)
  S.soil %>% select(sp, pair_name, pair.no) %>% arrange(pair.no) 
  
  # combine bare and soil NS dataframes
  colnames(med.phi)
  colnames(S.soil)
  med.phi$sp

  S.soil$sp 
  
  S.soil <- S.soil %>% filter(., pair_name !="NAHX") %>% filter(.,pair_name!="NARS")   # remove the NS.soil NAHX
  
  S.corr <-full_join(med.phi, S.soil, by = c("sp" = "pair_name")) %>% as.data.frame()
  #####
  
  ##### S only corr coef w/ all variables
  #####
  S.corr
  colnames(S.corr)
  chart.Correlation(S.corr[,c(3,9:33)], histogram= T, pch = 19) # numbers that it gives you are pearson correlation coefficients   
  
  bar.cor.S<-rcorr(as.matrix(S.corr[,c(3,8:32)]), type = "pearson")

  
  
  bar.coef.S<-data.frame(variable = row.names(bar.cor.S$r),
                       r = NA,
                       p = NA)
  
  for(i in 1:26){
  bar.coef.S$r[i] <- bar.cor.S$r[1,][[i]]
  bar.coef.S$p[i] <- bar.cor.S$P[1,][[i]]
  }
  bar.coef.S<-arrange(bar.coef.S, p)
  bar.coef.S
  adj.p <- rep(0.05,26)
  for (i in 1:26){
    adj.p[i] <- round(.05/(26-(i-1)),4)
  }
  adj.p
  bar.coef.S$adj.p <- adj.p
  bar.coef.S

  #####
  
  ## JUST NONSERPENTINE TAXA
  
  ##### make NS only df
  #####
  colnames(model1)
  # S =  1; NS = 2
  
  #extract all NS pops from model1
  NS.phi <- model1 %>% select(seq(17,45,2)) %>% gather("pair", "phi", 1:15) %>% group_by(pair) %>% summarise(med.phi = median(phi))
  NS.phi$soil <- "NS"  
  NS.phi
  NS.phi$pair.no <- c(1,10,11,12,13,14,15,2,3,4,5,6,7,8,9)  
  
  # mu.var goes w/ model1
  mu.var %>% select(pair_name, pair.no)
  
  NS.phi <- mu.var %>% select(pair_name, pair.no) %>% full_join(.,NS.phi)
  
  # make soil df for NS
  head(pairs) # soil data for all populations
  
  NS.soil <- pairs %>% filter(soil_type=="NS")
  colnames(NS.soil)
  NS.soil$pair.no <- as.numeric(NS.soil$pair.no)
  NS.soil %>% select(sp, pair_name, pair.no) %>% arrange(pair.no) # onyl 14 bc mgut/nahx/cosp aren't duplicated
    # pair numbers in this dataframe are totally differet
  
  # combine bare and soil NS dataframes
  colnames(NS.phi)
  colnames(NS.soil)
  NS.phi$pair_name # this doesn't have the NARS/NAHX pair, but has duplicated COSP/MGUT/NAHX
  NS.soil$pair_name # this doesn't have the COGR, MNUD or NARS pairs (bc isnt duplicated)
  
  # remove the NS.phi COGR and MNUD
  NS.phi <-NS.phi %>% filter(., pair_name!="MNUD") %>% filter(.,pair_name!="COGR")
    # remove the NS.soil NAHX
  NS.soil <- NS.soil %>% filter(., pair_name !="NAHX")
  
  NS.corr <-full_join(NS.phi, NS.soil, by = "pair_name") %>% as.data.frame()
  #####
  
  ##### NS only corr coef w/ all variables
  #####
  NS.corr
  colnames(NS.corr)
  chart.Correlation(NS.corr[,c(4,10:34)], histogram= T, pch = 19) # numbers that it gives you are pearson correlation coefficients   
  
  bar.cor.NS<-rcorr(as.matrix(NS.corr[,c(4,10:34)]), type = "pearson")
  str(bar.cor.NS)
  
  
  bar.coef.NS<-data.frame(variable = row.names(bar.cor.NS$r),
                       r = NA,
                       p = NA)
  
  for(i in 1:26){
  bar.coef.NS$r[i] <- bar.cor.NS$r[1,][[i]]
  bar.coef.NS$p[i] <- bar.cor.NS$P[1,][[i]]
  }
  bar.coef.NS<-arrange(bar.coef.NS, p)
  bar.coef.NS
  adj.p <- rep(0.05,26)
  for (i in 1:26){
    adj.p[i] <- round(.05/(26-(i-1)),4)
  }
  adj.p
  bar.coef.NS$adj.p <- adj.p
  bar.coef.NS

  #####
    
  ## All TAXA
  
  ##### make combined df
  #####
  NS.corr
  S.corr
  colnames(S.corr)
  S.corr1 <- S.corr %>% select(sp, phi.med,8:32)
  colnames(NS.corr)
  NS.corr1 <- NS.corr %>% select(sp,med.phi,10:34 )
  
  colnames(S.corr1) [2] <- "med.phi"
  colnames(NS.corr1)
  
  all.taxa.corr <- bind_rows(S.corr1, NS.corr1)
  #####
  
  ##### Corr coefs for all taxa phi (% bare ground) vs all soil variation
  #####
  all.taxa.corr
  colnames(all.taxa.corr)
  #chart.Correlation(all.taxa.corr[,c(3,9:33)], histogram= T, pch = 19) # numbers that it gives you are pearson correlation coefficients   
  
  bar.cor.ALL<-rcorr(as.matrix(all.taxa.corr[,c(2:27)]), type = "pearson")

  
  
  bar.coef.ALL<-data.frame(variable = row.names(bar.cor.ALL$r),
                       r = NA,
                       p = NA)
  nrow(bar.cor.ALL$r)
  for(i in 1:26){
  bar.coef.ALL$r[i] <- bar.cor.ALL$r[1,][[i]]
  bar.coef.ALL$p[i] <- bar.cor.ALL$P[1,][[i]]
  }
  bar.coef.ALL<-arrange(bar.coef.ALL, p)
  bar.coef.ALL
  adj.p <- rep(0.05,26)
  for (i in 1:26){
    adj.p[i] <- round(.05/(26-(i-1)),4)
  }
  adj.p
  bar.coef.ALL$adj.p <- adj.p
  bar.coef.ALL
  #####
  