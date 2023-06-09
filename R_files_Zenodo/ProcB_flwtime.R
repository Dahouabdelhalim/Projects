#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## Parallel evolution of phenological isolation across the speciation spectrum in serpentine-adapted plants. 
#####---------------------#####-----------------------######

## README - NOTES
  # all treatments are written in the form of seed_soil: S_S indicates serpentine seed in serpentine soil

  # to aid in readability of this script, any time you see a heading that looks like:
        #  ##### "heading title"
        #  #####
   # you can click on the down-facing arrow next to the line numbers on the left-hand side of the RStudio-script-window. Clicking this down-facing arrow will collapse the code in that given section. Clicking on the arrow again will expand the code.  


##### libraries
#####
library(rjags)
library(phylotate)
library(ape)
library(phytools)
library(mvtnorm)
library(geiger)
library(nlme)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggrepel)
library(ggplot2)
library(ggfortify)
library(wesanderson)
library(broom)
library(rr2)
library(extrafont)
font_import()
#loadfonts
loadfonts(device = "pdf")
#####

##### functions
#####
sem<-function(x){sd(x,na.rm=T)/sqrt(length(x))}

Euc_dist_19_PCs.clim<-function(x){
    sqrt(
      abs(x[x$soil=="S","PC1"]-x[x$soil=="NS","PC1"])^2 +
        abs(x[x$soil=="S","PC2"]-x[x$soil=="NS","PC2"])^2 +
        abs(x[x$soil=="S","PC3"]-x[x$soil=="NS","PC3"])^2 +
        abs(x[x$soil=="S","PC4"]-x[x$soil=="NS","PC4"])^2 +
        abs(x[x$soil=="S","PC5"]-x[x$soil=="NS","PC5"])^2 +
        abs(x[x$soil=="S","PC6"]-x[x$soil=="NS","PC6"])^2 +
        abs(x[x$soil=="S","PC7"]-x[x$soil=="NS","PC7"])^2 +
        abs(x[x$soil=="S","PC8"]-x[x$soil=="NS","PC8"])^2 +
        abs(x[x$soil=="S","PC9"]-x[x$soil=="NS","PC9"])^2 +
        abs(x[x$soil=="S","PC10"]-x[x$soil=="NS","PC10"])^2 +
        abs(x[x$soil=="S","PC11"]-x[x$soil=="NS","PC11"])^2 +
        abs(x[x$soil=="S","PC12"]-x[x$soil=="NS","PC12"])^2 +
        abs(x[x$soil=="S","PC13"]-x[x$soil=="NS","PC13"])^2 +
        abs(x[x$soil=="S","PC14"]-x[x$soil=="NS","PC14"])^2 +
        abs(x[x$soil=="S","PC15"]-x[x$soil=="NS","PC15"])^2 +
        abs(x[x$soil=="S","PC16"]-x[x$soil=="NS","PC16"])^2 +
        abs(x[x$soil=="S","PC17"]-x[x$soil=="NS","PC17"])^2 +
        abs(x[x$soil=="S","PC18"]-x[x$soil=="NS","PC18"])^2 +
        abs(x[x$soil=="S","PC19"]-x[x$soil=="NS","PC19"])^2 
    )
  }  
#####

getwd()

setwd("your_directory")

#### input dataframes  
#####


  # ITS genetic distance
    its <- read.csv("ITS_genetic_distance.csv")
  
  # Phylogenetic tree, from Sianta and Kay, 2019
    serp.tree <- read.tree("phylo_tree.tre")
    plot(serp.tree)
    
  # multivariate soil chemistry, from Sianta and Kay 2019
    soil.distance <- read.csv("multivariate_soil_distance.csv")
    
  # GPS points for all pops used to extract climatic data
    points <- read.csv("TAXA_GPS_google_earth.csv") 
    
##### 

##### greenhouse dataframe "flw" 
#####
    master1 <- read.csv("greenhouse_data.csv")
    master1$pair_name <- factor(master1$pair_name, levels = c("MGUT","CABR","CACO","COHT","COSP","NAHX","NAPB","PLER","TWILD","CABE_CAST", "CAGT_CAGA", "CLDV_CLHT" ,"COGR_COSP" ,"LADI_LAGL","MNUD_MGUT" ,"NAJP_NAHN" ,"NARS_NAHX"))
    master1$pair_type <- as.character(master1$pair_type)
    
    # tolerator NS taxa that serve as the NS sister for both an endemic and tolerator pair. Duplicate shared treatment (the nonserpentine seeds in nonserpentine soil) from the endemic pair, and give it the correct information for the corresponding tolerator pair. Then concatentate all samples together
    mgut<-master1 %>% filter(sp == "MGUT" & trt=="NS_NS")
    mgut[,"pair_name"] <- "MGUT"
    mgut[,"pair"] <- 14
    mgut[,"pair_type"] <- "T"
    
    cosp<-master1 %>% filter(sp == "COSP" & trt=="NS_NS")
    cosp[,"pair_name"] <- "COSP"
    cosp[,"pair"] <- 12
    cosp[,"pair_type"] <- "T"
   
    nahx<-master1 %>% filter(sp == "NAHX"& trt=="NS_NS")
    nahx[,"pair_name"] <- "NAHX"
    nahx[,"pair"] <- 17
    nahx[,"pair_type"] <- "T"
    
    flw<-rbind(master1,mgut)
    flw <- rbind(flw,cosp)
    flw <- rbind(flw,nahx)
    
    flw<-flw %>% filter(!is.na(days_to_flw)) # subset out all individuals that did not survive to flower
    flw
####


     
#####

    
### per pair (rows) and treatment (seed_soil, columns), sampe size of individuals that survived to flower
flw %>%
  group_by(pair_name, trt) %>%
  summarise(n = n()) %>%
  spread(trt, n)
    

#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## ITS genetic distances and phylogenetic tree
#####---------------------#####-----------------------######

# mean and standard deviation
its %>%
  group_by(pair_type) %>%
  summarise(avg.nt.diff = mean(nucleotide_differences),
            sd.nt.diff = sd(nucleotide_differences))


#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## Are flowering onset shifts common following adaptation to serpentine? and are they correlated with greater shifts in the edaphic and/or climatic environment?
#####---------------------#####-----------------------######

##### Average magnitude of flowering time shift across all pairs
#####
flw %>%
  filter(trt=="S_S"|trt=="NS_NS") %>%
  group_by(pair_name, trt) %>%
  summarise(avg.flw = mean(days_to_flw)) %>%
  spread(trt, avg.flw ) %>% mutate(shift = S_S-NS_NS) %>%
  filter(shift > 0) %>% ungroup() %>% summarise(avg = mean(shift)) # mean shift for when S later than NS
    
flw %>%
  filter(trt=="S_S"|trt=="NS_NS") %>%
  group_by(pair_name, trt) %>%
  summarise(avg.flw = mean(days_to_flw)) %>%
  spread(trt, avg.flw ) %>% 
  mutate(shift = abs(S_S-NS_NS)) %>% ungroup() %>%summarise(avg = mean(shift), sem = sem(shift))
#####

##### T.test for every pair (S_S vs NS_NS), Table S2
#####

# base data frame
tab.s1<-flw %>%
  filter(trt=="S_S"|trt=="NS_NS") %>%
  group_by(pair_name, trt) %>%
  summarise(avg.flw = mean(days_to_flw)) %>%
  spread(trt, avg.flw ) %>% # each column average flowering onset for each taxon in soil soil
  mutate(dif = S_S - NS_NS) %>% as.data.frame() # difference in the averages within each pair

# attach t.test results for all pairs

tab.s1$t <- NA
tab.s1$df <- NA
tab.s1$p.val <- NA

tab.s1[1,5:7] <- c(round(t.test(data = filter(flw, pair_name == "MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[2,5:7] <- c(round(t.test(data = filter(flw, pair_name == "CABR" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "CABR" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "CABR" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[3,5:7] <- c(round(t.test(data = filter(flw, pair_name == "CACO" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "CACO" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "CACO" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[4,5:7] <- c(round(t.test(data = filter(flw, pair_name == "COHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "COHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "COHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[5,5:7] <- c(round(t.test(data = filter(flw, pair_name == "COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[6,5:7] <- c(round(t.test(data = filter(flw, pair_name == "NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[7,5:7] <- c(round(t.test(data = filter(flw, pair_name == "NAPB" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "NAPB" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "NAPB" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[8,5:7] <- c(round(t.test(data = filter(flw, pair_name == "PLER" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "PLER" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "PLER" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[9,5:7] <- c(round(t.test(data = filter(flw, pair_name == "TWILD" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "TWILD" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "TWILD" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[10,5:7] <- c(round(t.test(data = filter(flw, pair_name == "CABE_CAST" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "CABE_CAST" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "CABE_CAST" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[11,5:7] <- c(round(t.test(data = filter(flw, pair_name == "CAGT_CAGA" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "CAGT_CAGA" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "CAGT_CAGA" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[12,5:7] <- c(round(t.test(data = filter(flw, pair_name == "CLDV_CLHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "CLDV_CLHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "CLDV_CLHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[13,5:7] <- c(round(t.test(data = filter(flw, pair_name == "COGR_COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "COGR_COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "COGR_COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[14,5:7] <- c(round(t.test(data = filter(flw, pair_name == "LADI_LAGL" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "LADI_LAGL" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "LADI_LAGL" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[15,5:7] <- c(round(t.test(data = filter(flw, pair_name == "MNUD_MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "MNUD_MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "MNUD_MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))

tab.s1[16,5:7] <- c(round(t.test(data = filter(flw, pair_name == "NAJP_NAHN" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "NAJP_NAHN" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "NAJP_NAHN" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))


tab.s1[17,5:7] <- c(round(t.test(data = filter(flw, pair_name == "NARS_NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "NARS_NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "NARS_NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw ~ trt)[[3]][[1]], 10))


  tab.s1
  tab.s1<- tab.s1 %>% arrange(p.val)
  
    # sequential bonferonni correction
  bon <- rep(0,17)
  for ( i in 1:17){
    bon[i]<- round(0.05/(17 - i + 1),5)
  }
  bon
  tab.s1$bon <- bon

  tab.s1
#####
  
  
  ## Note - to perform analysis of absolute shifts in flowering time vs multivariate soil/climatic difference, need the mean absolute shift from the posterior distribution of Bayesian Model 1. Run Bayesian Models 1 &2 in next section, then return to the relationship b/t the magnitude of onset shifts and edaphic/climatic divergence. 
  
#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## Do endemic sister taxa pairs have greater flowering onset shifts than tolerator sister taxa pairs?
#####---------------------#####-----------------------######  
  
##### Set up the coancestry matrix
#####
plot(serp.tree)
G1 <- vcv(serp.tree, corr = T)
G1
Ginv1<-solve(G1)
Ginv1
##### 
  
##### Make JAGS data frame
#####
# Make new columns w/ numeric values for pair type, soil and seed
    flw$E_T<-flw$pair_type
    flw$E_T<-gsub("\\\\<E\\\\>",1,flw$E_T)
    flw$E_T<-gsub("\\\\<T\\\\>",0,flw$E_T)
    flw$E_T <- as.numeric(flw$E_T)
      
    flw$soil.type<-flw$soil
    flw$soil.type<-gsub("\\\\<S\\\\>",1,flw$soil.type)
    flw$soil.type<-gsub("\\\\<NS\\\\>",2,flw$soil.type)
    flw$soil.type <- as.numeric(flw$soil.type)
    
    flw$seed.type<-flw$seed
    flw$seed.type<-gsub("\\\\<S\\\\>",1,flw$seed.type)
    flw$seed.type<-gsub("\\\\<NS\\\\>",2,flw$seed.type)
    flw$seed.type <- as.numeric(flw$seed.type)

# JAGS dataframe
    data = list(
      pair.no = as.double(flw$pair),
      soil.type = as.double(flw$soil.type),
      seed.type = as.double(flw$seed.type),
      n = nrow(flw),
      y = as.double(flw$days_to_flw)
    )
#####
  
##### JAGS model 1
#####
sink("flw_onset_m1.R")
cat("
    model{
    
    # lambda (mean of poisson distribution) prior for all pairs
    for (k in 1:17){
    for (m in 1:2){
    for (j in 1:2){
    lambda[j,m,k] ~ dgamma(0.001, 0.001)
    }
    }
    }

    # every pair/trt gets its own lambda
    for (i in 1:n){
    y[i] ~ dpois(lambda[seed.type[i], soil.type[i], pair.no[i]])
    }
    
    # derived value - adj.dif = S_S - NS_NS
    for (k in 1:17){
    adj.dif[k] <- abs(lambda[1,1,k] - lambda[2,2,k])
    adj.dif.dir[k] <- lambda[1,1,k] - lambda[2,2,k]
    }
    
    }
    ",fill = TRUE)
sink()

n.adapt = 5000
n.update = 10000
n.iter = 10000
#####
  
##### Run JAGS model 1
#####
m1 = jags.model("flw_onset_m1.R", data = data, n.chains = 3, n.adapt = n.adapt)
update(m1, n.iter = n.update)
m1c = coda.samples(m1, variable.names = c("lambda","adj.dif", "adj.dif.dir"), n.iter = n.iter, n.thin = 1)
#####
  
    ##### diagnostics model 1 - lambda subsets are [seed, soil, pair] - lambdas w/ convergence diagnostics > 1 are for 2 nonserpentine seed in serpentine soil treatments, which are not used the calculation of flw onset shifts. 
    #####
    gelman.diag(m1c, multivariate = F)
    #plot(m1c)
    #####
    
    ##### summary of all three chains model 1
    #####
    model1 <- as.data.frame(rbind(m1c[[1]], m1c[[2]], m1c[[3]]))
    dim(model1)
    #####

    ##### Figure 2 - flowering time divergence
    #####
    colnames(model1)
    fig2<-model1 %>%
    dplyr::select(18:34) %>%
    gather("pair","adj.dif.dir",1:17) %>%
    group_by(pair) %>%
    summarise(avg = mean(adj.dif.dir),
              low.95 = quantile(adj.dif.dir,0.025),
              high.95 = quantile(adj.dif.dir, 0.975)) 
    fig2$pair <- c(1, 10, 11, 12, 13, 14, 15, 16, 17, 2, 3, 4, 5, 6, 7, 8, 9)
    fig2 <- flw %>% group_by(pair_name, pair_type, pair) %>% summarise(n = n()) %>% dplyr::select(-n) %>% full_join(.,fig2, by = "pair")
    
  fig2 %>%  
  ggplot(.,aes(pair_name, avg, fill = pair_type))+
  geom_errorbar(aes(ymin = low.95, ymax = high.95), width = .4)+
    geom_point(size = 3.5, shape = 21) +
geom_hline(yintercept = 0, lty = 2) +
theme_minimal() + 
  scale_y_continuous(breaks = seq(-40,40,10)) +
  scale_fill_manual(values = c("#018571","#a6611a"), name = "Pair type", labels = c("Endemic","Tolerator"))+
  scale_color_manual(values = c("#018571","#a6611a"), name = "Pair type", labels = c("Endemic","Tolerator"))+
  theme(#axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  xlab("Pair name") + ylab("Average flowering onset shift (days) \\n ")
    
    
    #####
    
    
    
    
##### Summarise mean and variance of posteriors of the onset shifts from model 1 ("adj.dif" data frame) and make data for Model 2 - modeling whether Endemic vs tolerator pairs explains adj.dif
#####
  colnames(model1)
  adj.dif<-model1 %>%
        dplyr::select(1:17) %>%
        gather("pair","adj.dif",1:17) %>%
        group_by(pair) %>%
        summarise(D = mean(adj.dif), V = var(adj.dif)) %>%
        as.data.frame()
          
  adj.dif$pair <- c(1,10,11,12,13,14,15,16,17,2,3,4,5,6,7,8,9) # change pair column to actual pair numbers
  adj.dif<-adj.dif %>% arrange(pair)
  
  adj.dif<-flw %>%
    group_by(pair_name, pair, pair_type, E_T, year) %>%
    summarise(n()) %>%
    dplyr::select(-`n()`) %>%
    arrange(pair) %>%
    full_join(adj.dif) %>% 
    as.data.frame()
  
  adj.dif <- adj.dif %>%
    mutate(year1 = case_when(year == 2017 ~ 0,
                             year == 2018 ~ 1))
  
  colnames(Ginv1) # double check order of column names in coancestry matrix match pair # order

  
  ### JAGS data
  data.model2 = list(
      D = as.double(adj.dif$D),
      #V = as.double(adj.dif$V),
      E_T = as.double(adj.dif$E_T),
      Ginv = Ginv1,
      mu.0 = as.double(rep(0,17)),
      year = as.double(adj.dif$year1)
    )
#####

##### JAGS model 2 
#####
sink("flw_onset_adj_dif_m2.R")
cat("
    model{
    
    # data structure - every k is a row  - one column for D (=mean adj.dif[k]) and V (= variance adj.dif[k])
    
    for (k in 1:17){
    D[k] ~ dgamma(alpha1[k], beta1[k]) # incorporate sampling error
    
    alpha1[k] = z[k]^2 / sigma  
    beta1[k] = z[k] / sigma
    
    z[k] ~ dgamma(alpha2[k], beta2[k]) # incorporate process error
    
    alpha2[k] = mu[k]^2 / sigma2
    beta2[k] = mu[k] / sigma2
    
    log(mu[k]) <- b0[k] + b1 * E_T[k] + b2 + b3 * year[k]
    
    }

    # Priors on sigma (sampling error) 
    sigma ~ dunif(0,100)

    # Priors on sigma2 (process error)
    sigma2 ~ dunif(0, 1000)
    
    # Priors on regression coefficients
    b0 ~ dmnorm(mu.0, Ginv)
    b1 ~ dnorm(0,tau.b1)
    tau.b1 <- 1 / sigma2b
    sigma2b ~ dunif(0,100)
    b2 ~ dnorm(0,tau.b2)
    tau.b2 <- 1/sigma2c
    sigma2c ~ dunif(0,100)
    b3 ~ dnorm(0,tau.b3)
    tau.b3 <- 1/sigma3
    sigma3 ~ dunif(0,100)
    }
    
    ",fill = TRUE)
sink()

n.adapt = 10000
n.update = 300000
n.iter = 500000

inits = list(
  list(b0 = c(1, 0.8, 0.5, 0.74, 0.15, 3.5, 3.8, 0.9, 0.01, -0.01, 0.32, 3.3, 4, 0.52, 0.23, 0.4, 0.5),
       b1 = 2,
       b2 = 3,
       sigma = 10,
       sigma2b = 5,
       sigma2 = 2),
  list(b0 = c(1, 0.8, 0.5, 0.74, 0.15, 3.5, 3.8, 0.9, 0.01, -0.01, 0.32, 3.3, 4, 0.52, 0.23, 0.4, 0.5),
       b1 = 1.2,
       b2 = .5,
       sigma = 2,
       sigma2b = 1,
       sigma2 = 10),
  list(b0 = c(1, 0.8, 0.5, 0.74, 0.15, 3.5, 3.8, 0.9, 0.01, -0.01, 0.32, 3.3, 4, 0.52, 0.23, 0.4, 0.5),
       b1 = -0.8,
       b2 = 5,
       sigma = 100,
       sigma2b = 10,
       sigma2 = 1)
)
#####

    ##### Run JAGS Model 
    #####
    m2 = jags.model("flw_onset_adj_dif_m2.R", data = data.model2, n.chains = 3, n.adapt = n.adapt)
    update(m2, n.iter = n.update)
    m2c = coda.samples(m2, variable.names = c("z","mu","sigma2","sigma2b","b0","b1","b2","b3"), n.iter = n.iter, n.thin = 100)
    #####


    ##### diagnostics model 2
    #####
    gelman.diag(m2c, multivariate = F)
    #####

    ##### Summary of all 3 chains
    #####
    model2 <- as.data.frame(rbind(m2c[[1]], m2c[[2]], m2c[[3]]))
    dim(model2) #1.5 million rows. Thin out every 100
    model2thin <- model2[seq(1,nrow(model2),100),]
    dim(model2thin)
    #####
    
    ##### b1 posterior and summary statistics - pair type effect
    #####
    hist(model2thin[,"b1"], freq = FALSE, main = "B1 coefficient", xlab = "b1", ylim = c(0,2))
    lines(density(model2thin[,"b1"]), col = "blue", lwd = 2)
    abline(v = quantile(model2thin[,"b1"],c(0.025, .975)), col = "red", lty = 3, lwd = 2) #95% credible interval
    
    quantile(model2thin[, "b1"],c(0.025, .5, .975))
    1-ecdf(model2thin$b1)(0) # 92.3% of the distribution is greater than 0
    
    mean(model2thin[,"b1"]) #untransformed mean of b1 postier
    exp(mean(model2thin[,"b1"])) #transformed mean
    exp(mean(model2thin[,"b1"])) - exp(0) # actual effect, in # of days, that endemic pairs have on increasing onset shift compared to tolerator pairs
    #####

    ##### b2 posterior (intercept) - average shift of tolerators
    #####
    hist(model2thin[,"b2"], freq = FALSE, main = "B2 coefficient", xlab = "b2", ylim = c(0,2))
    lines(density(model2thin[,"b2"]), col = "blue", lwd = 2)
    mean(model2thin[,"b2"])
    exp(mean(model2thin[,"b2"])) # average tolerator shift
    quantile(model2thin[,"b2"],c(0.025, .5, .975))
    #####
    
    ##### b3 posterior (coefficient on year-grown-in-greenhouse effect) - no year effect
    #####
    hist(model2[,"b3"], freq = FALSE, main = "B3 coefficient", xlab = "b2", ylim = c(0,2))
    lines(density(model2[,"b3"]), col = "blue", lwd = 2)
    mean(model2[,"b3"])
    exp(mean(model2[,"b3"]))
    quantile(model2[,"b3"],c(0.025, .5, .975))
    #####

#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## Are flowering onset shifts correlated with greater shifts in the edaphic and/or climatic environment?
#####---------------------#####-----------------------######  
  
##### Relationship between flowering onset shits and shifts in soil habitat

  ## combine data frames
  #####
  soil.distance # multivariate soil distance data frame
  adj.dif # data frame w/ magnitude of flowering time shifts for each pair (adj.dif$D)
  flw.v.soil<-left_join(adj.dif, soil.distance, by = c("pair_name", "pair_type"))
  #####
  
  ##### make sure row names of flw.v.soil match the tip names in the serp.tree
  #####
  flw.v.soil$name <- c("CABE","CABR","CACO","CAGT","TWILD","CLDV","NAJP","NAPB","LADI","COGR","COHT","COSP","MNUD","MGUT","PLER","NARS","NAHX")
  rownames(flw.v.soil) <- flw.v.soil$name
  name.check(serp.tree, flw.v.soil)
  #####
  
  ##### PGLS
  #####
  summary(gls(D ~ distance + year, correlation = corBrownian(phy = serp.tree), data = flw.v.soil, method = "ML")) 
  
  #R2
  m.full <- gls(D ~ distance + year, correlation = corBrownian(phy = serp.tree), data = flw.v.soil, method = "ML")
  m.reduced <- gls(D ~ 1, correlation = corBrownian(phy = serp.tree), data = flw.v.soil, method = "ML")
  
  R2.lik(m.full, m.reduced)
  #####
  
  ##### Fig S1
  #####
  flw.v.soil %>%
    ggplot(.,aes(distance, D, label = pair_name))+
    geom_point(size = 4,aes(shape = pair_type)) +
    #geom_text_repel() +
    #xlab("Multivariate soil distance (Euclidean distance across all PCs)") +
    xlab("Multivariate soil distance") +
    ylab("Absolute shift in flowering onset (days)") +
    scale_shape_manual(values = c(9,16), name = "Pair type", labels = c("Endemic","Tolerator"))+
    theme_classic()+
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14)) 
  #####
  

##### Relationship between flowering onset shits and shifts in climatic environment
  
  ##### get world clim data
  #####
  library(raster)
  r <- getData("worldclim", var = "bio",res = 2.5)

  points<-points  %>% separate(latlon, c("lat","lon"), sep = ",")
  points$lat <- as.numeric(points$lat)
  points$lon <- as.numeric(points$lon) 
  
  coords <- data.frame(x = points$lon, y = points$lat)
  coords<-SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))  

  clim.values <- extract(r,coords)

  clim.df <- cbind.data.frame(coordinates(coords),clim.values)
  clim.df     
  colnames(clim.df)
  colnames(points)
  clim.df <- bind_cols(points[,1:2], clim.df)
  #####
  
  ##### spruce up clim.df
  #####
    head(clim.df)

    clim.df$sp
    clim.df$pair_name <- c("PLER","PLER","COSP","COSP","COGR_COSP","NAPB","NAPB","CACO","CACO","TWILD","TWILD","MGUT","MGUT","MNUD_MGUT","NAJP_NAHN","NAJP_NAHN","CAGT_CAGA","CAGT_CAGA","CLDV_CLHT","CLDV_CLHT","NAHX","NARS_NAHX","NAHX","LADI_LAGL","LADI_LAGL","CABE_CAST","CABE_CAST","CABR","CABR","COHT","COHT")  

    clim.df$pair<-clim.df$pair_name
    clim.df$pair<-gsub("\\\\<CABE_CAST\\\\>","1",clim.df$pair)
    clim.df$pair<-gsub("\\\\<CABR\\\\>","2",clim.df$pair)
    clim.df$pair<-gsub("\\\\<CACO\\\\>","3",clim.df$pair)
    clim.df$pair<-gsub("\\\\<CAGT_CAGA\\\\>","4",clim.df$pair)
    clim.df$pair<-gsub("\\\\<CLDV_CLHT\\\\>","6",clim.df$pair)
    clim.df$pair<-gsub("\\\\<COGR_COSP\\\\>","10",clim.df$pair)
    clim.df$pair<-gsub("\\\\<COHT\\\\>","11",clim.df$pair)
    clim.df$pair<-gsub("\\\\<COSP\\\\>","12",clim.df$pair)
    clim.df$pair<-gsub("\\\\<LADI_LAGL\\\\>","9",clim.df$pair)
    clim.df$pair<-gsub("\\\\<MGUT\\\\>","14",clim.df$pair)
    clim.df$pair<-gsub("\\\\<MNUD_MGUT\\\\>","13",clim.df$pair)
    clim.df$pair<-gsub("\\\\<NAJP_NAHN\\\\>","7",clim.df$pair)
    clim.df$pair<-gsub("\\\\<NAPB\\\\>","8",clim.df$pair)
    clim.df$pair<-gsub("\\\\<PLER\\\\>","15",clim.df$pair)
    clim.df$pair<-gsub("\\\\<NARS_NAHX\\\\>","16",clim.df$pair)
    clim.df$pair<-gsub("\\\\<NAHX\\\\>","17",clim.df$pair)
    clim.df$pair<-gsub("\\\\<TWILD\\\\>","5",clim.df$pair)
    clim.df$pair <- as.numeric(clim.df$pair)
    clim.df<-clim.df %>% arrange(pair)
    
    clim.df$pair_type<-clim.df$pair_name
    clim.df$pair_type<-gsub("\\\\<CABE_CAST\\\\>","E",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<CABR\\\\>","T",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<CACO\\\\>","T",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<CAGT_CAGA\\\\>","E",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<CLDV_CLHT\\\\>","E",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<COGR_COSP\\\\>","E",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<COHT\\\\>","T",clim.df$pair_type) ;	
    clim.df$pair_type<-gsub("\\\\<COSP\\\\>","T",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<LADI_LAGL\\\\>","E",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<MGUT\\\\>","T",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<MNUD_MGUT\\\\>","E",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<NAJP_NAHN\\\\>","E",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<NAPB\\\\>","T",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<PLER\\\\>","T",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<TWILD\\\\>","T",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<NARS_NAHX\\\\>","E",clim.df$pair_type);
    clim.df$pair_type<-gsub("\\\\<NAHX\\\\>","T",clim.df$pair_type)
  #####
  
  ##### generate PCA - for just climate
  #####
    colnames(clim.pca)
    clim.pca <- prcomp(clim.df[,c(5:23)], center = T, scale. = T)
    summary(clim.pca) 
    clim.pca$x
    
    autoplot(clim.pca,data=clim.df,colour="soil",label = T, label.label = "sp" ) +
      theme_classic()+
      theme(axis.text = element_text(size=11, color="black"),plot.title = element_text(hjust = 0.5),legend.key.size = unit(1,"line"),legend.spacing = unit(1,"line"),legend.box.spacing = unit(1,"line"), legend.text = element_text(size = 10), legend.title = element_text(size = 11), axis.title = element_text(size = 12)) +
      scale_linetype_discrete(guide = FALSE)+
      geom_line(aes(group = pair_name, lty = pair_type), lwd = .5)
    
  #####
    
  ##### Caculate Euclidean distances - for just climate 
  #####
    # Make dataframe with ind PC values - 19 total bioclim variables
    ind_values<-as.data.frame(clim.pca$x)

    colnames(clim.df)
    id<-clim.df[,c(1,2,24,25,26)]
    values<-bind_cols(id,ind_values)
    
    # Duplicate rows for COSP and MGUT 
    cosp<-filter(values,sp=="COSP"&soil=="NS")
    cosp$pair<-10
    cosp$pair_name <- "COGR_COSP"
    cosp$pair_type <- "E"

    mgut<-filter(values,sp=="MGUT"&soil=="NS")  
    mgut$pair<-13
    mgut$pair_name <- "MNUD_MGUT"
    mgut$pair_type <- "E"    
    
    nahx<-filter(values,sp=="NAHX"&soil=="NS")  
    nahx$pair<-16
    nahx$pair_name <- "NARS_NAHX"
    nahx$pair_type <- "E"  
    
    values<-bind_rows(values,cosp,mgut,nahx)
   
    # Create empty distance dataframe to store values  
    distances_clim<-data.frame (pair=c(1:17), distance=NA)
    

    # Loop to calculate distances for each pair

    for (i in 1:17){
      a<-values[values$pair==i,]
      distances_clim[i,2]<-Euc_dist_19_PCs.clim(a)
      i=i+1
     }
    
    # spruce up distance_clim dataframe w/ id information
    distances_clim <- clim.df %>%
      group_by(pair, pair_name, pair_type) %>% 
      summarise(n = n()) %>%
      dplyr::select(-n) %>%
      left_join(., distances_clim, by = "pair")
    distances_clim$pair_name <- as.factor(distances_clim$pair_name)
    
    # add in average flowering onset shift
    distances_clim <- adj.dif %>%
      dplyr::select(pair_name, D) %>%
      full_join(.,distances_clim, by = "pair_name")
  
    # add in year-grown-in-GH covariate
    distances_clim <- flw %>% filter(seed == "S") %>% group_by(pair_name,sp, year) %>% summarise(n = n()) %>% dplyr::select(-n) %>% left_join(distances_clim, ., by = "pair_name")
  
    rownames(distances_clim) <- distances_clim$sp # set up row names to match serp tree for pgls
  ######
    
  ##### PGLS for climatic dist vs adj shifts in onset
  #####

  name.check(serp.tree,distances_clim)

  summary(gls(D ~ distance + year, data = distances_clim, correlation = corBrownian(phy = serp.tree), method = "ML" ))
  #####
  
  ##### Fig S2
  #####
  distances_clim %>%
    ggplot(.,aes(distance, D, label = pair_name))+
    geom_point(size = 3, aes(shape = pair_type)) +
    scale_shape_manual(values = c(9,16), name = "Pair type", labels = c("Endemic","Tolerator"))+
    #geom_label_repel()+
    #geom_text(vjust = 1.2, hjust = -.1) +
    theme_classic() +
    xlab("Multivariate climatic distance") + ylab("Absolute difference in flowering \\nonset (days)") +
    #theme(axis.text = element_text(size = rel(1)),
    #      axis.title = element_text(size = rel(1.2))) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16)) 
  ##### 
  
#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## Are flowering onset shifts more genetically-based in endemic pairs than tolerator pairs?
#####---------------------#####-----------------------######    
  
    
##### Linear models of seed (genetic) and soil (plastic) effects on flowering onset shifts. Table S2 
#####
    flw$seed<-relevel(as.factor(flw$seed),ref='NS')
    flw$soil<-relevel(as.factor(flw$soil),ref='NS')
    
    # make dataframe to input coefficients
    
    flw.coef <- data.frame(pair_name = levels(flw$pair_name),
                           coef.int = NA,
                           coef.seed = NA,
                           coef.soil = NA,
                           Var.seed = NA,
                           Var.soil = NA,
                           Var.resid = NA,
                           F.seed = NA,
                           F.soil = NA,
                           p.seed = NA,
                           p.soil = NA)
      
    # make functionto work on a pair
    coef.anal<- function(pair){
      pair1 <- as.character(pair)
      data1 <- filter(flw, pair_name == pair1 & soil!="NS_C" & trt!="NS_S")%>% dplyr::select(pot, sp, pair_name, seed, soil, trt, days_to_flw)
      
      data1$days.std<-scale(data1$days_to_flw,center=T,scale=T)[1:nrow(data1),]
      
      mod<-aov(data = data1, days.std~seed+soil)
      
      SumSq<-anova(mod)[,2]
      VarComp<-SumSq/(length(data1$days.std)-1)
      
      c(coef(mod)[[1]], coef(mod)[[2]],coef(mod)[[3]], VarComp, anova(mod)[,4][1:2],anova(mod)[,5][1:2] )
      
    }
    
    # make loop to run through all of the pairs
    pairs <- levels(flw$pair_name)
    
    for (i in pairs){
    flw.coef[flw.coef$pair_name == i,2:11] <- coef.anal(i)
    }

flw.coef<-flw %>%
  group_by(pair_name, pair_type, year) %>%
  summarise(n = n()) %>%
  dplyr::select(-n) %>%
  left_join(flw.coef, ., by = "pair_name")
        
flw.coef$pair_type <- factor(flw.coef$pair_type, levels = c("T","E"))    
#####

##### Fig 3
#####
lab.key <- c(T = "Tolerator pairs",E = "Endemic pairs")
lab.key
flw.coef %>%
  dplyr::select(pair_name, Var.seed, Var.soil, Var.resid, pair_type) %>% 
  gather("Var.comp","value", 2:4) %>%
  ggplot(.,aes(pair_name, value, fill = Var.comp))+
  geom_bar(stat = "identity") +
  facet_grid (.~pair_type, scales = "free_x", labeller = labeller(pair_type = lab.key)) + 
  ylab( "Percent variance explained" ) +
  theme_classic()+
  scale_fill_manual(values = wes_palette(3, name = "FantasticFox1"), name = "Variance \\ncomponent", labels = c("Residual","Genetic","Plastic")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14),
        #axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) 
#####  

##### PGLS's on variance components
#####

flw.coef
rownames(flw.coef) <- c("MGUT","CABR","CACO","COHT","COSP","NAHX","NAPB","PLER","TWILD","CABE","CAGT","CLDV","COGR","LADI","MNUD","NAJP","NARS")
name.check(serp.tree, flw.coef)


# pgls of variance components on their own
summary(gls(Var.seed ~ pair_type + year,correlation = corBrownian(phy = serp.tree), data = flw.coef, method = "ML"))
summary(gls(Var.soil ~ pair_type + year,correlation = corBrownian(phy = serp.tree), data = flw.coef, method = "ML"))
summary(gls(Var.resid ~ pair_type,correlation = corBrownian(phy = serp.tree), data = flw.coef, method = "ML"))

# pgls of proportion of nonresidual variance that is due to soil or seed
summary(gls(Var.soil/(Var.seed+Var.soil) ~ pair_type + year,correlation = corBrownian(phy = serp.tree), data = flw.coef, method = "ML"))

flw.coef %>%
  mutate(perc_soil = Var.soil/(Var.seed+Var.soil)) %>%
  group_by(pair_type) %>%
  summarise(mean = mean(perc_soil), sem = sem(perc_soil))

summary(gls(Var.seed/(Var.seed+Var.soil) ~ pair_type + year,correlation = corBrownian(phy = serp.tree), data = flw.coef, method = "ML"))

flw.coef %>%
  mutate(perc_soil = Var.seed/(Var.seed+Var.soil)) %>%
  group_by(pair_type) %>%
  summarise(mean = mean(perc_soil), sem = sem(perc_soil))
#####

##### Fig S3
#####
ab<-flw %>%
  mutate(days_to_flw = flw_day - germ_day) %>%
  filter(!c(pair_name=="COSP"&year==2018)) %>%
  filter(trt=="NS_NS"|trt=="S_S"|trt=="S_NS")

ab$trt <- droplevels(ab$trt)  

ggplot(ab,aes(trt, days_to_flw, fill = trt)) +
  geom_boxplot() +
  facet_wrap(~pair_name, scales = "free") + 
  theme(#axis.text.x = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 16),
        legend.position = c(.75,.12),
        legend.direction = "vertical",
        legend.key.size = unit(1,"cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 12))+
  ylab("Flowering onset (days)") +
    #scale_x_discrete(labels=c("NS_NS" = "NS seed\\nNS soil", "S_NS" = "S seed\\nNS soil", "S_S" = "S seed\\nS soil")) +
 scale_fill_manual(values = wes_palette(3, name = "IsleofDogs1"), name = "Seed/soil treatment", labels = c("NS seed, NS soil","S seed, NS soil","S seed, S soil"))
#####



#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## Do endemic sister taxa pairs have stronger and more permanent phenological isolation than tolerator sister taxa pairs?
#####---------------------#####-----------------------######    

    ##### Make census df - column for pot, census day and total open - used to calculate RI
    #####
      cen.nos <- c("cen1","cen2","cen3","cen4","cen5","cen6","cen7","cen8","cen9","cen10","cen11","cen12","cen13","cen14","cen15","cen16","cen17","cen18","cen19","cen20","cen21","cen22","cen23")
      
      census <- data.frame(pot = NA,
                           pair_name = NA,
                           trt = NA,
                           census_day = NA,
                           open_flw = NA
                             )
      
      for(i in 1:length(cen.nos)){
        cen <- cen.nos[i]
        a<-dplyr::select(flw,"pot","pair_name","trt",starts_with(paste(cen, "_",sep=""))) %>%
          dplyr::select(-contains("new")) %>%
          gather("day","tot",5) %>% as.data.frame() %>% dplyr::select(-5)
        day<-a[1,4] # way to make sure that all the day #s are the same for the census
        a[,4]<- day
        colnames(a)<- c("pot","pair_name","trt","census_day","open_flw")
        census<-bind_rows(census,a)
      }


      census<-census %>%
        filter( !c(is.na(census_day))) # remove rows where a pot was not censused on a given census day
      
      
      census[is.na(census$open_flw),"open_flw"] <- 0
      
      census$pair_name <- factor (census$pair_name, levels = c("MGUT","CABR","CACO","COHT","COSP","NAHX","NAPB","PLER","TWILD","CABE_CAST","CAGT_CAGA","CLDV_CLHT","COGR_COSP","LADI_LAGL","MNUD_MGUT","NAJP_NAHN","NARS_NAHX"))
      
    #####
      
    ##### Figure S4-6 - Flowering time distributions
    #####

    # in home soils
    census %>%
        group_by(pair_name, trt, census_day) %>%
        summarise(open_flw_day = sum(open_flw)) %>%
        group_by(pair_name, trt) %>%
        mutate(total = sum(open_flw_day)) %>% 
        filter(trt == "S_S" | trt == "NS_NS") %>%
        ggplot(.,aes(census_day, open_flw_day, color = trt)) +
        geom_line(aes(group = trt), lwd = 0.75) +
        ylab("# open flowers") +
        xlab("census day") +
        facet_wrap(~pair_name, scale = "free_y", ncol = 4) +
        scale_color_manual(values = c("#009988","#EE7733"), name = "Treatment",labels = c("NS taxon, NS soil", "S taxon, S soil")) +
        theme_bw(base_size = 16) 
      
      # in NS soils
      census %>%
        group_by(pair_name, trt, census_day) %>%
        summarise(open_flw_day = sum(open_flw)) %>%
        group_by(pair_name, trt) %>%
        mutate(total = sum(open_flw_day)) %>% 
        filter(trt == "S_NS" | trt == "NS_NS") %>%
        ggplot(.,aes(census_day, open_flw_day, color = trt)) +
        geom_line(aes(group = trt), lwd = .75) +
        ylab("# open flowers") +
        xlab("census day") +
        facet_wrap(~pair_name, scale = "free_y", ncol = 4) +
        scale_color_manual(values = c("#009988","#EE3377"), name = "Treatment",labels = c("NS taxon, NS soil", "S taxon, NS soil")) +
        theme_bw(base_size = 16) 
      
      # in S soils
      census %>%
        group_by(pair_name, trt, census_day) %>%
        summarise(open_flw_day = sum(open_flw)) %>%
        group_by(pair_name, trt) %>%
        mutate(total = sum(open_flw_day)) %>% 
        filter(trt == "S_S" | trt == "NS_S") %>%
        ggplot(.,aes(census_day, open_flw_day, color = trt)) +
        geom_line(aes(group = trt), lwd = 0.75) +
        ylab("# open flowers") +
        xlab("census day") +
        facet_wrap(~pair_name, scale = "free_y", ncol = 4) +
        scale_color_manual(values = c("#33BBEE","#EE7733"), name = "Treatment",labels = c("NS taxon, S soil", "S taxon, S soil")) +
        theme_bw(base_size = 16) 
      
    #####
    

    ##### Calculate RI when pairs are in their home ("adjacent") soil 
    #####
      
        ## Using Sobel/Chen RI-4S1, but # open on day i are the % of total, and the total = 1. Symmetric RI for S and NS
        pGF4S1 <- function (ai,at,bi){ 
          (ai/at) * (bi/(ai+bi))
          }
    
        ## Summary table (="adj) prep - each row is a pair/census day with a column for total # flws in each trt
          adj <- census %>% 
            filter(trt =="S_S" | trt == "NS_NS") %>% 
            group_by(pair_name, trt, census_day) %>% 
            dplyr::summarise(tot_open = sum(open_flw)) %>%
            spread(trt, tot_open)
          head(adj)
    
          # # replace any NAs in the flower counts with 0s
          adj[is.na(adj$S_S),"S_S"]<-0
          adj[is.na(adj$NS_NS),"NS_NS"]<-0
      adj %>% filter(pair_name == "COSP")
 
      # Calculate RI --> RI_adj
      RI_adj<-data.frame(pair_name = unique(adj$pair_name),
                 RI = rep(NA,17))
                
      
      pairs <- unique(adj$pair_name)
      
      for (i in 1:length(pairs)){
        a <- adj %>% filter(pair_name == pairs[i])
        a <- a %>% mutate(NS_tot = sum(NS_NS), S_tot = sum(S_S), 
                          S_rel = S_S/S_tot, NS_rel = NS_NS/NS_tot) %>%
          mutate(RI = pGF4S1(S_rel, sum(S_rel), NS_rel))
        
        a[is.na(a$RI),"RI"] <-0
        
        RI_adj[RI_adj$pair_name == pairs[i],"RI"] <- 1 - 2*(sum(a$RI, na.rm = T))
      }
      
      RI_adj
 
      RI_adj <- flw %>% 
        group_by(pair_name, pair_type) %>%
        dplyr::select(pair_name, pair_type, tot_flw) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::select(-n) %>%
        left_join(.,RI_adj, by = ("pair_name")) %>% as.data.frame()
    #####
    
      
    ##### Bootstap RI when pairs are in their home ("adjacent") soil
    #####
      
      adj_boot_RI_rel <- data.frame(pair_name = rep(unique(census$pair_name),each = 10000),
                                    RI = rep(NA,170000)
      )
    
      pairs <- unique(census$pair_name)
    
    
      for (i in 1:17){
        for (j in 1:10000){
          pair <- pairs[i]
          xyz<- census %>% filter(pair_name == pair)
    
          ss<-xyz %>% filter(trt=="S_S")
          a<-ss %>%
            group_by(census_day) %>%
            mutate(open.bs=sample(open_flw,length(unique(pot)),replace=T)) %>%
            summarise(open.s = sum(open.bs)) %>% as.data.frame()
          a$tot.s <- sum(a$open.s, na.rm = T)
          a<-a %>% mutate(S_rel = open.s / tot.s)
    
          nsns<-xyz %>% filter(trt=="NS_NS")
          b<-nsns %>%
            group_by(census_day) %>%
            mutate(open.bs=sample(open_flw,length(unique(pot)),replace=T)) %>%
            summarise(open.ns = sum(open.bs)) %>% as.data.frame()
          b$tot.ns <- sum(b$open.ns, na.rm = T)
          b<-b %>% mutate(NS_rel = open.ns / tot.ns)
    
          c<-full_join(a,b, by = "census_day")
    
          c<-c %>%
            mutate(RI = pGF4S1(S_rel, sum(S_rel), NS_rel))
          c[is.na(c$RI),"RI"] <- 0
    
          adj_boot_RI_rel[(10000 * i - 10000 + j),2] <- c(1-2*sum(c$RI))
        }
      }


      # add in pair types
      adj_boot_RI_rel<-flw %>%
        group_by(pair_name,  pair_type) %>%
        summarise(n= n()) %>% dplyr::select(-n) %>%
        left_join(adj_boot_RI_rel,., by = "pair_name")
   
    #####
      
      
    ##### Calculate RI when pairs are in the common NS soil
    #####
        ### Summary table (=in.ns) prep - each row is a pair/census day with a column for total # flws in each trt
        in.ns <- census %>% 
          filter(trt =="S_NS" | trt == "NS_NS") %>% 
          group_by(pair_name, trt, census_day) %>% 
          summarise(tot_open = sum(open_flw))%>%
          spread(trt, tot_open)

        in.ns[is.na(in.ns$S_NS),"S_NS"]<-0
        in.ns[is.na(in.ns$NS_NS),"NS_NS"]<-0
  
        ## Calculate RI   = RI_in.ns df

        RI_in.ns<-data.frame(pair_name = unique(adj$pair_name),
                           RI = rep(NA,17))
        
        pairs <- unique(adj$pair_name)
        
        for (i in 1:length(pairs)){
          a <- in.ns %>% filter(pair_name == pairs[i])
          a <- a %>% mutate(NS_tot = sum(NS_NS), S_tot = sum(S_NS), 
                            S_rel = S_NS/S_tot, NS_rel = NS_NS/NS_tot) %>%
            mutate(RI = pGF4S1(S_rel, sum(S_rel), NS_rel))
          a[is.na(a$RI),"RI"] <-0
          
          RI_in.ns[RI_in.ns$pair_name == pairs[i],"RI"] <- 1 - 2*(sum(a$RI, na.rm = T))
        }
        
        RI_in.ns
        
        RI_in.ns <- flw %>% 
          group_by(pair_name, pair_type) %>%
          dplyr::select(pair_name, pair_type, tot_flw) %>%
          summarise(n = n()) %>%
          dplyr::select(-n) %>%
          left_join(.,RI_in.ns, by = ("pair_name")) %>% as.data.frame()
    #####
     
    ##### Bootstap RI when pairs are in NS soil
    #####
      
      ns_boot_RI_rel <- data.frame(pair_name = rep(unique(census$pair_name),each = 10000),
                                    RI = rep(NA,170000)
      )
    
      pairs <- unique(census$pair_name)
    
    
      for (i in 1:17){
        for (j in 1:10000){
          pair <- pairs[i]
          xyz<- census %>% filter(pair_name == pair)
    
          ss<-xyz %>% filter(trt=="S_NS")
          a<-ss %>%
            group_by(census_day) %>%
            mutate(open.bs=sample(open_flw,length(unique(pot)),replace=T)) %>%
            summarise(open.s = sum(open.bs)) %>% as.data.frame()
          a$tot.s <- sum(a$open.s, na.rm = T)
          a<-a %>% mutate(S_rel = open.s / tot.s)
    
          nsns<-xyz %>% filter(trt=="NS_NS")
          b<-nsns %>%
            group_by(census_day) %>%
            mutate(open.bs=sample(open_flw,length(unique(pot)),replace=T)) %>%
            summarise(open.ns = sum(open.bs)) %>% as.data.frame()
          b$tot.ns <- sum(b$open.ns, na.rm = T)
          b<-b %>% mutate(NS_rel = open.ns / tot.ns)
    
          c<-full_join(a,b, by = "census_day")
    
          c<-c %>%
            mutate(RI = pGF4S1(NS_rel, sum(NS_rel), S_rel))
          c[is.na(c$RI),"RI"] <- 0
    
          ns_boot_RI_rel[(10000 * i - 10000 + j),2] <- c(1-2*sum(c$RI))
        }
      }


      # add in pair types
      ns_boot_RI_rel<-flw %>%
        group_by(pair_name,  pair_type) %>%
        summarise(n= n()) %>% dplyr::select(-n) %>%
        left_join(ns_boot_RI_rel,., by = "pair_name")
   
    #####
           
    ##### Calculate permanence of RI between two ecological contexts
    #####

      ### Making dataframe - adj.v.ns
      adj.v.ns <- RI_adj %>% 
        full_join(.,RI_in.ns, by = c("pair_name", "pair_type")) %>%
        mutate(loss.RI = RI.x - RI.y)
      
      colnames(adj.v.ns)[3] <- "In_home_soil"
      colnames(adj.v.ns)[4] <- "In_NS_soil"
      
      # adding year
      colnames(flw)
      adj.v.ns <- flw %>%
        group_by(pair_name, pair, pair_type, year) %>%
        summarise(n = n()) %>%
        ungroup() %>%
        dplyr::select(pair_name, year) %>%
        left_join(adj.v.ns, ., by = "pair_name")
      adj.v.ns
    #####
      
    ##### Summary stats of RI at home, RI in common NS soil and permanence of RI
    #####
    adj.v.ns %>%
      gather("context", "value", 3:5) %>%
      group_by(context) %>%
      summarise(min = min(value),
                max = max(value),
                avg = mean(value),
                se = sem(value))
    
    adj.v.ns %>%
      gather("context", "value", 3:4) %>%
      summarise(avg = mean(value)) # average RI across all scenarios and pairs
      
    #####
      
    ##### Figure S7
    #####
    adj.v.ns
    
    ns_boot_RI_rel
    adj_boot_RI_rel
    
    ggplot(adj.v.ns,aes(pair_name, In_home_soil))+
      geom_violin(data = adj_boot_RI_rel, aes( x = pair_name, y = RI))+
      geom_point() 
    
    ggplot(adj.v.ns,aes(pair_name, In_NS_soil))+
      geom_violin(data = ns_boot_RI_rel, aes( x = pair_name, y = RI))+
      geom_point() 
    
    adj.v.ns <- ns_boot_RI_rel %>% group_by(pair_name) %>% summarise( ns_min = quantile(RI, 0.025), ns_max = quantile(RI,0.975)) %>% full_join(adj.v.ns,.)
    
    adj.v.ns <- adj_boot_RI_rel %>% group_by(pair_name) %>% summarise(adj_min = quantile(RI, 0.025), adj_max = quantile(RI,0.975)) %>% full_join(adj.v.ns,.)
    
   
    ggplot(adj.v.ns, aes(In_NS_soil, In_home_soil,fill = pair_type)) +
      geom_abline(intercept = 0, slope = 1, lty = 2)+
      geom_errorbar( aes(xmin = ns_min, xmax = ns_max), width = .025, alpha = 0.7)+
      geom_errorbar( aes(ymin = adj_min, ymax = adj_max), width = 0.025, alpha = 0.7)+
      geom_point(size = 4, shape = 21)+
      scale_fill_manual(values = c("#018571","#a6611a"), name = "Pair type", labels = c("Endemic","Tolerator"))+
      xlim(0,1) +
      ylim(0,1) +
      xlab("\\nPhenological isolation in pair's NS soil") +
      ylab("Phenological isolation in home soils\\n") +
      theme_minimal(base_size = 14)#+
      #geom_label_repel()
    #####
      
    ##### PGLS's
    #####
    adj.v.ns
    adj.v.ns$name <- c("MGUT","CABR","CACO","COHT","COSP","NAHX","NAPB","PLER","TWILD","CABE","CAGT","CLDV","COGR","LADI","MNUD","NAJP","NARS")
    rownames(adj.v.ns)<-adj.v.ns$name
    name.check(serp.tree, adj.v.ns )
    
    # in home soils
    summary(gls(In_home_soil ~ pair_type + year, correlation = corBrownian(phy = serp.tree), data = adj.v.ns, method = "ML"))  
           
    # in nonserp soils
    summary(gls(In_NS_soil ~ pair_type + year, correlation = corBrownian(phy = serp.tree), data = adj.v.ns, method = "ML"))       
    
    # loss of RI
    summary(gls(loss.RI ~ pair_type + year, correlation = corBrownian(phy = serp.tree), data = adj.v.ns, method = "ML")) 
    #####


#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## Has plasticity in flowering onset evolved following adaptation to serpentine and are plastic responses in an adaptive direction?
#####---------------------#####-----------------------######    


    ##### T-test of maternal family reaction norm slopes within each pair
    #####
    flw<-unite(flw,"fam_ID",c("seed","fam"),sep="_",remove=F) 
    
    ## create dataframe (fam.df) of difference flowering onset between individuals of the same maternal family that are in different soil treatments
    fam.df <- flw %>%
      filter(soil!="NS_C") %>%
      group_by(pair_name, pair_type, soil, fam_ID) %>%
      summarise(days_to_flw = mean(days_to_flw)) %>%# there are 8 families of the 1337 independent families that have 2 individuals from the same maternal familiy in the same soil -- everything else just has one. Here, I avg the flowering onset for those 8 families to use in downstream analyses. 
      spread(soil, days_to_flw) %>%
      mutate(slopeS_NS = S - NS) %>%
      separate(fam_ID, sep = "_", c("seed", "fam"))   %>% ungroup()

    fam.df$seed <- as.factor(fam.df$seed)
 
    ## sample size of the number of families for which I can calculate reaction norm slope - any NAs are seed sources that had no maternal families survive in both the serpentine and nonserpentine soils (typically because of death of all nonserpentine seeds in serpentine soils)
   fam.df %>%
     filter(!is.na(slopeS_NS)) %>%
     group_by(pair_name, seed) %>%
     summarise(n = n()) %>% 
     spread(seed, n)

    ## table of t-tests results w/ the tidy function
    fam.ttest <- fam.df %>%
      filter(!(pair_name == "NAHX" | pair_name == "CABE_CAST" | pair_name == "NARS_NAHX" | pair_name == "CLDV_CLHT" | pair_name == "NAJP_NAHN")) %>% # exclude any pair with a sample size < 2
      filter(!is.na(slopeS_NS)) %>%
      group_by(pair_name) %>%
      do(tidy(t.test(slopeS_NS~seed, data = .))) # "statistic" is t-test, parameter is df, estimate is S seed mean, estimate1 is NS seed mean
    #####
    
    ##### Table S4
    #####
  
    ## add in sample sizes
    fam.ttest <- fam.df %>%
       filter(!is.na(slopeS_NS)) %>%
       group_by(pair_name, seed) %>%
       summarise(n = n()) %>% 
       spread(seed, n) %>%
        left_join(.,fam.ttest, by = "pair_name") %>%
        dplyr::select(pair_name, NS, S, estimate, estimate2, statistic, parameter, p.value)

    colnames(fam.ttest)[4:8] <- c("NS_mean", "S_mean", "t","df","p") 
    
    # add in NS seed or S seed mean family reaction norm slopes for families for which t.tests could not be done
    fam.ttest[fam.ttest$pair_name == "NAHX", "S_mean"] <- fam.df %>% filter(pair_name == "NAHX" & seed == "S") %>% summarise(mean_slope = mean(slopeS_NS, na.rm = T)) %>% as.data.frame()
    
    fam.ttest[fam.ttest$pair_name == "CABE_CAST", "S_mean"] <- fam.df %>% filter(pair_name == "CABE_CAST" & seed == "S") %>% summarise(mean_slope = mean(slopeS_NS, na.rm = T)) %>% as.data.frame()
    
    fam.ttest[fam.ttest$pair_name == "CLDV_CLHT", "S_mean"] <- fam.df %>% filter(pair_name == "CLDV_CLHT" & seed == "S") %>% summarise(mean_slope = mean(slopeS_NS, na.rm = T)) %>% as.data.frame()
    
    fam.ttest[fam.ttest$pair_name == "NAJP_NAHN", "S_mean"] <- fam.df %>% filter(pair_name == "NAJP_NAHN" & seed == "S") %>% summarise(mean_slope = mean(slopeS_NS, na.rm = T)) %>% as.data.frame()
    
    fam.ttest[fam.ttest$pair_name == "NAJP_NAHN", "NS_mean"] <- fam.df %>% filter(pair_name == "NAJP_NAHN" & seed == "NS") %>% summarise(mean_slope = mean(slopeS_NS, na.rm = T)) %>% as.data.frame()
    
    fam.ttest[fam.ttest$pair_name == "NARS_NAHX", "S_mean"] <- fam.df %>% filter(pair_name == "NARS_NAHX" & seed == "S") %>% summarise(mean_slope = mean(slopeS_NS, na.rm = T)) %>% as.data.frame()   
        
    ## create columns that have mean and sample size combined in format of " mean (sample size) "
    for (i in 1:17){
     fam.ttest$NS.mean.n[i] <- paste(round(fam.ttest[i,"NS_mean"],3)," (",fam.ttest[i,"NS"],")")
     fam.ttest$S.mean.n[i] <- paste(round(fam.ttest[i,"S_mean"],3)," (",fam.ttest[i,"S"],")")
    }

    ## sequential bonferonni correction
    bon <- rep(0,12)
      for ( i in 1:12){
        bon[i]<- round(0.05/(12 - i + 1),5)
      }
      bon
    fam.ttest <- fam.ttest  %>% arrange(p)
  
    fam.ttest$bon <- c(bon, rep (NA,5))
    fam.ttest

    fam.ttest <- dplyr::select(fam.ttest, c(pair_name, NS.mean.n, S.mean.n, t, df, p, bon))
    fam.ttest
    #####
     
    ##### Figure S8
    #####
    fam.df %>%
        filter(!is.na(slopeS_NS)) %>%
        ggplot(.,aes(pair_name, slopeS_NS, color = seed)) +
        geom_boxplot(width = .5, position = position_dodge(width = .6), lwd = .75)+
        ylab("Slope of maternal reaction norms\\n(onset in S soil - onset in NS soil)")+
        geom_hline(yintercept = 0, lty = 2)+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 16),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 16),
              axis.title.x = element_blank(),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 14))+
        scale_color_manual(values = c( "#00A08A" , "#F2AD00" ), name = "Taxon", labels = c("Nonserpentine", "Serpentine"))+
        scale_y_continuous(limits = c(-60,80))+
        annotate(geom = "text", x = 5, y = 70, label = "sib in S soil flowers later than sib in NS soil", size = 5, fontface = "italic")+
        annotate(geom = "text", x = 5, y = -50, label = "sib in S soil flowers earlier than sib in NS soil", size = 5, fontface ="italic")
    #####

    ##### Phenotypic selection analysis - linear gradients - Table S5
    #####
        ## create dataframe to measure selection in serpentine soils 
        flw.sel <- flw %>% filter(trt == "S_S" | trt == "NS_S")
        
        ## function to scale
        scale_this <- function(x){
          (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
        }
        
        ## within each pair and seed source, scale the day of first flower
        flw.sel <- flw.sel %>%
          #mutate(flw.std = scale_this(flw_day),
          #       flw.std2 = flw.std^2) %>%
          group_by(pair_name,seed) %>% 
          mutate(flw.std = scale_this(flw_day),
                 flw.std2 = flw.std^2, # quadratic term 
                 mean.flw = mean(tot_flw),
                 w.rel = tot_flw/mean.flw) %>%
          dplyr::select(year, pair_name, seed, trt, flw_day, flw.std, tot_flw, flw.std2, mean.flw, w.rel) 

        
        # sample sizes
        flw.sel %>%
          group_by(pair_name, seed) %>%
          summarise(n = n()) %>%
          spread(seed, n) 
    
        # coefficient table - extract info for slopes
        flw.sel.lm <- flw.sel %>%
          filter(!(seed == "NS" & (pair_name == "NAHX" | pair_name == "CABE_CAST" |pair_name == "NARS_NAHX" | pair_name == "CLDV_CLHT" | pair_name == "COSP"| pair_name == "MGUT" | pair_name == "NAJP_NAHN"))) %>% # remove seed sources with less than 5 individuals
          group_by(pair_name, seed) %>%
          do(tidy(lm(w.rel ~ flw.std, data = .))) # stat is the t-value
          
      
        # add in the r2
        flw.sel.lm <- flw.sel %>%
          filter(!(seed == "NS" & (pair_name == "NAHX" | pair_name == "CABE_CAST" |pair_name == "NARS_NAHX" | pair_name == "CLDV_CLHT" | pair_name == "COSP"| pair_name == "MGUT" | pair_name == "NAJP_NAHN"))) %>%
          group_by(pair_name, seed) %>%
          do(glance(lm(w.rel ~ flw.std, data = .))) %>%
          dplyr::select(pair_name, seed, r.squared) %>%
          left_join(flw.sel.lm, .)
    
        # add in the sample sizes
        flw.sel.lm.TableS4 <- flw.sel %>%
          group_by(pair_name, seed) %>%
          summarise(n = n()) %>%
          spread(seed, n) %>%
          gather("seed", "n", 2:3) %>%
          left_join(., filter(flw.sel.lm, term == "flw.std")) 
    
        flw.sel.lm.TableS4 <- flw.sel.lm.TableS4 %>% 
          arrange(pair_name, seed)
    
 
        bon <- rep(0,27)
          for ( i in 1:27){
            bon[i]<- round(0.05/(27 - i + 1),5)
          }
          bon
        flw.sel.lm.TableS4 <- flw.sel.lm.TableS4  %>% arrange(p.value)
      
        flw.sel.lm.TableS4$bon <- c(bon, rep (NA,7))
        flw.sel.lm.TableS4
        #write.csv(flw.sel.lm.TableS4, "/Users/shelley/Dropbox/Desktop/Flowering time/ProcB/Reviews/new data/linear selection.csv")
    #####

    ##### Fig S9
    #####
        
    xyz <- flw.sel.lm %>%
      select(pair_name, seed, term, estimate) %>%
      spread(term, estimate) %>%
      full_join(., select(filter(flw.sel.lm, term == "flw.std"), p.value))    
    colnames(xyz)[4] <- "flw.std.estimate"
    
    full_join(flw.sel, xyz, by = c("pair_name", "seed")) %>%
      mutate(y = `(Intercept)` + flw.std*flw.std.estimate) %>%
      filter(w.rel < 4) %>%
      ggplot(.,aes(flw.std, w.rel))+
      geom_point(aes(color = seed), alpha = .7)+
      geom_line(aes(flw.std, y, color = seed), data = filter(full_join(flw.sel, xyz, by = c("pair_name", "seed")) %>% mutate( y = `(Intercept)` + flw.std*flw.std.estimate) , p.value < 0.05), lwd = 1)+
      #geom_line(aes(flw.std, y, color = seed), data = filter(full_join(flw.sel, xyz, by = c("pair_name", "seed")) %>% mutate( y = `(Intercept)` + flw.std*flw.std.estimate) , p.value > 0.05 | p.value > 0.1), lty=2, lwd = 1)+
      #facet_wrap(~pair_name, scales = "free") +
      facet_wrap(~pair_name) +
      scale_color_manual(values = c( "#00A08A" , "#F2AD00" ), name = "Taxon", labels = c("Nonserpentine", "Serpentine"))+
      ylab("Relative fitness") +
      xlab("Standardized day of first flower") +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            strip.text.x = element_text(size = 12))
    #####
        
    ##### Phenotypic selection analyses - quadratic and linear factors - Table S1
    #####
    
        flw.sel.quad <- flw.sel %>%
          filter(!(seed == "NS" & (pair_name == "NAHX" | pair_name == "CABE_CAST" |pair_name == "NARS_NAHX" | pair_name == "CLDV_CLHT" | pair_name == "COSP"| pair_name == "MGUT" | pair_name == "NAJP_NAHN"))) %>% # remove seed sources with less than 5 individuals
          group_by(pair_name, seed) %>%
          do(tidy(lm(w.rel ~ flw.std + flw.std2, data = .)))
        
        ## look at p-values for the linear and quad term... what's
        flw.sel.quad %>%
          filter(term == "flw.std" | term == "flw.std2") %>%
          select(-std.error) #%>%
          #write.csv("quadselanalysis.csv")

        bon
        flw.sel.quad2 <- flw.sel.quad %>% filter(term == "flw.std2") %>% arrange(p.value)
        flw.sel.quad2$bon <- c(bon)
   
         
        filter(p.value < 0.05) %>%
        group_by(pair_name, seed) %>%
        spread(term, p.value) %>%
        arrange(pair_name, seed)
        
    #####
        
    ##### Fig S9 - with quadratic equation
    #####
  
    abc <- flw.sel.quad %>%
          select(pair_name, seed, term, estimate) %>%
          spread(term, estimate) %>%
          full_join(., select(filter(flw.sel.quad, term == "flw.std2"), p.value))
     colnames(abc)[4:5] <- c("flw.std.est","flw.std2.est")
     
    full_join(flw.sel, abc, by = c("pair_name", "seed")) %>%   
        mutate(y = `(Intercept)` + flw.std.est*flw.std + flw.std2.est*flw.std2) %>% 
      ggplot(.,aes(flw.std, w.rel, color = seed))+
      geom_point()+
      geom_line(aes(flw.std, y, color = seed), data = filter(full_join(flw.sel, abc, by = c("pair_name", "seed")) %>% mutate(y = `(Intercept)` + flw.std.est*flw.std + flw.std2.est*flw.std2) , p.value < 0.05), lwd = 1)+
      #geom_line(aes(flw.std, y, color = seed), data = filter(full_join(flw.sel, abc, by = c("pair_name", "seed")) %>% mutate(y = `(Intercept)` + flw.std.est*flw.std + flw.std2.est*flw.std2) , p.value > 0.05), lwd = .5, lty = 2)+
      facet_wrap(~pair_name, scales = "free") +
      scale_color_manual(values = c( "#00A08A" , "#F2AD00" ), name = "Taxon", labels = c("Nonserpentine", "Serpentine"))+
      ylab("Flower number") +
      xlab("Standardized day of first flower") +
      theme_classic() +
      theme(axis.title = element_text(size = 16),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            strip.text.x = element_text(size = 12))
    #####

        
##### Fig S10
#####
flw %>%
  filter(trt != "S_NS") %>%
  ggplot(.,aes(trt, days_to_flw, fill = trt)) +
  geom_boxplot() +
  facet_wrap(~pair_name, scales = "free") + 
  theme(#axis.text.x = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 16),
        legend.position = c(.75,.12),
        legend.direction = "vertical",
        legend.key.size = unit(1,"cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 12))+
  ylab("Flowering onset (days)") +
  scale_fill_manual(values = wes_palette(3, name = "Darjeeling1"), name = "Seed/soil treatment", labels = c("NS seed, NS soil","NS seed, S soil","S seed, S soil"))     
#####
    
    
    
    
#*****************  FROM HERE BELOW ARE ALL ANALYSES THAT INVOLVE INCLUDING VARIATION IN GERMINATION TIMING
    
    
    
#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## Variation in germination timing
#####---------------------#####-----------------------######    
    
    
##### Days to germination from planting
#####
    
flw <- flw %>%
      mutate(date_planted2 = case_when(date_planted == 161117 ~ -4, # for 2016-2017 experiement, translate the dates (YYMMDD) to day-of-experiment
                                  date_planted == 161118 ~ -3,
                                  date_planted == 161119 ~ -2,
                                  date_planted == 161120 ~ -1,
                                  date_planted == 161121 ~ 0,
                                  date_planted == 161122 ~ 1,
                                  date_planted == 161123 ~ 2,
                                  date_planted == 161206 ~ 15,
                                  date_planted == 161207 ~ 16,
                                  date_planted == 161214 ~ 23,
                                  date_planted == 161215 ~ 24,
                                  date_planted == 161216 ~ 25,
                                  date_planted == 161219 ~ 28,
                                  )) %>%
      mutate(days_to_germ = case_when(year == 2017 ~ as.numeric(germ_day - date_planted2),
                                     year == 2018 ~ as.numeric(germ_day - date_planted)))
      
flw %>%    
    ggplot(.,aes(soil, days_to_germ, fill = seed)) +
    geom_boxplot()+
    facet_wrap(~pair_name) +
    ylab("Days to germination")+
    xlab("Soil")+
    labs(fill = "Taxon") +
    theme_bw(base_size = 16)


# Testing for differences in germination among pairs. 

## for pairs that have individuals survigin to flower in each of 4 treatments, run this model:
seed.soil.lm <- function(df) {
  model<-summary(lm(days_to_germ ~ seed*soil, data = df))
  coefs <- as.data.frame(round(model$coefficients,3))
  coefs$pair_name <- unique(as.character(df$pair_name))
  write.table(coefs, "/Users/shelley/Dropbox/Desktop/Flowering time/ProcB/DATA/germination_coefs.csv", append = T, sep = ",")
  return(coefs)
}

flw %>%
  filter(!c(pair_name == "NAHX" | pair_name =="CABE_CAST" | pair_name == "NARS_NAHX")) %>%
  by(., .$pair_name, seed.soil.lm)

# for pairs that have no NS individuals surviving to flower in S soils, run this model:
trt.lm <- function(df) {
  model<-summary(lm(days_to_germ ~ trt, data = df))
  coefs <- as.data.frame(round(model$coefficients,3))
  coefs$pair_name <- unique(as.character(df$pair_name))
  write.table(coefs, "/Users/shelley/Dropbox/Desktop/Flowering time/ProcB/DATA/germination_coefs.csv", append = T, sep = ",")
  return(coefs)
}     

flw %>%
  filter(c(pair_name == "NAHX" | pair_name =="CABE_CAST" | pair_name == "NARS_NAHX")) %>%
  by(., .$pair_name, trt.lm)


#####    

#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## Rerun flw time divergence analyses w/ flw onset calc'd from planting date (incorporate germination variation)
#####---------------------#####-----------------------######  

# 'days_to_flw2' is variable that incorporates differences in germination
flw <- flw %>%
  mutate(days_to_flw2 = case_when(year == 2017 ~ as.numeric(flw_day - date_planted2),
                                     year == 2018 ~ as.numeric(flw_day - date_planted)))


##### Average magnitude of flowering time shift across all pairs
#####
flw %>%
  filter(trt=="S_S"|trt=="NS_NS") %>%
  group_by(pair_name, trt) %>%
  summarise(avg.flw = mean(days_to_flw2)) %>%
  spread(trt, avg.flw ) %>%
  mutate(shift = abs(NS_NS - S_S)) %>% ungroup() %>%summarise(avg = mean(shift), sem = sem(shift))
#####

##### T.test for every pair (S_S vs NS_NS), Table S2
#####

# base data frame
tab.s1.germ<-flw %>%
  filter(trt=="S_S"|trt=="NS_NS") %>%
  group_by(pair_name, trt) %>%
  summarise(avg.flw = mean(days_to_flw2)) %>%
  spread(trt, avg.flw ) %>% # each column average flowering onset for each taxon in soil soil
  mutate(dif = S_S - NS_NS) %>% as.data.frame() # difference in the averages within each pair

# attach t.test results for all pairs

tab.s1.germ$t <- NA
tab.s1.germ$df <- NA
tab.s1.germ$p.val <- NA

tab.s1.germ[1,5:7] <- c(round(t.test(data = filter(flw, pair_name == "MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[2,5:7] <- c(round(t.test(data = filter(flw, pair_name == "CABR" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "CABR" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "CABR" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[3,5:7] <- c(round(t.test(data = filter(flw, pair_name == "CACO" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "CACO" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "CACO" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[4,5:7] <- c(round(t.test(data = filter(flw, pair_name == "COHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "COHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "COHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[5,5:7] <- c(round(t.test(data = filter(flw, pair_name == "COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[6,5:7] <- c(round(t.test(data = filter(flw, pair_name == "NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[7,5:7] <- c(round(t.test(data = filter(flw, pair_name == "NAPB" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "NAPB" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "NAPB" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[8,5:7] <- c(round(t.test(data = filter(flw, pair_name == "PLER" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "PLER" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "PLER" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[9,5:7] <- c(round(t.test(data = filter(flw, pair_name == "TWILD" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "TWILD" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "TWILD" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[10,5:7] <- c(round(t.test(data = filter(flw, pair_name == "CABE_CAST" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "CABE_CAST" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "CABE_CAST" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[11,5:7] <- c(round(t.test(data = filter(flw, pair_name == "CAGT_CAGA" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "CAGT_CAGA" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "CAGT_CAGA" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[12,5:7] <- c(round(t.test(data = filter(flw, pair_name == "CLDV_CLHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "CLDV_CLHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "CLDV_CLHT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[13,5:7] <- c(round(t.test(data = filter(flw, pair_name == "COGR_COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "COGR_COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "COGR_COSP" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[14,5:7] <- c(round(t.test(data = filter(flw, pair_name == "LADI_LAGL" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "LADI_LAGL" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "LADI_LAGL" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[15,5:7] <- c(round(t.test(data = filter(flw, pair_name == "MNUD_MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "MNUD_MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "MNUD_MGUT" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))

tab.s1.germ[16,5:7] <- c(round(t.test(data = filter(flw, pair_name == "NAJP_NAHN" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "NAJP_NAHN" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "NAJP_NAHN" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))


tab.s1.germ[17,5:7] <- c(round(t.test(data = filter(flw, pair_name == "NARS_NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[1]][[1]],10), round(t.test(data = filter(flw, pair_name == "NARS_NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[2]][[1]], 10), round(t.test(data = filter(flw, pair_name == "NARS_NAHX" & (trt == "S_S"|trt == "NS_NS")), days_to_flw2 ~ trt)[[3]][[1]], 10))


  tab.s1.germ
  tab.s1.germ<- tab.s1.germ %>% arrange(p.val)
  
    # sequential bonferonni correction
  bon <- rep(0,17)
  for ( i in 1:17){
    bon[i]<- round(0.05/(17 - i + 1),5)
  }
  bon
  tab.s1.germ$bon <- bon

  arrange(tab.s1.germ, pair_name)
#####

#####---------------------#####-----------------------######
#####---------------------#####-----------------------######
## Bayesian analysis - Flw onset w/ germination includes vs. E/T
#####---------------------#####-----------------------######  
  
##### Make JAGS data frame
#####

# JAGS dataframe
    data = list(
      pair.no = as.double(flw$pair),
      soil.type = as.double(flw$soil.type),
      seed.type = as.double(flw$seed.type),
      n = nrow(flw),
      y = as.double(flw$days_to_flw2)
    )
#####
  
##### JAGS model 1
#####
sink("flw_onset_m1.R")
cat("
    model{
    
    # lambda (mean of poisson distribution) prior for all pairs
    for (k in 1:17){
    for (m in 1:2){
    for (j in 1:2){
    lambda[j,m,k] ~ dgamma(0.001, 0.001)
    }
    }
    }

    # every pair/trt gets its own lambda
    for (i in 1:n){
    y[i] ~ dpois(lambda[seed.type[i], soil.type[i], pair.no[i]])
    }
    
    # derived value - adj.dif = S_S - NS_NS
    for (k in 1:17){
    adj.dif[k] <- abs(lambda[1,1,k] - lambda[2,2,k])
    adj.dif.dir[k] <- lambda[1,1,k] - lambda[2,2,k]
    }
    
    }
    ",fill = TRUE)
sink()

n.adapt = 5000
n.update = 10000
n.iter = 10000
#####
  
##### Run JAGS model 1
#####
m1.germ = jags.model("flw_onset_m1.R", data = data, n.chains = 3, n.adapt = n.adapt)
update(m1.germ, n.iter = n.update)
m1c.germ = coda.samples(m1.germ, variable.names = c("lambda","adj.dif", "adj.dif.dir"), n.iter = n.iter, n.thin = 1)
#####
  
    ##### diagnostics model 1 - lambda subsets are [seed, soil, pair] - lambdas w/ convergence diagnostics > 1 are for 2 nonserpentine seed in serpentine soil treatments, which are not used the calculation of flw onset shifts. 
    #####
    gelman.diag(m1c.germ, multivariate = F)
    #plot(m1c)
    #####
    
    ##### summary of all three chains model 1
    #####
    model1.germ <- as.data.frame(rbind(m1c.germ[[1]], m1c.germ[[2]], m1c.germ[[3]]))
    dim(model1.germ)
    #####

    ##### Figure 2 - flowering time divergence
    #####
    colnames(model1.germ)
    fig2.germ<-model1.germ %>%
    dplyr::select(18:34) %>%
    gather("pair","adj.dif.dir",1:17) %>%
    group_by(pair) %>%
    summarise(avg = mean(adj.dif.dir),
              low.95 = quantile(adj.dif.dir,0.025),
              high.95 = quantile(adj.dif.dir, 0.975)) 
    fig2.germ$pair <- c(1, 10, 11, 12, 13, 14, 15, 16, 17, 2, 3, 4, 5, 6, 7, 8, 9)
    fig2.germ <- flw %>% group_by(pair_name, pair_type, pair) %>% summarise(n = n()) %>% dplyr::select(-n) %>% full_join(.,fig2.germ, by = "pair")
    
    ## original fig2, calculated w/ germination included
    fig2.germ %>%  
    ggplot(.,aes(pair_name, avg, fill = pair_type))+
    geom_errorbar(aes(ymin = low.95, ymax = high.95), width = .4)+
      geom_point(size = 3.5, shape = 21) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_minimal() + 
    scale_y_continuous(breaks = seq(-40,40,10)) +
    scale_fill_manual(values = c("#018571","#a6611a"), name = "Pair type", labels = c("Endemic","Tolerator"))+
    scale_color_manual(values = c("#018571","#a6611a"), name = "Pair type", labels = c("Endemic","Tolerator"))+
    theme(#axis.text.x = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12)) +
    xlab("Pair name") + ylab("Average flowering onset shift (days) \\n ")
    
    ##comparison of flw onset shift w/ vs w/out germination included 
    fig2$calculated.from <- "germination"
    fig2.germ$calculated.from <- "planting"
    bind_rows(fig2, fig2.germ) %>%
        ggplot(.,aes(pair_name, avg, fill = pair_type))+
        geom_errorbar(aes(ymin = low.95, ymax = high.95), width = .4)+
        geom_point(size = 3.5, aes(shape = calculated.from)) +
        geom_hline(yintercept = 0, lty = 2) +
        theme_minimal() + 
        scale_y_continuous(breaks = seq(-40,40,10)) +
        scale_shape_manual(values = c(21, 25), name = "Calculated from:", labels = c("germination","planting")) +
        scale_fill_manual(values = c("#018571","#a6611a"), name = "Pair type", labels = c("Endemic","Tolerator"))+
        scale_color_manual(values = c("#018571","#a6611a"), name = "Pair type", labels = c("Endemic","Tolerator"))+
        theme_bw(base_family='Times New Roman') +
        theme(#axis.text.x = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12)) +
        guides(fill = guide_legend(override.aes = list(shape = 21) ))+
        xlab("Pair name") + 
        ylab("Flowering time divergence (days)\\n(S taxon - NS taxon) ")
      
    
    #####
    
    
    
    
##### Summarise mean and variance of posteriors of the onset shifts from model 1 ("adj.dif" data frame) and make data for Model 2 - modeling whether Endemic vs tolerator pairs explains adj.dif
#####
  colnames(model1.germ)
  adj.dif.germ<-model1.germ %>%
        dplyr::select(1:17) %>%
        gather("pair","adj.dif",1:17) %>%
        group_by(pair) %>%
        summarise(D = mean(adj.dif), V = var(adj.dif)) %>%
        as.data.frame()
          
  adj.dif.germ$pair <- c(1,10,11,12,13,14,15,16,17,2,3,4,5,6,7,8,9) # change pair column to actual pair numbers
  adj.dif.germ<-adj.dif %>% arrange(pair)
  
  adj.dif.germ<-flw %>%
    group_by(pair_name, pair, pair_type, E_T, year) %>%
    summarise(n()) %>%
    dplyr::select(-`n()`) %>%
    arrange(pair) %>%
    full_join(adj.dif.germ) %>% 
    as.data.frame()
  
  adj.dif.germ <- adj.dif.germ %>%
    mutate(year1 = case_when(year == 2017 ~ 0,
                             year == 2018 ~ 1))
  
  colnames(Ginv1) # double check order of column names in coancestry matrix match pair # order

  
  ### JAGS data
  data.model2.germ = list(
      D = as.double(adj.dif.germ$D),
      #V = as.double(adj.dif.germ$V),
      E_T = as.double(adj.dif.germ$E_T),
      Ginv = Ginv1,
      mu.0 = as.double(rep(0,17)),
      year = as.double(adj.dif.germ$year1)
    )
#####

##### JAGS model 2 
#####
sink("flw_onset_adj_dif_m2.R")
cat("
    model{
    
    # data structure - every k is a row  - one column for D (=mean adj.dif[k]) and V (= variance adj.dif[k])
    
    for (k in 1:17){
    D[k] ~ dgamma(alpha1[k], beta1[k]) # incorporate sampling error
    
    alpha1[k] = z[k]^2 / sigma  
    beta1[k] = z[k] / sigma
    
    z[k] ~ dgamma(alpha2[k], beta2[k]) # incorporate process error
    
    alpha2[k] = mu[k]^2 / sigma2
    beta2[k] = mu[k] / sigma2
    
    log(mu[k]) <- b0[k] + b1 * E_T[k] + b2 + b3 * year[k]
    
    }

    # Priors on sigma (sampling error) 
    sigma ~ dunif(0,100)

    # Priors on sigma2 (process error)
    sigma2 ~ dunif(0, 1000)
    
    # Priors on regression coefficients
    b0 ~ dmnorm(mu.0, Ginv)
    b1 ~ dnorm(0,tau.b1)
    tau.b1 <- 1 / sigma2b
    sigma2b ~ dunif(0,100)
    b2 ~ dnorm(0,tau.b2)
    tau.b2 <- 1/sigma2c
    sigma2c ~ dunif(0,100)
    b3 ~ dnorm(0,tau.b3)
    tau.b3 <- 1/sigma3
    sigma3 ~ dunif(0,100)
    }
    
    ",fill = TRUE)
sink()

n.adapt = 10000
n.update = 300000
n.iter = 500000

inits = list(
  list(b0 = c(1, 0.8, 0.5, 0.74, 0.15, 3.5, 3.8, 0.9, 0.01, -0.01, 0.32, 3.3, 4, 0.52, 0.23, 0.4, 0.5),
       b1 = 2,
       b2 = 3,
       sigma = 10,
       sigma2b = 5,
       sigma2 = 2),
  list(b0 = c(1, 0.8, 0.5, 0.74, 0.15, 3.5, 3.8, 0.9, 0.01, -0.01, 0.32, 3.3, 4, 0.52, 0.23, 0.4, 0.5),
       b1 = 1.2,
       b2 = .5,
       sigma = 2,
       sigma2b = 1,
       sigma2 = 10),
  list(b0 = c(1, 0.8, 0.5, 0.74, 0.15, 3.5, 3.8, 0.9, 0.01, -0.01, 0.32, 3.3, 4, 0.52, 0.23, 0.4, 0.5),
       b1 = -0.8,
       b2 = 5,
       sigma = 100,
       sigma2b = 10,
       sigma2 = 1)
)
#####

    ##### Run JAGS Model 
    #####
    m2.germ = jags.model("flw_onset_adj_dif_m2.R", data = data.model2.germ, n.chains = 3, n.adapt = n.adapt)
    update(m2.germ, n.iter = n.update)
    m2c.germ = coda.samples(m2.germ, variable.names = c("z","mu","sigma2","sigma2b","b0","b1","b2","b3"), n.iter = n.iter, n.thin = 100)
    #####


    ##### diagnostics model 2
    #####
    gelman.diag(m2c.germ, multivariate = F)
    #####

    ##### Summary of all 3 chains
    #####
    model2.germ <- as.data.frame(rbind(m2c.germ[[1]], m2c.germ[[2]], m2c.germ[[3]]))
    dim(model2.germ) #1.5 million rows. Thin out every 100
    model2thin.germ <- model2.germ[seq(1,nrow(model2.germ),100),]
    dim(model2thin.germ)
    #####
    
    ##### b1 posterior and summary statistics - pair type effect
    #####
    hist(model2thin.germ[,"b1"], freq = FALSE, main = "B1 coefficient", xlab = "b1", ylim = c(0,2))
    lines(density(model2thin.germ[,"b1"]), col = "blue", lwd = 2)
    abline(v = quantile(model2thin.germ[,"b1"],c(0.025, .975)), col = "red", lty = 3, lwd = 2) #95% credible interval
    
    quantile(model2thin.germ[, "b1"],c(0.025, .5, .975))
    1-ecdf(model2thin.germ$b1)(0) # 92.3% of the distribution is greater than 0
    
    mean(model2thin.germ[,"b1"]) #untransformed mean of b1 postier
    exp(mean(model2thin.germ[,"b1"])) #transformed mean
    exp(mean(model2thin.germ[,"b1"])) - exp(0) # actual effect, in # of days, that endemic pairs have on increasing onset shift compared to tolerator pairs
    #####

    ##### b2 posterior (intercept) - average shift of tolerators
    #####
    hist(model2thin.germ[,"b2"], freq = FALSE, main = "B2 coefficient", xlab = "b2", ylim = c(0,2))
    lines(density(model2thin.germ[,"b2"]), col = "blue", lwd = 2)
    mean(model2thin.germ[,"b2"])
    exp(mean(model2thin.germ[,"b2"])) # average tolerator shift
    quantile(model2thin.germ[,"b2"],c(0.025, .5, .975))
    #####
    
    ##### b3 posterior (coefficient on year-grown-in-greenhouse effect) - no year effect
    #####
    hist(model2.germ[,"b3"], freq = FALSE, main = "B3 coefficient", xlab = "b2", ylim = c(0,2))
    lines(density(model2.germ[,"b3"]), col = "blue", lwd = 2)
    mean(model2.germ[,"b3"])
    exp(mean(model2.germ[,"b3"]))
    quantile(model2.germ[,"b3"],c(0.025, .5, .975))
    #####


