rm(list=ls())

library(metafor)
library(ggtree)
library(tidyverse)
library(readxl)
library(ape)
library(phytools)
library(MCMCglmm)
library(cowplot)

# Importing data ----------------------------------------------------------


#Set the working directory
setwd("...")

study = read_excel("data.xlsx", sheet="study")
species = read_excel("data.xlsx", sheet="species")
observation = read_excel("data.xlsx", sheet="observation")
data = full_join(species, observation, by='commonName') %>% 
  full_join(study, by='citation')

#BRINGING IN PHYLOGENETIC TREE
tree = read.nexus("output.nex")
tree = averageTree(tree, method='symmetric.difference')

#Replacing old species names with modern ones
tree$tip.label = gsub('Dendroica_caerulescens', 'Setophaga_caerulescens', tree$tip.label)
tree$tip.label = gsub('Dendroica_chrysoparia', 'Setophaga_chrysoparia', tree$tip.label)
tree$tip.label = gsub('Dendroica_kirtlandii', 'Setophaga_kirtlandii', tree$tip.label)
tree$tip.label = gsub('Dendroica_petechia', 'Setophaga_petechia', tree$tip.label)
tree$tip.label = gsub('Dendroica_virens', 'Setophaga_virens', tree$tip.label)


#Inverting relatedness matrix for phylogenetic meta-analysis
invB = inverseA(tree, nodes='TIPS')

#Separating occupancy and density/abundance effects
dens = data %>% filter(responseCategory=='Density/Abundance')
occ = data %>% filter(responseCategory=='Occupancy')


#-------------------------------------------------------------------------
#Functions for collecting relevant data from models
#-------------------------------------------------------------------------

collectSamples = function(model,var,covType){
  if(covType=='categorical'){
    tmp1 = unique(tmpDens[,c("commonName", var)])
    tmp1 = data.frame(table(tmp1[[var]]))
    colnames(tmp1) = c("covariate", "nSpecies")

    tmp2 = unique(tmpDens[,c("citation", var)])
    tmp2 = data.frame(table(tmp2[[var]]))
    colnames(tmp2) = c("covariate", "nStudies")

    tmp3 = data.frame(table(tmpDens[[var]]))
    colnames(tmp3) = c("covariate", "n")

    res = merge(merge(tmp3, tmp1, by="covariate"), tmp2, by="covariate")
    res$response = "Density"

    tmp1 = unique(tmpOcc[,c("commonName", var)])
    tmp1 = data.frame(table(tmp1[[var]]))
    colnames(tmp1) = c("covariate", "nSpecies")

    tmp2 = unique(tmpOcc[,c("citation", var)])
    tmp2 = data.frame(table(tmp2[[var]]))
    colnames(tmp2) = c("covariate", "nStudies")

    tmp3 = data.frame(table(tmpOcc[[var]]))
    colnames(tmp3) = c("covariate", "n")

    tmp = merge(merge(tmp3, tmp1, by="covariate"), tmp2, by="covariate")
    tmp$response = "Occupancy"
    res = rbind(res, tmp)
    res$model = model
  } else{

    res = data.frame("covariate" = var,
                     "nSpecies" = length(unique(tmpDens$commonName)),
                     "nStudies" = length(unique(tmpDens$Authors)),
                     "n" = nrow(tmpDens),
                     "response" = "Density")

    tmp = data.frame("covariate" = var,
                     "nSpecies" = length(unique(tmpOcc$commonName)),
                     "nStudies" = length(unique(tmpOcc$Authors)),
                     "n" = nrow(tmpOcc),
                     "response" = "Occupancy")
    res = rbind(res, tmp)
    res$model = model
  }
  return(res)
}

collectResults = function(d,o,model,var, covType){
  tmp1 = data.frame(rbind(summary(d)$Gcovariances, summary(d)$Rcovariances))
  tmp1$pMCMC = NA
  tmp1$covariate = "random"
  if(covType=="categorical"){
    tmp1$n = paste("n = ", sum(table(tmpDens[[var]])), sep="")
  } else{
    tmp1$n = paste("n = ", nrow(tmpDens), sep="")
  }
  tmp2 = data.frame(summary(d)$solutions)
  tmp2$covariate = "fixed"
  if(covType=="categorical"){
    tmp2$n = paste("n = ", as.numeric(table(tmpDens[[var]])), sep="")
  } else{
    tmp2$n = paste("n = ", nrow(tmpDens), sep="")
  }
  tmp1 = rbind(tmp1, tmp2)
  tmp1$response = "density"
  tmp1$variable = row.names(tmp1)
  
  tmp3 = data.frame(rbind(summary(o)$Gcovariances, summary(o)$Rcovariances))
  tmp3$pMCMC = NA
  tmp3$covariate = "random"
  if(covType=="categorical"){
    tmp3$n = paste("n = ", sum(table(tmpOcc[[var]])), sep="")
  } else{
    tmp3$n = paste("n = ", nrow(tmpOcc), sep="")
  }
  tmp4 = data.frame(summary(o)$solutions)
  tmp4$covariate = "fixed"
  if(covType=="categorical"){
    tmp4$n = paste("n = ", as.numeric(table(tmpOcc[[var]])), sep="")
  } else{
    tmp4$n = paste("n = ", nrow(tmpOcc), sep="")
  }
  tmp3 = rbind(tmp3, tmp4)
  tmp3$response = "occupancy"
  tmp3$variable = row.names(tmp3)
  
  tmp = rbind(tmp1, tmp3)
  tmp$model = model
  
  tmp$variable = sub(var, "", tmp$variable)
  return(tmp)
}

comparisons = function(d, o, model, var){
  tmp1 = tmp2 = matrix(NA, nrow=d$Sol, ncol=choose(ncol(d$Sol),2))
  for(i in 1:(ncol(d$Sol)-1)){
    for(j in (i+1):ncol(d$Sol)){
      tmp1 = d$Sol[,i] - d$Sol[,j]
      tmp2 = o$Sol[,i] - o$Sol[,j]
      
      tmp3 = data.frame(cbind(mean(tmp1), HPDinterval(tmp1)))
      tmp4 = data.frame(cbind(mean(tmp2), HPDinterval(tmp2)))
      tmp3$response = "density"
      tmp4$response = "occupancy"
      tmp3$model = tmp4$model = model
      tmp3$var1 = tmp4$var1 = colnames(d$Sol)[i]
      tmp3$var2 = tmp4$var2 = colnames(d$Sol)[j]
      tmp3$var1 = sub(var, "", tmp3$var1)
      tmp3$var2 = sub(var, "", tmp3$var2)
      tmp4$var1 = sub(var, "", tmp4$var1)
      tmp4$var2 = sub(var, "", tmp4$var2)
      
      tmp3$n1 = paste("n = ", as.numeric(table(tmpDens[[var]]))[i], sep="")
      tmp3$n2 = paste("n = ", as.numeric(table(tmpDens[[var]]))[j], sep="")
      tmp4$n1 = paste("n = ", as.numeric(table(tmpOcc[[var]]))[i], sep="")
      tmp4$n2 = paste("n = ", as.numeric(table(tmpOcc[[var]]))[j], sep="")
      
      if (i==1 & j==2){
        tmp = rbind(tmp3, tmp4)
      } else{
        tmp = rbind(rbind(tmp, tmp3), tmp4)
      }
    }
  }
  tmp$sig = ifelse(tmp$lower <= 0 & tmp$upper >= 0, 0, 1)
  return(tmp)
}


#-------------------------------------------------------------------------
#Species models
#-------------------------------------------------------------------------

tmpDens = data.frame(dens)
tmpOcc = data.frame(occ)

prior = list(R=list(V=1, nu=1),
             G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.v=1000),
                    G2=list(V=1, fix=1)))

densSpecies = MCMCglmm(yi~commonName-1,
                       data=tmpDens,
                       random = ~citation + idh(SE):units,
                       prior=prior,
                       nitt=400000,
                       burnin=200000,
                       thin=100)

occSpecies = MCMCglmm(yi~commonName-1,
                      data=tmpOcc,
                      random = ~citation + idh(SE):units,
                      prior=prior,
                      nitt=400000,
                      burnin=200000,
                      thin=100)

#For checking convergence and effective sample sizes
# summary(densSpecies)
# plot(densSpecies$Sol)
# plot(densSpecies$VCV)
# 
# summary(occSpecies)
# plot(occSpecies$Sol)
# plot(occSpecies$VCV)

results = collectResults(densSpecies, occSpecies, "species", "commonName", "categorical")

resultsCount = collectSamples("species", "commonName", "categorical")


#-----------------------------------------------------------
#Creating phylogenetic effect size figure
#-----------------------------------------------------------

tmp = unique(data[,c("commonName", "scientificName"),])
tmp$commonName = as.character(tmp$commonName)
tree2 = tree
tree2$tip.label = tmp$commonName[match(tree2$tip.label, tmp$scientificName)]

tmp1 = results %>% filter(response=='density')
tmp2 = results %>% filter(response=='occupancy')

colnames(tmp1)[9] = colnames(tmp2)[9] = "commonName"

tmp1 = tmp %>% 
  full_join(tmp1, by='commonName') %>% 
  mutate(scientificName = gsub("_", " ", scientificName))
tmp2 = tmp %>% 
  full_join(tmp2, by='commonName') %>% 
  mutate(scientificName = gsub("_", " ", scientificName))


colnames(tmp1)[1] = colnames(tmp2)[1] = "id"

tmp1 = tmp1 %>% 
  mutate(commonName = id) %>% 
  mutate(Effects = as.integer(gsub("n = ", "", n)))

tmp2 = tmp2 %>% 
  mutate(commonName = id) %>% 
  mutate(Effects = as.integer(gsub("n = ", "", n)))


p = ggtree(tree2)
p2 = facet_plot(p, panel="Species", data=tmp1, geom=geom_text,
                aes(x=0, label=commonName), size=2.6)+
  theme_tree2()
p3 = facet_plot(p2, panel="Density effect", data=tmp1, geom=geom_errorbarh, aes(x=post.mean, xmin=l.95..CI, xmax=u.95..CI))
p4 = facet_plot(p3, panel="Density effect", data=tmp1, geom=geom_point, aes(x=post.mean, size=Effects))
p5 = facet_plot(p4, panel="Density effect", data=tmp1, geom=geom_vline, aes(xintercept=0), linetype='dashed', alpha=1)
p6 = facet_plot(p5, panel="Presence effect", data=tmp2, geom=geom_errorbarh, aes(x=post.mean, xmin=l.95..CI, xmax=u.95..CI))
p7 = facet_plot(p6, panel="Presence effect", data=tmp2, geom=geom_point, (aes(x=post.mean, size=Effects)))
p8 = facet_plot(p7, panel="Presence effect", data=tmp2, geom=geom_vline, aes(xintercept=0), linetype='dashed', alpha=1)
p9 = facet_labeller(p8, c(Tree='Phylogenetic tree', Species='Species'))+
  theme(strip.background=element_rect(fill='white', color='black'))+
  xlim_expand(c(-1,1), panel='Species')+
  theme(legend.position='right')+
  labs(size = "Sample size")+
  theme(legend.title=element_text(size=8))+
  theme(strip.text=element_text(size=8))

p9

#--------------------------------------------------------------
#Intercept models
#--------------------------------------------------------------

tmpDens = data.frame(dens)
tmpOcc = data.frame(occ)

prior = list(R=list(V=1, nu=1),
             G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.v=1000),
                    G2=list(V=1, nu=1, alpha.mu=0, alpha.v=1000),
                    G3=list(V=1, nu=1, alpha.mu=0, alpha.v=1000),
                    G4=list(V=1, fix=1)))

densIntercept = MCMCglmm(yi~1,
                         data=tmpDens,
                         random = ~scientificName + commonName + citation + idh(SE):units,
                         ginverse = list(scientificName=invB$Ainv),
                         prior=prior,
                         nitt=400000,
                         burnin=200000,
                         thin=100,
                         pr=T)

occIntercept = MCMCglmm(yi~1,
                        data=tmpOcc,
                        random = ~scientificName + commonName + citation + idh(SE):units,
                        ginverse = list(scientificName=invB$Ainv),
                        prior=prior,
                        nitt=400000,
                        burnin=200000,
                        thin=100,
                        pr=T)

#For checking effective sample sizes and convergence
# summary(densIntercept)
# plot(densIntercept$Sol)
# plot(densIntercept$VCV)
# 
# summary(occIntercept)
# plot(occIntercept$Sol)
# plot(occIntercept$VCV)

results = collectResults(densIntercept, occIntercept, "intercept", "intercept", "continuous")

#---------------------------------------------------------------------
#What proportion of variance is explained by phylogenetic relatedness?
#---------------------------------------------------------------------

#Heritability - about 30% of the density response can be attributed to phylogeny
tmp = densIntercept$VCV[,1]/(apply(densIntercept$VCV[,c(1,2,3,5)], 1, FUN='sum'))
mean(tmp)
quantile(tmp, c(0.025, 0.975))

#About 26% of the occupancy response can be explained by phylogeny
tmp = occIntercept$VCV[,1]/(apply(occIntercept$VCV[,c(1,2,3,5)], 1, FUN='sum'))
mean(tmp)
quantile(tmp, c(0.025, 0.975))

#--------------------------------------------------------------
#Effect size and funnel plots for the intercept models
#--------------------------------------------------------------

#Dens, marginalizing out units
tmpDens$prediction = predict(densIntercept, marginal=~idh(SE):units)
tmpDens$resid = tmpDens$yi - tmpDens$prediction

#Testing for publication bias in density/abundance studies
tmp1 = rma(yi=resid, vi=vi, data=tmpDens)
regtest(tmp1)

#Occ, marginalizing out units
tmpOcc$prediction = predict(occIntercept, marginal=~idh(SE):units)
tmpOcc$resid = tmpOcc$yi - tmpOcc$prediction

#Testing for publication bias in presence-absence studies
tmp2 = rma(yi=resid, vi=vi, data=tmpOcc)
regtest(tmp2)

tmp = results[which(results$covariate=='fixed'),]
tmp$response = ifelse(tmp$response=='occupancy', 'Pres-abs', 'Dens/abund')


a = ggplot(tmp, aes(x=response, y=post.mean))+
  geom_point()+
  theme_bw()+
  theme(panel.grid=element_blank())+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI))+
  geom_abline(intercept=0, slope=0, linetype='dashed')+
  coord_flip()+
  geom_text(data=tmp, aes(label=n), nudge_x=0.2)+
  ylab("Effect of social information treatments")+
  xlab("")+
  theme(axis.text.y = element_text(angle=90, vjust=1, hjust=0.5))+
  theme(axis.title = element_text(size=9))

tmpDens = tmpDens %>% 
  mutate(se = seq(0, 1.2, length.out=nrow(.))) %>% 
  mutate(lcl = coef(tmp1)-1.96*se) %>% 
  mutate(ucl = coef(tmp1)+1.96*se) %>% 
  mutate(residMean = rep(coef(tmp1), nrow(.)))


b = ggplot(tmpDens, aes(x=resid, y=sqrt(vi)))+
  geom_point()+
  theme_bw()+
  theme(panel.grid=element_blank())+
  scale_y_reverse(breaks=seq(0,1.2,0.2))+
  geom_line(aes(x=lcl, y=se), linetype='solid')+
  geom_line(aes(x=ucl, y=se), linetype='solid')+
  geom_line(aes(x=residMean, y=se), linetype='dashed')+
  xlab("Density/abundance model residuals")+
  ylab("Effect size SE")+
  theme(axis.title=element_text(size=9))

tmpOcc = tmpOcc %>% 
  mutate(se = seq(0, 2.2, length.out=nrow(.))) %>% 
  mutate(lcl = coef(tmp2)-1.96*se) %>% 
  mutate(ucl = coef(tmp2)+1.96*se) %>% 
  mutate(residMean = rep(coef(tmp2), nrow(.)))

c = ggplot(tmpOcc, aes(x=resid, y=sqrt(vi)))+
  geom_point()+
  theme_bw()+
  theme(panel.grid=element_blank())+
  scale_y_reverse(breaks=seq(0,2.4,0.4))+
  xlim(-4.5, 4.5)+
  geom_line(aes(x=lcl, y=se), linetype='solid')+
  geom_line(aes(x=ucl, y=se), linetype='solid')+
  geom_line(aes(x=residMean, y=se), linetype='dashed')+
  xlab("Presence-absence model residuals")+
  ylab("Effect size SE")+
  theme(axis.title=element_text(size=9))

plot_grid(a, b, c, labels=c("A", "B", "C"), ncol=1)


#------------------------------------------------------------------
#Moderator models
#------------------------------------------------------------------

prior = list(R=list(V=1, nu=1),
             G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.v=1000),
                    G2=list(V=1, nu=1, alpha.mu=0, alpha.v=1000),
                    G3=list(V=1, nu=1, alpha.mu=0, alpha.v=1000),
                    G4=list(V=1, fix=1)))

####
#Site quality
####
tmpDens = dens %>% 
  mutate(siteQuality = as.character(siteQuality)) %>% 
  mutate(siteQuality = ifelse(siteQuality=='Gradient', 'Gradient', 'High')) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(siteQuality = as.character(siteQuality)) %>% 
  mutate(siteQuality = ifelse(siteQuality=='Gradient', 'Gradient', 'High')) %>% 
  data.frame()


densQuality = MCMCglmm(yi~siteQuality-1,
                       data=tmpDens,
                       random = ~scientificName + commonName + citation + idh(SE):units,
                       ginverse = list(scientificName=invB$Ainv),
                       prior=prior,
                       nitt=400000,
                       burnin=200000,
                       thin=100)

occQuality = MCMCglmm(yi~siteQuality-1,
                      data=tmpOcc,
                      random = ~scientificName + commonName + citation + idh(SE):units,
                      ginverse = list(scientificName=invB$Ainv),
                      prior=prior,
                      nitt=400000,
                      burnin=200000,
                      thin=100)

#Checking sample sizes and convergence
# summary(densQuality)
# plot(densQuality$Sol)
# plot(densQuality$VCV)
# 
# summary(occQuality)
# plot(occQuality$Sol)
# plot(occQuality$VCV)


results = collectResults(densQuality, occQuality,
                         "Site quality", "siteQuality", "categorical")

resultsComp = comparisons(densQuality, occQuality, "Site quality", "siteQuality")



####
#Previous occupancy
####
tmpDens = dens %>% 
  mutate(previousOccupancy = as.character(previousOccupancy)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(previousOccupancy = as.character(previousOccupancy)) %>% 
  data.frame()


densPrevOcc = MCMCglmm(yi~previousOccupancy-1,
                       data=tmpDens,
                       random = ~scientificName + commonName + citation + idh(SE):units,
                       ginverse = list(scientificName=invB$Ainv),
                       prior=prior,
                       nitt=400000,
                       burnin=200000,
                       thin=100)

occPrevOcc = MCMCglmm(yi~previousOccupancy-1,
                      data=tmpOcc,
                      random = ~scientificName + commonName + citation + idh(SE):units,
                      ginverse = list(scientificName=invB$Ainv),
                      prior=prior,
                      nitt=400000,
                      burnin=200000,
                      thin=100)

#Checking sample sizes and convergence
# summary(densPrevOcc)
# plot(densPrevOcc$Sol)
# plot(densPrevOcc$VCV)
# 
# summary(occPrevOcc)
# plot(occPrevOcc$Sol)
# plot(occPrevOcc$VCV)


results = rbind(results, collectResults(densPrevOcc, occPrevOcc,
                                        "Previous occ.", "previousOccupancy", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densPrevOcc, occPrevOcc, "Previous occ.",
                                             "previousOccupancy"))




####
#Site size
####
tmpDens = dens %>% 
  mutate(siteSize = scale(as.numeric(siteSize), center=T, scale=T)) %>% 
  filter(!is.na(siteSize)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(siteSize = scale(as.numeric(siteSize), center=T, scale=T)) %>% 
  filter(!is.na(siteSize)) %>% 
  data.frame()

densSize = MCMCglmm(yi~siteSize,
                    data=tmpDens,
                    random = ~scientificName + commonName + citation + idh(SE):units,
                    ginverse = list(scientificName=invB$Ainv),
                    prior=prior,
                    nitt=400000,
                    burnin=200000,
                    thin=100)

occSize = MCMCglmm(yi~siteSize,
                   data=tmpOcc,
                   random = ~scientificName + commonName + citation + idh(SE):units,
                   ginverse = list(scientificName=invB$Ainv),
                   prior=prior,
                   nitt=400000,
                   burnin=200000,
                   thin=100)

#Checking sample sizes and convergence
# summary(densSize)
# plot(densSize$Sol)
# plot(densSize$VCV)
# 
# summary(occSize)
# plot(occSize$Sol)
# plot(occSize$VCV)

results = rbind(results, collectResults(densSize, occSize,
                                        "Size", "siteSize", "continuous"))




####
#Site latitude
####
tmpDens = dens %>% 
  mutate(latitude = scale(latitude, center=T, scale=T)) %>% 
  filter(!is.na(latitude)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(latitude = scale(latitude, center=T, scale=T)) %>% 
  filter(!is.na(latitude)) %>% 
  data.frame()


densLat = MCMCglmm(yi~latitude,
                   data=tmpDens,
                   random = ~scientificName + commonName + citation + idh(SE):units,
                   ginverse = list(scientificName=invB$Ainv),
                   prior=prior,
                   nitt=400000,
                   burnin=200000,
                   thin=100)

occLat = MCMCglmm(yi~latitude,
                  data=tmpOcc,
                  random = ~scientificName + commonName + citation + idh(SE):units,
                  ginverse = list(scientificName=invB$Ainv),
                  prior=prior,
                  nitt=400000,
                  burnin=200000,
                  thin=100)

#Checking sample sizes and convergence
# summary(densLat)
# plot(densLat$Sol)
# plot(densLat$VCV)
# 
# summary(occLat)
# plot(occLat$Sol)
# plot(occLat$VCV)


results = rbind(results, collectResults(densLat, occLat,
                                        "Latitude", "latitude", "continuous"))


####
#Conspecific density
####
tmpDens = dens %>% 
  mutate(eBirdFreq = scale(eBirdFreq, center=T, scale=T)) %>% 
  filter(!is.na(eBirdFreq)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(eBirdFreq = scale(eBirdFreq, center=T, scale=T)) %>% 
  filter(!is.na(eBirdFreq)) %>% 
  data.frame()

densDist = MCMCglmm(yi~eBirdFreq,
                    data=tmpDens,
                    random = ~scientificName + commonName + citation + idh(SE):units,
                    ginverse = list(scientificName=invB$Ainv),
                    prior=prior,
                    nitt=400000,
                    burnin=200000,
                    thin=100)

occDist = MCMCglmm(yi~eBirdFreq,
                   data=tmpOcc,
                   random = ~scientificName + commonName + citation + idh(SE):units,
                   ginverse = list(scientificName=invB$Ainv),
                   prior=prior,
                   nitt=400000,
                   burnin=200000,
                   thin=100)

#Checking sample sizes and convergence
# summary(densDist)
# plot(densDist$Sol)
# plot(densDist$VCV)
# 
# summary(occDist)
# plot(occDist$Sol)
# plot(occDist$VCV)


results = rbind(results, collectResults(densDist, occDist,
                                        "Local dist.", "eBirdFreq", "continuous"))


####
#Response measured
####
####
tmpDens = dens %>% 
  mutate(responseMeasured = as.character(responseMeasured)) %>% 
  filter(!is.na(responseMeasured)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(responseMeasured = as.character(responseMeasured)) %>% 
  filter(!is.na(responseMeasured)) %>% 
  data.frame()


densResponseType = MCMCglmm(yi~responseMeasured-1,
                            data=tmpDens,
                            random = ~scientificName + commonName + citation + idh(SE):units,
                            ginverse = list(scientificName=invB$Ainv),
                            prior=prior,
                            nitt=400000,
                            burnin=200000,
                            thin=100)

occResponseType = MCMCglmm(yi~responseMeasured-1,
                           data=tmpOcc,
                           random = ~scientificName + commonName + citation + idh(SE):units,
                           ginverse = list(scientificName=invB$Ainv),
                           prior=prior,
                           nitt=400000,
                           burnin=200000,
                           thin=100)

#Checking sample sizes and convergence
# summary(densResponseType)
# plot(densResponseType$Sol)
# plot(densResponseType$VCV)
# 
# summary(occResponseType)
# plot(occResponseType$Sol)
# plot(occResponseType$VCV)


results = rbind(results, collectResults(densResponseType, occResponseType,
                                        "Resp. measured", "responseMeasured", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densResponseType, occResponseType,
                                             "Resp. measured", "responseMeasured"))





####
#Monitoring duration
####
tmpDens = dens %>% 
  mutate(monitoringDuration = scale(as.numeric(monitoringDuration), center=T, scale=T)) %>% 
  filter(!is.na(monitoringDuration)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(monitoringDuration = scale(as.numeric(monitoringDuration), center=T, scale=T)) %>% 
  filter(!is.na(monitoringDuration)) %>% 
  data.frame()

densMonitoringDur = MCMCglmm(yi~monitoringDuration,
                             data=tmpDens,
                             random = ~scientificName + commonName + citation + idh(SE):units,
                             ginverse = list(scientificName=invB$Ainv),
                             prior=prior,
                             nitt=400000,
                             burnin=200000,
                             thin=100)

occMonitoringDur = MCMCglmm(yi~monitoringDuration,
                            data=tmpOcc,
                            random = ~scientificName + commonName + citation + idh(SE):units,
                            ginverse = list(scientificName=invB$Ainv),
                            prior=prior,
                            nitt=400000,
                            burnin=200000,
                            thin=100)

#Checking sample sizes and convergence
# summary(densMonitoringDur)
# plot(densMonitoringDur$Sol)
# plot(densMonitoringDur$VCV)
# 
# summary(occMonitoringDur)
# plot(occMonitoringDur$Sol)
# plot(occMonitoringDur$VCV)

results = rbind(results, collectResults(densMonitoringDur, occMonitoringDur,
                                        "Monit. duration", "monitoringDuration", "continuous"))


####
#Playback density
####
tmpDens = dens %>% 
  mutate(playbackDensity = scale(as.numeric(playbackDensity), center=T, scale=T)) %>% 
  filter(!is.na(playbackDensity)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(playbackDensity = scale(as.numeric(playbackDensity), center=T, scale=T)) %>% 
  filter(!is.na(playbackDensity)) %>% 
  data.frame()


densPBdens = MCMCglmm(yi~playbackDensity,
                      data=tmpDens,
                      random = ~scientificName + commonName + citation + idh(SE):units,
                      ginverse = list(scientificName=invB$Ainv),
                      prior=prior,
                      nitt=400000,
                      burnin=200000,
                      thin=100)

occPBdens = MCMCglmm(yi~playbackDensity,
                     data=tmpOcc,
                     random = ~scientificName + commonName + citation + idh(SE):units,
                     ginverse = list(scientificName=invB$Ainv),
                     prior=prior,
                     nitt=400000,
                     burnin=200000,
                     thin=100)

#Checking sample sizes and convergence
#summary(densPBdens)
#plot(densPBdens$Sol)
#plot(densPBdens$VCV)

#summary(occPBdens)
#plot(occPBdens$Sol)
#plot(occPBdens$VCV)


results = rbind(results, collectResults(densPBdens, occPBdens,
                                        "Playback density", "playbackDensity", "continuous"))


####
#Playback duration
####
tmpDens = dens %>% 
  mutate(pbDuration = scale(as.numeric(pbDuration), center=T, scale=T)) %>% 
  filter(!is.na(pbDuration)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(pbDuration = scale(as.numeric(pbDuration), center=T, scale=T)) %>% 
  filter(!is.na(pbDuration)) %>% 
  data.frame()


densPBduration = MCMCglmm(yi~pbDuration,
                          data=tmpDens,
                          random = ~scientificName + commonName + citation + idh(SE):units,
                          ginverse = list(scientificName=invB$Ainv),
                          prior=prior,
                          nitt=400000,
                          burnin=200000,
                          thin=100)

occPBduration = MCMCglmm(yi~pbDuration,
                         data=tmpOcc,
                         random = ~scientificName + commonName + citation + idh(SE):units,
                         ginverse = list(scientificName=invB$Ainv),
                         prior=prior,
                         nitt=400000,
                         burnin=200000,
                         thin=100)

#Checking sample sizes and convergence
#summary(densPBduration)
#plot(densPBduration$Sol)
#plot(densPBduration$VCV)

#summary(occPBduration)
#plot(occPBduration$Sol)
#plot(occPBduration$VCV)


results = rbind(results, collectResults(densPBduration, occPBduration,
                                        "Playback duration", "pbDuration", "continuous"))


####
#Cue type
####
tmpDens = dens %>% 
  mutate(cueType = as.character(cueType)) %>% 
  filter(!is.na(cueType)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(cueType = as.character(cueType)) %>% 
  filter(!is.na(cueType)) %>% 
  data.frame()


densCueType = MCMCglmm(yi~cueType-1,
                       data=tmpDens,
                       random = ~scientificName + commonName + citation + idh(SE):units,
                       ginverse = list(scientificName=invB$Ainv),
                       prior=prior,
                       nitt=400000,
                       burnin=200000,
                       thin=100)

occCueType = MCMCglmm(yi~cueType-1,
                      data=tmpOcc,
                      random = ~scientificName + commonName + citation + idh(SE):units,
                      ginverse = list(scientificName=invB$Ainv),
                      prior=prior,
                      nitt=400000,
                      burnin=200000,
                      thin=100)

#Checking sample sizes and convergence
#summary(densCueType)
#plot(densCueType$Sol)
#plot(densCueType$VCV)

#summary(occCueType)
#plot(occCueType$Sol)
#plot(occCueType$VCV)


results = rbind(results, collectResults(densCueType, occCueType,
                                        "Cue type", "cueType", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densCueType, occCueType, "Cue type",
                                             "cueType"))


####
#Signal origin
####
tmpDens = dens %>% 
  mutate(geographicOrigin = as.character(geographicOrigin)) %>% 
  filter(!is.na(geographicOrigin)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(geographicOrigin = as.character(geographicOrigin)) %>% 
  filter(!is.na(geographicOrigin)) %>% 
  data.frame()


densOrigin = MCMCglmm(yi~geographicOrigin-1,
                      data=tmpDens,
                      random = ~scientificName + commonName + citation + idh(SE):units,
                      ginverse = list(scientificName=invB$Ainv),
                      prior=prior,
                      nitt=400000,
                      burnin=200000,
                      thin=100)

occOrigin = MCMCglmm(yi~geographicOrigin-1,
                     data=tmpOcc,
                     random = ~scientificName + commonName + citation + idh(SE):units,
                     ginverse = list(scientificName=invB$Ainv),
                     prior=prior,
                     nitt=400000,
                     burnin=200000,
                     thin=100)

#Checking sample sizes and convergence
#summary(densOrigin)
#plot(densOrigin$Sol)
#plot(densOrigin$VCV)

#summary(occOrigin)
#plot(occOrigin$Sol)
#plot(occOrigin$VCV)


results = rbind(results, collectResults(densOrigin, occOrigin,
                                        "Origin", "geographicOrigin", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densOrigin, occOrigin, "Origin",
                                             "geographicOrigin"))

####
#Time of year PB applied
####
tmpDens = dens %>% 
  mutate(timeOfYear = as.character(timeOfYear)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(timeOfYear = as.character(timeOfYear)) %>% 
  data.frame()


densTOY = MCMCglmm(yi~timeOfYear-1,
                   data=tmpDens,
                   random = ~scientificName + commonName + citation + idh(SE):units,
                   ginverse = list(scientificName=invB$Ainv),
                   prior=prior,
                   nitt=400000,
                   burnin=200000,
                   thin=100)

occTOY = MCMCglmm(yi~timeOfYear-1,
                  data=tmpOcc,
                  random = ~scientificName + commonName + citation + idh(SE):units,
                  ginverse = list(scientificName=invB$Ainv),
                  prior=prior,
                  nitt=400000,
                  burnin=200000,
                  thin=100)

#Checking sample sizes and convergence
#summary(densTOY)
#plot(densTOY$Sol)
#plot(densTOY$VCV)

#summary(occTOY)
#plot(occTOY$Sol)
#plot(occTOY$VCV)

results = rbind(results, collectResults(densTOY, occTOY,
                                        "Time of year", "timeOfYear", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densTOY, occTOY, "Time of year",
                                             "timeOfYear"))

####
#Used decoys
####
tmpDens = dens %>% 
  mutate(decoy = as.character(decoy)) %>% 
  filter(!is.na(decoy)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(decoy = as.character(decoy)) %>% 
  filter(!is.na(decoy)) %>% 
  data.frame()

densDecoy = MCMCglmm(yi~decoy-1,
                     data=tmpDens,
                     random = ~scientificName + commonName + citation + idh(SE):units,
                     ginverse = list(scientificName=invB$Ainv),
                     prior=prior,
                     nitt=400000,
                     burnin=200000,
                     thin=100)

occDecoy = MCMCglmm(yi~decoy-1,
                    data=tmpOcc,
                    random = ~scientificName + commonName + citation + idh(SE):units,
                    ginverse = list(scientificName=invB$Ainv),
                    prior=prior,
                    nitt=400000,
                    burnin=200000,
                    thin=100)

#Checking sample sizes and convergence
#summary(densDecoy)
#plot(densDecoy$Sol)
#plot(densDecoy$VCV)

#summary(occDecoy)
#plot(occDecoy$Sol)
#plot(occDecoy$VCV)

results = rbind(results, collectResults(densDecoy, occDecoy,
                                        "Decoys", "decoy", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densDecoy, occDecoy, "Decoys",
                                             "decoy"))

####
#Quality of conspecific producing signal
####
tmpDens = dens %>% 
  mutate(conspQuality = as.character(conspQuality)) %>% 
  filter(!is.na(conspQuality)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(conspQuality = as.character(conspQuality)) %>% 
  filter(!is.na(conspQuality)) %>% 
  data.frame()

densFitness = MCMCglmm(yi~conspQuality-1,
                       data=tmpDens,
                       random = ~scientificName + commonName + citation + idh(SE):units,
                       ginverse = list(scientificName=invB$Ainv),
                       prior=prior,
                       nitt=400000,
                       burnin=200000,
                       thin=100)

occFitness = MCMCglmm(yi~conspQuality-1,
                      data=tmpOcc,
                      random = ~scientificName + commonName + citation + idh(SE):units,
                      ginverse = list(scientificName=invB$Ainv),
                      prior=prior,
                      nitt=400000,
                      burnin=200000,
                      thin=100)

#Checking sample sizes and convergence
#summary(densFitness)
#plot(densFitness$Sol)
#plot(densFitness$VCV)

#summary(occFitness)
#plot(occFitness$Sol)
#plot(occFitness$VCV)

results = rbind(results, collectResults(densFitness, occFitness,
                                        "Consp. quality", "conspQuality", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densFitness, occFitness, "Consp. quality",
                                             "conspQuality"))


####
#Longevity
####
tmpDens = dens %>% 
  mutate(longevity = scale(as.numeric(longevity), center=T, scale=T)) %>% 
  filter(!is.na(longevity)) %>% 
  data.frame()
tmpOcc = occ %>% 
  mutate(longevity = scale(as.numeric(longevity), center=T, scale=T)) %>% 
  filter(!is.na(longevity)) %>% 
  data.frame()


densLong = MCMCglmm(yi~longevity,
                    data=tmpDens,
                    random = ~scientificName + commonName + citation + idh(SE):units,
                    ginverse = list(scientificName=invB$Ainv),
                    prior=prior,
                    nitt=400000,
                    burnin=200000,
                    thin=100)

occLong = MCMCglmm(yi~longevity,
                   data=tmpOcc,
                   random = ~scientificName + commonName + citation + idh(SE):units,
                   ginverse = list(scientificName=invB$Ainv),
                   prior=prior,
                   nitt=400000,
                   burnin=200000,
                   thin=100)

#Checking sample sizes and convergence
#summary(densLong)
#plot(densLong$Sol)
#plot(densLong$VCV)

#summary(occLong)
#plot(occLong$Sol)
#plot(occLong$VCV)


results = rbind(results, collectResults(densLong, occLong,
                                        "Longevity", "longevity", "continuous"))


####
#Migratory status
####
tmpDens = dens %>% 
  filter(!is.na(migratory)) %>% 
  data.frame()
tmpOcc = occ %>% 
  filter(!is.na(migratory)) %>% 
  data.frame()


densMig = MCMCglmm(yi~migratory-1,
                   data=tmpDens,
                   random = ~scientificName + commonName + citation + idh(SE):units,
                   ginverse = list(scientificName=invB$Ainv),
                   prior=prior,
                   nitt=400000,
                   burnin=200000,
                   thin=100)

occMig = MCMCglmm(yi~migratory-1,
                  data=tmpOcc,
                  random = ~scientificName + commonName + citation + idh(SE):units,
                  ginverse = list(scientificName=invB$Ainv),
                  prior=prior,
                  nitt=400000,
                  burnin=200000,
                  thin=100)

#Checking sample sizes and convergence
#summary(densMig)
#plot(densMig$Sol)
#plot(densMig$VCV)

#summary(occMig)
#plot(occMig$Sol)
#plot(occMig$VCV)


results = rbind(results, collectResults(densMig, occMig,
                                        "Migratory", "migratory", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densMig, occMig, "Migratory", "migratory"))

####
#Polygamy
####
tmpDens = dens %>% 
  filter(!is.na(polygamy)) %>% 
  data.frame()
tmpOcc = occ %>% 
  filter(!is.na(polygamy)) %>% 
  data.frame()


densPolygamy = MCMCglmm(yi~polygamy-1,
                        data=tmpDens,
                        random = ~scientificName + commonName + citation + idh(SE):units,
                        ginverse = list(scientificName=invB$Ainv),
                        prior=prior,
                        nitt=400000,
                        burnin=200000,
                        thin=100)

occPolygamy = MCMCglmm(yi~polygamy-1,
                       data=tmpOcc,
                       random = ~scientificName + commonName + citation + idh(SE):units,
                       ginverse = list(scientificName=invB$Ainv),
                       prior=prior,
                       nitt=400000,
                       burnin=200000,
                       thin=100)

#Checking sample sizes and convergence
#summary(densPolygamy)
#plot(densPolygamy$Sol)
#plot(densPolygamy$VCV)

#summary(occPolygamy)
#plot(occPolygamy$Sol)
#plot(occPolygamy$VCV)


results = rbind(results, collectResults(densPolygamy, occPolygamy,
                                        "Polygamy", "polygamy", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densPolygamy, occPolygamy, "Polygamy",
                                             "polygamy"))


####
#Clustered breeding
####
tmpDens = dens %>% 
  filter(!is.na(clusteredBreeder)) %>% 
  data.frame()
tmpOcc = occ %>% 
  filter(!is.na(clusteredBreeder)) %>% 
  data.frame()


densClustered = MCMCglmm(yi~clusteredBreeder-1,
                         data=tmpDens,
                         random = ~scientificName + commonName + citation + idh(SE):units,
                         ginverse = list(scientificName=invB$Ainv),
                         prior=prior,
                         nitt=400000,
                         burnin=200000,
                         thin=100)

occClustered = MCMCglmm(yi~clusteredBreeder-1,
                        data=tmpOcc,
                        random = ~scientificName + commonName + citation + idh(SE):units,
                        ginverse = list(scientificName=invB$Ainv),
                        prior=prior,
                        nitt=400000,
                        burnin=200000,
                        thin=100)

#Checking sample sizes and convergence
#summary(densClustered)
#plot(densClustered$Sol)
#plot(densClustered$VCV)

#summary(occClustered)
#plot(occClustered$Sol)
#plot(occClustered$VCV)


results = rbind(results, collectResults(densClustered, occClustered,
                                        "Clustered breeder", "clusteredBreeder", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densClustered, occClustered, "Clustered breeder",
                                             "clusteredBreeder"))

####
#Habitat preference
####
tmpDens = dens %>% 
  filter(!is.na(habitat)) %>% 
  data.frame()
tmpOcc = occ %>% 
  filter(!is.na(habitat)) %>% 
  data.frame()

densEphemeral = MCMCglmm(yi~habitat-1,
                         data=tmpDens,
                         random = ~scientificName + commonName + citation + idh(SE):units,
                         ginverse = list(scientificName=invB$Ainv),
                         prior=prior,
                         nitt=400000,
                         burnin=200000,
                         thin=100)

occEphemeral = MCMCglmm(yi~habitat-1,
                        data=tmpOcc,
                        random = ~scientificName + commonName + citation + idh(SE):units,
                        ginverse = list(scientificName=invB$Ainv),
                        prior=prior,
                        nitt=400000,
                        burnin=200000,
                        thin=100)

#Checking sample sizes and convergence
#summary(densEphemeral)
#plot(densEphemeral$Sol)
#plot(densEphemeral$VCV)

#summary(occEphemeral)
#plot(occEphemeral$Sol)
#plot(occEphemeral$VCV)


results = rbind(results, collectResults(densEphemeral, occEphemeral,
                                        "Habitat", "habitat", "categorical"))

resultsComp = rbind(resultsComp, comparisons(densEphemeral, occEphemeral, "Habitat",
                                             "habitat"))



#--------------------------------------------------------------------
#Making moderators figure
#--------------------------------------------------------------------
tmp = results %>%
  filter(covariate=='fixed' & variable!='(Intercept)')

#A little cleanup
tmp$variable[which(tmp$model=='Clustered breeder' & tmp$variable=='Clustered')] = 'Yes'
tmp$variable[which(tmp$model=='Clustered breeder' & tmp$variable=='NoEvidence')] = 'No'
tmp$variable[which(tmp$model=='Cue type' & tmp$variable=='PerformanceInformation')] = 'PI'
tmp$variable[which(tmp$model=='Resp. measured' & tmp$variable=='Detected')] = 'Det.'
tmp$variable[which(tmp$model=='Resp. measured' & tmp$variable=='Territorial')] = 'Territ.'
tmp$variable[which(tmp$model=='Resp. measured' & tmp$variable=='Breeding')] = 'Breed.'
tmp$variable[which(tmp$model=='Origin' & tmp$variable=='Unknown')] = 'Unk.'
tmp$variable[which(tmp$model=='Origin' & tmp$variable=='Regional')] = 'Reg.'
tmp$variable[which(tmp$model=='Previous occ.' & tmp$variable=='Unknown')] = 'Unk.'
tmp$variable[which(tmp$model=='Previous occ.' & tmp$variable=='Unoccupied')] = 'Unocc.'
tmp$variable[which(tmp$model=='Previous occ.' & tmp$variable=='Occupied')] = 'Occ.'
tmp$variable[which(tmp$model=='Decoys' & tmp$variable=='NotUsed')] = 'Not used'


tmp$model = ifelse(tmp$model=='Longevity', 'A) Longevity', tmp$model)
tmp$model = ifelse(tmp$model=='Migratory', 'B) Migratory', tmp$model)
tmp$model = ifelse(tmp$model=='Polygamy', 'C) Polygamy', tmp$model)
tmp$model = ifelse(tmp$model=='Clustered breeder', 'D) Clustered breeder', tmp$model)
tmp$model = ifelse(tmp$model=='Habitat', 'E) Habitat', tmp$model)
tmp$model = ifelse(tmp$model=='Size', 'F) Size', tmp$model)
tmp$model = ifelse(tmp$model=='Latitude', 'G) Latitude', tmp$model)
tmp$model = ifelse(tmp$model=='Local dist.', 'H) Local dist.', tmp$model)
tmp$model = ifelse(tmp$model=='Previous occ.', 'I) Previous occ.', tmp$model)
tmp$model = ifelse(tmp$model=='Site quality', 'J) Site quality', tmp$model)
tmp$model = ifelse(tmp$model=='Playback duration', 'K) Playback duration', tmp$model)
tmp$model = ifelse(tmp$model=='Playback density', 'L) Playback density', tmp$model)
tmp$model = ifelse(tmp$model=='Time of year', 'M) Time of year', tmp$model)
tmp$model = ifelse(tmp$model=='Origin', 'N) Origin', tmp$model)
tmp$model = ifelse(tmp$model=='Cue type', 'O) Cue type', tmp$model)
tmp$model = ifelse(tmp$model=='Decoys', 'P) Decoys', tmp$model)
tmp$model = ifelse(tmp$model=='Consp. quality', 'Q) Consp. quality', tmp$model)
tmp$model = ifelse(tmp$model=='Monit. duration', 'R) Monit. duration', tmp$model)
tmp$model = ifelse(tmp$model=='Resp. measured', 'S) Resp. measured', tmp$model)
tmp$response = ifelse(tmp$response=='occupancy', 'Presence-absence', 'Density/abundance')
tmp = tmp %>%
  filter(variable!='NA') %>% 
  mutate(variable = factor(variable, levels=c('', 'Migrant', 'Resident',
                                              'Gradient',
                                              'Low', 'High', 'No', 'Yes',
                                              'Ephemeral', 'Static', 'Unk.',
                                              'Unocc.', 'Occ.', 'Post',
                                              'Pre/breeding', 'Reg.',
                                              'Local', 'Location', 'PI',
                                              'Not used', 'Used', 'Det.',
                                              'Territ.', 'Breed.')))

a = ggplot(tmp, aes(x=variable, y=post.mean, shape=response))+
  facet_wrap(~model, scales='free_y')+
  geom_point(position=position_dodge(0.5), size=3)+
  theme_bw()+
  geom_abline(intercept=0, slope=0, linetype='dashed')+
  geom_errorbar(aes(ymin=l.95..CI, ymax=u.95..CI), position=position_dodge(0.5), width=0)+
  theme(panel.grid=element_blank())+
  coord_flip()+
  theme(axis.text.y=element_text(angle=90, hjust=0.5))+
  theme(axis.ticks.y=element_blank())+
  theme(axis.title.y=element_blank())+
  ylab("Parameter estimate (95% credible interval)")+
  theme(strip.background=element_blank())+
  theme(legend.position='none')+
  theme(legend.title=element_blank())

tmp1 = tmp %>% filter(model=='A) Longevity')
b = ggplot(tmp1, aes(x=variable, y=post.mean, shape=response))+
  theme_bw()+
  geom_point(position=position_dodge(0.5), size=3)+
  theme(legend.title=element_blank(), legend.background=element_blank())

ggdraw(a) +
  draw_plot(get_legend(b), 0.8, 0.05, 0.2, 0.2)



#---------------------------------------------------------------------
#Level comparisons for categorical moderators
#---------------------------------------------------------------------
tmp = resultsComp %>% 
  filter(var1!='NA')
colnames(tmp) = c("Estimate", "LCL", "UCL", "Response", "Moderator",
                  "Level 1", "Level 2", "n 1", "n 2", "Significant")
tmp = tmp %>%
  select("Moderator", "Response", "Level 1", "Level 2", "n 1", "n 2",
         "Estimate", "LCL", "UCL") %>% 
  mutate(Moderator = replace(Moderator, Moderator=='Migratory', 'Migratory status'),
         Moderator = replace(Moderator, Moderator=='Habitat', 'Habitat preference'),
         Moderator = replace(Moderator, Moderator=="Previous occ.", "Previous occupation"),
         Moderator = replace(Moderator, Moderator=="Origin", "Signal origin"),
         Moderator = replace(Moderator, Moderator=='Consp. quality', "Conspecific quality"),
         Moderator = replace(Moderator, Moderator=='Resp. measured', "Response measured")) %>% 
  mutate(Moderator = factor(Moderator, levels=c("Migratory status", "Polygamy",
                                                "Clustered breeder", "Habitat preference",
                                                "Previous occupation", "Site quality",
                                                "Time of year", "Signal origin",
                                                "Cue type", "Decoys",
                                                "Conspecific quality", "Response measured"))) %>% 
  arrange(Moderator, Response) %>% 
  mutate(Response = ifelse(Response=='density', 'Density/abundance', 'Presence-absence')) %>% 
  mutate(Estimate = round(Estimate, 2),
         LCL = round(LCL, 2),
         UCL = round(UCL, 2)) %>% 
  mutate_if(is.factor, as.character)

for(i in unique(tmp$Moderator)){
  a = tmp[1,]
  a[1,] = ''
  b = c = a
  b$Moderator = i
  c$Response = 'Density/abundance'
  tmp1 = tmp %>%
    filter(Moderator==i & Response=='Density/abundance') %>% 
    mutate(Moderator = '', Response = '')
  tmp1 = rbind(a, b, c, tmp1)
  tmp2 = tmp %>% 
    filter(Moderator==i & Response=='Presence-absence') %>% 
    mutate(Moderator = '', Response = '')
  c$Response = 'Presence-absence'
  tmp2 = rbind(c, tmp2)
  tmp1 = rbind(tmp1, tmp2)
  
  tmp3 = tmp %>%
    filter(Moderator!=i)
  
  tmp = rbind(tmp3, tmp1)
  
}

tmp
