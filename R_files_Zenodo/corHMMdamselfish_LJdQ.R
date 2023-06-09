library(corHMM)
library(dplyr)
library(readr)
library(parallel)
library(phytools)


setwd("")


#tree
phy <- read.tree("tree_damsels3.nwk")
#phy <- ladderize(phy, right=F)
plot(phy, cex=0.5)

phy.c <- read.tree("tree_damsels3_clad.nwk")
phy.c <- ladderize(phy.c, right=F)
plot(phy.c,cex=0.5)

phy.tips <- sort(phy$tip.label)
phy.c.tips <- sort(phy.c$tip.label)

tips <- cbind(phy.tips,phy.c.tips)

cophyloplot(phy, phy.c, tips, space = 100, gap = 100)


#traits
damsels <- read_delim("damsels3.txt", delim = "\\t", escape_double = FALSE, trim_ws = TRUE)
data <- data.frame(as.list(damsels))
head(data)
tail(data)

table(data$T2)

#data[is.na(data)] <- '?'
#View(damsels)

#Which names are in the tree but not in the data frame?
setdiff(phy$tip.label, data$Genus_sp) 
#Which names are in the dataframe but not in the tree?
setdiff(data$Genus_sp, phy$tip.label)

#Checking the number of states for each trait
unique(data$T1) # 2 states (0 no red, 1 red)
unique(data$T2)  # 3 states (0 zooplanktivore, 1 omnivore, 2 herbivore)
unique(data$T3)  # 4 states (lws expression, 0 <1, 1 1-10, 2 10-20, 3 >20)
unique(data$T5)  # 4 states (lws expression, 0 <1, 1 1-10, 2 10-20, 3 >20)


# #Checking for the  combinations observed between traits
# data$t1_t2 <- paste(data$T1,data$T2)
# data$t1_t2_t3 <- paste(data$T1, data$T2, data$T3)
# data$t1_t3 <- paste(data$T1,data$T3)
# data$t2_t3 <- paste(data$T2,data$T3)

#warning: by default, the method will take into consideration only the  combinations of characters observed in the dataset. See flag 'collapse' of the function 'corHMM'
#Here you see the combinations of traits present in this dataset:
# t1_t2 <- unique(data$t1_t2); length(t1_t2) #there are 9 combinations assuming traits 1 and 2
# t1_t2_t3 <- unique(data$t1_t2_t3); length(t1_t2_t3) #there are 15 combinations
# t1_t3 <- unique(data$t1_t3); length(t1_t3) #there are 10 combinations
# t2_t3 <- unique(data$t2_t3); length(t2_t3) #9 combinations


traits <- c("T1","T2","T3")
models <- c('ARD','SYM','ER')
n.iterations <- 100 #ideally at least 100

#creating a table to be used to store each analysis
base <- as.data.frame(
  matrix(
    data=NA,
    ncol=6,
    nrow=length(models)
  )
)

names(base) <- c('trait',
              'method',
              'loglik',
              'aic',
              'aicc',
              'rate.cat')


base

#reconstructing each trait independently
#no hidden rate categories

#T1

#Creatoing  list to store full results
corHMM.fullResults.T1 <- list()

i <- 1 #T1

t1 <- base

for (j in 1:length(models)) { #loop to reconstruct  each trait using three distinct methods

# running the model for trait i and model j:
  corHMM.fullResults.T1 <- append(
    corHMM.fullResults.T1, list(
        model <-  corHMM(phy = phy,
                    data = select(data, Genus_sp, traits[i]), 
                    rate.cat = 1,
                    nstarts = n.iterations, 
                    model=models[j],
                    get.tip.states = T,
                    n.cores = detectCores()-2,
                    node.states = "marginal")
      )
    )

#storing the results to be used later
#corHMM.fullResults <- append(corHMM.fullResults, list(MK_1trait$states))

#summarizing the results to help choosing the best reconstruction 
  t1[j,] <- c(traits[i], 
            models[j], 
            model$loglik, 
            model$AIC, 
            model$AICc,
            model$rate.cat)


pdf(paste0(traits[i], "_", models[j],".pdf"))
plotRECON(phy.c,
          corHMM.fullResults.T1[[j]]$states,
          cex=0.5,
          pie.cex=0.5,
          node.depth = 0.1,
          title = paste(traits[i], models[j]))

plotMKmodel(corHMM.fullResults.T1[[j]])

dev.off()


} 

#saving T1 results
write_rds(corHMM.fullResults.T1, "T1_reconstructions.rds")

#Summary of your analysis:
t1

#if you want to see the results of three analyses 

corHMM.fullResults.T1[[1]]
corHMM.fullResults.T1[[2]]
corHMM.fullResults.T1[[3]]

#tip states for T1, only best model

write.csv(corHMM.fullResults.T1[[2]]$tip.states, "tip_states_rec_T1_ER.csv")


######## T2

#Creating  list to store full results
corHMM.fullResults.T2 <- list()

i <- 2 #T2

traits[i]

t2 <- base

for (j in 1:length(models)) { #loop to reconstruct  each trait using three distinct methods
  
  # running the model for trait i and model j:
  corHMM.fullResults.T2 <- append(
    corHMM.fullResults.T2, list(
      model <-  corHMM(phy = phy,
                       data = select(data, Genus_sp, traits[i]), 
                       rate.cat = 1,
                       nstarts = n.iterations, 
                       model=models[j],
                       get.tip.states = T,
                       n.cores = detectCores()-2,
                       node.states = "marginal")
    )
  )
  
  #storing the results to be used later
  #corHMM.fullResults <- append(corHMM.fullResults, list(MK_1trait$states))
  
  #summarizing the results to help choosing the best reconstruction 
  t2[j,] <- c(traits[i], 
             models[j], 
             model$loglik, 
             model$AIC,
             model$AICc, 
             model$rate.cat)
  
  
  pdf(paste0(traits[i], "_", models[j],".pdf"))
  plotRECON(phy.c,
            corHMM.fullResults.T2[[j]]$states,
            cex=0.5,
            pie.cex=0.5,
            node.depth = 0.1,
            title = paste(traits[i], models[j]))
  
  plotMKmodel(corHMM.fullResults.T2[[j]])
  
  dev.off()
  
  
}

#saving T2 results
write_rds(corHMM.fullResults.T2, "T2_reconstructions.rds")

#Summary of your analysis:
t2

#tip states for T2, only best model

write.csv(corHMM.fullResults.T2[[3]]$tip.states, "tip_states_rec_T2_ER.csv")



##########T3

#Creatoing  list to store full results
corHMM.fullResults.T3 <- list()

i <- 3 #T3

traits[i]

t3 <- base

for (j in 1:length(models)) { #loop to reconstruct  each trait using three distinct methods
  
  # running the model for trait i and model j:
  corHMM.fullResults.T3 <- append(
    corHMM.fullResults.T3, list(
      model <-  corHMM(phy = phy,
                       data = select(data, Genus_sp, traits[i]), 
                       rate.cat = 1,
                       nstarts = n.iterations, 
                       model=models[j],
                       get.tip.states = T,
                       n.cores = detectCores()-2,
                       node.states = "marginal")
    )
  )
  
  #storing the results to be used later
  #corHMM.fullResults <- append(corHMM.fullResults, list(MK_1trait$states))
  
  #summarizing the results to help choosing the best reconstruction 
  t3[j,] <- c(traits[i], 
              models[j], 
              model$loglik, 
              model$AIC, 
              model$AICc, 
              model$rate.cat)
  
  
  pdf(paste0(traits[i], "_", models[j],".pdf"))
  plotRECON(phy.c,
            corHMM.fullResults.T3[[j]]$states,
            cex=0.5,
            pie.cex=0.5,
            node.depth = 0.1,
            title = paste(traits[i], models[j]))
  
  plotMKmodel(corHMM.fullResults.T3[[j]])
  
  dev.off()
  
}

#saving T3 results
write_rds(corHMM.fullResults.T3, "T3_reconstructions.rds")

#Summary of your analysis:
t3


#### Trait dependent models


######### T1 x T2

#Creating  list to store full results
corHMM.fullResults.T1.T2 <- list()


t1.2 <- base
i <- c(1,2)
traits[i]

for (j in 1:length(models)) { #loop to reconstruct  each trait using three distinct methods
  
  # running the model for trait i and model j:
  corHMM.fullResults.T1.T2 <- append(
    corHMM.fullResults.T1.T2, list(
      model <-  corHMM(phy = phy,
                       data = select(data, Genus_sp, traits[i]), 
                       rate.cat = 1,
                       nstarts = n.iterations, 
                       model=models[j],
                       get.tip.states = T,
                       n.cores = detectCores()-2,
                       node.states = "marginal")
    )
  )
  
  #storing the results to be used later
  #corHMM.fullResults <- append(corHMM.fullResults, list(MK_1trait$states))
  
  #summarizing the results to help choosing the best reconstruction 
  t1.2[j,] <- c(paste(traits[i],collapse=","), 
              models[j], 
              model$loglik, 
              model$AIC, 
              model$AICc,
              model$rate.cat)
  
  
  pdf(paste(paste(traits[i],collapse="_"), "_", models[j],".pdf"))
  plotRECON(phy.c,
            corHMM.fullResults.T1.T2[[j]]$states,
            cex=0.5,
            pie.cex=0.5,
            node.depth = 0.1,
            title = paste(traits[i], models[j]))
  
  plotMKmodel(corHMM.fullResults.T1.T2[[j]])
  
  dev.off()
  
}




#saving T1xT2 results
write_rds(corHMM.fullResults.T1.T2, "T1_T2_reconstructions.rds")

#Summary of your analysis:
t1.2


write.csv(corHMM.fullResults.T1.T2[[j]]$tip.states, file = "T1_T2_tip_reconstructions.csv")

corHMM.fullResults.T1.T2[[3]]$tip.states


######### T1 x T3

#Creating  list to store full results
corHMM.fullResults.T1.T3 <- list()


t1.3 <- base
i <- c(1,3)
traits[i]

for (j in 1:length(models)) { #loop to reconstruct  each trait using three distinct methods
  
  # running the model for trait i and model j:
  corHMM.fullResults.T1.T3 <- append(
    corHMM.fullResults.T1.T3, list(
      model <-  corHMM(phy = phy,
                       data = select(data, Genus_sp, traits[i]), 
                       rate.cat = 1,
                       nstarts = n.iterations, 
                       model=models[j],
                       get.tip.states = T,
                       n.cores = detectCores()-2,
                       node.states = "marginal")
    )
  )
  
  #storing the results to be used later
  #corHMM.fullResults <- append(corHMM.fullResults, list(MK_1trait$states))
  
  #summarizing the results to help choosing the best reconstruction 
  t1.3[j,] <- c(paste(traits[i],collapse=","), 
                models[j], 
                model$loglik, 
                model$AIC, 
                model$AICc,
                model$rate.cat)
  
  
  pdf(paste(paste(traits[i],collapse="_"), "_", models[j],".pdf"))
  plotRECON(phy.c,
            corHMM.fullResults.T1.T3[[j]]$states,
            cex=0.5,
            pie.cex=0.5,
            node.depth = 0.1,
            title = paste(traits[i], models[j]))
  
  plotMKmodel(corHMM.fullResults.T1.T3[[j]])
  
  dev.off()
  
}

#saving T1xT3 results
write_rds(corHMM.fullResults.T1.T3, "T1_T3_reconstructions.rds")

#Summary of your analysis:
t1.3


######### T2 x T3

#Creating  list to store full results
corHMM.fullResults.T2.T3 <- list()


t2.3 <- base
i <- c(2,3)
traits[i]

for (j in 1:length(models)) { #loop to reconstruct  each trait using three distinct methods
  
  # running the model for trait i and model j:
  corHMM.fullResults.T2.T3 <- append(
    corHMM.fullResults.T2.T3, list(
      model <-  corHMM(phy = phy,
                       data = select(data, Genus_sp, traits[i]), 
                       rate.cat = 1,
                       nstarts = n.iterations, 
                       model=models[j],
                       get.tip.states = T,
                       n.cores = detectCores()-2,
                       node.states = "marginal")
    )
  )
  
  #storing the results to be used later
  #corHMM.fullResults <- append(corHMM.fullResults, list(MK_1trait$states))
  
  #summarizing the results to help choosing the best reconstruction 
  t2.3[j,] <- c(paste(traits[i],collapse=","), 
                models[j], 
                model$loglik, 
                model$AIC, 
                model$AICc, 
                
                model$rate.cat)
  
  
  pdf(paste(paste(traits[i],collapse="_"), "_", models[j],".pdf"))
  plotRECON(phy.c,
            corHMM.fullResults.T2.T3[[j]]$states,
            cex=0.5,
            pie.cex=0.5,
            node.depth = 0.1,
            title = paste(traits[i], models[j]))
  
  plotMKmodel(corHMM.fullResults.T2.T3[[j]])
  
  dev.off()
  
}

#saving T2xT3 results
write_rds(corHMM.fullResults.T2.T3, "T2_T3_reconstructions.rds")

#Summary of your analysis:
t2.3




######### T1 x T2 x T3

#Creating  list to store full results
corHMM.fullResults.T1.T2.T3 <- list()

t1.2.3 <- base
i <- c(1,2,3)
traits[i]

for (j in 1:length(models)) { #loop to reconstruct  each trait using three distinct methods
  
  # running the model for trait i and model j:
  corHMM.fullResults.T1.T2.T3 <- append(
    corHMM.fullResults.T1.T2.T3, list(
      model <-  corHMM(phy = phy,
                       data = select(data, Genus_sp, traits[i]), 
                       rate.cat = 1,
                       nstarts = n.iterations, 
                       model=models[j],
                       get.tip.states = T,
                       n.cores = detectCores()-1,
                       node.states = "marginal")
    )
  )
  
  #storing the results to be used later
  #corHMM.fullResults <- append(corHMM.fullResults, list(MK_1trait$states))
  
  #summarizing the results to help choosing the best reconstruction 
  t1.2.3[j,] <- c(paste(traits[i],collapse=","), 
                models[j], 
                model$loglik, 
                model$AIC, 
                model$rate.cat)
  
  
  pdf(paste(paste(traits[i],collapse="_"), "_", models[j],".pdf"))
  plotRECON(phy.c,
            corHMM.fullResults.T1.T2.T3[[j]]$states,
            cex=0.5,
            pie.cex=0.5,
            node.depth = 0.1,
            title = paste(traits[i], models[j]))
  
  plotMKmodel(corHMM.fullResults.T1.T2.T3[[j]])
  
  dev.off()
  
}

#saving T1xT3 results
write_rds(corHMM.fullResults.T1.T2.T3, "T1_T2_T3_reconstructions.rds")

#Summary of your analysis:
t1.2.3





### Uploading res



### Uploading res
setwd("/Users/jardimlu/Dropbox/Eawag/sara_stieb_corHMM/newData_20062022/corHMM")

T1 <- readRDS('T1_reconstructions.rds')

T1[[1]]$AICc
T1[[2]]$AICc
T1[[3]]$AICc

damsels.df <- as.data.frame(damsels)
traits.T1 <- damsels.df[,1:2]
traits.T1$state <- NA

#traits.T1$state <- ifelse(traits.T1$T1 == 0, "white", "red")

# Specifying the colours to use
red <-  "red"
non_red = "white"



pdf("figures_to_SS/supp/T1_color_ER_SYM_ARD.pdf")


plotRECON.ljq(phy.c, 
              T1[[1]]$states,
              cex=0.35, cex.lab=0.5,
              pie.cex=0.6,
              lwd = 10,
              label.offset = 0.27, adj = 0,
              piecolors = c("white","red"), edge.width = 0.5, pch = 19)

title ("Color, ARD", cex.main = 0.75)

tiplabels(pie = T1[[1]]$tip.states,
          tip = traits$tip.order,
          piecol=c(non_red,red),
          adj = .5, cex = .45)

legend("topleft", legend = c("Non-red","Red"), fill = c("white","red"),
       box.lwd = -1, bg = "white", cex = 0.5)

plotRECON.ljq(phy.c,
              T1[[2]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = c("white","red"), edge.width = 0.5)

title ("Color, SYM", cex.main = 0.75)

legend("topleft", legend = c("Non-red","Red"), fill = c("white","red"),
       box.lwd = -1, bg = "white", cex = 0.5)

tiplabels(pie = T1[[2]]$tip.states,
          #tip = traits$tip.order,
          piecol=c(non_red,red),
          adj = .5, cex = .45)


plotRECON.ljq(phy.c,
              T1[[3]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = c("white","red"), edge.width = 0.5)

legend("topleft", legend = c("Non-red","Red"), fill = c("white","red"),
       box.lwd = -1, bg = "white", cex = 0.5)

title ("Color, ER", cex.main = 0.75)

tiplabels(pie = T1[[3]]$tip.states,
          #tip = traits$tip.order,
          piecol=c(non_red,red),
          adj = .5, cex = .45)


dev.off()

#T2

T2 <- readRDS('T2_reconstructions.rds')

T2[[1]]$AICc
T2[[2]]$AICc
T2[[3]]$AICc


# Specifying the colours to use
#colorblind friendly

library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)

colorblind <- brewer.pal(3,'Reds')

zooplank <-  "white"
omniv <-  "orange"
herb <- "red"

colors <- c('white','orange','red')

states <- c("Zooplanktivores","Omnivores", "Herbivores")



pdf("figures_to_SS/supp/T2_TrophicGroups_ER_SYM_ARD.pdf")


plotRECON.ljq(phy.c, 
              T2[[1]]$states,
              cex=0.35, cex.lab=0.5,
              pie.cex=0.6,
              lwd = 10,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5, pch = 19)

title ("Trophic groups, ARD", cex.main = 0.75)

T2[[1]]$tip.states[c(1:2),] <- NA
T2[[2]]$tip.states[c(1:2),] <- NA
T2[[3]]$tip.states[c(1:2),] <- NA

tiplabels(pie = T2[[1]]$tip.states,
          tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.5)


plotRECON.ljq(phy.c,
              T2[[2]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5)

title ("Trophic groups, SYM", cex.main = 0.75)
legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.5)
tiplabels(pie = T2[[2]]$tip.states,
          #tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)


plotRECON.ljq(phy.c,
              T2[[3]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5)
legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.5)
title ("Trophic groups, ER", cex.main = 0.75)

tiplabels(pie = T2[[3]]$tip.states,
          #tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)


dev.off()

#T3

T3 <- readRDS('T3_reconstructions.rds')

T3[[1]]$AICc
T3[[2]]$AICc
T3[[3]]$AICc


T3[[1]]$AIC
T3[[2]]$AIC
T3[[3]]$AIC


# Specifying the colours to use

colors <- c('white',"yellow",'red','darkred')
states <- c("0-0.9%", "1-9.9%", "10-19.9%","20-31%")



pdf("figures_to_SS/supp/T3_GenExpr_ER_SYM_ARD.pdf")


plotRECON.ljq(phy.c, 
              T3[[1]]$states,
              cex=0.35, cex.lab=0.5,
              pie.cex=0.6,
              lwd = 10,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5, pch = 19)

title ("LWS expression, ARD", cex.main = 0.75)

T3[[1]]$tip.states[c(1:2),] <- NA
T3[[2]]$tip.states[c(1:2),] <- NA
T3[[3]]$tip.states[c(1:2),] <- NA

tiplabels(pie = T3[[1]]$tip.states,
          tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.5)


plotRECON.ljq(phy.c,
              T3[[2]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5)

title ("LWS expression, SYM", cex.main = 0.75)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.5)

tiplabels(pie = T3[[2]]$tip.states,
          #tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)


plotRECON.ljq(phy.c,
              T3[[3]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5)

title ("LWS expression, ER", cex.main = 0.75)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.5)

tiplabels(pie = T3[[3]]$tip.states,
          #tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)


dev.off()


#T1.T2

T1.T2<- readRDS('T1_T2_reconstructions.rds')

T1.T2[[1]]$AICc
T1.T2[[2]]$AICc
T1.T2[[3]]$AICc


# Specifying the colours to use

colfunc <- colorRampPalette(c("white", "yellow", "red", "black"))
colfunc(6)
plot(rep(1,6), col=colfunc(6),pch=15,cex=5.6)

colors <- colfunc(6)

states <- c("Non-red / Zooplanktivores",
            "Non-red / Omnivores",
            "Non-red / Herbivores",
            "Red / Zooplanktivores",
            "Red / Omnivores",
            "Red / Herbivores")


pdf("figures_to_SS/supp/T1.T2_GenExpr_ER_SYM_ARD.pdf")


plotRECON.ljq(phy.c, 
              T1.T2[[1]]$states,
              cex=0.35, cex.lab=0.5,
              pie.cex=0.6,
              lwd = 10,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5, pch = 19)

title ("Color & Trophic groups, ARD", cex.main = 0.75)

T1.T2[[1]]$tip.states[c(1:2),] <- NA
T1.T2[[2]]$tip.states[c(1:2),] <- NA
T1.T2[[3]]$tip.states[c(1:2),] <- NA

tiplabels(pie = T1.T2[[1]]$tip.states,
          tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.4)


plotRECON.ljq(phy.c,
              T1.T2[[2]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5)

title ("Color & Trophic groups, SYM", cex.main = 0.75)


legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.4)

tiplabels(pie = T1.T2[[2]]$tip.states,
          #tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)


plotRECON.ljq(phy.c,
              T1.T2[[3]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5)

title ("Color & Trophic groups, ER", cex.main = 0.75)


legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.4)

tiplabels(pie = T1.T2[[3]]$tip.states,
          #tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)

dev.off()

#T1.T3

T1.T3<- readRDS('T1_T3_reconstructions.rds')

T1.T3[[1]]$AICc
T1.T3[[2]]$AICc
T1.T3[[3]]$AICc


# Specifying the colours to use

colfunc <- colorRampPalette(c("white", "yellow", "red", "black"))
colfunc(8)
plot(rep(1,8), col=colfunc(8),pch=15,cex=5.6)

colors <- colfunc(8)

states <- c("Non-red/ 0-0.9%",
            "Non-red / 1-9.9%",
            "Non-red / 10-19.9%",
            "Non-red / 20-31%",
            
            "Red / 0-0.9%",
            "Red / 1-9.9%",
            "Red / 10-19.9%",
            "Red / 20-31%")




pdf("figures_to_SS/supp/T1.T3_Color_GenExpr_ER_SYM_ARD.pdf")


plotRECON.ljq(phy.c, 
              T1.T3[[1]]$states,
              cex=0.35, cex.lab=0.5, 
              pie.cex=0.6,
              lwd = 10,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5, pch = 19)

title ("Color & LWS expression, ARD", cex.main = 0.75)

T1.T3[[1]]$tip.states[c(1:2),] <- NA
T1.T3[[2]]$tip.states[c(1:2),] <- NA
T1.T3[[3]]$tip.states[c(1:2),] <- NA

tiplabels(pie = T1.T3[[1]]$tip.states,
          tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.4)


plotRECON.ljq(phy.c,
              T1.T3[[2]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5)

title ("Color & LWS expressions, SYM", cex.main = 0.75)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.4)

tiplabels(pie = T1.T3[[2]]$tip.states,
          #tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)


plotRECON.ljq(phy.c,
              T1.T3[[3]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5)

title ("Color & LWS expression, ER", cex.main = 0.75)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.4)

tiplabels(pie = T1.T3[[3]]$tip.states,
          #tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)

dev.off()


#T2.T3

T2.T3<- readRDS('T2_T3_reconstructions.rds')

T2.T3[[1]]$AICc
T2.T3[[2]]$AICc
T2.T3[[3]]$AICc


# Specifying the colours to use

colfunc <- colorRampPalette(c("yellow", "red", "black"))
colfunc(12)
plot(rep(1,12), col=colfunc(12),pch=15,cex=5.6)

colors <- colfunc(12)

states <- c("Zooplanktivores / 0-0.9%",
            "Zooplanktivores / 1-9.9%",
            "Zooplanktivores / 10-19.9%",
            "Zooplanktivores / 20-31%",
            
            "Omnivores / 0-0.9%",
            "Omnivores / 1-9.9%",
            "Omnivores / 10-19.9%",
            "Omnivores / 20-31%",
            
            "Herbivores / 0-0.9%",
            "Herbivores / 1-9.9%",
            "Herbivores / 10-19.9%",
            "Herbivores / 20-31%"
)




pdf("figures_to_SS/supp/T2.T3_Trophic_GenExpr_ER_SYM_ARD.pdf")


plotRECON.ljq(phy.c, 
              T2.T3[[1]]$states,
              cex=0.35, cex.lab=0.5, 
              pie.cex=0.6,
              lwd = 10,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5, pch = 19)

title ("Trophic groups & LWS expression, ARD", cex.main = 0.75)

T2.T3[[1]]$tip.states[c(1:2),] <- NA
T2.T3[[2]]$tip.states[c(1:2),] <- NA
T2.T3[[3]]$tip.states[c(1:2),] <- NA

tiplabels(pie = T2.T3[[1]]$tip.states,
          tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.34)


plotRECON.ljq(phy.c,
              T2.T3[[2]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5)

title ("Trophic groups & LWS expressions, SYM", cex.main = 0.75)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.34)

tiplabels(pie = T2.T3[[2]]$tip.states,
          #tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)


plotRECON.ljq(phy.c,
              T2.T3[[3]]$states,
              cex=0.35,
              pie.cex=0.6,
              node.depth = 0.1,
              label.offset = 0.27, adj = 0,
              piecolors = colors, edge.width = 0.5)

title ("Trophic groups & LWS expression, ER", cex.main = 0.75)

legend("topleft",  legend = states, fill = colors,
       box.lwd = -1, bg = "white", cex = 0.34)
tiplabels(pie = T2.T3[[3]]$tip.states,
          #tip = traits$tip.order,
          piecol=colors,
          adj = .5, cex = .45)

dev.off()


#Generating table with the parameters that were estimated by the models
base$delta.aic <- NA
base$delta.aicc <- NA

p.T1 <- base
for (i in 1:3) {
  run <- T1[[i]]
  p.T1$trait[i] <- 'Color'
  p.T1$method[i] <- models[i]
  p.T1$loglik[i] <- run$loglik
  p.T1$aic[i] <- run$AIC
  p.T1$aicc[i] <- run$AICc
  p.T1$rate.cat[i] <- run$rate.cat
}

p.T1$delta.aic <- round(p.T1$aic-min(p.T1$aic), digits = 2)
p.T1$delta.aicc <- round(p.T1$aicc-min(p.T1$aicc), digits = 2)



p.T2 <- base
for (i in 1:3) {
  run <- T2[[i]]
  p.T2$trait[i] <- 'Trophic groups'
  p.T2$method[i] <- models[i]
  p.T2$loglik[i] <- run$loglik
  p.T2$aic[i] <- run$AIC
  p.T2$aicc[i] <- run$AICc
  p.T2$rate.cat[i] <- run$rate.cat
}

p.T2$delta.aic <- round(p.T2$aic-min(p.T2$aic), digits = 2)
p.T2$delta.aicc <- round(p.T2$aicc-min(p.T2$aicc), digits = 2)



p.T3 <- base
for (i in 1:3) {
  run <- T3[[i]]  
  p.T3$trait[i] <- 'LWS expression'
  p.T3$method[i] <- models[i]
  p.T3$loglik[i] <- run$loglik
  p.T3$aic[i] <- run$AIC
  p.T3$aicc[i] <- run$AICc
  p.T3$rate.cat[i] <- run$rate.cat
}
p.T3$delta.aic <- round(p.T3$aic-min(p.T3$aic), digits = 2)
p.T3$delta.aicc <- round(p.T3$aicc-min(p.T3$aicc), digits = 2)



p.T1.T2 <- base
for (i in 1:3) {
  run <- T1.T2[[i]]  
  p.T1.T2$trait[i] <- 'Color & Trophic groups'
  p.T1.T2$method[i] <- models[i]
  p.T1.T2$loglik[i] <- run$loglik
  p.T1.T2$aic[i] <- run$AIC
  p.T1.T2$aicc[i] <- run$AICc
  p.T1.T2$rate.cat[i] <- run$rate.cat
}
p.T1.T2$delta.aic <- round(p.T1.T2$aic-min(p.T1.T2$aic), digits = 2)
p.T1.T2$delta.aicc <- round(p.T1.T2$aicc-min(p.T1.T2$aicc), digits = 2)


p.T1.T3 <- base
for (i in 1:3) {
  run <- T1.T3[[i]]  
  p.T1.T3$trait[i] <- 'Color & LWS expression'
  p.T1.T3$method[i] <- models[i]
  p.T1.T3$loglik[i] <- run$loglik
  p.T1.T3$aic[i] <- run$AIC
  p.T1.T3$aicc[i] <- run$AICc
  p.T1.T3$rate.cat[i] <- run$rate.cat
}
p.T1.T3$delta.aic <- round(p.T1.T3$aic-min(p.T1.T3$aic), digits = 2)
p.T1.T3$delta.aicc <- round(p.T1.T3$aicc-min(p.T1.T3$aicc), digits = 2)


p.T2.T3 <- base
for (i in 1:3) {
  run <- T2.T3[[i]]  
  p.T2.T3$trait[i] <- 'Trophic groups & LWS expression'
  p.T2.T3$method[i] <- models[i]
  p.T2.T3$loglik[i] <- run$loglik
  p.T2.T3$aic[i] <- run$AIC
  p.T2.T3$aicc[i] <- run$AICc
  p.T2.T3$rate.cat[i] <- run$rate.cat
}
p.T2.T3$delta.aic <- round(p.T2.T3$aic-min(p.T2.T3$aic), digits = 2)
p.T2.T3$delta.aicc <- round(p.T2.T3$aicc-min(p.T2.T3$aicc), digits = 2)


overall.table <- rbind(p.T1, p.T2, p.T3, p.T1.T2, p.T1.T3, p.T2.T3)

write.csv(overall.table, "summary.models.csv")


######
#### Modified function to avoid plotting legend
plotRECON.ljq <- function(phy, likelihoods, piecolors=NULL, cex=0.5, pie.cex=0.25, file=NULL, height=11, width=8.5, show.tip.label=TRUE, title=NULL, ...){
  if(is.null(piecolors)){
    piecolors=c("white","black","red","yellow","forestgreen","blue","coral","aquamarine","darkorchid","gold","grey","yellow","#3288BD","#E31A1C")
  }
  if(!is.null(file)){
    pdf(file, height=height, width=width,useDingbats=FALSE)
  }
  plot(phy, cex=cex, show.tip.label=show.tip.label, ...)
  
  if(!is.null(title)){
    title(main=title)
  }
  nodelabels(pie=likelihoods,piecol=piecolors, cex=pie.cex)
  states <- colnames(likelihoods)
  #legend(x="topleft", states, cex=0.8, pt.bg=piecolors,col="black",pch=21);
  
  if(!is.null(file)){
    dev.off()
  }
}









