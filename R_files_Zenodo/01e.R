library(ape)
library(picante)
library(betapart)
library(picante)
library(cowplot)
#Phylo-Ordinations
#per this tutorial: http://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html


phy <- read.tree("../DATA/phylogeny.analyzed.2016-01-05b.tre") 
phy$tip.label    #140
phy$tip.label[which(phy$tip.label == "Symphyotrichum_novae-angliae")] <-
  "Symphyotrichum_novae.angliae"
tree <- drop.tip(phy, which(!phy$tip.label %in% colnames(comm))) 
#convert phylogenety to a distance matrix
phy.dist <- cophenetic(tree)

#calculate phylogenetic MNTD beta diversity
head(df.PA) #from 00a.load
DF1 <-as.matrix((subset(df.PA, select = Agalinis_tenuifolia:Zizia_aurea)))
comm.mntd.dist <- comdistnt(DF1, phy.dist, abundance.weighted = FALSE) #P/A
comm.mntd.mds <- monoMDS(comm.mntd.dist, k=3)

###PLOT ORDINATIONS
#subset ordination by year and by PD treatment for plotting
as.data.frame(scores(comm.mntd.mds)) -> Sc.mntd.mds
Sc.mntd.mds$year <- df.PA$year 
Sc.mntd.mds$phyD <- as.factor(df.PA$phyD)
Sc.mntd.mds$plot <- df.PA$plot

start <- subset(Sc.mntd.mds, (year == "2016"))
spring <- subset(Sc.mntd.mds, (year == "2019.06"))
end <- subset(Sc.mntd.mds, (year == "2019.09"))
StoE <- rbind(start, spring, end)

Fig2B <-  StoE %>%             #sort ascending so that 2019 is plotted 
  ggplot() +
  geom_point(aes(x=MDS1, y=MDS2, shape=year, color = phyD), size=3) +
  theme_bw() +
  labs(shape="Year", color="PD") +
  geom_path(aes(x = MDS1, y = MDS2, group = plot), arrow = arrow(length = unit(0.15, "cm")), color="grey") +
  stat_ellipse(data=start, aes(x=MDS1,y=MDS2, group = phyD, color = phyD), linetype = 2) + 
  
  #geom_mark_ellipse(data=start, aes(x=MDS1,y=MDS2, group = phyD, color = phyD), linetype = 2) + 
  #stat_ellipse(data=spring, aes(x=NMDS1,y=NMDS2, group = phyD), color= "red") + 
  stat_ellipse(data=end, aes(x=MDS1,y=MDS2, color= phyD)) +
  
  #geom_mark_ellipse(data=end, aes(x=MDS1,y=MDS2, color= phyD)) +
  scale_color_manual("PD", values = c("L" = "#332288",
                                      "M" = "#117733",
                                      "H" = "#44AA99"
  )) +
  annotate("text", x = 1.3, y = 1.7, color="grey",
           label = "stress=0.19") + labs(x="NMDS1", y="NMDS2")

##############################
#PERMANOVA and dispersion test
#subset  by year and by PD treatment
head(df.PA) #from 00a.load
sub.df.PA <- subset(df.PA,  (year == "2016" |  year == "2019.06" | year == "2019.09" ))

DF1.sub <-as.matrix((subset(sub.df.PA, select = Agalinis_tenuifolia:Zizia_aurea)))
mntd.dist.sub <- comdistnt(DF1.sub, phy.dist, abundance.weighted = FALSE) #P/A

dis.ord <- adonis2(mntd.dist.sub ~ year + phyD, data = group.dat, strata="plot") #, permutations = perm) 
mod <- betadisper(mntd.dist.sub, group.dat$year)
mod.pd <- betadisper(mntd.dist.sub, group.dat$phyD)

modPO <- permutest(mod, pairwise = TRUE, permutations = 999) ## Permutation test for F

group.dat <- subset(sub.df.PA, select=c(year, phyD, plot))

groupsPD <- sub.df.PA %>% 
  unite(yearPD, c("year", phyD))
as.factor(groupsPD$yearPD) -> groupsPD$yearPD
mod2 <- betadisper(mntd.dist.sub, groupsPD$yearPD)
mod2P <- permutest(mod2, pairwise = TRUE, permutations = 999) ## Permutation test for F
mod2P

####Plot 
mod.pd <- data.frame(group=mod$group, distances=mod$distances)
mod.pd$group <- gsub( "2016_L", "Low planted", mod.pd$group)
mod.pd$group <- gsub( "2019.06_L", "Low year 3 spring", mod.pd$group)
mod.pd$group <- gsub( "2019.09_L", "Low year 3 summer", mod.pd$group)
mod.pd$group <- gsub( "2016_M", "Medium planted", mod.pd$group)
mod.pd$group <- gsub( "2019.06_M", "Medium year 3 spring", mod.pd$group)
mod.pd$group <- gsub( "2019.09_M", "Medium year 3 summer", mod.pd$group)
mod.pd$group <- gsub( "2016_H", "High planted", mod.pd$group)
mod.pd$group <- gsub( "2019.06_H", "High year 3 spring", mod.pd$group)
mod.pd$group <- gsub( "2019.09_H", "High year 3 summer", mod.pd$group)

mod.pd$group <- factor(mod.pd$group , levels=c("Low planted", "Low year 3 spring","Low year 3 summer","Medium planted", "Medium year 3 spring",  "Medium year 3 summer", "High planted", "High year 3 spring","High year 3 summer"))

ggplot(mod.pd, aes(x=group, y=distances, fill=group)) + geom_boxplot() +
  xlab("") +
  ylab("Distance to centroid") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  labs(fill="Survey", x= "Survey") + theme_bw() +
  scale_fill_manual("group", values = c("Low planted" = "#332288",
                                        "High planted" = "#44AA99", 
                                        "Low year 3 spring" = "#332288",
                                        "Low year 3 summer" = "#332288", 
                                        "High year 3 summer" = "#44AA99",
                                        "High year 3 spring" = "#44AA99", 
                                        "Medium planted" = "#117733",
                                        "Medium year 3 spring" = "#117733",
                                        "Medium year 3 summer" = "#117733"
  )) 



##Dissimilarity
#########################Dissimilarity by PhyD Class###################
temp = lapply(split(df.PA, df.PA$year), function(x) split(x, x$phyD))
nested.mntd.dis <- lapply(temp, function(x) lapply(x, function(x) comdistnt(subset(x, select = Agalinis_tenuifolia:Zizia_aurea), phy.dist, abundance.weighted = FALSE)))

nested.mean <- lapply(nested.mntd.dis, function(x) unlist(lapply(x, function(x) mean(x)) ) )

as.data.frame(nested.mean) -> diss
names(diss) <- names(nested.mean)
rownames(diss) -> diss$phyd
factor(diss$phyd, levels = c('L', 'M', 'H')) -> diss$phyd
melt(diss) -> diss.long
#test.fall <- subset(test, select = -c(`2018.05`, `2019.06`))
#test.spring <- subset(test, select = c(`2018.05`, `2019.06`))
#melt(test.fall) -> test.fall

P <- diss.long %>%
  mutate(group2 = case_when(variable %in% c("2018.05", "2019.06") ~ "A",
                            TRUE ~ "B")) %>% 
  ggplot(aes(x = variable, 
             y = value, 
             color = phyd, 
             group = interaction(phyd, group2), 
             linetype = group2)) +
  geom_line(size = 1) +
  scale_linetype_manual("season",
                        values = c("A" = "dashed",
                                   "B" = "solid"),
                        labels = c("A" = "spring",
                                   "B" = "fall")) + 
  scale_color_manual("PD",
                     values = c("L" = "#332288",
                                "M" = "#117733",
                                "H" = "#44AA99"
                     )) +
  geom_point(size = 3, alpha = 0.5) # only for comprehension, remove it
P <- P + labs(x = "year",  y= "dissimilarity (Jaccard's)", color="PD") + theme_bw() +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1)) +
  

#labs(title = "Mean dissimilarity",subtitle = "Presence/absence data", x = "year",  y= "Dissimilarity (Jaccard's)") + theme_bw() +
#ggsave(P, file="../OUT/meanDissPA.png", width = 8, height = 6)

#Confidence Intervals for dissimilarity
#temp from above bootstrap 
##
myReps = 1000
myCI = 0.95
myError <- myObs <- vector('list', length = length(temp))
names(myError) <- names(myObs) <- names(temp)
for(i in names(myError)) myError[[i]] <- list(L = numeric(0), M = numeric(0), H = numeric(0))

for(i in names(temp)) {
  for(j in c('L', 'M', 'H')) {
    print(paste('doing year', i, 'with phyD', j))
    bittyTemp <- temp[[i]][[j]][which(names(temp[[i]][[j]]) == 'Agalinis_tenuifolia'):which(names(temp[[i]][[j]]) == 'Zizia_aurea')]
    bittyTemp.d <- sapply(seq(myReps), function(x) {
      comdistnt(bittyTemp[sample(1:24, 12, replace = F), ], phy.dist) %>% mean
    })
    myObs[[i]][[j]] <- mean(comdistnt(bittyTemp, phy.dist))
    myError[[i]][[j]] <- c(obs = myObs[[i]][[j]], quantile(bittyTemp.d, c((1-myCI) / 2, myCI+((1-myCI) / 2)))) #half jacknifing in
  }
}
myError -> err
#lapply(err$`2016`, function(x) as.data.frame(melt(x))) -> B
as.data.frame(melt(err)) -> ERR #need to keep the obs, 2.5, 
valueIDs <- c("obs","lowCI", "highCI" )
cbind(ERR, valueIDs) -> NEWdf.test

names(NEWdf.test) <- c("value", "phyd", "year", "valueID")

data2 <- dcast(NEWdf.test, year + phyd ~ valueID, value.var = c("value"))
factor(data2$phyd, levels = c('L', 'M', 'H')) -> data2$phyd

PP <- data2 %>%
  mutate(group2 = case_when(year %in% c("2018.05", "2019.06") ~ "A",
                            TRUE ~ "B")) %>% 
  ggplot(aes(x = year, 
             y = obs, 
             color = phyd, 
             group = interaction(phyd, group2), 
             linetype = group2)) +
  geom_line(size = 1) +
  scale_linetype_manual("season",
                        values = c("A" = "dashed",
                                   "B" = "solid"),
                        labels = c("A" = "spring",
                                   "B" = "fall")) + 
  scale_color_manual("PD",
                     values = c("L" = "#332288",
                                "M" =  "#117733",
                                "H" =  "#44AA99"
                     )) +
  geom_point(size = 3, alpha = 0.5)  +
  geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=.2,
                position=position_dodge(.1)) 

PP <- PP + labs(x = "year",  y= "dissimilarity (MNTD)", color="PD") + theme_bw() +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))




cowplot::plot_grid(Fig2B, PP)

#############################FD treatments on Phylo ordinations##########################################################
#plot ordinations by FD
#subset ordination by year and by FD treatment for plotting
as.data.frame(scores(comm.mntd.mds)) -> Sc.mntd.mds
Sc.mntd.mds$year <- df.PA$year 
Sc.mntd.mds$fd <- as.factor(df.PA$fd)
Sc.mntd.mds$plot <- df.PA$plot

start <- subset(Sc.mntd.mds, (year == "2016"))
spring <- subset(Sc.mntd.mds, (year == "2019.06"))
end <- subset(Sc.mntd.mds, (year == "2019.09"))
StoE <- rbind(start, spring, end)

FigSd <-  StoE %>%             #sort ascending so that 2019 is plotted 
  ggplot() +
  geom_point(aes(x=MDS1, y=MDS2, shape=year, color = fd), size=3) +
  theme_bw() +
  labs(shape="Year", color="FD") +
  geom_path(aes(x = MDS1, y = MDS2, group = plot), arrow = arrow(length = unit(0.15, "cm")), color="grey") +
  stat_ellipse(data=start, aes(x=MDS1,y=MDS2, group = fd, color = fd), linetype = 2) + 
  
  #geom_mark_ellipse(data=start, aes(x=MDS1,y=MDS2, group = phyD, color = phyD), linetype = 2) + 
  #stat_ellipse(data=spring, aes(x=NMDS1,y=NMDS2, group = phyD), color= "red") + 
  stat_ellipse(data=end, aes(x=MDS1,y=MDS2, color= fd)) +
  
  #geom_mark_ellipse(data=end, aes(x=MDS1,y=MDS2, color= phyD)) +
  scale_color_manual("PD", values = c("L" = "#332288",
                                      "M" = "#117733",
                                      "H" = "#44AA99"
  )) +
  annotate("text", x = 1.3, y = 1.7, color="grey",
           label = "stress=0.19") + labs(x="NMDS1", y="NMDS2")


####Dissimilarity by FD###############
temp = lapply(split(df.PA, df.PA$year), function(x) split(x, x$fd))
nested.mntd.dis <- lapply(temp, function(x) lapply(x, function(x) comdistnt(subset(x, select = Agalinis_tenuifolia:Zizia_aurea), phy.dist, abundance.weighted = FALSE)))

nested.mean <- lapply(nested.mntd.dis, function(x) unlist(lapply(x, function(x) mean(x)) ) )

as.data.frame(nested.mean) -> diss
names(diss) <- names(nested.mean)
rownames(diss) -> diss$fd
factor(diss$fd, levels = c('L', 'H')) -> diss$fd
melt(diss) -> diss.long

  #Confidence Intervals for dissimilarity
  #temp from above bootstrap 
  ##
myReps = 1000
myCI = 0.95
myError <- myObs <- vector('list', length = length(temp))
names(myError) <- names(myObs) <- names(temp)
for(i in names(myError)) myError[[i]] <- list(L = numeric(0), H = numeric(0))

for(i in names(temp)) {
  for(j in c('L', 'H')) {
    print(paste('doing year', i, 'with fd', j))
    bittyTemp <- temp[[i]][[j]][which(names(temp[[i]][[j]]) == 'Agalinis_tenuifolia'):which(names(temp[[i]][[j]]) == 'Zizia_aurea')]
    bittyTemp.d <- sapply(seq(myReps), function(x) {
      comdistnt(bittyTemp[sample(1:24, 12, replace = F), ], phy.dist) %>% mean
    })
    myObs[[i]][[j]] <- mean(comdistnt(bittyTemp, phy.dist))
    myError[[i]][[j]] <- c(obs = myObs[[i]][[j]], quantile(bittyTemp.d, c((1-myCI) / 2, myCI+((1-myCI) / 2)))) #half jacknifing in
  }
}

myError -> err
as.data.frame(melt(err)) -> ERR #need to keep the obs, 2.5, 
valueIDs <- c("obs","lowCI", "highCI" )
cbind(ERR, valueIDs) -> NEWdf.test
names(NEWdf.test) <- c("value", "fd", "year", "valueID")
data2 <- dcast(NEWdf.test, year + fd ~ valueID, value.var = c("value"))
factor(data2$fd, levels = c('L', 'H')) -> data2$fd

PP <- data2 %>%
  mutate(group2 = case_when(year %in% c("2018.05", "2019.06") ~ "A",
                            TRUE ~ "B")) %>% 
  ggplot(aes(x = year, 
             y = obs, 
             color = fd, 
             group = interaction(fd, group2), 
             linetype = group2)) +
  geom_line(size = 1) +
  scale_linetype_manual("season",
                        values = c("A" = "dashed",
                                   "B" = "solid"),
                        labels = c("A" = "spring",
                                   "B" = "fall")) + 
  scale_color_manual("FD",
                     values = c("L" = "#332288",
                                "H" =  "#44AA99"
                     )) +
  geom_point(size = 3, alpha = 0.5)  +
  geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=.2,
                position=position_dodge(.1)) 

PP <- PP + labs(x = "year",  y= "dissimilarity (MNTD)", color="FD") + theme_bw() +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))






#######Dispersion by FD###############
group.datFD <- subset(sub.df.PA, select=c(year, fd, plot))
dis.ord <- adonis2(mntd.dist.sub ~ year + fd, data = group.datFD, strata="plot") #, permutations = perm) 

groupsFD <- sub.df.PA %>% 
  unite(yearFD, c("year", fd))

as.factor(groupsFD$yearFD) -> groupsFD$yearFD
mod2 <- betadisper(mntd.dist.sub, groupsFD$yearFD)
mod2P <- permutest(mod2, pairwise = TRUE, permutations = 999) ## Permutation test for F
mod2P

####Plot 
mod2.fd <- data.frame(group=mod2$group, distances=mod2$distances)
mod2.fd$group <- gsub( "2016_L", "Low planted", mod2.fd$group)
mod2.fd$group <- gsub( "2019.06_L", "Low year 3 spring", mod2.fd$group)
mod2.fd$group <- gsub( "2019.09_L", "Low year 3 summer", mod2.fd$group)
mod2.fd$group <- gsub( "2019.09_M", "Medium year 3 summer", mod2.fd$group)
mod2.fd$group <- gsub( "2016_H", "High planted", mod2.fd$group)
mod2.fd$group <- gsub( "2019.06_H", "High year 3 spring", mod2.fd$group)
mod2.fd$group <- gsub( "2019.09_H", "High year 3 summer", mod2.fd$group)

mod2.fd$group <- factor(mod2.fd$group , levels=c("Low planted", "Low year 3 spring","Low year 3 summer","High planted", "High year 3 spring","High year 3 summer"))

ggplot(mod2.fd, aes(x=group, y=distances, fill=group)) + geom_boxplot() +
  xlab("") +
  ylab("Distance to centroid") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  labs(fill="Survey", x= "Survey") + theme_bw() +
  scale_fill_manual("group", values = c("Low planted" = "#332288",
                                        "High planted" = "#44AA99", 
                                        "Low year 3 spring" = "#332288",
                                        "Low year 3 summer" = "#332288", 
                                        "High year 3 summer" = "#44AA99",
                                        "High year 3 spring" = "#44AA99" 
                                    
  )) 


