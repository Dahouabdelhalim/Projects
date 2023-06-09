####
####
# Differences in social network structure and implications 
# for reproductive isolation in a pair of stickleback ecotypes - code for figures and analysis 
####
####

#Load in packages
library(ggplot2)
library(emmeans)
library(survcomp)
library(dplyr)
library(lme4)
library(car)
library(stats)
library(igraph)
library(lmerTest)
library(DirectedClustering)
library(ggpubr) 
library(wPerm)
library(ggpattern)

#Load in data files
# Nodes and edges for adjacency matrix - mixed groups
nodes_mixed <- read.csv("~/Desktop/Desktop/fall_shoaling/final data + code/nodes_mixed.csv")
edges_mixed <- read.csv("~/Desktop/Desktop/fall_shoaling/final data + code/edges_mixed.csv")
# Nodes and edges for adjacency matrix - white groups
nodes_white <- read.csv("~/Desktop/Desktop/fall_shoaling/final data + code/nodes_white.csv")
edges_white <- read.csv("~/Desktop/Desktop/fall_shoaling/final data + code/edges_white.csv")
# Nodes and edges for adjacency matrix - common groups
nodes_common <- read.csv("~/Desktop/Desktop/fall_shoaling/final data + code/nodes_common.csv")
edges_common <- read.csv("~/Desktop/Desktop/fall_shoaling/final data + code/edges_common.csv")

# load in all edges as one file for some later analysis
all_edges <- read.csv("~/Desktop/Desktop/fall_shoaling/final data + code/all_edges.csv")

# Individual level data / information
ind_shoaling_data <- read.csv("~/Desktop/Desktop/fall_shoaling/final data + code/ind_shoaling_data.csv")
ind_shoaling_data$fish_type <- as.factor(ind_shoaling_data$fish_type)
ind_shoaling_data$group_type <- as.factor(ind_shoaling_data$group_type)
ind_shoaling_data$group <- as.factor(ind_shoaling_data$group)
ind_shoaling_data$family <- as.factor(ind_shoaling_data$family)
ind_shoaling_data$holding_tank <- as.factor(ind_shoaling_data$holding_tank)

# Separate out nodes and edges by group for constructing networks for each group
split_nodes <- split(nodes_mixed,rep(1:8,each=6))
split_edges <- split(edges_mixed,rep(1:8,each=15))
split_nodes_c <- split(nodes_common,rep(1:8,each=6))
split_edges_c <- split(edges_common,rep(1:8,each=15))
split_nodes_w <- split(nodes_white,rep(1:8,each=6))
split_edges_w <- split(edges_white,rep(1:8,each=15))

###########
# Construct networks from data frame - igraph 
g_MA <- graph_from_data_frame(d=split_edges$'1', vertices=split_nodes$'1', directed = FALSE)
g_MB <- graph_from_data_frame(d=split_edges$'2', vertices=split_nodes$'2', directed = FALSE)
g_MC <- graph_from_data_frame(d=split_edges$'3', vertices=split_nodes$'3', directed = FALSE)
g_MD <- graph_from_data_frame(d=split_edges$'4', vertices=split_nodes$'4', directed = FALSE)
g_ME <- graph_from_data_frame(d=split_edges$'5', vertices=split_nodes$'5', directed = FALSE)
g_MF <- graph_from_data_frame(d=split_edges$'6', vertices=split_nodes$'6', directed = FALSE)
g_MG <- graph_from_data_frame(d=split_edges$'7', vertices=split_nodes$'7', directed = FALSE)
g_MH <- graph_from_data_frame(d=split_edges$'8', vertices=split_nodes$'8', directed = FALSE)

g_CA <- graph_from_data_frame(d=split_edges_c$'1', vertices=split_nodes_c$'1', directed = FALSE)
g_CB <- graph_from_data_frame(d=split_edges_c$'2', vertices=split_nodes_c$'2', directed = FALSE)
g_CC <- graph_from_data_frame(d=split_edges_c$'3', vertices=split_nodes_c$'3', directed = FALSE)
g_CD <- graph_from_data_frame(d=split_edges_c$'4', vertices=split_nodes_c$'4', directed = FALSE)
g_CE <- graph_from_data_frame(d=split_edges_c$'5', vertices=split_nodes_c$'5', directed = FALSE)
g_CF <- graph_from_data_frame(d=split_edges_c$'6', vertices=split_nodes_c$'6', directed = FALSE)
g_CG <- graph_from_data_frame(d=split_edges_c$'7', vertices=split_nodes_c$'7', directed = FALSE)
g_CH <- graph_from_data_frame(d=split_edges_c$'8', vertices=split_nodes_c$'8', directed = FALSE)

g_WA <- graph_from_data_frame(d=split_edges_w$'1', vertices=split_nodes_w$'1', directed = FALSE)
g_WB <- graph_from_data_frame(d=split_edges_w$'2', vertices=split_nodes_w$'2', directed = FALSE)
g_WC <- graph_from_data_frame(d=split_edges_w$'3', vertices=split_nodes_w$'3', directed = FALSE)
g_WD <- graph_from_data_frame(d=split_edges_w$'4', vertices=split_nodes_w$'4', directed = FALSE)
g_WE <- graph_from_data_frame(d=split_edges_w$'5', vertices=split_nodes_w$'5', directed = FALSE)
g_WF <- graph_from_data_frame(d=split_edges_w$'6', vertices=split_nodes_w$'6', directed = FALSE)
g_WG <- graph_from_data_frame(d=split_edges_w$'7', vertices=split_nodes_w$'7', directed = FALSE)
g_WH <- graph_from_data_frame(d=split_edges_w$'8', vertices=split_nodes_w$'8', directed = FALSE)

# Convert networks to adjaceny matrix  to compute clustering 
adj_CA = as_adjacency_matrix(g_CA, sparse=F, attr="weight")
adj_CB = as_adjacency_matrix(g_CB, sparse=F, attr="weight")
adj_CC = as_adjacency_matrix(g_CC, sparse=F, attr="weight")
adj_CD = as_adjacency_matrix(g_CD, sparse=F, attr="weight")
adj_CE = as_adjacency_matrix(g_CE, sparse=F, attr="weight")
adj_CF = as_adjacency_matrix(g_CF, sparse=F, attr="weight")
adj_CG = as_adjacency_matrix(g_CG, sparse=F, attr="weight")
adj_CH = as_adjacency_matrix(g_CH, sparse=F, attr="weight")
#
adj_WA = as_adjacency_matrix(g_WA, sparse=F, attr="weight")
adj_WB = as_adjacency_matrix(g_WB, sparse=F, attr="weight")
adj_WC = as_adjacency_matrix(g_WC, sparse=F, attr="weight")
adj_WD = as_adjacency_matrix(g_WD, sparse=F, attr="weight")
adj_WE = as_adjacency_matrix(g_WE, sparse=F, attr="weight")
adj_WF = as_adjacency_matrix(g_WF, sparse=F, attr="weight")
adj_WG = as_adjacency_matrix(g_WG, sparse=F, attr="weight")
adj_WH = as_adjacency_matrix(g_WH, sparse=F, attr="weight")
#
adj_MA = as_adjacency_matrix(g_MA, sparse=F, attr="weight")
adj_MB = as_adjacency_matrix(g_MB, sparse=F, attr="weight")
adj_MC = as_adjacency_matrix(g_MC, sparse=F, attr="weight")
adj_MD = as_adjacency_matrix(g_MD, sparse=F, attr="weight")
adj_ME = as_adjacency_matrix(g_ME, sparse=F, attr="weight")
adj_MF = as_adjacency_matrix(g_MF, sparse=F, attr="weight")
adj_MG = as_adjacency_matrix(g_MG, sparse=F, attr="weight")
adj_MH = as_adjacency_matrix(g_MH, sparse=F, attr="weight")

# Compute clustering coefficient for the networks
ca <- as.data.frame(ClustF(adj_CA, type = "undirected"))
cb <- as.data.frame(ClustF(adj_CB, type = "undirected"))
cc <- as.data.frame(ClustF(adj_CC, type = "undirected"))
cd <- as.data.frame(ClustF(adj_CD, type = "undirected"))
ce <- as.data.frame(ClustF(adj_CE, type = "undirected"))
cf <- as.data.frame(ClustF(adj_CF, type = "undirected"))
cg <- as.data.frame(ClustF(adj_CG, type = "undirected"))
ch <- as.data.frame(ClustF(adj_CH, type = "undirected"))
clustering_common <- rbind(ca,cb,cc,cd,ce,cf,cg,ch)
#
ma <- as.data.frame(ClustF(adj_MA, type = "undirected"))
mb <- as.data.frame(ClustF(adj_MB, type = "undirected"))
mc <- as.data.frame(ClustF(adj_MC, type = "undirected"))
md <- as.data.frame(ClustF(adj_MD, type = "undirected"))
me <- as.data.frame(ClustF(adj_ME, type = "undirected"))
mf <- as.data.frame(ClustF(adj_MF, type = "undirected"))
mg <- as.data.frame(ClustF(adj_MG, type = "undirected"))
mh <- as.data.frame(ClustF(adj_MH, type = "undirected"))
clustering_mixed <- rbind(ma,mb,mc,md,me,mf,mg,mh)
#
wa <- as.data.frame(ClustF(adj_WA, type = "undirected"))
wb <- as.data.frame(ClustF(adj_WB, type = "undirected"))
wc <- as.data.frame(ClustF(adj_WC, type = "undirected"))
wd <- as.data.frame(ClustF(adj_WD, type = "undirected"))
we <- as.data.frame(ClustF(adj_WE, type = "undirected"))
wf <- as.data.frame(ClustF(adj_WF, type = "undirected"))
wg <- as.data.frame(ClustF(adj_WG, type = "undirected"))
wh <- as.data.frame(ClustF(adj_WH, type = "undirected"))
clustering_white <- rbind(wa,wb,wc,wd,we,wf,wg,wh)
#
cl <- rbind(clustering_common,clustering_mixed,clustering_white)
ind_shoaling_data <- cbind(cl, ind_shoaling_data)

###
# DATA ANALYSIS 

# Boldness - linear model
boldness_lm <- lm(log(latency_sec) ~ length_mm + fish_type,
                       data = ind_shoaling_data)
Anova(boldness_lm, type = "III")

# Body size - Kruskal Wallis test
kruskal.test(ind_shoaling_data$length_mm ~ ind_shoaling_data$fish_type)

# Activity - linear mixed model
activity_lm <- lmer(log(distance_swam_m) ~ length_mm + fish_type + group_type +
                      (1|group), data = ind_shoaling_data)
anova(activity_lm)

# Clustering coefficient - group level linear model
# get mean values for clustering by group type
means_cl <- as.data.frame(aggregate(LocalCC ~ group + group_type, data=ind_shoaling_data, mean))
# run model
means_cl_lm <- lm(LocalCC ~ group_type, data = means_cl)
Anova(means_cl_lm, type = "III")
emmeans(means_cl_lm, pairwise ~ group_type)

# Clustering coefficient - white in mixed vs white in homogeneous and 
#   common in mixed vs common in homogeneous 

# separate data by ecotype
shoaling_data_white <- ind_shoaling_data[ind_shoaling_data$fish_type == "W",]
shoaling_data_common <- ind_shoaling_data[ind_shoaling_data$fish_type == "C",]

# run tests 
# white
wilcox.test(LocalCC ~ group_type, data = shoaling_data_white)
# common
wilcox.test(LocalCC ~ group_type, data = shoaling_data_common)

# Exact permutation test for comparing clustering of
#     white vs. common individuals within mixed groups 

clustering_mixed <- ind_shoaling_data[ind_shoaling_data$group_type == "M",]

# reorder so values for each group are together when we split data into separate 
#     dataframesfor each group 
clustering_mixed <- clustering_mixed[order(clustering_mixed$group),]
#
clustering_mixed_split <- split(clustering_mixed,rep(1:8,each=6))

#vector to store Kruskal-Wallis test statistics   
clust_K_MA <- vector()
clust_K_MB <- vector()
clust_K_MC <- vector()
clust_K_MD <- vector()
clust_K_ME <- vector()
clust_K_MF <- vector()
clust_K_MG <- vector()
clust_K_MH <- vector()

# MA
dfs <- vector("list",3)

for(i in 1:1000){
  clustering_mixed_split$'1'$new <- sample(clustering_mixed_split$'1'$LocalCC, replace = FALSE)
  dfs[[i]] <- clustering_mixed_split$'1'
  kruskal_test <- kruskal.test(new ~ fish_type, data = dfs[[i]])
  clust_K_MA <- append(clust_K_MA, kruskal_test$statistic)
}

# MB
dfs <- vector("list",3)

for(i in 1:1000){
  clustering_mixed_split$'2'$new <- sample(clustering_mixed_split$'2'$LocalCC, replace = FALSE)
  dfs[[i]] <- clustering_mixed_split$'2'
  kruskal_test <- kruskal.test(new ~ fish_type, data = dfs[[i]])
  clust_K_MB <- append(clust_K_MB, kruskal_test$statistic)
}

# MC
dfs <- vector("list",3)

for(i in 1:1000){
  clustering_mixed_split$'3'$new <- sample(clustering_mixed_split$'3'$LocalCC, replace = FALSE)
  dfs[[i]] <- clustering_mixed_split$'3'
  kruskal_test <- kruskal.test(new ~ fish_type, data = dfs[[i]])
  clust_K_MC <- append(clust_K_MC, kruskal_test$statistic)
}

# MD
dfs <- vector("list",3)

for(i in 1:1000){
  clustering_mixed_split$'4'$new <- sample(clustering_mixed_split$'4'$LocalCC, replace = FALSE)
  dfs[[i]] <- clustering_mixed_split$'4'
  kruskal_test <- kruskal.test(new ~ fish_type, data = dfs[[i]])
  clust_K_MD <- append(clust_K_MD, kruskal_test$statistic)
}

# ME
dfs <- vector("list",3)

for(i in 1:1000){
  clustering_mixed_split$'5'$new <- sample(clustering_mixed_split$'5'$LocalCC, replace = FALSE)
  dfs[[i]] <- clustering_mixed_split$'5'
  kruskal_test <- kruskal.test(new ~ fish_type, data = dfs[[i]])
  clust_K_ME <- append(clust_K_ME, kruskal_test$statistic)
}

# MF
dfs <- vector("list",3)

for(i in 1:1000){
  clustering_mixed_split$'6'$new <- sample(clustering_mixed_split$'6'$LocalCC, replace = FALSE)
  dfs[[i]] <- clustering_mixed_split$'6'
  kruskal_test <- kruskal.test(new ~ fish_type, data = dfs[[i]])
  clust_K_MF <- append(clust_K_MF, kruskal_test$statistic)
}

# MG
dfs <- vector("list",3)

for(i in 1:1000){
  clustering_mixed_split$'7'$new <- sample(clustering_mixed_split$'7'$LocalCC, replace = FALSE)
  dfs[[i]] <- clustering_mixed_split$'7'
  kruskal_test <- kruskal.test(new ~ fish_type, data = dfs[[i]])
  clust_K_MG <- append(clust_K_MG, kruskal_test$statistic)
}

# MH
dfs <- vector("list",3)

for(i in 1:1000){
  clustering_mixed_split$'8'$new <- sample(clustering_mixed_split$'8'$LocalCC, replace = FALSE)
  dfs[[i]] <- clustering_mixed_split$'8'
  kruskal_test <- kruskal.test(new ~ fish_type, data = dfs[[i]])
  clust_K_MH <- append(clust_K_MH, kruskal_test$statistic)
}

###

# actual k scores
kMA <- kruskal.test(LocalCC ~ fish_type, data = clustering_mixed_split$'1')
kMB <- kruskal.test(LocalCC ~ fish_type, data = clustering_mixed_split$'2')
kMC <- kruskal.test(LocalCC ~ fish_type, data = clustering_mixed_split$'3')
kMD <- kruskal.test(LocalCC ~ fish_type, data = clustering_mixed_split$'4')
kME <- kruskal.test(LocalCC ~ fish_type, data = clustering_mixed_split$'5')
kMF <- kruskal.test(LocalCC ~ fish_type, data = clustering_mixed_split$'6')
kMG <- kruskal.test(LocalCC ~ fish_type, data = clustering_mixed_split$'7')
kMH <- kruskal.test(LocalCC ~ fish_type, data = clustering_mixed_split$'8')

# log transform of the Kruskal wallis test statistic distribution is approximately normal, 
#    so we can use t table to see if observed value is significantly different from simulated values
mean(log(clust_K_MA + 0.0000000000000001))
sd(log(clust_K_MA + 0.0000000000000001))
log(kMA$statistic)
#
pval_MA_cl <- 1 - pnorm(0.85, mean = -0.87, sd = 1.57)

#
mean(log(clust_K_MB + 0.0000000000000001))
sd(log(clust_K_MB + 0.0000000000000001))
log(kMB$statistic)
#
pval_MB_cl <- 1 - pnorm(0.17, mean = -0.97, sd = 1.57)

#
mean(log(clust_K_MC + 0.0000000000000001))
sd(log(clust_K_MC + 0.0000000000000001))
log(kMC$statistic)
#
pval_MC_cl <- 1 - pnorm(-0.85, mean = -0.93, sd = 1.57)

#
mean(log(clust_K_MD + 0.0000000000000001))
sd(log(clust_K_MD + 0.0000000000000001))
log(kMD$statistic)
#
pval_MD_cl <- 1 - pnorm(0.17, mean = -0.99, sd = 1.58)

#
mean(log(clust_K_ME + 0.0000000000000001))
sd(log(clust_K_ME + 0.0000000000000001))
log(kME$statistic)
#
pval_ME_cl <- 1 - pnorm(1.35, mean = -0.91, sd = 1.56)

#
mean(log(clust_K_MF + 0.0000000000000001))
sd(log(clust_K_MF + 0.0000000000000001))
log(kMF$statistic)
#
pval_MF_cl <- 1 - pnorm(-0.84, mean = -0.85, sd = 1.52)

#
mean(log(clust_K_MG + 0.0000000000000001))
sd(log(clust_K_MG + 0.0000000000000001))
log(kMG$statistic)
#
pval_MG_cl <- 1 - pnorm(-0.84, mean = -0.86, sd = 1.57)

#
mean(log(clust_K_MH + 0.0000000000000001))
sd(log(clust_K_MH + 0.0000000000000001))
log(kMH$statistic)
#
pval_MH_cl <- 1 - pnorm(-0.84, mean = -0.96, sd = 1.54)

# Take mean of p values and get p value
pvals_cl <- as.data.frame(rbind(pval_MA_cl, pval_MB_cl, pval_MC_cl, pval_MD_cl,
                                pval_ME_cl, pval_MF_cl, pval_MG_cl, pval_MH_cl))
mean(pvals_cl$V1)
#

# Interaction rate - group level linear model 
# get mean values for interaction rate by group type
means_interaction_rate <- as.data.frame(aggregate(weight ~ group + group_type, data=all_edges, mean))
# run model
means_interaction_rate_lm <- lm(weight ~ group_type, data = means_interaction_rate)
Anova(means_interaction_rate_lm, type = "III")

# Exact permutation test for comparing interaction rate of 
#       white vs. common individuals within mixed groups 

# start with edges which represent all interaction rates 
View(split_edges)

# set up vectors to store test Kruskal-Wallis test statistics   
kscores_MA <- vector()
kscores_MB <- vector()
kscores_MC <- vector()
kscores_MD <- vector()
kscores_ME <- vector()
kscores_MF <- vector()
kscores_MG <- vector()
kscores_MH <- vector()

# MA - repeat this for loop for all groups 
dfs <- vector("list",3)

for(i in 1:1000){
  split_edges$'1'$new <- sample(split_edges$'1'$weight, replace = FALSE)
  dfs[[i]] <- split_edges$'1'
  kruskal_test <- kruskal.test(new ~ interaction, data = dfs[[i]])
  kscores_MA <- append(kscores_MA, kruskal_test$statistic)
}

# MB
dfs <- vector("list",3)

for(i in 1:1000){
  split_edges$'2'$new <- sample(split_edges$'2'$weight, replace = FALSE)
  dfs[[i]] <- split_edges$'2'
  kruskal_test <- kruskal.test(new ~ interaction, data = dfs[[i]])
  kscores_MB <- append(kscores_MB, kruskal_test$statistic)
}

# MC
dfs <- vector("list",3)

for(i in 1:1000){
  split_edges$'3'$new <- sample(split_edges$'3'$weight, replace = FALSE)
  dfs[[i]] <- split_edges$'3'
  kruskal_test <- kruskal.test(new ~ interaction, data = dfs[[i]])
  kscores_MC <- append(kscores_MC, kruskal_test$statistic)
}

# MD
dfs <- vector("list",3)

for(i in 1:1000){
  split_edges$'4'$new <- sample(split_edges$'4'$weight, replace = FALSE)
  dfs[[i]] <- split_edges$'4'
  kruskal_test <- kruskal.test(new ~ interaction, data = dfs[[i]])
  kscores_MD <- append(kscores_MD, kruskal_test$statistic)
}

# ME
dfs <- vector("list",3)

for(i in 1:1000){
  split_edges$'5'$new <- sample(split_edges$'5'$weight, replace = FALSE)
  dfs[[i]] <- split_edges$'5'
  kruskal_test <- kruskal.test(new ~ interaction, data = dfs[[i]])
  kscores_ME <- append(kscores_ME, kruskal_test$statistic)
}

# MF
dfs <- vector("list",3)

for(i in 1:1000){
  split_edges$'6'$new <- sample(split_edges$'6'$weight, replace = FALSE)
  dfs[[i]] <- split_edges$'6'
  kruskal_test <- kruskal.test(new ~ interaction, data = dfs[[i]])
  kscores_MF <- append(kscores_MF, kruskal_test$statistic)
}

# MG
dfs <- vector("list",3)

for(i in 1:1000){
  split_edges$'7'$new <- sample(split_edges$'7'$weight, replace = FALSE)
  dfs[[i]] <- split_edges$'7'
  kruskal_test <- kruskal.test(new ~ interaction, data = dfs[[i]])
  kscores_MG <- append(kscores_MG, kruskal_test$statistic)
}

# MH
dfs <- vector("list",3)

for(i in 1:1000){
  split_edges$'8'$new <- sample(split_edges$'8'$weight, replace = FALSE)
  dfs[[i]] <- split_edges$'8'
  kruskal_test <- kruskal.test(new ~ interaction, data = dfs[[i]])
  kscores_MH <- append(kscores_MH, kruskal_test$statistic)
}

###

# actual k scores
kMA <- kruskal.test(weight ~ interaction, data = split_edges$'1')
kMB <- kruskal.test(weight ~ interaction, data = split_edges$'2')
kMC <- kruskal.test(weight ~ interaction, data = split_edges$'3')
kMD <- kruskal.test(weight ~ interaction, data = split_edges$'4')
kME <- kruskal.test(weight ~ interaction, data = split_edges$'5')
kMF <- kruskal.test(weight ~ interaction, data = split_edges$'6')
kMG <- kruskal.test(weight ~ interaction, data = split_edges$'7')
kMH <- kruskal.test(weight ~ interaction, data = split_edges$'8')

# log transform of the Kruskal wallis test statistic distribution is approximately normal, 
#    so we can use t table to see if observed value is significantly different from simulated values
mean(log(kscores_MA + 0.0000000000000001))
sd(log(kscores_MA + 0.0000000000000001))
log(kMA$statistic)
#
pval_MA <- 1 - pnorm(0.588, mean = 0.117, sd = 2.03)

#
mean(log(kscores_MB + 0.0000000000000001))
sd(log(kscores_MB + 0.0000000000000001))
log(kMB$statistic)
#
pval_MB <- 1 - pnorm(0.857, mean = 0.104, sd = 2.60)

#
mean(log(kscores_MC + 0.0000000000000001))
sd(log(kscores_MC + 0.0000000000000001))
log(kMC$statistic)
#
pval_MC <- 1 - pnorm(0.857, mean = 0.252, sd = 1.15)

#
mean(log(kscores_MD + 0.0000000000000001))
sd(log(kscores_MD + 0.0000000000000001))
log(kMD$statistic)
#
pval_MD <- 1 - pnorm(0.951, mean = 0.06, sd = 2.86)

#
mean(log(kscores_ME + 0.0000000000000001))
sd(log(kscores_ME + 0.0000000000000001))
log(kME$statistic)
#
pval_ME <- 1 - pnorm(1.11, mean = 0.24, sd = 1.16)

#
mean(log(kscores_MF + 0.0000000000000001))
sd(log(kscores_MF + 0.0000000000000001))
log(kMF$statistic)
#
pval_MF <- 1 - pnorm(1.63, mean = 0.18, sd = 2.00)

#
mean(log(kscores_MG + 0.0000000000000001))
sd(log(kscores_MG + 0.0000000000000001))
log(kMG$statistic)
#
pval_MG <- 1 - pnorm(-0.86, mean = 0.14, sd = 2.01)

#
mean(log(kscores_MH + 0.0000000000000001))
sd(log(kscores_MH + 0.0000000000000001))
log(kMH$statistic)
#
pval_MH <- 1 - pnorm(1.77, mean = 0.168, sd = 1.19)

# Take mean of p values and get p value
pvals_interaction_rate <- as.data.frame(rbind(pval_MA, pval_MB, pval_MC, pval_MD,
                                              pval_ME, pval_MF, pval_MG, pval_MH))
mean(pvals_interaction_rate$V1)
#

########

# Correlations between behaviors 

#Compute correlations between all behaviors for each ecotype 
perm.relation(shoaling_data_white$latency_sec, shoaling_data_white$distance_swam_m,
              method = "spearman", R = 10000)
perm.relation(shoaling_data_white$latency_sec, shoaling_data_white$LocalCC,
              method = "spearman", R = 10000)
perm.relation(shoaling_data_white$latency_sec, shoaling_data_white$total_interaction_time,
              method = "spearman", R = 10000)
perm.relation(shoaling_data_white$total_interaction_time, shoaling_data_white$distance_swam_m,
              method = "spearman", R = 10000)
#
perm.relation(shoaling_data_common$latency_sec, shoaling_data_common$distance_swam_m,
              method = "spearman", R = 10000)
perm.relation(shoaling_data_common$latency_sec, shoaling_data_common$LocalCC,
              method = "spearman", R = 10000)
perm.relation(shoaling_data_common$latency_sec, shoaling_data_common$total_interaction_time,
              method = "spearman", R = 10000)
perm.relation(shoaling_common$total_interaction_time, shoaling_data_common$distance_swam_m,
              method = "spearman", R = 10000)
