
remove(list = ls())
start_time = Sys.time()
setwd("~/GitHub/utlp-spm-mux")
library(igraph)
library(LinkPrediction)
library(lattice)
source("form_net.R")
source('structural_similarity.R')
# physicians_l2
# physicians_l1
# 'celegans_l1'
# 'celegans_l3'
network_layer1_name <- 'physicians_l1'
network_layer2_name <- 'physicians_l2'

net_l1_fpath <- paste0('mux-data/', network_layer1_name, '.txt')
net_l2_fpath <- paste0('mux-data/', network_layer2_name, '.txt')

gT <- graph_from_adjacency_matrix(form_net(net_l1_fpath), mode = "undirected", diag = FALSE)
global_gA <<- graph_from_adjacency_matrix(form_net(net_l2_fpath), mode = "undirected", diag = FALSE)

NCGT = nrow(as_adjacency_matrix(gT, type = "both"))
print(NCGT)
DENSITY_GT = sum(as_adjacency_matrix(gT, type = "both"))/(NCGT*NCGT-NCGT)
print(DENSITY_GT)
NCGA = nrow(as_adjacency_matrix(global_gA, type = "both"))
print(NCGA)
DENSITY_GA = sum(as_adjacency_matrix(global_gA, type = "both"))/(NCGA*NCGA-NCGA)
print(DENSITY_GA)

gT <<-erdos.renyi.game(nrow(as_adjacency_matrix(gT, type = "both")), DENSITY_GT, type = c("gnp"), directed = FALSE, loops = FALSE)
global_gA <<-erdos.renyi.game(nrow(as_adjacency_matrix(gT, type = "both")), DENSITY_GA, type = c("gnp"), directed = FALSE, loops = FALSE)


SS <- structural_similarity(gT)

KK = nrow(SS)
K_ToShow = KK*0.1
SS_M = SS[1:KK,1:KK]
SS_MUM = SS_M
#str(SS)
#str(SS_M)
index_Max = (which(SS_M == max(SS_M), arr.ind = TRUE))
index_Max = c(index_Max[1],index_Max[2])
index_Max_C = index_Max
print(index_Max_C)
SS_Diag = SS_M[index_Max]
# str(index_Max)
for (i in 2:KK) {
  SS_M[index_Max_C[1],]=-1
  SS_M[,index_Max_C[2]]=-1
  index_Max_C= which(SS_M == max(SS_M), arr.ind = TRUE)
  index_Max_C = c(index_Max_C[1],index_Max_C[2])
  print(index_Max_C)
  #Sys.sleep(0.5)
  SS_Diag = cbind(SS_Diag, SS_M[index_Max_C])
  index_Max = rbind(index_Max,index_Max_C)
}

PMR = matrix(0, nrow=KK, ncol=KK)
for(iPM in 1:KK){
  ToSetOne = index_Max[iPM,1]
  PMR[ToSetOne,iPM]=1
}
SS_PR = t(PMR) %*%  SS_MUM 
levelplot(t(SS_PR))
# Sys.sleep(1)
# levelplot(t(SS_PR))
# Sys.sleep(1)

PMC = matrix(0, nrow=KK, ncol=KK)
for(iPM in 1:KK){
  ToSetOne = index_Max[iPM,2]
  PMC[iPM,ToSetOne]=1
}

SS_PC = SS_MUM %*% t(PMC) 
# levelplot(t(SS_MUM))
# Sys.sleep(1)
# levelplot(t(SS_PC))
# Sys.sleep(1)

SS_P = t(PMR) %*%  SS_MUM %*% t(PMC)
# levelplot(t(SS_MUM))
# Sys.sleep(1)

StartOfShow = 3;
SS_P_ToShow = SS_P[StartOfShow:(K_ToShow),StartOfShow:(K_ToShow)]


# levelplot(t(SS_P_ToShow))
# p.strip <- list(cex=1.5, lines=2, fontface='bold') 
# ckey <- list(labels=list(cex=1.5), height=0.5) 
# x.scale <- list(cex=1, alternating=1) 
# y.scale <- list(cex=1, alternating=1) 

# levelplot(s, colorkey=ckey, par.strip.text=p.strip, 
#          scales=list(x=x.scale, y=y.scale)) 
# levelplot(s, colorkey=ckey, par.strip.text=p.strip, 
#           ) 
levelplot(t(SS_P_ToShow), 
          at = seq(0, 1, by=0.01),
          xlab = # 'Eigenvectors of alpha',
            'Eigenvectors of g(alpha)',
            #'Eigenvectors of  Randomized Physicians-Advice Network',
            # 'Eigen Vector Physicians-Advice Network',
          ylab = # 'Eigenvectors of beta'
            'Eigenvectors of g(beta)'
            #'Ei of  Randomized Physicians-Advice Network'
          )
            # 'Eigen Vector Physicians-Discuss Network')
# Sys.sleep(1)

print(sum(diag(SS_P_ToShow))/K_ToShow)
end_time = Sys.time()
end_time - start_time

## index_Max = cbind(index_Max,index_Max_C)
  


