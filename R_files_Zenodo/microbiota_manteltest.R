##Microbiota analysis 2019

require("ape")
require("vegan")
require("reshape2")
require("ade4")

#Import distance matrix of fish guts & water from Apoyo, Xiloa, Nicaragua & Managua

#Either weighted/unweighted unifrac or bray-curtis data
distmat_complete <- read.csv("distance-matrix_wuf.tsv", dec = ",", sep = "\\t", as.is = T)
#distmat_complete <- read.csv("distance-matrix_uwuf.tsv", dec = ",", sep = "\\t", as.is = T)
#distmat_complete <- read.csv("distance-matrix_bc.tsv", dec = ",", sep = "\\t", as.is = T)

rownames(distmat_complete) <- distmat_complete[,1]
distmat_complete <- distmat_complete[,-1]
names_vec <- rownames(distmat_complete)
colnames(distmat_complete) <- names_vec

#Import metadata file
ids_complete <- read.csv("samples.txt", dec = ",", sep = "\\t", as.is = T)

#Subset: only crater lakes

distmat_crater <- distmat_complete[c(2:14,16,17,19:21,23:57,59,61:65,67:90,92:95,98,100:103,107,110,111,114:116,118:127,129,130),c(2:14,16,17,19:21,23:57,59,61:65,67:90,92:95,98,100:103,107,110,111,114:116,118:127,129,130)]
ids_crater <- ids_complete[c(2:14,16,17,19:21,23:57,59,61:65,67:90,92:95,98,100:103,107,110,111,114:116,118:127,129,130),]

#Xiloa distmat

distmat_xiloa <- distmat_crater[c(1,5,8,10:12,14,17,24,26,27,29:31,33,35:37,40:45,47,48,50:52,54:56,58,60:63,65,67:69,71,72,75:77,80:82,89,91,95,98,101,106,110),c(1,5,8,10:12,14,17,24,26,27,29:31,33,35:37,40:45,47,48,50:52,54:56,58,60:63,65,67:69,71,72,75:77,80:82,89,91,95,98,101,106,110)]
ids_xiloa <- ids_crater[c(1,5,8,10:12,14,17,24,26,27,29:31,33,35:37,40:45,47,48,50:52,54:56,58,60:63,65,67:69,71,72,75:77,80:82,89,91,95,98,101,106,110),]

#Xiloa distmat, only SIA samples

distmat_xiloa_sia <- distmat_xiloa[c(1:5,7:38,40,42,44:47,49:51,53:56),c(1:5,7:38,40,42,44:47,49:51,53:56)]
ids_xiloa_sia <- ids_xiloa[c(1:5,7:38,40,42,44:47,49:51,53:56),]

ids_xiloa_sia$species <- as.factor(ids_xiloa_sia$species)

#Apoyo distmat

distmat_apoyo <- distmat_crater[c(2:4,6,7,9,13,15,16,18:23,25,28,32,34,38,39,46,49,53,57,59,64,66,70,73,74,78,79,83:88,90,92:94,96,97,99,100,102:105,107:109),c(2:4,6,7,9,13,15,16,18:23,25,28,32,34,38,39,46,49,53,57,59,64,66,70,73,74,78,79,83:88,90,92:94,96,97,99,100,102:105,107:109)]
ids_apoyo <- ids_crater[c(2:4,6,7,9,13,15,16,18:23,25,28,32,34,38,39,46,49,53,57,59,64,66,70,73,74,78,79,83:88,90,92:94,96,97,99,100,102:105,107:109),]

#Apoyo distmat, only SIA samples

distmat_apoyo_sia <- distmat_apoyo[c(1:45, 47:54),c(1:45, 47:54)]

ids_apoyo_sia <- ids_apoyo[-46,]
ids_apoyo_sia$species <- as.factor(ids_apoyo_sia$species)

#Xiloa: Pairwise distances between individuals, only SIA samples

D <- as.matrix(distmat_xiloa_sia)
D[upper.tri(D)] <- NA
x <- subset(melt(D), value!=0)
x <- as.data.frame(x)

row.names(x) <- NULL

#Nitrogen
n_dist <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  dif <- abs(ids_xiloa_sia$n[which(ids_xiloa_sia$id==id1)]-ids_xiloa_sia$n[which(ids_xiloa_sia$id==id2)])
  n_dist[g] <- dif
}

x$n_dist <- n_dist 

#carbon
c_dist <- NA

for (g in 1:length(x$Var1)){
  id1 <- x$Var1[g]
  id2 <- x$Var2[g]
  dif <- abs(ids_xiloa_sia$c[which(ids_xiloa_sia$id==id1)]-ids_xiloa_sia$c[which(ids_xiloa_sia$id==id2)])
  c_dist[g] <- dif
}

x$c_dist <- c_dist 

#Apoyo: Pairwise distances for between individuals, only SIA samples

D2 <- as.matrix(distmat_apoyo_sia)
D2[upper.tri(D2)] <- NA
x2 <- subset(melt(D2), value!=0)
x2 <- as.data.frame(x2)

row.names(x2) <- NULL

#Nitrogen
n_dist <- NA

for (g in 1:length(x2$Var1)){
  id1 <- x2$Var1[g]
  id2 <- x2$Var2[g]
  dif <- abs(ids_apoyo_sia$n[which(ids_apoyo_sia$id==id1)]-ids_apoyo_sia$n[which(ids_apoyo_sia$id==id2)])
  n_dist[g] <- dif
}

x2$n_dist <- n_dist 

#carbon
c_dist <- NA

for (g in 1:length(x2$Var1)){
  id1 <- x2$Var1[g]
  id2 <- x2$Var2[g]
  dif <- abs(ids_apoyo_sia$c[which(ids_apoyo_sia$id==id1)]-ids_apoyo_sia$c[which(ids_apoyo_sia$id==id2)])
  c_dist[g] <- dif
}

x2$c_dist <- c_dist 

#Mantel test

x_wuf <- x[,c(1,2,3)]
distmat_mb <- acast(x_wuf, Var1 ~ Var2)

x_n <- x[,c(1,2,4)]
distmat_n <- acast(x_n, Var1 ~ Var2)

x_c <- x[,c(1,2,5)]
distmat_c <- acast(x_c, Var1 ~ Var2)

x2_wuf <- x2[,c(1,2,3)]
distmat_mb2 <- acast(x2_wuf, Var1 ~ Var2)

x2_n <- x2[,c(1,2,4)]
distmat_n2 <- acast(x2_n, Var1 ~ Var2)

x2_c <- x2[,c(1,2,5)]
distmat_c2 <- acast(x2_c, Var1 ~ Var2)

#Remove 2 samples that occur only in rows or columns 
distmat_mb <- distmat_mb[-18,-25]
distmat_n <- distmat_n[-18,-25]
distmat_c <- distmat_c[-18,-25]

distmat_mb2 <- distmat_mb2[-4,-36]
distmat_n2 <- distmat_n2[-4,-36]
distmat_c2 <- distmat_c2[-4,-36]

diag(distmat_n) <- diag(distmat_c) <- diag(distmat_mb) <- diag(distmat_n2) <- diag(distmat_c2) <- diag(distmat_mb2) <- 0

for (i in 1:length(colnames(distmat_mb))){
  for (g in 1:length(rownames(distmat_mb))){
    if (is.na(distmat_mb[i,g])){
      distmat_mb[i,g] <- distmat_mb[g,i]
    }
  }
}

for (i in 1:length(colnames(distmat_n))){
  for (g in 1:length(rownames(distmat_n))){
    if (is.na(distmat_n[i,g])){
      distmat_n[i,g] <- distmat_n[g,i]
    }
  }
}

for (i in 1:length(colnames(distmat_c))){
  for (g in 1:length(rownames(distmat_c))){
    if (is.na(distmat_c[i,g])){
      distmat_c[i,g] <- distmat_c[g,i]
    }
  }
}

for (i in 1:length(colnames(distmat_mb2))){
  for (g in 1:length(rownames(distmat_mb2))){
    if (is.na(distmat_mb2[i,g])){
      distmat_mb2[i,g] <- distmat_mb2[g,i]
    }
  }
}

for (i in 1:length(colnames(distmat_n2))){
  for (g in 1:length(rownames(distmat_n2))){
    if (is.na(distmat_n2[i,g])){
      distmat_n2[i,g] <- distmat_n2[g,i]
    }
  }
}

for (i in 1:length(colnames(distmat_c2))){
  for (g in 1:length(rownames(distmat_c2))){
    if (is.na(distmat_c2[i,g])){
      distmat_c2[i,g] <- distmat_c2[g,i]
    }
  }
}

#mantel.rtest

dist_mb <- dist(distmat_mb)
dist_c <- dist(distmat_c)
dist_n <- dist(distmat_n)

dist_mb2 <- dist(distmat_mb2)
dist_c2 <- dist(distmat_c2)
dist_n2 <- dist(distmat_n2)

#Mantel test for Xiloa
mantel.rtest(dist_mb, dist_n, nrepet = 999)
mantel.rtest(dist_mb, dist_c, nrepet = 999)

#Mantel test for Apoyo
mantel.rtest(dist_mb2, dist_n2, nrepet = 999)
mantel.rtest(dist_mb2, dist_c2, nrepet = 999)