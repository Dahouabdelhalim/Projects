# run a PGLS of head shape and buccal shape on feeding and reproduction

library(geomorph)
library(phytools)

# load geomorph dataframe with procrustes coordinates
int_gmdf <- readRDS("data/int_gmdf.RDS")
ext_gmdf <- readRDS("data/ext_gmdf.RDS")

# PGLS of buccal cavity shape on feeding and reproduction:
int_fit <- procD.lm(data = int_gmdf,
                    coords ~ reproduction*feeding,
                    iter = 999)
summary(int_fit) # both feeding and reproduction are significant

# vs overall head shape:
ext_fit <- procD.pgls(data = ext_gmdf,
                    coords ~ reproduction*feeding,
                    phy = phy, iter = 9999)
summary(ext_fit) # only feeding is significant




# average shapes
fit_inxn <- procD.lm(data = int_gmdf, 
                   coords ~ inxn, 
                   iter = 999)

for (i in 1:length(unique(int_gmdf$inxn))) {
  
  combo <- unique(int_gmdf$inxn)[i]
  predicted.coordinates <- fit_inxn$fitted[which(int_gmdf$inxn == combo), ]
  predicted.3d <- arrayspecs(predicted.coordinates, 
                             p = int.gm.shapes$p,
                             k = int.gm.shapes$k)
  
  new.shape <- data.frame(x = predicted.3d[ , 1, 1],
                          y = predicted.3d[ , 2, 1],
                          group = rep(combo, dim(predicted.3d)[1]))
  
  pt.order <- c(5, 15:28, 6, 9, 29:42, 10, 11, 43:46, 12, 5)
  
  shape.points <- new.shape[pt.order, ]
  shape.poly <- sf::st_polygon(list(as.matrix(shape.points[ , 1:2])))
  plot(shape.poly, main = combo)
  points(shape.points[-1, ])
  
  shape.obj <- list(group = combo,
                    coords = new.shape[, 1:2],
                    points = shape.points,
                    polygon = shape.poly)
  if (i == 1) {
    predicted.shapes <- list(shape.obj)
  } else {
    predicted.shapes[[i]] <- shape.obj
  }
  names(predicted.shapes)[i] <- as.character(combo)
}

# get procrustes distances between shapes
shape_dist <- function(shape_list) {
  
  combos <- t(combn(1:length(shape_list), 2))
  combo_names <- t(combn(names(shape_list), 2))
  df_out <- data.frame(combo_names, procrustes_dist = NA)
  colnames(df_out)[1:2] <- c("shape1", "shape2")
  
  for (i in 1:nrow(df_out)) {
    
    s1 <- predicted.shapes[[combos[i, 1]]]$coords
    s2 <- predicted.shapes[[combos[i, 2]]]$coords
    
    d <- c()
    for (j in 1:nrow(s1)) {
      d <- c(d, dist(rbind(s1[j, ], s2[j, ])))
    }
    
    df_out$procrustes_dist[i] <- sqrt(sum(d))
    
  }
  df_out <- df_out[order(df_out$procrustes_dist), ]
  
  
  return(df_out)
}

shape_dist(predicted.shapes)
