# PCA intervals' merhod
# powered by jordi.marce.nogue@uni-hamburg.de

devtools::install_github("kassambara/factoextra") 


library(FactoMineR)
library(factoextra)

setwd("F:/Uni stuff/MSc/Bristol return - updated version/Jordi folder/intervals-method-sam/R-part")

stress.distrib = read.csv("intervals-models-sam-nowhales.csv", row.names=1, header = TRUE, sep = ",")


# Multivariate analysis PCA 
#--------------------------------------------------------------------------------------------------------

# ho fem de sb.stress

col.number = ncol(stress.distrib)
PCA.stress <- PCA(stress.distrib[,1:col.number],  graph = FALSE)

# colors by group
group.colors = row.names(stress.distrib)
#group.palette = c("chartreuse","cadetblue4","darkseagreen4","rosybrown3")

#colors by variable
interval.number = nrow(PCA.stress$var$coord)
interval.vector = seq(from = 1, to = interval.number, by=1 )
interval.colors = interval.vector
interval.palette = c("blue","cyan","green","chartreuse","yellow","gold","orange","red")

# Biplot: variables coloured by contribution to PCs

library(ggpubr)

fviz_pca_biplot(PCA.stress, 
                mean.point=F,                     # no vui la mitja de cada esp??cie
                axes.linetype = "solid",
    
                # Fill individuals by groups
                geom.ind=c("point", "text"),
                pointshape = 21, 
                pointsize = 5,            
                col.ind= "black",                 # que no hi hagi el cercle exterior
                fill.ind = group.colors,          # fill amb els grups
                alpha.ind = 1,                    # transparency, no = 1
                
                # Color variable by intervals
                geom.var = "arrow",
                col.var = interval.colors,
                arrowsize = 0.5,
                
                repel = TRUE) +                   # Aboid label overplotting
  
  
  #fill_palette(group.palette) +
  gradient_color(interval.palette)+
  labs(x = "PC1", y = "PC2")+
  
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

summary(PCA.stress)





