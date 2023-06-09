 Mtree <- read.nexus()
 Mtreedichotomy <- multi2di(Mtree, random = TRUE)
 Tshape <- as.treeshape(Mtreedichotomy)
 Ttips <- Ntip(Mtreedichotomy)
 shapevalue <- colless(Tshape)
 maxvalue <- ((Ttips-1)*(Ttips-2))/2
 N.shapevalue <- shapevalue/maxvalue
 N.shapevalue
