library(reshape2)
library(plyr)
library(dplyr)
library(rgl)
library(ggpubr)
library(cowplot)


# 1) WIREFRAME PLOT FUNCTION
# plot wireframe from an n x 3 matrix A and a W matrix
plot.coords <- function(A, W, points.col="black", points.cex=1, lines.col="black", lines.wd=2, bg.col=NULL, 
                        main=NULL, main.line=2, main.cex=2, legend=NULL, legend.pos="topright", legend.title="", 
                        legend.col=NULL, legend.cex=1.2, legend.lwd=2, legend.bty="n", params=NULL, add=FALSE) {
  if (!is.null(params)) {par3d(params)}
  points.col <- rep(points.col, length.out=nrow(A))
  points.cex <- rep(points.cex, length.out=nrow(A))
  lines.col <- rep(lines.col, length.out=nrow(W))
  lines.wd <- rep(lines.wd, length.out=nrow(W))
  if (!is.null(bg.col)) rgl.bg(sphere=TRUE, color=bg.col, lit=FALSE, back="fill")
  plot3d(A, type="s", col=points.col, xlab="", ylab="", zlab="", size=points.cex, aspect=FALSE, box=F, axes=F, add=add)
  if (!is.null(main) | !is.null(legend)) {
    if (!is.null(legend) & is.null(legend.col)) stop("must supply legend colors")
    bgplot3d({plot.new()
      if (!is.null(main)) title(main=main, line=main.line, cex.main=main.cex)
      if (!is.null(legend)) legend(legend.pos, title=legend.title, legend=legend, col=legend.col, lwd=legend.lwd, cex=legend.cex, bty=legend.bty)})}
  for (i in 1:nrow(W)) {
    segments3d(rbind(A[W[i,1],], A[W[i,2],]), lwd=lines.wd[i], col=lines.col[i])
  }
}



# 2) JAW LEVER GRAPHS 3D
#
# 2.1) Tangent PCA Op4.3dC, population means 
Op4.3dC.tanspace<-plotTangentSpace(Op4.3dC.res.geodf$coords, axis1 = 1, axis2 = 2, warpgrids = T, 
                                   groups = ldmk.factors$Population)

Op4.tan.group.1 <- Op4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC1 = mean(PC1))
Op4.tan.group.2 <- Op4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC2 = mean(PC2))
Op4.tan.group.3 <- Op4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC1max = max(PC1))
Op4.tan.group.4 <- Op4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC1min = min(PC1))
Op4.tan.group.5 <- Op4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC2max = max(PC2))
Op4.tan.group.6 <- Op4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC2min = min(PC2))
Op4.tan.group<-cbind(Op4.tan.group.1,Op4.tan.group.2[2],Op4.tan.group.3[2],Op4.tan.group.4[2],Op4.tan.group.5[2],Op4.tan.group.6[2])
Op4.tan.group<-distinct(base::merge(Op4.tan.dat[c("Habitat", "Watershed", "Population")],Op4.tan.group, by = "Population", all.x = F))

# 2.2) PCA means plot - Figure 2A (right)
Op4.tan.mean.points<-ggplot(data = Op4.tan.group, aes(x = PC1, y = PC2)) + 
  aes(colour=Habitat, shape=Watershed, stroke = 2)+ 
  scale_shape_manual(values = c(17,16,15,3,4,18)) +
  scale_color_manual(values = c("slategray2","khaki4"))+
  geom_line(aes(x = PC1, y = PC2,group = Watershed, linetype = "solid"), size = .5, color = "black") +
  scale_linetype(guide = F) + 
  geom_point(size = 1.75, fill = "white") + xlab("")+ylab("PC2")
Op4.tan.mean.points

# 2.3) PCA x LD - Figure 2A (left)
Op4.3dC.WxLDplot<- ggplot(data = Op4.dat<-data.frame(Op4.3dC.tanspace$pc.scores,ldmk.factors), 
                          aes(x=LD1, y=PC1))+geom_point(alpha = .3, stroke = 2)+
  xlab("LD1")+ylab("PC1") #+ #geom_density2d(alpha=.10) +
#stat_ellipse(type = "norm", level = .95)
Op4.group <- Op4.dat %>% group_by(Population) %>%dplyr::summarize(PC1 = mean(PC1))
Op4.group<-distinct(base::merge(Op4.dat[c("LD1","Habitat", "Watershed", "Population")],Op4.group, by = "Population", all.x = F))
Op4.3dC.WxLDplot<- Op4.3dC.WxLDplot+ aes(colour=Habitat, shape=Watershed, stroke = 1)+ 
  scale_shape_manual(values = c("Beaver" =17,"Boot"=16,
                                "Misty"=15,"Pye"=3,"Roberts"=4,"VillageBay"=18)) + xlab("") +
  scale_color_manual(values = c("slategray2","khaki4"))
Op4.3dC.WxLDplot<-Op4.3dC.WxLDplot +  
  geom_line(data = Op4.group, aes(x = LD1, y = PC1, group = Watershed), size = .5, color = "black") +
  geom_point(data = Op4.group, size = 3, stroke = 2.5)  
Op4.3dC.WxLDplot


####

###



# 2.4) Tangent PCA Dent.3dC, population means
Dent.3dC.tanspace<-plotTangentSpace(Dent.3dC.res.geodf$coords, axis1 = 1, axis2 = 2, warpgrids = T, 
                                    groups = ldmk.factors$Population)

Dent.tan.group.1 <- Dent.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC1 = mean(PC1))
Dent.tan.group.2 <- Dent.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC2 = mean(PC2))
Dent.tan.group.3 <- Dent.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC1max = max(PC1))
Dent.tan.group.4 <- Dent.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC1min = min(PC1))
Dent.tan.group.5 <- Dent.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC2max = max(PC2))
Dent.tan.group.6 <- Dent.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC2min = min(PC2))
Dent.tan.group<-cbind(Dent.tan.group.1,Dent.tan.group.2[2],Dent.tan.group.3[2],Dent.tan.group.4[2],Dent.tan.group.5[2],Dent.tan.group.6[2])
Dent.tan.group<-distinct(base::merge(Dent.tan.dat[c("Habitat", "Watershed", "Population")],Dent.tan.group, by = "Population", all.x = F))

# 2.5) PCA means plot - Figure 2B (right)
Dent.tan.mean.points<-ggplot(data = Dent.tan.group, aes(x = PC1, y = PC2)) + 
  aes(colour=Habitat, shape=Watershed, stroke = 2)+ 
  scale_shape_manual(values = c(17,16,15,3,4,18)) +
  scale_color_manual(values = c("slategray2","khaki4"))+
  geom_line(aes(x = PC1, y = PC2,group = Watershed, linetype = "solid"), size = .5, color = "black") +
  scale_linetype(guide = F) + 
  geom_point(size = 1.75, fill = "white") + xlab("")+ylab("PC2")
Dent.tan.mean.points

# 2.6) Figure 2B (left)
Dent.3dC.WxLDplot<- ggplot(data = Dent.dat<-data.frame(Dent.3dC.tanspace$pc.scores,ldmk.factors), 
                           aes(x=LD1, y=PC1))+geom_point(alpha = .3, stroke = 2)+
  xlab("LD1")+ylab("PC1") #+ #geom_density2d(alpha=.10) +
#stat_ellipse(type = "norm", level = .95)
Dent.group <- Dent.dat %>% group_by(Population) %>%dplyr::summarize(PC1 = mean(PC1))
Dent.group<-distinct(base::merge(Dent.dat[c("LD1","Habitat", "Watershed", "Population")],Dent.group, by = "Population", all.x = F))
Dent.3dC.WxLDplot<- Dent.3dC.WxLDplot+ aes(colour=Habitat, shape=Watershed, stroke = 1)+ 
  scale_shape_manual(values = c("Beaver" =17,"Boot"=16,
                                "Misty"=15,"Pye"=3,"Roberts"=4,"VillageBay"=18)) + xlab("") +
  scale_color_manual(values = c("slategray2","khaki4"))
Dent.3dC.WxLDplot<-Dent.3dC.WxLDplot +  
  geom_line(data = Dent.group, aes(x = LD1, y = PC1, group = Watershed), size = .5, color = "black") +
  geom_point(data = Dent.group, size = 3, stroke = 2.5)  
Dent.3dC.WxLDplot

###

# 2.7) Tangent PCA Max4.3dC, population means
Max4.3dC.tanspace<-plotTangentSpace(Max4.3dC.res.geodf$coords, axis1 = 1, axis2 = 2, warpgrids = T, 
                                    groups = ldmk.factors$Population)

Max4.tan.group.1 <- Max4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC1 = mean(PC1))
Max4.tan.group.2 <- Max4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC2 = mean(PC2))
Max4.tan.group.3 <- Max4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC1max = max(PC1))
Max4.tan.group.4 <- Max4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC1min = min(PC1))
Max4.tan.group.5 <- Max4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC2max = max(PC2))
Max4.tan.group.6 <- Max4.tan.dat %>% group_by(Population) %>%dplyr::summarize(PC2min = min(PC2))
Max4.tan.group<-cbind(Max4.tan.group.1,Max4.tan.group.2[2],Max4.tan.group.3[2],Max4.tan.group.4[2],Max4.tan.group.5[2],Max4.tan.group.6[2])
Max4.tan.group<-distinct(base::merge(Max4.tan.dat[c("Habitat", "Watershed", "Population")],Max4.tan.group, by = "Population", all.x = F))

# 2.8) PCA means plot - Figure 2C (right)
Max4.tan.mean.points<-ggplot(data = Max4.tan.group, aes(x = PC1, y = PC2)) + 
  aes(colour=Habitat, shape=Watershed, stroke = 2)+ 
  scale_shape_manual(values = c(17,16,15,3,4,18)) +
  scale_color_manual(values = c("slategray2","khaki4"))+
  geom_line(aes(x = PC1, y = PC2,group = Watershed, linetype = "solid"), size = .5, color = "black") +
  scale_linetype(guide = F) + 
  geom_point(size = 1.75, fill = "white") + xlab("PC1")+ylab("PC2")
Max4.tan.mean.points

# 2.9) Figure 2C (left)
Max4.3dC.WxLDplot<- ggplot(data = Max4.dat<-data.frame(Max4.3dC.tanspace$pc.scores,ldmk.factors), 
                           aes(x=LD1, y=PC1))+geom_point(alpha = .3, stroke = 2)+
  xlab("LD")+ylab("PC1") #+ #geom_density2d(alpha=.10) +
#stat_ellipse(type = "norm", level = .95)
Max4.group <- Max4.dat %>% group_by(Population) %>%dplyr::summarize(PC1 = mean(PC1))
Max4.group<-distinct(base::merge(Max4.dat[c("LD1","Habitat", "Watershed", "Population")],Max4.group, by = "Population", all.x = F))
Max4.3dC.WxLDplot<- Max4.3dC.WxLDplot+ aes(colour=Habitat, shape=Watershed, stroke = 1)+ 
  scale_shape_manual(values = c("Beaver" =17,"Boot"=16,
                                "Misty"=15,"Pye"=3,"Roberts"=4,"VillageBay"=18)) + xlab("Diet LD") +
  scale_color_manual(values = c("slategray2","khaki4"))
Max4.3dC.WxLDplot<-Max4.3dC.WxLDplot +  
  geom_line(data = Max4.group, aes(x = LD1, y = PC1, group = Watershed), size = .5, color = "black") +
  geom_point(data = Max4.group, size = 3, stroke = 2.5)  
Max4.3dC.WxLDplot










# 3) Reference wireframes for graphs
W.sym.tot<-cbind(c("vOpr","vOpr","Opr","12r","11r","11r","11r","6r","10r","vOpl","vOpl","Opl","12l","11l","11l","11l","6l","10l","3r","2r","11r"),
                 c("Opr","12r","11r","11r","2r","6r","10r","3r","3r","Opl","12l","11l","11l","2l","6l","10l","3l","3l","3l","2l","11l"))
W.Op4<-rbind(W.sym.tot[1:4,],W.sym.tot[10:13,])
W.Dent<-rbind(W.sym.tot[4:6,], W.sym.tot[13:15,])
W.Max4<-rbind(W.sym.tot[6:9,], W.sym.tot[15:18,])

# 3.1) Op4 3D min and max PC references
row.names(Op4bar.3d.gpa$consensus)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
Op4.3d.pc1minpos<-Op4bar.3d.gpa$consensus[1:4,]+Op4.3dC.tanspace$pc.shapes$PC1min
Op4.3d.pc1minneg<-Op4.3d.pc1minpos
Op4.3d.pc1minneg[,3]<-Op4.3d.pc1minpos[,3]*-1
Op4.3d.pc1min<-abind(Op4.3d.pc1minpos,Op4.3d.pc1minneg, along = 1)
Op4.3d.pc1maxpos<-Op4bar.3d.gpa$consensus[1:4,]+Op4.3dC.tanspace$pc.shapes$PC1max
Op4.3d.pc1maxneg<-Op4.3d.pc1maxpos
Op4.3d.pc1maxneg[,3]<-Op4.3d.pc1maxpos[,3]*-1
Op4.3d.pc1max<-abind(Op4.3d.pc1maxpos,Op4.3d.pc1maxneg, along = 1)
Op4.3d.pc2minpos<-Op4bar.3d.gpa$consensus[1:4,]+Op4.3dC.tanspace$pc.shapes$PC2min
Op4.3d.pc2minneg<-Op4.3d.pc2minpos
Op4.3d.pc2minneg[,3]<-Op4.3d.pc2minpos[,3]*-1
Op4.3d.pc2min<-abind(Op4.3d.pc2minpos,Op4.3d.pc2minneg, along = 1)
Op4.3d.pc2maxpos<-Op4bar.3d.gpa$consensus[1:4,]+Op4.3dC.tanspace$pc.shapes$PC2max
Op4.3d.pc2maxneg<-Op4.3d.pc2maxpos
Op4.3d.pc2maxneg[,3]<-Op4.3d.pc2maxpos[,3]*-1
Op4.3d.pc2max<-abind(Op4.3d.pc2maxpos,Op4.3d.pc2maxneg, along = 1)
row.names(Op4.3d.pc1min)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
row.names(Op4.3d.pc1max)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
row.names(Op4.3d.pc2min)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
row.names(Op4.3d.pc2max)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")


# 3.2) calibrate view
plot.coords(Op4bar.3d.gpa$consensus, W.Op4[,], points.col="black", lines.col="black", add = T)
view.Op<-par3d()
view2.Op<-par3d()

layout3d(matrix(1:4, nrow=2), sharedMouse = TRUE)
# OP4.3D PC axis 1
plot.coords(Op4.3d.pc1min, W.Op4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Op4bar.3d.gpa$consensus, W.Op4[], points.col="black", lines.col="black", add = T)
plot.coords(Op4.3d.pc1max, W.Op4[], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Op4bar.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view.Op$userMatrix)        
# OP4.3D PC axis 2
plot.coords(Op4.3d.pc2min, W.Op4[], points.col="slategray2", lines.col="slategray2")
plot.coords(Op4bar.3d.gpa$consensus, W.Op4[], points.col="black", lines.col="black", add = T)
plot.coords(Op4.3d.pc2max, W.Op4[], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Op4bar.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view.Op$userMatrix)        
# OP4.3D PC axis 1
plot.coords(Op4.3d.pc1min, W.Op4[], points.col="slategray2", lines.col="slategray2")
plot.coords(Op4bar.3d.gpa$consensus, W.Op4[], points.col="black", lines.col="black", add = T)
plot.coords(Op4.3d.pc1max, W.Op4[], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Op4bar.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view2.Op$userMatrix)        
# OP4.3D PC axis 2
plot.coords(Op4.3d.pc2min, W.Op4[], points.col="slategray2", lines.col="slategray2")
plot.coords(Op4bar.3d.gpa$consensus, W.Op4[], points.col="black", lines.col="black", add = T)
plot.coords(Op4.3d.pc2max, W.Op4[], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Op4bar.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view2.Op$userMatrix)  

# 3.3) Dent 3D PC reference wireframes
rownames(Dent.3d.gpa$consensus)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
Dent.3d.pc1minpos<-Dent.3d.gpa$consensus[1:4,]+Dent.3dC.tanspace$pc.shapes$PC1min
Dent.3d.pc1minneg<-Dent.3d.pc1minpos
Dent.3d.pc1minneg[,3]<-Dent.3d.pc1minpos[,3]*-1
Dent.3d.pc1min<-abind(Dent.3d.pc1minpos,Dent.3d.pc1minneg, along = 1)
Dent.3d.pc1maxpos<-Dent.3d.gpa$consensus[1:4,]+Dent.3dC.tanspace$pc.shapes$PC1max
Dent.3d.pc1maxneg<-Dent.3d.pc1maxpos
Dent.3d.pc1maxneg[,3]<-Dent.3d.pc1maxpos[,3]*-1
Dent.3d.pc1max<-abind(Dent.3d.pc1maxpos,Dent.3d.pc1maxneg, along = 1)
Dent.3d.pc2minpos<-Dent.3d.gpa$consensus[1:4,]+Dent.3dC.tanspace$pc.shapes$PC2min
Dent.3d.pc2minneg<-Dent.3d.pc2minpos
Dent.3d.pc2minneg[,3]<-Dent.3d.pc2minpos[,3]*-1
Dent.3d.pc2min<-abind(Dent.3d.pc2minpos,Dent.3d.pc2minneg, along = 1)
Dent.3d.pc2maxpos<-Dent.3d.gpa$consensus[1:4,]+Dent.3dC.tanspace$pc.shapes$PC2max
Dent.3d.pc2maxneg<-Dent.3d.pc2maxpos
Dent.3d.pc2maxneg[,3]<-Dent.3d.pc2maxpos[,3]*-1
Dent.3d.pc2max<-abind(Dent.3d.pc2maxpos,Dent.3d.pc2maxneg, along = 1)
row.names(Dent.3d.pc1min)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
row.names(Dent.3d.pc1max)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
row.names(Dent.3d.pc2min)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
row.names(Dent.3d.pc2max)<-c("12r","11r","6r","2r","12l","11l","6l","2l")

# 3.4) calibrate view
plot.coords(Dent.3d.gpa$consensus, W.Dent[], points.col="black", lines.col="black", add = T)
text3d(Dent.3d.gpa$consensus,adj = c(2,2),text = c("vOp","Op-H","An","A-Q"), col = "black", alpha = 1)
view.Dent<-par3d()
view2.Dent<-par3d()

layout3d(matrix(1:4, nrow=2), sharedMouse = TRUE)
# Dent.3D PC axis 1
plot.coords(Dent.3d.pc1min, W.Dent, points.col="slategray2", lines.col="slategray2")
plot.coords(Dent.3d.gpa$consensus, W.Dent[,], points.col="black", lines.col="black", add = T)
plot.coords(Dent.3d.pc1max, W.Dent[,], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Dent.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view.Dent$userMatrix)        
# Dent.3D PC axis 2
plot.coords(Dent.3d.pc2min, W.Dent[], points.col="slategray2", lines.col="slategray2")
plot.coords(Dent.3d.gpa$consensus, W.Dent[], points.col="black", lines.col="black", add = T)
plot.coords(Dent.3d.pc2max, W.Dent[], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Dent.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view.Dent$userMatrix)        
# Dent.3D PC axis 1
plot.coords(Dent.3d.pc1min, W.Dent, points.col="slategray2", lines.col="slategray2")
plot.coords(Dent.3d.gpa$consensus, W.Dent[,], points.col="black", lines.col="black", add = T)
plot.coords(Dent.3d.pc1max, W.Dent[,], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Dent.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view2.Dent$userMatrix)        
# Dent.3D PC axis 2
plot.coords(Dent.3d.pc2min, W.Dent[], points.col="slategray2", lines.col="slategray2")
plot.coords(Dent.3d.gpa$consensus, W.Dent[], points.col="black", lines.col="black", add = T)
plot.coords(Dent.3d.pc2max, W.Dent[], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Dent.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view2.Dent$userMatrix) 


# 3.5) Max4 3D PC reference wireframes
rownames(Max4bar.3d.gpa$consensus)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
Max4.3d.pc1minpos<-Max4bar.3d.gpa$consensus[1:4,]+Max4.3dC.tanspace$pc.shapes$PC1min
Max4.3d.pc1minneg<-Max4.3d.pc1minpos
Max4.3d.pc1minneg[,3]<-Max4.3d.pc1minpos[,3]*-1
Max4.3d.pc1min<-abind(Max4.3d.pc1minpos,Max4.3d.pc1minneg, along = 1)
Max4.3d.pc1maxpos<-Max4bar.3d.gpa$consensus[1:4,]+Max4.3dC.tanspace$pc.shapes$PC1max
Max4.3d.pc1maxneg<-Max4.3d.pc1maxpos
Max4.3d.pc1maxneg[,3]<-Max4.3d.pc1maxpos[,3]*-1
Max4.3d.pc1max<-abind(Max4.3d.pc1maxpos,Max4.3d.pc1maxneg, along = 1)
Max4.3d.pc2minpos<-Max4bar.3d.gpa$consensus[1:4,]+Max4.3dC.tanspace$pc.shapes$PC2min
Max4.3d.pc2minneg<-Max4.3d.pc2minpos
Max4.3d.pc2minneg[,3]<-Max4.3d.pc2minpos[,3]*-1
Max4.3d.pc2min<-abind(Max4.3d.pc2minpos,Max4.3d.pc2minneg, along = 1)
Max4.3d.pc2maxpos<-Max4bar.3d.gpa$consensus[1:4,]+Max4.3dC.tanspace$pc.shapes$PC2max
Max4.3d.pc2maxneg<-Max4.3d.pc2maxpos
Max4.3d.pc2maxneg[,3]<-Max4.3d.pc2maxpos[,3]*-1
Max4.3d.pc2max<-abind(Max4.3d.pc2maxpos,Max4.3d.pc2maxneg, along = 1)
row.names(Max4.3d.pc1min)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
row.names(Max4.3d.pc1max)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
row.names(Max4.3d.pc2min)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
row.names(Max4.3d.pc2max)<-c("11r","6r","3r","10r","11l","6l","3l","10l")

# 3.6) calibrate view
plot.coords(Max4bar.3d.gpa$consensus, W.Max4[,], points.col="black", lines.col="black", add = T)
text3d(Max4bar.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
view.Max<-par3d()
view2.Max<-par3d()

layout3d(matrix(1:4, nrow=2), sharedMouse = TRUE)
# Max4.3D PC axis 1
plot.coords(Max4.3d.pc1min, W.Max4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Max4bar.3d.gpa$consensus, W.Max4[,], points.col="black", lines.col="black", add = T)
plot.coords(Max4.3d.pc1max, W.Max4[,], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Max4bar.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view.Max$userMatrix)           
# Max4.3D PC axis 2
plot.coords(Max4.3d.pc2min, W.Max4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Max4bar.3d.gpa$consensus, W.Max4[,], points.col="black", lines.col="black", add = T)
plot.coords(Max4.3d.pc2max, W.Max4[,], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Max4bar.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view.Max$userMatrix)           
# Max4.3D PC axis 1
plot.coords(Max4.3d.pc1min, W.Max4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Max4bar.3d.gpa$consensus, W.Max4[,], points.col="black", lines.col="black", add = T)
plot.coords(Max4.3d.pc1max, W.Max4[,], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Max4bar.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view2.Max$userMatrix)           
# Max4.3D PC axis 2
plot.coords(Max4.3d.pc2min, W.Max4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Max4bar.3d.gpa$consensus, W.Max4[,], points.col="black", lines.col="black", add = T)
plot.coords(Max4.3d.pc2max, W.Max4[,], points.col="sienna1", lines.col="sienna1", add = T)
text3d(Max4bar.3d.gpa$consensus,adj = c(2,2), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = view2.Max$userMatrix)        


# 3.7) Figure 2
ggarrange(Op4.3dC.WxLDplot,Op4.tan.mean.points, Dent.3dC.WxLDplot,
          Dent.tan.mean.points, Max4.3dC.WxLDplot,Max4.tan.mean.points, 
          labels = c("A","", "B","", "C",""),
          ncol = 2, nrow = 3, common.legend = T, legend = "right")






#
#
#
# 4) Plots of pairwise angles 
# 4.1) Figure 3A
Op4.angle<-na.exclude(data.frame(Op4.3dC.angle.melted, Op4.2dC.angle.melted[3]))
Op4.angle.plot<-ggplot(data = Op4.angle, 
                       aes(x=value, y=value.1)) + geom_abline(slope = 1, intercept = 30, color = "red1") + 
  geom_abline(slope = 1, intercept = -30, color = "red1") +geom_point(aes(shape=X1), stroke = 1, color = "steelblue4", size = 2.5)+
  geom_point(aes(shape=X2), stroke = .75, color = "sienna1", size = 1.5)+
  xlab("")+ylab("2D Angle") + geom_abline(slope = 1, intercept = 0)
Op4.angle.plot<- Op4.angle.plot+ scale_shape_manual(values = c("Beaver" =8,"Boot"=3,
                                                               "Misty"=4,"Pye"=16,"Roberts"=17,"VillageBay"=15)) +
  scale_shape_manual(values = c("Beaver" =8,"Boot"=3,
                                "Misty"=4,"Pye"=16,"Roberts"=17,"VillageBay"=15)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,180)) + scale_y_continuous(expand = c(0, 0), limits = c(0,180)) +
  guides(colour = guide_legend(override.aes = list(shape = 19))) + 
  labs(shape = "Watershed")
Op4.angle.plot
# 4.2) Figure 3B
Dent.angle<-na.exclude(data.frame(Dent.3dC.angle.melted, Dent.2dC.angle.melted[3]))
Dent.angle.plot<-ggplot(data = Dent.angle, 
                        aes(x=value, y=value.1)) + geom_abline(slope = 1, intercept = 30, color = "red1") + 
  geom_abline(slope = 1, intercept = -30, color = "red1") +geom_point(aes(shape=X1), stroke = 1, color = "steelblue4", size = 2.5)+
  geom_point(aes(shape=X2), stroke = .75, color = "sienna1", size = 1.5)+
  xlab("")+ylab("2D Angle") + geom_abline(slope = 1, intercept = 0)
Dent.angle.plot<- Dent.angle.plot+ scale_shape_manual(values = c("Beaver" =8,"Boot"=3,
                                                                 "Misty"=4,"Pye"=16,"Roberts"=17,"VillageBay"=15)) +
  scale_shape_manual(values = c("Beaver" =8,"Boot"=3,
                                "Misty"=4,"Pye"=16,"Roberts"=17,"VillageBay"=15)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,180)) + scale_y_continuous(expand = c(0, 0), limits = c(0,180)) +
  guides(colour = guide_legend(override.aes = list(shape = 19))) + 
  labs(shape = "Watershed")
Dent.angle.plot
# 4.3) Figure 3C
Max4bar.angle<-na.exclude(data.frame(Max4bar.3dC.angle.melted, Max4bar.2dC.angle.melted[3]))
Max4bar.angle.plot<-ggplot(data = Max4bar.angle, 
                           aes(x=value, y=value.1)) + geom_abline(slope = 1, intercept = 30, color = "red1") + 
  geom_abline(slope = 1, intercept = -30, color = "red1") +geom_point(aes(shape=X1), stroke = 1, color = "steelblue4", size = 2.5)+
  geom_point(aes(shape=X2), stroke = .75, color = "sienna1", size = 1.5)+
  xlab("3D Angle")+ylab("2D Angle") + geom_abline(slope = 1, intercept = 0)
Max4bar.angle.plot<- Max4bar.angle.plot+ scale_shape_manual(values = c("Beaver" =8,"Boot"=3,
                                                                       "Misty"=4,"Pye"=16,"Roberts"=17,"VillageBay"=15)) +
  scale_shape_manual(values = c("Beaver" =8,"Boot"=3,
                                "Misty"=4,"Pye"=16,"Roberts"=17,"VillageBay"=15)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,180)) + scale_y_continuous(expand = c(0, 0), limits = c(0,180)) +
  guides(colour = guide_legend(override.aes = list(shape = 19))) + 
  labs(shape = "Watershed")
Max4bar.angle.plot

# 4.4) Figure 3
ggarrange(Op4.angle.plot, Dent.angle.plot, Max4bar.angle.plot,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3, common.legend = T, legend = "right")

#
#
#
#
#
#
# 5) Correlations between 2d and 3d procdist for B/W
# 5.1) Figure S7A
Op4.procdist
Op4.procdist.plot.BW<-ggplot(data = Op4.procdist, 
                             aes(x=procdist.3d, y=procdist.2d))+geom_point(aes(shape= watershed, color = habitat), stroke = 1, size = 2.5)+
  xlab("")+ylab("2D Proc. Distance") + geom_abline(slope = 1, intercept = 0) 
Op4.procdist.plot.BW<- Op4.procdist.plot.BW+
  scale_shape_manual(values = c("Beaver" =1,"Boot"=2,"Misty"=0,"Pye"=16,"Roberts"=17,"VillageBay"=15)) +
  scale_color_manual(values = c("steelblue4","sienna1")) 
Op4.procdist.plot.BW
# 5.2) Figure S7B
Dent.procdist
Dent.procdist.plot.BW<-ggplot(data = Dent.procdist, 
                              aes(x=procdist.3d, y=procdist.2d))+geom_point(aes(shape= watershed, color = habitat), stroke = 1, size = 2.5)+
  xlab("")+ylab("2D Proc. Distance") + geom_abline(slope = 1, intercept = 0) 
Dent.procdist.plot.BW<- Dent.procdist.plot.BW+ 
  scale_shape_manual(values = c("Beaver" =1,"Boot"=2,"Misty"=0,"Pye"=16,"Roberts"=17,"VillageBay"=15)) +
  scale_color_manual(values = c("steelblue4","sienna1")) 
Dent.procdist.plot.BW
# 5.3) Figure S7C
Max4bar.procdist
Max4bar.procdist.plot.BW<-ggplot(data = Max4bar.procdist, 
                                 aes(x=procdist.3d, y=procdist.2d))+geom_point(aes(shape= watershed, color = habitat), stroke = 1, size = 2.5)+
  xlab("3D Proc. Distance")+ylab("2D Proc. Distance") + geom_abline(slope = 1, intercept = 0) 
Max4bar.procdist.plot.BW<- Max4bar.procdist.plot.BW+ 
  scale_shape_manual(values = c("Beaver" =1,"Boot"=2,"Misty"=0,"Pye"=16,"Roberts"=17,"VillageBay"=15)) +
  scale_color_manual(values = c("steelblue4","sienna1")) 
Max4bar.procdist.plot.BW

# 5.4) Figure S7
ggarrange(Op4.procdist.plot.BW, Dent.procdist.plot.BW, Max4bar.procdist.plot.BW,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3, common.legend = T, legend = "right")


# 6) Figure S8
ggarrange(
  ggarrange(
    ggplot(data =  kinem.df.res, aes(x = Watershed, y = Op4.2d.KT.res, fill = Habitat)) +
      geom_boxplot(width = 0.4, outlier.shape = NA) + 
      ylab("Op4 KT Res.")+ xlab(NULL)+
      scale_color_manual(values = c("steelblue4","red3")) + scale_fill_manual(values = c("slategray2","darkolivegreen3")) + 
      theme(axis.ticks.x = element_blank()),
    ggplot(data =  kinem.df.res, aes(x = Watershed, y = Max4.2d.KT.res, fill = Habitat)) +
      geom_boxplot(width = 0.4, outlier.shape = NA) + 
      ylab("Max4 KT Res.")+ xlab(NULL)+
      scale_color_manual(values = c("steelblue4","red3")) + scale_fill_manual(values = c("slategray2","darkolivegreen3")) + 
      theme(axis.ticks.x = element_blank()),
    ggplot(data =  kinem.df.res, aes(x = Watershed, y = Mand.2d.closingLR.res, fill = Habitat)) +
      geom_boxplot(width = 0.4, outlier.shape = NA) + 
      ylab("Mand Closing LR")+ xlab(NULL)+
      scale_color_manual(values = c("steelblue4","red3")) + scale_fill_manual(values = c("slategray2","darkolivegreen3")) + 
      theme(axis.ticks.x = element_blank()),
    ggplot(data =  kinem.df.res, aes(x = Watershed, y = Mand.2d.openingLR.res, fill = Habitat)) +
      geom_boxplot(width = 0.4, outlier.shape = NA) + 
      ylab("Mand Opening LR")+ xlab(NULL)+
      scale_color_manual(values = c("steelblue4","red3")) + scale_fill_manual(values = c("slategray2","darkolivegreen3")) + 
      theme(axis.ticks.x = element_blank()),
    ncol = 1, nrow = 4, labels = c("A","B","C","D"), legend = "right", common.legend = T),
  ggarrange(
    ggplot(data =  kinem.df.res, aes(x = Watershed, y = Op4.3d.KT.30.res, fill = Habitat)) +
      geom_boxplot(width = 0.4, outlier.shape = NA) + 
      ylab(NULL)+ xlab(NULL)+
      scale_color_manual(values = c("steelblue4","red3")) + scale_fill_manual(values = c("slategray2","darkolivegreen3")) + 
      theme(axis.ticks.x = element_blank()),
    ggplot(data =  kinem.df.res, aes(x = Watershed, y = Max4.3d.KT.20.res, fill = Habitat)) +
      geom_boxplot(width = 0.4, outlier.shape = NA) + 
      ylab(NULL)+ xlab(NULL)+
      scale_color_manual(values = c("steelblue4","red3")) + scale_fill_manual(values = c("slategray2","darkolivegreen3")) + 
      theme(axis.ticks.x = element_blank()),
    ggplot(data =  kinem.df.res, aes(x = Watershed, y = Mand.3d.closingLR.res, fill = Habitat)) +
      geom_boxplot(width = 0.4, outlier.shape = NA) + 
      ylab(NULL)+ xlab(NULL)+
      scale_color_manual(values = c("steelblue4","red3")) + scale_fill_manual(values = c("slategray2","darkolivegreen3")) + 
      theme(axis.ticks.x = element_blank()),
    ggplot(data =  kinem.df.res, aes(x = Watershed, y = Mand.3d.openingLR.res, fill = Habitat)) +
      geom_boxplot(width = 0.4, outlier.shape = NA) + 
      ylab(NULL)+ xlab(NULL)+
      scale_color_manual(values = c("steelblue4","red3")) + scale_fill_manual(values = c("slategray2","darkolivegreen3")) + 
      theme(axis.ticks.x = element_blank()),
    ncol = 1, nrow = 4, labels = c("A","B","C","D"), legend = "right", common.legend = T),
  ncol = 2, nrow = 1, legend = "right", common.legend = T)





