library(rgl)

# WIREFRAMES FOR FIGURE 4


# Inherited from "1_SeparateShapeSpaces.R"
# Op4bar.3d.gpa
# Dent.3d.gpa
# Max4bar.3d.gpa
# spec.sym.3d.gpa


# Inherited from "2-CommonShapeSpace.R"
# Op4.2d3d.tanpca
# Dent.2d3d.tanpca
# Max4.2d3d.tanpca
# tot.2d3d.tanpca


# 1) wireframe function (for more info and explanation of commands within function, 
#    see http://www.randigriffin.com/2017/11/10/plotting-shape-changes-geomorph.html)
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

# 1.1) "W" commands define which points are connected by lines
W.sym.tot<-cbind(c("vOpr","vOpr","Opr","12r","11r","11r","11r","6r","10r","vOpl","vOpl","Opl","12l","11l","11l","11l","6l","10l","3r","2r","11r"),
                 c("Opr","12r","11r","11r","2r","6r","10r","3r","3r","Opl","12l","11l","11l","2l","6l","10l","3l","3l","3l","2l","11l"))
W.tot<-cbind(c("Op4bar","Op4bar" ,"Op4bar" ,"Mand" ,"Mand" ,"Mand" ,"Max4" ,"Max4" ,"Max4" ),c("vOp","vOp","Op","12","11","11","11","6","10"),c("Op","12","11","11","2","6","10","3","3"))
W.Op4<-rbind(W.sym.tot[1:4,],W.sym.tot[10:13,])
W.Dent<-rbind(W.sym.tot[4:6,], W.sym.tot[13:15,])
W.Max4<-rbind(W.sym.tot[6:9,], W.sym.tot[15:18,])

# 1.2) Op4 - defining coordinates for max and min along each PC
row.names(Op4bar.3d.gpa$consensus)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
Op4.2d3d.pc1minpos<-Op4bar.3d.gpa$consensus[1:4,]+Op4.2d3d.tanpca$pc.shapes$PC1min
Op4.2d3d.pc1minneg<-Op4.2d3d.pc1minpos
Op4.2d3d.pc1minneg[,3]<-Op4.2d3d.pc1minpos[,3]*-1
Op4.2d3d.pc1min<-abind(Op4.2d3d.pc1minpos,Op4.2d3d.pc1minneg, along = 1)
Op4.2d3d.pc1maxpos<-Op4bar.3d.gpa$consensus[1:4,]+Op4.2d3d.tanpca$pc.shapes$PC1max
Op4.2d3d.pc1maxneg<-Op4.2d3d.pc1maxpos
Op4.2d3d.pc1maxneg[,3]<-Op4.2d3d.pc1maxpos[,3]*-1
Op4.2d3d.pc1max<-abind(Op4.2d3d.pc1maxpos,Op4.2d3d.pc1maxneg, along = 1)
Op4.2d3d.pc2minpos<-Op4bar.3d.gpa$consensus[1:4,]+Op4.2d3d.tanpca$pc.shapes$PC2min
Op4.2d3d.pc2minneg<-Op4.2d3d.pc2minpos
Op4.2d3d.pc2minneg[,3]<-Op4.2d3d.pc2minpos[,3]*-1
Op4.2d3d.pc2min<-abind(Op4.2d3d.pc2minpos,Op4.2d3d.pc2minneg, along = 1)
Op4.2d3d.pc2maxpos<-Op4bar.3d.gpa$consensus[1:4,]+Op4.2d3d.tanpca$pc.shapes$PC2max
Op4.2d3d.pc2maxneg<-Op4.2d3d.pc2maxpos
Op4.2d3d.pc2maxneg[,3]<-Op4.2d3d.pc2maxpos[,3]*-1
Op4.2d3d.pc2max<-abind(Op4.2d3d.pc2maxpos,Op4.2d3d.pc2maxneg, along = 1)
row.names(Op4.2d3d.pc1min)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
row.names(Op4.2d3d.pc1max)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
row.names(Op4.2d3d.pc2min)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")
row.names(Op4.2d3d.pc2max)<-c("vOpr","Opr","12r","11r","vOpl","Opl","12l","11l")

# 1.3) Same as 1.2 for mandibular levers
row.names(Dent.3d.gpa$consensus)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
Dent.2d3d.pc1minpos<-Dent.3d.gpa$consensus[1:4,]+Dent.2d3d.tanpca$pc.shapes$PC1min
Dent.2d3d.pc1minneg<-Dent.2d3d.pc1minpos
Dent.2d3d.pc1minneg[,3]<-Dent.2d3d.pc1minpos[,3]*-1
Dent.2d3d.pc1min<-abind(Dent.2d3d.pc1minpos,Dent.2d3d.pc1minneg, along = 1)
Dent.2d3d.pc1maxpos<-Dent.3d.gpa$consensus[1:4,]+Dent.2d3d.tanpca$pc.shapes$PC1max
Dent.2d3d.pc1maxneg<-Dent.2d3d.pc1maxpos
Dent.2d3d.pc1maxneg[,3]<-Dent.2d3d.pc1maxpos[,3]*-1
Dent.2d3d.pc1max<-abind(Dent.2d3d.pc1maxpos,Dent.2d3d.pc1maxneg, along = 1)
Dent.2d3d.pc2minpos<-Dent.3d.gpa$consensus[1:4,]+Dent.2d3d.tanpca$pc.shapes$PC2min
Dent.2d3d.pc2minneg<-Dent.2d3d.pc2minpos
Dent.2d3d.pc2minneg[,3]<-Dent.2d3d.pc2minpos[,3]*-1
Dent.2d3d.pc2min<-abind(Dent.2d3d.pc2minpos,Dent.2d3d.pc2minneg, along = 1)
Dent.2d3d.pc2maxpos<-Dent.3d.gpa$consensus[1:4,]+Dent.2d3d.tanpca$pc.shapes$PC2max
Dent.2d3d.pc2maxneg<-Dent.2d3d.pc2maxpos
Dent.2d3d.pc2maxneg[,3]<-Dent.2d3d.pc2maxpos[,3]*-1
Dent.2d3d.pc2max<-abind(Dent.2d3d.pc2maxpos,Dent.2d3d.pc2maxneg, along = 1)
row.names(Dent.2d3d.pc1min)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
row.names(Dent.2d3d.pc1max)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
row.names(Dent.2d3d.pc2min)<-c("12r","11r","6r","2r","12l","11l","6l","2l")
row.names(Dent.2d3d.pc2max)<-c("12r","11r","6r","2r","12l","11l","6l","2l")

# 1.4) Same as above for Max4
row.names(Max4bar.3d.gpa$consensus)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
Max4.2d3d.pc1minpos<-Max4bar.3d.gpa$consensus[1:4,]+Max4.2d3d.tanpca$pc.shapes$PC1min
Max4.2d3d.pc1minneg<-Max4.2d3d.pc1minpos
Max4.2d3d.pc1minneg[,3]<-Max4.2d3d.pc1minpos[,3]*-1
Max4.2d3d.pc1min<-abind(Max4.2d3d.pc1minpos,Max4.2d3d.pc1minneg, along = 1)
Max4.2d3d.pc1maxpos<-Max4bar.3d.gpa$consensus[1:4,]+Max4.2d3d.tanpca$pc.shapes$PC1max
Max4.2d3d.pc1maxneg<-Max4.2d3d.pc1maxpos
Max4.2d3d.pc1maxneg[,3]<-Max4.2d3d.pc1maxpos[,3]*-1
Max4.2d3d.pc1max<-abind(Max4.2d3d.pc1maxpos,Max4.2d3d.pc1maxneg, along = 1)
Max4.2d3d.pc2minpos<-Max4bar.3d.gpa$consensus[1:4,]+Max4.2d3d.tanpca$pc.shapes$PC2min
Max4.2d3d.pc2minneg<-Max4.2d3d.pc2minpos
Max4.2d3d.pc2minneg[,3]<-Max4.2d3d.pc2minpos[,3]*-1
Max4.2d3d.pc2min<-abind(Max4.2d3d.pc2minpos,Max4.2d3d.pc2minneg, along = 1)
Max4.2d3d.pc2maxpos<-Max4bar.3d.gpa$consensus[1:4,]+Max4.2d3d.tanpca$pc.shapes$PC2max
Max4.2d3d.pc2maxneg<-Max4.2d3d.pc2maxpos
Max4.2d3d.pc2maxneg[,3]<-Max4.2d3d.pc2maxpos[,3]*-1
Max4.2d3d.pc2max<-abind(Max4.2d3d.pc2maxpos,Max4.2d3d.pc2maxneg, along = 1)
row.names(Max4.2d3d.pc1min)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
row.names(Max4.2d3d.pc1max)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
row.names(Max4.2d3d.pc2min)<-c("11r","6r","3r","10r","11l","6l","3l","10l")
row.names(Max4.2d3d.pc2max)<-c("11r","6r","3r","10r","11l","6l","3l","10l")

# 1.5) Same as above for total landmark set
row.names(spec.sym.3d.gpa$consensus)<-c("vOpr","Opr","12r","11r","6r","2r","3r","10r","vOpl","Opl","12l","11l","6l","2l","3l","10l")
tot.2d3d.pc1minpos<-spec.sym.3d.gpa$consensus[1:8,]+tot.2d3d.tanpca$pc.shapes$PC1min
tot.2d3d.pc1minneg<-tot.2d3d.pc1minpos
tot.2d3d.pc1minneg[,3]<-tot.2d3d.pc1minpos[,3]*-1
tot.2d3d.pc1min<-abind(tot.2d3d.pc1minpos,tot.2d3d.pc1minneg, along = 1)
tot.2d3d.pc1maxpos<-spec.sym.3d.gpa$consensus[1:8,]+tot.2d3d.tanpca$pc.shapes$PC1max
tot.2d3d.pc1maxneg<-tot.2d3d.pc1maxpos
tot.2d3d.pc1maxneg[,3]<-tot.2d3d.pc1maxpos[,3]*-1
tot.2d3d.pc1max<-abind(tot.2d3d.pc1maxpos,tot.2d3d.pc1maxneg, along = 1)
tot.2d3d.pc2minpos<-spec.sym.3d.gpa$consensus[1:8,]+tot.2d3d.tanpca$pc.shapes$PC2min
tot.2d3d.pc2minneg<-tot.2d3d.pc2minpos
tot.2d3d.pc2minneg[,3]<-tot.2d3d.pc2minpos[,3]*-1
tot.2d3d.pc2min<-abind(tot.2d3d.pc2minpos,tot.2d3d.pc2minneg, along = 1)
tot.2d3d.pc2maxpos<-spec.sym.3d.gpa$consensus[1:8,]+tot.2d3d.tanpca$pc.shapes$PC2max
tot.2d3d.pc2maxneg<-tot.2d3d.pc2maxpos
tot.2d3d.pc2maxneg[,3]<-tot.2d3d.pc2maxpos[,3]*-1
tot.2d3d.pc2max<-abind(tot.2d3d.pc2maxpos,tot.2d3d.pc2maxneg, along = 1)
row.names(tot.2d3d.pc1min)<-c("vOpr","Opr","12r","11r","6r","2r","3r","10r","vOpl","Opl","12l","11l","6l","2l","3l","10l")
row.names(tot.2d3d.pc1max)<-c("vOpr","Opr","12r","11r","6r","2r","3r","10r","vOpl","Opl","12l","11l","6l","2l","3l","10l")
row.names(tot.2d3d.pc2min)<-c("vOpr","Opr","12r","11r","6r","2r","3r","10r","vOpl","Opl","12l","11l","6l","2l","3l","10l")
row.names(tot.2d3d.pc2max)<-c("vOpr","Opr","12r","11r","6r","2r","3r","10r","vOpl","Opl","12l","11l","6l","2l","3l","10l")

#
#
#
#
# 2) calibrate views 
plot.coords(Op4bar.3d.gpa$consensus, W.Op4[,], points.col="black", lines.col="black", add = T) #plots lever system, can rotate in RGL window
text3d(Op4bar.3d.gpa$consensus,adj = c(2,2),text = c("vOp","Op-H","An","A-Q"), col = "black", alpha = 1) #Adds joint labels
Opview<-par3d() # fixes parameters for first view
Opview2<-par3d() # fixes parameters for second view

plot.coords(Dent.3d.gpa$consensus, W.Dent[,], points.col="black", lines.col="black", add = T)
text3d(Dent.3d.gpa$consensus,adj = c(2,2),text = c("An","A-Q","Art","D"), col = "black", alpha = 1)
Dentview<-par3d()
Dentview2<-par3d()

plot.coords(Max4bar.3d.gpa$consensus, W.Max4[,], points.col="black", lines.col="black", add = T)
text3d(Max4bar.3d.gpa$consensus,adj = c(2,2),text = c("A-Q","ES","Art","Max"), col = "black", alpha = 1)
Max4view<-par3d()
Max4view2<-par3d()

plot.coords(spec.sym.3d.gpa$consensus, W.sym.tot[,], points.col="black", lines.col="black", add = T)
text3d(tot.2d3d.tanpca$consensus,adj = c(2,2),text = c("A-Q","ES","Art","Max"), col = "black", alpha = 1)
totview<-par3d()
totview2<-par3d()

#
#
#
# 3) plots wireframes with min, max, and consensus on PC1 and PC2 for each landmark set
layout3d(matrix(1:16, nrow=2), sharedMouse = TRUE) #puts subsequent wireframes in grid, sharedMouse allows simultaneous rotation
# 3.1)  OP4.3D PC axis 1
plot.coords(Op4.2d3d.pc1min, W.Op4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Op4bar.3d.gpa$consensus, W.Op4[,], points.col="black", lines.col="black", add = T)
plot.coords(Op4.2d3d.pc1max, W.Op4[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Op4bar.3d.gpa$consensus,adj = c(2,2),text = c("vOp","Op-H","An","A-Q"), col = "black", alpha = 1) #this line adds joint labels
rgl.viewpoint(userMatrix = Opview$userMatrix)  #  Orients wireframe to perspective defined by par3d() above    
# OP4.3D PC axis 2
plot.coords(Op4.2d3d.pc2min, W.Op4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Op4bar.3d.gpa$consensus, W.Op4[,], points.col="black", lines.col="black", add = T)
plot.coords(Op4.2d3d.pc2max, W.Op4[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Op4bar.3d.gpa$consensus,adj = c(2,2),text = c("vOp","Op-H","An","A-Q"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Opview$userMatrix)        
# OP4.3D PC axis 1
plot.coords(Op4.2d3d.pc1min, W.Op4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Op4bar.3d.gpa$consensus, W.Op4[,], points.col="black", lines.col="black", add = T)
plot.coords(Op4.2d3d.pc1max, W.Op4[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Op4bar.3d.gpa$consensus,adj = c(2,2),text = c("vOp","Op-H","An","A-Q"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Opview2$userMatrix)        
# OP4.3D PC axis 2
plot.coords(Op4.2d3d.pc2min, W.Op4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Op4bar.3d.gpa$consensus, W.Op4[,], points.col="black", lines.col="black", add = T)
plot.coords(Op4.2d3d.pc2max, W.Op4[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Op4bar.3d.gpa$consensus,adj = c(2,2),text = c("vOp","Op-H","An","A-Q"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Opview2$userMatrix)         

# 3.2) Dent.3D PC axis 1
plot.coords(Dent.2d3d.pc1min, W.Dent[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Dent.3d.gpa$consensus, W.Dent[,], points.col="black", lines.col="black", add = T)
plot.coords(Dent.2d3d.pc1max, W.Dent[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Dent.3d.gpa$consensus,adj = c(2,2),text = c("An","A-Q","Art","D"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Dentview$userMatrix)        
# Dent.3D PC axis 2
plot.coords(Dent.2d3d.pc2min, W.Dent[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Dent.3d.gpa$consensus, W.Dent[,], points.col="black", lines.col="black", add = T)
plot.coords(Dent.2d3d.pc2max, W.Dent[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Dent.3d.gpa$consensus,adj = c(2,2),text = c("An","A-Q","Art","D"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Dentview$userMatrix)        
# Dent.3D PC axis 1
plot.coords(Dent.2d3d.pc1min, W.Dent[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Dent.3d.gpa$consensus, W.Dent[,], points.col="black", lines.col="black", add = T)
plot.coords(Dent.2d3d.pc1max, W.Dent[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Dent.3d.gpa$consensus,adj = c(2,2),text = c("An","A-Q","Art","D"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Dentview2$userMatrix)        
# Dent.3D PC axis 2
plot.coords(Dent.2d3d.pc2min, W.Dent[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Dent.3d.gpa$consensus, W.Dent[,], points.col="black", lines.col="black", add = T)
plot.coords(Dent.2d3d.pc2max, W.Dent[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Dent.3d.gpa$consensus,adj = c(2,2),text = c("An","A-Q","Art","D"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Dentview2$userMatrix)        

# 3.3) Max4.3D PC axis 1
plot.coords(Max4.2d3d.pc1min, W.Max4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Max4bar.3d.gpa$consensus, W.Max4[,], points.col="black", lines.col="black", add = T)
plot.coords(Max4.2d3d.pc1max, W.Max4[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Max4bar.3d.gpa$consensus,adj = c(2,2),text = c("A-Q","ES","Art","Max"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Max4view$userMatrix)           
# Max4.3D PC axis 2
plot.coords(Max4.2d3d.pc2min, W.Max4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Max4bar.3d.gpa$consensus, W.Max4[,], points.col="black", lines.col="black", add = T)
plot.coords(Max4.2d3d.pc2max, W.Max4[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Max4bar.3d.gpa$consensus,adj = c(2,2),text = c("A-Q","ES","Art","Max"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Max4view$userMatrix)           
# Max4.3D PC axis 1
plot.coords(Max4.2d3d.pc1min, W.Max4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Max4bar.3d.gpa$consensus, W.Max4[,], points.col="black", lines.col="black", add = T)
plot.coords(Max4.2d3d.pc1max, W.Max4[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Max4bar.3d.gpa$consensus,adj = c(2,2),text = c("A-Q","ES","Art","Max"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Max4view2$userMatrix)           
# Max4.3D PC axis 2
plot.coords(Max4.2d3d.pc2min, W.Max4[,], points.col="slategray2", lines.col="slategray2")
plot.coords(Max4bar.3d.gpa$consensus, W.Max4[,], points.col="black", lines.col="black", add = T)
plot.coords(Max4.2d3d.pc2max, W.Max4[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Max4bar.3d.gpa$consensus,adj = c(2,2),text = c("A-Q","ES","Art","Max"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = Max4view2$userMatrix)     

# 3.3) Tot.3D PC axis 1
plot.coords(tot.2d3d.pc1min, W.sym.tot[,], points.col="slategray2", lines.col="slategray2")
plot.coords(spec.sym.3d.gpa$consensus, W.sym.tot[,], points.col="black", lines.col="black", add = T)
plot.coords(tot.2d3d.pc1max, W.sym.tot[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Max4bar.3d.gpa$consensus,adj = c(2,2),text = c("A-Q","ES","Art","Max"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = totview$userMatrix)           
# Tot.3D PC axis 2
plot.coords(tot.2d3d.pc2min, W.sym.tot[,], points.col="slategray2", lines.col="slategray2")
plot.coords(spec.sym.3d.gpa$consensus, W.sym.tot[,], points.col="black", lines.col="black", add = T)
plot.coords(tot.2d3d.pc2max, W.sym.tot[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Max4bar.3d.gpa$consensus,adj = c(2,2),text = c("A-Q","ES","Art","Max"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = totview$userMatrix)           
# Tot.3D PC axis 1
plot.coords(tot.2d3d.pc1min, W.sym.tot[,], points.col="slategray2", lines.col="slategray2")
plot.coords(spec.sym.3d.gpa$consensus, W.sym.tot[,], points.col="black", lines.col="black", add = T)
plot.coords(tot.2d3d.pc1max, W.sym.tot[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Max4bar.3d.gpa$consensus,adj = c(2,2),text = c("A-Q","ES","Art","Max"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = totview2$userMatrix)           
# Tot.3D PC axis 2
plot.coords(tot.2d3d.pc2min, W.sym.tot[,], points.col="slategray2", lines.col="slategray2")
plot.coords(spec.sym.3d.gpa$consensus, W.sym.tot[,], points.col="black", lines.col="black", add = T)
plot.coords(tot.2d3d.pc2max, W.sym.tot[,], points.col="sienna1", lines.col="sienna1", add = T)
#text3d(Max4bar.3d.gpa$consensus,adj = c(2,2),text = c("A-Q","ES","Art","Max"), col = "black", alpha = 1)
rgl.viewpoint(userMatrix = totview2$userMatrix)     
