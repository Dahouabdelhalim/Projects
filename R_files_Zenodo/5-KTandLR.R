#download linkR, svgViewR, and matools from  https://github.com/aaronolsen, 
install.packages("~/My Cloud/Desktop/science!/grad school/McGill/Geometric Morphometrics for Biologists/linkR/linkR_1.2", repos = NULL, type="source")
install.packages("~/My Cloud/Desktop/science!/grad school/McGill/Geometric Morphometrics for Biologists/linkR/matools", repos = NULL, type="source")
install.packages("~/My Cloud/Desktop/science!/grad school/McGill/Geometric Morphometrics for Biologists/linkR/svgViewR_1.3.1", repos = NULL, type="source")
install.packages("GUniFrac")
install.packages("jmv")
install.packages("LearnGeom")
install.packages("maptools")
install.packages("ggforce")
install.packages("multcomp")
library(linkR)
library(matools)
library(svgViewR)
library(LearnGeom)
library(abind)
library(mvnormtest)
library(jmv)
library(vegan)
library(geomorph)
library(ggforce)
library(multcomp)
library(cowplot)
library(ggplot2)
library(car)
library(ggpubr)

# INHERITED FROM "1-SeparateShapeSpaces.R"
# Op4.2d.sym, Op4.3d.sym, Max4.2d.sym, Max4.3d.sym, Mand.2d.sym, Mand.3d.sym

#  KINEMATIC MODELS

anova<-stats::anova
abind<-abind::abind


#   1) function to calculate angle
angle.2d<-function(vertex, pA, pB){
  v1<- c(pA[1]-vertex[1],pA[2]-vertex[2])
  v2<- c(pB[1]-vertex[1],pB[2]-vertex[2])
  
  v1mag<-sqrt(v1[1]^2+v1[2]^2)
  v1norm<-c(v1[1]/v1mag,v1[2]/v1mag)
  v2mag<-sqrt(v2[1]^2+v2[2]^2)
  v2norm<-c(v2[1]/v2mag,v2[2]/v2mag)
  
  dot.prod<-(v1norm[1]*v2norm[1]+v1norm[2]*v2norm[2])
  angle<-acos(dot.prod)*(180/pi)
  angle
}


#
#   2) 2D OP 4-BAR 
#
count<-1
Op4.2d.KT<-vector("double",59) # These lines create vectors which the following loop will vill with KT values
Op4.2d.KT.15<-vector("double",59)
Op4.2d.KT.30<-vector("double",59)
repeat{
  
  ## Define the joint coordinates
  joint.coor<-cbind(Op4.2d.sym[,,count],c(0,0,0,0))
  
  # Define the joint types
  joint.types <- c("S", "R", "S", "R")
  
  # Define joint constraints
  joint.cons <- list(diag(3), c(0,0,1), diag(3), c(0,0,1))
  joint.cons.15<- list(diag(3), c(1,0,3.73206), diag(3), c(0,0,1))
  joint.cons.30<- list(diag(3), c(1,0,1.73205), diag(3), c(0,0,1))
  
  
  # Define two links connected by each joint
  joint.conn <- rbind( c('Op2','Op1'),c('Op1','Op3'), c('Op2','OpM'), c('OpM','Op3'))
  
  # Define ground link
  ground.link <- "Op3"
  
  # Define linkage
  mechanism <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                               joint.cons=joint.cons, body.conn=joint.conn, fixed=ground.link, input.joint= "Op", 
                               print.progress=FALSE)
  mechanism.15 <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                                  joint.cons=joint.cons.15, body.conn=joint.conn, fixed=ground.link, input.joint= "Op", 
                                  print.progress=FALSE)
  mechanism.30 <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                                  joint.cons=joint.cons.30, body.conn=joint.conn, fixed=ground.link, input.joint= "Op", 
                                  print.progress=FALSE)
  
  # Set input parameters
  input.param <- seq(0,-pi/18,length=30)
  
  # These lines create simple rectangles representing each body
  # This is necessary to visualize the bodies being animated
  mechanism <- associatePointShape(mechanism, shape='rect', body='Op1', ends=joint.coor[c('vOp', 'Op'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='Op3', ends=joint.coor[c('Op', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='Op2', ends=joint.coor[c('vOp', '12'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='OpM', ends=joint.coor[c('12', '11'), ], nvector=c(0,0,1), width=0.01)
  
  mechanism.15 <- associatePointShape(mechanism.15, shape='rect', body='Op1', ends=joint.coor[c('vOp', 'Op'), ], nvector=c(0,0,1), width=0.01)
  mechanism.15 <- associatePointShape(mechanism.15, shape='rect', body='Op3', ends=joint.coor[c('Op', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism.15 <- associatePointShape(mechanism.15, shape='rect', body='Op2', ends=joint.coor[c('vOp', '12'), ], nvector=c(0,0,1), width=0.01)
  mechanism.15 <- associatePointShape(mechanism.15, shape='rect', body='OpM', ends=joint.coor[c('12', '11'), ], nvector=c(0,0,1), width=0.01)
  
  mechanism.30 <- associatePointShape(mechanism.30, shape='rect', body='Op1', ends=joint.coor[c('vOp', 'Op'), ], nvector=c(0,0,1), width=0.01)
  mechanism.30 <- associatePointShape(mechanism.30, shape='rect', body='Op3', ends=joint.coor[c('Op', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism.30 <- associatePointShape(mechanism.30, shape='rect', body='Op2', ends=joint.coor[c('vOp', '12'), ], nvector=c(0,0,1), width=0.01)
  mechanism.30 <- associatePointShape(mechanism.30, shape='rect', body='OpM', ends=joint.coor[c('12', '11'), ], nvector=c(0,0,1), width=0.01)
  
  
  # Animate mechanism
  Op4.2d.animated <- animateMechanism(mechanism, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  Op4.2d.animated.15 <- animateMechanism(mechanism.15, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  Op4.2d.animated.30 <- animateMechanism(mechanism.30, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  
  # Draw mechanism
  Op4.2d.animate_mechanism <- drawMechanism(Op4.2d.animated, file='Op4.2D.html')
  Op4.2d.animate_mechanism.15 <- drawMechanism(Op4.2d.animated.15, file='Op4.2D.15.html')
  Op4.2d.animate_mechanism.30 <- drawMechanism(Op4.2d.animated.30, file='Op4.2D.30.html')
  
  theta.out<-angle.2d(Op4.2d.animated$joint.coor.anim["11",,1,1],Op4.2d.animated$joint.coor.anim["12",,1,1],
                      Op4.2d.animated$joint.coor.anim["12",,30,1] )
  theta.out.15<-angle.2d(Op4.2d.animated.15$joint.coor.anim["11",,1,1],Op4.2d.animated.15$joint.coor.anim["12",,1,1],
                         Op4.2d.animated.15$joint.coor.anim["12",,30,1] )
  theta.out.30<-angle.2d(Op4.2d.animated.30$joint.coor.anim["11",,1,1],Op4.2d.animated.30$joint.coor.anim["12",,1,1],
                         Op4.2d.animated.30$joint.coor.anim["12",,30,1] )
  Op4.2d.KT[count]<-theta.out/(pi/18*(180/pi))
  Op4.2d.KT.15[count]<-theta.out.15/(pi/18*(180/pi))
  Op4.2d.KT.30[count]<-theta.out.30/(pi/18*(180/pi))
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}


#
#   3) 3D OP 4-BAR
#
count<-1
Op4.3d.KT<-vector("double",59)
Op4.3d.KT.15<-vector("double",59)
Op4.3d.KT.30<-vector("double",59)
repeat{
  
  ## Define the joint coordinates
  joint.coor<-Op4.3d.sym[,,count]
  
  # Define the joint types
  joint.types <- c("S", "R", "S", "R")
  
  # Define joint constraints
  joint.cons <- list(diag(3), c(0,0,1), diag(3), c(0,0,1))
  joint.cons.15<- list(diag(3), c(1,0,3.73206), diag(3), c(0,0,1))
  joint.cons.30<- list(diag(3), c(1,0,1.73205), diag(3), c(0,0,1))
  
  
  # Define two links connected by each joint
  joint.conn <- rbind( c('Op2','Op1'),c('Op1','Op3'), c('Op2','OpM'), c('OpM','Op3'))
  
  # Define ground link
  ground.link <- "Op3"
  
  # Define linkage
  mechanism <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                               joint.cons=joint.cons, body.conn=joint.conn, fixed=ground.link, input.joint= "Op", 
                               print.progress=FALSE)
  mechanism.15 <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                                  joint.cons=joint.cons.15, body.conn=joint.conn, fixed=ground.link, input.joint= "Op", 
                                  print.progress=FALSE)
  mechanism.30 <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                                  joint.cons=joint.cons.30, body.conn=joint.conn, fixed=ground.link, input.joint= "Op", 
                                  print.progress=FALSE)
  
  # Set input parameters
  input.param <- seq(0,-pi/18,length=30)
  
  # These lines create simple rectangles representing each body
  # This is necessary to visualize the bodies being animated
  mechanism <- associatePointShape(mechanism, shape='rect', body='Op1', ends=joint.coor[c('vOp', 'Op'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='Op3', ends=joint.coor[c('Op', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='Op2', ends=joint.coor[c('vOp', '12'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='OpM', ends=joint.coor[c('12', '11'), ], nvector=c(0,0,1), width=0.01)
  
  mechanism.15 <- associatePointShape(mechanism.15, shape='rect', body='Op1', ends=joint.coor[c('vOp', 'Op'), ], nvector=c(0,0,1), width=0.01)
  mechanism.15 <- associatePointShape(mechanism.15, shape='rect', body='Op3', ends=joint.coor[c('Op', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism.15 <- associatePointShape(mechanism.15, shape='rect', body='Op2', ends=joint.coor[c('vOp', '12'), ], nvector=c(0,0,1), width=0.01)
  mechanism.15 <- associatePointShape(mechanism.15, shape='rect', body='OpM', ends=joint.coor[c('12', '11'), ], nvector=c(0,0,1), width=0.01)
  
  mechanism.30 <- associatePointShape(mechanism.30, shape='rect', body='Op1', ends=joint.coor[c('vOp', 'Op'), ], nvector=c(0,0,1), width=0.01)
  mechanism.30 <- associatePointShape(mechanism.30, shape='rect', body='Op3', ends=joint.coor[c('Op', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism.30 <- associatePointShape(mechanism.30, shape='rect', body='Op2', ends=joint.coor[c('vOp', '12'), ], nvector=c(0,0,1), width=0.01)
  mechanism.30 <- associatePointShape(mechanism.30, shape='rect', body='OpM', ends=joint.coor[c('12', '11'), ], nvector=c(0,0,1), width=0.01)
  
  # Animate mechanism
  Op4.3d.animated <- animateMechanism(mechanism, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  Op4.3d.animated.15 <- animateMechanism(mechanism.15, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  Op4.3d.animated.30 <- animateMechanism(mechanism.30, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  
  # Draw mechanism
  Op4.3d.animate_mechanism <- drawMechanism(Op4.3d.animated, file='Op4.3D.html')
  Op4.3d.animate_mechanism.15 <- drawMechanism(Op4.3d.animated.15, file='Op4.3D.15.html')
  Op4.3d.animate_mechanism.30 <- drawMechanism(Op4.3d.animated.30, file='Op4.3D.30.html')
  
  
  theta.out<-angle.2d(Op4.3d.animated$joint.coor.anim["11",,1,1],Op4.3d.animated$joint.coor.anim["12",,1,1],
                      Op4.3d.animated$joint.coor.anim["12",,30,1] )
  theta.out.15<-angle.2d(Op4.3d.animated.15$joint.coor.anim["11",,1,1],Op4.3d.animated.15$joint.coor.anim["12",,1,1],
                         Op4.3d.animated.15$joint.coor.anim["12",,30,1] )
  theta.out.30<-angle.2d(Op4.3d.animated.30$joint.coor.anim["11",,1,1],Op4.3d.animated.30$joint.coor.anim["12",,1,1],
                         Op4.3d.animated.30$joint.coor.anim["12",,30,1] )
  Op4.3d.KT[count]<-theta.out/(pi/18*(180/pi))
  Op4.3d.KT.15[count]<-theta.out.15/(pi/18*(180/pi))
  Op4.3d.KT.30[count]<-theta.out.30/(pi/18*(180/pi))
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}




#
#   4) 2D MAX 4-BAR
#
Max4.2d.sym
count<-1
Max4.2d.KT<-vector("double",59)
Max4.2d.KT.10<-vector("double",59)
Max4.2d.KT.20<-vector("double",59)
repeat{
  
  ## Define the joint coordinates
  joint.coor<-cbind(Max4.2d.sym[,,count],c(0,0,0,0))
  
  # Define the joint types
  joint.types <- c("R", "S", "S", "R")
  
  # Define joint constraints
  joint.cons <- list( c(0,0,1), diag(3), diag(3), c(0,0,1))
  joint.cons.10 <- list( c(0,0,1), diag(3), diag(3), c(0,1,5.67127))
  joint.cons.20 <- list( c(0,0,1), diag(3), diag(3), c(0,2,5.49496))
  
  
  # Define two links connected by each joint
  joint.conn <- rbind( c('Max3','MaxM'), c('MaxM','Max1'),c('Max2','Max1'), c('Max2','Max3'))
  
  # Define ground link
  ground.link <- "Max3"
  
  # Define linkage
  mechanism <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                               joint.cons=joint.cons, body.conn=joint.conn, fixed=ground.link, input.joint= "11", 
                               print.progress=FALSE)
  mechanism.10 <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                                  joint.cons=joint.cons.10, body.conn=joint.conn, fixed=ground.link, input.joint= "11", 
                                  print.progress=FALSE)
  mechanism.20 <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                                  joint.cons=joint.cons.20, body.conn=joint.conn, fixed=ground.link, input.joint= "11", 
                                  print.progress=FALSE)
  
  # Set input parameters
  input.param <- seq(0,-pi/18,length=30)
  
  # These lines create simple rectangles representing each body
  # This is necessary to visualize the bodies being animated
  mechanism <- associatePointShape(mechanism, shape='rect', body='Max1', ends=joint.coor[c('6', '3'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='Max2', ends=joint.coor[c('3', '10'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='Max3', ends=joint.coor[c('10', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='MaxM', ends=joint.coor[c('6', '11'), ], nvector=c(0,0,1), width=0.01)
  
  mechanism.10 <- associatePointShape(mechanism.10, shape='rect', body='Max1', ends=joint.coor[c('6', '3'), ], nvector=c(0,0,1), width=0.01)
  mechanism.10 <- associatePointShape(mechanism.10, shape='rect', body='Max2', ends=joint.coor[c('3', '10'), ], nvector=c(0,0,1), width=0.01)
  mechanism.10 <- associatePointShape(mechanism.10, shape='rect', body='Max3', ends=joint.coor[c('10', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism.10 <- associatePointShape(mechanism.10, shape='rect', body='MaxM', ends=joint.coor[c('6', '11'), ], nvector=c(0,0,1), width=0.01)
  
  mechanism.20 <- associatePointShape(mechanism.20, shape='rect', body='Max1', ends=joint.coor[c('6', '3'), ], nvector=c(0,0,1), width=0.01)
  mechanism.20 <- associatePointShape(mechanism.20, shape='rect', body='Max2', ends=joint.coor[c('3', '10'), ], nvector=c(0,0,1), width=0.01)
  mechanism.20 <- associatePointShape(mechanism.20, shape='rect', body='Max3', ends=joint.coor[c('10', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism.20 <- associatePointShape(mechanism.20, shape='rect', body='MaxM', ends=joint.coor[c('6', '11'), ], nvector=c(0,0,1), width=0.01)
  
  
  # Animate mechanism
  Max4.2d.animated <- animateMechanism(mechanism, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  Max4.2d.animated.10<- animateMechanism(mechanism.10, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  Max4.2d.animated.20<- animateMechanism(mechanism.20, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  
  # Draw mechanism
  Max4.2d.animate_mechanism <- drawMechanism(Max4.2d.animated, file='Max4.2D.html')
  Max4.2d.animate_mechanism.10 <- drawMechanism(Max4.2d.animated.10, file='Max4.2D.10.html')
  Max4.2d.animate_mechanism.20 <- drawMechanism(Max4.2d.animated.20, file='Max4.2D.20.html')
  
  theta.angle.10<-rbind(as.vector((rotationMatrixZYX(0,0,pi/18)%*%Max4.2d.animated.10$joint.coor.anim["10",,1,1])),
                        as.vector((rotationMatrixZYX(0,0,pi/18)%*%Max4.2d.animated.10$joint.coor.anim["3",,1,1])),
                        as.vector((rotationMatrixZYX(0,0,pi/18)%*%Max4.2d.animated.10$joint.coor.anim["3",,30,1])))
  theta.angle.20<-rbind(as.vector((rotationMatrixZYX(0,0,pi/9)%*%Max4.2d.animated.20$joint.coor.anim["10",,1,1])),
                        as.vector((rotationMatrixZYX(0,0,pi/9)%*%Max4.2d.animated.20$joint.coor.anim["3",,1,1])),
                        as.vector((rotationMatrixZYX(0,0,pi/9)%*%Max4.2d.animated.20$joint.coor.anim["3",,30,1])))
  
  theta.out<-angle.2d(Max4.2d.animated$joint.coor.anim["10",,1,1],Max4.2d.animated$joint.coor.anim["3",,1,1],
                      Max4.2d.animated$joint.coor.anim["3",,30,1] )
  theta.out.10<-angle.2d(theta.angle.10[1,],theta.angle.10[2,],theta.angle.10[3,])
  theta.out.20<-angle.2d(theta.angle.20[1,],theta.angle.20[2,],theta.angle.20[3,])
  
  Max4.2d.KT[count]<-theta.out/(pi/18*(180/pi))  
  Max4.2d.KT.10[count]<-theta.out.10/(pi/18*(180/pi)) 
  Max4.2d.KT.20[count]<-theta.out.20/(pi/18*(180/pi))  
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}


#
#   5) 3D MAX 4-BAR
#
Max4.3d.sym
count<-1
Max4.3d.KT<-vector("double",59)
Max4.3d.KT.10<-vector("double",59)
Max4.3d.KT.20<-vector("double",59)
repeat{
  
  ## Define the joint coordinates
  joint.coor<-Max4.3d.sym[,,count]
  
  # Define the joint types
  joint.types <- c("R", "S", "S", "R")
  
  # Define joint constraints
  joint.cons <- list( c(0,0,1), diag(3), diag(3), c(0,0,1))
  joint.cons.10 <- list( c(0,0,1), diag(3), diag(3), c(0,1,5.67127))
  joint.cons.20 <- list( c(0,0,1), diag(3), diag(3), c(0,2,5.49496))
  
  # Define two links connected by each joint
  joint.conn <- rbind( c('Max3','MaxM'), c('MaxM','Max1'),c('Max2','Max1'), c('Max2','Max3'))
  
  # Define ground link
  ground.link <- "Max3"
  
  # Define linkage
  mechanism <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                               joint.cons=joint.cons, body.conn=joint.conn, fixed=ground.link, input.joint= "11", 
                               print.progress=FALSE)
  mechanism.10 <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                                  joint.cons=joint.cons.10, body.conn=joint.conn, fixed=ground.link, input.joint= "11", 
                                  print.progress=FALSE)
  mechanism.20 <- defineMechanism(joint.coor=joint.coor, joint.types=joint.types, 
                                  joint.cons=joint.cons.20, body.conn=joint.conn, fixed=ground.link, input.joint= "11", 
                                  print.progress=FALSE)
  
  # Set input parameters
  input.param <- seq(0,-pi/18,length=30)
  
  # These lines create simple rectangles representing each body
  # This is necessary to visualize the bodies being animated
  mechanism <- associatePointShape(mechanism, shape='rect', body='Max1', ends=joint.coor[c('6', '3'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='Max2', ends=joint.coor[c('3', '10'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='Max3', ends=joint.coor[c('10', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism <- associatePointShape(mechanism, shape='rect', body='MaxM', ends=joint.coor[c('6', '11'), ], nvector=c(0,0,1), width=0.01)
  
  mechanism.10 <- associatePointShape(mechanism.10, shape='rect', body='Max1', ends=joint.coor[c('6', '3'), ], nvector=c(0,0,1), width=0.01)
  mechanism.10 <- associatePointShape(mechanism.10, shape='rect', body='Max2', ends=joint.coor[c('3', '10'), ], nvector=c(0,0,1), width=0.01)
  mechanism.10 <- associatePointShape(mechanism.10, shape='rect', body='Max3', ends=joint.coor[c('10', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism.10 <- associatePointShape(mechanism.10, shape='rect', body='MaxM', ends=joint.coor[c('6', '11'), ], nvector=c(0,0,1), width=0.01)
  
  mechanism.20 <- associatePointShape(mechanism.20, shape='rect', body='Max1', ends=joint.coor[c('6', '3'), ], nvector=c(0,0,1), width=0.01)
  mechanism.20 <- associatePointShape(mechanism.20, shape='rect', body='Max2', ends=joint.coor[c('3', '10'), ], nvector=c(0,0,1), width=0.01)
  mechanism.20 <- associatePointShape(mechanism.20, shape='rect', body='Max3', ends=joint.coor[c('10', '11'), ], nvector=c(0,0,1), width=0.01)
  mechanism.20 <- associatePointShape(mechanism.20, shape='rect', body='MaxM', ends=joint.coor[c('6', '11'), ], nvector=c(0,0,1), width=0.01)
  
  # Animate mechanism
  Max4.3d.animated<- animateMechanism(mechanism, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  Max4.3d.animated.10<- animateMechanism(mechanism.10, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  Max4.3d.animated.20<- animateMechanism(mechanism.20, input.param=input.param, print.progress=FALSE, check.inter.joint.dist=FALSE)
  
  # Draw mechanism
  Max4.3d.animate_mechanism <- drawMechanism(Max4.3d.animated, file='Max4.3D.html')
  Max4.3d.animate_mechanism.10 <- drawMechanism(Max4.3d.animated.10, file='Max4.3D.10.html')
  Max4.3d.animate_mechanism.20 <- drawMechanism(Max4.3d.animated.20, file='Max4.3D.20.html')
  
  theta.angle.10<-rbind(as.vector((rotationMatrixZYX(0,0,pi/18)%*%Max4.3d.animated.10$joint.coor.anim["10",,1,1])),
                        as.vector((rotationMatrixZYX(0,0,pi/18)%*%Max4.3d.animated.10$joint.coor.anim["3",,1,1])),
                        as.vector((rotationMatrixZYX(0,0,pi/18)%*%Max4.3d.animated.10$joint.coor.anim["3",,30,1])))
  theta.angle.20<-rbind(as.vector((rotationMatrixZYX(0,0,pi/9)%*%Max4.3d.animated.20$joint.coor.anim["10",,1,1])),
                        as.vector((rotationMatrixZYX(0,0,pi/9)%*%Max4.3d.animated.20$joint.coor.anim["3",,1,1])),
                        as.vector((rotationMatrixZYX(0,0,pi/9)%*%Max4.3d.animated.20$joint.coor.anim["3",,30,1])))
  
  theta.out<-angle.2d(Max4.3d.animated.20$joint.coor.anim["10",,1,1],Max4.3d.animated.20$joint.coor.anim["3",,1,1],
                      Max4.3d.animated.20$joint.coor.anim["3",,30,1])
  theta.out.10<-angle.2d(theta.angle.10[1,],theta.angle.10[2,],theta.angle.10[3,])
  theta.out.20<-angle.2d(theta.angle.20[1,],theta.angle.20[2,],theta.angle.20[3,])
  
  Max4.3d.KT[count]<-theta.out/(pi/18*(180/pi))  
  Max4.3d.KT.10[count]<-theta.out.10/(pi/18*(180/pi)) 
  Max4.3d.KT.20[count]<-theta.out.20/(pi/18*(180/pi)) 
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}



#
#   6) 2D MAND
#
Mand.2d.sym
count<-1
Mand.2d.openingLR<-vector("double",59)
Mand.2d.closingLR<-vector("double",59)
repeat{
  
  OpM.length<-sqrt((Mand.2d.sym["11",1,count]-Mand.2d.sym["12",1,count])^2+
                     (Mand.2d.sym["11",2,count]-Mand.2d.sym["12",2,count])^2)
  Mand.length<-sqrt((Mand.2d.sym["11",1,count]-Mand.2d.sym["2",1,count])^2+
                      (Mand.2d.sym["11",2,count]-Mand.2d.sym["2",2,count])^2)
  MaxM.length<-sqrt((Mand.2d.sym["11",1,count]-Mand.2d.sym["6",1,count])^2+
                      (Mand.2d.sym["11",2,count]-Mand.2d.sym["6",2,count])^2)
  
  Mand.2d.openingLR[count]<-(OpeningLR<-(Mand.length/OpM.length))
  Mand.2d.closingLR[count]<-(ClosingLR<-(Mand.length/MaxM.length))
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}


#
#   7) 3D Mand
#
Mand.3d.sym
count<-1
Mand.3d.openingLR<-vector("double",59)
Mand.3d.closingLR<-vector("double",59)
Mand.3d.RotAxis<-matrix(nrow = 3, ncol = 3)
repeat{
  
  Mand.3d.RotAxis<-matrix(rep(Mand.3d.sym["11",,count], times=3),nrow = 3, ncol = 3, byrow = T,dimnames = list(c("OpM","MaxM","Mand")))
  Mand.3d.RotAxis["OpM",3]<-Mand.3d.sym["12",3,count]
  Mand.3d.RotAxis["MaxM",3]<-Mand.3d.sym["6",3,count]
  Mand.3d.RotAxis["Mand",3]<-Mand.3d.sym["2",3,count]
  
  
  OpM.length<-sqrt((Mand.3d.RotAxis["OpM",1]-Mand.3d.sym["12",1,count])^2+
                     (Mand.3d.RotAxis["OpM",2]-Mand.3d.sym["12",2,count])^2+
                     (Mand.3d.RotAxis["OpM",3]-Mand.3d.sym["12",3,count])^2)
  Mand.length<-sqrt((Mand.3d.RotAxis["Mand",1]-Mand.3d.sym["2",1,count])^2+
                      (Mand.3d.RotAxis["Mand",2]-Mand.3d.sym["2",2,count])^2+
                      (Mand.3d.RotAxis["Mand",3]-Mand.3d.sym["2",3,count])^2)
  MaxM.length<-sqrt((Mand.3d.RotAxis["MaxM",1]-Mand.3d.sym["6",1,count])^2+
                      (Mand.3d.RotAxis["MaxM",2]-Mand.3d.sym["6",2,count])^2+
                      (Mand.3d.RotAxis["MaxM",3]-Mand.3d.sym["6",3,count])^2)
  
  Mand.3d.openingLR[count]<-(OpeningLR<-(Mand.length/OpM.length))
  Mand.3d.closingLR[count]<-(ClosingLR<-(Mand.length/MaxM.length))
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}

#  8) creates df of KT and LR caluclations
kinem.df<-cbind(spec.sym.2d.gpa$Csize,spec.sym.3d.gpa$Csize,Op4.2d.KT,Op4.2d.KT.15,Op4.2d.KT.30,Op4.3d.KT,Op4.3d.KT.15,Op4.3d.KT.30,
                Max4.2d.KT,Max4.2d.KT.10,Max4.2d.KT.20,Max4.3d.KT,Max4.3d.KT.10,Max4.3d.KT.20,
                Mand.2d.openingLR, Mand.2d.closingLR, Mand.3d.openingLR, Mand.3d.closingLR)
#  8.1) calculates centroid size residuals of KT and LR 
kinem.df.res<-cbind(Op4.2d.KT.res=lm(log10(spec.sym.2d.gpa$Csize) ~ Op4.2d.KT)$residuals,
                    Op4.2d.KT.15.res=lm(log10(spec.sym.2d.gpa$Csize) ~ Op4.2d.KT.15)$residuals,
                    Op4.2d.KT.30.res=lm(log10(spec.sym.2d.gpa$Csize) ~ Op4.2d.KT.30)$residuals,
                    Op4.3d.KT.res=lm(log10(spec.sym.3d.gpa$Csize) ~ Op4.3d.KT)$residuals,
                    Op4.3d.KT.15.res=lm(log10(spec.sym.3d.gpa$Csize) ~ Op4.3d.KT.15)$residuals,
                    Op4.3d.KT.30.res=lm(log10(spec.sym.3d.gpa$Csize) ~ Op4.3d.KT.30)$residuals,
                    
                    Max4.2d.KT.res=lm(log10(spec.sym.2d.gpa$Csize) ~ Max4.2d.KT)$residuals,
                    Max4.2d.KT.10.res=lm(log10(spec.sym.2d.gpa$Csize) ~ Max4.2d.KT.10)$residuals,
                    Max4.2d.KT.20.res=lm(log10(spec.sym.2d.gpa$Csize) ~ Max4.2d.KT.20)$residuals,
                    Max4.3d.KT.res=lm(log10(spec.sym.3d.gpa$Csize) ~ Max4.3d.KT)$residuals,
                    Max4.3d.KT.10.res=lm(log10(spec.sym.3d.gpa$Csize) ~ Max4.3d.KT.10)$residuals,
                    Max4.3d.KT.20.res=lm(log10(spec.sym.3d.gpa$Csize) ~ Max4.3d.KT.20)$residuals,
                    
                    Mand.2d.openingLR.res=lm(log10(spec.sym.2d.gpa$Csize) ~ Mand.2d.openingLR)$residuals,
                    Mand.2d.closingLR.res=lm(log10(spec.sym.2d.gpa$Csize) ~ Mand.2d.closingLR)$residuals,
                    Mand.3d.openingLR.res=lm(log10(spec.sym.3d.gpa$Csize) ~ Mand.3d.openingLR)$residuals,
                    Mand.3d.closingLR.res=lm(log10(spec.sym.3d.gpa$Csize) ~ Mand.3d.closingLR)$residuals,
                    ldmk.factors)
kinem.df.res.rrpp<-rrpp.data.frame(kinem.df.res)

#  9) Op4 Watershed*dietLD linear models
Op4.2d.KT.resxLD.lm<-(lm(Op4.2d.KT.res~Watershed*LD1, data = kinem.df.res))
Op4.2d.KT.15.resxLD.lm<-(lm(Op4.2d.KT.15.res~Watershed*LD1, data = kinem.df.res))
Op4.2d.KT.30.resxLD.lm<-(lm(Op4.2d.KT.30.res~Watershed*LD1, data = kinem.df.res))
Op4.3d.KT.resxLD.lm<-(lm(Op4.3d.KT.res~Watershed*LD1, data = kinem.df.res))
Op4.3d.KT.15.resxLD.lm<-(lm(Op4.3d.KT.15.res~Watershed*LD1, data = kinem.df.res))
Op4.3d.KT.30.resxLD.lm<-(lm(Op4.3d.KT.30.res~Watershed*LD1, data = kinem.df.res))
Op4.2d.KT.resxLD.lm
Anova(Op4.2d.KT.resxLD.lm, type = "II")
Op4.2d.KT.15.resxLD.lm
Anova(Op4.2d.KT.15.resxLD.lm, type = "II")
Op4.2d.KT.30.resxLD.lm
Anova(Op4.2d.KT.30.resxLD.lm, type = "II")
Op4.3d.KT.resxLD.lm
Anova(Op4.3d.KT.resxLD.lm, type = "II")
Op4.3d.KT.15.resxLD.lm
Anova(Op4.3d.KT.15.resxLD.lm, type = "II")
Op4.3d.KT.30.resxLD.lm
Anova(Op4.3d.KT.30.resxLD.lm, type = "II")
# 9.1) Op4 Watershed*simpson linear models
Op4.2d.KT.resxSimp.lm<-lm(Op4.2d.KT.res~Watershed*simpson, data = kinem.df.res)
Op4.2d.KT.15.resxSimp.lm<-lm(Op4.2d.KT.15.res~Watershed*simpson, data = kinem.df.res)
Op4.2d.KT.30.resxSimp.lm<-lm(Op4.2d.KT.30.res~Watershed*simpson, data = kinem.df.res)
Op4.3d.KT.resxSimp.lm<-lm(Op4.3d.KT.res~Watershed*simpson, data = kinem.df.res)
Op4.3d.KT.15.resxSimp.lm<-lm(Op4.3d.KT.15.res~Watershed*simpson, data = kinem.df.res)
Op4.3d.KT.30.resxSimp.lm<-lm(Op4.3d.KT.30.res~Watershed*simpson, data = kinem.df.res)
Op4.2d.KT.resxSimp.lm
Anova(Op4.2d.KT.resxSimp.lm, type = "II")
Op4.2d.KT.15.resxSimp.lm
Anova(Op4.2d.KT.15.resxSimp.lm, type = "II")
Op4.2d.KT.30.resxSimp.lm
Anova(Op4.2d.KT.30.resxSimp.lm, type = "II")
Op4.3d.KT.resxSimp.lm
Anova(Op4.3d.KT.resxSimp.lm, type = "II")
Op4.3d.KT.15.resxSimp.lm
Anova(Op4.3d.KT.15.resxSimp.lm, type = "II")
Op4.3d.KT.30.resxSimp.lm
Anova(Op4.3d.KT.30.resxSimp.lm, type = "II")
# 9.2) Op4 Watershed*habitat linear models
Op4.2d.KT.resxHabitat.lm<-(lm(Op4.2d.KT.res~Watershed*Habitat, data = kinem.df.res))
Op4.2d.KT.15.resxHabitat.lm<-(lm(Op4.2d.KT.15.res~Watershed*Habitat, data = kinem.df.res))
Op4.2d.KT.30.resxHabitat.lm<-(lm(Op4.2d.KT.30.res~Watershed*Habitat, data = kinem.df.res))
Op4.3d.KT.resxHabitat.lm<-(lm(Op4.3d.KT.res~Watershed*Habitat, data = kinem.df.res))
Op4.3d.KT.15.resxHabitat.lm<-(lm(Op4.3d.KT.15.res~Watershed*Habitat, data = kinem.df.res))
Op4.3d.KT.30.resxHabitat.lm<-(lm(Op4.3d.KT.30.res~Watershed*Habitat, data = kinem.df.res))
Op4.2d.KT.resxHabitat.lm
Anova(Op4.2d.KT.resxHabitat.lm, type = "II")
Op4.2d.KT.15.resxHabitat.lm
Anova(Op4.2d.KT.15.resxHabitat.lm, type = "II")
Op4.2d.KT.30.resxHabitat.lm
Anova(Op4.2d.KT.30.resxHabitat.lm, type = "II")
Op4.3d.KT.resxHabitat.lm
Anova(Op4.3d.KT.resxHabitat.lm, type = "II")
Op4.3d.KT.15.resxHabitat.lm
Anova(Op4.3d.KT.15.resxHabitat.lm, type = "II")
Op4.3d.KT.30.resxHabitat.lm
Anova(Op4.3d.KT.30.resxHabitat.lm, type = "II")

#  10) Max4KT Watershed*dietLD linear models
Max4.2d.KT.resxLD.lm<-(lm(Max4.2d.KT.res~Watershed*LD1, data = kinem.df.res))
Max4.2d.KT.10.resxLD.lm<-(lm(Max4.2d.KT.10.res~Watershed*LD1, data = kinem.df.res))
Max4.2d.KT.20.resxLD.lm<-(lm(Max4.2d.KT.20.res~Watershed*LD1, data = kinem.df.res))
Max4.3d.KT.resxLD.lm<-(lm(Max4.3d.KT.res~Watershed*LD1, data = kinem.df.res))
Max4.3d.KT.10.resxLD.lm<-(lm(Max4.3d.KT.10.res~Watershed*LD1, data = kinem.df.res))
Max4.3d.KT.20.resxLD.lm<-(lm(Max4.3d.KT.20.res~Watershed*LD1, data = kinem.df.res))
Max4.2d.KT.resxLD.lm
Anova(Max4.2d.KT.resxLD.lm, type = "II")
Max4.2d.KT.10.resxLD.lm
Anova(Max4.2d.KT.10.resxLD.lm, type = "II")
Max4.2d.KT.20.resxLD.lm
Anova(Max4.2d.KT.20.resxLD.lm, type = "II")
Max4.3d.KT.resxLD.lm
Anova(Max4.3d.KT.resxLD.lm, type = "II")
Max4.3d.KT.10.resxLD.lm
Anova(Max4.3d.KT.10.resxLD.lm, type = "II")
Max4.3d.KT.20.resxLD.lm
Anova(Max4.3d.KT.20.resxLD.lm, type = "II")
#  10.1) Max4KT Watershed*simpson linear models
Max4.2d.KT.resxSimp.lm<-(lm(Max4.2d.KT.res~Watershed*simpson, data = kinem.df.res))
Max4.2d.KT.10.resxSimp.lm<-(lm(Max4.2d.KT.10.res~Watershed*simpson, data = kinem.df.res))
Max4.2d.KT.20.resxSimp.lm<-(lm(Max4.2d.KT.20.res~Watershed*simpson, data = kinem.df.res))
Max4.3d.KT.resxSimp.lm<-(lm(Max4.3d.KT.res~Watershed*simpson, data = kinem.df.res))
Max4.3d.KT.10.resxSimp.lm<-(lm(Max4.3d.KT.10.res~Watershed*simpson, data = kinem.df.res))
Max4.3d.KT.20.resxSimp.lm<-(lm(Max4.3d.KT.20.res~Watershed*simpson, data = kinem.df.res))
Max4.2d.KT.resxSimp.lm
Anova(Max4.2d.KT.resxSimp.lm, type = "II")
Max4.2d.KT.10.resxSimp.lm
Anova(Max4.2d.KT.10.resxSimp.lm, type = "II")
Max4.2d.KT.20.resxSimp.lm
Anova(Max4.2d.KT.20.resxSimp.lm, type = "II")
Max4.3d.KT.resxSimp.lm
Anova(Max4.3d.KT.resxSimp.lm, type = "II")
Max4.3d.KT.10.resxSimp.lm
Anova(Max4.3d.KT.10.resxSimp.lm, type = "II")
Max4.3d.KT.20.resxSimp.lm
Anova(Max4.3d.KT.20.resxSimp.lm, type = "II")

#  10.2) Max4KT Watershed*Habitat linear models
Max4.2d.KT.resxHabitat.lm<-(lm(Max4.2d.KT.res~Watershed*Habitat, data = kinem.df.res))
Max4.2d.KT.10.resxHabitat.lm<-(lm(Max4.2d.KT.10.res~Watershed*Habitat, data = kinem.df.res))
Max4.2d.KT.20.resxHabitat.lm<-(lm(Max4.2d.KT.20.res~Watershed*Habitat, data = kinem.df.res))
Max4.3d.KT.resxHabitat.lm<-(lm(Max4.3d.KT.res~Watershed*Habitat, data = kinem.df.res))
Max4.3d.KT.10.resxHabitat.lm<-(lm(Max4.3d.KT.10.res~Watershed*Habitat, data = kinem.df.res))
Max4.3d.KT.20.resxHabitat.lm<-(lm(Max4.3d.KT.20.res~Watershed*Habitat, data = kinem.df.res))
Max4.2d.KT.resxHabitat.lm
Anova(Max4.2d.KT.resxHabitat.lm, type = "II")
Max4.2d.KT.10.resxHabitat.lm
Anova(Max4.2d.KT.10.resxHabitat.lm, type = "II")
Max4.2d.KT.20.resxHabitat.lm
Anova(Max4.2d.KT.20.resxHabitat.lm, type = "II")
Max4.3d.KT.resxHabitat.lm
Anova(Max4.3d.KT.resxHabitat.lm, type = "II")
Max4.3d.KT.10.resxHabitat.lm
Anova(Max4.3d.KT.10.resxHabitat.lm, type = "II")
Max4.3d.KT.20.resxHabitat.lm
Anova(Max4.3d.KT.20.resxHabitat.lm, type = "II")

#  11) Mandibular levers Watershed*dietLD linear models
Mand.2d.openingLR.resxLD.lm<-(lm(Mand.2d.openingLR.res~Watershed*LD1, data = kinem.df.res))
Mand.2d.closingLR.resxLD.lm<-(lm(Mand.2d.closingLR.res~Watershed*LD1, data = kinem.df.res))
Mand.3d.openingLR.resxLD.lm<-(lm(Mand.3d.openingLR.res~Watershed*LD1, data = kinem.df.res))
Mand.3d.closingLR.resxLD.lm<-(lm(Mand.3d.closingLR.res~Watershed*LD1, data = kinem.df.res))
Mand.2d.openingLR.resxLD.lm
Anova(Mand.2d.openingLR.resxLD.lm, type = "II")
Mand.2d.closingLR.resxLD.lm
Anova(Mand.2d.closingLR.resxLD.lm, type = "II")
Mand.3d.openingLR.resxLD.lm
Anova(Mand.3d.openingLR.resxLD.lm, type = "II")
Mand.3d.closingLR.resxLD.lm
Anova(Mand.3d.closingLR.resxLD.lm, type = "II")

#  11.1) Mandibular levers Watershed*simpson linear models
Mand.2d.openingLR.resxSimp.lm<-(lm(Mand.2d.openingLR.res~Watershed*simpson, data = kinem.df.res))
Mand.2d.closingLR.resxSimp.lm<-(lm(Mand.2d.closingLR.res~Watershed*simpson, data = kinem.df.res))
Mand.3d.openingLR.resxSimp.lm<-(lm(Mand.3d.openingLR.res~Watershed*simpson, data = kinem.df.res))
Mand.3d.closingLR.resxSimp.lm<-(lm(Mand.3d.closingLR.res~Watershed*simpson, data = kinem.df.res))
Mand.2d.openingLR.resxSimp.lm
Anova(Mand.2d.openingLR.resxSimp.lm, type = "II")
Mand.2d.closingLR.resxSimp.lm
Anova(Mand.2d.closingLR.resxSimp.lm, type = "II")
Mand.3d.openingLR.resxSimp.lm
Anova(Mand.3d.openingLR.resxSimp.lm, type = "II")
Mand.3d.closingLR.resxSimp.lm
Anova(Mand.3d.closingLR.resxSimp.lm, type = "II")

#  11.2) Mandibular levers Watershed*Habitat linear models
Mand.2d.openingLR.resxHabitat.lm<-(lm(Mand.2d.openingLR.res~Watershed*Habitat, data = kinem.df.res))
Mand.2d.closingLR.resxHabitat.lm<-(lm(Mand.2d.closingLR.res~Watershed*Habitat, data = kinem.df.res))
Mand.3d.openingLR.resxHabitat.lm<-(lm(Mand.3d.openingLR.res~Watershed*Habitat, data = kinem.df.res))
Mand.3d.closingLR.resxHabitat.lm<-(lm(Mand.3d.closingLR.res~Watershed*Habitat, data = kinem.df.res))
Mand.2d.openingLR.resxHabitat.lm
Anova(Mand.2d.openingLR.resxHabitat.lm, type = "II")
Mand.2d.closingLR.resxHabitat.lm
Anova(Mand.2d.closingLR.resxHabitat.lm, type = "II")
Mand.3d.openingLR.resxHabitat.lm
Anova(Mand.3d.openingLR.resxHabitat.lm, type = "II")
Mand.3d.closingLR.resxHabitat.lm
Anova(Mand.3d.closingLR.resxHabitat.lm, type = "II")

#  12) Biomechanical PERMANOVAs
KT.LR.2D.res.0.0<-cbind(kinem.df.res$Op4.2d.KT.res,kinem.df.res$Max4.2d.KT.res,
                        kinem.df.res$Mand.2d.openingLR.res,kinem.df.res$Mand.2d.closingLR.res)
KT.LR.2D.res.15.10<-cbind(kinem.df.res$Op4.2d.KT.15.res,kinem.df.res$Max4.2d.KT.10.res,
                          kinem.df.res$Mand.2d.openingLR.res,kinem.df.res$Mand.2d.closingLR.res)
KT.LR.2D.res.30.20<-cbind(kinem.df.res$Op4.2d.KT.30.res,kinem.df.res$Max4.2d.KT.20.res,
                          kinem.df.res$Mand.2d.openingLR.res,kinem.df.res$Mand.2d.closingLR.res)
KT.LR.3D.res.0.0<-cbind(kinem.df.res$Op4.3d.KT.res,kinem.df.res$Max4.3d.KT.res,
                        kinem.df.res$Mand.3d.openingLR.res,kinem.df.res$Mand.3d.closingLR.res)
KT.LR.3D.res.15.10<-cbind(kinem.df.res$Op4.3d.KT.15.res,kinem.df.res$Max4.3d.KT.10.res,
                          kinem.df.res$Mand.3d.openingLR.res,kinem.df.res$Mand.3d.closingLR.res)
KT.LR.3D.res.30.20<-cbind(kinem.df.res$Op4.3d.KT.30.res,kinem.df.res$Max4.3d.KT.20.res,
                          kinem.df.res$Mand.3d.openingLR.res,kinem.df.res$Mand.3d.closingLR.res)
# 2D Watershed*dietLD
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.2D.res.0.0~Watershed*LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.2D.res.15.10~ Watershed * LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.2D.res.30.20~ Watershed * LD1, SS.type = "II")))
# 3D Watershed*dietLD
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.3D.res.0.0~ Watershed * LD1,SS.type = "II")))
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.3D.res.15.10~ Watershed * LD1,SS.type = "II")))
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.3D.res.30.20~ Watershed * LD1,SS.type = "II")))
# 2D Watershed*simpson
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.2D.res.0.0~ Watershed * simpson, SS.type = "II")))
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.2D.res.15.10~ Watershed * simpson, SS.type = "II")))
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.2D.res.30.20~ Watershed * simpson, SS.type = "II")))
# 3D Watershed*simpson
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.3D.res.0.0~ Watershed * simpson, SS.type = "II")))
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.3D.res.15.10~ Watershed * simpson, SS.type = "II")))
summary(anova(lm.rrpp(data = kinem.df.res, KT.LR.3D.res.30.20~ Watershed * simpson, SS.type = "II")))


#  13) Dataframe of biomechanical disparity (3D values-2D values)
mech.disparity<-data.frame(kinem.df.res[,c("Population","Watershed", "Habitat","Spec", "LD1", "shannon", "simpson")],
                           Op4.KT.res.disp = kinem.df.res$Op4.3d.KT.res-kinem.df.res$Op4.2d.KT.res,
                           Op4.KT.15.res.disp = kinem.df.res$Op4.3d.KT.15.res-kinem.df.res$Op4.2d.KT.15.res,
                           Op4.KT.30.res.disp = kinem.df.res$Op4.3d.KT.30.res-kinem.df.res$Op4.2d.KT.30.res,
                           Op4.KT.3d15_2d.res.disp = kinem.df.res$Op4.3d.KT.15.res-kinem.df.res$Op4.2d.KT.res,
                           Op4.KT.3d30_2d.res.disp = kinem.df.res$Op4.3d.KT.30.res-kinem.df.res$Op4.2d.KT.res,
                           Max4.KT.res.disp = kinem.df.res$Max4.3d.KT.res-kinem.df.res$Max4.2d.KT.res,
                           Max4.KT.10.res.disp = kinem.df.res$Max4.3d.KT.10.res-kinem.df.res$Max4.2d.KT.10.res,
                           Max4.KT.20.res.disp = kinem.df.res$Max4.3d.KT.20.res-kinem.df.res$Max4.2d.KT.20.res,
                           Max4.KT.3d10_2d.res.disp = kinem.df.res$Max4.3d.KT.10.res-kinem.df.res$Max4.2d.KT.res,
                           Max4.KT.3d20_2d.res.disp = kinem.df.res$Max4.3d.KT.20.res-kinem.df.res$Max4.2d.KT.res,
                           Mand.openingLR.res.disp = kinem.df.res$Mand.3d.openingLR.res-kinem.df.res$Mand.2d.openingLR.res,
                           Mand.closingLR.res.disp = kinem.df.res$Mand.3d.closingLR.res-kinem.df.res$Mand.2d.closingLR.res
)

#  14) Disparity models analogous to 9-9.1
Op4.KT.res.disp.lm<-(lm(data = mech.disparity, Op4.KT.res.disp ~ Watershed * LD1))
Op4.KT.15.res.disp.lm<-(lm(data = mech.disparity, Op4.KT.15.res.disp ~ Watershed * LD1))
Op4.KT.30.res.disp.lm<-(lm(data = mech.disparity, Op4.KT.30.res.disp ~ Watershed * LD1))
Op4.KT.3d15_2d.res.disp.lm<-(lm(data = mech.disparity, Op4.KT.3d15_2d.res.disp ~ Watershed * LD1))
Op4.KT.3d30_2d.res.disp.lm<-(lm(data = mech.disparity, Op4.KT.3d30_2d.res.disp ~ Watershed * LD1))
Op4.KT.res.simpsondisp.lm<-(lm(data = mech.disparity, Op4.KT.res.disp ~ Watershed * simpson))
Op4.KT.15.res.simpsondisp.lm<-(lm(data = mech.disparity, Op4.KT.15.res.disp ~ Watershed * simpson))
Op4.KT.30.res.simpsondisp.lm<-(lm(data = mech.disparity, Op4.KT.30.res.disp ~ Watershed * simpson))
Op4.KT.3d15_2d.res.simpsondisp.lm<-(lm(data = mech.disparity, Op4.KT.3d15_2d.res.disp ~ Watershed * simpson))
Op4.KT.3d30_2d.res.simpsondisp.lm<-(lm(data = mech.disparity, Op4.KT.3d30_2d.res.disp ~ Watershed * simpson))
Op4.KT.res.disp.lm
Anova(Op4.KT.res.disp.lm, type = "II")
Op4.KT.15.res.disp.lm
Anova(Op4.KT.15.res.disp.lm, type = "II")
Op4.KT.30.res.disp.lm
Anova(Op4.KT.30.res.disp.lm, type = "II")
Op4.KT.3d15_2d.res.disp.lm
Anova(Op4.KT.3d15_2d.res.disp.lm, type = "II")
Op4.KT.3d30_2d.res.disp.lm
Anova(Op4.KT.3d30_2d.res.disp.lm, type = "II")
Op4.KT.res.simpsondisp.lm
Anova(Op4.KT.res.simpsondisp.lm, type = "II")
Op4.KT.15.res.simpsondisp.lm
Anova(Op4.KT.15.res.simpsondisp.lm, type = "II")
Op4.KT.30.res.simpsondisp.lm
Anova(Op4.KT.30.res.simpsondisp.lm, type = "II")
Op4.KT.3d15_2d.res.simpsondisp.lm
Anova(Op4.KT.3d15_2d.res.simpsondisp.lm, type = "II")
Op4.KT.3d30_2d.res.simpsondisp.lm
Anova(Op4.KT.3d30_2d.res.simpsondisp.lm, type = "II")

#  15) Disparity models analogous to 10-10.1
Max4.KT.res.disp.lm<-(lm(data = mech.disparity, Max4.KT.res.disp ~ Watershed * LD1))
Max4.KT.10.res.disp.lm<-(lm(data = mech.disparity, Max4.KT.10.res.disp ~ Watershed * LD1))
Max4.KT.20.res.disp.lm<-(lm(data = mech.disparity, Max4.KT.20.res.disp ~ Watershed * LD1))
Max4.KT.3d10_2d.res.disp.lm<-(lm(data = mech.disparity, Max4.KT.3d10_2d.res.disp ~ Watershed * LD1))
Max4.KT.3d20_2d.res.disp.lm<-(lm(data = mech.disparity, Max4.KT.3d20_2d.res.disp ~ Watershed * LD1))
Max4.KT.res.simpsondisp.lm<-(lm(data = mech.disparity, Max4.KT.res.disp ~ Watershed * simpson))
Max4.KT.10.res.simpsondisp.lm<-(lm(data = mech.disparity, Max4.KT.10.res.disp ~ Watershed * simpson))
Max4.KT.20.res.simpsondisp.lm<-(lm(data = mech.disparity, Max4.KT.20.res.disp ~ Watershed * simpson))
Max4.KT.3d10_2d.res.simpsondisp.lm<-(lm(data = mech.disparity, Max4.KT.3d10_2d.res.disp ~ Watershed * simpson))
Max4.KT.3d20_2d.res.simpsondisp.lm<-(lm(data = mech.disparity, Max4.KT.3d20_2d.res.disp ~ Watershed * simpson))
Max4.KT.res.disp.lm
Anova(Max4.KT.res.disp.lm, type = "II")
Max4.KT.10.res.disp.lm
Anova(Max4.KT.10.res.disp.lm, type = "II")
Max4.KT.20.res.disp.lm
Anova(Max4.KT.20.res.disp.lm, type = "II")
Max4.KT.3d10_2d.res.disp.lm
Anova(Max4.KT.3d10_2d.res.disp.lm, type = "II")
Max4.KT.3d20_2d.res.disp.lm
Anova(Max4.KT.3d20_2d.res.disp.lm, type = "II")
Max4.KT.res.simpsondisp.lm
Anova(Max4.KT.res.simpsondisp.lm, type = "II")
Max4.KT.10.res.simpsondisp.lm
Anova(Max4.KT.10.res.simpsondisp.lm, type = "II")
Max4.KT.20.res.simpsondisp.lm
Anova(Max4.KT.20.res.simpsondisp.lm, type = "II")
Max4.KT.3d10_2d.res.simpsondisp.lm
Anova(Max4.KT.3d10_2d.res.simpsondisp.lm, type = "II")
Max4.KT.3d20_2d.res.simpsondisp.lm
Anova(Max4.KT.3d20_2d.res.simpsondisp.lm, type = "II")

#  16) Disparity models analogous to 11-11.1
Mand.openingLR.res.disp.lm<-(lm(data = mech.disparity, Mand.openingLR.res.disp ~ Watershed * LD1))
Mand.closingLR.res.disp.lm<-(lm(data = mech.disparity, Mand.closingLR.res.disp ~ Watershed * LD1))
Mand.openingLR.res.simpsondisp.lm<-(lm(data = mech.disparity, Mand.openingLR.res.disp ~ Watershed * simpson))
Mand.closingLR.res.simpsondisp.lm<-(lm(data = mech.disparity, Mand.closingLR.res.disp ~ Watershed * simpson))
Mand.openingLR.res.disp.lm
Anova(Mand.openingLR.res.disp.lm, type = "II")
Mand.closingLR.res.disp.lm
Anova(Mand.closingLR.res.disp.lm, type = "II")
Mand.openingLR.res.simpsondisp.lm
Anova(Mand.openingLR.res.simpsondisp.lm, type = "II")
Mand.closingLR.res.simpsondisp.lm
Anova(Mand.closingLR.res.simpsondisp.lm, type = "II")

#  17) KT/LR DISPARITY PERMANOVAs, (dietLD, simpson, watershed*dietLD, watershed*simpson)
# 17.1)no tilted axes for Op4 and KT
disp.0.0<-cbind(mech.disparity$Op4.KT.res.disp,mech.disparity$Max4.KT.res.disp,
                mech.disparity$Mand.openingLR.res.disp,mech.disparity$Mand.closingLR.res.disp)
summary(anova(lm.rrpp(data = mech.disparity, disp.0.0~ LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.0.0~ simpson, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.0.0~ Watershed * LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.0.0~ Watershed * simpson, SS.type = "II")))
# 17.2) Op4 tilted 15º, Max4 tilted 10º
disp.15.10<-cbind(mech.disparity$Op4.KT.15.res.disp,mech.disparity$Max4.KT.10.res.disp,
                  mech.disparity$Mand.openingLR.res.disp,mech.disparity$Mand.closingLR.res.disp)
summary(anova(lm.rrpp(data = mech.disparity, disp.15.10~ LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.15.10~ simpson, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.15.10~ Watershed * LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.15.10~ Watershed * simpson, SS.type = "II")))
# 17.3) Op4 tilted 30º, Max4 tilted 20º
disp.30.20<-cbind(mech.disparity$Op4.KT.30.res.disp,mech.disparity$Max4.KT.20.res.disp,
                  mech.disparity$Mand.openingLR.res.disp,mech.disparity$Mand.closingLR.res.disp)
summary(anova(lm.rrpp(data = mech.disparity, disp.30.20~ LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.30.20~ simpson, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.30.20~ Watershed * LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.30.20~ Watershed * simpson, SS.type = "II")))
# 17.4) Op4 tilted 15º in 3D, planar in 2D; Max4 tilted 10º in 3D, planar in 2D
disp.15_0.10_0<-cbind(mech.disparity$Op4.KT.3d15_2d.res.disp,mech.disparity$Max4.KT.3d10_2d.res.disp,
                      mech.disparity$Mand.openingLR.res.disp,mech.disparity$Mand.closingLR.res.disp)
summary(anova(lm.rrpp(data = mech.disparity, disp.15_0.10_0~ LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.15_0.10_0~ simpson, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.15_0.10_0~ Watershed * LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.15_0.10_0~ Watershed * simpson, SS.type = "II")))
# 17.5) Op4 tileted 30º in 3D, planar in 2D; Max4 tilted 20º in 3D, planar in 2D
disp.30_0.20_0<-cbind(mech.disparity$Op4.KT.3d30_2d.res.disp,mech.disparity$Max4.KT.3d20_2d.res.disp,
                      mech.disparity$Mand.openingLR.res.disp,mech.disparity$Mand.closingLR.res.disp)
summary(anova(lm.rrpp(data = mech.disparity, disp.30_0.20_0~ LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.30_0.20_0~ simpson, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.30_0.20_0~ Watershed * LD1, SS.type = "II")))
summary(anova(lm.rrpp(data = mech.disparity, disp.30_0.20_0~ Watershed * simpson, SS.type = "II")))


## 18) Tests if disparity in calculated kinematic measures varies as consequence of diet LD  #####
summary(lm((Op4.3d.KT-Op4.2d.KT)~ldmk.factors$LD1))
summary(lm((Op4.3d.KT.15-Op4.2d.KT.15)~ldmk.factors$LD1))
summary(lm((Op4.3d.KT.15-Op4.2d.KT)~ldmk.factors$LD1))
summary(lm((Op4.3d.KT.30-Op4.2d.KT.30)~ldmk.factors$LD1))
summary(lm((Op4.3d.KT.30-Op4.2d.KT)~ldmk.factors$LD1))
summary(lm((Max4.3d.KT-Max4.2d.KT)~ldmk.factors$LD1))
summary(lm((Max4.3d.KT.10-Max4.2d.KT.10)~ldmk.factors$LD1))
summary(lm((Max4.3d.KT.10-Max4.2d.KT)~ldmk.factors$LD1))
summary(lm((Max4.3d.KT.20-Max4.2d.KT.20)~ldmk.factors$LD1))
summary(lm((Max4.3d.KT.20-Max4.2d.KT)~ldmk.factors$LD1))