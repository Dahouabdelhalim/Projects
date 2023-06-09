library(ggpubr)
library(data.table)

# PAIRWISE VECTOR COMPARISONS FOR SEPARATE SHAPE SPACES

# Inherited from "1-SeparateShapeSpaces.R"
# Op4.2dC.PW, Dent.2dC.PW, Max4.2dC.PW


# MORPHOLOGY VECTOR ANGLES
# 1) 2D
# 1.1) Op4
Op4.2dC.PW.angles<-summary(Op4.2dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi) #LakexLake pairwise angles (degrees)
Op4.2dC.PW.p<-summary(Op4.2dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P #LakexLake pairwise p-values
#     Figure S4, top left upper triangle (angles)
get_upper_tri<-function(Op4.2dC.PW.angles){
  Op4.2dC.PW.angles[upper.tri(Op4.2dC.PW.angles)] <- NA
  return(Op4.2dC.PW.angles)
}
Op4.2dC.angle.melted<-melt(get_upper_tri(Op4.2dC.PW.angles))
Op4.2dC.angle.melted[Op4.2dC.angle.melted==0]<-NA
Op4.2d.angle.matrix<-ggplot(data = Op4.2dC.angle.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 90, 
                                                    na.value = "transparent", limit = c(0,180), guide_colorbar(title = "Angle (deg)")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 3)
#     Figure S4, top left lower triangle (p-values)
Op4.2dC.PW.p
get_lower_tri<-function(Op4.2dC.PW.p){
  Op4.2dC.PW.p[lower.tri(Op4.2dC.PW.p)] <- NA
  return(Op4.2dC.PW.p)
}
Op4.2dC.p.melted<-melt(get_lower_tri(Op4.2dC.PW.p))
Op4.2dC.p.melted[Op4.2dC.p.melted==0]<-NA
Op4.2d.p.matrix<-ggplot(data = Op4.2dC.p.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

# 1.2) Dent
Dent.2dC.PW.angles<-summary(Dent.2dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Dent.2dC.PW.p<-summary(Dent.2dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
#     Figure S4, top center upper triangle (angles)
get_upper_tri<-function(Dent.2dC.PW.angles){
  Dent.2dC.PW.angles[upper.tri(Dent.2dC.PW.angles)] <- NA
  return(Dent.2dC.PW.angles)
}
Dent.2dC.angle.melted<-melt(get_upper_tri(Dent.2dC.PW.angles))
Dent.2dC.angle.melted[Dent.2dC.angle.melted==0]<-NA
Dent.2d.angle.matrix<-ggplot(data = Dent.2dC.angle.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 90, 
                                                    na.value = "transparent", limit = c(0,180), guide_colorbar(title = "Angle (deg)")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 3)
#     Figure S4, top center lower triangle (p-values)
Dent.2dC.PW.p
get_lower_tri<-function(Dent.2dC.PW.p){
  Dent.2dC.PW.p[lower.tri(Dent.2dC.PW.p)] <- NA
  return(Dent.2dC.PW.p)
}
Dent.2dC.p.melted<-melt(get_lower_tri(Dent.2dC.PW.p))
Dent.2dC.p.melted[Dent.2dC.p.melted==0]<-NA
Dent.2d.p.matrix<-ggplot(data = Dent.2dC.p.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

# 1.3) Max4
Max4bar.2dC.PW.angles<-summary(Max4.2dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Max4bar.2dC.PW.p<-summary(Max4.2dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
#     Figure S4, top right upper triangle (angles)
get_upper_tri<-function(Max4bar.2dC.PW.angles){
  Max4bar.2dC.PW.angles[upper.tri(Max4bar.2dC.PW.angles)] <- NA
  return(Max4bar.2dC.PW.angles)
}
Max4bar.2dC.angle.melted<-melt(get_upper_tri(Max4bar.2dC.PW.angles))
Max4bar.2dC.angle.melted[Max4bar.2dC.angle.melted==0]<-NA
Max4bar.2d.angle.matrix<-ggplot(data = Max4bar.2dC.angle.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 90, 
                                                    na.value = "transparent", limit = c(0,180), guide_colorbar(title = "Angle (deg)")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 3)
#     Figure S4, top right lower triangle (p-values)
Max4bar.2dC.PW.p
get_lower_tri<-function(Max4bar.2dC.PW.p){
  Max4bar.2dC.PW.p[lower.tri(Max4bar.2dC.PW.p)] <- NA
  return(Max4bar.2dC.PW.p)
}
Max4bar.2dC.p.melted<-melt(get_lower_tri(Max4bar.2dC.PW.p))
Max4bar.2dC.p.melted[Max4bar.2dC.p.melted==0]<-NA
Max4bar.2d.p.matrix<-ggplot(data = Max4bar.2dC.p.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

#
#
# 2) 3D
# 2.1)Op4
Op4.3dC.PW.angles<-summary(Op4.3dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Op4.3dC.PW.p<-summary(Op4.3dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
#     Figure S4 bottom left, upper triangle (angles)
get_upper_tri<-function(Op4.3dC.PW.angles){
  Op4.3dC.PW.angles[upper.tri(Op4.3dC.PW.angles)] <- NA
  return(Op4.3dC.PW.angles)
}
Op4.3dC.angle.melted<-melt(get_upper_tri(Op4.3dC.PW.angles))
Op4.3dC.angle.melted[Op4.3dC.angle.melted==0]<-NA
Op4.3d.angle.matrix<-ggplot(data = Op4.3dC.angle.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 90, 
                                                    na.value = "transparent", limit = c(0,180), guide_colorbar(title = "Angle (deg)")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 3)
#     Figure S4 bottom left, lower trangle (p-values)
Op4.3dC.PW.p
get_lower_tri<-function(Op4.3dC.PW.p){
  Op4.3dC.PW.p[lower.tri(Op4.3dC.PW.p)] <- NA
  return(Op4.3dC.PW.p)
}
Op4.3dC.p.melted<-melt(get_lower_tri(Op4.3dC.PW.p))
Op4.3dC.p.melted[Op4.3dC.p.melted==0]<-NA
Op4.3d.p.matrix<-ggplot(data = Op4.3dC.p.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)


# 2.2) Dent
Dent.3dC.PW.angles<-summary(Dent.3dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Dent.3dC.PW.p<-summary(Dent.3dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
#     Figure S4 bottom center, upper triangle (angles)
get_upper_tri<-function(Dent.3dC.PW.angles){
  Dent.3dC.PW.angles[upper.tri(Dent.3dC.PW.angles)] <- NA
  return(Dent.3dC.PW.angles)
}
Dent.3dC.angle.melted<-melt(get_upper_tri(Dent.3dC.PW.angles))
Dent.3dC.angle.melted[Dent.3dC.angle.melted==0]<-NA
Dent.3d.angle.matrix<-ggplot(data = Dent.3dC.angle.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 90, 
                                                    na.value = "transparent", limit = c(0,180), guide_colorbar(title = "Angle (deg)")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 3)
#     Figure S4 bottom center, lower triangle (p-values)
Dent.3dC.PW.p
get_lower_tri<-function(Dent.3dC.PW.p){
  Dent.3dC.PW.p[lower.tri(Dent.3dC.PW.p)] <- NA
  return(Dent.3dC.PW.p)
}
Dent.3dC.p.melted<-melt(get_lower_tri(Dent.3dC.PW.p))
Dent.3dC.p.melted[Dent.3dC.p.melted==0]<-NA
Dent.3d.p.matrix<-ggplot(data = Dent.3dC.p.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

# 2.3 Max4
Max4bar.3dC.PW.angles<-summary(Max4.3dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$angle*(180/pi)
Max4bar.3dC.PW.p<-summary(Max4.3dC.PW, confidence = 0.95, test.type = "VC", angle.type = "deg")$pairwise.tables$P
#     Figure S4 bottom right, upper triangle (angles)
get_upper_tri<-function(Max4bar.3dC.PW.angles){
  Max4bar.3dC.PW.angles[upper.tri(Max4bar.3dC.PW.angles)] <- NA
  return(Max4bar.3dC.PW.angles)
}
Max4bar.3dC.angle.melted<-melt(get_upper_tri(Max4bar.3dC.PW.angles))
Max4bar.3dC.angle.melted[Max4bar.3dC.angle.melted==0]<-NA
Max4bar.3d.angle.matrix<-ggplot(data = Max4bar.3dC.angle.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 90, 
                                                    na.value = "transparent", limit = c(0,180), guide_colorbar(title = "Angle (deg)")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,1)), color = "black", size = 3)
#     Figure S4 bottom right, lower triangle (p-values)
Max4bar.3dC.PW.p
get_lower_tri<-function(Max4bar.3dC.p){
  Max4bar.3dC.PW.p[lower.tri(Max4bar.3dC.PW.p)] <- NA
  return(Max4bar.3dC.PW.p)
}
Max4bar.3dC.p.melted<-melt(get_lower_tri(Max4bar.3dC.PW.p))
Max4bar.3dC.p.melted[Max4bar.3dC.p.melted==0]<-NA
Max4bar.3d.p.matrix<-ggplot(data = Max4bar.3dC.p.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

# 3) Arranges angle and p-value matrices into 3Var2 grid
angle.matrix<-ggarrange(Op4.2d.angle.matrix, Dent.2d.angle.matrix, Max4bar.2d.angle.matrix,
                        Op4.3d.angle.matrix, Dent.3d.angle.matrix, Max4bar.3d.angle.matrix,
                        ncol = 3, nrow = 2, common.legend = T, legend = "left")
pval.matrix<-ggarrange(Op4.2d.p.matrix, Dent.2d.p.matrix, Max4bar.2d.p.matrix,
                       Op4.3d.p.matrix, Dent.3d.p.matrix, Max4bar.3d.p.matrix,
                       ncol = 3, nrow = 2, common.legend = T, legend = "right")
#
#
#
#
#
# MORPHOLOGY VECTOR LENGTHS
# 4) 2D
# 4.1) OP4
Op4.2dC.PW.D<-summary(Op4.2dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Op4.2dC.PW.Dp<-summary(Op4.2dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$P
#     Figure S5, top left upper triangle (angles)
get_upper_tri<-function(Op4.2dC.PW.D){
  Op4.2dC.PW.D[upper.tri(Op4.2dC.PW.D)] <- NA
  return(Op4.2dC.PW.D)
}
Op4.2dC.D.melted<-melt(get_upper_tri(Op4.2dC.PW.D))
Op4.2dC.D.melted[Op4.2dC.D.melted==0]<-NA
Op4.2d.D.matrix<-ggplot(data = Op4.2dC.D.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red",midpoint = 0.04,
                                                    na.value = "transparent", limit = c(0,0.1), guide_colorbar(title = "Vector Length diff.")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)
#     Figure S5, top left lower triangle (p-values)
Op4.2dC.PW.Dp
get_lower_tri<-function(Op4.2dC.PW.Dp){
  Op4.2dC.PW.Dp[lower.tri(Op4.2dC.PW.Dp)] <- NA
  return(Op4.2dC.PW.Dp)
}
Op4.2dC.Dp.melted<-melt(get_lower_tri(Op4.2dC.PW.Dp))
Op4.2dC.Dp.melted[Op4.2dC.Dp.melted==0]<-NA
Op4.2d.Dp.matrix<-ggplot(data = Op4.2dC.Dp.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

# 4.2) Dent
Dent.2dC.PW.D<-summary(Dent.2dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Dent.2dC.PW.Dp<-summary(Dent.2dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$P
#     Figure S5, top center upper triangle (angles)
get_upper_tri<-function(Dent.2dC.PW.D){
  Dent.2dC.PW.D[upper.tri(Dent.2dC.PW.D)] <- NA
  return(Dent.2dC.PW.D)
}
Dent.2dC.D.melted<-melt(get_upper_tri(Dent.2dC.PW.D))
Dent.2dC.D.melted[Dent.2dC.D.melted==0]<-NA
Dent.2d.D.matrix<-ggplot(data = Dent.2dC.D.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 0.04,
                                                    na.value = "transparent", limit = c(0,0.1), guide_colorbar(title = "Vector Length diff.")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)
#     Figure S5, top center lower triangle (p-values)
Dent.2dC.PW.Dp
get_lower_tri<-function(Dent.2dC.PW.Dp){
  Dent.2dC.PW.Dp[lower.tri(Dent.2dC.PW.Dp)] <- NA
  return(Dent.2dC.PW.Dp)
}
Dent.2dC.Dp.melted<-melt(get_lower_tri(Dent.2dC.PW.Dp))
Dent.2dC.Dp.melted[Dent.2dC.Dp.melted==0]<-NA
Dent.2d.Dp.matrix<-ggplot(data = Dent.2dC.Dp.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

# 4.3) Max4
Max4bar.2dC.PW.D<-summary(Max4.2dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Max4bar.2dC.PW.Dp<-summary(Max4.2dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$P
#     Figure S5, top right upper triangle (angles)
get_upper_tri<-function(Max4bar.2dC.PW.D){
  Max4bar.2dC.PW.D[upper.tri(Max4bar.2dC.PW.D)] <- NA
  return(Max4bar.2dC.PW.D)
}
Max4bar.2dC.D.melted<-melt(get_upper_tri(Max4bar.2dC.PW.D))
Max4bar.2dC.D.melted[Max4bar.2dC.D.melted==0]<-NA
Max4bar.2d.D.matrix<-ggplot(data = Max4bar.2dC.D.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red",midpoint = 0.04,
                                                    na.value = "transparent", limit = c(0,0.1), guide_colorbar(title = "Vector Length diff.")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)
#     Figure S5, top right lower triangle (p-values)
Max4bar.2dC.PW.Dp
get_lower_tri<-function(Max4bar.2dC.PW.Dp){
  Max4bar.2dC.PW.Dp[lower.tri(Max4bar.2dC.PW.Dp)] <- NA
  return(Max4bar.2dC.PW.Dp)
}
Max4bar.2dC.Dp.melted<-melt(get_lower_tri(Max4bar.2dC.PW.Dp))
Max4bar.2dC.Dp.melted[Max4bar.2dC.Dp.melted==0]<-NA
Max4bar.2d.Dp.matrix<-ggplot(data = Max4bar.2dC.Dp.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)


# 5) 3D
# 5.1) Op4
Op4.3dC.PW.D<-summary(Op4.3dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Op4.3dC.PW.Dp<-summary(Op4.3dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$P
#     Figure S5, bottom left upper triangle (angles)
get_upper_tri<-function(Op4.3dC.PW.D){
  Op4.3dC.PW.D[upper.tri(Op4.3dC.PW.D)] <- NA
  return(Op4.3dC.PW.D)
}
Op4.3dC.D.melted<-melt(get_upper_tri(Op4.3dC.PW.D))
Op4.3dC.D.melted[Op4.3dC.D.melted==0]<-NA
Op4.3d.D.matrix<-ggplot(data = Op4.3dC.D.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 0.04,
                                                    na.value = "transparent", limit = c(0,0.1), guide_colorbar(title = "Vector Length diff.")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text.y = element_blank(),axis.text.x = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)
#     Figure S5, bottom left upper triangle (p-values)
Op4.3dC.PW.Dp
get_lower_tri<-function(Op4.3dC.PW.Dp){
  Op4.3dC.PW.Dp[lower.tri(Op4.3dC.PW.Dp)] <- NA
  return(Op4.3dC.PW.Dp)
}
Op4.3dC.Dp.melted<-melt(get_lower_tri(Op4.3dC.PW.Dp))
Op4.3dC.Dp.melted[Op4.3dC.Dp.melted==0]<-NA
Op4.3d.Dp.matrix<-ggplot(data = Op4.3dC.Dp.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)


# 5.2) Dent
Dent.3dC.PW.D<-summary(Dent.3dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Dent.3dC.PW.Dp<-summary(Dent.3dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$P
#     Figure S5, bottom center upper triangle (angles)
get_upper_tri<-function(Dent.3dC.PW.D){
  Dent.3dC.PW.D[upper.tri(Dent.3dC.PW.D)] <- NA
  return(Dent.3dC.PW.D)
}
Dent.3dC.D.melted<-melt(get_upper_tri(Dent.3dC.PW.D))
Dent.3dC.D.melted[Dent.3dC.D.melted==0]<-NA
Dent.3d.D.matrix<-ggplot(data = Dent.3dC.D.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 0.04,
                                                    na.value = "transparent", limit = c(0,0.1), guide_colorbar(title = "Vector Length diff.")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)
#     Figure S5, bottom center lower triangle (p-values)
Dent.3dC.PW.Dp
get_lower_tri<-function(Dent.3dC.PW.Dp){
  Dent.3dC.PW.Dp[lower.tri(Dent.3dC.PW.Dp)] <- NA
  return(Dent.3dC.PW.Dp)
}
Dent.3dC.Dp.melted<-melt(get_lower_tri(Dent.3dC.PW.Dp))
Dent.3dC.Dp.melted[Dent.3dC.p.melted==0]<-NA
Dent.3d.Dp.matrix<-ggplot(data = Dent.3dC.Dp.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

# 5.3) Max4
Max4bar.3dC.PW.D<-summary(Max4.3dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$D
Max4bar.3dC.PW.Dp<-summary(Max4.3dC.PW, confidence = 0.95, test.type = "dist")$pairwise.tables$P
#     Figure S5, bottom right upper triangle (angles)
get_upper_tri<-function(Max4bar.3dC.PW.D){
  Max4bar.3dC.PW.D[upper.tri(Max4bar.3dC.PW.D)] <- NA
  return(Max4bar.3dC.PW.D)
}
Max4bar.3dC.D.melted<-melt(get_upper_tri(Max4bar.3dC.PW.D))
Max4bar.3dC.D.melted[Max4bar.3dC.D.melted==0]<-NA
Max4bar.3d.D.matrix<-ggplot(data = Max4bar.3dC.D.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "top") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "steelblue",
                                                    mid = "white", high = "red", midpoint = 0.04,
                                                    na.value = "transparent", limit = c(0,0.1), guide_colorbar(title = "Vector Length diff.")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "top")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)
#     Figure S5, bottom right lower triangle (p-values)
Max4bar.3dC.PW.Dp
get_lower_tri<-function(Max4bar.3dC.PW.Dp){
  Max4bar.3dC.PW.Dp[lower.tri(Max4bar.3dC.PW.Dp)] <- NA
  return(Max4bar.3dC.PW.Dp)
}
Max4bar.3dC.Dp.melted<-melt(get_lower_tri(Max4bar.3dC.PW.Dp))
Max4bar.3dC.Dp.melted[Max4bar.3dC.Dp.melted==0]<-NA
Max4bar.3d.Dp.matrix<-ggplot(data = Max4bar.3dC.Dp.melted, aes(Var2, Var1, fill = value))+scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(low = "red",
                                                    mid = "white", high = "steelblue", midpoint = .3, 
                                                    na.value = "transparent", limit = c(0,1), guide_legend(title = "p Value")) +
  theme_minimal() + xlab(NULL)+ylab(NULL)+ theme(panel.grid.major = element_blank(),axis.text = element_blank(),legend.position = "bottom")+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "black", size = 3)

# 6) Arranges angle and p-value matrices into 3Var2 grid
dist.matrix<-ggarrange(Op4.2d.D.matrix, Dent.2d.D.matrix, Max4bar.2d.D.matrix,
                       Op4.3d.D.matrix, Dent.3d.D.matrix, Max4bar.3d.D.matrix,
                       ncol = 3, nrow = 2, common.legend = T, legend = "left")
distpval.matrix<-ggarrange(Op4.2d.Dp.matrix, Dent.2d.Dp.matrix, Max4bar.2d.Dp.matrix,
                           Op4.3d.Dp.matrix, Dent.3d.Dp.matrix, Max4bar.3d.Dp.matrix,
                           ncol = 3, nrow = 2, common.legend = T, legend = "right")


