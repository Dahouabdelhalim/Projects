#Script for project
#This is sample code; would be repeated for each body part in RGB Values
primate_color_tree=read.nexus("consensusTree_10kTrees_Primates_Version3.nex")
RGB_traits=read.csv("HaircolorRGB.csv")

#Get rid of excess blank rows
RGB_traits[-c(95:100),]
RGB_traits_total=RGB_traits[-c(95:100),]
rownames(RGB_traits_total)=RGB_traits_total$ï..species

#align species names in dataset with tip labels in tree
name.check(primate_color_tree, RGB_traits_total)

#Organize data to be in same order as tree
RGB_traits_total <- RGB_traits_total[match(primate_color_tree$tip.label,rownames(RGB_traits_total)),]
RGB_traits_m=data.matrix(RGB_traits_total)
#Specify which variables you want to look at, cut out missing data
RGB_traits_head=RGB_traits_m[,c(3:8,39:41)]

#Need to log transform data because of presence of zeros
RGB_traits_head=log((RGB_traits_head)+0.01)

#Need to omit rows with NA
RGB_traits_head <- na.omit(RGB_traits_head)

#Find correlation coefficients for related body region variables
round(cor(RGB_traits_head),2)

#Excise any correlated data over .9
RGB_traits_head=RGB_traits_head[,c(2,5,7,9)]

#Run the PCA
RGB_traits_head_pca=prcomp(RGB_traits_head,scale=TRUE)

#Do diagnostics
summary(RGB_traits_head_pca)                 # variance explained by the components

plot(RGB_traits_head_pca,type="lines")       # "scree" plot of eigenvalues

screeplot(RGB_traits_head_pca, type="lines") # same

RGB_traits_head_pca$rotation                 # eigenvectors for each component(with trait loadings)

RGB_traits_head_pca$sdev^2                   # eigenvalues for each component(variances)

predict(RGB_traits_head_pca)                 # PCA scores

#export PCA scores to csv file

write.csv(predict(RGB_traits_head_pca), "RGB_traits_head_scores.csv")

#run biplot
biplot(RGB_traits_head_pca, cex=0.5)

#Assess phylogenetic signal
RGB_traits_head_scores=read.csv("RGB_traits_head_scores.csv")
View(RGB_traits_head_scores)
RGB_traits_head_scores2=RGB_traits_head_scores$PC1
RGB_traits_head_scores2=matrix(RGB_traits_head_scores2)
rownames(RGB_traits_head_scores2)=RGB_traits_head_scores$X
View(RGB_traits_head_scores2)

#But before that, since you got rid of some rows you need to align new tree with data
name.check(primate_color_tree, RGB_traits_head_scores2)
primate_color_tree_head = drop.tip(primate_color_tree, c("Leontopithecus_rosalia","Macaca_nemestrina","Miopithecus_talapoin","Papio_anubis","Pithecia_monachus","Semnopithecus_entellus","Trachypithecus_phayrei","Trachypithecus_pileatus"))
name.check(primate_color_tree_head, RGB_traits_head_scores2)   #should say "OK" now

#RUn K and lambda analyses
phylosig(primate_color_tree_head, RGB_traits_head_scores2, method="K", test=TRUE, nsim=1000)
phylosig(primate_color_tree_head, RGB_traits_head_scores2, method="lambda", test=TRUE, nsim=1000)

#Plot values on Contmap
primate_color_tree_head=read.tree("headcolor_tree.txt")
RGB_traits_head_scores = read.csv("RGB_traits_head_scores.csv")

RGB_traits_head_scores2=RGB_traits_head_scores$PC1
names(RGB_traits_head_scores2)=RGB_traits_head_scores$X
RGBHeadPCA<-contMap(primate_color_tree_head, RGB_traits_head_scores2)
RGBHeadPCA<-setMap(RGBHeadPCA,
                       c("gray10","white"))
plot(RGBHeadPCA,fsize=c(0.5,0.7),
     leg.txt="Trait Value")

RGB_traits_head_scores2=RGB_traits_head_scores$PC1
names(RGB_traits_head_scores2)=RGB_traits_head_scores$X
RGBHeadPCA<-contMap(primate_color_tree_head, RGB_traits_head_scores2)
RGBHeadPCA<-setMap(RGBHeadPCA,
                   c("gray10","white"))
plot(RGBHeadPCA,fsize=c(0.5,0.7),
     leg.txt="Trait Value")

############Model Trait Evolution############
#First make sure you have no polytomies in your tree
primate_color_tree_head2 <- multi2di(primate_color_tree_head, random = TRUE)
primate_color_tree_head2 <- multi2di(primate_color_tree_head, random = FALSE)
is.binary.tree(primate_color_tree_head2)

#Then force the tree to be ultrametric
library(phangorn)
primate_color_tree_head3<-force.ultrametric(primate_color_tree_head2) ## default method
is.ultrametric(primate_color_tree_head3)

#Model BM
modelBM_RGB_head<-fitContinuous(primate_color_tree_head3, RGB_traits_head_scores2, model = "BM")
#Get Diagnostics
modelBM_RGB_head$opt

#Model OU
modelOU_RGB_head<-fitContinuous(primate_color_tree_head3, RGB_traits_head_scores2, model = "OU")
#Get Diagnostics
modelOU_RGB_head$opt

#Model EB
modelEB_RGB_head<-fitContinuous(primate_color_tree_head3, RGB_traits_head_scores2, model = "EB")
#Get Diagnostics
modelEB_RGB_head$opt

#Model "white" aka null model
modelwhite_RGB_head<-fitContinuous(primate_color_tree_head3, RGB_traits_head_scores2, model = "white")
#Get Diagnostics
modelwhite_RGB_head$opt


##To make useful biplots of PCA Results:
RGB_traits_torso_scores=read.csv("RGB_traits_torso_scores.csv")
#Get the percent of variance explained for PC1 and 2:
summary(RGB_traits_torso_pca)
#This line and the next one gives a title to the legend
scale_label_Clade = "Clade"
pdf("TorsoColorBiplot.pdf", height=8, width=11) # Tells R to save graph as pdf in your working directory
qplot(data=RGB_traits_torso_scores, x=RGB_traits_torso_scores$PC1, y=RGB_traits_torso_scores$PC2, color=RGB_traits_torso_scores$Clade)+stat_ellipse()+ 
  #Replace PC_dataset with dataset containing your PC scores and other relevant information. Stat_ellipse draws a 95%
  # confidence interval around the groups specified by color
  geom_point(size = 5) +  
  #Give Clade Groupings Color
  scale_color_manual(scale_label_Clade, values=c("#d55e00", "#cc79a7", "#0072b2", "#009e73"))+
  #Label your axes
  xlab("PC1 (66.34% of variance)")+
  ylab("PC2 (21.12% of variance)")+  
  # Format biplot:
  theme(axis.title.x=element_text(size=20))+
  theme(axis.title.y=element_text(size=20))+
  theme(legend.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+
  theme(axis.text.x=element_text(colour="black", size=18))+
  theme(axis.text.y=element_text(colour="black", size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")) # Changes background into white/black instead of default gray that ggplot has 
dev.off()
