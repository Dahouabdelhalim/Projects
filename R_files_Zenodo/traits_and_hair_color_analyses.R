##############Trait and hair color ANALYSES###########################
# Would be repeated for each body region, as well as the sets of discrete character datasets
##############FIRST STEP: SINGLE VS MULTIPLE OU#############################
# IMPORT DATASET 
RGB_fore = read.csv("RGB_traits_fore_scores.csv")
#or
RGB_head = read.csv("RGB_traits_head_scores.csv")
#or
RGB_hind = read.csv("RGB_traits_hind_scores.csv")
#or
RGB_tail = read.csv("RGB_traits_tail_scores.csv")
#or
RGB_torso = read.csv("RGB_traits_torso_scores.csv")

# ASSIGN SPECIES AS ROWNAMES TO YOUR DATA FRAME
rownames(RGB_fore) = RGB_fore$X
#or
rownames(RGB_head) = RGB_head$X
#or
rownames(RGB_hind) = RGB_hind$X
#or
rownames(RGB_tail) = RGB_tail$X
#or
rownames(RGB_torso) = RGB_torso$X

#Import tree and align new tree with data--MAY BE REDUDANT, SKIP
primate_color_tree=read.nexus("consensusTree_10kTrees_Primates_Version3.nex")
name.check(primate_color_tree, RGB_fore)
primate_color_tree_fore = drop.tip(primate_color_tree, c("whichever primate species don't match"))
name.check(primate_color_tree_fore, RGB_fore)

#Make sure you have no polytomies in your tree
primate_color_tree_fore2 <- multi2di(primate_color_tree_fore, random = TRUE)
primate_color_tree_fore2 <- multi2di(primate_color_tree_fore, random = FALSE)
is.binary.tree(primate_color_tree_fore2)

#Write the tree
write.tree(primate_color_tree_fore2, "forecolor_tree.txt")

###Note: these trees already exist so you likely don't have to do previous lines but should double check. Otherwise just import your tree:
primate_color_tree_fore2=read.tree("forecolor_tree.txt")
name.check(primate_color_tree_fore2, RGB_fore)
#or
primate_color_tree_head2=read.tree("headcolor_tree2.txt")
name.check(primate_color_tree_head2, RGB_head)
#or
primate_color_tree_hind2=read.tree("hindcolor_tree.txt")
name.check(primate_color_tree_hind2, RGB_hind)
#or
primate_color_tree_tail2=read.tree("tailcolor_tree.txt")
name.check(primate_color_tree_tail2, RGB_tail)
#or
primate_color_tree_torso2=read.tree("torsocolor_tree.txt")
name.check(primate_color_tree_torso2, RGB_torso)

# MATCH ORDER OF TAXA IN DATASET AND TREE ***IMPORTANT****
match(rownames(RGB_fore),primate_color_tree_fore2$tip.label)
RGB_fore_sort <- RGB_fore[match(primate_color_tree_fore2$tip.label,rownames(RGB_fore)),]
match(rownames(RGB_fore_sort),primate_color_tree_fore2$tip.label)
#or
match(rownames(RGB_head),primate_color_tree_head2$tip.label)
RGB_head_sort <- RGB_head[match(primate_color_tree_head2$tip.label,rownames(RGB_head)),]
match(rownames(RGB_head_sort),primate_color_tree_head2$tip.label)
#or
match(rownames(RGB_hind),primate_color_tree_hind2$tip.label)
RGB_hind_sort <- RGB_hind[match(primate_color_tree_hind2$tip.label,rownames(RGB_hind)),]
match(rownames(RGB_hind_sort),primate_color_tree_hind2$tip.label)
#or
match(rownames(RGB_tail),primate_color_tree_tail2$tip.label)
RGB_tail_sort <- RGB_tail[match(primate_color_tree_tail2$tip.label,rownames(RGB_tail)),]
match(rownames(RGB_tail_sort),primate_color_tree_tail2$tip.label)
#or
match(rownames(RGB_torso),primate_color_tree_torso2$tip.label)
RGB_torso_sort <- RGB_torso[match(primate_color_tree_torso2$tip.label,rownames(RGB_torso)),]
match(rownames(RGB_torso_sort),primate_color_tree_torso2$tip.label)

# CREATE A VECTOR WITH YOUR DISCRETE VARIABLE ONLY, Vision, Habitat, Clade, Activity Pattern
RGB_fore_vision = RGB_fore_sort$Vision
#or
RGB_fore_activity = RGB_fore_sort$Activity
#or
RGB_fore_habitat = RGB_fore_sort$Habitat_final
#or
RGB_fore_clade = RGB_fore_sort$Clade
#or
RGB_fore_s_activity = RGB_fore_sort$Simple_Activity

#or
RGB_head_vision = RGB_head_sort$Vision
#or
RGB_head_activity = RGB_head_sort$Activity
#or
RGB_head_habitat = RGB_head_sort$Habitat_final
#or
RGB_head_clade = RGB_head_sort$Clade
#or
RGB_head_s_activity = RGB_head_sort$Simple_Activity

#or
RGB_hind_vision = RGB_hind_sort$Vision
#or
RGB_hind_activity = RGB_hind_sort$Activity
#or
RGB_hind_habitat = RGB_hind_sort$Habitat_final
#or 
RGB_hind_clade=RGB_hind_sort$Clade
#or
RGB_hind_s_activity = RGB_hind_sort$Simple_Activity

#or
RGB_tail_vision = RGB_tail_sort$Vision
#or
RGB_tail_activity = RGB_tail_sort$Activity
#or
RGB_tail_habitat = RGB_tail_sort$Habitat_final
#or
RGB_tail_clade = RGB_tail_sort$Clade
#or
RGB_tail_s_activity = RGB_tail_sort$Simple_Activity

#or
RGB_torso_vision = RGB_torso_sort$Vision
#or
RGB_torso_activity = RGB_torso_sort$Activity
#or
RGB_torso_habitat = RGB_torso_sort$Habitat_final
#or
RGB_torso_clade = RGB_torso_sort$Clade
#or
RGB_torso_s_activity = RGB_torso_sort$Simple_Activity

#ASSIGN SPECIES NAMES AS ATTRIBUTES TO THE VECTOR
names(RGB_fore_vision) = RGB_fore_sort$X
#or
names(RGB_fore_activity) = RGB_fore_sort$X
#or
names(RGB_fore_habitat) = RGB_fore_sort$X
#or
names(RGB_fore_clade) = RGB_fore_sort$X
#or
names(RGB_fore_s_activity) = RGB_fore_sort$X


#or
names(RGB_head_vision) = RGB_head_sort$X
#or
names(RGB_head_activity) = RGB_head_sort$X
#or
names(RGB_head_habitat) = RGB_head_sort$X
#or
names(RGB_head_clade) = RGB_head_sort$X
#or
names(RGB_head_s_activity) = RGB_head_sort$X

#or
names(RGB_hind_vision) = RGB_hind_sort$X
#or
names(RGB_hind_activity) = RGB_hind_sort$X
#or
names(RGB_hind_habitat) = RGB_hind_sort$X
#or
names(RGB_hind_clade) = RGB_hind_sort$X
#or
names(RGB_hind_s_activity) = RGB_hind_sort$X

#or
names(RGB_tail_vision) = RGB_tail_sort$X
#or
names(RGB_tail_activity) = RGB_tail_sort$X
#or
names(RGB_tail_habitat) = RGB_tail_sort$X
#or
names(RGB_tail_clade) = RGB_tail_sort$X
#or
names(RGB_tail_s_activity) = RGB_tail_sort$X

#or
names(RGB_torso_vision) = RGB_torso_sort$X
#or
names(RGB_torso_activity) = RGB_torso_sort$X
#or
names(RGB_torso_habitat) = RGB_torso_sort$X
#or
names(RGB_torso_clade) = RGB_torso_sort$X
#or
names(RGB_torso_s_activity) = RGB_torso_sort$X

###Make your simmap for analysis
RGB_fore_tree_simmap<-make.simmap(primate_color_tree_fore2, RGB_fore_vision, model="ER", nsim=1)
#or
RGB_fore_tree_simmap<-make.simmap(primate_color_tree_fore2, RGB_fore_activity, model="ER", nsim=1)
#or
RGB_fore_tree_simmap<-make.simmap(primate_color_tree_fore2, RGB_fore_habitat, model="ER", nsim=1)
#or
RGB_fore_tree_simmap<-make.simmap(primate_color_tree_fore2, RGB_fore_clade, model="ER", nsim=1)
#or
RGB_fore_tree_simmap<-make.simmap(primate_color_tree_fore2, RGB_fore_s_activity, model="ER", nsim=1)

#or
RGB_head_tree_simmap<-make.simmap(primate_color_tree_head2, RGB_head_vision, model="ER", nsim=1)
#or
RGB_head_tree_simmap<-make.simmap(primate_color_tree_head2, RGB_head_activity, model="ER", nsim=1)
#or
RGB_head_tree_simmap<-make.simmap(primate_color_tree_head2, RGB_head_habitat, model="ER", nsim=1)
#or
RGB_head_tree_simmap<-make.simmap(primate_color_tree_head2, RGB_head_clade, model="ER", nsim=1)
#or
RGB_head_tree_simmap<-make.simmap(primate_color_tree_head2, RGB_head_s_activity, model="ER", nsim=1)


#or
RGB_hind_tree_simmap<-make.simmap(primate_color_tree_hind2, RGB_hind_vision, model="ER", nsim=1)
#or
RGB_hind_tree_simmap<-make.simmap(primate_color_tree_hind2, RGB_hind_activity, model="ER", nsim=1)
#or
RGB_hind_tree_simmap<-make.simmap(primate_color_tree_hind2, RGB_hind_habitat, model="ER", nsim=1)
#or
RGB_hind_tree_simmap<-make.simmap(primate_color_tree_hind2, RGB_hind_clade, model="ER", nsim=1)
#or
RGB_hind_tree_simmap<-make.simmap(primate_color_tree_hind2, RGB_hind_s_activity, model="ER", nsim=1)

#or
RGB_tail_tree_simmap<-make.simmap(primate_color_tree_tail2, RGB_tail_vision, model="ER", nsim=1)
#or
RGB_tail_tree_simmap<-make.simmap(primate_color_tree_tail2, RGB_tail_activity, model="ER", nsim=1)
#or
RGB_tail_tree_simmap<-make.simmap(primate_color_tree_tail2, RGB_tail_habitat, model="ER", nsim=1)
#or
RGB_tail_tree_simmap<-make.simmap(primate_color_tree_tail2, RGB_tail_clade, model="ER", nsim=1)
#or
RGB_tail_tree_simmap<-make.simmap(primate_color_tree_tail2, RGB_tail_s_activity, model="ER", nsim=1)


#or
RGB_torso_tree_simmap<-make.simmap(primate_color_tree_torso2, RGB_torso_vision, model="ER", nsim=1)
#or
RGB_torso_tree_simmap<-make.simmap(primate_color_tree_torso2, RGB_torso_activity, model="ER", nsim=1)
#or
RGB_torso_tree_simmap<-make.simmap(primate_color_tree_torso2, RGB_torso_habitat, model="ER", nsim=1)
#or
RGB_torso_tree_simmap<-make.simmap(primate_color_tree_torso2, RGB_torso_clade, model="ER", nsim=1)
#or
RGB_torso_tree_simmap<-make.simmap(primate_color_tree_torso2, RGB_torso_s_activity, model="ER", nsim=1)


plot(RGB_fore_tree_simmap, fsize=0.5)
summary(RGB_fore_tree_simmap)
#or
plot(RGB_head_tree_simmap, fsize=0.5)
summary(RGB_head_tree_simmap)
#or
plot(RGB_hind_tree_simmap, fsize=0.5)
summary(RGB_hind_tree_simmap)
#or
plot(RGB_tail_tree_simmap, fsize=0.5)
summary(RGB_tail_tree_simmap)
#or
plot(RGB_torso_tree_simmap, fsize=0.5)
summary(RGB_torso_tree_simmap)

###Test single OU model
RGB_fore_OU1<- mvOU(RGB_fore_tree_simmap, RGB_fore_sort$PC1, model="OU1", diagnostic=FALSE, echo=FALSE)
RGB_fore_OU1
#or
RGB_head_OU1<- mvOU(RGB_head_tree_simmap, RGB_head_sort$PC1, model="OU1", diagnostic=FALSE, echo=FALSE)
RGB_head_OU1
#or
RGB_hind_OU1<- mvOU(RGB_hind_tree_simmap, RGB_hind_sort$PC1, model="OU1", diagnostic=FALSE, echo=FALSE)
RGB_hind_OU1
#or
RGB_tail_OU1<- mvOU(RGB_tail_tree_simmap, RGB_tail_sort$PC1, model="OU1", diagnostic=FALSE, echo=FALSE)
RGB_tail_OU1
#or
RGB_torso_OU1<- mvOU(RGB_torso_tree_simmap, RGB_torso_sort$PC1, model="OU1", diagnostic=FALSE, echo=FALSE)
RGB_torso_OU1

###Now test multiple OU model
RGB_fore_OUM<- mvOU(RGB_fore_tree_simmap, RGB_fore_sort$PC1, model="OUM", diagnostic=FALSE, echo=FALSE)
RGB_fore_OUM
#or
RGB_head_OUM<- mvOU(RGB_head_tree_simmap, RGB_head_sort$PC1, model="OUM", diagnostic=FALSE, echo=FALSE)
RGB_head_OUM
#or
RGB_hind_OUM<- mvOU(RGB_hind_tree_simmap, RGB_hind_sort$PC1, model="OUM", diagnostic=FALSE, echo=FALSE)
RGB_hind_OUM
#or
RGB_tail_OUM<- mvOU(RGB_tail_tree_simmap, RGB_tail_sort$PC1, model="OUM", diagnostic=FALSE, echo=FALSE)
RGB_tail_OUM
#or
RGB_torso_OUM<- mvOU(RGB_torso_tree_simmap, RGB_torso_sort$PC1, model="OUM", diagnostic=FALSE, echo=FALSE)
RGB_torso_OUM


######################SECOND STEP: TESTING FOR DIFFERENT RATES OF EVOLUTION AMONG CLADES###############
#LOAD PHYTOOLS
# TESTING FOR DIFFERENT RATES OF EVOLUTION AMONG CLADES
brownie_RGB_fore = brownie.lite(RGB_fore_tree_simmap, RGB_fore_sort$PC1, maxit=2000, test="simulation", nsim=100)
brownie_RGB_fore
#or
brownie_RGB_head = brownie.lite(RGB_head_tree_simmap, RGB_head_sort$PC1, maxit=2000, test="simulation", nsim=100)
brownie_RGB_head
#or
brownie_RGB_hind = brownie.lite(RGB_hind_tree_simmap, RGB_hind_sort$PC1, maxit=2000, test="simulation", nsim=100)
brownie_RGB_hind
#or
brownie_RGB_tail = brownie.lite(RGB_tail_tree_simmap, RGB_tail_sort$PC1, maxit=2000, test="simulation", nsim=100)
brownie_RGB_tail
#or
brownie_RGB_torso = brownie.lite(RGB_torso_tree_simmap, RGB_torso_sort$PC1, maxit=2000, test="simulation", nsim=100)
brownie_RGB_torso
