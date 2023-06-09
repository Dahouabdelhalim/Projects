# CORN, MARTINEZ, BURRESS, and WAINWRIGHT at Systematic Biology
#doi: https://doi.org/10.1093/sysbio/syaa091

########################################
#         LOADING IN THE DATA          #
########################################

#This section is code from C. Martinez, modified by K. Corn

#load packages
require(geomorph)
require(tidyverse)

#load data table
species_table= read.csv(".../Corn_et_al_SysBio_scales.csv", header=T) #this is the table that contains the scale info
names(species_table)
species_table <- species_table[-c(37:42),] #take out dactylopus because it's wacky and screws everything up
species_table[85,2] 

require(car) #just fixing some species whose names are spelled wrong
species_table$species<-recode(species_table$species,"'Paracentropodon_rubripinnis'='Paracentropogon_rubripinnis'")
species_table$species<-recode(species_table$species,"'Siganus_sp_black'='Siganus_uspi'")
species_table$species<-recode(species_table$species,"'Siganus_sp_white'='Siganus_virgatus'")
species_table$species<-recode(species_table$species,"'Inermia_vittata'='Haemulon_vittatum'")
species_table$species<-recode(species_table$species,"'Centrogenys_species'='Centrogenys_vaigiensis'")
species_table$species<-recode(species_table$species,"'Chromileptes_altivelis'='Cromileptes_altivelis'")
species_table$species<-recode(species_table$species,"'Pomacanthus_xanthometapon'='Pomacanthus_xanthometopon'")
species_table$species<-recode(species_table$species,"'Parupeneus_cyclostoma'='Parupeneus_cyclostomus'")
species_table$species_specimen<-recode(species_table$species_specimen,"'Emmelichthyops_atlanticus_18'='Emmelichthyops_atlanticus_8'")
#if these start giving you problems, replace "recode" with "car::recode". recode can conflict with dplyr recode function

#make a vector w order that everything is in, then go through every object and compare and make sure that object is in the right order matching that vector

species_table$video_id <- paste0(species_table$species_specimen,"_",seq_along(species_table$species))

nrow(species_table)
n_strike= 175 #should be the same number as nrow(species_table)

counter= c(1:n_strike)
strike_labels= as.data.frame(matrix(nrow=(n_strike*10), ncol=1)) # table that complements data of 1750 aligned shapes
for(i in 1:n_strike){
  num= counter[i]
  strike_labels[c(1:10) + (i-1)*(10),]= rep(num,10)
}

#read in landmarks
coord_data= readland.tps(file=".../Corn_et_al_SysBio_LMs.TPS",specID = "imageID") 
minus_dactylopus <- coord_data[,,-c(361:420)]
#read in data table of the semilandmarks
semiland= as.matrix(read.table(".../Corn_et_al_SysBio_sliding_semi.txt", header=T))

unaligned_2d= two.d.array(minus_dactylopus) 



n_specimen=175

scaled= as.data.frame(matrix(nrow=1750, ncol=36))
for(i in 1:n_specimen){
  sc = species_table[i,4] #same as "video labels"
  shapes=  unaligned_2d[c(1:10) + (i-1)*(10),]
  scaled_shape= shapes/sc
  scaled[c(1:10) + (i-1)*(10),]= scaled_shape
}

scaled
rownames(scaled) <- rownames(unaligned_2d)

p=18 #number of landmarks
k=2 #number of lms
scaled_3d= arrayspecs(scaled, p=p, k=k)

#<---------------- GPA then put coords in 2D table ---------------->#
y= gpagen(scaled_3d, curves=semiland, ProcD=T) # $coords and $Csize  
y2= two.d.array(y$coords) 

tang = gm.prcomp(A=y$coords)
plot(tang, axis1 = 1, axis2 = 2)

n_specimen=175

## table of Procrustes distances b/t successive motion points
kinetic_procrustes= as.data.frame(matrix(ncol=n_specimen+1, nrow=12))
for(f in 1:n_specimen){
  spec= y2[10*(f-1) +c(1:10),]
  for(i in 1:9){	
    proc= sqrt(sum((spec[i,]-spec[i+1,])^2))
    kinetic_procrustes[i,f+1]= proc # successive procrustes distances b/t steps
    kinetic_procrustes[10,f+1]= sum(kinetic_procrustes[1:9,f+1]) # total path length
    kinetic_procrustes[11,f+1]= sqrt(sum((spec[1,]-spec[10,])^2)) # procrustes first to last shapes (linear vector)
    kinetic_procrustes[12,f+1]= kinetic_procrustes[11,f+1]/kinetic_procrustes[10,f+1] # ratio of linear to nonlinear  
  }}
kinetic_procrustes[,1]= c(1:9, "total", "linear","ratio")
colnames(kinetic_procrustes) <- c("dist",species_table$video_id)
# kinetic_procrustes	# rows are distances, columns are specimens
kinetic_procrustes2 = as.data.frame(t(kinetic_procrustes[,2:(n_specimen+1)])) # rows are specimens, columns are distances

kinetic_procrustes_total= as.data.frame(kinetic_procrustes2[,10])
efficiency_ratio_total <- as.data.frame(kinetic_procrustes2[,12])

require(ggplot2)

l= density(kinetic_procrustes_total[,1])
t=ggplot(kinetic_procrustes_total, aes(x=kinetic_procrustes_total[,1])) + geom_density(size=1.2, adjust=1.2) + labs(x="kinesis",y="density") + xlim(range(l$x)) + theme_classic(base_size=14)
t


indiv_kinesis <- as.vector(kinetic_procrustes2[,10])


#check dfs order
kin_order_check <- data.frame(species_table$species,rownames(kinetic_procrustes2))

#<---------------- MAKE DATA FRAME ---------------->#


data <- as.data.frame(species_table$species)
data[,2] <- species_table$species_specimen
data[,3] <- species_table$feeding
data[,4] <- kinetic_procrustes_total
colnames(data) <- c("species","specimen","type","kinesis")

#<---------------- JAW PROTRUSION ---------------->#

LM4_x= y2[,7] # head of max LM
LM4_y= y2[,8]

LM5_x= y2[,9] # premax LM
LM5_y= y2[,10]

premax_protrusion= as.data.frame(sqrt((( LM5_x - LM4_x)^2) + (( LM5_y - LM4_y)^2))) # premax protrusion for each image (n=1360)

premax_by_strike= as.data.frame(matrix(nrow=n_specimen, ncol=1)) # max protrusion within each strike
for(i in 1:n_specimen){
  spec= premax_protrusion[10*(i-1) +c(1:10),]
  maximum_premax= max(spec)
  minimum_premax= min(spec)
  premax_by_strike[i,]= maximum_premax - minimum_premax
}

data[,5] <- premax_by_strike[,1]

#<---------------- JAW LENGTH ---------------->#

require(data.table)
LM7_x= y2[,13] # jaw joint LM
LM7_y= y2[,14]

LM18_x= y2[,35] # distal tip of dentary LM
LM18_y= y2[,36]

jaw_length= as.data.frame(sqrt((( LM18_x - LM7_x)^2) + (( LM18_y - LM7_y)^2))) # jaw length for each image (n=n_specimen)

length_by_strike= as.data.frame(matrix(nrow=n_specimen, ncol=1)) # max protrusion within each strike
for(i in 1:n_specimen){
  spec= jaw_length[10*(i-1) +c(1:10),]
  starting_length= first(spec)
  length_by_strike[i,]= starting_length
}

data[,6] <- length_by_strike[,1]

#<---------------- MAXIMUM GAPE ---------------->#

LM5_x= y2[,9] # tip of premaxilla LM
LM5_y= y2[,10]

LM18_x= y2[,35] # distal tip of dentary LM
LM18_y= y2[,36]

max_gape= as.data.frame(sqrt((( LM18_x - LM5_x)^2) + (( LM18_y - LM5_y)^2))) # gape for each image (n=n_specimen0)

gape_by_strike= as.data.frame(matrix(nrow=n_specimen, ncol=1)) # max protrusion within each strike
for(i in 1:n_specimen){
  spec= max_gape[10*(i-1) +c(1:10),]
  maximum_gape= max(spec)
  minimum_gape= min(spec)
  gape_by_strike[i,]= maximum_gape - minimum_gape
}

data[,7] <- gape_by_strike[,1]

#<---------------- MAXILLARY ROTATION ---------------->#

LM4_x= y2[,7] # head of max LM  (VERTEX)
LM4_y= y2[,8]

LM3_x= y2[,5] # eye LM
LM3_y= y2[,6]

LM6_x= y2[,11] # distal max LM
LM6_y= y2[,12]

#distances
Dist_vert_2= sqrt((( LM3_x - LM4_x)^2) + (( LM3_y - LM4_y)^2)) # vertex and second points 
Dist_vert_3= sqrt((( LM6_x - LM4_x)^2) + (( LM6_y - LM4_y)^2)) # vertex and third points
Dist_2_3= sqrt((( LM6_x - LM3_x)^2) + (( LM6_y - LM3_y)^2)) # second and third points

angle= acos(((Dist_vert_2^2) + (Dist_vert_3^2) - (Dist_2_3^2))/(2*Dist_vert_2*Dist_vert_3))
maxilla= as.data.frame(angle * (180/pi))


maxilla_by_strike= as.data.frame(matrix(nrow=n_specimen, ncol=1)) # maxillary rotation within each strike
for(i in 1:n_specimen){
  spec= maxilla[10*(i-1) +c(1:10),]
  maximum_maxilla= max(spec)
  minimum_maxilla= min(spec)
  maxilla_by_strike[i,]= maximum_maxilla - minimum_maxilla
}

data[,8] <- maxilla_by_strike[,1]

#<---------------- LOWER JAW ROTATION ---------------->#

LM7_x= y2[,13] # quadrate-articular joint LM  (VERTEX)
LM7_y= y2[,14]

LM3_x= y2[,5] # eye LM
LM3_y= y2[,6]

LM18_x= y2[,35] # distal lower jaw LM
LM18_y= y2[,36]

#distances
Dist_vert_2= sqrt((( LM3_x - LM7_x)^2) + (( LM3_y - LM7_y)^2)) # vertex and second points 
Dist_vert_3= sqrt((( LM18_x - LM7_x)^2) + (( LM18_y - LM7_y)^2)) # vertex and third points
Dist_2_3= sqrt((( LM18_x - LM3_x)^2) + (( LM18_y - LM3_y)^2)) # second and third points

angle= acos(((Dist_vert_2^2) + (Dist_vert_3^2) - (Dist_2_3^2))/(2*Dist_vert_2*Dist_vert_3))
lowerjaw= as.data.frame(angle * (180/pi))


lowerjaw_by_strike= as.data.frame(matrix(nrow=n_specimen, ncol=1)) # lower jaw rotation within each strike
for(u in 1:n_specimen){
  spec= lowerjaw[10*(u-1) +c(1:10),]
  maximum_lowerjaw= max(spec)
  minimum_lowerjaw= min(spec)
  lowerjaw_by_strike[u,]= maximum_lowerjaw - minimum_lowerjaw
}

data[,9] <- lowerjaw_by_strike[,1]

#<---------------- CRANIAL ELEVATION ---------------->#

LM2_x= y2[,3] # cleithrum-post temporal joint LM  (VERTEX)
LM2_y= y2[,4]

LM3_x= y2[,5] # eye LM
LM3_y= y2[,6]

LM8_x= y2[,15] # pectoral fin LM
LM8_y= y2[,16]

#distances
Dist_vert_2= sqrt((( LM3_x - LM2_x)^2) + (( LM3_y - LM2_y)^2)) # vertex and second points 
Dist_vert_3= sqrt((( LM8_x - LM2_x)^2) + (( LM8_y - LM2_y)^2)) # vertex and third points
Dist_2_3= sqrt((( LM8_x - LM3_x)^2) + (( LM8_y - LM3_y)^2)) # second and third points

angle= acos(((Dist_vert_2^2) + (Dist_vert_3^2) - (Dist_2_3^2))/(2*Dist_vert_2*Dist_vert_3))
headr= as.data.frame(angle * (180/pi))


headr_by_strike= as.data.frame(matrix(nrow=n_specimen, ncol=1)) # head rotation within each strike
for(u in 1:n_specimen){
  spec= headr[10*(u-1) +c(1:10),]
  maximum_headr= max(spec)
  minimum_headr= min(spec)
  headr_by_strike[u,]= maximum_headr - minimum_headr
}

data[,10] <- headr_by_strike[,1]

#<---------------- HYOID DEPRESSION ---------------->#

LM7_x= y2[,13] # quadrate-articular joint LM (base)
LM7_y= y2[,14]

LM9_x= y2[,17] # pelvic LM (base)
LM9_y= y2[,18]

LM13_x= y2[,25] # 5th ventral landmark, from pelvic (variable)
LM13_y= y2[,26]

#distances
Dist_B= sqrt((( LM7_x - LM9_x)^2) + (( LM7_y - LM9_y)^2)) # base length
Dist_A= sqrt((( LM9_x - LM13_x)^2) + (( LM9_y - LM13_y)^2)) # pelvic to 5th ventral
Dist_C= sqrt((( LM7_x - LM13_x)^2) + (( LM7_y - LM13_y)^2)) # quadrate-articular to 5th ventral

s= (Dist_A + Dist_B + Dist_C)/2 # half perimeter
area= sqrt(s*(s-Dist_A)*(s-Dist_B)*(s-Dist_C))# Heron's formula 
height= as.data.frame((2*area)/Dist_B) # re-arrangement of the formula, Area=1/2bh, to solve for h

hyoid_by_strike= as.data.frame(matrix(nrow=n_specimen, ncol=1)) # hyoid depressin within each strike
for(u in 1:n_specimen){
  spec= height[10*(u-1) +c(1:10),]
  maximum_hyoid= max(spec)
  minimum_hyoid= min(spec)
  hyoid_by_strike[u,]= maximum_hyoid - minimum_hyoid
}

data[,11] <- hyoid_by_strike[,1]

#<---------------- OTHER CALCULATIONS ---------------->#


data[,12] <- data[,8]/data[,9]
data[,13] <- data[,5]/data[,4]
data[,14] <- efficiency_ratio_total


#<---------------- MAXILLA LENGTH ---------------->#

LM4_x= y2[,7] # maxilla base LM
LM4_y= y2[,8]

LM6_x= y2[,11] # distal tip of maxilla LM
LM6_y= y2[,12]

maxilla_length= as.data.frame(sqrt((( LM6_x - LM4_x)^2) + (( LM6_y - LM4_y)^2))) # maxilla length for each image (n=n_specimen)

max_length_by_strike= as.data.frame(matrix(nrow=n_specimen, ncol=1)) # max protrusion within each strike
for(i in 1:n_specimen){
  spec= maxilla_length[10*(i-1) +c(1:10),]
  starting_length= first(spec)
  max_length_by_strike[i,]= starting_length
}

data[,15] <- max_length_by_strike
data[,16] <- (data[,8]/360)*2*pi*data[,15]

colnames(data) <- c("species","specimen","type","kinesis","protrusion","jaw_length","gape","maxillary_rotation","lower_jaw_rotation","head_rotation","hyoid_depression","KT","ratio","efficiency", "maxillary_length","maxillary_arc")

kinesis <- data %>% 
  dplyr::group_by(specimen) %>% 
  dplyr::mutate_at(c("kinesis","protrusion","jaw_length","gape","maxillary_rotation","lower_jaw_rotation","head_rotation","hyoid_depression","KT","ratio","efficiency","maxillary_length","maxillary_arc"),mean) %>%
  dplyr::group_by(species,type) %>% 
  dplyr::summarize_at(c("kinesis","protrusion","jaw_length","gape","maxillary_rotation","lower_jaw_rotation","head_rotation","hyoid_depression","KT","ratio","efficiency","maxillary_length","maxillary_arc"),mean) %>%
  ungroup() #%>% 


kinesis




########################################
#         ADDING THE PHYLOGENY         #
########################################

#This section written by K. Corn
require(ape)
require(phytools)
require(geiger)
require(ggrepel)

tree1 <- read.tree("~/Desktop/Ch 1 - Reef Fish Kinesis/code + source data/new_rabosky_tree.tre")

species<-as.vector(kinesis$species)

#fixing some tree labels that don't match the taxa I have
tree1$tip.label<-car::recode(tree1$tip.label,"'Choerodon_fasciatus'='Choerodon_cyanodus'")
tree1$tip.label<-car::recode(tree1$tip.label,"'Notocirrhitus_splendens'='Cyprinocirrhites_polyactis'")
tree1$tip.label<-car::recode(tree1$tip.label,"'Paracirrhites_hemistictus'='Oxycirrhites_typus'")
tree1$tip.label<-car::recode(tree1$tip.label,"'Bodianus_tanyokidus'='Terelabrus_flavocephalus'")

trim_tree<-drop.tip(tree1, setdiff(tree1$tip.label, species))
trim_tree$tip.label

kinesis_tree <- ladderize(trim_tree, right=FALSE)
plot(kinesis_tree)

phyl_kindata <- kinesis %>%
  arrange(match(species, kinesis_tree$tip.label))
table(phyl_kindata$species == kinesis_tree$tip.label) #this had better match



kinesis_data <- as.data.frame(phyl_kindata)





########################################
#         TANGENT SPACE (Fig. 1)       #
########################################

#This section written by C. Martinez and K. Corn

feeding_labels= as.data.frame(matrix(nrow=(n_strike*10), ncol=1)) # table that complements data of 1750 aligned shapes
for(i in 1:n_strike){
  rw2= as.character(species_table[i,6]) 
  feeding_labels[c(1:10) + (i-1)*(10),]= rep(rw2,10)
}

frame_point= rep(c("strt","mid","mid","mid","mid","mid","mid","mid","mid","end"), n_specimen) #for size
clr= c("gold3","turqoise4") #for coloring your groups

ggplot(as.data.frame(tang$x), aes(x=Comp1, y=Comp2, color = feeding_labels[,1], shape=frame_point)) + geom_point(aes(size=frame_point), alpha=0.5) + scale_colour_manual(values=c("gold3","turquoise4")) + scale_shape_manual(values=c(16,16,1)) + scale_size_manual(values=c(3,1,2))+ geom_path(aes(group=as.factor(strike_labels[,1])), size=0.4, alpha=0.5) + coord_equal() + theme_classic() + theme(legend.position = "none") #+ geom_text(aes(label=as.factor(1:1750)), vjust=-1)

#you can highlight a single strike at a time by just adjusting the dataset you color. make the alpha on everything else 0.2 and make your strike of interest red or something like that.





########################################
#         Warp Grids                   #
########################################

#This section written by C. Martinez

## warp grids of PC extremes
PC1_min= tang$shapes$shapes.comp1$min
PC1_max= tang$shapes$shapes.comp1$max
PC2_min= tang$shapes$shapes.comp2$min
PC2_max= tang$shapes$shapes.comp2$max
PC3_min= tang$shapes$shapes.comp3$min
PC3_max= tang$shapes$shapes.comp3$max
PC4_min= tang$shapes$shapes.comp4$min
PC4_max= tang$shapes$shapes.comp4$max

mean_shape= mshape(y$coords)

lnks= as.data.frame(matrix(nrow=16, ncol=2))
lnks[,1]= c(2,2,1,4,5,6,7,9,10,11,12,13,14,15,16,17)
lnks[,2]= c(8,3,4,5,6,4,18,10,11,12,13,14,15,16,17,18)

#plot pc extreme against mean shape
plotRefToTarget(M1=mean_shape , M2=PC1_max, mag=1, links=lnks, gridPar=list(n.col.cell=10, grid.col="gray60", grid.lwd=3.2, tar.link.lwd=6, tar.link.col="black", tar.pt.size=1.7, tar.pt.bg="mediumorchid2"))

#see the mean shape
plotRefToTarget(M1=mean_shape ,M2=mean_shape, mag=1, links=lnks, gridPar=list(n.col.cell=10, grid.col="gray60", grid.lwd=3.2, tar.link.lwd=6, tar.link.col="olivedrab", tar.pt.size=1.7, tar.pt.bg="dodgerblue4"))





########################################
#         Head Shape                   #
########################################


#This section written by K. Corn

#We make a separate alignment of starting head shapes

first_frames <- as.data.frame(scaled) %>% 
  rownames_to_column(., var = "video") %>% 
  filter(str_detect(video, "1$"))
first_frames_named <- as.data.frame(first_frames[,2:37],row.names = first_frames[,1])
first_frames_3d= arrayspecs(first_frames_named, p=18, k=2)
first_frames_aligned= gpagen(first_frames_3d, curves=semiland)
first_frames_2d= two.d.array(first_frames_aligned$coords)

#next, average by individual

mean_first_frames <- as_tibble(as.data.frame(first_frames_2d) %>%
                                 rownames_to_column(., var = "video") %>% 
                                 mutate(species = species_table$species,specimen = species_table$species_specimen, feeding = species_table$feeding) %>%
                                 dplyr::group_by(specimen) %>%
                                 mutate_at(., vars(`1.X`:`18.Y`),mean) %>%
                                 dplyr::group_by(species,feeding) %>%
                                 summarize_at(., vars(`1.X`:`18.Y`),mean) %>%
                                 ungroup())

#Order head shape data with phylogeny

tree_shape <- mean_first_frames %>%
  arrange(match(species, kinesis_tree$tip.label))
table(tree_shape$species == kinesis_tree$tip.label)

matrix_tree_shape <- data.frame(tree_shape[,3:38],row.names = tree_shape$species)

#PC scores
tree_shape_3d= arrayspecs(matrix_tree_shape, p=18, k=2)
tree_tang_pcs <- gm.prcomp(A = tree_shape_3d)
plot(tree_tang_pcs)

tree_tang_pcs$x #component scores for all specimens

#in order to make a pca
ggplot(as.data.frame(tree_tang_pcs$x), aes(x=Comp1, y=Comp2, color = tree_shape$feeding)) + geom_point(size=3,alpha=0.8) + scale_colour_manual(values=c("gold3","turquoise4"))  + theme_classic() + geom_text_repel(aes(label=tree_shape$species),show.legend = F)


#I *highly* recommend doing some checks to make sure that everything is running properly
coris_mat <- as.matrix(tree_shape_3d[,,"Coris_formosa"]) #set up a matrix with a species that you *know* what the species is
hemi_mat <- as.matrix(tree_shape_3d[,,"Hemitaurichthys_zoster"])

test_mat <- first_frames_aligned$coords[,,"Coris_formosa_01_09_1"] #set up a matrix with a head shape from that species

plotRefToTarget(coris_mat,test_mat) #make sure that the resemblance is 1:1
plotRefToTarget(hemi_mat,test_mat) #an example of what it looks like when it's not a match




########################################
#         PCA                          #
########################################

#This section written by K. Corn

reef_pca <- prcomp(~protrusion+gape+maxillary_rotation+lower_jaw_rotation+hyoid_depression+head_rotation, data=kinesis_data, scale=TRUE)




########################################
#      Phylogenetic ANOVA              #
########################################

#This section written by K. Corn

require(geomorph)

procD_test_func_out <- list()

pgls_kin_mat <- as.data.frame(
  phyl_kindata %>%
    dplyr::select(kinesis:hyoid_depression))
rownames(pgls_kin_mat) <- phyl_kindata$species

x = phyl_kindata$type
names(x) <- phyl_kindata$species
y <- pgls_kin_mat
for(i in 1:8){
  out <- procD.pgls(y[i] ~ x, phy = kinesis_tree, iter=100000, seed = NULL)
  procD_test_func_out[i] <- out
  i <- i + 1
}




########################################
#      convergent evolution            #
########################################

#This section written by K. Corn

require(convevol)

#convevol with kinesis & component data

biting_tibble <- phyl_kindata %>% 
  filter(type == "bite")
biting_list <- as.character(biting_tibble$species)

convtraits <- phyl_kindata %>%
  dplyr::select(protrusion,maxillary_rotation,lower_jaw_rotation,head_rotation,gape,hyoid_depression,kinesis)
convtraits_df <- as.matrix(convtraits)
rownames(convtraits_df) <- phyl_kindata$species

#convrat: distance based metrics of convergence
reef_convrat <- convrat(phyl = kinesis_tree, phendata = convtraits_df, convtips = biting_list)
reef_convratsig <- convratsig(phyl = kinesis_tree, phendata = convtraits_df, convtips = biting_list,nsim=500)






########################################
#      disparity analyses              #
########################################

#This section written by K. Corn

disp_kin_mat <- as.data.frame(
  phyl_kindata %>%
    dplyr::select(kinesis:hyoid_depression, - jaw_length))
rownames(disp_kin_mat) <- phyl_kindata$species

#set up df
disparity_out <- data.frame(rep(NA,7))
disparity_out$trait <- rep(NA,7)
disparity_out$var_biters <- rep(NA,7)
disparity_out$var_suction <- rep(NA,7)
disparity_out$props <- rep(NA,7)
disparity_out$p.dfc <- rep(NA,7)

var_names_vec <- c("kinesis", "protrusion", "gape", "maxillary_rotation", "lower_jaw_rotation", "head_rotation", "hyoid_depression")

iter = 1
feeding4disp <- phyl_kindata$type
names(feeding4disp) <- phyl_kindata$species
#disp_kin_mat
for(iter in 1:7){
  morphdisp_data <- as.vector(disp_kin_mat[,iter]) #make data into a format geomorph will talk to
  names(morphdisp_data) <- rownames(disp_kin_mat)
  out <- morphol.disparity(morphdisp_data ~ feeding4disp, groups= ~ feeding4disp, iter = 10000, print.progress = F)
  disparity_out$trait[iter] <- var_names_vec[iter]
  disparity_out$var_biters[iter] <- out$Procrustes.var[[1]]
  disparity_out$var_suction[iter] <- out$Procrustes.var[[2]]
  disparity_out$p.dfc[iter] <- out$PV.dist.Pval[1,2]
  disparity_out$props[iter] <- out$Procrustes.var[[2]]/out$Procrustes.var[[1]]
  iter <- iter + 1
}

#mean ratio of disparities
mean(disparity_out$props)


########################################
#      head shape rates                #
########################################

#This section written by K. Corn

matrix_tree_shape <- data.frame(tree_shape[,3:38],row.names = tree_shape$species)

feeding <- as.vector(tree_shape$feeding)
names(feeding) <- tree_shape$species

morph_rates <- compare.evol.rates(
  A = matrix_tree_shape,
  phy = kinesis_tree,
  gp = feeding,
  iter = 10000,
  seed = NULL,
  method = c("simulation"),
  print.progress = TRUE
)





########################################
#     stochastic character mapping     #
########################################

#This section written by K. Corn

# Make simmaps

feeding_type <- as.vector(kinesis_data$type)
names(feeding_type) <- kinesis_data$species

require(phytools)
rfk_simmaps <- make.simmap(tree = kinesis_tree, x = feeding_type, Q="empirical", pi = "estimated", nsim = 1000)


cols<-setNames(c("gold3","turquoise4"),c("bite", "suction"))
plotSimmap(rfk_simmaps[[1]], cols, ftype = "off", type = "fan", offset = 1)




########################################
#            running OUwie             #
########################################

#This section written by K. Corn


#Set up OUwie
require(OUwie)


#Set up data frame to hold output
best_model_data <- tibble(.rows = 35000)
best_model_data$trait <- rep(NA,35000)
best_model_data$simmap_count <- rep(NA,35000)
best_model_data$model <- rep(NA,35000)
best_model_data$loglik <- rep(NA,35000)
best_model_data$aicc <- rep(NA,35000)
best_model_data$eigval <- rep(NA, 35000)
best_model_data$sigma.sq_biters <- rep(NA, 35000)
best_model_data$sigma.sq_suction <- rep(NA, 35000)
best_model_data$alpha <- rep(NA,35000)
best_model_data$theta_biters <- rep(NA,35000)
best_model_data$theta_suction <- rep(NA,35000)
best_model_data$theta_biters_se <- rep(NA,35000)
best_model_data$theta_suction_se <- rep(NA,35000)
best_model_data$saddle <- rep(NA,35000)

OUwie.model<-function(model, phy, data) {
  print(paste("Now starting model",model))
  return(OUwie(phy, data, model, simmap.tree=T, diagn=T))
}
models <- c("BM1","BMS","OU1","OUM","OUMV")

OUwie_dataprep <- as.data.frame(phyl_kindata %>%
                                  dplyr::select(species, type, kinesis, protrusion, gape, maxillary_rotation, lower_jaw_rotation, head_rotation, hyoid_depression))

var_names_vec <- c("kinesis", "protrusion", "gape", "maxillary_rotation", "lower_jaw_rotation", "head_rotation", "hyoid_depression")


j = 1
i = 1
r = 1
q = 1

for (j in 1:7){
  print(paste("Now starting trait: ",var_names_vec[j]))
  OUwie_data <- data.frame(OUwie_dataprep[,c(1,2,j+2)])
  i = 1
  for (i in 1:1000){
    print(paste("Now starting SIMMAP: ",i))
    tree <- rfk_simmaps[[i]]
    results <- lapply(models, OUwie.model, phy=tree, data=OUwie_data)
    q = 1
    
    for (q in 1:5) {
      each_model <- results[[q]]
      
      best_model_data$trait[r] <- var_names_vec[j]
      best_model_data$simmap_count[r] <- i
      best_model_data$model[r] <- each_model$model
      best_model_data$loglik[r] <- each_model$loglik
      best_model_data$aicc[r] <- each_model$AICc
      best_model_data$eigval[r] <- list(each_model$eigval)
      best_model_data$sigma.sq_biters[r] <- each_model$solution[2,1]
      best_model_data$sigma.sq_suction[r] <- each_model$solution[2,2]
      best_model_data$alpha[r] <- each_model$solution[1,1]
      best_model_data$theta_biters[r] <- each_model$theta[1,1]
      best_model_data$theta_suction[r] <- each_model$theta[2,1]
      best_model_data$theta_biters_se[r] <- each_model$theta[1,2]
      best_model_data$theta_suction_se[r] <- each_model$theta[2,2]
      best_model_data$saddle[r] <- any(each_model$eigval < 0) #all should be over to 0 if well-fit healthy OU model/good run
      
      q <- q + 1
      r <- r + 1
    }
    i <- i + 1
  }
  j <- j + 1
}

best_model_data






########################################
#           OUwie simulations          #
########################################

#This section written by K. Corn

require(phytools)

feeding_type <- as.vector(kinesis_data$type)
names(feeding_type) <- kinesis_data$species

rfk_simmaps <- make.simmap(tree = kinesis_tree, x = feeding_type, Q="empirical", pi = "estimated", nsim = 101)
sim_sum <- summary(rfk_simmaps)

cols<-setNames(c("gold3","turquoise4"),c("bite", "suction"))
plotSimmap(rfk_simmaps[[1]], cols, ftype = "off", type = "fan", offset = 1)

#<-------------------------- MAKE DFS -------------------------->

#set up df to catch the output
sim_out_data <- tibble(.rows = 2500) #I think this will be of length 5 x 3 x 3 x 5 = 225
sim_out_data$sim_under_model <- rep(NA,2500)
sim_out_data$sim_count <- rep(NA,2500)
sim_out_data$simmap_count <- rep(NA,2500)
sim_out_data$model_run <- rep(NA,2500)
sim_out_data$loglik <- rep(NA,2500)
sim_out_data$aicc <- rep(NA,2500)
sim_out_data$eigval <- rep(NA, 2500)
sim_out_data$sigma.sq_biters <- rep(NA, 2500)
sim_out_data$sigma.sq_suction <- rep(NA, 2500)
sim_out_data$alpha <- rep(NA,2500)
sim_out_data$theta_biters <- rep(NA,2500)
sim_out_data$theta_suction <- rep(NA,2500)
sim_out_data$theta_biters_se <- rep(NA,2500)
sim_out_data$theta_suction_se <- rep(NA,2500)
sim_out_data$saddle <- rep(NA,2500)

#set up simulated data frame
sim_params <- tibble(.rows = 5) 
sim_params$model <- c("BM1","BMS","OU1","OUM","OUMV")
sim_params$alpha = rep(NA, 5)
sim_params$sigma.sq = rep(NA, 5)
sim_params$theta = rep(NA, 5)
sim_params$theta0 = rep(NA, 5)

sim_params$alpha[c(1,2)] <- list(c(1e-10,1e-10)) #alpha = 1e-10 for BM1 & BMS
sim_params$alpha[3] <- list(c(0.1,0.1)) #alpha = 0.1 for OU1
sim_params$alpha[c(4,5)] <- list(c(1.0,1.0)) #alpha = 1.0 for OUM and OUMV
sim_params$sigma.sq[c(1,4)] <- list(c(0.45,0.45))
sim_params$sigma.sq[3] <- list(c(0.9,0.9))
sim_params$sigma.sq[c(2,5)] <- list(c(0.45,0.9))
sim_params$theta0 <- 1.0
sim_params$theta <- c(1.0, 0.5, 1.0, 0.5, 0.5)
sim_params$theta[c(1,2)] <- list(c(0.0,0.0))
sim_params$theta[3] <- list(c(1.0,1.0))
sim_params$theta[c(4,5)] <- list(c(1.0,2.0))



#<-------------------------- RUN OUwie -------------------------->
require(OUwie)

OUwie.model<-function(model, phy, data) {
  print(paste("Now starting model",model))
  return(OUwie(phy, data, model, simmap.tree=T, diagn=T))
}
models <- c("BM1","BMS","OU1","OUM","OUMV")

i = 1
r = 1
q = 1
m = 1
for (i in 1:3){
  print(paste("Now starting SIMMAP: ",i))
  tree <- rfk_simmaps[[i]]
  m = 1 
  for (m in 1:5){
    print(paste("Now simulating dataset under model: ",(models[m])))
    
    #set up simulated data
    simulated_data <- OUwie.sim(phy = rfk_simmaps[[i]], simmap.tree = T,scaleHeight=FALSE, alpha=sim_params$alpha[[m]], sigma.sq = sim_params$sigma.sq[[m]], theta0=sim_params$theta0[m], theta=sim_params$theta[[m]])		
    sim_data_3cols <- simulated_data
    sim_data_3cols[,3] <- simulated_data[,2]
    sim_data_3cols[,2] <- phyl_kindata$type #fixes to errors in simulations
    
    results <- lapply(models, OUwie.model, phy=tree, data=sim_data_3cols)
    q = 1
    
    for (q in 1:5) {
      each_model <- results[[q]]
      
      sim_out_data$sim_under_model[r] <- models[m]
      sim_out_data$sim_count[r] <- m
      sim_out_data$simmap_count[r] <- i
      sim_out_data$model_run[r] <- each_model$model
      sim_out_data$loglik[r] <- each_model$loglik
      sim_out_data$aicc[r] <- each_model$AICc
      sim_out_data$eigval[r] <- list(each_model$eigval)
      sim_out_data$sigma.sq_biters[r] <- each_model$solution[2,1]
      sim_out_data$sigma.sq_suction[r] <- each_model$solution[2,2]
      sim_out_data$alpha[r] <- each_model$solution[1,1]
      sim_out_data$theta_biters[r] <- each_model$theta[1,1]
      sim_out_data$theta_suction[r] <- each_model$theta[2,1]
      sim_out_data$theta_biters_se[r] <- each_model$theta[1,2]
      sim_out_data$theta_suction_se[r] <- each_model$theta[2,2]
      sim_out_data$saddle[r] <- any(each_model$eigval < 0) #Should be over  0 if well-fit healthy OU model/good run
      
      q <- q + 1
      r <- r + 1
    }
    m <- m + 1
  }
  i <- i + 1
}

sim_out_data