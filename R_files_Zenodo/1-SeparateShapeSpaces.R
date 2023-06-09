library(geomorph)
library(Morpho)
library(vegan)
library(ggplot2)
library(MASS)
library(shapes)
library(abind)
library(RRPP)

abind<-abind::abind

#reproducibility 
sink('session-info.txt')
devtools::session_info()
sink()
git2r::repository()

# 0.1) Importing landmark data 
###### IF USING LANDMARK FILE SUPPLIED IN DRYAD FILES, SKIP 0.1 AND 2.0-2.1.
###### DOWNLOAD OF .TPS LANDMARK FILE IS AT 2.15.
ldmk.list<-file.path<-list.files(pattern = ".pp", recursive = F)#Reads filepaths for each specimen's landmark file into list
ldmk.array<-lapply(ldmk.list, read.mpp)#Reads landmarks from meshlab pickpoints files into list of PxK matrices (converted to array later)
      # note that these landmark files include landmarks that were not included in analysis

# 0.2) Factors for specs (ID, Habitat, Watershed, etc.)
ldmk.factors<-read.csv("ldmk_factors.csv")#Reads Non-landmark specimen data (ID, population, habitat, watershed) into DF.
ldmk.factors.long<-read.csv("ldmk_factors_ext.csv")#Same as above but double length and with factor for specimen side.

#Reads diet data into DF, and makes specs rownames
diet<-(read.csv("03.diet_CTmorph.csv"))
rownames(diet)<-diet[,1]
diet<-diet[,-1]



#
#   1) Diet
hellingerdiet<-data.frame(diet[1:3],sqrt(diet[,4:45])) #completes hellinger transformation on diet data (already proportional)

# 1.1) dietPCA
diet.pca<-rda(hellingerdiet[4:45], scale = F)
diet.pca.sum<-summary(diet.pca)
diet.pc.scores<-data.frame(diet[1:3],diet.pca.sum$sites)

# 1.2) dietLDA
logdiet.lda<-lda(data = hellingerdiet, hellingerdiet[4:45], grouping = hellingerdiet$habitat)
logdiet.lda.sum<-summary(logdiet.lda)
lda.predictions <- predict(logdiet.lda, diet[4:45])
mean(hellingerdiet$habitat==lda.predictions$class) #proportion of fish correctly assigned to habitat category
diet.ld.scores<-data.frame(diet[1:3],lda.predictions$x) #creates DF with fish data and LD scores

# 1.3) Figure S4
LD1Box<-ggplot(data =  diet.ld.scores, aes(x = watershed, y = LD1, color=habitat)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0))+
  aes(color = habitat) + ylab("Diet LD")+ xlab("Watershed")+
  scale_color_manual(values = c("steelblue","red")) + theme(axis.ticks.x = element_blank())
LD1Box

# 1.4) calculates mean diet diversity (Shannon & Simpson) for each population
diet.pop.BeaL<-subset(diet, pop == "BeaverL")[-c(1:3)] #Beaver Lake
diet.BeaL.tot<-t(as.data.frame(colSums(diet.pop.BeaL)))
diet.pop.BeaS<-subset(diet, pop == "BeaverS")[-c(1:3)] #Beaver Stream
diet.BeaS.tot<-t(as.data.frame(colSums(diet.pop.BeaS)))

diet.pop.BooL<-subset(diet, pop == "BootL")[-c(1:3)] #Boot Lake
diet.BooL.tot<-t(as.data.frame(colSums(diet.pop.BooL)))
diet.pop.BooS<-subset(diet, pop == "BootS")[-c(1:3)] #Boot Stream
diet.BooS.tot<-t(as.data.frame(colSums(diet.pop.BooS)))

diet.pop.MisL<-subset(diet, pop == "MistyL")[-c(1:3)] #Misty Lake
diet.MisL.tot<-t(as.data.frame(colSums(diet.pop.MisL)))
diet.pop.MisS<-subset(diet, pop == "MistyS")[-c(1:3)] #Misty Stream
diet.MisS.tot<-t(as.data.frame(colSums(diet.pop.MisS)))

diet.pop.PyeL<-subset(diet, pop == "PyeL")[-c(1:3)] #Pye Lake
diet.PyeL.tot<-t(as.data.frame(colSums(diet.pop.PyeL)))
diet.pop.PyeS<-subset(diet, pop == "PyeS")[-c(1:3)] #Pye Stream
diet.PyeS.tot<-t(as.data.frame(colSums(diet.pop.PyeS)))

diet.pop.RobL<-subset(diet, pop == "RobertsL")[-c(1:3)] #Roberts Lake
diet.RobL.tot<-t(as.data.frame(colSums(diet.pop.RobL)))
diet.pop.RobS<-subset(diet, pop == "RobertsS")[-c(1:3)] #Roberts Stream
diet.RobS.tot<-t(as.data.frame(colSums(diet.pop.RobS)))

diet.pop.VBL<-subset(diet, pop == "VillageBayL")[-c(1:3)] #Village Bay Lake
diet.VBL.tot<-t(as.data.frame(colSums(diet.pop.VBL)))
diet.pop.VBS<-subset(diet, pop == "VillageBayS")[-c(1:3)] #Village Bay Stream
diet.VBS.tot<-t(as.data.frame(colSums(diet.pop.VBS)))

# 1.5) Creates diet.pop.tot; DF of populations w/ shannon and simpson diversity
diet.pop.tot<-as.data.frame(rbind(diet.BeaL.tot, diet.BeaS.tot, diet.BooL.tot, diet.BooS.tot,
                                  diet.MisL.tot, diet.MisS.tot, diet.PyeL.tot, diet.PyeS.tot, 
                                  diet.RobL.tot, diet.RobS.tot, diet.VBL.tot, diet.VBS.tot))
rownames(diet.pop.tot)<-(diet.pop.tot$pop=c("BeaverL", "BeaverS", "BootL", "BootS",
                                            "MistyL", "MistyS", "PyeL", "PyeS", 
                                            "RobertsL", "RobertsS", "VillageBayL", "VillageBayS"))
diet.pop.tot$shannon<-diversity(diet.pop.tot[,-43], index = "shannon")
diet.pop.tot$simpson<-diversity(diet.pop.tot[,-c(43:45)], index = "simpson")

# 1.6) Creates DIET.TOTAL; DF of individuals used for gut contents w/ diet PCs, LD score, 
    #and population shannon and simpson diversity. Used only for inclusion of diet variables in ldmk.factors
diet.total<-merge(data.frame(diet.pc.scores,"LD1"=diet.ld.scores$LD1),
                  diet.pop.tot[,c("shannon","simpson","pop")], by = "pop")
rownames(diet.total)<-rownames(diet.pc.scores)

# 1.7) Adds diversity, PC1, and LD1 to ldmk.factors
ldmk.factors<-merge(ldmk.factors, diet.pop.tot[c("shannon","simpson","pop")], by.x= "Population", by.y = "pop")
PC1<-aggregate(diet.total$PC1~diet.total$pop, FUN=mean)
names(PC1)<-c("Population","PC1")
LD1<-aggregate(diet.total$LD1~diet.total$pop, FUN=mean)
names(LD1)<-c("Population","LD1")
ldmk.factors<-merge(ldmk.factors, PC1, by = "Population")
ldmk.factors<-merge(ldmk.factors, LD1, by = "Population")


#
#
#
#   2) Creating 3D arrays of specimen landmarks
#
# 2.0) FUNCTION to convert landmarks to 3D array including specimen and landmark IDs
clean.specs<-function(spec.list, spec.names)
{
  names(spec.list)<-spec.names
  spec.array<-simplify2array(spec.list)
  
  rownames(spec.array)<-c("R20","RvOp","ROp","R13","R12","R11","R10","R9","R8","R7",
                          "R6","R5","R3","R2","R1","C14","C19","C.3","C.4","C.5",
                          "L20","LvOp","LOp","L13","L12","L11","L10","L9","L8","L7","L6",
                          "L5","L3","L2","L1","R15","R16","R17","L15","L16","L17","R4",
                          "L4","R18","L18")
  spec.array<-spec.array[c("C14","C19","C.3","C.4","C.5","ROp","RvOp","R1","R2","R3","R4",
                           "R5","R6","R7","R8","R9","R10","R11","R12","R13",
                           "R15","R16","R17","R18","R20","LOp","LvOp","L1","L2","L3","L4",
                           "L5","L6","L7","L8","L9","L10","L11","L12","L13",
                           "L15","L16","L17","L18","L20"),,]
  spec.array
}


# 2.1) Uses clean.specs() to convert full landmark set to 3D array
spec.array<-clean.specs(ldmk.array,ldmk.factors$Spec)
spec.array[spec.array ==0]<-NA #replaces zeros (missing landmarks) w/ NAs#
##
#
writeland.tps(spec.array[c("C14","C19","C.3","C.4","C.5","LvOp","LOp","L12","L11",
                                        "L6","L2","L3","L10","L13","RvOp","ROp","R12","R11",
                                        "R6","R2","R3","R10","R13"),,],"VI_2D3D-land.tps", specID = T)

# 2.15) to import landmark set provided in publication data
spec.array<-readland.tps("VI_2D3D-land.tps",specID="ID") #produces warning, but can be ignored b/c coords already scaled
rownames(spec.array)<-c("C14","C19","C.3","C.4","C.5","LvOp","LOp","L12","L11",
                        "L6","L2","L3","L10","L13","RvOp","ROp","R12","R11",
                        "R6","R2","R3","R10","R13")


# 2.2) Replaces NAs in row 2 of three specimens w/ mean coords of other 
#4 ldmks used to define plane for rotation in lever shape landmarks
spec.array[2,,c("S130100","S130380","S130206")]<-colMeans(spec.array[1:5,,c("S130100","S130380","S130206")], na.rm = T)



#
#
#
#   3) FUNCTION FOR LANDMARK ROTATION
# SEE https://gist.github.com/davidckatz/3812f9270991b09eef3a1cd45913d3e1#file-preferred-rotationI FOR DETAILS
# Used to rotate specimens so that XY plane = mid-sagittal plane.
preferred.rotation <- function(specimen, mids.end, anterior, 
                               posterior, lowest)
{
  # INDEXING
  
  if(is.character(anterior)){anterior <- which(rownames(specimen)==anterior)}
  
  if(is.character(posterior)){posterior <- which(rownames(specimen)==posterior)}
  
  if(is.character(mids.end)){mids.end <- which(rownames(specimen)==mids.end)}
  
  if(is.character(lowest)){lowest <- which(rownames(specimen)==lowest)}
  
  # midline landmark subset
  midlines <- specimen[1:mids.end,]
  #########################################################################
  
  # REORIENTATIONS
  # Identity matrix  
  id.rm <- diag(ncol(specimen))
  
  # ROTATIONS AND CONDITIONAL ROTATIONS
  # 1. primary rotation: minimize z-axis variation for midline landmarks
  minz.rm <- eigen(var(midlines))$vectors
  minz.shape<-specimen%*%minz.rm
  
  # 2. 
  #reflection, reverse sign on z-axis values.
  if(sign(det(minz.rm))==-1) {reflectz.rm <- diag(c(1,1,-1))} else 
  {reflectz.rm<-id.rm}
  reflectz.shape <- minz.shape%*%reflectz.rm
  
  # 3. Ensure X-axis is anatomical antero-posterior.
  antpost.test <- (reflectz.shape[anterior,1]-
                     reflectz.shape[posterior,1])^2 >=
    (reflectz.shape[anterior,2]-reflectz.shape[posterior,2])^2
  
  if(!antpost.test){
    ap.rads <- 90*(pi/180)
    ap.rm <- matrix(c(cos(ap.rads), sin(ap.rads), 0,
                      -sin(ap.rads),cos(ap.rads), 0,
                      0,0,1),3,3)} else {ap.rm <- id.rm}
  ap.shape <- reflectz.shape%*%ap.rm
  
  # 4. Rotate so that ant and post are horizontal to the X-axis
  # compute angle to x-axis (in radians)
  hz.rads <- atan((ap.shape[anterior,2]-ap.shape[posterior,2])/
                    (ap.shape[anterior,1]-ap.shape[posterior,1]))
  # rotation matrix
  hz.rm <- matrix(c(cos(hz.rads), sin(hz.rads), 0,
                    -sin(hz.rads), cos(hz.rads), 0,
                    0,0,1), 3, 3)
  hz.shape <- ap.shape%*%hz.rm
  
  # 5. Ensure up is + and down is - on the Y-axis. 
  # if not, flip.
  updown.test <- hz.shape[lowest,2]<hz.shape[anterior,2]
  
  if(!updown.test){ud.rm<-diag(c(1,-1,-1))} else {ud.rm<-id.rm}
  ud.shape <- hz.shape%*%ud.rm
  
  # COMBINE ROTATION MATRICES INTO A SINGLE MATRIX
  flush.rm <- minz.rm%*%reflectz.rm%*%ap.rm%*%hz.rm%*%ud.rm
  flush<-specimen%*%flush.rm
  
  # check if specimen is pointing toward the right. if not, flip.
  flush[,1]<-
    if((flush[anterior,1]-flush[posterior,1])<0){flush[,1]*(-1)}else{flush[,1]}
  ##########################################################################
  
  # OUTPUT
  list(flush=flush, 
       flush.rm=flush.rm)
}


#   3.1) LOOP ROTATING LEVER ARRAYS
#for specs
rotation.count<-1
repeat {
  spec.array[,,rotation.count]<-preferred.rotation(specimen = spec.array[,,rotation.count], 
                                                   mids.end = 5, anterior = 4, posterior = 3, lowest = 5)$flush
  spec.array[,3,rotation.count]<-
    if(spec.array["LOp",3,rotation.count]-spec.array["ROp",3,rotation.count]>0){spec.array[,3,rotation.count]*(-1)}else{
      spec.array[,3,rotation.count]}
  
  rotation.count<-rotation.count+1
  
  if(rotation.count > 59) {
    break
  }
}

#   3.2) Reduces spec.array to only necessary landmarks
spec.array<-spec.array[c("LvOp","LOp","L12","L11",
                         "L6","L2","L3","L10","L13","C.3","RvOp","ROp","R12","R11",
                         "R6","R2","R3","R10","R13"),,]

#   3.3) Interpolates locations of missing points from other specimens in population
spec.array.boots<-spec.array[,,16:20]
spec.array.boots[c("L2","R2"),,"S130423"]<-NA
spec.array.boots<-estimate.missing(spec.array.boots, method = "TPS")
spec.array.beavers<-spec.array[,,6:10]
spec.array.beavers["ROp",,"S130470"]<-NA
spec.array.beavers<-estimate.missing(spec.array.beavers, method = "TPS")
spec.array.mistys<-spec.array[,,26:30]
spec.array.mistys[c("R2","L2"),,"S130537"]<-NA
spec.array.mistys<-estimate.missing(spec.array.mistys, method = "TPS")
spec.array.vbs<-spec.array[,,55:59]
spec.array.vbs[c("RvOp","R2","L2"),,"S131084"]<-NA #opercle deformity
spec.array.vbs<-estimate.missing(spec.array.vbs, method = "TPS")

spec.array[,,16:20]<-spec.array.boots
spec.array[,,6:10]<-spec.array.beavers
spec.array[,,26:30]<-spec.array.mistys
spec.array[,,55:59]<-spec.array.vbs


#   4)  PREP DATA FOR PROC ANALYSES ON 2D LEVER LANDMARKS
#
# 
spec.array.2d<-spec.array[,-3,] #makES spec.array 2D by dropping Z coords

#splits landmarks into L and R arrays, keeps only landmarks in lever/linkage systems
spec.array.2dL<-as.array(spec.array.2d[c("LvOp","LOp","L12","L11",
                                         "L6","L2","L3","L10"#,"L13","C.3"
),,])
spec.array.2dR<-as.array(spec.array.2d[c("RvOp","ROp","R12","R11",
                                         "R6","R2","R3","R10"#,"R13","C.3"
),,])
rownames(spec.array.2dL)<-c("vOp","Op","12","11","6","2","3","10"#,"13","C"
)
rownames(spec.array.2dR)<-c("vOp","Op","12","11","6","2","3","10"#,"13","C"
)
spec.array.2d.sides<-abind::abind(spec.array.2dL,spec.array.2dR, along = 3) #Binds L and R arrays (along N dimension)

#gpa and bilat symmetry 
proc.2d.sides<-gpagen(spec.array.2d.sides) #Procrustes registration

#plot(proc.2d.sides)
bilat.symmetry(proc.2d.sides, ind = ldmk.factors.long$Spec, side = ldmk.factors.long$Side,
                    object.sym = F, print.progress = T, iter = 999)# Test of bilateral symmetry


#   4.1)  CALCULATES AVERAGE WITHIN-SEPCIMEN LANDMARK COORDINATES FROM BOTH L & R SIDES
spec.array.2d.means<-(spec.array.2dL+spec.array.2dR)/2 #within-spec means
proc.2d.means<-gpagen(spec.array.2d.means) #Procrustes registration
rownames(proc.2d.means$coords)<-rownames(spec.array.2d.means) 
proc.2d.geodf<-geomorph.data.frame(proc.2d.means, pop = ldmk.factors$Population,
                                   watershed = ldmk.factors$Watershed, 
                                   habitat = ldmk.factors$Habitat, shannon = ldmk.factors$shannon,
                                   simpson = ldmk.factors$simpson, dietPC1 = ldmk.factors$PC1,
                                   dietLD1 = ldmk.factors$LD1) #Creates geomorph DF from within-spec mean coords


#   4.2) PREP OF 2D LANDMARKS FOR LINKAGES SCRIPT
# reflection across sag. plane not necessary w/ 2D coords
spec.sym.2d<-abind(spec.array.2d.means, spec.array.2d.means, along = 1) # binds symmetric landmarks from both sides into 1 array
spec.sym.2d.gpa<-gpagen(spec.sym.2d) #Proc. registration
rownames(spec.sym.2d.gpa$coords)<-c("vOp","Op","12","11","6","2","3","10","vOp","Op","12","11","6","2","3","10")
spec.sym.2d.res<-arrayspecs(procD.lm(data=spec.sym.2d.gpa, coords[-(9:16),,]~Csize)$residuals,8,2) #size correction using Csize residuals, 
                          #reduction to one side
rownames(spec.sym.2d.res)<-c("vOp","Op","12","11","6","2","3","10")

Op4.2d.sym<-spec.sym.2d.gpa$coords[1:4,,] #landmark sets for each linkage/lever system
Max4.2d.sym<-spec.sym.2d.gpa$coords[c(4:5,7:8),,]
Mand.2d.sym<-spec.sym.2d.gpa$coords[c(3:6),,]



# 5) creates geomorph df for each 2D lever system.  Proc. Registration performed on dataset with 
      #both sides of each specimen, Procrustes ANOVAs only use one side to avoid overfitting to symmetric landmarks.
#####
# 5.1 - Op4
Op4bar.2d.sym<-spec.sym.2d[c(1:4,9:12),,] #produces subset of landmarks including only Op4 linkages (symmetric)
Op4bar.2d.gpa<-gpagen(Op4bar.2d.sym, print.progress = F) #Procrustes registration
rownames(Op4bar.2d.gpa$consensus)<-rownames(spec.sym.2d[c(1:4,9:12),,]) 
Op4bar.2d.geodf<-geomorph.data.frame(coords =Op4bar.2d.gpa$coords[1:4,,], Csize = Op4bar.2d.gpa$Csize) # makes geomorph dataframe

Op4.2dC.res<-arrayspecs(procD.lm(coords~Csize, data = Op4bar.2d.geodf, logsz = T, print.progress = F)$residuals, 4, 2) #makes residuals of Csize into landmark array
Op4.2dC.res.geodf<-geomorph.data.frame(coords=Op4.2dC.res,pop = ldmk.factors$Population, # makes geomorph dataframe with residuals
                                       watershed = ldmk.factors$Watershed, Csize = Op4bar.2d.gpa$Csize,
                                       habitat = ldmk.factors$Habitat, GCsize = proc.2d.geodf$Csize, shannon = ldmk.factors$shannon,
                                       simpson = ldmk.factors$simpson, dietPC1 = ldmk.factors$PC1,
                                       dietLD1 = ldmk.factors$LD1)

#analysis for effects of shape based on Op4bar centroid size (logged)
summary(procD.lm(coords~watershed*habitat, SS.type = "II", data = Op4.2dC.res.geodf, logsz = T, print.progress = F))
summary(Op4.2dC.WatxLD<-procD.lm(coords~watershed*dietLD1, SS.type = "II", data = Op4.2dC.res.geodf, logsz = T, print.progress = F))
Op4.2dC.PW<-pairwise(Op4.2dC.WatxLD, groups = Op4.2dC.res.geodf$watershed, covariate = Op4.2dC.res.geodf$dietLD1) # Required for vector comparisons 
summary(procD.lm(coords~watershed*dietPC1, SS.type = "II", data = Op4.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*shannon, SS.type = "II", data = Op4.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*simpson, SS.type = "II", data = Op4.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*shannon, SS.type = "II", data = Op4.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*simpson, SS.type = "II", data = Op4.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*dietLD1*simpson, SS.type = "II", data = Op4.2dC.res.geodf, logsz = T, print.progress = F))


# 5.2) same as above, but for mandibular levers 
Dent.2d.sym<-spec.sym.2d[c(3:6,11:14),,]
Dent.2d.gpa<-gpagen(Dent.2d.sym, print.progress = F)
rownames(Dent.2d.gpa$consensus)<-rownames(spec.sym.2d[c(3:6,11:14),,])
Dent.2d.geodf<-geomorph.data.frame(coords=Dent.2d.gpa$coords[1:4,,],Csize = Dent.2d.gpa$Csize)

Dent.2dC.res<-arrayspecs(procD.lm(coords~Csize, data = Dent.2d.geodf, logsz = T, print.progress = F)$residuals, 4, 2)
Dent.2dC.res.geodf<-geomorph.data.frame(coords=Dent.2dC.res,pop = ldmk.factors$Population,
                                        watershed = ldmk.factors$Watershed, Csize = Dent.2d.gpa$Csize,
                                        habitat = ldmk.factors$Habitat, GCsize = proc.2d.geodf$Csize, shannon = ldmk.factors$shannon,
                                        simpson = ldmk.factors$simpson, dietPC1 = ldmk.factors$PC1,
                                        dietLD1 = ldmk.factors$LD1)

#analysis for effects of shape based on Dent centroid size (logged)
summary(procD.lm(coords~watershed*habitat, SS.type = "II",data = Dent.2dC.res.geodf, logsz = T, print.progress = F))
summary(Dent.2dC.WatxLD<-procD.lm(coords~watershed*dietLD1, SS.type = "II", data = Dent.2dC.res.geodf, logsz = T, print.progress = F))
Dent.2dC.PW<-pairwise(Dent.2dC.WatxLD, groups = Dent.2dC.res.geodf$watershed, covariate = Dent.2dC.res.geodf$dietLD1)
summary(procD.lm(coords~watershed*dietPC1, SS.type = "II", data = Dent.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*shannon, SS.type = "II", data = Dent.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*simpson, SS.type = "II", data = Dent.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*shannon, SS.type = "II", data = Dent.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*simpson, SS.type = "II", data = Dent.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*dietLD1*simpson, SS.type = "II", data = Dent.2dC.res.geodf, logsz = T, print.progress = F))


#   5.3) same as above, but for Max4 linkage 
Max4bar.2d.sym<-spec.sym.2d[c(4:5,7:8,12:13,15:16),,]
Max4bar.2d.gpa<-gpagen(Max4bar.2d.sym, print.progress = F)
rownames(Max4bar.2d.gpa$consensus)<-rownames(spec.sym.2d[c(4:5,7:8,12:13,15:16),,])
Max4bar.2d.geodf<-geomorph.data.frame(coords =Max4bar.2d.gpa$coords[c(1:4),,], Csize = Max4bar.2d.gpa$Csize)

#creates geomorph df for Maxillary 4bar
Max4.2dC.res<-arrayspecs(procD.lm(coords~Csize, data = Max4bar.2d.geodf, logsz = T, print.progress = F)$residuals, 4, 2)
Max4.2dC.res.geodf<-geomorph.data.frame(coords=Max4.2dC.res,pop = ldmk.factors$Population,
                                        watershed = ldmk.factors$Watershed, Csize = Max4bar.2d.gpa$Csize,
                                        habitat = ldmk.factors$Habitat, GCsize = proc.2d.geodf$Csize, shannon = ldmk.factors$shannon,
                                        simpson = ldmk.factors$simpson, dietPC1 = ldmk.factors$PC1,
                                        dietLD1 = ldmk.factors$LD1)

#analysis for effects of shape based on Max4bar centroid size (logged)
summary(procD.lm(coords~watershed*habitat, SS.type = "II",data = Max4.2dC.res.geodf, logsz = T, print.progress = F))
summary(Max4.2dC.WatxLD<-procD.lm(coords~watershed*dietLD1, SS.type = "II", data = Max4.2dC.res.geodf, logsz = T, print.progress = F))
Max4.2dC.PW<-pairwise(Max4.2dC.WatxLD, groups = Max4.2dC.res.geodf$watershed, covariate = Max4.2dC.res.geodf$dietLD1)
summary(Max4.2dC.PW, confidence = 0.95,SS.type = "II", test.type = "VC", angle.type = "deg")
summary(procD.lm(coords~watershed*dietPC1, SS.type = "II", data = Max4.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*shannon, SS.type = "II", data = Max4.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*simpson, SS.type = "II", data = Max4.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*shannon, SS.type = "II", data = Max4.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*simpson, SS.type = "II", data = Max4.2dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*dietLD1*simpson, SS.type = "II", data = Max4.2dC.res.geodf, logsz = T, print.progress = F))



#
#
#
#
#    6) 6-6.2 IS SAME AS 4-4.2, BUT FOR 3D DATA

spec.array.3dL<-spec.array[c("LvOp","LOp","L12","L11",
                             "L6","L2","L3","L10"#,"L13","C.3"
),,]
spec.array.3dR<-spec.array[c("RvOp","ROp","R12","R11",
                             "R6","R2","R3","R10"#,"R13","C.3"
),,]
spec.array.3dL[,3,]<-spec.array.3dL[,3,]*(-1) #this reflects 3dL across midline
rownames(spec.array.3dL)<-c("vOp","Op","12","11","6","2","3","10"#,"13","C"
)
rownames(spec.array.3dR)<-c("vOp","Op","12","11","6","2","3","10"#,"13","C"
)
spec.array.3d.sides<-abind::abind(spec.array.3dL,spec.array.3dR, along = 3)

#gpa and bilat symmetry 
proc.3d.sides<-gpagen(spec.array.3d.sides)

bilat.symmetry(proc.3d.sides, ind = ldmk.factors.long$Spec, side = ldmk.factors.long$Side,
               object.sym = F, print.progress = T, iter = 999)

# 6.1) CALCULATES AVERAGE WITHIN-SEPCIMEN LANDMARK COORDINATES FROM BOTH L & R SIDES
spec.array.3d.means<-(spec.array.3dL+spec.array.3dR)/2
proc.3d.means<-gpagen(spec.array.3d.means)
rownames(proc.3d.means$coords)<-rownames(spec.array.3d.means)
proc.3d.geodf<-geomorph.data.frame(proc.3d.means, pop = ldmk.factors$Population,
                                   watershed = ldmk.factors$Watershed, 
                                   habitat = ldmk.factors$Habitat, shannon = ldmk.factors$shannon,
                                   simpson = ldmk.factors$simpson, dietPC1 = ldmk.factors$PC1,
                                   dietLD1 = ldmk.factors$LD1)

# 6.2) 3D prep gpa alignment of landmark coordinates.  Coords represent both sides to maintain orientation relative to sag. plane,
#but sides have been averaged to remove asymmetry.
spec.array.3d.means.flipped<-spec.array.3d.means[c("vOp","Op","12","11","6","2","3","10"),,]
spec.array.3d.means.flipped[,3,]<-spec.array.3d.means[c("vOp","Op","12","11","6","2","3","10"),3,]*(-1)

spec.sym.3d<-abind(spec.array.3d.means[c("vOp","Op","12","11","6","2","3","10"),,], spec.array.3d.means.flipped, along = 1)
spec.sym.3d.gpa<-gpagen(spec.sym.3d)
rownames(spec.sym.3d.gpa$coords)<-c("vOp","Op","12","11","6","2","3","10","vOp","Op","12","11","6","2","3","10")

#size correction
spec.sym.3d.res<-arrayspecs(procD.lm(data = spec.sym.3d.gpa, coords[-(9:16),,]~Csize)$residuals, 8,3)
rownames(spec.sym.3d.res)<-c("vOp","Op","12","11","6","2","3","10")

Op4.3d.sym<-spec.sym.3d.gpa$coords[1:4,,]
Max4.3d.sym<-spec.sym.3d.gpa$coords[c(4:5,7:8),,]
Mand.3d.sym<-spec.sym.3d.gpa$coords[c(3:6),,]

# 7) creates geomorph df for each 3D lever system (same as 5.1-5.3).  Proc. Registration performed on dataset with 
#both sides of each specimen, Procrustes ANOVAs only use one side to avoid overfitting to symmetric landmarks.
#####
# 5.1 - Op4
Op4bar.2d.sym<-spec.sym.2d[c(1:4,9:12),,] #produces subset of landmarks including only Op4 linkages (symmetric)
Op4bar.2d.gpa<-gpagen(Op4bar.2d.sym, print.progress = F) #Procrustes registration
rownames(Op4bar.2d.gpa$consensus)<-rownames(spec.sym.2d[c(1:4,9:12),,]) 
Op4bar.2d.geodf<-geomorph.data.frame(coords =Op4bar.2d.gpa$coords[1:4,,], Csize = Op4bar.2d.gpa$Csize) # makes geomorph dataframe

Op4.2dC.res<-arrayspecs(procD.lm(coords~Csize, data = Op4bar.2d.geodf, logsz = T, print.progress = F)$residuals, 4, 2) #makes residuals of Csize into landmark array
Op4.2dC.res.geodf<-geomorph.data.frame(coords=Op4.2dC.res,pop = ldmk.factors$Population, # makes geomorph dataframe with residuals
                                       watershed = ldmk.factors$Watershed, Csize = Op4bar.2d.gpa$Csize,
                                       habitat = ldmk.factors$Habitat, GCsize = proc.2d.geodf$Csize, shannon = ldmk.factors$shannon,
                                       simpson = ldmk.factors$simpson, dietPC1 = ldmk.factors$PC1,
                                       dietLD1 = ldmk.factors$LD1)

# 7.1) Op4 
#creates geomorph df for opercular 4bar
Op4bar.3d.sym<-spec.sym.3d[c(1:4,9:12),,]
Op4bar.3d.gpa<-gpagen(Op4bar.3d.sym)
###  Next two lines switch Y and Z axes because gpa turned shape on its side#####
Op4bar.3d.gpa$coords[,c(2,3),]<-Op4bar.3d.gpa$coords[,c(3,2),]
Op4bar.3d.gpa$consensus[,c(2,3)]<-Op4bar.3d.gpa$consensus[,c(3,2)]
#
rownames(Op4bar.3d.gpa$consensus)<-rownames(Op4bar.3d.sym)
Op4bar.3d.geodf<-geomorph.data.frame(coords = Op4bar.3d.gpa$coords[1:4,,], Csize = Op4bar.3d.gpa$Csize)

Op4.3dC.res<-arrayspecs(procD.lm(coords~Csize, data = Op4bar.3d.geodf, logsz = T, print.progress = F)$residuals, 4, 3)
Op4.3dC.res.geodf<-geomorph.data.frame(coords=Op4.3dC.res,pop = ldmk.factors$Population,
                                       watershed = ldmk.factors$Watershed, Csize = Op4bar.3d.gpa$Csize,
                                       habitat = ldmk.factors$Habitat, shannon = ldmk.factors$shannon,
                                       simpson = ldmk.factors$simpson, dietPC1 = ldmk.factors$PC1,
                                       dietLD1 = ldmk.factors$LD1)

#analysis for effects of shape based on Op4bar centroid size (logged)
summary(procD.lm(coords~watershed*habitat, SS.type = "II", data = Op4.3dC.res.geodf, logsz = T, print.progress = F))
summary(Op4.3dC.WatxLD<-procD.lm(coords~watershed*dietLD1, SS.type = "II", data = Op4.3dC.res.geodf, logsz = T, print.progress = F))
Op4.3dC.PW<-pairwise(Op4.3dC.WatxLD, groups = Op4.3dC.res.geodf$watershed, covariate = Op4.3dC.res.geodf$dietLD1)
summary(procD.lm(coords~watershed*dietPC1, SS.type = "II", data = Op4.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*shannon, SS.type = "II", data = Op4.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*simpson, SS.type = "II", data = Op4.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*shannon, SS.type = "II", data = Op4.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*simpson, SS.type = "II", data = Op4.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*dietLD1*simpson, SS.type = "II", data = Op4.3dC.res.geodf, logsz = T, print.progress = F))



# 8.2) Mand
#creates geomorph df for Mandibular lever
Dent.3d.sym<-spec.sym.3d[c(3:6,11:14),,]
Dent.3d.gpa<-gpagen(Dent.3d.sym)
###  Switched Y and Z axes because gpa turned shape on its side#####
Dent.3d.gpa$coords[,c(2,3),]<-Dent.3d.gpa$coords[,c(3,2),]*-1
Dent.3d.gpa$consensus[,c(2,3)]<-Dent.3d.gpa$consensus[,c(3,2)]*-1
#
rownames(Dent.3d.gpa$consensus)<-rownames(Dent.3d.sym)
Dent.3d.geodf<-geomorph.data.frame(coords = Dent.3d.gpa$coords[1:4,,], Csize = Dent.3d.gpa$Csize)

Dent.3dC.res<-arrayspecs(procD.lm(coords~Csize, data = Dent.3d.geodf, logsz = T, print.progress = F)$residuals, 4, 3)
Dent.3dC.res.geodf<-geomorph.data.frame(coords=Dent.3dC.res,pop = ldmk.factors$Population,
                                        watershed = ldmk.factors$Watershed, Csize = Dent.3d.gpa$Csize,
                                        habitat = ldmk.factors$Habitat, shannon = ldmk.factors$shannon,
                                        simpson = ldmk.factors$simpson, dietPC1 = ldmk.factors$PC1,
                                        dietLD1 = ldmk.factors$LD1)

#analysis for effects of shape based on Dent centroid size (logged)
summary(procD.lm(coords~watershed*habitat, SS.type = "II",data = Dent.3dC.res.geodf, logsz = T, print.progress = F))
summary(Dent.3dC.WatxLD<-procD.lm(coords~watershed*dietLD1, SS.type = "II", data = Dent.3dC.res.geodf, logsz = T, print.progress = F))
Dent.3dC.PW<-pairwise(Dent.3dC.WatxLD, groups = Dent.3dC.res.geodf$watershed, covariate = Dent.3dC.res.geodf$dietLD1)
summary(procD.lm(coords~watershed*dietPC1, SS.type = "II", data = Dent.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*shannon, SS.type = "II", data = Dent.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*simpson, SS.type = "II", data = Dent.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*shannon, SS.type = "II", data = Dent.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*simpson, SS.type = "II", data = Dent.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*dietLD1*simpson, SS.type = "II", data = Dent.3dC.res.geodf, logsz = T, print.progress = F))



# 8.3) Max4
#creates geomorph df for Maxillary 4bar
Max4bar.3d.sym<-spec.sym.3d[c(4:5,7:8,12:13,15:16),,]
Max4bar.3d.gpa<-gpagen(Max4bar.3d.sym)
###  Switched Y and Z axes because gpa turned shape on its side#####
Max4bar.3d.gpa$coords[,c(2,3),]<-Max4bar.3d.gpa$coords[,c(3,2),]*-1
Max4bar.3d.gpa$consensus[,c(2,3)]<-Max4bar.3d.gpa$consensus[,c(3,2)]*-1
#
rownames(Max4bar.3d.gpa$consensus)<-rownames(Max4bar.3d.sym)
Max4bar.3d.geodf<-geomorph.data.frame(coords=Max4bar.3d.gpa$coords[1:4,,], Csize = Max4bar.3d.gpa$Csize)

Max4.3dC.res<-arrayspecs(procD.lm(coords~Csize, data = Max4bar.3d.geodf, logsz = T, print.progress = F)$residuals, 4, 3)
Max4.3dC.res.geodf<-geomorph.data.frame(coords=Max4.3dC.res,pop = ldmk.factors$Population,
                                        watershed = ldmk.factors$Watershed, Csize = Max4bar.3d.gpa$Csize,
                                        habitat = ldmk.factors$Habitat,shannon = ldmk.factors$shannon,
                                        simpson = ldmk.factors$simpson, dietPC1 = ldmk.factors$PC1,
                                        dietLD1 = ldmk.factors$LD1)

#analysis for effects of shape based on Max4bar centroid size (logged)
summary(procD.lm(coords~watershed*habitat, SS.type = "II",data = Max4.3dC.res.geodf, logsz = T, print.progress = F))
summary(Max4.3dC.WatxLD<-procD.lm(coords~watershed*dietLD1, SS.type = "II", data = Max4.3dC.res.geodf, logsz = T, print.progress = F))
Max4.3dC.PW<-pairwise(Max4.3dC.WatxLD, groups = Max4.3dC.res.geodf$watershed, covariate = Max4.3dC.res.geodf$dietLD1)
summary(procD.lm(coords~watershed*dietPC1, SS.type = "II", data = Max4.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*shannon, SS.type = "II", data = Max4.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*simpson, SS.type = "II", data = Max4.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*shannon, SS.type = "II", data = Max4.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~dietLD1*simpson, SS.type = "II", data = Max4.3dC.res.geodf, logsz = T, print.progress = F))
summary(procD.lm(coords~watershed*dietLD1*simpson, SS.type = "II", data = Max4.3dC.res.geodf, logsz = T, print.progress = F))


#
#
#
#
#
# 9) Correlations between 2D and 3D shape spaces
# 9.1) MANTEL TESTS, correlations between pairwise distance matrices
mantel(Op4bar.2d.gpa$procD,Op4bar.3d.gpa$procD, method = "pearson")
mantel(Dent.2d.gpa$procD,Dent.3d.gpa$procD, method = "pearson")
mantel(Max4bar.2d.gpa$procD,Max4bar.3d.gpa$procD, method = "pearson")


# 9.2)  2D v 3D correlations between distances from shape space centroids
# Op4
Op4.proc.dist3d<-1
count<-1
repeat {
  Op4.proc.dist3d[count]<-procdist(Op4bar.3d.gpa$coords[,,count],Op4bar.3d.gpa$consensus, type = "full")
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}
Op4.proc.dist2d<-1
count<-1
repeat {
  Op4.proc.dist2d[count]<-procdist(Op4bar.2d.gpa$coords[,,count],Op4bar.2d.gpa$consensus, type = "full")
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}
shapiro.test((Op4.proc.dist2d))
shapiro.test((Op4.proc.dist3d))
shapiro.test(log10(Op4.proc.dist2d))
shapiro.test(log10(Op4.proc.dist3d))
cor.test(log10(Op4.proc.dist3d),log10(Op4.proc.dist2d), method = "pearson")
Op4.procdist<-data.frame(procdist.3d = log10(Op4.proc.dist3d), procdist.2d = log10(Op4.proc.dist2d), pop = ldmk.factors$Population,
                         watershed = ldmk.factors$Watershed, habitat = ldmk.factors$Habitat)

# 9.3) 2D v 3D correlations
# Mandibular Levers
Dent.proc.dist3d<-1
count<-1
repeat {
  Dent.proc.dist3d[count]<-procdist(Dent.3d.gpa$coords[,,count],Dent.3d.gpa$consensus, type = "full")
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}
Dent.proc.dist2d<-1
count<-1
repeat {
  Dent.proc.dist2d[count]<-procdist(Dent.2d.gpa$coords[,,count],Dent.2d.gpa$consensus, type = "full")
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}
shapiro.test((Dent.proc.dist2d))
shapiro.test((Dent.proc.dist3d))
shapiro.test(log10(Dent.proc.dist2d))
shapiro.test(log10(Dent.proc.dist3d))
cor.test(log10(Dent.proc.dist3d),log10(Dent.proc.dist2d), method = "pearson")
Dent.procdist<-data.frame(procdist.3d = log10(Dent.proc.dist3d), procdist.2d = log10(Dent.proc.dist2d), pop = ldmk.factors$Population,
                          watershed = ldmk.factors$Watershed, habitat = ldmk.factors$Habitat)

# 9.4) 2D v 3D correlations
# Max4
Max4bar.proc.dist3d<-1
count<-1
repeat {
  Max4bar.proc.dist3d[count]<-procdist(Max4bar.3d.gpa$coords[,,count],Max4bar.3d.gpa$consensus, type = "full")
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}
Max4bar.proc.dist2d<-1
count<-1
repeat {
  Max4bar.proc.dist2d[count]<-procdist(Max4bar.2d.gpa$coords[,,count],Max4bar.2d.gpa$consensus, type = "full")
  
  count<-count+1
  
  if(count > 59) {
    break
  }
}
shapiro.test((Max4bar.proc.dist2d))
shapiro.test((Max4bar.proc.dist3d))
shapiro.test(log10(Max4bar.proc.dist2d))
shapiro.test(log10(Max4bar.proc.dist3d))
cor.test(log10(Max4bar.proc.dist3d),log10(Max4bar.proc.dist2d), method = "pearson")
Max4bar.procdist<-data.frame(procdist.3d = log10(Max4bar.proc.dist3d), procdist.2d = log10(Max4bar.proc.dist2d), pop = ldmk.factors$Population,
                             watershed = ldmk.factors$Watershed, habitat = ldmk.factors$Habitat)

