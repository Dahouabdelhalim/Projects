#gdm
library(gdm)

#load Fst distance matrix
ref<-read.table("pw_all_V10_fst.txt", header = TRUE)
j8105.mt<-read.table("pw_j8105_fst.txt", header = TRUE)

#Load table with Env variables and spatial coordinates
sp.env<-read.table("Pop_GeoEnv_gdm.txt", header = TRUE)
head(sp.env)
env.MT <- sp.env[,c(1:3,4)] # reference
env.PR <- sp.env[,c(1:3,5)]
env.AI <- sp.env[,c(1:3,6)]
head(env.MT)
?formatsitepair
sp.ref<-formatsitepair(ref,3,dist = "bray", 
                    XColumn="X", 
                    YColumn="Y", 
                    weightType = "equal", 
                    predData = env.MT,
                    siteColumn = "Site")



ref.gdm<-gdm(sp.ref, geo = TRUE)
#for custom plotting
write.table(ref.gdm, "ref.gdm.txt", quote=F)
#extract deviance
ref.gdm$explained

plot(ref.gdm)

sp.j8105.mt<-formatsitepair(j2510.mt,3,dist = "bray", 
                            XColumn="X", 
                            YColumn="Y", 
                            weightType = "equal", 
                            predData = env.MT,
                            siteColumn = "Site")

j8105.mt.gdm<-gdm(sp.j8105.mt, geo = TRUE)
#for custom plotting
write.table(j8105.mt.gdm, "spline_j8105.txt", quote=F)
#extract deviance
j8105.mt.gdm$explained
plot(j8105.mt.gdm)

