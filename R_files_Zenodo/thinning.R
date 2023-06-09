##Thinning species records by distance##################################################################################
# Thinning para Pithecopus ayeaye at KU - 01 Novembro 2018#############################################################

library(knitr)

# Package needed
if(!require(spThin)){
  install.packages("spThin")
  library(spThin)
}


# directory 
setwd("C:/Users/l756n005/Desktop/MODELLING/Thinning/P_ayeaye")

# simple example

## reading data
P_ayeaye<-read.csv("P_ayeaye.csv")


## thinning records
help(thin) # function's help

#Thinning usando NND = 17 Km

thinnedP_ayeaye_a <- thin( loc.data = P_ayeaye, lat.col = "lat", long.col = "long", 
        spec.col = "colecao", thin.par = 17, reps = 50, 
        locs.thinned.list.return = TRUE, write.files = TRUE, 
        max.files = 5, out.dir = "C:/Users/l756n005/Desktop/MODELLING/Thinning/P_ayeaye/results/thin", 
        out.base = "P_ayeaye_a_thined", write.log.file = TRUE,
        log.file = "P_ayeaye_a_thined_full_log_file.txt" )


#Thinning usando NND = 5 Km

thinnedP_ayeaye_a <- thin( loc.data = P_ayeaye, lat.col = "lat", long.col = "long", 
        spec.col = "colecao", thin.par = 5, reps = 50, 
        locs.thinned.list.return = TRUE, write.files = TRUE, 
        max.files = 5, out.dir = "C:/Users/l756n005/Desktop/MODELLING/Thinning/P_ayeaye/results/thin1", 
        out.base = "P_ayeaye_a_thined", write.log.file = TRUE,
        log.file = "P_ayeaye_a_thined_full_log_file.txt" )

#Thinning usando NND = 1 Km

thinnedP_ayeaye_a <- thin( loc.data = P_ayeaye, lat.col = "lat", long.col = "long", 
        spec.col = "colecao", thin.par = 1, reps = 50, 
        locs.thinned.list.return = TRUE, write.files = TRUE, 
        max.files = 5, out.dir = "C:/Users/l756n005/Desktop/MODELLING/Thinning/P_ayeaye/results/thin2", 
        out.base = "P_ayeaye_a_thined", write.log.file = TRUE,
        log.file = "P_ayeaye_a_thined_full_log_file.txt" )


#####
# example considering differences in environmental heterogeneity
# in area of records distribution

## categories of heterogeneity need to be addigned previously

## reading data
tal<-read.csv ("tal.csv")

## thinning records
help(thin)  # function's help

## class 1 of heterogeity distance of thinning 20 km

## thinning records
thinnedtal_a <- thin( loc.data = tal[ which(tal$Territorio == 1) , ], lat.col = "Latitude", 
        long.col = "Longitude", spec.col = "Species", thin.par = 20, reps = 50, 
        locs.thinned.list.return = TRUE, write.files = TRUE, max.files = 5, 
        out.dir = "D:/R/Peltophryne_cu/1a_spThin/TAL1", out.base = "tal_a_thined", 
        write.log.file = TRUE, log.file = "tala_thined_full_log_file.txt" )

## class 2 of heterogeity distance of thinning 10 km

## thinning records
thinnedtal_b <- thin( loc.data = tal[ which(tal$Territorio == 2) , ], lat.col = "Latitude", 
                      long.col = "Longitude", spec.col = "Species", thin.par = 10, reps = 50, 
                      locs.thinned.list.return = TRUE, write.files = TRUE, max.files = 5, 
                      out.dir = "D:/R/Peltophryne_cu/1a_spThin/TAL1", out.base = "tal_b_thined", 
                      write.log.file = TRUE, log.file = "talb_thined_full_log_file.txt" )
