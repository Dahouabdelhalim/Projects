dyn.load("/usr/local/R/gnu/9.1/3.6.3/site/lib/gdal/2.4.4/lib/libgdal.so")
library(raster)
library(ResistanceGA)
library(geosphere)
library(rgeos)

all_lc<-raster("~/ResistanceGA/both/all_lc_2020.asc")
all_dem<-raster("~/ResistanceGA/both/all_dem_2020.asc")
all_maps<-stack(all_lc,all_dem)

all_response1<-read.csv("~/ResistanceGA/both/all_response_45m_2.csv", row.names = 1)
all_info_2<-read.csv("~/ResistanceGA/both/all_info_2.csv")
all_pairs<-read.csv("~/ResistanceGA/both/all_pairs_45m_2.csv", row.names=1)
all_points<-all_info_2[,c(3,2)]
all_points[,1]<-all_points[,1]+0.0011
all_xy<-SpatialPoints(all_points)
all_response<-lower(all_response1)
geofilter<-lower(all_pairs)
#Test if ResistanceGA will run
GA.inputs <- GA.prep(ASCII.dir =  all_maps,  
                     Results.dir = "/users/PAS1504/martin2537/ResistanceGA/both/multi/45m/", 
                     method = "LL", 
                     seed = 1, 
                     parallel = 2,
                     select.trans=list(NA,"A"))

jl.inputs<-jl.prep(n.Pops = length(all_xy),
                   CS_Point.File = all_xy,
                   response = all_response,
                   pairs_to_include = geofilter)

test.jl <- Run_CS.jl(jl.inputs = jl.inputs,r=all_lc, full.mat = T)
-1 %in% lower(test.jl)[geofilter == 1] 

Single.Surface_optim<-SS_optim(jl.inputs = jl.inputs, GA.inputs = GA.inputs)

Multi.Surface_optim <- MS_optim(jl.inputs = jl.inputs, 
                                GA.inputs = GA.inputs, diagnostic_plots = FALSE)
