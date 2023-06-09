### Code to run primary trend analyses in Knapp et al. 2016
### models were fit using the package 'INLA' version 0.0-1420281647

dat <- read.csv("Data.csv")  ### set the working directory to the location of "Data.csv" and "W.dat"
library(reshape)
library(INLA)

### main trend model for adults


cDat <- data.frame(y = dat$adult, year = dat$year,year1 = dat$year,date = dat$date, cell = dat$cell, cell2 = dat$cell,site = dat$site)

zoutfullRE <- inla(y ~ year1+date+ f(year,model='ar1',replicate = site) + f(site, model = 'iid') + f(cell,model='besag', graph.file='W.dat',replicate=year),
                   family="zeroinflatednbinomial1",data=cDat,control.compute=list(dic=TRUE),  verbose = TRUE)

### trends models for subadults and tadpoles

cDat <- data.frame(y = dat$adult, year = dat$year-1992,year1 = dat$year-1992,date = dat$date, cell = dat$cell, cell2 = dat$cell,site = dat$site)

zoutfullREsubadult <- inla(y ~ year1+date+ f(year,model='ar1',replicate = site)+ f(site, model = 'iid') + f(cell,model='besag', graph.file='W.dat',replicate=year) ,
                           family="zeroinflatednbinomial1",data=cDat, verbose= TRUE)

cDat <- data.frame(y = dat$adult, year = dat$year,year1 = dat$year,date = dat$date, cell = dat$cell, cell2 = dat$cell,site = dat$site)
zoutfullREtadpole <- inla(y ~ year1+ date+f(year,model='ar1',replicate = site)+ f(site, model = 'iid') + f(cell,model='besag', graph.file='W.dat',replicate=year) ,
                          family="zeroinflatednbinomial1",data=cDat, verbose= TRUE)

### covariate model

year <- dat$year-1992
fish <- as.numeric(dat$fish!="No")
cDat <- data.frame(y = dat$adult, water = dat$w, year = year , year1 = year,date = dat$date, cell = dat$cell ,cell2 = dat$cell ,site = dat$site,fish = fish, depth = dat$depth, elev = dat$elev, basin = dat$basin, basinY = dat$basin*year, fishY = fish*year, depthY = dat$depth*year,elevY = dat$elev*year)
outCOV <- inla(y ~ year1 + date + water + fish*year1 + depth + elev + basin + depthY + elevY + basinY + f(year,model='ar1',replicate = cell2) + f(site, model = 'iid') + f(cell,model='besag', graph.file='W.dat',replicate=year),family="zeroinflatednbinomial1",data=cDat,  verbose = TRUE)


#### spatial trend model

cDat <- data.frame(y = dat[,"adult"],  year = year , year1 = year, year2 = year,date = dat$date, cell = dat$cell ,site = dat$site)

outSpatialTrend <- inla(y ~ year1 + date  + f(year,model='ar1') + f(site, model = 'iid') + 
                          f(cell,year2,model='besag', graph.file='W.dat'),family="zeroinflatednbinomial1",
                        data=cDat, verbose = TRUE)


