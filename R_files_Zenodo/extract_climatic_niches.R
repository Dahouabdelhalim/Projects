# Ignacio Quintero
# 28 04 2016

# t(-_-t)
#===============================================================================


################################################################################
# Supplementary Electronic Material for:
# Quintero, I., Suchard, M., Jetz, W. (2022) Macroevolutionary dynamics of climatic niche space. Proceedings of the Royal Society B. doi: 10.1098/rspb.2022.0091


# code to extract the climatic niches of birds and create polygons using 
# kernel density estimation for input to beast analysis.
################################################################################


################################################################################
# The following publicly available data sets are required:
# new_2002           = New, Mark, David Lister, Mike Hulme, and Ian Makin. 2002. “A High-Resolution Data Set of Surface Climate over Global Land Areas.” Climate Research 21 (1): 1–25.
# weigelt_2013       = Weigelt, Patrick, Walter Jetz, and Holger Kreft. 2013. “Bioclimatic and Physical Characterization of the World’s Islands.” Proceedings of the National Academy of Sciences 110 (38): 15307–12.
# bird_distributions = Hurlbert, Allen H., and Walter Jetz. 2007. “Species Richness, Hotspots, and the Scale Dependence of Range Maps in Ecology and Conservation.” Proceedings of the National Academy of Sciences 104 (33): 13384–89. Available through Map of life: www.mol.org
# BirdLife_maps      = BirdLife International and Handbook of the Birds of the World (2015) Bird species distribution maps of the world. Version 7.0. Available at http://datazone.birdlife.org/species/requestdis.
# bird_elevations    = Quintero, Ignacio, and Walter Jetz. 2018. “Global Elevational Diversity and Diversification Birds.” Nature 555: 246–50.
# mountains          = Quintero, Ignacio, and Walter Jetz. 2018. “Global Elevational Diversity and Diversification Birds.” Nature 555: 246–50.
# south_america      = Global Administrative Areas ( 2012). GADM database of Global Administrative Areas, version 2.0. [online] URL: www.gadm.org.
# chelsa             = Karger, Dirk Nikolaus, Olaf Conrad, Jürgen Böhner, Tobias Kawohl, Holger Kreft, Rodrigo Wilber Soria-Auza, Niklaus E Zimmermann, H Peter Linder, and Michael Kessler. 2017. “Climatologies at High Resolution for the Earth’s Land Surface Areas.” Scientific Data 4: 170122.
# dem                = Danielson, J. J., & Gesch, D. B. (2011). Global multi-resolution terrain elevation data 2010 (GMTED2010) (p. 26). Washington, DC, USA: US Department of the Interior, US Geological Survey. https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php
#
# in the script, these data sets are marked as `<dataset>`
################################################################################

# load packages
library(data.table)
library(maptools)
library(rgeos)
library(rgdal)
library(raster)

# read data from New et al, 2002 
tmp = fread(<new_2002 with temperature data>)
pre = fread(<new_2002 with precipitation data>)
ele = fread(<new_2002 with elevation data>)

# set names
setnames(tmp,c('lat', 'lon', 'm1_tmp', 'm2_tmp', 'm3_tmp', 'm4_tmp', 'm5_tmp',
				'm6_tmp', 'm7_tmp', 'm8_tmp', 'm9_tmp', 'm10_tmp', 'm11_tmp', 'm12_tmp'))
colnames(pre) = c('lat', 'lon', 'm1_pre', 'm2_pre', 'm3_pre', 'm4_pre', 'm5_pre',
				'm6_pre', 'm7_pre', 'm8_pre', 'm9_pre', 'm10_pre', 'm11_pre', 'm12_pre',
				'm1_pre_cv', 'm2_pre_cv', 'm3_pre_cv', 'm4_pre_cv', 'm5_pre_cv', 'm6_pre_cv', 
				'm7_pre_cv', 'm8_pre_cv', 'm9_pre_cv', 'm10_pre_cv', 'm11_pre_cv', 'm12_pre_cv')
colnames(ele) = c('lat','lon','elevation')


elpre = merge(ele,pre,by=c('lat','lon'))
telpre = merge(tmp,elpre,by=c('lat','lon'))

dt_sp = SpatialPointsDataFrame(telpre[,2:1,with=F], data = telpre[,3:27,with=F])
proj4string(dt_sp) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

ras = raster(ext=extent(dt_sp),resolution=c(1/6,1/5.9999))

rs = rasterize(dt_sp,ras)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# load bird ranges
br = readShapePoly(<bird_distributions>)
# 3: non-breeding only
# 1: breeding only
# 2: resident year-round (including breeding)
# 4: 1 or 2 – not distinguished

# remove non-breeding ranges
br = br[!br@data$OccCode == 3,]
proj4string(br) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# divide migrant species
# which residents
wres = br[br@data$OccCode == 2 | br@data$OccCode == 4,]

# extract all months for resident species
ores = extract(rs, wres, small = TRUE, cellnumbers = TRUE)

# organize data
oresdt = lapply(ores, as.data.frame)
dtr    = rbindlist(oresdt)

# make lat long from cells ID
lalon = xyFromCell(rs, dtr[,cell])
lalon = as.data.frame(lalon)
setDT(lalon)
setnames(lalon,c('lon', 'lat'))

# nrows per species
pcs = sapply(oresdt, nrow)

# insert species names
nams = data.table(spid=rep(wres@data$SpecID,pcs),sp=rep(wres@data$Latin,pcs))
dtr  = cbind(nams,dtr,lalon)
dres = na.omit(dtr)

# estimate average temperature and precipitation
dres[,avg_tmp := rowMeans(.SD),.SDcols=c(colnames(dres)[5:16])]
dres[,avg_pre := rowMeans(.SD),.SDcols=c(colnames(dres)[18:29])]

dres = dres[,.(spid, sp, lon, lat, elevation, avg_tmp, avg_pre)]
dres[,elevation := elevation*1000]


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# which migrant
wmig = br[br@data$OccCode == 1,]
omig = extract(rs, wmig, small = TRUE, cellnumbers = TRUE)

# organize data
omigdt = lapply(omig,as.data.frame)
dtm    = rbindlist(omigdt)

# make lat long from cells ID
lalon = xyFromCell(rs, dtm[,cell])
lalon = as.data.frame(lalon)
setDT(lalon)
setnames(lalon,c('lon', 'lat'))

# nrows per species
pcs    = sapply(omigdt, nrow)

nams = data.table(spid=rep(wmig@data$SpecID,pcs),sp=rep(wmig@data$Latin,pcs))
dt0  = cbind(nams,dtm,lalon)
dmig = na.omit(dt0)

# decide which species are south and which north
dmig[,migloc := NA_character_]

for (i in 1:length(wmig)){

	ext = extent(wmig[i,])
	avgext = (ext@ymin + ext@ymax)/2

	if (avgext > 0.0) {
		dmig[spid == wmig[i,]@data$SpecID, migloc := 'north']
	}	else {
		if (avgext < 0.0) {
			dmig[spid == wmig[i,]@data$SpecID, migloc := 'south']
		} else {
			print(paste("error in", wmig[i,]@data$SpecID,sep = ''))
		}
	}
	cat(i, '\\n')
}

dmig[,avg_tmp := NA_real_]
dmig[,avg_pre := NA_real_]

# for northerns
dmig[migloc == 'north',
		 avg_tmp := rowMeans(.SD),.SDcols=c(colnames(dmig)[8:13])]
dmig[migloc == 'north',
		 avg_pre := rowMeans(.SD),.SDcols=c(colnames(dmig)[21:26])]

# for southerns
dmig[migloc == 'south',
		 avg_tmp := rowMeans(.SD),.SDcols=c(colnames(dmig)[c(5:7,14:16)])]
dmig[migloc == 'south',
		 avg_pre := rowMeans(.SD),.SDcols=c(colnames(dmig)[c(18:20,27:29)])]


dmig = dmig[,.(spid,sp,lon,lat,elevation,avg_tmp,avg_pre)]
dmig[,elevation := elevation*1000]


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# join residents with migratory
dc = rbind(dres, dmig)
setorder(dc, cols = 'spid')


###########################
### 156 Remaining species

# species lacking environmental data
asp  = unique(br@data$SpecID)
ssp  = unique(dc[,spid])
nspp = asp[is.na(match(asp,ssp))]

# ranges of species lacking data
nbb = br[!is.na(match(br@data$SpecID,nspp)),]

# read Weigelt data: 
wd = fread(<weigelt_2013>)
wd = wd[,.(Long,Lat,Elev,Temp,Prec)]
wdpp = SpatialPointsDataFrame(wd[,.(Long,Lat)],data=wd[,1:5,with=F])
proj4string(wdpp) = proj4string(br)

ovr = over(nbb, wdpp, returnList=T)

names(ovr) = nbb@data$SpecID

dt = rbindlist(ovr, use.names=T)

pcs = sapply(ovr, nrow)
nams = data.table(spid=rep(nbb@data$SpecID,pcs),sp=rep(nbb@data$Latin,pcs))

dt0 = cbind(nams,dt)

setnames(dt0, c('spid','sp','lon','lat','elevation','mean_tmp','ann_pre'))
dt0[, avg_pre := ann_pre/12]
setnames(dt0, 'mean_tmp', 'avg_tmp')
dt0[,ann_pre := NULL]

dc = rbind(dc, dt0)

# length(unique(dc[,spid]))
# number of species 9955!


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Still 38 spp to go.

# species lacking environmental data
asp  = unique(br@data$SpecID)
ssp  = unique(dc[,spid])
nspp = asp[is.na(match(asp,ssp))]

# ranges of species lacking data
nbb = br[!is.na(match(br@data$SpecID,nspp)),]

# get centroids
cent = gCentroid(nbb, byid=T)

library(fossil)
cc = coordinates(cent)
cwd = coordinates(wdpp)

ndt = data.table("Long", "Lat", "Elev", "Temp", "Prec", "spid", "sp")
setnames(ndt,c("Long", "Lat", "Elev", "Temp", "Prec", "spid", "sp"))

# get environmental data from closest island
for (j in 1:dim(cc)[1]) {
	dds = numeric(dim(cwd)[1])
	for (i in 1:dim(cwd)[1]) dds[i] = earth.dist(rbind(cc[j,],cwd[i,]))
	xx = wd[which.min(dds),]
	xx[,spid:=nbb@data[j,'SpecID']][,sp:=nbb@data[j,'Latin']]
	ndt = rbind(ndt,xx)
	cat(j,'\\n')
}

ndt=ndt[-1,]

setnames(ndt,c('lon','lat','elevation','avg_tmp','ann_pre','spid','sp'))
ndt[,avg_pre := as.numeric(ann_pre)/12]

ndt = ndt[,c(6,7,1:4,8),with=F]

dc = rbind(dc,ndt)

# length(unique(dc[,spid]))
# number of species 9993!

#remove duplicated
dc = dc[!duplicated(dc)]
dc[,spid:=as.integer(spid)][,sp := as.character(sp)][,lon:= as.double(lon)][,lat:= as.double(lat)][,elevation:=as.double(elevation)][,avg_tmp := as.double(avg_tmp)]

# DONE!







#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# refine using elevational ranges
ele = fread(<bird_elevations>)

# remove duplicated ranges
ele = unique(ele,by=c('SpecID','elev_min','elev_max'))

#first, just join unique ranges
dspp = unique(ele[duplicated(ele[,SpecID])][,SpecID])

u.ele = ele[is.na(match(ele[,SpecID],dspp))] # unique species elevations
u.ele = u.ele[,.(SpecID,elev_min,elev_max)]
setnames(u.ele,c('spid','el_mn','el_mx'))

u.dat = dc[!is.na(match(dc[,spid],u.ele[,spid]))]

u.dt = merge(u.dat, u.ele, by='spid')

#filter
nsp0 = unique(u.dt[,spid])
fudt = u.dt[elevation<=el_mx & elevation>=el_mn]
nsp1 = unique(fudt[,spid])
nmsp = nsp0[is.na(match(nsp0,nsp1))]

# length(nmsp)
# 86 species don't match
### extract from CHELSA

# altitude from GMTED2010
alt = raster(<dem>)
# from CHELSA
tmp = raster(<chelsa_bio10_01.tif>)
pre = raster(<chelsa_bio10_12.tif>)

rs = stack(alt,tmp,pre)


# get range of terrestrial values
NAvalue(tmp) = 32768
mx = cellStats(tmp, stat = 'max', na.rm = TRUE)
NAvalue(tmp) = -32768
mn = cellStats(tmp, stat = 'min', na.rm = TRUE)


# loop for all species to get climatic information
csp = br[br@data$SpecID==nmsp[1],]
edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
setDT(edt)
xy = xyFromCell(rs, edt[,cell])
edt[,lon := xy[,'x']]
edt[,lat := xy[,'y']]
edt[,cell := NULL]
setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
edt[,avg_tmp:=avg_tmp/10]
edt[,avg_pre:=ann_pre/12]
c.spp = cbind(spid=nmsp[1],sp=ele[SpecID==nmsp[1],Sp], edt)
c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
c.spp = merge(c.spp,u.ele[spid==nmsp[1]],by='spid')
c.spp = c.spp[elevation<=el_mx & elevation>=el_mn]

r0 = c.spp

for (i in 2:length(nmsp)) {
	csp = br[br@data$SpecID==nmsp[i],]
	if (any(csp@data$OccCode == 1)) {
		print(nmsp[i])
		if (any(csp@data$OccCode == 2 | csp@data$OccCode == 4)) {
			print("saved")
			csp = csp[!csp@data$OccCode == 1,]
		}
	}
	edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
	setDT(edt)
	xy = xyFromCell(rs, edt[,cell])
	edt[,lon := xy[,'x']]
	edt[,lat := xy[,'y']]
	edt[,cell := NULL]
	setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
	edt[,avg_tmp:=avg_tmp/10]
	edt[,avg_pre:=ann_pre/12]
	c.spp = cbind(spid=nmsp[i],sp=ele[SpecID==nmsp[i],Sp], edt)
	c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
	c.spp = merge(c.spp,u.ele[spid==nmsp[i]],by='spid')
	c.spp = c.spp[elevation<=el_mx & elevation>=el_mn]

	if(nrow(c.spp) != 0) r0 = rbind(r0,c.spp)
	
	cat(i,'\\n')
}



# length(unique(r0[,spid]))
# 10 still don't match
u.dt1 = rbind(fudt,r0)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# which still doesn't match, solve manually

nsp0 = unique(u.dt[,spid])
nsp1 = unique(u.dt1[,spid])
nmsp = nsp0[is.na(match(nsp0,nsp1))]

### length(nmsp)
# 10 species without data still

	# "Tragopan melanocephalus"
	i = 1
	csp = readShapePoly(<BirdLife_maps/Tragopan_melanocephalus/species_22679147.shp>)
	edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
	setDT(edt)
	xy = xyFromCell(rs, edt[,cell])
	edt[,lon := xy[,'x']]
	edt[,lat := xy[,'y']]
	edt[,cell := NULL]
	setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
	edt[,avg_tmp:=avg_tmp/10]
	edt[,avg_pre:=ann_pre/12]
	c.spp = cbind(spid=nmsp[i],sp=ele[SpecID==nmsp[i],Sp], edt)
	c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
	c.spp = merge(c.spp,u.ele[spid==nmsp[i]],by='spid')
	c.spp = c.spp[elevation<=el_mx & elevation>=el_mn]
	u.dt1 = rbind(u.dt1, c.spp)

	# "Eriocnemis mirabilis"
	i = 2
	csp = readShapePoly(<BirdLife_maps/Eriocnemis_mirabilis/species_22687939.shp>)
	edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
	setDT(edt)
	xy = xyFromCell(rs, edt[,cell])
	edt[,lon := xy[,'x']]
	edt[,lat := xy[,'y']]
	edt[,cell := NULL]
	setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
	edt[,avg_tmp:=avg_tmp/10]
	edt[,avg_pre:=ann_pre/12]
	c.spp = cbind(spid=nmsp[i],sp=ele[SpecID==nmsp[i],Sp], edt)
	c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
	c.spp = merge(c.spp,u.ele[spid==nmsp[i]],by='spid')
	c.spp = c.spp[elevation<=el_mx & elevation>=el_mn]
	u.dt1 = rbind(u.dt1, c.spp)

	# "Calidris melanotos" (elevational range is wrong!)
	i = 3
	csp = br[br@data$SpecID==nmsp[i],]
	c.spp = u.dt[spid==nmsp[i]]
	u.dt1 = rbind(u.dt1,c.spp)

	# "Pterodroma hasitata" (elevational range is wrong!)
	i = 4
	csp = br[br@data$SpecID==nmsp[i],]
	edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
	setDT(edt)
	xy = xyFromCell(rs, edt[,cell])
	edt[,lon := xy[,'x']]
	edt[,lat := xy[,'y']]
	edt[,cell := NULL]
	setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
	edt[,avg_tmp:=avg_tmp/10]
	edt[,avg_pre:=ann_pre/12]
	c.spp = cbind(spid=nmsp[i],sp=ele[SpecID==nmsp[i],Sp], edt)
	c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
	c.spp = na.omit(c.spp)
	c.spp = merge(c.spp,u.ele[spid==nmsp[i]],by='spid')
	u.dt1 = rbind(u.dt1,c.spp)

	# "Nectarinia ludovicensis"
	i = 5
	csp = br[br@data$SpecID==nmsp[i],]
	edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
	setDT(edt)
	xy = xyFromCell(rs, edt[,cell])
	edt[,lon := xy[,'x']]
	edt[,lat := xy[,'y']]
	edt[,cell := NULL]
	setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
	edt[,avg_tmp:=avg_tmp/10]
	edt[,avg_pre:=ann_pre/12]
	c.spp = cbind(spid=nmsp[i],sp=ele[SpecID==nmsp[i],Sp], edt)
	c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
	c.spp = merge(c.spp,u.ele[spid==nmsp[i]],by='spid')
	# get data over 2200 (assuming an error of 300 m in ele database)
	c.spp = c.spp[elevation>=2200]

	u.dt1 = rbind(u.dt1,c.spp)


	## Nectarinia loveridgei 
#	i = 6
	spidi = 11674
	csp = readShapePoly(<BirdLife_maps/Nectarinia_loveridgei/species_22717931.shp>)
	edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
	setDT(edt)
	xy = xyFromCell(rs, edt[,cell])
	edt[,lon := xy[,'x']]
	edt[,lat := xy[,'y']]
	edt[,cell := NULL]
	setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
	edt[,avg_tmp:=avg_tmp/10]
	edt[,avg_pre:=ann_pre/12]
	c.spp = cbind(spid=spidi,sp=ele[SpecID==spidi,Sp], edt)
	c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
	c.spp = merge(c.spp,u.ele[spid==spidi],by='spid')
	c.spp = c.spp[elevation<=el_mx & elevation>=el_mn]
	u.dt1 = rbind(u.dt1,c.spp)

	## Lagonosticta virata
	#i = 7
	spidi = 12097 
	csp = readShapePoly(<BirdLife_maps/Lagonosticta_virata/species_22719461.shp>)
	edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
	setDT(edt)
	xy = xyFromCell(rs, edt[,cell])
	edt[,lon := xy[,'x']]
	edt[,lat := xy[,'y']]
	edt[,cell := NULL]
	setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
	edt[,avg_tmp:=avg_tmp/10]
	edt[,avg_pre:=ann_pre/12]
	c.spp = cbind(spid=spidi,sp=ele[SpecID==spidi,Sp], edt)
	c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
	c.spp = merge(c.spp,u.ele[spid==spidi],by='spid')
	# closest to elevation
	c.spp = c.spp[elevation>= 800]
	u.dt1 = rbind(u.dt1,c.spp)

	## Liocichla bugunorum
	i = 8
	c.spp = u.dt[spid==nmsp[i]][2,]
	u.dt1 = rbind(u.dt1, c.spp)

	## Pterodroma caribbaea (elevational range might be wrong!)
	i = 9
	c.spp = u.dt[spid==nmsp[i]]
	u.dt1 = rbind(u.dt1,c.spp)

	## Leucosticte sillemi (elevational range is wrong!)
	i = 10
	c.spp = u.dt[spid==nmsp[i]]
	u.dt1 = rbind(u.dt1,c.spp)

# done!
nsp0 = unique(u.dt[,spid])
nsp1 = unique(u.dt1[,spid])
nmsp = nsp0[is.na(match(nsp0,nsp1))]
length(nmsp)





#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
##################
# second, join species without elevation data
##################

usp    = unique(dc[,spid])
usp.el = unique(ele[,SpecID])

#which spp don't have elevational data
wna = usp[is.na(match(usp,usp.el))]
u.dt2 = dc[!is.na(match(dc[,spid],wna)),]
u.dt12 = rbind(u.dt1,u.dt2,fill = TRUE)





#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
##################
# third, approximate elevations with long/lat
# for species with multiple ranges
##################

usp = unique(u.dt12[,spid])
mdt = dc[is.na(match(dc[,spid],usp)),]

# read mountains
mon = readShapePoly(<mountains>)
proj4string(mon) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

library(rgeos)
mcent = coordinates(gCentroid(mon,byid=TRUE))

mcent = cbind(mon@data[,c('Id','Name')],mcent)
setDT(mcent)
setnames(mcent,c('moId','name','lon.M','lat.M'))

ele = merge(ele,mcent,by='moId')
setnames(ele,'SpecID','spid')

################
# add Andes Centroids

#read countries for andes

sa_files = list.files(<south_america>)

load(sa_files[1])
arg = gadm
load(sa_files[2])
bol = gadm
load(sa_files[4])
chi = gadm
load(sa_files[5])
col = gadm
load(sa_files[6])
ecu = gadm
load(sa_files[10])
per = gadm
load(sa_files[14])
ven = gadm

c2 = coordinates(gCentroid(arg))
ele[country=='ARG',lat.M:=c2[2]][country=='ARG',lon.M:=c2[1]]
c2 = coordinates(gCentroid(bol))
ele[country=='BOL',lat.M:=c2[2]][country=='BOL',lon.M:=c2[1]]
c2 = coordinates(gCentroid(chi))
ele[country=='CHL',lat.M:=c2[2]][country=='CHL',lon.M:=c2[1]]
c2 = coordinates(gCentroid(col))
ele[country=='COL',lat.M:=c2[2]][country=='COL',lon.M:=c2[1]]
c2 = coordinates(gCentroid(ecu))
ele[country=='ECU',lat.M:=c2[2]][country=='ECU',lon.M:=c2[1]]
c2 = coordinates(gCentroid(per))
ele[country=='PER',lat.M:=c2[2]][country=='PER',lon.M:=c2[1]]
c2 = coordinates(gCentroid(ven))
ele[country=='VEN',lat.M:=c2[2]][country=='VEN',lon.M:=c2[1]]



#############################
# create function to find
# best mountain ele ranges
# according to distance to 
# mountain centroids

usp = unique(mdt[,spid])

mdt[,poid:=seq(.N)][,el_mn := 0L][,el_mx := 0L]
setkeyv(mdt,c('spid','poid'))
setkeyv(ele,'spid')


#super fast loop thanks to data.table
n=0
for (i in usp) {
	d1 = mdt[J(i)]
	d2 = ele[J(i)]
	dd = merge(d1,d2,by='spid',allow.cartesian=TRUE)
	w.mon = dd[,which.min(sqrt((lon-lon.M)^2 + (lat-lat.M)^2)), by=.(poid)]
	r_el = d2[w.mon[,V1],.(elev_min,elev_max)]
	mdt[J(i),':='(el_mn=as.integer(r_el[,elev_min]),el_mx=as.integer(r_el[,elev_max]))]
	n = n+1
	cat(n,'\\n')
}

#filter by elevation
mdt[J(10774),el_mn:=0L]
fmdt = mdt[elevation<=el_mx & elevation>=el_mn]

# which species are left without data?
nsp0 = unique(mdt[,spid])
nsp1 = unique(fmdt[,spid])
nmsp = nsp0[is.na(match(nsp0,nsp1))]

# 4 species without data
length(nmsp)


# loop for all species to get climatic information
csp = br[br@data$SpecID==nmsp[1],]
edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
setDT(edt)
xy = xyFromCell(rs, edt[,cell])
edt[,':='(lon = xy[,'x'], lat = xy[,'y'])]
edt[,cell := NULL]
setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
edt[,':='(avg_tmp=avg_tmp/10, avg_pre=ann_pre/12)]
c.spp = cbind(spid=nmsp[1],sp=ele[spid==nmsp[1],Sp], edt)
c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
c.spp[,poid:=seq(.N)]
d2 = ele[J(nmsp[1])]
dd = merge(c.spp,d2,by='spid',allow.cartesian=TRUE)
w.mon = dd[,which.min(sqrt((lon-lon.M)^2 + (lat-lat.M)^2)), by=.(poid)]
r_el = d2[w.mon[,V1],.(elev_min,elev_max)]
c.spp[,':='(el_mn=as.integer(r_el[,elev_min]),el_mx=as.integer(r_el[,elev_max]))]
c.spp = c.spp[elevation<=el_mx & elevation>=el_mn]

r0 = c.spp

for (i in 2:length(nmsp)) {
	csp = br[br@data$SpecID==nmsp[i],]
	edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
	setDT(edt)
	xy = xyFromCell(rs, edt[,cell])
	edt[,':='(lon = xy[,'x'], lat = xy[,'y'])]
	edt[,cell := NULL]
	setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
	edt[,':='(avg_tmp=avg_tmp/10, avg_pre=ann_pre/12)]
	c.spp = cbind(spid=nmsp[i],sp=ele[spid==nmsp[i],Sp], edt)
	c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
	c.spp[,poid:=seq(.N)]
	d2 = ele[J(nmsp[i])]
	dd = merge(c.spp,d2,by='spid',allow.cartesian=TRUE)
	w.mon = dd[,which.min(sqrt((lon-lon.M)^2 + (lat-lat.M)^2)), by=.(poid)]
	r_el = d2[w.mon[,V1],.(elev_min,elev_max)]
	c.spp[,':='(el_mn=as.integer(r_el[,elev_min]),el_mx=as.integer(r_el[,elev_max]))]
	c.spp = c.spp[elevation<=el_mx & elevation>=el_mn]

	if(nrow(c.spp)!=0) r0 = rbind(r0,c.spp)

	cat(i,'\\n')
}

fmdt[,poid:=NULL]
r0[,poid:=NULL]

fmdt1 = rbind(fmdt,r0)

nsp0 = unique(mdt[,spid])
nsp1 = unique(fmdt1[,spid])
nmsp = nsp0[is.na(match(nsp0,nsp1))]


# ONLY one species remaining: check manually
i = 1

csp = readShapePoly(<BirdLife_maps/Artisornis_moreaui/species_103771879.shp>)
edt = as.data.frame(extract(rs,csp,cellnumbers=TRUE,small=TRUE)[[1]])
setDT(edt)
xy = xyFromCell(rs, edt[,cell])
edt[,':='(lon = xy[,'x'], lat = xy[,'y'])]
edt[,cell := NULL]
setnames(edt,c('elevation','avg_tmp','ann_pre', 'lon', 'lat'))
edt[,':='(avg_tmp=avg_tmp/10, avg_pre=ann_pre/12)]
c.spp = cbind(spid=nmsp[i],sp=ele[spid==nmsp[i],Sp], edt)
c.spp = c.spp[,c(1:2,6,7,3:4,8),with=F]
c.spp[,poid:=seq(.N)]
d2 = ele[J(nmsp[i])]
dd = merge(c.spp,d2,by='spid',allow.cartesian=TRUE)
w.mon = dd[,which.min(sqrt((lon-lon.M)^2 + (lat-lat.M)^2)), by=.(poid)]
r_el = d2[w.mon[,V1],.(elev_min,elev_max)]
c.spp[,':='(el_mn=as.integer(r_el[,elev_min]),el_mx=as.integer(r_el[,elev_max]))]
c.spp = c.spp[elevation<=el_mx & elevation>=el_mn]
c.spp[,poid:=NULL]

fmdt2 = rbind(fmdt1,c.spp)


######
# BIND three sources of refined
# niche envelopes

fdt = rbind(u.dt12,fmdt2)
fdt[,sp := as.character(sp)]

# check all species present
length(unique(fdt[,spid]))

fdt = fdt[,1:7,with = F]

fdt = na.omit(fdt)


# DONE REFINING!
















#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Create niche envelopes
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# estimates the 2 dimensional
# Euclidean distance
eq.dist = function(x,y) sqrt((x[1L] - y[1L])^2 + (x[2L] - y[2L])^2)

# Given a list of Polygons
# it dissolves the holes by 
# connecting a 0 area passage
# to each hole
dissolve_holes = function(pols) {

	require(SDMTools)
	require(rgeos)

	p1 = slot(pols,'Polygons')

	# which holes
	wh = sapply(p1,slot,'hole')

	p.h = p1[which(wh)]

	if (length(p.h)==0L) return(pols)

	p.nh = p1[which(!wh)]

	np = list()
	pwh = c()

	d.cor = double(length(p.h))

	# determine the order to make passages
	for (i in 1L:length(p.h)) {

		h.c = slot(p.h[[i]],'coords')

		#hole of which
		how = sapply(p.nh, 
			function(x) all(pnt.in.poly(h.c,slot(x,'coords'))[,3L]==1L))
		wwh = which(how)[1]
		pwh = c(pwh,wwh)

		pnhi = p.nh[[wwh]]

		pc = slot(pnhi,'coords')

		distm = sapply(1L:nrow(h.c),function(x) apply(pc,1,eq.dist,h.c[x,]))
		d.cor[i] = min(distm)
	}

	# make passages
	for(j in order(d.cor)) {

		h.c = slot(p.h[[j]],'coords')
		#hole of which
		how = sapply(p.nh,
			function(x) all(pnt.in.poly(h.c,slot(x,'coords'))[,3L]==1L))
		wwh = which(how)
		pwh = c(pwh,wwh)

		pnhi = p.nh[[wwh]]

		pc = slot(pnhi,'coords')

		distm = sapply(1L:nrow(h.c),function(x) apply(pc,1,eq.dist,h.c[x,]))

		wm = which.min(distm)
		mrow = wm%%nrow(pc)
		mcol = ceiling(wm/nrow(pc))
		
		#remove closure
		h.c = h.c[-nrow(h.c),]

		h.c = h.c[c(mcol:nrow(h.c),seq_len(mcol)),]

		np = rbind(pc[seq_len(mrow),],h.c,pc[mrow:nrow(pc),])
		poli = Polygon(np)

		p.nh[[wwh]] = poli
	}

	fp = Polygons(p.nh,ID=slot(pols,'ID'))
	return(fp)
}




library(MASS)

# transform to logarithmic space log(x + 1)
fdt[,logpre := log1p(avg_pre)]

rt = range(fdt[,avg_tmp])
rp = range(fdt[,logpre])

# 1 percent for smooth
hx = 1/100*diff(rt)
hy = 1/100*diff(rp)

# Silverman's rule of thumb
Si_bw = function(d,n,sigma) {
	4/(d+2)^(1/(d+4)) * n^(-1/(d+4)) * sigma
}

# Scott's rule of thumb
Sc_bw = function(d,n,sigma) {
	n^(-1/(d+4)) * sigma
}

sit = Si_bw(2,nrow(fdt),sd(fdt[,avg_tmp]))
sip = Si_bw(2,nrow(fdt),sd(fdt[,logpre]))

sct = Sc_bw(2,nrow(fdt),sd(fdt[,avg_tmp]))
scp = Sc_bw(2,nrow(fdt),sd(fdt[,logpre]))


setkey(fdt,'spid')
spp = unique(fdt[,spid])

## sample 3 species to visually see the result
# sam = sample(spp,6)

# par(mfrow = c(2,3))
# for (j in sam) {
  
#   sub = fdt[J(j)]
# 	plot(sub[,avg_tmp],sub[,logpre],pch=16,col=rgb(0,0,0,0.09),bty='n',
# 		xlab='Average annual temperature', ylab='Log average monthly precipitation')

# 	rx = range(sub[,avg_tmp])
# 	ry = range(sub[,logpre])
# 	k2=kde2d(sub[,avg_tmp],sub[,logpre],h=c(hx,hy),n=150,
# 		lims=c(rx[1]-1,rx[2]+1,ry[1]-0.1,ry[2]+0.1))

# 	c = contourLines(k2,nlevels=40)

# 	wls = unlist(lapply(c,'[',1))
# 	c = c[as.numeric(which(wls==min(wls)))]

# 	lapply(c,polygon,col=NULL,border='dodgerblue3',lwd=2)
# }



########################
## Create polygons and save 
# them for analysis

# transform to log(x+1)
fdt[,log1p_pre := log1p(avg_pre)]

rt = range(fdt[,avg_tmp])
rp = range(fdt[,log1p_pre])

# 1 percent for smooth
hx = 1/100*diff(rt)
hy = 1/100*diff(rp)

setkey(fdt,'spid')
spp = unique(fdt[,spid])

# get 95% density polygons
Pol = list()
ii = 1
for (j in spp) {
	sub = fdt[J(j)]

	rx = range(sub[,avg_tmp])
	ry = range(sub[,log1p_pre])

	k2 = kde2d(sub[,avg_tmp],sub[,log1p_pre],h=c(hx,hy),n=150,
		lims=c(rx[1]-1,rx[2]+1,ry[1]-0.1,ry[2]+0.1))

	c = contourLines(k2,nlevels=40)

	wls = unlist(lapply(c,'[',1))
	c = c[as.numeric(which(wls==min(wls)))]

	c = lapply(c,function(x) as.matrix(cbind(x$x,x$y)))
	pol = lapply(c, Polygon)
	Pol.j = Polygons(pol,ID=gsub(' ', '_',as.character(sub[1,sp])))

	Pol.j = checkPolygonsHoles(Pol.j)

	Pol[[ii]] = Pol.j

	ii = ii+1
	cat(ii,'\\n')
}



# incorporate holes with 0 area 'passages'
Pol.wh = list()

for (i in 1:length(Pol)) {
	Pol.wh[[i]] = dissolve_holes(Pol[[i]])
	cat(i,'\\n')
}

################################################################################
# Pol.wh holds the climatic niche polygons for input to BEAST (after 
# transforming to kml)
################################################################################

# DONE!



