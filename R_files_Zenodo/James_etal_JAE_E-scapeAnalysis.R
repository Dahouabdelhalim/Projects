#' """ Data analysis
#'     for James et al. Jounral of Animal Ecology
#'     @author: Ryan James
#'     date: 11/12/21""

library(tidyverse)
library(landscapemetrics)
library(raster)
library(sf)
library(fasterize)

# function to calculate area and edge habitat around sampling point
# r = base raster file of sampling area 
# data = data frame with sample lat lon and mixing model results
# radius = radius around point to calculate metrics
# site = vector used to designate plot_id
circIEI = function(r, data, radius, site){
  require(raster)
  require(landscapemetrics)
  require(tidyverse)
  require(sf)
  
  # convert into sf
  cds = st_as_sf(data, coords = c('Lon', 'Lat'))
  
  # set crs to WGS 84
  st_crs(cds) = 4326
  
  # convert crs to utm 15
  cds = st_transform(cds, 2027)
  
  # make empty data frame to calculate cover area at each site
  df = tibble(Site = site, f_marsh = NA, f_water = NA, f_edge = NA, f_tree = NA)
  
  # calculate the cover areas around sampling point in circle of size radius
  met = sample_lsm(r, y = cds, plot_id = cds$Site, size = radius, 
                   what = c('lsm_c_pland'),
                   shape = 'circle', return_raster = T, progress = T)
  
  # assign cover area and measure edge area
  for (i in 1:nrow(df)){
    
    # filter out single plot
    a = met %>% filter(plot_id == df$Site[i], metric == 'pland') 
    
    # initilize the areas of habitat types
    f_tree = 0
    f_marsh = 0
    f_water = 0

    # renames cover types from 1-5 to tree, marsh, water
    for (j in 1:nrow(a)){
     if (a$class[j] == 1){
         f_tree = a$value[j]/100
      }else if (a$class[j] == 2){
         f_marsh = a$value[j]/100
      }else if (a$class[j] == 3){
        f_water = a$value[j]/100
     }
    }
    
    # calculate the total edge between habitat classes
    # multiply the number of cells that are adjacent
    # and multiply by resolution of cells
    # rr = raster layer of habitat cover types
    per = get_adjacencies(a$raster_sample_plots[1])[[1]]*res(a$raster_sample_plots[[1]])[1]
    
    # calculate total area of raster 
    t = lsm_l_ta(a$raster_sample_plots[1])
    ta = t$value*10000
    
    # calculate the edge of marsh (class 2) and water (class 3)
    if (f_water > 0 & f_marsh > 0){
      mar_edge = per['2','3']*2
    }else{
      mar_edge = 0
    }
    
    # calculates the edge of magrove and water
    if (f_water > 0 & f_tree > 0){
      man_edge = per['1','3']*2
    }else{
      man_edge = 0
    }
    
    # calculate f_edge
    f_edge = (man_edge + mar_edge)/ta
    
    df$f_marsh[i] = f_marsh
    df$f_water[i] = f_water
    df$f_edge[i] = f_edge
    df$f_tree[i] = f_tree
    
  }
  
  #data$Site = as.factor(data$Site)
  
  d = full_join(data, df, by = 'Site')
  
  d$edgeiei = d$Algae/(d$f_edge)
  d$marshiei = d$SPART/(d$f_marsh)
  d$wateriei = d$POM/(d$f_water) 
  
  return(d)
}

# function to make an E-scape
# r = raster of area
# iei = data frame with the IEI values for 
# size = size of cell to generate E-scape
# rast = F; if T outputs as raster, default is sf
# prog = F; if T progress of sample_lsm is T
E_scape = function(r, iei, size, rast = F, prog = F){
  require(raster)
  require(landscapemetrics)
  require(tidyverse)
  require(sf)
  require(fasterize)
  
  # make grid over raster
  grid_geom = st_make_grid(r, cellsize = size)
  grid = st_sf(geom = grid_geom)
  
  # calculate the cover areas around sampling point in circle of size radius
  met = sample_lsm(r, grid, size = 0,
                   level = "class", metric = 'pland',
                   return_raster = T, progress = prog)
  
  # make empty data frame to calculate cover area at each site
  df = tibble(Site = unique(met$plot_id), f_marsh = NA, f_water = NA, 
              f_edge = NA, f_tree = NA, pi = NA)
  
  # assign cover area and measure edge area
  for (i in 1:nrow(df)){
    
    # filter out single plot
    a = met %>% filter(plot_id == df$Site[i], metric == 'pland') 
    
    # initilize the areas of habitat types
    f_tree = 0
    f_marsh = 0
    f_water = 0
    
    # renames cover types from 1-5 to tree, marsh, water
    for (j in 1:nrow(a)){
      if (a$class[j] == 1){
        f_tree = a$value[j]/100
      }else if (a$class[j] == 2){
        f_marsh = a$value[j]/100
      }else if (a$class[j] == 3){
        f_water = a$value[j]/100
      }
    }
    
    # calculate the total edge between habitat classes
    # multiply the number of cells that are adjacent
    # and multiply by resolution of cells
    # rr = raster layer of habitat cover types
    per = get_adjacencies(a$raster_sample_plots[1])[[1]]*res(a$raster_sample_plots[[1]])[1]
    
    # calculate total area of raster 
    t = lsm_l_ta(a$raster_sample_plots[1])
    ta = t$value*10000
    
    # calculate the edge of marsh (class 2) and water (class 3)
    if (f_water > 0 & f_marsh > 0){
      mar_edge = per['2','3']*2
    }else{
      mar_edge = 0
    }
    
    # calculates the edge of magrove and water
    if (f_water > 0 & f_tree > 0){
      man_edge = per['1','3']*2
    }else{
      man_edge = 0
    }
    
    # calculate f_edge
    f_edge = (man_edge + mar_edge)/ta
    
    df$f_marsh[i] = f_marsh
    df$f_water[i] = f_water
    df$f_edge[i] = f_edge
    df$f_tree[i] = f_tree
    df$pi[i] = a$percentage_inside
    
    if (i/nrow(df) == 0.05){
      cat('5% done \\n')
    } else if (i/nrow(df) == 0.25) {
      cat('25% done \\n')
    } else if (i/nrow(df) == 0.5) {
      cat('50% done \\n')
    } else if (i/nrow(df) == 0.75){
      cat('75% done \\n')
    } else if (i/nrow(df) == 0.95){
      cat('95% done \\n')
    }
  
  }
  
  # calculate HRI
  edgeiei = median(iei$edgeiei, na.rm = T)
  marshiei = median(iei$marshiei, na.rm = T)
  wateriei = median(iei$wateriei, na.rm = T)
  
  df$HRI = df$f_marsh*marshiei + df$f_edge*edgeiei + df$f_water*wateriei

  
  d = bind_cols(grid, df) %>% filter(pi > 90)
  
  if (rast == T){
    x = (extent(d)[2]-extent(d)[1]) %/% size
    y = (extent(d)[4]-extent(d)[3]) %/% size
    
    # make a raster of the points from the area
    rast = raster(ncol = x, nrow = y)
    extent(rast) = extent(grid)
    
    d = fasterize(grid, rast, field = 'HRI')
  }
  
  return(d)
} 

# calculate IEI values for white shrimp at multiple spatial scales----
# load data as a sf 
# calculate the HRI around each sampling point
data = read_csv('wsData.csv') 

# load raster for Fouchon
r = raster('PtFou1m.tif')

# radius length around sample locations
size = c(50, 75, 100, 150, 200, 250, 300, 400, 500, 750, 1000, 1500)

# create data frame to store information
IEI = data.frame(size = size, 
                 edgeiei = NA, edgesd = NA, 
                 wateriei = NA, watersd = NA,
                 marshiei = NA, marshsd = NA,
                 hri = NA, hrisd = NA, 
                 bioT = NA, bioP = NA, bAIC = NA,
                 calT = NA, calP = NA, cAIC = NA,
                 tcalT = NA, tcalP = NA, tAIC = NA,
                 abnT = NA, abnP = NA, aAIC = NA,
                 wgtT = NA, wgtP = NA, wAIC = NA)

for (i in 1:length(size)){
  df = circIEI(r, data, size[i], data$Site)
  df$size = size[i]
  write_csv(df, paste0('data/ws',size[i],'.csv'))
  assign(paste0("circ", size[i]),df)
}

circ = list(circ50, circ75, circ100, circ150, 
            circ200, circ250, circ300, circ400, 
            circ500, circ750, circ1000, circ1500)

# iterate for each radius size
for (i in 1:length(circ)){
  df = circ[[i]]
  
  # generate IEI values for each habitat/resource combination
  edgeIEI = median(df$edgeiei, na.rm = T)
  marshIEI = median(df$marshiei, na.rm = T)
  waterIEI = median(df$wateriei, na.rm = T)
  
  # generate hri within each buffer around sampling location
  # calculate total calories (tcal)
  # calculate average weight (wgt) at each sampling location
  df = df %>% 
    mutate(hri = f_edge*edgeIEI + f_marsh*marshIEI + f_water*waterIEI,
           tcal = cal * bio,
           wgt = bio/abundance)
  
  # fill dataframe for each IEI and HRI with standard deviations
  IEI$edgeiei[i] = median(df$edgeiei, na.rm = T)
  IEI$edgesd[i] = sd(df$edgeiei[is.finite(df$edgeiei)], na.rm = T)
  IEI$wateriei[i] = median(df$wateriei, na.rm = T)
  IEI$watersd[i] = sd(df$wateriei[is.finite(df$wateriei)], na.rm = T)
  IEI$marshiei[i] = median(df$marshiei, na.rm = T)
  IEI$marshsd[i] = sd(df$marshiei[is.finite(df$marshiei)], na.rm = T)
  IEI$hri[i] = mean(df$hri, na.rm = T)
  IEI$hrisd[i] = sd(df$hri, na.rm = T)
  
  # remove outliers and NAs for biomass and glm
  # store T value, P value, and AIC of model
  dq = subset(df, bio < 30)
  regm = glm(formula = dq$bio ~ dq$hri, family = Gamma(link = 'log'))
  rs = summary(regm)
  IEI$bioT[i] = rs$coefficients[2,3]
  IEI$bioP[i] = rs$coefficients[2,4]
  IEI$bAIC[i] = regm$aic
  
  # remove outliers and NAs for cal/g and glm
  # store T value, P value, and AIC of model
  c = subset(df, cal > 3000)
  cals = glm(formula = c$cal ~ c$hri, family = gaussian)
  cs = summary(cals)
  IEI$calT[i] = cs$coefficients[2,3]
  IEI$calP[i] = cs$coefficients[2,4]
  IEI$cAIC[i] = cals$aic
  
  # remove outliers and NAs for total calories and glm
  # store T value, P value, and AIC of model
  dm = subset(df, tcal < 100000)
  tcalm = glm(formula = dm$tcal ~ dm$hri, family = Gamma(link = 'log'))
  ts = summary(tcalm)
  IEI$tcalT[i] = ts$coefficients[2,3]
  IEI$tcalP[i] = ts$coefficients[2,4]
  IEI$tAIC[i] = tcalm$aic
  
  # remove outliers and NAs for abundance and glm
  # store T value, P value, and AIC of model
  da = df %>% filter(abundance < 250)
  abnm = glm(formula = da$abundance ~ da$hri, family = Gamma(link = 'log'))
  ind = summary(abnm)
  IEI$abnT[i] = ind$coefficients[2,3]
  IEI$abnP[i] = ind$coefficients[2,4]
  IEI$aAIC[i] = abnm$aic
  
  # remove outliers and NAs for body size and glm
  # store T value, P value, and AIC of model
  dw = df %>% filter(wgt < 0.5)
  wgtm = glm(formula = dw$wgt ~ dw$hri, family = Gamma(link = 'log'))
  ww = summary(wgtm)
  IEI$wgtT[i] = ww$coefficients[2,3]
  IEI$wgtP[i] = ww$coefficients[2,4]
  IEI$wAIC[i] = wgtm$aic
  
  write_csv(df, paste0('data/circ', size[i],'.csv'))
  assign(paste0("c", df$size[1]),df)
  
}
write_csv(IEI, 'WSieiSizes.csv')


# Make E-scape based on 200 m radius cirlce-----
# load raster
r = raster('PtFou1m.tif')
# load IEI data
d400 = read_csv('data/circ200.csv')
# generate E-scape
E400 = E_scape(r, d400, 400, rast = T)
# write as raster for plotting in QGIS
writeRaster(E400, 'E_scape400fou.tiff')

#Figures (besides map) and tables in Manuscript
# tables -----
T1 = IEI %>% 
  select(size, edgeiei, edge25, edge75,
         wateriei, water25, water75,
         marshiei, marsh25, marsh75,
         hri, hrisd)%>%
  mutate_if(is.numeric, round, digits = 2)
sub1 =  tibble(Size = T1$size,
               EdgeIEI = paste0(format(T1$edgeiei, nsmall = 2),' (',format(T1$edge25, nsmall =2),'-',format(T1$edge75),')'),
               WaterIEI =  paste0(format(T1$wateriei, nsmall = 2),' (',format(T1$water25, nsmall =2),'-',format(T1$water75),')'),
               MarshIEI =  paste0(format(T1$marshiei, nsmall = 2),' (',format(T1$marsh25, nsmall =2),'-',format(T1$marsh75),')'),
               HRI = paste0(format(T1$hri, nsmall = 1), ' \\u00B1 ', format(T1$hrisd,nsmall =1)))
write.csv(sub1, 'Tables/Table1.csv', row.names = F)

T2 = IEI %>%
  select(size, hri , hrisd , 
         bioT , bioP , bAIC ,
         calT , calP , cAIC ,
         tcalT , tcalP , tAIC ,
         abnT , abnP , aAIC ,
         wgtT , wgtP , wAIC ) %>%
  mutate_if(is.numeric, round, digits = 1)
t2 = tibble(Size = T2$size,
            HRI = paste0(format(T2$hri, nsmall = 1), ' \\u00B1 ', format(T2$hrisd,nsmall =1)),
            T2$bioT , IEI$bioP , T2$bAIC ,
            T2$calT , IEI$calP , T2$cAIC ,
            T2$tcalT , IEI$tcalP , T2$tAIC ,
            T2$abnT , IEI$abnP , T2$aAIC ,
            T2$wgtT , IEI$wgtP , T2$wAIC 
)
write.csv(t2, 'Tables/Table2.csv', row.names = F)

# Figures----
library(tidyverse)
df = read_csv('data/circ200.csv')

dq = subset(df, bio < 30)
ggplot(dq, aes(hri, bio)) + 
  geom_point()+
  geom_smooth(method = "glm", se = F, method.args = list(family = Gamma(link = 'log')))+
  labs(x = 'Habitat resource index', y = 'White shrimp biomass (g)')+
  #scale_x_continuous(limits = c(0, 2.5))+
  theme_classic()+ theme(axis.title = element_text(size = 20), 
                         plot.title = element_text(size = 20)) + theme(axis.text = element_text(size = 17))
ggsave("HRIbio.tiff", units="in", width=5, height=4.5, dpi=600,compression = 'lzw')

c = subset(df, cal > 3000)
cals = glm(formula = c$cal ~ c$hri, family = gaussian)

# dm = subset(df, tcal < 200000)
dm = subset(df, tcal < 100000)
tcalm = glm(formula = dm$tcal ~ dm$hri, family = Gamma(link = 'log'))

ggplot(dm, aes(hri, tcal))+
  geom_point()+
  geom_smooth(method = "glm", se = F, method.args = list(family = Gamma(link = 'log')))+
  labs(x = 'Habitat resource index', y = 'White shrimp total calories')+
  #scale_x_continuous(limits = c(0, 2.5))+
  theme_classic()+ theme(axis.title = element_text(size = 20), 
                         plot.title = element_text(size = 20)) + theme(axis.text = element_text(size = 17))
ggsave("HRItcal.tiff", units="in", width=5, height=4.5, dpi=600,compression = 'lzw')


da = df %>% filter(abundance < 250)
abnm = glm(formula = da$abundance ~ da$hri, family = Gamma(link = 'log'))

ggplot(da, aes(hri, abundance))+
  geom_point()+
  geom_smooth(method = "glm", se = F, method.args = list(family = Gamma(link = 'log')))+
  labs(x = 'Habitat resource index', y = 'White shrimp abundance')+
  #scale_x_continuous(limits = c(0, 2.5))+
  theme_classic()+ theme(axis.title = element_text(size = 20), 
                         plot.title = element_text(size = 20)) + theme(axis.text = element_text(size = 17))
ggsave("HRIabund.tiff", units="in", width=5, height=4.5, dpi=600,compression = 'lzw')



dw = df %>% filter(wgt < 0.5)
wgtm = glm(formula = dw$wgt ~ dw$hri, family = Gamma(link = 'log'))

ggplot(dw, aes(hri, wgt))+
  geom_point()+
  geom_smooth(method = "glm", se = F, method.args = list(family = Gamma(link = 'log')))+
  labs(x = 'Habitat resource index', y = 'White shrimp mean weight (g)')+
  #scale_x_continuous(limits = c(0, 2.5))+
  theme_classic()+ theme(axis.title = element_text(size = 20), 
                         plot.title = element_text(size = 20)) + theme(axis.text = element_text(size = 17))
ggsave("HRIwgt.tiff", units="in", width=5, height=4.5, dpi=600,compression = 'lzw')


