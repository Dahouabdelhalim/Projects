## 03 - Point process analysis

# Libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sp)
library(spatstat)
library(raster)
library(reshape2)
library(ggplot2)
library(scico)
library(rnaturalearth)

# Read data and make folders ----------------------------------------------

flr.data <- read.csv(file = "./Data/Flickr_results_clean.csv",
                     stringsAsFactors = FALSE)

loc.outputs <- "./Outputs/"
loc.figure <- "./Figures/"
loc.gbifdata <- "./Data/GBIFData/"
gbif.data <- read.csv(file = paste0(loc.gbifdata, "7_compiled_GBIFData_coordcleanComplete.csv"),
                      stringsAsFactors = FALSE)

# Land area download / loading --------------------------------------------

ref <- rnaturalearth::ne_download(scale = 110, type = "land",
                                  category = "physical", load = TRUE)
saveRDS(file = "landmass.Rds", object = land)
land <- readRDS(file = "landmass.Rds")
# create window based on land, have to modify the rnaturalearth shapefile, to
# make it conform to spatstat requirements
land <- crop(land, extent(-180, 180, -60, 83.6341))
land <- as(land, "SpatialPolygons") 
f.land <- fortify(land)
table(f.land$group)
f.landlist <- split(f.land, f = f.land$group)
pol.landlist <- lapply(f.landlist, function(x){data.frame("x" = x$long, "y" = x$lat)} )
class(pol.landlist)
# have to reverse the order for holes
i <- 0
for(p in f.landlist){
  i <- i + 1
  # p <- f.list[1]
  df <- p
  if(!df$hole[1] == TRUE){
    df <- df[dim(df)[1]:1,]
    pol.landlist[[i]] <- data.frame("x" = df$long, "y" = df$lat)
  } else {
    pol.landlist[[i]] <- data.frame("x" = df$long, "y" = df$lat)
  } # end of if statement
} # end of loop  

# then convert the polygon to a owin format for spatstat
spatstat.options(checkpolygons=FALSE)
land.owin <- owin(poly = pol.landlist)
spatstat.options(checkpolygons=TRUE)
# plot(land.owin)


# Create ppp object with Flickr data --------------------------------------
lon <- flr.data$longitude
lat <- flr.data$latitude
xrange <- range(lon, na.rm=T)
yrange <- range(lat, na.rm=T)
flr.ppp <- ppp(lon, lat, xrange, yrange, window = land.owin)
# plot(flr.ppp)

# Create ppp object with GBIF data --------------------------------------
lon <- gbif.data$lon
lat <- gbif.data$lat
xrange <- range(lon, na.rm=T)
yrange <- range(lat, na.rm=T)
gbif.ppp <- ppp(lon, lat, xrange, yrange, window = land.owin)
# plot(flr.ppp)

# grid based on 10 degrees by 10 degrees (roughly)
q.res <- quadrat.test(flr.ppp,
                      xbreaks = seq(land.owin$xrange[1], land.owin$xrange[2], 10),
                      ny = 14)
saveRDS(file = paste0(loc.outputs, "/qt_results_flickr.Rds"), object = q.res)

q.res <- quadrat.test(gbif.ppp,
                      xbreaks = seq(land.owin$xrange[1], land.owin$xrange[2], 10),
                      ny = 14)
saveRDS(file = paste0(loc.outputs, "/qt_results_gbif.Rds"), object = q.res)

# Nearest neighbor distances (G function)
gtest.f <- Gest(flr.ppp, correction = c("rs", "km", "han"))
gtest.g <- Gest(gbif.ppp, correction = c("rs", "km", "han"))

gtdf.f <- data.frame("r" = gtest.f$r, "theo" = gtest.f$theo, "han" = gtest.f$han,
           "rs" = gtest.f$rs, "km" = gtest.f$km, "dataset" = "Flickr")
gtdf.g <- data.frame("r" = gtest.g$r, "theo" = gtest.g$theo, "han" = gtest.g$han,
           "rs" = gtest.g$rs, "km" = gtest.g$km, "dataset" = "GBIF")
gtdf <- rbind(gtdf.g, gtdf.f)
gtdf.m <- melt(gtdf, c("r", "dataset"))

gtest.plot <- ggplot(gtdf.m) +
  geom_line(aes(x = r, y = value, colour = variable), size = 1) +
  facet_wrap(.~dataset, scales = "free") +
  theme_bw() +
  scale_colour_manual(values = scico(palette = "roma", n = 4),
                      labels = c(expression(G[italic(pois)](r)),
                                 expression(Ĝ[italic(han)](r)),
                                 expression(Ĝ[italic(rs)](r)),
                                 expression(Ĝ[italic(km)](r)))) +
  labs(x = "r", y = "G(r)", colour = "") +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = 4))

ggsave(filename = paste0(loc.figure, "SUPPLEMENTARY FIGURE 1 - gtest_plot.png"),
       plot = gtest.plot, width = 21, height = 12, units = "cm",
       dpi = 300, device = "png")
ggsave(filename = paste0(loc.figure, "SUPPLEMENTARY FIGURE 1 - gtest_plot.pdf"),
       plot = gtest.plot, width = 21, height = 12, units = "cm",
       device = "pdf")

# Empty space distances (F function)
ftest.f <- Fest(flr.ppp, correction=c("rs", "km", "cs"))
ftest.g <- Fest(gbif.ppp, correction=c("rs", "km", "cs"))

ftdf.f <- data.frame("r" = ftest.f$r, "theo" = ftest.f$theo, "han" = ftest.f$cs,
                     "rs" = ftest.f$rs, "km" = ftest.f$km, "dataset" = "Flickr")
ftdf.g <- data.frame("r" = ftest.g$r, "theo" = ftest.g$theo, "han" = ftest.g$cs,
                     "rs" = ftest.g$rs, "km" = ftest.g$km, "dataset" = "GBIF")
ftdf <- rbind(ftdf.g, ftdf.f)
ftdf.m <- melt(ftdf, c("r", "dataset"))

ftest.plot <- ggplot(ftdf.m) +
  geom_line(aes(x = r, y = value, colour = variable), size = 1) +
  facet_wrap(.~dataset, scales = "free") +
  theme_bw() +
  scale_colour_manual(values = scico(palette = "roma", n = 4),
                      labels = c(expression(F[italic(pois)](r)),
                                 expression(F[italic(cs)](r)),
                                 expression(F[italic(bord)](r)),
                                 expression(F[italic(km)](r)))) +
  labs(x = "r", y = "F(r)", colour = "") +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = 4))

ggsave(filename = paste0(loc.figure, "SUPPLEMENTARY FIGURE 2 - ftest_plot.png"),
       plot = ftest.plot, width = 21, height = 12, units = "cm",
       dpi = 300, device = "png")
ggsave(filename = paste0(loc.figure, "SUPPLEMENTARY FIGURE 2 - ftest_plot.pdf"),
       plot = ftest.plot, width = 21, height = 12, units = "cm",
       device = "pdf")

# Pairwise distance (K function), or the distance between all pairs of points.
# due to the size of the dataset we cannot compute this with any correction
ktest.f <- Kest(flr.ppp, correction = "none",
              nlarge = 3000)
ktest.g <- Kest(gbif.ppp, correction = "none",
              nlarge = 3000)

ktdf.f <- data.frame("r" = ktest.f$r, "theo" = ktest.f$theo, "un" = ktest.f$un,
                     "dataset" = "Flickr")
ktdf.g <- data.frame("r" = ktest.g$r, "theo" = ktest.g$theo, "un" = ktest.g$un,
                     "dataset" = "GBIF")
ktdf <- rbind(ktdf.g, ktdf.f)
ktdf.m <- melt(ktdf, c("r", "dataset"))

ktest.plot <- ggplot(ktdf.m) +
  geom_line(aes(x = r, y = value, colour = variable), size = 1) +
  facet_wrap(.~dataset, scales = "free") +
  theme_bw() +
  scale_colour_manual(values = scico(palette = "roma", n = 4),
                      labels = c(expression(K[italic(pois)](r)),
                                 expression(K[italic(un)](r)))) +
  labs(x = "r", y = "K(r)", colour = "") +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = 4))

ktest.plot

ggsave(filename = paste0(loc.figure, "SUPPLEMENTARY FIGURE 3 - ktest_plot.png"),
       plot = ktest.plot, width = 21, height = 12, units = "cm",
       dpi = 300, device = "png")
ggsave(filename = paste0(loc.figure, "SUPPLEMENTARY FIGURE 3 - ktest_plot.pdf"),
       plot = ktest.plot, width = 21, height = 12, units = "cm",
       device = "pdf")

