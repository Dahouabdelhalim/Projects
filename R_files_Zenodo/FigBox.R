library(data.table)
library(weights)
library(readxl)
library(rasterVis)
library(pals)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls(all=TRUE))


#### taken from here: https://gist.github.com/scbrown86/2779137a9378df7b60afd23e0c45c188 ####
# The function that produces the colour matrix
colmat <- function(nquantiles = 3, upperleft = "#0096EB", upperright = "#820050", 
                   bottomleft = "#BEBEBE", bottomright = "#FFE60F",
                   xlab = "x label", ylab = "y label", plotLeg = TRUE,
                   saveLeg = TRUE) {
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(magrittr)
  require(classInt)
  # The colours can be changed by changing the HEX codes for:
  # upperleft, upperright, bottomleft, bottomright
  # From http://www.joshuastevens.net/cartography/make-a-bivariate-choropleth-map/
  # upperleft = "#64ACBE", upperright = "#574249", bottomleft = "#E8E8E8", bottomright = "#C85A5A",
  # upperleft = "#BE64AC", upperright = "#3B4994", bottomleft = "#E8E8E8", bottomright = "#5AC8C8",
  # upperleft = "#73AE80", upperright = "#2A5A5B", bottomleft = "#E8E8E8", bottomright = "#6C83B5", 
  # upperleft = "#9972AF", upperright = "#804D36", bottomleft = "#E8E8E8", bottomright = "#C8B35A",
  # upperleft = "#DA8DC8", upperright = "#697AA2", bottomleft = "#E8E8E8", bottomright = "#73BCA0",
  # Similar to Teuling, Stockli, Seneviratnea (2011) [https://doi.org/10.1002/joc.2153]
  # upperleft = "#F7900A", upperright = "#993A65", bottomleft = "#44B360", bottomright = "#3A88B5",
  # Viridis style
  # upperleft = "#FEF287", upperright = "#21908D", bottomleft = "#E8F4F3", bottomright = "#9874A1",
  # Similar to Fjeldsa, Bowie, Rahbek 2012
  # upperleft = "#34C21B", upperright = "#FFFFFF", bottomleft = "#595757",  bottomright = "#A874B8",
  # Default from original source
  # upperleft = "#0096EB", upperright = "#820050", bottomleft= "#BEBEBE", bottomright = "#FFE60F",
  my.data <- seq(0, 1, .01)
  # Default uses terciles (Lucchesi and Wikle [2017] doi: 10.1002/sta4.150)
  my.class <- classInt::classIntervals(my.data,
                                       n = nquantiles,
                                       style = "quantile" )
  my.pal.1 <- findColours(my.class, c(upperleft, bottomleft))
  my.pal.2 <- findColours(my.class, c(upperright, bottomright))
  col.matrix <- matrix(nrow = 101, ncol = 101, NA)
  for (i in 1:101) {
    my.col <- c(paste(my.pal.1[i]), paste(my.pal.2[i]))
    col.matrix[102 - i, ] <- findColours(my.class, my.col)
  }
  col.matrix.plot <- col.matrix %>%
    as.data.frame(.) %>% 
    mutate("Y" = row_number()) %>%
    mutate_at(.tbl = ., .vars = vars(starts_with("V")), .funs = list(as.character)) %>% 
    pivot_longer(data = ., cols = -Y, names_to = "X", values_to = "HEXCode") %>% 
    mutate("X" = as.integer(sub("V", "", .$X))) %>%
    distinct(as.factor(HEXCode), .keep_all = TRUE) %>%
    mutate(Y = rev(.$Y)) %>% 
    dplyr::select(-c(4)) %>%
    mutate("Y" = rep(seq(from = 1, to = nquantiles, by = 1), each = nquantiles),
           "X" = rep(seq(from = 1, to = nquantiles, by = 1), times = nquantiles)) %>%
    mutate("UID" = row_number())
  # Use plotLeg if you want a preview of the legend
  if (plotLeg) {
    p <- ggplot(col.matrix.plot, aes(X, Y, fill = HEXCode)) +
      geom_raster() +
      scale_fill_identity() +
      coord_equal(expand = FALSE) +
      theme_void() +
      theme(aspect.ratio = 1,
            axis.title = element_text(size = 12, colour = "black",hjust = 0.5, 
                                      vjust = 1),
            axis.title.y = element_text(angle = 90, hjust = 0.5)) +
      xlab(bquote(.(xlab) ~  symbol("\\256"))) +
      ylab(bquote(.(ylab) ~  symbol("\\256")))
    print(p)
    assign(
      x = "BivLegend",
      value = p,
      pos = .GlobalEnv
    )
  }
  # Use saveLeg if you want to save a copy of the legend
  if (saveLeg) {
    ggsave(filename = "bivLegend.pdf", plot = p, device = "pdf",
           path = "./", width = 4, height = 4, units = "in",
           dpi = 300)
  }
  seqs <- seq(0, 100, (100 / nquantiles))
  seqs[1] <- 1
  col.matrix <- col.matrix[c(seqs), c(seqs)]
}

# Function to assign colour-codes to raster data
# As before, by default assign tercile breaks
bivariate.map <- function(rasterx, rastery, colormatrix = col.matrix,
                          nquantiles = 3, export.colour.matrix = TRUE,
                          outname = paste0("colMatrix_rasValues", names(rasterx))) {
  # export.colour.matrix will export a data.frame of rastervalues and RGB codes 
  # to the global environment outname defines the name of the data.frame
  quanmean <- getValues(rasterx)
  temp <- data.frame(quanmean, quantile = rep(NA, length(quanmean)))
  brks <- with(temp, quantile(temp,
                              na.rm = TRUE,
                              probs = c(seq(0, 1, 1 / nquantiles))
  ))
  ## Add (very) small amount of noise to all but the first break
  ## https://stackoverflow.com/a/19846365/1710632
  brks[-1] <- brks[-1] + seq_along(brks[-1]) * .Machine$double.eps
  r1 <- within(temp, quantile <- cut(quanmean,
                                     breaks = brks,
                                     labels = 2:length(brks),
                                     include.lowest = TRUE
  ))
  quantr <- data.frame(r1[, 2])
  quanvar <- getValues(rastery)
  temp <- data.frame(quanvar, quantile = rep(NA, length(quanvar)))
  brks <- with(temp, quantile(temp,
                              na.rm = TRUE,
                              probs = c(seq(0, 1, 1 / nquantiles))
  ))
  brks[-1] <- brks[-1] + seq_along(brks[-1]) * .Machine$double.eps
  r2 <- within(temp, quantile <- cut(quanvar,
                                     breaks = brks,
                                     labels = 2:length(brks),
                                     include.lowest = TRUE
  ))
  quantr2 <- data.frame(r2[, 2])
  as.numeric.factor <- function(x) {
    as.numeric(levels(x))[x]
  }
  col.matrix2 <- colormatrix
  cn <- unique(colormatrix)
  for (i in 1:length(col.matrix2)) {
    ifelse(is.na(col.matrix2[i]),
           col.matrix2[i] <- 1, col.matrix2[i] <- which(
             col.matrix2[i] == cn
           )[1]
    )
  }
  # Export the colour.matrix to data.frame() in the global env
  # Can then save with write.table() and use in ArcMap/QGIS
  # Need to save the output raster as integer data-type
  if (export.colour.matrix) {
    # create a dataframe of colours corresponding to raster values
    exportCols <- as.data.frame(cbind(
      as.vector(col.matrix2), as.vector(colormatrix),
      t(col2rgb(as.vector(colormatrix)))
    ))
    # rename columns of data.frame()
    colnames(exportCols)[1:2] <- c("rasValue", "HEX")
    # Export to the global environment
    assign(
      x = outname,
      value = exportCols,
      pos = .GlobalEnv
    )
  }
  cols <- numeric(length(quantr[, 1]))
  for (i in 1:length(quantr[, 1])) {
    a <- as.numeric.factor(quantr[i, 1])
    b <- as.numeric.factor(quantr2[i, 1])
    cols[i] <- as.numeric(col.matrix2[b, a])
  }
  r <- rasterx
  r[1:length(r)] <- cols
  return(r)
}


#### read data ####
crop.diversity.data <- read_xlsx("Data_supplement_Kastner_et_al_2021.xlsx","data_crop_diversity_map", skip = 1, col_names = F)
crop.diversity.data <- as.data.table(crop.diversity.data)
setnames(crop.diversity.data, c("x",	"y", "crop.diversity","area.domestic", "area.export"))

crop.diversity.data[, export.share := area.export / (area.domestic + area.export)]
crop.diversity.data[export.share>1, export.share := 1][export.share<0, export.share := 0]


#### histograms ####
png("FigBox_histogram_domestic.png", width=6, height=4, units="in", res=600)
wtd.hist(crop.diversity.data[,crop.diversity], xlim = c(1,12), 
         weight = crop.diversity.data[,area.domestic] / 1000 / 1000, col = "beige",
         xlab = "",ylab = "", main = NULL, 
         cex.axis = 2, cex.lab = 2)
dev.off()

png("FigBox_histogram_export.png", width=6, height=4, units="in", res=600)
wtd.hist(crop.diversity.data[,crop.diversity], xlim = c(1,12), 
         weight = crop.diversity.data[,area.export] / 1000 / 1000, col = "lightblue",
         xlab = "",ylab = "", main = NULL, 
         cex.axis = 2, cex.lab = 2)
dev.off()


#### bivariate map and legend ####
export.share.map <- rasterFromXYZ(crop.diversity.data[,.(x,y,export.share)])
crop.diversity.map <- rasterFromXYZ(crop.diversity.data[,.(x,y,crop.diversity)])
nBreaks <- 3
# Create the colour matrix
col.matrix <- colmat(nquantiles = nBreaks, xlab = "Share for export production", ylab = "Crop diversity", 
                     ## non default colours
                     upperleft = brewer.seqseq2(9)[7], upperright = brewer.seqseq2(9)[9], 
                     bottomleft = brewer.seqseq2(9)[1], bottomright = brewer.seqseq2(9)[3],
                     saveLeg = FALSE, plotLeg = TRUE)

# create the bivariate raster
bivmap <- bivariate.map(rasterx = export.share.map, rastery = crop.diversity.map,
                        export.colour.matrix = TRUE, outname = "bivMapCols",
                        colormatrix = col.matrix, nquantiles = nBreaks)

# Convert to dataframe for plotting with ggplot
bivMapDF <- as.data.frame(bivmap, xy = TRUE)
bivMapDF <- as.data.table(bivMapDF)
setnames(bivMapDF,c("x","y","bivVal"))
bivMapDF <- bivMapDF[!is.na(bivVal)]

col <- unique(as.vector(col.matrix))
myTheme <- rasterTheme(region = col)
myTheme$axis.line$col <- "transparent"
my.at <- sort(unique(bivMapDF[,bivVal]))
my.at <- my.at - 0.5
my.at <- union(my.at, my.at[length(my.at)]+1)
biv.map <- rasterFromXYZ(bivMapDF)

lat.lon_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
eck4_crs<-"+proj=eck4 +ellps=WGS84"
projected_crs <- "+proj=natearth +ellps=WGS84 +over"

target_extent <- bbox(extent(-180, 180, -60, 90))
# shape files from https://github.com/nvkelso/natural-earth-vector
world.shp <- rgdal::readOGR("110m_physical/ne_110m_coastline.shp")
bb.shp <- rgdal::readOGR("110m_physical/ne_110m_graticules_all/ne_110m_wgs84_bounding_box.shp")

world.shp <- crop(world.shp, target_extent)
bb.shp <- crop(bb.shp, target_extent)
projection(world.shp) <- lat.lon_crs
projection(bb.shp) <- lat.lon_crs
world.shp_proj <- spTransform(world.shp, projected_crs) 
bb.shp_proj <- spTransform(bb.shp, projected_crs)

tf <- tempfile(fileext = '.tif')
tf2 <- tempfile(fileext = '.tif')
tf3 <- tempfile(fileext = '.tif')

projection(biv.map) <- eck4_crs
writeRaster(biv.map, tf,overwrite = TRUE)
gdalUtilities::gdalwarp(tf, tf2, t_srs = lat.lon_crs, r='bilinear', 
                        te = c(target_extent), multi=TRUE, overwrite = TRUE)
gdalUtilities::gdalwarp(tf2, tf3, t_srs = projected_crs, r='bilinear', 
                        multi=TRUE, overwrite = TRUE)

biv.map_proj <- raster(tf3)

png(paste0("FigBox_Map.png"), width=8, height=4, units="in", res=600)
plot(bb.shp_proj, col = "aliceblue", border = NA,
     xlim=extent(bb.shp_proj)[1:2], ylim=extent(bb.shp_proj)[3:4])
levelplot(biv.map_proj,margin = FALSE, at = my.at,par.settings = myTheme,maxpixels = 1e7,
          scales=list(draw=F), colorkey = F,
          main = "", add = T)+
  latticeExtra::layer(sp.lines(world.shp_proj, col="darkgray", lwd=0.5)) + 
  latticeExtra::layer(sp.lines(bb.shp_proj, col="darkgray", lwd=0.5))

dev.off()

png(paste0("FigBox_legend.png"), width=4, height=4, units="in", res=600)
BivLegend
dev.off()

unlink(c(tf,tf2,tf3))