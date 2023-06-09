library('gdm')
library('tidyverse')
library('phyloseq')
library('cowplot')
library('wesanderson')
library('metagMisc')
library('vegan')
library('ggvegan')

annadelHood <- readRDS(file = "phyloseq.rds")
sampleData <- read.csv(file = "sampleData.csv")

annadelHood <- phyloseq_standardize_otu_abundance(annadelHood, method = "total") %>%
  subset_samples(is.na(cPerc) == FALSE) %>%
  subset_samples(is.na(pH) == FALSE)

sampleData <- sampleData %>% filter(is.na(cPerc) == FALSE, is.na(pH) == FALSE)

sppTab <- annadelHood@otu_table@.Data %>% as.data.frame()
sppTab$sample <- rownames(sppTab)

gdmPredData <- sampleData %>%
  dplyr::select(sample, ecosystem, burned, ABS_COVER, damagePerc, lon, lat, nPerc, cPerc, pH)

gdmPredData$ecosystem <- (gdmPredData$ecosystem == "DF") %>% as.numeric()
gdmPredData$burned <- (gdmPredData$burned == "B") %>% as.numeric()

#compute site pairwise GDM matrix
sitePairMatrix <- formatsitepair(bioData = sppTab, bioFormat = 1, XColumn = "lon", YColumn = "lat", 
                                 predData = gdmPredData, 
                                 siteColumn = "sample")

gdm.1 <- gdm(sitePairMatrix, geo = T)

#variable importance
set.seed(900)
varImp <- gdm.varImp(sitePairMatrix, geo = T, parallel = TRUE)

gdm.1 %>% summary()

#extract the splines
gdm.1.splineDat <- isplineExtract(gdm.1)

#now we need to transform the splines for nice plotting
gdmSplinesX <- gdm.1.splineDat$x %>% data.frame()

#scale the variables so they can all go on the same X axis
gdmSplinesX <- lapply(gdmSplinesX, scale) %>% data.frame()

gdmSplinesX$val <- 1:length(gdmSplinesX$Geographic)
gdmSplinesX <- gdmSplinesX %>% pivot_longer(-val, names_to = "var")
gdmSplinesX$val <- paste(gdmSplinesX$val, gdmSplinesX$var)

gdmSplinesY <- gdm.1.splineDat$y %>% data.frame()
gdmSplinesY$val <- 1:length(gdmSplinesY$Geographic) 
gdmSplinesY <- gdmSplinesY %>% pivot_longer(-val, names_to = "var")
gdmSplinesY$val <- paste(gdmSplinesY$val, gdmSplinesY$var)

gdmSplines <- dplyr::select(gdmSplinesY, -var) %>% 
  full_join(gdmSplinesX, ., by = "val")

#rename them properly
gdmSplines$var <- gdmSplines$var %>%
  str_replace_all("ecosystem", "Ecosystem type") %>%
  str_replace_all("cPerc", "Soil C (%)") %>%
  str_replace_all("burned", "Fire") %>%
  str_replace_all("ABS_COVER", "Canopy cover (%)") %>%
  str_replace_all("Geographic", "Geography") %>%
  str_replace_all("pH", "Soil pH")

#plot
splinePlot <- ggplot(filter(gdmSplines, !(var %in% c("damagePerc", "nPerc"))),
       aes(x = value.x, 
           y = value.y, color = var)) +
  geom_line(size = 0.8) +
  directlabels::geom_dl(aes(label= var),method="angled.boxes") +
  scale_color_manual(values = c(
    wes_palette("Rushmore1")[3],
    wes_palette("Rushmore1")[1],
    wes_palette("Rushmore1")[5],
    wes_palette("Rushmore1")[2],
    "#bababa",
    wes_palette("Rushmore1")[4],
    wes_palette("Rushmore1")[6])
    ) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Value (scaled)") +
  ylab("Partial ecological distance")

#make a prediction plot
#first, extract predictions
gdmPredictions <- predict(gdm.1, sitePairMatrix)

predictPlot <- cbind(dplyr::select(sitePairMatrix, distance), gdmPredictions) %>%
  ggplot() +
  geom_point(aes(x = distance, y = gdmPredictions), alpha = 0.6,
             size = 2,
             color = wes_palette("Rushmore1")[3]) +
  geom_abline(aes(intercept = 0, slope = 1), 
              linetype = "dashed",
              color = wes_palette("Rushmore1")[4]) +
  ylim(0.55,1) +
  xlim(0.55,1) +
  xlab("Observed ecological distance") +
  ylab("Predicted ecological distance") +
  theme_minimal()

gdmPlot <- plot_grid(splinePlot, predictPlot,
                     labels = c("a", "b"),
          ncol = 2)

save_plot(filename = "Figure4.svg",
          plot = gdmPlot,
          base_width = 10,
          base_height = 4)