### Analyze mortality data from drought experiment with
### Structural Equation Models (SEM)

### 00. import mortality data --------------------------------------------------
setwd("~/Desktop/Dryad/")

dat <- read.csv(
      file = "mortality_data.csv" ,
      header = TRUE,
      stringsAsFactors = FALSE)

# remove doubled planted and unplanted cells
dat <- subset(dat, Remove != "Remove")

# add column for days to death (stopped watering on June 28th = Julian 179)
dat$days_to_death <- dat$Jdate_dead - 179

# remove resurrection populations
dat2 <-  dat[!grepl(pattern = "herb", x = dat$Population),]
nrow(dat2)


### 01. import climate data ----------------------------------------------------
clim_dat <- read.csv(file = "Plantago_climate_data_simple.csv", header = T, stringsAsFactors = F)

# merge data
mort_dat <- merge(dat2, clim_dat, by.x = "Population", by.y = "Population")
str(mort_dat)

# add days to germinate
mort_dat$days_to_germ <- mort_dat$Jdate_germ - mort_dat$Jdate_planted

# change names of columns for figure
names(mort_dat)[names(mort_dat) == "AI"] <- "Source.Aridity"
names(mort_dat)[names(mort_dat) == "Jdate_germ"]  <- "Germ.Time"
names(mort_dat)[names(mort_dat) == "days_to_death"]  <- "Days.To.Death"
names(mort_dat)[names(mort_dat) == "plant_height_mm_on_06_29"]  <- "Plant.Height"

# remove NAs
nrow(mort_dat)
mort_dat <- mort_dat[!is.na(mort_dat$Germ.Time), ]
mort_dat <- mort_dat[!is.na(mort_dat$Days.To.Death), ]
mort_dat <- mort_dat[!is.na(mort_dat$Plant.Height), ]


### 02. Structural Equation Model ----------------------------------------------
library(lavaan)
library(semPlot)

# lavaan DEFINITIONS
      # ~ regression (as normal R syntax)
      # =~ latent variable definition
      # ~~ covariance operator (unexplained covariance)
      
### specify direct and indirect effects in model
model1_specs <- 
      "# A path    aridity -> days to death
      Days.To.Death ~ A * Source.Aridity
      
      # B path    aridity -> germ time
      Germ.Time ~ B * Source.Aridity
      
      # C path    aridity -> height
      Plant.Height ~ C * Source.Aridity
      
      # D path    germ time -> height
      Plant.Height ~ D * Germ.Time
      
      # E path    height - > days to death
      Days.To.Death ~ E * Plant.Height 
      
      ### indirect and total effects
      CE := C * E
      BDE := B * D * E
      DE := D *E
      total_Aridity := A + CE + BDE
      "

# NOTE: 10,000 bootstrap samples takes a while
# m1 <- sem(model1_specs, data = mort_dat, se = "bootstrap", bootstrap = 10000, test = "Bollen.Stine")
m1 <- sem(model1_specs, data = mort_dat, se = "bootstrap", bootstrap = 1000, test = "Bollen.Stine")

summary(m1, standardized = T)
fitMeasures(m1, c("cfi", "rmsea", "srmr"))
parameterestimates(m1)

# guidelines for model fit
      # chi-squared p < 0.05 (larger p-values here indicate good fit)
            # chi-square (sensitive to large sample sizes)
      # CFI Comparative Fit Index > 0.95
      # RMSEA Root Mean Square Error of Approximation < 0.06
      # SRMR Standardized Root Mean Square Residual < 0.08


# plot model -------------------------------------------------------------------
p1 <-  semPaths(m1, what = "std", 
         nCharNodes = 0, 
         sizeMan = 24, 
         sizeMan2 = 8,
         edge.label.color = "black",
         edge.label.cex = 1.5, 
         fade = F, 
         layout = "tree3", 
         curvePivot = T, 
         curvePivotShape = 1,
         theme = "colorblind",
         rotation = 3, 
         residuals = F,
         color = "lightgray",
         asize = 5,
         border.width = 2,
         label.font = 4,
         edge.label.font = 2)

# add significance values to paths
library(semptools)

p2  <-  mark_sig(p1, m1)
plot(p2)


### 03. Structural Equation Model for SLA dataset ------------------------------
sla_dat <- read.csv(file = "SLA_data.csv", header = T, stringsAsFactors = F)

# add SLA column
sla_dat$SLA <- sla_dat$leaf_area_mm_squared / sla_dat$leaf_weight_mg

# remove records without SLA data
sla_dat <- sla_dat[!is.na(sla_dat$SLA),]

# remove resurrection populations
sla_dat <- sla_dat[!grepl(pattern = "herb", x = sla_dat$Population),]
nrow(sla_dat)

# merge SLA data and climate data
sla_dat <- merge(sla_dat, clim_dat, by.x = "Population", by.y = "Population")
nrow(sla_dat)

# add column for days to death (stopped watering on June 28th = Julian 179)
sla_dat$days_to_death <- sla_dat$Jdate_dead - 179

# change names of columns for figure
names(sla_dat)[names(sla_dat) == "AI"] <- "Source.Aridity"
names(sla_dat)[names(sla_dat) == "Jdate_germ"]  <- "Germ.Time"
names(sla_dat)[names(sla_dat) == "days_to_death"]  <- "Days.To.Death"
names(sla_dat)[names(sla_dat) == "plant_height_mm_on_06_29"]  <- "Plant.Height"
names(sla_dat)[names(sla_dat) == "SLA"]  <- "Specific.Leaf.Area"

# remove NAs
sla_dat <- sla_dat[!is.na(sla_dat$Germ.Time), ]
sla_dat <- sla_dat[!is.na(sla_dat$Days.To.Death), ]
sla_dat <- sla_dat[!is.na(sla_dat$Plant.Height), ]

nrow(sla_dat)

### build model
sla_m3_specs <- 
      " # path A, aridity -> days to death
      Days.To.Death ~ A * Source.Aridity
      
      # path B, aridity -> sla
      Specific.Leaf.Area ~ B * Source.Aridity
      
      # path C, aridity -> height
      Plant.Height ~ C * Source.Aridity
      
      # path D, aridity -> germ time
      Germ.Time ~ D * Source.Aridity
      
      # path E, germ time -> height
      Plant.Height ~ E * Germ.Time
      
      # path F germ time -> sla
      Specific.Leaf.Area ~~ F * Germ.Time
      
      # path G, sla -> height
      Plant.Height ~~ G * Specific.Leaf.Area
      
      # path H, sla -> days to death
      Days.To.Death ~ H * Specific.Leaf.Area
      
      # path I, plant height -> days to death 
      Days.To.Death ~ I * Plant.Height 

      
      ### indirect and total effects
      BH := B * H
      CI := C * I
      DEI := D * E * I
      EI := E * I
      
      GI := G * I
      FGI := F * G * I
      "

# NOTE, add 10,000 bootstrap samples, but model takes a while to run.
m3_sla <- sem(sla_m3_specs, data = sla_dat, se = "bootstrap", bootstrap = 1000, test = "Bollen.Stine")

summary(m3_sla, standardized = T)
fitMeasures(m3_sla, c("pvalue", "cfi", "rmsea", "srmr"))


# plot model
p3 <- semPaths(m3_sla, what = "std", 
               nCharNodes = 0, 
               sizeMan = 12, 
               sizeMan2 = 8,
               edge.label.color = "black",
               edge.label.cex = 1.15, 
               label.cex = 1.125,
               fade = F, 
               layout = "tree3", 
               curvePivot = T, 
               curvePivotShape = 1,
               theme = "colorblind",
               rotation = 3, 
               residuals = F,
               color = "lightgray",
               asize = 5,
               border.width = 2,
               label.font = 4,
               edge.label.font = 2)

p4  <-  mark_sig(p3, m3_sla)
plot(p4)


### 04. plot multi-panel Figure 6 ----------------------------------------------
dev.new()

par(mfrow = c(1,2))
plot(p2)
plot(p4)

# quartz.save(file = "Figure6_SEM.jpg", dpi = 300, type = "jpg")

