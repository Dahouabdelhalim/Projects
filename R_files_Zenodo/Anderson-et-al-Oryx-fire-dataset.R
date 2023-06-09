############################################################################
# 
# Author: Daniel M. Griffith & T. Michael Anderson
# Date: Feb 2018
# Description: Read and plot the Serengeti National Park fire burn scar 
#   dataset as reported in:
#
# Anderson, T.M., P.M. Ngoti, M.L. Nzunda, D.M. Griffith, J.D.M. Speed, 
#   F. Fossøy, E. Røskaft and B.J. Graae. The burning question: does fire 
#   affect habitat selection and forage preference of black rhinos (Diceros 
#   bicornis) in East African savannas? In press, Oryx.
#
############################################################################

# SET WD TO SCRIPT LOCATION
  setwd(".")

# LOAD LIBRARIES
  library(raster)

# READ RASTER DATA
  ffreq <- raster("Anderson-et-al-Oryx-fire-frequency-2000-2016.tif")
  fires <- stack("Anderson-et-al-Oryx-annual-fire-bands-2000-2016.tif")

  # ffreq (some pixels burned more than once per year)
    # class       : RasterLayer 
    # dimensions  : 1287, 995, 1280565, 1  (nrow, ncol, ncell, nlayers)
    # resolution  : 231.6564, 231.6564  (x, y)
    # extent      : 3763026, 3993524, -403082.1, -104940.3  (xmin, xmax, ymin, ymax)
    # coord. ref. : +proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs 
    # names       : Anderson.et.al.Oryx.fire.frequency.2000.2016 
    # min values  :                                            0 
    # max values  :                                           26 
  # fires (binary flags for burn scars in burn years)
    # class       : RasterStack 
    # dimensions  : 1287, 995, 1280565, 16  (nrow, ncol, ncell, nlayers)
    # resolution  : 231.6564, 231.6564  (x, y)
    # extent      : 3763026, 3993524, -403082.1, -104940.3  (xmin, xmax, ymin, ymax)
    # coord. ref. : +proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs 

# Label the burn years in 'fires' (burn years are assessed from May 1st to 
# April 30th of the following year)
  names(fires) <- paste("Burn_year_", 2000:2015, "_", 2001:2016, sep = "")

# Simple plotting for fire data
  plot(ffreq)
  plot(fires)



