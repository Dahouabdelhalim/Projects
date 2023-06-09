# A script to extract the R-packages used during analyses for citation in the manuscript. 

library(knitr)

pkgs <- c("tidyverse",
          "sf",
          "raster",
          "dssatr",
          "readxl",
          "writexl",
          "fasterize",
          "elevatr",
          "cowplot",
          "missForest",
          "ape",
          "ggtree",
          "viridis",
          "phytools",
          "ggnewscale",
          "wesanderson",
          "forcats",
          "gridExtra",
          "deeptime",
          "rnaturalearth",
          "sp",
          "car",
          "olsrr",
          "lcvplants",
          "LCVP",
          "picante")

knitr::write_bib(pkgs, "output/Rpackages_for_citation.bib")
