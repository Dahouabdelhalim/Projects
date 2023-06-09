
#-----------------------------------
# This script runs all computational
# notebooks used in the article:
#
# Fine-scale heterogeneity in 
# Schistosoma mansoni force of 
# infection measured through 
# serologic surveillance
#-----------------------------------
library(here)
here::here()

#-----------------------------------
# Raw study data processing
# to create public data files
#
# (not run for public replication)
#-----------------------------------
rmarkdown::render(here::here("R/1-mbita-schisto-format-data.Rmd"),
                  output_file = here::here("R/1-mbita-schisto-format-data.html"))

#-----------------------------------
# Table 1
#-----------------------------------
rmarkdown::render(here::here("R/2-mbita-schisto-sample-summary.Rmd"),
                  output_file = here::here("R/2-mbita-schisto-sample-summary.html"))

#-----------------------------------
# Figure 1
# Figure S2
# Figure S3
#-----------------------------------
rmarkdown::render(here::here("R/3-mbita-schisto-map.Rmd"),
                  output_file = here::here("R/3-mbita-schisto-map.html"))

#-----------------------------------
# Figure 4c, 4d
# Fig S5
#-----------------------------------
rmarkdown::render(here::here("R/4-mbita-schisto-by-age.Rmd"),
                  output_file = here::here("R/4-mbita-schisto-by-age.html"))

#-----------------------------------
# Figure 4a, 4b
# Figure S6
#-----------------------------------
rmarkdown::render(here::here("R/5-mbita-schisto-distance-from-lake.Rmd"),
                  output_file = here::here("R/5-mbita-schisto-distance-from-lake.html"))


#-----------------------------------
# Figure 2
# Figure 3
# Figure S4
#-----------------------------------
rmarkdown::render(here::here("R/6-mbita-schisto-SEA-vs-KK.Rmd"),
                             output_file = here::here("R/6-mbita-schisto-SEA-vs-KK.html"))

#-----------------------------------
# Figure S1
#-----------------------------------
rmarkdown::render(here::here("R/7-mbita-schisto-antibody-distributions.Rmd"),
                  output_file = here::here("R/7-mbita-schisto-antibody-distributions.html"))

