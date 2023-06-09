
setwd("C:/Users/pedro/OneDrive/Documentos/2019/INPA/ModelagemCoronata/regioes_biogeograficas/kuenm")

# Installing and loading packages
if(!require(devtools)){
  install.packages("devtools")
}

if(!require(kuenm)){
  devtools::install_github("marlonecobos/kuenm")
}

library(kuenm)

# Preparing variables to be used in arguments
lcoronata <- "lcoronata_enm_process"
kuenm_start(file.name = lcoronata)
