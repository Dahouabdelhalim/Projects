## RScript to fit light response curves
## Code written by Eleinis Avila-Lovera
## Last modified on 5/5/21

#Set working directory
setwd()

#Read in data
data_10_500<-read.csv("Chen_LICOR_10-500.csv", header=T)

## Using package "photosynthesis" (https://github.com/cdmuir/photosynthesis)
## Non-rectangular hyperbolic model of light response from  Marshall B, Biscoe P. 1980. A model for C3 leaves 
## describing the dependence of net photosynthesis on irradiance. J Ex Bot 31:29-39
library(photosynthesis)

#Fit curves using 10-500 PAR data
fits2<-fit_many(
  data=data_10_500, 
  varnames=list(A_net="photo", PPFD="par", group="id"
  ), 
  funct=fit_aq_response, group="id"
)

#Explore results
summary(fits2[[3]][[1]]) #print model summary for list 3
fits2[[3]][[2]] #print fitted parameters (second sublist) of list 3
fits2[[3]][[3]] #print graph (third sublist) of list 3

#Compile parameters into dataframe for analysis
fits2_pars <- compile_data(fits2,
                          output_type = "dataframe",
                          list_element = 2
)

write.csv(fits2_pars, file="A_Q_results2_10-500.csv", row.names=FALSE)
