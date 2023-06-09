#This source code reproduces the results shown in subsection
#3.2 "Comparison between contract classes"
#for Kolmogorov-Smirnov Test section

#author: Andrea Falcon-Cortes

#Function necessary to indicate in the graph the maximum difference
#between the CDFs
KS_interpolate <- function(x,x0,y0){
  
  y = vector("numeric", length = length(x))  
  
  for (i in 1:length(x)){
    j = 1
    if (x[i] <= x0[1]){
      y[i] = y0[1]
    } else if(x[i] >= x0[length(x0)]){
      y[i] = y0[length(y0)]
    } else{
      while(j < length(x0)){
        if(x[i] == x0[j]){
          y[i] = y0[j]
          j = length(x0)
        } else if(x0[j]<x[i] & x[i]<x0[j+1]){
          y[i] = y0[j]
          j = length(x0)
        }
        j = j+1 
      } 
    }
  }
  
  return(y)
}

#Necessary packages
library(dplyr)
library(tidyverse)
library(Rmisc)
library(scales)

#Data relative to contracts made from 2013 to 2020

#Upload the source list
data <- read.csv("~/RPS_Cont_2013-2020.csv", 
                 stringsAsFactors=FALSE)

#From the source list take only the non-dummy variables
#of type i) and ii)
fdata = subset(data, select = c(Status, DepID, ProvContID,
                                Year, BeginningWeek, EBWeeks,
                                Spending, T.Cont.Max, T.Spending.Max))

#To compute all the results shown in section 3.2 (KS-Test) it
#is necessary to separate the fdata list by periods and
#classes. The user needs to remember that the test is for
#different classes (by pairs) from the same period. 

#Here we show an example for compare EFOS vs PCS in the
#1st period (2013-2018). The same code works for the 
#remaining comparisons. 

#Separate fdata by periods (Only 1st period is shown)
dataEPN=fdata %>% filter(Year==2013 | Year==2014 | 
                           Year==2015 | Year==2016 | 
                           Year==2017 | Year==2018 )
dataEPN = subset(dataEPN,select=-c(Year))

#Separate each period by classes 
#(Only EFOS and PCS classes are shown)
dataEPN.efos = dataEPN %>% filter(Status == "EFOS")
dataEPN.pcs = dataEPN %>% filter(Status == "PCS")

#Since variables of type i) are related to the fraction of contracts
#that satisfy certain property, and variables of type ii) are related
#to the fraction of buyers with certain features, it is necessary
#to make a list of the buyers and their features per class
dataEPN.efosDC = subset(dataEPN.efos, 
                        select = c(DepID,T.Cont.Max, T.Spending.Max))
dataEPN.pcsDC = subset(dataEPN.pcs, 
                       select = c(DepID,T.Cont.Max, T.Spending.Max))
dataEPN.efosDC = distinct(dataEPN.efosDC)
dataEPN.pcsDC = distinct(dataEPN.pcsDC)

#Run the KS Test for the period and classes selected

#Array to save the statistical D for each non-dummy variable
ks.D = vector("numeric", length = ncol(dataEPN)-3)
#Array to save the name of each non-dummy variable
ks.names = vector("character", length = ncol(dataEPN)-3)
#Array to save the computed p-vale for each non-dummy variable
ks.pv = vector("numeric", length = ncol(dataEPN)-3)

for(j in 4:ncol(dataEPN)){
  if(4 <= j & j <= 6){
    ks.D[j-3] = ks.test(dataEPN.efos[,j],
                        dataEPN.pcs[,j], simulate.p.value = T)$statistic
    ks.names[j-3] = names(dataEPN)[j]
    ks.pv[j-3] = ks.test(dataEPN.efos[,j],
                         dataEPN.pcs[,j], simulate.p.value = T)$p.value
  } else{
      ks.D[j-3] = ks.test(dataEPN.efosDC[,names(dataEPN)[j]],
                          dataEPN.pcsDC[,names(dataEPN)[j]], simulate.p.value = T)$statistic
      ks.names[j-3] = names(dataEPN)[j]
      ks.pv[j-3] = ks.test(dataEPN.efosDC[,names(dataEPN)[j]],
                           dataEPN.pcsDC[,names(dataEPN)[j]], simulate.p.value = T)$p.value
  }
}

#Make a data frame with the obtained results
ks.EPN.e_p = data.frame(ks.names, ks.D, ks.pv, stringsAsFactors = F) 
#Arranged from the biggest D to the smaller
ks.EPN.e_p = ks.EPN.e_p %>% arrange(desc(ks.D))
#Take only those dummy non-variables from which the 
#statistical D is larger than 0.1
#and p_value < 0.05
ks.EPN.ep.diff  = which(ks.EPN.e_p$ks.D >= 0.1 & ks.EPN.e_p$ks.pv < 0.05)

ks.EPN.e_p[ks.EPN.ep.diff,]

#The next section code is for generate the figures related to 
#those variables which passed the KS-Test.

for (l in ks.EPN.ep.diff){
  
      l.main = paste("EFOS vs PCS", "\\n D=", scientific(ks.EPN.e_p $ks.D[l], digits = 3),
                 sep = "") #main title in the graph
    
      # axis' range
      x.min = min(environment(ecdf(dataEPN.efos[,ks.EPN.e_p $ks.names[l]]))$x, 
                  environment(ecdf(dataEPN.pcs[,ks.EPN.e_p $ks.names[l]]))$x)
      x.max = max(environment(ecdf(dataEPN.efos[,ks.EPN.e_p $ks.names[l]]))$x, 
                  environment(ecdf(dataEPN.pcs[,ks.EPN.e_p $ks.names[l]]))$x)
      
      y.min = min(environment(ecdf(dataEPN.efos[,ks.EPN.e_p $ks.names[l]]))$y, 
                  environment(ecdf(dataEPN.pcs[,ks.EPN.e_p $ks.names[l]]))$y)
      y.max = max(environment(ecdf(dataEPN.efos[,ks.EPN.e_p $ks.names[l]]))$y, 
                  environment(ecdf(dataEPN.pcs[,ks.EPN.e_p $ks.names[l]]))$y)
      
      # constant interpolation of the minor sample to compare it with the biggest
      y.comp = KS_interpolate(environment(ecdf(dataEPN.pcs[,ks.EPN.e_p $ks.names[l]]))$x,
                              environment(ecdf(dataEPN.efos[,ks.EPN.e_p $ks.names[l]]))$x,
                              environment(ecdf(dataEPN.efos[,ks.EPN.e_p $ks.names[l]]))$y)
      
      #differences between CDFs
      Dif = abs(y.comp-environment(ecdf(dataEPN.pcs[,ks.EPN.e_p $ks.names[l]]))$y)
      #maximum difference between CDFs
      Dif.max = max(Dif)
      #x value related to the maximum difference
      Dif.x = which(Dif == Dif.max)
      
      #plot
      name.g = paste("EFOSvsPCS-",ks.EPN.e_p$ks.names[l],"-CDF-1st.png", sep = "")
      png(name.g, width = 5, height = 5, units = "in", res = 300)
      
      #This re-scale is necessary to a better look of statistical D
      if(ks.EPN.e_p $ks.names[l] == "Spending"){
        x.max = 1e+05
      }
      if(ks.EPN.e_p $ks.names[l] == "EBWeeks"){
        x.max = 100
      }
      
      plot(ecdf(dataEPN.efos[,ks.EPN.e_p $ks.names[l]]), verticals = T, col = "red", 
           main = l.main, xlim = c(x.min, x.max),
           ylim = c(y.min, y.max), lty = 1, pch = 1, cex.main = 1.25,
           xlab = ks.EPN.e_p $ks.names[l], ylab = "CDF", col.01line = NULL, do.points = T,
           las = 1, cex.lab = 1.25, cex.axis = 1.25)
      lines(ecdf(dataEPN.pcs[,ks.EPN.e_p $ks.names[l]]), verticals = T,
            col = "blue",lty = 2, pch = 2, cex = 0.75, col.01line = NULL, do.points = T)
      #paint the statistical D (maximum difference between CDFs)
      abline(v = environment(ecdf(dataEPN.pcs[,ks.EPN.e_p $ks.names[l]]))$x[Dif.x], col = "grey",
             lty = 4, lwd = 2)
      abline(h = y.comp[Dif.x], col = "grey",
             lty = 4, lwd = 2)
      abline(h = environment(ecdf(dataEPN.pcs[,ks.EPN.e_p $ks.names[l]]))$y[Dif.x], col = "grey",
             lty = 4, lwd = 2)
      
      if(ks.EPN.e_p $ks.names[l] == "Spending"){
        legend(1e04,1.0, legend = c("EFOS", "PCS"), col = c("red", "blue"),
               lty = c(1,2), pch = c(1,2), bty = "n", cex = 1.0)
      }
      
      dev.off()
}
