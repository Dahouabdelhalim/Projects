install.packages('iNEXT')

library(iNEXT)
library(ggplot2)


### birds in DX
##################################################
birddata<-read.csv("DXBirdData.csv")
head(birddata)
out <- iNEXT(birddata,q=c(0,1), endpoint=NULL, nboot=50, conf = 0.95, datatype="abundance")
out

out$iNextEst

p1=ggiNEXT(out, se=T)  ## very messy as it plot three indices at 4 sites together
p1 + 
  theme_bw()+
  theme(text = element_text(size=18,face="bold"))+
  theme(panel.grid.major= element_blank(),
        panel.grid.minor = element_blank())

p2<-ggiNEXT(out, type=1, facet.var="order", color.var="site") ## a little bit better but still messy
p2

## The easiest and neatest way is to report the asymptotic diversity estimates and the statistics 
benthosdiv.site<-out$AsyEst

## give you a data frame of 7 variables:  "Site", "Diversity","Observed", "Estimator", "s.e.", "LCL", and "UCL"
write.csv(div.site, "bird-alpha-V2.csv")



########################################################################################################

###### fish in DX

fishdata<-read.csv("DXFishData.csv")
head(fishdata)
out <- iNEXT(fishdata,q=c(0,1), endpoint=NULL, nboot=50, conf = 0.95, datatype="abundance")


p1=ggiNEXT(out, se=T)  ## very messy as it plot three indices at 4 sites together
p1 + 
  theme_bw()+
  theme(text = element_text(size=18,face="bold"))+
  theme(panel.grid.major= element_blank(),
        panel.grid.minor = element_blank())

p2<-ggiNEXT(out, type=1, facet.var="order", color.var="site") ## a little bit better but still messy
p2

## The easiest and neatest way is to report the asymptotic diversity estimates and the statistics 
fish.site<-out$AsyEst

## give you a data frame of 7 variables:  "Site", "Diversity","Observed", "Estimator", "s.e.", "LCL", and "UCL"
write.csv(fishdiv.site, "fish-alpha-V2.csv")


##################################################################################################################

#### macrobenthos in DX

benthosdata<-read.csv("DXBenthosData.csv")
head(benthosdata)
out <- iNEXT(benthosdata,q=c(0,1), endpoint=NULL, nboot=50, conf = 0.95, datatype="abundance")## endpoint=NULL 默认是数量的2倍外推


p1=ggiNEXT(out, se=T)  ## very messy as it plot three indices at 4 sites together
p1 + 
  theme_bw()+
  theme(text = element_text(size=18,face="bold"))+
  theme(panel.grid.major= element_blank(),
        panel.grid.minor = element_blank())

p2<-ggiNEXT(out, type=1, facet.var="order", color.var="site") ## a little bit better but still messy
p2



## The easiest and neatest way is to report the asymptotic diversity estimates and the statistics 
benthosdiv.site<-out$AsyEst

## give you a data frame of 7 variables:  "Site", "Diversity","Observed", "Estimator", "s.e.", "LCL", and "UCL"
write.csv(benthosdiv.site, "benthos-alpha-V2.csv")


##################################################################################################################

#### mangrove in DX

mangrovedata<-read.csv("DXMangroveData.csv")
head(mangrovedata)
out <- iNEXT(mangrovedata,q=c(0,1), endpoint=NULL, nboot=50, conf = 0.95, datatype="abundance")

p1 + 
  theme_bw()+
  theme(text = element_text(size=18,face="bold"))+
  theme(panel.grid.major= element_blank(),
        panel.grid.minor = element_blank())

p2<-ggiNEXT(out, type=1, facet.var="order", color.var="site") ## a little bit better but still messy
p2




## The easiest and neatest way is to report the asymptotic diversity estimates and the statistics 
mangrovediv.site<-out$AsyEst

## give you a data frame of 7 variables:  "Site", "Diversity","Observed", "Estimator", "s.e.", "LCL", and "UCL"
write.csv(mangrovediv.site, "mangrove-alpha-V2.csv")

