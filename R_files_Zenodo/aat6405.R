######################################################################
###
### Data analysis for article "Strong Impacts of Biodiversity
### in a Large-Scale Subtropical Forest Experiment" by Huang et al.
###

### Mixed models were fitted using ASReml
### (https://www.vsni.co.uk/software/asreml-r/)
### and summarized using function test.asreml provided by one
### of the corresponding authors.
###
### Since ASReml is commercial software, example output is provided
### for the first analysis in the form of comments.

library(asreml)

if(!require(pascal)) {
    library(devtools)
    install_github("pascal-niklaus/pascal/pascal")
    library(pascal)
}

######################################################################
###
### Read data file

d <- read.csv("aat6405_data_plot.csv") # taken from .xlsx file

d$Y <- as.factor(d$year)
d$logSR <- log2(d$div)

######################################################################
###
### Table 1: Basal area

ba.asr <- asreml(ba ~ site
                    + logSR
                    + Y
                    + Y:site
                    + Y:logSR,
                 random = ~comm:idh(site)
                    + mu4.plot:idh(site)
                    + siteplot
                    + Y:comm:idh(site)
                    + Y:mu4.plot:idh(site),
                 keep.order = TRUE,
                 data = d)
test.asreml(ba.asr)

## ---- Wald tests:
##             Df denDF F.inc      Pr
## (Intercept)  1   120 305.7 < 2e-16 ***
## site         1   120  14.3 0.00024 ***
## logSR        1   112   7.4 0.00739 **
## Y            4   489 309.0 < 2e-16 ***
## Y:site       4   488   7.7 4.7e-06 ***
## Y:logSR      4   456  15.2 1.1e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
##
## ---- Variance components:
##                        gamma component std.error z.ratio constraint
## comm:site!site.A       8.675    38.317     9.304    4.12   Positive
## comm:site!site.B       2.474    10.928     5.240    2.09   Positive
## siteplot!siteplot.var  3.669    16.204     1.715    9.45   Positive
## mu4.plot:site!site.A   0.963     4.253     3.805    1.12   Positive
## mu4.plot:site!site.B   1.744     7.701     4.467    1.72   Positive
## Y:comm:site!site.A     0.998     4.409     0.657    6.71   Positive
## Y:comm:site!site.B     0.563     2.486     0.545    4.56   Positive
## Y:mu4.plot:site!site.A 0.105     0.462     0.399    1.16   Positive
## Y:mu4.plot:site!site.B 0.199     0.881     0.444    1.98   Positive
## R!variance             1.000     4.417     0.221   20.01   Positive
##
## ---- Dispersion:
## 2.1

######################################################################
###
### Table 1: Basal area increment

dba.asr <- asreml.nvc(deltaba ~ site
                         + logSR
                         + Y
                         + Y:site
                         + Y:logSR,
                      random =~ comm:idh(site)
                         + mu4.plot:idh(site)
                         + siteplot
                         + Y:comm:idh(site)
                         + Y:mu4.plot:idh(site),
                      keep.order = TRUE,
                      ginits=c(.8,.4,.4,-.1,.2,.1,.1,.1,-.1,1),
                      control = asreml.control(maxiter = 100),
                      data = d)
test.asreml(dba.asr)

######################################################################
###
### Analysis including shrub effects and plot size

d2 <- subset(d, pool %in% c("A1","B1"))

d2$div_shrub <- factor(d2$div.shrub)
d2$year_scale <- d2$year-2015
d2$lgshrub <- ifelse(d2$div.shrub > 0, log2(d2$div.shrub), 0)

d2$shrub <- factor(ifelse(d2$div.shrub == 0, "noshrub", "shrub"))

d2.asr<-asreml(volume ~ site
                  + logSR
                  + shrub
                  + plot.size
                  + lgshrub
                  + year_scale
                  + logSR:(shrub + plot.size + lgshrub)
                  + year_scale:(logSR + shrub+plot.size + lgshrub),
               random =~ comm:idh(site)
                  + mu4.plot:idh(site)
                  + siteplot
                  + year_scale:comm:idh(site)
                  + year_scale:mu4.plot:idh(site)
                  + year_scale:siteplot,
               keep.order = TRUE,
               data=d2)
test.asreml(d2.asr)

######################################################################
###
### Analysis of Fig. 2 data: Net effects (NE),
### complementarity effects (CE), and selection effects (SE)

d3 <- subset(d, div > 1)
d3$year_scale <- d3$year - mean(d3$year)
d3$CE <- sqrt(abs(d3$CE)) * sign(d3$CE)
d3$SE <- sqrt(abs(d3$SE)) * sign(d3$SE)
d3$NE <- sqrt(abs(d3$NE)) * sign(d3$NE)

set.seed(54321)
d3.asr <- asreml.nvc(CE ~ site
                        + logSR
                        + year_scale
                        + year_scale:logSR,
                     random =~ comm:idh(site)
                        + mu4.plot:idh(site)
                        + siteplot
                        + year_scale:(comm:idh(site) + mu4.plot:idh(site) + siteplot),
                     keep.order = TRUE,
                     control = asreml.control(maxiter = 100),
                     data = d3)
test.asreml(d3.asr)

d3.asr <- asreml.nvc(SE ~ site
                        + logSR
                        + year_scale
                        + year_scale:logSR,
                     random =~ comm:idh(site)
                        + mu4.plot:idh(site)
                        + siteplot
                        + year_scale:(comm:idh(site) + mu4.plot:idh(site) + siteplot),
                     keep.order = TRUE,
                     control = asreml.control(maxiter = 100),
                     data=d3)
test.asreml(d3.asr)

d3.asr <- asreml.nvc(NE ~ site
                        + logSR
                        + year_scale
                        + year_scale:logSR,
                     random =~ comm:idh(site)
                        + mu4.plot:idh(site)
                        + siteplot
                        + year_scale:(comm:idh(site) + mu4.plot:idh(site) + siteplot),
                     keep.order = TRUE,
                     control = asreml.control(maxiter = 100),
                     data=d3)
test.asreml(d3.asr)

### End of file
###
######################################################################
