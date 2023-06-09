## Covid clinical trial meta-analysis
## cyp v-2022

library("ggplot2")
library("waffle")
library("dplyr")

# XXX set working directory
wd <- "~/Desktop/oct_dataset"
setwd(wd)
# load data
data <- read.csv("./oct-dataset_for-manuscript_2022-08-03.csv",header=T)

### Paper t-tests
## ICMJE member journals vs not for citations
# t = -5.3182, df = 104.35, p-value = 6.007e-07 (mean no: 41.98152, mean yes: 617.79048 ) - S
t.test(as.numeric(citations_dimensions)~icmje_member,data=data)
## ICMJE member journals vs not for altmetrics
# t = -6.3867, df = 106.79, p-value = 4.51e-09 (mean no: 306.1039, mean yes: 3419.8491 ) - S
t.test(as.numeric(altmetric_score)~icmje_member,data=data)
## ICMJE recommendations (i.e. follow?) vs not for citations
# t = -4.7315, df = 275.35, p-value = 3.567e-06 (mean no: 45.38148, mean yes: 264.15299) - S
t.test(as.numeric(citations_dimensions)~icmje_follow,data=data)
## ICMJE recommendations (i.e. follow?) vs not for altmetrics
# t = -5.37, df = 306.75, p-value = 1.558e-07 (mean no: 317.3444, mean yes: 1521.7993) - S
t.test(as.numeric(altmetric_score)~icmje_follow,data=data)
## Preprints found vs not for citations
# t = -0.55581, df = 198.6, p-value = 0.579 (mean no: 146.8297, mean yes: 178.7323) - NS
t.test(as.numeric(citations_dimensions)~preprint_found,data=data)
## Preprints found vs not for altmetrics
# t = -0.41217, df = 345.58, p-value = 0.6805 (mean no: 897.6156, mean yes: 985.3672) - NS
t.test(as.numeric(altmetric_score)~preprint_found,data=data)


## Paper figures
data_avail_request <- paste(data$data_avail_statement,data$data_avail_request_author,sep="_")
# A) to photoshop legend
# availability statement? no, yes --> available by request? no, yes
waffle(table(data_avail_request),colors=c("#D3D3D3","#708090","#DA614E"),xlab="1 square = 1 publication",rows=20,size=0.5)

# B)
data_avail_repo <- paste(data$data_avail_statement,data$data_avail_repo,sep="_")
# availability statement? no, yes --> available by repo? no, yes
waffle(table(data_avail_repo),colors=c("#D3D3D3","#708090","#4EC7DA"),xlab="1 square = 1 publication",rows=20,size=0.5)

