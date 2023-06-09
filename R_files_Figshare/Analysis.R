library('MASS')
library('pscl')
library('dummies')
library('lmtest')
library('stargazer')
library('gdata')
library('dplyr')
library('rsq')
library('Hmisc')
library('car')
library(ggplot2)
library("readxl")
library('caret')
## Variables Operationalization
##

setwd("~/Dropbox/CADS/Instagram Data")
Data = readRDS('Data.rds')

color = read.csv("color_variation.csv",header = FALSE, sep = ",")
#luminance2 = read.csv("luminance.csv",header = FALSE, sep = ",")
luminance =read.csv("luminance_entropy.csv",header=FALSE, sep = ",")
lc = t(luminance)
se = read.csv("subband_entropy.csv",header=FALSE, sep = ",")
fc = read.csv("feature_congestion.csv",header=FALSE, sep = ",")

luminance2 = read.csv("/Users/GijsOvergoor/Downloads/luminance2.csv", header =FALSE, sep =",")

image_name = read.table("image_name_new.txt")

lumi = cbind(image_name,luminance2)
colnames(lumi) = c("imageId","lumi2")


index = rep(0, dim(Data)[1])

for (i in 1: dim(Data)[1]){
  print(i)
  index[i] = which(Data$imageId[i] == newdata$imageId)
}

df = cbind(Data,newdata[index,-1])

#saveRDS(df,'Data_20180131.rds')

newdata = cbind(image_name,color[,-1],lc,se[,-1],fc[,-1])

colnames(newdata) = c("imageId","color","luminance","subband entropy","feature congestion")



df = left_join(Data,lumi,by = "imageId")

df = readRDS('Data_20180131.rds')
im_features=read_excel('/Users/gijso/Downloads/VCMatlab/image_features.xlsx')


for (i in 1:10){
  t = paste0('/Users/gijso/Downloads/Instagram/complexity/Objects/',i,".csv")
  tmp = read.csv(t,header=FALSE)
  tmp[,1]=gsub(".jpg","",tmp[,1])
  if (i == 1) {objects =tmp}
  else { 
  objects = rbind(objects,tmp)
}}

controls1 = read_excel('/Users/gijso/Downloads/VCMatlab/Concept Controls/anp.xlsx')
controls1$id = gsub("-biconcept.mat","",controls1$id)

controls2 = read_excel('/Users/gijso/Downloads/VCMatlab/Concept Controls/places_faces.xlsx')
controls2$id = gsub(".jpg","",controls2$id)


colnames(im_features)[c(32,33,35)] = c('CE','SC','id')

im_features$id = gsub(".jpg","",im_features$id)

N = dim(im_features)[1]

indicator = rep(0,N)

for (i in 1:N){
  tmp = which(df$imageId==im_features$id[i])
  if (length(tmp)>0){
  indicator[i] = length(tmp)}
}

im_features = im_features[-which(indicator==0),]
controls2 = controls2[-which(indicator==0),]


N = dim(df)[1]
m = rep(0,N)

for (i in 1:N){
  m[i] = which(controls1$id==df$imageId[i])
}
controls1 = controls1[m,]

m2 = rep(0,N)

for (i in 1:N){
  m2[i] = which(objects[,1]==df$imageId[i])
}
new_objects = objects[m2,]





im_features$m2[which(is.na(im_features$m2))]=0
im_features$m11[which(is.na(im_features$m11))]=0

 
controls= cbind(controls1[,-11],controls2[,-12])


allfeatures = cbind(df[,c("Image Size","threshold.0.01","Edges","color","luminance","subband entropy","feature congestion")],im_features[,-c(34,35)],new_objects[,2])
colnames(allfeatures)[41] = "object_rcnn"
allfeatures$fc2 = allfeatures$`feature congestion`-allfeatures$irv/0.0269

allfeatures$luminance=allfeatures$luminance/-10000000
allfeatures$`Image Size`=allfeatures$`Image Size`/1000000

allfeatures$lcv=allfeatures$lcv*-1 
allfeatures$lcv2 = log(allfeatures$lcv2+1)




# df$co
# ## delete zeros ## 
# zeros = which(df$color==0)
# 
# df = df[-zeros,]
# 
# zeros2 = which(df$followedByCount==0)
# 
# df  = df [-zeros2,]
# 
# zeros3 = which(df$Edges==0)
# df = df[-zeros3,]
# 
# na4 = which(is.na(df$imageTagCount))
# df = df[-na4,]

Data = df
n = dim(Data)[1]

Data$luminance = Data$luminance/-1000000


## Cut off all fully black or fully white images (manually inspected) ##
tmp = sort(Data$`feature congestion`)[140]
selecter = which(Data$`feature congestion` <= tmp)

Data = Data[-selecter,]
allfeatures = allfeatures[-selecter,]
controls = controls[-selecter,]

selecter2 = which(allfeatures$SC>3)
Data= Data[-selecter2,]
allfeatures = allfeatures[-selecter2,]
controls =controls[-selecter2,]

selecter3 = which(allfeatures$m11>7)
Data= Data[-selecter3,]
allfeatures = allfeatures[-selecter3,]
controls =controls[-selecter3,]

## exclude brands with less than 52 posts ## 
brands = unique(Data$brandId)



for (i in 1:length(brands)){
  tmp =which(brands[i]==Data$brandId)
  if (length(tmp)<52) {
   Data = Data[-tmp,]
 allfeatures = allfeatures[-tmp,]
   controls = controls[-tmp,]
  
    }
}


brands = unique(Data$brandId)

# brand investigation of luminance

#brand_lumi = matrix(0,length(brands),4)



#for (i in 1:length(brands)){
#  tmp =which(brands[i]==Data$brandId)
#brand_lumi[i,1] = mean(Data$luminance[tmp])
#  brand_lumi[i,2] = sd(Data$luminance[tmp])
#  brand_lumi[i,3] = min(Data$luminance[tmp])
#  brand_lumi[i,4] = max(Data$luminance[tmp])
#}

#mean(Data$luminance)
#sd(Data$luminance)
#min(Data$luminance)
#max(Data$luminance)



#########


n = dim(Data)[1]
## Filter or No Filter, Binary variable ##
filter = rep(0,n)
tmps = Data$imageFilter[1]
tmp = which(Data$imageFilter!=tmps)
filter[tmp] = 1

## tags ##
tags = Data$imageTagCount

## likes ## 
likes = Data$imageLikeCount

## loglikes ##
loglikes = log(Data$imageLikeCount+1)

## comments ##
comments = Data$imageCommentCount

## log comments ##
logcomments = log(Data$imageCommentCount+1)

## textual sentiment ##
positive = Data$Positive
negative = Data$Negative*-1

## image dimensions ##
H = Data$imageHeight
W = Data$imageWidth

H612 = which(H == 612)
H480 = which(H == 480)
H320 = which(H ==320)

dimensions = rep(1,n)
dimensions[H612] = 2
dimensions[H480] = 3
dimensions[H320] = 4

dimension = dummy(dimensions)
colnames(dimension) = c('640x640','612x612','480x480','320x320')

## followers ##
followers = log(Data$followedByCount)

## number of posts ##
posts = log(Data$postedMedia)

## weekend ##
weekend = Data$weekend

## season ## 
seasons = Data$season


for (i in 1:n){
  if (Data$month[i] %% 3 == 0 ){
    if (Data$day[i] < 22 ){
      if (seasons[i] > 1){
        seasons[i] = seasons[i]-1
      }else {
        seasons[i] = 4
      }
    }
  }
}

season = dummy(seasons)
colnames(season) = c('winter','spring','summer','fall')


## time of day ##
times = Data$timeofday
time = dummy(times)  

colnames(time) = c('morning','afternoon','evening','night')

## rush hour ## 
rush = rep(0,dim(Data)[1])
for (i in 1:dim(Data)[1]){
  if (Data$hour[i]<20 && Data$hour[i]>=17){
    if(Data$weekend[i]==0){
      rush[i]=1
    }
  }
}


## Feature Complexity ## 
jpeg = Data$`Image Size`/1000

## Feature Complexity ##
edges = Data$Edges*100
color = Data$color
luminance = Data$luminance

## High-Level complexity ##
quantities = Data$threshold.0.01
#tmp = which(quantities>12)
#quantity = rep(0,n)
#quantity[tmp] = 1
quantity = quantities

dissimilarities = (1- Data$Similarity)*100
#tmps = which(dissimilarities>40.37)
#dissimilarity=rep(0,n)
#dissimilarity[tmps]= 1
dissimilarity = dissimilarities

fc = Data$`feature congestion`
se = Data$`subband entropy`

## Normalization
normalization <- function(x) {
  (x - mean(x))/sd(x)
}
norm= function(x){ (x-min(x))/(max(x)-min(x))}







allfeatures$unique=allfeatures[,2]*Data$Similarity




for (i in 1:dim(allfeatures)[2]){
  allfeatures[,i] = norm(allfeatures[,i])
}

res = round(cor(allfeatures),2)

write.csv(res,'correlation_matrix.csv')


allfeatures=readRDS("allfeatures.rds")
likes=readRDS("likes.rds")

og_complexity = allfeatures[,c(2:4,7,34,35)]
og_complexity =cbind(og_complexity,og_complexity$Edges^2,og_complexity$color^2)



object_count = df[,c(3,31:34,43)]
write.csv(object_count,file="object_count.csv")

photography = allfeatures[,c(8:21)]
photography = photography[,-10] #brightness and clarity are too strongly correlated 

feature_complexity = allfeatures[,c(3,4,5,34,36,37,38)]
design_complexity = allfeatures[c(39,26:29,35,42,41)]                                 


new_complexity = allfeatures[,c(7,23:26,34,35)]
new_complexity = cbind(new_complexity,new_complexity$lcv2^2,new_complexity$edv^2,new_complexity$ccv^2,new_complexity$`feature congestion`^2)
M_complexity = allfeatures[,c(1,3,30:38)]
colnames(M_complexity)[c(6)]= c("m7")
for (i in 1:dim(M_complexity)[2]){
  M_complexity = cbind(M_complexity,M_complexity[,i]^2)
  colnames(M_complexity)[i+11] = paste("m",i,"-squared",sep="")
}



colnames(M_complexity)= c("M7 - Jpeg File Size" , "M6 - Edge Density",  "M1 - Contrast" , "M2 - Correlation"  ,"M3 - Energy"  ,"M4 - Homogeneity" , "M5 - Frequency Factor" , "M8 - Region Count" , "M9 - Colorfulness" , "M10 - Number of Colors" ,"M11 - Color Harmony")


M_complexity = M_complexity[,-c(3:4,6:8)]




fc =cbind(feature_complexity[,-c(2,6,7)],feature_complexity$Edges^2,feature_complexity$luminance^2, feature_complexity$m5^2,feature_complexity$m9^2)
assym = (design_complexity$ahv+design_complexity$avv)/2
dc = cbind(design_complexity[,c(5,6,8)],assym, assym^2,design_complexity$irv^2, design_complexity$object_rcnn^2,design_complexity$m8^2)



ALLN = as.data.frame(cbind(likes,fc,dc,photography,controls,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))

ALLN2 = as.data.frame(cbind(loglikes,fc,dc,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))

ALL2 = as.data.frame(cbind(comments,fc,dc,photography,controls,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))


saveRDS(ALLN,'allvars_regression.rds')


#,feature_complexity[,-c(2,5)] set up for feature complexity ( Edges, luminance, number of colors or colorfulness)
#fc = cbind(feature_complexity[,-c(2,5)],feature_complexity$Edges^2,feature_complexity$luminance^2, feature_complexity$m5^2, feature_complexity$m10^2, feature_complexity$m11^2)

fc = (fc$Edges+fc$luminance+fc$m9)/3
fc2 = fc^2
dc = (dc$object_rcnn+dc$irv+dc$assym)/3
dc2 = dc^2


negb= glm.nb(likes~.,data=ALLN[l,],maxit = 1000) 
summary(negb)
rsq(negb,adj=TRUE)



stargazer(negb,summary=TRUE,single.row=TRUE,keep.stat =c("n","adj.rsq"), initial.zero = FALSE )


# correlation matrix complexity variables #

M_complexity = M_complexity[c(3,4,5,6,7,1,2,8,9,10,11)]
colnames(dc) = c("Irregularity of OA", "Region Count", "Objects","Asymmilarity of OA")
colnames(fc) = c("Edge Density", "Luminance", "Frequency Factor","Color")
clutter = allfeatures[,c(6,7)]
colnames(clutter) = c("Clutter - feature congestion", "Clutter - subband entropy")

dc = dc

corplotter(as.data.frame(cbind(fc,dc,M_complexity,clutter)))

## Visualizations ##

beta_color = coef(negb)[c("m9", "`feature_complexity$m9^2`")]
beta_luminance = coef(negb)[c("luminance", "`feature_complexity$luminance^2`")]
beta_edge = coef(negb)[c("Edges", "`feature_complexity$Edges^2`")]

beta_objects = coef(negb)[c("object_rcnn", "`design_complexity$object_rcnn^2`")]
beta_ir = coef(negb)[c("irv", "`design_complexity$irv^2`")]
beta_as = coef(negb)[c("assym", "`assym^2`")]

fun_color = function(x)  beta_color[1]*x+beta_color[2]*x^2
fun_luminance =function(x)  beta_luminance[1]*x+beta_luminance[2]*x^2
fun_edges =function(x) beta_edge[1]*x+beta_edge[2]*x^2
fun_objects =function(x) beta_objects[1]*x+beta_objects[2]*x^2
fun_irv =function(x) beta_ir[1]*x+beta_ir[2]*x^2
fun_as =function(x) beta_as[1]*x+beta_as[2]*x^2


max_colorx= beta_color[1]/(-2*beta_color[2])
max_colory = fun_color(max_colorx)


plot(fun_color,yaxt='n', ylab="",xlab = "")
plot(fun_luminance,yaxt='n', ylab="",xlab = "")
plot(fun_edges,yaxt='n', ylab="",xlab = "")
plot(fun_objects,yaxt='n', ylab="",xlab = "")
plot(fun_irv,yaxt='n', ylab="",xlab = "")
plot(fun_as,yaxt='n', ylab="",xlab = "")



###

  
p = ggplot(data=data.frame(ALLN[,c(5,3,2,13,12,14)]), mapping = aes(x= object_rcnn))  

p+stat_function(fun = fun_objects) + xlim(0,1) + ylim(-0.1,0.1)# + geom_histogram( aes(y = stat(count/sum(count))), bins=30)






colnames(ALLN[,c(5,3,2,13,12,14)])


poisson = glm(likes~.,data=ALLN, family="poisson",maxit=1000)
summary(poisson)



logLik(poisson)


feature_complex = allfeatures$`Image Size`
feature_complex_squared =feature_complex^2
assym = (design_complexity$ahv+design_complexity$avv)/2
design_complex = (dc$object_rcnn+dc$irv+dc$assym)/3
design_complex_squared = design_complex^2
data_m1 = as.data.frame(cbind(likes,feature_complex,design_complex,photography,controls,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))
data_m2 = as.data.frame(cbind(likes,feature_complex,feature_complex_squared,design_complex,design_complex_squared,photography,controls,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))

fc =cbind(feature_complexity[,-c(2,6)])
assym = (design_complexity$ahv+design_complexity$avv)/2
dc = cbind(design_complexity[,c(5,6,8)],assym)

data_m3= as.data.frame(cbind(likes,fc,dc,photography,controls,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))

fc =cbind(feature_complexity[,-c(2,6)],feature_complexity$Edges^2,feature_complexity$luminance^2, feature_complexity$m5^2,feature_complexity$m9^2, feature_complexity$m11^2)
assym = (design_complexity$ahv+design_complexity$avv)/2
dc = cbind(design_complexity[,c(5,6,8)],assym, assym^2,design_complexity$irv^2, design_complexity$object_rcnn^2, design_complexity$m8^2)

data_m4 = as.data.frame(cbind(likes,fc,dc,photography,controls,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))

m1= glm.nb(likes~.,data=data_m1,maxit = 1000) 
m2= glm.nb(likes~.,data=data_m2,maxit = 1000) 
m3= glm.nb(likes~.,data=data_m3,maxit = 1000) 
m4= glm.nb(likes~.,data=data_m4,maxit = 1000) 

summary(m3)

stargazer(m1,m3,m2,m4,summary=TRUE,single.row=TRUE,keep.stat =c("n","adj.rsq"),omit=c("posts","tags","afternoon","evening","night","weekend","spring","summer","fall"  ), initial.zero = FALSE )

data_m11 = as.data.frame(cbind(likes,fc,dc,followers, posts,positive,negative,tags))
data_m12 = as.data.frame(cbind(likes,fc,dc,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))
data_m13 = as.data.frame(cbind(likes,fc,dc,photography,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))
data_m14 =as.data.frame(cbind(likes,fc,dc,photography,controls,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))

m11= glm.nb(likes~.,data=data_m11,maxit = 1000) 
m12= glm.nb(likes~.,data=data_m12,maxit = 1000)   
m13= glm.nb(likes~.,data=data_m13,maxit = 1000) 
m14= glm.nb(likes~.,data=data_m14,maxit = 1000) 


rsq(negb,adj=TRUE)

stargazer(m11,m12,m13,m14,summary=TRUE,single.row=TRUE,keep.stat =c("n","adj.rsq"), initial.zero = FALSE )  


## Correlation Table ## 

fc =cbind(feature_complexity[,-c(2,6,7)])
assym = (design_complexity$ahv+design_complexity$avv)/2
dc = cbind(design_complexity[,c(5,6,8)],assym)


ALLN = as.data.frame(cbind(likes,fc,dc,followers, positive,negative))

cnames = c( "Likes"  ,"Edge Density","Luminance","Frequency Factor","Color","Edges2","Luminance2","Frequency Factor2","Color2","Irregularity", "Region Count",        "Objects" ,    "Asymmetry" , "Asymmetry2",   "Irregularity2", "Region Count2",        "Objects2" ,     "Diagonal Dominance"     ,  "Rule of Thirds"    ,        "Vertical Physical Dominance"      ,     "Horizontal Physical Dominance"    ,       "Horizontal Color Balance"   ,  "Vertical Color Balance"   ,    "FG Size Difference"   ,    "FG Color Difference"      ,  "FG Texture Difference"    ,   "Saturation" ,    "Contrast"     ,  "Clarity"  ,      "Warmth"     ,    "Crazy_car"    ,  "classic_castle",   "hot_girls"     ,  "outdoor_party"  ,"busy_office"  ,   "amazing_food"  ,  "hot_cup"   ,      "cute_animals",   "outdoor_wedding" ,"favorite_team" ,  "art_studio"  ,    "bakery/shop" ,    "beach"   ,        "clean_room"  ,    "coffee_shop" ,    "desert/sand"  ,   "museum/indoor" ,  "nursery"  ,       "ocean"         ,  "playroom"     ,  "face"      ,      "Followers" ,      "Posts"     ,      "Text Sentiment Positive"  ,      "Text Sentiment Negative"    ,    "Hashtags"       ,    "afternoon"  ,     "evening" ,        "night"   ,        "weekend"    ,     "spring"       ,   "summer"    ,     "fall")
colnames(ALLN)=cnames


source("~/Dropbox/CADS/Instagram Data/correlation_plot.R")
cnames =c("Likes","Edge Density","Luminance","Frequency Factor","Color","Irregularity of OA","Region Count","Objects","Asymmetry of OA","Followers","Text Sentiment Positive","Text Sentiment Negative")
ALLN2 = ALLN[,c(1,5,3,2,4,8,6,9,7,10,11,12)]
corplotter(ALLN2)

medium = ALLN[,c(2:22)]

corplotter(ALLN)



stargazer(res1, type = "latex")

brand = Data[,c("brandName","industry")]
brand = brand[!duplicated(brand$brandName),]


write.csv(brand,'/Users/gijso/Downloads/brand.csv')
## Specified Datasets ##

## Subsets ## 

industries = unique(Data$industry)

AM = which(Data$industry==industries[2]) # automotive
FP = which(Data$industry==industries[7]) # Fashion & Personal Care
FS = which(Data$industry==industries[9]) # retail

SLT = which(Data$industry==industries[4]) # sports, leisure, travel


negb= glm.nb(likes~.,data=ALLN[AM,],maxit = 1000) 
summary(negb)


for (i in 1:length(industries)){
  print(i)
  print(length(which(Data$industry==industries[i])))
}


## 
filters = Data %>% select(imageFilter) %>% count(imageFilter)

sorted = sort(as.numeric(filters$n),decreasing = T)

filters_top6 = filters[which(filters$n>sorted[7]),]

filters_top = unique(Data$imageFilter)

filters_top = filters_top[c(1,6,8,11,3,4)]

filter_avg = matrix(NA,6,8)
vc = cbind(fc[,c(1:5)],dc[,c(1:3)])

for (i in 1:6){
  tmp = which(filters_top[i]==Data$imageFilter)
  for (j in 1:8){
        filter_avg[i,j] = mean(vc[tmp,j])
        
  }
}






## optimal value based on the coefficients  ##
optimal_lumi = -full2$coefficient[[2]]/(2*full2$coefficient[[3]])
optimal_edge = -full2$coefficient[[4]]/(2*full2$coefficient[[5]])


newALLN = ALLN

## Set optimal filter values ## 
l_up = 1.15
l_down = 0.75
e_up = 1.12
e_down = 1.1

## replace values ## 
for (i in 1:148066){
  if(l_up* ALLN[i,2]<= optimal_lumi){
    newALLN[i,2]= ALLN[i,2]*l_up}
  else if(l_down* ALLN[i,2] >= optimal_lumi){
    newALLN[i,2]= ALLN[i,2]*l_down
    print(2)
  } else { newALLN[i,2] = optimal_lumi
  print(3)}
}




for (i in 1:148066){
  if(e_up* ALLN[i,4]<= optimal_edge){
    newALLN[i,4]= ALLN[i,4]*e_up}
  else if(e_down* ALLN[i,4] >= optimal_edge){
    newALLN[i,4]= ALLN[i,4]*e_down
  } else { newALLN[i,4] = optimal_edge}
  print(i)
}


newALLN[,3] = newALLN[,2]^2
newALLN[,5] = newALLN[,4]^2


X = newALLN[,-1]
yhat_opti =exp( as.matrix(cbind(rep(1,148066),X)) %*% full2$coefficients)
yhat = full2$fitted.values

avg_improvement =sum(yhat_opti-yhat)/sum(yhat)





## Prediction ## 


fc =cbind(feature_complexity[,-c(2,6)],feature_complexity$Edges^2,feature_complexity$luminance^2, feature_complexity$m5^2,feature_complexity$m9^2, feature_complexity$m11^2)
assym = (design_complexity$ahv+design_complexity$avv)/2
dc = cbind(design_complexity[,c(5,8)],assym, assym^2,design_complexity$irv^2, design_complexity$object_rcnn^2)

fc = allfeatures$`Image Size`
dc = (dc$object_rcnn+dc$irv+dc$assym)/3

dc=dc$object_rcnn


newALLN = as.data.frame(cbind(likes,fc,dc,photography,controls,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))


#,feature_complexity[,-c(2,5)] set up for feature complexity ( Edges, luminance, number of colors or colorfulness)
#fc = cbind(feature_complexity[,-c(2,5)],feature_complexity$Edges^2,feature_complexity$luminance^2, feature_complexity$m5^2, feature_complexity$m10^2, feature_complexity$m11^2)


n=dim(newALLN)[1]



set.seed(123)
for (i in 1:5){

train_ind = sample(seq_len(n), size = floor(0.7*n))

train = newALLN[train_ind,]
test = newALLN[-train_ind,]


model = glm.nb(likes~.,data=train,maxit = 1000) 
summary(model)

yhat = exp(predict(model,test))

print(sqrt(mean((yhat-test[,1])^2)))
print(cor(yhat,test[,1],method = "spearman"))
print(sqrt(mean((yhat[-which(test[,1]>10000)]-test[-which(test[,1]>10000),1])^2)))
}





##############

## create optimal feature complexity data ## 

constant = rep(1,n)
opti_lumi = rep(optimal_lumi,n)
opti_lumi2 = rep(optimal_lumi^2,n)
opti_edge = rep(optimal_edge,n)
opti_edge2 = rep(optimal_edge^2,n)  

X = cbind(constant,opti_lumi,opti_lumi2, opti_edge,opti_edge2, ALLN[,-c(1:5)])

yhat_opti =exp( as.matrix(X) %*% full$coefficients)
yhat = full$fitted.values

avg_improvement =sum(yhat_opti-yhat)/sum(yhat)


#################

## brand fixed-effects instead of activity and followers (too much for my computer)

ALLN = as.data.frame(cbind(likes,nluminance,nluminance2, nedges, nedges2,nquantity,ndissimilarity,ninteraction,nfc,followers, posts,positive,negative,tags,time[,-1],weekend,season[,-1]))

brands = Data$brandId
brand_names = unique(Data$brandName)

brand_mat = dummy(brands)

colnames(brand_mat) = brand_names

newALLN = as.data.frame(cbind(likes,ncolor,ncolor2,nluminance,nluminance2, nedges, nedges2,nquantity,ndissimilarity,ninteraction,nfc,positive,negative,tags,time[,-1],weekend,season[,-1],brand_mat[,-1]))

full= glm.nb(likes~.,data=newALLN,maxit = 100) 
summary(full)
rsq(full,adj=TRUE)


#################

fm <- list("Basic Model" = basic, "Extended Model" = extended , "Full Model" = full)
round(sapply(fm, function(x) coef(x)[1:10]), digits = 8)

## testing poisson vs. neg binomial m3 ##

P = glm(comments~.,data=,family = poisson,maxit = 1000)
NB = glm.nb(likes~.,data=ALL[tmp,],maxit = 1000)
NB2 = glm.nb(comments~.,data= ALL2[tmp,], maxit= 1000)
## H = hurdle(likes ~ ., data = LAH1, dist = "negbin",maxit = 1000)
## ZI = zeroinfl(comments ~ ., data = CSM, dist = "negbin", EM = TRUE, maxit=1000)

## All full models ##
NBLL = glm.nb(likes~.,data=LAL,maxit = 1000)
NBLM = glm.nb(likes~.,data=LAM,maxit = 1000)
NBLH1 = glm.nb(likes~.,data=LAH1,maxit = 1000)
NBLH2 = glm.nb(likes~.,data=LAH2,maxit = 1000)

NBCL = glm.nb(comments~.,data=CAL,maxit = 1000)
NBCM = glm.nb(comments~.,data=CAM,maxit = 1000)
NBCH1 = glm.nb(comments~.,data=CAH1,maxit = 1000)
NBCH2 = glm.nb(comments~.,data=CAH2,maxit = 1000)

stargazer(full,title = "Influence of Low-Level Complexity on consumer engagement",align= TRUE,dep.var.labels = c("Likes","Comments"),covariate.labels = c("File Size","File Size Squared","Followers","Posts","Textual Sentiment - Positive","Textual Sentiment - Negative","Tags","Filter","Afternoon","Evening","Night","Weekend","Spring","Summer","Fall"),type="latex",digits = 3,no.space =TRUE)

stargazer(NBLM,NBCM,title = "Influence of Mid-Level Complexity on consumer engagement",align= TRUE,dep.var.labels = c("Likes","Comments"),covariate.labels = c("Edge Percentage","Followers","Posts","Textual Sentiment - Positive","Textual Sentiment - Negative","Tags","Filter","Afternoon","Evening","Night","Weekend","Spring","Summer","Fall"),type="latex",digits = 3,no.space =TRUE)

stargazer(NBLH1,NBCH1,title = "Influence of High-Level Complexity on consumer engagement",align= TRUE,dep.var.labels = c("Likes","Comments"),covariate.labels = c("Object Quantity","Followers","Posts","Textual Sentiment - Positive","Textual Sentiment - Negative","Tags","Filter","Afternoon","Evening","Night","Weekend","Spring","Summer","Fall"),type="latex",digits = 3,no.space =TRUE)

stargazer(NBLH2,NBCH2,title = "Influence of High-Level Complexity on consumer engagement",align= TRUE,dep.var.labels = c("Likes","Comments"),covariate.labels = c("Object Dissimilarity","Followers","Posts","Textual Sentiment - Positive","Textual Sentiment - Negative","Tags","Filter","Afternoon","Evening","Night","Weekend","Spring","Summer","Fall"),type="latex",digits = 3,no.space =TRUE)

## Show comparison of 3 datasets ##
NBLL1  = glm.nb(likes~.,data=LSL,maxit = 1000)
NBLL2 = glm.nb(likes~.,data=LEL,maxit = 1000)


NB = glm.nb(likes~.,data=ALL,maxit = 1000)
NB2 = glm.nb(comments~.,data= ALL2, maxit= 1000)



stargazer(NB,NB2,title = "All measures of complexity in one single model",dep.var.labels = c("Likes","Comments"),covariate.labels = c("File Size","File Size Squared","Edges","Quantity","Dissimilarity","Followers","Posts","Caption Positive","Caption Negative","Tags","Filter","Afternoon","Evening","Night","Weekend","Spring","Summer","Fall"),type="text",no.space =TRUE,keep.stat=c("rsq"))
NB_SLT = glm.nb(likes~.,data=ALL[SLT,],maxit = 1000)
NB2_SLT = glm.nb(comments~.,data= ALL2[SLT,], maxit= 1000)


NB_FS = glm.nb(likes~.,data=ALL[FS,],maxit = 1000)
NB2_FS = glm.nb(comments~.,data= ALL2[FS,], maxit= 1000)
?'glm.nb'

NB_FP = glm.nb(likes~.,data=ALL[FP,],maxit = 1000)
NB2_FP = glm.nb(comments~.,data= ALL2[FP,], maxit= 1000)

stargazer(NB_SLT,NB2_SLT,title = "Sports, Leisure and Travel",dep.var.labels = c("Likes","Comments"),covariate.labels = c("File Size","File Size Squared","Edges","Quantity","Dissimilarity","Followers","Posts","Caption Positive","Caption Negative","Tags","Filter","Afternoon","Evening","Night","Weekend","Spring","Summer","Fall"),type="latex",no.space =TRUE,keep.stat=c("rsq"))
stargazer(NB_FS,NB2_FS,title = "Food",dep.var.labels = c("Likes","Comments"),covariate.labels = c("File Size","File Size Squared","Edges","Quantity","Dissimilarity","Followers","Posts","Caption Positive","Caption Negative","Tags","Filter","Afternoon","Evening","Night","Weekend","Spring","Summer","Fall"),type="latex",no.space =TRUE,keep.stat=c("rsq"))
stargazer(NB_FP,NB2_FP,title = "Fashion",dep.var.labels = c("Likes","Comments"),covariate.labels = c("File Size","File Size Squared","Edges","Quantity","Dissimilarity","Followers","Posts","Caption Positive","Caption Negative","Tags","Filter","Afternoon","Evening","Night","Weekend","Spring","Summer","Fall"),type="latex",no.space =TRUE,keep.stat=c("rsq"))

1 - NB_SLT$deviance/NB_SLT$null.deviance
1 - NB2_SLT$deviance/NB2_SLT$null.deviance

1 - NB_FS$deviance/NB_FS$null.deviance
1 - NB2_FS$deviance/NB2_FS$null.deviance


1 - NB_FP$deviance/NB_FP$null.deviance
1 - NB2_FP$deviance/NB2_FP$null.deviance

brands = unique(Data$brandName)
posts = rep(0,n)

for (i in 1:n){
  tmp = which(Data$brandName==Data$brandName[i])
  posts[i] = length(tmp)
}

i = 20
### Descriptive Statistics ### 

descriptives = cbind(sapply(ALL, mean, na.rm=TRUE),sapply(ALL,sd,na.rm=TRUE),sapply(ALL,min,na.rm=TRUE),sapply(ALL,max,na.rm=TRUE))
descriptives

brands = Data[!duplicated(Data$brandId),]
mean(brands$followedByCount)
min(brands$followedByCount)
max(brands$followedByCount)
sd(brands$followedByCount)

## visualization


v = df$threshold.0.01*(1- df$Similarity)
l =which(v== sort(v)[25000])                # lowest, in increasing order.
h = which(v == sort(v, decreasing=TRUE)[25000])
df$imageId[l]
df$imageId[h]

ids = Data$imageId[which(Data$imageFilter=="Normal")]
