### JAE: Episodes of high recruitment buffer against climate-driven mass mortality events in a North Pacific seabird population
### 4/11/2020
## Mike Johns

require(reshape)
require(RMark)
require(ggplot2)
require(devtools)
require(scales)
require(dplyr)
require(gridExtra)

capt.histories <- read.csv("enter file location here")
covariates <- read.csv("enter file location here")

##### set up the prdel model ===========================================================================
sub.red <- subset(capt.histories, year1 > 1988)
#sub.red <- subset(sub, year > 1987) # was 82 (37) 88 = 31

junk<-melt(sub.red,id.var=c("id","year"),measure.var="recap")
y=cast(junk,id ~ year, mean)
#fill in all the days when the animal wasn't seen with zeros
y[is.na(y)]=0
y

x <- y
x$sum <- rowSums(x)
x[x$sum==0,]

# across all columns: check to see if there are any issues with the data
y %>% filter_all(any_vars(. %in% c(0.5,1.5,2.5,3.5,4)))

#function to read the matrix and create the capture history strings
pasty<-function(x) 
{
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {
    y<-(x[i,])
    out[i]<-paste(y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10],y[11],
                  y[12],y[13],y[14],y[15],y[16],y[17],y[18],y[19],y[20],y[21],y[22],y[23],
                  y[24],y[25],y[26],y[27],y[28],y[29],y[30],y[31],sep="")
  }
  return(out)
}

#capture history data frame for RMark  
capt.hist<-data.frame(ch=pasty(y[,2:32]))  

pradel.processed=process.data(capt.hist,model="Pradrec",begin.time = 1989)

pradel.ddl=make.design.data(pradel.processed)

### add covariates ###
pradel.ddl$Phi=merge_design.covariates(pradel.ddl$Phi,covariates)
pradel.ddl$p=merge_design.covariates(pradel.ddl$p,covariates)
pradel.ddl$f=merge_design.covariates(pradel.ddl$f,covariates)

#pradel survival and recruitment
run.ms=function(){
  # phi and p time variant 
  Phi.t=list(formula= ~ -1 + time)# model survival SST
  
  p.t=list(formula= ~ -1 + time)

  ## null models for f
  f.t=list(formula= ~ -1 + time)
  f.null=list(formula= ~ -1)
  
  ## lag 3 =====================================================================================================
  # no environment
  f.lag3.1 = list(formula = ~ occ)
  f.lag3.2 = list(formula = ~ db.lag3)
  f.lag3.3 = list(formula = ~ occ + db.lag3)
  f.lag3.4 = list(formula = ~ occ * db.lag3)
  
  # with sst
  f.lag3.5 = list(formula = ~ occ*db.lag3 + sst.sp)
  f.lag3.6 = list(formula = ~ occ + db.lag3*sst.sp)
  f.lag3.7 = list(formula = ~ occ*sst.sp + db.lag3)
  f.lag3.8 = list(formula = ~ occ + db.lag3 + sst.sp)
  f.lag3.9 = list(formula = ~ db.lag3 + sst.sp)
  f.lag3.10 = list(formula = ~ occ + sst.sp)
  f.lag3.11 = list(formula = ~ db.lag3*sst.sp)
  f.lag3.12 = list(formula = ~ occ*sst.sp)
  f.lag3.13 = list(formula = ~ sst.sp)
  
  # with beuti
  f.lag3.14 = list(formula = ~ occ*db.lag3 + beuti.sp)
  f.lag3.15 = list(formula = ~ occ + db.lag3*beuti.sp)
  f.lag3.16 = list(formula = ~ occ*beuti.sp + db.lag3)
  f.lag3.17 = list(formula = ~ occ + db.lag3 + beuti.sp)
  f.lag3.18 = list(formula = ~ db.lag3 + beuti.sp)
  f.lag3.19 = list(formula = ~ occ + beuti.sp)
  f.lag3.20 = list(formula = ~ db.lag3*beuti.sp)
  f.lag3.21 = list(formula = ~ occ*beuti.sp)
  f.lag3.22 = list(formula = ~ beuti.sp)
  
  # with soi.sp
  f.lag3.41 = list(formula = ~ occ*db.lag3 + soi.sp)
  f.lag3.42 = list(formula = ~ occ + db.lag3*soi.sp)
  f.lag3.43 = list(formula = ~ occ*soi.sp + db.lag3)
  f.lag3.44 = list(formula = ~ occ + db.lag3 + soi.sp)
  f.lag3.45 = list(formula = ~ db.lag3 + soi.sp)
  f.lag3.46 = list(formula = ~ db.lag3*soi.sp)
  f.lag3.46.5 = list(formula = ~ soi.sp*occ)
  f.lag3.46.6 = list(formula = ~ soi.sp+occ)
  f.lag3.47 = list(formula = ~ soi.sp)
  
  
  ## lag 4 ============================================================================
  # no environment
  f.lag4.1 = list(formula = ~ db.lag4)
  f.lag4.2 = list(formula = ~ occ + db.lag4)
  f.lag4.3 = list(formula = ~ occ * db.lag4)
  
  # with sst
  f.lag4.4 = list(formula = ~ occ*db.lag4 + sst.sp)
  f.lag4.5 = list(formula = ~ occ + db.lag4*sst.sp)
  f.lag4.6 = list(formula = ~ occ*sst.sp + db.lag4)
  f.lag4.7 = list(formula = ~ occ + db.lag4 + sst.sp)
  f.lag4.8 = list(formula = ~ db.lag4 + sst.sp)
  f.lag4.9 = list(formula = ~ db.lag4*sst.sp)
  
  # with beuti
  f.lag4.11 = list(formula = ~ occ*db.lag4 + beuti.sp)
  f.lag4.12 = list(formula = ~ occ + db.lag4*beuti.sp)
  f.lag4.13 = list(formula = ~ occ*beuti.sp + db.lag4)
  f.lag4.14 = list(formula = ~ occ + db.lag4 + beuti.sp)
  f.lag4.15 = list(formula = ~ db.lag4 + beuti.sp)
  f.lag4.16 = list(formula = ~ db.lag4*beuti.sp)
  
  # with soi.sp
  f.lag4.32 = list(formula = ~ occ*db.lag4 + soi.sp)
  f.lag4.33 = list(formula = ~ occ + db.lag4*soi.sp)
  f.lag4.34 = list(formula = ~ occ*soi.sp + db.lag4)
  f.lag4.35 = list(formula = ~ occ + db.lag4 + soi.sp)
  f.lag4.36 = list(formula = ~ db.lag4 + soi.sp)
  f.lag4.37 = list(formula = ~ db.lag4*soi.sp)
  
  
  
  ms.model.list=create.model.list("Pradrec")
  ms.results=mark.wrapper(ms.model.list, data=pradel.processed,ddl=pradel.ddl)
  
}

## best model
#run.ms=function(){
  # phi and p time variant 
#  Phi.t=list(formula= ~ -1 + time)# model survival SST
#  p.t=list(formula= ~ -1 + time)
#  f.lag4.12 = list(formula = ~ occ + db.lag4*beuti.sp)
  
 # ms.model.list=create.model.list("Pradrec")
#  ms.results=mark.wrapper(ms.model.list, data=pradel.processed,ddl=pradel.ddl)
  
#}

ms.results=run.ms()

# without c-hat correction
ms.results

## code for adjusting model ranking by increments of c-hat
adjust.chat(chat=2,ms.results) # AIC selection corrected

## top model information
best.model <- "Phi.t.p.t.f.lag4.12"
ms.results$Phi.t.p.t.f.lag4.12$results$beta

# derived lambda
l.t <- c(1,ms.results$Phi.t.p.t.f.lag4.12$results$derived$`Lambda Population Change`$estimate)
up.l <- c(1,ms.results$Phi.t.p.t.f.lag4.12$results$derived$`Lambda Population Change`$ucl)
down.l <- c(1,ms.results$Phi.t.p.t.f.lag4.12$results$derived$`Lambda Population Change`$lcl)

## geometric mean of lambda
vals <- l.t[2:31]
steps <- 1/length(vals)
geo <- sum(vals) * steps

### psuedo R2 for top ranking models #########################################
run.ms.var=function(){
  # phi and p time variant 
  #Phi.t=list(formula= ~ time)# model survival SST
  Phi.t=list(formula= ~ -1 + time)# model survival SST
  p.t=list(formula= ~ -1 + time)
  f.lag4.dot = list(formula = ~ 1)
  f.lag4.1 = list(formula = ~ occ + db.lag4*beuti.sp)
  f.lag4.2 = list(formula = ~ occ + db.lag3 + beuti.sp)
  f.lag4.3 = list(formula = ~ occ + db.lag3*beuti.sp)
  f.lag4.4 = list(formula = ~ occ*db.lag3 + beuti.sp)
  f.lag4.5 = list(formula = ~ occ*beuti.sp + db.lag3)
  f.lag4.6 = list(formula = ~ occ + db.lag3)
  f.lag4.7 = list(formula = ~ occ + beuti.sp)
  f.lag4.8 = list(formula = ~ occ)
  f.lag4.9 = list(formula = ~ beuti.sp)
  f.lag4.10 = list(formula = ~ -1 + time)
  f.lag4.11 = list(formula = ~ db.lag4)
  
  ms.model.list=create.model.list("Pradrec")
  ms.results=mark.wrapper(ms.model.list, data=pradel.processed,ddl=pradel.ddl)
  
}

var.results=run.ms.var()

null.p <- var.results$Phi.t.p.t.f.lag4.dot$results$deviance
v1 <- var.results$Phi.t.p.t.f.lag4.1$results$deviance
v2 <- var.results$Phi.t.p.t.f.lag4.2$results$deviance
v3 <- var.results$Phi.t.p.t.f.lag4.3$results$deviance
v4 <- var.results$Phi.t.p.t.f.lag4.4$results$deviance
v5 <- var.results$Phi.t.p.t.f.lag4.5$results$deviance
v6 <- var.results$Phi.t.p.t.f.lag4.6$results$deviance
v7 <- var.results$Phi.t.p.t.f.lag4.7$results$deviance
v8 <- var.results$Phi.t.p.t.f.lag4.8$results$deviance
v9 <- var.results$Phi.t.p.t.f.lag4.9$results$deviance
time <- var.results$Phi.t.p.t.f.lag4.10$results$deviance
v11 <- var.results$Phi.t.p.t.f.lag4.11$results$deviance

r1 <- (null.p - v1) / (null.p - time)
r2 <- (null.p - v2) / (null.p - time)
r3 <- (null.p - v3) / (null.p - time)
r4 <- (null.p - v4) / (null.p - time)
r5 <- (null.p - v5) / (null.p - time)
r6 <- (null.p - v6) / (null.p - time)
r7 <- (null.p - v7) / (null.p - time)
r8 <- (null.p - v8) / (null.p - time)
r9 <- (null.p - v9) / (null.p - time)
r10 <- (null.p - time) / (null.p - time)
r11 <- (null.p - v11) / (null.p - time)

## variance explained by each parameter in top model ###########################
run.ms=function(){
  # phi and p time variant 
  #Phi.t=list(formula= ~ time)# model survival SST
  Phi.t=list(formula= ~ -1 + time)# model survival SST
  
  p.t=list(formula= ~ -1 + time)
  #p.1=list(formula= ~ soi.sp.thres + Time + I(Time^2))
  
  ## null models for f
  f.t=list(formula= ~ time)
  f.dot=list(formula= ~ 1)
  f.occ= list(formula = ~ occ)
  f.db= list(formula = ~ db.lag4)
  f.beuti= list(formula = ~ beuti.sp)
  f.x= list(formula = ~ beuti.sp*db.lag4)
  
  ms.model.list=create.model.list("Pradrec")
  ms.results=mark.wrapper(ms.model.list, data=pradel.processed,ddl=pradel.ddl)
  
}
ms.results=suppressWarnings(run.ms())
ms.results

beuti <- ms.results$Phi.t.p.t.f.beuti$results$deviance
occ <- ms.results$Phi.t.p.t.f.occ$results$deviance
db <- ms.results$Phi.t.p.t.f.db$results$deviance
x <- ms.results$Phi.t.p.t.f.x$results$deviance
null <- ms.results$Phi.t.p.t.f.dot$results$deviance
t <- ms.results$Phi.t.p.t.f.t$results$deviance

beuti.dev <- (null - beuti) / (null - t)
occ.dev <- (null - occ) / (null - t)
db.dev <- (null - db) / (null - t)
x.dev <- (null - x) / (null - t)

variances <- data.frame(value = c(beuti.dev,occ.dev,db.dev,x.dev),
                        covariate = c("upwell","occ","db","upwell X db"))

### model selection and R2 squared for survival and recruitment
run.surv=function(){
  Phi.t=list(formula= ~  time)# 
  Phi.dot=list(formula= ~ 1 )# 
  Phi.soi=list(formula= ~ soi.win )# 
  Phi.soi2=list(formula= ~ soi.win + I(soi.win^2) )# 
  Phi.sst=list(formula= ~ sst.winter )# 
  Phi.sst2=list(formula= ~ sst.winter + I(sst.winter^2) )#
  Phi.mei=list(formula= ~ mei.winter )# 
  Phi.mei2=list(formula= ~ mei.winter + I(mei.winter^2) )#
  Phi.beuti=list(formula= ~ beuti.winter )# 
  Phi.beuti2=list(formula= ~ beuti.winter + I(beuti.winter^2) )#
  
  p.t=list(formula= ~  time)
  #p.1=list(formula= ~ soi.sp.thres + Time + I(Time^2))
  
  ## null models for f
  f.t=list(formula= ~ time)
  
  ms.model.list=create.model.list("Pradrec")
  ms.results=mark.wrapper(ms.model.list, data=pradel.processed,ddl=pradel.ddl)
  
}
surv.results=suppressWarnings(run.surv())
surv.results

null.surv <- surv.results$Phi.dot.p.t.f.t$results$deviance
t.surv <- surv.results$Phi.t.p.t.f.t$results$deviance
beuti.surv <- surv.results$Phi.beuti.p.t.f.t$results$deviance
beuti2.surv <- surv.results$Phi.beuti2.p.t.f.t$results$deviance
soi.surv <- surv.results$Phi.soi.p.t.f.t$results$deviance
soi2.surv <- surv.results$Phi.soi2.p.t.f.t$results$deviance
sst.surv <- surv.results$Phi.sst.p.t.f.t$results$deviance
sst2.surv <- surv.results$Phi.sst2.p.t.f.t$results$deviance
mei.surv <- surv.results$Phi.mei.p.t.f.t$results$deviance
mei2.surv <- surv.results$Phi.mei2.p.t.f.t$results$deviance

beuti.dev.surv <- (null.surv - beuti.surv) / (null.surv - t.surv)
beuti2.dev.surv <- (null.surv - beuti2.surv) / (null.surv - t.surv)
soi.dev.surv <- (null.surv - soi.surv) / (null.surv - t.surv)
soi2.dev.surv <- (null.surv - soi2.surv) / (null.surv - t.surv)
sst.dev.surv <- (null.surv - sst.surv) / (null.surv - t.surv)
sst2.dev.surv <- (null.surv - sst2.surv) / (null.surv - t.surv)
mei.dev.surv <- (null.surv - mei.surv) / (null.surv - t.surv)
mei2.dev.surv <- (null.surv - mei2.surv) / (null.surv - t.surv)

#####
run.p=function(){
  Phi.t=list(formula= ~ time)# 
  
  p.t=list(formula= ~  time)
  p.dot=list(formula= ~ 1 )# 
  p.soi=list(formula= ~ soi.win )# 
  p.soi2=list(formula= ~ soi.win + I(soi.win^2) )# 
  p.sst=list(formula= ~ sst.winter )# 
  p.sst2=list(formula= ~ sst.winter + I(sst.winter^2) )#
  p.mei=list(formula= ~ mei.winter )# 
  p.mei2=list(formula= ~ mei.winter + I(mei.winter^2) )#
  p.beuti=list(formula= ~ beuti.winter )# 
  p.beuti2=list(formula= ~ beuti.winter + I(beuti.winter^2) )#
  
  ## null models for f
  f.t=list(formula= ~ time)
  
  ms.model.list=create.model.list("Pradrec")
  ms.results=mark.wrapper(ms.model.list, data=pradel.processed,ddl=pradel.ddl)
  
}
p.results=suppressWarnings(run.p())
p.results

null.p <- p.results$Phi.t.p.dot.f.t$results$deviance
t.p <- p.results$Phi.t.p.t.f.t$results$deviance
beuti.p <- p.results$Phi.t.p.beuti.f.t$results$deviance
beuti2.p <- p.results$Phi.t.p.beuti2.f.t$results$deviance
soi.p <- p.results$Phi.t.p.soi.f.t$results$deviance
soi2.p <- p.results$Phi.t.p.soi2.f.t$results$deviance
sst.p <- p.results$Phi.t.p.sst.f.t$results$deviance
sst2.p <- p.results$Phi.t.p.sst2.f.t$results$deviance
mei.p <- p.results$Phi.t.p.mei.f.t$results$deviance
mei2.p <- p.results$Phi.t.p.mei2.f.t$results$deviance

beuti.dev.p <- (null.p - beuti.p) / (null.p - t.p)
beuti2.dev.p <- (null.p - beuti2.p) / (null.p - t.p)
soi.dev.p <- (null.p - soi.p) / (null.p - t.p)
soi2.dev.p <- (null.p - soi2.p) / (null.p - t.p)
sst.dev.p <- (null.p - sst.p) / (null.p - t.p)
sst2.dev.p <- (null.p - sst2.p) / (null.p - t.p)
mei.dev.p <- (null.p - mei.p) / (null.p - t.p)
mei2.dev.p <- (null.p - mei2.p) / (null.p - t.p)

sandp <- data.frame(value =c(beuti.dev.p,beuti2.dev.p,soi.dev.p,soi2.dev.p,sst.dev.p,sst2.dev.p,mei.dev.p,mei2.dev.p,
                             beuti.dev.surv,beuti2.dev.surv,soi.dev.surv,soi2.dev.surv,sst.dev.surv,sst2.dev.surv,mei.dev.surv,mei2.dev.surv),
                    covariate = rep(c("upwell","upwell2","soi","soi2","sst","sst2","mei","mei2"),2),
                    parameter = c(rep("detection",8),rep("survival",8)))
