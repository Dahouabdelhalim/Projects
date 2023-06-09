###Analysis R Script.R
##File associated with
##EM Caves, PA Green, MN Zipple, D Bharath, S Peters, S Johnsen, and S Nowicki “Comparison of categorical color perception in two Estrildid finches.” The American Naturalist.
##Author(s) for correspondence: Stephen Nowicki (snowicki@duke.edu) and Eleanor Caves (eleanor.caves@gmail.com)


library(dplyr)
library(tidyr)

agg_mean<-read.csv("Zebra Finch Luminance Training Data.csv")

####Generate Model of Luminance Dataset Performance (i.e. definitely did not cross 5/6 boundary)####

expected<-lm(I(pass.freq-1/15)~0+contrast+chrom.dist,data=subset(agg_mean,col1>5|col2<6))
expected_cont_only<-lm(I(pass.freq-1/15)~0+contrast,data=subset(agg_mean,col1>5|col2<6))
expected_chrom_only<-lm(I(pass.freq-1/15)~0+chrom.dist,data=subset(agg_mean,col1>5|col2<6))
summary(expected)
AIC(expected,expected_chrom_only,expected_cont_only)

agg_mean$exp<-(1/15+.534*agg_mean$contrast+.0168*agg_mean$chrom.dist)

####How well does the training model predict pass frequency over a different range?####

agg_caves<-read.csv("Zebra Finch Beak Data.csv")

agg_caves$exp<-(1/15+.534*agg_caves$contrast+.0167*agg_caves$chrom.dist)
agg_caves$across5<-ifelse(agg_caves$col1<6&agg_caves$col2>=6,1,0)

summary(lm(pass.freq~exp*across5,agg_caves))

###########Plot Zebra Finch Results#######

plot(pass.freq~exp,agg_caves,xlim=c(0,0.90),ylim=c(0,0.90),type="n",xlab="Expected Pass Frequency",ylab="Observed Pass Frequency",main="Zebra Finches")
points(pass.freq~exp,subset(agg_caves,across5==1),xlim=c(0,0.8),ylim=c(0,0.8),pch=19,cex=1)
points(pass.freq~exp,subset(agg_caves,across5==0),xlim=c(0,0.8),ylim=c(0,0.8),pch=19,col="blue",cex=1)
curve(expr=1*x,1/15,1,add = T,lty=2)
curve(expr=0.09+0.44*x,0.1,0.4,lty=2,lwd=3,col="blue",add=T)
curve(expr=.05+0.993*x,0.25,0.85,lty=2,lwd=3,col="black",add=T)


####Bengalese Finches dataset####
beng_agg<-read.csv("Bengalese Finch Data.csv")

####Bengalese Finch Analysis####

expected<-lm(I(pass.freq-1/15)~0+contrast+chrom.dist,subset(beng_agg,training==1))
expected_chrom_only<-lm(I(pass.freq-1/15)~0+chrom.dist,subset(beng_agg,training==1))
expected_cont_only<-lm(I(pass.freq-1/15)~0+contrast,subset(beng_agg,training==1))
AIC(expected,expected_chrom_only,expected_cont_only)

beng_agg$exp<-1/15+0.97*beng_agg$contrast+0.022*beng_agg$chrom.dist

summary(lm(pass.freq~exp,subset(beng_agg,training==1)))

beng_agg$across5<-ifelse((beng_agg$col1<6&beng_agg$col2<6)|(beng_agg$col1>5&beng_agg$col2<9&beng_agg$col1<9),0,
                        ifelse(beng_agg$col1<6&beng_agg$col2>5&beng_agg$col2<9,1,NA))

summary(lm(pass.freq~exp,subset(beng_agg,across5==1)))
summary(lm(pass.freq~exp,subset(beng_agg,across5==0)))
summary(lm(pass.freq~exp*across5,beng_agg))

#combined Figure
par(mfrow=c(1,2))
plot(pass.freq~exp,subset(beng_agg),type="n",xlim=c(0,1),ylim=c(0,1),xlab="Expected Pass Frequency",ylab="Observed Pass Frequency",main="Bengalese Finches")
points(pass.freq~exp,subset(beng_agg,col1<6&col2>5&col2<9),pch=19,cex=1)
points(pass.freq~exp,subset(beng_agg,(col1<6&col2<6)|(col1>5&col2<9&col1<9)),pch=19,col="blue",cex=1)

curve(expr=-0.04+0.86*x,0.1,0.4,lty=2,lwd=3,col="blue",add=T)
curve(expr=-0.03+0.85*x,0.3,0.9,lty=2,lwd=3,col="black",add=T)
curve(expr=1*x,1/15,0.9,lty=2,lwd=1,col="black",add=T)

plot(pass.freq~exp,agg_caves,xlim=c(0,0.90),ylim=c(0,0.90),type="n",xlab="Expected Pass Frequency",ylab="Observed Pass Frequency",main="Zebra Finches")
points(pass.freq~exp,subset(agg_caves,across5==1),xlim=c(0,0.8),ylim=c(0,0.8),pch=19,cex=1)
points(pass.freq~exp,subset(agg_caves,across5==0),xlim=c(0,0.8),ylim=c(0,0.8),pch=19,col="blue",cex=1)
curve(expr=1*x,1/15,1,add = T,lty=2)
curve(expr=0.09+0.44*x,0.1,0.4,lty=2,lwd=3,col="blue",add=T)
curve(expr=.05+0.993*x,0.25,0.85,lty=2,lwd=3,col="black",add=T)

##############################



######The code below generates the results in Table 1###
beng_agg_price<-subset(beng_agg,col2<9&col1<9)
beng_agg_price$across5<-ifelse(beng_agg_price$col1<6&beng_agg_price$col2>5,1,0)

bf_cat_mod1<-lm(pass.freq~across5,subset(beng_agg_price,col2-col1==1))
bf_con_mod1<-lm(pass.freq~contrast,subset(beng_agg_price,col2-col1==1))

AIC(bf_cat_mod1,bf_con_mod1)
summary(bf_cat_mod1)
summary(bf_con_mod1)

bf_cat_mod2<-lm(pass.freq~across5,subset(beng_agg_price,col2-col1==2))
bf_con_mod2<-lm(pass.freq~contrast,subset(beng_agg_price,col2-col1==2))

AIC(bf_cat_mod2,bf_con_mod2)
summary(bf_cat_mod2)
summary(bf_con_mod2)


bf_cat_mod3<-lm(pass.freq~across5,subset(beng_agg_price,col2-col1==3))
bf_con_mod3<-lm(pass.freq~contrast,subset(beng_agg_price,col2-col1==3))

AIC(bf_cat_mod3,bf_con_mod3)
summary(bf_cat_mod3)
summary(bf_con_mod3)

bf_cat_mod_1_3<-lm(pass.freq~across5+chrom.dist,subset(beng_agg_price,col2-col1<=3))
bf_con_mod_1_3<-lm(pass.freq~contrast+chrom.dist,subset(beng_agg_price,col2-col1<=3))
summary(bf_cat_mod_1_3)
summary(bf_con_mod_1_3)
AIC(bf_cat_mod_1_3,bf_con_mod_1_3)

###
zf_agg_price<-agg_caves

zf_cat_mod1<-lm(pass.freq~across5,subset(zf_agg_price,col2-col1==1))
zf_con_mod1<-lm(pass.freq~contrast,subset(zf_agg_price,col2-col1==1))

AIC(zf_cat_mod1,zf_con_mod1)
summary(zf_cat_mod1)
summary(zf_con_mod1)

zf_cat_mod2<-lm(pass.freq~across5,subset(zf_agg_price,col2-col1==2))
zf_con_mod2<-lm(pass.freq~contrast,subset(zf_agg_price,col2-col1==2))

AIC(zf_cat_mod2,zf_con_mod2)
summary(zf_cat_mod2)
summary(zf_con_mod2)


zf_cat_mod3<-lm(pass.freq~across5,subset(zf_agg_price,col2-col1==3))
zf_con_mod3<-lm(pass.freq~contrast,subset(zf_agg_price,col2-col1==3))

AIC(zf_cat_mod3,zf_con_mod3)
summary(zf_cat_mod3)
summary(zf_con_mod3)

zf_cat_mod_1_3<-lm(pass.freq~across5+chrom.dist,subset(zf_agg_price,col2-col1<=3))
zf_con_mod_1_3<-lm(pass.freq~contrast+chrom.dist,subset(zf_agg_price,col2-col1<=3))
summary(zf_cat_mod_1_3)
summary(zf_con_mod_1_3)
AIC(zf_cat_mod_1_3,zf_con_mod_1_3)

