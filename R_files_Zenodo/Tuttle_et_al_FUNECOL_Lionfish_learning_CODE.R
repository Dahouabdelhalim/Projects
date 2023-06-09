#######################################################################################################
##
## Analysis code for
## "Differential learning by native versus invasive predators to avoid distasteful cleaning mutualists"
## by Lillian J. Tuttle (tuttlel@hawaii.edu), Robert W. Lamb, and Allison Stringer
## in Functional Ecology 2021
##
#######################################################################################################



## Question: Is there a difference among the predator species in the proportion 
## that ate a cleaner goby?

## Answer (see below): No, there was no difference among the predator species 
## (lionfish, graysby, coney) in the proportion that ate a cleaner goby (p=0.760, 
## Fisher's exact test).

## Hungry PTVO: n=31 out of 42, hungry CECR: n=23 of 32, hungry CEFU: n=12 of 30
Input =(
  "Predator  AteELGE  DidntEatELGE
  PTVO       14        17
  CECR       11        12
  CEFU        4         8
  ")
Matriz = as.matrix(read.table(textConnection(Input),
                              header=TRUE, 
                              row.names=1))
Matriz
fisher.test(Matriz, alternative="two.sided")
#p-value = 0.7599
#alternative hypothesis: two.sided



## Question: Is the likelihood of a LIONFISH eating the cleaner goby affected by 
## the order of exposure to the prey species, and/or the ratio of lionfish TL:
## cleaner goby TL? (I did not include the region where the trial was conducted 
## because PTVO were hungry only at LSI. I did not include whether or not the 
## lionfish ate the COGL bc all hungry lionfish that ate ELGE also ate COGL.)

## Answer (see below): A lionfish was 8.3 times more likely to eat a cleaner goby for 
## every 1 unit increase in the ratio of lionfish to cleaner goby total lengths (95% 
## CI 1.4, 50.4, p=0.02, binary logistic regression).  There was no effect of the
## order of exposure to goby species on the likelihood of a lionfish eating a cleaner
## goby (p=0.43).

hungry.lf<-read.csv('/FILEPATH/Tuttle_et_al_FUNECOL_cleaner_hungry_PTVO.csv',header=T)
hungry.lf
ateelge.lf<-ifelse(hungry.lf$AteELGE=='y',1,0)
order.cogl.lf<-ifelse(hungry.lf$Order=='COGL',1,0)
ratio.TL.PTVO_ELGE<-as.vector(hungry.lf$TL/hungry.lf$ELGETL)
## Binary logistic regression model
mod.hungry.lf<-glm(ateelge.lf ~ order.cogl.lf + ratio.TL.PTVO_ELGE, data=hungry.lf,
                   family=binomial(link="logit")) 
summary(mod.hungry.lf)
#Coefficients:
#                  Estimate Std. Error z value Pr(>|z|)  
#(Intercept)         -7.4210     3.2942  -2.253   0.0243 *
#order.cogl.lf       -0.6647     0.8361  -0.795   0.4266  
#ratio.TL.PTVO_ELGE   2.1116     0.9223   2.289   0.0221 *
#Null deviance: 42.684  on 30  degrees of freedom
#Residual deviance: 34.470  on 28  degrees of freedom
#AIC: 40.47
## Calculating 95% CI for effect of TL ratio on likelihood of a lionfish eating ELGE
ci.lf<-c((2.1116-1.96*0.9223),2.1116,(2.1116+1.96*0.9223))
exp(ci.lf)
#[1]  1.355123  8.261449 50.365580
## Playing with empirical logit graphs
n=tapply(ateelge.lf, ratio.TL.PTVO_ELGE, length)
n
y=tapply(ateelge.lf, ratio.TL.PTVO_ELGE, sum)
elogits=log((y+0.5)/(n-y+0.5))
idx=c(2.5, 2.8421052631579, 2.88235294117647, 3.04761904761905, 3.07142857142857, 
      3.11111111111111, 3.2, 3.21052631578947, 3.26086956521739, 3.28571428571429, 
      3.33333333333333, 3.38461538461538, 3.4, 3.42857142857143, 3.47058823529412, 
      3.47368421052632, 3.52631578947368, 3.59090909090909, 3.85714285714286, 
      3.93333333333333, 3.9375, 3.94444444444444, 4, 4.15384615384615,
      4.41666666666667, 4.93333333333333, 5.38461538461538)
plot(idx,elogits, xlab="Ratio of lionfish to cleaner goby total lengths", 
     ylab="Empirical logit")



## Question: Is the likelihood of a GRAYSBY eating the cleaner goby affected by 
## the order of exposure to the prey species, the ratio of graysby TL:cleaner 
## goby TL, or the region where the trial was conducted? (I did not include 
## whether or not the graysby ate the COGL bc nearly all hungry graysby that ate 
## ELGE also ate COGL.)

## Answer (see below): There was no evidence to suggest that the likelihood of a 
## graysby eating a cleaner goby was affected by the order of exposure to goby
## species, the ratio of graysby to cleaner goby total lengths, or the region
## where the trial was conducted (all p>0.1, binary logistic regression).

hungry.cecr<-read.csv('/FILEPATH/Tuttle_et_al_FUNECOL_cleaner_hungry_CECR.csv',header=T)
hungry.cecr
ateelge.cecr<-ifelse(hungry.cecr$AteELGE=='y',1,0)
order.cogl.cecr<-ifelse(hungry.cecr$Order=='COGL',1,0)
loc<-ifelse(hungry.cecr$Locality=='LSI',1,0)
ratio.TL.CECR_ELGE<-as.vector(hungry.cecr$TL/hungry.cecr$ELGETL)
## Binary logistic regression model
mod.hungry.cecr<-glm(ateelge.cecr ~ order.cogl.cecr + ratio.TL.CECR_ELGE + loc, 
                     data=hungry.cecr, family=binomial(link="logit")) 
summary(mod.hungry.cecr)
#Coefficients:
#                   Estimate Std. Error z value Pr(>|z|)
#(Intercept)          3.0712     2.2079   1.391    0.164
#order.cogl.cecr     -1.3004     0.9805  -1.326    0.185
#ratio.TL.CECR_ELGE  -0.5545     0.4190  -1.323    0.186
#loc                  0.5851     0.9938   0.589    0.556
#Null deviance: 31.841  on 22  degrees of freedom
#Residual deviance: 28.845  on 19  degrees of freedom
#AIC: 36.845



## Question: Is the likelihood of a CONEY eating the cleaner goby affected by 
## the order of exposure to the prey species, the ratio of coney TL:cleaner 
## goby TL, or the region where the trial was conducted? (I did not include 
## whether or not the coney ate the COGL bc nearly all hungry coney that ate 
## ELGE also ate COGL.)

## Answer (see below): There was no evidence to suggest that the likelihood of a 
## coney eating a cleaner goby was affected by the order of exposure to goby
## species, the ratio of coney to cleaner goby total lengths, or the region
## where the trial was conducted (all p>0.1, binary logistic regression).

hungry.cefu<-read.csv('/FILEPATH/Tuttle_et_al_FUNECOL_cleaner_hungry_CEFU.csv',header=T)
hungry.cefu
ateelge.cefu<-ifelse(hungry.cefu$AteELGE=='y',1,0)
order.cogl.cefu<-ifelse(hungry.cefu$Order=='COGL',1,0)
loc2<-ifelse(hungry.cefu$Locality=='LSI',1,0)
ratio.TL.CEFU_ELGE<-as.vector(hungry.cefu$TL/hungry.cefu$ELGETL)
## Binary logistic regression
mod.hungry.cefu<-glm(ateelge.cefu ~ order.cogl.cefu + ratio.TL.CEFU_ELGE + loc2, 
                     data=hungry.cefu, family=binomial(link="logit")) 
summary(mod.hungry.cefu)
#Coefficients:
#                    Estimate Std. Error z value Pr(>|z|)
#(Intercept)           5.1738     3.8136   1.357    0.175
#order.cogl.cefu       0.3023     1.9765   0.153    0.878
#ratio.TL.CEFU_ELGE   -1.1453     0.9022  -1.269    0.204
#loc2                -21.5417  6304.4673  -0.003    0.997
#Null deviance: 15.2763  on 11  degrees of freedom
#Residual deviance:  6.8122  on  8  degrees of freedom
#AIC: 14.812



## Question: Is there a difference in the proportion of trials in which a LIONFISH 
## struck at a cleaner goby in the first trial versus the fourth trial? In the 
## first versus the third? In the first versus the second?

## Answer (see below): The proportion of trials in which a lionfish struck at a 
## cleaner goby declined significantly after only one trial (p<0.001, McNemar’s 
## exact test), and continued to decline over the course of 4 trials.

T1vT4 = as.table(rbind(c(5, 14),c(0, 0)))
colnames(T1vT4) <- rownames(T1vT4) <- c("Yes", "No")
names(dimnames(T1vT4)) = c("T1", "T4")
T1vT4
margin.table(T1vT4, 1)
margin.table(T1vT4, 2)
sum(T1vT4)
mcnemar.test(T1vT4, correct=FALSE)
#McNemar's chi-squared = 14, df = 1, p-value = 0.0001828
## If I didn't take the within-subjects nature of my data into account:
prop.test(rbind(margin.table(T1vT4, 1), margin.table(T1vT4, 2)), correct=FALSE)

T1vT3 = as.table(rbind(c(10, 11),c(0, 0)))
colnames(T1vT3) <- rownames(T1vT3) <- c("Yes", "No")
names(dimnames(T1vT3)) = c("T1", "T4")
T1vT3
margin.table(T1vT3, 1)
margin.table(T1vT3, 2)
sum(T1vT3)
mcnemar.test(T1vT3, correct=FALSE)
#McNemar's chi-squared = 13, df = 1, p-value = 0.0003115
## If I didn't take the within-subjects nature of my data into account:
prop.test(rbind(margin.table(T1vT3, 1), margin.table(T1vT3, 2)), correct=FALSE)

T1vT2 = as.table(rbind(c(14, 13),c(0, 0)))
colnames(T1vT2) <- rownames(T1vT2) <- c("Yes", "No")
names(dimnames(T1vT2)) = c("T1", "T4")
T1vT2
margin.table(T1vT2, 1)
margin.table(T1vT2, 2)
sum(T1vT2)
mcnemar.test(T1vT2, correct=FALSE)
#McNemar's chi-squared = 13, df = 1, p-value = 0.0003115
## If I didn't take the within-subjects nature of my data into account:
prop.test(rbind(margin.table(T1vT2, 1), margin.table(T1vT2, 2)), correct=FALSE)

#Alternative that allows for multiple comparisons among timepoints while 
#accounting for ind'l predator (n=19 lionfish that had 4 consec trials):
mydata <- data.frame(
  outcome = c(1,1,0,0, 1,0,0,0, 1,0,1,0, 1,0,1,0, 1,1,0,0, 
              1,1,1,1, 1,1,0,0, 1,0,1,1, 1,1,1,0, 1,1,1,1, 
              1,0,0,0, 1,0,0,0, 1,0,0,0, 1,0,0,0, 1,0,0,0, 
              1,1,1,1, 1,1,0,1, 1,1,0,0, 1,0,1,0),
  time = gl(4,1,76,labels=1:4),
  lionfish = gl(19,4,labels=letters[1:19])
)
mydata$outcome <- factor(
  mydata$outcome, levels = c(1, 0),
  labels = c("strike", "no strike")
)
# Cross-tabulation
xtabs(~outcome + time, mydata)
#           time
#outcome      1  2  3  4
#  strike    19  9  8  5
#  no strike  0 10 11 14

library(rstatix)
# Compare the proportion of strikes across timepoints
cochran_qtest(mydata, outcome ~ time|lionfish)
# A tibble: 1 x 6
#  .y.         n statistic    df         p  method          
#* <chr>   <int>     <dbl> <dbl>     <dbl>  <chr>           
#1 outcome    19      24.2     3 0.0000231  Cochran's Q test

# pairwise comparisons between timepoints
pairwise_mcnemar_test(mydata, outcome ~ time|lionfish)
# A tibble: 6 x 6
#  group1 group2        p   p.adj p.adj.signif method      
#* <chr>  <chr>     <dbl>   <dbl> <chr>        <chr>       
#1 1      2      0.00443  0.0266  *            McNemar test
#2 1      3      0.00257  0.0154  *            McNemar test
#3 1      4      0.000512 0.00307 **           McNemar test
#4 2      3      1        1       ns           McNemar test
#5 2      4      0.221    1       ns           McNemar test
#6 3      4      0.371    1       ns           McNemar test



## Question: Is there a difference in the proportion of GRAYSBY that struck at a
## cleaner goby in the first trial versus the fourth trial? In the first versus 
## the third? In the first versus the second?

## Answer (see below): The proportion of trials in which a graysby struck at a 
## cleaner goby declined over time (Fig. 3), but less precipitously than with 
## lionfish, and the decline was non-significant (p>0.05, McNemar’s exact test).

T1vT4.cecr = as.table(rbind(c(5, 3),c(0, 0)))
colnames(T1vT4.cecr) <- rownames(T1vT4.cecr) <- c("Yes", "No")
names(dimnames(T1vT4.cecr)) = c("T1", "T4")
T1vT4.cecr
margin.table(T1vT4.cecr, 1)
margin.table(T1vT4.cecr, 2)
sum(T1vT4.cecr)
mcnemar.test(T1vT4.cecr, correct=FALSE)
#McNemar's chi-squared = 3, df = 1, p-value = 0.08326
## If I didn't take the within-subjects nature of my data into account:
prop.test(rbind(margin.table(T1vT4.cecr, 1), margin.table(T1vT4.cecr, 2)), correct=FALSE)

T1vT3.cecr = as.table(rbind(c(7,1),c(0, 0)))
colnames(T1vT3.cecr) <- rownames(T1vT3.cecr) <- c("Yes", "No")
names(dimnames(T1vT3.cecr)) = c("T1", "T4")
T1vT3.cecr
margin.table(T1vT3.cecr, 1)
margin.table(T1vT3.cecr, 2)
sum(T1vT3.cecr)
mcnemar.test(T1vT3.cecr, correct=FALSE)
#McNemar's chi-squared = 1, df = 1, p-value = 0.3173
## If I didn't take the within-subjects nature of my data into account:
prop.test(rbind(margin.table(T1vT3.cecr, 1), margin.table(T1vT3.cecr, 2)), correct=FALSE)

T1vT2.cecr = as.table(rbind(c(8,1),c(0, 0)))
colnames(T1vT2.cecr) <- rownames(T1vT2.cecr) <- c("Yes", "No")
names(dimnames(T1vT2.cecr)) = c("T1", "T4")
T1vT2.cecr
margin.table(T1vT2.cecr, 1)
margin.table(T1vT2.cecr, 2)
sum(T1vT2.cecr)
mcnemar.test(T1vT2.cecr, correct=FALSE)
#McNemar's chi-squared = 1, df = 1, p-value = 0.3173
## If I didn't take the within-subjects nature of my data into account:
prop.test(rbind(margin.table(T1vT2.cecr, 1), margin.table(T1vT2.cecr, 2)), correct=FALSE)

#Alternative that allows for multiple comparisons among timepoints while 
#accounting for ind'l predator:
mydata_cecr <- data.frame(
  outcome = c(1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,0, 
              1,1,1,1, 1,0,0,0, 1,1,1,1, 1,1,1,0),
  time = gl(4,1,32,labels=1:4),
  cecr = gl(8,4,labels=letters[1:8])
)
mydata_cecr$outcome <- factor(
  mydata_cecr$outcome, levels = c(1, 0),
  labels = c("strike", "no strike")
)
# Cross-tabulation
xtabs(~outcome + time, mydata_cecr)
#           time
#outcome     1 2 3 4
#  strike    8 7 7 5
#  no strike 0 1 1 3

# Compare the proportion of strikes across timepoints
cochran_qtest(mydata_cecr, outcome ~ time|cecr)
# A tibble: 1 x 6
#  .y.         n statistic    df      p  method          
#* <chr>   <int>     <dbl> <dbl>  <dbl>  <chr>           
#1 outcome     8      6.33     3 0.0965  Cochran's Q test

# pairwise comparisons between timepoints not necessary bc Q test p>0.05
#pairwise_mcnemar_test(mydata_cecr, outcome ~ time|cecr)



## Question: 2007 DATA: Is there a difference in the proportion of trials in 
## which a LIONFISH struck at a cleaner goby in the first trial versus the 
## eighth trial? In the first versus the 2nd, 3rd, 4th, 5th, 6th, 7th?

## Answer (see below): The proportion of trials in which a lionfish struck at a 
## cleaner goby declined significantly after only one trial (p<0.001, McNemar’s 
## exact test), and continued to decline over the course of 4 trials.

T1vT8 = as.table(rbind(c(1, 7),c(0, 0)))
colnames(T1vT8) <- rownames(T1vT4) <- c("Yes", "No")
names(dimnames(T1vT8)) = c("T1", "T8")
T1vT8
margin.table(T1vT8, 1)
margin.table(T1vT8, 2)
sum(T1vT8)
mcnemar.test(T1vT8, correct=FALSE)
#McNemar's chi-squared = 7, df = 1, p-value = 0.008151

T1vT5.6.7 = as.table(rbind(c(1, 8),c(0, 0)))
colnames(T1vT5.6.7) <- rownames(T1vT5.6.7) <- c("Yes", "No")
names(dimnames(T1vT5.6.7)) = c("T1", "T5, 6, and 7")
T1vT5.6.7
margin.table(T1vT5.6.7, 1)
margin.table(T1vT5.6.7, 2)
sum(T1vT5.6.7)
mcnemar.test(T1vT5.6.7, correct=FALSE)
#McNemar's chi-squared = 8, df = 1, p-value = 0.004678

T1vT4 = as.table(rbind(c(2, 7),c(0, 0)))
colnames(T1vT4) <- rownames(T1vT4) <- c("Yes", "No")
names(dimnames(T1vT4)) = c("T1", "T4")
T1vT4
margin.table(T1vT4, 1)
margin.table(T1vT4, 2)
sum(T1vT4)
mcnemar.test(T1vT4, correct=FALSE)
#McNemar's chi-squared = 7, df = 1, p-value = 0.008151

T1vT2.3 = as.table(rbind(c(3, 6),c(0, 0)))
colnames(T1vT2.3) <- rownames(T1vT2.3) <- c("Yes", "No")
names(dimnames(T1vT2.3)) = c("T1", "T2 and 3")
T1vT2.3
margin.table(T1vT2.3, 1)
margin.table(T1vT2.3, 2)
sum(T1vT2.3)
mcnemar.test(T1vT2.3, correct=FALSE)
#McNemar's chi-squared = 6, df = 1, p-value = 0.01431

#Alternative that allows for multiple comparisons among timepoints while 
#accounting for ind'l predator:
mydata_lf2007 <- data.frame(
  outcome = c(1,0,0,0,0,0,0,0, 1,1,0,0,0,0,0,0, 1,1,0,0,0,0,0,0, 
              1,0,0,1,0,0,0,0, 1,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1, 
              1,0,1,0,0,0,0,0, 1,0,0,0,0,0,0,0),
  time = gl(8,1,64,labels=1:8),
  ptvo = gl(8,8,labels=letters[1:8])
)
mydata_lf2007$outcome <- factor(
  mydata_lf2007$outcome, levels = c(1, 0),
  labels = c("strike", "no strike")
)
# Cross-tabulation
xtabs(~outcome + time, mydata_lf2007)
#           time
#outcome     1 2 3 4
#  strike    8 7 7 5
#  no strike 0 1 1 3

# Compare the proportion of strikes across timepoints
cochran_qtest(mydata_lf2007, outcome ~ time|ptvo)
# A tibble: 1 x 6
#  .y.         n statistic    df         p  method          
#* <chr>   <int>     <dbl> <dbl>     <dbl>  <chr>           
#1 outcome     8      32.4     7 0.0000348  Cochran's Q test

# pairwise comparisons between timepoints
lf2007 <- pairwise_mcnemar_test(mydata_lf2007, outcome ~ time|ptvo, correct = F)
#overly conservative adjustments bc I do not want all paired comparisons -- so 
#focus on (unadjusted) comparisons between timepoints 1 and 2 (1 and 3 also?)
# A tibble: 28 x 6
#  group1 group2      p  p.adj p.adj.signif method      
#* <chr>  <chr>   <dbl> <dbl> <chr>        <chr>       
#1 1      2      0.0253  0.557 ns           McNemar test
#2 1      3      0.0143  0.315 ns           McNemar test



## Question: Is there a difference in the gill ventilation rate experienced by 
## LIONFISH after consuming a cleaner goby versus a bridled goby, in terms of 
## (1) number of gill opercular beats per minute, and (2) minutes until gill
## ventilation rate returns to normal (within 6 beats per minute)?

## Answer (see below): Yes, there is a difference in the number of gill
## opercular beats per minute (p<0.0001) and the time until gill ventilation
## rates returned to normal (p<0.0001) in lionfish after eating a cleaner
## goby versus a bridled goby.

gmm<-read.csv('/FILEPATH/Tuttle_et_al_FUNECOL_gmm.csv', header=T)
head(gmm)
ptvo<-subset(gmm, Pred_sp=="PTVO")
cecr<-subset(gmm, Pred_sp=="CECR")
wilcox.test(X0min_gmm~Prey_sp, data=ptvo, exact=FALSE, correct=FALSE)
#W = 0, p-value = 1.494e-06
wilcox.test(Time_to_normal~Prey_sp, data=ptvo, exact=FALSE, correct=FALSE)
#W = 0, p-value = 1.429e-06
wilcox.test(BL_gmm~Prey_sp, data=ptvo, exact=FALSE, correct=FALSE)
#W = 86.5, p-value = 0.01211
wilcox.test(BL_gmm~Prey_sp, data=cecr, exact=FALSE, correct=FALSE)
#W = 72, p-value = 0.5693
wilcox.test(BL_gmm~Prey_sp, data=gmm, exact=FALSE, correct=FALSE)
#W = 467.5, p-value = 0.9195



## Question: Is there a difference in the gill ventilation rate experienced by 
## GRAYSBY after consuming a cleaner goby versus a bridled goby, in terms of 
## (1) number of gill opercular beats per minute, and (2) minutes until gill
## ventilation rate returns to normal (within 6 beats per minute)?

## Answer (see below): Yes, there is a difference in the number of gill
## opercular beats per minute (p<0.0001) and the time until gill ventilation
## rates returned to normal (p<0.0001) in graysby after eating a cleaner
## goby versus a bridled goby.

cecr<-subset(gmm, Pred_sp=="CECR")
wilcox.test(X0min_gmm~Prey_sp, data=cecr, exact=FALSE, correct=FALSE)
#W = 11.5, p-value = 0.001116
wilcox.test(Time_to_normal~Prey_sp, data=cecr, exact=FALSE, correct=FALSE)
#W = 18, p-value = 0.004175



## Question: After eating a CLEANER goby, is there a difference between 
## lionfish and graysby in the (1) number of gill opercular beats per minute, 
## and (2) minutes until gill ventilation rate returns to normal (within 6 
## beats per minute)?

## Answer (see below): When compared to graysby, lionfish had significantly 
## higher gill ventilation rates after eating a cleaner goby (p=0.0002, rank
## sum test), and it also took significantly longer for lionfish gill 
## ventilation rates to return to normal after eating a cleaner goby than it 
## did for graysby (p=0.019, rank-sum test).

elge<-subset(gmm, Prey_sp=="ELGE")
wilcox.test(X0min_gmm~Pred_sp, data=elge, exact=FALSE, correct=FALSE)
#W = 0, p-value = 0.0002297
wilcox.test(Time_to_normal~Pred_sp, data=elge, exact=FALSE, correct=FALSE)
#W = 16.5, p-value = 0.01868



## Question: After eating a BRIDLED goby, is there a difference between 
## lionfish and graysby in the (1) number of gill opercular beats per minute, 
## and (2) minutes until gill ventilation rate returns to normal (within 6 
## beats per minute)?

## Answer (see below): When compared to graysby, lionfish had significantly 
## higher gill ventilation rates after eating a bridled goby (p<0.0001, rank
## sum test). However, there was no difference between graysby and lionfish 
## in the time it took for their gill ventilation rates to return to normal 
## after eating a bridled goby (0.887, rank-sum test).

cogl<-subset(gmm, Prey_sp=="COGL")
wilcox.test(X0min_gmm~Pred_sp, data=cogl, exact=FALSE, correct=FALSE)
#W = 41.5, p-value = 5.031e-06
wilcox.test(Time_to_normal~Pred_sp, data=cogl, exact=FALSE, correct=FALSE)
#W = 258.5, p-value = 0.8869


## FIELD EXPERIMENT
## Plotting the change in abundance of Elacatinus spp., averaged by treatment
change_ab<-read.csv('/FILEPATH/Tuttle_et_al_FUNECOL_change_ab_plot.csv',header=T)
change_ab
high <- change_ab[change_ab$tx=="high",]
low <- change_ab[change_ab$tx=="low",]
highcol <- "red"
lowcol <- "black"
plot(mean_change_ab~date_block, data=change_ab, type='n', 
     ylab='Change in abundance', xlab='Day',
     ylim=c(-16,0), xlim=c(0,80))
abline(h=0, lty="dotted")
points(mean_change_ab~date_block, data=change_ab, 
       subset=tx=="high", col = highcol, pch = 16)
points(mean_change_ab~date_block, data=change_ab, 
       subset=tx=="low", col = lowcol, pch = 15)
lines(mean_change_ab~date_block, data=change_ab, 
      subset=tx=="high", col = highcol, lty=2)
lines(mean_change_ab~date_block, data=change_ab, 
      subset=tx=="low", col = lowcol, lty=1)
segments(high$date_block, high$mean_change_ab - high$sem_change_ab,
         high$date_block, high$mean_change_ab + high$sem_change_ab, col=highcol)
segments(low$date_block, low$mean_change_ab - low$sem_change_ab,
         low$date_block, low$mean_change_ab + low$sem_change_ab, col=lowcol)
legend(2,-13,c('high-lionfish','low-lionfish'), 
       pch=c(16,15), 
       col=c(highcol, lowcol))

## Plotting the change in density of Elacatinus spp., averaged by treatment
change.de<-read.csv('/FILEPATH/Tuttle_et_al_FUNECOL_change_ab_plot.csv',header=T)
change.de
high <- change.de[change.de$tx=="high",]
low <- change.de[change.de$tx=="low",]
plot(mean_change_de~date_block, data=change.de, type='n',
     ylab='', xlab='', bty='l',
     ylim=c(-.07,.01), xlim=c(0,80))
points(mean_change_de~date_block, data=change.de, 
       subset=tx=="high", pch = 15, cex=1.5)
points(mean_change_de~date_block, data=change.de, 
       subset=tx=="low", pch = 0, cex=1.5)
lines(mean_change_de~date_block, data=change.de, 
      subset=tx=="high", lty=1, lwd=3)
lines(mean_change_de~date_block, data=change.de, 
      subset=tx=="low", lty=2, lwd=3)
segments(high$date_block, high$mean_change_de - high$sem_change_de,
         high$date_block, high$mean_change_de + high$sem_change_de)
segments(low$date_block, low$mean_change_de - low$sem_change_de,
         low$date_block, low$mean_change_de + low$sem_change_de)
abline(h=0, lty="dotted", lwd=2)
legend(0,-0.05,c('',''), bty='n',
       pch=c(15,0), cex=c(1.5,1.5), lty=c(1,2), lwd=c(3,3))



## Question: Did lionfish cause a change in the abundance of Elacatinus spp.
## (E. genie and E. evelynae; NOT E. atronasus nor E. horsti)?

## Answer (see below): There was a decline in the number of Elacatinus gobies
## on all study reefs over the course of 10 weeks. However, the presence of 
## lionfish did not affect Elacatinus spp. abundance (LRT, p=0.696)

## Statistics for change in abundance of Elacatinus spp., by treatment
cleaner<-read.csv('/FILEPATH/Tuttle_et_al_FUNECOL_cleaner_abund_subset.csv',header=T)
head(cleaner)
## create new dataframe without baseline "0" rows (for "change in" analyses) 
dim(cleaner)
noBL <- cleaner[!cleaner$date_block=="0",]
dim(noBL)

library(nlme)
library(lattice)

## look at the data
xyplot(change_ab ~ day | treatment,
       type = "p", data = noBL, groups = reef)
xyplot(change_ab ~ day | treatment, data = noBL,
       panel = function(x, y) {
         panel.xyplot(x, y)
         panel.grid(h = -1, v = -1)
         panel.abline(h = 0, lty = 2)
         panel.loess(x, y, col = 1, lwd = 2) })

plot(change_ab ~ ptvo_density_m2, data=noBL)
plot(change_ab ~ cephalopholis_density_m2, data=noBL)

## fit a linear model with gls and full fixed effects using noBL
f.full <- formula(change_ab ~ treatment * day + ptvo_density_m2 
                  + cephalopholis_density_m2)
M.gls <- gls(f.full, data = noBL)
anova(M.gls)
summary(M.gls)

## now we take a detailed look at the residuals
E.gls <- resid(M.gls, type="normalized")
F.gls <- fitted(M.gls)
op <- par(mfrow = c(3,3), mar = c(2,2,4,1))
plot(F.gls, E.gls, col = as.numeric(noBL$reef), main = "Fitted GLS")
abline(0,0)
plot(E.gls ~ noBL$day, main = "Time")
abline(0,0)
boxplot(E.gls ~ noBL$reef, main = "Reef")
abline(0,0)
boxplot(E.gls ~ noBL$treatment, main = "Treatment")
abline(0,0)
boxplot(E.gls ~ noBL$pair, main = "Pair")
abline(0,0)
plot(E.gls ~ noBL$cephalopholis_density_m2, main = "Cephalopholis")
abline(0,0)
plot(E.gls ~ noBL$ptvo_density_m2, main = "Lionfish")
abline(0,0)
par(op)
## patterns:
## 1.  normality and homoscedasticity seem like fair assumptions overall
## 2.  residuals for two reefs are higher (T08, T09)
## 3.  residuals for one reef pair is higher (T08, T09)

## We'll try adding a temporal correlation term just in case...
M.gls.c <- update(M.gls, .~.,
                     correlation = corCAR1(form = ~ day | reef))
anova(M.gls, M.gls.c) # No need to include temporal autocorrelation

## next, we'll look for the optimal model in terms of the fixed predictors:
## day, treatment (high v low lionfish density), day:tx, lionfish density, 
## and Cephalopholis spp. density
qqnorm(M.gls, abline=c(0,1)) #looks approx normal
hist(E.gls)
## we'll use stepwise selection and likelihood ratio tests to
## decide which variables should be in the model
## therefore, we have to refit using "ML"
M1.full <- update(M.gls, method = "ML",
                  control = lmeControl(maxIter = 200, msMaxIter = 200))
library(MASS)
step <- stepAIC(M1.full, direction="both")
step$anova # display results -- I want to keep treatment in the model bc 
## it is the source of the research question
M1.redA <- update(M1.full, .~. -treatment:day)
M1.redB <- update(M1.redA, .~. -cephalopholis_density_m2)
M1.reduced <- update(M1.redB, .~. -ptvo_density_m2)
summary(M1.reduced)
anova(M1.reduced, M1.full) #no evidence against the reduced model

## now we'll refit using REML
M.final <- update(M1.reduced, method = "REML")
summary(M.final)
#Coefficients:
#                  Value Std.Error    t-value p-value
#(Intercept)   2.5035555  4.590976  0.5453210  0.5897
#treatmentlow -1.1875000  3.010490 -0.3944541  0.6961
#day          -0.1998202  0.081761 -2.4439486  0.0208
intervals(M.final)
#Coefficients:
#                  lower       est.       upper
#(Intercept)  -6.8860446  2.5035555 11.89315552
#treatmentlow -7.3446427 -1.1875000  4.96964271
#day          -0.3670407 -0.1998202 -0.03259975



## Question: Is there a difference in Cephalopholis spp. density on 
## low- vs. high-lionfish treatment reefs?

## Answer (see below): There is no evidence that Cephalopholis spp.
## densities differ on low- and high-lionfish treatment reefs (p>0.1,
## rank-sum test).

## Statistics for change in abundance of Elacatinus spp., by treatment
ceph<-read.csv('/FILEPATH/Tuttle_et_al_FUNECOL_mean_cephalopholis_density.csv',header=T)
ceph
wilcox.test(mean_ceph_dens~treatment, data=ceph, exact=FALSE, correct=FALSE)
#W = 8, p-value = 1
wilcox.test(mean_ceph_juv_dens~treatment, data=ceph, exact=FALSE, correct=FALSE)
#W = 5, p-value = 0.3865
wilcox.test(mean_ceph_adult_dens~treatment, data=ceph, exact=FALSE, correct=FALSE)
#W = 10, p-value = 0.5637



## Question: Did lionfish cause a change in the DENSITY of Elacatinus spp.
## (E. genie and E. evelynae; NOT E. atronasus nor E. horsti)?

## Answer (see below): The density of cleaning gobies on experimental reefs 
## declined significantly over time (GLMM: t = -12.636, df = 22, p < 0.001) 
## at a rate of -0.00075 fish per m2 per day (95% CI -0.00062, -0.00087), 
## but lionfish had no effect on the change in density (GLMM: p = 0.696).

## Statistics for change in DENSITY of Elacatinus spp., by treatment
cleaner<-read.csv('/FILEPATH/Tuttle_et_al_FUNECOL_cleaner_dens_subset.csv',header=T)
head(cleaner)
## create new dataframe without baseline "0" rows (for "change in" analyses) 
dim(cleaner)
noBL <- cleaner[!cleaner$date_block=="0",]
dim(noBL)

library(nlme)
library(lattice)

## look at the data
xyplot(change_de ~ day | treatment,
       type = "p", data = noBL, groups = reef)
xyplot(change_de ~ day | treatment, data = noBL,
       panel = function(x, y) {
         panel.xyplot(x, y)
         panel.grid(h = -1, v = -1)
         panel.abline(h = 0, lty = 2)
         panel.loess(x, y, col = 1, lwd = 2) })
plot(change_de ~ ptvo_density_m2, data=noBL)
plot(change_de ~ cephalopholis_density_m2, data=noBL)
plot(change_de ~ ceph_juv_de, data=noBL)
plot(change_de ~ ceph_adult_de, data=noBL)

## fit a linear model with gls and full fixed effects using noBL
f.full <- formula(change_de ~ treatment * day)
M.gls <- gls(f.full, data = noBL)
anova(M.gls)
summary(M.gls)
#Generalized least squares fit by REML
#      AIC      BIC   logLik
#-73.70792 -67.0469 41.85396
#Coefficients:
#                        Value  Std.Error    t-value p-value
#(Intercept)      -0.000084880 0.02604705 -0.0032587  0.9974
#treatmentlow     -0.008105306 0.03683609 -0.2200371  0.8274
#day              -0.000590570 0.00049102 -1.2027496  0.2391
#treatmentlow:day  0.000115759 0.00069440  0.1667031  0.8688

## now we take a detailed look at the residuals
E.gls <- resid(M.gls, type="normalized")
F.gls <- fitted(M.gls)
op <- par(mfrow = c(3,3), mar = c(2,2,4,1))
plot(F.gls, E.gls, col = as.numeric(noBL$reef), main = "Fitted GLS")
abline(0,0)
plot(E.gls ~ noBL$day, main = "Time")
abline(0,0)
boxplot(E.gls ~ noBL$reef, main = "Reef")
abline(0,0)
boxplot(E.gls ~ noBL$treatment, main = "Treatment")
abline(0,0)
boxplot(E.gls ~ noBL$pair, main = "Pair")
abline(0,0)
plot(E.gls ~ noBL$cephalopholis_density_m2, main = "Cephalopholis")
abline(0,0)
plot(E.gls ~ noBL$ptvo_density_m2, main = "Lionfish")
abline(0,0)
plot(E.gls ~ noBL$ceph_juv_de, main = "Cephalopholis juv")
abline(0,0)
plot(E.gls ~ noBL$ceph_adult_de, main = "Cephalopholis adult")
abline(0,0)
par(op)
## patterns:
## 1.  heteroscedasticity with respect to reef, pair
## 2.  residuals for some reefs/pairs are high or low
## 3.  normality and heteroscedasticity look good w respect to Time and Tx
##     (although independence an issue bc of repeated measures...)
## 4.  one data point with much higher Ceph density, which might drive any 
##     correlation bw Ceph density and change in Elac density

#fit a model with a random reef effect
M.rrs <- lme(fixed = f.full, data = noBL,
             random = ~1 | reef)
a1 <- anova(M.gls, M.rrs)
a1
#      Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#M.gls     1  5 -73.70792 -67.04690 41.85396                        
#M.rrs     2  6 -97.59231 -89.59908 54.79615 1 vs 2 25.88438  <.0001
0.5 * (1 - pchisq(a1$L.Ratio[2], 1))
## definitely best to include reefs as random effects structure

## graphic methods to check for departures from assumptions, again
E.rrs <- resid(M.rrs, type = "normalized")
F.rrs <- fitted(M.rrs)
op <- par(mfrow = c(3,3), mar = c(2,2,4,1))
plot(F.rrs, E.rrs, col = as.numeric(noBL$reef), main = "Fitted LME.rrs")
abline(0,0)
plot(E.rrs ~ noBL$day, main = "Time")
abline(0,0)
boxplot(E.rrs ~ noBL$reef, main = "Reef")
abline(0,0)
boxplot(E.rrs ~ noBL$treatment, main = "Treatment")
abline(0,0)
boxplot(E.rrs ~ noBL$pair, main = "Pair")
abline(0,0)
plot(E.rrs ~ noBL$cephalopholis_density_m2, main = "Cephalopholis")
abline(0,0)
plot(E.rrs ~ noBL$ptvo_density_m2, main = "Lionfish")
abline(0,0)
plot(E.rrs ~ noBL$ceph_juv_de, main = "Cephalopholis juv")
abline(0,0)
plot(E.rrs ~ noBL$ceph_adult_de, main = "Cephalopholis adult")
abline(0,0)
par(op)
#some issues still evident:
#1.  heteroscedasticity with respect to reef
#2.  middle time points have smaller range of residuals than 1st and last

## might want to add Cephalopholis density as random effect
M.rrs2 <- lme(fixed = f.full, data = noBL,
             random = list(~1|reef, ~1|cephalopholis_density_m2))
summary(M.rrs2)
a2 <- anova(M.gls, M.rrs, M.rrs2)
a2
#       Model df       AIC       BIC   logLik   Test   L.Ratio p-value
#M.gls      1  5 -73.70792 -67.04690 41.85396                         
#M.rrs      2  6 -97.59231 -89.59908 54.79615 1 vs 2 25.884382  <.0001
#M.rrs2     3  7 -95.60209 -86.27666 54.80105 2 vs 3  0.009785  0.9212
## no evidence to suggest Ceph density necessary as random effect

## might want to add Cephalopholis JUVENILE density as random effect
M.rrs3 <- lme(fixed = f.full, data = noBL,
              random = list(~1|reef, ~1|ceph_juv_de))
summary(M.rrs3)
a3 <- anova(M.gls, M.rrs, M.rrs3)
a3
## no evidence to suggest Ceph juv density necessary as random effect

## might want to add Cephalopholis ADULT density as random effect
M.rrs4 <- lme(fixed = f.full, data = noBL,
              random = list(~1|reef, ~1|ceph_adult_de))
summary(M.rrs4)
a4 <- anova(M.gls, M.rrs, M.rrs4)
a4
## no evidence to suggest Ceph adult density necessary as random effect

#next we'll try models with different variances for different groups
M.rrs.vr <- update(M.rrs, weights = varIdent(form = ~1 | reef))
a5 <- anova(M.rrs, M.rrs.vr)
a5 # looks like variance structure helps!
#         Model df        AIC       BIC   logLik   Test  L.Ratio p-value
#M.rrs        1  6  -97.59231 -89.59908 54.79615                        
#M.rrs.vr     2 13 -101.52118 -84.20252 63.76059 1 vs 2 17.92887  0.0123

#we'll check the residual patterns from our best model so far...
E.rrs.vr <- resid(M.rrs.vr, type = "normalized")
F.rrs.vr <- fitted(M.rrs.vr)
op <- par(mfrow = c(3,3), mar = c(2,2,4,1))
plot(F.rrs.vr, E.rrs.vr, col = as.numeric(noBL$reef), main = "Fitted LME.rrs")
abline(0,0)
plot(E.rrs.vr ~ noBL$day, main = "Time")
abline(0,0)
boxplot(E.rrs.vr ~ noBL$reef, main = "Reef")
abline(0,0)
boxplot(E.rrs.vr ~ noBL$treatment, main = "Treatment")
abline(0,0)
boxplot(E.rrs.vr ~ noBL$pair, main = "Pair")
abline(0,0)
plot(E.rrs.vr ~ noBL$cephalopholis_density_m2, main = "Cephalopholis")
abline(0,0)
plot(E.rrs.vr ~ noBL$ptvo_density_m2, main = "Lionfish")
abline(0,0)
plot(E.rrs.vr ~ noBL$ceph_juv_de, main = "Cephalopholis juv")
abline(0,0)
plot(E.rrs.vr ~ noBL$ceph_adult_de, main = "Cephalopholis adult")
abline(0,0)
par(op)
#looks WAY better! heteroscedasticity and normality no longer an issue
#with respect to reef, pair

## We'll try adding a temporal correlation term just in case...
M.rrs.vr.c <- update(M.rrs.vr, .~.,
                  correlation = corCAR1(form = ~ day | reef))
anova(M.gls, M.rrs, M.rrs2, M.rrs.vr, M.rrs.vr.c) # Model W/O autocorr factor
## performs marginally better (i.e., no evidence against model 4)
#           Model df        AIC       BIC   logLik   Test   L.Ratio p-value
#M.gls          1  5  -73.70792 -67.04690 41.85396                        
#M.rrs          2  6  -97.59231 -89.59908 54.79615 1 vs 2  25.88438  <.0001
#M.rrs2         3  7  -95.60209 -86.27666 54.80105 2 vs 3  0.009785  0.9212
#M.rrs.vr       4 13 -101.52118 -84.20252 63.76059 3 vs 4 17.919084  0.0064
#M.rrs.vr.c     5 14  -99.52118 -80.87031 63.76059 4 vs 5  0.000000  1.0000
## next, we'll look for the optimal model in terms of the fixed predictors:
## day, treatment (high v low lionfish density), day:tx
qqnorm(M.rrs.vr, abline=c(0,1))
hist(E.rrs.vr)

## we'll use likelihood ratio tests to decide which variables should be 
## in the model. therefore, we have to refit using "ML"
M1.full <- update(M.rrs.vr, method = "ML",
                  control = lmeControl(maxIter = 200, msMaxIter = 200))
anova(M1.full) #day and tx:day interaction are significant, but not tx
#              numDF denDF   F-value p-value
#(Intercept)       1    22   6.58132  0.0176
#treatment         1     6   0.00915  0.9269
#day               1    22 162.85886  <.0001
#treatment:day     1    22  43.36880  <.0001

M1.A <- update(M1.full, .~. -treatment:day)
anova(M1.A, M1.full) #keep interaction
#        Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#M1.A        1 12 -140.8431 -123.2543 82.42157                        
#M1.full     2 13 -148.9995 -129.9450 87.49978 1 vs 2 10.15642  0.0014

M2.A <- update(M1.A, .~. -treatment)
M2.B <- update(M1.A, .~. -day)
anova(M1.A, M2.A)
anova(M1.A, M2.B)

#now we'll refit using REML
M.final <- update(M1.full, method = "REML")
summary(M.final)
#      AIC       BIC   logLik
#-101.5212 -84.20252 63.76059
#Fixed effects: list(f.full) 
#                       Value   Std.Error DF    t-value p-value
#(Intercept)       0.00800772 0.018336046 22   0.436720  0.6666
#treatmentlow     -0.03337886 0.026120200  6  -1.277894  0.2485
#day              -0.00075150 0.000059472 22 -12.636216  0.0000
#treatmentlow:day  0.00062460 0.000105017 22   5.947623  0.0000
intervals(M.final)
#Approximate 95% confidence intervals
#Fixed effects:
#                         lower          est.         upper
#(Intercept)      -0.0300189120  0.0080077194  0.0460343508
#treatmentlow     -0.0972926885 -0.0333788603  0.0305349678
#day              -0.0008748427 -0.0007515047 -0.0006281667
#treatmentlow:day  0.0004068096  0.0006246016  0.0008423935
