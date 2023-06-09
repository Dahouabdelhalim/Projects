# D. pavonius experiment

########## load packages ##########
library(readr)    # read data
library(lme4)     # random effects models
library(ggplot2)  # plots

########## load data ##########
dp_data <- read_delim("dp.data.txt", "\\t",    escape_double = FALSE, trim_ws = TRUE)
#test_dp <- read_delim("test.dp.txt", "\\t",  escape_double = FALSE, trim_ws = TRUE)
dp.data<- data.frame(dp_data)
head(dp.data)

# each row in the data frame is a single combinaiton of family*cohort*pathogen_chamber

######## Clean up data: Remove combinations with <10 individuals ########

 sub.dp.data<-subset(dp.data, tab.N>9)
 length(dp.data$tab.D)          # original number of rows
 length(sub.dp.data$tab.D)      # number of rows in subset
 head(sub.dp.data);

###########  Test effect of dew chamber ############
 
 library(lme4)
 
# start with a full model, that fits both path and chamber as random effects, with chamber nested within path
  mod1<- glmer(cbind(tab.D, tab.N-tab.D)~1+(1|Path/Chamber.rep), binomial, data=sub.dp.data)   
  summary(mod1)
  
  mod1.2<- glmer(cbind(tab.D, tab.N-tab.D)~1+(1|Path), binomial, data=sub.dp.data)   
  summary(mod1.2)
  anova(mod1,mod1.2) # test effect of omitting chamber
  
  mod1.3<- glmer(cbind(tab.D, tab.N-tab.D)~1+(1|Chamber.rep), binomial, data=sub.dp.data)   
  anova(mod1,mod1.3) # test effect of omitting pathogen
  
  ##### associated figure ()
  library(ggplot2)
  p1 <- ggplot(sub.dp.data, aes(x=Path, fill=as.factor(Chamber.rep), y=tab.D/tab.N)) + 
        xlab("Pathogen genotype")+
        ylab("Proportion infected")+
        ylim(0,1)+
        geom_boxplot(outlier.shape=NA)+
        theme_classic()+
        geom_jitter(shape=16, position=position_jitterdodge(.2), aes(shape=Path), col='black',alpha=0.7)
  p1  
  
##########################################################################################  
#####################  Main analysis: Glm models testing effects of path, age, and family   #####################
##########################################################################################  
  
  ## make wowrking data set:  each row is a single family*age*path
  
  big.combo<-paste(dp.data$Cohort, dp.data$Fam, dp.data$Path)
  newdata<-data.frame(dp.data, big.combo)
  tab.N<-tapply(newdata$tab.N,newdata$big.combo,sum)
  tab.D<-as.vector(tapply(newdata$tab.D,newdata$big.combo,sum))
  data.newdata<-data.frame(cbind(tab.N, tab.D))
  ID<-c(row.names(data.newdata))
  data.newdata<-data.frame(cbind(ID,data.newdata[,]))
  library(dplyr)
  library(tidyr)
  new.cols<-data.newdata %>% separate(ID, c("Cohort", "Fam","Path"))  
  head(new.cols)
  cohort.name<-ifelse(new.cols$Cohort=="A1B1", "11-month",ifelse(new.cols$Cohort=="A1B2","10-month",
                                                                 ifelse(new.cols$Cohort=="A1B3","9-month","1-month" )))
  days.old<-ifelse(new.cols$Cohort=="A1B1", 322,ifelse(new.cols$Cohort=="A1B2",260,
                                                       ifelse(new.cols$Cohort=="A1B3",247,30)))
  new.cols<-data.frame(new.cols, cohort.name, days.old)
  sub.new.cols<-subset(new.cols, tab.N>9)
  head(sub.new.cols); length(sub.new.cols$tab.D)
############################################################################
  ##################### Test effect of cohorts ######################
  
  ############ main effects of cohort ############
  
      cohort.mod<-  glm(cbind(tab.D, tab.N-tab.D)~Cohort,  quasibinomial, data=sub.new.cols)
      summary(cohort.mod)
  
  ############### Which cohorts differ?  Can any be combined? #############################
  #### Where is the difference in cohorts?  Test pairwise, and use Bonferroni to compare ###
   
  crit.q=0.05/6   #  0.008333333
  
  subB1.B2<-subset(sub.new.cols, Cohort=="A1B1" | Cohort=="A1B2")
  subB1.B3<-subset(sub.new.cols, Cohort=="A1B1" | Cohort=="A1B3")
  subB1.A3<-subset(sub.new.cols, Cohort=="A1B1" | Cohort=="A3")
  subB2.B3<-subset(sub.new.cols, Cohort=="A1B2" | Cohort=="A1B3")
  subB2.A3<-subset(sub.new.cols, Cohort=="A1B2" | Cohort=="A3")
  subB3.A3<-subset(sub.new.cols, Cohort=="A1B3" | Cohort=="A3")
  
  test1<-  glm(cbind(tab.D, tab.N-tab.D)~Cohort, quasibinomial, data=subB1.B2); anova(test1, test="Chi")
  test2<-  glm(cbind(tab.D, tab.N-tab.D)~Cohort, quasibinomial, data=subB1.B3); anova(test2, test="Chi")
  test3<-  glm(cbind(tab.D, tab.N-tab.D)~Cohort, quasibinomial, data=subB1.A3); anova(test3, test="Chi")
  test4<-  glm(cbind(tab.D, tab.N-tab.D)~Cohort, quasibinomial, data=subB2.B3); anova(test4, test="Chi")  # n.s.
  test5<-  glm(cbind(tab.D, tab.N-tab.D)~Cohort, quasibinomial, data=subB2.A3); anova(test5, test="Chi")
  test6<-  glm(cbind(tab.D, tab.N-tab.D)~Cohort, quasibinomial, data=subB3.A3); anova(test6, test="Chi")  # n.s. (Does not pass crit q)

  ## make figure of cohorts:   
  FigS2.2<- ggplot(sub.new.cols, aes(x=days.old, y=tab.D/tab.N, color=Cohort))+
            geom_jitter(shape=16, position = position_jitter(width = 15), aes(size=tab.N, alpha=.5))+
            xlab("Days old")+
            ylab("Proportion infected")+
            theme(text = element_text(size = 14))  
  FigS2.2

  ### supports combining A1B2 (8.5 months) and A1B3 (8 months) together ###

#############################################################################
############  combine A1B2 and A1B3 into a single age group #################
#############################################################################
  unique(dp.data$Cohort)
  new.cohorts<-ifelse(dp.data$Cohort=="A1B1", 'A1', ifelse(dp.data$Cohort=="A1B2", "A2",
                                                             ifelse(dp.data$Cohort=="A1B3", "A2", "A3")))
  dp.data<-data.frame(dp.data, new.cohorts)

  
  #### make combos of new cohort, path and family 
  big.combo<-paste(dp.data$new.cohorts, dp.data$Fam, dp.data$Path)
  newdata3<-data.frame(dp.data, big.combo)
  tab.N<-tapply(newdata3$tab.N,newdata3$big.combo,sum)
  tab.D<-as.vector(tapply(newdata3$tab.D,newdata3$big.combo,sum))
  data.newdata3<-data.frame(cbind(tab.N, tab.D))
  ID<-c(row.names(data.newdata3))
  data.newdata3<-data.frame(cbind(ID,data.newdata3[,]))
  library(dplyr)
  library(tidyr)
  new.cols3<-data.newdata3 %>% separate(ID, c("Cohort",  "Fam","Path"))  
  head(new.cols3)
  unique(new.cols3$Cohort)
  sub.new.cols3<-subset(new.cols3, tab.N>9)
  head(sub.new.cols3); length(sub.new.cols3$tab.D)
  
  ### Final data set:  each row is one combination of family, age (a1, a2, a3), and path ###
  
################################################################################################################      
  ######## overall glm with cohort as a categorical predictor with 4 categorial predictors ############
#########################################################################################################
 
  #### this uses sub.new.cols3 where each row is one combination of family,path and age (3 gropus)  ###
  
  range(sub.new.cols3$tab.N) # check that each row has 10 or more
  
###### series of models fit with glm.  AIC scores summarized in Table S2. ###
  
  #additive (model A in Fig S3)
  mod2.A<- glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam, 
               binomial, data=sub.new.cols3)
  summary(mod2.A)
  
  ### test pairwise interactions 
  # model B, Table S3
  mod2.B<- glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam
                                        +Path:Cohort,
                                        binomial, data=sub.new.cols3)
  summary(mod2.B); 
  anova(mod2.B, test="Chi")
  anova(mod2.A, mod2.B, test="Chi")
  
  # model C, Table S3
  mod2.C<- glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+ 
                                         Path:Fam,
                                         binomial, data=sub.new.cols3)
  summary(mod2.C); 
  anova(mod2.C, test="Chi")
  anova(mod2.A, mod2.C, test="Chi")
  
  # model D, Fig S3
  mod2.D<- glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+ 
                                         Cohort:Fam,
                                         binomial, data=sub.new.cols3)
  summary(mod2.D); 
  anova(mod2.D, test="Chi")
  anova(mod2.A, mod2.D, test="Chi")
  
  ######## test 2 pairwise interactions ########
  
  mod2.E<- glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+ 
                 Path:Cohort+ Path:Fam,
               binomial, data=sub.new.cols3)
  summary(mod2.E); 
  anova(mod2.E, test="Chi")
  anova(mod2.E, mod2.B, test="Chi")
  anova(mod2.E, mod2.C, test="Chi")
  
  ## model F in Table S3 ## 
  mod2.F<- glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+ 
                                         Path:Cohort+ Cohort:Fam,
                                         binomial, data=sub.new.cols3)
  summary(mod2.F); 
  anova(mod2.F, test="Chi")  # ** best fit -table 1 in ms
  anova(mod2.B, mod2.F, test="Chi")
  anova(mod2.D, mod2.F, test="Chi")  ## test cohort*path
  
  ### same one but with quasibinomial
  mod2.F.q<- glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+ 
                                           Path:Cohort+ Cohort:Fam,
                                          quasibinomial, data=sub.new.cols3)
  summary(mod2.F.q)
  anova(mod2.F.q, test='Chi')
 
  ## Model G in Table S3
  mod2.G<- glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+ 
                                         Path:Fam+ Cohort:Fam,
                                         binomial, data=sub.new.cols3)
  summary(mod2.G); 
  anova(mod2.G, test="Chi")
  anova(mod2.C, mod2.G, test="Chi")
  anova(mod2.D, mod2.G, test="Chi")
  
  
  ## Model H (Not shown; Fully saturated, no power to test interaciton effect
  mod2.G<- glm(cbind(tab.D, tab.N-tab.D)~Path*Cohort*Fam, 
                                         binomial, data=sub.new.cols3)
  summary(mod2.G) # does not converge.
  anova(mod2.G, test="Chi")
  y.pred = predict(mod2.G, sub.new.cols3, type="response")
  hist(y.pred)
  
#### Try with subsets of the data ####
  
# compare 1 month and 8 months
  sub.2.3<-subset(sub.new.cols3, Cohort!="A1"); length(sub.2.3$Cohort)
  
  mod2.G.sub.2.3<- glm(cbind(tab.D, tab.N-tab.D)~Path*Cohort*Fam, binomial, data=sub.2.3)
               summary(mod2.G.sub.2.3)
               anova(mod2.G.sub.2.3, test="Chi")
               
# compare 1 month and 10 months
    sub.1.3<-subset(sub.new.cols3, Cohort!="A2"); length(sub.1.3$Cohort)
               
    mod2.G.sub.1.3<- glm(cbind(tab.D, tab.N-tab.D)~Path*Cohort*Fam, binomial, data=sub.1.3)
               summary(mod2.G.sub.1.3)
               anova(mod2.G.sub.1.3, test="Chi")
               
# compare 8 months and 10 months
    sub.1.2<-subset(sub.new.cols3, Cohort!="A3"); length(sub.1.2$Cohort)
               
    mod2.G.sub.1.2<- glm(cbind(tab.D, tab.N-tab.D)~Path*Cohort*Fam, binomial, data=sub.1.2)
               summary(mod2.G.sub.1.2)
               anova(mod2.G.sub.1.2, test="Chi")
               
               
############################################################################
  ## figure showing main effects and interaciton of age and pathogen
############################################################################  
  
  # rename age groups with useful labels (number of months)
   new.age.cat<-ifelse(sub.new.cols3$Cohort=="A1",10, ifelse(sub.new.cols3$Cohort=="A2",8,1))
   sub.new.cols3<-data.frame(sub.new.cols3, new.age.cat)
   
   Fig3 <- ggplot(sub.new.cols3, aes(x=as.factor(new.age.cat), y=tab.D/tab.N, fill=Path))+ 
           xlab("Host age (months)")+
           ylab("Proportion infected")+
           geom_boxplot(outlier.shape=NA)+
           theme_classic()+
           theme(text = element_text(size = 14))+
           geom_jitter(position=position_jitterdodge(.1))#+
           scale_x_discrete(labels = c('Juv.','Adult'))
   Fig3
  
############################################################################### 
###############################################################################
   #### Analyze correlations between resistance at different ages ###
   ## this code will:
      # 1) generate residuals of model that accounts for pathogen affects
      # 1) generate a wide-form data set for each age comparions with these residuals
      # 2) test the correlation 
      # 3) Make the scatterplot figures for Fig 4A-C
      # 4) Run the chi-squared test of indepdendnce to test sig of age*family interactions
      # 5) Reformat the data to long-form to make the interaciton plots (Fig 4D-F)
   
#2a: run glm with only path and path*fam (but not fam)
   
   mod.p<-glm(cbind(tab.D, tab.N-tab.D)~Path+Path:Cohort,  binomial, data=sub.new.cols3 )
   summary(mod.p)
   anova(mod.p, test="Chi")
   #plot(mod)  # a little skewed on the qq plot
   
   # (an aside, check significane with quasibinomial)
   mod.p2<-glm(cbind(tab.D, tab.N-tab.D)~Path+Path:Cohort,  quasibinomial, data=sub.new.cols3 )
   summary(mod.p2)
   anova(mod.p2, test="Chi")
   
   ##### 2b: now save the residuals  
   res<-mod.p$residuals  # save resids!

   ##### 2Cneed to add these residuals to the data set 
   
   newdata<-cbind(sub.new.cols3, res)
   ######################################################################33
   ######################subset and make plots!!! ########################
   head(newdata) # subset(new.cols2.3, Path=='LP1')
   
   path.fam<-paste(newdata$Path, newdata$Fam )  # first make path.fam name
   newdata<-cbind(path.fam,newdata)
   head(newdata)
   
   #####################################
   ####### A3 vs A2 ##########
   #####################################
   sub.adult<-subset(newdata, Cohort=="A2", select=c(1,5,6,8))   # set adult age
   colnames(sub.adult)<-c("Path.Fam", "Num.A","Dis.A","res.A")
   head(sub.adult)
   
   sub.juv<-subset(newdata, Cohort=="A3",select=c(1,5,6,8))
   colnames(sub.juv)<-c("Path.Fam", "Num.J","Dis.J","res.J"); 
   head(sub.juv)
   
   wide.3.2<-merge(sub.juv, sub.adult, by.x='Path.Fam', sort = T)  # merge by path.fam instead of fam
   head(wide.3.2)
   colnames(wide.3.2)<-c('Path-Fam', 'jN', 'jD','res.J','aN','aD','res.A')
   mins<-apply(wide.3.2[,c(2,5)], 1, FUN=min)   # update mins column
   wide.3.2<-data.frame(wide.3.2,mins)
   head(wide.3.2)
   length(wide.3.2$jN)
   
   range(wide.3.2$mins)  # data set of A3 vs A2 adjusted for path.
   
   res2<- cor.test(wide.3.2$res.J, wide.3.2$res.A,  method = "pearson")
   res2
   
   modA<-lm(wide.3.2$res.A~wide.3.2$res.J)
   summary(modA)
   
   modA.w<-lm(wide.3.2$res.A~wide.3.2$res.J, weights=wide.3.2$mins)
   summary(modA.w)
   
   library(dplyr)
   library(tidyr)
   wide.3.2.path<-wide.3.2 %>% separate(Path.Fam, c("Path","Fam"))  
   head(wide.3.2.path)
   
   library(ggplot2)
   
corr.plotA<-  ggplot(wide.3.2.path, aes(x=res.J, y=res.A))+
              geom_point(shape=1, stroke=1, aes(size=mins, color=Path))+
              theme_classic()+
              theme(text = element_text(size = 14))+
              labs(title="", x="Resid(1-month prop. infect.)", y = " Resid (8-month prop. infect.)")+
              geom_abline(intercept= modA.w$coefficients[1] , slope = modA$coefficients[2] , col='black', lty=3)+
              theme(legend.position = "none") 

corr.plotA  
###### turn data to long form for interaction plot

### now make long for correlation plot

wide<- wide.3.2

age<-rep('Juv',length(wide$jN))
age2<-rep('Adult',length(wide$jN))

long.all<-data.frame(wide[,1], age, wide[,c(4,8)]) 
colnames(long.all)<-c("Path.Fam",'age','P','mins')

long.all.1<-data.frame(wide[,1], age2, wide[c(7,8)])
colnames(long.all.1)<-c("Path.Fam",'age','P','mins')
long.all<-rbind(long.all, long.all.1)
head(long.all)    
long.3.2<- long.all


long.3.2.path<-long.3.2 %>% separate(Path.Fam, c("Path","Fam"))  
long.3.2.path<-data.frame(long.3.2.path, long.all$Path.Fam)
head(long.3.2.path)

## Interaction figure 
sum.plotA <-  ggplot(long.3.2.path, aes(x=age, y=P,group=as.factor(long.all.Path.Fam)) )+ 
              geom_line( aes(color=Path)) + #alpha=0.5,
              geom_point(aes(color=Path))+
              labs(title="", x="Host Age", y = "Resid (Prop. infect.)")+
              theme_classic()+
              theme(legend.position = 'none')+
              theme(text = element_text(size = 14))+
              scale_x_discrete(labels = c('1 month','8 months'))

sum.plotA

## stats test for interaction effects
## subset to just two age comparisons, test with full data set

  af.data<-subset(sub.new.cols3, Cohort!="A1") 

  af.test<-glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+Path:Cohort +Cohort:Fam, binomial, data=af.data)
          summary(af.test)
          anova(af.test, test="Chi")

  af.test.quasi<-glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+Path:Cohort +Cohort:Fam, quasibinomial, data=af.data)
          anova(af.test.quasi, test="Chi")
          #summary(af.test.quasi) # dispersion 1.281876
    af.test.quasi.2<-glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+Path:Cohort , quasibinomial, data=af.data)
    anova(af.test.quasi, af.test.quasi.2, test="Chi")

  #####################################
    ####### 1 month  vs 10 months ##########
  #####################################
    sub.adult<-subset(newdata, Cohort=="A1", select=c(1,5,6,8))   # set adult age
    colnames(sub.adult)<-c("Path.Fam", "Num.A","Dis.A","res.A")
    head(sub.adult)
    
    sub.juv<-subset(newdata, Cohort=="A3",select=c(1,5,6,8))
    colnames(sub.juv)<-c("Path.Fam", "Num.J","Dis.J","res.J"); 
    head(sub.juv)
    
    wide.3.1<-merge(sub.juv, sub.adult, by.x='Path.Fam', sort = T)  # merge by path.fam instead of fam
    head(wide.3.1)
    colnames(wide.3.1)<-c('Path-Fam', 'jN', 'jD','res.J','aN','aD','res.A')
    mins<-apply(wide.3.1[,c(2,5)], 1, FUN=min)   # update mins column
    wide.3.1<-data.frame(wide.3.1,mins)
    head(wide.3.1)
    length(wide.3.1$jN)
    
    range(wide.3.1$mins)  # data set of A3 vs A2 adjusted for path.
    
    res2<- cor.test(wide.3.1$res.J, wide.3.1$res.A,  method = "pearson")
    res2
    
    modA<-lm(wide.3.1$res.A~wide.3.1$res.J)
    summary(modA)
    
    modA.w<-lm(wide.3.1$res.A~wide.3.1$res.J, weights=wide.3.1$mins)
    summary(modA.w)
    
    library(dplyr)
    library(tidyr)
    wide.3.1.path<-wide.3.1 %>% separate(Path.Fam, c("Path","Fam"))  
    head(wide.3.1.path)
    
corr.plotB<-  ggplot(wide.3.1.path, aes(x=res.J, y=res.A))+
              geom_point( shape=1, stroke=1,aes(size=mins, col=Path))+
              theme_classic()+
              theme(text = element_text(size = 14))+
              labs(title="", x="Resid(1-month prop. infect.)", y = " Resid(10-month prop. infect.)")+
              geom_abline(intercept= modA$coefficients[1] , slope = modA$coefficients[2] , col='black', lty=3)+
              theme(legend.position = "none") 
corr.plotB    
    ### now make long for correlation plot
    
    wide<- wide.3.1
    
    age<-rep('Juv',length(wide$jN))
    age2<-rep('Adult',length(wide$jN))
    
    long.all<-data.frame(wide[,1], age, wide[,c(4,8)]) 
    colnames(long.all)<-c("Path.Fam",'age','P','mins')
    
    long.all.1<-data.frame(wide[,1], age2, wide[c(7,8)])
    colnames(long.all.1)<-c("Path.Fam",'age','P','mins')
    long.all<-rbind(long.all, long.all.1)
    head(long.all)    
    long.3.1<- long.all
    
    
    long.3.1.path<-long.3.1 %>% separate(Path.Fam, c("Path","Fam"))  
    long.3.1.path<-data.frame(long.3.1.path, long.all$Path.Fam)
    head(long.3.1.path)
    
## make interaction plot
sum.plotB <-  ggplot(long.3.1.path, aes(x=age, y=P, group=as.factor(long.all.Path.Fam)) )+ 
              geom_line(  aes(color=Path)) +
              geom_point(aes(color=Path))+
              labs(title="", x="Host Age", y = "Resid (Prop. infect.)")+
              theme_classic()+
              theme(legend.position = 'none')+
              theme(text = element_text(size = 14))+
              scale_x_discrete(labels = c('1 month','10 months'))
sum.plotB   
    ## stats test for interaction effects
    ## subset to just two age comparisons, test with full data set
    
    af.data<-subset(sub.new.cols3, Cohort!="A2")
    
    af.test<-glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+Path:Cohort +Cohort:Fam, binomial, data=af.data)
    summary(af.test)
    anova(af.test, test="Chi")
    
    af.test.2<-glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+Path:Cohort, binomial, data=af.data)
    summary(af.test)
    anova(af.test, af.test.2,test="Chi")
    
#####################################
## 8 month vs 10 months
#####################################
    sub.adult<-subset(newdata, Cohort=="A1", select=c(1,5,6,8))   # set adult age
    colnames(sub.adult)<-c("Path.Fam", "Num.A","Dis.A","res.A")
    head(sub.adult)
    
    sub.juv<-subset(newdata, Cohort=="A2",select=c(1,5,6,8))
    colnames(sub.juv)<-c("Path.Fam", "Num.J","Dis.J","res.J"); 
    head(sub.juv)
    
    wide.2.1<-merge(sub.juv, sub.adult, by.x='Path.Fam', sort = T)  # merge by path.fam instead of fam
    head(wide.2.1)
    colnames(wide.2.1)<-c('Path-Fam', 'jN', 'jD','res.J','aN','aD','res.A')
    mins<-apply(wide.2.1[,c(2,5)], 1, FUN=min)   # update mins column
    wide.2.1<-data.frame(wide.2.1,mins)
    head(wide.2.1)
    length(wide.2.1$jN)
    
    range(wide.2.1$mins)  # data set of A3 vs A2 adjusted for path.
    
    res2<- cor.test(wide.2.1$res.J, wide.2.1$res.A,  method = "pearson")
    res2
    
    modA<-lm(wide.2.1$res.A~wide.2.1$res.J)
    summary(modA)
    
    modA.w<-lm(wide.2.1$res.A~wide.2.1$res.J, weights=wide.2.1$mins)
    summary(modA.w)
    
    library(dplyr)
    library(tidyr)
    wide.2.1.path<-wide.2.1 %>% separate(Path.Fam, c("Path","Fam"))  
    head(wide.2.1.path)
    
corr.plotC<-  ggplot(wide.2.1.path, aes(x=res.J, y=res.A))+
              geom_point(shape=1, stroke=1,aes(size=mins, color=Path))+
              theme_classic()+
              theme(text = element_text(size = 14))+
              labs(title="", x="Resid(8-months prop. infect.)", y = " Resid(10-month prop. infect.)")+
              geom_abline(intercept= modA$coefficients[1] , slope = modA$coefficients[2] , col='black', lty=3)+
              theme(legend.position = "none") 
corr.plotC    
    ### now make long for correlation plot
    
    wide<- wide.2.1
    
    age<-rep('Juv',length(wide$jN))
    age2<-rep('Adult',length(wide$jN))
    
    long.all<-data.frame(wide[,1], age, wide[,c(4,8)]) 
    colnames(long.all)<-c("Path.Fam",'age','P','mins')
    
    long.all.1<-data.frame(wide[,1], age2, wide[c(7,8)])
    colnames(long.all.1)<-c("Path.Fam",'age','P','mins')
    long.all<-rbind(long.all, long.all.1)
    head(long.all)    
    long.2.3<- long.all
    
    long.2.3.path<-long.2.3 %>% separate(Path.Fam, c("Path","Fam"))  
    long.2.3.path<-data.frame(long.2.3.path, long.all$Path.Fam)
    head(long.2.3.path)
## Make interaction plot

sum.plotC <-  ggplot(long.2.3.path, aes(x=age, y=P, group=as.factor(long.all.Path.Fam)) )+ 
              geom_line(  aes(color=Path)) +
              geom_point(aes(color=Path))+
              labs(title="", x="Host Age", y = "Resid (Prop. infect.)")+
              theme_classic()+
              theme(legend.position = 'none')+
              theme(text = element_text(size = 14))+
              scale_x_discrete(labels = c('8 months','10 months'))
sum.plotC    
    
    ## stats test for interaction effects
    # # subset to just two age comparisons, test with full data set
    
    af.data<-subset(sub.new.cols3, Cohort!="A3")
    
    af.test<-glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+Path:Cohort +Cohort:Fam, binomial, data=af.data)
    summary(af.test)
    anova(af.test, test="Chi")
    
      test.cohort<-glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+Path:Cohort, binomial, data=af.data)
      anova(af.test, test.cohort, test="Chi")
      
      test.path<-glm(cbind(tab.D, tab.N-tab.D)~Path+Cohort+Fam+Cohort:Fam, binomial, data=af.data)
      anova(af.test, test.path, test="Chi")
    
    library("ggpubr")    
    figure4 <- ggarrange(corr.plotA, corr.plotB, corr.plotC,
                         sum.plotA, sum.plotB, sum.plotC,
                         labels = c("A", "B","C", "D","E","F"))#, nrow=3)
    
    figure4
    
    
    