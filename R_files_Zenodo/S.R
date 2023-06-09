################################################3
##  Silene vulgaris 

## save all pdfs as 9 x 6

## Load pakages
library(readr)    # reads data in
library(ggplot2)  # for graphics
library(ggpubr)   # for arranging graphics into multiframe panel figures

# load data
library(readr)
mh_sv <- read_delim("mh-sv.txt", "\\t", escape_double = FALSE, trim_ws = TRUE)
sv.data<-as.data.frame(mh_sv)
prop<-sv.data$Dis/sv.data$N   # calculate proportion diseased
sv.data<-data.frame(sv.data,prop)
head(sv.data)

## Calculate overall infection rate for both ages
tab.d<-tapply(sv.data$Dis, sv.data$age, sum); tab.d
tab.n<-tapply(sv.data$N, sv.data$age, sum); tab.n
tab.d/tab.n

############## Logistic regression

# Saturated model with interaction
modB<- glm( cbind(Dis, N-Dis)~ as.factor(age)*Family, binomial, data=sv.data )
summary(modB)
anova(modB, test="Chi")  # but no DF for testing interaction effect

## Addititive model
  modB.2<- glm( cbind(Dis, N-Dis)~ as.factor(age)+Family, binomial, data=sv.data )
  summary(modB.2)
  anova(modB.2, test="Chi")

# Additive model with Family fit first
  modB.2<- glm( cbind(Dis, N-Dis)~ Family+as.factor(age), binomial, data=sv.data )
  summary(modB.2)
  anova(modB.2, test="Chi")

## Slight overdispersion. Test with a quasibinomial model 
  modB.3<- glm( cbind(Dis, N-Dis)~ as.factor(age)+Family, quasibinomial, data=sv.data )
  summary(modB.3)
  anova(modB.3, test="Chi") # age and family are still significant
  

####### make boxplot Fig 2a ##############

  #calcu sizing
  max(sv.data$N)
  size.scale3<-sv.data$N/100

boxplotB <- ggplot(sv.data, aes(x=as.factor(age), y=Dis/N))+ 
            xlab("Host Age")+
            ylab("Proportion infected")+
            geom_boxplot()+#(outlier.shape=NA)+
            theme_classic()+
            ylim(0,1)+
            theme(text = element_text(size = 14))+
            geom_jitter(width=.2, color='blue',shape=20, alpha=0.5, aes(size=size.scale3))+#alpha=0.4,
            theme(legend.position="none")+         
            scale_x_discrete(labels = c('Juvenile','Adult'))
boxplotB 

############ Correlation tests and Figure #########

## step 1, put the data in wide form and cut families <10 
sv2.s<-subset(sv.data, age==1, select=c(1,3,4))
sv2.a<-subset(sv.data, age==2, select=c(1,3,4))
sv.wide<-cbind(sv2.s, sv2.a[,c(2:3)])
names(sv.wide)<-c("Family","N.s","Dis.s","N.a","Dis.a")
sv.wide
prop.s<-sv.wide$Dis.s/sv.wide$N.s
prop.a<-sv.wide$Dis.a/sv.wide$N.a
sv.wide<-data.frame(sv.wide, prop.s, prop.a)
sv.wide
mins<-apply(sv.wide[,c(2,4)], 1, FUN=min)
sv.wide<-data.frame(sv.wide,mins)
head(sv.wide)

##### Correlation test ###
res2<- cor.test(sv.wide$prop.s, sv.wide$prop.a,  method = "pearson")
res2


### Make correlation plot Fig 2b ####
max(sv.wide$mins) #55

modB<-lm(prop.a~prop.s, weights=mins, data=sv.wide)
summary(modB)

size.scale4<-sv.wide$mins/100

corr.plotB<-  ggplot(sv.wide, aes(x=prop.s, y=prop.a))+
              geom_point(col='blue', shape=1, stroke=1, aes(size=size.scale4))+
              theme_classic()+
              theme(text = element_text(size = 14))+
              labs(title="", x="Proportion infected (Juvenile)", y = "Proportion infected (Adult)")+
              theme(legend.position="none")+    
              geom_abline(intercept= modB$coefficients[1], slope = modB$coefficients[2] , col='black', lty=1)

corr.plotB

########## test for contingency:  Can age * juv family predict adult family? ###########
  mean.j<-mean(sv.wide$prop.s); mean.j
  mean.a<-mean(sv.wide$prop.a); mean.a
  age.factor <- mean.a/mean.j; age.factor
  pred.aD    <- sv.wide$prop.s * age.factor *sv.wide$N.a

# calculate chi square (omit cases where juveniles are 0)
  Calc.chi<- ((pred.aD - sv.wide$Dis.a)^2)/pred.aD; Calc.chi
  Calc.chi[is.infinite(Calc.chi)] <- NA
  Calc.chi<-as.vector(na.omit(Calc.chi))
  Calc.chi

  cbind(sv.wide[,c(1,4,5)], pred.aD, Calc.chi) 

  # calc df, chi-sq, p
  df.num<- length(Calc.chi); df.num 
  chi.num<- sum(Calc.chi)  ; chi.num
  pchisq(chi.num, df=df.num, lower.tail=FALSE)

############### make interaction plot

## First resahpe wide trimmed data into long form
age<-rep('Juv',length(sv.wide$Family)) 
age2<-rep('Adult',length(sv.wide$Family)) # adult
long.all<-data.frame(sv.wide[,1], age, sv.wide[c(6,8)])
colnames(long.all)<-c("Fams",'age','P','mins')
long.all.1<-data.frame(sv.wide[,1], age2, sv.wide[c(7,8)])
colnames(long.all.1)<-c("Fams",'age','P','mins')
long.all<-rbind(long.all, long.all.1)
head(long.all)


interact.plotB <-  ggplot(long.all, aes(x=age, y=P, group=as.factor(Fams)) )+ 
                   geom_line(color='blue') +
                   geom_point()+
                   ylim(0,.75)+
                   labs(title="", x="Host Age", y = "Proportion infected")+
                   theme_classic()+
                   theme(legend.position = 'none')+
                   theme(text = element_text(size = 14))+
                   scale_x_discrete(labels = c('Juv.','Adult'))

interact.plotB 

library("ggpubr")
figure2 <- ggarrange(boxplotB, corr.plotB, interact.plotB, labels = c("A", "B","C"))#, nrow=1)
figure2


