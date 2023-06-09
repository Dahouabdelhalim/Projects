################ Dating Violence Questionnaire (DVQ, ex-Cuvino) 
############## variable computations

# carico le librerie
library(psych)
library(lavaan)


############### loading test-retest datafile
detach(retest.cuvino.data)
retest.cuvino.data <- read.table("total_sample_test_retest_PLOSONE.csv", sep=",", header=T)
attach(retest.cuvino.data)
head(retest.cuvino.data) # visualizza il file dati
names(retest.cuvino.data)


### variable labels 
# "cod" = participant code                   
# "eta" = age of participants                 
# "gender" = gender of participant (value label: 1'M' 2'F')                
# "anno" = year of school/university                
#  "psicologia" = university 
# "redditomensile" = income 
# "lavororetribuito" = type of job    
# "sessopartner" =          
# "etapartner" =             
# "cuv1a" to "cuv61" = DVQ items at the first observation                
# "ey1" to "ey48" = EPQ items                 
# "st1" to  "es7" = BRS items
# "source01" = participants who agreed to partecipate in the retest
# "source02" = participants with test-retest measures             
# "cuv1a_r" to "cuv61_r" = DVQ items at the second observation              



##################### Analysis 

############# Paragraph: "Participants"
### First observation: TEST
### distribution of gender
### value label gender 1'M' 2'F'.
table(gender)

### Mean age
describe(eta)
### Mean age as function of gender
by(eta, gender, describe)

### counting how many missing values at the DVQ for each participant
item.nr<-10
retest.cuvino.data$miss.count<-0
table(retest.cuvino.data$miss.count)
for (i in 1:42) {
  print(i)
  retest.cuvino.data$miss.count<-retest.cuvino.data$miss.count+(is.na(retest.cuvino.data[,item.nr]==T)*1)
  item.nr<-item.nr+2
}
### number of subjects with less than 10 missing responses on DVQ
table(retest.cuvino.data$miss.count<=10)
### excluding participants with more than 10 missing responses on DVQ
retest.cuvino.data <- subset(retest.cuvino.data, retest.cuvino.data$miss.count<=10) 

### number of participants who agree to participate to the retest
table(retest.cuvino.data$source01)

### number of participants for which both observations are available:Test-Retest
table(retest.cuvino.data$source02)

### Age Mean and SD for participants included in the retest sample
by(eta, source01, describe)
table(gender, source01)
table(gender, source02)

################### paragraph: "RESULTS"
### Shapiro-Wilk test for normality of distribution 
item.nr<-10
for (i in 1:42) {
  print(i)
  print(shapiro.test(retest.cuvino.data[,item.nr])) 
  item.nr<-item.nr+2
}

### frequencies 
# item 43: Have you ever been afraid of your partner?
table(retest.cuvino.data[,'cuv43'])

# Item 44: Did you ever feel the sensation of being trapped in your romantic relationship?
table(retest.cuvino.data[,'cuv44'])

# Item 45: Did you ever feel maltreated by your partner?
table(retest.cuvino.data[,'cuv45'])

# item 46: Do you know someone (a friend of yours or a neighborhood that was maltreated by the partner during their romantic relationship?
table(retest.cuvino.data[,'cuv46'])

# average relationship duration
describe(retest.cuvino.data[,'cuv48'])

# average age of onset of the relationship
describe(retest.cuvino.data[,'cuv49'])

# number of attempts to break the romantic relationship with the partner
table(retest.cuvino.data[,'cuv56']>=1)

# help from someone in order to break his/her relationships
table(retest.cuvino.data[,'cuv58'])

# how long it takes to break the relationship
describe(retest.cuvino.data[,'cuv57'])



######################################## CFA - BEGIN
### NOTE TO THE EDITOR: results concerning CFA and reported in the article 
### in the paragraph "Confirmatory Factor Analysis" was performed with M-PLUS. 
### The code is not included.

####################### Paragraph: Reliabilities and Construct Validity of the 8-factors model
### FIRST observations
### DESAPEGO
names(retest.cuvino.data) # only first observation
data.desapego<-cbind(retest.cuvino.data$cuv6a,retest.cuvino.data$cuv14a,retest.cuvino.data$cuv22a,
                     retest.cuvino.data$cuv30a,retest.cuvino.data$cuv32a,retest.cuvino.data$cuv33a,
                     retest.cuvino.data$cuv37a)
data.desapego<-as.matrix(data.frame(data.desapego))
alpha(data.desapego)
retest.cuvino.data$desapego<- rowSums(data.desapego) 
describe(retest.cuvino.data$desapego)

### HUMILIACION
data.humiliacion <- cbind(retest.cuvino.data$cuv7a,retest.cuvino.data$cuv15a,retest.cuvino.data$cuv23a,
                          retest.cuvino.data$cuv31a,retest.cuvino.data$cuv36a,retest.cuvino.data$cuv40a,
                          retest.cuvino.data$cuv41a)
data.humiliacion<-as.matrix(data.frame(data.humiliacion))
alpha(data.humiliacion)
retest.cuvino.data$humiliacion<- rowSums(data.humiliacion) 
describe(retest.cuvino.data$humiliacion)

### SEXUAL 
data.sexual <- cbind(retest.cuvino.data$cuv2a,retest.cuvino.data$cuv10a,retest.cuvino.data$cuv18a,
                     retest.cuvino.data$cuv26a,retest.cuvino.data$cuv34a,retest.cuvino.data$cuv39a) 
data.sexual<-as.matrix(data.frame(data.sexual))
alpha(data.sexual)
retest.cuvino.data$sexual<- rowSums(data.sexual) 
describe(retest.cuvino.data$sexual)

### COERCION 
data.coercion <-cbind(retest.cuvino.data$cuv1a, retest.cuvino.data$cuv9a, retest.cuvino.data$cuv17a, 
                      retest.cuvino.data$cuv25a, retest.cuvino.data$cuv38a, retest.cuvino.data$cuv42a)
data.coercion<-as.matrix(data.frame(data.coercion))
alpha(data.coercion)
retest.cuvino.data$coercion<- rowSums(data.coercion) 
describe(retest.cuvino.data$coercion)

### PHYSIC
data.fisic <- cbind(retest.cuvino.data$cuv5a, retest.cuvino.data$cuv13a, retest.cuvino.data$cuv20a, 
                    retest.cuvino.data$cuv21a, retest.cuvino.data$cuv29a)
data.fisic<-as.matrix(data.frame(data.fisic))
alpha(data.fisic)
retest.cuvino.data$fisic<- rowSums(data.fisic) 
describe(retest.cuvino.data$fisic)

### GENERO
data.genero <-cbind(retest.cuvino.data$cuv3a,retest.cuvino.data$cuv11a,retest.cuvino.data$cuv19a,
                    retest.cuvino.data$cuv27a,retest.cuvino.data$cuv35a)
data.genero<-as.matrix(data.frame(data.genero))
alpha(data.genero)
retest.cuvino.data$genero<- rowSums(data.genero) 
describe(retest.cuvino.data$genero)

### CASTIGO EMOTIVO
data.castigo_em <-cbind(retest.cuvino.data$cuv8a,retest.cuvino.data$cuv16a,retest.cuvino.data$cuv24a)
data.castigo_em<-as.matrix(data.frame(data.castigo_em))
alpha(data.castigo_em)
retest.cuvino.data$castigo_em<- rowSums(data.castigo_em) 
describe(retest.cuvino.data$castigo_em)

### INSTRUMENTAL
data.instrumental <- cbind(retest.cuvino.data$cuv12a,retest.cuvino.data$cuv4a,retest.cuvino.data$cuv28a)
data.instrumental<-as.matrix(data.frame(data.instrumental))
alpha(data.instrumental)
retest.cuvino.data$instrumental<- rowSums(data.instrumental) 
describe(retest.cuvino.data$instrumental)

### Spearman correlations 
construct.validity.data<-cbind(retest.cuvino.data$desapego,
                               retest.cuvino.data$humiliacion,
                               retest.cuvino.data$sexual,
                               retest.cuvino.data$coercion,
                               retest.cuvino.data$fisic,
                               retest.cuvino.data$genero,
                               retest.cuvino.data$castigo_em,
                               retest.cuvino.data$instrumental)

round(cor(construct.validity.data, method=c("spearman")),3)



############## SECOND OBSERVATION: RETEST
### selecting those participants that have both measures
table(retest.cuvino.data$source02)
test.retest.data <- subset(retest.cuvino.data, retest.cuvino.data$source02==1)
table(test.retest.data$source02)

### DESAPEGO
names(test.retest.data) # first and second observations
data.desapego_r<-cbind(test.retest.data$cuv6a_r,
                       test.retest.data$cuv14a_r,
                       test.retest.data$cuv22a_r,
                       test.retest.data$cuv30a_r,
                       test.retest.data$cuv32a_r,
                       test.retest.data$cuv33a_r,
                       test.retest.data$cuv37a_r)
#data.desapego<-cbind(cfa.data$cuv6a,cfa.data$cuv14a,cfa.data$cuv22a,cfa.data$cuv30a,cfa.data$cuv32a,cfa.data$cuv33a,cfa.data$cuv37a)
data.desapego_r<-as.matrix(data.frame(data.desapego_r))
alpha(data.desapego_r)
test.retest.data$desapego_r<- rowSums(data.desapego_r) 
describe(test.retest.data$desapego_r)

### HUMILIACION
data.humiliacion_r <- cbind(test.retest.data$cuv7a_r,
                            test.retest.data$cuv15a_r,
                            test.retest.data$cuv23a_r,
                            test.retest.data$cuv31a_r,
                            test.retest.data$cuv36a_r,
                            test.retest.data$cuv40a_r,
                            test.retest.data$cuv41a_r)
data.humiliacion_r<-as.matrix(data.frame(data.humiliacion_r))
alpha(data.humiliacion_r)
test.retest.data$humiliacion_r<- rowSums(data.humiliacion_r) 
describe(test.retest.data$humiliacion_r)

### SEXUAL 
data.sexual_r <- cbind(test.retest.data$cuv2a_r,
                       test.retest.data$cuv10a_r,
                       test.retest.data$cuv18a_r,
                       test.retest.data$cuv26a_r,
                       test.retest.data$cuv34a_r,
                       test.retest.data$cuv39a_r) 
data.sexual_r<-as.matrix(data.frame(data.sexual_r))
alpha(data.sexual_r)
test.retest.data$sexual_r<- rowSums(data.sexual_r) 
describe(test.retest.data$sexual_r)

### COERCION 
data.coercion_r <-cbind(test.retest.data$cuv1a_r, 
                        test.retest.data$cuv9a_r, 
                        test.retest.data$cuv17a_r, 
                        test.retest.data$cuv25a_r, 
                        test.retest.data$cuv38a_r, 
                        test.retest.data$cuv42a_r)
data.coercion_r<-as.matrix(data.frame(data.coercion_r))
alpha(data.coercion_r)
test.retest.data$coercion_r<- rowSums(data.coercion_r) 
describe(test.retest.data$coercion_r)

### PHYSIC
data.fisic_r <- cbind(test.retest.data$cuv5a_r, 
                      test.retest.data$cuv13a_r, 
                      test.retest.data$cuv20a_r, 
                      test.retest.data$cuv21a_r, 
                      test.retest.data$cuv29a_r)
data.fisic_r<-as.matrix(data.frame(data.fisic_r))
alpha(data.fisic_r)
test.retest.data$fisic_r<- rowSums(data.fisic_r) 
describe(test.retest.data$fisic_r)

### GENERO
data.genero_r <-cbind(test.retest.data$cuv3a_r,
                      test.retest.data$cuv11a_r,
                      test.retest.data$cuv19a_r,
                      test.retest.data$cuv27a_r,
                      test.retest.data$cuv35a_r)
data.genero_r<-as.matrix(data.frame(data.genero_r))
alpha(data.genero_r)
test.retest.data$genero_r<- rowSums(data.genero_r) 
describe(test.retest.data$genero_r)

### CASTIGO EMOTIVO
data.castigo_em_r <-cbind(test.retest.data$cuv8a_r,
                          test.retest.data$cuv16a_r,
                          test.retest.data$cuv24a_r)
data.castigo_em_r<-as.matrix(data.frame(data.castigo_em_r))
alpha(data.castigo_em_r)
test.retest.data$castigo_em_r<- rowSums(data.castigo_em_r) 
describe(test.retest.data$castigo_em_r)

### INSTRUMENTAL
data.instrumental_r <- cbind(test.retest.data$cuv12a_r,
                             test.retest.data$cuv4a_r,
                             test.retest.data$cuv28a_r)
data.instrumental_r<-as.matrix(data.frame(data.instrumental_r))
alpha(data.instrumental_r)
test.retest.data$instrumental_r<- rowSums(data.instrumental_r) 
describe(test.retest.data$instrumental_r)

### Spearman correlations: SECOND OBSERVATION
construct.validity.data.r<-cbind(test.retest.data$desapego_r,
                                 test.retest.data$humiliacion_r,
                                 test.retest.data$sexual_r,
                                 test.retest.data$coercion_r,
                                 test.retest.data$fisic_r,
                                 test.retest.data$genero_r,
                                 test.retest.data$castigo_em_r,
                                 test.retest.data$instrumental_r)

round(cor(construct.validity.data.r, method=c("spearman"), use="complete.obs"),3)

### test-retest reliability
construct.validity.data.tt1<-cbind(test.retest.data$desapego,
                                 test.retest.data$humiliacion,
                                 test.retest.data$sexual,
                                 test.retest.data$coercion,
                                 test.retest.data$fisic,
                                 test.retest.data$genero,
                                 test.retest.data$castigo_em,
                                 test.retest.data$instrumental)

round(cor(construct.validity.data.tt1, construct.validity.data.r, method=c("spearman"), use="complete.obs"),3)

# confidence interval for Bonferroni correction (alpha=0.05; number of comparisons=28; p-corrected=0.002)
### first wave
rho<-0.158
n=383
p.corrected=0.998
r.con(rho,n,p.corrected,twotailed=TRUE)

### second wave
# confidence interval for Bonferroni correction (alpha=0.05; number of comparisons=28; p-corrected=0.002)
rho<-0.26
n=139
p.corrected=0.998
r.con(rho,n,p.corrected,twotailed=TRUE)


################ MEANS and SD as function of gender and age groups: table 4
# two age groups: less than 21 and over 21
names(retest.cuvino.data)
retest.cuvino.data$age.group<-(retest.cuvino.data$eta<=21)*1+(retest.cuvino.data$eta>21)*2
table(retest.cuvino.data$age.group)

table(retest.cuvino.data$gender,retest.cuvino.data$age.group)
by(construct.validity.data, list(retest.cuvino.data$gender,retest.cuvino.data$age.group), describe)



########################## CONSTRUCT VALIDITY OF CUVINO 

########### EPQ 
# Compute EPQ_PSY = ey10+ey14+ey22+ey31+ey39+(1-ey2)+(1-ey6)+(1-ey18)+(1-ey26)+(1-ey28)+(1-ey35)+(1-ey43).
# Compute EPQ_EXT = ey3+ey7+ey11+ey15+ey19+ey23+ey32+ey36+ey44+ey48+(1-ey27)+(1-ey41).
# Compute EPQ_NEV= ey1+ey5+ey9+ey13+ey17+ey21+ey25+ey30+ey34+ey38+ey42+ey46.
# Compute EPQ_LIE= ey4+ey16+ey45+(1-ey8)+(1-ey12)+(1-ey20)+(1-ey24)+(1-ey29)+(1-ey33)+(1-ey37)+(1-ey40)+(1-ey47).
# CORRELATIONS EPQ_PSY EPQ_EXT EPQ_NEV EPQ_LIE.

retest.cuvino.data$epq.psy <- retest.cuvino.data[, c('ey10')] + 
  retest.cuvino.data[, c('ey14')] + 
  retest.cuvino.data[, c('ey22')] + 
  retest.cuvino.data[, c('ey31')] + 
  retest.cuvino.data[, c('ey39')] + 
  (1-retest.cuvino.data[, c('ey2')]) + 
  (1-retest.cuvino.data[, c('ey6')]) + 
  (1-retest.cuvino.data[, c('ey18')]) + 
  (1-retest.cuvino.data[, c('ey26')]) + 
  (1-retest.cuvino.data[, c('ey28')]) + 
  (1-retest.cuvino.data[, c('ey35')]) + 
  (1-retest.cuvino.data[, c('ey43')])
describe(retest.cuvino.data$epq.psy)

retest.cuvino.data$epq.ext <- retest.cuvino.data[, c('ey3')] + 
  retest.cuvino.data[, c('ey7')] + 
  retest.cuvino.data[, c('ey11')] + 
  retest.cuvino.data[, c('ey15')] + 
  retest.cuvino.data[, c('ey19')] + 
  retest.cuvino.data[, c('ey23')] + 
  retest.cuvino.data[, c('ey32')] + 
  retest.cuvino.data[, c('ey36')] + 
  retest.cuvino.data[, c('ey44')] + 
  retest.cuvino.data[, c('ey48')] + 
  (1-retest.cuvino.data[, c('ey27')]) + 
  (1-retest.cuvino.data[, c('ey41')])
describe(retest.cuvino.data$epq.ext)

retest.cuvino.data$epq.nev <- retest.cuvino.data[, c('ey1')] + 
  retest.cuvino.data[, c('ey5')] + 
  retest.cuvino.data[, c('ey9')] + 
  retest.cuvino.data[, c('ey13')] + 
  retest.cuvino.data[, c('ey17')] + 
  retest.cuvino.data[, c('ey21')] + 
  retest.cuvino.data[, c('ey25')] + 
  retest.cuvino.data[, c('ey30')] + 
  retest.cuvino.data[, c('ey34')] + 
  retest.cuvino.data[, c('ey38')] + 
  retest.cuvino.data[, c('ey42')] + 
  retest.cuvino.data[, c('ey46')]
describe(retest.cuvino.data$epq.nev)

retest.cuvino.data$epq.lie <- retest.cuvino.data[, c('ey4')] + 
  retest.cuvino.data[, c('ey16')] + 
  retest.cuvino.data[, c('ey45')] + 
  (1-retest.cuvino.data[, c('ey8')])+
  (1-retest.cuvino.data[, c('ey12')])+
  (1-retest.cuvino.data[, c('ey20')])+
  (1-retest.cuvino.data[, c('ey24')])+
  (1-retest.cuvino.data[, c('ey29')])+
  (1-retest.cuvino.data[, c('ey33')])+
  (1-retest.cuvino.data[, c('ey37')])+
  (1-retest.cuvino.data[, c('ey40')])+
  (1-retest.cuvino.data[, c('ey47')])
describe(retest.cuvino.data$epq.lie)

epq.data<-cbind(retest.cuvino.data$epq.psy,
                retest.cuvino.data$epq.ext,
                retest.cuvino.data$epq.nev,
                retest.cuvino.data$epq.lie)
round(cor(epq.data, method=c("spearman"), use="complete.obs"),3)

round(cor(construct.validity.data, epq.data, method=c("spearman"), use="complete.obs"),3)

# confidence interval for Bonferroni correction (alpha=0.05; number of comparisons=28; p-corrected=0.002)
rho<-0.140
n=383
p=0.99
r.con(rho,n,p,twotailed=TRUE)

#######  Breaking Relationship Scale (BRS)
data.stalk <-cbind(retest.cuvino.data[,'brs1'],
                     retest.cuvino.data[,'brs2'],
                     retest.cuvino.data[,'brs3'],
                     retest.cuvino.data[,'brs4'],
                     retest.cuvino.data[,'brs5'],
                     retest.cuvino.data[,'brs6'],
                     retest.cuvino.data[,'brs7'])
data.stalk<-as.matrix(data.frame(data.stalk))
alpha(data.stalk)

stalk.ind<-rowSums(data.stalk)
describe(stalk.ind)

fit.stalk.ind.cuv.0 <- lm(stalk.ind ~epq.psy+epq.ext+epq.nev, 
                          data=retest.cuvino.data)
summary(fit.stalk.ind.cuv.0) # show results
fit.stalk.ind.cuv.1 <- lm(stalk.ind ~epq.psy+epq.ext+epq.nev+ 
                            desapego+ humiliacion+ sexual+ coercion+ fisic+ genero+ castigo_em+ instrumental, 
                          data=retest.cuvino.data)
summary(fit.stalk.ind.cuv.1) # show results

fit.stalk.ind.cuv.1.std <- lm(scale(stalk.ind) ~scale(epq.psy)+scale(epq.ext)+scale(epq.nev)+ 
                                scale(desapego)+ scale(humiliacion)+ scale(sexual)+ scale(coercion)+ scale(fisic)+ scale(genero)+ scale(castigo_em)+ scale(instrumental), 
                              data=retest.cuvino.data)
summary(fit.stalk.ind.cuv.1.std) # show results
