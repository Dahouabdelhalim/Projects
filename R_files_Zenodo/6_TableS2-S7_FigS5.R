library(foreign)
library(stargazer)
library(sandwich)
library(lmtest)
#library(tidyverse)
#library(caret)
library(regclass)

### 1. PREPARE DATA TABLE
dataCore=read.table("MASZK_dataCore.csv", sep=",", header = T)


dataCore$pf_deny=0
dataCore$pf_deny[dataCore$PfizerRating==1]=1
dataCore$mo_deny=0
dataCore$mo_deny[dataCore$ModernaRating==1]=1
dataCore$as_deny=0
dataCore$as_deny[dataCore$AstraRating==1]=1
dataCore$sp_deny=0
dataCore$sp_deny[dataCore$SputRating==1]=1
dataCore$si_deny=0
dataCore$si_deny[dataCore$SinoRating==1]=1


dataCore$hedu=0
dataCore$hedu[dataCore$Schooling==4]=1

# handle NA
dataCore$Covid[is.na(dataCore$Covid)]=0

dataCore$SerCovid[is.na(dataCore$SerCovid)]=0

# add month factor ???

dataCore$Vaccine1Date

# standardize

data_reg=dataCore[,c("pf_deny", "mo_deny", "as_deny", "sp_deny", "si_deny", "VaccineDoctor", "VaccineScientist", "VaccineSceptical", "VaccinePolitician", "VaccineFamily", "VaccineFriends",
                     "VaccineJournalist", "VaccineCelebrity", "Age", "IsFemale", "hedu", "CitySize", "WealthPreCovid", "Smoking",
                       "ChronicIll", "Covid",  "SerCovid")]
sc_data_reg=as.data.frame(scale(data_reg))

### 2. Prepare Empirics

# LPM on unacceptable

md_pf=lm(data=dataCore, pf_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         #+ WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)
md_mo=lm(data=dataCore, mo_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         #+ WealthPreCovid 
         + Smoking + ChronicIll  
         + Covid + SerCovid
)
md_as=lm(data=dataCore, as_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         #+ WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)
md_sp=lm(data=dataCore, sp_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         #+ WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)
md_si=lm(data=dataCore, si_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         #+ WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)

cov_pf         <- vcovHC(md_pf, type = "HC1")
robust_pf    <- sqrt(diag(cov_pf))

cov_mo         <- vcovHC(md_mo, type = "HC1")
robust_mo    <- sqrt(diag(cov_mo))

cov_as         <- vcovHC(md_as, type = "HC1")
robust_as    <- sqrt(diag(cov_as))

cov_sp         <- vcovHC(md_sp, type = "HC1")
robust_sp    <- sqrt(diag(cov_sp))

cov_si         <- vcovHC(md_si, type = "HC1")
robust_si    <- sqrt(diag(cov_si))


# Stargazer output (with and with RSE)
stargazer(md_pf, md_mo, md_as, md_sp, md_si, 
          se=list(robust_pf, robust_mo, robust_as, robust_sp, robust_si),
          dep.var.labels="Vaccine Denial, linear probability models",
          column.labels = c("Pfizer","Moderna", "Astrazeneca", "Sputnik", "Sinopharm"),
          covariate.labels=c("Advice from Doctors", "Advice from Scientists", "Advice from Anti-Vaccine Propagators", "Advice from Politicians", 
                             "Advice from Family", "Advice from Friends", "Advice from Journalists", "Advice from Celebrities",
                             "Age", "Female", "University", "City Category by Size", "Smoking", "Chronic Illness", "Covid Previously",
                             "Serious Covid Previously"),
          omit.stat=c("LL","ser","f"), ci=F,  single.row=F)

# Standardized table
st_pf=lm(data=sc_data_reg, pf_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         #+ WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)
st_mo=lm(data=sc_data_reg, mo_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         #+ WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)
st_as=lm(data=sc_data_reg, as_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         #+ WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)
st_sp=lm(data=sc_data_reg, sp_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         #+ WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)
st_si=lm(data=sc_data_reg, si_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         #+ WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)

cov_st_pf         <- vcovHC(st_pf, type = "HC1")
robust_st_pf    <- sqrt(diag(cov_st_pf))

cov_st_mo         <- vcovHC(st_mo, type = "HC1")
robust_st_mo    <- sqrt(diag(cov_st_mo))

cov_st_as         <- vcovHC(st_as, type = "HC1")
robust_st_as    <- sqrt(diag(cov_st_as))

cov_st_sp         <- vcovHC(st_sp, type = "HC1")
robust_st_sp    <- sqrt(diag(cov_st_sp))

cov_st_si         <- vcovHC(st_si, type = "HC1")
robust_st_si    <- sqrt(diag(cov_st_si))

stargazer(st_pf, st_mo, st_as, st_sp, st_si, 
          dep.var.labels="Vaccine Denial, linear probability models",
          se=list(robust_st_pf, robust_st_mo, robust_st_as, robust_st_sp, robust_st_si),
          column.labels = c("Pfizer","Moderna", "Astrazeneca", "Sputnik", "Sinopharm"),
          covariate.labels=c("Advice from Doctors", "Advice from Scientists", "Advice from Anti-Vaccine Propagators", "Advice from Politicians", 
                             "Advice from Family", "Advice from Friends", "Advice from Journalists", "Advice from Celebrities",
                             "Age", "Female", "University", "City Category by Size", "Smoking", "Chronic Illness", "Covid Previously", "Serious Covid Previously"),
          omit.stat=c("LL","ser","f"), ci=F,  single.row=F)

# Logit on unacceptable

d_r$pf_deny_d=as.factor(d_r$pf_deny)

table(d_r$pf_deny)

md_log_pf=glm(data=data_reg, pf_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
                Age + IsFemale + hedu + CitySize 
              #+ WealthPreCovid 
              + Smoking + ChronicIll 
              + Covid + SerCovid, family = "binomial"
)
md_log_mo=glm(data=data_reg, mo_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
                Age + IsFemale + hedu + CitySize 
              #+ WealthPreCovid 
              + Smoking + ChronicIll 
              + Covid + SerCovid, family = "binomial"
)
md_log_as=glm(data=data_reg, as_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
                Age + IsFemale + hedu + CitySize 
              #+ WealthPreCovid 
              + Smoking + ChronicIll 
              + Covid + SerCovid, family = "binomial"
)
md_log_sp=glm(data=data_reg, sp_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
                Age + IsFemale + hedu + CitySize 
              #+ WealthPreCovid 
              + Smoking + ChronicIll 
              + Covid + SerCovid, family = "binomial"
)
md_log_si=glm(data=data_reg, si_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
                Age + IsFemale + hedu + CitySize 
              #+ WealthPreCovid 
              + Smoking + ChronicIll 
              + Covid + SerCovid, family = "binomial"
)

cov_md_log_pf         <- vcovHC(md_log_pf, type = "HC1")
robust_log_pf    <- sqrt(diag(cov_md_log_pf))

cov_md_log_mo         <- vcovHC(md_log_mo, type = "HC1")
robust_log_mo    <- sqrt(diag(cov_md_log_mo))

cov_md_log_as         <- vcovHC(md_log_as, type = "HC1")
robust_log_as    <- sqrt(diag(cov_md_log_as))

cov_md_log_sp         <- vcovHC(md_log_sp, type = "HC1")
robust_log_sp    <- sqrt(diag(cov_md_log_sp))

cov_md_log_si         <- vcovHC(md_log_si, type = "HC1")
robust_log_si    <- sqrt(diag(cov_md_log_si))



stargazer(md_log_pf, md_log_mo, md_log_as, md_log_sp, md_log_si, 
          se=list(robust_log_pf, robust_log_mo, robust_log_as, robust_log_sp, robust_log_si),
          dep.var.labels="",
          column.labels = c("Pfizer","Moderna", "Astrazeneca", "Sputnik", "Sinopharm"),
          covariate.labels=c("Advice from Doctors", "Advice from Scientists", "Advice from Anti-Vaccine Propagators", "Advice from Politicians", 
                             "Advice from Family", "Advice from Friends", "Advice from Journalists", "Advice from Celebrities",
                             "Age", "Female", "University", "City Category by Size", "Smoking", "Chronic Illness", "Covid Previously", "Serious Covid Previously"),
          omit.stat=c("LL","ser","f"), ci=F,  single.row=F)



## ROBUSTNESS

# those who are stable in assessment

md_pf_ch=lm(data=dataCore[dataCore$PfizerChange==0,], pf_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
              Age + IsFemale + hedu + CitySize 
            #+ WealthPreCovid 
            + Smoking + ChronicIll 
            + Covid + SerCovid
)
md_mo_ch=lm(data=dataCore[dataCore$ModernaChange==0,], mo_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
              Age + IsFemale + hedu + CitySize 
            #+ WealthPreCovid 
            + Smoking + ChronicIll 
            + Covid + SerCovid
)
md_as_ch=lm(data=dataCore[dataCore$AstraZenecaChange==0,], as_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
              Age + IsFemale + hedu + CitySize 
            #+ WealthPreCovid 
            + Smoking + ChronicIll 
            + Covid + SerCovid
)
md_sp_ch=lm(data=dataCore[dataCore$SputnikChange==0,], sp_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
              Age + IsFemale + hedu + CitySize 
            #+ WealthPreCovid 
            + Smoking + ChronicIll 
            + Covid + SerCovid
)
md_si_ch=lm(data=dataCore[dataCore$SinopharmChange==0,], si_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
              Age + IsFemale + hedu + CitySize 
            #+ WealthPreCovid 
            + Smoking + ChronicIll 
            + Covid + SerCovid
)

cov_ch_pf         <- vcovHC(md_pf_ch, type = "HC1")
robust_ch_pf    <- sqrt(diag(cov_ch_pf))

cov_ch_mo         <- vcovHC(md_mo_ch, type = "HC1")
robust_ch_mo    <- sqrt(diag(cov_ch_mo))

cov_ch_as         <- vcovHC(md_as_ch, type = "HC1")
robust_ch_as    <- sqrt(diag(cov_ch_as))

cov_ch_sp         <- vcovHC(md_sp_ch, type = "HC1")
robust_ch_sp    <- sqrt(diag(cov_ch_sp))

cov_ch_si         <- vcovHC(md_si_ch, type = "HC1")
robust_ch_si    <- sqrt(diag(cov_ch_si))

stargazer(md_pf_ch, md_mo_ch, md_as_ch, md_sp_ch, md_si_ch, 
          se=list(robust_ch_pf, robust_ch_mo, robust_ch_as, robust_ch_sp, robust_ch_si),
          dep.var.labels="",
          column.labels = c("Pfizer","Moderna", "Astrazeneca", "Sputnik", "Sinopharm"),
          covariate.labels=c("Advice from Doctors", "Advice from Scientists", "Advice from Anti-Vaccine Propagators", "Advice from Politicians", 
                             "Advice from Family", "Advice from Friends", "Advice from Journalists", "Advice from Celebrities",
                             "Age", "Female", "University", "City Category by Size", "Smoking", "Chronic Illness", "Covid Previously", "Serious Covid Previously"),         
          omit.stat=c("LL","ser","f"), ci=F,  single.row=F)

# Wealth
w_pf=lm(data=dataCore, pf_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         + WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)
w_mo=lm(data=dataCore, mo_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         + WealthPreCovid 
         + Smoking + ChronicIll  
         + Covid + SerCovid
)
w_as=lm(data=dataCore, as_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         + WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)
w_sp=lm(data=dataCore, sp_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         + WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)
w_si=lm(data=dataCore, si_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         + WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)

cov_w_pf         <- vcovHC(w_pf, type = "HC1")
robust_w_pf    <- sqrt(diag(cov_w_pf))

cov_w_mo         <- vcovHC(w_mo, type = "HC1")
robust_w_mo    <- sqrt(diag(cov_w_mo))

cov_w_as         <- vcovHC(w_as, type = "HC1")
robust_w_as    <- sqrt(diag(cov_w_as))

cov_w_sp         <- vcovHC(w_sp, type = "HC1")
robust_w_sp    <- sqrt(diag(cov_w_sp))

cov_w_si         <- vcovHC(w_si, type = "HC1")
robust_w_si    <- sqrt(diag(cov_w_si))

stargazer(w_pf, w_mo, w_as, w_sp, w_si, 
          se=list(robust_w_pf, robust_w_mo, robust_w_as, robust_w_sp, robust_w_si),
          dep.var.labels="Vaccine Denial, linear probability models",
          column.labels = c("Pfizer","Moderna", "Astrazeneca", "Sputnik", "Sinopharm"),
          covariate.labels=c("Advice from Doctors", "Advice from Scientists", "Advice from Anti-Vaccine Propagators", "Advice from Politicians", 
                             "Advice from Family", "Advice from Friends", "Advice from Journalists", "Advice from Celebrities",
                             "Age", "Female", "University", "City Category by Size", "Wealth pre COVID", "Smoking", "Chronic Illness", 
                             "Covid Previously", "Serious Covid Previously"), 
          omit.stat=c("LL","ser","f"), ci=F,  single.row=F)


# Multicollinearity is not present - all VIF values are between 1 and 2

VIF(md_pf)


# Heteroskedasticity tests: 

#Pfizer, Moderna works with robust standard errors - significant does not change for most coefficients
coeftest(md_pf, vcov = vcovHC(md_pf, "HC1"))
coeftest(md_mo, vcov = vcovHC(md_mo, "HC1"))
coeftest(md_as, vcov = vcovHC(md_as, "HC1"))
coeftest(md_sp, vcov = vcovHC(md_sp, "HC1")) # Celebrities are not significant and smoking is significant positive with robust standard errors
coeftest(md_si, vcov = vcovHC(md_si, "HC1"))


# Prepare summary statistics, VIF, and correlation
md_pf=lm(data=dataCore, pf_deny ~ VaccineDoctor + VaccineScientist + VaccineSceptical + VaccinePolitician + VaccineFamily + VaccineFriends + VaccineJournalist + VaccineCelebrity +
           Age + IsFemale + hedu + CitySize 
         + WealthPreCovid 
         + Smoking + ChronicIll 
         + Covid + SerCovid
)

d=dataCore[,c("VaccineDoctor", "VaccineScientist", "VaccineSceptical", "VaccinePolitician", "VaccineFamily", 
              "VaccineFriends", "VaccineJournalist", "VaccineCelebrity", "Age", "IsFemale", "hedu", "CitySize", 
              "WealthPreCovid", "Smoking", "ChronicIll", "Covid", "SerCovid")]
names(d)=c("Adv. from Doct.", "Adv. from Sci.", "Adv. from Anti-vacc. Prop.", "Adv. from Pol.", "Adv. from Family", 
           "Adv. from Friends", "Adv. from Journ.", "Adv. from Celeb.", "Age", "Female", "University", "City Cat. by Size", 
           "Wealth Pre Covid", "Smoking", "Chronic Illness", "Covid Previously", "Serious Covid Previously")
sum.d=summary(d)
VIF(md_pf)

#install.packages("xtable")
library(xtable)

print(xtable(t(sum.d), type = "latex"), file = "var_sum.tex", include.rownames=FALSE)

#install.packages("corrplot")
library(corrplot)

?cor

M=cor(d, use="pairwise.complete.obs")

png("correlation.png", width=600, height=600)
corrplot(M, method = "color")
dev.off()


