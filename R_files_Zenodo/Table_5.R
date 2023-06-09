
# Script to reproduce information in Table 5 from:
#
# Cowling BJ, Ng DMW, Ip DKM, Liao Q, Lam WWT, Wu JT, et al.
# Community Psychological and Behavioral Responses through the First Wave
# of the 2009 Influenza A(H1N1) Pandemic in Hong Kong
# JID 2010; 202(6):867-876.


source("dataframe.r")

require(MASS)

# avoid eating out
fit.m <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
  fit.m[[i]] <- glmmPQL(pm10_1~sex+factor(age_gp)+factor(edu)+factor(ax1_score)+
  factor(ph1)+factor(bf1)+factor(bf2)+factor(bf5)+factor(bf4)+bf3_1+bf3_3+bf3_4,
  random = ~ 1 | wave, family = binomial,data=subset.nomiss[[i]])
}
round(combine.mi(fit.m,10),2)

# avoid using public transport
fit.m <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
  fit.m[[i]] <- glmmPQL(pm10_2~sex+factor(age_gp)+factor(edu)+factor(ax1_score)+
  factor(ph1)+factor(bf1)+factor(bf2)+factor(bf5)+factor(bf4)+bf3_1+bf3_3+bf3_4,
  random = ~ 1 | wave, family = binomial,data=subset.nomiss[[i]])
}
round(combine.mi(fit.m,10),2)

# avoid crowded places due to swine flu
fit.m <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
  fit.m[[i]] <- glmmPQL(pm10_3~sex+factor(age_gp)+factor(edu)+factor(ax1_score)+
  factor(ph1)+factor(bf1)+factor(bf2)+factor(bf5)+factor(bf4)+bf3_1+bf3_3+bf3_4,
  random = ~ 1 | wave, family = binomial,data=subset.nomiss[[i]])
}
round(combine.mi(fit.m,10),2)

# reschedule travel plan
fit.m <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
  fit.m[[i]] <- glmmPQL(pm10_5~sex+factor(age_gp)+factor(edu)+factor(ax1_score)+
  factor(ph1)+factor(bf1)+factor(bf2)+factor(bf5)+factor(bf4)+bf3_1+bf3_3+bf3_4,
  random = ~ 1 | wave, family = binomial,data=subset.nomiss[[i]])
}
round(combine.mi(fit.m,10),2)

# End of script
