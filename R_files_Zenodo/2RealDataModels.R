###Run models on real data 

options(stringsAsFactors = FALSE)

modeldata<-read.csv("modeldata.csv")

modeldatasexed<-modeldata[-which(is.na(modeldata$sex)),]

modelfocal<-modeldatasexed[-which(is.na(modeldatasexed$proportion_groups)),]

library(logistf)

#effect on degree on disease
degree_log<-logistf(formula = tsd ~ avg_degree + sex, data = modeldatasexed, pl = TRUE,
                    alpha = 0.05, firth = TRUE)

summary(degree_log)

infdegree_log<-logistf(formula = tsd ~ proportion_infected + sex, data = modeldatasexed, pl = TRUE,
                       alpha = 0.05, firth = TRUE)

summary(infdegree_log)

social_log<-logistf(formula = tsd ~ proportion_social + sex, data = modelfocal, pl = TRUE,
                    alpha = 0.05, firth = TRUE)

summary(social_log)

groups_log<-logistf(formula = tsd ~ proportion_groups + sex, data = modelfocal, pl = TRUE,
                    alpha = 0.05, firth = TRUE)

summary(groups_log)
