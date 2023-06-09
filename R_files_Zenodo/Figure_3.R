### FIGURE 3A ###

#Kaplan-Meier, file: survival.csv
library(survival)
library(survminer)

larvalongev<-Surv(time=data_survival$time,event=data_survival$dead)
larvalongev

KM_larvalongev<-survfit(larvalongev~data_survival$diet)
KM_larvalongev
plot(KM_larvalongev,conf.int = "none",xlab="days",ylab="prop. alive",lty=1:4)

# Log-Rank 
survdiff(larvalongev~data_survival$diet) 
# pairwise comparisons
paircomp<- pairwise_survdiff(Surv(time, dead) ~ diet,data = data_survival) 
paircomp

# Larval masses, file: weight.csv
kruskal.test(data_weight$weight~data_weight$diet,data=data_weight)
# pairwise comparisons
library("dunn.test")
library(FSA) 
dunnTest(data_weight$weight~data_weight$diet,data=data_weight,
         method="bonferroni")



### FIGURE 3B ###

# Proportion surviving, file: caterpillar_mass&survival.csv
larvaldata$died_1_alive_0<-as.numeric(larvaldata$died_1_alive_0)
chisq.test(larvaldata$died_1_alive_0 , larvaldata$diet, correct = FALSE)

# Larval masses, file: caterpillar_mass&survival.csv
kruskal.test(final_fw_mg ~ diet,data=larvaldata)
dunnTest(final_fw_mg ~ diet,data=larvaldata,
         method="bonferroni")



### FIGURE 3C ###

# Hyles 1 larva
#Kaplan-Meier, file: survival1.csv
library(survival)
library(survminer)
larvalongev1<-Surv(time=data_survival1$day,event=data_survival1$dead)
larvalongev1

KM_larvalongev1<-survfit(larvalongev1~data_survival1$diet)
KM_larvalongev1

#Log-Rank 
survdiff(larvalongev1~data_survival1$diet)
# pairwise comparisons
larva1<- pairwise_survdiff(Surv(day, dead) ~ diet,data = data_survival1) 
larva1

# Larval masses,file: weight_1cat.csv
kruskal.test(data_weight1$weight~data_weight1$diet,data=data_weight1)
# pairwise comparisons
library("dunn.test")
library(FSA) 
dunnTest(data_weight1$weight~data_weight1$diet,data=data_weight1,
         method="bonferroni")



# Hyles 4 larvae
#Kaplan-Meier, file: survival4.csv
library(survival)
library(survminer)
larvalongev4<-Surv(time=data_survival4$week,event=data_survival4$dead)
larvalongev4

KM_larvalongev4<-survfit(larvalongev4~data_survival1$diet)
KM_larvalongev4

#Log-Rank 
survdiff(larvalongev4~data_survival1$diet)
# pairwise comparisons
larva1<- pairwise_survdiff(Surv(day, dead) ~ diet,data = data_survival4) 
larva1


# Larval masses,file: weight_4cat.csv
kruskal.test(data_weight4$weight~data_weight4$diet,data=data_weight4)
# pairwise comparisons
library("dunn.test")
library(FSA) 
dunnTest(data_weight4$weight~data_weight4$diet,data=data_weight4,
         method="bonferroni")