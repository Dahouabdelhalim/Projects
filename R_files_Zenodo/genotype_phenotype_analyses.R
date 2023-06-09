##genotype-phenotype graphs and analyses

library(ggplot2)
library(lubridate)
library(rsq)
library(plyr)
library(dplyr)
library(gridExtra)
rm(list = ls())
setwd("")
data <- read.csv("lifehistoryfish.csv")
levels(data$sex)<- c("Female","Male")

##GLM between life history ecotype, genotype, and sex
data$lifehistory_binom<- 0
data$lifehistory_binom[data$lifehistory == "migratory"]<- 1

##include only fish that migrated (were detected at stationary antenna) Feb-May
data$antennadate <- mdy(data$antennadate)
data$month <- as.numeric(month(data$antennadate))
data <- data[(data$month %in% c(2,3,4,5)|is.na(data$month)),]

migratory <- subset(data, lifehistory=="migratory")
resident <- subset(data, lifehistory=="resident")

##summary tables
table(data$sex,data$lifehistory)
table(data$lifehistory)

##resident:
table(data$lifehistory, data$sex)
table(resident$genotype)/(sum(!(is.na(resident$genotype))))*100
table(migratory$genotype)/(sum(!(is.na(migratory$genotype))))*100

##binomial test for resident fish (predict more males):
sum(resident$sex == "F",na.rm = T)## 13 resident females
sum(resident$sex == "M",na.rm = T)## 45 resident males
binom.test(45,(13+45),0.5,"greater")
binom.test(45,(13+45),0.5,"two.sided")##CI of 65.3 to 87.7% ##report two sided test (more conservative)
45/(13+45) ##77.6% were male

#migratory: (predict less males)
sum(migratory$sex == "F",na.rm=T) ##96 female
sum(migratory$sex == "M",na.rm=T) ##59 male
binom.test(59,(96+59),0.5,"less")
binom.test(59,(96+59),0.5,"two.sided") ##CI from 30.2% to 46.1% ##report two sided teset (more conservative)
59/(96+59) ##38.1% male
96/(96+59)##68.9% female


###########
#####GLMS
###########
levels(data$genotype)
#data$genotype<- factor(data$genotype, levels=levels(data$genotype)[c(2,1,3)])
fitsex <- glm(lifehistory_binom~sex, data = data, family = "binomial")
summary(fitsex)
rsq(fitsex) ## 0.31

fitgeno <- glm(lifehistory_binom~genotype, data = data, family = "binomial")
summary(fitgeno)
rsq(fitgeno) # 0.20 (0.20)

fit1 <- glm(lifehistory_binom~genotype + sex, data = data, family = "binomial")
summary(fit1)
rsq(fit1) #0.45

##including interaction
fit2 <- glm(lifehistory_binom~genotype * sex, data = data,family="binomial")
summary(fit2)

fit3 <- glm(lifehistory_binom~genotype + sex + location,data=data, family = "binomial")
fit4 <- glm(lifehistory_binom~genotype * sex + location,data=data, family = "binomial")

BIC(fit1);BIC(fit2);BIC(fit3);BIC(fit4);BIC(fitsex);BIC(fitgeno)

###Plot coefficients for model with no interaction
coefficients <- as.data.frame(summary(fit1)$coefficients)
coefficients$parameter <- as.factor(rownames(coefficients))
coefficients<- droplevels(subset(coefficients, parameter != "(Intercept)"))
levels(coefficients$parameter)<- c("Migratory \\n Genotype", "Resident \\n Genotype", "Male")
coefficients$xmin <- coefficients$Estimate - coefficients$`Std. Error`
coefficients$xmax <- coefficients$Estimate + coefficients$`Std. Error`

##Figure 5b
fig5b<- ggplot(data=coefficients)+geom_point(aes(y = parameter, x = Estimate), size = 3)+
  geom_vline(xintercept = 0, color = "gray51", linetype=2)+
  geom_errorbarh(aes(x=Estimate, y = parameter, xmin = xmin, xmax = xmax,height=0.1))+
  theme_classic(11)+labs(y = "Parameter")+
  geom_text(x=-2.2,y=.6,label="B)",size=3)+
  geom_text(x=-2.0,y=3.2,label="-2.0",size=3)+
  geom_text(x=-1.3,y=2.2,label="-1.3",size=3)+
  geom_text(x=1.85,y=1.2,label="1.9",size=3)
fig5b

##plotting probabilities for model with no interaction
newdata <- data.frame(genotype = rep(c("Heterozygote","Migratory","Resident"),2), sex = c(rep("Male",3), rep("Female",3)))
pred <- data.frame(predict(fit1, newdata, type = "response", se.fit = T))
pred <- cbind(pred, newdata)
pred$prob_min <- pred$fit - pred$se.fit
pred$prob_max <- pred$fit + pred$se.fit

pred$genotype<- factor(pred$genotype, levels=levels(pred$genotype)[c(2,1,3)])
pred$genotype2<- pred$genotype
levels(pred$genotype2) <- c("Migratory\\nGenotype","Heterozygote\\nGenotype","Resident\\nGenotype")

##Figure 5a
fig5a <- ggplot(data=pred)+geom_errorbar(aes(x=genotype,ymin = prob_min, ymax = prob_max,group=sex),width=0.1)+
  geom_point(aes(x=genotype, y = fit, shape = sex,color=sex),size=3)+scale_color_manual(values=c("gray20","gray60"))+scale_shape_manual(values = c(15,17))+
  labs(x = "Genotype", y = "Probability of Out-Migrating")+theme_classic(11)+
  theme(legend.position = c(0.3,.15), legend.title=element_blank(),legend.background = element_blank(),
        legend.key.height = unit(.8,"line"))+
  geom_text(x=.6,y=0.25,label="A)",size=3)
fig5a

fig5 <- grid.arrange(fig5a,fig5b,nrow=1)

#ggsave("fig5.jpg", plot = fig5,
  #     scale = 1, width = 6, height = 2.5, units = c("in"),
   #    dpi = 350 )

#######FREQUNECY PLOTS##########
tables1 <- read.csv("tables1.csv")
juvfish_table <- subset(tables1, age == 0)
table_juv <- data.frame(Freq = c(sum(juvfish_table$migratory),sum(juvfish_table$heterozygote),sum(juvfish_table$resident)),
                         genotype = c("Migratory","Heterozygote","Resident"))


##copy table to make alterations for figure
table_juv2 <- table_juv
table_juv2$percent <- round(table_juv2$Freq/sum(table_juv2$Freq)*100,1)
table_juv2$percent<- table_juv2$percent/2
table_juv2$Freq<- table_juv2$Freq/2

##artificially doubling: assuming 50-50 for Juvenile fish
table_juv2<- rbind(table_juv2,table_juv2)
table_juv2$sex<- "Female"
table_juv2$sex[4:6]<- "Male"
table_juv2$sex<- as.factor(table_juv2$sex)
#table_juv$genotypesex <- as.factor(paste(table_juv$genotype,table_juv$sex, sep=" "))

table_juv2<- table_juv2[order(table_juv2$genotype),]
table_juv2$lifehistory <- "Juvenile Fish"

##make a Frequency plots from life history fish

data$genotype<- factor(data$genotype,levels=levels(data$genotype)[c(2,1,3)])


#migratory frequency table
table(migratory$genotype)/(dim(migratory)[1])*100
table_mig <- data.frame(table(migratory$genotype, migratory$sex))
colnames(table_mig)[1:2]<- c("genotype","sex")
#table_mig$genotypesex <- as.factor(paste(table_mig$genotype,table_mig$sex, sep=" "))
table_mig<- table_mig[order(table_mig$genotype),]
table_mig$percent <- round(table_mig$Freq/sum(table_mig$Freq)*100,1)

#residentfrequency table
table(resident$genotype)/(dim(resident)[1])*100
table_res <- data.frame(table(resident$genotype, resident$sex))
colnames(table_res)[1:2]<- c("genotype","sex")
#table_res$genotypesex <- as.factor(paste(table_res$genotype,table_res$sex, sep=" "))
table_res<- table_res[order(table_res$genotype),]
table_res$percent <- round(table_res$Freq/sum(table_res$Freq)*100,1)

table_mig$lifehistory <- "Migratory"
table_res$lifehistory <- "Resident"

table_all <- rbind(table_juv2, table_mig, table_res)
table_all$sexgenotype <- as.factor(paste(table_all$sex, table_all$genotype, sep = " "))
table_all$sexgenotype<- factor(table_all$sexgenotype, levels=levels(table_all$sexgenotype)[c(2,1,3,5,4,6)])

##re order levels so they are consistent for chi sq test
resident$genotype<- factor(resident$genotype, levels = levels(resident$genotype)[c(2,1,3)])
migratory$genotype<- factor(migratory$genotype, levels = levels(migratory$genotype)[c(2,1,3)])
##chi square tests ##using genotype frequencies and juvenile genotype frequencies as the null expectation
chisq.test(table(resident$genotype), p = table_juv$Freq/sum(table_juv$Freq))
chisq.test(table(migratory$genotype), p = table_juv$Freq/sum(table_juv$Freq))


###Bar plots for Figure 3###
###including juvenile fish
table_all$percent[table_all$Freq==0]<-NA
colors <- c("salmon1","tomato3","dark red","skyblue2","steelblue3","dark blue")
bp <- ggplot(data=table_all, aes(x = lifehistory,y=percent,fill=sexgenotype))+geom_bar(stat="identity",color="black")+
  scale_fill_manual(values = colors)+theme_classic(16)+theme(legend.title = element_blank())+
  labs(x="Ecotype",y="Percent")+
 geom_text(aes(label=paste(percent, "%",sep="")), position = position_stack(vjust = 0.5),color="white",fontface="bold")
bp             

##black and white graph
table_all$percent[table_all$Freq==0]<-NA
colors <- c("gray90","gray60","gray30","gray90","gray60","gray30")
bp <- ggplot(data=table_all, aes(x = lifehistory,y=percent,fill=sexgenotype))+geom_bar(stat="identity",color="black")+
  scale_fill_manual(values = colors)+theme_classic(11)+theme(legend.title = element_blank())+
  labs(x="Ecotype",y="Percent")+
  geom_text(aes(label=paste(round(percent,1), "%",sep="")), position = position_stack(vjust = 0.5),color="white",size=2.5)
bp

fig3 <- ggplot(data=table_all, aes(x = lifehistory,y=percent,fill=sexgenotype))+geom_bar(stat="identity",color="black")+
  scale_fill_manual(values = colors)+theme_classic(11)+theme(legend.title = element_blank())+
  labs(x="Ecotype",y="Percent")
fig3
#ggsave("fig3.jpg",plot=fig3,scale = 1, width = 4.5, height = 3.5, units = c("in"),
  #    dpi = 350)


###Changing genotype frequencies through age classes, and migratory vs non-migratory (Fig 4)

####
#year_summary_test <- ddply(tables1, .(year_capture),summarize,num_mig = sum(migratory),num_het = sum(heterozygote),num_res = sum(resident))
colnames(tables1)[3:5]<- c("num_mig","num_het","num_res")

tables1$num_fish = tables1$num_mig + tables1$num_res+tables1$num_het
tables1$num_mig_alleles = tables1$num_mig*2+tables1$num_het
tables1$num_res_alleles = tables1$num_res*2+tables1$num_het
tables1$prop_mig_alleles=tables1$num_mig_alleles/(tables1$num_fish*2)
tables1$prop_res_alleles=tables1$num_res_alleles/(tables1$num_fish*2)
tables1$year_age_0 <- tables1$year-tables1$age
tables1$Age <- as.factor(as.character(tables1$age))
levels(tables1$Age)<- c("Age-0", "Age-1", "Age-2")
tables1$year_age_0 <- as.factor(as.character(tables1$year_age_0))


bw2 <- c("black", "gray15","gray30","gray45","gray60","gray75")
p1 <- ggplot(data=tables1, aes(x=Age,y=prop_mig_alleles,group=year_age_0))+
  geom_point(aes(color=year_age_0,shape=year_age_0),size=3.4)+
  geom_line(aes(color=year_age_0),size=0.8)+
  labs(x="",y="Proportion Migratory Alleles",color="Year Age-0",shape="Year Age-0")+
  theme_classic(11)+
  scale_color_manual(values=bw2)+
  scale_shape_manual(values = c(15:19,15))+
  scale_x_discrete(expand=c(0.1,0.1))+
  theme(legend.position = c(0.8,0.7),legend.background = element_blank())+ylim(0,0.9)+
  geom_text(x=1,y=.02,label="A)",size=3)+
  theme(plot.margin = unit(c(1,0,1,1), "line"))
p1

##part 2
summary_migratory <- ddply(migratory,.(year),summarize,
                 num_mig = sum(genotype == "Migratory",na.rm=T),
                 num_het = sum(genotype == "Heterozygote",na.rm=T),
                 num_res = sum(genotype == "Resident",na.rm=T),
                 num_fish = num_mig + num_res+num_het,
                 num_mig_alleles = num_mig*2+num_het,
                 num_res_alleles = num_res*2+num_het,
                 prop_mig_alleles=num_mig_alleles/(num_fish*2),
                 prop_res_alleles=num_res_alleles/(num_fish*2))
summary_migratory$type <- as.factor("Migratory Fish")
tables1_byyear <- ddply(tables1,.(year_capture),summarize,
  num_mig = sum(num_mig),num_het = sum(num_het),num_res = sum(num_res),num_fish=sum(num_fish))
  
tables1_byyear$num_mig_alleles= tables1_byyear$num_mig*2+tables1_byyear$num_het
tables1_byyear$num_res_alleles=tables1_byyear$num_res*2+tables1_byyear$num_het
tables1_byyear$prop_mig_alleles = tables1_byyear$num_mig_alleles/(tables1_byyear$num_fish*2)
tables1_byyear$prop_res_alleles=tables1_byyear$num_res_alleles/(tables1_byyear$num_fish*2)
tables1_byyear$type <- as.factor("All Fish")
colnames(tables1_byyear)[1]<- "year"
colnames(summary_migratory)

year_summary_combo = rbind(tables1_byyear,summary_migratory)
year_summary_combo$year<- as.factor(as.character(year_summary_combo$year))

bw <- c("gray30","gray45","gray60","gray75")
p2 <- ggplot(data=year_summary_combo, aes(x=type,y=prop_mig_alleles,group=year))+
  geom_point(aes(color=year,shape=year),size=3.4)+
  geom_line(aes(color=year),size=0.8)+
  scale_color_manual(values=bw)+
  scale_shape_manual(values = c(17,18,19,15))+
  labs(x = "",y=" ",color="Year Captured",shape="Year Captured")+
  theme_classic(11)+
  scale_x_discrete(expand=c(0.2,0.2))+
  theme(legend.position = c(0.3,0.83), legend.background = element_blank())+ylim(0,0.9)+
  geom_text(x=1,y=0.02,label="B)",size=3)+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  theme(plot.margin = unit(c(1,1,1,0),"line"))
p2

fig4 <- grid.arrange(p1,p2,nrow=1)
#ggsave("fig4.jpg",plot=fig4,scale = 1, width = 6, height = 3.5, units = c("in"),
 #    dpi = 350)

