### Otter memory test ###
require("lme4")
require("MuMIn")
library("ggplot2")
#install.packages("Rmisc")
library("Rmisc")

###TASK SOLVING MEMORY###

##Read Data into R##
solvemem<-read.csv(file.choose(), fileEncoding="UTF-8-BOM")
names(solvemem)
summary(solvemem)
str(solvemem)

#Remove NAs##
solvemem<-subset(solvemem,!is.na(SolveI)) 
solvemem$SolveI

##Make R read discrete variables as discrete instead of continuous##
solvemem$Presentation<-as.factor(solvemem$Presentation) 
solvemem$Presentation
solvemem$Difficulty<-as.factor(solvemem$Difficulty)
solvemem$Difficulty
solvemem$Sex<-as.factor(solvemem$Sex)
solvemem$Sex

##Create Model##
smem <- glmer(SolveI~ solvemem$Presentation  + solvemem$Sex + solvemem$Difficulty + AgeCon + (1|OtterID),
                     data = solvemem, family = Gamma(link=log))

##Check Model assumptions##
summary(smem)
plot(resid(smem) ~ fitted(smem))
options(scipen = 999)
hist(resid(smem))
qqnorm(resid(smem))
qqline(resid(smem))


##Need this to set what to do if encounter an NA, will break otherwise##
options(na.action = "na.fail")

##Rank all possible combinations of variables##
AIC_Models <- dredge(smem, rank="AIC")
AIC_Models

##Subset by delta AIC cut-off (<2; can also be set to <6 to be more lenient)##
AIC_ModelsNEW <- subset(AIC_Models, delta<=2)
AIC_ModelsNEW 

##Apply nesting rule to eliminate more complex models from the top model set##
AIC_ModelsSUB <- subset(AIC_ModelsNEW, !nested(.))
AIC_ModelsSUB

##Run best model
smemBEST <- glmer(SolveI~ solvemem$Presentation  + solvemem$Sex + solvemem$Difficulty + (1|OtterID),
                     data = solvemem, family = Gamma(link=log)) #Gamma(link=log)as data positively skewed

##Check Assumptions again##
summary(smemBEST)
plot(resid(smemBEST) ~ fitted(smemBEST))
hist(resid(smemBEST))
qqnorm(resid(smemBEST))
qqline(resid(smemBEST))

##Extract fitted data and attach to database##
smemPRD<-cbind(solvemem, Fit=predict(smemBEST, newdata=solvemem, type='response'))

##Subset fitted values for variables of interest to report in manuscript##
smemRound1<-subset(smemPRD, Presentation=="1")
summary(smemRound1$Fit)
sd(smemRound1$Fit)
smemRound2<-subset(smemPRD, Presentation=="2")
summary(smemRound2$Fit)
sd(smemRound2$Fit)
smemTask1<-subset(smemPRD, Difficulty=="1")
summary(smemTask1$Fit)
sd(smemTask1$Fit)
smemTask2<-subset(smemPRD, Difficulty=="2")
summary(smemTask2$Fit)
sd(smemTask2$Fit)
smemTask3<-subset(smemPRD, Difficulty=="3")
summary(smemTask3$Fit)
sd(smemTask3$Fit)
smemTask4<-subset(smemPRD, Difficulty=="4")
summary(smemTask4$Fit)
sd(smemTask4$Fit)
smemTask5<-subset(smemPRD, Difficulty=="5")
summary(smemTask5$Fit)
sd(smemTask5$Fit)
smemMale<-subset(smemPRD, Sex=="0")
summary(smemMale$Fit)
sd(smemMale$Fit)
smemFemale<-subset(smemPRD, Sex=="1")
summary(smemFemale$Fit)
sd(smemFemale$Fit)

##Plot fitted values##
P1<-ggplot(smemPRD, aes(x = Presentation, y = Fit)) + 
  #geom_violin(width = 1)+
  geom_boxplot()+
  scale_y_continuous(limits=c(0,800),breaks=c(0,200,400,600,800))+
  theme_bw()+
        theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("Task Presentation Round") + ylab("Solve Latency (sec)") +
  ggtitle("(a)")+
        theme(plot.title = element_text(vjust = - 10))+
        theme(plot.title = element_text(hjust = "0.96"))+     
        theme(plot.title = element_text(color = "black", size = 13, face = "bold"))

P2<-ggplot(smemPRD, aes(x = Difficulty, y = Fit)) + 
  geom_boxplot()+
  scale_y_continuous(limits=c(0,800),breaks=c(0,200,400,600,800))+
   theme_bw()+
        theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("Task Type") +
  theme(axis.title.y=element_blank())+
  #theme(axis.text.y = element_text( color="White"))+
  ggtitle("(b)")+
        theme(plot.title = element_text(vjust = - 10))+
        theme(plot.title = element_text(hjust = "0.96"))+     
        theme(plot.title = element_text(color = "black", size = 13, face = "bold"))

P3<-ggplot(smemPRD, aes(x = Sex, y = Fit)) + 
  geom_boxplot()+
  scale_y_continuous(limits=c(0,800),breaks=c(0,200,400,600,800))+
   theme_bw()+
        theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_discrete(labels=c("0" = "Male", "1" = "Female"))+
  xlab("Sex") +
  theme(axis.title.y=element_blank())+
  #theme(axis.text.y = element_text( color="White"))+
  ggtitle("(c)")+
        theme(plot.title = element_text(vjust = - 10))+
        theme(plot.title = element_text(hjust = "0.96"))+     
        theme(plot.title = element_text(color = "black", size = 13, face = "bold"))

multiplot(P1,P2,P3,cols=3)