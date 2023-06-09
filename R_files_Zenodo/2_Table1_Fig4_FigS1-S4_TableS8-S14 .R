# 11 Oct 2021
# Kristof Kutasi
# Running time: 1 minute
packages <- c("dplyr","lubridate","xtable","tidyr","qdapTools","gghighlight","hexbin",
              "ggpubr","censReg","VGAM","stargazer","gridExtra","reshape2","svglite",
              "forcats","grid")
lapply(packages, require, character.only = TRUE)


dataCore=read.table("MASZK_dataCore.csv", sep=",", header = T)

# Table 1
# Sumary table "Table 1"
t1 = dataCore[!is.na(dataCore$VaccineType),] %>% select(c(Id, IsFemale, Schooling, ChronicIll, RejectedAny, Age, VaccineType,
                           PfizerRating,ModernaRating,AstraRating,SputRating,SinoRating))
t1$Category = ifelse(t1$VaccineType == "NoVaccine", 3,
                     ifelse(t1$RejectedAny == 0, 1, 2))
t1$RejectedAny = NULL
t1$AgeGroup = ifelse(t1$Age >= 20 & t1$Age < 40, 1, 
                     ifelse(t1$Age >= 40 & t1$Age < 60,2,
                            ifelse(t1$Age >= 60 & t1$Age < 80,3,
                                   ifelse(t1$Age >= 80,4,NA))))
t1$Education = ifelse(t1$Schooling == 1, 1,
                      ifelse(t1$Schooling == 2 | t1$Schooling == 3, 2,3))

t12 = t1[!is.na(t1$PfizerRating) & t1$PfizerRating == 1 & !is.na(t1$ModernaRating) & 
           t1$ModernaRating == 1 & !is.na(t1$AstraRating) & t1$AstraRating == 1 & 
           !is.na(t1$SputRating) & t1$SputRating == 1 & !is.na(t1$SinoRating) & 
           t1$SinoRating == 1,]

t13 = t1[(!is.na(t1$PfizerRating) & t1$PfizerRating == 1) | (!is.na(t1$ModernaRating) & 
                                                               t1$ModernaRating == 1) | (!is.na(t1$AstraRating) & t1$AstraRating == 1) | 
           (!is.na(t1$SputRating) & t1$SputRating == 1) | (!is.na(t1$SinoRating) & 
                                                             t1$SinoRating == 1),]

t1 = t1 %>% select(c(Id,IsFemale,Education,ChronicIll,AgeGroup,VaccineType,Category)) %>%
  group_by(Category)
t1 = t1 %>% summarise(Men = sum(IsFemale == 0), Women = sum(IsFemale == 1),
                      University = sum(Education ==3),"High-school" = sum(Education ==2), Elementary = sum(Education ==1),
                      "Not Chronic Illness" = sum(ChronicIll == 0), "Chronic Illness" = sum(ChronicIll == 1),
                      "Age 20-39" = sum(AgeGroup ==1, na.rm = T), "Age 40-59" = sum(AgeGroup ==2, na.rm = T), 
                      "Age 60-79" = sum(AgeGroup ==3, na.rm = T), "Age 80+" = sum(AgeGroup ==4, na.rm = T),
                      Pfizer = sum(VaccineType == "Pfizer"), Moderna = sum(VaccineType == "Moderna"), AstraZeneca = sum(VaccineType == "Astrazeneca"),
                      Sputnik = sum(VaccineType == "Sputnik"), Sinopharm = sum(VaccineType == "Sinopharm"))
t1 = as.data.frame(t(t1[2:17])); colnames(t1) = c("Vaccinated without rejection","Vaccinated after rejection", "Non-vaccinated")

t12 = t12 %>% summarise(Men = sum(IsFemale == 0), Women = sum(IsFemale == 1),
                        University = sum(Education ==3),"High-school" = sum(Education ==2), Elementary = sum(Education ==1),
                        "Not Chronic Illness" = sum(ChronicIll == 0), "Chronic Illness" = sum(ChronicIll == 1),
                        "Age 20-39" = sum(AgeGroup ==1, na.rm = T), "Age 40-59" = sum(AgeGroup ==2, na.rm = T), 
                        "Age 60-79" = sum(AgeGroup ==3, na.rm = T), "Age 80+" = sum(AgeGroup ==4, na.rm = T),
                        Pfizer = sum(VaccineType == "Pfizer"), Moderna = sum(VaccineType == "Moderna"), AstraZeneca = sum(VaccineType == "Astrazeneca"),
                        Sputnik = sum(VaccineType == "Sputnik"), Sinopharm = sum(VaccineType == "Sinopharm"))
t12 = as.data.frame(t(t12)); colnames(t12) = c("Rates all vaccines unacceptable")

t13 = t13 %>% summarise(Men = sum(IsFemale == 0), Women = sum(IsFemale == 1),
                        University = sum(Education ==3),"High-school" = sum(Education ==2), Elementary = sum(Education ==1),
                        "Not Chronic Illness" = sum(ChronicIll == 0), "Chronic Illness" = sum(ChronicIll == 1),
                        "Age 20-39" = sum(AgeGroup ==1, na.rm = T), "Age 40-59" = sum(AgeGroup ==2, na.rm = T), 
                        "Age 60-79" = sum(AgeGroup ==3, na.rm = T), "Age 80+" = sum(AgeGroup ==4, na.rm = T),
                        Pfizer = sum(VaccineType == "Pfizer"), Moderna = sum(VaccineType == "Moderna"), AstraZeneca = sum(VaccineType == "Astrazeneca"),
                        Sputnik = sum(VaccineType == "Sputnik"), Sinopharm = sum(VaccineType == "Sinopharm"))
t13 = as.data.frame(t(t13)); colnames(t13) = c("Rates one vaccine acceptable but another unacceptable")
Table1 = cbind(t12,t13,t1)
Table1 = round(Table1/nrow(dataCore)*100,1)
print(xtable(Table1,digits=1), type="latex", file="table1.tex")

# Figure 4
fig4Data = dataCore[dataCore$Vaccinated == 1,]
# Figure 4A
fig4NoReject = fig4Data[fig4Data$RejectedAny == 0,] %>% 
  filter(!is.na(fig4Data[fig4Data$RejectedAny == 0,]$VaccineType)) %>% 
  select(c(VaccineType,PfizerRating,ModernaRating,AstraRating,SputRating,SinoRating)) %>% 
  group_by(VaccineType) %>% summarise(
    Obs = n(),
    Pfizer = mean(PfizerRating, na.rm = TRUE),
    Moderna = mean(ModernaRating, na.rm = TRUE),
    Astrazeneca = mean(AstraRating, na.rm = TRUE),
    Sputnik = mean(SputRating, na.rm = TRUE),
    Sinopharm = mean(SinoRating, na.rm = TRUE)
  )
fig4NoReject[,3:7] = round(fig4NoReject[,3:7], digits = 2)
MeltedFig4NoReject <- melt(fig4NoReject %>% select(c(VaccineType,Pfizer,Moderna,Astrazeneca,Sputnik,Sinopharm)))
MatrixHeatMeanNoReject = MeltedFig4NoReject %>% mutate(VaccineType = fct_relevel(VaccineType,
                                                                                 "Pfizer","Moderna","Astrazeneca","Sputnik","Sinopharm"))  %>%
  ggplot(aes(x=VaccineType, y=variable, fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 3, limit = c(1,5), space = "Lab", 
                       name="Value") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 7, hjust = 1),
        axis.text.y=element_blank(),legend.position="none",
        plot.title = element_text(hjust = 0.5,size = 10))+
  coord_fixed() +
  ylab("Average vaccine rating") +
  xlab("Accepted vaccine") + 
  ggtitle("Rejected None")+ 
  geom_text(aes(VaccineType, variable, label = value), color = "black", size = 3)

# Figure 4B
fig4RejectedAny = fig4Data[fig4Data$RejectedAny == 1,] %>% 
  filter(!is.na(fig4Data[fig4Data$RejectedAny == 1,]$VaccineType)) %>% 
  select(c(VaccineType,PfizerRating,ModernaRating,AstraRating,SputRating,SinoRating)) %>% 
  group_by(VaccineType) %>% summarise(
    Obs = n(),
    Pfizer = mean(PfizerRating, na.rm = TRUE),
    Moderna = mean(ModernaRating, na.rm = TRUE),
    Astrazeneca = mean(AstraRating, na.rm = TRUE),
    Sputnik = mean(SputRating, na.rm = TRUE),
    Sinopharm = mean(SinoRating, na.rm = TRUE)
  )
fig4RejectedAny[,3:7] = round(fig4RejectedAny[,3:7], digits = 2)
MeltedFig4RejectedAny <- melt(fig4RejectedAny %>% select(c(VaccineType,Pfizer,Moderna,Astrazeneca,Sputnik,Sinopharm)))
MatrixHeatMeanRejectedAny = MeltedFig4RejectedAny %>% mutate(VaccineType = fct_relevel(VaccineType,
                                                                                       "Pfizer","Moderna","Astrazeneca","Sputnik","Sinopharm"))  %>%
  ggplot(aes(x=VaccineType, y=variable, fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 3, limit = c(1,5), space = "Lab", 
                       name="Value") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 7, hjust = 1),
        axis.text.y=element_blank(),legend.position="none",
        plot.title = element_text(hjust = 0.5,size = 10))+
  coord_fixed() +
  ylab("Average vaccine rating") +
  xlab("Accepted vaccine") + 
  ggtitle("Rejected Any")+ 
  geom_text(aes(VaccineType, variable, label = value), color = "black", size = 3)

# Figure 4C
fig4Central = fig4Data[fig4Data$RejectedAny == 1,] %>% 
  filter(!is.na(fig4Data[fig4Data$RejectedAny == 1,]$CentralVaccineAllocation)) %>% 
  select(c(CentralVaccineAllocation,PfizerRating,ModernaRating,AstraRating,SputRating,SinoRating)) %>% 
  group_by(CentralVaccineAllocation) %>% summarise(
    Obs = n(),
    Pfizer = mean(PfizerRating, na.rm = TRUE),
    Moderna = mean(ModernaRating, na.rm = TRUE),
    Astrazeneca = mean(AstraRating, na.rm = TRUE),
    Sputnik = mean(SputRating, na.rm = TRUE),
    Sinopharm = mean(SinoRating, na.rm = TRUE)
  )
fig4Central[,3:7] = round(fig4Central[,3:7], digits = 2)
MeltedFig4Central <- melt(fig4Central %>% select(c(CentralVaccineAllocation,Pfizer,Moderna,Astrazeneca,Sputnik,Sinopharm)))
MatrixHeatMeanCentral = MeltedFig4Central %>% mutate(CentralVaccineAllocation = fct_relevel(CentralVaccineAllocation,
                                                                                            "Pfizer","Moderna","Astrazeneca","Sputnik","Sinopharm"))  %>%
  ggplot(aes(x=CentralVaccineAllocation, y=variable, fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 3, limit = c(1,5), space = "Lab", 
                       name="Value") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 7, hjust = 1),
        axis.text.y=element_blank(),legend.position="none",
        plot.title = element_text(hjust = 0.5,size = 10))+
  coord_fixed() +
  ggtitle("Assigned but Rejected")+ 
  geom_text(aes(CentralVaccineAllocation, variable, label = value), color = "black", size = 3)

fig4 = ggarrange(MatrixHeatMeanNoReject + rremove("ylab") + rremove("xlab"), 
                 MatrixHeatMeanRejectedAny + rremove("ylab") + rremove("xlab"), 
                 MatrixHeatMeanCentral +rremove("ylab") + rremove("xlab"),
                 labels = c("A","B","C"),
                 ncol = 3, nrow = 1,
                 common.legend = TRUE, legend = "right",
                 align = "hv", 
                 font.label = list(size = 5, color = "black", face = "bold", family = NULL, position = "top"))
ggsave(file="Output/fig4.svg")

# Figure S1
# Vaccination date, age, vaccination type, chronic illness and rejection status distribution
figS1Data = dataCore %>% select(c(Id,VaccineType,Vaccine1Date,Age,Group)) %>% drop_na()
figS1 = ggplot(figS1Data, aes(x = Vaccine1Date, y = Age, col = Group, shape = VaccineType)) + 
  geom_point() +
  scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  xlab("Vaccination Date (1st dose)")
ggsave(file="Output/figS1.svg")

# Figure S2
# Vaccine type distribution conditionally on chronic illness and rejection status
figS2Data = dataCore %>% select(c(Id,VaccineType,Vaccine1Date,Age,Group)) %>% drop_na() %>% 
  group_by(VaccineType,Group) %>% summarise(Obs = n(), MeanAge = mean(Age), SdAge = sd(Age), 
                                            MeanV1Date = mean(Vaccine1Date), SdV1Date = sd(Vaccine1Date))
tmp = dataCore %>% select(c(Id,VaccineType,Vaccine1Date,Age,Group)) %>% drop_na() %>% 
  ungroup() %>% group_by(VaccineType) %>% summarise(VaccObs = n())
figS2Data = figS2Data %>% left_join(tmp, by="VaccineType"); tmp =NULL
figS2Data$WithinVaccRatio = figS2Data$Obs/figS2Data$VaccObs
figS2Data$MeanAge = round(figS2Data$MeanAge,2); figS2Data$SdAge = round(figS2Data$SdAge,2)
figS2Data$WeekV1 = as.numeric(strftime(figS2Data$MeanV1Date, format = "%V"))
figS2Data$SdWeekV1 = figS2Data$SdV1Date/7

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plotS2_1 = ggplot(figS2Data, aes(fill=Group, x=VaccineType, y=Obs)) + 
  geom_bar(position="dodge", stat = "identity") +
  theme(text = element_text(size=10),
        plot.title = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1,size = 10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        legend.position="bottom") + 
  guides(fill=guide_legend(title = "Groups: ",nrow=2,byrow=TRUE)) +
  ggtitle("") +
  ylab("Number of Observations") + 
  xlab("")
plotS2_2 = ggplot(figS2Data, aes(fill=Group, x=VaccineType, y=WithinVaccRatio)) + 
  geom_bar(stat = "identity", width = 0.7) +
  theme(text = element_text(size=10),
        plot.title = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1,size = 10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) +
  ylim(0,1) + 
  ggtitle("") +
  ylab("Group Ratio") +
  xlab("")
mylegend<-g_legend(plotS2_1)
figS2 = grid.arrange(arrangeGrob(plotS2_1 + theme(legend.position="none"), plotS2_2 + 
                                     theme(legend.position="none"), nrow = 1),mylegend, nrow=2,heights=c(5, 1),
                       top = textGrob("Within Vaccine Group Distribution",gp=gpar(fontsize=15)))
ggsave(file="Output/figS2.svg")

# Figure S3
# Average age distribution conditionally on vaccine type, chronic illness and rejection status
figS3 = ggplot(figS2Data, aes(fill=Group, x=VaccineType, y=MeanAge)) + 
  geom_bar(position="dodge", stat = "identity") +
  theme(text = element_text(size=10),
        plot.title = element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1,size = 10),        
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) + 
  ggtitle("Avearage Age Within Vaccines and Groups\\n") +
  ylab("Average Age") +
  xlab("")  + 
  geom_errorbar(
    aes(ymin = MeanAge - SdAge, 
        ymax = MeanAge + SdAge), 
    width=.2,
    position=position_dodge(.9),
    color = "black"
  )
ggsave(file="Output/figS3.svg")

# Figure S4
# Average vaccination week distribution in 2021 conditionally on vaccine type, chronic illness and rejection status.
figS4 = ggplot(figS2Data, aes(fill=Group, x=VaccineType, y=WeekV1)) + 
  geom_bar(position="dodge", stat = "identity") +
  theme(text = element_text(size=10),
        plot.title = element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1,size = 10),        
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) + 
  ggtitle("Average vaccination week in 2021 (1st Dose)\\n") +
  ylab("Week") +
  xlab("") + 
  geom_errorbar(
    aes(ymin = WeekV1 - SdWeekV1, 
        ymax = WeekV1 + SdWeekV1), 
    width=.2,
    position=position_dodge(.9),
    color = "black"
  )
ggsave(file="Output/figS4.svg")

# Table S8
fig4NoRejectSd = fig4Data[fig4Data$RejectedAny == 0,] %>% 
  filter(!is.na(fig4Data[fig4Data$RejectedAny == 0,]$VaccineType)) %>% 
  select(c(VaccineType,PfizerRating,ModernaRating,AstraRating,SputRating,SinoRating)) %>% 
  group_by(VaccineType) %>% summarise(
    Pfizer = sd(PfizerRating, na.rm = TRUE),
    Moderna = sd(ModernaRating, na.rm = TRUE),
    Astrazeneca = sd(AstraRating, na.rm = TRUE),
    Sputnik = sd(SputRating, na.rm = TRUE),
    Sinopharm = sd(SinoRating, na.rm = TRUE)
  )
fig4NoRejectSd[,2:6] = round(fig4NoRejectSd[,2:6], digits = 2)

fig4NoRejectObs = fig4Data[fig4Data$RejectedAny == 0,] %>% 
  filter(!is.na(fig4Data[fig4Data$RejectedAny == 0,]$VaccineType)) %>% 
  group_by(VaccineType) %>% summarise(
    Obs = n()
  )

MeltedFig4NoReject <- melt(fig4NoReject %>% select(-c(Obs))); colnames(MeltedFig4NoReject) = c("VaccineType","variable","Mean")
Meltedfig4NoRejectSd <- melt(fig4NoRejectSd); colnames(Meltedfig4NoRejectSd) = c("VaccineType","variable","SD")

MeltedFig4NoReject = MeltedFig4NoReject %>% left_join(Meltedfig4NoRejectSd, by=c("VaccineType","variable")) %>% 
  left_join(fig4NoRejectObs, by=c("VaccineType"))

colnames(MeltedFig4NoReject) = c("Accepted Vaccine","Rated Vaccine","Mean","Standard Deviation","Number of Observations")
MeltedFig4NoReject <- MeltedFig4NoReject[, c(1,2,5,3,4)]
MeltedFig4NoReject = MeltedFig4NoReject[order(MeltedFig4NoReject$`Accepted Vaccine`),]
print.xtable(xtable(MeltedFig4NoReject), file="Output/HeatMapTableNotRejected.tex")
write.table(MeltedFig4NoReject,"Output/HeatMapTableNotRejected.txt")

# Table S9
fig4RejectedAnySd = fig4Data[fig4Data$RejectedAny == 1,] %>% 
  filter(!is.na(fig4Data[fig4Data$RejectedAny == 1,]$VaccineType)) %>% 
  select(c(VaccineType,PfizerRating,ModernaRating,AstraRating,SputRating,SinoRating)) %>% 
  group_by(VaccineType) %>% summarise(
    Pfizer = sd(PfizerRating, na.rm = TRUE),
    Moderna = sd(ModernaRating, na.rm = TRUE),
    Astrazeneca = sd(AstraRating, na.rm = TRUE),
    Sputnik = sd(SputRating, na.rm = TRUE),
    Sinopharm = sd(SinoRating, na.rm = TRUE)
  )
fig4RejectedAnySd[,2:6] = round(fig4RejectedAnySd[,2:6], digits = 2)

fig4RejectedAnyObs = fig4Data[fig4Data$RejectedAny == 1,] %>% 
  filter(!is.na(fig4Data[fig4Data$RejectedAny == 1,]$VaccineType)) %>% 
  group_by(VaccineType) %>% summarise(
    Obs = n()
  )

MeltedFig4RejectedAny <- melt(fig4RejectedAny %>% select(-c(Obs))); colnames(MeltedFig4RejectedAny) = c("VaccineType","variable","Mean")
Meltedfig4RejectedAnySd <- melt(fig4RejectedAnySd); colnames(Meltedfig4RejectedAnySd) = c("VaccineType","variable","SD")

MeltedFig4RejectedAny = MeltedFig4RejectedAny %>% left_join(Meltedfig4RejectedAnySd, by=c("VaccineType","variable")) %>% 
  left_join(fig4RejectedAnyObs, by=c("VaccineType"))

colnames(MeltedFig4RejectedAny) = c("Accepted Vaccine","Rated Vaccine","Mean","Standard Deviation","Number of Observations")
MeltedFig4RejectedAny <- MeltedFig4RejectedAny[, c(1,2,5,3,4)]
MeltedFig4RejectedAny = MeltedFig4RejectedAny[order(MeltedFig4RejectedAny$`Accepted Vaccine`),]
print.xtable(xtable(MeltedFig4RejectedAny), file="Output/HeatMapTableRejectedAny.tex")
write.table(MeltedFig4RejectedAny,"Output/HeatMapTableRejectedAny.txt")

# Table S10
fig4CentralSd = fig4Data[fig4Data$RejectedAny == 1,] %>% 
  filter(!is.na(fig4Data[fig4Data$RejectedAny == 1,]$CentralVaccineAllocation)) %>% 
  select(c(CentralVaccineAllocation,PfizerRating,ModernaRating,AstraRating,SputRating,SinoRating)) %>% 
  group_by(CentralVaccineAllocation) %>% summarise(
    Pfizer = sd(PfizerRating, na.rm = TRUE),
    Moderna = sd(ModernaRating, na.rm = TRUE),
    Astrazeneca = sd(AstraRating, na.rm = TRUE),
    Sputnik = sd(SputRating, na.rm = TRUE),
    Sinopharm = sd(SinoRating, na.rm = TRUE)
  )
fig4CentralSd[is.na(fig4CentralSd)] <- 0
fig4CentralSd[,2:6] = round(fig4CentralSd[,2:6], digits = 2)

fig4CentralObs = fig4Data[fig4Data$RejectedAny == 1,] %>% 
  filter(!is.na(fig4Data[fig4Data$RejectedAny == 1,]$CentralVaccineAllocation)) %>% 
  group_by(CentralVaccineAllocation) %>% summarise(
    Obs = n()
  )

MeltedFig4Central <- melt(fig4Central %>% select(-c(Obs))); colnames(MeltedFig4Central) = c("CentralVaccineAllocation","variable","Mean")
MeltedFig4CentralSd <- melt(fig4CentralSd); colnames(MeltedFig4CentralSd) = c("CentralVaccineAllocation","variable","SD")

MeltedFig4Central = MeltedFig4Central %>% left_join(MeltedFig4CentralSd, by=c("CentralVaccineAllocation","variable")) %>% 
  left_join(fig4CentralObs, by=c("CentralVaccineAllocation"))
colnames(MeltedFig4Central) = c("Assigned Vaccine","Rated Vaccine","Mean","Standard Deviation","Number of Observations")
MeltedFig4Central <- MeltedFig4Central[, c(1,2,5,3,4)]
MeltedFig4Central = MeltedFig4Central[order(MeltedFig4Central$`Assigned Vaccine`),]
print.xtable(xtable(MeltedFig4Central), file="Output/HeatMapTableCentral.tex")
write.table(MeltedFig4Central,"Output/HeatMapTableCentral.txt")

# Table S11
# Change attitudes towards vaccines
t11 = table(dataCore$PfizerChange)
t11 = rbind(t11,table(dataCore$ModernaChange))
t11 = rbind(t11,table(dataCore$AstraZenecaChange))
t11 = rbind(t11,table(dataCore$SputnikChange))
t11 = rbind(t11,table(dataCore$SinopharmChange))
t11 = cbind(t11[,1],t11[,3])
rownames(t11) = c("Pfizer","Moderna","AstraZeneca","Sputnik","Sinopharm")
colnames(t11) = c("Negative change","Positive change")
print(xtable(t11), type="latex", file="table11.tex")

# Table S12
dataWait = dataCore %>% select(c(Id,Vaccinated,VaccineType,Vaccine1Date,CentralVaccineAllocation,
                                 EarliestRejectionDate,Wait,PfizerRating,ModernaRating,AstraRating,
                                 SputRating,SinoRating,RejectedAny,RejectedPfizer,RejectedPfizer,
                                 RejectedModerna,RejectedAstra,RejectedSputnik,RejectedSinoph))
dataRating = dataCore[dataCore$Vaccinated == 1,] %>% 
  select(c(PfizerRating,ModernaRating,AstraRating,SputRating,SinoRating,
           JanssenRating,VaccineType,Id))
dataRating[1:6][is.na(dataRating[1:6])] = 0
dataRating = dataRating %>% mutate(MaxRating = pmax(PfizerRating,ModernaRating,AstraRating,SputRating,SinoRating,JanssenRating))

dataRating[is.na(dataRating$VaccineType),]$VaccineType = "6"
dataRating[dataRating$VaccineType == "Pfizer",]$VaccineType = "1"
dataRating[dataRating$VaccineType == "Moderna",]$VaccineType = "2"
dataRating[dataRating$VaccineType == "Astrazeneca",]$VaccineType = "3"
dataRating[dataRating$VaccineType == "Sputnik",]$VaccineType = "4"
dataRating[dataRating$VaccineType == "Sinopharm",]$VaccineType = "5"
dataRating$VaccineType = as.numeric(dataRating$VaccineType)

for (i in 1:739) {
  dataRating[i,10] = dataRating[i,dataRating[i,]$VaccineType]
}
names(dataRating)[10] = "ObtainedRating"
dataRating$GotHighestRating = ifelse(dataRating$MaxRating == dataRating$ObtainedRating,1,0)
dataNoBest = dataRating[dataRating$GotHighestRating == 0,]
for (i in 1:nrow(dataNoBest)) {
  for (j in 1:6) {
    dataNoBest[i,j+11] = ifelse(dataNoBest[i,9] == dataNoBest[i,j],1,0)
  }
}
names(dataNoBest)[12:17] = c("PfizerBest","ModernaBest","AstraBest","SputnikBest","SinoBest","JanssenBest")
table12 = dataNoBest[dataNoBest$VaccineType != "6",] %>% 
  group_by(VaccineType) %>% summarise(obsNum = n(),
                                      meanPfizerBest = mean(PfizerBest),
                                      meanModernaBest = mean(ModernaBest),
                                      meanAstraBest = mean(AstraBest),
                                      meanSputnikBest = mean(SputnikBest),
                                      meanSinoBest = mean(SinoBest))
table12[,3:7] = round(table12[,3:7],2)
table12$VaccineType = as.character(table12$VaccineType)
table12[table12$VaccineType == "1",]$VaccineType = "Pfizer"
table12[table12$VaccineType == "2",]$VaccineType = "Moderna"
table12[table12$VaccineType == "3",]$VaccineType = "Astrazeneca"
table12[table12$VaccineType == "4",]$VaccineType = "Sputnik"
table12[table12$VaccineType == "5",]$VaccineType = "Sinopharm"
colnames(table12) = c("Accepted vaccine", "Number of Observations", "Pfizer Best",
                      "Moderna Best", "Astrazeneca Best", "Sputnik Best","Sinopharm Best")
print(xtable(table12), type="latex", file="table12.tex")

# Table S13
dataRating = dataRating %>% 
  left_join(dataNoBest %>% select(c("PfizerBest","ModernaBest","AstraBest","SputnikBest","SinoBest","JanssenBest","Id")), by = "Id")
dataWait = dataWait %>% left_join(dataRating %>% select(c("Id","MaxRating","ObtainedRating","GotHighestRating",
                                                          "PfizerBest","ModernaBest","AstraBest","SputnikBest","SinoBest","JanssenBest")), by = "Id")
dataWait[is.na(dataWait$GotHighestRating),]$GotHighestRating = -1
dataWait$Diff_in_weeks = round(difftime(dataWait$Vaccine1Date, dataWait$EarliestRejectionDate, units = "weeks"),0) # weeks
dataRejectedBestWait = dataWait[dataWait$VaccineType != "NoVaccine" & !is.na(dataWait$VaccineType),] %>% 
  select(c(Id,VaccineType,RejectedAny,GotHighestRating,CentralVaccineAllocation,
           Vaccine1Date,EarliestRejectionDate,Wait,Diff_in_weeks,ObtainedRating,
           MaxRating,RejectedPfizer,RejectedPfizer,RejectedModerna,RejectedAstra,
           RejectedSputnik,RejectedSinoph))
table13 = dataRejectedBestWait %>% ungroup() %>%
  group_by(RejectedAny, GotHighestRating) %>% 
  summarise(obsNum = n(),
            meanExpectedWait = mean(Wait, na.rm = TRUE),
            meanActualWait = mean(Diff_in_weeks, na.rm = TRUE),
            sdExpectedWait = sd(Wait, na.rm = TRUE),
            sdActualWait = sd(Diff_in_weeks, na.rm = TRUE))  %>%
  mutate(lowerCIExpected = meanExpectedWait - qt(1 - (0.1 / 2), obsNum - 1) * sdExpectedWait/sqrt(obsNum),
         upperCIExpected = meanExpectedWait + qt(1 - (0.1 / 2), obsNum - 1) * sdExpectedWait/sqrt(obsNum),
         lowerCIActual = meanActualWait - qt(1 - (0.1 / 2), obsNum - 1) * sdActualWait/sqrt(obsNum),
         upperCIActual = meanActualWait + qt(1 - (0.1 / 2), obsNum - 1) * sdActualWait/sqrt(obsNum))
table13[sapply(table13, is.numeric)] <- lapply(table13[sapply(table13, is.numeric)], round, 2)
print(xtable(table13), type="latex", file="table13.tex")

# Table S14
dataNoBestSummary = dataRejectedBestWait[dataRejectedBestWait$GotHighestRating == 0,] %>% group_by(VaccineType) %>% 
  filter(if_any(everything(), ~ !is.na(.))) %>% 
  summarise(obsNumAN = n(),
            meanExpected = mean(Wait, na.rm = TRUE),
            sdExpectedAN = sd(Wait, na.rm = TRUE),
            lowerCIExpected =  meanExpected - qt(1 - (0.1 / 2), obsNumAN - 1) * sdExpectedAN/sqrt(obsNumAN),
            upperCIExpected =  meanExpected + qt(1 - (0.1 / 2), obsNumAN - 1) * sdExpectedAN/sqrt(obsNumAN))
dataBestRejectPf = dataRejectedBestWait[dataRejectedBestWait$GotHighestRating == 1 & dataRejectedBestWait$RejectedAny == 1 & dataRejectedBestWait$CentralVaccineAllocation == "Pfizer",]  %>% 
  filter(if_any(everything(), ~ !is.na(.))) %>% 
  summarise(VaccineType = "Pfizer",
            obsNumRB = n(),
            meanActual = mean(Diff_in_weeks, na.rm = TRUE),
            sdActual = sd(Diff_in_weeks, na.rm = TRUE),
            lowerCIActual =  meanActual - qt(1 - (0.1 / 2), obsNumRB - 1) * sdActual/sqrt(obsNumRB),
            upperCIActual =  meanActual + qt(1 - (0.1 / 2), obsNumRB - 1) * sdActual/sqrt(obsNumRB))
dataBestRejectMo = dataRejectedBestWait[dataRejectedBestWait$GotHighestRating == 1 & dataRejectedBestWait$RejectedAny == 1 & dataRejectedBestWait$CentralVaccineAllocation == "Moderna",]  %>% 
  filter(if_any(everything(), ~ !is.na(.))) %>% 
  summarise(VaccineType = "Moderna",
            obsNumRB = n(),
            meanActual = mean(Diff_in_weeks, na.rm = TRUE),
            sdActual = sd(Diff_in_weeks, na.rm = TRUE),
            lowerCIActual =  meanActual - qt(1 - (0.1 / 2), obsNumRB - 1) * sdActual/sqrt(obsNumRB),
            upperCIActual =  meanActual + qt(1 - (0.1 / 2), obsNumRB - 1) * sdActual/sqrt(obsNumRB))
dataBestRejectAs = dataRejectedBestWait[dataRejectedBestWait$GotHighestRating == 1 & dataRejectedBestWait$RejectedAny == 1 & dataRejectedBestWait$CentralVaccineAllocation == "Astrazeneca",]  %>% 
  filter(if_any(everything(), ~ !is.na(.))) %>% 
  summarise(VaccineType = "Astrazeneca",
            obsNumRB = n(),
            meanActual = mean(Diff_in_weeks, na.rm = TRUE),
            sdActual = sd(Diff_in_weeks, na.rm = TRUE),
            lowerCIActual =  meanActual - qt(1 - (0.1 / 2), obsNumRB - 1) * sdActual/sqrt(obsNumRB),
            upperCIActual =  meanActual + qt(1 - (0.1 / 2), obsNumRB - 1) * sdActual/sqrt(obsNumRB))
dataBestRejectSp = dataRejectedBestWait[dataRejectedBestWait$GotHighestRating == 1 & dataRejectedBestWait$RejectedAny == 1 & dataRejectedBestWait$CentralVaccineAllocation == "Sputnik",]  %>% 
  filter(if_any(everything(), ~ !is.na(.))) %>% 
  summarise(VaccineType = "Sputnik",
            obsNumRB = n(),
            meanActual = mean(Diff_in_weeks, na.rm = TRUE),
            sdActual = sd(Diff_in_weeks, na.rm = TRUE),
            lowerCIActual =  meanActual - qt(1 - (0.1 / 2), obsNumRB - 1) * sdActual/sqrt(obsNumRB),
            upperCIActual =  meanActual + qt(1 - (0.1 / 2), obsNumRB - 1) * sdActual/sqrt(obsNumRB))
dataBestRejectSi = dataRejectedBestWait[dataRejectedBestWait$GotHighestRating == 1 & dataRejectedBestWait$RejectedAny == 1 & dataRejectedBestWait$CentralVaccineAllocation == "Sinopharm",]  %>% 
  filter(if_any(everything(), ~ !is.na(.))) %>% 
  summarise(VaccineType = "Sinopharm",
            obsNumRB = n(),
            meanActual = mean(Diff_in_weeks, na.rm = TRUE),
            sdActual = sd(Diff_in_weeks, na.rm = TRUE),
            lowerCIActual =  meanActual - qt(1 - (0.1 / 2), obsNumRB - 1) * sdActual/sqrt(obsNumRB),
            upperCIActual =  meanActual + qt(1 - (0.1 / 2), obsNumRB - 1) * sdActual/sqrt(obsNumRB))
dataBestReject = rbind(dataBestRejectPf,dataBestRejectMo,dataBestRejectAs,dataBestRejectSp,dataBestRejectSi)

table14 = dataNoBestSummary %>% left_join(dataBestReject, by = "VaccineType")
table14[sapply(table14, is.numeric)] <- lapply(table14[sapply(table14, is.numeric)], round, 2)
table14 <- table14[, c("VaccineType", "obsNumAN", "obsNumRB","meanExpected","meanActual",
                                 "lowerCIExpected","upperCIExpected","lowerCIActual","upperCIActual")]
print(xtable(table14), type="latex", file="table14.tex")



