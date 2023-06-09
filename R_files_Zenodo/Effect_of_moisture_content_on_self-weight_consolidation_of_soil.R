## Importing data

pc15_replicate_1 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 15pc replicate 1.csv")
pc15_replicate_2 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 15pc replicate 2.csv")
pc15_replicate_3 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 15pc replicate 3.csv")
pc15_replicate_4 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 15pc replicate 4.csv")
pc15_replicate_5 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 15pc replicate 5.csv")
pc15_replicate_6 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 15pc replicate 6.csv")

pc20_replicate_1 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 20pc replicate 1.csv")
pc20_replicate_2 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 20pc replicate 2.csv")
pc20_replicate_3 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 20pc replicate 3.csv")
pc20_replicate_4 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 20pc replicate 4.csv")
pc20_replicate_5 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 20pc replicate 5.csv")
pc20_replicate_6 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 20pc replicate 6.csv")

pc30_replicate_1 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 30pc replicate 1.csv")
pc30_replicate_2 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 30pc replicate 2.csv")
pc30_replicate_3 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 30pc replicate 3.csv")
pc30_replicate_4 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 30pc replicate 4.csv")
pc30_replicate_5 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 30pc replicate 5.csv")

pc40_replicate_1 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 40pc replicate 1.csv")
pc40_replicate_2 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 40pc replicate 2.csv")
pc40_replicate_3 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 40pc replicate 3.csv")
pc40_replicate_4 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 40pc replicate 4.csv")
pc40_replicate_5 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 40pc replicate 5.csv")

pc50_replicate_1 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 50pc replicate 1.csv")
pc50_replicate_2 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 50pc replicate 2.csv")
pc50_replicate_3 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 50pc replicate 3.csv")
pc50_replicate_4 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 50pc replicate 4.csv")
pc50_replicate_5 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 50pc replicate 5.csv")
pc50_replicate_6 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 50pc replicate 6.csv")

pc60_replicate_1 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 60pc replicate 1.csv")
pc60_replicate_2 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 60pc replicate 2.csv")
pc60_replicate_3 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 60pc replicate 3.csv")
pc60_replicate_4 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Moisture 60pc replicate 4.csv")

# Step1: Eliminating unnecessary columns and rows and reassigning row numbers
# Step2: Creating new columns representing size of samples and binding all data sets into a single data set
# Step 3: Declaring stress and strain as numeric and character vectors
# Step 4: Creating new columns representing size of samples and binding all data sets into a single data set


pc15_replicate_1 <- pc15_replicate_1[3:nrow(pc15_replicate_1),1:2]
row.names(pc15_replicate_1)<- 1:nrow(pc15_replicate_1)

a<- rep("15% moisture replicate 1",nrow(pc15_replicate_1))
a=data.frame(a)
pc15_replicate_1 = cbind(pc15_replicate_1, a)
colnames(pc15_replicate_1)=c("stress","strain","location")# naming column names
rm(a)

pc15_replicate_1$stress<- as.numeric(as.character(pc15_replicate_1$stress))
pc15_replicate_1$strain<- as.numeric(as.character(pc15_replicate_1$strain))

##########

pc15_replicate_2 <- pc15_replicate_2[3:nrow(pc15_replicate_2),1:2]
row.names(pc15_replicate_2)<- 1:nrow(pc15_replicate_2)

a<- rep("15% moisture replicate 2",nrow(pc15_replicate_2))
a=data.frame(a)
pc15_replicate_2 = cbind(pc15_replicate_2, a)
colnames(pc15_replicate_2)=c("stress","strain","location")# naming column names
rm(a)

pc15_replicate_2$stress<- as.numeric(as.character(pc15_replicate_2$stress))
pc15_replicate_2$strain<- as.numeric(as.character(pc15_replicate_2$strain))

##########


pc15_replicate_3 <- pc15_replicate_3[3:nrow(pc15_replicate_3),1:2]
row.names(pc15_replicate_3)<- 1:nrow(pc15_replicate_3)

a<- rep("15% moisture replicate 3",nrow(pc15_replicate_3))
a=data.frame(a)
pc15_replicate_3 = cbind(pc15_replicate_3, a)
colnames(pc15_replicate_3)=c("stress","strain","location")# naming column names
rm(a)

pc15_replicate_3$stress<- as.numeric(as.character(pc15_replicate_3$stress))
pc15_replicate_3$strain<- as.numeric(as.character(pc15_replicate_3$strain))

##########

pc15_replicate_4 <- pc15_replicate_4[3:nrow(pc15_replicate_4),1:2]
row.names(pc15_replicate_4)<- 1:nrow(pc15_replicate_4)

a<- rep("15% moisture replicate 4",nrow(pc15_replicate_4))
a=data.frame(a)
pc15_replicate_4 = cbind(pc15_replicate_4, a)
colnames(pc15_replicate_4)=c("stress","strain","location")# naming column names
rm(a)

pc15_replicate_4$stress<- as.numeric(as.character(pc15_replicate_4$stress))
pc15_replicate_4$strain<- as.numeric(as.character(pc15_replicate_4$strain))

##########

pc15_replicate_5 <- pc15_replicate_5[3:nrow(pc15_replicate_5),1:2]
row.names(pc15_replicate_5)<- 1:nrow(pc15_replicate_5)

a<- rep("15% moisture replicate 5",nrow(pc15_replicate_5))
a=data.frame(a)
pc15_replicate_5 = cbind(pc15_replicate_5, a)
colnames(pc15_replicate_5)=c("stress","strain","location")# naming column names
rm(a)

pc15_replicate_5$stress<- as.numeric(as.character(pc15_replicate_5$stress))
pc15_replicate_5$strain<- as.numeric(as.character(pc15_replicate_5$strain))

##########

pc15_replicate_6 <- pc15_replicate_6[3:nrow(pc15_replicate_6),1:2]
row.names(pc15_replicate_6)<- 1:nrow(pc15_replicate_6)

a<- rep("15% moisture replicate 6",nrow(pc15_replicate_6))
a=data.frame(a)
pc15_replicate_6 = cbind(pc15_replicate_6, a)
colnames(pc15_replicate_6)=c("stress","strain","location")# naming column names
rm(a)

pc15_replicate_6$stress<- as.numeric(as.character(pc15_replicate_6$stress))
pc15_replicate_6$strain<- as.numeric(as.character(pc15_replicate_6$strain))

##########

pc20_replicate_1 <- pc20_replicate_1[3:nrow(pc20_replicate_1),1:2]
row.names(pc20_replicate_1)<- 1:nrow(pc20_replicate_1)

a<- rep("20% moisture replicate 1",nrow(pc20_replicate_1))
a=data.frame(a)
pc20_replicate_1 = cbind(pc20_replicate_1, a)
colnames(pc20_replicate_1)=c("stress","strain","location")# naming column names
rm(a)

pc20_replicate_1$stress<- as.numeric(as.character(pc20_replicate_1$stress))
pc20_replicate_1$strain<- as.numeric(as.character(pc20_replicate_1$strain))

##########

pc20_replicate_2 <- pc20_replicate_2[3:nrow(pc20_replicate_2),1:2]
row.names(pc20_replicate_2)<- 1:nrow(pc20_replicate_2)

a<- rep("20% moisture replicate 2",nrow(pc20_replicate_2))
a=data.frame(a)
pc20_replicate_2 = cbind(pc20_replicate_2, a)
colnames(pc20_replicate_2)=c("stress","strain","location")# naming column names
rm(a)

pc20_replicate_2$stress<- as.numeric(as.character(pc20_replicate_2$stress))
pc20_replicate_2$strain<- as.numeric(as.character(pc20_replicate_2$strain))

##########

pc20_replicate_3 <- pc20_replicate_3[3:nrow(pc20_replicate_3),1:2]
row.names(pc20_replicate_3)<- 1:nrow(pc20_replicate_3)

a<- rep("20% moisture replicate 3",nrow(pc20_replicate_3))
a=data.frame(a)
pc20_replicate_3 = cbind(pc20_replicate_3, a)
colnames(pc20_replicate_3)=c("stress","strain","location")# naming column names
rm(a)

pc20_replicate_3$stress<- as.numeric(as.character(pc20_replicate_3$stress))
pc20_replicate_3$strain<- as.numeric(as.character(pc20_replicate_3$strain))

##########

pc20_replicate_4 <- pc20_replicate_4[3:nrow(pc20_replicate_4),1:2]
row.names(pc20_replicate_4)<- 1:nrow(pc20_replicate_4)

a<- rep("20% moisture replicate 4",nrow(pc20_replicate_4))
a=data.frame(a)
pc20_replicate_4 = cbind(pc20_replicate_4, a)
colnames(pc20_replicate_4)=c("stress","strain","location")# naming column names
rm(a)

pc20_replicate_4$stress<- as.numeric(as.character(pc20_replicate_4$stress))
pc20_replicate_4$strain<- as.numeric(as.character(pc20_replicate_4$strain))

##########

pc20_replicate_5 <- pc20_replicate_5[3:nrow(pc20_replicate_5),1:2]
row.names(pc20_replicate_5)<- 1:nrow(pc20_replicate_5)

a<- rep("20% moisture replicate 5",nrow(pc20_replicate_5))
a=data.frame(a)
pc20_replicate_5 = cbind(pc20_replicate_5, a)
colnames(pc20_replicate_5)=c("stress","strain","location")# naming column names
rm(a)

pc20_replicate_5$stress<- as.numeric(as.character(pc20_replicate_5$stress))
pc20_replicate_5$strain<- as.numeric(as.character(pc20_replicate_5$strain))

##########

pc20_replicate_6 <- pc20_replicate_6[3:nrow(pc20_replicate_6),1:2]
row.names(pc20_replicate_6)<- 1:nrow(pc20_replicate_6)

a<- rep("20% moisture replicate 6",nrow(pc20_replicate_6))
a=data.frame(a)
pc20_replicate_6 = cbind(pc20_replicate_6, a)
colnames(pc20_replicate_6)=c("stress","strain","location")# naming column names
rm(a)

pc20_replicate_6$stress<- as.numeric(as.character(pc20_replicate_6$stress))
pc20_replicate_6$strain<- as.numeric(as.character(pc20_replicate_6$strain))

##########

pc30_replicate_1 <- pc30_replicate_1[3:nrow(pc30_replicate_1),1:2]
row.names(pc30_replicate_1)<- 1:nrow(pc30_replicate_1)

a<- rep("30% moisture replicate 1",nrow(pc30_replicate_1))
a=data.frame(a)
pc30_replicate_1 = cbind(pc30_replicate_1, a)
colnames(pc30_replicate_1)=c("stress","strain","location")# naming column names
rm(a)

pc30_replicate_1$stress<- as.numeric(as.character(pc30_replicate_1$stress))
pc30_replicate_1$strain<- as.numeric(as.character(pc30_replicate_1$strain))

##########

pc30_replicate_2 <- pc30_replicate_2[3:nrow(pc30_replicate_2),1:2]
row.names(pc30_replicate_2)<- 1:nrow(pc30_replicate_2)

a<- rep("30% moisture replicate 2",nrow(pc30_replicate_2))
a=data.frame(a)
pc30_replicate_2 = cbind(pc30_replicate_2, a)
colnames(pc30_replicate_2)=c("stress","strain","location")# naming column names
rm(a)

pc30_replicate_2$stress<- as.numeric(as.character(pc30_replicate_2$stress))
pc30_replicate_2$strain<- as.numeric(as.character(pc30_replicate_2$strain))

##########

pc30_replicate_3 <- pc30_replicate_3[3:nrow(pc30_replicate_3),1:2]
row.names(pc30_replicate_3)<- 1:nrow(pc30_replicate_3)

a<- rep("30% moisture replicate 3",nrow(pc30_replicate_3))
a=data.frame(a)
pc30_replicate_3 = cbind(pc30_replicate_3, a)
colnames(pc30_replicate_3)=c("stress","strain","location")# naming column names
rm(a)

pc30_replicate_3$stress<- as.numeric(as.character(pc30_replicate_3$stress))
pc30_replicate_3$strain<- as.numeric(as.character(pc30_replicate_3$strain))

##########

pc30_replicate_4 <- pc30_replicate_4[3:nrow(pc30_replicate_4),1:2]
row.names(pc30_replicate_4)<- 1:nrow(pc30_replicate_4)

a<- rep("30% moisture replicate 4",nrow(pc30_replicate_4))
a=data.frame(a)
pc30_replicate_4 = cbind(pc30_replicate_4, a)
colnames(pc30_replicate_4)=c("stress","strain","location")# naming column names
rm(a)

pc30_replicate_4$stress<- as.numeric(as.character(pc30_replicate_4$stress))
pc30_replicate_4$strain<- as.numeric(as.character(pc30_replicate_4$strain))

##########

pc30_replicate_5 <- pc30_replicate_5[3:nrow(pc30_replicate_5),1:2]
row.names(pc30_replicate_5)<- 1:nrow(pc30_replicate_5)

a<- rep("30% moisture replicate 5",nrow(pc30_replicate_5))
a=data.frame(a)
pc30_replicate_5 = cbind(pc30_replicate_5, a)
colnames(pc30_replicate_5)=c("stress","strain","location")# naming column names
rm(a)

pc30_replicate_5$stress<- as.numeric(as.character(pc30_replicate_5$stress))
pc30_replicate_5$strain<- as.numeric(as.character(pc30_replicate_5$strain))

##########

pc40_replicate_1 <- pc40_replicate_1[3:nrow(pc40_replicate_1),1:2]
row.names(pc40_replicate_1)<- 1:nrow(pc40_replicate_1)

a<- rep("40% moisture replicate 1",nrow(pc40_replicate_1))
a=data.frame(a)
pc40_replicate_1 = cbind(pc40_replicate_1, a)
colnames(pc40_replicate_1)=c("stress","strain","location")# naming column names
rm(a)

pc40_replicate_1$stress<- as.numeric(as.character(pc40_replicate_1$stress))
pc40_replicate_1$strain<- as.numeric(as.character(pc40_replicate_1$strain))

##########

pc40_replicate_2 <- pc40_replicate_2[3:nrow(pc40_replicate_2),1:2]
row.names(pc40_replicate_2)<- 1:nrow(pc40_replicate_2)

a<- rep("40% moisture replicate 2",nrow(pc40_replicate_2))
a=data.frame(a)
pc40_replicate_2 = cbind(pc40_replicate_2, a)
colnames(pc40_replicate_2)=c("stress","strain","location")# naming column names
rm(a)

pc40_replicate_2$stress<- as.numeric(as.character(pc40_replicate_2$stress))
pc40_replicate_2$strain<- as.numeric(as.character(pc40_replicate_2$strain))

##########

pc40_replicate_3 <- pc40_replicate_3[3:nrow(pc40_replicate_3),1:2]
row.names(pc40_replicate_3)<- 1:nrow(pc40_replicate_3)

a<- rep("40% moisture replicate 3",nrow(pc40_replicate_3))
a=data.frame(a)
pc40_replicate_3 = cbind(pc40_replicate_3, a)
colnames(pc40_replicate_3)=c("stress","strain","location")# naming column names
rm(a)

pc40_replicate_3$stress<- as.numeric(as.character(pc40_replicate_3$stress))
pc40_replicate_3$strain<- as.numeric(as.character(pc40_replicate_3$strain))

##########

pc40_replicate_4 <- pc40_replicate_4[3:nrow(pc40_replicate_4),1:2]
row.names(pc40_replicate_4)<- 1:nrow(pc40_replicate_4)

a<- rep("40% moisture replicate 4",nrow(pc40_replicate_4))
a=data.frame(a)
pc40_replicate_4 = cbind(pc40_replicate_4, a)
colnames(pc40_replicate_4)=c("stress","strain","location")# naming column names
rm(a)

pc40_replicate_4$stress<- as.numeric(as.character(pc40_replicate_4$stress))
pc40_replicate_4$strain<- as.numeric(as.character(pc40_replicate_4$strain))

##########

pc40_replicate_5 <- pc40_replicate_5[3:nrow(pc40_replicate_5),1:2]
row.names(pc40_replicate_5)<- 1:nrow(pc40_replicate_5)

a<- rep("40% moisture replicate 5",nrow(pc40_replicate_5))
a=data.frame(a)
pc40_replicate_5 = cbind(pc40_replicate_5, a)
colnames(pc40_replicate_5)=c("stress","strain","location")# naming column names
rm(a)

pc40_replicate_5$stress<- as.numeric(as.character(pc40_replicate_5$stress))
pc40_replicate_5$strain<- as.numeric(as.character(pc40_replicate_5$strain))

##########

pc50_replicate_1 <- pc50_replicate_1[3:nrow(pc50_replicate_1),1:2]
row.names(pc50_replicate_1)<- 1:nrow(pc50_replicate_1)

a<- rep("50% moisture replicate 1",nrow(pc50_replicate_1))
a=data.frame(a)
pc50_replicate_1 = cbind(pc50_replicate_1, a)
colnames(pc50_replicate_1)=c("stress","strain","location")# naming column names
rm(a)

pc50_replicate_1$stress<- as.numeric(as.character(pc50_replicate_1$stress))
pc50_replicate_1$strain<- as.numeric(as.character(pc50_replicate_1$strain))

##########

pc50_replicate_2 <- pc50_replicate_2[3:nrow(pc50_replicate_2),1:2]
row.names(pc50_replicate_2)<- 1:nrow(pc50_replicate_2)

a<- rep("50% moisture replicate 2",nrow(pc50_replicate_2))
a=data.frame(a)
pc50_replicate_2 = cbind(pc50_replicate_2, a)
colnames(pc50_replicate_2)=c("stress","strain","location")# naming column names
rm(a)

pc50_replicate_2$stress<- as.numeric(as.character(pc50_replicate_2$stress))
pc50_replicate_2$strain<- as.numeric(as.character(pc50_replicate_2$strain))

##########

pc50_replicate_3 <- pc50_replicate_3[3:nrow(pc50_replicate_3),1:2]
row.names(pc50_replicate_3)<- 1:nrow(pc50_replicate_3)

a<- rep("50% moisture replicate 3",nrow(pc50_replicate_3))
a=data.frame(a)
pc50_replicate_3 = cbind(pc50_replicate_3, a)
colnames(pc50_replicate_3)=c("stress","strain","location")# naming column names
rm(a)

pc50_replicate_3$stress<- as.numeric(as.character(pc50_replicate_3$stress))
pc50_replicate_3$strain<- as.numeric(as.character(pc50_replicate_3$strain))

##########

pc50_replicate_4 <- pc50_replicate_4[3:nrow(pc50_replicate_4),1:2]
row.names(pc50_replicate_4)<- 1:nrow(pc50_replicate_4)

a<- rep("50% moisture replicate 4",nrow(pc50_replicate_4))
a=data.frame(a)
pc50_replicate_4 = cbind(pc50_replicate_4, a)
colnames(pc50_replicate_4)=c("stress","strain","location")# naming column names
rm(a)

pc50_replicate_4$stress<- as.numeric(as.character(pc50_replicate_4$stress))
pc50_replicate_4$strain<- as.numeric(as.character(pc50_replicate_4$strain))

##########

pc50_replicate_5 <- pc50_replicate_5[3:nrow(pc50_replicate_5),1:2]
row.names(pc50_replicate_5)<- 1:nrow(pc50_replicate_5)

a<- rep("50% moisture replicate 5",nrow(pc50_replicate_5))
a=data.frame(a)
pc50_replicate_5 = cbind(pc50_replicate_5, a)
colnames(pc50_replicate_5)=c("stress","strain","location")# naming column names
rm(a)

pc50_replicate_5$stress<- as.numeric(as.character(pc50_replicate_5$stress))
pc50_replicate_5$strain<- as.numeric(as.character(pc50_replicate_5$strain))

##########

pc50_replicate_6 <- pc50_replicate_6[3:nrow(pc50_replicate_6),1:2]
row.names(pc50_replicate_6)<- 1:nrow(pc50_replicate_6)

a<- rep("50% moisture replicate 6",nrow(pc50_replicate_6))
a=data.frame(a)
pc50_replicate_6 = cbind(pc50_replicate_6, a)
colnames(pc50_replicate_6)=c("stress","strain","location")# naming column names
rm(a)

pc50_replicate_6$stress<- as.numeric(as.character(pc50_replicate_6$stress))
pc50_replicate_6$strain<- as.numeric(as.character(pc50_replicate_6$strain))

##########

pc60_replicate_1 <- pc60_replicate_1[3:nrow(pc60_replicate_1),1:2]
row.names(pc60_replicate_1)<- 1:nrow(pc60_replicate_1)

a<- rep("60% moisture replicate 1",nrow(pc60_replicate_1))
a=data.frame(a)
pc60_replicate_1 = cbind(pc60_replicate_1, a)
colnames(pc60_replicate_1)=c("stress","strain","location")# naming column names
rm(a)

pc60_replicate_1$stress<- as.numeric(as.character(pc60_replicate_1$stress))
pc60_replicate_1$strain<- as.numeric(as.character(pc60_replicate_1$strain))

##########

pc60_replicate_2 <- pc60_replicate_2[3:nrow(pc60_replicate_2),1:2]
row.names(pc60_replicate_2)<- 1:nrow(pc60_replicate_2)

a<- rep("60% moisture replicate 2",nrow(pc60_replicate_2))
a=data.frame(a)
pc60_replicate_2 = cbind(pc60_replicate_2, a)
colnames(pc60_replicate_2)=c("stress","strain","location")# naming column names
rm(a)

pc60_replicate_2$stress<- as.numeric(as.character(pc60_replicate_2$stress))
pc60_replicate_2$strain<- as.numeric(as.character(pc60_replicate_2$strain))

##########

pc60_replicate_3 <- pc60_replicate_3[3:nrow(pc60_replicate_3),1:2]
row.names(pc60_replicate_3)<- 1:nrow(pc60_replicate_3)

a<- rep("60% moisture replicate 3",nrow(pc60_replicate_3))
a=data.frame(a)
pc60_replicate_3 = cbind(pc60_replicate_3, a)
colnames(pc60_replicate_3)=c("stress","strain","location")# naming column names
rm(a)

pc60_replicate_3$stress<- as.numeric(as.character(pc60_replicate_3$stress))
pc60_replicate_3$strain<- as.numeric(as.character(pc60_replicate_3$strain))

##########

pc60_replicate_4 <- pc60_replicate_4[3:nrow(pc60_replicate_4),1:2]
row.names(pc60_replicate_4)<- 1:nrow(pc60_replicate_4)

a<- rep("60% moisture replicate 4",nrow(pc60_replicate_4))
a=data.frame(a)
pc60_replicate_4 = cbind(pc60_replicate_4, a)
colnames(pc60_replicate_4)=c("stress","strain","location")# naming column names
rm(a)

pc60_replicate_4$stress<- as.numeric(as.character(pc60_replicate_4$stress))
pc60_replicate_4$strain<- as.numeric(as.character(pc60_replicate_4$strain))

######################################################################
## Plotting each sample individually and removing seating stress

plot(pc15_replicate_1$strain, pc15_replicate_1$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="15% moisture replicate 1")

##########

plot(pc15_replicate_2$strain, pc15_replicate_2$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="15% moisture replicate 2")

##########

plot(pc15_replicate_3$strain, pc15_replicate_3$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="15% moisture replicate 3")

##########

plot(pc15_replicate_4$strain, pc15_replicate_4$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="15% moisture replicate 4")

##########

plot(pc15_replicate_5$strain, pc15_replicate_5$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="15% moisture replicate 5")

##########

plot(pc15_replicate_6$strain, pc15_replicate_6$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="15% moisture replicate 6")

##########

plot(pc20_replicate_1$strain, pc20_replicate_1$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="20% moisture replicate 1")

##########

plot(pc20_replicate_2$strain, pc20_replicate_2$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="20% moisture replicate 2")

##########

plot(pc20_replicate_3$strain, pc20_replicate_3$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="20% moisture replicate 3")

##########

plot(pc20_replicate_4$strain, pc20_replicate_4$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="20% moisture replicate 4")

##########

plot(pc20_replicate_5$strain, pc20_replicate_5$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="20% moisture replicate 5")

##########

plot(pc20_replicate_6$strain, pc20_replicate_6$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="20% moisture replicate 6")
##########

plot(pc30_replicate_1$strain, pc30_replicate_1$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="30% moisture replicate 1")

##########

plot(pc30_replicate_2$strain, pc30_replicate_2$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="30% moisture replicate 2")

##########

plot(pc30_replicate_3$strain, pc30_replicate_3$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="30% moisture replicate 3")

##########

plot(pc30_replicate_4$strain, pc30_replicate_4$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="30% moisture replicate 4")

##########

plot(pc30_replicate_5$strain, pc30_replicate_5$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="30% moisture replicate 5")

##########

plot(pc40_replicate_1$strain, pc40_replicate_1$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="40% moisture replicate 1")

##########

plot(pc40_replicate_2$strain, pc40_replicate_2$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="40% moisture replicate 2")

##########

plot(pc40_replicate_3$strain, pc40_replicate_3$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="40% moisture replicate 3")

##########

plot(pc40_replicate_4$strain, pc40_replicate_4$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="40% moisture replicate 4")

##########

plot(pc40_replicate_5$strain, pc40_replicate_5$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="40% moisture replicate 5")

##########

plot(pc50_replicate_1$strain, pc50_replicate_1$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="50% moisture replicate 1")

##########

plot(pc50_replicate_2$strain, pc50_replicate_2$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="50% moisture replicate 2")

##########
plot(pc50_replicate_3$strain, pc50_replicate_3$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="50% moisture replicate 3")

##########

plot(pc50_replicate_4$strain, pc50_replicate_4$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="50% moisture replicate 4")

##########

plot(pc50_replicate_5$strain, pc50_replicate_5$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="50% moisture replicate 5")

##########

plot(pc50_replicate_6$strain, pc50_replicate_6$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="50% moisture replicate 6")

##########

plot(pc60_replicate_1$strain, pc60_replicate_1$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="60% moisture replicate 1")

##########

plot(pc60_replicate_2$strain, pc60_replicate_2$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="60% moisture replicate 2")

##########

plot(pc60_replicate_3$strain, pc60_replicate_3$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="60% moisture replicate 3")

##########

plot(pc60_replicate_4$strain, pc60_replicate_4$stress,
     type="l", col="red",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="60% moisture replicate 4")


##################################################
### peak stress values plotting

pc15 <- c(max(pc15_replicate_1$stress),
          max(pc15_replicate_2$stress),
          max(pc15_replicate_3$stress),
          max(pc15_replicate_4$stress),
          max(pc15_replicate_5$stress),
          max(pc15_replicate_6$stress))

pc20 <- c(max(pc20_replicate_1$stress),
          max(pc20_replicate_2$stress),
          max(pc20_replicate_3$stress),
          max(pc20_replicate_4$stress),
          max(pc20_replicate_5$stress),
          max(pc20_replicate_6$stress))

pc30 <- c(max(pc30_replicate_1$stress),
          max(pc30_replicate_2$stress),
          max(pc30_replicate_3$stress),
          max(pc30_replicate_4$stress),
          max(pc30_replicate_5$stress))

pc40 <- c(max(pc40_replicate_1$stress),
          max(pc40_replicate_2$stress),
          max(pc40_replicate_3$stress),
          max(pc40_replicate_4$stress),
          max(pc40_replicate_5$stress))

pc50 <- c(max(pc50_replicate_1$stress),
          max(pc50_replicate_2$stress),
          max(pc50_replicate_3$stress),
          max(pc50_replicate_4$stress),
          max(pc50_replicate_5$stress),
          max(pc50_replicate_6$stress))

pc60 <- c(max(pc60_replicate_1$stress),
          max(pc60_replicate_2$stress),
          max(pc60_replicate_3$stress),
          max(pc60_replicate_4$stress))

##########
## Importing data

data <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 4 - Effect of moisture content on self-weight consolidation of soil/Effect of moisture content on self-weight consolidation of soil - statistics.csv")

mod1<- lm(Stress[1:20] ~ Moisture[1:20], data=data)
summary(mod1)
shapiro.test(residuals(mod1))

mod2<- aov(Stress[1:20] ~ Moisture[1:20], data=data)
summary(mod2)
TukeyHSD(mod2)

mod3<- lm(Stress[21:32] ~ Moisture[21:32], data=data)
summary(mod3)
shapiro.test(residuals(mod3))

mod4<- aov(Stress[21:32] ~ Moisture[21:32], data=data)
summary(mod4)
TukeyHSD(mod4)
