#######################################################
######## Importing data
Bolus_1 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Bolus replicate 1.csv")
Bolus_2 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Bolus replicate 2.csv")
Bolus_3 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Bolus replicate 3.csv")
Bolus_4 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Bolus replicate 4.csv")
Bolus_5 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Bolus replicate 5.csv")
Bolus_6 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Bolus replicate 6.csv")
Bolus_7 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Bolus replicate 7.csv")
Bolus_8 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Bolus replicate 8.csv")
Bolus_9 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Bolus replicate 9.csv")
Bolus_10 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Bolus replicate 10.csv")


Control_1 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 1.csv")
Control_2 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 2.csv")
Control_3 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 3.csv")
Control_4 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 4.csv")
Control_5 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 5.csv")
Control_6 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 6.csv")
Control_7 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 7.csv")
Control_8 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 8.csv")
Control_9 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 9.csv")
Control_10 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 10.csv")
Control_11 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 11.csv")
Control_12 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 12.csv")
Control_13 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 13.csv")
Control_14 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 14.csv")
Control_15 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 15.csv")
Control_16 <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/Control replicate 16.csv")



# Step1: Eliminating unnecessary columns and rows and reassigning row numbers
# Step2: Creating new columns representing size of samples and binding all data sets into a single data set
# Step 3: Declaring stress and strain as numeric and character vectors
# Step 4: Creating new columns representing size of samples and binding all data sets into a single data set
# Step 5: Graph plotting

Bolus_1 <- Bolus_1[4:nrow(Bolus_1),1:2]
row.names(Bolus_1)<- 1:nrow(Bolus_1)

a<- rep("Bolus_1",nrow(Bolus_1))
a=data.frame(a)
Bolus_1 = cbind(Bolus_1, a)
colnames(Bolus_1)=c("stress","strain","location")# naming column names
rm(a)

Bolus_1$stress<- as.numeric(as.character(Bolus_1$stress))
Bolus_1$strain<- as.numeric(as.character(Bolus_1$strain))

plot(Bolus_1$strain, Bolus_1$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Bolus_1")

####

Bolus_2 <- Bolus_2[4:nrow(Bolus_2),1:2]
row.names(Bolus_2)<- 1:nrow(Bolus_2)

a<- rep("Bolus_2",nrow(Bolus_2))
a=data.frame(a)
Bolus_2 = cbind(Bolus_2, a)
colnames(Bolus_2)=c("stress","strain","location")# naming column names
rm(a)

Bolus_2$stress<- as.numeric(as.character(Bolus_2$stress))
Bolus_2$strain<- as.numeric(as.character(Bolus_2$strain))

plot(Bolus_2$strain, Bolus_2$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Bolus_2")

####

Bolus_3 <- Bolus_3[4:nrow(Bolus_3),1:2]
row.names(Bolus_3)<- 1:nrow(Bolus_3)

a<- rep("Bolus_3",nrow(Bolus_3))
a=data.frame(a)
Bolus_3 = cbind(Bolus_3, a)
colnames(Bolus_3)=c("stress","strain","location")# naming column names
rm(a)

Bolus_3$stress<- as.numeric(as.character(Bolus_3$stress))
Bolus_3$strain<- as.numeric(as.character(Bolus_3$strain))

plot(Bolus_3$strain, Bolus_3$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Bolus_3")

####

Bolus_4 <- Bolus_4[4:nrow(Bolus_4),1:2]
row.names(Bolus_4)<- 1:nrow(Bolus_4)

a<- rep("Bolus_4",nrow(Bolus_4))
a=data.frame(a)
Bolus_4 = cbind(Bolus_4, a)
colnames(Bolus_4)=c("stress","strain","location")# naming column names
rm(a)

Bolus_4$stress<- as.numeric(as.character(Bolus_4$stress))
Bolus_4$strain<- as.numeric(as.character(Bolus_4$strain))

plot(Bolus_4$strain, Bolus_4$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Bolus_4")

####

Bolus_5 <- Bolus_5[4:nrow(Bolus_5),1:2]
row.names(Bolus_5)<- 1:nrow(Bolus_5)

a<- rep("Bolus_5",nrow(Bolus_5))
a=data.frame(a)
Bolus_5 = cbind(Bolus_5, a)
colnames(Bolus_5)=c("stress","strain","location")# naming column names
rm(a)

Bolus_5$stress<- as.numeric(as.character(Bolus_5$stress))
Bolus_5$strain<- as.numeric(as.character(Bolus_5$strain))

plot(Bolus_5$strain, Bolus_5$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Bolus_5")

####

Bolus_6 <- Bolus_6[4:nrow(Bolus_6),1:2]
row.names(Bolus_6)<- 1:nrow(Bolus_6)

a<- rep("Bolus_6",nrow(Bolus_6))
a=data.frame(a)
Bolus_6 = cbind(Bolus_6, a)
colnames(Bolus_6)=c("stress","strain","location")# naming column names
rm(a)

Bolus_6$stress<- as.numeric(as.character(Bolus_6$stress))
Bolus_6$strain<- as.numeric(as.character(Bolus_6$strain))

plot(Bolus_6$strain, Bolus_6$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Bolus_6")

####

Bolus_7 <- Bolus_7[4:nrow(Bolus_7),1:2]
row.names(Bolus_7)<- 1:nrow(Bolus_7)

a<- rep("Bolus_7",nrow(Bolus_7))
a=data.frame(a)
Bolus_7 = cbind(Bolus_7, a)
colnames(Bolus_7)=c("stress","strain","location")# naming column names
rm(a)

Bolus_7$stress<- as.numeric(as.character(Bolus_7$stress))
Bolus_7$strain<- as.numeric(as.character(Bolus_7$strain))

plot(Bolus_7$strain, Bolus_7$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Bolus_7")

####

Bolus_8 <- Bolus_8[4:nrow(Bolus_8),1:2]
row.names(Bolus_8)<- 1:nrow(Bolus_8)

a<- rep("Bolus_8",nrow(Bolus_8))
a=data.frame(a)
Bolus_8 = cbind(Bolus_8, a)
colnames(Bolus_8)=c("stress","strain","location")# naming column names
rm(a)

Bolus_8$stress<- as.numeric(as.character(Bolus_8$stress))
Bolus_8$strain<- as.numeric(as.character(Bolus_8$strain))

plot(Bolus_8$strain, Bolus_8$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Bolus_8")

####

Bolus_9 <- Bolus_9[4:nrow(Bolus_9),1:2]
row.names(Bolus_9)<- 1:nrow(Bolus_9)

a<- rep("Bolus_9",nrow(Bolus_9))
a=data.frame(a)
Bolus_9 = cbind(Bolus_9, a)
colnames(Bolus_9)=c("stress","strain","location")# naming column names
rm(a)

Bolus_9$stress<- as.numeric(as.character(Bolus_9$stress))
Bolus_9$strain<- as.numeric(as.character(Bolus_9$strain))

plot(Bolus_9$strain, Bolus_9$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Bolus_9")

####

Bolus_10 <- Bolus_10[4:nrow(Bolus_10),1:2]
row.names(Bolus_10)<- 1:nrow(Bolus_10)

a<- rep("Bolus_10",nrow(Bolus_10))
a=data.frame(a)
Bolus_10 = cbind(Bolus_10, a)
colnames(Bolus_10)=c("stress","strain","location")# naming column names
rm(a)

Bolus_10$stress<- as.numeric(as.character(Bolus_10$stress))
Bolus_10$strain<- as.numeric(as.character(Bolus_10$strain))

plot(Bolus_10$strain, Bolus_10$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Bolus_10")

####

Control_1 <- Control_1[4:nrow(Control_1),1:2]
row.names(Control_1)<- 1:nrow(Control_1)

a<- rep("Control_1",nrow(Control_1))
a=data.frame(a)
Control_1 = cbind(Control_1, a)
colnames(Control_1)=c("stress","strain","location")# naming column names
rm(a)

Control_1$stress<- as.numeric(as.character(Control_1$stress))
Control_1$strain<- as.numeric(as.character(Control_1$strain))

plot(Control_1$strain, Control_1$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_1")

####

Control_2 <- Control_2[4:nrow(Control_2),1:2]
row.names(Control_2)<- 1:nrow(Control_2)

a<- rep("Control_2",nrow(Control_2))
a=data.frame(a)
Control_2 = cbind(Control_2, a)
colnames(Control_2)=c("stress","strain","location")# naming column names
rm(a)

Control_2$stress<- as.numeric(as.character(Control_2$stress))
Control_2$strain<- as.numeric(as.character(Control_2$strain))

plot(Control_2$strain, Control_2$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_2")

####

Control_3 <- Control_3[4:nrow(Control_3),1:2]
row.names(Control_3)<- 1:nrow(Control_3)

a<- rep("Control_3",nrow(Control_3))
a=data.frame(a)
Control_3 = cbind(Control_3, a)
colnames(Control_3)=c("stress","strain","location")# naming column names
rm(a)

Control_3$stress<- as.numeric(as.character(Control_3$stress), na.omit=TRUE)
Control_3$strain<- as.numeric(as.character(Control_3$strain), na.omit=TRUE)

Control_3$stress[complete.cases(Control_3$stress)]

plot(Control_3$strain, Control_3$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_3")

####

Control_4 <- Control_4[4:nrow(Control_4),1:2]
row.names(Control_4)<- 1:nrow(Control_4)

a<- rep("Control_4",nrow(Control_4))
a=data.frame(a)
Control_4 = cbind(Control_4, a)
colnames(Control_4)=c("stress","strain","location")# naming column names
rm(a)

Control_4$stress<- as.numeric(as.character(Control_4$stress))
Control_4$strain<- as.numeric(as.character(Control_4$strain))

plot(Control_4$strain, Control_4$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_4")

####

Control_5 <- Control_5[4:nrow(Control_5),1:2]
row.names(Control_5)<- 1:nrow(Control_5)

a<- rep("Control_5",nrow(Control_5))
a=data.frame(a)
Control_5 = cbind(Control_5, a)
colnames(Control_5)=c("stress","strain","location")# naming column names
rm(a)

Control_5$stress<- as.numeric(as.character(Control_5$stress))
Control_5$strain<- as.numeric(as.character(Control_5$strain))

plot(Control_5$strain, Control_5$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_5")

####

Control_6 <- Control_6[4:nrow(Control_6),1:2]
row.names(Control_6)<- 1:nrow(Control_6)

a<- rep("Control_6",nrow(Control_6))
a=data.frame(a)
Control_6 = cbind(Control_6, a)
colnames(Control_6)=c("stress","strain","location")# naming column names
rm(a)

Control_6$stress<- as.numeric(as.character(Control_6$stress))
Control_6$strain<- as.numeric(as.character(Control_6$strain))

plot(Control_6$strain, Control_6$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_6")

####

Control_7 <- Control_7[4:nrow(Control_7),1:2]
row.names(Control_7)<- 1:nrow(Control_7)

a<- rep("Control_7",nrow(Control_7))
a=data.frame(a)
Control_7 = cbind(Control_7, a)
colnames(Control_7)=c("stress","strain","location")# naming column names
rm(a)

Control_7$stress<- as.numeric(as.character(Control_7$stress))
Control_7$strain<- as.numeric(as.character(Control_7$strain))

plot(Control_7$strain, Control_7$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_7")

####

Control_8 <- Control_8[4:nrow(Control_8),1:2]
row.names(Control_8)<- 1:nrow(Control_8)

a<- rep("Control_8",nrow(Control_8))
a=data.frame(a)
Control_8 = cbind(Control_8, a)
colnames(Control_8)=c("stress","strain","location")# naming column names
rm(a)

Control_8$stress<- as.numeric(as.character(Control_8$stress))
Control_8$strain<- as.numeric(as.character(Control_8$strain))

####

Control_9 <- Control_9[4:nrow(Control_9),1:2]
row.names(Control_9)<- 1:nrow(Control_9)

a<- rep("Control_9",nrow(Control_9))
a=data.frame(a)
Control_9 = cbind(Control_9, a)
colnames(Control_9)=c("stress","strain","location")# naming column names
rm(a)

Control_9$stress<- as.numeric(as.character(Control_9$stress))
Control_9$strain<- as.numeric(as.character(Control_9$strain))

plot(Control_9$strain, Control_9$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_9")

####

Control_10 <- Control_10[4:nrow(Control_10),1:2]
row.names(Control_10)<- 1:nrow(Control_10)

a<- rep("Control_10",nrow(Control_10))
a=data.frame(a)
Control_10 = cbind(Control_10, a)
colnames(Control_10)=c("stress","strain","location")# naming column names
rm(a)

Control_10$stress<- as.numeric(as.character(Control_10$stress))
Control_10$strain<- as.numeric(as.character(Control_10$strain))

plot(Control_10$strain, Control_10$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_10")

####

Control_11 <- Control_11[4:nrow(Control_11),1:2]
row.names(Control_11)<- 1:nrow(Control_11)

a<- rep("Control_11",nrow(Control_11))
a=data.frame(a)
Control_11 = cbind(Control_11, a)
colnames(Control_11)=c("stress","strain","location")# naming column names
rm(a)

Control_11$stress<- as.numeric(as.character(Control_11$stress))
Control_11$strain<- as.numeric(as.character(Control_11$strain))

plot(Control_11$strain, Control_11$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_11")

####

Control_12 <- Control_12[4:nrow(Control_12),1:2]
row.names(Control_12)<- 1:nrow(Control_12)

a<- rep("Control_12",nrow(Control_12))
a=data.frame(a)
Control_12 = cbind(Control_12, a)
colnames(Control_12)=c("stress","strain","location")# naming column names
rm(a)

Control_12$stress<- as.numeric(as.character(Control_12$stress))
Control_12$strain<- as.numeric(as.character(Control_12$strain))

plot(Control_12$strain, Control_12$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_12")

####

Control_13 <- Control_13[4:nrow(Control_13),1:2]
row.names(Control_13)<- 1:nrow(Control_13)

a<- rep("Control_13",nrow(Control_13))
a=data.frame(a)
Control_13 = cbind(Control_13, a)
colnames(Control_13)=c("stress","strain","location")# naming column names
rm(a)

Control_13$stress<- as.numeric(as.character(Control_13$stress))
Control_13$strain<- as.numeric(as.character(Control_13$strain))

plot(Control_13$strain, Control_13$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_13")

####

Control_14 <- Control_14[4:nrow(Control_14),1:2]
row.names(Control_14)<- 1:nrow(Control_14)

a<- rep("Control_14",nrow(Control_14))
a=data.frame(a)
Control_14 = cbind(Control_14, a)
colnames(Control_14)=c("stress","strain","location")# naming column names
rm(a)

Control_14$stress<- as.numeric(as.character(Control_14$stress))
Control_14$strain<- as.numeric(as.character(Control_14$strain))

plot(Control_14$strain, Control_14$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_14")

####

Control_15 <- Control_15[4:nrow(Control_15),1:2]
row.names(Control_15)<- 1:nrow(Control_15)

a<- rep("Control_15",nrow(Control_15))
a=data.frame(a)
Control_15 = cbind(Control_15, a)
colnames(Control_15)=c("stress","strain","location")# naming column names
rm(a)

Control_15$stress<- as.numeric(as.character(Control_15$stress))
Control_15$strain<- as.numeric(as.character(Control_15$strain))

plot(Control_15$strain, Control_15$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_15")

####

Control_16 <- Control_16[4:nrow(Control_16),1:2]
row.names(Control_16)<- 1:nrow(Control_16)

a<- rep("Control_16",nrow(Control_16))
a=data.frame(a)
Control_16 = cbind(Control_16, a)
colnames(Control_16)=c("stress","strain","location")# naming column names
rm(a)

Control_16$stress<- as.numeric(as.character(Control_16$stress))
Control_16$strain<- as.numeric(as.character(Control_16$strain))

plot(Control_16$strain, Control_16$stress,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="Control_16")

##################################################
### Getting peak stress values

max(Bolus_1$stress)
max(Bolus_2$stress)
max(Bolus_3$stress)
max(Bolus_4$stress)
max(Bolus_5$stress)
max(Bolus_6$stress)
max(Bolus_7$stress)
max(Bolus_8$stress)
max(Bolus_9$stress)
max(Bolus_10$stress)

max(Control_1$stress)
max(Control_2$stress)
max(Control_3$stress[complete.cases(Control_3$stress)])
max(Control_4$stress)
max(Control_5$stress)
max(Control_6$stress)
max(Control_7$stress)
max(Control_8$stress)
max(Control_9$stress)
max(Control_10$stress)
max(Control_11$stress)
max(Control_12$stress)
max(Control_13$stress)
max(Control_14$stress)
max(Control_15$stress)
max(Control_16$stress)


#################
#### All in-situ samples

All_in_situ_samples <- read.csv("~/Desktop/Strength and weathering resistance of termite mounds/Data for Figure 6 - Contribution of termite secretions towards strength of aggregated soil/All in-situ samples.csv")
head(All_in_situ_samples)

stress<- All_in_situ_samples$Peak.compressive.stress..kPa.
stress<- stress[complete.cases(stress)]
length(stress)
strain<- All_in_situ_samples$Strain....
strain<- strain[complete.cases(strain)]
length(strain)


## In-situ sample 1

stress_sample_1 <- stress[1:160280]
strain_sample_1 <- strain[1:160280]

plot(strain_sample_1, stress_sample_1,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="In-situ Sample 1")

max(stress_sample_1)  ## In-situ sample 1

## In-situ sample 2

stress_sample_2 <- stress[160281:311101]
strain_sample_2 <- strain[160281:311101]

plot(strain_sample_2, stress_sample_2,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="In-situ Sample 2")

max(stress_sample_2)  ## In-situ sample 2

## In-situ sample 3

stress_sample_3 <- stress[321702:421121]
strain_sample_3 <- strain[321702:421121]

plot(strain_sample_3, stress_sample_3,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="In-situ Sample 3")

max(stress_sample_3)  ## In-situ sample 3

## In-situ sample 4

stress_sample_4 <- stress[471121:583401]
strain_sample_4 <- strain[471121:583401]

plot(strain_sample_4, stress_sample_4,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="In-situ Sample 4")

max(stress_sample_4)  ## In-situ sample 4

## In-situ sample 5

stress_sample_5 <- stress[603401:718000]
strain_sample_5 <- strain[603401:718000]

plot(strain_sample_5, stress_sample_5,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="In-situ Sample 5")

max(stress_sample_5)  ## In-situ sample 5

## In-situ sample 6

stress_sample_6 <- stress[735000:845880]
strain_sample_6 <- strain[735000:845880]

plot(strain_sample_6, stress_sample_6,
     type="l", col="blue",
     xlab="Strain (%)", ylab="Stress (kPa)",
     main="In-situ Sample 6")

max(stress_sample_6)  ## In-situ sample 6



################
### Compiling values from the above data

Bolus<- c(1104.562,
          1024.609,
          759.6703,
          933.2659,
          1450.803,
          996.1583,
          1416.017,
          1063.847,
          1712.036,
          1417.967)

In_situ<- c(1073.576,
            1545.049,
            2222.666,
            1187.868,
            844.9758,
            814.3635)

Control<- c(1566.716,
            787.6314,
            1605.056,
            1371.458,
            1524.02,
            1162.912,
            2017.758,
            2509.146,
            2298.899,
            4342.053,
            2405.814,
            851.8726,
            1165.44,
            572.1303,
            789.4288,
            698.665)

boxplot(Bolus, In_situ, Control)
a<- c(Bolus, In_situ, Control)
shapiro.test(a)

wilcox.test(Bolus, In_situ, paired=FALSE)
wilcox.test(Bolus, Control, paired=FALSE)
wilcox.test(Control, In_situ, paired=FALSE)
