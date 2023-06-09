#### TITLE: Scripts used to create datasets, from the raw data, used for analysis  

# Paper: Neonate personality affects early-life resource acquisition in a large social mammal

# Authors: Bawan Amin, Dómhnall J. Jennings, Alison Norman, Andrew Ryan, Vasiliki Ioannidis, Alice Magee, Hayley-Anne Haughey, Amy Haigh, Simone Ciuti


# Description: The following script has the function of creating the 4 final datasets that are used for the bivariate models reported in Amin et al. - Neonate personality affects early-life resource acquisition in a large social mammal

# The script starts with the raw data and transforms it, step by step, into the final datasets. Steps are annotated and if still unclear, then any query can be sent to the main author, Bawan Amin (bawan.amin@ucdconnect.ie)



# Load raw dataset and summarize 
db = read.csv("Raw_data.csv")
summary(db)


# Change structure of certain columns
db$FawnID = as.factor(db$FawnID)
db$Variable = as.factor(db$Variable)
db$Beh = as.factor(db$Beh)
db$Season = as.factor(db$Season)
db$Sex = as.factor(db$Sex)
db$Capture = as.factor(db$Capture)
db$Year = as.factor(db$Year)
db$Group = as.factor(db$Group)

summary(db)

# Model complained about NAs in the fixed factors. NAs made sense: variables like meancap, Airtemp and Weight are only available for Heartrate, whereas the other variables in the model are only available for Time Budgets. Here I have filled out the NAs with the mean of the column, but in the model specify which variables should be used for which response. 
db$meancap = ifelse(is.na(db$meancap), 0.38, db$meancap  )
db$X..Deer = ifelse (is.na(db$X..Deer), mean(db$X..Deer, na.rm = T), db$X..Deer)
db$Air.temperature = ifelse (is.na(db$Air.temperature), mean(db$Air.temperature, na.rm = T), db$Air.temperature)
db$Weight = ifelse (is.na(db$Weight), mean(db$Weight, na.rm = T), db$Weight)
db$JD = ifelse (is.na(db$JD), mean(db$JD, na.rm = T), db$JD)
db$X..People = ifelse (is.na(db$X..People), mean(db$X..People, na.rm = T), db$X..People)
db$HRtime = ifelse (is.na(db$HRtime), mean(db$HRtime, na.rm = T), db$HRtime)
db$Capture = ifelse (is.na(db$Capture), 7, db$Capture)
db$X..Bucks = ifelse (is.na(db$X..Bucks), mean(db$X..Bucks, na.rm = T), db$X..Bucks)
db$BDAY = ifelse (is.na(db$BDAY), mean(db$BDAY, na.rm = T), db$BDAY)
db$Capture = as.factor(db$Capture)
db$Position = ifelse (is.na(db$Position), mean(db$Position, na.rm = T), db$Position)



#Divide datasets up into seasons and behaviours, and apply transformation on the response traits.  

# Focal observation data
Try = db[db$Response > 0,]
e = min(Try$Response)

logitTransform <- function(p) { log( (p+e) / ((1-p)+e) ) }

SummerVig = subset (db, Season == "Summer2018" | Season == "Summer2019")
SummerVig = subset (SummerVig, Beh == "Vigilance")

# Logits of 0 lead to -Inf. Therefore, we change the zero's into half of the minimum, as has been suggested when applying the logit transformation

SummerVig$Response = logitTransform(SummerVig$Response)
SummerVig$Response = scale(SummerVig$Response)


AutumnVig = subset (db, Season == "Autumnwinter2018" | Season == "Autumnwinter2019")
AutumnVig = subset (AutumnVig, Beh == "Vigilance")
Try = AutumnVig[AutumnVig$Response > 0,]
e = min(Try$Response)

AutumnVig$Response = logitTransform(AutumnVig$Response)
AutumnVig$Response = scale(AutumnVig$Response)



# Neonate capture data
CaptureHR = subset(db, Beh == "HRend")
CaptureHR$Response = log(CaptureHR$Response)
CaptureHR$Response = scale (CaptureHR$Response)


CaptureLat = subset(db, Beh == "Latency")
CaptureLat$Response = ifelse (CaptureLat$Response > 9, 10, CaptureLat$Response) #See Amin et al, Functional Ecology, 2021
CaptureLat$Response = log(CaptureLat$Response+1)
CaptureLat$Response = scale(CaptureLat$Response)



# Compiling seasons together

Trivar_Autumn = rbind(CaptureHR, CaptureLat, AutumnVig)
Trivar_Summer = rbind(CaptureHR, CaptureLat, SummerVig)

# Make a list of animals sampled in the scanning collections 
AutumnVig$Observed = 1
SummerVig$Observed = 1

Autumn = na.omit(as.data.frame(tapply(AutumnVig$Observed, AutumnVig$FawnID, min)))
Autumn$FawnID = rownames(Autumn)
colnames(Autumn)[1] = "Autumn"

Summer = na.omit(as.data.frame(tapply(SummerVig$Observed, SummerVig$FawnID, min)))
Summer$FawnID = rownames(Summer)
colnames(Summer)[1] = "Summer"


#Censor all data where animals weren't sampled in both 
names(Trivar_Autumn)

Trivar_Autumn = merge (Trivar_Autumn, Autumn, by = "FawnID")
Trivar_Autumn = subset (Trivar_Autumn, Autumn == 1)

Trivar_Autumn$Variable = Trivar_Autumn$Beh

Trivar_Autumn = Trivar_Autumn[, c("FawnID", "Response", "Variable", "meancap", "Weight", "Year", "Capture", "Season", "JD", "Time", "X..People", "X..Deer", "BDAY", "Sex", "Position", "X..Bucks", "Air.temperature", "Duration")]

Trivar_Autumn$Duration = ifelse(is.na(Trivar_Autumn$Duration), 1, Trivar_Autumn$Duration)

Trivar_Autumn = na.omit(Trivar_Autumn)


names(Trivar_Summer)
Trivar_Summer = merge (Trivar_Summer, Summer, by = "FawnID")
Trivar_Summer = subset (Trivar_Summer, Summer == 1)

Trivar_Summer$Variable = Trivar_Summer$Beh

Trivar_Summer = Trivar_Summer[, c("FawnID", "Response", "Variable", "meancap", "Weight", "Year", "Capture", "Season", "JD", "Time", "X..People", "X..Deer", "BDAY", "Sex", "Position", "X..Bucks", "Air.temperature", "Duration")]

Trivar_Summer$Duration = ifelse(is.na(Trivar_Summer$Duration), 1, Trivar_Summer$Duration)

Trivar_Summer = na.omit(Trivar_Summer)


# Towards final datasets
summary(Trivar_Autumn)
Trivar_Autumn$Variable = as.factor(Trivar_Autumn$Variable)
Trivar_Autumn$Season = as.factor(Trivar_Autumn$Season)
Trivar_Autumn$Sex = as.factor(Trivar_Autumn$Sex)
Trivar_Autumn$Year = as.factor(Trivar_Autumn$Year)
summary(Trivar_Autumn)

Emergence = read.csv("Herd_Emergence.csv")
Emergence = Emergence [, c("FawnID", "Emergence")]
Trivar_Autumn = merge (Trivar_Autumn, Emergence, by = "FawnID")
Trivar_Autumn$Days_in_herd = Trivar_Autumn$JD - Trivar_Autumn$Emergence

HR_Autumn_final = subset(Trivar_Autumn, Variable == "HRend" | Variable == "Vigilance")
Latency_Autumn_final = subset(Trivar_Autumn, Variable == "Latency" | Variable == "Vigilance")

summary(Trivar_Summer)
Trivar_Summer$Variable = as.factor(Trivar_Summer$Variable)
Trivar_Summer$Season = as.factor(Trivar_Summer$Season)
Trivar_Summer$Sex = as.factor(Trivar_Summer$Sex)
Trivar_Summer$Year = as.factor(Trivar_Summer$Year)
summary(Trivar_Summer)

Emergence = read.csv("Herd_Emergence.csv")
Emergence = Emergence [, c("FawnID", "Emergence")]
Trivar_Summer = merge (Trivar_Summer, Emergence, by = "FawnID", all.x = T)
Trivar_Summer$Days_in_herd = Trivar_Summer$JD - Trivar_Summer$Emergence

HR_Summer_final = subset(Trivar_Summer, Variable == "HRend" | Variable == "Vigilance")
Latency_Summer_final = subset(Trivar_Summer, Variable == "Latency" | Variable == "Vigilance")


# Uncomment and run lines below to write the datasets into csv files. 

write.csv(HR_Autumn_final, "HR_Autumn_final.csv", row.names = F)
write.csv(Latency_Autumn_final, "Latency_Autumn_final.csv", row.names = F)
write.csv(HR_Summer_final, "HR_Summer_final.csv", row.names = F)
write.csv(Latency_Summer_final, "Latency_Summer_final.csv", row.names = F)
