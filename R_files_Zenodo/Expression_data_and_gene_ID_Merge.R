#setwd
# Merge with Fisher annotated data
annotated<-read.csv(" E.carinata Annotated_Trimmed_3.8.21.csv", header = TRUE, sep=",")

# This code removes the repeated rows (everything in the row has to be the same)
annotated<-unique(annotated)
# selecting columns to keep
annotated.new<-annotated[,c(1, 5)]
# removes repeateed rows
annotated.new<-unique(annotated.new)

# Remove all rows that don't contain a FlyID - end up with 4646 annotated genes
# Contains genes that have been cross checked and selected -- see File named Fisher annotated Trimmed 3.8.21
# for details on the duplicates that remained in the dataset and why they were selected over other duplicate
# Also details which duplicates were removed altogether
annotated.remove.NA<-subset(annotated.new, annotated.new$GeneID != "")

#1 Merge Anterior-Legs by TrinityID, All
AL.expression<-read.csv("DE_Anterior-Legs.csv", header=T)
AL_dup<-AL.expression[duplicated(AL.expression$GeneID) | duplicated(AL.expression$GeneID, fromLast = TRUE),]
AL.expression.merge<-merge(AL.expression, annotated.remove.NA, by="TrinityID", sort=TRUE)

#2 Merge Posterior-Legs by TrinityID, All
PL.expression<-read.csv("DE_Posterior-Legs.csv", header=T)
PL_dup<-PL.expression[duplicated(PL.expression$GeneID) | duplicated(PL.expression$GeneID, fromLast = TRUE),]
PL.expression.merge<-merge(PL.expression, annotated.remove.NA, by="TrinityID", sort=TRUE)

#3 Merge Anterior-Posterior by TrinityID, All
AP.expression<-read.csv("DE_Anterior-Posterior.csv", header=T)
AP_dup<-AP.expression[duplicated(AP.expression$GeneID) | duplicated(AP.expression$GeneID, fromLast = TRUE),]
AP.expression.merge<-merge(AP.expression, annotated.remove.NA, by="TrinityID", sort=TRUE)

#4 Merge Anterior-Wing by TrinityID, All
AW.expression<-read.csv("DE_Anterior-Wings.csv", header=T)
AW_dup<-AW.expression[duplicated(AW.expression$TrinityID) | duplicated(AW.expression$GeneID, fromLast = TRUE),]
AW.expression.merge<-merge(AW.expression, annotated.remove.NA, by="TrinityID", sort=TRUE)

#5 Merge Posterior-Wing by TrinityID, All
PW.expression<-read.csv("DE_Posterior-Wings.csv", header=T)
PW_dup<-PW.expression[duplicated(PW.expression$TrinityID) | duplicated(PW.expression$TrinityID, fromLast = TRUE),]
PW.expression.merge<-merge(PW.expression, annotated.remove.NA, by="TrinityID", sort=TRUE)

#6 Merge Legs-Wing by TrinityID, All
LW.expression<-read.csv("DE_Wings-Legs.csv", header=T)
LW_dup<-LW.expression[duplicated(LW.expression$TrinityID) | duplicated(LW.expression$TrinityID, fromLast = TRUE),]
LW.expression.merge<-merge(LW.expression, annotated.remove.NA, by="TrinityID", sort=TRUE)

# Export All files
#1
write.csv(LW.expression.merge, "Leg-Wing_3.13.21.csv")
#2
write.csv(PW.expression.merge, "Posterior-Wing_3.13.21.csv")
#3
write.csv(AW.expression.merge, "Anterior-Wing_3.13.21.csv")
#4
write.csv(AP.expression.merge, "Anterior-Posterior_3.13.21.csv")
#5
write.csv(PL.expression.merge, "Posterior-Leg_3.13.21.csv")
#6
write.csv(AL.expression.merge, "Anterior-Leg_3.13.21.csv")
