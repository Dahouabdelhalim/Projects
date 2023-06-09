###########Home range max overlap determination

rm(list=ls()) 

library(plyr)
library(reshape2)

##import csvs
voles0 = read.csv("2016_Home Range Overlap_ALL.csv")
voleinfo = read.csv(file.choose())
head(voles0)
head(voleinfo)

##merge vole data to each column to get sexes
voles1 = merge(voles0, voleinfo,
              by = c("Vole.ID"))
head(voles1)
voles = merge(voles1, voleinfo,
              by.x = "Overlapping.ID",
              by.y = "Vole.ID")
head(voles)
###subset to necessary columns only & only M-F pairs
voles = subset(voles, select=c("HR.Period", "Enclosure.x", "Vole.ID", "Sex.x", 
                               "Overlapping.ID", "Sex.y", "Overlap"))
      
voles = subset(voles, !voles$Sex.x == voles$Sex.y)
head(voles, 20)

####now, convert to only max overlaps- but may have doubles####
Volemax = ddply( voles,. (HR.Period, Vole.ID),
                 summarise, MaxOver = max(Overlap),
                 Max2 = sort(Overlap, decreasing=T)[2],
                 sum = sum(Overlap))
head(Volemax, 20)
Volemax$sub = (Volemax$MaxOver - Volemax$Max2)
Volemax$Propo = (Volemax$MaxOver - (Volemax$sum-Volemax$MaxOver) )

hist(Volemax$Propo)
hist(Volemax$sub)
hist(Volemax$MaxOver)

summary(Volemax$Max2)
boxplot(Volemax$Max2)


write.csv(Volemax, "HR_Overlap trial.csv")
###add associated vole-but may be doubles
Max = merge(Volemax, voles, 
            by.x = c("HR.Period", "Vole.ID", "MaxOver"),
            by.y = c("HR.Period", "Vole.ID", "Overlap"), 
            all.x = TRUE)
head(Max)

###now to find/remove doubles (unlikely but just in case)
doubles = ddply( Max,. (HR.Period, Vole.ID),
                 summarise, Count = length(Vole.ID))
head(doubles)

###add counts to the voles file, exclude doubles
voles2 = merge(voles, doubles, 
               by = c("HR.Period", "Vole.ID"))
head(voles2, 30)
subset(voles2, Count >1)
voles3 = voles2[voles2$Count == 1,]


##now merge again, doesn't include doubles###
head(Volemax, 30); head(voles3)

Max2 = merge(Volemax, voles3, 
            by.x = c("HR.Period", "Vole.ID", "MaxOver"),
            by.y = c("HR.Period", "Vole.ID", "Overlap"), 
            all.x = TRUE)

Max2$Sex.y = NULL; Max2$Count = NULL; Max2$Sex.x = NULL;
head(Max2, 20)




##############Now, determine # of females overlapped!!!###########
volesover = voles[!voles$Overlap == 0,]  ##remove overlaps of 0
head(volesover)
#determine # of opposite sex overlapped for each vole
Countover = ddply(volesover,. (HR.Period, Enclosure.x,Vole.ID, Sex.x),
                  summarise, Count = length(unique(Overlapping.ID)))

##now add sex for overlapping voles again
Countover$Sex.y = ifelse(Countover$Sex.x == "M", "F", "M")
head(Countover)
summary(Countover)

##determine total potential overlaps for each vole
MFcount = ddply(voles,. (HR.Period, Enclosure.x, Sex.x),
                  summarise, Tot = length(unique(Vole.ID)))
head(MFcount)

#merge them together
Countover2 = merge(Countover, MFcount, 
                   by.x = c("HR.Period", "Enclosure.x", "Sex.y"),
                   by.y = c("HR.Period", "Enclosure.x", "Sex.x"),
                   all.x = TRUE)
head(Countover2)

###Now calculate proportion overlapped
Countover2$Prop = Countover2$Count/Countover2$Tot
head(Countover2)
Countover2$Sex.y = NULL; Countover2$Tot = NULL
    


###now merge, make csv
head(Max2); head(Countover2)
Final = merge(Max2, Countover2, by = c("HR.Period", "Vole.ID"), all = TRUE)
Final$Enclosure.x.x = NULL; Final$Sex.x = NULL; Final$Enclosure.x.y = NULL
summary(Final)
head(Final)

#make 0 if no overlaps
Final[c("Prop")][is.na(Final[c("Prop")])] = 0


write.csv(Final, "HR_Max Overlap 2016.csv")
