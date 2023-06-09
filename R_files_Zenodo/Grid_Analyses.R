#Clean the R environment - remove all previously defined objects
rm(list=ls()) 

# Read in data from a CSV text file
data0 = read.csv("Grid data_adults_R.csv",header=T)
head(data0)
voleinfo = read.csv(file.choose(), header = TRUE)
head(voleinfo)

#load libraries
library(plyr)

##then remove voles capturesd <2x for period 3
captures <- ddply(data0,.(TrapPeriod,
                          Enclosure,
                          Vole_ID), 
                  summarise, Tot=length(Vole_ID))
head(captures); head(data0)

mergecap <- merge(data0, captures, 
                  by = c("TrapPeriod", "Enclosure", "Vole_ID"),
                  all.x= TRUE); head(mergecap)
#data = subset(mergecap, subset = Tot > 1); head(data)  ##this removes <2x/every period
data = mergecap[!(mergecap$TrapPeriod == "Grid3" & mergecap$Tot < 2 ),]


##add sex to data
data = merge(data, voleinfo, by.x = c("Enclosure", "Vole_ID"), 
                                      by.y = c("Enclosure","Vole.ID"),
                                all.x = TRUE)




#####Find Highest space-use overlap w/female using grid data######################

w.data <- ddply(data,.(TrapPeriod,
                       Enclosure,
                       NR_num,
                       Sex, 
                       Vole_ID), 
                       summarise, freq=length(Vole_ID))
head(w.data)

#Drop frequency 
w.data <- data.frame(w.data[1:5]); head(w.data)

#Subset males 
males <- subset(w.data, subset = Sex == "M"); head(males)

#Subset females 
females <- subset(w.data, subset = Sex == "F"); head(females)
                       
#Merge males and females based on trap NR
merge.data <- merge(males, females, by = c("TrapPeriod", "Enclosure", "NR_num"),
                    all= TRUE)
head(merge.data)       
                
#Clean up datafram 

new.data <- data.frame(merge.data$TrapPeriod, 
                       merge.data$Enclosure, 
                       #merge.data$NR_num, 
                       merge.data$Vole_ID.x, 
                       merge.data$Vole_ID.y)                
colnames(new.data)[1:4] <- c("trap.period", 
                             "enclosure", 
                             #"nr.num", 
                             "vole.id.m", 
                             "vole.id.f")
head(new.data)

#build function 
count.duplicates <- function(new.data){
  x <- do.call('paste', c(new.data, sep = '\\r'))
  ox <- order(x)
  rl <- rle(x[ox])
  cbind(new.data[ox[cumsum(rl$lengths)],,drop=FALSE],pair.count = rl$lengths)
}                       
         
#Run function 
dup <- count.duplicates(new.data)

#Return to new data and find the totalnumber of times that males were captured 
gurgle <- ddply(new.data[1:3],.(trap.period,
                       enclosure,
                       vole.id.m), 
                       summarise, freq=length(vole.id.m))
#Drop NA values for males 
w.gurgle <- subset(gurgle, subset = vole.id.m != "<NA>")


#Merge total captures with paired captures
final.data <- merge(dup, w.gurgle, by = c("trap.period", "enclosure","vole.id.m"), 
                    all = TRUE)
head(final.data)

#Do some math
final.data$math <- final.data$pair.count / final.data$freq
head(final.data)

##Overlap of males w/females- this is proportion overlap with all females per male
 
final.data2 <- subset(final.data, vole.id.m != "<NA>")  ##remove male-less captures
head(final.data2)  ###lists males overlap w/every female (or no females)


#now, remove female-less captures UNLESS male shared space w/no F in a period
final.nonas = final.data2[!(final.data2$math < 1 & is.na(final.data2$vole.id.f)) ,]  

      #puts 0 if male overlapped w/no females
final.nonas$math = ifelse(is.na(final.nonas$vole.id.f)  
                     , 0, final.nonas$math) 
head(final.nonas) ###has all overlaps w/females



########Max overlap of each male########
MaxMales = ddply(final.nonas,.(trap.period, enclosure, vole.id.m), 
                               summarise, math = max(math))
head(MaxMales)


###Now to add associated female for each period
  #remove voles w/multiple "max overlap" females
    #merge to add associated females- includes multiples
Max2 = merge(MaxMales, final.nonas, by = c("trap.period", 
            "vole.id.m", "enclosure", "math"),
            all.x=TRUE)
head(Max2)

  # count finds voles w/multiple females listed
Mcounts = ddply(Max2,.(trap.period, enclosure, vole.id.m), 
                 summarise, freq = length(vole.id.m))  
head(Mcounts)
    #add counts to data to remove cases of 2+ females
Mcounts2 = merge(final.nonas, Mcounts, by = 
              c("trap.period",  "enclosure", "vole.id.m"),
             all.x=TRUE)
Mcounts3 = subset(Mcounts2, freq.y == 1)
head(Mcounts3)


###now, adds females ONLY if there is 1 primary female overlaped
Over.final = merge(MaxMales, Mcounts3, by = c("trap.period", "vole.id.m", "enclosure", "math"),
             all.x=TRUE)

Over.final$pair.count = NULL; Over.final$freq.x = NULL; Over.final$freq.y = NULL;
Over.final$maxover = Over.final$math; Over.final$math = NULL;
head(Over.final)



############Determining proportion of females overlapped#########
head(final.data2)
final.data3 = subset(final.data2, vole.id.m != "<NA>")  ##remove male-less captures 
head(final.data3)  ###final.data3 has male-less  captures removed
final.data3
Countfem = ddply(final.data3, .(trap.period, enclosure, ###counts number of females overlapped
                               vole.id.m), 
                      summarize, count = length(unique(
                        vole.id.f[!is.na(vole.id.f)])))
##Check- some males may have 0, NAs should not be counted as a female
head(Countfem)

###total captured females per encl per time period
final.data4 = subset(final.data, vole.id.f != "<NA>")  ##remove female-less captures ONLY
head(final.data4)
Totfem = ddply(final.data4, .(trap.period, enclosure), 
                 summarize, TotF = length(unique(vole.id.f)))
head(Totfem)


###merge data to include total females/encl/time
Merge.fem = merge(Countfem, Totfem, by = c("trap.period", "enclosure"),
                all=TRUE)
head(Merge.fem)
Merge.fem$prop = Merge.fem$count/Merge.fem$TotF
head(Merge.fem)

Merge.fem$count <- NULL; Merge.fem$TotF = NULL  ###void unnecessary columns
head(Merge.fem)   ####Has prop for each male for every time period!!!




#######Combine max overlap & proportion overlap into 1 nice file!######
head(Over.final)
Allgrid = merge(Merge.fem, Over.final, 
          by = c("trap.period", "enclosure","vole.id.m"), 
                all = TRUE)
Allgrid = subset(Allgrid, vole.id.m != "<NA>")
head(Allgrid)
write.csv(Allgrid, file = "AllGridlong.csv")

###reshape into nice wide form for easy viewing
Allcross <- reshape(Allgrid, 
                    timevar = "trap.period",
                    idvar = c("vole.id.m", "enclosure"),
                    direction = "wide")

head(Allcross) ####final product has all variables for each male at all 3 periods!!!!

write.csv(Allcross, file = "Allgridfinal.csv")





                     