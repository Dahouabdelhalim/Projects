    #This is a script to import a number of text
    #output files from Gen5 and count colonized wells

    #intended for thrips nectar dispersal
    #assays collected initial and post-treatment

    #written by Marshall McMunn - mmcmunn@gmail.com

#rm(list=ls())
library("plyr")
library("ggplot2")
library("reshape2")

#set work dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("plate reader exports/")

#list folder names
datesExpt<-list.files()

#a function to gather plate reader files (rowwise .txt exports) from a folder, threshold wells, and produce a data frame
#very specific naming of files required with strainID_replicate_date.txt (EC52_3_102820.txt)

gatherPlateFiles<-function(folderName, thresholdSD){
          #move to subfolder
          setwd(paste(folderName ,"/", sep=""))
          
          #within each folder
          #treatment from each file name
          fileNames<-list.files(pattern=".txt")
          
          
          #load all files with read.delim
          d <- lapply(fileNames, function(x) read.delim( x))
          #move back up working directory
          setwd("..")
          
          #get times of all plates, store
          times<-lapply( d, subset, Software.Version=="Time")
          names(times)<-c(fileNames)
          
          #get dates of all plates, store
          dates<-lapply( d, subset, Software.Version=="Date")
          names(dates)<-c(fileNames)
          
          
          #unlist the dataframes and take only rows that have well data
          names(d)<-fileNames
          d1<-do.call(rbind, d)
          #anything with format letter then number
          d1<-d1[ grepl(pattern ="[A-H][0-9]", d1[,1]) ,]
          
          #split the rownames and take out the useful info (1st element of 3)
          d1$sourceFile<-unlist( strsplit(rownames(d1), "\\\\."))[c(TRUE, FALSE, FALSE)]
          d1$sourceFile<-paste(d1$sourceFile , ".txt", sep = "")
          
          #add an index for plate
          colnames(d1) <- c("well" , "absorb600" , "sourceFile")
          newFile<-d1$sourceFile!=c(1,d1$sourceFile[-length(d1$sourceFile)])
          #a running sum of the unique plates from a single trial
          d1$FileCount<-as.factor(cumsum(newFile))
          
          #ID for plate within files
          newPlate<-d1$well=="A1"
          plateNumber<-cumsum(newPlate)

          #match up dates and times from files, plates 
          #(this uses the fact that sequencing is identical between the unlisted lists of dataframes)
          d1$date<-do.call(rbind, dates)[plateNumber ,2]
          d1$time<-do.call(rbind, times)[plateNumber ,2]
          
          #reset rownames to row numbers
          rownames(d1)<-1:nrow(d1)
          
        #this breaks up the file name to retreive microbe, then replicate within trial, then startdate
          #get microbe name
          d1$microbe<-unlist(strsplit(d1$sourceFile , "\\\\_"))[c(TRUE, FALSE, FALSE)]
          #take all numbers from the useful part of the name, this is replicate
          d1$replicate<-unlist(strsplit(d1$sourceFile , "\\\\_"))[c(FALSE, TRUE, FALSE)]
          #get date from file name
          d1$startDate<-unlist(strsplit(d1$sourceFile , "\\\\_|\\\\."))[c(FALSE, FALSE, TRUE, FALSE)]
          
          
          #split well into rows and columns
          d1$wellCol <- gsub("[^[:digit:]]", "", d1$well)
          #take all letters from the useful part, this is treatment
          d1$wellRow <- gsub("[[:digit:]]", "", d1$well)
          
          #this gets rid of wells that were masked during reads (makes NA)
          #there were a few early trials where a thrips fell in and that well was masked.
          #later trials this practice was dropped
          d1$absorb600<-as.numeric(as.character(d1$absorb600))
          
          #make date columns real dates in R
          d1$readDate<-as.Date(d1$date, format = "%m/%d/%Y")
          d1$startDate <-as.Date(d1$startDate, format = "%m%d%y")
          
          #unique wells within trial
          d1$uniqueWellID<-rownames(d1)
          
          #create a column for threhold using
          #deviance from the initial control plates cal C
          for(i in unique(d1$readDate)){
            #for each trial, subset
            temp<-d1[d1$readDate==i , ]
            
            #calculate mean OD on calibration control
            cMean<-mean(temp[temp$replicate=="cal"&temp$microbe=="C" , "absorb600"], na.rm = TRUE)
            
            #calculate sd for calibration control OD
            cSD<-sd(temp[temp$replicate=="cal"&temp$microbe=="C" , "absorb600"], na.rm = TRUE)
            
            #this parameter "thresholdSD" isset in the function call (~line 141)            
            threshold<-cMean+(thresholdSD*cSD)
            #coverts OD to thresholded elevated = 1, or not elevated = 0
            temp$elevatedOD<-ifelse(temp$absorb600>threshold, 1 , 0)
            
            #pastes the output into the trial level dataset "d1" using unique well ID's by trial
            d1[ match(temp$uniqueWellID , d1$uniqueWellID) , "thresholdOD"]<-threshold
            d1[ match(temp$uniqueWellID , d1$uniqueWellID) , "elevatedOD"]<-temp$elevatedOD
          }
          
          #if the plate reader had too high a reading, this was occupied well 
          d1$thresholdOD<-ifelse(d1$absorb600=="OVRFLW" , 1 ,d1$thresholdOD )
          
          #plate ID column
          d1$plateID<-unlist( strsplit(d1$sourceFile, "\\\\."))[c(TRUE, FALSE)]
          
          #calibration logical
          d1$isCal<-grepl(pattern= "[Cc][Aa][Ll]" , d1$plateID)
          
          #control logical
          d1$isControl<-d1$microbe=="C"|d1$microbe=="c"
          d1
}


#an empty list to catch dataframes on loop output
out<-list()
for(m in unique(datesExpt)){
  #applying our plate gather to each subfolder
d1<-gatherPlateFiles(folderName=m , thresholdSD = 6)
  #saving dataframe as list element
out[[m]]<-d1
}


#rbind all dataframes together
dAll<-do.call(rbind , out)
head(dAll)
setwd("..")
write.csv(dAll , "thripsPlatesBulked.csv")

#read in bulked data
#dAll<-read.csv("thripsPlatesBulked.csv")

sum(!dAll$isCal)
#create plate level dataframe
#start with elevated total within each plate
dPlate <- data.frame(total=with(dAll , tapply( elevatedOD,list(microbe, interaction(startDate, replicate)) ,sum, na.rm=TRUE )))
head(dAll)

dPlate <- melt(t(dPlate))
colnames(dPlate)<-c("code", "microbe" , "wellsElevated")

#match rows from the well-level dataframe for some covariates. cbind these rows
dPlate<-cbind(dPlate, 
dAll[
match(
  paste(dPlate$code,dPlate$microbe, sep = ".") ,
gsub("\\\\-", ".", paste("total" , dAll$startDate ,dAll$replicate, dAll$microbe, sep = "."))
),
])

dPlate
#adjusted non-control totals
dPlate$wellsElevated.post.adj<-with(dPlate, ifelse( !isControl&!isCal ,wellsElevated-32 ,wellsElevated  ))

#new column for proportion available wells inoculated
dPlate$propWellsInoc <- ifelse(dPlate$isControl, dPlate$wellsElevated.post.adj/96, dPlate$wellsElevated.post.adj/64)
dPlate<-dPlate[!duplicated(as.list(dPlate))]

dPlate$microbe <- droplevels(dPlate$microbe)
levels(dPlate$microbe)

#get rid of NA rows
dPlate<-na.omit(dPlate)

#assign fungi bacteria as factor
refKeyKingdom<-data.frame(microbe = sort(unique(dPlate$microbe)), kingdom = c("Fungi","Control","Fungi","Bacteria",
                                                               "Fungi","Fungi","Bacteria","Bacteria","Fungi",
                                                               "Bacteria","Bacteria","Bacteria"))

#add kingdom for color
dPlate$kingdom<-refKeyKingdom[match(dPlate$microbe, refKeyKingdom$microbe), "kingdom"]
dPlate$kingdom<-factor(dPlate$kingdom, levels = c("Bacteria", "Fungi","Control"))
sum(!dPlate$isCal)

#plot counts - adjusted only by subtracting 32 for the starting treatment wells
ggplot(data = dPlate[dPlate$isCal==FALSE , ] , aes(x = reorder(microbe, -propWellsInoc ) , y =wellsElevated.post.adj, fill = kingdom)) + 
  geom_bar(position="dodge", stat="summary", fun.y = "mean")+
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.5) + xlab("microbe")+
  ylab("number of wells \\ncolonized") + theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold")) + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
ggtitle("24 hour dispersal via thrips - number wells colonized")+scale_fill_manual(values=c("#56B4E9", "#E69F00","#999999"))+
ggsave("figures/rawCounts.pdf")


head(dPlate)
#stopped here
unique()
dev.off()
#prop counts - adjusted by dividing by the total number of uncolonized wells at outset
ggplot(data = na.omit(dPlate[dPlate$isCal==FALSE , ] ), aes(x = reorder(microbe, -propWellsInoc ) , y =propWellsInoc, fill = kingdom)) + 
  geom_bar(position="dodge", stat="summary", fun.y = "mean")+
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.5) + xlab("microbe")+
  ylab("proportion available wells \\ncolonized") + theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold")) +
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("24 hour dispersal via thrips - proportion wells colonized")+scale_fill_manual(values=c("#56B4E9", "#E69F00","#999999"))
ggsave("figures/propCounts.pdf")

#prop counts - adjusted by dividing by the total number of uncolonized wells at outset
ggplot(data = na.omit(dPlate[dPlate$isCal==FALSE , ] ), aes(x = reorder(microbe, -propWellsInoc , FUN=median ) , y =propWellsInoc, fill = kingdom)) + 
  geom_boxplot()+
  xlab("microbe")+
  ylab("proportion available wells \\ncolonized") + theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold")) +
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  ggtitle("24 hour dispersal via thrips - proportion wells colonized")+scale_fill_manual(values=c("#999999","#56B4E9", "#E69F00"))
ggsave("figures/propCounts.pdf")


#difference in proportion counts beyond control
#create new column (prop control elevated)
dateContam<-data.frame(with(dPlate[dPlate$isControl==TRUE,] , tapply(propWellsInoc,date , mean)))
colnames(dateContam)<-"propWellsInoc"
dateContam$date<-rownames(dateContam)

#contamination rate in controls by date
ggplot(data = dateContam)+aes(x = date, y = propWellsInoc)+geom_point()+
  ggtitle("contamination rate in Controls - proportion wells")

dPlate$trialContamRate<-dateContam[match(dPlate$date , dateContam$date) , "propWellsInoc"]

dPlate$propInocAboveControl<-dPlate$propWellsInoc-dPlate$trialContamRate

table(dPlate[dPlate$isCal==FALSE ,"microbe" ])
ggplot(data = na.omit(dPlate[dPlate$isCal==FALSE , ] ), aes(x = reorder(microbe, -propInocAboveControl ) , y =propInocAboveControl, fill = kingdom)) + 
  geom_bar(position="dodge", stat="summary", fun.y = "mean")+
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.5) + xlab("microbe")+
  ylab("proportion available wells \\ncolonized") + theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold")) +
  theme(axis.text.x = element_text(angle = 60, hjust=1))+scale_fill_manual(values=c("#56B4E9", "#E69F00","#999999"))+
  ggtitle("24 hour dispersal via thrips - proportion wells colonized beyond Control")
ggsave("figures/propCountsAboveControl.pdf")


head(dPlate)

dPlate$wellsAvailable <- ifelse(dPlate$isControl, 96, 64)
dPlate$kingdom<-factor(dPlate$kingdom , levels = c("Control" , "Bacteria" , "Fungi"))

#are bacteria better at dispersing than fungi?
dPlateNO_C<-dPlate[dPlate$kingdom!="Control" , ]
library(lme4)
Mmicrobe<-with(dPlateNO_C , glmer(cbind(wellsElevated.post.adj, 
                      wellsAvailable) ~ kingdom + (1|microbe), 
                family = "binomial" ))


avgKingdoms<-with(dPlate , tapply(propWellsInoc ,list(microbe,kingdom) ,mean ))[,-1]
t.test(avgKingdoms[,1], avgKingdoms[,2])


summary(Mmicrobe)
#yes

#how much better? average microbe strains then average those within kingdom
avgKingdoms<-colMeans(with(dPlate , tapply(propWellsInoc ,list(microbe,kingdom) ,mean )), na.rm=TRUE)
avgKingdoms["Bacteria"]/avgKingdoms["Fungi"] 
#mean 31% higher
medKingdoms<-colMeans(with(dPlate , tapply(propWellsInoc ,list(microbe,kingdom) ,median )), na.rm=TRUE)
medKingdoms["Bacteria"]/medKingdoms["Fungi"]
#median 3.6X higher. wow.
sum(!is.na(dat$Glucose_tot_mgml))

