#Script: percentiles.R
#Author: Diego Chavarro (dchavarro@gmail.com), Based on the paper Do Synthesis Centers Synthesize? Research Policy, forthcoming.
#Version: 0.3
#Note: this script is only provided as a guidance that needs to be adapted to your specific computing environment. It computes median and top 10% diversity adapted from Uzzi's paper: Uzzi, B., Mukherjee, S., Stringer, M., & Jones, B. (2013). Atypical Combinations and Scientific Impact. Science, 342(6157), 468-472. https://doi.org/10.1126/science.1240474





#require library tidyr (needed for this script to work)
 library(tidyr)

#set work-space
 setwd("PATH\\\\TO\\\\YOUR-WORKSPACE")

#Load topic distance matrix
 distData <-read.table("PATH\\\\TO\\\\distMatrix1.csv", header=TRUE, sep=",")
#convert distance matrix into long format for calculations 
 distDataLong <- gather(distData, topic, weight, T1:T152)
#apply column names distance matrix. First column is topic 1, second column topic 2, and third column is the distance between the 2.
 colnames(distDataLong) <- c("topic1","topic2","distance")
#remove the "T" and leave only the numbers of the topics 1 to 152
 distDataLong$topic1 <- gsub("T", "", distDataLong$topic1)
 distDataLong$topic2 <- gsub("T", "", distDataLong$topic2)
#ensure cells are numeric
 distDataLong$topic2 <- as.numeric(distDataLong$topic2)
 distDataLong$topic1 <- as.numeric(distDataLong$topic1)
#sort out matrix
 distDataLong <- distDataLong[order(distDataLong$topic2,distDataLong$topic1),]
#read papers database. In the repository there is one file per batch of approx. 50,000 papers. change file to analyze a different batch, or implement a loop.
 art <- read.table("PATH#\\\\TO\\\\Table 3 .csv", header=TRUE, sep=",", row.names="id")
#calculate proportions of each paper each of the 152 topics
 prop <- art / rowSums(art, na.rm=TRUE)
#count number of papers in batch
 count <- nrow(prop)

#calculations for each paper in batch
	for(i in 1:count){
		#paperTopicMatrix combinations
		 PaperTopicX <- as.numeric(prop[i,]) %*% t(as.numeric(prop[i,])) 
		#convert to data frame
		 PaperTopicX <- as.data.frame(PaperTopicX)
		#add topic numbers to columns
		 PaperTopicX <- cbind(X = seq(1:152),PaperTopicX)
		#convert to Long format
		 PaperTopicXLong <- gather(PaperTopicX,topic, weight, V1:V152)
		#name the columns of the long format table for papers and topics
		 colnames(PaperTopicXLong) <- c("topic1","topic2","weight")
		#remove V from variable names, so that the topics match those of the distance matrix
		 PaperTopicXLong$topic2 <- gsub("V", "", distDataLong$topic2)
		#ensure cells are numeric
		 PaperTopicXLong$topic2 <- as.numeric(PaperTopicXLong$topic2)
		 PaperTopicXLong$topic1 <- as.numeric(PaperTopicXLong$topic1)
		#merge topic distance and paper distribution of topics
		 paper <- merge(PaperTopicXLong,distDataLong,by=c("topic1","topic2"))
		#add identifiers to papers
		 paper <- cbind(Id = rownames(prop[i,]),paper)
		#sort papers 
		 OrderedPaper <- paper[order(paper$distance,paper$topic1),]
		#calculate the cumulative sum of sum of each paper's weights
		 OrderedPaper <- cbind(OrderedPaper, cumul = cumsum(OrderedPaper$weight))
		#find the values at which the cumulative sum reaches 0.5
		 paperPercentil50 <-head (OrderedPaper[OrderedPaper$cumul >=0.5,],1)
		#find the values at which the cumulative sum reaches 0.9
		 paperPercentil90 <-head (OrderedPaper[OrderedPaper$cumul >=0.9,],1)
		#merge results
		 paperPercentil <- cbind(paperPercentil50, paperPercentil90)
		#write results
		 write.table(paperPercentil, "percentil.txt", append=TRUE, quote = TRUE, sep = "\\t", na = "0", row.names=TRUE, col.names=FALSE)

		

	}
