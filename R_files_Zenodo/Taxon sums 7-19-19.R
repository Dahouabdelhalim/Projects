######## Taxon sums #####

data<-read.csv("taxon_sum_7-19-19.csv")
head(data)
length(data$taxon) # 44*9 = 396

table(data$taxon) # 44 in each of 9 taxa

class(data$taxon)
class(data$year44)
class(data$non_normalized) # numeric

# summarySE - http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/ 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


### St. dev ###

data$taxon <- factor(data$taxon, levels=c("urchin","seagrass","coral","molluscs","mammals","elasmo","turtle","decapod","fish")) # re-orders


data_summarize1<-summarySE(data,measurevar="year44", groupvars=c("taxon"))
head(data_summarize1)

#### TAXON NORMALIZED ###

library(ggplot2) # http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
ggplot(data_summarize1, aes(x=taxon, y=year44,#fill="red"
)) + 
  labs(title = "Sum of normalized disease reports by taxon")+
  ylab("Sum of normalized disease reports (1970-2013)")+
  xlab("Taxonomic group")+
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=year44-1*sd, ymax=year44+1*sd),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) 


#### TAXON NON-NORMALIZED ###
data$taxon <- factor(data$taxon, levels=c("urchin","seagrass","coral","molluscs","mammals","elasmo","turtle","decapod","fish")) # re-orders


data_summarize1<-summarySE(data,measurevar="non_normalized", groupvars=c("taxon"))


library(ggplot2) # http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
ggplot(data_summarize1, aes(x=taxon, y=non_normalized,#fill="red"
)) + 
  labs(title = "Sum of non-normalized disease reports by taxon")+
  ylab("Sum of non-normalized disease reports (1970-2013)")+
  xlab("Taxonomic group")+
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=non_normalized-1*sd, ymax=non_normalized+1*sd),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))

