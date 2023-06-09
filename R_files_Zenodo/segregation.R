
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("working directory")
seg <- read.csv("Segregation.csv")

head(seg)
str(seg)

# change Chromosome column to factor
seg$Chromosome <- as.factor(seg$Chromosome)


#####---------------------------#####
#####---------------------------#####
# Make the graph
#####---------------------------#####

unique(seg$Chromosome)
colnames(seg)

##### Make new column "pos.cm" w/ continuous positions
#####

seg$pos.cm <- seg$Position_cM

# testing
max(seg[seg$Chromosome == 1, "Position_cM"])
seg[seg$Chromosome == 2, "Position_cM"] + max(seg[seg$Chromosome == 1, "Position_cM"])

#
seg[seg$Chromosome == 2, "pos.cm"] <- seg[seg$Chromosome == 2, "Position_cM"] + max(seg[seg$Chromosome == 1, "Position_cM"])
seg[seg$Chromosome == 3, "pos.cm"] <- seg[seg$Chromosome == 3, "Position_cM"] + max(seg[seg$Chromosome == 2, "pos.cm"])
seg[seg$Chromosome == 4, "pos.cm"] <- seg[seg$Chromosome == 4, "Position_cM"] + max(seg[seg$Chromosome == 3, "pos.cm"])
seg[seg$Chromosome == 5, "pos.cm"] <- seg[seg$Chromosome == 5, "Position_cM"] + max(seg[seg$Chromosome == 4, "pos.cm"])
seg[seg$Chromosome == 6, "pos.cm"] <- seg[seg$Chromosome == 6, "Position_cM"] + max(seg[seg$Chromosome == 5, "pos.cm"])
seg[seg$Chromosome == 7, "pos.cm"] <- seg[seg$Chromosome == 7, "Position_cM"] + max(seg[seg$Chromosome == 6, "pos.cm"])
seg[seg$Chromosome == 8, "pos.cm"] <- seg[seg$Chromosome == 8, "Position_cM"] + max(seg[seg$Chromosome == 7, "pos.cm"])
seg[seg$Chromosome == 9, "pos.cm"] <- seg[seg$Chromosome == 9, "Position_cM"] + max(seg[seg$Chromosome == 8, "pos.cm"])
#####

# add in lines for chromo boundaries
  bounds <- c(max(seg[seg$Chromosome == 1, "pos.cm"]),
              max(seg[seg$Chromosome == 2, "pos.cm"]),
              max(seg[seg$Chromosome == 3, "pos.cm"]),
              max(seg[seg$Chromosome == 4, "pos.cm"]),
              max(seg[seg$Chromosome == 5, "pos.cm"]),
              max(seg[seg$Chromosome == 6, "pos.cm"]),
              max(seg[seg$Chromosome == 7, "pos.cm"]),
              max(seg[seg$Chromosome == 8, "pos.cm"]),
              max(seg[seg$Chromosome == 9, "pos.cm"]))

# graph
ggplot(seg, aes(x = pos.cm, y = f.P., color = Chromosome))+
  geom_point() +
  theme_classic()+
  geom_hline(yintercept = 0.5, lty = 2) +
  geom_vline(xintercept = bounds, col = "grey") +
  ylim(0.3 , 0.65) +
  xlab ("Position (cM)") +
  ylab ("Frequency of P allele") 





#####---------------------------#####
#####---------------------------#####
# Chi-squared test
#####---------------------------#####

# doing a chisq: https://rpubs.com/nmccurtin/chisquare
  # to run the test, input a vector of observed integers (x), and a vector of their expected proportions (p)
  # H0 is that observed values occur in expected proportions
  # HA is that observed values don't occur in expected proportions
  
  chisq.test(x = c(20,45,22), p = c(0.25, 0.5, 0.25)) # p = 0.9
  chisq.test(x = c(40,45,22), p = c(0.25, 0.5, 0.25)) # p = 0.01


# example of running chisq on one row:
  chisq.test(x = seg[1,c("PP","PS","SS")], p = c(0.25, 0.5, 0.25))
  
  
# example of assign test to a dummy variable and extract values  
  abc <- chisq.test(x = seg[1,c("PP","PS","SS")], p = c(0.25, 0.5, 0.25))
  str(abc)
  
  abc$statistic[[1]] # Chi sq value
  abc$parameter[[1]] # df - 2 for all of them, so won't put this in loop.
  abc$p.value[[1]] # p.value
    
  
# Make a summary table to hold stats
chi.summ <- data_frame(marker = seg$Marker, # marker names are in same order in this table as in original csv
                       chi.sq = NA,
                       pvalue = NA)
#data_frame deprecated
chi.summ <- tibble(marker = seg$Marker, # marker names are in same order in this table as in original csv
                       chi.sq = numeric(length=732),
                       pvalue = numeric(length=732))
  
# Make the loop
for( i in 1:nrow(seg)){
  abc <- chisq.test(x = seg[i , c("PP","PS","SS")], p = c(0.25, 0.5, 0.25)) # for row i, do test
  
  chi.summ[i, "chi.sq" ] <- abc$statistic[[1]] # assign chisq value to summary table
  
  chi.summ[i, "pvalue" ] <- abc$p.value[[1]] # assign pvalue to summary table
}


# Join the summary table with the initial table
seg <- full_join(seg, chi.summ, by = c("Marker" = "marker"))


# Make a categorical variable for the pvalues, so can visualize which are significant
 for (i in 1:nrow(seg)){
   ifelse( seg$pvalue[i] < 0.05, seg$orig.p[i] <- "sig", seg$orig.p[i] <- "ns")
   # syntax is 1) a conditional, 2) what to do if conditional is true, 3) what to do if conditional is false
 }


# plot, w/ color for p-values
 ggplot(seg, aes(x = pos.cm, y = f.P., color = orig.p))+
  geom_point() +
  theme_classic()+
  geom_hline(yintercept = 0.5, lty = 2) +
  geom_vline(xintercept = bounds, col = "grey") +
  ylim(0.3 , 0.65) +
  xlab ("Position (cM)") +
  ylab ("Frequency of P allele") 
 
 
#####---------------------------#####
#####---------------------------#####
# Adjusting p-values - tried w/ the strict and adjusted bonferroni... look the same on the graph
#####---------------------------##### 
 
 
# Strict bonferroni - only signifacant if less than 0.5/#tests
 0.05/nrow(seg)
 
 # make a column to say if a marker is significant or not w/ a strict bonferroni
 for (i in 1:nrow(seg)){
   ifelse( seg$pvalue[i] < 0.05/nrow(seg), seg$strict.b[i] <- "sig", seg$strict.b[i] <- "ns")
   # syntax is 1) a conditional, 2) what to do if conditional is true, 3) what to do if conditional is false
 }

  ggplot(seg, aes(x = pos.cm, y = f.P., color = strict.b))+
  geom_point() +
  theme_classic()+
  geom_hline(yintercept = 0.5, lty = 2) +
  geom_vline(xintercept = bounds, col = "grey") +
  ylim(0.3 , 0.65) +
  xlab ("Position (cM)") +
  ylab ("Frequency of P allele") 
  
#Adjusted bonferroni - order your pvalues in ascending order - smallest pvalue is assessed for sig against the most stringent adjusted pvalue (0.05/#tests), then the next smalles pvalue is assessed for sig against the adjusted pvalue of (0.05/(#tests-1)), and so on
  
  # make a vector of all the adjusted p-values to compare against
  bon <- rep(0, nrow(seg))
  
  for(i in 1:nrow(seg)){
    bon[i] <- 0.05/(nrow(seg) - i+1)
  }
  
  # rearrange dataframe so that p-values are in ascending order
  colnames(seg)
  seg <- arrange(seg, pvalue)
  
  # attach on the adjusted p-values
  seg$bon <- bon

  # determine, for each row, if the pvalue is sig after the adjustment
   for (i in 1:nrow(seg)){
   ifelse( seg$pvalue[i] < seg$bon, seg$adj.b[i] <- "sig", seg$adj.b[i] <- "ns")
   # syntax is 1) a conditional, 2) what to do if conditional is true, 3) what to do if conditional is false
   }
  
### Graphs
  
  ggplot(seg, aes(x = pos.cm, y = f.P., color = adj.b))+
  geom_point() +
  theme_classic()+
  geom_hline(yintercept = 0.5, lty = 2) +
  geom_vline(xintercept = bounds, col = "grey") +
  ylim(0.3 , 0.65) +
  xlab ("Position (cM)") +
  ylab ("Frequency of P allele") 

pdf(file="segregation_distortion.pdf", width=7, height=6, pointsize=10)  
  ggplot(seg, aes(x = pos.cm, y = f.P., color = Chromosome, shape = adj.b))+
  geom_point(size = 2, stroke = .75) +
  scale_shape_manual(values = c(16, 24)) + # can find all sorts of different ggplot shapes online
  theme_classic()+
  geom_hline(yintercept = 0.5, lty = 2) +
  geom_vline(xintercept = bounds, col = "grey") +
  ylim(0.45 , 0.65) +
  xlab ("Position (cM)") +
  ylab ("Frequency of P allele")+
    theme(legend.position = "top")
  
  dev.off()

  