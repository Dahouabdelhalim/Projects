#load the required libraries
library(ggplot2)
library(overlapping)

#read in IBS values for within- and among-population comparisons
IBS_invasive_within_pop <- read.csv(file="IBS_invasive_within_population.txt",header = FALSE, sep = ,)
IBS_invasive_among_pop <- read.csv(file="IBS_invasive_among_populations.txt",header = FALSE, sep = ,)

IBS_native_within_pop <- read.csv(file="IBS_native_within_population.txt",header = FALSE, sep = ,)
IBS_native_among_pop <- read.csv(file="IBS_native_among_populations.txt",header = FALSE, sep = ,)

#plot IBS values for the native range
ggplot(data = IBS_native_among_pop, aes(x=V1)) +
  geom_line(stat="density")+
  geom_line(stat="density",data=IBS_native_within_pop, aes(x=V1), linetype="longdash")+
  xlab("\\nPairwise IBS values")+
  ylab("Count\\n")+
  xlim(0.85,1)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#plot IBS values for the invasive range
ggplot(data = IBS_invasive_among_pop, aes(x=V1)) +
  geom_density(color=NA, fill="grey90")+
  geom_density(color=NA,data=IBS_invasive_within_pop, aes(x=V1),fill="grey90")+
  geom_line(stat="density")+
  geom_line(stat="density",data=IBS_invasive_within_pop, aes(x=V1), linetype="longdash")+
  xlab("\\nPairwise IBS values")+
  ylab("Count\\n")+
  xlim(0.85,1)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


#estimate overlap in the two distributions
native_within<-IBS_native_within_pop$V1
native_among<-IBS_native_among_pop$V1

native <- list(X1=native_within, X2=native_among)
out <- overlap(native, plot=TRUE)

invasive_within<-IBS_invasive_within_pop$V1
invasive_among<-IBS_invasive_among_pop$V1

invasive <- list(X1=invasive_within, X2=invasive_among)
out <- overlap(invasive, plot=TRUE)
