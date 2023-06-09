##########################################################
#   Gawel et al. RSOS Contrasting Ecological Roles      ##
#   of non-native ungulates in a novel ecosystem.       ##
#   Supp Information - justification for using          ##
#   species instead of time in detailed                 ##
#   seedling plot analysis (GLMM)                       ## 
#   supp info on time in ground for seedlings           ##
##########################################################

##load packages###
library(plyr)
library(ggplot2)

##load data
ungulate <- read.csv("seedlingplots_guam.csv", stringsAsFactors=F)

#assign class of each column
factor_cols <- c("island","site", "species", "trt", "deer.abundance", "pig.abundance")
numeric_cols <- c("numplant", "alive", "dead", "length.exposure")
ungulate[factor_cols] <- lapply(ungulate[factor_cols], as.factor)
ungulate[numeric_cols] <- lapply(ungulate[numeric_cols], as.numeric)

#look at data
str(ungulate)
with(ungulate, table(species, length.exposure))

##Does length  of time in the ground affect treatment effects on seedling survival?

unique(ungulate$length.exposure)
hist(ungulate$length.exposure)
with(ungulate, table(species, length.exposure))

ggplot(ungulate,aes(species, length.exposure))+
  geom_boxplot()+
  ylab("No. days in ground")+
  scale_x_discrete("Species", labels=c("A. mariannensis", "C. papaya", "M. citrifolia", "O. oppositifolia", "P. serratifolia", "P. mariana"))+
  theme_minimal()+
  theme(axis.title.y=element_text(size=10, face="bold"),
        axis.title.x  = element_text(size=10,face ="bold"),
        axis.text.x  = element_text(size=10,angle=90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"),
        strip.background = element_rect(colour = "white", fill="white"))

#Planting dates were staggered based on species. Because outplanting dates corresponded to species ready for outplanting, and since there are six species, we see that the time a seedling was in the ground can be grouped into six discrete groups.

#C. papaya: 183 days
#A. mariannensis: 230 - 232 days
#M. citrifolia: 372 - 378 days
#O. oppositifolia: 411 - 418 days
#P. serratifolia: 205 - 208 days
#P. mariana: 110 - 114 days

#Therefore, we group them into six groups by assigning them just one number instead of a range by replacing ranges of values with just one value picked from each group.

ungulate$length.exposure <- gsub("231", "230", ungulate$length.exposure)
ungulate$length.exposure <- gsub("232", "230", ungulate$length.exposure)
ungulate$length.exposure <- gsub("372", "378", ungulate$length.exposure)
ungulate$length.exposure <- gsub("379", "378", ungulate$length.exposure)
ungulate$length.exposure <- gsub("418", "414", ungulate$length.exposure)
ungulate$length.exposure <- gsub("415", "414", ungulate$length.exposure)
ungulate$length.exposure <- gsub("416", "414", ungulate$length.exposure)
ungulate$length.exposure <- gsub("411", "414", ungulate$length.exposure)
ungulate$length.exposure <- gsub("412", "414", ungulate$length.exposure)
ungulate$length.exposure <- gsub("208", "206", ungulate$length.exposure)
ungulate$length.exposure <- gsub("205", "206", ungulate$length.exposure)
ungulate$length.exposure <- gsub("113", "112", ungulate$length.exposure)
ungulate$length.exposure <- gsub("114", "112", ungulate$length.exposure)
ungulate$length.exposure <- gsub("111", "112", ungulate$length.exposure)
ungulate$length.exposure <- gsub("110", "112", ungulate$length.exposure)

unique(ungulate$length.exposure)

ungulate$length.exposure <- factor(ungulate$length.exposure, ordered = TRUE)

##Now, create plot similar to figure 1 in manuscript, except use our discrete lengths of time, corresponding to species, on the y axis, and order by length of time in ground:

newdf<-ddply(ungulate,.(trt, length.exposure),summarise,
             prop=sum(alive)/sum(numplant),
             low=prop.test(sum(alive),sum(numplant))$conf.int[1],
             upper=prop.test(sum(alive),sum(numplant))$conf.int[2])

newdf$length.exposure <- factor(newdf$length.exposure, levels = c("112", "183", "206", "230","378","414"))

#Relabel, matching values to species in x axis:
#C. papaya: 183 days
#A. mariannensis: 230 - 232 days
#M. citrifolia: 372 - 378 days
#O. oppositifolia: 411 - 418 days
#P. serratifolia: 205 - 208 days
#P. mariana: 110 - 114 days

ggplot(newdf,aes(length.exposure,prop, ymin=low, ymax=upper, shape=trt))+
  geom_point(stat="identity", position=position_dodge(width=0.3), size=3)+
  geom_errorbar(width=0.2, position=position_dodge(width=0.3))+
  ylab("Proportion seedling survival")+
  ylim(0,1.1)+
  scale_x_discrete("Species and no. days in ground", labels=c("P. mariana 110-114", "C. papaya ~183", "P. serratifolia 205-208", "A. mariannensis 230-232", "M. citrifolia 372-378","O. oppositifolia 411-418"))+
  scale_shape_manual(values=c(17,19), breaks=c("fenced","ungulate"), labels=c("No ungulates", "Ungulates"))+
  theme_minimal()+
  theme(axis.title.y=element_text(size=10, face="bold"),
        axis.text.x  = element_text(size=10,angle=90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        strip.background = element_rect(colour = "white", fill="white"))

