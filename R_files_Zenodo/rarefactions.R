#load required packages
library(plyr)
library(dplyr)
library(reshape)
library(reshape2)
library(tidyr)
library(vegan)
library(data.table)
library(ggplot2)
library(scales)

##############################################
#read in data
##############################################

network <- read.csv("Pollen networks.csv", header=TRUE)

##################################
#data wrangling for metrics----
##################################

#remove insects with no pollen
network$pollen.found[network$pollen.found == 0] <- NA
network <- na.omit(network)
network$Site <- factor(network$Site)
network.melt <- melt(network, id.vars=c(1,3,8), measure.vars=c(14:29))
network.melt$value[network.melt$value == 0] <- NA
network.melt <- na.omit(network.melt)
colnames(network.melt)[3] <- "pollinator"
colnames(network.melt)[4] <- "plant"
trim <- function (x) gsub("^\\\\s+|\\\\s+$", "", x)
network.melt$pollinator <- trim(network.melt$pollinator)
network.melt$plant <- trim(network.melt$plant)
network.melt$pollinator <- gsub(" ", "_",  network.melt$pollinator, fixed = TRUE)
network.melt$plant <- gsub("..", "_",  network.melt$plant, fixed = TRUE)
network.melt$plant <- gsub(".", "_",  network.melt$plant, fixed = TRUE)
network.melt$plant <- sub("_$", "", network.melt$plant)
network.melt$plant<- as.factor(as.character(network.melt$plant))
network.melt$pollinator <- as.factor(as.character(network.melt$pollinator))
network.melt <- as.data.frame(network.melt)
network.melt.2 <- untable(network.melt, num = network.melt$value)
network.melt.2$Site <- factor(network.melt.2$Site)

#remove singletons
network.melt.2<-ddply(network.melt.2,.(Site),function(x){
  if(nrow(x)==1){
    return(NULL)
  }
  else{
    return(x)
  }
})

#create empty list
acc.list <- list()

#run loop over each site
for (j in levels(network.melt.2[, 1])){
  web <- subset(network.melt.2, Site == j)#iterate over site
  id <- as.vector(web$Site[1])
  web <- web %>% mutate(int = group_indices(web))
  web <- tibble::rowid_to_column(web, "ID")
  web$ID <- as.factor(as.character(web$ID))
  web$bp <- paste(web$pollinator, web$plant, sep = "_")
  web$bp <- as.factor(web$bp)
  web <- as.data.frame(web)
  web.2 <- dcast(web, ID ~ bp, fun.aggregate = sum, na.rm =F, value.var="int", fill = 0)
  chao <- specpool(web.2[,-1])
  all <- specaccum(web.2[,-1], method = "random", permutations=10000)
  df <- all[c(3,4,5)]
  df <- data.frame(t(matrix(unlist(df), nrow=length(df), byrow=T)))
  colnames(df)[1] <- "interactions"
  colnames(df)[2] <- "unique"
  colnames(df)[3] <- "sd"
  df$id <- id
  chao <- specpool(web.2[,-1], smallsample = TRUE)
  chao$id <- id
  df <- merge(df,chao)
  df$per.chao <- df$Species/df$chao
  acc.list[[j]] <- df
}

acc.df <- rbind.fill(lapply(acc.list, as.data.frame))
chao.m <- acc.df %>% distinct(id, per.chao, .keep_all = T)
chao.m <- chao.m[c(1,5,6,7,14)]
colnames(chao.m)[1] <- "Site"
colnames(chao.m)[2] <- "Observed unique interactions"
colnames(chao.m)[5] <- "Percent of chao"
write.csv(chao.m, file="C:/Users/jstavert/Documents/Mareeba_Pollen_Network/results/chao-table.csv", row.names=FALSE)

p <- ggplot()
p <- p + xlab("Number of plant-pollinator interaction events") + ylab("Number of unique plant-pollinator interactions")
p <- p + geom_line(data=acc.df, aes(x=interactions, y=unique), size=0.6)
p <- p + geom_ribbon(data=acc.df, aes(ymax=unique+sd, ymin=unique-sd, x=interactions), alpha = 0.35,linetype=0)
p <- p + geom_hline(data = acc.df, aes(yintercept = chao), linetype="longdash", colour="red")
p <- p + geom_hline(data = acc.df, aes(yintercept = chao+chao.se), linetype="dotted")
p <- p + geom_hline(data = acc.df, aes(yintercept = chao-chao.se), linetype="dotted")
p <- p + scale_x_continuous(breaks= pretty_breaks())
p <- p + facet_wrap(~id, scales="free")
p <- p + theme(axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()) +
  theme(axis.text.x=element_text(angle = 360, vjust = 0.5, hjust = 0.5, size =12),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =12),
        axis.title.y=element_text(size=22, vjust = 1),
        axis.title.x=element_text(size=22, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(axis.ticks.length = unit(1.5, "mm"),
        axis.ticks = element_line(colour = 'black', size = 0.4))+
  theme(strip.background = element_rect(colour="NA", fill=NA),
        strip.text = element_text(size=12))
p <- p + theme(legend.key = element_rect(fill = NA, color = NA))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p <- p + theme(axis.title.x=element_text(margin=margin(20,0,0,0)))
p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.4))
p

ggsave(p,file="graphs/sampling-completeness.pdf", device = "pdf",dpi=320,width=14,height=11,units = c("in"))

###########################
#END
###########################