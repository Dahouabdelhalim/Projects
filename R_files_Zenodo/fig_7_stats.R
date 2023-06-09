require(tidyr)
require(ggplot2)
require(gridExtra)

######
## Recovery
######

reco <- read.csv('Processed_data/Valve_opening_recovery.csv')
reco <- subset(reco, !is.na(reco$Recovery_min))
reco$Interval[is.na(reco$Interval)] <- "NA"
reco$Interval <- as.factor(reco$Interval)
reco$Interval <- factor(reco$Interval, levels = c('NA', '5', '10', '20', '40', '80', '160'))

# we have to analyse the non-responders, the non-recoverers, and the others all separate
reco.sub <- subset(reco, reco$Max_open > .5 & reco$Mean_gape_bef_fst/reco$Max_open > .1)
reco.non.resp <- subset(reco.sub, reco.sub$Recovery_min == 0)
reco.non.rec <- subset(reco.sub, reco.sub$Recovery_min == 70)
reco.recovs <- subset(reco.sub, reco.sub$Recovery_min > 0 & reco.sub$Recovery_min < 70)

reco.cnt <- NULL
for(i in c('NA', '5', '10', '20', '40', '80', '160')){
  reco.cnt.nresp <- data.frame('Interval' = i, 'Type' = 'No closure', 'Count' = length(reco.non.resp$Mussel[reco.non.resp$Interval == i]))
  reco.cnt.rec <- data.frame('Interval' = i, 'Type' = 'Not recovered', 'Count' = length(reco.non.rec$Mussel[reco.non.rec$Interval == i]))
  reco.cnt.recovs <- data.frame('Interval' = i, 'Type' = 'Recovered', 'Count' = length(reco.recovs$Mussel[reco.recovs$Interval == i]))
  
  total.ind <- reco.cnt.nresp$Count + reco.cnt.rec$Count + reco.cnt.recovs$Count
  reco.cnt.nresp$Prop <- reco.cnt.nresp$Count / total.ind
  reco.cnt.rec$Prop <- reco.cnt.rec$Count / total.ind
  reco.cnt.recovs$Prop <- reco.cnt.recovs$Count / total.ind
  
  if(exists('reco.cnt') && is.data.frame(get('reco.cnt'))){
    reco.cnt <- rbind(reco.cnt, reco.cnt.nresp, reco.cnt.rec, reco.cnt.recovs)
  }else{
    reco.cnt <- rbind(reco.cnt.nresp, reco.cnt.rec, reco.cnt.recovs)
  }
}

reco.cnt$Interval <- as.factor(reco.cnt$Interval)
reco.cnt$Interval <- factor(reco.cnt$Interval, levels = c("NA", "5", "10", "20", "40", "80", "160"))
reco.cnt$Type <- as.factor(reco.cnt$Type)
reco.cnt$Type <- factor(reco.cnt$Type, levels = c("No closure", "Recovered", "Not recovered"))


p.types <- ggplot(data = reco.cnt, aes(x = Interval, y = Count, fill = Type))+
  geom_bar(stat="identity", position=position_dodge(), color = 'black')+
  ggtitle('Type of responses after the pulse train onset')+
  ylab("Number of individuals")+xlab("Time interval between pulses (s)")+labs(tag = 'a)')+
  scale_fill_manual(values=c("#C49A00", "#00C094", "#FB61D7"))+
  theme_bw()+theme(legend.position="bottom")
p.types

p.types.prop <- ggplot(data = reco.cnt, aes(x = Interval, y = Prop, fill = Type))+
  geom_bar(stat="identity", position=position_dodge(), color = 'black')+
  ggtitle('Type of responses after pulse train onset')+
  ylab("Proportion of individuals")+xlab("Time interval between pulses (s)")+labs(tag = 'a)')+
  scale_fill_manual(values=c("#C49A00", "#00C094", "#FB61D7"))+
  scale_x_discrete(labels=c(paste("Control\\n(n = ", sum(reco.cnt$Count[reco.cnt$Interval == "NA"]),")", sep = ""), 
                            paste("5\\n(n = ", sum(reco.cnt$Count[reco.cnt$Interval == "5"]),")", sep = ""),
                            paste("10\\n(n = ", sum(reco.cnt$Count[reco.cnt$Interval == "10"]),")", sep = ""),
                            paste("20\\n(n = ", sum(reco.cnt$Count[reco.cnt$Interval == "20"]),")", sep = ""),
                            paste("40\\n(n = ", sum(reco.cnt$Count[reco.cnt$Interval == "40"]),")", sep = ""), 
                            paste("80\\n(n = ", sum(reco.cnt$Count[reco.cnt$Interval == "80"]),")", sep = ""),
                            paste("160\\n(n = ", sum(reco.cnt$Count[reco.cnt$Interval == "160"]),")", sep = "")))+
  theme_bw()+theme(legend.position="bottom")
p.types.prop

# n of this plot
sum(reco.cnt$Count)

# Create trendline for plot B
trendl <- data.frame('Interval' = seq(5, 160, 1))
trendl$Interval_log <- log(trendl$Interval)
trendl$Modelled <- 4.591 + (trendl$Interval_log * 2.5)
trendl$Treatment <- 'exposure'

reco.p <- subset(reco, reco$Mean_gape_bef_fst/reco$Max_open > .1 & reco$Max_open > .5 &
                   reco$Recovery_min != 0 & reco$Recovery_min != 70)
p.recovs <- ggplot(data = reco.p, aes(x = Interval, y = Recovery_min, fill = Treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2), show.legend = FALSE, shape = 1)+
  ggtitle("Recovery time (mussels that recovered)")+
  ylab("Recovery time (min)")+xlab("Time interval between pulses (s)")+labs(tag = 'b)')+
  scale_fill_discrete(labels = c("Control", "Exposure"))+
  scale_x_discrete(labels=c(paste("Control\\n(n = ", reco.cnt$Count[reco.cnt$Interval == "NA" & reco.cnt$Type == "Recovered"],")", sep = ""), 
                            paste("5\\n(n = ", reco.cnt$Count[reco.cnt$Interval == "5" & reco.cnt$Type == "Recovered"],")", sep = ""),
                            paste("10\\n(n = ", reco.cnt$Count[reco.cnt$Interval == "10" & reco.cnt$Type == "Recovered"],")", sep = ""),
                            paste("20\\n(n = ", reco.cnt$Count[reco.cnt$Interval == "20" & reco.cnt$Type == "Recovered"],")", sep = ""),
                            paste("40\\n(n = ", reco.cnt$Count[reco.cnt$Interval == "40" & reco.cnt$Type == "Recovered"],")", sep = ""), 
                            paste("80\\n(n = ", reco.cnt$Count[reco.cnt$Interval == "80" & reco.cnt$Type == "Recovered"],")", sep = ""),
                            paste("160\\n(n = ", reco.cnt$Count[reco.cnt$Interval == "160" & reco.cnt$Type == "Recovered"],")", sep = "")))+
  theme_bw()+theme(legend.position="bottom")+
  geom_text(aes(label = '*', x = 7.3, y = 58), size = 9)
p.recovs

overview <- grid.arrange(p.types.prop, p.recovs, nrow = 1, widths = c(.5, .5))
ggsave('Figures/Fig_7ab_HQ.jpg', overview, png(), units = "in", width = 8.3, height = 4, dpi = 600)

## Stats on the recovery times
reco.stats <- subset(reco, reco$Max_open > .5 & reco$Mean_gape_bef_fst/reco$Max_open > .1 &
                       reco$Recovery_min != 0 & reco$Recovery_min != 70 & reco$Interval != 'NA')
reco.stats$Interval <- as.numeric(as.character(reco.stats$Interval))
hist(reco.stats$Recovery_min)
reco.stats$Interval_log <- log(reco.stats$Interval)
m1 <- glm(Recovery_min ~ Interval_log, data = reco.stats, family = 'gaussian')
m2 <- glm(Recovery_min ~ Interval, data = reco.stats, family = 'gaussian')
aov(m1,m2) # AICc of m1 is better
#plot(m1)
summary(m1)
require(effectsize)
standardize_parameters(m1)