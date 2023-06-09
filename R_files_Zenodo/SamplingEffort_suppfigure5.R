library(tidyverse)
library(viridisLite)

### 

#Generate supplementary figure 4
#Time course of Damaraland mole-rat activity recording via an RFID reader array

###

###set working directory
setwd('C:/Users/...')
###load data
table=read.csv('DMR_activity_table.csv')
df = table
###labels for durations
labeldays = c('24h-48h','48h-72h','72h-96h','96h-120h','120h-144h')

df$Date <- as.Date(df$Date)

plotdata <- df %>% 
  group_by(ReadingID) %>% 
  mutate(firstread = min(Date), read.days = max(day)) %>% 
  select(Colony, firstread, read.days, colonySize) %>% 
  distinct() 

plotdata <- plotdata %>% 
  group_by(Colony) %>% 
  summarise(groupfirstread = min(firstread)) %>% 
  arrange(groupfirstread) %>% 
  mutate(id = 1:n()) %>% 
  right_join(plotdata)
  
head(plotdata)
labelorder = plotdata %>% distinct(id,.keep_all = T)
labelorder = labelorder[order(labelorder$id),]
grouplabels = labelorder$Colony

suppfigure_4 = ggplot(plotdata, aes(x = firstread, y = id, group = id)) +
  geom_vline(xintercept = seq.Date(as.Date("2014-01-01"),
                                   as.Date("2015-09-01"), 
                                   by = "1 month"), 
             colour = "grey", linetype = "longdash", alpha = 0.6) +
  geom_line(col = "darkgrey") +
  geom_point(aes(size = colonySize, col = as.factor(read.days)), alpha = 0.8) +
  xlab("Date") + 
  ylab("Group Identity") + 
  theme_classic() + 
  theme(axis.title.x = element_text(size = rel(1.3)), 
        axis.title.y = element_text(size = rel(1.3)), 
        axis.text.x = element_text( colour = "black", size = rel(1.05),angle=45, hjust=1),
        axis.text.y = element_text( colour = "black", size = rel(1.05))) + 
  scale_y_continuous(breaks = 1:19, labels = grouplabels) +
  scale_x_date(date_breaks = "3 month", date_labels = "%b %Y") + 
  guides(col = guide_legend("Read Days"), 
         size = guide_legend("Group Size")) +
  scale_colour_viridis_d(labels = labeldays)
suppfigure_4
pdf('SupplementaryFigure5.pdf',width=8,height = 6)
suppfigure_4
dev.off()


