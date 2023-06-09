library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

pubs<-read.csv('Fig1-CEP_keywords_PubsPerYr.csv', sep=",", header=T)

# graph cumulative publications per year, by category
pubs.cumsum<-pubs %>%
  group_by(as.factor(Search)) %>%
  arrange(Publication_year) %>%
  mutate(cumsum=cumsum(Record_count))

pubs.cumsum<-as.data.frame(pubs.cumsum)
pubs.cumsum$label<-rep("", nrow(pubs.cumsum))
pubs.cumsum[pubs.cumsum$Publication_year==max(pubs.cumsum$Publication_year), "label"]<-pubs.cumsum[pubs.cumsum$Publication_year==max(pubs.cumsum$Publication_year), "cumsum"]

ggplot(pubs.cumsum, aes(x=Publication_year, y=cumsum, color=as.factor(Search)))+
  geom_step()+
  theme_bw()+
  labs(x="Year", y="Cumulative record count", colour="Search terms")+
  scale_color_manual(values=c("#b35806", "#d8daeb", "#b2abd2", "#8073ac", "#542788"))+
  theme(axis.text=element_text(size=11),
  	axis.title=element_text(size=11),
  	legend.position=c(.25, .65),
  	legend.background=element_rect(fill="white", colour="black", size=0.25))+
  	geom_text_repel(data=pubs.cumsum, aes(label=label), nudge_y=100, nudge_x=2, na.rm=T, show.legend=F)+
  	xlim(1965,2025)+
  	ylim(0,5250)
  
ggsave("pubs-per-year-cumulative.png", height=3, width=3.15, dpi=600)