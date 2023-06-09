exx = read.csv("/2021/data/data_with_plumx.csv")
attach(exx)
class(exx)
#cor.test(citation_count,Social_count,method = "spearm",exact=FALSE)
#cor.test(citation_count,Mention_count,method = "kendall")
model <- lm(citation_count ~ Tweet_count, data = exx)
model
summary(model)


##-------------scatter plot

exx = read.csv("/2021/data/data_with_plumx.csv")
library(reshape2)
library(ggplot2)
yt<-ggplot(data = extractedPapers, aes(x = citation_count, y =Tweet_count )) +
  geom_point(size=1)+ 
  labs(y="Tweets", x = "Citations")+
  scale_x_continuous(limits = c(1, 500))+
  scale_y_continuous(limits = c(1, 500))+
  theme_bw()+
  theme(axis.title = element_text(size=12),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=10),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.25, 0.80),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'),  
        legend.key.height = unit(0.4, 'cm'),  
        legend.key.width = unit(0.5, 'cm'),
        legend.text=element_text(size=7),
        panel.border = element_rect(size=0.2, linetype="solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_line(size=0.1)
  ) 
plot(yt)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE) +facet_wrap(~ Year)

#----------------------
#-- Twitter paper barplot
exx = read.csv("/Users/randalchkr/Dropbox/PHD/work/Yusra/papergraphs/2021/data/groupeddataT.csv")
library(ggplot2)
plot1<-ggplot(exx, aes(x=cat2, y=value,fill=cate)) + 
  geom_bar(stat="identity", position=position_dodge(),width = 0.8)+
  labs(y="Count", x = "Years")+ 
  scale_x_continuous(breaks = c(1:7),labels=c("2015", "2016", "2017","2018","2019","2020","2021"))+
  scale_y_continuous(limits = c(0, 210000),,breaks=seq(0,200000,25000))+
  theme_minimal()+
  theme( 
    panel.grid.major.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.x =  element_line(size=0.1), 
    panel.grid.minor.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white",color="black",size=0.2),
    axis.ticks.y =  element_line(size=0.1), 
    axis.title = element_text(size=10),
    axis.text.y= element_text(size=8),
    axis.text.x= element_text(size=8), 
    legend.position ="top",
    legend.key.size = unit(0.5, 'cm'), 
    legend.key.height = unit(0.4, 'cm'), 
    legend.key.width = unit(0.5, 'cm'),
    legend.text=element_text(size=8),
    legend.title = element_blank()
  ) +scale_fill_manual(values=c( "#41b6c4","#a1dab4","#0868ac","grey", "#253494","#e0e000"),
                       breaks=c("a", "b", "c","d","e","f"),
                       labels=c("Captures", "Usage counts", "Citations","Social media counts","Tweets","Mentions")
  )
plot1

library(reshape2)
library(ggplot2)
yt<-ggplot(data = groupeddataT, aes(x = citation_count, y =Tweet_count )) +
  geom_point(size=1)+ 
  labs(y="Tweets", x = "Citations")+
  scale_x_continuous(limits = c(1, 500))+
  scale_y_continuous(limits = c(1, 500))+
  theme_bw()+
  theme(axis.title = element_text(size=12),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=10),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.25, 0.80),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'),  
        legend.key.height = unit(0.4, 'cm'),  
        legend.key.width = unit(0.5, 'cm'),
        legend.text=element_text(size=7),
        panel.border = element_rect(size=0.2, linetype="solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_line(size=0.1)
  ) 
plot(yt)+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE) +facet_wrap(~ Year)
