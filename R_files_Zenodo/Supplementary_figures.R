#load data
D_cod <- read_excel("MASTER copy survey2020 covid19.xlsx", 
                    sheet = "Quant recoded data for analysis")


#[DESCRIPTIVE STAT FIGURE FOR Q22]
# Q22	How was communication with your research group affected by Covid-19?
# 1	Much better
# 2	Somewhat better
# 3 Unchanged
# 4 Somewhat worse
# 5 Much worse
# 6 Not applicable

print(
  D_cod%>%
    mutate(Q22_1=as.character(Q22_1))%>%
    count(Q22_1, sort = TRUE)%>%
    mutate(percent=(n/sum(n))*100)%>%
    mutate(Q22_1 = replace(Q22_1, Q22_1 == 1, "Much better"))%>%
    mutate(Q22_1 = replace(Q22_1, Q22_1 == 2, "Somewhat better"))%>%
    mutate(Q22_1 = replace(Q22_1, Q22_1 == 3, "Unchanged"))%>%
    mutate(Q22_1 = replace(Q22_1, Q22_1 == 4, "Somewhat worse"))%>%
    mutate(Q22_1 = replace(Q22_1, Q22_1 == 5, "Much worse"))%>%
    mutate(Q22_1 = replace(Q22_1, Q22_1 == 6, "Not applicable"))%>%
    ggplot(., aes(x="", y=percent, fill=factor(Q22_1, levels=c("Much better", "Somewhat better", "Unchanged", "Somewhat worse", "Much worse", "Not applicable"))))+
    geom_bar(stat = "identity", position="stack")+
    geom_text(aes(label=n), size=5, colour=OYSTERcolor, position = position_stack(vjust = 0.5))+
    labs(y="% Percentage",x="", title="Q22:	How was communication with your research group affected by Covid-19?")+
    scale_fill_manual(values=c("green", "darkolivegreen1", "gray", "#fc9272", "#de2d26", "cornsilk3"))+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_flip()+
    OYSTERtheme+
    theme(legend.title = element_blank(),
          title = element_text(size=9))+
    guides(fill=guide_legend(reverse=TRUE))
)
ggsave(filename = "output/FigS1.png", dpi=300, units="mm", width=170, height = 55)
shell.exec(paste(getwd(),"/output/FigS1.png", sep=""))

#[DESCRIPTIVE STAT FIGURE FOR Q23]
# Q23	Has your relationship with your supervisor changed since Covid-19?
# 1	Much better
# 2	Somewhat better
# 3 Unchanged
# 4 Somewhat worse
# 5 Much worse

print(
  D_cod%>%
    mutate(Q23_1=as.character(Q23_1))%>%
    count(Q23_1, sort = TRUE)%>%
    mutate(percent=(n/sum(n))*100)%>%
    mutate(Q23_1 = replace(Q23_1, Q23_1 == 1, "Much better"))%>%
    mutate(Q23_1 = replace(Q23_1, Q23_1 == 2, "Somewhat better"))%>%
    mutate(Q23_1 = replace(Q23_1, Q23_1 == 3, "Unchanged"))%>%
    mutate(Q23_1 = replace(Q23_1, Q23_1 == 4, "Somewhat worse"))%>%
    mutate(Q23_1 = replace(Q23_1, Q23_1 == 5, "Much worse"))%>%
    ggplot(., aes(x="", y=percent, fill=factor(Q23_1, levels=c("Much better", "Somewhat better", "Unchanged", "Somewhat worse", "Much worse"))))+
    geom_bar(stat = "identity", position="stack")+
    geom_text(aes(label=n), size=5, colour=OYSTERcolor, position = position_stack(vjust = 0.5))+
    labs(y="% Percentage",x="", title="Q23:	Has your relationship with your supervisor changed since Covid-19?")+
    scale_fill_manual(values=c("green", "darkolivegreen1", "gray", "#fc9272", "#de2d26"))+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_flip()+
    OYSTERtheme+
    theme(legend.title = element_blank())+
    guides(fill=guide_legend(reverse=TRUE))
)
ggsave(filename = "output/FigS2.png", dpi=300, units="mm", width=170, height = 50)
shell.exec(paste(getwd(),"/output/FigS2.png", sep=""))

#Q24 & Q25 plot
#What has your supervisor/employer done well/poorly during Covid-19?
#Data in "Q24_25_Supervisors_wellVSbad.xlsx" (local copy)

Supervisor = factor(c(1,2,1,2,1,2,1,2,1,2),labels=c("Well", "Poorly"))
Keyword = factor(c(1,1,2,2,3,3,4,4,5,5),labels=c("Communication","Support","Flexibility", "Management","Availability"))
number_answers = c(148,85,59,41,45,26,41,24,46,41)
df=data.frame(Supervisor,Keyword,number_answers) 


print(
  df%>%group_by(Supervisor)%>%
    mutate(perc=number_answers/sum(number_answers))%>%
    ggplot(., aes(y=Supervisor, x=-perc*100, fill=Keyword)) +
    geom_bar(stat="identity")+
    geom_text(aes(label=number_answers), position = position_stack(vjust = 0.5), size=5, colour=OYSTERcolor)+
    scale_x_continuous(breaks = c(-100, -75, -50, -25, 0), labels= c("0", "25", "50", "75", "100"), expand = c(0,0))+
    scale_y_discrete(expand=c(0,0))+
    scale_fill_brewer(palette="Pastel1")+
    labs(x="% Percentage", y="", title="Q24 & Q25: What has your supervisor/employer done well/poorly during Covid-19?")+
    theme_minimal()+
    OYSTERtheme+
    theme(
          legend.title = element_blank(),
          axis.text.y = element_text(colour = OYSTERcolor, face="bold"),
          plot.title = element_text(size=10, face="bold")))
ggsave(filename = "output/FigS3.png", dpi=300, units="mm", width=170, height = 50)
shell.exec(paste(getwd(),"/output/FigS3.png", sep=""))


#[DESCRIPTIVE STAT FIGURE FOR Q28]
# Q28	Overall, how would you rate the effect of Covid-19 and the related measures on your career in the long term?
# 1	Very negative
# 2	Somewhat negative
# 3 No effect
# 4 Somewhat positive
# 5 Very positive

print(
  D_cod%>%
    mutate(Q28_1=as.character(Q28_1))%>%
    count(Q28_1, sort = TRUE)%>%
    mutate(percent=(n/sum(n))*100)%>%
    mutate(Q28_1 = replace(Q28_1, Q28_1 == 1, "Very negative"))%>%
    mutate(Q28_1 = replace(Q28_1, Q28_1 == 2, "Somewhat negative"))%>%
    mutate(Q28_1 = replace(Q28_1, Q28_1 == 3, "No effect"))%>%
    mutate(Q28_1 = replace(Q28_1, Q28_1 == 4, "Somewhat positive"))%>%
    mutate(Q28_1 = replace(Q28_1, Q28_1 == 5, "Very positive"))%>%
    ggplot(., aes(x="", y=percent, fill=factor(Q28_1, levels=c("Very negative", "Somewhat negative", "No effect", "Somewhat positive", "Very positive"))))+
    geom_bar(stat = "identity", position="stack")+
    geom_text(aes(label=n), size=5, colour=OYSTERcolor, position = position_stack(vjust = 0.5))+
    labs(y="% Percentage",x="", title="Q28:	How would you rate the effect of Covid-19 and the related measures on your career in the long term?")+
    scale_fill_manual(values=c("#de2d26", "#fc9272", "gray", "darkolivegreen1","green"  ))+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_flip()+
    OYSTERtheme+
    theme(legend.title = element_blank(),
          title = element_text(size=7))+
    guides(fill=guide_legend(reverse=TRUE))
)
ggsave(filename = "output/FigS4.png", dpi=300, units="mm", width=170, height = 50)
shell.exec(paste(getwd(),"/output/FigS4.png", sep=""))

#[DESCRIPTIVE STAT FIGURE FOR Q29]
# Q29	Since the Covid-19 crisis, have you considered changing your career path?
# 0	It has made me more dedicated to my career
# 1	I don't want to change career
# 2 I have not thought about it 
# 3 I would like to change career
# 4 I have taken steps towards changing my career

print(
  D%>%
    mutate(Q29=as.character(Q29))%>%
    count(Q29, sort = TRUE)%>%
    mutate(percent=(n/sum(n))*100)%>%
    mutate(Q29 = replace(Q29, Q29 == 0, "It has made me more dedicated to my career"))%>%
    mutate(Q29 = replace(Q29, Q29 == 1, "I don't want to change career"))%>%
    mutate(Q29 = replace(Q29, Q29 == 2, "I have not thought about it"))%>%
    mutate(Q29 = replace(Q29, Q29 == 3, "I would like to change career"))%>%
    mutate(Q29 = replace(Q29, Q29 == 4, "I have taken steps towards changing my career"))%>%
    ggplot(., aes(x="", y=percent, fill=factor(Q29, levels=c("It has made me more dedicated to my career", "I don't want to change career", "I have not thought about it", "I would like to change career", "I have taken steps towards changing my career"))))+
    geom_bar(stat = "identity", position="stack")+
    geom_text(aes(label=n), size=4, colour=OYSTERcolor, position = position_stack(vjust = 0.5))+
    labs(y="% Percentage",x="", title="Q29:	Since the COVID-19 crisis, have you considered changing your career path?")+
    scale_fill_brewer(palette="Pastel1")+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_flip()+
    OYSTERtheme+
    theme(legend.title = element_blank(),
          title = element_text(size=9))+
    guides(fill=guide_legend(reverse=TRUE))
)
ggsave(filename = "output/FigS5.png", dpi=300, units="mm", width=170, height = 50)
shell.exec(paste(getwd(),"/output/FigS5.png", sep=""))

#[DESCRIPTIVE STAT FIGURE FOR Q27]
# Q27	Do you think the way we do marine research will be affected by COVID-19 over the next 10 years?
# 0	Yes
# 1	No

print(
  D%>%
    mutate(Q27=as.character(Q27))%>%
    count(Q27, sort = TRUE)%>%
    mutate(percent=(n/sum(n))*100)%>%
    mutate(Q27 = replace(Q27, Q27 == 0, "Yes"))%>%
    mutate(Q27 = replace(Q27, Q27 == 1, "No"))%>%
    ggplot(., aes(x="", y=percent, fill=factor(Q27, levels=c("Yes", "No"))))+
    geom_bar(stat = "identity", position="stack")+
    geom_text(aes(label=n), size=6, colour=OYSTERcolor, position = position_stack(vjust = 0.5))+
    labs(y="% Percentage",x="", title="Q27: Do you think the way we do marine research will be affected by COVID-19 over the next 10 years?")+
    scale_fill_manual(values=c("#fed9a6", "#b3cde3"))+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_flip()+
    OYSTERtheme+
    theme(legend.title = element_blank(),
          title = element_text(size=7))+
    guides(fill=guide_legend(reverse=TRUE))
)
ggsave(filename = "output/FigS6.png", dpi=300, units="mm", width=170, height = 40)
shell.exec(paste(getwd(),"/output/FigS6.png", sep=""))

#[DESCRIPTIVE STAT FIGURE FOR Q26]
# Q26	In a post Covid-19 world, do you forsee permanent changes to your day-to-day research practices?
# 0	Yes
# 1	No

print(
  D%>%
    mutate(Q26=as.character(Q26))%>%
    count(Q26, sort = TRUE)%>%
    mutate(percent=(n/sum(n))*100)%>%
    mutate(Q26 = replace(Q26, Q26 == 0, "Yes"))%>%
    mutate(Q26 = replace(Q26, Q26 == 1, "No"))%>%
    ggplot(., aes(x="", y=percent, fill=factor(Q26, levels=c("Yes", "No"))))+
    geom_bar(stat = "identity", position="stack")+
    geom_text(aes(label=n), size=7, colour=OYSTERcolor, position = position_stack(vjust = 0.5))+
    labs(y="% Percentage",x="", title="Q26:	In a post COVID-19 world, do you foresee permanent changes to your day-to-day research practices?")+
    scale_fill_manual(values=c("#fed9a6", "#b3cde3"))+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_flip()+
    OYSTERtheme+
    theme(legend.title = element_blank(),
          title = element_text(size=7))+
    guides(fill=guide_legend(reverse=TRUE))
  )

ggsave(filename = "output/FigS7.png", dpi=300, units="mm", width=170, height = 40)
shell.exec(paste(getwd(),"/output/FigS7.png", sep=""))

#[DESCRIPTIVE STAT FIGURE FOR Q16]
# Q16	Has your funding been extended?
# 0	Yes
# 1	No
# 2 Don't know

print(
  D%>%
    mutate(Q16=as.character(Q16))%>%
    count(Q16, sort = TRUE)%>%
    mutate(percent=(n/sum(n))*100)%>%
    mutate(Q16 = replace(Q16, Q16 == 0, "Yes"))%>%
    mutate(Q16 = replace(Q16, Q16 == 1, "No"))%>%
    mutate(Q16 = replace(Q16, Q16 == 2, "Don't know"))%>%
    ggplot(., aes(x="", y=percent, fill=factor(Q16, levels=c("Yes", "No", "Don't know"))))+
    geom_bar(stat = "identity", position="stack")+
    geom_text(aes(label=n), size=7, colour=OYSTERcolor, position = position_stack(vjust = 0.5))+
    labs(y="% Percentage",x="", title="Q16: Has your funding been extended?")+
    scale_fill_manual(values=c("#fed9a6", "#b3cde3", "gray"))+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_flip()+
    OYSTERtheme+
    theme(legend.title = element_blank())+
    guides(fill=guide_legend(reverse=TRUE))
)

ggsave(filename = "output/FigS8.png", dpi=300, units="mm", width=170, height = 40)
shell.exec(paste(getwd(),"/output/FigS8.png", sep=""))


#[DESCRIPTIVE STAT FIGURE FOR Q20]
# Q20	Would you prefer to work from home more often than you did before Covid-19 measures were introduced?
# 0	Yes, sometimes
# 1	Yes, part time
# 2 Yes, full time
# 3 No

print(
  D%>%
    mutate(Q20=as.character(Q20))%>%
    count(Q20, sort = TRUE)%>%
    mutate(percent=(n/sum(n))*100)%>%
    mutate(Q20 = replace(Q20, Q20 == 0, "Yes, sometimes"))%>%
    mutate(Q20 = replace(Q20, Q20 == 1, "Yes, part time"))%>%
    mutate(Q20 = replace(Q20, Q20 == 2, "Yes, full time"))%>%
    mutate(Q20 = replace(Q20, Q20 == 3, "No"))%>%
    ggplot(., aes(x="", y=percent, fill=factor(Q20, levels=c("Yes, full time", "Yes, part time", "Yes, sometimes", "No"))))+
    geom_bar(stat = "identity", position="stack")+
    geom_text(aes(label=n), size=5, colour=OYSTERcolor, position = position_stack(vjust = 0.5))+
    labs(y="% Percentage",x="", title="Q20: Would you prefer to work from home more often than you did before COVID-19?")+
    scale_fill_brewer(palette="Pastel1")+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_flip()+
    OYSTERtheme+
    theme(legend.title = element_blank(),
          title=element_text(size=8))+
    guides(fill=guide_legend(reverse=TRUE))
)
ggsave(filename = "output/FigS9.png", dpi=300, units="mm", width=170, height = 40)
shell.exec(paste(getwd(),"/output/FigS9.png", sep=""))


#[DESCRIPTIVE STAT FIGURE FOR Q11] - FigS10
# Q11	Did you have any care duties during the Covid-19 pandemic?
# 0	Yes
# 1	No

print(
  D%>%
  mutate(Q11=as.character(Q11))%>%
  count(Q11, sort = TRUE)%>%
  mutate(percent=(n/sum(n))*100)%>%
  mutate(Q11 = replace(Q11, Q11 == 0, "Yes"))%>%
  mutate(Q11 = replace(Q11, Q11 == 1, "No"))%>%
  ggplot(., aes(x="", y=percent, fill=factor(Q11, levels=c("Yes", "No"))))+
  geom_bar(stat = "identity", position="stack")+
  geom_text(aes(label=n), size=7, colour=OYSTERcolor, position = position_stack(vjust = 0.5))+
  labs(y="% Percentage",x="", title="Q11: Did you have any care duties during the Covid-19 pandemic?")+
  scale_fill_manual(values=c("#fed9a6", "#b3cde3"))+
  scale_x_discrete(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  coord_flip()+
  OYSTERtheme+
  theme(legend.title = element_blank())+
    guides(fill=guide_legend(reverse=TRUE))
)

ggsave(filename = "output/FigS10.png", dpi=300, units="mm", width=170, height = 45)
shell.exec(paste(getwd(),"/output/FigS10.png", sep=""))





