#Career stage vs gender
print(
  D %>% 
  count(stage, gender) %>% 
  mutate(perc = n / nrow(D),
         stage = factor(stage, 
                        levels = c("Research assistant/Researcher/other ECR", "Postdoctoral researcher", "PhD student")))%>%
  ggplot(., aes(y = stage, x = perc, fill=gender)) +
  geom_bar(position="stack", stat = "identity")+
  OYSTERtheme+
  scale_fill_manual(values=OYSTERpalette)+
  # coord_flip() +
  labs(y="Career stage", x="Percentage of respondents (%)", fill="Gender")+
  scale_x_continuous(labels = c("0", "20", "40", "60"))+
  scale_y_discrete(labels=NULL)+
  geom_text(data=function(x) subset(x,stage=="PhD student"), aes(x=.1,label=stage), hjust=-2.6, colour="black", size=7)+
geom_text(data=function(x) subset(x,stage=="Postdoctoral researcher"), aes(x=.1,label=stage), hjust=-1.35, colour="black", size=7)+
geom_text(data=function(x) subset(x,stage=="Research assistant/Researcher/other ECR"), aes(x=.1,label=stage),hjust=-.75,  colour="black",  size=7)
)
