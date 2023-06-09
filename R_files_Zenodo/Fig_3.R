##Figure 3- DBH by Trt

#make EMM a data frame
DBH_trt_diff<-as.data.frame(emm_DBHe_late) #from EMM by Trt levels
DBH_trt_diff

DBH_trt_plot<-ggplot(DBH_trt_diff,aes(x=Trt,y=response,fill=Trt, color=Trt))+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL), width=0,colour="black") +
  geom_point(aes(fill=Trt),size=1.5)+
  geom_point(shape = 1,size = 1.5,colour = "black")+
  scale_color_manual(values = c("blue","darkgreen"))+  
  scale_y_continuous()+
  theme_classic()+
  labs(x="Restoration treatment", y = 'DBH (cm yr-1)', size = 4)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 10,color='black'))+
  theme(axis.text = element_text(size = 10,color='black'))+
  scale_y_continuous(expand = c(0,0), limits = c(0,.45))

DBH_trt_plot
