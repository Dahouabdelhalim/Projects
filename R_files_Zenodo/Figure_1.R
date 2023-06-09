###Figure 1- Predicted probs for treatments- early and late succ seedlings and all saplings- based on EMMs and their CIs

##Figure 1a:
#Make EMMs into data frame
surv_early_trt<-as.data.frame(emm_ICC)
surv_earlty_trt

#Figure 1a- early-succ seedlings
seed_dot_early<-ggplot(surv_early_trt,aes(x=Trt,y=prob,fill=Trt, color=Trt))+
  geom_errorbar(aes(ymax=asymp.UCL,ymin=asymp.LCL), width=0,colour="black") +
  geom_point(aes(fill=Trt),size=1)+
  geom_point(shape = 1,size = 1,colour = "black")+
  scale_color_manual(values = c("orange","blue","darkgreen"))+  
  scale_y_continuous()+
  theme_classic()+
  labs(x="Restoration treatment", y = "Predicted survival probability", size = 4)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 10,color='black'))+
  theme(axis.text = element_text(size = 10,color='black'))+
  scale_y_continuous(expand = c(0,0), limits = c(.6,1))

seed_dot_early

##Figure 1b
#Model that includes natural regeneration for proper alignment on plot- remove in post
seedlings_late_C<-seedlings_main%>%
  filter(!succ_stage=="Early")

surv_late_C<-glmer(survival~Trt*year_restsc+htsc+(1|Species),family=binomial, seedlings_late_C)
eml_ICC_C <- emmeans(surv_late_C,"Trt",type='response',bias.adj=T,sigma=reSD_late)

surv_late_trt<-as.data.frame(eml_ICC_C)
surv_late_trt

#Figure 1b- latey-succ seedlings
seed_dot_late<-ggplot(surv_late_trt,aes(x=Trt,y=prob,fill=Trt, color=Trt))+
  geom_errorbar(aes(ymax=asymp.UCL,ymin=asymp.LCL), width=0,colour="black") +
  geom_point(aes(fill=Trt),size=1)+
  geom_point(shape = 1,size = 1,colour = "black")+
  scale_color_manual(values = c("orange","blue","darkgreen"))+  
  scale_y_continuous()+
  theme_classic()+
  labs(x="Restoration treatment", y = "Predicted survival probability", size = 4)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 10,color='black'))+
  theme(axis.text = element_text(size = 10,color='black'))+
  scale_y_continuous(expand = c(0,0), limits = c(.6,1))
seed_dot_late

##Figure 1c- sap survival all recruits
#EEM as data frame
sap_by_trt<-as.data.frame(emm_sap)
sap_by_trt

#Fig 1c
sap_dot_all<-ggplot(sap_by_trt,aes(x=Trt,y=prob,fill=Trt, color=Trt))+
  geom_errorbar(aes(ymax=asymp.UCL,ymin=asymp.LCL), width=0,colour="black") +
  geom_point(aes(fill=Trt),size=1)+
  geom_point(shape = 1,size = 1,colour = "black")+
  scale_color_manual(values = c("orange","blue","darkgreen"))+  
  scale_y_continuous()+
  theme_classic()+
  labs(x="Restoration treatment", y = "Predicted survival probability", size = 4)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 10,color='black'))+
  theme(axis.text = element_text(size = 10,color='black'))+
  scale_y_continuous(expand = c(0,0), limits = c(.6,1))

sap_dot_all

#Combine panels and finish in post
grid.arrange(seed_dot_early,seed_dot_late,sap_dot_all)
```