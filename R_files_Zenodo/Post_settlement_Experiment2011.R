#Post-mortality oyster experiment conducted in Tomales Bay with outplanted juvenile oysters


data = read.csv("postsettlementM2011.csv")

# subsetting the data for each site, then use glm (binomial error, logistic regression) to test for procedural artifcats within each site.
  
#site 1
site1 = subset(data, site == "1")
site1.model = glm(cbind(live,dead)~treatment, data = site1, family ="binomial")
anova(site1.model, test="LR")

#Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
#NULL                         16     145.44              
#treatment  2   55.931        14      89.51 7.157e-13 ***
 
#get multcomp package
summary(glht(site1.model, mcp(treatment="Tukey"))) 

#Linear Hypotheses:
#  Estimate Std. Error z value Pr(>|z|)    
#cage.control - cage == 0     -0.8394     0.2768  -3.033  0.00683 ** 
#  control - cage == 0          -2.3099     0.3428  -6.738  < 0.001 ***
#  control - cage.control == 0  -1.4705     0.3338  -4.405  < 0.001 ***

# predation present (control differs from cage) and procedural artifact present (cage-control and control differ)
  
#site 2
site2 = subset(data, site == "2")
site2.model = glm(cbind(live,dead)~treatment, data = site2, family ="binomial")
anova(site2.model, test="LR")

#Df Deviance Resid. Df Resid. Dev Pr(>Chi)
#NULL                         17     99.114         
#treatment  2  0.15567        15     98.958   0.9251

summary(glht(site2.model, mcp(treatment="Tukey"))) 

#Linear Hypotheses:
#  Estimate Std. Error z value Pr(>|z|)
#cage.control - cage == 0    -0.06669    0.25827  -0.258    0.964
#control - cage == 0          0.03339    0.25843   0.129    0.991
#control - cage.control == 0  0.10008    0.25836   0.387    0.921

#no predation and no procedural artifact

site3 = subset(data, site == "3")
site3.model = glm(cbind(live,dead)~treatment, data = site3, family ="binomial")
anova(site3.model, test="LR")

#Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
#NULL                         16    107.484              
#treatment  2   34.221        14     73.264 3.707e-08 ***

summary(glht(site3.model, mcp(treatment="Tukey")))


#Linear Hypotheses:
#  Estimate Std. Error z value Pr(>|z|)    
#cage.control - cage == 0     -1.5906     0.3203  -4.966   <1e-05 ***
#  control - cage == 0          -1.4706     0.3127  -4.703   <1e-05 ***
#  control - cage.control == 0   0.1201     0.3468   0.346    0.936 

#predation present without a procedural artifact

site3.2 = subset(data, site == "3.2")
site3.2.model = glm(cbind(live,dead)~treatment, data = site3.2, family ="binomial")
anova(site3.2.model, test="LR")

#Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
#NULL                         17     70.932              
#treatment  2   18.692        15     52.240 8.731e-05 ***

summary(glht(site3.2.model, mcp(treatment="Tukey")))

#Estimate Std. Error z value Pr(>|z|)    
#cage.control - cage == 0     -0.4709     0.2606  -1.807   0.1671    
#control - cage == 0          -1.1400     0.2702  -4.219   <0.001 ***
#  control - cage.control == 0  -0.6690     0.2682  -2.495   0.0337 *  

#predation present and a procedural artifact


site3.5 = subset(data, site == "3.5")
site3.5.model = glm(cbind(live,dead)~treatment, data = site3.5, family ="binomial")
anova(site3.5.model, test="LR")

#Df Deviance Resid. Df Resid. Dev Pr(>Chi)   
#NULL                         17     72.876            
#treatment  2   10.789        15     62.087 0.004541 **

summary(glht(site3.5.model, mcp(treatment="Tukey")))

#Linear Hypotheses:
#  Estimate Std. Error z value Pr(>|z|)   
#cage.control - cage == 0     -0.8358     0.2681  -3.118   0.0051 **
#  control - cage == 0          -0.6132     0.2632  -2.330   0.0517 . 
#control - cage.control == 0   0.2226     0.2728   0.816   0.6930   

#predation present, but no artifact


site3.75 = subset(data, site == "3.75")
site3.75.model = glm(cbind(live,dead)~treatment, data = site3.75, family ="binomial")
anova(site3.75.model, test="LR")

#Df Deviance Resid. Df Resid. Dev Pr(>Chi)
#NULL                         17     138.39         
#treatment  2    1.634        15     136.75   0.4418

#no predation or artifact

site4 = subset(data, site == "4")
site4.model = glm(cbind(live,dead)~treatment, data = site4, family ="binomial")
anova(site4.model, test="LR")

#Df Deviance Resid. Df Resid. Dev Pr(>Chi)   
#NULL                         16    106.235            
#treatment  2   10.761        14     95.473 0.004605 **

summary(glht(site4.model, mcp(treatment="Tukey")))

#Estimate Std. Error z value Pr(>|z|)   
#cage.control - cage == 0      0.6871     0.2816   2.439  0.03900 * 
#  control - cage == 0          -0.1335     0.2938  -0.455  0.89227   
#control - cage.control == 0  -0.8206     0.2715  -3.023  0.00708 **

#no predation but there is a procedural control

# Remove procedural control sites from the data set by removing sites, 1, 3.2, and 4 ####
  
data2 = data[data$site != '1',] 
data3 = data2[data2$site != '3.2',] 
data4 = data3[data3$site != '4',] 
final.data = data4[data4$treatment != 'cage.control',] # remove cage-control because already evaluated procedural artifacts above. Now, focus solely on controls and cages

clean.model = glm(cbind(live,dead)~ treatment*distance, data = final.data, family = "binomial")
anova(clean.model, test="LR")


Analysis of Deviance Table

Model: binomial, link: logit

Response: cbind(live, dead)

Terms added sequentially (first to last)


Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
NULL                                  46     292.19              
treatment           1  11.5794        45     280.61 0.0006669 ***
  distance            1   4.1786        44     276.43 0.0409365 *  
  treatment:distance  1   0.4571        43     275.97 0.4989629     

filtered=ddply(final.data, c("treatment"), summarize,
               n=length(which(!is.na(mort))),
               mean.mort=mean(mort, na.rm=T),
               sd=sd(mort, na.rm=T),
               sterr=sd/sqrt(n))

treatment  n mean.surv        sd      sterr
1      cage 23 0.4847826 0.2781574 0.05799983
2   control 24 0.3750000 0.2236068 0.04564355

filtered
treatment  n mean.mort        sd      sterr
1      cage 23 0.5152174 0.2781574 0.05799983
2   control 24 0.6250000 0.2236068 0.04564355


final.data$treatment=factor(final.data$treatment, levels = c("cage", "control"))
levels(final.data$treatment)[levels(final.data$treatment)=="cage"] <- "Cage"
levels(final.data$treatment)[levels(final.data$treatment)=="control"] <- "Control"

#boxplot for Figure 3A
Treatment.fig <-ggplot(final.data, aes(factor(treatment), y=mort)) + 
  geom_boxplot()+
  #facet_grid( ~ estuary) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  #scale_fill_grey(start = 0, end = 1.0, na.value = "red")+
  
#geom_errorbar(aes(ymin=mean.mort-sterr, ymax=mean.mort+sterr), width=0.5, size=0.5, position=position_dodge(1.0)) +
ylim(0,1)+
  xlab("Experimental treatment") + ylab("Post-settlement mortality\\n (proportional)") +
  guides(fill=FALSE) + # To hide the treatment legend
  theme_classic(base_size = 20) +
  theme(strip.text.x = element_text(face="bold", size=18),
        strip.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size = 18),
        axis.title.y = element_text(face="bold", size=18))


clean.model1a <- glm(cbind(dead,live)~treatment*distance, data = final.data, family ="binomial")

preds <- predict(clean.model1a, type="response", se=T)
preds.df <- data.frame(predictions = preds$fit, 
                       ymin = preds$fit - 2*preds$se.fit,
                       ymax = preds$fit + 2*preds$se.fit,
                       distance=final.data$distance,
                       treatment = final.data$treatment)


Distance.fig<-ggplot(final.data, aes(x=distance, y=mort)) +  
  geom_ribbon(data = preds.df, aes(color=treatment, y = predictions, ymin =
                                     ymin, ymax =
                                     ymax, linetype=NA), fill="grey70",
              alpha = 0.25, show_guide=F) +
  geom_line(data = preds.df, aes(linetype=treatment, y = predictions)) +
  geom_point(aes(shape=treatment, fill=treatment), cex = 5,
             position=position_jitter(width=0.2, height=0), alpha=0.7, color="black") +
  xlab("Distance (km) from ocean") + ylab("Post-settlement mortality\\n(proportional)")  +
  scale_fill_manual(values=c("black", NA)) +
  scale_color_manual(values=c("black", "gray80")) +
  scale_shape_manual(values=c(21, 24)) +
  guides(fill = guide_legend("treatment"), linetype =
           guide_legend("treatment"),
         shape = guide_legend("treatment")) +
  #geom_smooth(method = "lm", colour = "black") + 
  ylim(0,1.25)+
  coord_cartesian()+
  #scale_color_gradient() +
  theme_classic(base_size = 20)  +
  theme(axis.title.x = element_text(face="bold", size = 18),
        axis.title.y = element_text(face="bold", size=18))+
  theme(legend.position="None") 


install.packages("gridExtra")
library(gridExtra)
figure3 = grid.arrange(Treatment.fig, Distance.fig, ncol=1)






