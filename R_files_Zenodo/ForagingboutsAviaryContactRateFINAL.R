##RFID data analysis of foraging bout lengths for Hawley et al. High virulence is associated with pathogen spreadability in a songbird-bacterial system  
rm(list=ls())

##first import the RFID bout data 
library(readr)
feeding_bouts_4sRule_2017 <- read_csv("Google Drive/My Drive/Rcode/DanaHawley/feeding_bouts_4sRule_2017.csv")
View(feeding_bouts_4sRule_2017)

allbouts<-feeding_bouts_4sRule_2017 
head(allbouts)
dim(allbouts) #157496 9
names(allbouts)

##subset bouts to those that are at least 3 seconds in length (as defined in manuscript)
allboutsmin=subset(allbouts, Bout_duration_s > 3)


##Match the PIT tags to bird ID the crude way
index <- c('1129', '1136', '1147', '1153', '1162', '1168', '1169', '1170', '1171')
index_PIT <- c('0700EDFF14', '0700EDD671', '0700EE1002', '0700EE1D5A', '0700EE16FC',
               '0700EDB1D6', '0700EDABD6', '0700EE0656', '0700EDBC2A')


#match flock treatment to PIT
#not the most efficient way to do this, but what worked the quickest

#treatment
CON <- c('0700EDED3D','0700EE1002','0700EDBB52','0700ED8E73','0700EDAB2E',
         '0700EE3293','0700EE2EA4','0700EDC42F','0700EDB1D6','0700EE1974',
         '0700EDBED7','0700EDBDFD','0700EDABD6','0700EDA5A5','0700EE34D3')
LOW <- c('0700EE3042','0700ED9E69','0700EE2AB0','0700EDA26D','0700EE0656',
         '0700ED8748','0700EDB342','0700EDA95A','0700EE352B','0700EDBC2A',
         '0700EDFF14','0700EDA993','0700EE1BA4','0700ED9E70','0700ED99D2') #1165, 01101775AE = 0700ED99D2
HIGH <- c('0700ED9462','0700EDBAB1','0110176895','0700EE2081','0700EE16FC',
          '0700ED89B9','0700EDD671','0700EE216A','0700ED9A67','0110176B0D',
          '0700ED8B7A','0700EDA2B5','0700EE2810','0700EE1D5A','0700EDD4BC')

#NOTE: double checked, 1165 was listed as PIT 01101775AE and as 0700ED99D2. 2nd PIT 0700ED99D2 is correct.
# Exactly 4854 rows of data missing when using 01101 PIT and exactly 4854 rows of data added back in when
# using 0700ED PIT instead, indicating that 0700ED is the correct PIT Tag to use for bird 1165.

conbouts <- allboutsmin[allboutsmin$PIT_Tag %in% CON,]
conbouts$Treatment <- 'Control'

lowbouts <- allboutsmin[allboutsmin$PIT_Tag %in% LOW,]
lowbouts$Treatment <- 'Low'

highbouts <- allboutsmin[allboutsmin$PIT_Tag %in% HIGH,]
highbouts$Treatment <- 'High'

#groups
CON1 <- c('0700EDED3D','0700EE1002','0700EDBB52','0700ED8E73','0700EDAB2E')
CON2 <- c('0700EE3293','0700EE2EA4','0700EDC42F','0700EDB1D6','0700EE1974')
CON3 <- c('0700EDBED7','0700EDBDFD','0700EDABD6','0700EDA5A5','0700EE34D3')
LOW1 <- c('0700EE3042','0700ED9E69','0700EE2AB0','0700EDA26D','0700EE0656')
LOW2 <- c('0700ED8748','0700EDB342','0700EDA95A','0700EE352B','0700EDBC2A')
LOW3 <- c('0700EDFF14','0700EDA993','0700EE1BA4','0700ED9E70','0700ED99D2') #1165, 01101775AE = 0700ED99D2
HIGH1 <- c('0700ED9462','0700EDBAB1','0110176895','0700EE2081','0700EE16FC')
HIGH2 <- c('0700ED89B9','0700EDD671','0700EE216A','0700ED9A67','0110176B0D')
HIGH3 <- c('0700ED8B7A','0700EDA2B5','0700EE2810','0700EE1D5A','0700EDD4BC')

con1b <- conbouts[conbouts$PIT_Tag %in% CON1,]
con1b$Group <- 'CON1'
con2b <- conbouts[conbouts$PIT_Tag %in% CON2,]
con2b$Group <- 'CON2'
con3b <- conbouts[conbouts$PIT_Tag %in% CON3,]
con3b$Group <- 'CON3'

conbouts2 <- rbind(con1b, con2b, con3b)

low1b <- lowbouts[lowbouts$PIT_Tag %in% LOW1,]
low1b$Group <- 'LOW1'
low2b <- lowbouts[lowbouts$PIT_Tag %in% LOW2,]
low2b$Group <- 'LOW2'
low3b <- lowbouts[lowbouts$PIT_Tag %in% LOW3,]
low3b$Group <- 'LOW3'

lowbouts2 <- rbind(low1b, low2b, low3b)

high1b <- highbouts[highbouts$PIT_Tag %in% HIGH1,]
high1b$Group <- 'HIGH1'
high2b <- highbouts[highbouts$PIT_Tag %in% HIGH2,]
high2b$Group <- 'HIGH2'
high3b <- highbouts[highbouts$PIT_Tag %in% HIGH3,]
high3b$Group <- 'HIGH3'

highbouts2 <- rbind(high1b, high2b, high3b)

allbouts2 <- rbind(conbouts2, lowbouts2, highbouts2)

library(doBy)

##make a dataframe for each categorical time interval; pre-inoculation is 10/8 only; other categories only include days that birds were not caught and sampled
pre_inf<- c('10/08/17')
early_inf<-c('10/10/17', '10/11/17')
peak_inf<- c('10/17/17', '10/18/17', '10/21/17')

##link dates to category 
pre <- allbouts2[allbouts2$Date %in% pre_inf,]
pre$inf_period <- 'PRE'
early<- allbouts2[allbouts2$Date %in% early_inf,]
early$inf_period <- 'EARLY'
peak <- allbouts2[allbouts2$Date %in% peak_inf,]
peak$inf_period <- 'PEAK'


allbouts3 <- rbind(pre, early, peak)


##distinguish between index birds and flockmates
allbouts3$Status <- ifelse(allbouts3$PIT_Tag %in% index_PIT, 'index', 'flockmate')


##Relevel bout data (it's not treating Treatment as a "factor" first make it a factor
allbouts3$Treatment = as.factor(allbouts3$Treatment)
levels(allbouts3$Treatment)

allbouts3$inf_period = as.factor(allbouts3$inf_period)
levels(allbouts3$inf_period)

library(tidyverse)
allbouts3 <- mutate(allbouts3, 
                    Treatment=factor (Treatment, 
                                      levels=c("Control", "Low", "High")))
levels(allbouts3$Treatment)

library(tidyverse)
allbouts3 <- mutate(allbouts3, 
                    inf_period=factor (inf_period, 
                                       levels=c("PRE", "EARLY", "PEAK")))
levels(allbouts3$inf_period)

##add column to log transform bout length data 
allbouts3$ln_bout_duration <- log(allbouts3$Bout_duration_s)

##subset to index birds only
allbouts3.in=subset(allbouts3,Status=="index")

##subset to peak bouts
allbouts3.peak=subset(allbouts3, inf_period=="PEAK")

##subset to flockmates only
allbouts3.fm=subset(allbouts3,Status=="flockmate")

##now match index sex to PIT tag
Male <- c('0700EDFF14', '0700EDD671', '0700EE1002', '0700EE1D5A',
          '0700EDB1D6', '0700EDABD6', '0700EDBC2A')
Female <- c('0700EE16FC', '0700EE0656')

##distinguish sex for index birds
allbouts3.in$Index_Sex <- ifelse(allbouts3.in$PIT_Tag %in% Male, 'M', 'F')


##now begin analyses
hist(allbouts3$Bout_duration_s)

##data are very right skewed so a GLMM is needed

##first look at flockmates vs. index birds by treatment at peak infection
##tried poisson family first given the right-skewed nature of the untransformed data (though using log link function so that may not matter)
glm_status=glmer(Bout_duration_s ~Treatment * Status + (1|PIT_Tag), data = allbouts3.peak, family = poisson (link="log"))
summary(glm_status)
Anova(glm_status)
emmeans(glm_status, pairwise~Treatment * Status, type="response")
summary(allEffects(glm_status), type="response")
##the interaction between treatment and status is significant

##but need to check for overdispersion using Ben Bolkers function and yes, the poisson is very overdispersed so I will use a gamma 
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(glm_status)

##now using gamma where overdispersion is irrelevant
glm_status=glmer(Bout_duration_s ~Treatment * Status + (1|PIT_Tag), data = allbouts3.peak, family = Gamma (link="log"))
summary(glm_status)
Anova(glm_status)
plot(glm_status)
resid(glm_status)
plot(resid(glm_status))
hist(resid(glm_status))
shapiro.test(resid(glm_status))
emmeans(glm_status, pairwise~Treatment * Status, type="response")
summary(allEffects(glm_status), type="response")
status_bouts_emm = emmeans(glm_status, ~Treatment * Status, type="response")
##export post-hoc pairwise comparisons for easy display in table format (these are the tables in supplementary material)
pm = pwpm(status_bouts_emm)
clipr::write_clip(pm)
pm
##the interaction between treatment and status is still very significant

##now look at just index birds over time (IF I include Index_Sex, it is not significant, so was removed from model)
glm_index=glmer(Bout_duration_s~Treatment * inf_period + (1|PIT_Tag), data = allbouts3.in, family = Gamma (link="log"))
summary(glm_index)
Anova(glm_index)
emmeans(glm_index, pairwise~Treatment * inf_period, type="response")
summary(allEffects(glm_index), type="response")
indextime_bouts_emm = emmeans(glm_index, ~Treatment * inf_period, type="response")
##export post-hoc pairwise comparisons for easy display in table format (these are the tables in supplementary material)
pm = pwpm(indextime_bouts_emm)
clipr::write_clip(pm)
pm



library(dplyr)

##Change name of facet label levels (time) for graphing purposes
levels(allbouts3.in$inf_period)<-c("PRE-INFECTION", "EARLY INFECTION", "PEAK INFECTION")
levels(allbouts3.in$inf_period)

##re-level flockmate and index for graph
allbouts3.peak$Status = as.factor(allbouts3.peak$Status)
levels(allbouts3.peak$Status)

library(tidyverse)
allbouts3.peak <- mutate(allbouts3.peak, 
                         Status=factor (Status, 
                                              levels=c("index", "flockmate")))
levels(allbouts3.peak$Status)

##Change name of facet label levels (status) for graphing purposes
levels(allbouts3.peak$Status)<-c("index birds", "flockmates")
levels(allbouts3.peak$Status)
table(allbouts3.peak$Status)


##plotting raw data on foraging bout length at peak infection - this was included in manuscript as Fig 4; used log10 transformed data because the raw data are very hard to visualize
Fig=ggplot(data=allbouts3.peak,aes(x=Treatment,y=ln_bout_duration))+
  facet_wrap(~Status,ncol=1,nrow=2)+ #this is creating multiple "panels" for site
  geom_boxplot()+
  geom_point(aes(), alpha = 0.03, size = 1, position=position_jitter(height=.05, width=.05))+
  ##ylab(expression(paste("Bout Duration (Log", [10], "(sec))"+ #change y axes
  ylab(expression(paste("Bout Duration ",log[10], "(sec)")))+
  xlab("Virulence Treatment")+
  theme_bw()+
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=18),
        axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 0, l = 0)),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        axis.text.x = element_text(face="italic"),
        strip.text = element_text(size = 20))
Fig



##now a graph for index birds over time- this graph was added to supplement as Fix S2 upon revision
Fig=ggplot(data=allbouts3.in,aes(x=Treatment,y=ln_bout_duration))+
  facet_wrap(~inf_period,ncol=1,nrow=3)+ #this is creating multiple "panels" for in period
  geom_boxplot()+
  geom_point(aes(), alpha = 0.03, size = 1, position=position_jitter(height=.05, width=.05))+
  ##ylab(expression(paste("Bout Duration (Log", [10], "(sec))"+ #change y axes
  ylab(expression(paste("Bout Duration ",log[10], "(sec)")))+
  xlab("Virulence Treatment")+
  theme_bw()+
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=18),
        axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 0, l = 0)),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        axis.text.x = element_text(face="italic"),
        strip.text = element_text(size = 20))
Fig


