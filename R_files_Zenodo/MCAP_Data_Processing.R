################################## Formatting ############################################## #####
library(dplyr)
library(tidyverse)
library(janitor)
library(readxl)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)
library(rlang)
library(glmm)
library(lme4)
library(lsmeans)
library(effects)
library(emmeans)
library(car)
library(ResourceSelection)

setwd("~/Desktop")                                                  ####set your working directory
rawdata<-read_excel("MCAP_Symbiont_Data.xls")                                   ####import file
rawdata<-rawdata[-(1:6),]   #before comma = rows, and after = columns         ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them
#colnames --> on column 7 change to ct (also a rownames function)

##set arguments 
copy.n.C<-33   #copy number  --> 33 Actin gene sequence in C for this assay
copy.n.D<-3
copy.n.Mcap<-1
copy.n.ratioCD<-(copy.n.C/copy.n.D)
copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827  #differences in flur based on the reporter
fluo.D<-0
fluo.Mcap<-0.84815

data<-format_data%>%  #tidyverse
  dplyr::select(sample_name,target_name,ct_mean, ct_sd)%>%                           ####take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%  #way to create or change something in the data sheet w/out moving or deleting it (make columns numeric)    ####set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%   #take avg. sd. and mean and group it together                                  ####group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>% #has to be passed^ (. = all the columns) (no.omit = dont throw out NAs)  ####take mean of technical replicates
  filter(sample_name!="Control")%>%    #! = not equal
  filter(sample_name!="Control+")%>%
  filter(sample_name!="Control-")%>%
  filter(sample_name!="Contol-")%>%
  filter(sample_name!="control +")%>%
  filter(sample_name!="control -")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
Mcap<-filter(data,target_name=="Mcap")  #cbind = column bind (taking data and bringing together)      ####make separate sets of columns for C and D so that each sample has all its data by row

int<-left_join(left_join(C,D,by='sample_name'),Mcap,by="sample_name")
target_adults <- c("11B","11C","203B","209B","213C","219A","219B","219C","221B","3C","209B.2")
target_bundles <- c("11A","11A.2","12A","12B.2","19B","19B.2","1B","1B.3","201A","201A.2","201B","201B.2","203B","203B.2","209A","209A.2","211B","211B.2",
                    "213A","213A.2","213B","213B.2","214B","214B.2","219A","219A.2","219B","221A","221A.2","221B","2A","2A2","2B","2B.2","3A","3B","3B.2",
                    "4A","4A.2","4B","4B.2","2A.2","12B")
final<-int%>%
  select(-target_name.x,-target_name.y,-target_name)%>%                                                                                    ###final formatting, get rid of redundant sample name column from previous step
  filter(sample_name!="Control")%>%                                                                                                      ###remove control
  rename(d_mean=ct_mean.y)%>%rename(d_sd=ct_sd.y)%>%rename(c_mean=ct_mean.x)%>%rename(c_sd=ct_sd.x)%>%rename(mcap_mean=ct_mean)%>%rename(mcap_sd=ct_sd)%>%                                          ###rename redudntant names to C,D columns - all lowercase
  mutate(c_mean=c_mean-fluo.C)%>%mutate(d_mean=d_mean-fluo.D)%>%
  mutate(presence=case_when(d_mean>0~"CD",is.na(d_mean)~"C"))%>% ###set presence/absence in new column based on ct values
  mutate(dm_ratio=(2^(mcap_mean-d_mean))/copy.n.ratioDM)%>%    
  mutate(cm_ratio=(2^(mcap_mean-c_mean))/copy.n.ratioCM)%>%    
  mutate(cd_ratio=(cm_ratio/dm_ratio))%>%                                                                                               ###calculate ratio, log of 2^(difference in ct values)
  mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>% 
  mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%     #take out A,B,C on the sample_name --> regular expressions (numbers)                                                               ###extract numbers from sample name, set as colony number
  mutate(phenotype=case_when(colony%%2==0~"Nonbleached",colony%%2==1~"Bleached"))%>%    #%% = print the remainder                                                  ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%  #ross --> prop. of C
  mutate(prop_c=case_when(is.na(prop_c)~"1",TRUE ~ as.character(prop_c)))%>%
  mutate(prop_c=as.numeric(prop_c))%>%
  mutate(prop_d=1-prop_c)%>%   #prob of D (1-prop of C)
  mutate(abundance=case_when(prop_d>prop_c~"D>C",prop_c>prop_d~"C>D"))%>%
  mutate(sd_warning=case_when(c_sd>1~"c*",d_sd>1~"d*"))%>%                  #sd warning                                                             ###set warnings where standard deviation of tech replicates is >1
  mutate(ct_warning=case_when(c_mean>38~"c*",d_mean>38~"d*"))%>%           #ct warning                                                              ###set warnings where ct values is later than 34
  mutate(rep=case_when((is.na(d_sd)&d_mean>0)|(is.na(c_sd)&c_mean>0)~"1replicate"))%>%  ###set warning if only one replicate processed (no standard deviation but mean present)
  mutate(dataset="Larvae")%>% # or mutate(dataset="Parents")%>% or  mutate(dataset="Bundles")%>%
  select(sample_name,sd_warning,ct_warning,rep, colony,phenotype, everything())%>%
  filter(!(sample_name %in% target_bundles))

#can use select to re-order things ^ 
parents <- final 
bundles <- final 
larvae <- final

grand_final <- rbind(parents, bundles, larvae) #combine the rows into the final dataset 

saveRDS(grand_final, "final_data_MASTER")
write.table(grand_final,"final_data_MASTER.txt",sep="\\t",quote=FALSE,row.names=FALSE) #save as a text file

################################## Symbiont Chi-squared Analysis ############################################## 

data <- readRDS("final_data_MASTER")

#Bundles --> proportion Durusdinium
format.1<-data%>%ungroup()%>%
  filter(phenotype=="Bleached")%>%
  filter(dataset=="Bundles")%>%
  select(phenotype,prop_d,-sample_name,colony)

counts_bl<-as.data.frame(table(cut(format.1$prop_d, breaks=seq(0,1, by=.1)))) #counts of prop c from bleached colonies 

format.2<-data%>%ungroup()%>%
  filter(phenotype=="Nonbleached")%>%
  filter(dataset=="Bundles")%>%
  select(phenotype,prop_d,-sample_name, colony)

counts_nbl<-as.data.frame(table(cut(format.2$prop_d, breaks=seq(0,1, by=.1)))) #counts of prop c from nonbleached colonies

chi_data<-as.data.frame(t(cbind(counts_bl[,2],counts_nbl[,2])))
row.names(chi_data)<-c("nonbleached","bleached")

(Xsq <- chisq.test(chi_data)) #small sample size, so fisher.exact test will help
fisher.test(chi_data)

#Parents --> proportion Durusdinium

format.3<-data%>%ungroup()%>%
  filter(phenotype=="Bleached")%>%
  filter(dataset=="Parents")%>%
  select(phenotype,prop_d,-sample_name)

counts_bl<-as.data.frame(table(cut(format.3$prop_d, breaks=seq(0,1, by=.1)))) #counts of prop c from bleached colonies 

format.4<-data%>%ungroup()%>%
  filter(phenotype=="Nonbleached")%>%
  filter(dataset=="Parents")%>%
  select(phenotype,prop_d,-sample_name)

counts_nbl<-as.data.frame(table(cut(format.4$prop_d, breaks=seq(0,1, by=.1)))) #counts of prop c from nonbleached colonies

chi_data<-as.data.frame(t(cbind(counts_bl[,2],counts_nbl[,2])))
row.names(chi_data)<-c("nonbleached","bleached")

(Xsq <- chisq.test(chi_data)) #small sample size, so fisher.exact test will help
fisher.test(chi_data)

#Figure 2 historgram - ggolot
#Bundles
hist(format.1$prop_d,main="bleached")
hist(format.2$prop_d,main="nonbleached")
hist.chart <- rbind(format.1,format.2)
str(hist.chart)
hist.chart$phenotype=as.factor(hist.chart$phenotype)

b <- ggplot(hist.chart)+
  geom_histogram(aes(prop_d,fill=phenotype),bins = 7.5, color="black", alpha=0.7)+
  theme_classic()+
  scale_fill_manual(values=c("#F39C12","#0000FF"))+
  xlab(expression(paste("Proportion ", italic("Durusdinium "), "in Gamete Bundles")))+
  #xlab("Proportion of D Symbionts in Gamete Bundles")+
  ylab("Count")+
  scale_y_continuous(breaks=c(seq(0,50,5)))+
  scale_x_continuous(breaks=c(seq(0,1,.1)))+
  guides(fill=guide_legend(title="Phenotype"))
b 

#Parents
hist(format.1$prop_d,main="Bleached")
hist(format.2$prop_d,main="Nonbleached")
hist.chart <- rbind(format.3,format.4)
str(hist.chart)
hist.chart$phenotype=as.factor(hist.chart$phenotype)

a <- ggplot(hist.chart)+
  geom_histogram(aes(prop_d,fill=phenotype),bins = 7.5, color="black", alpha=0.7)+
  theme_classic()+
  scale_fill_manual(values=c("#F39C12","#0000FF"))+
  xlab(expression(paste("Proportion ", italic("Durusdinium "), "in Parent Colonies")))+
  #xlab("Proportion of D Symbionts in Parent Colonies")+
  ylab("Count")+
  scale_y_continuous(breaks=c(seq(0,35,5)))+
  scale_x_continuous(breaks=c(seq(0,1,.1)))+
  guides(fill=guide_legend(title="Phenotype"))
a

#larvae at t0 
Larvae_t0 <- data%>%filter(dataset=="Larvae")%>%
  filter(.,grepl("0", timepoint))%>%
  select(prop_d,treatment,phenotype, timepoint,cd_ratio)

L <- ggplot(Larvae_t0)+
  geom_histogram(aes(prop_d,fill=phenotype),bins = 7.5, color="black", alpha=0.7)+
  theme_classic()+
  scale_fill_manual(values=c("#F39C12", "gray","#0000FF"))+
  xlab(expression(paste("Proportion ", italic("Durusdinium "), "in Larvae 12 h Post Fertilization")))+
  #xlab("Proportion of D Symbionts in Larvae 12 h Post-fertilization")+
  ylab("Count")+
  scale_y_continuous(breaks=c(seq(0,40,5)))+
  scale_x_continuous(breaks=c(seq(0,1,.1)))+
  guides(fill=guide_legend(title="Phenotype"))
L 

#combining all 3 histograms: parents, bundles and larvae
prow <- plot_grid(a + theme(legend.position = "none"),
                  b + theme(legend.position = "none"),
                  L + theme(legend.position = "none"), align = 'vh', 
                  labels = c("a","b","c"),  hjust = -1, nrow = 3) #take out the legent in both 
legend <- get_legend(L) #add the legend in manually 
p <- plot_grid(prow, legend , rel_widths = c(1.5, .3))
p

################################ Symbiont GLM############################################

#Bundle samples from night 2 of spawning collections (07/13/2018)
bundles_b.d <- data%>%filter(dataset=="Bundles")%>%
  filter(.,grepl("B",sample_name))%>%
  group_by(colony)%>%
  summarise(bundles_prop_d=mean(prop_d))%>%
  mutate(phenotype=case_when(colony%%2==0~"Nonbleached",colony%%2==1~"Bleached"))%>%
  select(colony,prop_d,abundance,dataset,phenotype,presence)

#Adult samples
adults.d <- data%>%filter(dataset=="Parents")%>%
  group_by(colony)%>%
  summarise(adult_prop_d=mean(prop_d))%>%
  mutate(phenotype=case_when(colony%%2==0~"Nonbleached",colony%%2==1~"Bleached"))

#Both nonbleached and bleached phenotypes
LR_data_all<- left_join(adults.d,bundles_b.d,by="colony")%>% 
  select(-sample_name)%>%
  filter(adult_prop_d<1.0)%>%
  filter(prop_d<1.0)

LR_data_all$phenotype<-as.factor(LR_data_all$phenotype)

#GLM 
fit<-glm(LR_data_all$adult_prop_d~LR_data_all$prop_d,family=quasibinomial("logit"))
summary(fit)
plot(fit)

#nonbleached colonies only
LR_data_nonbleached<- left_join(adults.d,bundles_b.d,by="colony")%>% 
  select(-sample_name)%>%
  filter(phenotype=="Nonbleached")

#GLM
fit<-glm(LR_data_nonbleached$adult_prop_d~LR_data$prop_d,family=quasibinomial("logit"))
summary(fit)
plot(fit)

#Bleached colonies only
LR_data_bleached<- left_join(adults.d,bundles_b.d,by="colony")%>% 
  select(-sample_name)%>%
  filter(phenotype=="Bleached")

#GLM
fit<-glm(LR_data_bleached$adult_prop_d~LR_data$prop_d,family=quasibinomial("logit"))
summary(fit)
plot(fit)

#Symbiont:Host 

#Parents s:h 
adults.sh <- data%>%filter(dataset=="Parents")%>%
  group_by(colony)%>%
  summarise(adult_logsh=mean(logSH))

#Bundles (night 2) s:h
bundles_b.sh <- data%>%filter(dataset=="Bundles")%>%
  filter(.,grepl("B",sample_name))%>%
  select(colony,logSH)

SH_data_all<- left_join(adults.sh,bundles_b.sh,by="colony")%>% 
  select(-sample_name)%>%

#GLM
fit<-glm(SH_data_all$adult_logSH~SH_data_all$logSH,family=quasibinomial("logit"))
summary(fit)
plot(fit)
  
#Figure 3
#Nonbleached and bleached 
c <- ggplot(LR_data_all)+
  geom_point(aes(adult_prop_d,prop_d,fill=phenotype.x),size = 3,shape=21,color="black",alpha=0.7)+
  scale_fill_manual(values = c("#F39C12","#0000FF"))+
  theme_classic()+
  theme(axis.title = element_blank())+ #gets rid of axis labels
  theme(axis.title.y=element_blank())+
  #xlab("Proportion of D Symbionts in Parent Colonies")+
  #ylab("Proportion of D Symbionts in Gamete Bundles")+
  geom_smooth(aes(adult_prop_d,prop_d),method='lm',se = FALSE,color="black")+
  guides(fill=guide_legend(title="Phenotype"))+
  scale_x_continuous(breaks=c(seq(0,1,.25)))
c

#Nonbleached only 
d <- ggplot(LR_data_nonbleached)+
  geom_point(aes(adult_prop_d,prop_d,fill=phenotype),size = 3,shape=21,color="black",alpha=0.7)+
  scale_fill_manual(values = c("#0000FF"))+
  theme_classic()+
  theme(axis.title = element_blank())+ #gets rid of axis labels
  theme(axis.title.y=element_blank())+
  #xlab("Proportion of D Symbionts in Parent Colonies")+
  #ylab("Proportion of D Symbionts in Gamete Bundles")+
  geom_smooth(aes(adult_prop_d,prop_d),method='lm',se = FALSE,color="black")+
  guides(fill=guide_legend(title="Phenotype"))+
  scale_x_continuous(breaks=(seq(0, 1, .25)), limits = c(0, 1))
#scale_x_continuous(breaks=(seq(0,1,.25))) 
d

#Bleached only
e <- ggplot(LR_data_bleached)+
  geom_point(aes(adult_prop_d,prop_d,fill=phenotype),size = 3,shape=21,color="black",alpha=0.7)+
  scale_fill_manual(values = c("#F39C12"))+
  theme_classic()+
  theme(axis.title = element_blank())+ #gets rid of axis labels
  theme(axis.title.y=element_blank())+
  #xlab("Proportion of D Symbionts in Parent Colonies")+
  #ylab("Proportion of D Symbionts in Gamete Bundles")+
  geom_smooth(aes(adult_prop_d,prop_d),method='lm',se = FALSE,color="black")+
  guides(fill=guide_legend(title="Phenotype"))+
  scale_x_continuous(breaks=(seq(0, 1, .25)), limits = c(0, 1))+
  scale_y_continuous(breaks=(seq(0, 1, .25)), limits = c(0, 1))
e

#All GLM plots combined 
prow <- plot_grid(
  e + theme(legend.position="none"), #take out the legend 
  d + theme(legend.position="none"),
  c + theme(legend.position="none"),
  align = 'v',
  labels = c("a", "b", "c"),
  hjust = 1,
  nrow = 1
)
prow

legend <- get_legend(c + theme(legend.box.margin = margin(0, 0, 0, 12))
)#add the legend in manually
f <- plot_grid(prow, legend, rel_widths = c(3, .4))
f

y.grob <- textGrob("Proportion Durusdinium in Bundles", 
                   gp=gpar(fontface="bold", fontsize=10), rot=90)

x.grob <- textGrob("Proportion Durusdinium in Parents", 
                   gp=gpar(fontface="bold", fontsize=10))

grid.arrange(arrangeGrob(f, left = y.grob, bottom = x.grob))

################################ Symbiont ANOVA############################################

data$logSH = log(data$sh_ratio) #log of S:H and adds it to the file

aov1 <- data%>% 
  select(colony,phenotype,sample_name,prop_c, prop_d, logSH,dataset)%>%
  filter(.,grepl("B",sample_name)&dataset=="Bundles")

aov2 <- data%>% 
  select(colony,phenotype,sample_name,prop_c, prop_d,logSH,dataset)%>%
  filter(dataset=="Parents")
aovdata<-rbind(aov2,aov1) #combine adults and bundles by row

#two-way ANOVA for prop D
aov.d <- aov(prop_d ~ phenotype + dataset + phenotype*dataset, data = combined)
summary(aov.d)

#ANOVA for sh 
#log transform sh data

aovdata <- aovdata[-c(50,55,56),] #remove outliers from sh --> colonies 222 (a, and b) and 220b from adults 

aov.sh <- aov(logSH ~ phenotype + dataset + phenotype*dataset, data = aovdata)
summary(aov.sh)

#Figure 4: Boxplot

data.1 <- data%>%
  select(phenotype, prop_d, sample_name, logSH, sh_ratio, dataset)%>%
  filter(dataset=="Bundles")%>%
  filter(.,grepl("B",sample_name)&dataset=="Bundles")

data.2 <- data%>% 
  select(phenotype,prop_d, sample_name, logSH, sh_ratio, dataset)%>%
  filter(dataset=="Parents")
data.2 <- data.2[-c(50,55,56),]

data.3 <- data%>%
  select(phenotype,timepoint,prop_d,logSH, sh_ratio, dataset)%>%
  filter(timepoint=="0")

data.box.final <- rbind(data.1,data.2,data.3)

str(data.box.final)
data.box.final$phenotype=as.factor(data.box.final$phenotype)
data.box.final$dataset=as.factor(data.box.final$dataset)

plot1 <- ggplot(data = data.box.final, aes(x=dataset, y=prop_d))+
  geom_boxplot(aes(fill=phenotype),alpha=0.7)+
  theme_classic()+
  scale_fill_manual(values=c("#F39C12","gray", "#0000FF"))+
  xlab("Life History Stage")+
  ylab(expression(paste("Proportion ", italic("Durusdinium"))))+
  # ylab("Proportion of Durusdinium")+
  guides(fill=guide_legend(title="Phenotype"))
plot1

plot2 <- ggplot(data=data.box.final, aes(x=dataset, y=logSH)) +
  geom_boxplot(aes(fill=phenotype),alpha=0.7) +
  theme_classic()+
  scale_fill_manual(values=c("#F39C12","gray", "#0000FF"))+
  xlab("Life History Stage")+
  ylab("Log S:H Cell Ratio")+  
  #theme(axis.title = element_blank())+ #gets rid of axis labels
  #theme(axis.title.y=element_blank())+
  guides(fill=guide_legend(title="Phenotype"))
plot2

plot3 <- ggplot(data = data.box, aes(x= phenotype, y=prop_d))+
  geom_boxplot(aes(fill=phenotype),alpha=0.7)+
  theme_classic()+
  scale_fill_manual(values=c("#FFFFFF","#999999","#996633"))+
  #xlab("Life History Stage")+
  #ylab("Proportion of D Symbionts")+
  guides(fill=guide_legend(title="Phenotype"))
plot3

#Combined boxplots:proportion Durusdnium and symbiont densities
prow <- plot_grid(plot1 + theme(legend.position = "none") , 
                  plot2 + theme(legend.position = "none"), 
                  plot3 + theme(legend.position = "none"),
                  align = "v",
                  labels = c("a","b"),  hjust = -1, nrow = 1) #take out the legend in both 
legend <- get_legend(plot1) #add the legend in manually 
p <- plot_grid(prow, legend , rel_widths = c(1.5, .3))
p

x.grob <- textGrob("Life History Stage", 
                   gp=gpar(fontface="bold", fontsize=12))

grid.arrange(arrangeGrob(p, left = y.grob, bottom = x.grob))

################################ Larval Size ############################################

setwd("~/Desktop") #set working directory
size<-read_excel("MCAP_Size_Data.xlsx",sheet="Measurements") #to read your file
str(size) 

#Two-way ANOVA
two_way <- size %>% 
  separate(Phenotype, c("Temp", "Phenotype"))%>%
  select(Temp,Phenotype,Area)

anova2<-aov(Area~as.factor(Temp)*as.factor(Phenotype),data=two_way)
leveneTest(Area~as.factor(Temp)*as.factor(Phenotype),data=two_way)     #p=0.00398, reject null. Center does not equal median
summary(anova2) 

TukeyHSD(anova2)
#High cross diff than all others
#High:Bleached-Ambient:Bleached       0.0006039
#High:Nonbleached-High:Bleached        0.0000001
#Ambient:Cross-High:Bleached            0.0000000
#High:Cross-High:Bleached              0.0000000
#Ambient:Nonbleached-High:Bleached     0.0000045

#Only with the bleached are the phenotypic effects seen
#High:Nonbleached-Ambient:Nonbleached  0.9971686
#High:Cross-Ambient:Cross               1.0000000
#High:Bleached-Ambient:Bleached       0.0006039

#Figure 5: Bar graph

long_size<-size%>%
  select(Phenotype,Area)%>%
  group_by(Phenotype)%>%
  summarise(avg=mean(Area), #(Percentage[!is.na(Percentage)],na.rm=TRUE) <- this is to take the mean even with blank cells in the data
            stdev=sd(Area),
            se=(sd(Area)/sqrt(length(Area))))

str(long_size)

s <- ggplot(long_size, aes(x=Phenotype, y=avg, fill=Phenotype))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se, width=0))+
  theme_classic()+
  annotate("text", x = "Ambient Bleached", y = 0.10, label = "a")+
  annotate("text", x = "Ambient Cross", y = 0.11, label = "b")+
  annotate("text", x = "Ambient Nonbleached", y = 0.105, label = "ab")+
  annotate("text", x = "High Bleached", y = 0.10, label = "c")+
  annotate("text", x = "High Cross", y = 0.11, label = "b")+
  annotate("text", x = "High Nonbleached", y = 0.105, label = "ab")+
  ylab(expression(Size (mm^2)))+ 
  xlab("Treatment and Phenotype")+
  scale_y_continuous(breaks=seq(0,0.11,0.02),limits=c(0,0.11))+
  theme(legend.position="none")+
  scale_fill_manual(values=c("#66CCFF", "#0000CC", "#000099","#FF9999","#FF0000", "#CC0000"),
                    labels=c("Ambient Bleached","Ambient Cross","Ambient Nonbleached","High Bleached","High Cross","High Nonbleached"))+
  theme(axis.text.x = element_text(color="black",size=7))
s

################################ Temperature Data ############################################

setwd("~/Desktop") 
larvae<-read_excel("MCAP_Temperature_Data.xlsx",sheet="Sheet1") #to read your file
str(larvae)

temp_larvae<-larvae%>%       
  filter(Running_Time=="1") #Keeping data that was actually durnig the experiment duration

graph_larvae<-temp_larvae%>%
  select(Hours,High_Bleached,High_Cross,High_Nonbleached,Ambient_Bleach,Ambient_Cross,Ambient_Nonbleached)%>%
  gather(Treatment,"temperature",-Hours)
