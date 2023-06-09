rm(list=ls())
setwd("~/Desktop/DATA/Piper Partitioning 2018")
#setwd("~/Lauren/PPP JTE 2018/data")
library(tidyr) #tidy data
library(ggplot2) #graphing
library(dplyr) #data wrangling
library(MASS) #modern applied stats with s
library(ggthemes) #graphing
library(lme4) #mixed model library
library(multcomp) #Simultaneous Inference in General Parametric Models
library(vegan) #MDS
library(MuMIn)
library(plyr)
library(EcoSimR)

poop <- read.csv(file="Maynard_etal_CarolliaDiet.csv",head=TRUE)

poop<-poop[-c(56),]#bat doesn't have age, sex, rs info


##Objective 1A
#all Piper=1, no or some Piper=0
poop$Diet.Mix <- NA
for(i in 1:length(poop$prop_totpiper)){
	if(poop$prop_totpiper[i]==0){poop$Diet.Mix[i]="0"}
	if(poop$prop_totpiper[i]==1){poop$Diet.Mix[i]="1"}
	if(poop$prop_totpiper[i]>0 & poop$prop_totpiper[i]<1){poop$Diet.Mix[i]="0"}
}
poop$Diet.Mix <- as.numeric(poop$Diet.Mix)


##fixing age categories
levels(poop$age)
poop$age[which(poop$age=="SA")] <- "J"
poop$age[which(poop$age=="")] <- NA
poop$age <- droplevels(poop$age)

##fixing rs categories
levels(poop$rs)
poop$rs[which(poop$rs=="")] <- NA
poop$rs <- droplevels(poop$rs)



#glm on binomial diet data
#create global model
poop$Diet.Mix <- as.numeric(poop$Diet.Mix)
m1 <- glmer(Diet.Mix ~ species + age + sex + rs+ (1|j_date), data=poop, family=binomial, na.action="na.fail")
summary(m1)
drop1(m1, test="Chisq")
#age p-value=0.939932


#Model checking
cor(model.matrix(m1)[,-1])
plot(m1)

#simplified model with non-significant terms deleted
m2 <- glmer(Diet.Mix ~ species + sex + rs + (1|j_date), data=poop, family=binomial, na.action="na.fail")
summary(m2)
drop1(m2, test="Chisq")
##pvals: species=0.004703, sex=0.029890, rs=0.080403

#So there is a significant effect of species 
#and a significant effect of sex and a marginal effect of rs

#tukey test to assess differences among species
#Cp-Cc p=0.06, Cs-Cc p=0.0397, Cs-Cp p=0.982
summary(glht(m2, linfct=mcp(species="Tukey")))

#stats
aggregate(Diet.Mix~sex, data=poop, FUN=sum) 
summary(poop$sex)

#plots
poop$bat <- NA
for(i in 1:length(poop$species)){
	if(poop$species[i]=="Cc"){poop$bat[i]="C. castanea"}
	if(poop$species[i]=="Cp"){poop$bat[i]="C. perspicillata"}
	if(poop$species[i]=="Cs"){poop$bat[i]="C. sowelli"}
}
#categories for graph
poop$Diet.Mix.Cat <- NA
for(i in 1:length(poop$prop_totpiper)){
	if(poop$prop_totpiper[i]==0){poop$Diet.Mix.Cat[i]="Some or No Piper"}
	if(poop$prop_totpiper[i]==1){poop$Diet.Mix.Cat[i]="All Piper"}
	if(poop$prop_totpiper[i]>0 & poop$prop_totpiper[i]<1){poop$Diet.Mix.Cat[i]="Some or No Piper"}
}
poop$Sex <- NA
for(i in 1:length(poop$sex)){
	if(poop$sex[i]=="M"){poop$Sex[i]="Male"}
	if(poop$sex[i]=="F"){poop$Sex[i]="Female"}
}

leg_piper <- c(expression(paste("Exclusively", italic(" Piper "))), 
			   expression(paste("Mixed Genera or No ",italic("Piper"))))

sexdiet<-ggplot(poop[order(poop$Diet.Mix.Cat,decreasing=T),],aes(x = Sex,fill = factor(Diet.Mix.Cat))) + 
	geom_bar(position = position_fill(reverse = TRUE)) + coord_flip()+
	labs(y = "Proportion of samples", x = "Bat sex", fill = "")+
	theme_minimal()+
	scale_fill_manual(values=c('#88419d','#fdbb84'),labels = leg_piper )+
	guides(fill=guide_legend(title=''))+
	theme(legend.position="bottom")+
	scale_x_discrete(limits=c("Female","Male"))

speciesdiet<-ggplot(poop,aes(x = bat,fill = factor(Diet.Mix.Cat))) + 
	geom_bar(position = position_fill(reverse = TRUE)) + coord_flip()+
	labs(y = "", x = "Bat species", fill = "")+
	theme_minimal()+
	theme(axis.text.y = element_text(face = "italic"))+
	scale_fill_manual(values=c('#88419d','#fdbb84'))+ guides(fill=FALSE)+
	scale_x_discrete(limits=c("C. sowelli","C. perspicillata", "C. castanea"))

library(gridExtra)
library(ggpubr)

tiff('allpiper.tiff', units="in", width=5, height=5, res=300)
ggarrange(speciesdiet, sexdiet, 
		  labels = c("a", "b"),heights = c(2, 2),
		  ncol = 1, nrow = 2, align = "v")
dev.off()

sexdiet_bw<-ggplot(poop[order(poop$Diet.Mix.Cat,decreasing=T),],aes(x = Sex,fill = factor(Diet.Mix.Cat))) + 
	geom_bar(position = position_fill(reverse = TRUE)) + coord_flip()+
	labs(y = "Proportion of samples", x = "Bat sex", fill = "")+
	theme_minimal()+
	scale_fill_manual(values=c('#bdbdbd','#252525'),labels = leg_piper )+
	guides(fill=guide_legend(title=''))+
	theme(legend.position="bottom")+
	scale_x_discrete(limits=c("Female","Male"))

speciesdiet_bw<-ggplot(poop,aes(x = bat,fill = factor(Diet.Mix.Cat))) + 
	geom_bar(position = position_fill(reverse = TRUE)) + coord_flip()+
	labs(y = "", x = "Bat species", fill = "")+
	theme_minimal()+
	theme(axis.text.y = element_text(face = "italic"))+
	scale_fill_manual(values=c('#bdbdbd','#252525'))+ guides(fill=FALSE)+
	scale_x_discrete(limits=c("C. sowelli","C. perspicillata", "C. castanea"))

tiff('allpiper_bw.tiff', units="in", width=5, height=5, res=300)
ggarrange(speciesdiet_bw, sexdiet_bw, 
		  labels = c("a", "b"),heights = c(2, 2),
		  ncol = 1, nrow = 2, align = "v")
dev.off()


##Objective 1B
#group early species p/a
poop <- poop %>% 
	mutate(early_pa = pa_sf + pa_multi + pa_umbricola + pa_peltatum + pa_auritum)
#group mid and late species p/a
poop <- poop %>% 
	mutate(mid_late_pa = pa_urostachyum + pa_hold  + pa_cyn +pa_colonense + pa_ret + pa_pera + pa_nudifolium + pa_evasum)
sum(poop$early_pa)
sum(poop$mid_late_pa)



response.succ <- cbind(poop$early_pa, poop$mid_late_pa)
m3 <- glmer(response.succ ~ species + age + sex + rs+ (1|j_date), data=poop, family=binomial, na.action="na.fail")
summary(m3)
drop1(m3, test="Chisq") #use this p-value for rs, 0.75559

m4 <- glmer(response.succ ~ species + age + sex + (1|j_date), data=poop, family=binomial, na.action="na.fail")
summary(m4)
drop1(m4, test="Chisq") #use this p-value for sex, p=0.86511 (warning, failed to converge)

m5 <- glmer(response.succ ~ species + age + (1|j_date), data=poop, family=binomial, na.action="na.fail")
summary(m5)
drop1(m5, test="Chisq") #use this p-value for species, p=0.8215 (warning, failed to converge)

m6 <- glmer(response.succ ~ age + (1|j_date), data=poop, family=binomial, na.action="na.fail")
summary(m6)
drop1(m6, test="Chisq") #age p=0.0385

##significant effect of age only



##setting up the plot
pa.succ <- aggregate(response.succ~age, data=poop, FUN=sum) 
age = c("Adult","Adult", "Juvenile", "Juvenile") 
s = c("Early-successional", "Mid-/late-successional",
	  "Early-successional", "Mid-/late-successional") 
p = c(150/185, 35/185, 32/49, 17/49) 
df = data.frame(age, s, p) 

#load italics in legends
my_y_title_succ <- expression(paste(italic("Carollia "), " Age Class"))

#plot
succplot<-ggplot(df, aes(x=age, y=p, fill = factor(s))) +
	geom_bar(stat = 'identity',  position = position_fill(reverse = TRUE)) +
	theme_minimal() +
	theme(legend.position="right")+
	scale_fill_manual(values=c('#006837','#d9f0a3'))+
	coord_flip()+
	labs(y='Proportion of occurrences')+
	guides(fill=guide_legend(title=''))+
	theme(legend.position="right",legend.text.align = 0)+
	labs(x="Bat age class")

tiff('succplot.tiff', units="in", width=6, height=2, res=300)
succplot
dev.off()

succplot_bw<-ggplot(df, aes(x=age, y=p, fill = factor(s))) +
	geom_bar(stat = 'identity',  position = position_fill(reverse = TRUE)) +
	theme_minimal() +
	theme(legend.position="right")+
	scale_fill_manual(values=c('#252525','#bdbdbd'))+
	coord_flip()+
	labs(y='Proportion of occurrences')+
	guides(fill=guide_legend(title=''))+
	theme(legend.position="right",legend.text.align = 0)+
	labs(x="Bat age class")

tiff('succplot_bw.tiff', units="in", width=5, height=2, res=300)
succplot_bw
dev.off()

##Objective 1C
sp <- poop[,24:30]#5 common pipers, uncommon, non-piper
bat <- poop[,4:7]


#Run the NMDS 
mds1 <- metaMDS(sp, k=3,try=200)  #run the NMDS
mds1 #print the NMDS results
cor(vegdist(sp), dist(mds1$points))^2  #calculate R^2 for NMDS
##no convergence warnings

#model checking
library(goeveg)
dimcheckMDS(sp, distance = "bray", trymax = 20,
			autotransform = TRUE)

par(mfrow = c(1, 2))
stressplot(mds1, main = "Shepard plot")
gof <- goodness(mds1)
plot(mds1, type = "t", main = "Goodness of fit", display = 'sites')
points(mds1, cex = 2*gof/mean(gof))


##plots from stackoverflow
## with metaMDS plot:
plot(mds1, display="si", las=1, type = "n") # for an empty plot
points(mds1, pch = as.numeric(bat$species), col= as.numeric(bat$species))
## or with generic plot:
plot(mds1$points[,1], mds1$points[,2], pch = as.numeric(bat$species),
	 col= as.numeric(bat$species),
	 xlab = "NMDS1", ylab= "NMDS2",
	 asp = 1, las = 1) # this is new

#plot
#basic
plot(mds1$points[,1], mds1$points[,2], pch = as.numeric(bat$species),
	 col= as.numeric(bat$species),
	 xlab = "NMDS1", ylab= "NMDS2")

leg_mds <- c(expression(paste(italic("Carollia castanea"))),
			 expression(paste(italic("Carollia perspicillata"))), 
			 expression(paste(italic("Carollia sowelli"))))
leg_tit_mds <-expression(paste(italic("Carollia "), "Species"))	



#plot w/ labels, etc.
tiff('NMDSplots.tiff', units="in", width=8, height=7, res=300)
par(mfrow = c(2, 2))

nmdssp<-{plot(mds1$points[,1], mds1$points[,2], pch = c(1, 2, 0)[as.numeric(bat$species)],
	 col=c('#4d004b','#8c6bb1','#9ebcda')[as.numeric(bat$species)], 
	 xlab = "NMDS1", ylab= "NMDS2", las = 1, asp = 1)
legend('topright', legend = leg_mds, cex = .95, 
	   pch = c(1, 2, 0) [unique(bat$species)],
	   col = c('#4d004b','#8c6bb1','#9ebcda')[unique(bat$species)])

ordiellipse(mds1, bat$species, display = "sites", kind = "sd", label = F, lwd=2,
			col = c('#4d004b','#8c6bb1','#9ebcda')[unique(bat$species)])}



#plot with sex
nmdssex<-{plot(mds1$points[,1], mds1$points[,2], pch = c(3,5)[as.numeric(bat$sex)],
	 col=c('pink','navy')[as.numeric(bat$sex)], 
	 xlab = "NMDS1", ylab= "NMDS2", las = 1, asp = 1)
legend('topright', legend = c("Female", "Male"), cex = .95, 
	   pch = c(3,5) [unique(bat$sex)],
	   col = c('pink','navy')[unique(bat$sex)])

ordiellipse(mds1, bat$sex, display = "sites", kind = "sd", label = F, lwd=2,
			col = c('pink','navy')[unique(bat$sex)])}

#plot with age
nmdsage<-{plot(mds1$points[,1], mds1$points[,2], pch = c(4, 6)[as.numeric(bat$age)],
	 col=c('grey','black')[as.numeric(bat$age)], 
	 xlab = "NMDS1", ylab= "NMDS2", las = 1, asp = 1)
legend('topright', legend = c("Adult", "Juvenile"), cex = .95, 
	   pch = c(4, 6) [unique(bat$age)],
	   col = c('grey','black')[unique(bat$age)])

ordiellipse(mds1, bat$age, display = "sites", kind = "sd", label = F, lwd=2,
			col = c('grey','black')[unique(bat$age)])}

#plot with rs
nmdsrs<-{plot(mds1$points[,1], mds1$points[,2], pch = c(8, 9)[as.numeric(bat$rs)],
			   col=c('darkgreen','lightgreen')[as.numeric(bat$rs)], 
			   xlab = "NMDS1", ylab= "NMDS2", las = 1, asp = 1)
	legend('topright', legend = c("Non-reproductive", "Reproductive"), cex = .95, 
		   pch = c(8, 9) [unique(bat$rs)],
		   col = c('darkgreen','lightgreen')[unique(bat$rs)])
	
	ordiellipse(mds1, bat$rs, display = "sites", kind = "sd", label = F, lwd=2,
				col = c('darkgreen','lightgreen')[unique(bat$rs)])}
dev.off()



library(gridExtra)
library(ggpubr)


species = c("Cc", "Cp", "Cs")
s = c("1","1","1","2","2","2", "3", "3", "3",
	  "4","4","4", "5", "5", "5",
	  "6","6","6", "7", "7", "7")
p = c(38, 31, 69, 5, 4, 6, 4, 10, 10, 4, 5, 1, 11, 7, 15, 5, 4, 10, 1, 4 ,9) 
pie = data.frame(species, s, p) 
pie$bat <- NA
for(i in 1:length(pie$species)){
	if(pie$species[i]=="Cc"){pie$bat[i]="C. castanea"}
	if(pie$species[i]=="Cp"){pie$bat[i]="C. perspicillata"}
	if(pie$species[i]=="Cs"){pie$bat[i]="C. sowelli"}
}

#load italics in legends
my_y_title <- expression(paste(italic("Carollia "), "Species"))
my_leg_title <-expression(paste(italic("Piper "), "Species"))
leg_labs <- c(expression(paste(italic("P. sancti-felicis"))),expression(paste(italic("P. reticulatum"))), expression(paste(italic("P. colonense"))),
			  expression(paste(italic("P. multiplinervium"))), expression(paste(italic("P. umbricola"))),
			  expression(paste("Uncommon", italic(" Piper "), "spp.")), expression(paste("Unknown/non-",italic("Piper "), "spp.")))

dietplot<-ggplot(pie, aes(x=bat, y=p, fill = factor(s))) +
	geom_bar(stat = 'identity',  position = position_fill(reverse = TRUE)) +
	theme_minimal() +
	theme(legend.position="bottom")+
	coord_flip()+
	labs(y='Proportion of occurrences', x='')+
	guides(fill=guide_legend(title=''))+
	theme(legend.position="right", axis.text.y = element_text(face = "italic"),legend.text.align = 0)+
	scale_x_discrete(limits=c("C. sowelli","C. perspicillata", "C. castanea"))+
	scale_fill_manual(values=c('#4d004b','#810f7c','#8c6bb1','#8c96c6','#9ebcda','#bfd3e6','#fc8d62'), labels=leg_labs)

tiff('dietplot.tiff', units="in", width=6, height=2, res=300)
dietplot
dev.off()
	

#rarefaction curve
sp <- poop[,28:34]#5 common pipers, uncommon, non-piper
bat <- poop[,4:7]
sac <- specaccum(sp,method ="rarefaction", gamma = "chao1")

plot(sac, add= F, random = F, ci = 2, ci.type = c("polygon"), 
	 col = par("fg"), xvar = c( "individuals"))


#pianka index
poop_niche <- poop[,c(4,24:28,31:39)]  #species
sex_niche <-poop[,c(5,24:28,31:39)]  #sex
age_niche <-poop[,c(6,24:28,31:39)]  #age
rs_niche <-poop[,c(7,24:28,31:39)]  #rs

niche_sum <- poop_niche %>%
	group_by(species) %>%
	summarise_all(sum)

sex_sum <- sex_niche %>%
	group_by(sex) %>%
	summarise_all(sum)

age_sum <- age_niche %>%
	group_by(age) %>%
	summarise_all(sum)

rs_sum <- rs_niche %>%
	group_by(rs) %>%
	summarise_all(sum)


#Cc vs Cp
pianka(niche_sum[1:2,2:14])
#Cc vs Cp
pianka(niche_sum[c(1, 3),2:14])
#Cs vs Cp
pianka(niche_sum[2:3,2:14])
#sex
pianka(sex_sum[1:2,2:14])
#age
pianka(age_sum[1:2,2:14])
#rs
pianka(rs_sum[1:2,2:14])

##preference trials
pref <- read.csv(file="Maynard_etal_CarolliaPreference.csv",head=TRUE)
pref2 <- dplyr::select(pref, choice, bat_sp, piper, sex, age, rs, weight.bat, BatID, date) 
pref2$age[which(pref2$age=="SA")] <- "J"

m10=glmer(data=pref2, choice ~ bat_sp * piper + sex + age + rs + weight.bat + (1|BatID), 
         family=binomial, na.action = "na.fail")
summary(m10)
drop1(m10 ,test="Chisq")#report p-vals for all non-sig factors
#sex p=0.88, age p=0.165. rs p=0.9076, weight p=0.773

m20=glmer(data=pref2, choice ~ bat_sp * piper + (1|BatID), 
         family=binomial, na.action = "na.fail")
summary(m20)
drop1(m20, test="Chisq")  ##bat x piper p=0.01761
##Errors, believe the model is just a bit complex for this dataset


###Dividing data by Piper species

pref2_Pcol <- pref2[which(pref2$piper=="P. colonense"),]
pref2_Ppelt <- pref2[which(pref2$piper=="P. peltatum"),]
pref2_Psilv <- pref2[which(pref2$piper=="P. silvivagum"),]
pref2_Pret <- pref2[which(pref2$piper=="P. reticulatum"),]
pref2_Ppera <- pref2[which(pref2$piper=="P. peracuminatum"),]
pref2_Pumb <- pref2[which(pref2$piper=="P. umbricola"),]


##Running separate models for each
m1.Pcol <- glm(choice ~ bat_sp,family=binomial, 
               na.action = "na.fail", data=pref2_Pcol)
drop1(m1.Pcol, test="Chisq")
#Significant effect of species, p=0.03582


m1.Ppelt <- glm(choice ~ bat_sp, family=binomial, na.action = "na.fail", 
                data=pref2_Ppelt)
drop1(m1.Ppelt, test="Chisq")
#marginal effect of species, p=0.06627


m1.Psilv <- glm(choice ~ bat_sp, family=binomial, na.action = "na.fail", 
                data=pref2_Psilv)
drop1(m1.Psilv, test="Chisq")

m1.Pret <- glm(choice ~ bat_sp, family=binomial, na.action = "na.fail", 
               data=pref2_Pret)
drop1(m1.Pret, test="Chisq")

m1.Ppera <- glm(choice ~ bat_sp, family=binomial, na.action = "na.fail", 
                data=pref2_Ppera)
drop1(m1.Ppera, test="Chisq")

m1.Pumb <- glm(choice ~ bat_sp, family=binomial, na.action = "na.fail", 
               data=pref2_Pumb)
drop1(m1.Pumb, test="Chisq")

##No effects for other species

#Tukey test on two significant species

summary(glht(m1.Pcol, linfct=mcp(bat_sp="Tukey")))#marginal Cs-Cc, but overall nothing
summary(glht(m1.Ppelt, linfct=mcp(bat_sp="Tukey")))#nothing
