##Omy5 genotypes from all years, relating to distance upstream and above vs. below barriers in Elder Creek
##The first part of the code looks at how upstream distance relates to proportion of migratory alleles per sample pool
##The second part of the code explores the change in density of migratory alleles for juvenile fish above and below barriers in Elder Creek
##This code makes Table 1, Fig 3, Fig 5, and Table S1 in the manuscript and associated analyses
##This code was written and run on a PC
##COntact Suzanne Kelson, skelson@berkeley.edu, with questions

library(ggplot2)
library(plyr);library(dplyr)
library(lme4);library(lmerTest)
library(rsq)
library(reshape2)
library(gridExtra)
library(ggpubr)
rm(list =ls())

##SET WORKING DIRECTORY
setwd("")

location_colors <- RColorBrewer::brewer.pal(6, "Set1") ##2-5 are Elder Creek, 6 is South Fork, 1 is Fox Creek
elder_colors <- location_colors[2:5]
fox_colors <- location_colors[1]

##read in data files
habitat <- read.csv("habitat_2014_2015_2016_2017.csv")
omy5 <- read.csv("platemaps_omy5genotypes.csv")
omy5 <- select(omy5, sample_ID, reason_included, FID, year, date, stream,location, pool,pool.unique, FL, wt, PIT, age_class, group, PCA.pc.1, PCA.pc.2, genotype,missing_data)
omy5 <- droplevels(subset(omy5, reason_included == "electrofishing survey"))
omy5$genotype[omy5$missing_data > 250]<- NA ##remove genotypes for fish that were missing data at over 250 SNPs


#merge with habitat data so that we only pools that were part of the survey
habitat <- select(habitat,location,pool,year,distance,mouthdistance,pool_run,length_m,SA_m2,vol_m3,maxD_m)
omy5_efish <- merge(omy5, habitat, by = c("location","pool","year"),all.y=T) ##keep all the rows in habitat frame - there are some pools where we caught no fish
omy5_efish$year <- as.factor(as.character(omy5_efish$year))
##remove pools with no genetic data, since we subsetted in 2015-2017 years
omy5_efish <- subset(omy5_efish, !(is.na(FID)))
omy5_efish$location<- factor(omy5_efish$location, levels=levels(omy5_efish$location)[c(3,2,1,4,5)])##re order levels from upstream to downstream


##Summarize data (Table 1)
data_summary <- ddply(omy5_efish, .(year, location),summarize,
                      num_pools = sum(!(duplicated(pool))),
                      num_fish = sum(!(is.na(FL))),
                      num_yoy = sum(age_class == "YOY"))
#write.csv(data_summary, "table1.csv",row.names = F)

#########################################################################
### PART ONE - Proportion of migratory alleles by distance upstream ######
##########################################################################

#function to summarize data frame
genotypes_by_distance_fun <- function(df){
  df_summary <- ddply(df, .(year,location,pool,mouthdistance,distance,SA_m2,length_m),summarize,
        num_resident = sum(genotype == "Resident",na.rm=T),
        num_het = sum(genotype== "Heterozygote",na.rm=T),
        num_migratory = sum(genotype == "Migratory",na.rm = T),
        num_fish = num_resident+num_het+num_migratory,
        prop_res_fish = round(num_resident/num_fish,2),
        prop_het_fish = round(num_het/num_fish,2),
        prop_mig_fish = round(num_migratory/num_fish,2),
        num_res_alleles = 2*num_resident+num_het,
        num_mig_alleles = 2*num_migratory + num_het,
        num_alleles = 2*num_fish,
        prop_res_alleles = round(num_res_alleles/num_alleles,2),
        prop_mig_alleles = round(num_mig_alleles/num_alleles,2))
        df_summary$num_migratory_m=df_summary$num_migratory/df_summary$length_m
        df_summary$num_resident_m=df_summary$num_resident/df_summary$length_m
        df_summary$num_het_m=df_summary$num_het/df_summary$length_m
        df_summary$num_migratory_m2=df_summary$num_migratory/df_summary$SA_m2
        df_summary$num_resident_m2=df_summary$num_resident/df_summary$SA_m2
        df_summary$num_het_m2=df_summary$num_het/df_summary$SA_m2
        df_summary$num_fish_m=df_summary$num_fish/df_summary$SA_m
        df_summary$num_fish_m2=df_summary$num_fish/df_summary$SA_m2
        df_summary$num_mig_alleles_m2 = df_summary$num_mig_alleles/df_summary$SA_m2
  return(df_summary)
}

##################
##Elder Creek ####
##################
genotype_summary_elder <- genotypes_by_distance_fun(subset(omy5_efish, stream == "Elder"))

###Plot of proportion migratory alleles by distance from mouth in Elder Creek (Fig 3a)
mig_alleles <- ggplot(data=genotype_summary_elder,aes(y = prop_mig_alleles, x = distance))+
  geom_smooth(data=subset(genotype_summary_elder,year=="2014"),method="lm",color="gray15",linetype=1,se=F,size=.6)+
  geom_smooth(data=subset(genotype_summary_elder,year=="2015"),method="lm",color="gray35",linetype=1,se=F,size=.6)+
  geom_smooth(data=subset(genotype_summary_elder,year=="2016"),method="lm",color="gray55",linetype=1,se=F,size=.6)+
  geom_smooth(data=subset(genotype_summary_elder,year=="2017"),method="lm",color="gray71",linetype=1,se=F,size=.6)+
  geom_point(data=genotype_summary_elder, aes(size=num_fish/5,color=location,pch=year),alpha=.6,size=1.1)+
  labs(y = "Proportion Migratory Alleles", x = "Distance from Mouth (m)", title = "Elder Creek")+
  scale_color_manual(values=elder_colors,name="Location",labels=c("Below Waterfall","Above Waterfall","Misery","Paralyze"))+
  theme_classic(9)+
  scale_shape_manual(values=c(15,16,17,18),name="Year")+
  scale_size_continuous(name="Num. Fish")+
  scale_linetype_manual(name="Year",values=c(1,2,3,4),labels=c("2014","2015","2016","2017"))+
  theme(plot.title=element_text(size=10),legend.key.size = unit(0.5, "lines"))
mig_alleles


###model as a proportion
fit14 <- glm(prop_mig_alleles~distance,family="quasibinomial",data=subset(genotype_summary_elder, year == "2014"))
fit15 <- glm(prop_mig_alleles~distance,family="quasibinomial",data=subset(genotype_summary_elder, year == "2015"))
fit16 <- glm(prop_mig_alleles~distance,family="quasibinomial",data=subset(genotype_summary_elder, year == "2016"))
fit17 <- glm(prop_mig_alleles~distance,family="quasibinomial",data=subset(genotype_summary_elder, year == "2017"))

attach(genotype_summary_elder) ##attach so that the package rsq works with binomial model
fit.table.glm <- (data.frame(year=c("2014","2015","2016","2017")))
fit.table.glm$r.squared <- c(rsq(fit14),rsq(fit15),rsq(fit16),rsq(fit17))
fit.table.glm$coeff <- c(summary(fit14)$coefficients[2,1],summary(fit15)$coefficients[2,1],summary(fit16)$coefficients[2,1],summary(fit17)$coefficients[2,1])
fit.table.glm$stderr <- c(summary(fit14)$coefficients[2,2],summary(fit15)$coefficients[2,2],summary(fit16)$coefficients[2,2],summary(fit17)$coefficients[2,2])
fit.table.glm$zval <- c(summary(fit14)$coefficients[2,3],summary(fit15)$coefficients[2,3],summary(fit16)$coefficients[2,3],summary(fit17)$coefficients[2,3])
fit.table.glm$pval <- c(summary(fit14)$coefficients[2,4],summary(fit15)$coefficients[2,4],summary(fit16)$coefficients[2,4],summary(fit17)$coefficients[2,4])
detach(genotype_summary_elder)

#############
##Fox Creek##
############
genotype_summary_fox <- genotypes_by_distance_fun(subset(omy5_efish, stream == "Fox"))

#plots with lines for each year when relationship is significant (2015 and 2017)
mig_alleles_fox <- ggplot(data=genotype_summary_fox, aes(y = prop_mig_alleles, x = distance))+
  geom_smooth(data=subset(genotype_summary_fox,year=="2015"),method="lm",col="gray35",se=F,linetype=1,size=0.6)+
  geom_smooth(data=subset(genotype_summary_fox,year=="2017"),method="lm",col="gray71",se=F,linetype=1,size=0.6)+
  geom_point(aes(size=num_fish/5,pch=year,color=location),alpha=.6,size=1.1)+
  scale_color_manual(values=fox_colors,guide=F)+theme_classic(9)+
  scale_shape_manual(values=c(15,16,17,18),name="Year")+
  scale_size_continuous(name="Num. Fish")+
  labs(y = "Proportion of Migratory Alleles", x = "Distance from Mouth (m)", title = "Fox Creek")+
  theme(plot.title=element_text(size=10),legend.key.size = unit(0.5, "lines"))
mig_alleles_fox


##FIGURE 3
proportion_plot <- ggarrange(mig_alleles, mig_alleles_fox, nrow = 2, labels = c("A","B"), font.label = list(size=11))
ggsave("proportion_plot.pdf", plot = proportion_plot,
       scale = 1, width = 3.3, height = 5, units = c("in"),dpi = 350 )

###model as a proportion
fit14 <- glm(prop_mig_alleles~distance,family="quasibinomial", data=subset(genotype_summary_fox, year == "2014"))
fit15 <- glm(prop_mig_alleles~distance,family="quasibinomial", data=subset(genotype_summary_fox, year == "2015"))
fit16 <- glm(prop_mig_alleles~distance,family="quasibinomial", data=subset(genotype_summary_fox, year == "2016"))
fit17 <- glm(prop_mig_alleles~distance,family="quasibinomial", data=subset(genotype_summary_fox, year == "2017"))

attach(genotype_summary_fox)
fit.table.glm.fox <- (data.frame(year=c("2014","2015","2016","2017")))
fit.table.glm.fox$r.squared <- c(rsq(fit14),rsq(fit15),rsq(fit16),rsq(fit17))
fit.table.glm.fox$coeff <- c(summary(fit14)$coefficients[2,1],summary(fit15)$coefficients[2,1],summary(fit16)$coefficients[2,1],summary(fit17)$coefficients[2,1])
fit.table.glm.fox$stderr <- c(summary(fit14)$coefficients[2,2],summary(fit15)$coefficients[2,2],summary(fit16)$coefficients[2,2],summary(fit17)$coefficients[2,2])
fit.table.glm.fox$zval <- c(summary(fit14)$coefficients[2,3],summary(fit15)$coefficients[2,3],summary(fit16)$coefficients[2,3],summary(fit17)$coefficients[2,3])
fit.table.glm.fox$pval <- c(summary(fit14)$coefficients[2,4],summary(fit15)$coefficients[2,4],summary(fit16)$coefficients[2,4],summary(fit17)$coefficients[2,4])
detach(genotype_summary_fox)

###################################################
##### PART TWO - A - FOX JUVENILE FISH BY YEAR ####
###################################################
fox_efish_yoy <- subset(omy5_efish, stream=="Fox" & FL<85)
fox_yoy_summary<- genotypes_by_distance_fun(fox_efish_yoy)

fit <- glm(num_mig_alleles~SA_m2+year,data=fox_yoy_summary)
summary(fit)
anova(fit)

############################################################################
### PART TWO - B - JUVENILE FISH- ABOVE VS BELOW BARRIERS IN ELDER CREEK ####
############################################################################

#######################
#######Statistics - comparing above vs below each barrier using glms ####
###############
elder_yoy_efish <- subset(omy5_efish, stream == "Elder" & FL <85)
elder_yoy_summary <- genotypes_by_distance_fun(elder_yoy_efish)
elder_yoy_summary<- droplevels(elder_yoy_summary)
elder_yoy_summary$num_mig_alleles_m2 <- elder_yoy_summary$num_mig_alleles/elder_yoy_summary$SA_m2

##################################
### ABOVE VS BELOW WATERFALL #####
##################################

##look for effect of year on above vs. below comparison
elder_yoy_summary$above_below <- "Above"
elder_yoy_summary$above_below[elder_yoy_summary$location == "Below"]<- "Below"
elder_yoy_summary$above_below<- as.factor(elder_yoy_summary$above_below)
elder_yoy_summary$above_below<- factor(elder_yoy_summary$above_below, levels=levels(elder_yoy_summary$above_below)[c(2,1)])
fit_no_interaction<-glm(num_mig_alleles ~ SA_m2+above_below+year,data=elder_yoy_summary,family="poisson")
fit <- glm(num_mig_alleles ~ SA_m2+above_below*year,data=elder_yoy_summary,family="poisson")
summary(fit)
anova<-anova(fit, fit_no_interaction) ###compare full model and reduced model
##p-value for including an interaction effect: 
1-pchisq(abs(anova$Deviance[2]), abs(anova$Df[2]))
##choose the full model with the interaction effect, save in a dataframe for reference
abovebelow_table <- data.frame(summary(fit)$coefficients)

##since there is an interaction, plot each year seperately
fit14<-glm(num_mig_alleles ~ SA_m2+above_below,data=subset(elder_yoy_summary,year == "2014"),family="poisson")
fit15<-glm(num_mig_alleles ~ SA_m2+above_below,data=subset(elder_yoy_summary,year == "2015"),family="poisson")
fit16<-glm(num_mig_alleles ~ SA_m2+above_below,data=subset(elder_yoy_summary,year == "2016"),family="poisson")
fit17<-glm(num_mig_alleles ~ SA_m2+above_below,data=subset(elder_yoy_summary,year == "2017"),family="poisson")

##table S1 for above vs. below waterfall comparison
fit.table.abovebelow.yearly <- data.frame(rbind(summary(fit14)$coefficients,summary(fit15)$coefficients,summary(fit16)$coefficients,summary(fit17)$coefficients))
fit.table.abovebelow.yearly$barrier <- as.factor("Elder waterfall")
fit.table.abovebelow.yearly$year <- c(rep("2014",3),rep("2015",3),rep("2016",3),rep("2017",3))
fit.table.abovebelow.yearly$coefficient <- rep(c("Intercept","SA_m2","Upstream"),4)


##calculating confidence intervals for difference in means in each year
confint <- data.frame(year = c("2014","2015","2016","2017"), lower = NA, upper = NA,pvalue = NA)
elder_yoy_summary_subset<- subset(elder_yoy_summary, year == "2014")
ttest<- t.test(num_mig_alleles_m2~above_below, data = elder_yoy_summary_subset)
confint[1,2]<- ttest$conf.int[1];confint[1,3] = ttest$conf.int[2];confint[1,4]=ttest$p.value
elder_yoy_summary_subset<- subset(elder_yoy_summary, year == "2015")
ttest <- t.test(num_mig_alleles_m2~above_below, data = elder_yoy_summary_subset)
confint[2,2]<- ttest$conf.int[1];confint[2,3] = ttest$conf.int[2];confint[2,4]=ttest$p.value
elder_yoy_summary_subset<- subset(elder_yoy_summary, year == "2016")
ttest <- t.test(num_mig_alleles_m2~above_below, data = elder_yoy_summary_subset)
confint[3,2]<- ttest$conf.int[1];confint[3,3] = ttest$conf.int[2];confint[3,4]=ttest$p.value
elder_yoy_summary_subset<- subset(elder_yoy_summary, year == "2017")
ttest <- t.test(num_mig_alleles_m2~above_below, data = elder_yoy_summary_subset)
confint[4,2]<- ttest$conf.int[1];confint[4,3] = ttest$conf.int[2];confint[4,4]=ttest$p.valu

#######################
###Look at above the waterfall vs. Misery
###############
above_vs_misery <- droplevels(subset(elder_yoy_summary, location == "Misery" | location == "Above"))
levels(above_vs_misery$location)
fit_no_interaction <- glm(num_mig_alleles~SA_m2+location+year, data=above_vs_misery,family="poisson")
fit <- glm(num_mig_alleles~SA_m2+location*year, data=above_vs_misery,family="poisson")
#summary(fit)
anova<-anova(fit, fit_no_interaction) 
##p-value for including an interaction effect: 
1-pchisq(abs(anova$Deviance[2]), abs(anova$Df[2]))
#anova(fit)
##save the coefficients for the full model
abovemis_table <- data.frame(summary(fit)$coefficients)
fit14<-glm(num_mig_alleles ~ SA_m2+location,data=subset(above_vs_misery,year == "2014"),family="poisson")
fit15<-glm(num_mig_alleles ~ SA_m2+location,data=subset(above_vs_misery,year == "2015"),family="poisson")
fit16<-glm(num_mig_alleles ~ SA_m2+location,data=subset(above_vs_misery,year == "2016"),family="poisson")
fit17<-glm(num_mig_alleles ~ SA_m2+location,data=subset(above_vs_misery,year == "2017"),family="poisson")

##part of Table S1
fit.table.abovemisery.yearly <- data.frame(rbind(summary(fit14)$coefficients,summary(fit15)$coefficients,summary(fit16)$coefficients,summary(fit17)$coefficients))
fit.table.abovemisery.yearly$barrier <- as.factor("Misery confluence")
fit.table.abovemisery.yearly$year <- c(rep("2014",3),rep("2015",3),rep("2016",3),rep("2017",3))
fit.table.abovemisery.yearly$coefficient <- rep(c("Intercept","SA_m2","Upstream"),4)

##calculating confidence intervals for difference in means in each year
confint <- data.frame(year = c("2014","2015","2016","2017"), lower = NA, upper = NA,pvalue = NA)
above_vs_misery_subset<- subset(above_vs_misery, year == "2014")
ttest<- t.test(num_mig_alleles_m2~location, data = above_vs_misery_subset)
confint[1,2]<- ttest$conf.int[1];confint[1,3] = ttest$conf.int[2];confint[1,4]=ttest$p.value
above_vs_misery_subset<- subset(above_vs_misery, year == "2015")
ttest <- t.test(num_mig_alleles_m2~location, data = above_vs_misery_subset)
confint[2,2]<- ttest$conf.int[1];confint[2,3] = ttest$conf.int[2];confint[2,4]=ttest$p.value
above_vs_misery_subset<- subset(above_vs_misery, year == "2016")
ttest <- t.test(num_mig_alleles_m2~location, data = above_vs_misery_subset)
confint[3,2]<- ttest$conf.int[1];confint[3,3] = ttest$conf.int[2];confint[3,4]=ttest$p.value
above_vs_misery_subset<- subset(above_vs_misery, year == "2017")
ttest <- t.test(num_mig_alleles_m2~location, data = above_vs_misery_subset)
confint[4,2]<- ttest$conf.int[1];confint[4,3] = ttest$conf.int[2];confint[4,4]=ttest$p.value


###########
###Look at above waterfall vs. Paralyze
#########
above_vs_paralyze <- droplevels(subset(elder_yoy_summary, location == "Paralyze"|location == "Above"))
fit_no_interaction <- glm(num_mig_alleles~SA_m2+location+year, data=above_vs_paralyze,family="poisson")
fit <- glm(num_mig_alleles~SA_m2+location*year, data=above_vs_paralyze,family="poisson")
#summary(fit)
#anova(fit)
anova<-anova(fit, fit_no_interaction) 
##p-value for including an interaction effect: 
1-pchisq(abs(anova$Deviance[2]), abs(anova$Df[2]))

aboveparal_table <- data.frame(summary(fit)$coefficients)
##table and fits by year
fit14<-glm(num_mig_alleles ~ SA_m2+location,data=subset(above_vs_paralyze,year == "2014"),family="poisson")
fit15<-glm(num_mig_alleles ~ SA_m2+location,data=subset(above_vs_paralyze,year == "2015"),family="poisson")
fit16<-glm(num_mig_alleles ~ SA_m2+location,data=subset(above_vs_paralyze,year == "2016"),family="poisson")
fit17<-glm(num_mig_alleles ~ SA_m2+location,data=subset(above_vs_paralyze,year == "2017"),family="poisson")
##part of Table S1
fit.table.aboveparal.yearly <- data.frame(rbind(summary(fit14)$coefficients,summary(fit15)$coefficients,summary(fit16)$coefficients,summary(fit17)$coefficients))
fit.table.aboveparal.yearly$barrier <- as.factor("Paralyze confluence")
fit.table.aboveparal.yearly$year <- c(rep("2014",3),rep("2015",3),rep("2016",3),rep("2017",3))
fit.table.aboveparal.yearly$coefficient <- rep(c("Intercept","SA_m2","Upstream"),4)

confint <- data.frame(year = c("2014","2015","2016","2017"), lower = NA, upper = NA,pvalue = NA)
above_vs_paralyze_subset<- subset(above_vs_paralyze, year == "2014")
ttest<- t.test(num_mig_alleles_m2~location, data = above_vs_paralyze_subset)
confint[1,2]<- ttest$conf.int[1];confint[1,3] = ttest$conf.int[2];confint[1,4]=ttest$p.value
above_vs_paralyze_subset<- subset(above_vs_paralyze, year == "2015")
ttest <- t.test(num_mig_alleles_m2~location, data = above_vs_paralyze_subset)
confint[2,2]<- ttest$conf.int[1];confint[2,3] = ttest$conf.int[2];confint[2,4]=ttest$p.value
above_vs_paralyze_subset<- subset(above_vs_paralyze, year == "2016")
ttest <- t.test(num_mig_alleles_m2~location, data = above_vs_paralyze_subset)
confint[3,2]<- ttest$conf.int[1];confint[3,3] = ttest$conf.int[2];confint[3,4]=ttest$p.value
above_vs_paralyze_subset<- subset(above_vs_paralyze, year == "2017")
ttest <- t.test(num_mig_alleles_m2~location, data = above_vs_paralyze_subset)
confint[4,2]<- ttest$conf.int[1];confint[4,3] = ttest$conf.int[2];confint[4,4]=ttest$p.value

#############
###combined plot for downstream vs upstream comparison - Fig 5
#########

baseplot <- ggplot()+scale_fill_manual(values=c("gray30","gray80"),name="",labels=c("Downstream","Upstream"))+ theme_classic(9)

abovebelow <- baseplot+geom_violin(data=elder_yoy_summary,aes(x=year,y=num_mig_alleles_m2,fill=above_below))+
  labs(x="",y="Mig. Alleles/m2",title="Elder Waterfall")+
  annotate("text", x=1, y=2, label="*")+annotate("text", x=2, y=2, label="*")+annotate("text", x=3, y=2, label="*")+annotate("text", x=4, y=2, label = "*")
abovemis <- baseplot+geom_violin(data=above_vs_misery,aes(x=year,y=num_mig_alleles_m2,fill=location))+
  labs(x="",y="Mig. Alleles/m2",title="Misery Confluence")+
  annotate("text", x=1, y=.6, label="*")+annotate("text", x=2, y=.6, label="*")
aboveparal <- baseplot+geom_violin(data=above_vs_paralyze,aes(x=year,y=num_mig_alleles_m2,fill=location))+
  labs(x="",y="Mig. Alleles/m2",title="Paralyze Confluence")+
  annotate("text", x=2, y=.8, label="*")+annotate("text", x=3, y=.8, label="*")+annotate("text", x=4, y=.8, label = "*")
above_below_boxpl<- ggarrange(abovebelow,abovemis,aboveparal,nrow=3,common.legend = T, labels = c("A","B","C"), font.label = list(size=10))
ggsave("Fig5.pdf", plot = above_below_boxpl,
              scale = 1, width = 3, height = 4.5, units = c("in"),dpi = 350 )

##COMBINE FIT TABLES FOR YEARLY MODELS - Table S1
fit.tables.yearly <- rbind(fit.table.abovebelow.yearly, fit.table.abovemisery.yearly, fit.table.aboveparal.yearly)
fit.tables.yearly$coeff_min <- fit.tables.yearly$Estimate-fit.tables.yearly$Std..Error
fit.tables.yearly$coeff_max <- fit.tables.yearly$Estimate+fit.tables.yearly$Std..Error
fit.tables.yearly.subset <- subset(fit.tables.yearly, coefficient == "Upstream"& Pr...z..<=0.05)


############################################################################## 
### PART THREE - SUMMARIZING JUVENILE DATA BY LOCATION WITHIN ELDER CREEK ##### 
###############################################################################
#this data is summarized in maps in the manuscript

###summarize percent of alleles in each location
percent_alleles_by_location <- ddply(elder_yoy_efish, .(year, location),summarize,
                                     num_mig_fish = sum(genotype == "Migratory",na.rm=T),
                                    num_het_fish = sum(genotype == "Heterozygote",na.rm = T),
                                   num_res_fish = sum(genotype == "Resident",na.rm = T),
                                  num_fish = num_mig_fish+num_het_fish + num_res_fish,
                                 num_mig_alleles = 2*sum(genotype == "Migratory",na.rm=T)+sum(genotype=="Heterozygote",na.rm=T),
                                num_res_alleles = 2*sum(genotype == "Resident",na.rm=T)+sum(genotype=="Heterozygote",na.rm=T),
                               num_alleles = 2*sum(genotype == "Migratory",na.rm=T)+2*sum(genotype=="Heterozygote",na.rm=T)+2*sum(genotype == "Resident",na.rm=T))

alleles_fish_by_year<- ddply(percent_alleles_by_location, .(year),summarize,
                         num_mig_alleles_year_total = sum(num_mig_alleles),
                        num_res_alleles_year_total = sum(num_res_alleles),
                       num_fish_year_total = sum(num_fish))

percent_alleles_by_location<- merge(alleles_fish_by_year,percent_alleles_by_location, by = "year")
percent_alleles_by_location$percent_mig_for_region <- round(percent_alleles_by_location$num_mig_alleles/percent_alleles_by_location$num_mig_alleles_year_total*100,1)
percent_alleles_by_location$percent_res_for_region <- round(percent_alleles_by_location$num_res_alleles/percent_alleles_by_location$num_res_alleles_year_total*100,1)

###Summarize Fox data
percent_alleles_by_year_fox <- ddply(fox_efish_yoy, .(year),summarize,
                                     num_mig_fish = sum(genotype == "Migratory",na.rm=T),
                                     num_het_fish = sum(genotype == "Heterozygote",na.rm = T),
                                     num_res_fish = sum(genotype == "Resident",na.rm = T),
                                     num_fish = num_mig_fish+num_het_fish + num_res_fish,
                                     num_mig_alleles = 2*sum(genotype == "Migratory",na.rm=T)+sum(genotype=="Heterozygote",na.rm=T),
                                     num_res_alleles = 2*sum(genotype == "Resident",na.rm=T)+sum(genotype=="Heterozygote",na.rm=T),
                                     num_alleles = 2*sum(genotype == "Migratory",na.rm=T)+2*sum(genotype=="Heterozygote",na.rm=T)+2*sum(genotype == "Resident",na.rm=T))

percent_alleles_by_year_fox$percent_mig <- round(percent_alleles_by_year_fox$num_mig_alleles/percent_alleles_by_year_fox$num_alleles*100,1)
percent_alleles_by_year_fox$percent_res <- round(percent_alleles_by_year_fox$num_res_alleles/percent_alleles_by_year_fox$num_alleles*100,1)
