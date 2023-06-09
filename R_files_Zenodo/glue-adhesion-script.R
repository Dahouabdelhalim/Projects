###Script associated with the article "The glue produced by Drosophila melanogaster for pupa adhesion is universal"

rm(list=ls(all=T))
setwd("/perso/borne/Documents/Projet_these/Pull_off_force_measurements/Kiel_experiments_August_2019/Article_redaction/Submitted_JEB")
#install.packages("gridExtra")
#install.packages("ggplot2")
library(ggplot2)
library(gridExtra)
library(reshape2)



#Upload the table
data <- read.csv("Borne_2020_glue_table_S1.csv", sep=",", fileEncoding = "UTF-7")
#n=328
#Correct maximal force with initial force
data$corr_exact_max_force <- data$max_force - data$mean_initial_force
#Select data which are good for analysis (data belongs to one of the four groups)
data_good <- subset(data, data$group == "peeling_tape" |
                      data$group == "peeling_pupa" |
                      data$group == "peeling_tape_and_peeling_pupa" |
                      data$group == "ok")
#n=235



#***********************************METHODS*************************************************
#Test if group and substrate impact on forces
anova1 <- aov(corr_exact_max_force~substrate*group, data= data_good)
summary(anova1)
#Substrate but not group impact forces



#***********************************RESULTS*************************************************


#######Pupa adhesion on different substrates################################################
#Test which substrates give significant difference in force
anova2 <- aov(corr_exact_max_force~substrate, data= data_good)
summary(anova2)
#post-hoc test to find pairwise comparison between substrates and correct for multiple tests
results_substrate <- TukeyHSD(anova2)
results_substrate
#Only comparisons with teflon are significant


#******FIGURE4******************************************
data_good$substrate <- factor(data$substrate, levels = c('glass-ROTH','glass-superfrost','glass-plasma','glass-PLL','glass-PLL-PEG', 'teflon', 'smooth-resin', '1MIC-resin','9MIC-resin','P1000-resin','P80-resin'), ordered = TRUE)
n_fun <- function(data_good){
  return(data.frame(y = 500, label = paste0("n = ",length(data_good))))
}

p1 <- ggplot2::ggplot(data= data_good, aes(x= substrate, y= corr_exact_max_force, fill = substrate))
p1 <- p1 + geom_boxplot(width= 0.4, outlier.colour = "white", fill = "grey") + xlab("substrate") + ylab("Force (mN)") 
p1 <- p1 + geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 1, position=position_jitter(0.1))
p1 <- p1 + stat_summary(fun.data = n_fun, geom = "text", size=6) 
p1 <- p1 + theme_classic() + theme(plot.title = element_text(color="black", size=20, face = "bold"), legend.position="none", axis.title = element_text(size = 20), axis.text = element_text(color="black", size=15))
p1

#*******************************************************


#######Contact areas########################################################################
#Select glass-type substrates
data_glass <- subset(data_good, data_good$substrate == "glass-superfrost" | 
                       data_good$substrate == "glass-plasma" | 
                       data_good$substrate == "glass-ROTH" | 
                       data_good$substrate == "glass-PLL" | 
                       data_good$substrate == "glass-PLL-PEG")
#Test if print types and substrate impact on forces
anova3 <- aov(corr_exact_max_force~substrate*rupture_types, data= data_glass)
summary(anova3)
#print types but not substrate impact force
#Find which comparison are significant
anova4 <- aov(corr_exact_max_force~rupture_types, data= data_glass)
summary(anova4)
#post-hoc test to find pairwise comparison between substrates and correct for multiple tests
results_print <- TukeyHSD(anova4)
results_print
#Only comparison between in_between and stays_on_pupa
#Values for these 2 groups
stays_pupa <- subset(data_glass, data_glass$rupture_types == "stays_on_pupa")
median(stays_pupa$corr_exact_max_force)
sd(stays_pupa$corr_exact_max_force)
in_between <- subset(data_glass, data_glass$rupture_types == "in_between")
median(in_between$corr_exact_max_force)
sd(in_between$corr_exact_max_force)

#******FIGURE5D******************************************
data_glass_print <- subset(data_glass, data_glass$rupture_types != "NA")
data_glass_print$rupture_types <- factor(data_glass_print$rupture_types, levels = c('stays_on_substrate','in_between','stays_on_pupa', ordered = TRUE))
p2 <- p2 <- ggplot(data_glass_print, aes(x= rupture_types, y=corr_exact_max_force, fill))
p2 <- p2 + geom_boxplot(outlier.colour = "white", width= 0.3, fill="grey") 
p2 <- p2 + geom_point(aes(group=rupture_types), shape = 21, colour = "black", fill = "white", size = 1.8, stroke = 1, position = position_jitter(width=0.1)) +
  xlab("print") + ylab("force")
p2 <- p2 + theme_classic() + theme(plot.title = element_text(color="black", size=20, face = "bold"), legend.position=c(0.9, 0.9), axis.title = element_text(size = 20), axis.text = element_text(color="black", size=15))
p2


#*******************************************************


#Convert px to mm for areas
data_glass$area_green_px <- as.numeric(as.character(data_glass$area_green_px))
data_glass$area_green_mm <- (data_glass$area_green_px)*(1/(data_glass$scale_print_px/data_glass$scale_print_mm)^2)
data_glass$area_black_px <- as.numeric(as.character(data_glass$area_black_px))
data_glass$area_black_mm <- (data_glass$area_black_px)*(1/(data_glass$scale_print_px/data_glass$scale_print_mm)^2)
data_glass$area_red_px <- as.numeric(as.character(data_glass$area_red_px))
data_glass$area_red_mm <- (data_glass$area_red_px)*(1/(data_glass$scale_print_px/data_glass$scale_print_mm)^2)

#Find if substrates impact on black area
anova5 <- aov(area_black_mm~substrate, data= data_glass)
summary(anova5)
#not significant

#Find if substrates impact on green area
anova5 <- aov(area_green_mm~substrate, data= data_glass)
summary(anova5)
#not significant

#Find if substrates impact on red area
anova5 <- aov(area_red_mm~substrate, data= data_glass)
summary(anova5)
#not significant


#******FIGURE5F******************************************
my_color <- c("black", "darkgreen", "red")
area_data <- melt(data_glass,id.vars='substrate', measure.vars=c('area_black_mm','area_green_mm','area_red_mm'))
area_data$substrate <- factor(area_data$substrate, levels = c('glass-ROTH','glass-superfrost','glass-plasma','glass-PLL','glass-PLL-PEG'), ordered = TRUE)
p3 <-  ggplot(area_data) 
p3 <-  p3 + geom_boxplot(outlier.colour = "white", width= 0.6, aes(x=substrate, y=value, col=variable)) + 
  geom_point(shape = 21, fill= "white", position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.4), aes(x=substrate, y=value, color=variable)) +
  ylab("area (mm²)")
p3 <- p3 + theme_classic() + theme(plot.title = element_text(color="black", size=20, face = "bold"), legend.position=c(0.9,0.9), axis.title = element_text(size = 20), axis.text = element_text(color="black", size=15))
p3 <- p3 + scale_colour_manual(name= "Print areas", values=my_color)
p3

#combined substrates
p3b <- ggplot(area_data)
p3b <- p3b + geom_boxplot(aes(x=variable, y=value))
p3b

#Minimum and maximum of red area
min(data_glass$area_red_mm, na.rm =T)
max(data_glass$area_red_mm, na.rm =T)
#Minimu and maximum of green area
min(data_glass$area_green_mm, na.rm =T)
max(data_glass$area_green_mm, na.rm =T)



#*******************************************************
#Find correlation between force and areas
cor.test( ~ corr_exact_max_force + area_black_mm,
         data=data_glass,
         method = "pearson",
         conf.level = 0.95)
lm_black <- lm(corr_exact_max_force~area_black_mm, data=data_glass)
summary(lm_black)

cor.test( ~ corr_exact_max_force + area_green_mm,
         data=data_glass,
         method = "pearson",
         conf.level = 0.95)
lm_green <- lm(corr_exact_max_force~area_green_mm, data=data_glass)
summary(lm_green)

cor.test( ~ corr_exact_max_force + area_red_mm,
         data=data_glass,
         method = "pearson",
        conf.level = 0.95)
lm_red <- lm(corr_exact_max_force~area_red_mm, data=data_glass)
summary(lm_red)



#***********************************SUPPLEMENTARY FIGURES*************************************************

#******SUP_FIGURE1******************************************
par(mfrow=c(2,1))
p4 <- ggplot(data_glass, aes(x=corr_exact_max_force, y=area_black_mm)) + geom_point() + xlab("Force") + ylab("black area")
#ajouter les droites de régression
p4 <- p4 + stat_smooth(method="lm", se=F, col="darkgrey")
p4 <- p4 + theme_classic() + theme(plot.title = element_text(color="black", size=20, face = "bold"), legend.position=c(0.9,0.9), 
                                   axis.title = element_text(size = 20), axis.text = element_text(color="black", size=15))
p4

p5 <- ggplot(data_glass, aes(x=corr_exact_max_force, y=area_green_mm)) + geom_point(col="darkgreen") + xlab("Force") + ylab("green area")
#ajouter les droites de régression
p5 <- p5 + stat_smooth(method="lm", se=F, col="darkgrey")
p5 <- p5 + theme_classic() + theme(plot.title = element_text(color="black", size=20, face = "bold"), legend.position=c(0.9,0.9), 
                                   axis.title = element_text(size = 20), axis.text = element_text(color="black", size=15))
p5

p6 <- ggplot(data_glass, aes(x=corr_exact_max_force, y=area_red_mm)) + geom_point(col="red") + xlab("Force") + ylab("red area")
#ajouter les droites de régression
p6 <- p6 + stat_smooth(method="lm", se=F, col="darkgrey")
p6 <- p6 + theme_classic() + theme(plot.title = element_text(color="black", size=20, face = "bold"), legend.position=c(0.9,0.9), 
                                   axis.title = element_text(size = 20), axis.text = element_text(color="black", size=15))
p6

grid.arrange(p4, p5, p6, nrow = 1)






###########TEST CORRELATION WITH OTHER PARAMETERS
#ranges for these parameters
min(data_good$temperature, na.rm =T)
max(data_good$temperature, na.rm =T)

min(data_good$humidity, na.rm =T)
max(data_good$humidity, na.rm =T)

min(data_good$atmospheric_pressure, na.rm =T)
max(data_good$atmospheric_pressure, na.rm =T)

min(data_good$max_age_hr, na.rm =T)
max(data_good$max_age_hr, na.rm =T)

#Test if humidity/Temperature/Pressure impact on forces
lm1 <- lm(corr_exact_max_force~substrate 
           + temperature
           + humidity 
           + atmospheric_pressure
           + max_age_hr, data=data_good)
summary(lm1)
#Only substrates is significant


###########DISTRIBUTION GROUPS BY SUBSTRATES
#Distribution of group depending on substrates
p9 <- ggplot(data_good, aes(x= substrate, fill= group))
p9 <- p9 + geom_bar()
p9




#Find correlation between force and areas for each glass substrates: 
data_glass <- subset(data_good, data_good$substrate == "glass-superfrost" | 
                       data_good$substrate == "glass-plasma" | 
                       data_good$substrate == "glass-ROTH" | 
                       data_good$substrate == "glass-PLL" | 
                       data_good$substrate == "glass-PLL-PEG")
# Glass-ROTH
data_substrate <- subset(data_glass, data_glass$substrate == "glass-ROTH")
lm_black <- lm(corr_exact_max_force~area_black_mm, data=data_substrate)
summary(lm_black)
lm_green <- lm(corr_exact_max_force~area_green_mm, data=data_substrate)
summary(lm_green)
lm_red <- lm(corr_exact_max_force~area_red_mm, data=data_substrate)
summary(lm_red)
#Not significant

# Glass-superfrost
data_substrate <- subset(data_glass, data_glass$substrate == "glass-superfrost")
lm_black <- lm(corr_exact_max_force~area_black_mm, data=data_substrate)
summary(lm_black)
lm_green <- lm(corr_exact_max_force~area_green_mm, data=data_substrate)
summary(lm_green)
lm_red <- lm(corr_exact_max_force~area_red_mm, data=data_substrate)
summary(lm_red)
#Not significant

# Glass-PLL
data_substrate <- subset(data_glass, data_glass$substrate == "glass-PLL")
lm_black <- lm(corr_exact_max_force~area_black_mm, data=data_substrate)
summary(lm_black)
lm_green <- lm(corr_exact_max_force~area_green_mm, data=data_substrate)
summary(lm_green)
lm_red <- lm(corr_exact_max_force~area_red_mm, data=data_substrate)
summary(lm_red)
#Not significant

# Glass-PLL-PEG
data_substrate <- subset(data_glass, data_glass$substrate == "glass-PLL-PEG")
lm_black <- lm(corr_exact_max_force~area_black_mm, data=data_substrate)
summary(lm_black)
lm_green <- lm(corr_exact_max_force~area_green_mm, data=data_substrate)
summary(lm_green)
lm_red <- lm(corr_exact_max_force~area_red_mm, data=data_substrate)
summary(lm_red)
#Not significant

# Glass-plama
data_substrate <- subset(data_glass, data_glass$substrate == "glass-plasma")
lm_black <- lm(corr_exact_max_force~area_black_mm, data=data_substrate)
summary(lm_black)
lm_green <- lm(corr_exact_max_force~area_green_mm, data=data_substrate)
summary(lm_green)
lm_red <- lm(corr_exact_max_force~area_red_mm, data=data_substrate)
summary(lm_red)
#Not significant

