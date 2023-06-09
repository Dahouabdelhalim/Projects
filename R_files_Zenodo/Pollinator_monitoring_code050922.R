obs <- read.csv(file = "MDO_data071321.csv", sep = ",", header = TRUE)
climate_dat <-read.csv(file = "climate_data.csv", sep = ",", header = TRUE)

dates <- unique(mg_dat$Date)

library(janitor)
library(ggplot2)
library(reshape2)
library(car)
library(bbmle)
library(multcomp)
library(emmeans)
library(vegan)
library(IDPmisc)
library(multcompView)
library(effects)
library(tidyr)
library(sars)
library(chron)
library(splitstackshape)
library(dplyr)
library(sjPlot)
library(jtools)
library(SiZer)

##format climate data
climate <- climate_dat
str(climate)
colnames(climate)[1] <- "Date"      
colnames(climate)[2] <- "TempF"      
colnames(climate)[3] <- "RelHumper"   
colnames(climate)[4] <- "StationPressure"      
colnames(climate)[5] <- "WindSpeedmph"      
colnames(climate)[6] <- "WindDirectiondeg"      
colnames(climate)[7] <- "Windgustmph"
colnames(climate)[8] <- "Percpitation1hr(in)"      

climate2 <- climate[c(1:3, 5, 8)]
climate2$TempF <- as.numeric(climate2$TempF)
climate2$RelHumper <- as.numeric(climate2$RelHumper)
climate2$WindSpeedmph <- as.numeric(climate2$WindSpeedmph)
climate2$`Percpitation1hr(in)` <- as.numeric(climate2$`Percpitation1hr(in)`)
climate2$Date <- as.character(climate2$Date)
climatedate <- data.frame(concat.split(climate2, split.col = 1, sep = " "))
climate_final <- climatedate[c(2:7)]
colnames(climate_final)[5] <- "Date"
colnames(climate_final)[6] <- "Time"


###format time
climate_final$Time<- paste(climate_final$Time, ":00", sep = "")
climate_final$Time <-  chron(times = climate_final$Time, format = c("h:m:s"))


##Check Portulaca 

port <- droplevels(subset(obs, Cultivar == "Portulaca 'Sundial Light Pink'"))
port$count <- 1
str(port)
port_ag <- aggregate(list(length = port$count), by = list(Time = port$Date), sum)

###############Observational data
###########Total Abundance
########


obs$total_viz <- rowSums(obs[c(9:14)])

total_in_30 <- aggregate(list(complete = obs$total_viz), by = list(Date = obs$Date, A_P = obs$A_P, StartTime = obs$StartTime, 
                                                                   Cultivar = obs$Cultivar, Replicate = obs$Replicate, Floral_Area_inch2 = obs$Floral_Area_inch2), sum)
obs2 <- merge(obs, total_in_30, by = c("Date", "StartTime", "Cultivar", "Replicate", "Floral_Area_inch2", "A_P")) 


obs2$Time <- as.character(obs2$Time)

obs2$Time2 <- chron(times = obs2$Time, format = c("h:m:s"))

obs2$mins <- minutes(obs2$Time2) + 1

IDs <- unique(obs2[c(1:7, 17)])
IDs2 <- IDs [rep(seq_len(nrow(IDs)), each = 31), ]
IDs2$mins <- as.numeric(rep(1:31))

obs3 <- unique(merge(obs2,IDs2, by = c("mins", "Cultivar", "Replicate", "Date","StartTime", "Floral_Area_inch2", "A_P", "InfloresenceNum", "complete"), all = TRUE))

obs3$total_viz <- ifelse(is.na(obs3$total_viz), 0, obs3$total_viz)

obs3$accumulated_viz <- ave(obs3$total_viz, obs3$Cultivar, obs3$Replicate, obs3$Date, FUN=cumsum)

###find the max vsitors/minute
accumulated_viz_fin <- aggregate(list(viz_min = obs3$accumulated_viz), by = list(Date = obs3$Date, A_P = obs3$A_P, StartTime = obs3$StartTime, 
                                                                                 Cultivar = obs3$Cultivar, Replicate = obs3$Replicate, Floral_Area_inch2 = obs3$Floral_Area_inch2, mins = obs3$mins, complete = obs3$complete), max)


accumulated_viz_fin$Floral_Area_m2 <- accumulated_viz_fin$Floral_Area_inch2*0.00064516

accumulated_viz_fin$vis_rate_by_min <- (accumulated_viz_fin$viz_min/accumulated_viz_fin$mins/accumulated_viz_fin$Floral_Area_m2)

accumulated_viz_fin$Plant_Name <- paste(accumulated_viz_fin$Cultivar, accumulated_viz_fin$Replicate, sep = "")

obs3.5 <- accumulated_viz_fin

###run piecewise regression for each cultivar/time
obs3.5$StartTime <- paste(obs3.5$StartTime, ":00", sep = "")
obs3.5$StartTime <- chron(times = obs3.5$StartTime, format = c("h:m:s"))
# obs3.5$StartTime <- as.numeric(obs3.5$StartTime)

###Add time of day variable
check <- unique(obs3.5[c("Cultivar", "Date", "StartTime", "complete")])
check2 <- check[order(check$StartTime),]
check3 <- check2[order(check2$Date),]
check3$t_o_d <- rep(c("Early", "Early", "Early", "Mid", "Mid", "Mid", "Late", "Late", "Late"))

##add time information back to original DF
Asymp2 <- merge(obs3.5, check3, by = c("Cultivar", "Date", "StartTime"))
Asymp3 <- droplevels(subset(Asymp2, vis_rate_by_min > 0))
piece_rep <- Asymp3 %>% group_by(Cultivar, Date, Replicate) %>% 
  do(model = piecewise.linear(.$mins, .$vis_rate_by_min, middle = 1))


piece_sum <- piece_rep$model
names(piece_sum) <- paste(piece_rep$Cultivar, piece_rep$Date, piece_rep$Replicate, sep = ",")
change.p <- data.frame(sapply(piece_sum, "[[", 1))
colnames(change.p)[1] <-  "change_point_mins"
change.p$plantID <- rownames(change.p)
change.p

##create a function to extract predicted values

fun_pred <- function(x)
{predict(x, c(1:31))}

pred_vals <- data.frame(list(sapply(piece_sum, fun_pred)), check.names = F)
pred_vals$mins <- c(1:31)


###format data for merging
pred_melt <- melt(pred_vals, id = "mins", measure.vars = 1:101)
colnames(pred_melt)[3] <- "predicted"
colnames(pred_melt)[2] <- "plantID"


Asymp2$plantID <- paste(Asymp2$Cultivar, Asymp2$Date, Asymp2$Replicate, sep = ",")

plot_fin <- merge(Asymp2, pred_melt, by = c("mins", "plantID"))

plot_fin$predicted <- ifelse(plot_fin$vis_rate_by_min == 0, NA, plot_fin$predicted)
plot_fin_cult <- droplevels(subset(plot_fin, Cultivar == "Digitalis 'Panther Pink'" ))
cults <- unique(plot_fin$Cultivar)
cults
ggplot(plot_fin_cult, aes( x = mins, y = vis_rate_by_min, group = plantID)) + geom_point(aes(fill = Cultivar), color = "grey") + 
  geom_line(aes(y = predicted, color = Cultivar), size = 1.5) + theme(legend.position = "none") + theme_classic() + scale_color_manual(values = '#66a182')

#yellow
#scale_color_manual(values = '#d29b43')
#scale_color_manual(values = '#66a182')

###convert to DF for modelling, etc. 
MDO1 <- cbind(rownames(change.p), data.frame(change.p, row.names=NULL))
MDO2 <- data.frame(concat.split(MDO1, split.col = 1, sep = ","))
MDO3 <- MDO2[c(2:6)]
colnames(MDO3)[3] <- "Cultivar"
colnames(MDO3)[4] <- "Date"
colnames(MDO3)[5] <- "Replicate"

MDO4 <- merge(MDO3, check3, by = c("Cultivar", "Date"))

MDO4$Date <- as.character(MDO4$Date)

###plot

##first calculate the mean and SE

MDO_means <- aggregate(change_point_mins ~ Cultivar, MDO4, FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))

cultivars <- unique(obs3.5[c(2,4)])

MDO_means2 <- merge(MDO_means, cultivars, by = "Cultivar")

###then plot!
MDO_lims <- aes(ymax = change_point_mins[,1] + change_point_mins[,2], 
                ymin = change_point_mins[,1] - change_point_mins[,2])

MDO_means2$Cultivar <- factor(MDO_means2$Cultivar, levels = unique(MDO_means2$Cultivar [order(MDO_means2$change_point_mins[,1])]))

mean_MDO <- ggplot(MDO_means2, aes(x = Cultivar, y = change_point_mins[,1], fill = A_P)) + geom_bar(stat = "identity", position = "stack") + 
  theme_classic() + geom_errorbar(MDO_lims, width = 0) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_text(aes(label = round(change_point_mins[,1], digits = 1)), stat = "identity", nudge_y = 5,  color = "darkgrey") +  scale_fill_manual(values = c('#d29b43' ,'#66a182'))
mean_MDO


##
###Master Line
###Merge change point with DF
cults
plot_fin_cult2 <- droplevels(subset(plot_fin, Cultivar == "Verbena 'Meteor Shower'"))
plot_fin_cp <- merge(plot_fin_cult2, MDO_means2, by = c("Cultivar"))

ggplot(plot_fin_cp, aes(x = mins, y = vis_rate_by_min, group = Cultivar)) + geom_smooth(se = FALSE, color = "darkgrey", size = 3) + 
  theme_classic() + geom_point(aes(x = change_point_mins[,1]), y = 37, size = 6) + 
  geom_errorbarh(aes(xmin = change_point_mins[,1] - change_point_mins[,2], 
                     xmax = change_point_mins[,1] + change_point_mins[,2], y = 37, height = .5)) + ylim(0,150)


######Modelling MDO

#check distribution
hist(MDO4$change_point_mins)

##Merge with hourly climate data
MDO4$Hour <- hours(MDO4$StartTime)
climate_final$Hour <- hours(climate_final$Time)
MDO_clim <- merge(climate_final, MDO4, by = c("Date", "Hour"))

MDO_clim_fin <- merge(MDO_clim, obs3.5, by = c("Cultivar", "Date", "Replicate", "StartTime", "complete"))
MDO_clim_fin2 <- unique(MDO_clim_fin[c(1:10, 12:16, 19, 21)])

MDO_clim_fin2$vis_rate <- MDO_clim_fin2$complete/MDO_clim_fin2$Floral_Area_m2/30

###Assign time of day variable
MDO_clim_fin2$time_num <- paste(hours(MDO_clim_fin2$StartTime), (minutes(MDO_clim_fin2$StartTime)*(10/6)), sep = ".")

MDO_clim_fin2$t_o_d <- ifelse(MDO_clim_fin2$time_num <= 11.5, "AM", ifelse(MDO_clim_fin2$time_num >= 13.00, "PM", "MID"))


###Analysis

MDO_lm6.5 <- glmer(change_point_mins ~ scale(TempF) +  scale(Floral_Area_m2) +
                     WindSpeedmph + t_o_d + (1|Plant_Name) + Cultivar, MDO_clim_fin2, family = "Gamma"(link = "log"), 
                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

###GVIF1/2 should be less that 2 (Square it for VIF)
vif(MDO_lm6.5)

summary(MDO_lm6.5)

##Residuals plot#
##Verify that the model is a good fit
plot(MDO_lm6.5)

r.squaredGLMM(MDO_lm6.5)

car::Anova(MDO_lm6.5, type = "III")

#######
####Extract the slopes from MDO
slopes <- piece_sum
slopes_2 <- sapply(slopes, "[", 2)
slopes_2.5 <- sapply(slopes_2, "[", 1)
slopes_2.5 <- data.frame(slopes_2.5)
slopes_2.5

cult_names <- data.frame(names(slopes_2))
cult_names$order <- 1:101


slopes_3 <- cbind(rownames(slopes_2.5), data.frame(slopes_2.5, row.names=NULL))
slopes_5 <- melt(slopes_3, id.vars = 1, measure.vars = 2:102)

slopes_6 <- slopes_5 %>% 
  mutate(order = rep(1:101, each=3))

slopes_7 <- merge(cult_names, slopes_6, by = "order")

slopes_8 <- data.frame(concat.split(slopes_7, split.col = 2, sep = ","))
slopes_8 <- slopes_8[c(3, 5, 6, 7, 8)]
colnames(slopes_8 )[1] <- "variable"
colnames(slopes_8 )[3] <- "Cultivar"
colnames(slopes_8 )[4] <- "Date"
colnames(slopes_8 )[5] <- "Replicate"

slopes_a <- droplevels(subset(slopes_8, variable == "x"))

slopes_a$binned <- ifelse(slopes_a$value > 0, 0, 1)

###Run binomial model of binned slopes 

try_bin <- glm(binned ~ Cultivar, data = slopes_a, family = binomial)
post <- data.frame(emmeans(try_bin, ~ Cultivar , type = "response"))
###############

time_viz <- unique(merge(obs2,IDs2, by = c("mins", "Cultivar", "Replicate", "Date","StartTime", "Floral_Area_inch2", "A_P", "InfloresenceNum"), all = TRUE))


time_viz$Bee <- ifelse(is.na(time_viz$Bee), 0, time_viz$Bee)
time_viz$Fly <- ifelse(is.na(time_viz$Fly), 0, time_viz$Fly)
time_viz$Butterfly <- ifelse(is.na(time_viz$Butterfly), 0, time_viz$Butterfly)

time_viz_c <- droplevels(subset(time_viz, Cultivar == "Veronica 'Purple Candles'"))

time_plot <- ggplot(time_viz_c, aes(x = Time, y = Butterfly)) + geom_point() + facet_grid(~Replicate)
time_plot


time_viz$Bee_ac <- ave(time_viz$Bee, time_viz$Cultivar, time_viz$Replicate, time_viz$Date, FUN=cumsum)


#######Minimum time to determine taxonomic groups
###############
############

####Bees
bees_in_30 <- aggregate(list(bee_complete = obs$Bee), by = list(Date = obs$Date, A_P = obs$A_P, StartTime = obs$StartTime, 
                                                                Cultivar = obs$Cultivar, Replicate = obs$Replicate, Floral_Area_inch2 = obs$Floral_Area_inch2), sum)
bee_obs2 <- merge(obs, bees_in_30, by = c("Date", "StartTime", "Cultivar", "Replicate", "Floral_Area_inch2", "A_P")) 


###Flies

flies_in_30 <- aggregate(list(fly_complete = obs$Fly), by = list(Date = obs$Date, A_P = obs$A_P, StartTime = obs$StartTime, 
                                                                 Cultivar = obs$Cultivar, Replicate = obs$Replicate, Floral_Area_inch2 = obs$Floral_Area_inch2), sum)
fly_obs2 <- merge(obs, flies_in_30, by = c("Date", "StartTime", "Cultivar", "Replicate", "Floral_Area_inch2", "A_P")) 


###Leps

leps_in_30 <- aggregate(list(lep_complete = obs$Butterfly), by = list(Date = obs$Date, A_P = obs$A_P, StartTime = obs$StartTime, 
                                                                      Cultivar = obs$Cultivar, Replicate = obs$Replicate, Floral_Area_inch2 = obs$Floral_Area_inch2), sum)
lep_obs2 <- merge(obs, leps_in_30, by = c("Date", "StartTime", "Cultivar", "Replicate", "Floral_Area_inch2", "A_P")) 



####Merge it all together

taxa_obs3 <- unique(merge(bees_in_30,IDs2, by = c("Cultivar", "Replicate", "Date","StartTime", "Floral_Area_inch2", "A_P"), all = TRUE))
taxa_obs4 <- unique(merge(flies_in_30, taxa_obs3, by = c("Cultivar", "Replicate", "Date","StartTime", "Floral_Area_inch2", "A_P"), all = TRUE))
taxa_obs5 <- unique(merge(leps_in_30, taxa_obs4, by = c("Cultivar", "Replicate", "Date","StartTime", "Floral_Area_inch2", "A_P"), all = TRUE))


#########
######Define plants as "bee"  "fly" "butterfly" visited 

prop_tax <- data.frame(prop.table(as.matrix(taxa_obs5[c(7:9)]), margin = 1))
colnames(prop_tax)[1] <- "prop_lep"
colnames(prop_tax)[2] <- "prop_fly"
colnames(prop_tax)[3] <- "prop_bee"

taxa_obs_prop <- cbind(prop_tax, taxa_obs5)

############
#Classify based on proportion

taxa_obs_prop$Main_Taxa <- ifelse(c(taxa_obs_prop$prop_bee >= 0.3 & c(taxa_obs_prop$prop_fly < 0.3 & taxa_obs_prop$prop_lep < 0.3)), "Bee",
                                  ifelse(c(taxa_obs_prop$prop_fly >= 0.3 & c(taxa_obs_prop$prop_bee < 0.3 & taxa_obs_prop$prop_lep < 0.3)), "Fly", 
                                         ifelse(c(taxa_obs_prop$prop_lep >= 0.3 & c(taxa_obs_prop$prop_bee < 0.3 & taxa_obs_prop$prop_fly < 0.3)),"Lep", 
                                                ifelse(c(taxa_obs_prop$prop_bee >= 0.3 & taxa_obs_prop$prop_fly >= 0.3), "Bee_Fly", 
                                                       ifelse(c(taxa_obs_prop$prop_bee >= 0.3 & taxa_obs_prop$prop_lep >= 0.3), "Bee_Lep", "Lep_Fly")))))



#####Check with fig 1

####
test <- taxa_obs5

taxa_fig.5 <- aggregate(list(Bee =  test$bee_complete, Lep = test$lep_complete, Fly= test$fly_complete), by = list(
  Cultivar = test$Cultivar), FUN = mean)

prop_tax2 <- data.frame(prop.table(as.matrix(taxa_fig.5[c(2:4)]), margin = 1))
colnames(prop_tax2)[1] <- "prop_bee"
colnames(prop_tax2)[2] <- "prop_lep"
colnames(prop_tax2)[3] <- "prop_fly"

taxa_fig1 <- cbind(prop_tax2, taxa_fig.5)


taxa_fig1$Main_Taxa <- ifelse(c(taxa_fig1$prop_bee >= 0.3 & c(taxa_fig1$prop_fly < 0.3 & taxa_fig1$prop_lep < 0.3)), "Bee",
                              ifelse(c(taxa_fig1$prop_fly >= 0.3 & c(taxa_fig1$prop_bee < 0.3 & taxa_fig1$prop_lep < 0.3)), "Fly", 
                                     ifelse(c(taxa_fig1$prop_lep >= 0.3 & c(taxa_fig1$prop_bee < 0.3 & taxa_fig1$prop_fly < 0.3)),"Lep", 
                                            ifelse(c(taxa_fig1$prop_bee >= 0.3 & taxa_fig1$prop_fly >= 0.3), "Bee_Fly", 
                                                   ifelse(c(taxa_fig1$prop_bee >= 0.3 & taxa_fig1$prop_lep >= 0.3), "Bee_Lep", "Lep_Fly")))))







###Model 

taxa_mod <- taxa_obs_prop

taxa_mod$total_abund <- taxa_mod$lep_complete + taxa_mod$fly_complete + taxa_mod$bee_complete

taxa_mod1 <- unique(taxa_mod[c(1:13, 16:17)])

taxa_mod1$ID <- paste(taxa_mod1$Cultivar, taxa_mod1$Replicate, sep = "")

taxa_model1 <- glmer(total_abund ~ Main_Taxa + (1|ID) + Floral_Area_inch2, data = taxa_mod1, family = "poisson" (link = "log"))

summary(taxa_model1)

plot(taxa_model1)

r.squaredGLMM(taxa_model1)

car::Anova(taxa_model1, type = "III")


######Plot lines to represent waffle plot size 

taxa_ab_means <- data.frame(emmeans(taxa_model1, ~ Main_Taxa, type = "response"))


taxa_prop_means <- data.frame(emmeans(taxa_model1_prop, ~ Main_Taxa + variable, type = "response"))

taxa_ab_means$sqrt.means <- sqrt(taxa_ab_means$rate)


sqrt_plot <- ggplot(taxa_ab_means, aes(x = Main_Taxa, y = sqrt.means)) + geom_bar(stat = "identity")

sqrt_plot

####Taxa_proportion
taxa_melt_prop2 <-  na.omit(melt(taxa_mod1, id = c("Cultivar", "Replicate", "Date", "Main_Taxa", "Floral_Area_inch2", "ID"), measure.vars = c("prop_lep", "prop_fly", "prop_bee")))


tax_means <- aggregate(value ~ Main_Taxa + variable, taxa_melt_prop2, mean)

bees <- droplevels(subset(tax_means, Main_Taxa == "Bee"))
bees$per <- bees$value*100
bees2 <- bees[order(-bees$per),]
waffle(bees2$per)

taxa_model1_prop <- lm(value ~ Main_Taxa + variable + Main_Taxa:variable + ID, data = taxa_melt_prop2)

summary(taxa_model1_prop)



###plot

taxa_ab_means <- emmeans(taxa_model1, ~ Main_Taxa, type = "response", adjust = "tukey")


taxa_prop_means <- data.frame(emmeans(taxa_model1_prop, ~ Main_Taxa + variable, type = "response"))

###Make a waffle chart
library(waffle)

bees <- droplevels(subset(taxa_prop_means, Main_Taxa == "Lep"))
bees$sum <- sum(bees$emmean)
bees$per <- round(bees$emmean*100)
bees2 <- bees[order(-bees$per),]
waffle(bees2$per, rows = 10, size = 0.33)


######################
###find proportion of dominant taxa by minute
###################
#####################

##Merge taxa proportions with original minute/minute data

obs4 <- obs3

obs4$Bee <- ifelse(is.na(obs4$Bee), 0, obs4$Bee)
obs4$Fly <- ifelse(is.na(obs4$Fly), 0, obs4$Fly)
obs4$Butterfly <- ifelse(is.na(obs4$Butterfly), 0, obs4$Butterfly)

###order DF correctly for cumsum function
obs5<-obs4[order(obs4$Time),]
obs5.5 <-obs5[order(obs5$mins),]
obs6 <- obs5.5[order(obs5.5$Floral_Area_inch2),]

obs6$bee_acc_min <- ave(obs6$Bee, obs6$Cultivar, obs6$Replicate, obs6$Date, FUN=cumsum)
obs6$fly_acc_min <- ave(obs6$Fly, obs6$Cultivar, obs6$Replicate, obs6$Date, FUN=cumsum)
obs6$lep_acc_min <- ave(obs6$Butterfly, obs6$Cultivar, obs6$Replicate, obs6$Date, FUN=cumsum)

accumulated_taxa_viz_fin <- aggregate(list(bee_acc_min_fin = obs6$bee_acc_min, fly_acc_min_fin = obs6$fly_acc_min, 
                                           lep_acc_min_fin = obs6$lep_acc_min),
                                      by = list(Date = obs6$Date, A_P = obs6$A_P, StartTime = obs6$StartTime,
                                                Cultivar = obs6$Cultivar, Replicate = obs6$Replicate, 
                                                Floral_Area_inch2 = obs6$Floral_Area_inch2, 
                                                mins = obs6$mins), max)

##Proportion by minute

prop_tax_min <- data.frame(prop.table(as.matrix(accumulated_taxa_viz_fin[c(8:10)]), margin = 1))
colnames(prop_tax_min)[1] <- "prop_bee_min"
colnames(prop_tax_min)[2] <- "prop_fly_min"
colnames(prop_tax_min)[3] <- "prop_lep_min"

taxa_obs_prop_fin <- cbind(prop_tax_min, accumulated_taxa_viz_fin)

taxa_obs_prop_fin$prop_bee_min <- ifelse(is.na(taxa_obs_prop_fin$prop_bee_min), 0, taxa_obs_prop_fin$prop_bee_min)
taxa_obs_prop_fin$prop_fly_min <- ifelse(is.na(taxa_obs_prop_fin$prop_fly_min), 0, taxa_obs_prop_fin$prop_fly_min)
taxa_obs_prop_fin$prop_lep_min <- ifelse(is.na(taxa_obs_prop_fin$prop_lep_min), 0, taxa_obs_prop_fin$prop_lep_min)

####Merge with original Taxa DF

Taxa_final <- merge(taxa_obs_prop_fin, taxa_obs_prop, by = c("mins", "Cultivar", "Replicate", "Date", "Floral_Area_inch2", "A_P", "StartTime"), all = TRUE)



#####Plot taxa

taxa_melt <- melt(Taxa_final, id = c("Cultivar", "Date", "Replicate", "StartTime", "Floral_Area_inch2", "A_P", "Main_Taxa", "mins"), 
                  measure.vars = c("prop_bee_min", "prop_fly_min", "prop_lep_min"))


taxa_plot <- ggplot(taxa_melt, aes(x = mins, y = value, color = variable)) + geom_smooth(method =lm, formula = y ~ log(x)) + facet_grid(~Cultivar) 

taxa_plot

Taxa_min <- taxa_melt %>% group_by(Cultivar, Date, Replicate, variable)


Taxa_final$MT_det <- ifelse(c(Taxa_final$Main_Taxa == "Bee" & Taxa_final$prop_bee_min >= Taxa_final$prop_bee), "TRUE", 
                            ifelse(c(Taxa_final$Main_Taxa == "Fly" & Taxa_final$prop_fly_min >= Taxa_final$prop_fly), "TRUE",
                                   ifelse(c(Taxa_final$Main_Taxa == "Lep" & Taxa_final$prop_lep_min >= Taxa_final$prop_lep), "TRUE", 
                                          ifelse(c(Taxa_final$Main_Taxa == "Bee_Fly" & c(Taxa_final$prop_bee_min >= Taxa_final$prop_bee & Taxa_final$prop_fly_min >= Taxa_final$prop_fly)), "TRUE",
                                                 ifelse(c(Taxa_final$Main_Taxa == "Bee_Lep" & c(Taxa_final$prop_bee_min >= Taxa_final$prop_bee & Taxa_final$prop_lep_min >= Taxa_final$prop_lep)), "TRUE",
                                                        ifelse(c(Taxa_final$Main_Taxa == "Lep_Fly" & c(Taxa_final$prop_fly_min >= Taxa_final$prop_fly & Taxa_final$prop_lep_min >= Taxa_final$prop_lep)), "TRUE","FALSE"))))))


taxa_times <- droplevels(subset(Taxa_final, MT_det == "TRUE"))

Taxa_time_min <- aggregate(list(minutes  = taxa_times$min), by = list(Date = taxa_times$Date, A_P = taxa_times$A_P, StartTime = taxa_times$StartTime, 
                                                                      Cultivar = taxa_times$Cultivar, Replicate = taxa_times$Replicate, Floral_Area_inch2 = taxa_times$Floral_Area_inch2,
                                                                      Main_Taxa = taxa_times$Main_Taxa), min)

Taxa_time_min$plantID <- paste(Taxa_time_min$Cultivar, Taxa_time_min$Replicate, sep = " ")

hist(Taxa_time_min$minutes)


####Model and extract means 
###Add climate data

Taxa_time_min$StartTime <- paste(Taxa_time_min$StartTime, ":00", sep = "")
Taxa_time_min$StartTime <-  chron(times = Taxa_time_min$StartTime, format = c("h:m:s"))

####set hour for both
climate_final$Hour <- hours(climate_final$Time)
Taxa_time_min$Hour <- hours(Taxa_time_min$StartTime)

Taxa_time_min_clim <- merge(Taxa_time_min, climate_final, by = c("Date", "Hour"))
Taxa_time_min_clim$TempF_sc <- scale(Taxa_time_min_clim$TempF)
Taxa_time_min_clim$FA_sc <- scale(Taxa_time_min_clim$Floral_Area_inch2)

###Assign time of day variable
Taxa_time_min_clim$time_num <- paste(hours(Taxa_time_min_clim$StartTime), (minutes(Taxa_time_min_clim$StartTime)*(10/6)), sep = ".")

Taxa_time_min_clim$t_o_d <- ifelse(Taxa_time_min_clim$time_num <= 11.5, "AM", ifelse(Taxa_time_min_clim$time_num >= 13.00, "PM", "MID"))

##Remove outlier
Taxa_time_min_clim2 <- droplevels(subset(Taxa_time_min_clim, Floral_Area_inch2 != 14.38))
Taxa_time_min_clim3 <- droplevels(subset(Taxa_time_min_clim2, Main_Taxa != "Lep_Fly"))

###Use this model
Minutes_lm7 <- glmer(minutes ~ Cultivar + Main_Taxa + TempF_sc + t_o_d + FA_sc + 
                       (1|plantID), data = Taxa_time_min_clim3, family = "poisson",
                     glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

#Taxa_time_min_clim$residuals <- resid(Minutes_lm7)

library(MuMIn)

summary(Minutes_lm7)
plot(Minutes_lm7)

##Check r2
r.squaredGLMM(Minutes_lm7)

plot(predict(Minutes_lm7), Taxa_time_min_clim$minutes)
car::Anova(Minutes_lm7, typ = "III")

vif(Minutes_lm7)

##plot interactions

summary(Minutes_lm7)

Min_means <- emmeans(Minutes_lm7, pairwise ~ Main_Taxa, type = "response", adjust = "Tukey")
Min_means

###plot

Min_means_plot <- na.omit(data.frame(emmeans(Minutes_lm7,  ~ Main_Taxa, type = "response")))
Min_means_plot

c <- emmeans(Minutes_lm7,  ~ Main_Taxa, type = "response")

letters <- data.frame(cld(c, Letters = LETTERS))
letters

Min_means_plot$Main_Taxa <- factor(Min_means_plot$Main_Taxa, levels = Min_means_plot$Main_Taxa[order(Min_means_plot$rate)])

taxa_mins_plot <- ggplot(Min_means_plot, aes(x = Main_Taxa, y = rate, fill = Main_Taxa)) + geom_bar(stat = "identity", position = "dodge") +
  theme_classic() + geom_errorbar(aes(ymax = rate + SE,
                                      ymin = rate - SE), width = 0) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

taxa_mins_plot


############33
########Model overall attractiveness 
#############
########33

obs6 <- taxa_obs5
obs6$Floral_Area_m2 <- obs6$Floral_Area_inch2*0.00064516

obs6$Bee_rate <- obs6$bee_complete/obs6$Floral_Area_m2/30
obs6$Fly_rate <- obs6$fly_complete/obs6$Floral_Area_m2/30
obs6$Lep_rate <- obs6$lep_complete/obs6$Floral_Area_m2/30

obs_6_melt <- melt(obs6, id = c("Cultivar", "Date", "Replicate", "StartTime", "Floral_Area_m2", "A_P"), 
                   measure.vars = c("bee_complete", "fly_complete", "lep_complete"))


###Add climate data

obs_6_melt$StartTime <- paste(obs_6_melt$StartTime, ":00", sep = "")
obs_6_melt$StartTime <-  chron(times = obs_6_melt$StartTime, format = c("h:m:s"))

####set hour for both
climate_final$Hour <- hours(climate_final$Time)
obs_6_melt$Hour <- hours(obs_6_melt$StartTime)

obs6_clim <- merge(obs_6_melt, climate_final, by = c("Date", "Hour"))
obs6_clim$Tempc <- (obs6_clim$TempF - 32) * (5/9)

obs6_clim$WindSpeedkmph <- (obs6_clim$WindSpeedmph * 5) / 3.1
obs6_clim$TempC_sc <- scale(obs6_clim$Tempc)
obs6_clim$FA_sc <- scale(obs6_clim$Floral_Area_m2)
obs6_clim$Wind_Sc <- scale(obs6_clim$WindSpeedkmph)

obs6_clim$plantID <- paste(obs6_clim$Cultivar, obs6_clim$Replicate, sep = " ")
obs6_clim$value <- ifelse(is.na(obs6_clim$value), 0, obs6_clim$value)
hist(obs6_clim$value)

att_lm4 <- lmer(value ~ Cultivar + variable + TempC_sc + variable:TempC_sc + FA_sc +  variable:Wind_Sc + Wind_Sc+ variable:FA_sc + 
                  variable:Cultivar +
                  (1|plantID), data = obs6_clim)

summary(att_lm4)

###Remove bc/ of VIF
att_lm4_hum <- lmer(value ~ Cultivar + variable + TempC_sc + variable:TempC_sc + FA_sc +  variable:Wind_Sc + Wind_Sc+ variable:FA_sc + 
                      variable:Cultivar + RelHumper + variable:RelHumper + 
                      (1|plantID), data = obs6_clim)

att_lm4_full <- lmer(value ~ Cultivar + variable + Tempc + variable:Tempc + Floral_Area_m2 +  variable:WindSpeedkmph + WindSpeedkmph + variable:Floral_Area_m2 + 
                       variable:Cultivar +
                       (1|plantID), data = obs6_clim)

plot(predict(att_lm4), obs6_clim$value)

r.squaredGLMM(att_lm4)

vif(att_lm4_hum)
plot(att_lm5)

check_means <- data.frame(emmeans(att_lm4_full, ~ variable + Cultivar))
check_means

summary(att_lm5)
str(obs6_clim)

car::Anova(att_lm4_hum, type = "III")

my_cols <- colors = c("#c3d5e6", "#354d61", "#6089ad")

plot_model(att_lm4, type = "pred", terms = c("Wind_Sc", "variable"), 
           group.terms = c(1,2,3), colors = c("#6089ad", "#c3d5e6", "#354d61")) + theme_classic() + ylim(-10,50) 
plot_model(att_lm4_full, type = "pred", terms = c("WindSpeedkmph", "variable"), 
           group.terms = c(1,2,3), colors = c("#6089ad", "#c3d5e6", "#354d61")) + theme_classic() + ylim(-10,50)

plot_model(att_lm4, type = "pred", terms = c("FA_sc", "variable"), 
           group.terms = c(1,2,3), colors = c("#6089ad", "#c3d5e6", "#354d61")) + theme_classic() + ylim(-10,50)
plot_model(att_lm4_full, type = "pred", terms = c("Floral_Area_m2", "variable"), 
           group.terms = c(1,2,3), colors = c("#6089ad", "#c3d5e6", "#354d61")) + theme_classic() + ylim(-10,50)

plot_model(att_lm4, type = "pred", terms = c("TempF_sc", "variable"), 
           group.terms = c(1,2,3), colors = c("#6089ad", "#c3d5e6", "#354d61")) + theme_classic() +  ylim(-10,50)
plot_model(att_lm4_full, type = "pred", terms = c("Tempc", "variable"), 
           group.terms = c(1,2,3), colors = c("#6089ad", "#c3d5e6", "#354d61")) + theme_classic() + ylim(-10,50)


plot_model(att_lm4_hum, type = "pred", terms = c("RelHumper", "variable"), 
           group.terms = c(1,2,3), colors = c("#6089ad", "#c3d5e6", "#354d61")) + theme_classic() +  ylim(-10,50)

plot_model(att_lm4, type = "est", terms = c("Wind_Sc", "FA_sc", "Tempc_sc"))

############
####Correlate temp, wind speed, time of day
###########


Temp <- obs6_clim$TempF
Hour <- obs6_clim$Hour
Hum <- obs6_clim$RelHumper
cor.test(Temp, Hour, data = obs6_clim, method = "pearson")
# 

Wind <- obs6_clim$WindSpeedmph
cor.test(Wind, Hour, data = obs6_clim, method = "pearson")

#cor.test(Temp, Wind, data = obs6_clim, method = "pearson")

cor.test(Temp, Hum, data = obs6_clim, method = "pearson")

cor.test(Wind, Hum, data = obs6_clim, method = "pearson")

bee <- droplevels(subset(obs6_clim, variable == "bee_complete"))
vis <- bee$value
temp <- bee$Floral_Area_m2
cor.test(temp, vis, data = bee, method = "pearson")


##########3
########Rank attractiveness to each group by minute######
#######

Taxa_Ranks <- melt(Taxa_final, id = c("Cultivar", "Date", "Replicate", "StartTime", "Floral_Area_inch2", "A_P", "Main_Taxa", "mins"), 
                   measure.vars = c("bee_acc_min_fin", "lep_acc_min_fin", "fly_acc_min_fin"))

Taxa_Ranks$Floral_Area_m2 <- Taxa_Ranks$Floral_Area_inch2*0.00064516

Taxa_Ranks$value_stand <- Taxa_Ranks$value/ Taxa_Ranks$Floral_Area_m2

###Find Rank for each taxa group by minute

Ranks_per_min <- data.frame(Taxa_Ranks %>% group_by(Date, variable, mins) %>% 
                              mutate(rank_taxa_minutes = rank(value_stand)))

Ranks_per_min$Taxa_Group <- ifelse(Ranks_per_min$variable == "bee_acc_min_fin", "Bee", 
                                   ifelse(Ranks_per_min$variable == "fly_acc_min_fin", "Fly", "Lep"))

###Rank taxa total
rank_comp <- unique(obs_6_melt)
rank_comp$value_stand <- rank_comp$value / rank_comp$Floral_Area_m2


Ranks_comp1 <- data.frame(rank_comp %>% group_by(Date, variable) %>% 
                            mutate(rank_taxa_complete = rank(value_stand)))

Ranks_comp1$Taxa_Group <- ifelse(Ranks_comp1$variable == "bee_complete", "Bee", ifelse(Ranks_comp1$variable == "fly_complete", "Fly", "Lep"))

###Merge

Ranks_final <- merge(Ranks_per_min, Ranks_comp1, by = c("Date", "Cultivar", "Replicate", 
                                                        "A_P", "Taxa_Group"))

Ranks_final
####Do ICC for each minute 

icc_by_minute <- Ranks_final %>% group_by(Date, Taxa_Group, mins) %>% 
  do(Rank_cor = icc(.[c(14, 21)], type = "agreement")$value)


###Plot ICC to determine threshold 

icc_df <- data.frame(icc_by_minute)

rank_df_2 <- icc_df$Rank_cor

test <- lapply(rank_df_2, "first")

test2 <- data.frame(sapply(test, "[[", 1))

colnames(test2)[1] <- "icc_value"

icc_df2 <- cbind(icc_df, test2)


######Find the minimum time at which icc > 90

rank_times <- droplevels(subset(icc_df2, icc_value >= 0.90))

rank_times_min <- aggregate(list(mins= rank_times$mins), by = list(Date = rank_times$Date,
                                                                   Taxa_Group = rank_times$Taxa_Group), min)

rank_times_means <- aggregate(mins ~ Taxa_Group, rank_times_min, FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))

rank_times_means

icc_plot <- ggplot(icc_df2, aes(x = mins, y = icc_value, color = Taxa_Group)) + 
  geom_smooth(se = FALSE) + facet_grid(~Taxa_Group, scales = "free") + theme_classic()

icc_plot


#####Overall abundance
########3

Abund_Ranks <- melt(Asymp3 , id = c("Cultivar", "Date", "Replicate", "StartTime", "Floral_Area_inch2", "A_P", "mins"), 
                    measure.vars = c("vis_rate_by_min"))

###Find Rank for each taxa group by minute

Ab_Ranks_per_min <- data.frame(Abund_Ranks %>% group_by(Date, mins) %>% 
                                 mutate(rank_ab_minutes = rank(value)))

###Rank taxa total
rank_comp_ab <- unique(Asymp3[c(1:5, 8, 10)])

rank_comp_ab$comp_rate <- rank_comp_ab$complete.x/ rank_comp_ab$Floral_Area_m2

ab_Ranks_comp1 <- data.frame(rank_comp_ab %>% group_by(Date) %>% 
                               mutate(rank_ab_complete = rank(comp_rate)))

###Merge

ab_Ranks_final <- merge(Ab_Ranks_per_min, ab_Ranks_comp1, by = c("Date", "Cultivar", "Replicate", 
                                                                 "A_P"))

ab_Ranks_final
####Do ICC for each minute 

ab_icc_by_minute <- ab_Ranks_final %>% group_by(Date, mins) %>% 
  do(Rank_cor = icc(.[c(10,15)], type = "agreement")$value)


###Plot ICC to determine threshold 

ab_icc_df <- data.frame(ab_icc_by_minute)

ab_rank_df_2 <- ab_icc_df$Rank_cor

ab_test <- lapply(ab_rank_df_2, "first")

ab_test2 <- data.frame(sapply(ab_test, "[[", 1))

colnames(ab_test2)[1] <- "icc_value"

ab_icc_df2 <- cbind(ab_icc_df, ab_test2)


######Find the minimum time at which icc > 90

ab_rank_times <- droplevels(subset(ab_icc_df2, icc_value >= 0.90))

ab_rank_times_min <- aggregate(list(mins= ab_rank_times$mins), by = list(Date = ab_rank_times$Date), min)
ab_rank_times_min$x <- "x"
ab_rank_times_means <- aggregate(mins ~ x, ab_rank_times_min, FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))

ab_rank_times_means 
############################
###############
#MASTER GARDENER DATA
##############
##########################

mg_dat <- read.csv(file = "/Users/apple/Desktop/forHRI/MG_orders.csv", sep = ",", header = TRUE)
mg_dat2 <- read.csv(file = "/Users/apple/Desktop/forHRI/MG_descriptive.csv", sep = ",", header = TRUE)
mg_FA <- read.csv(file = "/Users/apple/Desktop/forHRI/MGobservations_FA.csv", sep = ",", header = TRUE)

library(ggplot2)
library(chron)
library(dplyr)
library(IDPmisc)
library(reshape2)
library(emmeans)
library(lme4)
library(car)
library(splitstackshape)
library(bbmle)

##simplify descriptive data set
mg_dat2$Bee <- mg_dat2$BumbleBee + mg_dat2$HoneyBee + mg_dat2$HairyLeg_LongAntennaeBee + mg_dat2$SmallGreenBee + mg_dat2$SmallDarkBee + mg_dat2$OtherBee
mg_dat2$Fly <- mg_dat2$BeeorWaspLikeFly + mg_dat2$HouseFly + mg_dat2$OtherFly
mg_dat2$Butterfly <- mg_dat2$WhiteButterfly + mg_dat2$BrownButterfly + mg_dat2$OtherButterfly
mg_dat2$Beetle <- mg_dat2$OrangeBeetle + mg_dat2$OtherBeetle

#Simplify DF for merging
mg_dat3 <- mg_dat2[c(1:5, 22:25, 20)]

#create variable for dataset type
mg_dat3$coll_type <- "Descriptive"
mg_dat3$ExperienceLevel <- "Novice"
colnames(mg_dat3)[10] <- "Unknown"
colnames(mg_dat)[11] <- "Unknown"
mg_dat$coll_type <- "Simple"

#Merge datasets

mg_merged <- merge(mg_dat, mg_dat3, by = c("Date", "Observer", "Cultivar", "Replicate", "ExperienceLevel", "StartTime", 
                                           "Bee", "Fly", "Butterfly", "Beetle", "Unknown", "coll_type"), all = TRUE)

dates <- unique(mg_merged$Date)

###fix duplicated observer IDS
mg_merged $Observer <- ifelse(c(mg_merged $Observer == "N" & mg_merged $coll_type == "Simple"), "X", mg_merged $Observer)
###Add Floral Area 

mg_FA_merged <- merge(mg_merged, mg_FA, by = c("Cultivar", "Replicate"))

##melt for plotting
mg_melt <- melt(mg_FA_merged, id = c("Cultivar", "Date", "Observer", "Replicate", "ExperienceLevel", "StartTime", "coll_type", "Area_inch2", "A_P"), 
                measure.vars = c("Bee", "Fly", "Butterfly", "Beetle", "Unknown"))
mg_melt$value <- as.numeric(mg_melt$value)
mg_melt$viz_rate <- mg_melt$value / mg_melt$Area_inch2

##Calculate mean and SE (for visualzing)

MG_obsmean <- aggregate(viz_rate ~ ExperienceLevel + coll_type + variable + A_P, mg_melt, FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))

# mg_anova <- aov(viz_rate ~ ExperienceLevel + variable + ExperienceLevel:variable + coll_type + variable:coll_type + A_P + ExperienceLevel:A_P, mg_melt)   
# summary(mg_anova)


###plot
limits <- aes(ymax = MG_obsmean$viz_rate[,1] + MG_obsmean$viz_rate[,2], 
              ymin = MG_obsmean$viz_rate[,1] - MG_obsmean$viz_rate[,2])

mean_obs <- ggplot(MG_obsmean, aes(x = variable, y = MG_obsmean$viz_rate[,1], fill = variable)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_classic() + geom_errorbar(limits, width = 0) +facet_grid(~ExperienceLevel + coll_type + A_P)
mean_obs


##convert to rank data (overall attraction)

ranks_cor <- mg_melt %>% group_by(Date, Observer, ExperienceLevel, variable, coll_type) %>% 
  mutate(rank_att_vr = rank(viz_rate))


###test correlation between rank attractiveness (Inter-rater reliability)
library(irr)

##Format dataframes (split in-ground and containers, make dfs wide)
##annuals
ans_mg <- data.frame(droplevels(subset(ranks_cor, A_P == "A")))

##subset taxa###
##BEES##
bees_ans_mg <- droplevels(subset(ans_mg, variable == "Bee"))
bees_ans_mg_novice <- droplevels(subset(bees_ans_mg, ExperienceLevel == "Novice"))
bee_ans_cast_n <- dcast(bees_ans_mg_novice, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Interclass correllation
#kappam.fleiss(bee_ans_cast_n[c(2:ncol(bee_ans_cast_n))])

icc(bee_ans_cast_n[c(2:ncol(bee_ans_cast_n))], type = "agreement")

###for experts
bees_ans_mg_exp <- droplevels(subset(bees_ans_mg, ExperienceLevel == "Expert"))
bee_ans_cast_n_e <- dcast(bees_ans_mg_exp, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Interclass correllation
#kappam.fleiss(bee_ans_cast_n_e[c(2:ncol(bee_ans_cast_n_e))])

icc(bee_ans_cast_n_e[c(2:ncol(bee_ans_cast_n_e))], type = "agreement")

#############
##Flies###
#########
fly_ans_mg <- droplevels(subset(ans_mg, variable == "Fly"))
fly_ans_mg_novice <- droplevels(subset(fly_ans_mg, ExperienceLevel == "Novice"))
fly_ans_cast_n <- dcast(fly_ans_mg_novice, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Intracalss correlation
#kappam.fleiss(fly_ans_cast_n[c(2:ncol(fly_ans_cast_n))])
icc(fly_ans_cast_n[c(2:ncol(fly_ans_cast_n))], type = "agreement")
###for experts
fly_ans_mg_exp <- droplevels(subset(fly_ans_mg, ExperienceLevel == "Expert"))
fly_ans_cast_n_e <- dcast(fly_ans_mg_exp, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Intraclass correlation
#kappam.fleiss(fly_ans_cast_n_e[c(2:ncol(fly_ans_cast_n_e))])
icc(fly_ans_cast_n_e[c(2:ncol(fly_ans_cast_n_e))], type = "agreement")

#############
##Butterfly###
#########
lep_ans_mg <- droplevels(subset(ans_mg, variable == "Butterfly"))
lep_ans_mg_novice <- droplevels(subset(lep_ans_mg, ExperienceLevel == "Novice"))
lep_ans_cast_n <- dcast(lep_ans_mg_novice, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Intraclass correlation
# kappam.fleiss(lep_ans_cast_n[c(2:ncol(lep_ans_cast_n))])
icc(lep_ans_cast_n[c(2:ncol(lep_ans_cast_n))], type = "agreement")

###for experts
lep_ans_mg_exp <- droplevels(subset(lep_ans_mg, ExperienceLevel == "Expert"))
lep_ans_cast_n_e <- dcast(lep_ans_mg_exp, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Intraclass correlation
#kappam.fleiss(lep_ans_cast_n_e[c(2:ncol(lep_ans_cast_n_e))])
icc(lep_ans_cast_n_e[c(2:ncol(lep_ans_cast_n_e))], type = "agreement")


#############
##Beetle###
#########
beet_ans_mg <- droplevels(subset(ans_mg, variable == "Beetle"))
beet_ans_mg_novice <- droplevels(subset(beet_ans_mg, ExperienceLevel == "Novice"))
beet_ans_cast_n <- dcast(beet_ans_mg_novice, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
kappam.fleiss(beet_ans_cast_n[c(2:ncol(beet_ans_cast_n))])
icc(beet_ans_cast_n[c(2:ncol(beet_ans_cast_n))], type = "agreement")

###for experts
beet_ans_mg_exp <- droplevels(subset(beet_ans_mg, ExperienceLevel == "Expert"))
beet_ans_cast_n_e <- dcast(beet_ans_mg_exp, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
kappam.fleiss(beet_ans_cast_n_e[c(2:ncol(beet_ans_cast_n_e))])
icc(beet_ans_cast_n_e[c(2:ncol(beet_ans_cast_n_e))], type = "agreement")


####Total abundance
####

tot_ab <- mg_melt %>% group_by(Date, Observer, ExperienceLevel, Cultivar) %>% 
  mutate(sum_abun = sum(value))
tot_ab <- data.frame(tot_ab)
tot_ab$abund_rate<- tot_ab$sum_abun/tot_ab$Area_inch2
tot_ab2 <- unique(tot_ab[c(1:9, 13:14)])

ranks_cor_ab <- tot_ab2 %>% group_by(Date, Observer, ExperienceLevel, coll_type) %>% 
  mutate(rank_att_vr = rank(abund_rate))

ranks_cor_ab_a <- droplevels(subset(ranks_cor_ab, A_P == "A"))
##Novice
tot_ans_mg_novice <- droplevels(subset(ranks_cor_ab_a, ExperienceLevel == "Novice"))
tot_ans_cast_n <- dcast(tot_ans_mg_novice, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
kappam.fleiss(tot_ans_cast_n[c(2:ncol(tot_ans_cast_n))])
icc(tot_ans_cast_n[c(2:ncol(tot_ans_cast_n))], type = "agreement")

##Experienced
tot_ans_mg_exp <- droplevels(subset(ranks_cor_ab_a, ExperienceLevel == "Expert"))
tot_ans_cast_e <- dcast(tot_ans_mg_exp, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
kappam.fleiss(tot_ans_cast_e[c(2:ncol(tot_ans_cast_e))])
icc(tot_ans_cast_e[c(2:ncol(tot_ans_cast_e))], type = "agreement")


#########
##Perennials
###################
###########
pers_mg <- droplevels(subset(ranks_cor, A_P == "P"))

##subset taxa###
##BEES##
bees_pers_mg <- droplevels(subset(pers_mg, variable == "Bee"))
bees_pers_mg_novice <- droplevels(subset(bees_pers_mg, ExperienceLevel == "Novice"))
bee_pers_cast_n <- dcast(bees_pers_mg_novice, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
icc(bee_pers_cast_n[c(2:ncol(bee_pers_cast_n))], type = "agreement")

###for experts
bees_pers_mg_exp <- droplevels(subset(bees_pers_mg, ExperienceLevel == "Expert"))
bee_pers_cast_n_e <- dcast(bees_pers_mg_exp, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
icc(bee_pers_cast_n_e[c(2:ncol(bee_pers_cast_n_e))], type = "agreement")

#############
##Flies###
#########
fly_pers_mg <- droplevels(subset(pers_mg, variable == "Fly"))
fly_pers_mg_novice <- droplevels(subset(fly_pers_mg, ExperienceLevel == "Novice"))
fly_pers_cast_n <- dcast(fly_pers_mg_novice, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
icc(fly_pers_cast_n[c(2:ncol(fly_pers_cast_n))], type = "agreement")

###for experts
fly_pers_mg_exp <- droplevels(subset(fly_pers_mg, ExperienceLevel == "Expert"))
fly_pers_cast_n_e <- dcast(fly_pers_mg_exp, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
icc(fly_pers_cast_n_e[c(2:ncol(fly_pers_cast_n_e))], type = "agreement")

#############
##Butterfly###
#########
lep_pers_mg <- droplevels(subset(pers_mg, variable == "Butterfly"))
lep_pers_mg_novice <- droplevels(subset(lep_pers_mg, ExperienceLevel == "Novice"))
lep_pers_cast_n <- dcast(lep_pers_mg_novice, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
icc(lep_pers_cast_n[c(2:ncol(lep_pers_cast_n))], type = "agreement")

###for experts
lep_pers_mg_exp <- droplevels(subset(lep_pers_mg, ExperienceLevel == "Expert"))
lep_pers_cast_n_e <- dcast(lep_pers_mg_exp, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
icc(lep_pers_cast_n_e[c(2:ncol(lep_pers_cast_n_e))], type = "agreement")


#############
##Beetle###
#########
beet_pers_mg <- droplevels(subset(pers_mg, variable == "Beetle"))
beet_pers_mg_novice <- droplevels(subset(beet_pers_mg, ExperienceLevel == "Novice"))
beet_pers_cast_n <- dcast(beet_pers_mg_novice, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

icc(beet_pers_cast_n[c(2:ncol(beet_pers_cast_n))], type = "agreement")

###for experts
beet_pers_mg_exp <- droplevels(subset(beet_pers_mg, ExperienceLevel == "Expert"))
beet_pers_cast_n_e <- dcast(beet_pers_mg_exp, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
icc(beet_pers_cast_n_e[c(2:ncol(beet_pers_cast_n_e))], type = "agreement")


####Total abundance
####

ranks_cor_ab_p <- droplevels(subset(ranks_cor_ab, A_P == "P"))
##Novice
tot_pers_mg_novice <- droplevels(subset(ranks_cor_ab_p, ExperienceLevel == "Novice"))
tot_pers_cast_n <- dcast(tot_pers_mg_novice, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
icc(tot_pers_cast_n[c(2:ncol(tot_pers_cast_n))], type = "agreement")

##Experienced
tot_pers_mg_exp <- droplevels(subset(ranks_cor_ab_p, ExperienceLevel == "Expert"))
tot_pers_cast_e <- dcast(tot_pers_mg_exp, Cultivar ~ Observer, value.var = "rank_att_vr", fun.aggregate = sum)

##Kappam.fleiss correlation
icc(tot_pers_cast_e[c(2:ncol(tot_pers_cast_e))], type = "agreement")


#########
########
#Is ICC more accurate with broader groups?##

ranks_cor_ab$hml <- ifelse(ranks_cor_ab$abund_rate > 1, 3, ifelse(c(ranks_cor_ab$abund_rate <=1 & ranks_cor_ab$abund_rate > .1), 2, 1))

ranks_cor_ab$hml <- ifelse(ranks_cor_ab$abund_rate == 0, 0, ranks_cor_ab$hml)


###Subset annuals                                                                                                                                  
ranks_hml_a <- droplevels(subset(ranks_cor_ab, A_P == "A"))
##Novice
ranks_hml_a_n <- droplevels(subset(ranks_hml_a , ExperienceLevel == "Novice"))
ranks_hml_a_c <- dcast(ranks_hml_a_n, Cultivar ~ Observer, value.var = "hml", fun.aggregate = sum)

## correlation
icc(ranks_hml_a_c [c(2:ncol(ranks_hml_a_c ))], type = "agreement")

##Experienced
ranks_hml_a_e <- droplevels(subset(ranks_hml_a , ExperienceLevel == "Expert"))
ranks_hml_a_ec <- dcast(ranks_hml_a_e, Cultivar ~ Observer, value.var = "hml", fun.aggregate = sum)

## correlation
icc(ranks_hml_a_ec [c(2:ncol(ranks_hml_a_ec ))], type = "agreement")

#################
####modelling
#######
#######
mg_FA_merged$StartTime <- paste(mg_FA_merged$StartTime, ":00", sep = "")
mg_FA_merged$StartTime <-  chron(times = mg_FA_merged$StartTime, format = c("h:m:s"))

####set hour for both
climate_final$Hour <- hours(climate_final$Time)
mg_FA_merged$Hour <- hours(mg_FA_merged$StartTime)

##merged
mg_climate_merge <- merge(mg_FA_merged, climate_final, by = c("Date", "Hour"), all = FALSE)

###calculate total abundance for modelling
mg_climate_merge$tot_abund <- rowSums(mg_climate_merge[c(8:12)])

##############################
##################
##Mixed Effects Model
#################
##############################


####Check coll type 

mg_nov <- droplevels(subset(mg_climate_merge, ExperienceLevel == "Novice"))
coll_mod <- glmer(tot_abund ~ offset(log(Area_inch2)) + Cultivar + scale(TempF) + coll_type+ 
                    (1|plantID), family = poisson, data = mg_nov, 
                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


# ##make plantIDvariable
mg_climate_merge$plantID <- paste(mg_climate_merge$Cultivar,mg_climate_merge$Replicate, sep = " ")
#
mg_climate_merge$Cultivar
mg_mod <- glmer(tot_abund ~ offset(log(Area_inch2)) + Cultivar + scale(TempF) + ExperienceLevel + 
                  (1|plantID), family = poisson, data = mg_climate_merge, 
                glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(mg_mod)
vif(mg_mod)

plot(mg_mod)

r.squaredGLMM(mg_mod)

car::Anova(mg_mod)

####Extract emmeans
experience_means <- emmeans(mg_mod, pairwise ~ ExperienceLevel, adjust = "tukey", type = "response")
experience_means
#
obs_means <- data.frame(emmeans(mg_mod, ~ ExperienceLevel))
#
limits <- aes(ymax = emmean + SE,
              ymin = emmean - SE)
#
mean_obs2 <- ggplot(obs_means, aes(x = ExperienceLevel, y = emmean, fill = ExperienceLevel)) + geom_bar(stat = "identity", position = "dodge") +
  theme_classic() + geom_errorbar(aes(ymax = emmean + SE,
                                      ymin = emmean - SE), width = 0) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
mean_obs2

################################
###################
########By Taxonomic Groups
################
#####################
#  

###make a proportional data table to identify dominant visitors 

taxa_prop <- data.frame(prop.table(as.matrix(mg_climate_merge[c(8:12)]), margin = 1))
colnames(taxa_prop)[1] <- "prop_bee"
colnames(taxa_prop)[2] <- "prop_fly"
colnames(taxa_prop)[3] <- "prop_lep"
colnames(taxa_prop)[4] <- "prop_beetle"
colnames(taxa_prop)[5] <- "prop_unk"

taxa_prop2 <- cbind(mg_climate_merge, taxa_prop)
taxa_prop2$prop_bee <- ifelse(is.na(taxa_prop2$prop_bee), 0, taxa_prop2$prop_bee)
taxa_prop2$prop_fly <- ifelse(is.na(taxa_prop2$prop_fly), 0, taxa_prop2$prop_fly)
taxa_prop2$prop_lep <- ifelse(is.na(taxa_prop2$prop_lep), 0, taxa_prop2$prop_lep)
taxa_prop2$prop_beetle <- ifelse(is.na(taxa_prop2$prop_beetle), 0, taxa_prop2$prop_beetle)
taxa_prop2$prop_unk <- ifelse(is.na(taxa_prop2$prop_unk), 0, taxa_prop2$prop_unk)

####

taxa_prop_melt <-melt(taxa_prop2, id = c("Cultivar", "Date", "Observer", "Replicate", "ExperienceLevel", "StartTime", "coll_type", "Area_inch2", "A_P", "TempF", "plantID"), 
                      measure.vars = c("prop_bee", "prop_fly", "prop_lep", "prop_beetle", "prop_unk"))

##Check coll_type
mg_nov_tax <- droplevels(subset(taxa_prop_melt, ExperienceLevel == "Novice"))
coll_mod_tax <- lmer(value ~ offset(log(Area_inch2)) + Cultivar + scale(TempF) + 
                       coll_type + (1|plantID) + variable + coll_type:variable + 
                       Cultivar:variable ,data = mg_nov_tax)
plot(coll_mod_tax)

vif(coll_mod_tax)

r.squaredGLMM(coll_mod_tax)
car::Anova(coll_mod_tax)

taxa_prop <- emmeans(coll_mod_tax, pairwise ~ coll_type + variable, adjust = "tukey", type = "response")
taxa_prop2 <- data.frame(emmeans(coll_mod_tax, ~ coll_type + variable, adjust = "tukey", type = "response"))


taxa_coll_plot <- ggplot(taxa_prop2, aes(x = variable, y = emmean, fill = coll_type)) + geom_bar(stat = "identity", position = "dodge") +
  theme_classic() + 
  geom_errorbar(aes( ymax = taxa_prop2$emmean + taxa_prop2$SE, ymin = taxa_prop2$emmean - taxa_prop2$SE), width = 0, position = position_dodge(width = .9))

taxa_coll_plot


##use this one
tax_mod4_prop <- lmer(value ~ offset(log(Area_inch2)) + Cultivar + scale(TempF) + ExperienceLevel + (1|plantID) + variable + ExperienceLevel:variable + 
                        Cultivar:variable ,data = taxa_prop_melt)


summary(tax_mod4_prop)

r.squaredGLMM(tax_mod4_prop)

car::Anova(tax_mod5_prop)

vif(tax_mod4_prop)

###Taxa observation error
taxa_exp_means <- emmeans(tax_mod5_prop, pairwise ~ ExperienceLevel + variable, adjust = "tukey", type = "response")
taxa_exp_means
#
taxa_means <- data.frame(emmeans(tax_mod5_prop, ~ variable + ExperienceLevel))
#
#
taxa_mean <- ggplot(taxa_means, aes(x = ExperienceLevel, y = emmean, fill = variable)) + geom_bar(stat = "identity", position = "fill") +
  theme_classic() + ylim (0,1)

#+ geom_errorbar(aes(ymax = emmean + SE,
# ymin = emmean - SE), width = 0, position=position_dodge(width=0.9)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
taxa_mean

####Taxa by cultivar

taxa_means2 <- emmeans(tax_mod5_prop, pairwise ~ Cultivar , adjust = "tukey", type = "response")
taxa_means2
#
taxa_means3 <- data.frame(emmeans(tax_mod5_prop, ~ Cultivar + variable))
#
#
taxa_mean_plot2 <- ggplot(taxa_means3, aes(x = Cultivar, y = emmean, fill = variable)) + geom_bar(stat = "identity", position = "fill") +
  theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,1)
taxa_mean_plot2



####Get main taxa for fig 1
taxa_fig1.2 <- mg_climate_merge
taxa_fig.3 <-melt(taxa_fig1.2, id = c("Cultivar", "Date", "Observer", "Replicate", "ExperienceLevel", "StartTime", "coll_type", "Area_inch2", "A_P", "TempF", "plantID"), 
                  measure.vars = c("Bee", "Fly", "Butterfly"))

####


taxa_fig.5 <- aggregate(value ~ variable + Cultivar, taxa_fig.3, FUN = mean)

taxa_fig1.5 <- dcast(taxa_fig.5, Cultivar ~ variable, value.var = "value", fun.aggregate = NULL)


taxa_fig <- data.frame(prop.table(as.matrix(taxa_fig1.5[c(2:4)]), margin = 1))
colnames(taxa_fig)[1] <- "prop_bee"
colnames(taxa_fig)[2] <- "prop_fly"
colnames(taxa_fig)[3] <- "prop_lep"

taxa_fig.2 <- cbind(taxa_fig1.5, taxa_fig)
taxa_fig.2$prop_bee <- ifelse(is.na(taxa_fig.2$prop_bee), 0, taxa_fig.2$prop_bee)
taxa_fig.2$prop_fly <- ifelse(is.na(taxa_fig.2$prop_fly), 0, taxa_fig.2$prop_fly)
taxa_fig.2$prop_lep <- ifelse(is.na(taxa_fig.2$prop_lep), 0, taxa_fig.2$prop_lep)

taxa_fig1 <- taxa_fig.2

taxa_fig1$Main_Taxa <- ifelse(c(taxa_fig1$prop_bee >= 0.3 & c(taxa_fig1$prop_fly < 0.3 & taxa_fig1$prop_lep < 0.3)), "Bee",
                              ifelse(c(taxa_fig1$prop_fly >= 0.3 & c(taxa_fig1$prop_bee < 0.3 & taxa_fig1$prop_lep < 0.3)), "Fly", 
                                     ifelse(c(taxa_fig1$prop_lep >= 0.3 & c(taxa_fig1$prop_bee < 0.3 & taxa_fig1$prop_fly < 0.3)),"Lep", 
                                            ifelse(c(taxa_fig1$prop_bee >= 0.3 & taxa_fig1$prop_fly >= 0.3), "Bee_Fly", 
                                                   ifelse(c(taxa_fig1$prop_bee >= 0.3 & taxa_fig1$prop_lep >= 0.3), "Bee_Lep", "Lep_Fly")))))


