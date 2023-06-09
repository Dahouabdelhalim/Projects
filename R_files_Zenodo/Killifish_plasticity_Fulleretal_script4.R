#read in the packages
library(ggplot2) #for graphics
library(dplyr) #for summarizing
library(gridExtra) #for making a big grid
library(binom) #for testing proportions from an expected ratio

#user sets the working directory
setwd()

#read in the data
#read in the data
census <- read.csv("Killifish_plasticity.csv")
census

#get the total number of males with red,yellow, or orange
census$total_number_color <- census$totalanyorange + census$totalanyred + census$totalanyyellow

#summarize 'any yellow' and 'any red' by dam-sire combinations, sum across water treatments
#calculate the proportion of any red or any yellow offspring pooled across the
#clear and tea-stained families for a give dam-sire combo
census1_sire_dam <- census %>%
  group_by(sire,dam,cross,malecolorpattern) %>%
  summarise(total_anyred = sum(totalanyred),
            total_anyyel=sum(totalanyyellow),
            total_anyorange=sum(totalanyorange),
            total_number_color = sum(total_number_color),
            total_solidyellow=sum(totalsolidyellow),
            total_solidred=sum(totalsolidred))

census1_sire_dam$prop_anyred = census1_sire_dam$total_anyred/census1_sire_dam$total_number_color
census1_sire_dam$prop_anyyel= census1_sire_dam$total_anyyel/census1_sire_dam$total_number_color
census1_sire_dam$prop_solidred = census1_sire_dam$total_solidred/census1_sire_dam$total_number_color
census1_sire_dam$prop_solidyellow = census1_sire_dam$total_solidyellow/census1_sire_dam$total_number_color

#analyze the red sires
#subset the red sires
Red_sires <- subset(census1_sire_dam, (census1_sire_dam$malecolorpattern=="rr"|census1_sire_dam$malecolorpattern=="rb"))
Red_sires #print them out

#create vectors to hold the results
#for each clutch from a red father, we are going to ask if the proportion of sons 
#with any red on them differed from a null expectation of 0%, 50% or 100%
#we will do this by calculating the 95% confidence intervals and asking if they differ from
#0, 0.5 or 1 
consistent_zero_red <- vector() #declare holder vector for test proportion against zero
consistent_fifty_red <- vector() #holder for test proportion against 50%
consistent_all_red <- vector() #holder for test proportion against 100%

#we will also calculate the upper and lower 95% confidence interval for each proportion
upper_CL_red <- vector() #holder for the upper confidence limits
lower_CL_red <- vector() #holder for the lower confidence limits

#declare the same vectors for the analysis of the yellow frequencies for red sires
consistent_zero_yel <- vector() #declare holder vector for test proportion against zero
consistent_fifty_yel <- vector() #holder for test proportion against 50%
consistent_all_yel <- vector() #holder for test proportion against 100%
upper_CL_yel <- vector() #holder for upper 95% confidence limit
lower_CL_yel <- vector() #holder for lower 95% confidence limit

#cycle through the data set for Red_sires
for(i in 1:nrow(Red_sires)){
  holder_red<- NULL #object to hold the result of the binomial test for red
  holder_red <- binom.confint(x=Red_sires$total_anyred[i], n = Red_sires$total_number_color[i], method = "ac") # get the proportion and the confidence intervals
  upper_CL_red[i] <- holder_red$upper #save upper confidence limit
  lower_CL_red[i] <- holder_red$lower #save lower confidence limit
  #test whether the 95% confidence intervals including zero. Is the lower
  #confidence interval lower than zero? 
  consistent_zero_red[i] <- ifelse(lower_CL_red[i]<0, "yes","no") #test zero result
  #test if the 95% confidence interval includes 0.5
  consistent_fifty_red[i] <- ifelse((lower_CL_red[i]<0.5&upper_CL_red[i]>0.5),"yes","no") #test 50% result
  #test if the 95% confidence interval inclues 1
  consistent_all_red[i] <- ifelse(upper_CL_red[i]>1, "yes","no") #test 100% result
  
  holder_yel<- NULL#object to hold the result of the binomial test for yellow
  holder_yel <- binom.confint(x=Red_sires$total_anyyel[i], n = Red_sires$total_number_color[i], method = "ac") #get the proportion of yellow with upper and lower CLs
  upper_CL_yel[i] <- holder_yel$upper #save upper CL
  lower_CL_yel[i] <- holder_yel$lower #save lower CL
  consistent_zero_yel[i] <- ifelse(lower_CL_yel[i]<0, "yes","no") #test if zero in 95% confidence limits
  consistent_fifty_yel[i] <- ifelse((lower_CL_yel[i]<0.5&upper_CL_yel[i]>0.5),"yes","no") #test if 0.5 in 95% confidence limits
  consistent_all_yel[i] <- ifelse(upper_CL_yel[i]>1, "yes","no") #test if 1 in 95% confidence limits
}
#combine all the data together 
Red_sires_sum <- cbind.data.frame(Red_sires, upper_CL_red, lower_CL_red, consistent_zero_red, consistent_fifty_red, consistent_all_red, Red_sires$prop_anyyel, upper_CL_yel, lower_CL_yel, consistent_zero_yel, consistent_fifty_yel, consistent_all_yel)
Red_sires_sum #look at them

#create a term to determine if they deviated from the expected proportions
#combine results to see if any deviations for red proportions
#the logic here is that we ask whether the ratios were consistent with either 0, 0.5, or 1
#if it is consistent with any one of these, then deviate is no
#if it is inconsistent with all three, then deviate is yes. 
Red_sires_sum$deviate_red <- ifelse(((as.numeric(Red_sires_sum$consistent_fifty_red=="yes")+as.numeric(Red_sires_sum$consistent_zero_red=="yes") + as.numeric(Red_sires_sum$consistent_all_red=='yes'))>0),"no","yes")
#combine results to see if any deviations for yellow proportions
#same logic as for red
Red_sires_sum$deviate_yel <- ifelse(((as.numeric(Red_sires_sum$consistent_fifty_yel=="yes")+as.numeric(Red_sires_sum$consistent_zero_yel=="yes") + as.numeric(Red_sires_sum$consistent_all_yel=='yes'))>0),"no","yes")
Red_sires_sum

#add the results to the Red_sires data frame
Red_sires$deviate_red <- ifelse(((as.numeric(Red_sires_sum$consistent_fifty_red=="yes")+as.numeric(Red_sires_sum$consistent_zero_red=="yes") + as.numeric(Red_sires_sum$consistent_all_red=='yes'))>0),"no","yes")
Red_sires$deviate_yel <- ifelse(((as.numeric(Red_sires_sum$consistent_fifty_yel=="yes")+as.numeric(Red_sires_sum$consistent_zero_yel=="yes") + as.numeric(Red_sires_sum$consistent_all_yel=='yes'))>0),"no","yes")
Red_sires$up_cl_red <- Red_sires_sum$upper_CL_red
Red_sires$lower_cl_red <- Red_sires_sum$lower_CL_red
Red_sires$up_cl_yel <- Red_sires_sum$upper_CL_yel
Red_sires$lower_cl_yel <- Red_sires_sum$lower_CL_yel

#----------------Analyze results for the yellow sires-------------------------

#subset the yellow sires
Yellow_sires <- subset(census1_sire_dam, (census1_sire_dam$malecolorpattern=="yy"|census1_sire_dam$malecolorpattern=="yb"))
Yellow_sires

#create vectors to hold the results
#for each clutch from a yellow father, we are going to ask if the proportion of sons 
#with any red on them differed from a null expectation of 0%, 25% or 50%
#we will do this by calculating the 95% confidence intervals and asking if they differ from
#0, 0.25 or 1 
consistent_fifty_red <- vector() # test % red offspring from 50%
consistent_25_red <- vector()  #test % red offspring from 25%
consistent_zero_red <- vector() #test % red offspring from 0
upper_CL_red <- vector() #upper CL on % red
lower_CL_red <- vector() #lower CL on % red

#for each clutch from a yellow father, we are going to ask if the proportion of sons 
#with any red on them differed from a null expectation of 100%, 75% or 50%
#we will do this by calculating the 95% confidence intervals and asking if they differ from
#1, 0.75, or 0.5
consistent_fifty_yel <- vector() #test if yellow offspring from 50%
consistent_75_yel <- vector()  #test if yellow offspring differ from 75%
consistent_all_yel <- vector() #test if yellow offspring differ from 100%
upper_CL_yel <- vector() #hold upper CL on prop yellow
lower_CL_yel <- vector() #hold lower CL on prop yellow

#test all yellow sires for deviations from expected
for(i in 1:nrow(Yellow_sires)){
  holder_red<- NULL
  holder_red <- binom.confint(x=Yellow_sires$total_anyred[i], n = Yellow_sires$total_number_color[i], method = "ac") # calclate prop red
  upper_CL_red[i] <- holder_red$upper #hold upper CL on prop red
  lower_CL_red[i] <- holder_red$lower #hold lower CL on prop red
  consistent_fifty_red[i] <- ifelse((lower_CL_red[i]<0.5&upper_CL_red[i]>0.5),"yes","no") #test if it differs from 50%
  consistent_25_red[i]<- ifelse((lower_CL_red[i]<0.25&upper_CL_red[i]>0.25),"yes","no") #test if it differs from 25%
  consistent_zero_red[i] <- ifelse(lower_CL_red[i]<0, "yes","no") #test if it differs from 0%
  
  holder_yel<- NULL
  holder_yel <- binom.confint(x=Yellow_sires$total_anyyel[i], n = Yellow_sires$total_number_color[i], method = "ac") #calculate proportion of yellow & CLs
  upper_CL_yel[i] <- holder_yel$upper #hold upper CL on yellow
  lower_CL_yel[i] <- holder_yel$lower #hold lower CL on yellow
  consistent_fifty_yel[i] <- ifelse((lower_CL_yel[i]<0.5&upper_CL_yel[i]>0.5),"yes","no") #test if differs from 50%
  consistent_75_yel[i]<- ifelse((lower_CL_yel[i]<0.75&upper_CL_yel[i]>0.75),"yes","no") #test if differs from 75%
  consistent_all_yel[i] <- ifelse(upper_CL_yel[i]>1, "yes","no") #test if differs from 1
  
}

#combine all the data
Yellow_sires_sum <- cbind.data.frame(Yellow_sires[,c(1:5)], Yellow_sires$prop_anyred, upper_CL_red, lower_CL_red,  consistent_fifty_red, consistent_25_red,consistent_zero_red,Yellow_sires$prop_anyyel, upper_CL_yel, lower_CL_yel,  consistent_fifty_yel, consistent_75_yel,consistent_all_yel )

#combine results to see if any deviations for red proportions
Yellow_sires_sum$deviate_red <- ifelse(((as.numeric(Yellow_sires_sum$consistent_fifty_red=="yes")+as.numeric(Yellow_sires_sum$consistent_25_red=="yes") + as.numeric(Yellow_sires_sum$consistent_zero_red=='yes'))>0),"no","yes")
#combine results to see if any deviations for yellow proportions
Yellow_sires_sum$deviate_yel <- ifelse(((as.numeric(Yellow_sires_sum$consistent_fifty_yel=="yes")+as.numeric(Yellow_sires_sum$consistent_75_yel=="yes") + as.numeric(Yellow_sires_sum$consistent_all_yel=='yes'))>0),"no","yes")

#add the summarized results for red to the yellow_sires data frame
Yellow_sires$deviate_red <- ifelse(((as.numeric(Yellow_sires_sum$consistent_fifty_red=="yes")+as.numeric(Yellow_sires_sum$consistent_25_red=="yes") + as.numeric(Yellow_sires_sum$consistent_zero_red=='yes'))>0),"no","yes")
#add the summarized results for yellow to the yellow_sires data frame
Yellow_sires$deviate_yel <- ifelse(((as.numeric(Yellow_sires_sum$consistent_fifty_yel=="yes")+as.numeric(Yellow_sires_sum$consistent_75_yel=="yes") + as.numeric(Yellow_sires_sum$consistent_all_yel=='yes'))>0),"no","yes")
Yellow_sires #look at the data frame
Yellow_sires$up_cl_red <- Yellow_sires_sum$upper_CL_red
Yellow_sires$lower_cl_red <- Yellow_sires_sum$lower_CL_red
Yellow_sires$up_cl_yel <- Yellow_sires_sum$upper_CL_yel
Yellow_sires$lower_cl_yel <- Yellow_sires_sum$lower_CL_yel


#combine the two data frames
combo_data_frame <- rbind.data.frame(Yellow_sires, Red_sires)
combo_data_frame  #look at it

#-----------Make Supplemental table 3-------------------------

Suppl_table_3 <- combo_data_frame[,c(1,2,3,4,9,13,14,11,10,15,16,12)]
Suppl_table_3 


#-----------Look for patterns in deviations from expected with regards to cross-------
#There is no evidence of having more deviations in some cross types than others
#These results are not in the manuscript. 
#make the table
red_result<- table(combo_data_frame$cross, combo_data_frame$deviate_red)
red_result #look at it
#run Chi-squared
chisq.test(red_result, simulate.p.value=TRUE)

#make the table
yel_result<- table(combo_data_frame$cross, combo_data_frame$deviate_yel)
yel_result #look at it
#run Chi-square
chisq.test(yel_result, simulate.p.value=TRUE)

#----------make the graphs--------------------------------------

#create a theme for the graphs
mytheme <- theme(
  axis.text=element_text(size=12),
  axis.title = element_text(size=14),
  legend.text=element_text(size=12),
  legend.position="none"
)

Any_Red_Distribution <- ggplot(combo_data_frame, aes(x=malecolorpattern, y=prop_anyred, color=deviate_red)) + 
  geom_dotplot(binaxis = "y", stackdir="center", stroke=2.5) +
  ylab("Proportion Sons Any Red") + 
  guides(color = "none") +
  scale_color_manual(values=c("black","red")) + 
  scale_x_discrete(labels=c("R-B","R-R","Y-B","Y-Y")) + 
  scale_y_continuous(limits=c(-0.05,1.05), breaks=seq(0,1,0.25)) + 
  xlab("Male Color Pattern") + 
  theme_bw() + mytheme 


Any_Yellow_Distribution <- ggplot(combo_data_frame, aes(x=malecolorpattern, y=prop_anyyel, color=deviate_yel)) + 
  geom_dotplot(binaxis = "y", stackdir="center", stroke=2.5) +
  guides(color = "none") +
  scale_y_continuous(breaks=seq(0,1,0.25)) +  scale_color_manual(values=c("black","yellow")) +
  scale_x_discrete(labels=c("R-B","R-R","Y-B","Y-Y")) + 
  ylab("Proportion Sons Any Yellow") + 
  xlab("Male Color Pattern") + 
  theme_bw() + mytheme

#the graph below is figure 4
Fig4 <- grid.arrange(ncol=1, nrow=2, Any_Red_Distribution, Any_Yellow_Distribution)
Fig4


#-----The next bit of code was used to make supplemental graphs 3 through 6
mytheme <- theme(
  axis.text=element_text(size=12),
  axis.title = element_text(size=12),
  legend.text=element_text(size=12)
)

#get cross means to put on the graph
crossmeans_anyred <- summarise(group_by(census1_sire_dam, cross),
                               mean_any_red = mean(prop_anyred))

#put new order and labels on color pattern
census1_sire_dam$malecolorpattern<- factor(census1_sire_dam$malecolorpattern, levels=c("rb","rr","yb","yy"), labels=c("R-B","R-R","Y-B","Y-Y"))

#put new labels on crosses
labels <- c(spmxspf = "Spring", swmxswf = "Swamp", spmxswf="Spm x Swf", swmxspf="Swm x Spf")


#make supplemental figure 3 Proportion solid red graph--------------------
#get the means to put on the graph 
crossmeans_solidred <- summarise(group_by(census1_sire_dam, cross),
                                 mean_solid_red = mean(prop_solidred))


Suppl_fig_3 <-ggplot(census1_sire_dam, aes(x=malecolorpattern, y=prop_solidred)) + 
  geom_dotplot(aes(x=malecolorpattern,group=malecolorpattern),
               binaxis = "y", stackdir="center", dotsize=0.8) + 
  stat_summary(aes(x=malecolorpattern,group=malecolorpattern),fun="mean", geom='point', shape=95, size=5, position=position_nudge(0.2)) + 
  stat_summary(aes(x=malecolorpattern,group=malecolorpattern), fun.data = mean_se, geom = "errorbar", size=0.5, width=0.05, position=position_nudge(0.2)) +
  ylab("Proportion Sons Solid Red") + 
  facet_wrap(~cross, labeller = labeller(cross=labels))+
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.02,1.02)) + 
  xlab("Male Color Pattern") + 
  geom_text(
    data    = crossmeans_solidred,
    mapping = aes(x = 3, y = 0.9, label = paste("mean =",format(round(mean_solid_red,2), nsmall=2))),
  )+
  theme_bw() + 
  mytheme

Suppl_fig_3


#-----------Supplemental figure 4 Proportion solid yellow---------------------------------------------------
#get cross means to put on the graphs
crossmeans_solidyellow <- summarise(group_by(census1_sire_dam, cross),
                                    mean_solid_yellow = mean(prop_solidyellow))

Suppl_fig_4 <-ggplot(census1_sire_dam, aes(x=malecolorpattern, y=prop_solidyellow)) + 
  geom_dotplot(aes(x=malecolorpattern,group=malecolorpattern),
               binaxis = "y", stackdir="center", dotsize=0.8) + 
  stat_summary(aes(x=as.numeric(malecolorpattern)+0.2,group=malecolorpattern),fun="mean", geom='point', shape=95, size=5, position=position_nudge(0.2)) + 
  stat_summary(aes(x=as.numeric(malecolorpattern)+0.2,group=malecolorpattern),fun.data = mean_se, geom = "errorbar", size=0.5, width=0.05, position=position_nudge(0.2)) +
  ylab("Proportion Sons Solid Yellow") + 
  facet_wrap(~cross, labeller = labeller(cross=labels))+
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.02,1.02)) + 
  xlab("Male Color Pattern") + 
  geom_text(
    data    = crossmeans_solidyellow,
    mapping = aes(x = 3, y = 0.9, label = paste("mean =",round(mean_solid_yellow,2))))+
  theme_bw() + 
  mytheme

Suppl_fig_4 


#---------- Make Supplemental figure 5 -------- any red graph --------------------
Suppl_fig_5 <-ggplot(census1_sire_dam, aes(x=malecolorpattern, y=prop_anyred)) + 
  geom_dotplot(aes(x=malecolorpattern,group=malecolorpattern),
               binaxis = "y", stackdir="center", dotsize=0.8) + 
  stat_summary(aes(x=malecolorpattern,group=malecolorpattern),fun="mean", geom='point', shape=95, size=5, position=position_nudge(0.2)) + 
  stat_summary(aes(x=malecolorpattern,group=malecolorpattern),fun.data = mean_se, geom = "errorbar", size=0.5, width=0.05, position=position_nudge(0.2)) +
  ylab("Proportion Sons Any Red") + 
  facet_wrap(~cross, labeller = labeller(cross=labels))+
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.02,1.02)) + 
  xlab("Male Color Pattern") + 
  geom_text(
    data    = crossmeans_anyred,
    mapping = aes(x = 3, y = 0.9, label = paste("mean =",round(mean_any_red,2)))
  )+
  
  theme_bw() + 
  mytheme

Suppl_fig_5 # note that the placement of the means varies depending on size of graph


#----------Make Supplemental figure 6-----Any yellow graph------------------------
#get means to put on the graph
crossmeans_yel <- summarise(group_by(census1_sire_dam, cross),
                            mean_any_yel = mean(prop_anyyel))


Suppl_fig_6 <- ggplot(census1_sire_dam, aes(x=malecolorpattern, y=prop_anyyel)) + 
  geom_dotplot(aes(x=malecolorpattern,group=malecolorpattern),
               binaxis = "y", stackdir="center", dotsize=0.8) + 
  stat_summary(aes(x=malecolorpattern,group=malecolorpattern),fun="mean", geom='point', shape=95, size=5, position=position_nudge(0.2)) + 
  stat_summary(aes(x=malecolorpattern,group=malecolorpattern),fun.data = mean_se, geom = "errorbar", size=0.5, width=0.05, position=position_nudge(0.2)) +
  ylab("Proportion Sons Any Yellow") + 
  facet_wrap(~cross, labeller = labeller(cross=labels))+
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(-0.02,1.02)) + 
  xlab("Male Color Pattern") + 
  theme_bw() + 
  geom_text(
    data    = crossmeans_yel,
    mapping = aes(x = 3, y = 0.1, label = paste("mean =",format(round(mean_any_yel,2), nsmall=2)))
  )+
  
  mytheme

Suppl_fig_6


