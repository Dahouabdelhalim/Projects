#Ms Cc Early and Developmental Heat Shock Experiments

#Figures for Manuscript

#load libraries
library(readr)
library(ggplot2)
library(tidyr)
library(Rmisc)
library(dplyr)
library(viridis)
library(cowplot)
library(extrafont)


#load data----------

ehs <- read_csv("data/Ms+Cc_EHS_complete_clean.csv")



ehs_lng <- read_csv("data/Ms+Cc_EHS_incomplete_clean_long.csv", 
                    col_types = cols(hs.num = col_factor(levels = c("0", 
                                                                    "1", "2", "3", "4"))))


ehs_wowdis <- read_csv("data/EHS_mongos_for_diss.csv", 
                       col_types = cols(hs.temp = col_factor(levels = c("40",  "42")), 
                                        hs.num = col_factor(levels = c("1", "2", "3", "4"))))


rhs <- read_csv("data/Ms+Cc_RHS_complete_clean.csv")



rhs_lng <- read_csv("data/Ms+Cc_RHS_comp_clean_long.csv")



theme_set(theme_classic())

#data cleaning-------------------

#EHS data cleaning

#change hs.temp 0 to 35
ehs$hs.temp <- ifelse(ehs$hs.temp==0, 35, ehs$hs.temp)

#make hs.temp a factor for colors when plotting
ehs$hs.temp <- factor(ehs$hs.temp)

#remove wanderers
ehs <- subset(ehs, class!="wand")


#remove individuals that were not dissected (lost in freezer or with missing data)
ehs$num.unem[is.na(ehs$num.unem)] <- 0
ehs <- subset(ehs, num.unem >= 0)


#create an all treatment column that combines hs.temp and hs.num
ehs <- unite(ehs, hs_treat, hs.temp, hs.num, sep="_", remove = FALSE)



#EHS dissected WOWE data cleaning

#remove individual that couldn't be found to be dissected
ehs_wowdis$diss.wasp[is.na(ehs_wowdis$diss.wasp)] <- 0.5
ehs_wowdis <- subset(ehs_wowdis, diss.wasp!=0.5)




#RHS data cleaning

#remove low late heat shocks
rhs <- subset(rhs, hs.cal!="low")

#remove individuals that died
rhs <- subset(rhs, died.bf5==0)

#drop rows with NAs in tot.load--1 that did not have num_em recorded, 2 that got lost in freezer and
#not dissected
rhs <- drop_na(rhs, tot.load)

#remove individuals with load greater than 300
rhs <- subset(rhs, tot.load<300)

#remove wanderers
rhs$date.wand.j[is.na(rhs$date.wand.j)] <- 0
rhs <- subset(rhs, date.wand.j==0)


#creating column "class"--em or mongo
rhs$date.em.j[is.na(rhs$date.em.j)]<-0
rhs$date.cull.j[is.na(rhs$date.cull.j)]<-0
rhs$num.em[is.na(rhs$num.em)]<-0


rhs$class<-ifelse(rhs$date.em.j>0 | rhs$num.em>0, "em",
                  ifelse(rhs$date.cull.j>0, "mongo", "unk"))

#subset out any with class "unk"
rhs <- subset(rhs, class!="unk")


#order shock stage as a factor
rhs$shock.stage <- factor(rhs$shock.stage, levels = c("control", "early", "mid", "late"))


#make long dataframe of wasp mass, number and sex
rhs_wam <- rhs %>% gather(sex, ad.mass, ind.fem.mass, ind.male.mass)
rhs_wam$sex<-gsub("ind.fem.mass", "female", rhs_wam$sex)
rhs_wam$sex<-gsub("ind.male.mass", "male", rhs_wam$sex)

rhs_wnum <- rhs %>% gather(sex, num_wasp, fem.ecl, male.ecl)
rhs_wnum$sex <- gsub("fem.ecl", "female", rhs_wnum$sex)
rhs_wnum$sex <- gsub("male.ecl", "male", rhs_wnum$sex)

rhs_wnum <- select(rhs_wnum, id, shock.stage, sex, num_wasp)

rhs_wlng <- merge(rhs_wam, rhs_wnum, by=c("id", "shock.stage", "sex"))


#Fig 1: proportion of WOWE/em for each treatment; both expts----------------------

#RHS table

#determining num of each class in each treatment, and the total number in each treatment
n.table<-table(rhs$shock.stage, rhs$class)

#converting tables to data frames
n.table <- as.data.frame.matrix(n.table) 

#make a shock.stage column
n.table <- data.frame(shock.stage = row.names(n.table), n.table)

#create tot.n column
n.table$tot.n <- n.table$em + n.table$mongo

n.table <- gather(n.table, class, n, em, mongo)

#create proportion column
n.table$prop <- n.table$n / n.table$tot.n

#level shock stage correctly
n.table$shock.stage <- factor(n.table$shock.stage, levels=c("control", "early", "mid", "late"))


#plot proportion of each class for each shock stage
rhs_class_plot <- ggplot(n.table, aes(x = shock.stage, y = prop, fill = class))
rhs_class_plot <- rhs_class_plot + geom_bar(position="fill", stat="identity",
                                            width = .7
) + scale_fill_manual(values = c("#95D840", "#1F9F88"),
                    breaks = c("em", "mongo"),
                    labels = c("Emergence", "WOWE"),
                    name="Outcome"
) + labs(x="Shock Stage", y="Proportion"
) + scale_x_discrete(labels=c("Control", "Early", "Middle", "Late")                   
) + theme(axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2.5, "mm"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none")

rhs_class_plot





#EHS plot

#create a new class description: WOWE with wasp larvae--no emergence, but wasp larvae found upon dissection
#may not keep this in final figure, since only a subset of WOWEs were dissected
ehs$tot.unem[is.na(ehs$tot.unem)] <- 0

ehs$class <- ifelse(ehs$class=="mongo" & ehs$tot.unem > 0, "mwwl", ehs$class)

#make a column that includes mwwl as mongos--since only a subset were dissected
ehs$class2 <- ifelse(ehs$class=="mongo" | ehs$class=="mwwl", "mongo", "em")


#determining num of each class in each treatment, and the total number in each treatment
ehs.n.table<-table(ehs$hs_treat, ehs$class2)


#converting tables to data frames
ehs.n.table <- as.data.frame.matrix(ehs.n.table) 

#make a shock.stage column
ehs.n.table <- data.frame(hs_treat = row.names(ehs.n.table), ehs.n.table)

#create tot.n column
ehs.n.table$tot.n <- ehs.n.table$em + ehs.n.table$mongo 

ehs.n.table <- gather(ehs.n.table, class, n, em, mongo)

#create proportion column
ehs.n.table$prop <- ehs.n.table$n / ehs.n.table$tot.n

#create separate hs.temp and hs.treat columns
ehs.n.table <- ehs.n.table %>% separate(hs_treat, c("hs.temp", "hs.num"), sep = "_")

#plot 35 and hs treatments separately
ehsnt_35 <- subset(ehs.n.table, hs.temp==35)
ehsnt_hs <- subset(ehs.n.table, hs.temp!=35)


#plot proportion of each class for each heat shock temp and number
ehs.class.plot35 <- ggplot(ehsnt_35, aes(x=hs.num, y=prop, fill=class))
ehs.class.plot35 <- ehs.class.plot35 + geom_bar(position="fill", stat="identity",
                                                width = .9
) + scale_fill_manual(values=c("#95D840", "#1F9F88"),
                    breaks=c("em", "mongo"),
                    labels=c("Emergence", "WOWE"),
                    name="Outcome"
) + labs(x="Days in Heat Shock", y="Proportion"
) + facet_wrap(~hs.temp
) + theme(strip.text = element_text(size=15),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          legend.position = "none")

ehs.class.plot35


#plot for legend
for_legend <- ggplot(ehsnt_hs, aes(x=hs.num, y=prop, fill=class))
for_legend <- for_legend + geom_bar(position="fill", stat="identity"
) + scale_fill_manual(values=c("#95D840", "#1F9F88"),
                      breaks=c("em", "mongo"),
                      labels=c("WE", "WOWE"),
                      name="Emergence"
) + labs(x="Days in Heat Shock", y="Proportion"
) + facet_wrap(~hs.temp
) + theme(axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_blank(),
          axis.ticks.x = element_line(colour = 'black', size = 1),
          axis.ticks.length.x = unit(2.5, "mm"),
          axis.ticks.y = element_blank(), 
          axis.text.x = element_text(size = 15),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_blank(),
          legend.text=element_text(size=12),
          legend.title=element_text(size=16),
          legend.background = element_rect(color="black",linetype="solid",size=1))


#get legend for ehs class plots
ehs_legend <- get_legend(for_legend)



ehs.class.ploths <- ggplot(ehsnt_hs, aes(x=hs.num, y=prop, fill=class))
ehs.class.ploths <- ehs.class.ploths + geom_bar(position="fill", stat="identity",
                                                width = .9
) + scale_fill_manual(values=c("#95D840", "#1F9F88"),
                      breaks=c("em", "mongo"),
                      labels=c("Emergence", "WOWE"),
                      name="Outcome"
) + labs(x="Days in Heat Shock", y="Proportion"
) + facet_wrap(~hs.temp
) + theme(strip.text = element_text(size=15),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_blank(),
          axis.ticks.x = element_line(colour = 'black', size = 1),
          axis.ticks.length.x = unit(2.5, "mm"),
          axis.ticks.y = element_blank(), 
          axis.text.x = element_text(size = 15),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_blank(),
          legend.position = "none")


ehs.class.ploths





#combine the ehs plots
ehs.class.plot <- plot_grid(ehs.class.plot35, ehs.class.ploths, align = "h", rel_widths = c(1 , 3.5))

ehs.class.plot


#combine rhs and legend
rhs_class_plot_leg <- plot_grid(rhs_class_plot, ehs_legend, rel_widths = c(2, 1))

rhs_class_plot_leg



class.plot <- plot_grid(rhs_class_plot_leg, ehs.class.plot,
                        nrow = 2, labels = c("A", "B"))

class.plot








#Fig 2: Bar plot of the mn number of wasp larvae by stage: immat 2nd, mat 2nd, emerged, eclosed-----

#EHS select to columns that are useful for this plot
ehs_ws <- select(ehs, id, hs.temp, hs.num, num.mat2, num.em, num.ecl, load)

#make all NAs = 0
ehs_ws[is.na(ehs_ws)] <- 0

#calculate the number that died between emergence and eclosion
ehs_ws$died_em <- ehs_ws$num.em - ehs_ws$num.ecl

#calculate the number of wasps accounted for by mat2, em and ecl
ehs_ws$calc <- ehs_ws$num.mat2 + ehs_ws$died_em + ehs_ws$num.ecl

#calculate the immat 2 load
ehs_ws$num.immat2 <- ehs_ws$load - ehs_ws$calc

#make a long data set
ehs_lws <- gather(ehs_ws, stage, num, num.immat2, num.mat2, died_em, num.ecl)

#make stage a factor so that the levels are correct
ehs_lws$stage <- factor(ehs_lws$stage, levels = c("num.ecl", "died_em", "num.mat2", "num.immat2"))



#RHS--select to only necessary columns
rhs_ws <- select(rhs, id, shock.stage, num.unem.mat, num.unem.im, num.em, num.ecl, tot.load)

#make all NAs = 0
rhs_ws[is.na(rhs_ws)] <- 0

#calculate the number that died between emergence and eclosion
rhs_ws$died_em <- rhs_ws$num.em - rhs_ws$num.ecl


#create long dataset
rhs_lws <- gather(rhs_ws, stage, num, num.unem.mat, num.unem.im, died_em, num.ecl)

#make stage a factor so the levels are in the correct order
rhs_lws$stage <- factor(rhs_lws$stage, levels = c("num.ecl", "died_em", "num.unem.mat", "num.unem.im"))



#EHS PLOT

#calculate means and variances
stage_sum <- summarySE(ehs_lws, measurevar = "num",
                       groupvars = c("hs.temp", "hs.num", "stage"),
                       na.rm = TRUE)

stage_sum


nui_stmn_35 <- stage_sum[4, 5]
num_stmn_35 <- stage_sum[3, 5] + nui_stmn_35
de_stmn_35 <- stage_sum[2, 5] + num_stmn_35
ne_stmn_35 <- stage_sum[1, 5] + de_stmn_35

f40.1_nui_stmn <- stage_sum[8, 5]
f40.1_num_stmn <- stage_sum[7, 5] + f40.1_nui_stmn
f40.1_de_stmn <- stage_sum[6, 5] + f40.1_num_stmn
f40.1_ne_stmn <- stage_sum[5, 5] + f40.1_de_stmn

f40.2_nui_stmn <- stage_sum[12, 5]
f40.2_num_stmn <- stage_sum[11, 5] + f40.2_nui_stmn
f40.2_de_stmn <- stage_sum[10, 5] + f40.2_num_stmn
f40.2_ne_stmn <- stage_sum[9, 5] + f40.2_de_stmn

f40.3_nui_stmn <- stage_sum[16, 5]
f40.3_num_stmn <- stage_sum[15, 5] + f40.3_nui_stmn
f40.3_de_stmn <- stage_sum[14, 5] + f40.3_num_stmn
f40.3_ne_stmn <- stage_sum[13, 5] + f40.3_de_stmn

f40.4_nui_stmn <- stage_sum[20, 5]
f40.4_num_stmn <- stage_sum[19, 5] + f40.4_nui_stmn
f40.4_de_stmn <- stage_sum[18, 5] + f40.4_num_stmn
f40.4_ne_stmn <- stage_sum[17, 5] + f40.4_de_stmn

f42.1_nui_stmn <- stage_sum[24, 5]
f42.1_num_stmn <- stage_sum[23, 5] + f42.1_nui_stmn
f42.1_de_stmn <- stage_sum[22, 5] + f42.1_num_stmn
f42.1_ne_stmn <- stage_sum[21, 5] + f42.1_de_stmn

f42.2_nui_stmn <- stage_sum[28, 5]
f42.2_num_stmn <- stage_sum[27, 5] + f42.2_nui_stmn
f42.2_de_stmn <- stage_sum[26, 5] + f42.2_num_stmn
f42.2_ne_stmn <- stage_sum[25, 5] + f42.2_de_stmn

f42.3_nui_stmn <- stage_sum[32, 5]
f42.3_num_stmn <- stage_sum[31, 5] + f42.3_nui_stmn
f42.3_de_stmn <- stage_sum[30, 5] + f42.3_num_stmn
f42.3_ne_stmn <- stage_sum[29, 5] + f42.3_de_stmn

f42.4_nui_stmn <- stage_sum[36, 5]
f42.4_num_stmn <- stage_sum[35, 5] + f42.4_nui_stmn
f42.4_de_stmn <- stage_sum[34, 5] + f42.4_num_stmn
f42.4_ne_stmn <- stage_sum[33, 5] + f42.4_de_stmn


stage_sum$stckd_mn <- c(ne_stmn_35, de_stmn_35, num_stmn_35, nui_stmn_35,
                        f40.1_ne_stmn, f40.1_de_stmn, f40.1_num_stmn, f40.1_nui_stmn,
                        f40.2_ne_stmn, f40.2_de_stmn, f40.2_num_stmn, f40.2_nui_stmn,
                        f40.3_ne_stmn, f40.3_de_stmn, f40.3_num_stmn, f40.3_nui_stmn,
                        f40.4_ne_stmn, f40.4_de_stmn, f40.4_num_stmn, f40.4_nui_stmn,
                        f42.1_ne_stmn, f42.1_de_stmn, f42.1_num_stmn, f42.1_nui_stmn,
                        f42.2_ne_stmn, f42.2_de_stmn, f42.2_num_stmn, f42.2_nui_stmn,
                        f42.3_ne_stmn, f42.3_de_stmn, f42.3_num_stmn, f42.3_nui_stmn,
                        f42.4_ne_stmn, f42.4_de_stmn, f42.4_num_stmn, f42.4_nui_stmn)


#subset summary into control and hs treatments
stage_sum35 <- subset(stage_sum, hs.temp==35)
stage_sum_hw <- subset(stage_sum, hs.temp!=35)


#bar plot of stage--control treatment
waspst_plot35 <- ggplot(stage_sum35, aes(x=as.character(hs.num), y=num, fill=stage))
waspst_plot35 <- waspst_plot35 + geom_bar(stat = "identity"
) + geom_errorbar(aes(ymin = stckd_mn - se, ymax = stckd_mn + se),
                  width = .5
) + scale_fill_manual(values = c("#FDCD31", "#F17020", "#E8226F", "#6F22A8"),
                      breaks = c("num.ecl", "died_em", "num.mat2", "num.immat2"),
                      labels = c("Eclosed", "Emerged", "Mature 2nd", "Immature 2nd"),
                      name = "Stage Survived"
) + scale_y_continuous(limits = c(0, 130),
                       breaks = seq(0, 130, by=20)
) + labs(x="Days in Heat Wave", y="Number of Parasitoids"
) + facet_wrap(~hs.temp
) + theme(strip.text.x = element_text(size = 15),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          legend.position = "none")

waspst_plot35


#bar plot of stage--heat wave treatments
waspst_plot_hw <- ggplot(stage_sum_hw, aes(x=as.character(hs.num), y=num, fill=stage))
waspst_plot_hw <- waspst_plot_hw + geom_bar(stat = "identity"
) + geom_errorbar(aes(ymin = stckd_mn - se, ymax = stckd_mn + se),
                  width = .5
) + scale_fill_manual(values = c("#FDCD31", "#F17020", "#E8226F", "#6F22A8"),
                      breaks = c("num.ecl", "died_em", "num.mat2", "num.immat2"),
                      labels = c("Eclosed", "Emerged", "Mature 2nd", "Immature 2nd"),
                      name = "Stage Survived"
) + scale_y_continuous(limits = c(0, 130),
                       breaks = seq(0, 130, by=20)
) + labs(x="Days in Early Heat Shock", y="Number of Parasitoids"
) + facet_wrap(~hs.temp
) + theme(strip.text.x = element_text(size = 15),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_blank(),
          legend.position = "none")

waspst_plot_hw


#combine control and heat wave plots
ehs_wasp_stage_fig <- plot_grid(waspst_plot35, waspst_plot_hw,
                                nrow = 1, rel_widths = c(1, 3),
                                align = "h")
ehs_wasp_stage_fig



#RHS PLOT

rhs_stage_sum <- summarySE(rhs_lws, measurevar = "num",
                           groupvars = c("shock.stage", "stage"),
                           na.rm = TRUE)

#make column adding means to se so that error bars appear in the correct place on stacked bars
con_nui_stmn <- rhs_stage_sum[4, 4]
con_num_stmn <- rhs_stage_sum[3, 4] + con_nui_stmn
con_de_stmn <- rhs_stage_sum[2, 4] + con_num_stmn
con_ne_stmn <- rhs_stage_sum[1, 4] + con_de_stmn

erly_nui_stmn <- rhs_stage_sum[8, 4]
erly_num_stmn <- rhs_stage_sum[7, 4] + erly_nui_stmn
erly_de_stmn <- rhs_stage_sum[6, 4] + erly_num_stmn
erly_ne_stmn <- rhs_stage_sum[5, 4] + erly_de_stmn

mid_nui_stmn <- rhs_stage_sum[12, 4]
mid_num_stmn <- rhs_stage_sum[11, 4] + mid_nui_stmn
mid_de_stmn <- rhs_stage_sum[10, 4] + mid_num_stmn
mid_ne_stmn <- rhs_stage_sum[9, 4] + mid_de_stmn

lte_nui_stmn <- rhs_stage_sum[16, 4]
lte_num_stmn <- rhs_stage_sum[15, 4] + lte_nui_stmn
lte_de_stmn <- rhs_stage_sum[14, 4] + lte_num_stmn
lte_ne_stmn <- rhs_stage_sum[13, 4] + lte_de_stmn


rhs_stage_sum$stckd_mn <- c(con_ne_stmn, con_de_stmn, con_num_stmn, con_nui_stmn,
                            erly_ne_stmn, erly_de_stmn, erly_num_stmn, erly_nui_stmn,
                            mid_ne_stmn, mid_de_stmn, mid_num_stmn, mid_nui_stmn,
                            lte_ne_stmn, lte_de_stmn, lte_nui_stmn, lte_num_stmn)



rhs_waspst_plot <- ggplot(rhs_stage_sum, aes(x=shock.stage, y=num, fill=stage))
rhs_waspst_plot <- rhs_waspst_plot + geom_bar(stat = "identity",
                                              width = .7
) + geom_errorbar(aes(ymin=stckd_mn - se, ymax=stckd_mn + se),
                  width=.5
) + scale_fill_manual(values = c("#FDCD31", "#F17020", "#E8226F", "#6F22A8"),
                      breaks = c("num.ecl", "died_em", "num.unem.mat", "num.unem.im"),
                      labels = c("Eclosed", "Emerged", "Mature 2nd", "Immature 2nd"),
                      name = "Final Stage"
) + scale_x_discrete(breaks = c("control", "early", "mid", "late"),
                     labels = c("Control", "Early", "Middle", "Late")
) + scale_y_continuous(limits = c(0, 130),
                       breaks = seq(0, 130, by=20)
) + labs(x="Heat Shock Stage", y="Number of Parasitoids"
) + theme(axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size=15))

rhs_waspst_plot


#combine rhs and ehs wasp stage bar plots into one figure
wasp_stage_fig <- plot_grid(rhs_waspst_plot, ehs_wasp_stage_fig,
                            nrow = 2, align = "h",
                            labels = c("A", "B"),
                            rel_widths = c(1, 3),
                            rel_heights = c(1, 1.2))

wasp_stage_fig



#Fig 3: RHS figures of wasp survival and mass------------------

#RHS wasp survival to eclosion


#calculate mean proportion of wasps survivng to eclosion by shock stage and sex
prop_ecl_sum <- summarySE(rhs, measurevar = "ps.ecl",
                          groupvars = "shock.stage",
                          na.rm = TRUE)
prop_ecl_sum

#make group for connecting lines
prop_ecl_sum$group <- ifelse(prop_ecl_sum$shock.stage=="control", 1, 2)


#plot mean proportion of wasp survival to eclosion by shock stage 

#plot for combining
propecl_plot <- ggplot(prop_ecl_sum, aes(x=shock.stage, y=ps.ecl, group=group))
propecl_plot <- propecl_plot + geom_point(size=5, stroke=2, shape=15, color="#666666"
) + geom_line(size=2, color="#666666"
) + geom_errorbar(aes(ymin = ps.ecl - se, ymax = ps.ecl + se),
                  width=.5, size=1.2, color="#666666"
) + scale_x_discrete(labels = c("Control", "Early", "Middle", "Late")
) + labs(x="Shock Stage", y="Proportion Surviving to Eclosion"
) + theme(axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          legend.position = "none")

propecl_plot



#calculate mean wasp mass by shock stage and sex
wadm_sum <- summarySE(rhs_wlng, measurevar = "ad.mass",
                      groupvars = c("shock.stage", "sex"),
                      na.rm = TRUE)




#plot mean proportion of wasp survival to eclosion by shock stage and sex
wadmass_plot <- ggplot(wadm_sum, aes(x=shock.stage, y=ad.mass, group=sex, color=sex))
wadmass_plot <- wadmass_plot + geom_point(aes(shape=sex, size=N),
                               stroke=2
) + geom_line(size=2
) + geom_errorbar(aes(ymin = ad.mass - se, ymax = ad.mass + se),
                  width=.5, size=1.2
) + scale_color_manual(values = c("#000000", "#E69F00"),
                       labels = c("Female", "Male"),
                       name = "Sex"
) + scale_shape_manual(values = c(16, 17),
                       labels = c("Female", "Male"),
                       name = "Sex"
) + scale_x_discrete(labels = c("Control", "Early", "Middle", "Late")
) + labs(x="Shock Stage", y="Adult Wasp Mass [mg]"
) + theme(axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          legend.position = "none")

wadmass_plot



#for legend
wadmass_plot_leg <- ggplot(wadm_sum, aes(x=shock.stage, y=ad.mass, group=sex, color=sex))
wadmass_plot_leg <- wadmass_plot_leg + geom_point(aes(shape=sex),
                                          stroke=2, size=5
) + geom_line(size=2
) + geom_errorbar(aes(ymin = ad.mass - se, ymax = ad.mass + se),
                  width=.5, size=1.2
) + scale_color_manual(values = c("#000000", "#E69F00"),
                       labels = c("Female", "Male"),
                       name = "Sex"
) + scale_shape_manual(values = c(16, 17),
                       labels = c("Female", "Male"),
                       name = "Sex"
) + scale_x_discrete(labels = c("Control", "Early", "Middle", "Late")
) + labs(x="Shock Stage", y="Adult Wasp Mass [mg]"
) + theme(axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          legend.background = element_rect(color="black",linetype="solid",size=1))

wadmass_plot_leg


wadmass_leg <- get_legend(wadmass_plot_leg)


#combine RHS wasp survival and mass plots
rhs_wasp_survmass <- plot_grid(propecl_plot, wadmass_plot, 
                               align = "vh", nrow = 2)

rhs_wasp_survmass


rhs_wasp_survmass_leg <- plot_grid(rhs_wasp_survmass, wadmass_leg, 
                                   rel_widths = c(1, .2))
rhs_wasp_survmass_leg







#Suppl: RHS wasp adult mass by shock stage, sex and load----------------

#remove early treatment to make the plot look nicer
rhs_wwlng <- subset(rhs_wlng, shock.stage!="early")

#create a labeller for the facet wrap titles
shock_labs <- c("control" = "Control",
                "mid" = "Middle",
                "late" = "Late")


admass_plot <- ggplot(rhs_wwlng, aes(x=tot.load, y=ad.mass, group=sex, color=sex))
admass_plot + geom_point(aes(shape=sex),
                         size=4
) + geom_smooth(method = lm, se=FALSE, size=1,
) + scale_color_manual(values = c("#000000", "#E69F00"),
                       labels = c("Female", "Male"),
                       name = "Sex"
) + scale_shape_manual(values = c(16, 17),
                       labels = c("Female", "Male"),
                       name = "Sex"
) + labs(x="Parasitoid Load", y="Adult Mass [mg]"
) + facet_wrap(~shock.stage, labeller = as_labeller(shock_labs))


#Suppl: EHS dissected WOWE figures-----

#make table of how many dissected WOWEs have wasp larvae

#Table of number of WOWEs by treatment that had wasp larvae or not
wlf <- dplyr::count(ehs_wowdis, hs.temp, hs.num, diss.wasp)


#make data frame for plotting
wlf$tot.n <- c(6, 6, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7)
wlf$prop <- wlf$n / wlf$tot.n
wlf$diss.wasp <- factor(wlf$diss.wasp)


#Bar plot of the proportion of dissected WOWEs that had wasp larvae
prop_wwasp_plot <- ggplot(wlf, aes(x=hs.num, y=prop, fill=diss.wasp))
prop_wwasp_plot <- prop_wwasp_plot + geom_bar(position = "fill", stat = "identity"
) + scale_fill_manual(values = c("#1F9F88", "#453275"),
                      breaks = c(0, 1),
                      name = "Wasp Larvae"
) + labs(x="Days in Heat Shock", y="Proportion"
) + facet_wrap(~hs.temp
) + theme(strip.text = element_text(size = 18),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          legend.position = "none")

prop_wwasp_plot



#plot of mean mass at culling by days in heat shock and whether wasp larvae were found

#make diss.wasp a factor
ehs_wowdis$diss.wasp <- factor(ehs_wowdis$diss.wasp)

#convert mass.end to grams
ehs_wowdis$mass.end.g <- ehs_wowdis$mass.end / 1000


#find mean mass.end for WOWEs with and without wasp larvae
wlf_sum <- summarySE(ehs_wowdis, measurevar = "mass.end.g",
                     groupvars = c("hs.temp", "hs.num", "diss.wasp"),
                     na.rm = TRUE)

wlf_sum




#for combining
mn_mssend_disw_plot <- ggplot(wlf_sum, aes(x=hs.num, y=mass.end.g, group=diss.wasp, color=diss.wasp))
mn_mssend_disw_plot <- mn_mssend_disw_plot + geom_point(aes(shape=diss.wasp),
                                                        size=6
) + geom_line(size=1.5
) + geom_errorbar(aes(ymin=mass.end.g - sd, ymax=mass.end.g + sd),
                  size=1, width=.5
) + scale_color_manual(values = c("#1F9F88", "#453275"),
                       breaks = c(0, 1),
                       name = "Wasp Larvae"
) + scale_shape_manual(values = c(16, 17),
                       breaks = c(0, 1),
                       name = "Wasp Larvae"
) + labs(x="Days in Heat Shock", y="Mass at Cull [g]"
) + facet_wrap(~hs.temp
) + theme(strip.text = element_text(size = 18),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size = 18),
          legend.position = "none")

mn_mssend_disw_plot



#for legend
mn_mssend_disw_plot_leg <- ggplot(wlf_sum, aes(x=hs.num, y=mass.end.g, group=diss.wasp, color=diss.wasp))
mn_mssend_disw_plot_leg <- mn_mssend_disw_plot_leg + geom_point(aes(shape=diss.wasp),
                                                        size=6
) + geom_line(size=1.5
) + geom_errorbar(aes(ymin=mass.end.g - sd, ymax=mass.end.g + sd),
                  size=1, width=.5
) + scale_color_manual(values = c("#1F9F88", "#453275"),
                       breaks = c(0, 1),
                       name = "Wasp Larvae \\nPresent"
) + scale_shape_manual(values = c(16, 17),
                       breaks = c(0, 1),
                       name = "Wasp Larvae \\nPresent"
) + labs(x="Days in Heat Shock", y="Proportion"
) + facet_wrap(~hs.temp
) + theme(strip.text = element_text(size = 18),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.key.width = unit(1, "cm"),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15))


diss_wasp_leg <- get_legend(mn_mssend_disw_plot_leg)



#Violin and scatter plot of WOWE mass at end, highlighting the ones that I dissected

#subset ehs to only wowes
ehs_wowe <- subset(ehs, class=="mongo")

ehs_wowe$hs.num <- factor(ehs_wowe$hs.num)

#subset dissected wowes to make plotting and color easier
ehs_wowdis_0 <- subset(ehs_wowdis, diss.wasp==0)
ehs_wowdis_1 <- subset(ehs_wowdis, diss.wasp==1)

#calculate mass.end in grams
ehs_wowe$mass.end.g <- ehs_wowe$mass.end / 1000
ehs_wowdis$mass.end.g <- ehs_wowdis$mass.end / 1000
ehs_wowdis_0$mass.end.g <- ehs_wowdis_0$mass.end / 1000
ehs_wowdis_1$mass.end.g <- ehs_wowdis_1$mass.end / 1000


#For combining
wowe_mssend_viol_plot <- ggplot(ehs_wowe, aes(x=hs.num, y=mass.end.g, color=hs.num))
wowe_mssend_viol_plot <- wowe_mssend_viol_plot + geom_violin(aes(fill=hs.num),
                                     alpha=.4
) + scale_fill_manual(values=c("#9105A5", "#CD4479", "#EC9221", "#E8D900"),
                    breaks=c(1,2,3,4),
                    name="Days in Heat Shock"
) + scale_color_manual(values=c("#9105A5", "#CD4479", "#EC9221", "#E8D900"),
                     breaks=c(1,2,3,4),
                     name="Days in Heat Shock"
) + geom_jitter(size=3, width = .2
) + geom_point(data = ehs_wowdis_0, aes(x=hs.num, y=mass.end.g),
               size = 4, shape = 8, stroke = 1.75, color = "#1F9F88" 
) + geom_point(data = ehs_wowdis_1, aes(x=hs.num, y=mass.end.g),
               size = 4, shape = 8, stroke = 1.75, color = "#453275" 
) + labs(x="Days in Heat Shock", y="Mass at Cull [g]"
) + facet_wrap(~hs.temp
) + theme(strip.text = element_text(size = 18),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.position = "none")

wowe_mssend_viol_plot



#for legend
wowe_mssend_viol_plot_leg <- ggplot(ehs_wowe, aes(x=hs.num, y=mass.end.g, color=hs.num))
wowe_mssend_viol_plot_leg <- wowe_mssend_viol_plot_leg + geom_violin(aes(fill=hs.num),
                                                             alpha=.4
) + scale_fill_manual(values=c("#9105A5", "#CD4479", "#EC9221", "#E8D900"),
                      breaks=c(1,2,3,4),
                      name="Days in Heat Shock"
) + scale_color_manual(values=c("#9105A5", "#CD4479", "#EC9221", "#E8D900"),
                       breaks=c(1,2,3,4),
                       name="Days in Heat Shock"
) + geom_jitter(size=3, width = .2
) + geom_point(data = ehs_wowdis_0, aes(x=hs.num, y=mass.end.g),
               size = 4, shape = 8, stroke = 1.75, color = "#1F9F88" 
) + geom_point(data = ehs_wowdis_1, aes(x=hs.num, y=mass.end.g),
               size = 4, shape = 8, stroke = 1.75, color = "#453275" 
) + labs(x="Days in Heat Shock", y="Mass at Cull [g]"
) + facet_wrap(~hs.temp
) + theme(strip.text = element_text(size = 18),
          axis.line.x=element_line(colour = 'black', size = 1),
          axis.line.y=element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = 'black', size = 1),
          axis.ticks.length = unit(2.5, "mm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15))

wowe_mssend_viol_plot_leg

viol_leg <- get_legend(wowe_mssend_viol_plot_leg)



#get legend for dissected individuals
diss_wowe_leg <- ggplot(ehs_wowdis, aes(x=hs.num, y=mass.end.g, color=factor(diss.wasp)))
diss_wowe_leg + geom_point(size = 4, shape = 8, stroke = 1.75
) + scale_color_manual(values = c("#1F9F88", "#453275"),
                       breaks = c(0, 1),
                       labels = c("No", "Yes"),
                       name = "Wasp Larvae Found"
) + theme(legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))



#combine all plots into 3 panels
suppl_wowe_fig <- plot_grid(prop_wwasp_plot, mn_mssend_disw_plot, wowe_mssend_viol_plot,
                            nrow = 3, labels = c("A", "B", "C"),
                            align = "vh")

suppl_wowe_fig



#combine legends
suppl_leg <- plot_grid(diss_wasp_leg, viol_leg, 
                       nrow = 2, align = "v")

suppl_leg


#add legends to figure
suppl_wowe_fig <- plot_grid(suppl_wowe_fig, suppl_leg, 
                            ncol = 2, rel_widths = c(5, 1))

suppl_wowe_fig

