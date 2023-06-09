
#### Small mammals reduce distance-dependence and increase seed predation risk in tropical rainforest fragments ####
#### Biotropica
#### Authors: Aparna KRISHNAN*, Anand M OSURI, Meghna KRISHNADAS

#### Description: This file contains the R codes for Figure 2 and Figure S1. The codes generates a graphical visualisation of the proportion of seeds predated over time (in weeks) and the proportion of seeds consumed by mammal and insect predators in different treatment categories

library(dplyr)
library(ggplot2)
library(reshape2)
library(patchwork)

#### FIGURE 2 ####

seed_plots <- read.csv("seed_plots_raw.csv", header = T, stringsAsFactors = F)

# Removing all rows after the week 7 for the graphs
seed_plots <- seed_plots[-which(seed_plots$Week > 6),]

# proportion seeds predated
seed_plots$pct_pred <- seed_plots$Predated_Seeds*10

seed_pred <- seed_plots %>% group_by(Species_Name, 
                                    Landscape_Type,
                                    Plot_Location, 
                                    Treatment,
                                    Week) %>% summarise(
                                      mean_pct_pred = mean(pct_pred),
                                      se_pct_pred = sd(pct_pred)/sqrt(length(pct_pred)),
                                      ss = length(pct_pred))

seed_pred$Landscape_Type <- ifelse(seed_pred$Landscape_Type == "Contiguous", "Contiguous forest", "Forest fragments")

seed_nexl <- seed_pred[which(seed_pred$Treatment == "No_Exclosure"),]
seed_exl <- seed_pred[which(seed_pred$Treatment == "Exclosure"),]

#### Acronychia pedunculata ###

# Open Plots (No Exclosure)
seed_nexl_ap <- seed_nexl[which(seed_nexl$Species_Name == "AP"),]


AP_NE <- ggplot(data = seed_nexl_ap, aes(x = Week, y = mean_pct_pred))+
  geom_point(position = position_dodge(width = 0.5), size = 1.55, stroke = 0.5, aes(colour = Plot_Location, shape = Plot_Location))+
  geom_line(aes(colour = Plot_Location),position = position_dodge(width = 0.5))+
  geom_errorbar(data = seed_nexl_ap, 
                position = position_dodge(width = 0.5), width = 0.2,
                aes(ymin = mean_pct_pred - se_pct_pred, ymax =  mean_pct_pred + se_pct_pred, 
                    colour = Plot_Location))+
  xlab("Time (weeks)")+ ylab("% Predated")+
  scale_x_continuous(limits = c(-0.2,6.5),
                     breaks = round(0:6))+
  scale_y_continuous(limits = c(0,100))+
  ggtitle("Open plots","Acronychia pedunculata")+
  scale_color_manual("Plot Location", labels = c( "Far","Near"), values=c("#E7B800", "#00AFBB"))+
  scale_shape_manual("Plot Location", labels = c( "Far","Near"), values=c(19,17))+
  theme(axis.line = element_line(colour = "black"), 
        legend.position = "none",
        panel.background = element_blank(),
        plot.subtitle = element_text(face = "italic"))+
  scale_alpha(guide = F)+
  scale_linetype_manual(guide = guide_legend(reverse = T))+
  facet_wrap(.~factor(Landscape_Type, levels = c( "Forest fragments", "Contiguous forest")))


# Exclosure plots
seed_exl_ap <- seed_exl[which(seed_exl$Species_Name == "AP"),]
AP_E <- ggplot(data = seed_exl_ap, aes(x = Week, y = mean_pct_pred))+
  geom_point(position = position_dodge(width = 0.5), size = 1.55, stroke = 0.5, aes(colour = Plot_Location, shape = Plot_Location))+
  geom_line(aes(colour = Plot_Location),position = position_dodge(width = 0.5))+
  geom_errorbar(data = seed_exl_ap, 
                position = position_dodge(width = 0.5), width = 0.2,
                aes(ymin = mean_pct_pred - se_pct_pred, ymax =  mean_pct_pred + se_pct_pred, 
                    colour = Plot_Location))+
  xlab("Time (weeks)")+ ylab("% Predated")+
  scale_x_continuous(limits = c(-0.2,6.5),
                     breaks = round(0:6))+
  scale_y_continuous(limits = c(0,100))+
  ggtitle("Mammal exclosure plots","Acronychia pedunculata")+
  scale_color_manual("Landscape type", labels = c( "Far","Near"), values=c("#E7B800", "#00AFBB"))+
  scale_shape_manual("Landscape type", labels = c( "Far","Near"), values=c(19,17))+
  theme(axis.line = element_line(colour = "black"), 
        legend.position = "none",
        panel.background = element_blank(),
        plot.subtitle = element_text(face = "italic"))+
  scale_alpha(guide = F)+
  scale_linetype_manual(guide = guide_legend(reverse = T))+
  facet_wrap(.~factor(Landscape_Type, levels = c( "Forest fragments", "Contiguous forest")))

AP_NE + AP_E

#### Cullenia exarillata ###

# Open Plots (No Exclosure)
seed_nexl_ce <- seed_nexl[which(seed_nexl$Species_Name == "CE"),]

CE_NE <- ggplot(data = seed_nexl_ce, aes(x = Week, y = mean_pct_pred))+
  geom_point(position = position_dodge(width = 0.5), size = 1.55, stroke = 0.5, aes(colour = Plot_Location, shape = Plot_Location))+
  geom_line(aes(colour = Plot_Location),position = position_dodge(width = 0.5))+
  geom_errorbar(data = seed_nexl_ce, 
                position = position_dodge(width = 0.5), width = 0.2,
                aes(ymin = mean_pct_pred - se_pct_pred, ymax =  mean_pct_pred + se_pct_pred, 
                    colour = Plot_Location))+
  xlab("Time (weeks)")+ ylab("% Predated")+
  scale_x_continuous(limits = c(-0.2,6.5),
                     breaks = round(0:6))+
  scale_y_continuous(limits = c(0,100))+
  ggtitle("","Cullenia exarillata")+
  scale_color_manual("Landscape type", labels = c( "Far","Near"), values=c("#E7B800", "#00AFBB"))+
  scale_shape_manual("Landscape type", labels = c( "Far","Near"), values=c(19,17))+
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_blank(),
        legend.position = "none",
        plot.subtitle = element_text(face = "italic"))+
  scale_alpha(guide = F)+
  scale_linetype_manual(guide = guide_legend(reverse = T))+
  facet_wrap(.~factor(Landscape_Type, levels = c( "Forest fragments", "Contiguous forest")))


# Exclosure plots
seed_exl_ce <- seed_exl[which(seed_exl$Species_Name == "CE"),]
seed_exl_ce <- seed_exl_ce[-which(seed_exl_ce$Week > 5),]

CE_E <- ggplot(data = seed_exl_ce, aes(x = Week, y = mean_pct_pred))+
  geom_point(position = position_dodge(width = 0.5), size = 1.55, stroke = 0.5, aes(colour = Plot_Location, shape = Plot_Location))+
  geom_line(aes(colour = Plot_Location),position = position_dodge(width = 0.5))+
  geom_errorbar(data = seed_exl_ce, 
                position = position_dodge(width = 0.5), width = 0.2,
                aes(ymin = mean_pct_pred - se_pct_pred, ymax =  mean_pct_pred + se_pct_pred, 
                    colour = Plot_Location))+
  xlab("Time (weeks)")+ ylab("% Predated")+
  scale_x_continuous(limits = c(-0.2,6.5),
                     breaks = round(0:5))+
  scale_y_continuous(limits = c(0,100))+
  ggtitle("","Cullenia exarillata")+
  scale_color_manual("Landscape type", labels = c( "Far","Near"), values=c("#E7B800", "#00AFBB"))+
  scale_shape_manual("Landscape type", labels = c( "Far","Near"), values=c(19,17))+
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_blank(),
        legend.position = "none",
        plot.subtitle = element_text(face = "italic"))+
  scale_alpha(guide = F)+
  scale_linetype_manual(guide = guide_legend(reverse = T))+
  facet_wrap(.~factor(Landscape_Type, levels = c( "Forest fragments", "Contiguous forest")))

CE_NE + CE_E

#### Syzygium rubicundum ###

# Open Plots (No Exclosure)
seed_nexl_sr <- seed_nexl[which(seed_nexl$Species_Name == "SR"),]

SR_NE <- ggplot(data = seed_nexl_sr, aes(x = Week, y = mean_pct_pred))+
  geom_point(position = position_dodge(width = 0.5), size = 1.55, stroke = 0.5, aes(colour = Plot_Location, shape = Plot_Location))+
  geom_line(aes(colour = Plot_Location),position = position_dodge(width = 0.5))+
  geom_errorbar(data = seed_nexl_sr, 
                position = position_dodge(width = 0.5), width = 0.2,
                aes(ymin = mean_pct_pred - se_pct_pred, ymax =  mean_pct_pred + se_pct_pred, 
                    colour = Plot_Location))+
  xlab("Time (weeks)")+ ylab("% Predated")+
  scale_x_continuous(limits = c(-0.2,6.5),
                     breaks = round(0:6))+
  scale_y_continuous(limits = c(0,100))+
  ggtitle("","Syzigium rubicundum")+
  scale_color_manual("Landscape type", labels = c( "Far","Near"), values=c("#E7B800", "#00AFBB"))+
  scale_shape_manual("Landscape type", labels = c( "Far","Near"), values=c(19,17))+
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_blank(),
        legend.position = "none",
        plot.subtitle = element_text(face = "italic"))+
  scale_alpha(guide = F)+
  scale_linetype_manual(guide = guide_legend(reverse = T))+
  facet_wrap(.~factor(Landscape_Type, levels = c( "Forest fragments", "Contiguous forest")))


# Exclosure plots
seed_exl_sr <- seed_exl[which(seed_exl$Species_Name == "SR"),]

SR_E <- ggplot(data = seed_exl_sr, aes(x = Week, y = mean_pct_pred))+
  geom_point(position = position_dodge(width = 0.5), size = 1.55, stroke = 0.5, aes(colour = Plot_Location, shape = Plot_Location))+
  geom_line(aes(colour = Plot_Location),position = position_dodge(width = 0.5))+
  geom_errorbar(data = seed_exl_sr, 
                position = position_dodge(width = 0.5), width = 0.2,
                aes(ymin = mean_pct_pred - se_pct_pred, ymax =  mean_pct_pred + se_pct_pred, 
                    colour = Plot_Location))+
  xlab("Time (weeks)")+ ylab("% Predated")+
  scale_x_continuous(limits = c(-0.2,6.5),
                     breaks = round(0:6))+
  scale_y_continuous(limits = c(0,100))+
  ggtitle("","Syzigium rubicundum")+
  scale_color_manual("Landscape type", labels = c( "Far","Near"), values=c("#E7B800", "#00AFBB"))+
  scale_shape_manual("Landscape type", labels = c( "Far","Near"), values=c(19,17))+
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_blank(),
        legend.position = "none",
        plot.subtitle = element_text(face = "italic"))+
  scale_alpha(guide = F)+
  scale_linetype_manual(guide = guide_legend(reverse = T))+
  facet_wrap(.~factor(Landscape_Type, levels = c( "Forest fragments", "Contiguous forest")))

SR_NE + SR_E

#### Ormosia travancorica ###

# Open Plots (No Exclosure)
seed_nexl_ot <- seed_nexl[which(seed_nexl$Species_Name == "OT"),]
unique(seed_nexl_ot$Landscape_Type)


OT_NE <- ggplot(data = seed_nexl_ot, aes(x = Week, y = mean_pct_pred))+
  geom_point( position = position_dodge(width = 0.5), size = 1.55, stroke = 0.5, aes(colour = Plot_Location, shape = Plot_Location))+
  geom_line(aes(colour = Plot_Location),position = position_dodge(width = 0.5))+
  geom_errorbar(data = seed_nexl_ot, 
                position = position_dodge(width = 0.5), width = 0.2,
                aes(ymin = mean_pct_pred - se_pct_pred, ymax =  mean_pct_pred + se_pct_pred, 
                    colour = Plot_Location))+
  xlab("Time (weeks)")+ ylab("% Predated")+
  scale_x_continuous(limits = c(-0.2,6.5),
                     breaks = round(0:6))+
  scale_y_continuous(limits = c(0,100))+
  ggtitle("","Ormosia travancorica")+
  scale_color_manual("Plot Location", labels = c( "Far","Near"), values=c("#E7B800", "#00AFBB"))+
  scale_shape_manual("Plot Location", labels = c( "Far","Near"), values=c(19,17))+
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_blank(),
        plot.subtitle = element_text(face = "italic"),
        legend.position = "none",
        legend.key=element_blank())+
  scale_alpha(guide = F)+
  scale_linetype_manual(guide = guide_legend(reverse = T))+
  facet_wrap(.~factor(Landscape_Type, levels = c( "Forest fragments", "Contiguous forest")))

Figure_2 <- (AP_NE + AP_E) / (CE_NE + CE_E) /(SR_NE + SR_E) / (OT_NE + ggplot()) 
ggsave(filename = "Figure_2", path = "~/Desktop", width = 8.5, height = 11, device='tiff', dpi=700)


## FIGURE S3 ####

library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

seed_plots <- read.csv("seed_plots_raw.csv", header = T, stringsAsFactors = F)

# Extracting the 6th Week
seed_plots_6 <- seed_plots %>% filter(Week == 6)

seed_pred <- pivot_longer(seed_plots_6, cols = c(11,12),  names_to = "Predator_Species", values_to = "No_Predated")
seed_pred$tot_seeds <- 10
who_predated <- seed_pred %>% group_by(Species_Name, Landscape_Type, Treatment, Predator_Species) %>% summarise(Pct_Pred = 100*sum(No_Predated)/sum(tot_seeds))
lab_treatment <- c('Exclosure' = "Mammal exclosure plots",
                   'No_Exclosure' = "Open plots")

# Acronychia pedunculata
ap <- who_predated %>% filter(Species_Name == "AP") 
AP <- ggplot(ap, aes(x = Landscape_Type, y = Pct_Pred))+
  geom_bar(stat = "identity", aes(fill = Predator_Species), width = 0.6)+
  xlab("Treatment")+ylab("% Predation")+
  ggtitle("", "Acronychia pedunculata")+
  scale_fill_manual("Species Group", labels = c("Insect Predated", "Mammal Predated"), values=c("#E7B800","#00AFBB"))+
  theme(axis.line = element_line(colour = "black"), 
        plot.subtitle = element_text(face = "italic"),
        panel.background = element_blank())+
  scale_y_continuous(expand = c(0.01,0), limits = c(0,100))+
  facet_wrap(.~Treatment, labeller = as_labeller(lab_treatment))
  
# Cullenia exarillata
ce <- who_predated %>% filter(Species_Name == "CE") 
CE <- ggplot(ce, aes(x = Landscape_Type, y = Pct_Pred))+
  geom_bar(stat = "identity", aes(fill = Predator_Species), width = 0.6)+
  xlab("Treatment")+ylab("% Predation")+
  ggtitle("", "Cullenia exarillata")+
  scale_fill_manual("Species Group", labels = c("Insect Predated", "Mammal Predated"), values=c("#E7B800","#00AFBB"))+
  theme(axis.line = element_line(colour = "black"), 
        plot.subtitle = element_text(face = "italic"),
        panel.background = element_blank())+
  scale_y_continuous(expand = c(0.01,0), limits = c(0,100))+
  facet_wrap(.~Treatment, labeller = as_labeller(lab_treatment))

# Ormosia travancorica
ot <- who_predated %>% filter(Species_Name == "OT")

#creating an empty row for the exclosure plot
temp <- ot[1,] 
temp["Treatment"] <- "Exclosure"
temp["Pct_Pred"] <- 0
ot[5,] <- temp

OT <- ggplot(ot, aes(x = Landscape_Type, y = Pct_Pred))+
  geom_bar(stat = "identity", aes(fill = Predator_Species), width = 0.6)+
  xlab("Treatment")+ylab("% Predation")+
  ggtitle("", "Ormosia travancorica")+
  scale_fill_manual("Species Group", labels = c("Insect Predated", "Mammal Predated"), values=c("#E7B800","#00AFBB"))+
  theme(axis.line = element_line(colour = "black"), 
        plot.subtitle = element_text(face = "italic"),
        panel.background = element_blank())+
  scale_y_continuous(expand = c(0.01,0), limits = c(0,100))+
  facet_wrap(.~Treatment, labeller = as_labeller(lab_treatment))

# Syzygium rubicundum
sr <- who_predated %>% filter(Species_Name == "SR") 
SR <- ggplot(sr, aes(x = Landscape_Type, y = Pct_Pred))+
  geom_bar(stat = "identity", aes(fill = Predator_Species), width = 0.6)+
  xlab("Treatment")+ylab("% Predation")+
  ggtitle("", "Syzygium rubicundum")+
  scale_fill_manual("Species Group", labels = c("Insect Predated", "Mammal Predated"), values=c("#E7B800","#00AFBB"))+
  theme(axis.line = element_line(colour = "black"), 
        plot.subtitle = element_text(face = "italic"),
        panel.background = element_blank())+
  scale_y_continuous(expand = c(0.01,0), limits = c(0,100))+
  facet_wrap(.~Treatment, labeller = as_labeller(lab_treatment))

Figure_S1 <- AP / CE / SR / OT
ggsave(filename = "Figure_S1", path = "~/Desktop", width = 6, height = 11, device='pdf', dpi=600)

