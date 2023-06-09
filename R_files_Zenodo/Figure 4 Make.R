## R code to produce Figure 4 in:
## Johnston EC, Cunning R, Burgess SC. Cophylogeny and specificity between cryptic coral species (Pocillopora spp.) at Moâ€™orea and their symbionts (Symbiodiniaceae).
## Code written by Erika Johnston. February 2022. Send comments or corrections to ejohnston@bio.fsu.edu
## R version 4.1.2 (2022-03-01)

# initialization ----------------------------------------------------------
rm(list=ls()) #clear all variables

#setwd("")

library(ggplot2)
library(ggstance)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(tidyr)
library(ggh4x)
library(forcats)


# Import data
dat <- read.csv("Figure 4 and 6 data.csv")

# Import colors for ITS2 sequence data
my_colors <- read.csv("Figure 4 colors.csv")

#### Plot ITS2 type profiles ####

#Order ITS2 type profiles
dat$Type_profile = factor(dat$Type_profile,levels=c(
                       "C1.C42.2",
                       "C1.C3.C42.2.C1b",
                       "C1.C3.C42.2.C1b.C3cg",
                       "C42.2.C1.C1b.C3.C1au",
                       "C42.2.C1.C42i.C3.C1b.C1au",
                       "C1.C42.2.C3.C1b.C1au.C115n",
                       "C1.C42.2.C3cg.C1b.C115k.C45c",
                       "C1.C42.2.C1b.C1au.C3.C3cg.C1j",
                       "C1.C1ag.C1ah.C42.2.C3cg.C1b.C3",
                       "C1.C42.2.C1b.C3cg.C3.C1au.C41p",
                       "C1.C1m.C1ag.C42.2.C3cg.C1b.C3cw.C1au",
                       "C1d.C42.2.C1.C1k.C1b",
                       "C1d.C42.2.C1.C3cg.C1b.C1au",
                       "C1d.C1.C42.2.C1b.C3cg.C3",
                       "C1d.C1.C42.2.C1k.C1b.C3cg.C3cw",
                       "C1.C42.2.C1d.C1b.C3.C3cg.C115n.C1au",
                       "C1d.C1.C42.2.C3cg.C1b.C45c.C115k.C3.C1au",
                       "C1.C1ah.C1ag.C42.2.C3cg.C1b.C1dj.C115k.C45c",
                       "C1d.C1.C42.2.C3.C1b.C3cg.C115k.C45c.C1au.C41p",
                       "C42a.C42.2.C1.C1b",
                       "C42.2.C42a.C1b.C1.C1j",
                       "C42a.C42.2.C1b.C1.C42b",
                       "C42a.C42.2.C1.C3.C1b.C1au",
                       "C42a.C42.2.C1.C1b.C3.C1j.C1au",
                       "C42a.C42.2.C1b.C1.C42ao.C1au",
                       "C42.2.C42a.C1.C1j.C1b.C1ci.C1au",
                       "C42g.C1.C42.2.C42a.C1b.C42h.C3.C1au",
                       "C1.C42.2.C42a.C1b.C1az.C1au.C115d.C3.C1j",
                       "C42a.C1.C42.2.C1j.C1b.C1au.C3.C115l.C115k",
                       "C1.C42.2.C42g.C42a.C1b.C3.C1az.C1au.C42az",
                       "C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C115k.C115d.C41p",
                       "C3.C116.C3dc",
                       "C116.C1.C116q",
                       "C15.C15ev.C15dt",
                       "C116.C116i.C116q",
                       "C42.2.C116.C1.C1b",
                       "A1.A1fi",
                       "D1",
                       "D1ch",
                       "D1.D6",
                       "D1.D6.D4.D2.2.D1h.D2c.D2f"))

#Set ITS2 profile colors
my_colors_profile <- c("C1.C42.2"="#ffedfe",
                       "C1.C3.C42.2.C1b"="#fbccff",
                       "C1.C3.C42.2.C1b.C3cg"="#f2b3ff",
                       "C42.2.C1.C1b.C3.C1au"="#fa7de1",
                       "C42.2.C1.C42i.C3.C1b.C1au"="#d152eb",
                       "C1.C42.2.C3.C1b.C1au.C115n"="#bf4e92",
                       "C1.C42.2.C3cg.C1b.C115k.C45c"="#e35675",
                       "C1.C42.2.C1b.C1au.C3.C3cg.C1j"="#bd2848",
                       "C1.C1ag.C1ah.C42.2.C3cg.C1b.C3"="#8a017a",
                       "C1.C42.2.C1b.C3cg.C3.C1au.C41p"="#4a0475",
                       "C1.C1m.C1ag.C42.2.C3cg.C1b.C3cw.C1au"="#230436",
                       "C1d.C42.2.C1.C1k.C1b"="#fff5b5",
                       "C1d.C42.2.C1.C3cg.C1b.C1au"="#fcf5ca",
                       "C1d.C1.C42.2.C1b.C3cg.C3"="#f5e171",
                       "C1d.C1.C42.2.C1k.C1b.C3cg.C3cw"="#fcbf74",
                       "C1.C42.2.C1d.C1b.C3.C3cg.C115n.C1au"="#ed9f18",
                       "C1d.C1.C42.2.C3cg.C1b.C45c.C115k.C3.C1au"="#de7309",
                       "C1.C1ah.C1ag.C42.2.C3cg.C1b.C1dj.C115k.C45c"="#e85702",
                       "C1d.C1.C42.2.C3.C1b.C3cg.C115k.C45c.C1au.C41p"="#b85707",
                       "C42a.C42.2.C1.C1b"="#c4fffc",
                       "C42.2.C42a.C1b.C1.C1j"="#9cfffa",
                       "C42a.C42.2.C1b.C1.C42b"="#5cf1ff",
                       "C42a.C42.2.C1.C3.C1b.C1au"="#74f1fc",
                       "C42a.C42.2.C1.C1b.C3.C1j.C1au"="#21dced",
                       "C42a.C42.2.C1b.C1.C42ao.C1au"="#259dd9",
                       "C42.2.C42a.C1.C1j.C1b.C1ci.C1au"="#52a0fa",
                       "C42g.C1.C42.2.C42a.C1b.C42h.C3.C1au"="#286bb8",
                       "C1.C42.2.C42a.C1b.C1az.C1au.C115d.C3.C1j"="#0353ad",
                       "C42a.C1.C42.2.C1j.C1b.C1au.C3.C115l.C115k"="#143ec9",
                       "C1.C42.2.C42g.C42a.C1b.C3.C1az.C1au.C42az"="#01068f",
                       "C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C115k.C115d.C41p"="#05165c",
                       "C3.C116.C3dc"="#edebeb",
                       "C116.C1.C116q"="#b8b6b6",
                       "C15.C15ev.C15dt"="#adadad",
                       "C116.C116i.C116q"="#858585",
                       "C42.2.C116.C1.C1b"="#5e5d5d",
                       "A1.A1fi"="#ffee00",
                       "D1"="#eafc5d",
                       "D1ch"="#c0d428",
                       "D1.D6"="#a1c91e",
                       "D1.D6.D4.D2.2.D1h.D2c.D2f"="#018a26")

#Order Pocillopora species and haplotypes
dat$Species.haplotype = factor(dat$Species.haplotype,levels=c(
  "Haplotype 10",
  "P. verrucosa",
  "P. meandrina",
  "Haplotype 8a",
  "P. grandis",
  "Haplotype 11",
  "Haplotype 2"))

#Order depths
dat$Depth.m = factor(dat$Depth.m,levels=c("5","10","20"))

# Shorten Haplotype 2 to Hap 2 to plot
Hap.labs <- c("Haplotype 10","P. meandrina","Haplotype 11","Hap 2","P. grandis","P. verrucosa","Haplotype 8a")
names(Hap.labs) <- c("Haplotype 10","P. meandrina","Haplotype 11","Haplotype 2","P. grandis","P. verrucosa","Haplotype 8a")

# Add "m" to depths to plot
Depth.labs <- c("5m","10m","20m")
names(Depth.labs) <- c("5", "10", "20")


####  Plot ITS2 type profiles
ggplot(dat, aes(x=Coral.ID, y=Type_profile_Prop))+
  geom_col(aes(fill = Type_profile), colour="grey", size=0.005)+
  labs(x="Colony", y="Proportion")+
  theme_bw(base_size = 10)+
  scale_y_continuous(labels=function(Prop)Prop+1)+
  scale_fill_manual(values=my_colors_profile, name="ITS2 Type Profile")+
  guides(color = guide_legend(override.aes = list(size=1.5)), fill=guide_legend(ncol=1, title.theme = element_text(angle = 90)))+
  theme(panel.spacing=unit(0.2,"lines"),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), text=element_text(size=7),panel.grid = element_blank(), strip.text = element_text(size = 7), legend.position = "bottom", legend.key.size = unit(0.3, "cm"))+
  facet_nested(.~Species.haplotype+Depth.m, scales="free", space="free",
               labeller = labeller(Species.haplotype = Hap.labs, Depth.m = Depth.labs))


#### Plot ITS2 sequences ####

# Convert to long format to plot by ITS2 sequence
dat_long <- gather(dat,Symbio.clade,Prop,7:628)

# Remove redundant ITS2 type profile data - Plot based on majority ITS2 type profile
dat_long <- dat_long %>%
  filter(Type_profile_Prop > 0.43)

#Order ITS2 type sequences by the order in which they appear in dat_long
dat_long$Symbio.clade <- fct_inorder(dat_long$Symbio.clade)

####  Plot ITS2 sequences 
ggplot(dat_long, aes(x=Coral.ID,y=desc(Prop)))+
  geom_col(aes(fill = Symbio.clade), colour="grey", size=0.005)+
  labs(x="Colony", y="Proportion")+
  theme_bw(base_size = 10)+
  scale_y_continuous(labels=function(Prop)Prop+1)+
  scale_fill_manual(values=my_colors$my_colors, name="ITS2 Sequences")+
  guides(color = guide_legend(override.aes = list(size=0.25)), fill=guide_legend(ncol=15, title.theme = element_text(angle = 90)))+
  theme(panel.spacing=unit(0.2,"lines"),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), text=element_text(size=7),panel.grid = element_blank(), strip.text = element_text(size = 7), legend.position = "bottom", legend.key.size = unit(0.3, "cm"))+
  facet_nested(.~Species.haplotype+Depth.m, scales="free", space="free",
             labeller = labeller(Species.haplotype = Hap.labs,Depth.m = Depth.labs))

