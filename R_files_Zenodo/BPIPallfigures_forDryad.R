####Code for all figures in BPIP paper####
#Note that these figures were edited in Adobe Illustrator after generation in R

#All packages used
library(ggplot2)
library(scales)
library(cowplot)
library(ggpubr)
library(Cairo)

#Adding colorblind and printer friendly palette
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

########################FIG 1##############################
cfudata <- read.table("BPIPCFUforR_forDryad.txt", header = TRUE, fill = TRUE)

#Setting Week variable as factor
cfudata$Week <- as.factor(cfudata$Week)

#Separating MgCl2 control from sample data
cfudata <- subset(cfudata, Treatment != "MGCL2")

#Removing Phage only data
cfudatanoP <- subset(cfudata, Treatment != "P")

#Making new column with log10CFU
cfudatanoP["log10CFU"] <- NA
cfudatanoP$log10CFU <- log10(cfudatanoP$CFUperml)

#Facet labels
Envirolabels <- c(
  `Invitro` = "In vitro lines",
  `Plant` = "In planta lines")

#Creating dataframe with threshold values for plant and in vitro data
thresholdsfig1 <- data.frame(Environment = c("Invitro", "Plant"), yval = c(log10(166.5), log10(.3)))

#Plot code for figure 1
mixplot <- ggplot(cfudatanoP, aes(x=Week, col= Treatment, linetype=Treatment)) + 
  facet_wrap(~ Environment, labeller = as_labeller(Envirolabels)) +
  geom_line(aes(x= Week, y=log10CFU, col = Treatment, group = Sample, linetype=Treatment)) + 
  labs(y="Bacterial density (log10 CFU per ml culture or 1mm2 leaf tissue)", 
       x= "Passage number",
       color=NULL) +
  scale_color_manual(values=cbPalette, 
                     name="Treatment",
                     breaks=c("B", "BP", "Bpanc", "BancP"),
                     labels=c("Bacteria only", "Coevolving", "Evolving B, ancestral P", "Ancestral B, evolving P")) +
  scale_linetype_manual(name="Treatment", values=c("solid", "dashed", "dashed", "solid"), breaks=c("B", "BP", "Bpanc", "BancP"),
                        labels=c("Bacteria only", "Coevolving", "Evolving B, ancestral P", "Ancestral B, evolving P"))+
  scale_y_continuous(limits = c(0, 11.5)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(data=thresholdsfig1, aes(yintercept=yval))

#Saving plot as .eps file
ggsave("Figure_1.eps", plot=mixplot, width= 300, height = 200, units=c("mm"), dpi=600)

########################FIG 2############################
res <- read.table("BPIPresfreqonly_forDryad.txt", header = TRUE, fill = TRUE)

#Setting variables as factors, and levels of treatment
res$Bweek <- as.factor(res$Bweek)
res$Line <- as.factor(res$Line)
res$Treatment <- factor(res$Treatment, levels = c("B", "BP", "Bpanc", "BancP"))

#Adding labels for facets
weeklabels <- c('3' = "Week 3 bacteria", '6' = "Week 6 bacteria")

#Subsetting plant and invitro data
plantonly <- subset(res, Environment != "Invitro")
invitroonly <- subset(res, Environment != "Plant")

#Plant plot
plantplot <- ggplot(plantonly, aes(x=Treatment, y=Prop_notS, fill = Phage)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.45,
               position=position_dodge(0.8)) +
  scale_fill_manual(values = c("#000000", "#B3B3B3"), name= "Phage type", breaks = c("Ancestral", "Contemporary"), labels=c("Ancestral","Contemporary")) +
  facet_wrap(~ Bweek, labeller=labeller(Bweek = weeklabels)) +
  ylab("Proportion not sensitive") +
  labs(subtitle = "In planta lines") +
  scale_x_discrete(breaks=c("B","BP","BancP", "Bpanc"),
                   labels=c("Bacteria only", "Coevolving", "Ancestral B, evolving P", "Evolving B, ancestral P")) +
  ylim(0,1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1))
#Note: plant plot values are incorrectly plotted in ggplot due to issues with resolution
#Values were adjusted in illustrator to all line up at 0
#Specific values are listed in figure legend, and available in Dryad

#In vitro plot
invitroplot <-  ggplot(invitroonly, aes(Treatment, Prop_notS, fill = Phage)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.45,
               position=position_dodge(0.8)) +
  labs(subtitle = "In vitro lines") +
  scale_fill_manual(values = c("#000000", "#B3B3B3"), name= "Phage type") +
  facet_wrap(~ Bweek, labeller=labeller(Bweek = weeklabels)) +
  ylab("Proportion not sensitive") +
  scale_x_discrete(breaks=c("B","BP","BancP", "Bpanc"),
                   labels=c("Bacteria only", "Coevolving", "Ancestral B, evolving P", "Evolving B, ancestral P")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1))

#Combined plot for publication
combinedplot <- plot_grid(invitroplot, plantplot, labels = c("A", "B"))

#Saving plot as .eps file
ggsave("Figure_2.eps", plot=combinedplot, width= 400, height = 200, units=c("mm"), dpi=600)

########################FIG 3############################
res3avg <- read.table("Avgdleaf_res3_feb2018_forDryad.txt", header = TRUE, fill = TRUE)

#Setting Time variable as factor
res3avg$Time <- as.factor(res3avg$Time)

#Removing Time 0 due to too many points below limit of detection
res3avgno0 <- subset(res3avg, Time!= "0")

#Removing ancestral PT23 from plot (no replicates)
res3avgnoPT23or0 <- subset(res3avgno0, Treatment != "PT23")

#Facet labels
facetlabels <- c('B' = "Bacteria only", 'BP' = "Coevolving", 'Bpanc' = "Evolving B, ancestral P")

#Final plot
avgcfuplot3 <- ggplot(subset(res3avgnoPT23or0), aes(Time, log10(AvgCFUperml), fill = Phage)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6,
               position=position_dodge(0.8)) +
  facet_wrap (~ Treatment, labeller=labeller(Treatment = facetlabels)) +
  labs(y= "Bacterial density (log10 CFU per 1mm2 leaf tissue)", x= "Time (Hours post-inoculation)") +
  scale_fill_manual(values=c("#000000", "#B3B3B3"), 
                    name="Phage",
                    breaks=c("No", "Yes")) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 7.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Saving plot as .eps file
ggsave("Figure_3.eps", plot=avgcfuplot3, width= 225, height = 150, units=c("mm"), dpi=600)

########################FIG 4############################
d <- read.table("resinvivo2_forR_forDryad.txt", header = TRUE, fill = TRUE)

#Setting Time as factor
d$Time <- as.factor(d$Time)

#Subsetting samples
dsamples <- subset(d, Treatment != "Phage_only")
dsamplesno0 <- subset(dsamples, Time != "0")
dsamplesnoPT <- subset(dsamples, Treatment != "Ancestral")

#Want samples with phage only to plot phage over time
dphage <- subset(d, Treatment != "Ancestral")
dphageyes <- subset(dphage, Phage != "No")

#Treatment facet labeller
facetlabels <- c('Bonly' = "Bacteria only", 'Bpanc' = "Evolving bacteria with ancestral phage")

#Plot over time bacterial density
Bovertime <- ggplot(dsamplesnoPT, aes(x=Time, y=log10(AvgB), fill=Phage)) +
  scale_fill_manual(values=c("#000000", "#B3B3B3"), labels=c("No","Yes")) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.45,
               position=position_dodge(0.8)) +
  facet_wrap(~Treatment, labeller=labeller(Treatment = facetlabels)) +
  labs(y= "Bacterial density (log10 copies per 1mm2 leaf tissue)", x= "Time (Hours post-inoculation)", subtitle = "Bacterial density over time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0, 7.5)

#Subsetting samples with phage applied
dphagesamples <- subset(dsamplesnoPT, Phage != "No")

#Plot of phage over time
Povertime <- ggplot(dphageyes, aes(x=Time, y=log10(AvgP), fill=Treatment)) +
  scale_fill_manual(values=c("#000000", "#B3B3B3", "#FFFFFF"), labels=c("B only + phage","B evolved + phage", "P only")) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.45,
               position=position_dodge(0.8)) +
  labs(y= "Phage density (log10 copies per 1mm2 leaf tissue)", x= "Time (Hours post-inoculation)", subtitle = "Phage density over time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0, 5.25)

#Combining all into one plot
combinedplot <- plot_grid(Bovertime, Povertime, labels = c("A", "B"))

#Saving plot as .eps file
ggsave("Figure_4.eps", plot=combinedplot, width= 400, height = 200, units=c("mm"), dpi=600)

########################FIG 5############################
pdata <- read.table("Fall 2017 Phage cocktail plant data for R_forDryad.txt", header = TRUE, fill = TRUE)

#Setting Week, Phage, Treatment as factors
pdata$Week <- as.factor(pdata$Week)
pdata$Phage <- factor(pdata$Phage, levels=c("None", "FRS", "SHL"))
pdata$Treatment <- factor(pdata$Treatment, levels = c("B", "BFRS", "BSHL", "BFRSSHL"))

#Setting up facet labels
Envirolabels <- c(
  `Invitro` = "In vitro lines",
  `Plant` = "In planta lines")

##Want to subset data that does not include phage only or MGCL2 controls
pcnoFRSonly <- subset(pdata, Treatment != "FRS")
pcnophageonly <- subset(pcnoFRSonly, Treatment !="SHL")
pcsamples <- subset(pcnophageonly, Treatment != "MGCL2")

#Removing PT23
pcsamplesnoPT <- subset(pcsamples, Treatment != "PT23")

#Subsetting pcsamplesnoPT, removing phage plates
pcsamplesnoPTorFRS <- subset(pcsamplesnoPT, Phage != "FRS")
pcsamplesonly <- subset(pcsamplesnoPTorFRS, Phage != "SHL")

#Calculating log10CFU data
pcsamplesonly$log10CFU <- log10(pcsamplesonly$CFUperml)

#Creating dataframe with threshold values for plant and in vitro data
thresholdsfig5 <- data.frame(Environment = c("Invitro", "Plant"), yval = c(log10(166.5), log10(.295)))

#Making graphs of each line over time separately
mixplot <- ggplot(pcsamplesonly, aes(x=Week, col= Treatment, linetype=Treatment)) + 
  facet_wrap(~ Environment, labeller = as_labeller(Envirolabels)) +
  geom_line(aes(x= Week, y=log10CFU, col = Treatment, group = Sample, linetype=Treatment)) + 
  labs(y="Bacterial density (log10 CFU per ml culture or 1mm2 leaf tissue)", 
       x= "Passage number",
       color=NULL) +
  scale_color_manual(values=c("#000000", "#56B4E9", "#009E73", "#E69F00"), 
                     name="Treatment",
                     breaks=c("B", "BFRS", "BSHL", "BFRSSHL"),
                     labels=c("Bacteria only", "Bacteria + FRS", "Bacteria + SHL", "Bacteria + FRS + SHL")) +
  scale_linetype_manual(name="Treatment", values=c("solid", "dashed", "solid", "dashed"), breaks=c("B", "BFRS", "BSHL", "BFRSSHL"),
                        labels=c("Bacteria only", "Bacteria + FRS", "Bacteria + SHL", "Bacteria + FRS + SHL"))+
  scale_y_continuous(limits = c(0, 9.5)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(subtitle = "Bacterial densities after each passage") +
  geom_hline(data=thresholdsfig5, aes(yintercept=yval))

#Separating week 3 data
pcsamplesnoPTwithphage <- subset(pcsamplesnoPT, Phage != "None")
pcsamplesnoPTwithphagenowk1 <- subset(pcsamplesnoPTwithphage, Week != "1")
pcsamplesnoPTwithphage3 <- subset(pcsamplesnoPTwithphagenowk1, Week != "2")

#Setting levels of Treatment
pcsamplesnoPTwithphage3$Treatment <- factor(pcsamplesnoPTwithphage3$Treatment, levels = c("B", "BFRS", "BSHL", "BFRSSHL"))

#Propotion resistance plot
combinedresplot <- ggplot(pcsamplesnoPTwithphage3, aes(x=Treatment, y=Propres_minusPT23prop, fill = Phage)) +
  scale_fill_manual(values = c("#000000", "#B3B3B3"),
                    name="Phage type") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.45,
               position=position_dodge(0.8)) +
  facet_wrap(~ Environment, labeller=labeller(Environment = Envirolabels)) +
  ylab("Proportion resistant") +
  scale_x_discrete(breaks=c("B","BFRS","BSHL", "BFRSSHL"),
                   labels=c("Bacteria only", "Bacteria + FRS", "Bacteria + SHL", "Bacteria + FRS + SHL")) +
  ylim(0,1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(subtitle = "Proportion of resistant colonies after third passage")

#Combined densities over time and proportion resistance plot
combinedplot <- plot_grid(mixplot, combinedresplot, labels=c("A", "B"))

#Saving plot as .eps file
ggsave("Figure_5.eps", plot=combinedplot, width= 400, height = 200, units=c("mm"), dpi=600)

########################FIG 6############################
d <- read.table("RS_FRS_invivoinvitro_cleanedJan2019_forRforDryad.txt", header = TRUE, fill = TRUE)

#Setting Time as factor
d$Time <- as.factor(d$Time)

#Subsetting non-control samples
dsamples <- subset(d, Res_category != "NA")

#Subsetting plant and invitro data
dplant <- subset(dsamples, Environment != "Media")
dmedia <- subset(dsamples, Environment != "Plant")

#Adding labels
facetlabels <- c('R' = "Resistant", 'S' = "Sensitive")
phagelabels <- c('N' = "No", 'Y' = "Yes")

#Plotting plant data
plantdata <- ggplot(dplant, aes(x=Time, y=log10(AvgB_limitadjusted), fill=Phage)) +
  scale_fill_manual(values=c("#000000", "#B3B3B3"), labels=c("No","Yes")) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.45,
               position=position_dodge(0.8)) +
  facet_wrap(~ Res_category, labeller=labeller(Res_category = facetlabels)) +
  labs(y= "Bacterial density (log10 copies per 1mm2 leaf tissue)", x= "Time (Hours post-inoculation)", subtitle ="In planta lines") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0, 7.5) +
  scale_x_discrete(breaks=c("0","24","72"),
                   labels=c("1", "24", "72")) +
  labs(subtitle = "In planta growth") +
  geom_hline(aes(yintercept=log10(2.93)))

#Plotting in vitro data
mediadata <- ggplot(dmedia, aes(x=Time, y=log10(AvgB_limitadjusted), fill=Phage)) +
  scale_fill_manual(values=c("#000000", "#B3B3B3"), labels=c("No","Yes")) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.45,
               position=position_dodge(0.8)) +
  facet_wrap(~ Res_category, labeller=labeller(Res_category = facetlabels)) +
  labs(y= "Bacterial density (log10 copies per ml)", x= "Time (Hours post-inoculation)", subtitle ="In vitro lines") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0, 10) +
  scale_x_discrete(breaks=c("0","24","72"),
                   labels=c("1", "24", "72")) +
  labs(subtitle = "In vitro growth") +
  geom_hline(aes(yintercept=log10(824.175)))

#Combined plot both panels
combinedplot <- plot_grid(mediadata, plantdata, labels = c("A", "B"))

#Saving plot as .eps file
ggsave("Figure_6.eps", plot=combinedplot, width= 400, height = 200, units=c("mm"), dpi=600)

########################SUPP FIG 1############################
pfudata <- read.table("AmplifiedPFUsforR_forDryad.txt", header = TRUE, fill = TRUE)

#Setting Week as factor
pfudata$Week <- as.factor(pfudata$Week)

#Separating MgCl2 from sample data
pfudata <- subset(pfudata, Treatment != "MGCL2")

#Adding log10 column
pfudata["log10PFU"] <- NA
pfudata$log10PFU <- log10(pfudata$PFUperml_limitadjusted)

#Renaming treatments
labels <- c(B = "Bacteria only", BP = "Coevolving B and P",  BANCP = "Ancestral B, evolving P", BPANC = "Ancestral P, evolving B")

#Renaming environments for facet labels
Envirolabels <- c(
  `Invitro` = "In vitro lines",
  `Plant` = "In planta lines")

#Setting thresholds
thresholdssupfig1 <- data.frame(Environment = c("Invitro", "Plant"), yval = c(log10(13320), log10(24)))

#Making plot of densities over time for each line individually, color for treatments, facet by environment
mixplotPFU <- ggplot(pfudata, aes(x=Week, col= Treatment, linetype = Treatment)) + 
  facet_wrap(~ Environment, labeller = as_labeller(Envirolabels)) +
  geom_line(aes(x= Week, y=log10PFU, col = Treatment, group = Line, linetype = Treatment)) + 
  labs(y="Amplified phage density (log10 PFU per ml culture or 1mm2 leaf tissue)", 
       x= "Passage number",
       color=NULL) +
  scale_color_manual(values=cbPalette, 
                     name="Treatment",
                     breaks=c("B", "BP", "BPANC", "BANCP", "P"),
                     labels=c("Bacteria only", "Coevolving", "Evolving B, ancestral P", "Ancestral B, evolving P", "Phage only")) +
  scale_y_continuous(limits = c(0, 14)) +
  scale_linetype_manual(name="Treatment", values=c("solid", "dashed", "dashed", "solid", "solid"), breaks=c("B", "BP", "BPANC", "BANCP", "P"),
                        labels=c("Bacteria only", "Coevolving", "Evolving B, ancestral P", "Ancestral B, evolving P", "Phage only"))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(data=thresholdssupfig1, aes(yintercept=yval))

#Saving plot as .eps file
ggsave("Supp_Figure_1.eps", plot=mixplotPFU, width= 300, height = 200, units=c("mm"), dpi=600)

########################SUPP FIG 2############################
d <- read.table("RS_FRS_invivoinvitro_cleanedJan2019_forRforDryad.txt", header = TRUE, fill = TRUE)

#Setting Time as factor
d$Time <- as.factor(d$Time)

#Adding 1 to CFU values due to values being too low to take log (below 1)
d$CFUperml_limitadjusted <- (d$CFUperml_limitadjusted + 1)

#Subsetting experimental samples without controls
dsamples <- subset(d, Res_category != "NA")

#Subsetting plant and invitro data
dplant <- subset(dsamples, Environment != "Media")
dmedia <- subset(dsamples, Environment != "Plant")

#Adding labeller
facetlabels <- c('R' = "Resistant", 'S' = "Sensitive")
phagelabels <- c('N' = "No", 'Y' = "Yes")

#Plotting plant data
plantdataCFU <- ggplot(dplant, aes(x=Time, y=log10(CFUperml_limitadjusted), fill=Phage)) +
  scale_fill_manual(values=c("#000000", "#B3B3B3"), labels=c("No","Yes")) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.45,
               position=position_dodge(0.8)) +
  facet_wrap(~ Res_category, labeller=labeller(Res_category = facetlabels)) +
  labs(y= "Bacterial density (log10 CFU per 1mm2 leaf tissue)", x= "Time (Hours post-inoculation)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0,7) +
  scale_x_discrete(breaks=c("0","24","72"),
                   labels=c("1", "24", "72")) +
  labs(subtitle = "In planta growth") +
  geom_hline(aes(yintercept=log10(1.59)))

#Plotting in vitro data
mediadataCFU <- ggplot(dmedia, aes(x=Time, y=log10(CFUperml_limitadjusted), fill=Phage)) +
  scale_fill_manual(values=c("#000000", "#B3B3B3"), labels=c("No","Yes")) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.45,
               position=position_dodge(0.8)) +
  facet_wrap(~ Res_category, labeller=labeller(Res_category = facetlabels)) +
  labs(y= "Bacterial density (log10 CFU per ml)", x= "Time (Hours post-inoculation)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0,10) +
  scale_x_discrete(breaks=c("0","24","72"),
                   labels=c("1", "24", "72")) +
  labs(subtitle= "In vitro growth") +
  geom_hline(aes(yintercept=log10(167.5)))

#Note: +1 was added to value of threshold for geom_hline in both plots, to match data

#Combined figure with both panels
combinedCFU <- plot_grid(mediadataCFU, plantdataCFU, labels = c("A", "B"))

#Saving plot as .eps file
ggsave("Supp_Figure_2.eps", plot=combinedCFU, width= 400, height = 200, units=c("mm"), dpi=600)

########################SUPP FIG 3############################
d <- read.table("RS_FRS_invivoinvitro_cleanedJan2019_forRforDryad.txt", header = TRUE, fill = TRUE)

#Setting Time as factor
d$Time <- as.factor(d$Time)

#Adding 1 to CFU per ml, limit adjusted values
d$CFUperml_limitadjusted <- (d$CFUperml_limitadjusted + 1)

#Adding 1 to ddPCR values
d$AvgB_limitadjusted <- (d$AvgB_limitadjusted + 1)

#Subsetting experimental samples separate from controls
dsamples <- subset(d, Res_category != "NA")

#Removing any samples that have a CFU or AvgB value at the limit of detection
dsamplesforcorr <- subset(dsamples, CFUorAvgBatlimit_forcorr != "Y")

#Correlation plot of ddPCR and CFU values
corr <- ggplot(dsamplesforcorr, aes(x=log10(CFUperml_limitadjusted), y=log10(AvgB_limitadjusted))) +
  geom_point(aes(x=log10(CFUperml_limitadjusted), y=log10(AvgB_limitadjusted), fill=Phage), size=3, pch=21, colour="black") +
  scale_fill_manual(values = c("#000000", "#B3B3B3"),
                    name="Phage")  +           
  geom_smooth(method="lm", se=F, fullrange=T, colour="black", linetype="dashed") +
  geom_abline(intercept=0, slope=1) +
  stat_cor(method="pearson") +
  labs(y= "ddPCR concentration (log10 copies per ml culture or 1mm2 leaf tissue)", x= "CFU concentration (log10 CFU per ml culture or 1mm2 leaf tissue)") +
  scale_y_continuous(limits = c(0, 10)) +
  scale_x_continuous(limits = c(0, 10))

#Saving plot as .eps file
ggsave("Supp_Figure_3.eps", plot=corr, width= 200, height = 200, units=c("mm"), dpi=600, device = cairo_ps)
