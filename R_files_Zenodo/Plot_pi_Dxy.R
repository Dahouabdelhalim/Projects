#load the required library 
library(ggplot2)
library(tidyr)

#read in pi estimates along the X chromosome for each population (these were obtained using the "pi_dxy.sh" script)
SOR<-read.csv(file="per_population_pi/SOR.pi.txt", header = TRUE, sep=' ')
MAR<-read.csv(file="per_population_pi/MAR.pi.txt", header = TRUE, sep=' ')
LHA<-read.csv(file="per_population_pi/LHA.pi.txt", header = TRUE, sep=' ')
LHB<-read.csv(file="per_population_pi/LHB.pi.txt", header = TRUE, sep=' ')
GUA<-read.csv(file="per_population_pi/GUA.pi.txt", header = TRUE, sep=' ')
JIC<-read.csv(file="per_population_pi/JIC.pi.txt", header = TRUE, sep=' ')
CAB<-read.csv(file="per_population_pi/CAB.pi.txt", header = TRUE, sep=' ')
ESM<-read.csv(file="per_population_pi/ESM.pi.txt", header = TRUE, sep=' ')
POR<-read.csv(file="per_population_pi/POR.pi.txt", header = TRUE, sep=' ')
BAH<-read.csv(file="per_population_pi/BAH.pi.txt", header = TRUE, sep=' ')
ALA<-read.csv(file="per_population_pi/ALA.pi.txt", header = TRUE, sep=' ')
BRE<-read.csv(file="per_population_pi/BRE.pi.txt", header = TRUE, sep=' ')
BRO<-read.csv(file="per_population_pi/BRO.pi.txt", header = TRUE, sep=' ')
CIT<-read.csv(file="per_population_pi/CIT.pi.txt", header = TRUE, sep=' ')
CMB<-read.csv(file="per_population_pi/CMB.pi.txt", header = TRUE, sep=' ')
COL<-read.csv(file="per_population_pi/COL.pi.txt", header = TRUE, sep=' ')
DUV<-read.csv(file="per_population_pi/DUV.pi.txt", header = TRUE, sep=' ')
FLA<-read.csv(file="per_population_pi/FLA.pi.txt", header = TRUE, sep=' ')
GLA<-read.csv(file="per_population_pi/GLA.pi.txt", header = TRUE, sep=' ')
HAM<-read.csv(file="per_population_pi/HAM.pi.txt", header = TRUE, sep=' ')
HEN<-read.csv(file="per_population_pi/HEN.pi.txt", header = TRUE, sep=' ')
HIG<-read.csv(file="per_population_pi/HIG.pi.txt", header = TRUE, sep=' ')
LAK<-read.csv(file="per_population_pi/LAK.pi.txt", header = TRUE, sep=' ')
LEE<-read.csv(file="per_population_pi/LEE.pi.txt", header = TRUE, sep=' ')
LEV<-read.csv(file="per_population_pi/LEV.pi.txt", header = TRUE, sep=' ')
LOW<-read.csv(file="per_population_pi/LOW.pi.txt", header = TRUE, sep=' ')
MAN<-read.csv(file="per_population_pi/MAN.pi.txt", header = TRUE, sep=' ')
MARi<-read.csv(file="per_population_pi/MARi.pi.txt", header = TRUE, sep=' ')
MIA<-read.csv(file="per_population_pi/MIA.pi.txt", header = TRUE, sep=' ')
MON<-read.csv(file="per_population_pi/MON.pi.txt", header = TRUE, sep=' ')
ORA<-read.csv(file="per_population_pi/ORA.pi.txt", header = TRUE, sep=' ')
PAL<-read.csv(file="per_population_pi/PAL.pi.txt", header = TRUE, sep=' ')
POL<-read.csv(file="per_population_pi/POL.pi.txt", header = TRUE, sep=' ')
SAR<-read.csv(file="per_population_pi/SAR.pi.txt", header = TRUE, sep=' ')
STJ<-read.csv(file="per_population_pi/STJ.pi.txt", header = TRUE, sep=' ')
STL<-read.csv(file="per_population_pi/STL.pi.txt", header = TRUE, sep=' ')
STP<-read.csv(file="per_population_pi/STP.pi.txt", header = TRUE, sep=' ')
TAM<-read.csv(file="per_population_pi/TAM.pi.txt", header = TRUE, sep=' ')
TIF<-read.csv(file="per_population_pi/TIF.pi.txt", header = TRUE, sep=' ')
VOL<-read.csv(file="per_population_pi/VOL.pi.txt", header = TRUE, sep=' ')

###plot pi along chr7 for each population; the VOL population is plotted here as an example
ggplot(data = VOL,
       aes(x = mid, 
           y = pi_VOL))+
  annotate("rect",xmin=70000000, xmax=88000000, ymin=0, ymax=0.03, alpha=0.1, fill="black", size=0)+
  geom_point(size=1.5, colour="black", shape=21)+
  geom_smooth(span=0.2,se = FALSE, colour="red", size=0.8, method = "loess")+
  geom_vline(xintercept = 21000000, linetype="dashed")+
  geom_vline(xintercept = 92000000, linetype="dashed")+
  labs(x = "Chromosome 7 position", y = "Nucleotide diversity") +
  #only plot pi values up to 0.03, to better visualize trends at the low end of nucleotide diversity distribution (as per Figs. S19-20)
  ylim(0,0.03)+
  theme_bw()+
  theme( 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))

#read in Dxy estimates calculated along the X chromosome for each genome scan population pair (these were obtained using the "pi_dxy.sh" script)
LHA_LHB<-read.csv(file="Dxy_unique_population_pairs/LHA_LHB.dxy.txt", header = TRUE, sep=' ')
MAR_GUA<-read.csv(file="Dxy_unique_population_pairs/MAR_GUA.dxy.txt", header = TRUE, sep=' ')
CAB_JIC<-read.csv(file="Dxy_unique_population_pairs/CAB_JIC.dxy.txt", header = TRUE, sep=' ')
ESM_POR<-read.csv(file="Dxy_unique_population_pairs/ESM_POR.dxy.txt", header = TRUE, sep=' ')
CMB_ORA<-read.csv(file="Dxy_unique_population_pairs/CMB_ORA.dxy.txt", header = TRUE, sep=' ')
BRE_LEE<-read.csv(file="Dxy_unique_population_pairs/BRE_LEE.dxy.txt", header = TRUE, sep=' ')
CIT_SAR<-read.csv(file="Dxy_unique_population_pairs/CIT_SAR.dxy.txt", header = TRUE, sep=' ')
ALA_POL<-read.csv(file="Dxy_unique_population_pairs/ALA_POL.dxy.txt", header = TRUE, sep=' ')
FLA_LAK<-read.csv(file="Dxy_unique_population_pairs/FLA_LAK.dxy.txt", header = TRUE, sep=' ')
LEV_TAM<-read.csv(file="Dxy_unique_population_pairs/LEV_TAM.dxy.txt", header = TRUE, sep=' ')
COL_DUV<-read.csv(file="Dxy_unique_population_pairs/COL_DUV.dxy.txt", header = TRUE, sep=' ')
HIG_MIA<-read.csv(file="Dxy_unique_population_pairs/HIG_MIA.dxy.txt", header = TRUE, sep=' ')
GLA_PAL<-read.csv(file="Dxy_unique_population_pairs/GLA_PAL.dxy.txt", header = TRUE, sep=' ')
HEN_STJ<-read.csv(file="Dxy_unique_population_pairs/HEN_STJ.dxy.txt", header = TRUE, sep=' ')
HAM_MARi<-read.csv(file="Dxy_unique_population_pairs/HAM_MARi.dxy.txt", header = TRUE, sep=' ')
MAN_STP<-read.csv(file="Dxy_unique_population_pairs/MAN_STP.dxy.txt", header = TRUE, sep=' ')
MON_VOL<-read.csv(file="Dxy_unique_population_pairs/MON_VOL.dxy.txt", header = TRUE, sep=' ')
LOW_STL<-read.csv(file="Dxy_unique_population_pairs/LOW_STL.dxy.txt", header = TRUE, sep=' ')
BRO_TIF<-read.csv(file="Dxy_unique_population_pairs/BRO_TIF.dxy.txt", header = TRUE, sep=' ')

###plot Dxy along chr7 for each population pair; the BRO_TIF population pair is plotted here as an example
ggplot(data = BRO_TIF,
       aes(x = mid, 
           y = dxy_BRO_TIF))+
  annotate("rect",xmin=70000000, xmax=88000000, ymin=0, ymax=0.03, alpha=0.1, fill="black", size=0)+
  geom_point(size=1.5, colour="black", shape=21)+
  geom_smooth(span=0.2,se = FALSE, colour="red", size=0.8, method = "loess")+
  geom_vline(xintercept = 21000000, linetype="dashed")+
  geom_vline(xintercept = 92000000, linetype="dashed")+
  labs(x = "Chromosome 7 position", y = "Dxy") +
  #only plot Dxy values up to 0.03, to better visualize trends at the low end of nucleotide diversity distribution (as per Figs. S21-22)
  ylim(0,0.03)+
  theme_bw()+
  theme( 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))


#analyses for the X chromosome divergent locus vs. a candidate neutral locus of the same size, for each range
pi_all<-read.csv(file="pi_divergent_locus_vs_rest_of_X_no_PAR/pi_combined.txt", header = TRUE, sep=' ')
pi_all$locus = factor(pi_all$locus, levels=c("background", "divergent"))
pi_all$range = factor(pi_all$range, levels=c("native", "invasive"))

#also get means for each category, used for plotting
means<-as.data.frame(aggregate(x=pi_all$pi,by=list(pi_all$locus,pi_all$range),FUN=mean, na.rm=TRUE,na.action=NULL))
names(means)[names(means) == "Group.2"] <- "range"
names(means)[names(means) == "Group.1"] <- "locus"
names(means)[names(means) == "x"] <- "mean_pi"


#plot all windows and averages per category
ggplot(data = pi_all,
       aes(x = locus, 
           y = pi, colour=range))+
  geom_jitter(width=0.2, shape=1, alpha=0.2)+
  geom_point(data=means, aes(x=locus, y=mean_pi), colour="lightgrey", size=3.5)+
  ylim(0,0.1)+
  labs(x = NULL, y = NULL) +
  theme_bw()+
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank())+
  facet_wrap(~range)+
  scale_colour_manual(values = c("red","black"))


#test for differences in pi for each range, divergent locus vs. background
#also test for differences between the two ranges, at the divergent locus categories
#input files were obtained using the "pi_dxy.sh" (average_pi_native.txt and average_pi_invasive.txt); we annotated these file to add population IDs  

Native<-read.csv(file="Native_pi_neutral_vs_divergent.csv",header = TRUE, sep = ,)
Invasive<-read.csv(file="Invasive_pi_neutral_vs_divergent.csv",header = TRUE, sep = ,)

native_divergent <- Native[Native$locus == 'divergent',]
native_divergent_pi<-native_divergent$Pi
native_neutral <- Native[Native$locus == 'background',]
native_neutral_pi<-native_neutral$Pi

invasive_divergent <- Invasive[Invasive$locus == 'divergent',]
invasive_divergent_pi<-invasive_divergent$Pi
invasive_neutral <- Invasive[Invasive$locus == 'background',]
invasive_neutral_pi<-invasive_neutral$Pi

#test whether divergent locus pi averages are smaller than neutral locus pi averages, for each range
wilcox.test(native_divergent_pi, native_neutral_pi,alternative = "less", paired = TRUE)
wilcox.test(invasive_divergent_pi, invasive_neutral_pi,alternative = "less", paired = TRUE)

#test whether native range pi averages are smaller than invasive range pi averages, for each locus category
wilcox.test(native_divergent_pi, invasive_divergent_pi,alternative = "less")
wilcox.test(native_neutral_pi, invasive_neutral_pi,alternative = "less")

#adjust P values for multiple comparisons using Bonferroni
pvalues<-c("0.0009766","9.313e-10","2.359e-09","1.999e-05")
p.adjust(pvalues, method="bonferroni")


