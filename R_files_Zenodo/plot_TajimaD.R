#load the required library 
library(ggplot2)

##read in data files (these are obtained using the "TajimaD_windowed.sh" script)
MAR<-read.csv(file="MAR.Tajima.D.txt", header = TRUE, sep='\\t')
MAR_CHR7<-MAR[MAR$CHROM == '7',]
SOR<-read.csv(file="SOR.Tajima.D.txt", header = TRUE, sep='\\t')
SOR_CHR7<-SOR[SOR$CHROM == '7',]
LHA<-read.csv(file="LHA.Tajima.D.txt", header = TRUE, sep='\\t')
LHA_CHR7<-LHA[LHA$CHROM == '7',]
LHB<-read.csv(file="LHB.Tajima.D.txt", header = TRUE, sep='\\t')
LHB_CHR7<-LHB[LHB$CHROM == '7',]
GUA<-read.csv(file="GUA.Tajima.D.txt", header = TRUE, sep='\\t')
GUA_CHR7<-GUA[GUA$CHROM == '7',]
JIC<-read.csv(file="JIC.Tajima.D.txt", header = TRUE, sep='\\t')
JIC_CHR7<-JIC[JIC$CHROM == '7',]
CAB<-read.csv(file="CAB.Tajima.D.txt", header = TRUE, sep='\\t')
CAB_CHR7<-CAB[CAB$CHROM == '7',]
ESM<-read.csv(file="ESM.Tajima.D.txt", header = TRUE, sep='\\t')
ESM_CHR7<-ESM[ESM$CHROM == '7',]
POR<-read.csv(file="POR.Tajima.D.txt", header = TRUE, sep='\\t')
POR_CHR7<-POR[POR$CHROM == '7',]
BAH<-read.csv(file="BAH.Tajima.D.txt", header = TRUE, sep='\\t')
BAH_CHR7<-BAH[BAH$CHROM == '7',]


ALA<-read.csv(file="ALA.Tajima.D.txt", header = TRUE, sep='\\t')
ALA_CHR7<-ALA[ALA$CHROM == '7',]
BRE<-read.csv(file="BRE.Tajima.D.txt", header = TRUE, sep='\\t')
BRE_CHR7<-BRE[BRE$CHROM == '7',]
BRO<-read.csv(file="BRO.Tajima.D.txt", header = TRUE, sep='\\t')
BRO_CHR7<-BRO[BRO$CHROM == '7',]
CIT<-read.csv(file="CIT.Tajima.D.txt", header = TRUE, sep='\\t')
CIT_CHR7<-CIT[CIT$CHROM == '7',]
CMB<-read.csv(file="CMB.Tajima.D.txt", header = TRUE, sep='\\t')
CMB_CHR7<-CMB[CMB$CHROM == '7',]
COL<-read.csv(file="COL.Tajima.D.txt", header = TRUE, sep='\\t')
COL_CHR7<-COL[COL$CHROM == '7',]
DUV<-read.csv(file="DUV.Tajima.D.txt", header = TRUE, sep='\\t')
DUV_CHR7<-DUV[DUV$CHROM == '7',]
FLA<-read.csv(file="FLA.Tajima.D.txt", header = TRUE, sep='\\t')
FLA_CHR7<-FLA[FLA$CHROM == '7',]
GLA<-read.csv(file="GLA.Tajima.D.txt", header = TRUE, sep='\\t')
GLA_CHR7<-GLA[GLA$CHROM == '7',]
HAM<-read.csv(file="HAM.Tajima.D.txt", header = TRUE, sep='\\t')
HAM_CHR7<-HAM[HAM$CHROM == '7',]
HEN<-read.csv(file="HEN.Tajima.D.txt", header = TRUE, sep='\\t')
HEN_CHR7<-HEN[HEN$CHROM == '7',]
HIG<-read.csv(file="HIG.Tajima.D.txt", header = TRUE, sep='\\t')
HIG_CHR7<-HIG[HIG$CHROM == '7',]
LAK<-read.csv(file="LAK.Tajima.D.txt", header = TRUE, sep='\\t')
LAK_CHR7<-LAK[LAK$CHROM == '7',]
LEE<-read.csv(file="LEE.Tajima.D.txt", header = TRUE, sep='\\t')
LEE_CHR7<-LEE[LEE$CHROM == '7',]
LEV<-read.csv(file="LEV.Tajima.D.txt", header = TRUE, sep='\\t')
LEV_CHR7<-LEV[LEV$CHROM == '7',]
LOW<-read.csv(file="LOW.Tajima.D.txt", header = TRUE, sep='\\t')
LOW_CHR7<-LOW[LOW$CHROM == '7',]
MAN<-read.csv(file="MAN.Tajima.D.txt", header = TRUE, sep='\\t')
MAN_CHR7<-MAN[MAN$CHROM == '7',]
MAR<-read.csv(file="MAR.Tajima.D.txt", header = TRUE, sep='\\t')
MAR_CHR7<-MAR[MAR$CHROM == '7',]
MIA<-read.csv(file="MIA.Tajima.D.txt", header = TRUE, sep='\\t')
MIA_CHR7<-MIA[MIA$CHROM == '7',]
MON<-read.csv(file="MON.Tajima.D.txt", header = TRUE, sep='\\t')
MON_CHR7<-MON[MON$CHROM == '7',]
ORA<-read.csv(file="ORA.Tajima.D.txt", header = TRUE, sep='\\t')
ORA_CHR7<-ORA[ORA$CHROM == '7',]
PAL<-read.csv(file="PAL.Tajima.D.txt", header = TRUE, sep='\\t')
PAL_CHR7<-PAL[PAL$CHROM == '7',]
POL<-read.csv(file="POL.Tajima.D.txt", header = TRUE, sep='\\t')
POL_CHR7<-POL[POL$CHROM == '7',]
SAR<-read.csv(file="SAR.Tajima.D.txt", header = TRUE, sep='\\t')
SAR_CHR7<-SAR[SAR$CHROM == '7',]
STJ<-read.csv(file="STJ.Tajima.D.txt", header = TRUE, sep='\\t')
STJ_CHR7<-STJ[STJ$CHROM == '7',]
STL<-read.csv(file="STL.Tajima.D.txt", header = TRUE, sep='\\t')
STL_CHR7<-STL[STL$CHROM == '7',]
STP<-read.csv(file="STP.Tajima.D.txt", header = TRUE, sep='\\t')
STP_CHR7<-STP[STP$CHROM == '7',]
TAM<-read.csv(file="TAM.Tajima.D.txt", header = TRUE, sep='\\t')
TAM_CHR7<-TAM[TAM$CHROM == '7',]
TIF<-read.csv(file="TIF.Tajima.D.txt", header = TRUE, sep='\\t')
TIF_CHR7<-TIF[TIF$CHROM == '7',]
VOL<-read.csv(file="VOL.Tajima.D.txt", header = TRUE, sep='\\t')
VOL_CHR7<-VOL[VOL$CHROM == '7',]

###plot TajimaD along chr7 for each population; the SOR population is plotted here as an example
ggplot(data = SOR_CHR7,
       aes(x = BIN_START, 
           y = TajimaD))+
  annotate("rect",xmin=70000000, xmax=88000000, ymin=-2.5, ymax=2.5, alpha=0.1, fill="black", size=0)+
  geom_point(size=1.5, colour="black", shape=21)+
  geom_smooth(span=0.2,se = FALSE, colour="red", size=0.8, method = "loess")+
  geom_vline(xintercept = 21000000, linetype="dashed")+
  geom_vline(xintercept = 92000000, linetype="dashed")+
  labs(x = "\\nChromosome 7 position", y = "Tajima D") +
  theme_bw()+
  theme( 
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))

###plot average Tajima D for the native and invasive ranges along chromosome 7 (Fig. 2C)
Cuba_average_TajimaD_CHR7<-read.csv(file="native.chr7.average_D.50d.with_markers.txt", header = TRUE, sep=' ')
Florida_average_TajimaD_CHR7<-read.csv(file="invasive.chr7.average_D.50d.with_markers.txt", header = TRUE, sep=' ')


###average Tajima D for chromosome 7
ggplot(data = Cuba_average_TajimaD_CHR7,
       aes(x = BP, 
           y = average_D))+
  annotate("rect",xmin=70000000, xmax=88000000, ymin=-1.2, ymax=1.2, alpha=0.1, fill="black", size=0)+
  geom_point(colour="#576cadff",shape=21, size=1)+
  geom_point(data=Florida_average_TajimaD_CHR7, colour="#576cadff", size=1)+
  geom_vline(xintercept = 21000000, linetype="dashed")+
  geom_vline(xintercept = 92000000, linetype="dashed")+
  geom_smooth(data=Florida_average_TajimaD_CHR7, span=0.2,se = FALSE, colour="black", size=1, method = "loess")+
  geom_smooth(span=0.2,se = FALSE, colour="red", size=1, method = "loess")+
  labs(x = NULL, y = NULL) +
  theme_bw()+
  theme(legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())


