# R script to run the metabolic network model for larger G. barretti sponges (Leys et al. 2018, doi: 10.1002/lno.10623)
# Date: June 2020
# This script belongs to De Kluijver et al. 2020 doi: 10.3389/fmars.2020.596251

#load packages
require(LIM)
require(dplyr)
require(ggplot2)
require(reshape2)
require(Hmisc)

#read the input file
file    <- "MN_Gbarretti_larger_fit.input"     				  # LIM inputfile for larger sponges
readLIM     <- Read(file)                                     # this reads the LIM file, but doesn't set it up

LIM<-testsetup(readLIM)                                       # first run LIM_additionalMN_setup.R, that contains the function testsetup
  x0      <- Ldei(LIM)                                        # least distance solution
  Ldei(LIM)$X                                                 # show least distance solution


# sensitivity analyses on parameters and data
Xss<-NULL                                                     # create a file to store modeloutput
RCs<-NULL                                                     # create a file to store parameters

#read parameters from LIM file
#stoichiometries
CN_Spo<-readLIM$pars$val[which(readLIM$pars$name == "CN_Spo")] 
CN_DOM<-readLIM$pars$val[which(readLIM$pars$name == "CN_DOM")] 
CN_Bac<-readLIM$pars$val[which(readLIM$pars$name == "CN_Bac")]

#efficiencies
fp_N_DOM<-readLIM$pars$val[which(readLIM$pars$name == "fp_N_DOM")]
fp_N_Bac<-readLIM$pars$val[which(readLIM$pars$name == "fp_N_Bac")]
Cfix_nit1<-readLIM$pars$val[which(readLIM$pars$name == "Cfix_nit1")]
Cfix_nit2<-readLIM$pars$val[which(readLIM$pars$name == "Cfix_nit2")]
Cfix_ana<-readLIM$pars$val[which(readLIM$pars$name == "Cfix_ana")]
fp_N_Denit <-readLIM$pars$val[which(readLIM$pars$name == "fp_N_Denit")]
fp_N_DNRA <-readLIM$pars$val[which(readLIM$pars$name == "fp_N_DNRA")]
RQ<-readLIM$pars$val[which(readLIM$pars$name=="RQ")]

#measured fluxes
Oxyupt<-readLIM$pars$val[which(readLIM$pars$name == "Oxyupt")]
DOCupt<-readLIM$pars$val[which(readLIM$pars$name == "DOCupt")]
BacCupt<-readLIM$pars$val[which(readLIM$pars$name == "BacCupt")]
NH4upt<-readLIM$pars$val[which(readLIM$pars$name == "NH4upt")]
NO2upt<-readLIM$pars$val[which(readLIM$pars$name == "NO2upt")]
NO3rel<-readLIM$pars$val[which(readLIM$pars$name == "NO3rel")]


#uncertainty - variability analysis with ranges
  
for (i in 1:1000000){ 
  #create ranges
    CN_Spo<-c(CN_Spo, runif(1, min = 3.9, max = 4.5))       # measured range (mean +/- sd)
    CN_DOM<-c(CN_DOM, runif(1, min = 6, max = 11))          # range from labile to refractory
    CN_Bac<-c(CN_Bac, runif(1,min = 4, max = 6 ))           # range for marine bacteria
    
    fp_N_DOM<-c(fp_N_DOM, runif(1,min=0.2,max=0.8))         # range 
    fp_N_Bac<-c(fp_N_Bac, runif(1,min=0.2,max=0.8))         # range
    fp_N_Denit<-c(fp_N_Denit,runif(1,min=0.3,max=0.6))      # range
    fp_N_DNRA<-c(fp_N_DNRA,runif(1,min=0.3,max=0.6))        # range
    
    Cfix_nit1<-c(Cfix_nit1,runif(1,min=0.04,max=0.11))      # literature range
    Cfix_nit2<-c(Cfix_nit2,runif(1,min=0.01,max=0.03))      # literature range
    Cfix_ana<-c(Cfix_ana,runif(1,min=0.03,max=0.1))         # literature range
    RQ<-c(RQ,runif(1,min=1,max=1.16))                       # range from C only to redfield OM
    
    Oxyupt<-c(Oxyupt, runif(1, min = 1.9, max = 13.5))      # measured Oxyupt (mean +/- sd)
    DOCupt<-c(DOCupt, runif(1, min = 3.2, max = 7.3))       # measured DOCupt (mean +/- sd)
    BacCupt<-c(BacCupt, runif(1, min = 0.19, max = 0.41))   # measured Bacupt (mean +/- sd)
    NH4upt<-c(NH4upt, runif(1, min = 0.0064, max = 0.013))  # measured NH4upt (mean +/- sd)
    NO2upt<-c(NO2upt, runif(1, min = 0.0055, max = 0.031))  # measured NO2upt (mean +/- sd)
    NO3rel<-c(NO3rel, runif(1, min = -1.35, max = -0.28))   # measured NO3rel (mean +/- sd)
    
    ##read values from ranges into LIM input 
    readLIM$pars$val[which(readLIM$pars$name == "CN_Spo")] <- CN_Spo[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "NC_Spo")] <- 1/CN_Spo[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "CN_DOM")] <- CN_DOM[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "CN_Bac")] <- CN_Bac[i+1]
  
    readLIM$pars$val[which(readLIM$pars$name == "Oxyupt")] <- Oxyupt[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "DOCupt")] <- DOCupt[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "BacCupt")] <- BacCupt[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "NH4upt")] <- NH4upt[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "NO2upt")] <- NO2upt[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "NO3rel")] <- NO3rel[i+1]
    
    readLIM$pars$val[which(readLIM$pars$name =="fp_N_DOM")] <- fp_N_DOM[i+1]
    readLIM$pars$val[which(readLIM$pars$name =="fp_N_Bac")] <- fp_N_Bac[i+1]
    readLIM$pars$val[which(readLIM$pars$name =="fex_N_DOM")] <- 1-fp_N_DOM[i+1]
    readLIM$pars$val[which(readLIM$pars$name =="fex_N_Bac")] <- 1-fp_N_Bac[i+1]
    
    readLIM$pars$val[which(readLIM$pars$name =="fp_N_Denit")] <- fp_N_Denit[i+1]
    readLIM$pars$val[which(readLIM$pars$name =="fp_N_DNRA")] <- fp_N_DNRA[i+1]
    readLIM$pars$val[which(readLIM$pars$name =="fex_N_Denit")] <- 1-fp_N_Denit[i+1]
    readLIM$pars$val[which(readLIM$pars$name =="fex_N_DNRA")] <- 1-fp_N_DNRA[i+1]
    
    readLIM$pars$val[which(readLIM$pars$name == "Cfix_nit1")]<- Cfix_nit1[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "Cfix_nit2")]<- Cfix_nit2[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "Cfix_ana")]<-  Cfix_ana[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "RQ")]<-RQ[i+1]
    
    fr_C_Bac<-1- CN_Spo[i+1]*fp_N_Bac[i+1]/CN_Bac[i+1]
    if(fr_C_Bac>0){  readLIM$pars$val[which(readLIM$pars$name =="fr_C_Bac")] <-fr_C_Bac}
    else{ next}
    readLIM$pars$val[which(readLIM$pars$name =="fr_C_DOM")] <- 1 - CN_Spo[i+1]*fp_N_DOM[i+1]/CN_DOM[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "fr_C_Denit")] <- 1 - CN_Spo[i+1]*fp_N_Denit[i+1]/CN_DOM[i+1]
    readLIM$pars$val[which(readLIM$pars$name == "fr_C_DNRA")] <- 1 - CN_Spo[i+1]*fp_N_DNRA[i+1]/CN_DOM[i+1]
    
  LIM<-testsetup(readLIM)                                       #setup LIM
  if(Ldei(LIM)$IsError==TRUE) next                              #skip when combination is not feasible
  RCs<-rbind(RCs,LIM$Parameters)                                #write parameters in file
  Xss<-rbind(Xss, Xsample(LIM, iter=500))                       #write feasible solutions in file
  
  NGE_C<-Xss[,"r10"]*CN_Spo[i+1]/(DOCupt[i+1]+BacCupt[i+1])     #calculate and add CNGE 
  Xss2<-cbind(Xss,NGE_C)

  }

# analyse parameters from uncertainty analysis
write.csv(RCs,"allRCs.csv")
RCs_unique<-  RCs%>% group_by(name) %>% distinct(val)%>% group_by(name) %>% filter(n()>1)     # select variable parameters

# calculate carbon production efficiencies (1- respired fractions)
fp_C_DOM<-data.frame(name="fp_C_DOM",val=1-RCs_unique$val[RCs_uniqu$name=="fr_C_DOM"])
fp_C_Bac<-data_frame(name="fp_C_Bac",val=1-RCs_unique$val[RCs_unique$name=="fr_C_Bac"])
fp_C_Denit<-data.frame(name="fp_C_Denit",val=1-RCs_unique$val[RCs_unique$name=="fr_C_Denit"])
fp_C_DNRA<-data.frame(name="fp_C_DNRA",val=1-RCs_unique$val[RCs_unique$name=="fr_C_DNRA"])
RCs_unique2<-bind_rows(RCs_unique,fp_C_DOM,fp_C_Bac,fp_C_Denit,fp_C_DNRA)                     # add production efficiencies

# read parameters that are later used for post-proces calculations and plots
CNDOM<-RCs_unique[RCs_unique$name=="CN_DOM",]
CNSpo<-RCs_unique[RCs_unique$name=="CN_Spo",]
r11_tDOCupt<-RCs_unique[RCs_unique$name=="DOCupt",]
fr_C_Denit<-RCs_unique[RCs_unique$name=="fr_C_Denit",]
fr_C_DNRA<-RCs_unique[RCs_unique$name=="fr_C_DNRA",]
r1_BacCupt<-RCs_unique[RCs_unique$name=="BacCupt",]
fp_N_DOM<-RCs_unique[RCs_unique$name == "fp_N_DOM",]
fp_N_Denit<-RCs_unique[RCs_unique$name == "fp_N_Denit",]
fp_N_DNRA<-RCs_unique[RCs_unique$name == "fp_N_DNRA",]

##statistics on parameters
RCstat<-RCs_unique2 %>% group_by(name) %>% summarise (mean = mean(val),sd=sd(val),min = quantile(val,probs=0.05),median=median(val),
                                                              max = quantile(val,probs=0.95),n=n())
RCstat <- RCstat %>%  bind_rows(RCstat %>%
              filter(name %in% c("BacCupt", "DOCupt")) %>%
              summarise_if(is.numeric, sum) %>%
              mutate(name = "OCupt")) %>%
              bind_rows(RCstatL %>%
              filter(name %in% c("NO2upt", "NH4upt")) %>%
              summarise_if(is.numeric, sum) %>%
              mutate(name = "NO2-NH4upt"))               

RCstat

##plot parameters
theme_set(theme_bw(15))
ggplot(RCs_unique2,aes(x=val))+geom_histogram()+facet_wrap(~name,scales="free")       #all parameters

#create figure 5 (histograms of parameters)
RCplot<-ggplot(RCs_unique2,aes(x=val))+geom_histogram()+facet_wrap(~name,scales="free")+ggtitle("larger sponges") #,fill=factor(name)
RCplot %+% subset(RCs_uniqueL2,name %in% c("CN_DOM","fp_N_DOM","fp_C_DOM")) #save with 1000-400
ggsave("DOMparameters.jpg",width=18, height=7,units="cm", dpi=300,scale=1.5)

# analyse model solutions of fluxes

colnames(Xss2)<-c("BacNupt","DONupt","denit","DNRA","nitr_NH4","nitr_NO2","biomassprod_nitrif","anammox","biomassprod_anammox","Prod","tDONupt","cOxyupt","cNH4upt","cNO3rel","cNO2upt","CO2rel","N2rel","NGE_C")
output<-as.data.frame(Xss2)

# average per 500 (xsample) and calculate with reaction coefficients
output2<-sapply(output[1:18], function(x) colMeans(matrix(x, nrow=500)))
output2<-data.frame(output2)
output2$CNDOM<-CNDOM$val
output2$r11_tDOCupt<-r11_tDOCupt$val
output2$CNSpo<-CNSpo$val
output2$fr_C_Denit<-fr_C_Denit$val
output2$fr_C_DNRA<-fr_C_DNRA$val
output2$fp_C_DOM<-fp_C_DOM$val
output2$r1_BacCupt<-r1_BacCupt$val

#postproces calculations
output3<-output2 %>% 
  mutate(RQ = abs(cOxyupt/CO2rel)) %>%
  mutate(ONupt = BacNupt + tDONupt) %>%
  mutate(Nfix = biomassprod_nitrif+biomassprod_anammox) %>%
  mutate(fracfix = Nfix / Prod) %>%
  mutate(r10_totCprod = Prod * CNSpo) %>%
  mutate(r7_r9_Cfix = Nfix * CNSpo ) %>%
  mutate(r3_denit_N = denit * 0.8 * CNDOM*fr_C_Denit) %>%
  mutate(r4_DNRA_N = DNRA * 0.56 * CNDOM*fr_C_DNRA) %>%
  mutate(r17_Ngas_loss = N2rel * 2) %>%
  mutate(NGEN = Prod/ONupt) %>%
  mutate(fracDenitNit= r3_denit_N/nitr_NO2) %>%
  mutate(fracDNRANit= r4_DNRA_N/nitr_NO2) %>%
  mutate(fracN2TotN = r17_Ngas_loss/(ONupt)) %>%
  mutate(fracOxDOM = DONupt/tDONupt) %>%
  mutate(fracDNRADOM = DNRA/tDONupt) %>%
  mutate(fracDenitDOM = denit/tDONupt) %>%
  mutate(OxyNO3 = cOxyupt/cNO3rel) %>%
  mutate(r2_DOCupt = DONupt*CNDOM) %>%
  mutate(r3_denitDOCupt=denit*CNDOM) %>%
  mutate(r4_DNRADOCupt=DNRA*CNDOM)  %>%
  mutate(r2_CProd_DOC = DONupt*CNDOM*fp_C_DOM)  %>%
  mutate(r3_CProd_denit = denit*CNDOM*fp_C_Denit$val) %>%
  mutate(r4_CProd_DNRA = DNRA*CNDOM*fp_C_DNRA$val) %>%
  mutate(fracOxDOM = DONupt/tDONupt)  %>%
  mutate(fracNH4DNRADenit = r4_DNRA_N/nitr_NH4) %>%
  mutate(OxyNit = (nitr_NH4*1.5+nitr_NO2*0.5)/cOxyupt) %>%
  mutate(r8_anammox = anammox * 2) %>% 
  mutate(Nsymb_Prod = (r3_CProd_denit + r4_CProd_DNRA+r7_r9_Cfix)/r10_totCprod) %>%
  mutate(DOCprodeff = (r2_CProd_DOC+r3_CProd_denit+r4_CProd_DNRA)/(tDONupt*CNDOM)) %>%
  mutate(DONprodeff = (DONupt*fp_N_DOM$val+denit*fp_N_Denit$val+DNRA*fp_N_Denit$val)/tDONupt)

write.csv(output3,"allmodeloutput.csv")           

#look at correlations between modelled fluxes
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) ) #redefine y-axis; x-axis stays the same
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

ss<-subset(output2,select=c(BacNupt,DONupt,Prod,denit,DNRA,nitr_NH3,nitr_NO2,anammox,cOxyupt,cNH4upt,cNO3rel,cNO2upt,CO2rel,N2rel))
ss<-ss[c(1:5000),]
pairs(ss, upper.panel = NULL, diag.panel = panel.hist,
      pch = ".", cex = 2,  cex.labels = 1.3,main = "Reaction network large sponges")

# statistics on modeloutput
output4<-data.frame(colMeans(output3),apply(output3,2,sd),apply(output3,2,quantile, probs = 0.05),apply(output3,2,quantile,probs=0.95))
colnames(output4)<-c("mean","sd","min","max")
write.csv(output4,file="outputstats.csv")

#bind statistics RCs and output in one file
output4$inoutput<-"calculated"
output4$name<-row.names(output4)
RCstat_SL$inoutput<-"measured"
dfall<-bind_rows(output4,RCstat) 
dfall
write.csv(dfall_SL,file="outputGeodia_all.csv")

# plot C fluxes (Fig. 3)
output4$group[output4$name %in% c("r10_totCprod","r2_CProd_DOC","r3_CProd_denit","r4_CProd_DNRA","r7_r9_Cfix")]<-"Cproduction"
output4$group[output4$name %in% c("r11_tDOCupt","r1_BacCupt","r2_DOCupt","r3_denitDOCupt","r4_DNRADOCupt")]<-"Cassimilation"
output4$name[which(output4$name=="CO2rel")]<-"r16_CO2rel"
output4$group[output4$name %in% c("r16_CO2rel")]<-"Crelease"

Cbudget2<-output4[output4$group %in% c("Cproduction","Cassimilation","Crelease"),]

Cbudget2$name <-factor(Cbudget2$name, levels = c("r1_BacCupt","r2_DOCupt","r2_CProd_DOC","r3_denitDOCupt","r3_CProd_denit",
                                                 "r4_DNRADOCupt","r4_CProd_DNRA","r7_r9_Cfix","r10_totCprod","r11_tDOCupt","r16_CO2rel"))

ggplot(Cbudget2,aes(x=name,y=abs(mean),fill=group))+geom_bar(stat="identity")+
  geom_errorbar(aes(ymin = abs(mean) - sd, ymax = abs(mean) + sd), width=0.2)+
  scale_fill_manual(values=c("darkorange1","#009E73","blue"))+
                    labs(x="",y = expression(paste(mu,"mol ",cm^{-3},d^{-1})),fill="",title="larger sponges")+
  theme(axis.text.x = element_text(angle = 90))

ggsave("Cbudget.jpg",width=18,height=9,units="cm", dpi=300,scale=1.5)

# plot N fluxes (Fig. 4)
output4$name[which(output4$name=="nitr_NH4")]<-"r5_NH4_oxidation"
output4$name[which(output4$name=="nitr_NO2")]<-"r6_NO2_oxidation"
output4$group[output4$name %in% c("r3_denit_N","r4_DNRA_N","r8_anammox")]<-"anoxic"
output4$group[output4$name %in% c("r5_NH4_oxidation","r6_NO2_oxidation")]<-"oxic"
Nbudget2<-output4[output4$name %in% c("r3_denit_N","r4_DNRA_N","r5_NH4_oxidation","r6_NO2_oxidation","r8_anammox"),]
Nbudget2$name <-factor(Nbudget2$name, levels = c("r3_denit_N","r4_DNRA_N","r5_NH4_oxidation","r6_NO2_oxidation","r8_anammox"))
ggplot(Nbudget2,aes(x=name,y=abs(mean),fill=group))+geom_bar(stat="identity")+
  geom_errorbar(aes(ymin = abs(mean) - sd, ymax = abs(mean) + sd), width=0.2)+
  scale_fill_manual(values=c("darkorange1", "blue2"))+
  labs(x="",y = expression(paste(mu,"mol ",cm^{-3},d^{-1})),fill="",title="larger sponges")+
  theme(axis.text.x = element_text(angle = 90))

ggsave("Nbudget.jpg",width=10.5,height=9,units="cm", dpi=300,scale=1.5)
   