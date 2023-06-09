FUNCANIMAL <- function(id, mat, vecteur, mydata, listparam, Chemin){
  
  
  tablo = mydata
  j=1
  for (i in (1:length(vecteur))){
    if (vecteur[i]==TRUE){
      tablo[i, "Valeur"] = paste(rep(mat[id,listparam[j]], tablo[i,"nbcols"]), collapse=" ")
      j = j+1
    }
  }

  write.table(tablo,"fichier_entrees.csv",sep=";", col.names = TRUE,
              row.names = FALSE)
  

  shell("C:/Users/acadero/AppData/Local/Programs/Python/Python36/python.exe Main_python.py") # Windows version
  # system("python3.4 Main_python.py")               # Linux version
  
  performances <- read.table("Performances_bilan.txt", sep ="", header = TRUE)
  Aliment_cout <- read.table("Aliment_cout.txt", sep ="", header = TRUE)
  


# Indicateurs economiques ---------------------------------------------------  

  tauxCC_fr=0.765       # rendement carcasse froid
  tauxCC_ch=0.79        # rendement carcasse chaud
  
  Coeffcor_TMP = 61.082/57.2648277921286
  #   La référence GTE 2015 la plus proche des caractéristiques de ta 
  #   simulation est celle des élevages naisseurs engraisseurs conduits en 7 
  #   bandes (toutes les 3 semaines). Si on prend les 213 élevages correspondant 
  #   au tiers supérieur de cette référence (les 33% meilleurs élevages triés selon la marge),
  #   on se rapproche des caractéristiques de ta simulation (moyenne de 232 truies présentes, engraissement entre 31 et 120kg avec GMQ 
  #   observé de 826g/j).
  #   On a TMP=61.082 et écart type 0,752.
  #   Sachant comme on le disait qu'une majorité des élevages allote au poids, 
  #   distribue un aliment biphase et rationne en engraissement (même si on ne 
  #   gère pas ces codifications en GTE).
  
  performances$TMPcor= performances$TMP*Coeffcor_TMP
  
  GDP1 <- read.csv(paste(Chemin,"/GDP_2015.csv", sep =""), header=TRUE, sep=";")
  
  performances$plusvalue<-0
  for (j in 3:17) {
    for (i in 3:15) {
      if (nrow(performances[(performances$TMPcor>GDP1$TMPmin[j])&
                            (performances$TMPcor<=GDP1$TMPmax[j])&
                            (performances$PVFin>(as.numeric(GDP1[1,i])/tauxCC_ch))&
                            (performances$PVFin<=(as.numeric(GDP1[2,i])/tauxCC_ch)),])>0) {  
        
        performances[(performances$TMPcor>GDP1$TMPmin[j])&
                       (performances$TMPcor<=GDP1$TMPmax[j])&
                       (performances$PVFin>(as.numeric(GDP1[1,i])/tauxCC_ch))&
                       (performances$PVFin<=(as.numeric(GDP1[2,i])/tauxCC_ch)),]$plusvalue<-GDP1[j,i]
      }
    }
  }
  
  performances$plusvalue[performances$Depart=="mort"]=0
  
  
  # Creation colonnes des indicateurs economiques
  prixbase=1.314                     # euros/kg en f?vrier 2014
  
  performances$prixvente <- ifelse (performances$Depart=="mort",
                                    0,
                                    performances$plusvalue+prixbase+0.02)
  # euros/ kg carcasse
  # 0.02 plus-value tracabilit?
  
  performances$produit<- (performances$prixvente*performances$PVFin*tauxCC_fr)   # euros/ carcasse de porc
  
  performances$marge_Eng <- (performances$produit - (2.35*25 + (performances$PVInit-25)*0.67)- performances$CoutAlimTot)
  # euros/ porc
  # Seules les charges alimentaires en engraissement 
  # et le prix d'achat du porcelet à la fin du post-sevrage sont déduits 
  # Cout de revient du porcelet de 25 kg en fevrier 2014 : 2,35 €/kg
  # Cout du kg vif supplémentaire : 0,67 €/kg
  
# Indicateurs environnementaux ----------------------------------------------  

  
#----------------------------------------------------------------------------
  #FATTENING
#----------------------------------------------------------------------------
  
  
  # N COMPOUNDS
  
  # Housing Emissions (kg manure N ex-animal *emission factor) 
  
  performances$N_NH3EmitH= (performances$RejeteN / 1000) * 0.24          # kg N-NH3 / pig
  # Rigolot et al.(2010b): EF is 0.24 kg N-NH3/kg N excreted
  
  performances$N_N2OEmitH= (performances$RejeteN / 1000) * 0.002         # kg N-N2O / pig
  # IPCC (2006): TABLE 10.21-Pit storage bellow animal confinements
  
  performances$N_NOxEmitH= (performances$RejeteN / 1000) * 0.002         # kg N-NOx / pig
  # Dammgen & Hutchings cited by IPCC (2006)
  
  performances$N_NO3EmitH=(performances$RejeteN / 1000) * 0              # kg N-NO3 / pig
  
  # N leaving after housing emissions (kg N / pig)
  performances$N_St = (performances$RejeteN/1000 - performances$N_NH3EmitH 
                       - performances$N_N2OEmitH - performances$N_NOxEmitH)
  
  
  # Storage Emissions (kg manure N ex-housing *emission factor) 
  
  performances$N_NH3EmitSt= performances$N_St * 0.1                      # kg N-NH3 / pig
  # Rigolot et al.(2010b): Table 4
  
  performances$N_N2OEmitSt= performances$N_St * 0.001                    # kg N-N2O / pig
  # IPCC (2006)
  
  performances$N_NOxEmitSt= performances$N_St * 0.005                    # kg N-NOx / pig
  # Dämmgen & Hutchings(2008)cited by IPCC(2006)
  
  # N leaving after storage emissions (kg N / pig)
  performances$N_field= (performances$N_St - performances$N_NH3EmitSt
                         - performances$N_N2OEmitSt - performances$N_NOxEmitSt)
  
  
  # Emissions during field application (kg manure N ex-storage *emission factor) 
  
  performances$N_NH3EmitFA= performances$N_field * 0.07                  # kg N-NH3 / pig
  # Andersen et al. (2001)
  
  performances$N_N2OEmitFA= performances$N_field * 0.01                  # kg N-N2O / pig
  # IPCC (2006)
  
  performances$N_NOxEmitFA= performances$N_field * 0.001                 # kg N-NOx / pig
  # Nemecek and Kägi (2007)
  
  # Content of manure after field application (kg N / pig)
  performances$N_AfterFA= (performances$N_field - performances$N_NH3EmitFA
                           - performances$N_N2OEmitFA - performances$N_NOxEmitFA)
  
  
  # Emissions (No3 leaching) after field application
  performances$N_NO3EmitFA= performances$N_AfterFA * 0.05                # kg N / pig
  # Garcia-Launay et al.(2014)
  
  # Extraloss as nitrates compared to mineral fertilizer
  performances$NO3_leach_increase = 0.05 * performances$N_field          # kg N / pig
  # (IPCC, 2006 ; eventuellement revoir/ publi Danemark)
  
  
  # Total emissions of N compounds (sum of all emissions, kg)
  performances$N_NH3Emit= ( performances$N_NH3EmitH
                            + performances$N_NH3EmitSt
                            + performances$N_NH3EmitFA )  
  performances$N_N2OEmit= ( performances$N_N2OEmitH
                            + performances$N_N2OEmitSt
                            + performances$N_N2OEmitFA )
  performances$N_NOxEmit= ( performances$N_NOxEmitH
                            + performances$N_NOxEmitSt
                            + performances$N_NOxEmitFA )
  performances$N_NO3Emit= ( performances$N_NO3EmitFA 
                            + performances$NO3_leach_increase)  
  
  # Total (kg N / pig)
  performances$N_Emit= (performances$N_NH3Emit+performances$N_N2OEmit 
                        +performances$N_NOxEmit+performances$N_NO3Emit)
  
  
  # Fertilizer N substitution - nutrient balance
  performances$Nsubs= 0.75 * performances$N_field                       # kg N / pig
  # (Nguyen et al., 2010)
  
  
  # -------------------------------------------------------------------------
  
  # CH4
  
  # CH4 emissions by enteric fermentation (Digested fibre intake)
  
  # Residu = MO - MAT - MG - Amidon - Sucres 
  # Resdig = Residu * coeff de digestibilit?
 
  performances$CH4_Enteric= (performances$Resdig * 1000        # g / pig
                             * 670                             # 670 J / g
                             / 1000000) / 55.65                # MJ / pig
  # enteric. Rigolot et al.(2010a)
  
  
  # CH4 emissions at storage 
  
  # performances$MOexc : kg of OM excreted
  performances$VS= performances$MOexc                    # kg, VS equals volatile solids cf equation IPCC 2006
  # Rigolot et al.(2010b)
  
  B0=0.42
  # IPCC (2006) gives a value of 0.45 m3CH4/kg VS for Western Europe (and and 0.29 for Latin América)
  
  TempMoy = 2/3*12+1/3*12.0
  # 2/3 indoor slurry storage at 20?C and 1/3 outdoor st at 12?C  (Garcia-Launay et al., 2014)
  
  performances$MCF=(0.1517*(TempMoy^2)-2.1569*TempMoy+24.634)/100
  # Garcia-Launay et al. (2014): MCFlisier=0.1517*(average temp^2)-2.1569*(average temp)+24.634
  # MCF : methan conversion factor
  
  performances$CH4_St=performances$VS*(B0*0.67)*performances$MCF
  # kg of CH4 storage. Rigolot et al. (2010b) from IPCC(2006); 0.67 for convert m3 to kg
  
  
  # Total of CH4 emissions
  performances$CH4Emit=performances$CH4_Enteric+performances$CH4_St 
  # kg of CH4 emitted by enteric fermentation and manure storage
  
  # ------------------------------------------------------------------------
  
  # K for spreading - without treatment
  
  performances$K_field= (performances$RejeteK                     # g / pig
                         )/ 1000                                  # kg of potassium excreted / pig
  

  # ------------------------------------------------------------------------
  
  # Cu, Zn and Se for spreading - without treatment - To calculate Ecotoxicity
  
  # performances$Cu_field= (performances$RejeteCu                   # mg / pig
  # )/ 1000 / 1000                          # kg of cooper excreted / pig
  # 
  # performances$Zn_field= (performances$RejeteZn                   # mg / pig
  # )/ 1000 / 1000                          # kg of zinc excreted / pig
  # 
  # performances$Se_int= (performances$IC * performances$GainPV 
  #                       * 30 * 0.05 / 1000)    # g of selenium intake
  # # 30 mg/kg of Se in Premix ingested to 5% (Garcia-Launay et al., 2014)
  # 
  # performances$Se_ret= (performances$GainPV) * 0.888 / 1000 
  # # g of selenium retained
  # 
  # performances$Se_field= (performances$Se_int - performances$Se_ret) / 1000 
  # # kg of selenium excreted
  
  # ------------------------------------------------------------------------
  
  # Potential of PO4 leaching (IPCC, 2006)
  performances$P_field= performances$RejeteP / 1000             # kg of phosphorus excreted / pig
  
  performances$Psubs= performances$P_field
  # kg fertilizer P substitution - nutrient balance (Sommer et al., 2008); 100% substitution
  
  performances$PO4_EmitH= performances$P_field * 0
  # kg PO4 emitted housing
  
  performances$PO4_leach= performances$P_field - performances$Psubs 
  # kg PO4 leaching
  
  # ------------------------------------------------------------------------
  
  # Manure management
  
  Amount_slurry=403             # kg slurry / fattening pig
  # Amount of manure produced. FR=403; BR=3.9(kg/pig/d); durée = 105 jours
  
  # ------------------------------------------------------------------------
  
  # LCA pre-Calculations
  
  # Feed production (in kg of feed / pig)
  
  performances$Electricity_feed= performances$CumulAlim * 41 / 1000    # kWh
  # 41 kWh of electricity/ton of feed (Garcia-Launay et al.(2014)
  
  performances$NaturalGas_feed= performances$CumulAlim *20.5 / 1000   # kWh
  # 20.5 kWh of natural gas/ton of feed (Garcia-Launay et al.(2014)
  
  performances$RoadTr_feed= 30 * performances$CumulAlim / 1000        # t.km
  # for 30 km of feed factory until farm (Garcia-Launay et al.(2014)
  # 30 t.km = (30 t sur 1 km) OR (1 t sur 30 km)
  
  
  # Housing
  
  Outdoor_housing=0     # m2 used in pig production in outdoor system    # m2 / pig
  
  Electricity_housing= 11.5   # 11.5 kWh of electricity / pig (Garcia-Launay et al.(2014))
  
  
  # Slurry
  
  Electricity_slurry=0     # kWh used for slurry treatment
  
  Transport_slurry= (Amount_slurry / 1000) * 10         # t.km
  # for 10 km of distance of farm to spreading
  
  Spreading_slurry= (Amount_slurry / 1000)              # t manure spread
  # (t) slurry to spreading
  
  
  # Other calculations
  
  performances$IC= ( performances$CumulAlim 
                     / (performances$PVFin - performances$PVInit) )  # kg de nourriture / kg gain de poids
  # IC de tous les porcs (m?me ceux qui sont morts avant d'?tre vendus ? l'abattoir)
  
  performances$GainPV=(performances$PVFin - performances$PVInit)        # kg pig (during all fattening)
  # IC de tous les porcs (m?me ceux qui sont morts avant d'?tre vendus ? l'abattoir)
  
  # ------------------------------------------------------------------------
  
  # LCA Calculations
  # Soja moyen du Bresil = 70% de la zone Centre Ouest (<1% en déforestation)
  # 30% de la zone Sud (0% en deforestation)
  
  # ------------------1.Climate Change------------------------
  
  ges_alimA = 0.6958          # kg CO2 / kg aliment
  ges_alimB = 0.4871          # kg CO2 / kg aliment
  
  Aliment_cout$ges <- ifelse (Aliment_cout$NomAlim=="AlimentA",
                              Aliment_cout$PctAli/100*Aliment_cout$QuantiteRation*ges_alimA,
                              Aliment_cout$PctAli/100*Aliment_cout$QuantiteRation*ges_alimB) 
  # kg CO2-eq/porc/aliment/jour
  
  vecteur2<-data.frame(aggregate(Aliment_cout$ges ~ Aliment_cout$Porc, FUN=sum))
  colnames(vecteur2)<-c("Porc","GESAlim") # kg CO2-eq/porc
  
  if (any(performances$AgeFin==performances$AgeInit)) {
    matrice <- data.frame(performances$Porc[performances$AgeFin==performances$AgeInit], 0)
    colnames(matrice)<-c("Porc","GESAlim")
    vecteur2 <- rbind(vecteur2, matrice)
  }
  
  vecteur2 <- vecteur2[order(as.character(vecteur2$Porc)),]
  performances <- performances[order(performances$Porc),] 
  row.names(performances) <- 1:nrow(performances)
  
  performances$GESAlim <- ifelse (performances$Porc==vecteur2$Porc, 
                                  vecteur2$GESAlim, 
                                  NA)  # kg CO2-eq/porc
                                  
  
  # LCA Climate Change (kg CO2-eq / pig) TOTAL
  
  performances$CC= (performances$GESAlim                                     # kg CO2-eq / pig
                    + performances$CH4Emit * 25                              # kg CH4 / pig             -> kg CO2-eq / pig 
                    + performances$N_N2OEmit * 44 / 28 * 298                 # kg N_N2O / pig           -> kg CO2-eq / pig
                    + (performances$N_NH3Emit*0.01)* 44 / 28 * 298           # kg N_NH3 / pig           -> kg CO2-eq / pig
                    # NH3 slurry (kg) indirect N2O emission (IPCC 2006)
                    # 1% du N-NH3 passe en N-N2O
                    + (performances$N_NOxEmit * 0.01)* 44/28 * 298           # kg N_NOx / pig           -> kg CO2-eq / pig
                    # NOx slurry (kg) indirect N2O emission  (IPCC 2006)
                    + (performances$N_NO3Emit * 0.0075)* 44/28 * 298         # kg N_NO3 / pig           -> kg CO2-eq / pig
                    # NO3 slurry (kg) indirect N2O emission  (IPCC 2006)
                    - performances$Nsubs * 6.576494                            # kg N / pig               -> kg CO2-eq / pig
                    # Substitution N mineral
                    - performances$Psubs * 4.053827                              # kg P / pig               -> kg CO2-eq / pig
                    # Substitution P mineral
                    - performances$K_field * 0.732281                          # kg K / pig               -> kg CO2-eq / pig
                    # Substitution K mineral  
                    + Electricity_housing * 0.00871855                # kWh / pig    -> kg CO2-eq / pig
                    + performances$Electricity_feed * 0.00871855      # kWh / pig    -> kg CO2-eq / pig
                    + Transport_slurry * 0.39964633                          # kg CO2 / t.km            -> kg CO2-eq / pig
                    # 10 km de transport
                    + Spreading_slurry * 0.01739462                          # kg CO2 / t manure spread -> kg CO2-eq / pig
                    # Débit d'épandage = 140 kg/s 
                    # Quantité de lisier à épandre = 403 kg
                    # Donc temps d'épandage  = 403/140 = 2.88 secondes
                    # Coeffient de conversion pour 1h d'épandage = 21.743272
                    # Donc Coeff pour 1h /3600 *2.88 = 0.01739462
                    + performances$NaturalGas_feed * 0.02060937              # kWh / kg feed            -> kg CO2-eq / pig
                    # natural gas for feed production 
                    + performances$RoadTr_feed * 0.13552743 )                # kgCO2 / t.km             -> kg CO2-eq / pig
                    # road transport for feed production
                 
  
  
  # Climate change ( kg CO2-eq / pig) FEED
  
  performances$CC_F= (performances$GESAlim                                            # kg CO2-eq / pig
                      + performances$NaturalGas_feed * 0.02060937   # kWh / kg feed  -> kg CO2-eq / pig
                      + performances$RoadTr_feed * 0.13552743       # kgCO2 / t.km   -> kg CO2-eq / pig
                      + performances$Electricity_feed * 0.00871855   # kWh / pig      -> kg CO2-eq / pig
                      )
  
  
  # Climate change ( kg CO2-eq / pig) HOUSING
  
  performances$CC_H= ( (performances$CH4_Enteric + performances$CH4_St*0.33)*25               # kg CO2-eq / pig
                       # 33% de temps de stockage en b?timent sous les animaux
                       + performances$N_N2OEmitH * 44 / 28 * 298             # kg N_N2O / pig -> kg CO2-eq / pig
                       + (performances$N_NH3EmitH * 0.01)* 44 / 28 * 298      # kg N_NH3 / pig -> kg CO2-eq / pig
                       + (performances$N_NOxEmitH * 0.01)* 44/28 * 298        # kg N_NOx / pig -> kg CO2-eq / pig
                       + (performances$N_NO3EmitH * 0.0075)* 44/28 * 298      # kg N_NO3 / pig -> kg CO2-eq / pig
                       + Electricity_housing * 0.00871855                     # kWh / pig      -> kg CO2-eq / pig
  )
  
  
  # Climate change ( kg CO2-eq / pig) MANURE
  
  performances$CC_M=performances$CC-performances$CC_F-performances$CC_H
  
  
  # ----------------------2.Acidification------------------------
  
  
  acid_alimA = 7.19         # g SO2-eq / kg aliment
  acid_alimB = 6.62         # g SO2-eq / kg aliment
  
  Aliment_cout$Acid <- ifelse (Aliment_cout$NomAlim=="AlimentA", 
                               Aliment_cout$PctAli/100*Aliment_cout$QuantiteRation*acid_alimA, 
                               Aliment_cout$PctAli/100*Aliment_cout$QuantiteRation*acid_alimB) 
  # g SO2-eq/porc/aliment/jour
  
  vecteur3<-data.frame(aggregate(Aliment_cout$Acid ~ Aliment_cout$Porc, FUN=sum))
  colnames(vecteur3)<-c("Porc","AcidAlim") # kg SO2-eq/porc
  
  if (any(performances$AgeFin==performances$AgeInit)) {
    colnames(matrice)<-c("Porc","AcidAlim")
    vecteur3 <- rbind(vecteur3, matrice)}
  
  vecteur3 <- vecteur3[order(as.character(vecteur3$Porc)),]
  performances <- performances[order(performances$Porc),] 
  row.names(performances) <- 1:nrow(performances)
  
  performances$AcidAlim <- ifelse (performances$Porc==vecteur3$Porc, 
                                   vecteur3$AcidAlim/1000, 
                                   NA)
  # performances$AcidAlim in kg SO2 / pig
  
  
  # LCA Acidification (g SO2-eq / pig)
  
  performances$AC= ( (performances$AcidAlim                          # kg SO2-eq / pig
                      + performances$N_NH3Emit * 17 / 14 * 1.6       # kg N_NH3 / pig -> kg SO2-eq / pig
                      - performances$Nsubs * 0.019816                  # kg N / pig     -> kg SO2-eq / pig
                      - performances$Psubs * 0.083455                  # kg P / pig     -> kg SO2-eq / pig
                      - performances$K_field * 0.005344                # kg K / pig     -> kg SO2-eq / pig
                      + Transport_slurry * 0.0024778163                # t.km           -> kg SO2-eq / pig
                      # 10 km de transport
                      + Spreading_slurry * 0.0001179684                # t              -> kg SO2-eq / pig
                      # D?bit d'épandage = 140 kg/s 
                      # Quantité de lisier à épandre = 403 kg
                      # Donc temps d'épandage  = 403/140 = 2.88 secondes
                      # Coeffient de conversion pour 1h d'épandage = 0.14746051
                      # Donc Coeff pour 1h /3600 *2.88 = 0.0001179684
                      + performances$RoadTr_feed * 0.00071396        # t.km           -> kg SO2-eq / pig
                      + Electricity_housing * 0.00004142              # kWh / pig      -> kg SO2-eq / pig
                      + performances$Electricity_feed * 0.00004142    # kWh / pig      -> kg SO2-eq / pig
                      + performances$NaturalGas_feed * 0.00007305    # kWh / kg feed  -> kg CO2-eq / pig
                      ) * 1000 )                                     # Convert kg to g SO2-eq

  
  # Acidification of feed( g SO2-eq / pig) FEED
  performances$AC_F= ( (performances$AcidAlim                          # kg SO2-eq / kg pig
                        + performances$RoadTr_feed * 0.00071396        # t.km -> kg SO2-eq / pig
                        + performances$NaturalGas_feed * 0.00007305    # kWh / kg feed  -> kg CO2-eq / pig
                        + performances$Electricity_feed * 0.00004142
                        ) * 1000 )                                     # Convert kg to g SO2-eq

  
  # Acidification of feed( g SO2-eq / pig) HOUSING
  performances$AC_H= ( (performances$N_NH3EmitH * 17 / 14 * 1.6       # kg SO2-eq / pig
                        + Electricity_housing * 0.00004142             # kWh / pig      -> kg SO2-eq / pig
                        )*1000 )     

  
  # Acidification of feed( g SO2-eq / pig) MANURE                   
  performances$AC_M= performances$AC - performances$AC_F - performances$AC_H            
  
  
  # -------------------3.Eutrophisation------------------------        
  
  
  EUTR_alimA = 3.827         # g PO4 / kg aliment
  EUTR_alimB = 3.345         # g PO4 / kg aliment
  
  Aliment_cout$EUTR <- ifelse (Aliment_cout$NomAlim=="AlimentA", 
                               Aliment_cout$PctAli/100*Aliment_cout$QuantiteRation*EUTR_alimA, 
                               Aliment_cout$PctAli/100*Aliment_cout$QuantiteRation*EUTR_alimB) 
  # g PO4-eq/porc/aliment/jour
  
  vecteur4<-data.frame(aggregate(Aliment_cout$EUTR ~ Aliment_cout$Porc, FUN=sum))
  colnames(vecteur4)<-c("Porc","EUTRAlim")  # g PO4-eq/porc
  
  if (any(performances$AgeFin==performances$AgeInit)) {
    colnames(matrice)<-c("Porc","EUTRAlim")
    vecteur4 <- rbind(vecteur4, matrice)
  }
  
  vecteur4 <- vecteur4[order(as.character(vecteur4$Porc)),]
  performances <- performances[order(performances$Porc),] 
  row.names(performances) <- 1:nrow(performances)
  
  performances$EUTRAlim <- ifelse (performances$Porc==vecteur4$Porc, 
                                   vecteur4$EUTRAlim/1000, 
                                   NA)
  # kg PO4-eq / pig
  
  # Eutrophication of diet( g PO4-eq / pig)  TOTAL
  
  performances$EU=( (performances$EUTRAlim                            # kg PO4-eq / pig
                     + performances$N_NH3Emit * 17 / 14 * 0.35        # kg N_NH3/pig -> kg PO4-eq / pig
                     + performances$N_NO3Emit * 62 / 14 * 0.1         # kg N_NO3/pig -> kg PO4-eq / pig
                     - performances$Nsubs * 0.008919                    # kg N/pig     -> kg PO4-eq / pig
                     - performances$Psubs * 0.096651                    # kg P/pig     -> kg PO4-eq / pig
                     - performances$K_field * 0.001153                  # kg K/pig     -> kg PO4-eq / pig
                     + Transport_slurry * 0.00066817162                  # t.km         -> kg PO4-eq / pig
                     # 10 km de transport
                     + Spreading_slurry * 0.00003017293                # t            -> kg PO4-eq / pig
                     # D?bit d'épandage = 140 kg/s 
                     # Quantité de lisier à épandre = 403 kg
                     # Donc temps d'épandage  = 403/140 = 2.88 secondes
                     # Coeffient de conversion pour 1h d'épandage = 0.037716159
                     # Donc Coeff pour 1h /3600 *2.88 = 0.00003017293
                     + performances$RoadTr_feed * 0.00013245                # t.km         -> kg PO4-eq / pig
                     + Electricity_housing * 0.00000341              # kWh / pig    -> kg SO2-eq / pig
                     + performances$Electricity_feed * 0.00000341    # kWh / pig    -> kg SO2-eq / pig
                     + performances$NaturalGas_feed * 0.00000215    # kWh / kg feed  -> kg CO2-eq / pig
                     ) * 1000)                                        # Convert kg to g PO4-eq
  
  
  # Eutrophication of diet( g PO4-eq / pig)  FEED
  performances$EU_F=( (performances$EUTRAlim                        # kg PO4-eq / pig
                       + performances$RoadTr_feed * 0.00013245      # t.km         -> kg PO4-eq / pig
                       + performances$Electricity_feed * 0.00000341               # kWh / pig    -> kg SO2-eq / pig
                       + performances$NaturalGas_feed * 0.00000215  # kWh / kg feed  -> kg CO2-eq / pig
                       )* 1000)                                          # Convert kg to g PO4-eq
  
  
  # Eutrophication of diet( g PO4-eq / pig)  HOUSING
  performances$EU_H=( (performances$N_NH3EmitH * 17 / 14 * 0.35       # kg N_NH3/pig -> kg PO4-eq / pig
                       + performances$N_NO3EmitH * 62 / 14 * 0.1      # kg N_NO3/pig -> kg PO4-eq / pig
                       + Electricity_housing * 0.00000341              # kWh / pig    -> kg SO2-eq / pig 
                       ) * 1000)              # Convert kg to g PO4-eq
  
  
  # Eutrophication of diet( g PO4-eq / pig)  MANURE
  performances$EU_M= performances$EU - performances$EU_F - performances$EU_H
  
  
  # -----------------4.Cumulative Energy Demand------------------
  
  
  Energie_alimA = 6.8472        # MJ / kg aliment
  Energie_alimB = 5.049         # MJ / kg aliment
  
  Aliment_cout$Energie <- ifelse (Aliment_cout$NomAlim=="AlimentA", 
                                  Aliment_cout$PctAli/100*Aliment_cout$QuantiteRation*Energie_alimA, 
                                  Aliment_cout$PctAli/100*Aliment_cout$QuantiteRation*Energie_alimB) 
  # MJ/porc/aliment/jour
  
  vecteur5<-data.frame(aggregate(Aliment_cout$Energie ~ Aliment_cout$Porc, FUN=sum))
  colnames(vecteur5)<-c("Porc","EnergieAlim") # MJ/porc
  
  if (any(performances$AgeFin==performances$AgeInit)) {
    colnames(matrice)<-c("Porc","EnergieAlim")
    vecteur5 <- rbind(vecteur5, matrice)}
  
  vecteur5 <- vecteur5[order(as.character(vecteur5$Porc)),]
  performances <- performances[order(performances$Porc),] 
  row.names(performances) <- 1:nrow(performances)
  
  performances$EnergieAlim <- ifelse (performances$Porc==vecteur5$Porc, 
                                      vecteur5$EnergieAlim, 
                                      NA)
  # MJ / pig
  
  
  # LCA Cumulative Energy Demand (MJ / pig)
  
  performances$CED= (performances$EnergieAlim                                          # MJ / pig
                     - performances$Nsubs * 59.858776                # kg N / pig        -> MJ / pig
                     - performances$Psubs * 62.806341                # kg N / pig        -> MJ / pig
                     - performances$K_field * 12.6701042               # kg N / pig        -> MJ / pig
                     + Electricity_housing * 0.90513108                 # kWh / pig         -> MJ / pig
                     + performances$Electricity_feed * 0.90513108       # kWh / pig         -> MJ / pig
                     + Transport_slurry * 5.2284446               # t.km              -> MJ / pig
                     # 10 km de transport
                     + Spreading_slurry * 0.2907724               # t                 -> MJ / pig
                     # Débit d'épandage = 140 kg/s 
                     # Quantité de lisier à épandre = 403 kg
                     # Donc temps d'épandage  = 403/140 = 2.88 secondes
                     # Coeffient de conversion pour 1h d'épandage = 363.46499
                     # Donc Coeff pour 1h /3600 *2.88 = 0.290772
                     + performances$NaturalGas_feed * 0.33626306   # kWh / ton of feed -> MJ / pig
                     + performances$RoadTr_feed * 2.22760640)      # t.km              -> MJ / pig
  
  
  # cumulative energy demand of diet( MJ/ pig)   FEED
  performances$CED_F= (performances$EnergieAlim                                          # MJ / pig
                       + performances$NaturalGas_feed * 0.33626306      # kWh / ton of feed -> MJ / pig
                       + performances$RoadTr_feed * 2.22760640          # t.km              -> MJ / pig
                       + performances$Electricity_feed * 0.90513108                   # kWh / pig         -> MJ / pig
                       )         
  
  
  # cumulative energy demand of diet( MJ/ pig) HOUSING
  performances$CED_H= (Electricity_housing * 0.90513108)                 # kWh / pig         -> MJ / pig
  
  
  # cumulative energy demand of diet( MJ/ pig) MANURE
  performances$CED_M= performances$CED - performances$CED_F - performances$CED_H 


# ------------------5.Land occupation----------------------


  Surface_alimA = 1.3442           # m2 year / kg aliment
  Surface_alimB = 1.1679           # m2 year / kg aliment
  
  Aliment_cout$Surface <- ifelse (Aliment_cout$NomAlim=="AlimentA", 
                                  Aliment_cout$PctAli/100*Aliment_cout$QuantiteRation*Surface_alimA, 
                                  Aliment_cout$PctAli/100*Aliment_cout$QuantiteRation*Surface_alimB) 
  # m2 year/porc/aliment/jour
  
  vecteur6<-data.frame(aggregate(Aliment_cout$Surface ~ Aliment_cout$Porc, FUN=sum))
  colnames(vecteur6)<-c("Porc","SurfaceAlim") # m2 year/porc
  
  if (any(performances$AgeFin==performances$AgeInit)) {
    colnames(matrice)<-c("Porc","SurfaceAlim")
    vecteur6 <- rbind(vecteur6, matrice)
  }
  
  vecteur6 <- vecteur6[order(as.character(vecteur6$Porc)),]
  performances <- performances[order(performances$Porc),] 
  row.names(performances) <- 1:nrow(performances)
  
  performances$SurfaceAlim <- ifelse (performances$Porc==vecteur6$Porc, 
                                      vecteur6$SurfaceAlim, 
                                      NA)
  # m2 year / pig
  
  
  # Land occupation of diet(m2 year / pig)  TOTAL
  performances$Land= (performances$SurfaceAlim                           # m2 year/ pig
                      + Electricity_housing * 0.00045690               # kWh / pig -> m2 year/ pig
                      + performances$Electricity_feed * 0.00045690     # kWh / pig -> m2 year/ pig
                      + Outdoor_housing * 1                              # m2 year/ pig  
                      + performances$RoadTr_feed * 0.01129125            # t.km -> m2 year/ pig
                      + performances$NaturalGas_feed * 0.00006732        # kWh / kg feed -> kg CO2-eq / pig
                      + Transport_slurry * 0.074090177                    # t.km    -> m2 year / pig
                      # 10 km de transport
                      + Spreading_slurry * 0.0002984623 )                 # t / pig
                      # Débit d'épandage = 140 kg/s 
                      # Quantité de lisier à épandre = 403 kg
                      # Donc temps d'épandage  = 403/140 = 2.88 secondes
                      # Coeffient de conversion pour 1h d'épandage = 0.37307785
                      # Donc Coeff pour 1h /3600 *2.88 = 0.0002984623
                      
  
  # Land occupation of diet(m2 year / pig) FEED  
  performances$Land_F= (performances$SurfaceAlim                         # m2 year/ pig             
                        + performances$RoadTr_feed * 0.01129125          # t.km      -> m2 year/ pig
                        + performances$Electricity_feed * 0.00045690                   # kWh / pig -> m2 year/ pig
                        + performances$NaturalGas_feed * 0.00006732      # kWh / kg feed -> kg CO2-eq / pig
                        )                     
  
  # Land occupation of diet(m2 year / pig) HOUSING
  performances$Land_H= (Electricity_housing * 0.00045690                  # kWh / pig -> m2 year/ pig
                        + Outdoor_housing * 1)                           # m2 year/ pig
  
  
  # Land occupation of diet(m2 year / pig) MANURE
  performances$Land_M= performances$Land - performances$Land_F - performances$Land_H
  


#----------------------------------------------------------------------------
  #SOW
#----------------------------------------------------------------------------
  
  # A l'echelle de la truie
  
  RejeteN_sow = 21479       # g N / pig
  RejeteK_sow = 7481        # g K / pig
  RejeteP_sow = 4877        # g P / pig
  MOex_sow = 364.93         # kg of OM excreted during a year
  CumulAlim_sow = 1207      # kg feed / animal / year
  Amount_slurry_sow = 4600   # kg slurry / animal / year 
  
  # A l'echelle de l'aliment
  
  Feed_gestation_sow = 0.78   # %
  Feed_lactation_sow = 0.22   # %
  MAT_AlimGest = 140      # g / kg MS
  MAT_AlimLact = 165      # g / kg MS
  ResDt_AlimGest = 101.3    # g / kg feed
  ResDt_AlimLact = 108.8    # g / kg feed
  GES_AlimGest = 0.39     # kg CO2-eq / kg gestation feed
  GES_AlimLact = 0.59     # kg CO2-eq / kg lactation feed
  Acid_AlimGest = 4.73    # g SO2-eq / kg feed
  Acid_AlimLact = 6.57    # g SO2-eq / kg feed
  EUTR_AlimGest = 3.22       # g PO4-eq / kg feed
  EUTR_AlimLact = 3.64       # g PO4-eq / kg feed
  Energie_AlimGest = 3.49      # MJ / kg feed
  Energie_AlimLact = 5.35      # MJ / kg feed
  Surface_AlimGest = 1.34      # m2 year/ kg feed
  Surface_AlimLact = 1.31      # m2 year/ kg feed
  Cout_AlimGest = 206.35       # euros / feed ton 
  Cout_AlimLact = 247.53       # euros / tonne d'aliment pour un porcelet
  
  
  # param_feed_sow <- matrix(list(MAT_AlimGest, MAT_AlimLact,
  #                          ResDt_AlimGest, ResDt_AlimLact,
  #                          GES_AlimGest, GES_AlimLact,
  #                          Acid_AlimGest, Acid_AlimLact,
  #                          EUTR_AlimGest, EUTR_AlimLact,
  #                          Energie_AlimGest, Energie_AlimLact,
  #                          Surface_AlimGest, Surface_AlimLact,
  #                          Cout_AlimGest, Cout_AlimLact,
  #                          Feed_gestation_sow, Feed_lactation_sow),
  #                          nrow=2, ncol=9, byrow=FALSE)
  # colnames(param_feed_sow)<-c("MAT","ResDt","GES","Acid",
  #                             "EUTR","Energie","Surface",
  #                             "CoutAlim","Proportion")
  # rownames(param_feed_sow)<-c("Gestation","Lactation")
  
  
  # Calculations
  
  CoutAlim_sow = round((Cout_AlimGest*Feed_gestation_sow
                   + Cout_AlimLact*Feed_lactation_sow), digits=2)
  # Co?t moyen de l'aliment truie en euros par tonne d'aliments
  
  
# N COMPOUNDS ------------------------
  
  
  # Housing Emissions (kg manure N ex-animal *emission factor) 
  
  N_NH3EmitH_sow = (RejeteN_sow / 1000) * 0.24          # kg N-NH3 / pig
  # Rigolot et al.(2010b): EF is 0.24 kg N-NH3/kg N excreted
  
  N_N2OEmitH_sow = (RejeteN_sow / 1000) * 0.002         # kg N-N2O / pig
  # IPCC (2006): TABLE 10.21-Pit sotorage bellow animal confinements
  
  N_NOxEmitH_sow = (RejeteN_sow / 1000) * 0.002         # kg N-NOx / pig
  # Dämmgen & Hutchings(2008)cited by IPCC(2006)
  
  N_NO3EmitH_sow =(RejeteN_sow / 1000) * 0              # kg N-NO3 / pig
  
  # N leaving after housing emissions (kg N / pig)
  N_St_sow = (RejeteN_sow/1000 - N_NH3EmitH_sow 
              - N_N2OEmitH_sow - N_NOxEmitH_sow)
  
  
  # Storage Emissions (kg manure N ex-housing *emission factor) 
  
  N_NH3EmitSt_sow = N_St_sow * 0.1                      # kg N-NH3 / pig
  # Rigolot et al.(2010b): Table 4
  
  N_N2OEmitSt_sow = N_St_sow * 0.001                    # kg N-N2O / pig
  # IPCC (2006)
  
  N_NOxEmitSt_sow = N_St_sow * 0.005                    # kg N-NOx / pig
  # Dämmgen & Hutchings(2008)cited by IPCC(2006)
  
  # N leaving after storage emissions (kg N / pig)
  N_field_sow = (N_St_sow - N_NH3EmitSt_sow
                 - N_N2OEmitSt_sow - N_NOxEmitSt_sow)
  
  
  # Emissions during field application (kg manure N ex-storage *emission factor) 
  
  N_NH3EmitFA_sow = N_field_sow * 0.07                  # kg N-NH3 / pig
  # Andersen et al. (2001)
  
  N_N2OEmitFA_sow = N_field_sow * 0.01                  # kg N-N2O / pig
  # IPCC (2006)
  
  N_NOxEmitFA_sow = N_field_sow * 0.001                 # kg N-NOx / pig
  # Nemecek and Kägi (2007)
  
  # Content of manure after field application (kg N / pig)
  N_AfterFA_sow = (N_field_sow - N_NH3EmitFA_sow
                   - N_N2OEmitFA_sow - N_NOxEmitFA_sow)
  
  
  # Emissions (No3 leaching) after field application
  N_NO3EmitFA_sow = N_AfterFA_sow * 0.05                # kg N / pig
  # Garcia-Launay et al.(2014)
  
  # Extraloss as nitrates compared to mineral fertilizer
  NO3_leach_increase_sow = 0.05 * N_field_sow          # kg N / pig
  # (IPCC, 2006 ; eventuellement revoir/ publi Danemark)
  
  
  # Total emissions of N compounds (sum of all emissions, kg)
  N_NH3Emit_sow = ( N_NH3EmitH_sow
                    + N_NH3EmitSt_sow 
                    + N_NH3EmitFA_sow )   
  N_N2OEmit_sow = ( N_N2OEmitH_sow 
                    + N_N2OEmitSt_sow 
                    + N_N2OEmitFA_sow )
  N_NOxEmit_sow = ( N_NOxEmitH_sow 
                    + N_NOxEmitSt_sow 
                    + N_NOxEmitFA_sow )
  N_NO3Emit_sow = ( N_NO3EmitFA_sow 
                    + NO3_leach_increase_sow )
  
  
  # Total (kg N / pig)
  N_Emit_sow = ( N_NH3Emit_sow 
                + N_N2OEmit_sow 
                + N_NOxEmit_sow 
                + N_NO3Emit_sow )
  
  
  # Fertilizer N substitution - nutrient balance
  Nsubs_sow = 0.75 * N_field_sow                       # kg N / pig
  # (Nguyen et al., 2010)
  

  
# CH4 EMISSIONS ------------------------
  
  ResDt_sow = ( ((ResDt_AlimGest * Feed_gestation_sow 
                + ResDt_AlimLact * Feed_lactation_sow)/1000)
                * CumulAlim_sow )
  # kg / sow

  # CH4 emissions by enteric fermentation (Digested fibre intake)

  CH4_Enteric_sow = (ResDt_sow * 1000                       # g / piglet  
                     * 1340                                 # 1340 J / g
                     / 1000000) / 55.65                     # MJ / pig
  # enteric. Rigolot et al.(2010a)
  
  
  # CH4 emissions at storage 
  
  # performances$MOexc : kg of OM excreted
  VS_sow = MOex_sow                    # kg, VS equals volatile solids cf equation IPCC 2006
  # Rigolot et al.(2010b)
  
  B0=0.42
  # IPCC (2006) gives a value of 0.45 m3 CH4/kg VS, for Western Europe (and and 0.29 for Latin América)
  
  TempMoy_sow = 2/3*12+1/3*12
  # 2/3 indoor slurry storage at 20?C and 1/3 outdoor st. at 12?C  (Garcia-Launay et al., 2014)
  
  MCF_sow=(0.1517*(TempMoy_sow^2)-2.1569*TempMoy_sow+24.634)/100    
  # Garcia-Launay et al. (2014): MCFlisier=0.1517*(average temp^2)-2.1569*(average temp)+24.634
  # MCF : methane conversion factor
  
  CH4_St_sow = VS_sow*(B0*0.67)*MCF_sow                              
  # kg of CH4 storage / sow
  # Rigolot et al. (2010b) from IPCC(2006)
  # 0.67 for convert m3 to kg
  
  # Total of CH4 emissions
  CH4Emit_sow = ( CH4_Enteric_sow 
                  + CH4_St_sow ) 
  # kg of CH4 emitted by enteric fermentation and manure storage

  
# K EMISSIONS ------------------------

  K_field_sow = (RejeteK_sow               # g / pig
                 )/ 1000                   # kg of potassium excreted / pig
  
  
# PO4 LEACHING ------------------------

  P_field_sow = RejeteP_sow/1000        # kg P / pig
  Psubs_sow = P_field_sow               # kg P / pig (100% substitution P mineral)
  
  
# LCA pre-Calculations ----------------
  
  # Feed production (kg of feed / pig)
  
  Electricity_feed_sow = CumulAlim_sow * 41 / 1000         # kWh
  NaturalGas_feed_sow = CumulAlim_sow * 20.5 / 1000        # kWh
  RoadTr_feed_sow = 30 * CumulAlim_sow / 1000              # t.km
  
  # Housing
  
  Outdoor_housing_sow = 0     # m2 used in pig production in outdoor system
  Electricity_housing_sow = 11.5   # kWh
  
  # Slurry
  
  Transport_slurry_sow = (Amount_slurry_sow / 1000) * 10         # t.km
  Spreading_slurry_sow = (Amount_slurry_sow / 1000)              # t manure spread


  
# ------------------1.Climate Change------------------------
  
  GESAlim_sow = ( (GES_AlimGest * Feed_gestation_sow 
                  + GES_AlimLact * Feed_lactation_sow)
                  *CumulAlim_sow )
  # kg CO2-eq / sow
  
# TOTAL (kg CO2-eq / pig)
  CC_sow= (GESAlim_sow                                     # kg CO2-eq / pig
           + CH4Emit_sow * 25                              # kg CH4 / pig             -> kg CO2-eq / pig 
           + N_N2OEmit_sow * 44 / 28 * 298                 # kg N_N2O / pig           -> kg CO2-eq / pig
           + (N_NH3Emit_sow * 0.01)* 44 / 28 * 298         # kg N_NH3 / pig           -> kg CO2-eq / pig
           # NH3 slurry (kg) indirect N2O emission (IPCC 2006)
           # 1% du N-NH3 passe en N-N2O
           + (N_NOxEmit_sow * 0.01)* 44/28 * 298           # kg N_NOx / pig           -> kg CO2-eq / pig
           # NOx slurry (kg) indirect N2O emission  (IPCC 2006)
           + (N_NO3Emit_sow * 0.0075)* 44/28 * 298         # kg N_NO3 / pig           -> kg CO2-eq / pig
           # NO3 slurry (kg) indirect N2O emission  (IPCC 2006)
           - Nsubs_sow * 6.576494                            # kg N / pig               -> kg CO2-eq / pig
           # Substitution N mineral
           - Psubs_sow * 4.053827                              # kg P / pig               -> kg CO2-eq / pig
           # Substitution P mineral
           - K_field_sow * 0.732281                          # kg K / pig               -> kg CO2-eq / pig
           # Substitution K mineral 
           + (Electricity_housing_sow+Electricity_feed_sow)* 0.00871855  # kWh / pig -> kg CO2-eq / pig
           + Transport_slurry_sow * 0.39964633                 # kg CO2 / t.km            -> kg CO2-eq / pig
           # 10 km de transport
           + Spreading_slurry_sow * 0.1984678                 # kg CO2 / t manure spread -> kg CO2-eq / pig
           # Débit d'épandage = 140 kg/s 
           # Quantité de lisier à épandre = 4600 kg par truie par an
           # Donc temps d'épandage  = 4600/140 = 32.86 secondes
           # Coeffient de conversion pour 1h d'épandage = 21.743272
           # Donc Coeff pour 1h /3600 *32.86 = 0.1984678
           + NaturalGas_feed_sow * 0.02060937              # kWh / kg feed            -> kg CO2-eq / pig
           # natural gas for feed production 
           + RoadTr_feed_sow * 0.13552743 )                # kgCO2 / t.km             -> kg CO2-eq / pig
           # road transport for feed production
  
  
# FEED (kg CO2-eq / pig)
  CC_F_sow = (GESAlim_sow                                  # kg CO2-eq / pig
              + NaturalGas_feed_sow * 0.02060937           # kWh / kg feed  -> kg CO2-eq / pig
              + RoadTr_feed_sow * 0.13552743               # kgCO2 / t.km   -> kg CO2-eq / pig
              + Electricity_feed_sow * 0.00871855 )         # kWh / pig      -> kg CO2-eq / pig


# HOUSING (kg CO2-eq / pig)
  CC_H_sow = ( (CH4_Enteric_sow + CH4_St_sow*0.33)*25               # kg CO2-eq / pig
                 # 33% de temps de stockage en batiment sous les animaux
                + N_N2OEmitH_sow * 44 / 28 * 298                     # kg N_N2O / pig -> kg CO2-eq / pig
                + (N_NH3EmitH_sow * 0.01)* 44/28 * 298               # kg N_NH3 / pig -> kg CO2-eq / pig
                + (N_NOxEmitH_sow * 0.01)* 44/28 * 298               # kg N_NOx / pig -> kg CO2-eq / pig
                + (N_NO3EmitH_sow * 0.0075)* 44/28 * 298             # kg N_NO3 / pig -> kg CO2-eq / pig
                + Electricity_housing * 0.00871855 )                 # kWh / pig      -> kg CO2-eq / pig

  
# MANURE (kg CO2-eq / pig)
  CC_M_sow = CC_sow - CC_F_sow - CC_H_sow
  
  
# ----------------------2.Acidification------------------------

  AcidAlim_sow = ( (Acid_AlimGest * Feed_gestation_sow 
                    + Acid_AlimLact * Feed_lactation_sow)
                   * CumulAlim_sow
                   / 1000 )                # Convert g to kg SO2-eq
  
# TOTAL (g SO2-eq / pig)
  AC_sow = ( (AcidAlim_sow                                        # kg SO2-eq / pig
              + N_NH3Emit_sow * 17 / 14 * 1.6                     # kg N_NH3 / pig -> kg SO2-eq / pig
              - Nsubs_sow * 0.019816                                # kg N / pig     -> kg SO2-eq / pig
              - Psubs_sow * 0.083455                                # kg P / pig     -> kg SO2-eq / pig
              - K_field_sow * 0.005344                              # kg K / pig     -> kg SO2-eq / pig
              + Transport_slurry_sow * 0.0024778163                     # t.km           -> kg SO2-eq / pig
              # 10 km de transport
              + Spreading_slurry_sow * 0.001345987                    # t              -> kg SO2-eq / pig
              # Débit d'épandage = 140 kg/s 
              # Quantité de lisier à épandre = 4600 kg par truie par an
              # Donc temps d'épandage  = 4600/140 = 32.86 secondes
              # Coeffient de conversion pour 1h d'épandage = 0.14746051
              # Donc Coeff pour 1h /3600 *32.86 = 0.001345987
              + RoadTr_feed_sow * 0.00071396                                 # t.km           -> kg SO2-eq / pig
              + (Electricity_housing_sow+Electricity_feed_sow) * 0.00004142    # kWh / sow      -> kg SO2-eq / pig
              + NaturalGas_feed_sow * 0.00007305                         # kWh / kg feed  -> kg CO2-eq / pig
              ) * 1000 )                                     # Convert kg to g SO2-eq
  
  
# FEED (g SO2-eq / pig)
  AC_F_sow = ( (AcidAlim_sow                                      # kg SO2-eq / kg pig
                + RoadTr_feed_sow * 0.00071396                    # t.km -> kg SO2-eq / pig
                + NaturalGas_feed_sow * 0.00007305                # kWh / kg feed  -> kg CO2-eq / pig
                + Electricity_feed_sow * 0.00004142
                ) * 1000 )                                        # Convert kg to g SO2-eq

  
# HOUSING (g SO2-eq / pig)
  AC_H_sow = ( (N_NH3EmitH_sow * 17 / 14 * 1.6                   # kg SO2-eq / pig
                + Electricity_housing_sow * 0.00004142            # kWh / pig  -> kg SO2-eq / pig
                )*1000 ) 
  

# MANURE (g SO2-eq / pig)                   
  AC_M_sow = AC_sow - AC_F_sow - AC_H_sow 
  
  
# -------------------3.Eutrophisation------------------------
  
  EUTRAlim_sow = ( (EUTR_AlimGest * Feed_gestation_sow
                  + EUTR_AlimLact * Feed_lactation_sow)
                  * CumulAlim_sow 
                  /1000 )              # Convert g to kg PO4-eq
  
  
  # TOTAL (g PO4-eq / pig)
  EU_sow =( (EUTRAlim_sow                            # kg PO4-eq / pig
             + N_NH3Emit_sow * 17 / 14 * 0.35        # kg N_NH3/pig -> kg PO4-eq / pig
             + N_NO3Emit_sow * 62 / 14 * 0.1         # kg N_NO3/pig -> kg PO4-eq / pig
             - Nsubs_sow * 0.008919                    # kg N/pig     -> kg PO4-eq / pig
             - Psubs_sow * 0.096651                    # kg P/pig     -> kg PO4-eq / pig
             - K_field_sow * 0.001153                  # kg K/pig     -> kg PO4-eq / pig
             + Transport_slurry_sow * 0.00066817162         # t.km         -> kg PO4-eq / pig
             # 10 km de transport ? une vitesse de 10km/h donc 1h de transport
             + Spreading_slurry_sow * 0.0003442647        # t            -> kg PO4-eq / pig             + RoadTr_feed_sow * 0.0002              # t.km         -> kg PO4-eq / pig
             # Débit d'épandage = 140 kg/s 
             # Quantité de lisier à épandre = 4600 kg par truie par an
             # Donc temps d'épandage  = 4600/140 = 32.86 secondes
             # Coeffient de conversion pour 1h d'épandage = 0.037716159
             # Donc Coeff pour 1h /3600 *32.86 = 0.0003442647
             + RoadTr_feed_sow * 0.00013245        # t.km         -> kg PO4-eq / pig
             + (Electricity_housing_sow+Electricity_feed_sow) * 0.00000341    # kWh / pig    -> kg SO2-eq / pig
             + NaturalGas_feed_sow * 0.00000215    # kWh / kg feed  -> kg CO2-eq / pig
             ) * 1000)                             # Convert kg to g PO4-eq
  

  # FEED (g PO4-eq / pig)
  EU_F_sow = ( (EUTRAlim_sow                            # kg PO4-eq / pig
                + RoadTr_feed_sow * 0.00013245          # t.km         -> kg PO4-eq / pig
                + Electricity_feed_sow * 0.00000341      # kWh / pig    -> kg SO2-eq / pig
                + NaturalGas_feed_sow * 0.00000215      # kWh / kg feed  -> kg CO2-eq / pig
                )* 1000)                                # Convert kg to g PO4-eq
  

  # HOUSING (g PO4-eq / pig)
  EU_H_sow = ( (N_NH3EmitH_sow * 17 / 14 * 0.35        # kg N_NH3/pig -> kg PO4-eq / pig
                + N_NO3EmitH_sow * 62 / 14 * 0.1       # kg N_NO3/pig -> kg PO4-eq / pig
                + Electricity_housing * 0.00000341      # kWh / pig    -> kg SO2-eq / pig 
                ) * 1000)                              # Convert kg to g PO4-eq
  

  # MANURE (g PO4-eq / pig)
  EU_M_sow = EU_sow - EU_F_sow - EU_H_sow
  
  
# -----------------4.Cumulative Energy Demand------------------
  
  EnergieAlim_sow = ( (Energie_AlimGest * Feed_gestation_sow 
                     + Energie_AlimGest * Feed_lactation_sow)
                     *CumulAlim_sow )
  # MJ / pig
  
  # TOTAL (MJ/ pig)
  CED_sow = (EnergieAlim_sow                         # MJ / pig
             - Nsubs_sow * 59.858776                   # kg N / pig        -> MJ / pig
             - Psubs_sow * 62.806341                   # kg N / pig        -> MJ / pig
             - K_field_sow * 12.670042                  # kg N / pig        -> MJ / pig
             + (Electricity_housing_sow+Electricity_feed_sow) * 0.90513108      # kWh / pig         -> MJ / pig
             + Transport_slurry_sow * 5.2284446         # t.km              -> MJ / pig
             # 10 km de transport
             + Spreading_slurry_sow * 3.317628         # t                 -> MJ / pig
             # Débit d'épandage = 140 kg/s 
             # Quantité de lisier à épandre = 4600 kg par truie et par an
             # Donc temps d'épandage  = 4600/140 = 32.86 secondes
             # Coeffient de conversion pour 1h d'épandage = 363.46499
             # Donc Coeff pour 1h /3600 *32.86 = 3.317628
             + NaturalGas_feed_sow * 0.33626306       # kWh / ton of feed -> MJ / pig
             + RoadTr_feed_sow * 2.22760640 )         # t.km              -> MJ / pig
  

  # FEED (MJ/ pig)
  CED_F_sow = (EnergieAlim_sow                            # MJ / pig
               + NaturalGas_feed_sow * 0.33626306         # kWh / ton of feed -> MJ / pig
               + RoadTr_feed_sow * 2.22760640             # t.km              -> MJ / pig
               + Electricity_feed_sow * 0.90513108 )       # kWh / pig         -> MJ / pig


  # HOUSING (MJ/ pig)
  CED_H_sow = (Electricity_housing_sow * 0.90513108)       # kWh / pig         -> MJ / pig

  
  # MANURE (MJ/ pig)
  CED_M_sow = CED_sow - CED_F_sow - CED_H_sow 
  

# ------------------5.Land occupation----------------------
  
  SurfaceAlim_sow = ( (Surface_AlimGest * Feed_gestation_sow 
                       + Surface_AlimLact * Feed_lactation_sow)
                       * CumulAlim_sow )
  # m2 year/ pig
  
  # TOTAL (m2 year / pig)
  Land_sow = (SurfaceAlim_sow                             # m2 year/ pig
              + (Electricity_housing_sow+Electricity_feed_sow)*0.00045690        # kWh / pig -> m2 year/ pig
              + Outdoor_housing_sow * 1                   # m2 year/ pig  
              + RoadTr_feed_sow * 0.01129125              # t.km -> m2 year/ pig
              + NaturalGas_feed_sow * 0.00006732          # kWh / kg feed -> kg CO2-eq / pig
              + Transport_slurry_sow * 0.074090177             # t.km              -> MJ / pig
              # 10 km de transport 
              + Spreading_slurry_sow * 0.003405372 )          # t      
              # Débit d'épandage = 140 kg/s 
              # Quantité de lisier à épandre = 4600 kg par truie par an
              # Donc temps d'épandage  = 403/140 = 32.86 secondes
              # Coeffient de conversion pour 1h d'épandage = 0.37307785
              # Donc Coeff pour 1h /3600 *32.86 = 0.003405372


  # FEED (m2 year / pig)  
  Land_F_sow = (SurfaceAlim_sow                           # m2 year/ pig             
                + RoadTr_feed_sow * 0.01129125            # t.km      -> m2 year/ pig
                + Electricity_feed_sow * 0.00045690             # kWh / pig -> m2 year/ pig
                + NaturalGas_feed_sow * 0.00006732 )      # kWh / kg feed -> kg CO2-eq / pig


  # HOUSING (m2 year / pig)
  Land_H_sow = (Electricity_housing_sow * 0.00045690       # kWh / pig -> m2 year/ pig
                + Outdoor_housing_sow * 1)                # m2 year/ pig
  
  
  # MANURE (m2 year / pig)
  Land_M_sow = Land_sow - Land_F_sow - Land_H_sow
  
  
#----------------------------------------------------------------------------
  # POST WEANING
#----------------------------------------------------------------------------
  
  # A l'echelle du porcelet
  
  RejeteN_piglet = 607      # g N / pig
  RejeteK_piglet = 312      # g K / piglet
  RejeteP_piglet = 139      # g P / piglet
  MOex_piglet = 12.09       # kg of OM excreted /piglet
  CumulAlim_piglet = 41.2    # kg feed / animal
  Amount_slurry_piglet = 46  # kg slurry / animal
  
  # A l'echelle de l'aliment
  
  Feed_phase1_piglet = 0.33   # %
  Feed_phase2_piglet = 0.67   # %
  MAT_AlimPS1 = 200      # g / kg MS
  MAT_AlimPS2 = 180      # g / kg MS
  ResDc_AlimPS1 = 195.7    # g / kg feed
  ResDc_AlimPS2 = 112.7    # g / kg feed
  GES_AlimPS1 = 1.03     # kg CO2-eq / kg phase1 feed
  GES_AlimPS2 = 0.77     # kg CO2-eq / kg phase2 feed
  Acid_AlimPS1 = 7.5     # g SO2-eq / kg feed
  Acid_AlimPS2 = 7.1     # g SO2-eq / kg feed
  EUTR_AlimPS1 = 4.1          # g PO4-eq / kg feed
  EUTR_AlimPS2 = 4.2          # g PO4-eq / kg feed
  Energie_AlimPS1 = 6.97      # MJ / kg feed
  Energie_AlimPS2 = 6.75      # MJ / kg feed
  Surface_AlimPS1 = 1.42      # m2 year/ kg feed
  Surface_AlimPS2 = 1.44      # m2 year/ kg feed
  Cout_AlimPS1 = 434.06       # euros / feed ton 
  Cout_AlimPS2 = 290.33       # euros / tonne d'aliment pour un porcelet
  
  
  # param_feed_piglet<- matrix(list(MAT_AlimPS1, MAT_AlimPS2,
  #                               ResDc_AlimPS1, ResDc_AlimPS2,
  #                               GES_AlimPS1, GES_AlimPS2,
  #                               Acid_AlimPS1, Acid_AlimPS2,
  #                               EUTR_AlimPS1, EUTR_AlimPS2,
  #                               Energie_AlimPS1, Energie_AlimPS2,
  #                               Surface_AlimPS1, Surface_AlimPS2,
  #                               Cout_AlimPS1, Cout_AlimPS2,
  #                               Feed_phase1_piglet, Feed_phase2_piglet),
  #                          nrow=2, ncol=9, byrow=FALSE)
  # colnames(param_feed_piglet)<-c("MAT","ResDt","GES","Acid",
  #                             "EUTR","Energie","Surface",
  #                             "CoutAlim","Proportion")
  # rownames(param_feed_piglet)<-c("Phase 1","Phase 2")
  
  
  # Calculations
  
  CoutAlim_piglet = round((Cout_AlimPS1*Feed_phase1_piglet
                        + Cout_AlimPS2*Feed_phase2_piglet), digits=2)
  # Co?t moyen de l'aliment porcelet en euros par tonne d'aliments
  
  
  
  # N COMPOUNDS ------------------------
  
  
  # Housing Emissions (kg manure N ex-animal *emission factor) 
  
  N_NH3EmitH_piglet = (RejeteN_piglet / 1000) * 0.24    # kg N-NH3 / piglet
  # Rigolot et al.(2010b): EF is 0.24 kg N-NH3/kg N excreted
  
  N_N2OEmitH_piglet = (RejeteN_piglet / 1000) * 0.002    # kg N-N2O / piglet
  # IPCC (2006): TABLE 10.21-Pit sotorage bellow animal confinements
  
  N_NOxEmitH_piglet = (RejeteN_piglet / 1000) * 0.002     # kg N-NOx / piglet
  # Dämmgen & Hutchings(2008)cited by IPCC(2006)
  
  N_NO3EmitH_piglet =(RejeteN_piglet / 1000) * 0          # kg N-NO3 / piglet
  
  # N leaving after housing emissions (kg N / pig)
  N_St_piglet = (RejeteN_piglet/1000 - N_NH3EmitH_piglet 
              - N_N2OEmitH_piglet - N_NOxEmitH_piglet)
  
  
  
  # Storage Emissions (kg manure N ex-housing *emission factor) 
  
  N_NH3EmitSt_piglet = N_St_piglet * 0.1                # kg N-NH3 / piglet
  # Rigolot et al.(2010b): Table 4
  
  N_N2OEmitSt_piglet = N_St_piglet * 0.001              # kg N-N2O / piglet
  # IPCC (2006)
  
  N_NOxEmitSt_piglet = N_St_piglet * 0.005              # kg N-NOx / piglet
  # Dämmgen & Hutchings(2008)cited by IPCC(2006)
  
  # N leaving after storage emissions (kg N / pig)
  N_field_piglet = (N_St_piglet - N_NH3EmitSt_piglet
                 - N_N2OEmitSt_piglet - N_NOxEmitSt_piglet)
  
  
  
  # Emissions during field application (kg manure N ex-storage *emission factor) 
  
  N_NH3EmitFA_piglet = N_field_piglet * 0.07            # kg N-NH3 / pig
  # Andersen et al. (2001)
  
  N_N2OEmitFA_piglet = N_field_piglet * 0.01            # kg N-N2O / pig
  # IPCC (2006)
  
  N_NOxEmitFA_piglet = N_field_piglet * 0.001           # kg N-NOx / pig
  # Nemecek and Kägi (2007)
  
  # Content of manure after field application (kg N / pig)
  N_AfterFA_piglet = (N_field_piglet - N_NH3EmitFA_piglet
                   - N_N2OEmitFA_piglet - N_NOxEmitFA_piglet)
  
  # Emissions (No3 leaching) after field application
  N_NO3EmitFA_piglet = N_AfterFA_piglet * 0.05                # kg N / pig
  # Garcia-Launay et al.(2014)
  
  # Extraloss as nitrates compared to mineral fertilizer
  NO3_leach_increase_piglet = 0.05 * N_field_piglet          # kg N / pig
  # (IPCC, 2006 ; eventuellement revoir/ publi Danemark)
  
  
  # Total emissions of N compounds (sum of all emissions, kg)
  N_NH3Emit_piglet = ( N_NH3EmitH_piglet
                    + N_NH3EmitSt_piglet 
                    + N_NH3EmitFA_piglet )   
  N_N2OEmit_piglet = ( N_N2OEmitH_piglet 
                    + N_N2OEmitSt_piglet 
                    + N_N2OEmitFA_piglet )
  N_NOxEmit_piglet = ( N_NOxEmitH_piglet 
                    + N_NOxEmitSt_piglet 
                    + N_NOxEmitFA_piglet )
  N_NO3Emit_piglet = ( N_NO3EmitFA_piglet 
                       + NO3_leach_increase_piglet)
  
  
  # Total (kg N / pig)
  N_Emit_piglet = ( N_NH3Emit_piglet 
                 + N_N2OEmit_piglet 
                 + N_NOxEmit_piglet 
                 + N_NO3Emit_piglet )
  
  
  # Fertilizer N substitution - nutrient balance
  Nsubs_piglet = 0.75 * N_field_piglet                       # kg N / pig
  # (Nguyen et al., 2010)
  
  
  # CH4 EMISSIONS ------------------------
  
  ResDc_piglet = ( ((ResDc_AlimPS1 * Feed_phase1_piglet 
                + ResDc_AlimPS2 * Feed_phase2_piglet)/1000)
                * CumulAlim_piglet )
  # kg / piglet

  # CH4 emissions by enteric fermentation (Digested fibre intake)
  
  CH4_Enteric_piglet = (ResDc_piglet * 1000                 # g / piglet  
                     * 670                                  # 670 J / g
                     / 1000000) / 55.65                     # MJ / pig
  # enteric. Rigolot et al.(2010a)
  
  
  # CH4 emissions at storage 
  
  # performances$MOexc : kg of OM excreted
  VS_piglet = MOex_piglet                    # kg, VS equals volatile solids cf equation IPCC 2006
  # Rigolot et al.(2010b)
  
  B0=0.42
  # IPCC (2006) gives a value of 0.45 m3 CH4/kg VS, for Western Europe (and and 0.29 for Latin América)
  
  TempMoy_piglet = 2/3*12+1/3*12
  # 2/3 indoor slurry storage at 20?C and 1/3 outdoor st. at 12?C  (Garcia-Launay et al., 2014)
  
  MCF_piglet = (0.1517*(TempMoy_piglet^2)-2.1569*TempMoy_piglet+24.634)/100    
  # Garcia-Launay et al. (2014): MCFlisier=0.1517*(average temp^2)-2.1569*(average temp)+24.634
  # MCF : methane conversion factor
  
  CH4_St_piglet = VS_piglet*(B0*0.67)*MCF_piglet                              
  # kg of CH4 storage / sow
  # Rigolot et al. (2010b) from IPCC(2006)
  # 0.67 for convert m3 to kg
  
  # Total of CH4 emissions
  CH4Emit_piglet = ( CH4_Enteric_piglet 
                   + CH4_St_piglet ) 
  # kg of CH4 emitted by enteric fermentation and manure storage

  
  # K EMISSIONS ------------------------
  
  K_field_piglet = (RejeteK_piglet            # g / piglet
                    )/ 1000                   # kg of potassium excreted / piglet
  
  
  # PO4 LEACHING ------------------------
  
  P_field_piglet = RejeteP_piglet/1000        # kg P / piglet
  Psubs_piglet = P_field_piglet               # kg P / piglet (100% substitution P mineral)
  
  
  
  # LCA pre-Calculations ----------------
  
  # Feed production (kg of feed / pig)
  
  Electricity_feed_piglet = CumulAlim_piglet * 41 / 1000         # kWh
  NaturalGas_feed_piglet = CumulAlim_piglet * 20.5 / 1000        # kWh 
  RoadTr_feed_piglet = 30 * CumulAlim_piglet / 1000              # t.km
  
  # Housing
  
  Outdoor_housing_piglet = 0     # m2 used in pig production in outdoor system
  Electricity_housing_piglet = 11.5   # kWh
  
  # Slurry
  
  Transport_slurry_piglet = (Amount_slurry_piglet / 1000) * 10         # t.km
  Spreading_slurry_piglet = (Amount_slurry_piglet / 1000)              # t manure spread
  
  
  
  # ------------------1.Climate Change------------------------
  
  GESAlim_piglet = ( (GES_AlimPS1 * Feed_phase1_piglet 
                      + GES_AlimPS2 * Feed_phase2_piglet)
                     *CumulAlim_piglet)
  # kg CO2-eq / pig
  
  # TOTAL (kg CO2-eq / pig)
  CC_piglet= (GESAlim_piglet                               # kg CO2-eq / pig
           + CH4Emit_piglet * 25                           # kg CH4 / pig             -> kg CO2-eq / pig 
           + N_N2OEmit_piglet * 44 / 28 * 298              # kg N_N2O / pig           -> kg CO2-eq / pig
           + (N_NH3Emit_piglet * 0.01)* 44 / 28 * 298      # kg N_NH3 / pig           -> kg CO2-eq / pig
           # NH3 slurry (kg) indirect N2O emission (IPCC 2006)
           # 1% du N-NH3 passe en N-N2O
           + (N_NOxEmit_piglet * 0.01)* 44/28 * 298        # kg N_NOx / pig           -> kg CO2-eq / pig
           # NOx slurry (kg) indirect N2O emission  (IPCC 2006)
           + (N_NO3Emit_piglet * 0.0075)* 44/28 * 298      # kg N_NO3 / pig           -> kg CO2-eq / pig
           # NO3 slurry (kg) indirect N2O emission  (IPCC 2006)
           - Nsubs_piglet * 6.576494                         # kg N / pig               -> kg CO2-eq / pig
           # Substitution N mineral
           - Psubs_piglet * 4.053827                           # kg P / pig               -> kg CO2-eq / pig
           # Substitution P mineral
           - K_field_piglet * 0.732281                       # kg K / pig               -> kg CO2-eq / pig
           # Substitution K mineral 
           + (Electricity_housing_piglet+Electricity_feed_piglet) * 0.00871855             # kWh / pig                -> kg CO2-eq / pig
           + Transport_slurry_piglet * 0.39964633                 # kg CO2 / t.km            -> kg CO2-eq / pig
           # 10 km de transport
           + Spreading_slurry_piglet *  0.001993133                 # kg CO2 / t manure spread -> kg CO2-eq / pig
           # Débit d'épandage = 140 kg/s 
           # Quantité de lisier à épandre = 46 kg par porcelet
           # Donc temps d'épandage  = 46/140 = 0.33 secondes
           # Coeffient de conversion pour 1h d'épandage = 21.743272
           # Donc Coeff pour 1h /3600 *0.33 =  0.001993133
           + NaturalGas_feed_piglet * 0.02060937              # kWh / kg feed            -> kg CO2-eq / pig
           # natural gas for feed production 
           + RoadTr_feed_piglet * 0.13552743 )                 # kgCO2 / t.km             -> kg CO2-eq / pig
           # road transport for feed production

  
  # FEED (kg CO2-eq / pig)
  CC_F_piglet = (GESAlim_piglet                            # kg CO2-eq / pig
              + NaturalGas_feed_piglet * 0.02060937        # kWh / kg feed  -> kg CO2-eq / pig
              + RoadTr_feed_piglet * 0.13552743            # kgCO2 / t.km   -> kg CO2-eq / pig
              + Electricity_feed_piglet * 0.00871855        # kWh / pig      -> kg CO2-eq / pig
              )
  

  # HOUSING (kg CO2-eq / pig)
  CC_H_piglet = ( (CH4_Enteric_piglet + CH4_St_piglet*0.33)*25        # kg CO2-eq / pig
               # 33% de temps de stockage en batiment sous les animaux
               + N_N2OEmitH_piglet * 44 / 28 * 298                     # kg N_N2O / pig -> kg CO2-eq / pig
               + (N_NH3EmitH_piglet * 0.01)* 44 / 28 * 298             # kg N_NH3 / pig -> kg CO2-eq / pig
               + (N_NOxEmitH_piglet * 0.01)* 44/28 * 298               # kg N_NOx / pig -> kg CO2-eq / pig
               + (N_NO3EmitH_piglet * 0.0075)* 44/28 * 298             # kg N_NO3 / pig -> kg CO2-eq / pig
               + Electricity_housing_piglet * 0.00871855 )              # kWh / pig      -> kg CO2-eq / pig

  
  # MANURE (kg CO2-eq / pig)
  CC_M_piglet = CC_piglet - CC_F_piglet - CC_H_piglet
  
  
  # ----------------------2.Acidification------------------------
  
  AcidAlim_piglet = ( (Acid_AlimPS1 * Feed_phase1_piglet 
                       + Acid_AlimPS2 * Feed_phase2_piglet)
                      *CumulAlim_piglet
                      / 1000)                                     # Convert g to kg SO2-eq
  
  # TOTAL (g SO2-eq / pig)
  AC_piglet = ( (AcidAlim_piglet                                  # kg SO2-eq / pig
              + N_NH3Emit_piglet * 17 / 14 * 1.6                  # kg N_NH3 / pig -> kg SO2-eq / pig
              - Nsubs_piglet * 0.019816                             # kg N / pig     -> kg SO2-eq / pig
              - Psubs_piglet * 0.083455                             # kg P / pig     -> kg SO2-eq / pig
              - K_field_piglet * 0.005344                           # kg K / pig     -> kg SO2-eq / pig
              + Transport_slurry_piglet * 0.0024778163                     # t.km           -> kg SO2-eq / pig
              # 10 km de transport
              + Spreading_slurry_piglet * 0.00001351721                    # t              -> kg SO2-eq / pig
              # Débit d'épandage = 140 kg/s 
              # Quantité de lisier à épandre = 46 kg par porcelet
              # Donc temps d'épandage  = 46/140 = 0.33 secondes
              # Coeffient de conversion pour 1h d'épandage = 0.14746051
              # Donc Coeff pour 1h /3600 *0.33 = 0.00001351721
              + RoadTr_feed_piglet * 0.00071396                       # t.km           -> kg SO2-eq / pig
              + (Electricity_housing_piglet+Electricity_feed_piglet) * 0.00004142        # kWh / pig      -> kg SO2-eq / pig
              + NaturalGas_feed_piglet * 0.00007305                   # kWh / kg feed  -> kg CO2-eq / pig
              ) * 1000 )                                              # Convert kg to g SO2-eq


  # FEED (g SO2-eq / pig)
  AC_F_piglet = ( (AcidAlim_piglet                                # kg SO2-eq / kg pig
                + RoadTr_feed_piglet * 0.00071396                 # t.km -> kg SO2-eq / pig
                + NaturalGas_feed_piglet * 0.00007305             # kWh / kg feed  -> kg CO2-eq / pig
                + Electricity_feed_piglet * 0.00004142
                ) * 1000 )                                        # Convert kg to g SO2-eq


  # HOUSING (g SO2-eq / pig)
  AC_H_piglet = ( (N_NH3EmitH_piglet * 17 / 14 * 1.6             # kg SO2-eq / pig
                  + Electricity_housing_piglet * 0.00004142       # kWh / pig      -> kg SO2-eq / pig
                  )*1000 )     

  
  # MANURE (g SO2-eq / pig)                   
  AC_M_piglet = AC_piglet - AC_F_piglet - AC_H_piglet
  
  # -------------------3.Eutrophisation------------------------
  
  EUTRAlim_piglet = ( (EUTR_AlimPS1 * Feed_phase1_piglet
                      + EUTR_AlimPS2 * Feed_phase2_piglet)
                      *CumulAlim_piglet
                      /1000 )         # Convert g to kg PO4-eq
  # kg PO4-eq / pig
  
  # TOTAL (g PO4-eq / pig)
  EU_piglet =( (EUTRAlim_piglet                      # kg PO4-eq / pig
             + N_NH3Emit_piglet * 17 / 14 * 0.35     # kg N_NH3/pig -> kg PO4-eq / pig
             + N_NO3Emit_piglet * 62 / 14 * 0.1      # kg N_NO3/pig -> kg PO4-eq / pig
             - Nsubs_piglet * 0.008919                 # kg N/pig     -> kg PO4-eq / pig
             - Psubs_piglet * 0.096651                 # kg P/pig     -> kg PO4-eq / pig
             - K_field_piglet * 0.001153               # kg K/pig     -> kg PO4-eq / pig
             + Transport_slurry_piglet * 0.00066817162         # t.km         -> kg PO4-eq / pig
             # 10 km de transport
             + Spreading_slurry_piglet * 0.000003457315       # t            -> kg PO4-eq / pig            
             # Débit d'épandage = 140 kg/s 
             # Quantité de lisier à épandre = 46 kg par porcelet
             # Donc temps d'épandage  = 46/140 = 0.33 secondes
             # Coeffient de conversion pour 1h d'épandage = 0.037716159
             # Donc Coeff pour 1h /3600 *0.33 = 0.000003457315
             + RoadTr_feed_piglet * 0.00013245                # t.km         -> kg PO4-eq / pig
             + (Electricity_housing_piglet+Electricity_feed_piglet) * 0.00000341    # kWh / pig    -> kg SO2-eq / pig
             + NaturalGas_feed_piglet * 0.00000215    # kWh / kg feed  -> kg CO2-eq / pig
             ) * 1000)                                          # Convert kg to g PO4-eq
  

  # FEED (g PO4-eq / pig)
  EU_F_piglet = ( (EUTRAlim_piglet                      # kg PO4-eq / pig
                + RoadTr_feed_piglet * 0.00013245       # t.km         -> kg PO4-eq / pig
                + Electricity_feed_piglet * 0.00000341      # kWh / pig    -> kg SO2-eq / pig
                + NaturalGas_feed_piglet * 0.00000215      # kWh / kg feed  -> kg CO2-eq / pig
                )* 1000)                                   # Convert kg to g PO4-eq
  
  
  # HOUSING (g PO4-eq / pig)
  EU_H_piglet = ( (N_NH3EmitH_piglet * 17 / 14 * 0.35     # kg N_NH3/pig -> kg PO4-eq / pig
                + N_NO3EmitH_piglet * 62 / 14 * 0.1       # kg N_NO3/pig -> kg PO4-eq / pig
                + Electricity_housing * 0.00000341              # kWh / pig    -> kg SO2-eq / pig 
                ) * 1000)                                      # Convert kg to g PO4-eq
  
  # MANURE (g PO4-eq / pig)
  EU_M_piglet = EU_piglet - EU_F_piglet - EU_H_piglet
  
  
  # -----------------4.Cumulative Energy Demand------------------
  
  EnergieAlim_piglet = ( (Energie_AlimPS1*Feed_phase1_piglet 
                          + Energie_AlimPS2*Feed_phase2_piglet)
                         * CumulAlim_piglet )
  # MJ / pig
  
  # TOTAL (MJ/ pig)
  CED_piglet = (EnergieAlim_piglet                      # MJ / pig
             - Nsubs_piglet * 59.858776                   # kg N / pig        -> MJ / pig
             - Psubs_piglet * 62.806341                   # kg N / pig        -> MJ / pig
             - K_field_piglet * 12.670042                  # kg N / pig        -> MJ / pig
             + (Electricity_housing_piglet+Electricity_feed_piglet) * 0.90513108      # kWh / pig         -> MJ / pig
             + Transport_slurry_piglet * 5.2284446            # t.km              -> MJ / pig
             # 10 km de transport
             + Spreading_slurry_piglet * 0.03331762            # t                 -> MJ / pig
             # Débit d'épandage = 140 kg/s 
             # Quantité de lisier à épandre = 46 kg par porcelet
             # Donc temps d'épandage  = 46/140 = 0.33 secondes
             # Coeffient de conversion pour 1h d'épandage = 363.46499
             # Donc Coeff pour 1h /3600 *0.33 = 0.03331762
             + NaturalGas_feed_piglet * 0.33626306          # kWh / ton of feed -> MJ / pig
             + RoadTr_feed_piglet * 2.22760640)             # t.km              -> MJ / pig
  

  # FEED (MJ/ pig)
  CED_F_piglet = (EnergieAlim_piglet                    # MJ / pig
               + NaturalGas_feed_piglet * 0.33626306    # kWh / ton of feed -> MJ / pig
               + RoadTr_feed_piglet * 2.22760640        # t.km              -> MJ / pig
               + Electricity_feed_piglet * 0.90513108 )         # kWh / pig         -> MJ / pig
               
  
  # HOUSING (MJ/ pig)
  CED_H_piglet = (Electricity_housing_piglet * 0.90513108)    # kWh / pig         -> MJ / pig
  
  
  # MANURE (MJ/ pig)
  CED_M_piglet = CED_piglet - CED_F_piglet - CED_H_piglet 
  
  
  # ------------------5.Land occupation----------------------
  
  SurfaceAlim_piglet = ( (Surface_AlimPS1 * Feed_phase1_piglet 
                          + Surface_AlimPS2 * Feed_phase2_piglet)
                         * CumulAlim_piglet )
  # m2 year / pig
  
  # TOTAL (m2 year / pig)
  Land_piglet = (SurfaceAlim_piglet                          # m2 year/ pig
              + (Electricity_housing_piglet+Electricity_feed_piglet) * 0.00045690                # kWh / pig -> m2 year/ pig
              + Outdoor_housing_piglet * 1                   # m2 year/ pig  
              + RoadTr_feed_piglet * 0.01129125              # t.km -> m2 year/ pig
              + NaturalGas_feed_piglet * 0.00006732          # kWh / kg feed -> kg CO2-eq / pig
              + Transport_slurry_piglet * 0.074090177                # t.km              -> MJ / pig
              # 10 km de transport
              + Spreading_slurry_piglet * 0.0000341988 )            # t      
              # Débit d'épandage = 140 kg/s 
              # Quantité de lisier à épandre = 46 kg
              # Donc temps d'épandage  = 46/140 = 0.33 secondes
              # Coeffient de conversion pour 1h d'épandage = 0.37307785
              # Donc Coeff pour 1h /3600 *0.33 = 0.0000341988
  

  # FEED (m2 year / pig)  
  Land_F_piglet = (SurfaceAlim_piglet                        # m2 year/ pig             
                + RoadTr_feed_piglet * 0.01129125            # t.km      -> m2 year/ pig
                + Electricity_feed_piglet * 0.00045690          # kWh / pig -> m2 year/ pig
                + NaturalGas_feed_piglet * 0.00006732 )        # kWh / kg feed -> kg CO2-eq / pig

  
  # HOUSING (m2 year / pig)
  Land_H_piglet = (Electricity_housing_piglet * 0.00045690    # kWh / pig -> m2 year/ pig
                + Outdoor_housing_piglet * 1)                # m2 year/ pig
  
  
  # MANURE (m2 year / pig)
  Land_M_piglet = Land_piglet - Land_F_piglet - Land_H_piglet
  
  
  #------------------------------------------------------------
  
  
  # Sur toute la p?riode (de la truie ? la vente du porc ? l'abattoir)
  
  # Number of piglets after post weaning produced by a sow by year  
  Nb_piglets = 28.9*(1-0.024)
  
  
  
  performances$CC_tot = ( (performances$CC + (CC_sow/Nb_piglets) + CC_piglet)
                / performances$PVFin )
  # kg CO2-eq / kg of pig live weight
  
  
  performances$CC_F_tot = ( (performances$CC_F + (CC_F_sow/Nb_piglets) + CC_F_piglet)
                 / performances$PVFin )
  # kg CO2-eq / kg of pig live weight
  
  
  performances$CC_H_tot = ( (performances$CC_H + (CC_H_sow/Nb_piglets) + CC_H_piglet)
                 / performances$PVFin ) 
  # kg CO2-eq / kg of pig live weight
  
  
  performances$CC_M_tot = ( (performances$CC_M + (CC_M_sow/Nb_piglets) + CC_M_piglet)
                 / performances$PVFin ) 
  # kg CO2-eq / kg of pig live weight
  
  
  performances$AC_tot = ( (performances$AC + (AC_sow/Nb_piglets) + AC_piglet)
               / performances$PVFin )
  # g SO2-eq / kg of pig live weight
  
  
  performances$AC_F_tot = ( (performances$AC_F + (AC_F_sow/Nb_piglets) + AC_F_piglet)
                 / performances$PVFin ) 
  # g SO2-eq / kg of pig live weight
  
  
  performances$AC_H_tot = ( (performances$AC_H + (AC_H_sow/Nb_piglets) + AC_H_piglet)
                 / performances$PVFin )
  # g SO2-eq / kg of pig live weight
  
  
  performances$AC_M_tot = ( (performances$AC_M + (AC_M_sow/Nb_piglets) + AC_M_piglet)
                 / performances$PVFin ) 
  # g SO2-eq / kg of pig live weight
  
  
  performances$EU_tot = ( (performances$EU + (EU_sow/Nb_piglets) + EU_piglet)
               / performances$PVFin ) 
  # g PO4-eq / kg of pig live weight
  
  
  performances$EU_F_tot = ( (performances$EU_F + (EU_F_sow/Nb_piglets) + EU_F_piglet)
                 / performances$PVFin ) 
  # g PO4-eq / kg of pig live weight
  
  
  performances$EU_H_tot = ( (performances$EU_H + (EU_H_sow/Nb_piglets) + EU_H_piglet)
                 / performances$PVFin ) 
  # g PO4-eq / kg of pig live weight
  
  
  performances$EU_M_tot = ( (performances$EU_M + (EU_M_sow/Nb_piglets) + EU_M_piglet)
                 / performances$PVFin ) 
  # g PO4-eq / kg of pig live weight
  
  
  performances$CED_tot = ( (performances$CED + (CED_sow/Nb_piglets) + CED_piglet)
                / performances$PVFin ) 
  # MJ / kg of pig live weight
  
  
  performances$CED_F_tot = ( (performances$CED_F + (CED_F_sow/Nb_piglets) + CED_F_piglet)
                  / performances$PVFin ) 
  # MJ / kg of pig live weight
  
  
  performances$CED_H_tot = ( (performances$CED_H + (CED_H_sow/Nb_piglets) + CED_H_piglet)
                  / performances$PVFin ) 
  # MJ / kg of pig live weight pig
  
  
  performances$CED_M_tot = ( (performances$CED_M + (CED_M_sow/Nb_piglets) + CED_M_piglet)
                  / performances$PVFin ) 
  # MJ / kg of pig live weight
  
  
  performances$Land_tot = ( (performances$Land + (Land_sow/Nb_piglets) + Land_piglet)
                 / performances$PVFin ) 
  # m2 year / kg of pig live weight
  
  
  performances$Land_F_tot = ( (performances$Land_F + (Land_F_sow/Nb_piglets) + Land_F_piglet)
                   / performances$PVFin ) 
  # m2 year / kg of pig live weight
  
  
  performances$Land_H_tot = ( (performances$Land_H + (Land_H_sow/Nb_piglets) + Land_H_piglet)
                   / performances$PVFin ) 
  # m2 year / kg of pig live weight
  
  
  performances$Land_M_tot = ( (performances$Land_M + (Land_M_sow/Nb_piglets) + Land_M_piglet)
                   / performances$PVFin ) 
  # m2 year / kg of pig live weight
  

  write.table(performances,paste("performances_",as.character(id),".txt",sep=""),sep=" ", col.names = TRUE,
              row.names = FALSE)
  
  
  
  
  
  
  
  
  
  

  # A L'ECHELLE DE L'ELEVAGE
  
  # IMPACTS ENVIRONNEMENTAUX (/kg de poids vif de porc)---------------------------------
  # En engraissement : tous les porcs (m?me ceux morts avant d'?tre vendus) sont pris en compte
  # En maternit? et post-sevrage : seuls les porcs vendus sont pris en compte
  
  # Number of pigs produced by a sow by year  
  # Nb_pigs = 28.9*(1-0.024)*(1-0.037)
  # Prolificacy = 28.9 piglets weaned / sow / year (GTTT and GTE 2014)
  # Post weaning mortality rate = 2.4%
  # Mortality rate fattening = 3.7%

# Number of piglets after post weaning produced by a sow by year  
  Nb_piglets = 28.9*(1-0.024)

# Nb_piglets = Nb_pigs
# 1 piglet = X pigs
# X = Nb_pigs / Nb_piglets

  
  # Sur toute la p?riode (de la truie ? la vente du porc ? l'abattoir)-------
  
  CC =  round(( sum(performances$CC + (CC_sow/Nb_piglets) + CC_piglet)
          / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 )
  # kg CO2-eq / kg of pig live weight
  

  CC_F = round(( sum(performances$CC_F + (CC_F_sow/Nb_piglets) + CC_F_piglet)
           / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 )
  # kg CO2-eq / kg of pig live weight
  
  
  CC_H = round(( sum(performances$CC_H + (CC_H_sow/Nb_piglets) + CC_H_piglet)
           / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # kg CO2-eq / kg of pig live weight
  
  
  CC_M = round(( sum(performances$CC_M + (CC_M_sow/Nb_piglets) + CC_M_piglet)
           / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # kg CO2-eq / kg of pig live weight
  
  
  AC = round(( sum(performances$AC + (AC_sow/Nb_piglets) + AC_piglet)
         / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 )
  # g SO2-eq / kg of pig live weight
  
  
  AC_F = round(( sum(performances$AC_F + (AC_F_sow/Nb_piglets) + AC_F_piglet)
           / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # g SO2-eq / kg of pig live weight
  
  
  AC_H = round(( sum(performances$AC_H + (AC_H_sow/Nb_piglets) + AC_H_piglet)
           / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 )
  # g SO2-eq / kg of pig live weight
  
  
  AC_M = round(( sum(performances$AC_M + (AC_M_sow/Nb_piglets) + AC_M_piglet)
           / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # g SO2-eq / kg of pig live weight
  
  
  EU = round(( sum(performances$EU + (EU_sow/Nb_piglets) + EU_piglet)
         / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # g PO4-eq / kg of pig live weight
  
  
  EU_F = round(( sum(performances$EU_F + (EU_F_sow/Nb_piglets) + EU_F_piglet)
           / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # g PO4-eq / kg of pig live weight
  
  
  EU_H = round(( sum(performances$EU_H + (EU_H_sow/Nb_piglets) + EU_H_piglet)
           / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # g PO4-eq / kg of pig live weight
  
  
  EU_M = round(( sum(performances$EU_M + (EU_M_sow/Nb_piglets) + EU_M_piglet)
           / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # g PO4-eq / kg of pig live weight
  
  
  CED = round(( sum(performances$CED + (CED_sow/Nb_piglets) + CED_piglet)
          / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # MJ / kg of pig live weight
  
  
  CED_F = round(( sum(performances$CED_F + (CED_F_sow/Nb_piglets) + CED_F_piglet)
            / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # MJ / kg of pig live weight
  
  
  CED_H = round(( sum(performances$CED_H + (CED_H_sow/Nb_piglets) + CED_H_piglet)
            / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # MJ / kg of pig live weight pig
  
  
  CED_M = round(( sum(performances$CED_M + (CED_M_sow/Nb_piglets) + CED_M_piglet)
            / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # MJ / kg of pig live weight
  
  
  Land = round(( sum(performances$Land + (Land_sow/Nb_piglets) + Land_piglet)
           / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # m2 year / kg of pig live weight
  
  
  Land_F = round(( sum(performances$Land_F + (Land_F_sow/Nb_piglets) + Land_F_piglet)
             / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # m2 year / kg of pig live weight
  
  
  Land_H = round(( sum(performances$Land_H + (Land_H_sow/Nb_piglets) + Land_H_piglet)
             / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # m2 year / kg of pig live weight
  
  
  Land_M = round(( sum(performances$Land_M + (Land_M_sow/Nb_piglets) + Land_M_piglet)
             / sum(performances$PVFin[performances$Depart!="mort"]) ), digits=3 ) 
  # m2 year / kg of pig live weight
  
  
  # Phase d'engraissement-----------------------------------------------------
  
  CC_Eng =  round(( sum(performances$CC)
          / sum(performances$PVFin[performances$Depart!="mort"]
                -performances$PVInit[performances$Depart!="mort"]) ), digits=3 )
  # kg CO2-eq / kg of live weight gain
  
  
  CC_F_Eng = round(( sum(performances$CC_F)
           / sum(performances$PVFin[performances$Depart!="mort"]
                 -performances$PVInit[performances$Depart!="mort"]) ), digits=3 )
  # kg CO2-eq / kg of live weight gain
  
  
  CC_H_Eng = round(( sum(performances$CC_H)
           / sum(performances$PVFin[performances$Depart!="mort"]
                 -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # kg CO2-eq / kg of live weight gain
  
  
  CC_M_Eng = round(( sum(performances$CC_M)
           / sum(performances$PVFin[performances$Depart!="mort"]
                 -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # kg CO2-eq / kg of live weight gain
  
  
  AC_Eng = round(( sum(performances$AC)
         / sum(performances$PVFin[performances$Depart!="mort"]
               -performances$PVInit[performances$Depart!="mort"]) ), digits=3 )
  # g SO2-eq / kg of live weight gain
  
  
  AC_F_Eng = round(( sum(performances$AC_F)
           / sum(performances$PVFin[performances$Depart!="mort"]
                 -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # g SO2-eq / kg of live weight gain
  
  
  AC_H_Eng = round(( sum(performances$AC_H)
           / sum(performances$PVFin[performances$Depart!="mort"]
                 -performances$PVInit[performances$Depart!="mort"]) ), digits=3 )
  # g SO2-eq / kg of live weight gain
  
  
  AC_M_Eng = round(( sum(performances$AC_M)
           / sum(performances$PVFin[performances$Depart!="mort"]
                 -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # g SO2-eq / kg of live weight gain
  
  
  EU_Eng = round(( sum(performances$EU)
         / sum(performances$PVFin[performances$Depart!="mort"]
               -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # g PO4-eq / kg of live weight gain
  
  
  EU_F_Eng = round(( sum(performances$EU_F)
           / sum(performances$PVFin[performances$Depart!="mort"]
                 -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # g PO4-eq / kg of live weight gain
  
  
  EU_H_Eng = round(( sum(performances$EU_H)
           / sum(performances$PVFin[performances$Depart!="mort"]
                 -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # g PO4-eq / kg of live weight gain
  
  
  EU_M_Eng = round(( sum(performances$EU_M)
           / sum(performances$PVFin[performances$Depart!="mort"]
                 -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # g PO4-eq / kg of live weight gain
  
  
  CED_Eng = round(( sum(performances$CED)
          / sum(performances$PVFin[performances$Depart!="mort"]
                -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # MJ / kg of live weight gain
  
  
  CED_F_Eng = round(( sum(performances$CED_F)
            / sum(performances$PVFin[performances$Depart!="mort"]
                  -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # MJ / kg of live weight gain
  
  
  CED_H_Eng = round(( sum(performances$CED_H)
            / sum(performances$PVFin[performances$Depart!="mort"]
                  -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # MJ / kg of live weight gain
  
  
  CED_M_Eng = round(( sum(performances$CED_M)
            / sum(performances$PVFin[performances$Depart!="mort"]
                  -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # MJ / kg of live weight gain
  
  
  Land_Eng = round(( sum(performances$Land)
           / sum(performances$PVFin[performances$Depart!="mort"]
                 -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # m2 year / kg of live weight gain
  
  
  Land_F_Eng = round(( sum(performances$Land_F)
             / sum(performances$PVFin[performances$Depart!="mort"]
                   -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # m2 year / kg of live weight gain
  
  
  Land_H_Eng = round(( sum(performances$Land_H)
             / sum(performances$PVFin[performances$Depart!="mort"]
                   -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # m2 year / kg of live weight gain
  
  
  Land_M_Eng = round(( sum(performances$Land_M)
             / sum(performances$PVFin[performances$Depart!="mort"]
                   -performances$PVInit[performances$Depart!="mort"]) ), digits=3 ) 
  # m2 year / kg of live weight gain
  
  
  
  # IMPACTS ECONOMIQUES----------------------------------------------------------
  # selon les cas, tous ou une partie des porcs sont pris en compte
  
  prixvente = round( sum(performances$prixvente)
                / nrow(performances[performances$Depart!="mort",]), digits=2 )
  # euros / kg of carcass pig
  # Seuls les porcs vendus sont pris en compte
  
  produit = round( sum(performances$produit)
              / nrow(performances[performances$Depart!="mort",]), digits=2 )
  # euros / carcass pig
  # Seuls les porcs vendus sont pris en compte
  

  marge_Eng_reelle = round( sum(performances$marge_Eng)
                     / nrow(performances[performances$Depart!="mort",]),
                     digits=2 )
  # euros / pig

  marge_Eng_estimee = round( mean(performances$marge_Eng[performances$Depart!="mort"]),
                             digits=2 )
  # euros / pig
  
  
  CoutAlimTot = round( sum(performances$CoutAlimTot
                           + (CoutAlim_sow/1000*CumulAlim_sow/Nb_piglets)
                           + (CoutAlim_piglet/1000*CumulAlim_piglet) )
                          / nrow(performances[performances$Depart!="mort",]), 
                          digits=2 )
  # euros / pig
  # Charges alimentaires : tous les porcs sont pris en compte (vendus et morts avant la vente)
  
  CoutAlim_Eng = round( sum(performances$CoutAlimTot)
                       / nrow(performances[performances$Depart!="mort",]), 
                       digits=2 )
  # euros / pig
  # Charges alimentaires : tous les porcs sont pris en compte (vendus et morts avant la vente)
  
  
  # RESULTATS TECHNIQUES--------------------------------------------------
  
  IC_Eng = round( mean(performances$IC[performances$Depart!="mort"]),
                  digits=3 )
  # kg feed / kg of live weight pig
  # Seuls les porcs en engraissement vendus sont pris en compte
  
  
  RejeteN_Eng = round( mean(performances$RejeteN[performances$Depart!="mort"]), 
               digits=3 )
  # g N / pig killed at the slaughterhouse

  
  RejeteP_Eng = round( mean(performances$RejeteP[performances$Depart!="mort"]), 
                            digits=3 )
  # g N / pig killed at the slaughterhouse

  
  MOexc_Eng = round( mean(performances$MOexc[performances$Depart!="mort"]),
                     digits=3 )
  # kg excreted MO / pig killed at the slaughterhouse

  
  Resdig_Eng = round( mean(performances$Resdig[performances$Depart!="mort"]),
                      digits=3 )
  # kg digestible residues / pig killed at the slaughterhouse
  
  

  PVFin = round( mean(performances$PVFin[performances$Depart!="mort"]),
                 digits=3 )
  
  AgeFin = round( mean(performances$AgeFin[performances$Depart!="mort"]),
                  digits=3)
  
  Gamme <- round( length(performances$PVFin[(performances$Depart!="mort")&
                                         (performances$PVFin*tauxCC_fr>=82)&
                                         (performances$PVFin*tauxCC_fr<105)]) 
           /length(performances$PVFin[performances$Depart!="mort"]),
           digits=4 )

  Coeur_gamme <- round( length(performances$PVFin[(performances$Depart!="mort")&
                                               (performances$PVFin*tauxCC_fr>=87)&
                                               (performances$PVFin*tauxCC_fr<99)]) 
                    /length(performances$PVFin[performances$Depart!="mort"]),
                    digits=4 )
  
  Lourds = round( length(performances$PVFin[(performances$Depart!="mort")&
                                         (performances$PVFin*tauxCC_fr>=105)]) 
                /length(performances$PVFin[performances$Depart!="mort"]),
                digits=4 )
  
  Legers = round( length(performances$PVFin[(performances$Depart!="mort")&
                                         (performances$PVFin*tauxCC_fr<82)]) 
             /length(performances$PVFin[performances$Depart!="mort"]),
             digits=4 )

  Perte = round( nrow(performances[performances$Depart=="mort",])
            /nrow(performances), 
            digits=4 )
     
  


  resultats = c(id, CC, CC_F, CC_H, CC_M,
                AC, AC_F, AC_H, AC_M,
                EU, EU_F, EU_H, EU_M,
                CED, CED_F, CED_H, CED_M,
                Land, Land_F, Land_H, Land_M, 
                CC_Eng, CC_F_Eng, CC_H_Eng, CC_M_Eng,
                AC_Eng, AC_F_Eng, AC_H_Eng, AC_M_Eng,
                EU_Eng, EU_F_Eng, EU_H_Eng, EU_M_Eng,
                CED_Eng, CED_F_Eng, CED_H_Eng, CED_M_Eng,
                Land_Eng, Land_F_Eng, Land_H_Eng, Land_M_Eng,
                prixvente, produit, 
                marge_Eng_reelle, marge_Eng_estimee,
                CoutAlimTot, CoutAlim_Eng,
                IC_Eng, RejeteN_Eng, RejeteP_Eng,
                MOexc_Eng, Resdig_Eng,
                PVFin, AgeFin, Gamme,
                Coeur_gamme, Lourds, Legers, Perte)


  return(resultats)
}