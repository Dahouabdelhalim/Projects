!Input file Linear Inverse Model Coral Reef Aqaba Jordan Season: Summer
!Used in van Hoytema et al. (2023) DOI 10.1007/s00338-022-02339-3		
!To be used with R package LIM https://cran.r-project.org/web/packages/LIM/index.html
										
## parameters											
!water column organic C composition											
minphy	=	0.004	!min	and	max	percentages	that pelagic compartments make up of total organic C in water column
maxphy	=	0.006									
minbac	=	0.007									
maxbac	=	0.014									
minpoc	=	0.032									
maxpoc	=	0.095									
mindoc	=	0.886									
maxdoc	=	0.957									

!phytoplankton											
mingppphy	=	0									
maxgppphy	=	16.1				! pelagic GPP					
minRE	=	0.05	          !min	and	max	fraction	of	phytoplankton	GPP	going to respiration	
maxRE	=	0.3									
minEE	=	0.05	!min	and	max	phytoplankton	excretion	efficiency	of	DOC	
maxEE	=	0.5									

!pelagic respiration											
minresppel	=	68.7									
maxresppel	=	90.8									

!pelagic bacteria											
minprodbac	=	1.0	!min	and	max	bacterial	productivity	based	on	biomass	
maxprodbac	=	18.4									
minBGE	=	0.1	      !min	and	max	bacterial	growth	efficiency 			
maxBGE	=	0.6									

!fluxes per benthic compartment below follow the following naming convention min or max contraint, compartment from, compartment to
!hard corals											
mindichco	=	102.1									
maxdichco	=	166.8									
minhcodic	=	98.4									
maxhcodic	=	153.8									
minhcodoc	=	0.0									
maxhcodoc	=	13.0																		
minhcopoc	=	1.5									
maxhcopoc	=	16.3																

!soft corals											
mindicsco	=	19.8									
maxdicsco	=	35.5									
minscodic	=	19.2									
maxscodic	=	34.6									
minscodoc	=	0.0									
maxscodoc	=	0.6																	
minscopoc	=	0.8									
maxscopoc	=	5.5									

!coral POC release fraction that dissolves into DOC											
minCDC = 0.56         											
maxCDC = 0.8											
											

!coral heterotrophic feeding											
minhet	=	0.2	!min	and	max	fraction	of	coral	respiration	contributed	
maxhet	=	0.5									
scsa = 0.257 !soft coral surface area of not-Xeniidae fraction											

!coral mucus to sediment fluxes percentage released coral POC stright to Sediment OC											
minhcosoc = 0.34											
maxhcosoc = 0.63											

!macroalgae											
mindicmac	=	18.3									
maxdicmac	=	25.5									
minmacdic	=	4.6									
maxmacdic	=	6.7									
minmacdoc	=	0.0									
maxmacdoc	=	1.8																	
minmacpoc	=	0.6									
maxmacpoc	=	1.0																	
expmac = 8.8											

!turfalgae											
mindictur	=	3.1									
maxdictur	=	4.3									
minturdic	=	1.0									
maxturdic	=	1.5									
minturdoc	=	0.0									
maxturdoc	=	0.4																	
minturpoc	=	0.2									
maxturpoc	=	0.4																	
exptur	=	0.8									

!coral rock											
minphyroc	=	1.1									
maxphyroc	=	1.6									
mindicroc	=	13.7									
maxdicroc	=	28.0									
minrocdic	=	9.1									
maxrocdic	=	23.5									
minrocdoc	=	0.0									
maxrocdoc	=	2.3																	
minrocpoc	=	0.5									
maxrocpoc	=	3.1									

!sediment organic carbon											
mindicsoc	=	4.9									
maxdicsoc	=	8.1									
minsocdic	=	3.3									
maxsocdic	=	4.7									
minAE = 0.40            !min and max assimilation efficiency											
maxAE = 0.80											
minPE = 0.30            !min and max productivity efficiency											
maxPE = 0.60											
minpocsoc = 0.12        !percentages that hard coral to sediment POC flow constitutes of bulk POC sedimentation naumann et al. 2012 											
maxpocsoc = 0.62     											
expsoc = 1.2										
##end parameters											




##equations

expmac = expomacr

exptur = expoturf

expsoc = exposoc

##end equations

##inequalities

inphy = [minphy,maxphy]*intotal            !fractional breakdown of inflow of C into the reef
inbac = [minbac,maxbac]*intotal
inpoc = [minpoc,maxpoc]*intotal
indoc = [mindoc,maxdoc]*intotal

ppphy = [mingppphy,maxgppphy]                        !pelagic gpp
respphy = [minRE,maxRE] * ppphy             !metabolic constraints phytoplankton
excphy = [minEE,maxEE] * ppphy

pelresp = [minresppel,maxresppel]           !pelagic r

prodbac = [minBGE,maxBGE] * uptabac         !metabolic constraints pelagic bacteria
prodbac = [minprodbac,maxprodbac]

!hard corals
pphc	=	[mindichco,maxdichco] 
rehc = [minhcodic,maxhcodic]
hetcor = [minhet,maxhet]*rehc
hcpcreleased = [minhcopoc,maxhcopoc]
hco_poc -> doc = [minCDC,maxCDC] * hcpcreleased
hco_poc -> soc = [minhcosoc,maxhcosoc] * hcpcreleased
hcdc = [minhcodoc,maxhcodoc]           

!soft corals
ppsc	=	[mindicsco,maxdicsco]
resc = [minscodic,maxscodic]
hetscor = [minhet,maxhet]*resc*scsa
scpcreleased = [minscopoc,maxscopoc]
sco_poc -> doc = [minCDC,maxCDC] * scpcreleased
scdc = [minscodoc,maxscodoc]

!macroalgae
ppma	=	[mindicmac,maxdicmac]
rema = [minmacdic,maxmacdic]
mapc = [minmacpoc,maxmacpoc]
madc = [minmacdoc,maxmacdoc]

!turf algae
ppta	=	[mindictur,maxdictur]
reta = [minturdic,maxturdic]
tapc = [minturpoc,maxturpoc]
tadc = [minturdoc,maxturdoc]

!coral rock
ppro	=	[mindicroc,maxdicroc]
rero = [minrocdic,maxrocdic]
ropc = [minrocpoc,maxrocpoc]
rodc = [minrocdoc,maxrocdoc]
phyroc = [minphyroc,maxphyroc]  

!sediment organic carbon
ppsoc	= [mindicsoc,maxdicsoc]
resoc = [minsocdic,maxsocdic]
assisoc = [minAE,maxAE]*uptasoc           !metabolic constraints sediment life
prodsoc = [minPE,maxPE]*assisoc
hco_poc -> soc = [minpocsoc,maxpocsoc]*total_depo
##end inequalities

##variables
inphy = inf -> phy                         !pelagic inflows from offshore
inbac = inf -> bac
inpoc = inf -> poc
indoc = inf -> doc

intotal = inphy + inbac + inpoc + indoc

pelresp = phy -> dic + bac -> dic        !water column respiration

!phytoplankton
respphy = phy -> dic                    !phytoplankton respiration
ppphy = dic -> phy                       !phy photosynthesis
excphy = phy -> doc         !excretion DOC by PHY

!pelagic bacteria
respbac = bac -> dic                    !pelagic bacteria respiration
uptabac = flowto(bac) - inf -> bac
prodbac = uptabac - respbac

pphc	=	dic -> hco   
rehc = hco -> dic
hetcor = poc -> hco + bac -> hco
hcpcreleased = hco_poc -> poc + hco_poc -> doc + hco_poc -> soc                      
hcdc = hco -> doc                       

ppsc	=	dic -> sco
resc = sco -> dic
hetscor = poc -> sco + bac -> sco
scpcreleased = sco_poc -> poc + sco_poc -> doc 
scdc = sco -> doc

ppma	=	dic -> mac
rema = mac -> dic
mapc = mac -> poc
madc = mac -> doc
expomacr = mac -> exp

ppta	=	dic -> tur
reta = tur -> dic
tapc = tur -> poc
tadc = tur -> doc
expoturf = tur -> exp 

ppro	=	dic -> roc
rero = roc -> dic
ropc = roc -> poc
rodc = roc -> doc
phyroc = phy -> roc   

ppsoc	=	dic -> soc
resoc = soc -> dic
pocsoc = poc -> soc
uptasoc = flowto(soc) - dic -> soc
assisoc = uptasoc - soc -> exp
prodsoc = assisoc - resoc
exposoc = soc -> exp
total_depo = poc -> soc + hco_poc -> soc
##end variables

##stocks##

hco              !hard corals
hco_poc          !hard coral poc released
sco              !soft corals
sco_poc          !soft coral poc released
mac              !macroalgae
tur              !turf algae
roc              !coralrock biota
soc              !sediment organic carbon

phy              !phytoplankton
bac              !bacteria
poc              !pelagic particulate organic carbon
doc              !pelagic dissolved organic carbon


##end stocks##

##externals##

dic
inf              !inflow from outside of reef
exp
##end externals##



##flows##

!inflow
inf -> phy         !inflows
inf -> poc
inf -> doc
inf -> bac

!pelagic compartments
dic -> phy        !pelagic GPP and respiration
phy -> dic
phy -> doc
phy -> poc

doc -> bac
bac -> doc
bac -> poc
bac -> dic

poc -> doc

!hard corals
dic -> hco       
bac -> hco       
poc -> hco      
doc -> hco       
hco -> dic 
hco -> hco_poc
hco_poc -> poc
hco_poc -> soc      !direct coral mucus string sedimentation following naumann et al. 2012
hco_poc -> doc
hco -> doc      


!soft corals
dic -> sco      
bac -> sco      
poc -> sco    
doc -> sco      
sco -> dic
sco -> sco_poc
sco_poc -> poc    
sco_poc -> doc
sco -> doc       

!macroalgae
dic -> mac      
mac -> dic     
mac -> poc    
mac -> doc       
doc -> mac
poc -> mac
mac -> exp

!turf algae
dic -> tur      
tur -> dic       
tur -> poc      
tur -> doc 
poc -> tur
doc -> tur
tur -> exp

!coralrock
dic -> roc       
phy -> roc       
bac -> roc    
poc -> roc      
doc -> roc    
roc -> dic     
roc -> poc      
roc -> doc

!sediment
dic -> soc
poc -> soc     
doc -> soc
soc -> dic     
soc -> doc      
soc -> exp      

##end flows##

