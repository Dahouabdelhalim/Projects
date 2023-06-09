#Load Claddis and other required packages
library(Claddis)
library(paleotree)
library(strap)

#Read in data file
nexus.data <- ReadMorphNexus("Disparity_analysis_input_data_FINAL.txt")

#Generate a distance matrix
dist.data <- MorphDistMatrix(nexus.data)

#Check for incalculable distances
any(dist.data$comp.char.matrix == 0)

#Remove taxa to avoid distance matrix without gaps
trimmed.max.data <- TrimMorphDistMatrix(dist.data$max.dist.matrix)
trimmed<-trimmed.max.data$removed.taxa

#Define our groups including all taxa 
pseudosuchians_minusphyto <-c("Gracilisuchus_stipanicicorum", "Turfanosuchus_dabanensis", "Ornithosuchus_longidens", "Riojasuchus_tenuisceps", "Revueltosaurus_callenderi", "Stagonolepis_robertsoni", "Aetosaurus_ferratus", "Longosuchus_meadei", "Ticinosuchus_ferox", "Qianosuchus_mixtus", "Arizonasaurus_babbitti", "Poposaurus_gracilis_yale", "Lotosaurus_adentus", "Sillosuchus_longicervix", "Effigia_okeeffeae", "Shuvosaurus_inexpectatus", "Combined_Prestosuchus", "Saurosuchus_galilei", "Batrachotomus_kuperferzellensis", "Fasolasuchus_tenax", "Rauisuchus_triradentes", "Postosuchus_kirkpatricki", "Postosuchus_alisonae", "CM_73372", "Hesperosuchus_agilis", "Hesperosuchus_agilis_", "Dromicosuchus_grallator", "Dibothrosuchus_elaphros", "Terrestrisuchus_gracilis", "Mandasuchus_total", "Nundasuchus", "Vivaron_haydeni","Teratosaurus_suevicus","Yonghesuchus_sangbiensis","Xilousuchus_sapingensis","Polonosuchus_silesiacus","Poposaurus_gracilis_holotype")
pseudosuchians_withphyto <-c("Diandongosuchus", "Parasuchus_hislopi","Smilosuchus_gregorii","Pseudopalatus_pristinus", "Gracilisuchus_stipanicicorum", "Turfanosuchus_dabanensis", "Ornithosuchus_longidens", "Riojasuchus_tenuisceps", "Revueltosaurus_callenderi", "Stagonolepis_robertsoni", "Aetosaurus_ferratus", "Longosuchus_meadei", "Ticinosuchus_ferox", "Qianosuchus_mixtus", "Arizonasaurus_babbitti", "Poposaurus_gracilis_yale", "Lotosaurus_adentus", "Sillosuchus_longicervix", "Effigia_okeeffeae", "Shuvosaurus_inexpectatus", "Combined_Prestosuchus", "Saurosuchus_galilei", "Batrachotomus_kuperferzellensis", "Fasolasuchus_tenax", "Rauisuchus_triradentes", "Postosuchus_kirkpatricki", "Postosuchus_alisonae", "CM_73372", "Hesperosuchus_agilis", "Hesperosuchus_agilis_", "Dromicosuchus_grallator", "Dibothrosuchus_elaphros", "Terrestrisuchus_gracilis", "Mandasuchus_total", "Nundasuchus", "Vivaron_haydeni","Teratosaurus_suevicus","Yonghesuchus_sangbiensis","Xilousuchus_sapingensis","Polonosuchus_silesiacus","Poposaurus_gracilis_holotype")
avemetatarsalians <- c("Eudimorphodon_ranzii", "Lagerpeton_chanarensis", "Marasuchus_lilloensis", "Eucoelophysis_baldwini", "Sacisaurus_agudoensis", "Lewisuchus_Pseudolagosuchus", "Eocursor_parvus", "Silesaurus_opolensis", "Pisanosaurus_mertii",  "Herrerasaurus_ischigualastensis", "Staurikosaurus_pricei", "Eoraptor_lunensis", "Saturnalia_tupiniquim", "Plateosaurus_engelhardti", "Efraasia_minor", "Tawa_hallae", "Coelophysis_bauri", "Asilisaurus_kongwe_combined", "Teleocrater_combined", "Yarasuchus_deccanensis", "Saltopus", "Diodoris_", "Scleromochlus", "Dongusuchus","Lutungutali_sitwensis","Spondylosoma","Dromomeron_romeri","Dromomeron_gregorii","Dromomeron_gigas","PVSJ_883","Ignotosaurus_fragilis","Technosaurus","Sanjuansaurus_gordilloi","Guaibasaurus_candelariensis","Pampadromaeus_barberenai","Panphagia_protos")
pseudosuchians_MidTr <- c("Turfanosuchus_dabanensis", "Ticinosuchus_ferox", "Arizonasaurus_babbitti", "Lotosaurus_adentus", "Mandasuchus_total", "Nundasuchus", "Gracilisuchus_stipanicicorum", "Batrachotomus_kuperferzellensis", "Qianosuchus_mixtus")     
pseudosuchians_MidTr_withphyto <- c("Diandongosuchus","Turfanosuchus_dabanensis", "Ticinosuchus_ferox", "Arizonasaurus_babbitti", "Lotosaurus_adentus", "Mandasuchus_total", "Nundasuchus", "Gracilisuchus_stipanicicorum", "Batrachotomus_kuperferzellensis", "Qianosuchus_mixtus")     
pseudosuchians_carnian <- c("Combined_Prestosuchus", "Sillosuchus_longicervix", "Saurosuchus_galilei", "Rauisuchus_triradentes","Polonosuchus_silesiacus","Yonghesuchus_sangbiensis")
pseudosuchians_carnian_withphyto <- c("Parasuchus_hislopi","Combined_Prestosuchus", "Sillosuchus_longicervix", "Saurosuchus_galilei", "Rauisuchus_triradentes","Polonosuchus_silesiacus","Yonghesuchus_sangbiensis")
pseudosuchians_EarlyNorian <- c("Ornithosuchus_longidens", "Stagonolepis_robertsoni", "Longosuchus_meadei", "Poposaurus_gracilis_yale", "Shuvosaurus_inexpectatus", "Postosuchus_kirkpatricki", "Postosuchus_alisonae", "Hesperosuchus_agilis", "Dromicosuchus_grallator","Poposaurus_gracilis_holotype") 
pseudosuchians_EarlyNorian_withphyto <- c("Smilosuchus_gregorii","Pseudopalatus_pristinus", "Ornithosuchus_longidens", "Stagonolepis_robertsoni", "Longosuchus_meadei", "Poposaurus_gracilis_yale", "Shuvosaurus_inexpectatus", "Postosuchus_kirkpatricki", "Postosuchus_alisonae", "Hesperosuchus_agilis", "Dromicosuchus_grallator","Poposaurus_gracilis_holotype") 
pseudosuchians_lateNorian <- c("Revueltosaurus_callenderi", "Aetosaurus_ferratus", "Fasolasuchus_tenax", "Terrestrisuchus_gracilis", "Vivaron_haydeni", "CM_73372", "Hesperosuchus_agilis_","Teratosaurus_suevicus")
#Note inclusion of Pseudopalatus pristinus in late Norian-Rhaetian bin to represent phytosaurs from that interval
pseudosuchians_lateNorian_withphyto <- c("Pseudopalatus_pristinus", "Revueltosaurus_callenderi", "Aetosaurus_ferratus", "Fasolasuchus_tenax", "Terrestrisuchus_gracilis", "Vivaron_haydeni", "CM_73372", "Hesperosuchus_agilis_","Teratosaurus_suevicus")
avemetatarsalians_MidTr <- c("Teleocrater_combined", "Yarasuchus_deccanensis", "Asilisaurus_kongwe_combined","Dongusuchus","Lutungutali_sitwensis")
avemetatarsalians_carnian <- c("Saturnalia_tupiniquim", "Lagerpeton_chanarensis", "Marasuchus_lilloensis", "Lewisuchus_Pseudolagosuchus", "Pisanosaurus_mertii", "Eoraptor_lunensis", "Herrerasaurus_ischigualastensis", "Staurikosaurus_pricei","Diodoris_","Spondylosoma","PVSJ_883","Ignotosaurus_fragilis","Sanjuansaurus_gordilloi","Guaibasaurus_candelariensis","Pampadromaeus_barberenai","Panphagia_protos") 
avemetatarsalians_earlyNorian<-c("Saltopus", "Scleromochlus", "Plateosaurus_engelhardti", "Sacisaurus_agudoensis", "Silesaurus_opolensis","Diodoris_","Dromomeron_gregorii","Technosaurus")
avemetatarsalians_lateNorian<-c("Efraasia_minor", "Tawa_hallae", "Coelophysis_bauri", "Eudimorphodon_ranzii", "Eucoelophysis_baldwini", "Eocursor_parvus","Dromomeron_romeri","Dromomeron_gigas")

#Define our groups excluding taxa trimmed by TrimMorphDistMatrix
pseudosuchians_minusphyto.trimmed <- setdiff(pseudosuchians_minusphyto,trimmed) 
pseudosuchians_withphyto.trimmed <- setdiff(pseudosuchians_withphyto,trimmed)
avemetatarsalians.trimmed <- setdiff(avemetatarsalians,trimmed)
pseudosuchians_MidTr.trimmed <- setdiff(pseudosuchians_MidTr,trimmed) 
pseudosuchians_MidTr_withphyto.trimmed <- setdiff(pseudosuchians_MidTr_withphyto,trimmed)
pseudosuchians_carnian.trimmed <- setdiff(pseudosuchians_carnian,trimmed) 
pseudosuchians_carnian_withphyto.trimmed <- setdiff(pseudosuchians_carnian_withphyto,trimmed) 
pseudosuchians_EarlyNorian.trimmed <- setdiff(pseudosuchians_EarlyNorian,trimmed) 
pseudosuchians_EarlyNorian_withphyto.trimmed <- setdiff(pseudosuchians_EarlyNorian_withphyto,trimmed) 
pseudosuchians_lateNorian.trimmed <- setdiff(pseudosuchians_lateNorian,trimmed)
pseudosuchians_lateNorian_withphyto.trimmed <- setdiff(pseudosuchians_lateNorian_withphyto,trimmed)
avemetatarsalians_MidTr.trimmed <- setdiff(avemetatarsalians_MidTr,trimmed) 
avemetatarsalians_carnian.trimmed <- setdiff(avemetatarsalians_carnian,trimmed)   
avemetatarsalians_earlyNorian.trimmed <- setdiff(avemetatarsalians_earlyNorian,trimmed) 
avemetatarsalians_lateNorian.trimmed <- setdiff(avemetatarsalians_lateNorian,trimmed) 

#Create a list containing all these groups
taxongroups<-list(pseudosuchians_minusphyto.trimmed, pseudosuchians_withphyto.trimmed, avemetatarsalians.trimmed, pseudosuchians_MidTr.trimmed, pseudosuchians_MidTr_withphyto.trimmed, pseudosuchians_carnian.trimmed, pseudosuchians_carnian_withphyto.trimmed, pseudosuchians_EarlyNorian.trimmed, pseudosuchians_EarlyNorian_withphyto.trimmed, pseudosuchians_lateNorian.trimmed, pseudosuchians_lateNorian_withphyto.trimmed, avemetatarsalians_MidTr.trimmed, avemetatarsalians_carnian.trimmed, avemetatarsalians_earlyNorian.trimmed, avemetatarsalians_lateNorian.trimmed)

#Calculate WMPD for each of these groups
results<-matrix(nrow=15, ncol=4,
                dimnames= list(c("pseudosuchians_minusphyto.trimmed", "pseudosuchians_withphyto.trimmed", "avemetatarsalians.trimmed", "pseudosuchians_MidTr.trimmed", "pseudosuchians_MidTr_withphyto.trimmed", "pseudosuchians_carnian.trimmed", "pseudosuchians_carnian_withphyto.trimmed", "pseudosuchians_EarlyNorian.trimmed", "pseudosuchians_EarlyNorian_withphyto.trimmed", "pseudosuchians_lateNorian.trimmed", "pseudosuchians_lateNorian_withphyto.trimmed", "avemetatarsalians_MidTr.trimmed", "avemetatarsalians_carnian.trimmed", "avemetatarsalians_earlyNorian.trimmed", "avemetatarsalians_lateNorian.trimmed"),
                               c("WMPD","Bootstrap_WMPD","Lower_bound","Upper_bound")))

for(i in 1:15) {mord.dist <- as.dist(trimmed.max.data$dist.matrix[taxongroups[[i]], taxongroups[[i]]])
                mord.comparable.char<- as.dist(dist.data$comp.char.matrix[taxongroups[[i]], taxongroups[[i]]])
                results[[i,1]]<-sum(mord.comparable.char * mord.dist) / sum(mord.comparable.char)
                }

# Bootstrapping
bootstrapped_distances <- list() # Create empty list to store output
resample_matrix <- function(nexus.data) { # Function to actually do the resampling:
  characters <- sample(1:ncol(nexus.data$matrix), replace = TRUE)
  nexus.data$matrix <- nexus.data$matrix[, characters]
  nexus.data$ordering <- nexus.data$ordering[characters]
  nexus.data$weights <- nexus.data$weights[characters]
  nexus.data$max.vals <- nexus.data$max.vals[characters]
  nexus.data$min.vals <- nexus.data$min.vals[characters]
  return(nexus.data)
}
for(i in 1:1000) bootstrapped_distances[[(length(bootstrapped_distances) + 1)]]<- MorphDistMatrixFast(resample_matrix(nexus.data))$max.dist.matrix
function(bootstrapped_distances){}


#Calculate WMPD for each bootstrapped matrix
for(i in 1:15) {
                mord.dist.boot <- list()
                for(y in 1:1000) {
                      mord.dist.boot [[(length(mord.dist.boot) + 1)]]<- as.dist(bootstrapped_distances[[y]][taxongroups[[i]], taxongroups[[i]]])
                                  }                
                
                mord.comparable.char<-as.dist(dist.data$comp.char.matrix[taxongroups[[i]], taxongroups[[i]]])
                WMPD<-c(1:1000)
                for(z in 1:1000) {
                  WMPD [z]<-sum(mord.comparable.char * mord.dist.boot[[z]]) / sum(mord.comparable.char)
                }  
                results[[i,2]]<-mean(sort(WMPD))
                #Calculate the boostrapped WMPD
                sorted.means<-sort(WMPD)
                #Calculate confidence intervals for the mean
                results[[i,3]]<-sorted.means[length(sorted.means)*0.025]
                results[[i,4]]<-sorted.means[length(sorted.means)*0.975]
          
}

write.csv(results,file="teleocrater_results.csv")
