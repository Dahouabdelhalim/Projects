# script to describe and assess DADs

load('DADs_data.RData')
load('DADs_metadata.RData')
load('DADSAD_data.RData')
ls() # dietary data is called 'diet', diet metadata is called 'metadata', data for environment-diet comparisons is called 'DADSAD'

library(stringr)


############ (1) Assigning body mass values to each consumer (note that these are already present in the metadata table so this step is unnecessary to repeat but reported here for clarity)

traits <- read.csv() # read in the Pantheria database (Jones et al. 2009 Ecology)

metadata$BodyMass <- traits$X5.1_AdultBodyMass_g[match(metadata$FocalSpecies, traits$MSW05_Binomial)]

unique(metadata$FocalSpecies[which(is.na(metadata$BodyMass)==TRUE & metadata$Class == 'Mammalia')])
metadata$BodyMass[which(metadata$FocalSpecies == 'Damaliscus_lunatus_topi')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Damaliscus_lunatus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Naemorhedus_goral_bedfordi')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Naemorhedus_goral')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Tragelaphus_scriptus_meneliki')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Tragelaphus_scriptus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Kobus_ellipsiprymnus_defassa')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Kobus_ellipsiprymnus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Odocoileus_virginianus_clavium')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Odocoileus_virginianus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Bos_indicus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Bos_taurus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Giraffa_camelopardalis_tippelskirchi')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Giraffa_camelopardalis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Cervus_canadensis')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Cervus_elaphus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Camelus_bactrianus_ferus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Camelus_bactrianus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Rangifer_tarandus_pearyi')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Rangifer_tarandus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Hydropotes_inermis_argyropus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Hydropotes_inermis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Damaliscus_dorcas')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Damaliscus_pygargus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Oryx_gazella_gazella')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Oryx_gazella')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Budorcas_taxicolor_whitei')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Budorcas_taxicolor')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Ovis_canadensis_canadensis')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Ovis_canadensis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Ovis_musimon')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Ovis_aries')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Redunca_fulvorufula_chnaleri')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Redunca_fulvorufula')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Bos_gaurus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Bos_frontalis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Cervus_timorensis_russa')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Rusa_timorensis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Cervus_unicolor')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Rusa_unicolor')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Poephagus_mutus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Bos_grunniens')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Tayassu_tajacu')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Pecari_tajacu')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Canis_familiaris')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Canis_lupus')] # because wikipedia gives it as lupus
metadata$BodyMass[which(metadata$FocalSpecies == 'Canis_familiaris_dingo')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Canis_lupus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Felis_pardalis')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Leopardus_pardalis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Felis_concolor')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Puma_concolor')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Herpestes_auropunctatus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Herpestes_javanicus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Panthera_uncia')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Uncia_uncia')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Felis_aurata')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Profelis_aurata')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Lagothrix_lagotricha_cana')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Lagothrix_cana')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Lagothrix_lagotricha_poeppigii')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Lagothrix_poeppigii')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Ateles_belzebuth_belzebuth')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Ateles_belzebuth')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Cacajao_calvus_calvus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Cacajao_calvus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Cacajao_calvus_ucayalii')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Cacajao_calvus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Lemur_rubriventer')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Eulemur_rubriventer')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Loris_lydekkerianus_lydekkerianus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Loris_lydekkerianus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Ateles_paniscus_paniscus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Ateles_paniscus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Erythrocebus_patas_pyrrhonotus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Erythrocebus_patas')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Trachypithecus_auratus_sondaicus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Trachypithecus_auratus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Macaca_fuscata_yakui')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Macaca_fuscata')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Cercopithecus_aethiops')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Chlorocebus_aethiops')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Cercopithecus_sabaeus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Chlorocebus_sabaeus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Cacajao_melanocephalus_melanocephalus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Cacajao_melanocephalus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Presbytis_pileata')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Trachypithecus_pileatus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Arctocephalus_pusillus_pusillus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Arctocephalus_pusillus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Arctocephalus_pusillus_doriferus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Arctocephalus_pusillus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Dipodomys_spectabilis_spectabilis')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Dipodomys_spectabilis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eptesicus_nilssoni')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Eptesicus_nilsoni')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eutomias_amoenus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Tamias_amoenus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Macaca_fuscata_yakui')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Macaca_fuscata')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Myocaster_coypus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Myocastor_coypus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Odobenus_rosmarus_rosmarus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Odobenus_rosmarus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Phoca_hispida')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Pusa_hispida')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Physeter_macrocephalus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Physeter_catodon')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Pipistrellus_mimus')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Pipistrellus_tenuis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Sciurus_ingrami')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Sciurus_aestuans')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Spermophilus_beecheyi_fisheri')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Spermophilus_beecheyi')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Cervus_albirostris')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Przewalskium_albirostris')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Cervus_duvauceli')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Rucervus_duvaucelii')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Capricornis_thar')] <- traits$X5.1_AdultBodyMass_g[which(traits$MSW05_Binomial == 'Capricornis_sumatraensis')] # genus level 
sort(unique(metadata$FocalSpecies[which(is.na(metadata$BodyMass) == TRUE & metadata$Class == 'Mammalia')]))

mam <- read.table() # read in EltonTraits database for mammals (Wilman et al. 2014 Ecology)
sort(unique(metadata$FocalSpecies[which(metadata$BodyMass == -999 & metadata$Class == 'Mammalia')]))
metadata$BodyMass[which(metadata$FocalSpecies == 'Chlorocebus_djamdjamensis')] <- mam$BodyMass.Value[which(mam$Scientific == 'Chlorocebus djamdjamensis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Conepatus_chinga')] <- mam$BodyMass.Value[which(mam$Scientific == 'Conepatus chinga')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Lagothrix_lagotricha_cana')] <- mam$BodyMass.Value[which(mam$Scientific == 'Lagothrix cana')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Lagothrix_lagotricha_poeppigii')] <- mam$BodyMass.Value[which(mam$Scientific == 'Lagothrix lagotricha')]
sort(unique(metadata$FocalSpecies[which(metadata$BodyMass == -999 | is.na(metadata$BodyMass) == TRUE & metadata$Class == 'Mammalia')]))

# add in birds
ave <- read.delim() # read in the EltonTraits database for birds (Wilman et al. 2014 Ecology)
ave$Scientific <- as.character(ave$Scientific)
ave$Scientific <- str_replace_all(ave$Scientific, ' ', '_')
for(i in 1:nrow(metadata)){
	if(metadata$Class[i] == 'Aves'){
		if(metadata$FocalSpecies[i] %in% ave$Scientific){
			metadata$BodyMass[i] <- ave$BodyMass.Value[which(ave$Scientific == metadata$FocalSpecies[i])]
		}
	}
}
sort(unique(metadata$FocalSpecies[which(is.na(metadata$BodyMass) == TRUE & metadata$Class == 'Aves')]))
metadata$BodyMass[which(metadata$FocalSpecies == 'Anas_discolor')] <- ave$BodyMass.Value[which(ave$Scientific == 'Anas_discors')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Asio_clamator')] <- ave$BodyMass.Value[which(ave$Scientific == 'Pseudoscops_clamator')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Anthropoides_virgo')] <- ave$BodyMass.Value[which(ave$Scientific == 'Grus_virgo')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Ara_macaro')] <- ave$BodyMass.Value[which(ave$Scientific == 'Ara_macao')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Bucorvus_leadbeateri')] <- ave$BodyMass.Value[which(ave$Scientific == 'Bucorvus_cafer')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Canachites_canadensis')] <- ave$BodyMass.Value[which(ave$Scientific == 'Dendragapus_canadensis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Diomedea_epomorpha_epomorpha')] <- ave$BodyMass.Value[which(ave$Scientific == 'Diomedea_epomophora')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Diomedea_epomorpha_sanfordi')] <- ave$BodyMass.Value[which(ave$Scientific == 'Diomedea_sanfordi')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Emberiza_calandra')] <- ave$BodyMass.Value[which(ave$Scientific == 'Miliaria_calandra')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Empidonax_traillii_extimus')] <- ave$BodyMass.Value[which(ave$Scientific == 'Empidonax_traillii')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Haliaeetus_leucopcephalus_alascanus')] <- ave$BodyMass.Value[which(ave$Scientific == 'Haliaeetus_leucocephalus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Hemiphaga_novaeseelandiae_chathamensis')] <- ave$BodyMass.Value[which(ave$Scientific == 'Hemiphaga_novaeseelandiae')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Hieraaetus_fasciatus')] <- ave$BodyMass.Value[which(ave$Scientific == 'Aquila_fasciatus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Lanius_meridionalis_koenigi')] <- ave$BodyMass.Value[which(ave$Scientific == 'Lanius_meridionalis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Musophaga_johnstoni')] <- ave$BodyMass.Value[which(ave$Scientific == 'Ruwenzorornis_johnstoni')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Onychoprion_fuscatus')] <- ave$BodyMass.Value[which(ave$Scientific == 'Sterna_fuscata')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Onychoprion_lunatus')] <- ave$BodyMass.Value[which(ave$Scientific == 'Sterna_lunata')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Phalacrocorax_auritis')] <- ave$BodyMass.Value[which(ave$Scientific == 'Phalacrocorax_auritus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Speotyto_cunicularia')] <- ave$BodyMass.Value[which(ave$Scientific == 'Athene_cunicularia')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Strix_occidentalis_lucida')] <- ave$BodyMass.Value[which(ave$Scientific == 'Strix_occidentalis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Tauraco_schuettii')] <- ave$BodyMass.Value[which(ave$Scientific == 'Tauraco_schuetti')]
sort(unique(metadata$FocalSpecies[which(is.na(metadata$BodyMass) == TRUE & metadata$Class == 'Aves')]))


# amphibians from AmphiBIO
amp <- read.csv() # read in AmphiBIO database (Oliveira et al. 2017 Scientific Data)
amp$Species <- as.character(amp$Species)
amp$Species <- str_replace_all(amp$Species, ' ', '_')
sort(unique(metadata$FocalSpecies[which(is.na(metadata$BodyMass) == TRUE & metadata$Class == 'Amphibia')]))
for(i in 1:nrow(metadata)){
	if(metadata$Class[i] == 'Amphibia'){
		if(metadata$FocalSpecies[i] %in% amp$Species){
			metadata$BodyMass[i] <- amp$Body_mass_g[which(amp$Species == metadata$FocalSpecies[i])]
		}
	}
}
notin <- sort(unique(metadata$FocalSpecies[which(is.na(metadata$BodyMass) == TRUE & metadata$Class == 'Amphibia')]))
notin[which(!notin %in% amp$Species)]
metadata$BodyMass[which(metadata$FocalSpecies == 'Bufo_coniferus')] <- amp$Body_mass_g[which(amp$Species == 'Incilius_coniferus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Bufo_haemititicus')] <- amp$Body_mass_g[which(amp$Species == 'Rhaebo_haematiticus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Bufo_typhonius')] <- amp$Body_mass_g[which(amp$Species == 'Rhinella_roqueana')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Colostethus_ingunialis')] <- amp$Body_mass_g[which(amp$Species == 'Colostethus_inguinalis')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Colostethus_nubicola')] <- amp$Body_mass_g[which(amp$Species == 'Silverstoneia_nubicola')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Colostethus_talamancae')] <- amp$Body_mass_g[which(amp$Species == 'Allobates_talamancae')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Dendrobates_fulguritus')] <- amp$Body_mass_g[which(amp$Species == 'Ranitomeya_fulgurita')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Dendrobates_minutus')] <- amp$Body_mass_g[which(amp$Species == 'Ranitomeya_minuta')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Desmognathus_fuscus_fuscus')] <- amp$Body_mass_g[which(amp$Species == 'Desmognathus_fuscus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_altae')] <- amp$Body_mass_g[which(amp$Species == 'Pristimantis_altae')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_biporcatus')] <- amp$Body_mass_g[which(amp$Species == 'Strabomantis_biporcatus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_bransfordi')] <- amp$Body_mass_g[which(amp$Species == 'Craugastor_bransfordii')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_cruentus')] <- amp$Body_mass_g[which(amp$Species == 'Pristimantis_cruentus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_frater')] <- amp$Body_mass_g[which(amp$Species == 'Pristimantis_frater')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_golmeri')] <- amp$Body_mass_g[which(amp$Species == 'Craugastor_gollmeri')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_longirostris')] <- amp$Body_mass_g[which(amp$Species == 'Craugastor_longirostris')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_talamancae')] <- amp$Body_mass_g[which(amp$Species == 'Craugastor_talamancae')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_vocator')] <- amp$Body_mass_g[which(amp$Species == 'Diasporus_vocator')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Leptodactylus_ocellatus')] <- amp$Body_mass_g[which(amp$Species == 'Leptodactylus_latrans')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Rana_berlandieri')] <- amp$Body_mass_g[which(amp$Species == 'Lithobates_berlandieri')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Rana_clamitans')] <- amp$Body_mass_g[which(amp$Species == 'Lithobates_clamitans')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Rana_grylio')] <- amp$Body_mass_g[which(amp$Species == 'Lithobates_grylio')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Rana_nigromaculata')] <- amp$Body_mass_g[which(amp$Species == 'Pelophylax_nigromaculatus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Rana_porosa_brevipoda')] <- amp$Body_mass_g[which(amp$Species == 'Pelophylax_porosus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Rana_ridibunda')] <- amp$Body_mass_g[which(amp$Species == 'Pelophylax_ridibundus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Schoutedenella_xenodactyloides')] <- amp$Body_mass_g[which(amp$Species == 'Arthroleptis_xenodactyloides')]

# lizards
rept <- read.csv() # read in amniote database (Myhrvold et al. 2015 Ecology; https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/15-0846R.1)
rept$Bin <- paste(rept$genus, rept$species, sep='_')
rept2 <- rept[which(rept$class == 'Reptilia' & rept$adult_body_mass_g != -999),]

sort(unique(metadata$FocalSpecies[which(is.na(metadata$BodyMass) == TRUE & metadata$Class == 'Reptilia')]))
for(i in 1:nrow(metadata)){
	if(metadata$Class[i] == 'Reptilia'){
		if(metadata$FocalSpecies[i] %in% rept2$Bin){
			metadata$BodyMass[i] <- rept2$adult_body_mass_g[which(rept2$Bin == metadata$FocalSpecies[i])]
		}
	}
}
sort(unique(metadata$FocalSpecies[which(is.na(metadata$BodyMass) == TRUE & metadata$Class == 'Reptilia')]))
metadata$BodyMass[which(metadata$FocalSpecies == 'Cnemidophorus_sexlineatus')] <- rept$adult_body_mass_g[which(rept$Bin == 'Aspidoscelis_sexlineata')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Crotalus_lepidus_klauberi')] <- rept$adult_body_mass_g[which(rept$Bin == 'Crotalus_lepidus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Leiocephalus_armouri')] <- rept$adult_body_mass_g[which(rept$Bin == 'Leiocephalus_carinatus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Leiocephalus_virescens')] <- rept$adult_body_mass_g[which(rept$Bin == 'Leiocephalus_carinatus')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Morelia_spilota_spilota')] <- rept$adult_body_mass_g[which(rept$Bin == 'Morelia_spilota')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Nerodia_harteri_paucimaculata')] <- rept$adult_body_mass_g[which(rept$Bin == 'Nerodia_harteri')]
metadata$BodyMass[which(metadata$FocalSpecies == 'Python_reticulatus')] <- rept$adult_body_mass_g[which(rept$Bin == 'Malayopython_reticulatus')]
sort(unique(metadata$FocalSpecies[which(is.na(metadata$BodyMass) == TRUE & metadata$Class == 'Reptilia')]))

length(which(metadata$BodyMass == -999 | is.na(metadata$BodyMass) == TRUE))
unique(metadata$FocalSpecies[which(metadata$BodyMass == -999 | is.na(metadata$BodyMass) == TRUE)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Morone_mississippiensis')] <- mean(c(6050,11000)) #https://eol.org/pages/207906/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Thuunus_maccoyii')] <- mean(c(143000, 260000)) # https://eol.org/pages/46577338/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Lepidogobius_lepidus')] <- 14.95 # https://bioone.org/journals/northwest-science/volume-88/issue-2/046.088.0208/The-Effect-of-Region-Body-Size-and-Sample-Size-on/10.3955/046.088.0208.full
metadata$BodyMass[which(metadata$FocalSpecies == 'Lophius_piscatorius')] <- mean(c(31735, 57700)) # https://eol.org/pages/46566106/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Rhomboplites_aurorubens')] <- mean(c(3170, 1743.5)) # https://eol.org/pages/46580805
metadata$BodyMass[which(metadata$FocalSpecies == 'Pterois_volitans')] <- 1442 # https://eol.org/pages/46568062
metadata$BodyMass[which(metadata$FocalSpecies == 'Galeorhinus_galeus')] <- mean(c(24585, 44700)) # https://eol.org/pages/46559964/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Mustelus_asterias')] <- 4760 # https://eol.org/pages/46559983/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Scyliorhinus_stellaris')] <- mean(c(mean(c(1320,726)))) # represents genus average https://eol.org/search?utf8=%E2%9C%93&q=scyliorhinus
metadata$BodyMass[which(metadata$FocalSpecies == 'Squatina_squatina')] <- 80000 # https://eol.org/pages/46560323/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Raja_brachyura')] <- mean(c(14300, 7865)) # https://eol.org/pages/46560573/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Raja_clavata')] <- mean(c(18000, 9900)) # https://eol.org/pages/46560547/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Scyliorhinus_canicula')] <- mean(c(1320,726)) # https://eol.org/pages/46559871/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Squalus_acanthias')] <- mean(c(9100, 5005)) # https://eol.org/pages/46560201/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Solea_solea')] <- mean(c(3000, 1650)) # https://eol.org/pages/46570285/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Solea_senegalensis')] <- mean(c(mean(c(3000,1650)))) # genus average https://eol.org/search?utf8=%E2%9C%93&q=solea+
metadata$BodyMass[which(metadata$FocalSpecies == 'Coryphaena_hippurus')] <- mean(c(40000, 22000)) # https://eol.org/pages/46578341/data
metadata$BodyMass[which(metadata$FocalSpecies == 'Carcharodon_carcharias')] <- mean(c(1870000, 3230000)) # https://eol.org/pages/46559751/data


length(which(!metadata$BodyMass == -999 | is.na(metadata$BodyMass) == FALSE))
unique(metadata$FocalSpecies[which(metadata$BodyMass == -999 | is.na(metadata$BodyMass) == TRUE)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Atelopus_varius')] <- mean(amp$Body_mass_g[grep('Atelopus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Bufo_coniferus')] <- mean(amp$Body_mass_g[grep('Bufo', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Bufo_haemititicus')] <- mean(amp$Body_mass_g[grep('Bufo', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Bufo_typhonius')] <- mean(amp$Body_mass_g[grep('Bufo', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_altae')] <- mean(amp$Body_mass_g[grep('Eleutherodactylus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_biporcatus')] <- mean(amp$Body_mass_g[grep('Eleutherodactylus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_cruentus')] <- mean(amp$Body_mass_g[grep('Eleutherodactylus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_frater')] <- mean(amp$Body_mass_g[grep('Eleutherodactylus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_longirostris')] <- mean(amp$Body_mass_g[grep('Eleutherodactylus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_talamancae')] <- mean(amp$Body_mass_g[grep('Eleutherodactylus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_vocator')] <- mean(amp$Body_mass_g[grep('Eleutherodactylus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Leiocephalus_loxogrammus')] <- mean(rept$adult_body_mass_g[grep('Leiocephalus', rept$Bin)][which(rept$adult_body_mass_g[grep('Leiocephalus', rept$Bin)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Leiocephalus_inaguae')] <- mean(rept$adult_body_mass_g[grep('Leiocephalus', rept$Bin)][which(rept$adult_body_mass_g[grep('Leiocephalus', rept$Bin)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Leiocephalus_greenwayi')] <- mean(rept$adult_body_mass_g[grep('Leiocephalus', rept$Bin)][which(rept$adult_body_mass_g[grep('Leiocephalus', rept$Bin)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Raja_montagui')] <- mean(c(mean(c(14300, 7865)), mean(c(18000,9900)), mean(c(4500,2475)), mean(c(97100,53405)))) # genus average https://eol.org/search?utf8=%E2%9C%93&q=raja++
metadata$BodyMass[which(metadata$FocalSpecies == 'Raja_naevus')] <- mean(c(mean(c(14300, 7865)), mean(c(18000,9900)), mean(c(4500,2475)), mean(c(97100,53405))))  # genus average https://eol.org/search?utf8=%E2%9C%93&q=raja++
metadata$BodyMass[which(metadata$FocalSpecies == 'Desmognathus_fuscus_fuscus')] <- mean(amp$Body_mass_g[grep('Desmognathus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Cercopithecus_sabaeus')] <- mean(traits$X5.1_AdultBodyMass_g[grep('Chlorocebus', traits$MSW05_Genus)][which(traits$X5.1_AdultBodyMass_g[grep('Chlorocebus', traits$MSW05_Genus)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Trachypithecus_delacouri')] <- mean(traits$X5.1_AdultBodyMass_g[grep('Trachypithecus', traits$MSW05_Genus)][which(traits$X5.1_AdultBodyMass_g[grep('Trachypithecus', traits$MSW05_Genus)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Psammophis_schokari')] <- mean(rept$adult_body_mass_g[grep('Psammophis', rept$Bin)][which(rept$adult_body_mass_g[grep('Psammophis', rept$Bin)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Helicops_infrataeniatus')] <- mean(rept$adult_body_mass_g[grep('Helicops', rept$Bin)][which(rept$adult_body_mass_g[grep('Helicops', rept$Bin)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Laticauda_semifasciata')] <- mean(rept$adult_body_mass_g[grep('Laticauda', rept$Bin)][which(rept$adult_body_mass_g[grep('Laticauda', rept$Bin)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Schoutedenella_xenodactyloides')] <- mean(amp$Body_mass_g[grep('Arthroleptis', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_coqui')] <- mean(amp$Body_mass_g[grep('Eleutherodactylus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Rana_grylio')] <- mean(amp$Body_mass_g[grep('Rana', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Rana_berlandieri')] <- mean(amp$Body_mass_g[grep('Rana', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Lama_guanicoe')] <- mean(traits$X5.1_AdultBodyMass_g[grep('Lama', traits$MSW05_Genus)][which(traits$X5.1_AdultBodyMass_g[grep('Lama', traits$MSW05_Genus)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_bransfordi')] <- mean(amp$Body_mass_g[grep('Eleutherodactylus', amp$Species)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Eleutherodactylus_golmeri')] <- mean(amp$Body_mass_g[grep('Eleutherodactylus', amp$Species)], na.rm=TRUE)

length(which(!metadata$BodyMass == -999 | is.na(metadata$BodyMass) == FALSE))
length(which(metadata$BodyMass == -999 | is.na(metadata$BodyMass) == TRUE))
unique(metadata$FocalSpecies[which(metadata$BodyMass == -999 | is.na(metadata$BodyMass) == TRUE)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Colostethus_ingunialis')] <- mean(amp$Body_mass_g[grep('Dendrobatidae', amp$Family)], na.rm=TRUE) 
metadata$BodyMass[which(metadata$FocalSpecies == 'Colostethus_nubicola')] <- mean(amp$Body_mass_g[grep('Dendrobatidae', amp$Family)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Colostethus_pratti')] <- mean(amp$Body_mass_g[grep('Dendrobatidae', amp$Family)], na.rm=TRUE) 
metadata$BodyMass[which(metadata$FocalSpecies == 'Colostethus_talamancae')] <- mean(amp$Body_mass_g[grep('Dendrobatidae', amp$Family)], na.rm=TRUE) 
metadata$BodyMass[which(metadata$FocalSpecies == 'Dendrobates_auratus')] <- mean(amp$Body_mass_g[grep('Dendrobatidae', amp$Family)], na.rm=TRUE) 
metadata$BodyMass[which(metadata$FocalSpecies == 'Dendrobates_fulguritus')] <- mean(amp$Body_mass_g[grep('Dendrobatidae', amp$Family)], na.rm=TRUE) 
metadata$BodyMass[which(metadata$FocalSpecies == 'Dendrobates_minutus')] <- mean(amp$Body_mass_g[grep('Dendrobatidae', amp$Family)], na.rm=TRUE) 
metadata$BodyMass[which(metadata$FocalSpecies == 'Eurycea_bislineata')] <- mean(amp$Body_mass_g[grep('Plethodontidae', amp$Family)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Gyrinophilus_porphyriticus')] <- mean(amp$Body_mass_g[grep('Plethodontidae', amp$Family)], na.rm=TRUE)
metadata$BodyMass[which(metadata$FocalSpecies == 'Uromacer_frenatus')] <- mean(rept$adult_body_mass_g[grep('Colubridae', rept$family)][which(rept$adult_body_mass_g[grep('Colubridae', rept$family)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Afronatrix_anoscopus')] <- mean(rept$adult_body_mass_g[grep('Colubridae', rept$family)][which(rept$adult_body_mass_g[grep('Colubridae', rept$family)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Antillopis_parvifrons')] <- mean(rept$adult_body_mass_g[grep('Colubridae', rept$family)][which(rept$adult_body_mass_g[grep('Colubridae', rept$family)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Symphimus_mayae')] <- mean(rept$adult_body_mass_g[grep('Colubridae', rept$family)][which(rept$adult_body_mass_g[grep('Colubridae', rept$family)] > 0)])
metadata$BodyMass[which(metadata$FocalSpecies == 'Ranitomeya_virolinensis')] <- mean(amp$Body_mass_g[grep('Dendrobatidae', amp$Family)], na.rm=TRUE) 
unique(metadata$FocalSpecies[which(metadata$BodyMass == -999 | is.na(metadata$BodyMass) == TRUE)])



########## (2) Assessing DAD shape

# data is the vector of diet contributions for foods that a population eats, aicweights is a logical (if TRUE, the aicweight of each model is returned alongside AICc)
DAD_shape <- function(data, aicweights){

	df <- data.frame(x=1:length(data), y=rev(sort(data/100)))

	if(nrow(df) > 3){
		ones <- rep(1, nrow(df)) # to fit horizontal line to data with nlsLM, which is used here rather than lm to ensure
		# consistency in fitting approach across DAD shapes
		fit_hor <- try(minpack.lm::nlsLM(y ~ a * ones, data=df, start=c(a=0.1)), silent=TRUE)
		if(class(fit_hor) != 'try-error'){
			Equal=AICcmodavg::AICc(fit_hor)
			}else{
				Equal=NA
			}

		if(nrow(df) > 4){
		fit_lin <- try(minpack.lm::nlsLM(y ~ a + b*x, data=df, start=c(a=0.5, b=-0.1), upper=c(Inf,Inf), lower=c(0,-Inf)), silent=TRUE)	
		if(class(fit_lin) != 'try-error'){
			LinHeir=AICcmodavg::AICc(fit_lin)
			}else{
				LinHeir=NA
			}


		if(nrow(df) > 5){
			fit_up <- try(minpack.lm::nlsLM(y ~ (c ^(x-d))+f, data=df, start=c(c=0.1, f=0.1, d=0.1), 
				control=minpack.lm::nls.lm.control(maxfev=10000, maxiter=1024), upper=c(1,Inf,Inf), lower=c(0,-Inf,-Inf)), silent=TRUE)
			if(class(fit_up) != 'try-error'){
				FewCManyRar=AICcmodavg::AICc(fit_up)
			}else{
				FewCManyRar=NA
			}

			fit_down <- try(minpack.lm::nlsLM(y ~ -((((x-1)^g)/h)-b), data=df, start=c(h= 10, i=0.5, g=2), 
				control=minpack.lm::nls.lm.control(maxfev=10000, maxiter=1024), upper=c(Inf,1,Inf), lower=c(0,0,2)), silent=TRUE)
			if(class(fit_down) != 'try-error'){
				WeakHeir=AICcmodavg::AICc(fit_down)
			}else{
				WeakHeir=NA
			}

			if(nrow(df) > 6){
				fit_log <- try(minpack.lm::nlsLM(y ~ (j/(1 + exp(k*(x-m))))+n, data=df, start=c(j=0, k=2, m=4,n=0), 
					control=minpack.lm::nls.lm.control(maxfev=10000, maxiter=1024), upper=c(1,Inf,Inf,1), lower=c(0,1,4,0)), silent=TRUE)
				if(class(fit_log) != 'try-error'){
					StapSupp=AICcmodavg::AICc(fit_log)
				}else{
					StapSupp=NA
				}
			}else{
				StapSupp=NA
			}	
		}else{
			FewCManyRar=NA
			WeakHeir=NA
			StapSupp=NA
		}
	}else{
		LinHeir=NA
		FewCManyRar=NA
		WeakHeir=NA
		StapSupp=NA
	}
	}else{
		Equal=NA
		LinHeir=NA
		FewCManyRar=NA
		WeakHeir=NA
		StapSupp=NA
	}

	if(aicweights == TRUE){
		aiccs <- c(Equal, LinHeir, FewCManyRar, WeakHeir, StapSupp)
		names(aiccs) <- c('Equal', 'LinHeir', 'FewCManyRar', 'WeakHeir', 'StapSupp')
		if(all(is.na(aiccs)) == TRUE){
			res <- data.frame(Equal, LinHeir, FewCManyRar, WeakHeir, StapSupp, NA, NA, NA, NA, NA)
		}else{
			if(any(is.na(aiccs)) == TRUE){
				aiccs <- aiccs[-which(is.na(aiccs) == TRUE)]
			}
			delAIC <- aiccs - min(aiccs)
			relLik <- exp(-0.5 * delAIC)
			aw <- relLik/sum(relLik)
			if(length(aw) != 5){
				mods <- c('Equal', 'LinHeir', 'FewCManyRar', 'WeakHeir', 'StapSupp')
				out <- mods[which(!mods %in% names(aw))]
				for(i in 1:length(out)){
					aw <- c(aw, NA)
					names(aw)[length(aw)] <- out[i]
				}
			}
			res <- data.frame(Equal, LinHeir, FewCManyRar, WeakHeir, StapSupp, aw[which(names(aw) == 'Equal')],
				aw[which(names(aw) == 'LinHeir')],aw[which(names(aw) == 'FewCManyRar')],aw[which(names(aw) == 'WeakHeir')],
				aw[which(names(aw) == 'StapSupp')])
		}
	}else{
		res <- data.frame(Equal, LinHeir, FewCManyRar, WeakHeir, StapSupp)
	}
	
	return(res)
}

# create data.frame for results
res <- data.frame(Diet=names(diet), Equal=NA, LinHeir=NA, FewCManyRar=NA, WeakHeir=NA, StapSupp=NA, 
	awEqual=NA, awLinHeir=NA, awFewCManyRar=NA, awWeakHeir=NA, awStapSupp=NA)
for(i in 1:length(diet)){
	g <- DAD_shape(diet[[i]], aicweights=TRUE)
	res[i,2:ncol(res)] <- g
}

	############ trim results to the diets that have enough food items to be fit 

	# first remove the ones too small to fit (because they have 1-3 food items)
not_fit <- res[which(apply(res[,2:6], MAR=1, FUN=function(x) all(is.na(x) == TRUE)) == TRUE),]
sapply(diet[which(apply(res[,2:6], MAR=1, FUN=function(x) all(is.na(x) == TRUE)) == TRUE)], length)
res2 <- res[-which(apply(res[,2:6], MAR=1, FUN=function(x) all(is.na(x) == TRUE)) == TRUE),]

	# assess which ones were only fit to one (because they only have four food items)
only_horizontal <- res2[which(apply(res2[,2:6], MAR=1, FUN=function(x) length(which(is.na(x) == TRUE))) == 4),]
sapply(diet[which(names(diet) %in% as.character(res2[which(apply(res2[,2:6], MAR=1, FUN=function(x) length(which(is.na(x) == TRUE))) == 4),1]))], length)
res3 <- res2[-which(apply(res2[,2:6], MAR=1, FUN=function(x) length(which(is.na(x) == TRUE))) == 4),]

nrow(res3) # n = 1130 diets remain
best <- apply(res3[,2:6], MAR=1, FUN=function(x) names(x)[which(x == min(x, na.rm=TRUE))]) # what was the best shape for each?
within_two <- apply(res3[,2:6], MAR=1, FUN=function(x) names(x)[which(x-2 < min(x, na.rm=TRUE))]) # what shapes were within 2 of the minimum AICc?
													
length(within_two[which(sapply(within_two, length) == 1)]) # how many diets had no other shapes within 2 of the minimum AICc
table(unlist(within_two[which(sapply(within_two, length) == 1)]))/length(within_two[which(sapply(within_two, length) == 1)]) # whats the breakdown by shape?






########## (3) Measuring sp50

# dat is a vector describing the diet contributions of foods to a diet, level is the cumulative diet proportion to estimate generalization at (e.g, 0.5 for sp50, 0.75 for sp75)
sp50_fun <- function(dat, level=0.5){
	dat <- rev(sort(dat)) # order diet items by percent abundance
	dat <- dat/100 # convert to proportional abundance
	d2 <- c(0, dat[1]) # the next five lines create a data frame that begins at 0,0 adding one food per row and the cumulative contribution of all preceding foods
	for(i in 2:length(dat)){
		d2 <- c(d2, sum(dat[1:i]))
	}
	df <- data.frame(x=0:length(dat), y=d2)
	sides <- df[c(which(df$y == max(df$y[which(df$y < level)])), which(df$y > level)[1]),] # isolate the two foods directly to the sides of 'level'
	fit <- lm(y ~ x, sides) # fit a line to those two points
	sp50 <- (level-coef(fit)[[1]])/coef(fit)[[2]] # using the line and desired level, solve for x (the number of foods required to get to the desired dietary proportion)
	return(sp50)
}

all(metadata$DietID == names(diet)) # confirm the data and metadata are in the same order
metadata$sp50 <- sapply(diet, function(x) sp50_fun(x, level=0.5)) # note that transformation of data from percent to proportional contributions is executed within sp50_fun

summary(metadata$sp50)
summary(metadata$sp50[which(metadata$SampleSize >= 20 & metadata$IntraAnnualExtent >= 180)])

########## (4) Evaluate models to predict sp50
library(glmulti)

	# first limit data to the diets reporting all variables
dat <- metadata[which(is.na(metadata$Latitude) == FALSE & 
	is.na(metadata$BodyMass) == FALSE &
	is.na(metadata$Environment) == FALSE &
	is.na(metadata$ProportionDietFruit) == FALSE & 
	is.na(metadata$ConsumerType) == FALSE &
	is.na(metadata$SampleSize) == FALSE &
	is.na(metadata$IntraAnnualExtent) == FALSE &
	is.na(metadata$InterAnnualExtent) == FALSE &
	is.na(metadata$TaxonomicResolution) == FALSE &
	metadata$SamplingMethod %in% c('Fecal', 'Stomach', 'Observation')),]

	# prepare variables for model
dat$FruitIn <- 0
dat$FruitIn[which(dat$ProportionDietFruit > 0)] <- 1
dat$BodyMass <- scale(log10(dat$BodyMass))
dat$Latitude <- scale(abs(dat$Latitude))
dat$SampleSize <- scale(log10(dat$SampleSize))
dat$InterAnnualExtent <- scale(log10(dat$InterAnnualExtent))
dat$IntraAnnualExtent <- scale(dat$IntraAnnualExtent)
dat$TaxonomicResolution <- scale(dat$TaxonomicResolution)
colnames(dat)[which(colnames(dat) == 'sp50')] <- 'Y'
dat <- dat[,which(colnames(dat) %in% c('Y', 'FruitIn', 'ConsumerType', 'BodyMass', 'Environment', 'Latitude', 
	'SampleSize', 'TaxonomicResolution', 'SamplingMethod', 'InterAnnualExtent', 'IntraAnnualExtent'))]
dat$Y <- log10(dat$Y)
nrow(dat)

	# generate the set of additive models with an intercept (n = 1024)
res <- glmulti(Y ~ FruitIn + factor(ConsumerType) + BodyMass + Latitude + factor(Environment) + SampleSize + TaxonomicResolution + factor(SamplingMethod) +
	IntraAnnualExtent + InterAnnualExtent, intercept=TRUE, level=1, data=dat, confsetsize=1100,
	fitfunction=lm, crit='aicc')
plot(res)
print(res)

	# isolate the models with the most support
top <- weightable(res)
top <- top[top$aicc <= min(top$aicc) + 2,]
top
summary(res@objects[[1]])
summary(res@objects[[which(sapply(res@objects, function(x) length(coef(x))) == max(sapply(res@objects, function(x) length(coef(x)))))]]) # with all vars

	# investigate the relative importance of each variable
plot(res, type="s")
x <- res
ww = exp(-(x@crits - x@crits[1])/2)
ww = ww/sum(ww)
clartou = function(x) {
    pieces <- sort(strsplit(x, ":")[[1]])
    if (length(pieces) > 1) 
        paste(pieces[1], ":", pieces[2], sep = "")
    else x
}
tet = lapply(x@formulas, function(x) sapply(attr(delete.response(terms(x)), "term.labels"), clartou))
allt <- unique(unlist(tet))
imp <- sapply(allt, function(x) sum(ww[sapply(tet, function(t) x %in% t)]))
imp



########### (5) Case studies of DAD-SAD relationships


	# compare SAD and DAD shapes; note p-values are not corrected for multiple comparisons (as they are in the manuscript with p.adjust and method 'BH')
names(DADSAD)
ks.test(rev(sort(DADSAD[[1]]$Environment)), rev(sort(DADSAD[[1]]$Diet)))
ks.test(rev(sort(DADSAD[[2]]$Environment)), rev(sort(DADSAD[[2]]$Diet)))
ks.test(rev(sort(DADSAD[[3]]$Environment)), rev(sort(DADSAD[[3]]$Diet)))
ks.test(rev(sort(DADSAD[[4]]$Environment)), rev(sort(DADSAD[[4]]$Diet)))
ks.test(rev(sort(DADSAD[[5]]$Environment)), rev(sort(DADSAD[[5]]$Diet)))
ks.test(rev(sort(DADSAD[[6]]$Environment)), rev(sort(DADSAD[[6]]$Diet)))
	

	# correlation of ranks
for(i in 1:length(DADSAD)){
	d <- DADSAD[[i]]
	d$EnvRank <- rank(-d$Environment, ties.method='random') # rank availability from most common to most rare
	d$DietRank <- rank(-d$Diet, ties.method='random') # rank diet contributions from most common to most rare
	corel <- cor.test(d$EnvRank, d$DietRank, method='spearman') # what is the correlation between rank abundance and rank diet contribution

	if(any(duplicated(d$Diet) == TRUE | duplicated(d$Environment) == TRUE)){ # when foods have equal proportions, repeat process randomly breaking ties each time to get a mean correlation
		set.seed(1984)
		corel <- list(rho=corel$estimate, p.value=corel$p.value)
		for(x in 1:9){
			d$EnvRank <- rank(-d$Environment, ties.method='random') 
			d$DietRank <- rank(-d$Diet, ties.method='random')
			cor <- cor.test(d$EnvRank, d$DietRank, method='spearman') 
			corel$rho <- c(corel$rho, cor$estimate)
			corel$p.value <- c(corel$p.value, cor$p.value)
		}
		print(names(DADSAD)[i])
		print(mean(corel$rho)) # because of the stochastic nature of sampling 10 different tie-breaks, the mean coefficient will likely differ quantitatively but not qualitatively if the loop is run more than once
		print(mean(corel$p.value))
	}else{
		print(names(DADSAD)[i])
		print(corel)
	}

}


	# diets composition compared to random sampling

	# sad is vector of proportional availability, n is number of food selections to make
rand_diet <- function(sad, n){
	diets <- matrix(0, ncol=length(sad), nrow=1)
	colnames(diets) <- names(sad)
	x <- table(sample(names(sad), size=n, replace=TRUE, prob=sad))
	ranab <- x/sum(x)
	diets[,match(names(ranab), colnames(diets))] <- ranab
	return(diets)
}

	# for each DAD-SAD comparison, compute the dissimilarity between diet and environment, then sample random diets based on availability and compare observed dissimilarity to random-diet dissimilarity with environment 
for(i in 1:length(DADSAD)){	
	obs <- vegan::vegdist(t(DADSAD[[i]][,c('Environment','Diet')]), binary=FALSE, method='jaccard')
	sad <- DADSAD[[i]]$Environment
	names(sad) <- DADSAD[[i]]$Resource

	n=1000
	null_jac <- c()
	for(x in 1:n){
		null_diet <- rand_diet(sad, n=100)
		null_jac <- c(null_jac, vegan::vegdist(rbind(t(DADSAD[[i]]$Environment), null_diet), binary=FALSE, method='jaccard'))
	}

	p <- (length(which(null_jac >= obs))+1)/(n+1) # compute Monte Carlo p-value for the hypothesis that observed diet-environment dissimilarity is greater than the null distribution of diet-environment dissimilarity
	print(names(DADSAD)[i])
	print(p)
}