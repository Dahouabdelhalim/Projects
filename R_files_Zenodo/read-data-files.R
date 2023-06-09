## Copyright (c) 2018 Frédéric Vanwindekens, Dóra Drexler
## Centre wallon de Recherches agronomiques (CRA-W)
## Hungarian Research Institute of Organic Agriculture (ÖMKi),
## Projet DiverIMPACTS (https://www.diverimpacts.net/)
## License cc-by

data <- read.csv("./survey_438634_R_data_file.csv",
                 quote = "'\\"",
                 na.strings = c("", "\\"\\""),
                 stringsAsFactors = FALSE,
                 skip = 5) ## the five-lines-license

# LimeSurvey Field type: F
data[, 1] <- as.numeric(data[, 1])
attributes(data)$variable.labels[1] <- "id"
names(data)[1] <- "id"
# LimeSurvey Field type: DATETIME23.2
data[, 2] <- as.character(data[, 2])
attributes(data)$variable.labels[2] <- "submitdate"
names(data)[2] <- "submitdate"
# LimeSurvey Field type: 
data[, 3] <- as.character(data[, 3])
attributes(data)$variable.labels[3] <- "lastpage"
names(data)[3] <- "lastpage"
# LimeSurvey Field type: A
data[, 4] <- as.character(data[, 4])
attributes(data)$variable.labels[4] <- "startlanguage"
names(data)[4] <- "startlanguage"
# LimeSurvey Field type: A
data[, 5] <- as.character(data[, 5])
attributes(data)$variable.labels[5] <- "token"
names(data)[5] <- "token"
# LimeSurvey Field type: DATETIME23.2
data[, 6] <- as.character(data[, 6])
attributes(data)$variable.labels[6] <- "startdate"
names(data)[6] <- "startdate"
# LimeSurvey Field type: DATETIME23.2
data[, 7] <- as.character(data[, 7])
attributes(data)$variable.labels[7] <- "datestamp"
names(data)[7] <- "datestamp"
# LimeSurvey Field type: A
data[, 8] <- as.character(data[, 8])
attributes(data)$variable.labels[8] <- "ipaddr"
names(data)[8] <- "ipaddr"
# LimeSurvey Field type: A
data[, 9] <- as.character(data[, 9])
attributes(data)$variable.labels[9] <- "1.1 Title/Name of diversification initiative"
names(data)[9] <- "A1"
# LimeSurvey Field type: A
data[, 10] <- as.character(data[, 10])
attributes(data)$variable.labels[10] <- "[Country] 1.2 Geographical location of diversification initiative"
names(data)[10] <- "A2_SQ001"
# LimeSurvey Field type: A
data[, 11] <- as.character(data[, 11])
attributes(data)$variable.labels[11] <- "[Region(s)] 1.2 Geographical location of diversification initiative"
names(data)[11] <- "A2_SQ002"
# LimeSurvey Field type: F
data[, 12] <- as.numeric(data[, 12])
attributes(data)$variable.labels[12] <- "[Rotation was applied (different crops in successive growing years) Please describe the typical plant order.] 1.3 Please describe the production practice before the diversification initiative. Use the comment box to complete your answer."
data[, 12] <- factor(data[, 12], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[12] <- "A3_SQ001"
# LimeSurvey Field type: A
data[, 13] <- as.character(data[, 13])
attributes(data)$variable.labels[13] <- "[Comment] 1.3 Please describe the production practice before the diversification initiative. Use the comment box to complete your answer."
names(data)[13] <- "A3_SQ001comment"
# LimeSurvey Field type: F
data[, 14] <- as.numeric(data[, 14])
attributes(data)$variable.labels[14] <- "[Multicropping was applied (different crops within one growing year) Please describe the typical plant order.] 1.3 Please describe the production practice before the diversification initiative. Use the comment box to complete your answer."
data[, 14] <- factor(data[, 14], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[14] <- "A3_SQ002"
# LimeSurvey Field type: A
data[, 15] <- as.character(data[, 15])
attributes(data)$variable.labels[15] <- "[Comment] 1.3 Please describe the production practice before the diversification initiative. Use the comment box to complete your answer."
names(data)[15] <- "A3_SQ002comment"
# LimeSurvey Field type: F
data[, 16] <- as.numeric(data[, 16])
attributes(data)$variable.labels[16] <- "[Intercropping was applied (growing different species or cultivars in proximity in the same field) Please describe the typical intercropping (e.g. mixed, row or strip intercropping).] 1.3 Please describe the production practice before the diversification initiative. Use the comment box to complete your answer."
data[, 16] <- factor(data[, 16], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[16] <- "A3_SQ003"
# LimeSurvey Field type: A
data[, 17] <- as.character(data[, 17])
attributes(data)$variable.labels[17] <- "[Comment] 1.3 Please describe the production practice before the diversification initiative. Use the comment box to complete your answer."
names(data)[17] <- "A3_SQ003comment"
# LimeSurvey Field type: F
data[, 18] <- as.numeric(data[, 18])
attributes(data)$variable.labels[18] <- "[None of these were applied] 1.3 Please describe the production practice before the diversification initiative. Use the comment box to complete your answer."
data[, 18] <- factor(data[, 18], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[18] <- "A3_SQ004"
# LimeSurvey Field type: A
data[, 19] <- as.character(data[, 19])
attributes(data)$variable.labels[19] <- "[Comment] 1.3 Please describe the production practice before the diversification initiative. Use the comment box to complete your answer."
names(data)[19] <- "A3_SQ004comment"
# LimeSurvey Field type: F
data[, 20] <- as.numeric(data[, 20])
attributes(data)$variable.labels[20] <- "[New cash or cover crop(s) were included in the rotation (these could be annual or perennial)] 	1.4 What is/was done? Please indicate in the comment box the starting date of the diversification initiative, and where possible, the end (Year). Please add this information before selecting the next relevant answer."
data[, 20] <- factor(data[, 20], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[20] <- "A4_SQ009"
# LimeSurvey Field type: A
data[, 21] <- as.character(data[, 21])
attributes(data)$variable.labels[21] <- "[Comment] 	1.4 What is/was done? Please indicate in the comment box the starting date of the diversification initiative, and where possible, the end (Year). Please add this information before selecting the next relevant answer."
names(data)[21] <- "A4_SQ009comment"
# LimeSurvey Field type: F
data[, 22] <- as.numeric(data[, 22])
attributes(data)$variable.labels[22] <- "[New cash or cover crop(s) were included within one season (multicropping)] 	1.4 What is/was done? Please indicate in the comment box the starting date of the diversification initiative, and where possible, the end (Year). Please add this information before selecting the next relevant answer."
data[, 22] <- factor(data[, 22], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[22] <- "A4_SQ010"
# LimeSurvey Field type: A
data[, 23] <- as.character(data[, 23])
attributes(data)$variable.labels[23] <- "[Comment] 	1.4 What is/was done? Please indicate in the comment box the starting date of the diversification initiative, and where possible, the end (Year). Please add this information before selecting the next relevant answer."
names(data)[23] <- "A4_SQ010comment"
# LimeSurvey Field type: F
data[, 24] <- as.numeric(data[, 24])
attributes(data)$variable.labels[24] <- "[New crop mixture(s) were applied (intercropping)] 	1.4 What is/was done? Please indicate in the comment box the starting date of the diversification initiative, and where possible, the end (Year). Please add this information before selecting the next relevant answer."
data[, 24] <- factor(data[, 24], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[24] <- "A4_SQ008"
# LimeSurvey Field type: A
data[, 25] <- as.character(data[, 25])
attributes(data)$variable.labels[25] <- "[Comment] 	1.4 What is/was done? Please indicate in the comment box the starting date of the diversification initiative, and where possible, the end (Year). Please add this information before selecting the next relevant answer."
names(data)[25] <- "A4_SQ008comment"
# LimeSurvey Field type: A
data[, 26] <- as.character(data[, 26])
attributes(data)$variable.labels[26] <- "1.4.1 Please list which crops were included in the rotation"
names(data)[26] <- "A41"
# LimeSurvey Field type: A
data[, 27] <- as.character(data[, 27])
attributes(data)$variable.labels[27] <- "1.4.2 Please indicate which cash or cover crop(s) were included within one season (multiple cropping)."
names(data)[27] <- "A42"
# LimeSurvey Field type: A
data[, 28] <- as.character(data[, 28])
attributes(data)$variable.labels[28] <- "[Which crops were sown?] 1.4.3 Please describe shortly the applied new intercropping"
names(data)[28] <- "A43_SQ001"
# LimeSurvey Field type: A
data[, 29] <- as.character(data[, 29])
attributes(data)$variable.labels[29] <- "[Which was the main crop?] 1.4.3 Please describe shortly the applied new intercropping"
names(data)[29] <- "A43_SQ002"
# LimeSurvey Field type: A
data[, 30] <- as.character(data[, 30])
attributes(data)$variable.labels[30] <- "[Which crops were harvested (grain/silage)?] 1.4.3 Please describe shortly the applied new intercropping"
names(data)[30] <- "A43_SQ003"
# LimeSurvey Field type: A
data[, 31] <- as.character(data[, 31])
attributes(data)$variable.labels[31] <- "[What was the sowing design (mixed, alternated row, strip)?] 1.4.3 Please describe shortly the applied new intercropping"
names(data)[31] <- "A43_SQ004"
# LimeSurvey Field type: A
data[, 32] <- as.character(data[, 32])
attributes(data)$variable.labels[32] <- "[Was the usual sowing density of the crops modified?] 1.4.3 Please describe shortly the applied new intercropping"
names(data)[32] <- "A43_SQ005"
# LimeSurvey Field type: F
data[, 33] <- as.numeric(data[, 33])
attributes(data)$variable.labels[33] <- "[New food product ] 1.5 What is/was the targeted outcome of the initiative?"
data[, 33] <- factor(data[, 33], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[33] <- "A5_SQ001"
# LimeSurvey Field type: F
data[, 34] <- as.numeric(data[, 34])
attributes(data)$variable.labels[34] <- "[New feed product ] 1.5 What is/was the targeted outcome of the initiative?"
data[, 34] <- factor(data[, 34], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[34] <- "A5_SQ002"
# LimeSurvey Field type: F
data[, 35] <- as.numeric(data[, 35])
attributes(data)$variable.labels[35] <- "[New industrial product (e.g. bioenergy carriers, biomaterials, biochemicals)] 1.5 What is/was the targeted outcome of the initiative?"
data[, 35] <- factor(data[, 35], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[35] <- "A5_SQ003"
# LimeSurvey Field type: F
data[, 36] <- as.numeric(data[, 36])
attributes(data)$variable.labels[36] <- "[Improved crop production stability ] 1.5 What is/was the targeted outcome of the initiative?"
data[, 36] <- factor(data[, 36], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[36] <- "A5_SQ004"
# LimeSurvey Field type: F
data[, 37] <- as.numeric(data[, 37])
attributes(data)$variable.labels[37] <- "[Better cash crop quality] 1.5 What is/was the targeted outcome of the initiative?"
data[, 37] <- factor(data[, 37], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[37] <- "A5_SQ005"
# LimeSurvey Field type: F
data[, 38] <- as.numeric(data[, 38])
attributes(data)$variable.labels[38] <- "[Higher yield levels] 1.5 What is/was the targeted outcome of the initiative?"
data[, 38] <- factor(data[, 38], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[38] <- "A5_SQ015"
# LimeSurvey Field type: F
data[, 39] <- as.numeric(data[, 39])
attributes(data)$variable.labels[39] <- "[Lower input levels] 1.5 What is/was the targeted outcome of the initiative?"
data[, 39] <- factor(data[, 39], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[39] <- "A5_SQ014"
# LimeSurvey Field type: F
data[, 40] <- as.numeric(data[, 40])
attributes(data)$variable.labels[40] <- "[Higher economic income] 1.5 What is/was the targeted outcome of the initiative?"
data[, 40] <- factor(data[, 40], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[40] <- "A5_SQ013"
# LimeSurvey Field type: F
data[, 41] <- as.numeric(data[, 41])
attributes(data)$variable.labels[41] <- "[Improved environmental preservation] 1.5 What is/was the targeted outcome of the initiative?"
data[, 41] <- factor(data[, 41], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[41] <- "A5_SQ006"
# LimeSurvey Field type: F
data[, 42] <- as.numeric(data[, 42])
attributes(data)$variable.labels[42] <- "[Landscape aesthetics] 1.5 What is/was the targeted outcome of the initiative?"
data[, 42] <- factor(data[, 42], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[42] <- "A5_SQ007"
# LimeSurvey Field type: F
data[, 43] <- as.numeric(data[, 43])
attributes(data)$variable.labels[43] <- "[New value chain] 1.5 What is/was the targeted outcome of the initiative?"
data[, 43] <- factor(data[, 43], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[43] <- "A5_SQ008"
# LimeSurvey Field type: F
data[, 44] <- as.numeric(data[, 44])
attributes(data)$variable.labels[44] <- "[New certification label] 1.5 What is/was the targeted outcome of the initiative?"
data[, 44] <- factor(data[, 44], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[44] <- "A5_SQ009"
# LimeSurvey Field type: F
data[, 45] <- as.numeric(data[, 45])
attributes(data)$variable.labels[45] <- "[New production cooperative] 1.5 What is/was the targeted outcome of the initiative?"
data[, 45] <- factor(data[, 45], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[45] <- "A5_SQ010"
# LimeSurvey Field type: F
data[, 46] <- as.numeric(data[, 46])
attributes(data)$variable.labels[46] <- "[Compliance to legal requirements] 1.5 What is/was the targeted outcome of the initiative?"
data[, 46] <- factor(data[, 46], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[46] <- "A5_SQ011"
# LimeSurvey Field type: F
data[, 47] <- as.numeric(data[, 47])
attributes(data)$variable.labels[47] <- "[New information support tools for professionals] 1.5 What is/was the targeted outcome of the initiative?"
data[, 47] <- factor(data[, 47], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[47] <- "A5_SQ012"
# LimeSurvey Field type: F
data[, 48] <- as.numeric(data[, 48])
attributes(data)$variable.labels[48] <- "[Other] 1.5 What is/was the targeted outcome of the initiative?"
data[, 48] <- factor(data[, 48], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[48] <- "A5_SQ016"
# LimeSurvey Field type: A
data[, 49] <- as.character(data[, 49])
attributes(data)$variable.labels[49] <- "[Targeted outcome] 1.5a You have selected \\'Other\\' as a targeted outcome of your diversification initiative. Please specify what this targeted outcome is/was."
names(data)[49] <- "A5a_SQ0001"
# LimeSurvey Field type: F
data[, 50] <- as.numeric(data[, 50])
attributes(data)$variable.labels[50] <- "[Farmers] 1.6 Who is/was involved in the initiative?"
data[, 50] <- factor(data[, 50], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[50] <- "A6_SQ001"
# LimeSurvey Field type: F
data[, 51] <- as.numeric(data[, 51])
attributes(data)$variable.labels[51] <- "[Advisors] 1.6 Who is/was involved in the initiative?"
data[, 51] <- factor(data[, 51], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[51] <- "A6_SQ002"
# LimeSurvey Field type: F
data[, 52] <- as.numeric(data[, 52])
attributes(data)$variable.labels[52] <- "[Researchers] 1.6 Who is/was involved in the initiative?"
data[, 52] <- factor(data[, 52], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[52] <- "A6_SQ003"
# LimeSurvey Field type: F
data[, 53] <- as.numeric(data[, 53])
attributes(data)$variable.labels[53] <- "[Commercial companies] 1.6 Who is/was involved in the initiative?"
data[, 53] <- factor(data[, 53], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[53] <- "A6_SQ004"
# LimeSurvey Field type: F
data[, 54] <- as.numeric(data[, 54])
attributes(data)$variable.labels[54] <- "[Consumers] 1.6 Who is/was involved in the initiative?"
data[, 54] <- factor(data[, 54], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[54] <- "A6_SQ005"
# LimeSurvey Field type: F
data[, 55] <- as.numeric(data[, 55])
attributes(data)$variable.labels[55] <- "[Authorities] 1.6 Who is/was involved in the initiative?"
data[, 55] <- factor(data[, 55], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[55] <- "A6_SQ006"
# LimeSurvey Field type: F
data[, 56] <- as.numeric(data[, 56])
attributes(data)$variable.labels[56] <- "[Certification organizations (e.g. for organic farming)] 1.6 Who is/was involved in the initiative?"
data[, 56] <- factor(data[, 56], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[56] <- "A6_SQ007"
# LimeSurvey Field type: F
data[, 57] <- as.numeric(data[, 57])
attributes(data)$variable.labels[57] <- "[NGOs (non governemental organisations that represent civil society)] 1.6 Who is/was involved in the initiative?"
data[, 57] <- factor(data[, 57], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[57] <- "A6_SQ008"
# LimeSurvey Field type: F
data[, 58] <- as.numeric(data[, 58])
attributes(data)$variable.labels[58] <- "[Other] 1.6 Who is/was involved in the initiative?"
data[, 58] <- factor(data[, 58], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[58] <- "A6_SQ099"
# LimeSurvey Field type: A
data[, 59] <- as.character(data[, 59])
attributes(data)$variable.labels[59] <- "[Involved actors] 1.6a You have selected \\'Other\\' for involved actors of your diversification initiative. Please specify who these actors are/were."
names(data)[59] <- "A6a_SQ0001"
# LimeSurvey Field type: F
data[, 60] <- as.numeric(data[, 60])
attributes(data)$variable.labels[60] <- "[As individual farmer(s)] 1.6.1 How were farmers involved in the diversification initiative?"
data[, 60] <- factor(data[, 60], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[60] <- "A61_SQ001"
# LimeSurvey Field type: F
data[, 61] <- as.numeric(data[, 61])
attributes(data)$variable.labels[61] <- "[Through informal farmer group(s)] 1.6.1 How were farmers involved in the diversification initiative?"
data[, 61] <- factor(data[, 61], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[61] <- "A61_SQ002"
# LimeSurvey Field type: F
data[, 62] <- as.numeric(data[, 62])
attributes(data)$variable.labels[62] <- "[Through farmer association(s)] 1.6.1 How were farmers involved in the diversification initiative?"
data[, 62] <- factor(data[, 62], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[62] <- "A61_SQ003"
# LimeSurvey Field type: A
data[, 63] <- as.character(data[, 63])
attributes(data)$variable.labels[63] <- "[Other] 1.6.1 How were farmers involved in the diversification initiative?"
names(data)[63] <- "A61_other"
# LimeSurvey Field type: F
data[, 64] <- as.numeric(data[, 64])
attributes(data)$variable.labels[64] <- "[As individual advisor(s)] 1.6.2 How were advisors involved in the diversification initiative?"
data[, 64] <- factor(data[, 64], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[64] <- "A62_SQ001"
# LimeSurvey Field type: F
data[, 65] <- as.numeric(data[, 65])
attributes(data)$variable.labels[65] <- "[Through (national or regional) state advisory service organization(s)] 1.6.2 How were advisors involved in the diversification initiative?"
data[, 65] <- factor(data[, 65], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[65] <- "A62_SQ002"
# LimeSurvey Field type: F
data[, 66] <- as.numeric(data[, 66])
attributes(data)$variable.labels[66] <- "[Through private advisory service organization(s)] 1.6.2 How were advisors involved in the diversification initiative?"
data[, 66] <- factor(data[, 66], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[66] <- "A62_SQ003"
# LimeSurvey Field type: A
data[, 67] <- as.character(data[, 67])
attributes(data)$variable.labels[67] <- "[Other] 1.6.2 How were advisors involved in the diversification initiative?"
names(data)[67] <- "A62_other"
# LimeSurvey Field type: F
data[, 68] <- as.numeric(data[, 68])
attributes(data)$variable.labels[68] <- "[As individual researcher(s)] 1.6.3 How were researchers involved in the diversification initiative?"
data[, 68] <- factor(data[, 68], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[68] <- "A63_SQ001"
# LimeSurvey Field type: F
data[, 69] <- as.numeric(data[, 69])
attributes(data)$variable.labels[69] <- "[Through University] 1.6.3 How were researchers involved in the diversification initiative?"
data[, 69] <- factor(data[, 69], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[69] <- "A63_SQ002"
# LimeSurvey Field type: F
data[, 70] <- as.numeric(data[, 70])
attributes(data)$variable.labels[70] <- "[Through state research organization(s)] 1.6.3 How were researchers involved in the diversification initiative?"
data[, 70] <- factor(data[, 70], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[70] <- "A63_SQ003"
# LimeSurvey Field type: F
data[, 71] <- as.numeric(data[, 71])
attributes(data)$variable.labels[71] <- "[Through private research organization(s)] 1.6.3 How were researchers involved in the diversification initiative?"
data[, 71] <- factor(data[, 71], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[71] <- "A63_SQ004"
# LimeSurvey Field type: A
data[, 72] <- as.character(data[, 72])
attributes(data)$variable.labels[72] <- "[Other] 1.6.3 How were researchers involved in the diversification initiative?"
names(data)[72] <- "A63_other"
# LimeSurvey Field type: F
data[, 73] <- as.numeric(data[, 73])
attributes(data)$variable.labels[73] <- "[Seed producer(s)] 1.6.4 What kind of commercial companies were involved in the diversification initiative?"
data[, 73] <- factor(data[, 73], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[73] <- "A64_SQ001"
# LimeSurvey Field type: F
data[, 74] <- as.numeric(data[, 74])
attributes(data)$variable.labels[74] <- "[Plant protection companie(s)] 1.6.4 What kind of commercial companies were involved in the diversification initiative?"
data[, 74] <- factor(data[, 74], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[74] <- "A64_SQ006"
# LimeSurvey Field type: F
data[, 75] <- as.numeric(data[, 75])
attributes(data)$variable.labels[75] <- "[Fertilizer producer(s)] 1.6.4 What kind of commercial companies were involved in the diversification initiative?"
data[, 75] <- factor(data[, 75], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[75] <- "A64_SQ008"
# LimeSurvey Field type: F
data[, 76] <- as.numeric(data[, 76])
attributes(data)$variable.labels[76] <- "[Machine producer(s)] 1.6.4 What kind of commercial companies were involved in the diversification initiative?"
data[, 76] <- factor(data[, 76], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[76] <- "A64_SQ007"
# LimeSurvey Field type: F
data[, 77] <- as.numeric(data[, 77])
attributes(data)$variable.labels[77] <- "[Processor(s)] 1.6.4 What kind of commercial companies were involved in the diversification initiative?"
data[, 77] <- factor(data[, 77], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[77] <- "A64_SQ005"
# LimeSurvey Field type: F
data[, 78] <- as.numeric(data[, 78])
attributes(data)$variable.labels[78] <- "[Trader(s)] 1.6.4 What kind of commercial companies were involved in the diversification initiative?"
data[, 78] <- factor(data[, 78], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[78] <- "A64_SQ002"
# LimeSurvey Field type: F
data[, 79] <- as.numeric(data[, 79])
attributes(data)$variable.labels[79] <- "[Retailer(s)] 1.6.4 What kind of commercial companies were involved in the diversification initiative?"
data[, 79] <- factor(data[, 79], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[79] <- "A64_SQ003"
# LimeSurvey Field type: A
data[, 80] <- as.character(data[, 80])
attributes(data)$variable.labels[80] <- "[Other] 1.6.4 What kind of commercial companies were involved in the diversification initiative?"
names(data)[80] <- "A64_other"
# LimeSurvey Field type: F
data[, 81] <- as.numeric(data[, 81])
attributes(data)$variable.labels[81] <- "[Ministry] 1.6.5 What kind of authorities were involved in the diversification initiative?"
data[, 81] <- factor(data[, 81], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[81] <- "A65_SQ001"
# LimeSurvey Field type: F
data[, 82] <- as.numeric(data[, 82])
attributes(data)$variable.labels[82] <- "[National authority] 1.6.5 What kind of authorities were involved in the diversification initiative?"
data[, 82] <- factor(data[, 82], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[82] <- "A65_SQ004"
# LimeSurvey Field type: F
data[, 83] <- as.numeric(data[, 83])
attributes(data)$variable.labels[83] <- "[Regional authority] 1.6.5 What kind of authorities were involved in the diversification initiative?"
data[, 83] <- factor(data[, 83], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[83] <- "A65_SQ002"
# LimeSurvey Field type: F
data[, 84] <- as.numeric(data[, 84])
attributes(data)$variable.labels[84] <- "[Local authority] 1.6.5 What kind of authorities were involved in the diversification initiative?"
data[, 84] <- factor(data[, 84], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[84] <- "A65_SQ003"
# LimeSurvey Field type: A
data[, 85] <- as.character(data[, 85])
attributes(data)$variable.labels[85] <- "[Other] 1.6.5 What kind of authorities were involved in the diversification initiative?"
names(data)[85] <- "A65_other"
# LimeSurvey Field type: A
data[, 86] <- as.character(data[, 86])
attributes(data)$variable.labels[86] <- "1.6.6 What kind of certification organizations were involved in the diversification initiative? Please specify."
names(data)[86] <- "A66"
# LimeSurvey Field type: F
data[, 87] <- as.numeric(data[, 87])
attributes(data)$variable.labels[87] <- "[Consumer organizations] 1.6.7 What kind of NGOs were involved in the diversification initiative?"
data[, 87] <- factor(data[, 87], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[87] <- "A67_SQ001"
# LimeSurvey Field type: F
data[, 88] <- as.numeric(data[, 88])
attributes(data)$variable.labels[88] <- "[Environmental organizations] 1.6.7 What kind of NGOs were involved in the diversification initiative?"
data[, 88] <- factor(data[, 88], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[88] <- "A67_SQ002"
# LimeSurvey Field type: F
data[, 89] <- as.numeric(data[, 89])
attributes(data)$variable.labels[89] <- "[Farmer organizations] 1.6.7 What kind of NGOs were involved in the diversification initiative?"
data[, 89] <- factor(data[, 89], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[89] <- "A67_SQ003"
# LimeSurvey Field type: A
data[, 90] <- as.character(data[, 90])
attributes(data)$variable.labels[90] <- "[Other] 1.6.7 What kind of NGOs were involved in the diversification initiative?"
names(data)[90] <- "A67_other"
# LimeSurvey Field type: F
data[, 91] <- as.numeric(data[, 91])
attributes(data)$variable.labels[91] <- "1.7 Were you involved?"
data[, 91] <- factor(data[, 91], levels=c(1,2),labels=c("Yes", "No"))
names(data)[91] <- "A7"
# LimeSurvey Field type: F
data[, 92] <- as.numeric(data[, 92])
attributes(data)$variable.labels[92] <- "[As a farmer] 1.7.1 Please specify how you were involved."
data[, 92] <- factor(data[, 92], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[92] <- "A71_SQ001"
# LimeSurvey Field type: F
data[, 93] <- as.numeric(data[, 93])
attributes(data)$variable.labels[93] <- "[As an advisor] 1.7.1 Please specify how you were involved."
data[, 93] <- factor(data[, 93], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[93] <- "A71_SQ002"
# LimeSurvey Field type: F
data[, 94] <- as.numeric(data[, 94])
attributes(data)$variable.labels[94] <- "[As a researcher] 1.7.1 Please specify how you were involved."
data[, 94] <- factor(data[, 94], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[94] <- "A71_SQ003"
# LimeSurvey Field type: F
data[, 95] <- as.numeric(data[, 95])
attributes(data)$variable.labels[95] <- "[As a commercial company] 1.7.1 Please specify how you were involved."
data[, 95] <- factor(data[, 95], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[95] <- "A71_SQ004"
# LimeSurvey Field type: F
data[, 96] <- as.numeric(data[, 96])
attributes(data)$variable.labels[96] <- "[As a consumer] 1.7.1 Please specify how you were involved."
data[, 96] <- factor(data[, 96], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[96] <- "A71_SQ005"
# LimeSurvey Field type: F
data[, 97] <- as.numeric(data[, 97])
attributes(data)$variable.labels[97] <- "[As an authority] 1.7.1 Please specify how you were involved."
data[, 97] <- factor(data[, 97], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[97] <- "A71_SQ006"
# LimeSurvey Field type: F
data[, 98] <- as.numeric(data[, 98])
attributes(data)$variable.labels[98] <- "[As a certification organization] 1.7.1 Please specify how you were involved."
data[, 98] <- factor(data[, 98], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[98] <- "A71_SQ007"
# LimeSurvey Field type: F
data[, 99] <- as.numeric(data[, 99])
attributes(data)$variable.labels[99] <- "[As an NGO] 1.7.1 Please specify how you were involved."
data[, 99] <- factor(data[, 99], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[99] <- "A71_SQ008"
# LimeSurvey Field type: F
data[, 100] <- as.numeric(data[, 100])
attributes(data)$variable.labels[100] <- "[{438634X176X2730SQ0001.shown}] 1.7.1 Please specify how you were involved."
data[, 100] <- factor(data[, 100], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[100] <- "A71_SQ0001"
# LimeSurvey Field type: A
data[, 101] <- as.character(data[, 101])
attributes(data)$variable.labels[101] <- "1.8 Did you initiate the diversification initiative?"
data[, 101] <- factor(data[, 101], levels=c("A1","A2"),labels=c("Yes", "No"))
names(data)[101] <- "A8"
# LimeSurvey Field type: A
data[, 102] <- as.character(data[, 102])
attributes(data)$variable.labels[102] <- "[Other] 1.8 Did you initiate the diversification initiative?"
names(data)[102] <- "A8_other"
# LimeSurvey Field type: F
data[, 103] <- as.numeric(data[, 103])
attributes(data)$variable.labels[103] <- "[Farmers] 1.9 Who initiated the diversification initiative"
data[, 103] <- factor(data[, 103], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[103] <- "A9_SQ001"
# LimeSurvey Field type: F
data[, 104] <- as.numeric(data[, 104])
attributes(data)$variable.labels[104] <- "[Advisors] 1.9 Who initiated the diversification initiative"
data[, 104] <- factor(data[, 104], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[104] <- "A9_SQ002"
# LimeSurvey Field type: F
data[, 105] <- as.numeric(data[, 105])
attributes(data)$variable.labels[105] <- "[Researchers] 1.9 Who initiated the diversification initiative"
data[, 105] <- factor(data[, 105], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[105] <- "A9_SQ003"
# LimeSurvey Field type: F
data[, 106] <- as.numeric(data[, 106])
attributes(data)$variable.labels[106] <- "[Commercial companies] 1.9 Who initiated the diversification initiative"
data[, 106] <- factor(data[, 106], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[106] <- "A9_SQ004"
# LimeSurvey Field type: F
data[, 107] <- as.numeric(data[, 107])
attributes(data)$variable.labels[107] <- "[Consumers] 1.9 Who initiated the diversification initiative"
data[, 107] <- factor(data[, 107], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[107] <- "A9_SQ005"
# LimeSurvey Field type: F
data[, 108] <- as.numeric(data[, 108])
attributes(data)$variable.labels[108] <- "[Authorities] 1.9 Who initiated the diversification initiative"
data[, 108] <- factor(data[, 108], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[108] <- "A9_SQ006"
# LimeSurvey Field type: F
data[, 109] <- as.numeric(data[, 109])
attributes(data)$variable.labels[109] <- "[Certification organizations (e.g. for organic farming)] 1.9 Who initiated the diversification initiative"
data[, 109] <- factor(data[, 109], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[109] <- "A9_SQ007"
# LimeSurvey Field type: F
data[, 110] <- as.numeric(data[, 110])
attributes(data)$variable.labels[110] <- "[NGOs (non governemental organisations that represent civil society)] 1.9 Who initiated the diversification initiative"
data[, 110] <- factor(data[, 110], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[110] <- "A9_SQ008"
# LimeSurvey Field type: F
data[, 111] <- as.numeric(data[, 111])
attributes(data)$variable.labels[111] <- "[{438634X176X2730SQ0001.shown}] 1.9 Who initiated the diversification initiative"
data[, 111] <- factor(data[, 111], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[111] <- "A9_SQ0001"
# LimeSurvey Field type: F
data[, 112] <- as.numeric(data[, 112])
attributes(data)$variable.labels[112] <- "[On commercial farms, by farmer(s)/farmers association(s)] 1.10 Where was the diversification first tested?"
data[, 112] <- factor(data[, 112], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[112] <- "A10_SQ001"
# LimeSurvey Field type: F
data[, 113] <- as.numeric(data[, 113])
attributes(data)$variable.labels[113] <- "[On private trial plots, by commercial firm(s)] 1.10 Where was the diversification first tested?"
data[, 113] <- factor(data[, 113], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[113] <- "A10_SQ002"
# LimeSurvey Field type: F
data[, 114] <- as.numeric(data[, 114])
attributes(data)$variable.labels[114] <- "[At scientific institutions (research centre, university)] 1.10 Where was the diversification first tested?"
data[, 114] <- factor(data[, 114], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[114] <- "A10_SQ003"
# LimeSurvey Field type: A
data[, 115] <- as.character(data[, 115])
attributes(data)$variable.labels[115] <- "[Other] 1.10 Where was the diversification first tested?"
names(data)[115] <- "A10_other"
# LimeSurvey Field type: F
data[, 116] <- as.numeric(data[, 116])
attributes(data)$variable.labels[116] <- "[Arable cropping systems] 1.11a What type(s) of farming systems were involved in the diversification initiative?"
data[, 116] <- factor(data[, 116], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[116] <- "A11a_SQ001"
# LimeSurvey Field type: F
data[, 117] <- as.numeric(data[, 117])
attributes(data)$variable.labels[117] <- "[Vegetable cropping systems] 1.11a What type(s) of farming systems were involved in the diversification initiative?"
data[, 117] <- factor(data[, 117], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[117] <- "A11a_SQ002"
# LimeSurvey Field type: F
data[, 118] <- as.numeric(data[, 118])
attributes(data)$variable.labels[118] <- "[Arable cropping systems with animals] 1.11a What type(s) of farming systems were involved in the diversification initiative?"
data[, 118] <- factor(data[, 118], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[118] <- "A11a_SQ003"
# LimeSurvey Field type: F
data[, 119] <- as.numeric(data[, 119])
attributes(data)$variable.labels[119] <- "[Vegetable cropping systems with animals] 1.11a What type(s) of farming systems were involved in the diversification initiative?"
data[, 119] <- factor(data[, 119], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[119] <- "A11a_SQ004"
# LimeSurvey Field type: A
data[, 120] <- as.character(data[, 120])
attributes(data)$variable.labels[120] <- "[Other] 1.11a What type(s) of farming systems were involved in the diversification initiative?"
names(data)[120] <- "A11a_other"
# LimeSurvey Field type: A
data[, 121] <- as.character(data[, 121])
attributes(data)$variable.labels[121] <- "1.11 Was the crop diversification initiative conducted on certified organic area?"
data[, 121] <- factor(data[, 121], levels=c("A1","A2","A3","A4"),labels=c("No, the area was not certified organic", "No, the area was not certified organic, but no pesticide and chemical fertilizers were applied", "Yes, the area was certified organic completely", "Both organic and non-organic areas were involved"))
names(data)[121] <- "A11"
# LimeSurvey Field type: F
data[, 122] <- as.numeric(data[, 122])
attributes(data)$variable.labels[122] <- "[EU project(s). Please specify project name(s).] 1.12 How is/was the initiative funded?"
data[, 122] <- factor(data[, 122], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[122] <- "A12_SQ001"
# LimeSurvey Field type: A
data[, 123] <- as.character(data[, 123])
attributes(data)$variable.labels[123] <- "[Comment] 1.12 How is/was the initiative funded?"
names(data)[123] <- "A12_SQ001comment"
# LimeSurvey Field type: F
data[, 124] <- as.numeric(data[, 124])
attributes(data)$variable.labels[124] <- "[National or regional project funding. Please specify project name(s).] 1.12 How is/was the initiative funded?"
data[, 124] <- factor(data[, 124], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[124] <- "A12_SQ002"
# LimeSurvey Field type: A
data[, 125] <- as.character(data[, 125])
attributes(data)$variable.labels[125] <- "[Comment] 1.12 How is/was the initiative funded?"
names(data)[125] <- "A12_SQ002comment"
# LimeSurvey Field type: F
data[, 126] <- as.numeric(data[, 126])
attributes(data)$variable.labels[126] <- "[Commercial company funding. Please specify company name(s).] 1.12 How is/was the initiative funded?"
data[, 126] <- factor(data[, 126], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[126] <- "A12_SQ003"
# LimeSurvey Field type: A
data[, 127] <- as.character(data[, 127])
attributes(data)$variable.labels[127] <- "[Comment] 1.12 How is/was the initiative funded?"
names(data)[127] <- "A12_SQ003comment"
# LimeSurvey Field type: F
data[, 128] <- as.numeric(data[, 128])
attributes(data)$variable.labels[128] <- "[Private non-commercial funding. Please specify e.g. Foundation name(s).] 1.12 How is/was the initiative funded?"
data[, 128] <- factor(data[, 128], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[128] <- "A12_SQ004"
# LimeSurvey Field type: A
data[, 129] <- as.character(data[, 129])
attributes(data)$variable.labels[129] <- "[Comment] 1.12 How is/was the initiative funded?"
names(data)[129] <- "A12_SQ004comment"
# LimeSurvey Field type: F
data[, 130] <- as.numeric(data[, 130])
attributes(data)$variable.labels[130] <- "[Self-funding. Please specify who contributed.] 1.12 How is/was the initiative funded?"
data[, 130] <- factor(data[, 130], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[130] <- "A12_SQ005"
# LimeSurvey Field type: A
data[, 131] <- as.character(data[, 131])
attributes(data)$variable.labels[131] <- "[Comment] 1.12 How is/was the initiative funded?"
names(data)[131] <- "A12_SQ005comment"
# LimeSurvey Field type: A
data[, 132] <- as.character(data[, 132])
attributes(data)$variable.labels[132] <- "[Other] 1.12 How is/was the initiative funded?"
names(data)[132] <- "A12_other"
# LimeSurvey Field type: A
data[, 133] <- as.character(data[, 133])
attributes(data)$variable.labels[133] <- "[Other comment] 1.12 How is/was the initiative funded?"
names(data)[133] <- "A12_othercomment"
# LimeSurvey Field type: A
data[, 134] <- as.character(data[, 134])
attributes(data)$variable.labels[134] <- "1.13 What was the magnitude of funding dedicated to the initiative? (in total in Euro)"
names(data)[134] <- "A13"
# LimeSurvey Field type: F
data[, 135] <- as.numeric(data[, 135])
attributes(data)$variable.labels[135] <- "[Breeding] 1.14 Did the diversification initiative include the following upstream value chain levels (levels that supply agricultural production)?"
data[, 135] <- factor(data[, 135], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[135] <- "A14_SQ001"
# LimeSurvey Field type: F
data[, 136] <- as.numeric(data[, 136])
attributes(data)$variable.labels[136] <- "[Seed production] 1.14 Did the diversification initiative include the following upstream value chain levels (levels that supply agricultural production)?"
data[, 136] <- factor(data[, 136], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[136] <- "A14_SQ002"
# LimeSurvey Field type: F
data[, 137] <- as.numeric(data[, 137])
attributes(data)$variable.labels[137] <- "[Development of inputs needed for production ] 1.14 Did the diversification initiative include the following upstream value chain levels (levels that supply agricultural production)?"
data[, 137] <- factor(data[, 137], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[137] <- "A14_SQ003"
# LimeSurvey Field type: F
data[, 138] <- as.numeric(data[, 138])
attributes(data)$variable.labels[138] <- "[Organization of inputs needed for production] 1.14 Did the diversification initiative include the following upstream value chain levels (levels that supply agricultural production)?"
data[, 138] <- factor(data[, 138], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[138] <- "A14_SQ004"
# LimeSurvey Field type: F
data[, 139] <- as.numeric(data[, 139])
attributes(data)$variable.labels[139] <- "[Machinery development] 1.14 Did the diversification initiative include the following upstream value chain levels (levels that supply agricultural production)?"
data[, 139] <- factor(data[, 139], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[139] <- "A14_SQ005"
# LimeSurvey Field type: F
data[, 140] <- as.numeric(data[, 140])
attributes(data)$variable.labels[140] <- "[Organization of machinery needed for production] 1.14 Did the diversification initiative include the following upstream value chain levels (levels that supply agricultural production)?"
data[, 140] <- factor(data[, 140], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[140] <- "A14_SQ006"
# LimeSurvey Field type: F
data[, 141] <- as.numeric(data[, 141])
attributes(data)$variable.labels[141] <- "[Other] 1.14 Did the diversification initiative include the following upstream value chain levels (levels that supply agricultural production)?"
data[, 141] <- factor(data[, 141], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[141] <- "A14_SQ008"
# LimeSurvey Field type: F
data[, 142] <- as.numeric(data[, 142])
attributes(data)$variable.labels[142] <- "[None] 1.14 Did the diversification initiative include the following upstream value chain levels (levels that supply agricultural production)?"
data[, 142] <- factor(data[, 142], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[142] <- "A14_SQ077"
# LimeSurvey Field type: A
data[, 143] <- as.character(data[, 143])
attributes(data)$variable.labels[143] <- "[Upstream value chain level] 1.4a You have selected \\'Other\\' as an upstream value chain level that was included in your diversification initiative. Please specify what this upstream value chain level is/was."
names(data)[143] <- "A14a_SQ0001"
# LimeSurvey Field type: F
data[, 144] <- as.numeric(data[, 144])
attributes(data)$variable.labels[144] <- "[Quality assurance] 1.15 Did the initiative include the following down stream value chain levels (levels that come after primary production)?"
data[, 144] <- factor(data[, 144], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[144] <- "A15_SQ01"
# LimeSurvey Field type: F
data[, 145] <- as.numeric(data[, 145])
attributes(data)$variable.labels[145] <- "[Transportation] 1.15 Did the initiative include the following down stream value chain levels (levels that come after primary production)?"
data[, 145] <- factor(data[, 145], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[145] <- "A15_SQ02"
# LimeSurvey Field type: F
data[, 146] <- as.numeric(data[, 146])
attributes(data)$variable.labels[146] <- "[Logistics] 1.15 Did the initiative include the following down stream value chain levels (levels that come after primary production)?"
data[, 146] <- factor(data[, 146], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[146] <- "A15_SQ03"
# LimeSurvey Field type: F
data[, 147] <- as.numeric(data[, 147])
attributes(data)$variable.labels[147] <- "[Processing] 1.15 Did the initiative include the following down stream value chain levels (levels that come after primary production)?"
data[, 147] <- factor(data[, 147], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[147] <- "A15_SQ04"
# LimeSurvey Field type: F
data[, 148] <- as.numeric(data[, 148])
attributes(data)$variable.labels[148] <- "[Marketing] 1.15 Did the initiative include the following down stream value chain levels (levels that come after primary production)?"
data[, 148] <- factor(data[, 148], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[148] <- "A15_SQ05"
# LimeSurvey Field type: F
data[, 149] <- as.numeric(data[, 149])
attributes(data)$variable.labels[149] <- "[Sales] 1.15 Did the initiative include the following down stream value chain levels (levels that come after primary production)?"
data[, 149] <- factor(data[, 149], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[149] <- "A15_SQ06"
# LimeSurvey Field type: F
data[, 150] <- as.numeric(data[, 150])
attributes(data)$variable.labels[150] <- "[Other] 1.15 Did the initiative include the following down stream value chain levels (levels that come after primary production)?"
data[, 150] <- factor(data[, 150], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[150] <- "A15_SQ08"
# LimeSurvey Field type: F
data[, 151] <- as.numeric(data[, 151])
attributes(data)$variable.labels[151] <- "[None] 1.15 Did the initiative include the following down stream value chain levels (levels that come after primary production)?"
data[, 151] <- factor(data[, 151], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[151] <- "A15_SQ77"
# LimeSurvey Field type: A
data[, 152] <- as.character(data[, 152])
attributes(data)$variable.labels[152] <- "[Downstream value chain level] 1.5a You have selected \\'Other\\' as a downstream value chain level that was included in your diversification initiative. Please specify what this downstream value chain level is/was."
names(data)[152] <- "A15a_SQ0002"
# LimeSurvey Field type: F
data[, 153] <- as.numeric(data[, 153])
attributes(data)$variable.labels[153] <- "[Internal communication among participants (e.g. sharing experiences) (please specify communication tools)] 1.16 What kind of communication approaches were used during the diversification initiative?"
data[, 153] <- factor(data[, 153], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[153] <- "A16_SQ001"
# LimeSurvey Field type: A
data[, 154] <- as.character(data[, 154])
attributes(data)$variable.labels[154] <- "[Comment] 1.16 What kind of communication approaches were used during the diversification initiative?"
names(data)[154] <- "A16_SQ001comment"
# LimeSurvey Field type: F
data[, 155] <- as.numeric(data[, 155])
attributes(data)$variable.labels[155] <- "[External communication for professionals (e.g. publishing results) (please specify communication tools)] 1.16 What kind of communication approaches were used during the diversification initiative?"
data[, 155] <- factor(data[, 155], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[155] <- "A16_SQ002"
# LimeSurvey Field type: A
data[, 156] <- as.character(data[, 156])
attributes(data)$variable.labels[156] <- "[Comment] 1.16 What kind of communication approaches were used during the diversification initiative?"
names(data)[156] <- "A16_SQ002comment"
# LimeSurvey Field type: F
data[, 157] <- as.numeric(data[, 157])
attributes(data)$variable.labels[157] <- "[External communication for the public (e.g. media appearances) (please specify communication tools)] 1.16 What kind of communication approaches were used during the diversification initiative?"
data[, 157] <- factor(data[, 157], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[157] <- "A16_SQ003"
# LimeSurvey Field type: A
data[, 158] <- as.character(data[, 158])
attributes(data)$variable.labels[158] <- "[Comment] 1.16 What kind of communication approaches were used during the diversification initiative?"
names(data)[158] <- "A16_SQ003comment"
# LimeSurvey Field type: A
data[, 159] <- as.character(data[, 159])
attributes(data)$variable.labels[159] <- "[Other] 1.16 What kind of communication approaches were used during the diversification initiative?"
names(data)[159] <- "A16_other"
# LimeSurvey Field type: A
data[, 160] <- as.character(data[, 160])
attributes(data)$variable.labels[160] <- "[Other comment] 1.16 What kind of communication approaches were used during the diversification initiative?"
names(data)[160] <- "A16_othercomment"
# LimeSurvey Field type: A
data[, 161] <- as.character(data[, 161])
attributes(data)$variable.labels[161] <- "1.17 Did the initiative use a specific co-innovation methodology? If so, please specify."
names(data)[161] <- "A17"
# LimeSurvey Field type: A
data[, 162] <- as.character(data[, 162])
attributes(data)$variable.labels[162] <- "[Website(s)] 1.18 Where can we learn more about the initiative?"
names(data)[162] <- "A18_SQ001"
# LimeSurvey Field type: A
data[, 163] <- as.character(data[, 163])
attributes(data)$variable.labels[163] <- "[Publication(s)] 1.18 Where can we learn more about the initiative?"
names(data)[163] <- "A18_SQ002"
# LimeSurvey Field type: A
data[, 164] <- as.character(data[, 164])
attributes(data)$variable.labels[164] <- "[Video(s)] 1.18 Where can we learn more about the initiative?"
names(data)[164] <- "A18_SQ003"
# LimeSurvey Field type: A
data[, 165] <- as.character(data[, 165])
attributes(data)$variable.labels[165] <- "[Other] 1.18 Where can we learn more about the initiative?"
names(data)[165] <- "A18_SQ005"
# LimeSurvey Field type: F
data[, 166] <- as.numeric(data[, 166])
attributes(data)$variable.labels[166] <- "[ The diversification initiative was] 2.1 In your opinion is/was the diversification initiative successful?"
data[, 166] <- factor(data[, 166], levels=c(1,2,3,4,5),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[166] <- "B1_SQ001"
# LimeSurvey Field type: A
data[, 167] <- as.character(data[, 167])
attributes(data)$variable.labels[167] <- "[New food product ] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 167] <- factor(data[, 167], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[167] <- "B2_SQ001"
# LimeSurvey Field type: A
data[, 168] <- as.character(data[, 168])
attributes(data)$variable.labels[168] <- "[New feed product ] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 168] <- factor(data[, 168], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[168] <- "B2_SQ002"
# LimeSurvey Field type: A
data[, 169] <- as.character(data[, 169])
attributes(data)$variable.labels[169] <- "[New industrial product (e.g. bioenergy carriers, biomaterials, biochemicals) ] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 169] <- factor(data[, 169], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[169] <- "B2_SQ003"
# LimeSurvey Field type: A
data[, 170] <- as.character(data[, 170])
attributes(data)$variable.labels[170] <- "[Improved crop production stability] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 170] <- factor(data[, 170], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[170] <- "B2_SQ004"
# LimeSurvey Field type: A
data[, 171] <- as.character(data[, 171])
attributes(data)$variable.labels[171] <- "[Better cash crop quality] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 171] <- factor(data[, 171], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[171] <- "B2_SQ005"
# LimeSurvey Field type: A
data[, 172] <- as.character(data[, 172])
attributes(data)$variable.labels[172] <- "[Higher yield levels] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 172] <- factor(data[, 172], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[172] <- "B2_SQ015"
# LimeSurvey Field type: A
data[, 173] <- as.character(data[, 173])
attributes(data)$variable.labels[173] <- "[Lower input levels] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 173] <- factor(data[, 173], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[173] <- "B2_SQ014"
# LimeSurvey Field type: A
data[, 174] <- as.character(data[, 174])
attributes(data)$variable.labels[174] <- "[Higher economic income] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 174] <- factor(data[, 174], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[174] <- "B2_SQ013"
# LimeSurvey Field type: A
data[, 175] <- as.character(data[, 175])
attributes(data)$variable.labels[175] <- "[Improved environmental preservation] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 175] <- factor(data[, 175], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[175] <- "B2_SQ006"
# LimeSurvey Field type: A
data[, 176] <- as.character(data[, 176])
attributes(data)$variable.labels[176] <- "[Landscape aesthetics] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 176] <- factor(data[, 176], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[176] <- "B2_SQ007"
# LimeSurvey Field type: A
data[, 177] <- as.character(data[, 177])
attributes(data)$variable.labels[177] <- "[New value chain ] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 177] <- factor(data[, 177], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[177] <- "B2_SQ008"
# LimeSurvey Field type: A
data[, 178] <- as.character(data[, 178])
attributes(data)$variable.labels[178] <- "[New certification label ] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 178] <- factor(data[, 178], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[178] <- "B2_SQ009"
# LimeSurvey Field type: A
data[, 179] <- as.character(data[, 179])
attributes(data)$variable.labels[179] <- "[New production cooperative] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 179] <- factor(data[, 179], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[179] <- "B2_SQ010"
# LimeSurvey Field type: A
data[, 180] <- as.character(data[, 180])
attributes(data)$variable.labels[180] <- "[Compliance to legal requirements] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 180] <- factor(data[, 180], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[180] <- "B2_SQ011"
# LimeSurvey Field type: A
data[, 181] <- as.character(data[, 181])
attributes(data)$variable.labels[181] <- "[New information support tools for professionals] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 181] <- factor(data[, 181], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[181] <- "B2_SQ012"
# LimeSurvey Field type: A
data[, 182] <- as.character(data[, 182])
attributes(data)$variable.labels[182] <- "[{438634X176X2718SQ0001.shown}] 2.2 How successful is/was the realization of the diversification initiative’s targeted outcomes?"
data[, 182] <- factor(data[, 182], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[182] <- "B2_SQ0001"
# LimeSurvey Field type: A
data[, 183] <- as.character(data[, 183])
attributes(data)$variable.labels[183] <- "[Breeding] 2.3 Was the (re-)organization of upstream value chain levels (levels supplying primary production) successful?"
data[, 183] <- factor(data[, 183], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[183] <- "B3_SQ001"
# LimeSurvey Field type: A
data[, 184] <- as.character(data[, 184])
attributes(data)$variable.labels[184] <- "[Seed production] 2.3 Was the (re-)organization of upstream value chain levels (levels supplying primary production) successful?"
data[, 184] <- factor(data[, 184], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[184] <- "B3_SQ002"
# LimeSurvey Field type: A
data[, 185] <- as.character(data[, 185])
attributes(data)$variable.labels[185] <- "[Development of inputs needed for production ] 2.3 Was the (re-)organization of upstream value chain levels (levels supplying primary production) successful?"
data[, 185] <- factor(data[, 185], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[185] <- "B3_SQ003"
# LimeSurvey Field type: A
data[, 186] <- as.character(data[, 186])
attributes(data)$variable.labels[186] <- "[Organization of inputs needed for production] 2.3 Was the (re-)organization of upstream value chain levels (levels supplying primary production) successful?"
data[, 186] <- factor(data[, 186], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[186] <- "B3_SQ004"
# LimeSurvey Field type: A
data[, 187] <- as.character(data[, 187])
attributes(data)$variable.labels[187] <- "[Machinery development] 2.3 Was the (re-)organization of upstream value chain levels (levels supplying primary production) successful?"
data[, 187] <- factor(data[, 187], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[187] <- "B3_SQ005"
# LimeSurvey Field type: A
data[, 188] <- as.character(data[, 188])
attributes(data)$variable.labels[188] <- "[Organization of machinery needed for production] 2.3 Was the (re-)organization of upstream value chain levels (levels supplying primary production) successful?"
data[, 188] <- factor(data[, 188], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[188] <- "B3_SQ006"
# LimeSurvey Field type: A
data[, 189] <- as.character(data[, 189])
attributes(data)$variable.labels[189] <- "[{438634X176X2820SQ0001.shown}] 2.3 Was the (re-)organization of upstream value chain levels (levels supplying primary production) successful?"
data[, 189] <- factor(data[, 189], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[189] <- "B3_SQ008"
# LimeSurvey Field type: A
data[, 190] <- as.character(data[, 190])
attributes(data)$variable.labels[190] <- "[Quality assurance] 2.4 Was the (re-)organization of downstream value chain levels (levels after primary production) successful?"
data[, 190] <- factor(data[, 190], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[190] <- "B4_SQ01"
# LimeSurvey Field type: A
data[, 191] <- as.character(data[, 191])
attributes(data)$variable.labels[191] <- "[Transportation] 2.4 Was the (re-)organization of downstream value chain levels (levels after primary production) successful?"
data[, 191] <- factor(data[, 191], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[191] <- "B4_SQ02"
# LimeSurvey Field type: A
data[, 192] <- as.character(data[, 192])
attributes(data)$variable.labels[192] <- "[Logistics] 2.4 Was the (re-)organization of downstream value chain levels (levels after primary production) successful?"
data[, 192] <- factor(data[, 192], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[192] <- "B4_SQ03"
# LimeSurvey Field type: A
data[, 193] <- as.character(data[, 193])
attributes(data)$variable.labels[193] <- "[Processing] 2.4 Was the (re-)organization of downstream value chain levels (levels after primary production) successful?"
data[, 193] <- factor(data[, 193], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[193] <- "B4_SQ04"
# LimeSurvey Field type: A
data[, 194] <- as.character(data[, 194])
attributes(data)$variable.labels[194] <- "[Marketing] 2.4 Was the (re-)organization of downstream value chain levels (levels after primary production) successful?"
data[, 194] <- factor(data[, 194], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[194] <- "B4_SQ05"
# LimeSurvey Field type: A
data[, 195] <- as.character(data[, 195])
attributes(data)$variable.labels[195] <- "[Sales] 2.4 Was the (re-)organization of downstream value chain levels (levels after primary production) successful?"
data[, 195] <- factor(data[, 195], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[195] <- "B4_SQ06"
# LimeSurvey Field type: A
data[, 196] <- as.character(data[, 196])
attributes(data)$variable.labels[196] <- "[{438634X176X2831SQ0002.shown}] 2.4 Was the (re-)organization of downstream value chain levels (levels after primary production) successful?"
data[, 196] <- factor(data[, 196], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[196] <- "B4_SQ08"
# LimeSurvey Field type: A
data[, 197] <- as.character(data[, 197])
attributes(data)$variable.labels[197] <- "[Internal communication among participants (e.g. sharing experiences)] 2.5 In your opinion, were the applied communication approaches successful?"
data[, 197] <- factor(data[, 197], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[197] <- "B5_SQ001"
# LimeSurvey Field type: A
data[, 198] <- as.character(data[, 198])
attributes(data)$variable.labels[198] <- "[External communication for professionals (e.g. publishing results)] 2.5 In your opinion, were the applied communication approaches successful?"
data[, 198] <- factor(data[, 198], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[198] <- "B5_SQ002"
# LimeSurvey Field type: A
data[, 199] <- as.character(data[, 199])
attributes(data)$variable.labels[199] <- "[External communication for the public (e.g. media appearances)] 2.5 In your opinion, were the applied communication approaches successful?"
data[, 199] <- factor(data[, 199], levels=c("L001","L002","L003","L004","L005"),labels=c("not at all successful", "slightly successful", "moderately successful", "successful", "overwhelmingly successful"))
names(data)[199] <- "B5_SQ003"
# LimeSurvey Field type: A
data[, 200] <- as.character(data[, 200])
attributes(data)$variable.labels[200] <- "[Professional expertise of involved actors] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 200] <- factor(data[, 200], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[200] <- "B6_1#0"
# LimeSurvey Field type: A
data[, 201] <- as.character(data[, 201])
attributes(data)$variable.labels[201] <- "[Professional expertise of involved actors] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 201] <- factor(data[, 201], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[201] <- "B6_1#1"
# LimeSurvey Field type: A
data[, 202] <- as.character(data[, 202])
attributes(data)$variable.labels[202] <- "[Technical solutions/tools applied] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 202] <- factor(data[, 202], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[202] <- "B6_2#0"
# LimeSurvey Field type: A
data[, 203] <- as.character(data[, 203])
attributes(data)$variable.labels[203] <- "[Technical solutions/tools applied] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 203] <- factor(data[, 203], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[203] <- "B6_2#1"
# LimeSurvey Field type: A
data[, 204] <- as.character(data[, 204])
attributes(data)$variable.labels[204] <- "[Availability of inputs (including seeds)] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 204] <- factor(data[, 204], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[204] <- "B6_3#0"
# LimeSurvey Field type: A
data[, 205] <- as.character(data[, 205])
attributes(data)$variable.labels[205] <- "[Availability of inputs (including seeds)] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 205] <- factor(data[, 205], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[205] <- "B6_3#1"
# LimeSurvey Field type: A
data[, 206] <- as.character(data[, 206])
attributes(data)$variable.labels[206] <- "[Market conditions] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 206] <- factor(data[, 206], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[206] <- "B6_4#0"
# LimeSurvey Field type: A
data[, 207] <- as.character(data[, 207])
attributes(data)$variable.labels[207] <- "[Market conditions] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 207] <- factor(data[, 207], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[207] <- "B6_4#1"
# LimeSurvey Field type: A
data[, 208] <- as.character(data[, 208])
attributes(data)$variable.labels[208] <- "[Amount of financial resources committed] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 208] <- factor(data[, 208], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[208] <- "B6_5#0"
# LimeSurvey Field type: A
data[, 209] <- as.character(data[, 209])
attributes(data)$variable.labels[209] <- "[Amount of financial resources committed] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 209] <- factor(data[, 209], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[209] <- "B6_5#1"
# LimeSurvey Field type: A
data[, 210] <- as.character(data[, 210])
attributes(data)$variable.labels[210] <- "[Available timeframe of the initiative ] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 210] <- factor(data[, 210], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[210] <- "B6_6#0"
# LimeSurvey Field type: A
data[, 211] <- as.character(data[, 211])
attributes(data)$variable.labels[211] <- "[Available timeframe of the initiative ] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 211] <- factor(data[, 211], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[211] <- "B6_6#1"
# LimeSurvey Field type: A
data[, 212] <- as.character(data[, 212])
attributes(data)$variable.labels[212] <- "[Commitment of the involved actors ] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 212] <- factor(data[, 212], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[212] <- "B6_7#0"
# LimeSurvey Field type: A
data[, 213] <- as.character(data[, 213])
attributes(data)$variable.labels[213] <- "[Commitment of the involved actors ] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 213] <- factor(data[, 213], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[213] <- "B6_7#1"
# LimeSurvey Field type: A
data[, 214] <- as.character(data[, 214])
attributes(data)$variable.labels[214] <- "[Organization of actors within the diversification initiative] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 214] <- factor(data[, 214], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[214] <- "B6_8#0"
# LimeSurvey Field type: A
data[, 215] <- as.character(data[, 215])
attributes(data)$variable.labels[215] <- "[Organization of actors within the diversification initiative] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 215] <- factor(data[, 215], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[215] <- "B6_8#1"
# LimeSurvey Field type: A
data[, 216] <- as.character(data[, 216])
attributes(data)$variable.labels[216] <- "[Communication activities] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 216] <- factor(data[, 216], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[216] <- "B6_9#0"
# LimeSurvey Field type: A
data[, 217] <- as.character(data[, 217])
attributes(data)$variable.labels[217] <- "[Communication activities] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 217] <- factor(data[, 217], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[217] <- "B6_9#1"
# LimeSurvey Field type: A
data[, 218] <- as.character(data[, 218])
attributes(data)$variable.labels[218] <- "[General public interest] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 218] <- factor(data[, 218], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[218] <- "B6_10#0"
# LimeSurvey Field type: A
data[, 219] <- as.character(data[, 219])
attributes(data)$variable.labels[219] <- "[General public interest] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 219] <- factor(data[, 219], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[219] <- "B6_10#1"
# LimeSurvey Field type: A
data[, 220] <- as.character(data[, 220])
attributes(data)$variable.labels[220] <- "[Other] [Scale 1] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 220] <- factor(data[, 220], levels=c("A1","A2","A3","A4","A5"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[220] <- "B6_11#0"
# LimeSurvey Field type: A
data[, 221] <- as.character(data[, 221])
attributes(data)$variable.labels[221] <- "[Other] [Scale 2] 	2.6 Which factors have contributed to the success respectively failure of your diversification initiative?"
data[, 221] <- factor(data[, 221], levels=c("A6","A7","A8","A9","A10"),labels=c("not contributed at all", "slighty contributed", "moderately contributed", "strongly contributed", "very strongly contributed"))
names(data)[221] <- "B6_11#1"
# LimeSurvey Field type: A
data[, 222] <- as.character(data[, 222])
attributes(data)$variable.labels[222] <- "2.6.1 You have rated a factor contributing to the failure of the diversification initiative in the \\'Other\\' category of question 2.6. Please specify what this factor is."
names(data)[222] <- "B61"
# LimeSurvey Field type: A
data[, 223] <- as.character(data[, 223])
attributes(data)$variable.labels[223] <- "2.6.1 You have rated a factor contributing to the success of the diversification initiative in the \\'Other\\' category of question 2.6. Please specify what this factor is."
names(data)[223] <- "B62"
# LimeSurvey Field type: F
data[, 224] <- as.numeric(data[, 224])
attributes(data)$variable.labels[224] <- "[Agronomic] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 224] <- factor(data[, 224], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[224] <- "B63_SQ001"
# LimeSurvey Field type: F
data[, 225] <- as.numeric(data[, 225])
attributes(data)$variable.labels[225] <- "[Economic] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 225] <- factor(data[, 225], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[225] <- "B63_SQ002"
# LimeSurvey Field type: F
data[, 226] <- as.numeric(data[, 226])
attributes(data)$variable.labels[226] <- "[Policy] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 226] <- factor(data[, 226], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[226] <- "B63_SQ003"
# LimeSurvey Field type: F
data[, 227] <- as.numeric(data[, 227])
attributes(data)$variable.labels[227] <- "[Social sciences] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 227] <- factor(data[, 227], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[227] <- "B63_SQ004"
# LimeSurvey Field type: F
data[, 228] <- as.numeric(data[, 228])
attributes(data)$variable.labels[228] <- "[Natural sciences] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 228] <- factor(data[, 228], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[228] <- "B63_SQ005"
# LimeSurvey Field type: F
data[, 229] <- as.numeric(data[, 229])
attributes(data)$variable.labels[229] <- "[Marketing] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 229] <- factor(data[, 229], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[229] <- "B63_SQ006"
# LimeSurvey Field type: F
data[, 230] <- as.numeric(data[, 230])
attributes(data)$variable.labels[230] <- "[Sales] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 230] <- factor(data[, 230], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[230] <- "B63_SQ007"
# LimeSurvey Field type: F
data[, 231] <- as.numeric(data[, 231])
attributes(data)$variable.labels[231] <- "[Communication] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 231] <- factor(data[, 231], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[231] <- "B63_SQ008"
# LimeSurvey Field type: A
data[, 232] <- as.character(data[, 232])
attributes(data)$variable.labels[232] <- "[Other] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please specify which expertise is/was most relevant?"
names(data)[232] <- "B63_other"
# LimeSurvey Field type: A
data[, 233] <- as.character(data[, 233])
attributes(data)$variable.labels[233] <- "2.6.4 You have rated technical solutions/tools applied as a relevant factor contributing to the failure of your diversification initiative. Please specify which tools were/are most relevant in your case?"
names(data)[233] <- "B64"
# LimeSurvey Field type: F
data[, 234] <- as.numeric(data[, 234])
attributes(data)$variable.labels[234] <- "[Seeds] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the failure of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
data[, 234] <- factor(data[, 234], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[234] <- "B65_SQ001"
# LimeSurvey Field type: F
data[, 235] <- as.numeric(data[, 235])
attributes(data)$variable.labels[235] <- "[Fertilizers] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the failure of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
data[, 235] <- factor(data[, 235], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[235] <- "B65_SQ002"
# LimeSurvey Field type: F
data[, 236] <- as.numeric(data[, 236])
attributes(data)$variable.labels[236] <- "[Herbicides] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the failure of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
data[, 236] <- factor(data[, 236], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[236] <- "B65_SQ003"
# LimeSurvey Field type: F
data[, 237] <- as.numeric(data[, 237])
attributes(data)$variable.labels[237] <- "[Pesticides] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the failure of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
data[, 237] <- factor(data[, 237], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[237] <- "B65_SQ004"
# LimeSurvey Field type: F
data[, 238] <- as.numeric(data[, 238])
attributes(data)$variable.labels[238] <- "[Plant growth promoters] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the failure of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
data[, 238] <- factor(data[, 238], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[238] <- "B65_SQ005"
# LimeSurvey Field type: F
data[, 239] <- as.numeric(data[, 239])
attributes(data)$variable.labels[239] <- "[Other] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the failure of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
names(data)[239] <- "B65_other"
# LimeSurvey Field type: A
data[, 240] <- as.character(data[, 240])
attributes(data)$variable.labels[240] <- "2.6.6 You have rated market conditions as a relevant factor contributing to the failure of your diversification initiative. Please specify what were/are the most relevant conditions in your case?"
names(data)[240] <- "B66"
# LimeSurvey Field type: A
data[, 241] <- as.character(data[, 241])
attributes(data)$variable.labels[241] <- "2.6.7 You have rated committed financial resources as a relevant factor contributing to the failure of your diversification initiative. Please specify how the funding influenced your initiative?"
names(data)[241] <- "B67"
# LimeSurvey Field type: A
data[, 242] <- as.character(data[, 242])
attributes(data)$variable.labels[242] <- "[Farmers] 2.6.8 You have rated commitment of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please indicate the commitment of involved actors."
data[, 242] <- factor(data[, 242], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[242] <- "B68_SQ001"
# LimeSurvey Field type: A
data[, 243] <- as.character(data[, 243])
attributes(data)$variable.labels[243] <- "[Advisors] 2.6.8 You have rated commitment of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please indicate the commitment of involved actors."
data[, 243] <- factor(data[, 243], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[243] <- "B68_SQ002"
# LimeSurvey Field type: A
data[, 244] <- as.character(data[, 244])
attributes(data)$variable.labels[244] <- "[Researchers] 2.6.8 You have rated commitment of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please indicate the commitment of involved actors."
data[, 244] <- factor(data[, 244], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[244] <- "B68_SQ003"
# LimeSurvey Field type: A
data[, 245] <- as.character(data[, 245])
attributes(data)$variable.labels[245] <- "[Commercial companies] 2.6.8 You have rated commitment of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please indicate the commitment of involved actors."
data[, 245] <- factor(data[, 245], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[245] <- "B68_SQ004"
# LimeSurvey Field type: A
data[, 246] <- as.character(data[, 246])
attributes(data)$variable.labels[246] <- "[Consumers] 2.6.8 You have rated commitment of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please indicate the commitment of involved actors."
data[, 246] <- factor(data[, 246], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[246] <- "B68_SQ005"
# LimeSurvey Field type: A
data[, 247] <- as.character(data[, 247])
attributes(data)$variable.labels[247] <- "[Authorities] 2.6.8 You have rated commitment of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please indicate the commitment of involved actors."
data[, 247] <- factor(data[, 247], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[247] <- "B68_SQ006"
# LimeSurvey Field type: A
data[, 248] <- as.character(data[, 248])
attributes(data)$variable.labels[248] <- "[Certification organizations] 2.6.8 You have rated commitment of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please indicate the commitment of involved actors."
data[, 248] <- factor(data[, 248], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[248] <- "B68_SQ007"
# LimeSurvey Field type: A
data[, 249] <- as.character(data[, 249])
attributes(data)$variable.labels[249] <- "[NGOs] 2.6.8 You have rated commitment of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please indicate the commitment of involved actors."
data[, 249] <- factor(data[, 249], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[249] <- "B68_SQ008"
# LimeSurvey Field type: A
data[, 250] <- as.character(data[, 250])
attributes(data)$variable.labels[250] <- "[{438634X176X2730SQ0001.shown}] 2.6.8 You have rated commitment of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please indicate the commitment of involved actors."
data[, 250] <- factor(data[, 250], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[250] <- "B68_SQ0001"
# LimeSurvey Field type: A
data[, 251] <- as.character(data[, 251])
attributes(data)$variable.labels[251] <- "2.6.10 You have rated internal organization of involved actors as a relevant factor contributing to the failure of your diversification initiative. Please specify how the actors were organized?"
names(data)[251] <- "B610"
# LimeSurvey Field type: A
data[, 252] <- as.character(data[, 252])
attributes(data)$variable.labels[252] <- "2.6.11 You have rated general public interest as a relevant factor contributing to the failure of your diversification initiative. Please specify in what way it influenced your initiative?"
names(data)[252] <- "B611"
# LimeSurvey Field type: F
data[, 253] <- as.numeric(data[, 253])
attributes(data)$variable.labels[253] <- "[Agronomic] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the success of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 253] <- factor(data[, 253], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[253] <- "B63s_SQ001"
# LimeSurvey Field type: F
data[, 254] <- as.numeric(data[, 254])
attributes(data)$variable.labels[254] <- "[Economic] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the success of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 254] <- factor(data[, 254], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[254] <- "B63s_SQ002"
# LimeSurvey Field type: F
data[, 255] <- as.numeric(data[, 255])
attributes(data)$variable.labels[255] <- "[Policy] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the success of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 255] <- factor(data[, 255], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[255] <- "B63s_SQ003"
# LimeSurvey Field type: F
data[, 256] <- as.numeric(data[, 256])
attributes(data)$variable.labels[256] <- "[Social sciences] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the success of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 256] <- factor(data[, 256], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[256] <- "B63s_SQ004"
# LimeSurvey Field type: F
data[, 257] <- as.numeric(data[, 257])
attributes(data)$variable.labels[257] <- "[Natural sciences] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the success of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 257] <- factor(data[, 257], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[257] <- "B63s_SQ005"
# LimeSurvey Field type: F
data[, 258] <- as.numeric(data[, 258])
attributes(data)$variable.labels[258] <- "[Marketing] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the success of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 258] <- factor(data[, 258], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[258] <- "B63s_SQ006"
# LimeSurvey Field type: F
data[, 259] <- as.numeric(data[, 259])
attributes(data)$variable.labels[259] <- "[Sales] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the success of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 259] <- factor(data[, 259], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[259] <- "B63s_SQ007"
# LimeSurvey Field type: F
data[, 260] <- as.numeric(data[, 260])
attributes(data)$variable.labels[260] <- "[Communication] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the success of your diversification initiative. Please specify which expertise is/was most relevant?"
data[, 260] <- factor(data[, 260], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[260] <- "B63s_SQ008"
# LimeSurvey Field type: A
data[, 261] <- as.character(data[, 261])
attributes(data)$variable.labels[261] <- "[Other] 2.6.3 You have rated professional expertise of involved actors as a relevant factor contributing to the success of your diversification initiative. Please specify which expertise is/was most relevant?"
names(data)[261] <- "B63s_other"
# LimeSurvey Field type: A
data[, 262] <- as.character(data[, 262])
attributes(data)$variable.labels[262] <- "2.6.4 You have rated technical solutions/tools applied as a relevant factor contributing to the success of your diversification initiative. Please specify which tools were/are most relevant in your case?"
names(data)[262] <- "B64s"
# LimeSurvey Field type: F
data[, 263] <- as.numeric(data[, 263])
attributes(data)$variable.labels[263] <- "[Seeds] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the success of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
data[, 263] <- factor(data[, 263], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[263] <- "B65s_SQ001"
# LimeSurvey Field type: F
data[, 264] <- as.numeric(data[, 264])
attributes(data)$variable.labels[264] <- "[Fertilizers] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the success of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
data[, 264] <- factor(data[, 264], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[264] <- "B65s_SQ002"
# LimeSurvey Field type: F
data[, 265] <- as.numeric(data[, 265])
attributes(data)$variable.labels[265] <- "[Herbicides] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the success of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
data[, 265] <- factor(data[, 265], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[265] <- "B65s_SQ003"
# LimeSurvey Field type: F
data[, 266] <- as.numeric(data[, 266])
attributes(data)$variable.labels[266] <- "[Pesticides] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the success of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
data[, 266] <- factor(data[, 266], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[266] <- "B65s_SQ004"
# LimeSurvey Field type: F
data[, 267] <- as.numeric(data[, 267])
attributes(data)$variable.labels[267] <- "[Plant growth promoters] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the success of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
data[, 267] <- factor(data[, 267], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[267] <- "B65s_SQ005"
# LimeSurvey Field type: A
data[, 268] <- as.character(data[, 268])
attributes(data)$variable.labels[268] <- "[Other] 2.6.5 You have rated availability of inputs (including seeds) as a relevant factor contributing to the success of your diversification initiative. Please specify which input availabilities were most relevant in your case?"
names(data)[268] <- "B65s_other"
# LimeSurvey Field type: A
data[, 269] <- as.character(data[, 269])
attributes(data)$variable.labels[269] <- "2.6.6 You have rated market conditions as a relevant factor contributing to the success of your diversification initiative. Please specify what were/are the most relevant conditions in your case?"
names(data)[269] <- "B66s"
# LimeSurvey Field type: A
data[, 270] <- as.character(data[, 270])
attributes(data)$variable.labels[270] <- "2.6.7 You have rated committed financial resources as a relevant factor contributing to the success of your diversification initiative. Please specify how the funding influenced your initiative?"
names(data)[270] <- "B67s"
# LimeSurvey Field type: A
data[, 271] <- as.character(data[, 271])
attributes(data)$variable.labels[271] <- "[Farmers] 2.6.9 You have rated commitment of involved actors as a relevant factor contributing to the success of your diversification initiative. Please indicate the commitment of participants."
data[, 271] <- factor(data[, 271], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[271] <- "B69_SQ001"
# LimeSurvey Field type: A
data[, 272] <- as.character(data[, 272])
attributes(data)$variable.labels[272] <- "[Advisors] 2.6.9 You have rated commitment of involved actors as a relevant factor contributing to the success of your diversification initiative. Please indicate the commitment of participants."
data[, 272] <- factor(data[, 272], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[272] <- "B69_SQ002"
# LimeSurvey Field type: A
data[, 273] <- as.character(data[, 273])
attributes(data)$variable.labels[273] <- "[Researchers] 2.6.9 You have rated commitment of involved actors as a relevant factor contributing to the success of your diversification initiative. Please indicate the commitment of participants."
data[, 273] <- factor(data[, 273], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[273] <- "B69_SQ003"
# LimeSurvey Field type: A
data[, 274] <- as.character(data[, 274])
attributes(data)$variable.labels[274] <- "[Commercial companies] 2.6.9 You have rated commitment of involved actors as a relevant factor contributing to the success of your diversification initiative. Please indicate the commitment of participants."
data[, 274] <- factor(data[, 274], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[274] <- "B69_SQ004"
# LimeSurvey Field type: A
data[, 275] <- as.character(data[, 275])
attributes(data)$variable.labels[275] <- "[Consumers] 2.6.9 You have rated commitment of involved actors as a relevant factor contributing to the success of your diversification initiative. Please indicate the commitment of participants."
data[, 275] <- factor(data[, 275], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[275] <- "B69_SQ005"
# LimeSurvey Field type: A
data[, 276] <- as.character(data[, 276])
attributes(data)$variable.labels[276] <- "[Authorities] 2.6.9 You have rated commitment of involved actors as a relevant factor contributing to the success of your diversification initiative. Please indicate the commitment of participants."
data[, 276] <- factor(data[, 276], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[276] <- "B69_SQ006"
# LimeSurvey Field type: A
data[, 277] <- as.character(data[, 277])
attributes(data)$variable.labels[277] <- "[Certification organizations] 2.6.9 You have rated commitment of involved actors as a relevant factor contributing to the success of your diversification initiative. Please indicate the commitment of participants."
data[, 277] <- factor(data[, 277], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[277] <- "B69_SQ007"
# LimeSurvey Field type: A
data[, 278] <- as.character(data[, 278])
attributes(data)$variable.labels[278] <- "[NGOs] 2.6.9 You have rated commitment of involved actors as a relevant factor contributing to the success of your diversification initiative. Please indicate the commitment of participants."
data[, 278] <- factor(data[, 278], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[278] <- "B69_SQ008"
# LimeSurvey Field type: A
data[, 279] <- as.character(data[, 279])
attributes(data)$variable.labels[279] <- "[{438634X176X2730SQ0001.shown}] 2.6.9 You have rated commitment of involved actors as a relevant factor contributing to the success of your diversification initiative. Please indicate the commitment of participants."
data[, 279] <- factor(data[, 279], levels=c("L001","L002","L003","L004","L005"),labels=c("Not at all commited", "Slightly commited", "Moderately commited", "Very commited", "Overwhelmingly commited"))
names(data)[279] <- "B69_SQ0001"
# LimeSurvey Field type: A
data[, 280] <- as.character(data[, 280])
attributes(data)$variable.labels[280] <- "2.6.10 You have rated internal organization of involved actors as a relevant factor contributing to the success of your diversification initiative. Please specify how the actors were organized?"
names(data)[280] <- "B610s"
# LimeSurvey Field type: A
data[, 281] <- as.character(data[, 281])
attributes(data)$variable.labels[281] <- "2.6.11 You have rated general public interest as a relevant factor contributing to the success of your diversification initiative. Please specify in what way it influenced your initiative?"
names(data)[281] <- "B611s"
# LimeSurvey Field type: F
data[, 282] <- as.numeric(data[, 282])
attributes(data)$variable.labels[282] <- "In the following picture the main steps of the life cycle of a project are modelled. In which state is your crop diversification initiative? Please indicate the number corresponding to this state, a number between 1 and 9."
names(data)[282] <- "C1"
# LimeSurvey Field type: F
data[, 283] <- as.numeric(data[, 283])
attributes(data)$variable.labels[283] <- "[No, it remained the same (please type 0 in the comment box)] 3.2 Did the size of the crop diversification area change during the lifetime of the diversification initiative?"
data[, 283] <- factor(data[, 283], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[283] <- "C2_SQ002"
# LimeSurvey Field type: A
data[, 284] <- as.character(data[, 284])
attributes(data)$variable.labels[284] <- "[Comment] 3.2 Did the size of the crop diversification area change during the lifetime of the diversification initiative?"
names(data)[284] <- "C2_SQ002comment"
# LimeSurvey Field type: F
data[, 285] <- as.numeric(data[, 285])
attributes(data)$variable.labels[285] <- "[Yes, it increased by app. … % (please specify)] 3.2 Did the size of the crop diversification area change during the lifetime of the diversification initiative?"
data[, 285] <- factor(data[, 285], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[285] <- "C2_SQ003"
# LimeSurvey Field type: A
data[, 286] <- as.character(data[, 286])
attributes(data)$variable.labels[286] <- "[Comment] 3.2 Did the size of the crop diversification area change during the lifetime of the diversification initiative?"
names(data)[286] <- "C2_SQ003comment"
# LimeSurvey Field type: F
data[, 287] <- as.numeric(data[, 287])
attributes(data)$variable.labels[287] <- "[Yes, it decreased by app. … % (please specify)] 3.2 Did the size of the crop diversification area change during the lifetime of the diversification initiative?"
data[, 287] <- factor(data[, 287], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[287] <- "C2_SQ004"
# LimeSurvey Field type: A
data[, 288] <- as.character(data[, 288])
attributes(data)$variable.labels[288] <- "[Comment] 3.2 Did the size of the crop diversification area change during the lifetime of the diversification initiative?"
names(data)[288] <- "C2_SQ004comment"
# LimeSurvey Field type: A
data[, 289] <- as.character(data[, 289])
attributes(data)$variable.labels[289] <- "[] [Year] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[289] <- "C3_SQ001_SQ001"
# LimeSurvey Field type: A
data[, 290] <- as.character(data[, 290])
attributes(data)$variable.labels[290] <- "[] [Number of participating farms] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[290] <- "C3_SQ001_SQ002"
# LimeSurvey Field type: A
data[, 291] <- as.character(data[, 291])
attributes(data)$variable.labels[291] <- "[] [Approximate acreage of diversification initiative (ha)] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[291] <- "C3_SQ001_SQ003"
# LimeSurvey Field type: A
data[, 292] <- as.character(data[, 292])
attributes(data)$variable.labels[292] <- "[] [Year] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[292] <- "C3_SQ002_SQ001"
# LimeSurvey Field type: A
data[, 293] <- as.character(data[, 293])
attributes(data)$variable.labels[293] <- "[] [Number of participating farms] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[293] <- "C3_SQ002_SQ002"
# LimeSurvey Field type: A
data[, 294] <- as.character(data[, 294])
attributes(data)$variable.labels[294] <- "[] [Approximate acreage of diversification initiative (ha)] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[294] <- "C3_SQ002_SQ003"
# LimeSurvey Field type: A
data[, 295] <- as.character(data[, 295])
attributes(data)$variable.labels[295] <- "[] [Year] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[295] <- "C3_SQ003_SQ001"
# LimeSurvey Field type: A
data[, 296] <- as.character(data[, 296])
attributes(data)$variable.labels[296] <- "[] [Number of participating farms] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[296] <- "C3_SQ003_SQ002"
# LimeSurvey Field type: A
data[, 297] <- as.character(data[, 297])
attributes(data)$variable.labels[297] <- "[] [Approximate acreage of diversification initiative (ha)] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[297] <- "C3_SQ003_SQ003"
# LimeSurvey Field type: A
data[, 298] <- as.character(data[, 298])
attributes(data)$variable.labels[298] <- "[] [Year] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[298] <- "C3_SQ004_SQ001"
# LimeSurvey Field type: A
data[, 299] <- as.character(data[, 299])
attributes(data)$variable.labels[299] <- "[] [Number of participating farms] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[299] <- "C3_SQ004_SQ002"
# LimeSurvey Field type: A
data[, 300] <- as.character(data[, 300])
attributes(data)$variable.labels[300] <- "[] [Approximate acreage of diversification initiative (ha)] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[300] <- "C3_SQ004_SQ003"
# LimeSurvey Field type: A
data[, 301] <- as.character(data[, 301])
attributes(data)$variable.labels[301] <- "[] [Year] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[301] <- "C3_SQ005_SQ001"
# LimeSurvey Field type: A
data[, 302] <- as.character(data[, 302])
attributes(data)$variable.labels[302] <- "[] [Number of participating farms] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[302] <- "C3_SQ005_SQ002"
# LimeSurvey Field type: A
data[, 303] <- as.character(data[, 303])
attributes(data)$variable.labels[303] <- "[] [Approximate acreage of diversification initiative (ha)] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[303] <- "C3_SQ005_SQ003"
# LimeSurvey Field type: A
data[, 304] <- as.character(data[, 304])
attributes(data)$variable.labels[304] <- "[] [Year] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[304] <- "C3_SQ006_SQ001"
# LimeSurvey Field type: A
data[, 305] <- as.character(data[, 305])
attributes(data)$variable.labels[305] <- "[] [Number of participating farms] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[305] <- "C3_SQ006_SQ002"
# LimeSurvey Field type: A
data[, 306] <- as.character(data[, 306])
attributes(data)$variable.labels[306] <- "[] [Approximate acreage of diversification initiative (ha)] 3.3 How many farmers take/took part in the diversification initiative? Please specify the dynamics according to years. You may also define intervals, if the project was longer than six years."
names(data)[306] <- "C3_SQ006_SQ003"
# LimeSurvey Field type: A
data[, 307] <- as.character(data[, 307])
attributes(data)$variable.labels[307] <- "[Agronomic (e.g. water availability)] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 307] <- factor(data[, 307], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[307] <- "C4_1#0"
# LimeSurvey Field type: A
data[, 308] <- as.character(data[, 308])
attributes(data)$variable.labels[308] <- "[Agronomic (e.g. water availability)] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 308] <- factor(data[, 308], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[308] <- "C4_1#1"
# LimeSurvey Field type: A
data[, 309] <- as.character(data[, 309])
attributes(data)$variable.labels[309] <- "[Economic (e.g. product price)] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 309] <- factor(data[, 309], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[309] <- "C4_2#0"
# LimeSurvey Field type: A
data[, 310] <- as.character(data[, 310])
attributes(data)$variable.labels[310] <- "[Economic (e.g. product price)] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 310] <- factor(data[, 310], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[310] <- "C4_2#1"
# LimeSurvey Field type: A
data[, 311] <- as.character(data[, 311])
attributes(data)$variable.labels[311] <- "[Public policy (e.g. regulations)] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 311] <- factor(data[, 311], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[311] <- "C4_3#0"
# LimeSurvey Field type: A
data[, 312] <- as.character(data[, 312])
attributes(data)$variable.labels[312] <- "[Public policy (e.g. regulations)] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 312] <- factor(data[, 312], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[312] <- "C4_3#1"
# LimeSurvey Field type: A
data[, 313] <- as.character(data[, 313])
attributes(data)$variable.labels[313] <- "[Personal interactions (e.g. team work)] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 313] <- factor(data[, 313], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[313] <- "C4_4#0"
# LimeSurvey Field type: A
data[, 314] <- as.character(data[, 314])
attributes(data)$variable.labels[314] <- "[Personal interactions (e.g. team work)] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 314] <- factor(data[, 314], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[314] <- "C4_4#1"
# LimeSurvey Field type: A
data[, 315] <- as.character(data[, 315])
attributes(data)$variable.labels[315] <- "[Breeding (e.g. new varieties)] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 315] <- factor(data[, 315], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[315] <- "C4_SQ001#0"
# LimeSurvey Field type: A
data[, 316] <- as.character(data[, 316])
attributes(data)$variable.labels[316] <- "[Breeding (e.g. new varieties)] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 316] <- factor(data[, 316], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[316] <- "C4_SQ001#1"
# LimeSurvey Field type: A
data[, 317] <- as.character(data[, 317])
attributes(data)$variable.labels[317] <- "[Seed production (e.g. seed availability)] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 317] <- factor(data[, 317], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[317] <- "C4_SQ002#0"
# LimeSurvey Field type: A
data[, 318] <- as.character(data[, 318])
attributes(data)$variable.labels[318] <- "[Seed production (e.g. seed availability)] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 318] <- factor(data[, 318], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[318] <- "C4_SQ002#1"
# LimeSurvey Field type: A
data[, 319] <- as.character(data[, 319])
attributes(data)$variable.labels[319] <- "[Development of inputs needed for production] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 319] <- factor(data[, 319], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[319] <- "C4_SQ003#0"
# LimeSurvey Field type: A
data[, 320] <- as.character(data[, 320])
attributes(data)$variable.labels[320] <- "[Development of inputs needed for production] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 320] <- factor(data[, 320], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[320] <- "C4_SQ003#1"
# LimeSurvey Field type: A
data[, 321] <- as.character(data[, 321])
attributes(data)$variable.labels[321] <- "[Organization of inputs needed for production] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 321] <- factor(data[, 321], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[321] <- "C4_SQ004#0"
# LimeSurvey Field type: A
data[, 322] <- as.character(data[, 322])
attributes(data)$variable.labels[322] <- "[Organization of inputs needed for production] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 322] <- factor(data[, 322], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[322] <- "C4_SQ004#1"
# LimeSurvey Field type: A
data[, 323] <- as.character(data[, 323])
attributes(data)$variable.labels[323] <- "[Machinery development] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 323] <- factor(data[, 323], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[323] <- "C4_SQ005#0"
# LimeSurvey Field type: A
data[, 324] <- as.character(data[, 324])
attributes(data)$variable.labels[324] <- "[Machinery development] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 324] <- factor(data[, 324], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[324] <- "C4_SQ005#1"
# LimeSurvey Field type: A
data[, 325] <- as.character(data[, 325])
attributes(data)$variable.labels[325] <- "[Organization of machinery needed for production] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 325] <- factor(data[, 325], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[325] <- "C4_SQ006#0"
# LimeSurvey Field type: A
data[, 326] <- as.character(data[, 326])
attributes(data)$variable.labels[326] <- "[Organization of machinery needed for production] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 326] <- factor(data[, 326], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[326] <- "C4_SQ006#1"
# LimeSurvey Field type: A
data[, 327] <- as.character(data[, 327])
attributes(data)$variable.labels[327] <- "[Quality assurance] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 327] <- factor(data[, 327], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[327] <- "C4_SQ01#0"
# LimeSurvey Field type: A
data[, 328] <- as.character(data[, 328])
attributes(data)$variable.labels[328] <- "[Quality assurance] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 328] <- factor(data[, 328], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[328] <- "C4_SQ01#1"
# LimeSurvey Field type: A
data[, 329] <- as.character(data[, 329])
attributes(data)$variable.labels[329] <- "[Transportation] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 329] <- factor(data[, 329], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[329] <- "C4_SQ02#0"
# LimeSurvey Field type: A
data[, 330] <- as.character(data[, 330])
attributes(data)$variable.labels[330] <- "[Transportation] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 330] <- factor(data[, 330], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[330] <- "C4_SQ02#1"
# LimeSurvey Field type: A
data[, 331] <- as.character(data[, 331])
attributes(data)$variable.labels[331] <- "[Logistics] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 331] <- factor(data[, 331], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[331] <- "C4_SQ03#0"
# LimeSurvey Field type: A
data[, 332] <- as.character(data[, 332])
attributes(data)$variable.labels[332] <- "[Logistics] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 332] <- factor(data[, 332], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[332] <- "C4_SQ03#1"
# LimeSurvey Field type: A
data[, 333] <- as.character(data[, 333])
attributes(data)$variable.labels[333] <- "[Processing] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 333] <- factor(data[, 333], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[333] <- "C4_SQ04#0"
# LimeSurvey Field type: A
data[, 334] <- as.character(data[, 334])
attributes(data)$variable.labels[334] <- "[Processing] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 334] <- factor(data[, 334], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[334] <- "C4_SQ04#1"
# LimeSurvey Field type: A
data[, 335] <- as.character(data[, 335])
attributes(data)$variable.labels[335] <- "[Marketing] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 335] <- factor(data[, 335], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[335] <- "C4_SQ05#0"
# LimeSurvey Field type: A
data[, 336] <- as.character(data[, 336])
attributes(data)$variable.labels[336] <- "[Marketing] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 336] <- factor(data[, 336], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[336] <- "C4_SQ05#1"
# LimeSurvey Field type: A
data[, 337] <- as.character(data[, 337])
attributes(data)$variable.labels[337] <- "[Sales] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 337] <- factor(data[, 337], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[337] <- "C4_SQ06#0"
# LimeSurvey Field type: A
data[, 338] <- as.character(data[, 338])
attributes(data)$variable.labels[338] <- "[Sales] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 338] <- factor(data[, 338], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[338] <- "C4_SQ06#1"
# LimeSurvey Field type: A
data[, 339] <- as.character(data[, 339])
attributes(data)$variable.labels[339] <- "[{438634X176X2820SQ0001.shown}] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 339] <- factor(data[, 339], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[339] <- "C4_SQ0001#0"
# LimeSurvey Field type: A
data[, 340] <- as.character(data[, 340])
attributes(data)$variable.labels[340] <- "[{438634X176X2820SQ0001.shown}] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 340] <- factor(data[, 340], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[340] <- "C4_SQ0001#1"
# LimeSurvey Field type: A
data[, 341] <- as.character(data[, 341])
attributes(data)$variable.labels[341] <- "[{438634X176X2831SQ0002.shown}] [Scale 1] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 341] <- factor(data[, 341], levels=c("A2","A3","A4","A5"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[341] <- "C4_SQ0002#0"
# LimeSurvey Field type: A
data[, 342] <- as.character(data[, 342])
attributes(data)$variable.labels[342] <- "[{438634X176X2831SQ0002.shown}] [Scale 2] 	3.4 Did the diversification initiative encounter any drawbacks or enablers during its lifetime? Please rate the following potential drawbacks/enablers on both scales."
data[, 342] <- factor(data[, 342], levels=c("A7","A8","A9","A10"),labels=c("slightly relevant", "moderately relevant", "very relevant", "overwhelmingly relevant"))
names(data)[342] <- "C4_SQ0002#1"
# LimeSurvey Field type: F
data[, 343] <- as.numeric(data[, 343])
attributes(data)$variable.labels[343] <- "[Water availability (for irrigation)] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 343] <- factor(data[, 343], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[343] <- "C41_1"
# LimeSurvey Field type: F
data[, 344] <- as.numeric(data[, 344])
attributes(data)$variable.labels[344] <- "[Crop protection] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 344] <- factor(data[, 344], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[344] <- "C41_2"
# LimeSurvey Field type: F
data[, 345] <- as.numeric(data[, 345])
attributes(data)$variable.labels[345] <- "[Weed management] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 345] <- factor(data[, 345], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[345] <- "C41_3"
# LimeSurvey Field type: F
data[, 346] <- as.numeric(data[, 346])
attributes(data)$variable.labels[346] <- "[Soil type] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 346] <- factor(data[, 346], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[346] <- "C41_4"
# LimeSurvey Field type: F
data[, 347] <- as.numeric(data[, 347])
attributes(data)$variable.labels[347] <- "[Climatic issues (e.g. drought)] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 347] <- factor(data[, 347], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[347] <- "C41_5"
# LimeSurvey Field type: F
data[, 348] <- as.numeric(data[, 348])
attributes(data)$variable.labels[348] <- "[Machinery] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 348] <- factor(data[, 348], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[348] <- "C41_6"
# LimeSurvey Field type: F
data[, 349] <- as.numeric(data[, 349])
attributes(data)$variable.labels[349] <- "[Quality of agricultural product] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 349] <- factor(data[, 349], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[349] <- "C41_7"
# LimeSurvey Field type: F
data[, 350] <- as.numeric(data[, 350])
attributes(data)$variable.labels[350] <- "[Yield of agricultural product] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 350] <- factor(data[, 350], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[350] <- "C41_8"
# LimeSurvey Field type: F
data[, 351] <- as.numeric(data[, 351])
attributes(data)$variable.labels[351] <- "[Storage] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 351] <- factor(data[, 351], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[351] <- "C41_9"
# LimeSurvey Field type: F
data[, 352] <- as.numeric(data[, 352])
attributes(data)$variable.labels[352] <- "[Expertise available on diversification] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 352] <- factor(data[, 352], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[352] <- "C41_10"
# LimeSurvey Field type: F
data[, 353] <- as.numeric(data[, 353])
attributes(data)$variable.labels[353] <- "[Availability of new crop propagation materials] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 353] <- factor(data[, 353], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[353] <- "C41_11"
# LimeSurvey Field type: F
data[, 354] <- as.numeric(data[, 354])
attributes(data)$variable.labels[354] <- "[Availability of inputs] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 354] <- factor(data[, 354], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[354] <- "C41_12"
# LimeSurvey Field type: A
data[, 355] <- as.character(data[, 355])
attributes(data)$variable.labels[355] <- "[Other] 3.4.1 You have evaluated agronomic drawbacks as very important for your initiative. Please specify which agronomic parameters were most important?"
names(data)[355] <- "C41_other"
# LimeSurvey Field type: F
data[, 356] <- as.numeric(data[, 356])
attributes(data)$variable.labels[356] <- "[Low price of new product(s)] 3.4.2 You have evaluated economic drawbacks as very important for your initiative. Please specify which economic parameters were most important?"
data[, 356] <- factor(data[, 356], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[356] <- "C42_1"
# LimeSurvey Field type: F
data[, 357] <- as.numeric(data[, 357])
attributes(data)$variable.labels[357] <- "[High cultivation costs] 3.4.2 You have evaluated economic drawbacks as very important for your initiative. Please specify which economic parameters were most important?"
data[, 357] <- factor(data[, 357], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[357] <- "C42_2"
# LimeSurvey Field type: F
data[, 358] <- as.numeric(data[, 358])
attributes(data)$variable.labels[358] <- "[Low yield] 3.4.2 You have evaluated economic drawbacks as very important for your initiative. Please specify which economic parameters were most important?"
data[, 358] <- factor(data[, 358], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[358] <- "C42_3"
# LimeSurvey Field type: F
data[, 359] <- as.numeric(data[, 359])
attributes(data)$variable.labels[359] <- "[Higher risk of crop loss] 3.4.2 You have evaluated economic drawbacks as very important for your initiative. Please specify which economic parameters were most important?"
data[, 359] <- factor(data[, 359], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[359] <- "C42_4"
# LimeSurvey Field type: F
data[, 360] <- as.numeric(data[, 360])
attributes(data)$variable.labels[360] <- "[Competition with mainstream producers] 3.4.2 You have evaluated economic drawbacks as very important for your initiative. Please specify which economic parameters were most important?"
data[, 360] <- factor(data[, 360], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[360] <- "C42_5"
# LimeSurvey Field type: F
data[, 361] <- as.numeric(data[, 361])
attributes(data)$variable.labels[361] <- "[Buyers’ focus on single crops] 3.4.2 You have evaluated economic drawbacks as very important for your initiative. Please specify which economic parameters were most important?"
data[, 361] <- factor(data[, 361], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[361] <- "C42_6"
# LimeSurvey Field type: A
data[, 362] <- as.character(data[, 362])
attributes(data)$variable.labels[362] <- "[Other] 3.4.2 You have evaluated economic drawbacks as very important for your initiative. Please specify which economic parameters were most important?"
names(data)[362] <- "C42_other"
# LimeSurvey Field type: F
data[, 363] <- as.numeric(data[, 363])
attributes(data)$variable.labels[363] <- "[CAP Pillar I. - Greening measure] 3.4.3 You have evaluated policy drawbacks as very important for your initiative. Please specify which policy instuments were most important?"
data[, 363] <- factor(data[, 363], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[363] <- "C43_1"
# LimeSurvey Field type: F
data[, 364] <- as.numeric(data[, 364])
attributes(data)$variable.labels[364] <- "[CAP Pillar II - Rural Development Program] 3.4.3 You have evaluated policy drawbacks as very important for your initiative. Please specify which policy instuments were most important?"
data[, 364] <- factor(data[, 364], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[364] <- "C43_2"
# LimeSurvey Field type: F
data[, 365] <- as.numeric(data[, 365])
attributes(data)$variable.labels[365] <- "[Nitrate Directive] 3.4.3 You have evaluated policy drawbacks as very important for your initiative. Please specify which policy instuments were most important?"
data[, 365] <- factor(data[, 365], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[365] <- "C43_3"
# LimeSurvey Field type: F
data[, 366] <- as.numeric(data[, 366])
attributes(data)$variable.labels[366] <- "[Water Framework Directive] 3.4.3 You have evaluated policy drawbacks as very important for your initiative. Please specify which policy instuments were most important?"
data[, 366] <- factor(data[, 366], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[366] <- "C43_4"
# LimeSurvey Field type: F
data[, 367] <- as.numeric(data[, 367])
attributes(data)$variable.labels[367] <- "[Sustainable Use of Pesticides Directive] 3.4.3 You have evaluated policy drawbacks as very important for your initiative. Please specify which policy instuments were most important?"
data[, 367] <- factor(data[, 367], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[367] <- "C43_7"
# LimeSurvey Field type: F
data[, 368] <- as.numeric(data[, 368])
attributes(data)$variable.labels[368] <- "[European Seed Law] 3.4.3 You have evaluated policy drawbacks as very important for your initiative. Please specify which policy instuments were most important?"
data[, 368] <- factor(data[, 368], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[368] <- "C43_6"
# LimeSurvey Field type: A
data[, 369] <- as.character(data[, 369])
attributes(data)$variable.labels[369] <- "[Other] 3.4.3 You have evaluated policy drawbacks as very important for your initiative. Please specify which policy instuments were most important?"
names(data)[369] <- "C43_other"
# LimeSurvey Field type: F
data[, 370] <- as.numeric(data[, 370])
attributes(data)$variable.labels[370] <- "[Knowledge exchange with the actors of the initiative] 3.4.4 You have evaluated personal interactions as very important drawbacks for your initiative. Please specify which personal interactions were most important?"
data[, 370] <- factor(data[, 370], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[370] <- "C44_1"
# LimeSurvey Field type: F
data[, 371] <- as.numeric(data[, 371])
attributes(data)$variable.labels[371] <- "[Public interaction/participation (e.g. field days)] 3.4.4 You have evaluated personal interactions as very important drawbacks for your initiative. Please specify which personal interactions were most important?"
data[, 371] <- factor(data[, 371], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[371] <- "C44_2"
# LimeSurvey Field type: F
data[, 372] <- as.numeric(data[, 372])
attributes(data)$variable.labels[372] <- "[New contacts established through the initiative] 3.4.4 You have evaluated personal interactions as very important drawbacks for your initiative. Please specify which personal interactions were most important?"
data[, 372] <- factor(data[, 372], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[372] <- "C44_3"
# LimeSurvey Field type: F
data[, 373] <- as.numeric(data[, 373])
attributes(data)$variable.labels[373] <- "[New cooperations established through the initiative] 3.4.4 You have evaluated personal interactions as very important drawbacks for your initiative. Please specify which personal interactions were most important?"
data[, 373] <- factor(data[, 373], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[373] <- "C44_4"
# LimeSurvey Field type: F
data[, 374] <- as.numeric(data[, 374])
attributes(data)$variable.labels[374] <- "[New skills acquired through the initiative] 3.4.4 You have evaluated personal interactions as very important drawbacks for your initiative. Please specify which personal interactions were most important?"
data[, 374] <- factor(data[, 374], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[374] <- "C44_5"
# LimeSurvey Field type: A
data[, 375] <- as.character(data[, 375])
attributes(data)$variable.labels[375] <- "[Other] 3.4.4 You have evaluated personal interactions as very important drawbacks for your initiative. Please specify which personal interactions were most important?"
names(data)[375] <- "C44_other"
# LimeSurvey Field type: F
data[, 376] <- as.numeric(data[, 376])
attributes(data)$variable.labels[376] <- "[Water availability (for irrigation)] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 376] <- factor(data[, 376], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[376] <- "C41e_1"
# LimeSurvey Field type: F
data[, 377] <- as.numeric(data[, 377])
attributes(data)$variable.labels[377] <- "[Crop protection] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 377] <- factor(data[, 377], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[377] <- "C41e_2"
# LimeSurvey Field type: F
data[, 378] <- as.numeric(data[, 378])
attributes(data)$variable.labels[378] <- "[Weed management] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 378] <- factor(data[, 378], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[378] <- "C41e_3"
# LimeSurvey Field type: F
data[, 379] <- as.numeric(data[, 379])
attributes(data)$variable.labels[379] <- "[Soil type] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 379] <- factor(data[, 379], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[379] <- "C41e_4"
# LimeSurvey Field type: F
data[, 380] <- as.numeric(data[, 380])
attributes(data)$variable.labels[380] <- "[Climatic issues (e.g. drought)] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 380] <- factor(data[, 380], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[380] <- "C41e_5"
# LimeSurvey Field type: F
data[, 381] <- as.numeric(data[, 381])
attributes(data)$variable.labels[381] <- "[Machinery] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 381] <- factor(data[, 381], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[381] <- "C41e_6"
# LimeSurvey Field type: F
data[, 382] <- as.numeric(data[, 382])
attributes(data)$variable.labels[382] <- "[Quality of agricultural product] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 382] <- factor(data[, 382], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[382] <- "C41e_7"
# LimeSurvey Field type: F
data[, 383] <- as.numeric(data[, 383])
attributes(data)$variable.labels[383] <- "[Yield of agricultural product] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 383] <- factor(data[, 383], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[383] <- "C41e_8"
# LimeSurvey Field type: F
data[, 384] <- as.numeric(data[, 384])
attributes(data)$variable.labels[384] <- "[Storage] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 384] <- factor(data[, 384], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[384] <- "C41e_9"
# LimeSurvey Field type: F
data[, 385] <- as.numeric(data[, 385])
attributes(data)$variable.labels[385] <- "[Expertise available on diversification] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 385] <- factor(data[, 385], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[385] <- "C41e_10"
# LimeSurvey Field type: F
data[, 386] <- as.numeric(data[, 386])
attributes(data)$variable.labels[386] <- "[Availability of new crop propagation materials] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 386] <- factor(data[, 386], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[386] <- "C41e_11"
# LimeSurvey Field type: F
data[, 387] <- as.numeric(data[, 387])
attributes(data)$variable.labels[387] <- "[Availability of inputs] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
data[, 387] <- factor(data[, 387], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[387] <- "C41e_12"
# LimeSurvey Field type: A
data[, 388] <- as.character(data[, 388])
attributes(data)$variable.labels[388] <- "[Other] 3.4.1 You have evaluated agronomic enablers as very important for your initiative. Please specify which agronomic parameters were most important?"
names(data)[388] <- "C41e_other"
# LimeSurvey Field type: F
data[, 389] <- as.numeric(data[, 389])
attributes(data)$variable.labels[389] <- "[Higher price of new product(s)] 3.4.2 You have evaluated economic enablers as very important for your initiative. Please specify which economic parameters were most important?"
data[, 389] <- factor(data[, 389], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[389] <- "C42e_1"
# LimeSurvey Field type: F
data[, 390] <- as.numeric(data[, 390])
attributes(data)$variable.labels[390] <- "[Lower cultivation costs] 3.4.2 You have evaluated economic enablers as very important for your initiative. Please specify which economic parameters were most important?"
data[, 390] <- factor(data[, 390], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[390] <- "C42e_2"
# LimeSurvey Field type: F
data[, 391] <- as.numeric(data[, 391])
attributes(data)$variable.labels[391] <- "[Higher yield] 3.4.2 You have evaluated economic enablers as very important for your initiative. Please specify which economic parameters were most important?"
data[, 391] <- factor(data[, 391], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[391] <- "C42e_3"
# LimeSurvey Field type: F
data[, 392] <- as.numeric(data[, 392])
attributes(data)$variable.labels[392] <- "[Lower risk of crop loss] 3.4.2 You have evaluated economic enablers as very important for your initiative. Please specify which economic parameters were most important?"
data[, 392] <- factor(data[, 392], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[392] <- "C42e_4"
# LimeSurvey Field type: F
data[, 393] <- as.numeric(data[, 393])
attributes(data)$variable.labels[393] <- "[Better competition with mainstream producers] 3.4.2 You have evaluated economic enablers as very important for your initiative. Please specify which economic parameters were most important?"
data[, 393] <- factor(data[, 393], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[393] <- "C42e_5"
# LimeSurvey Field type: F
data[, 394] <- as.numeric(data[, 394])
attributes(data)$variable.labels[394] <- "[Buyers’ willingness to purchase mixed crops] 3.4.2 You have evaluated economic enablers as very important for your initiative. Please specify which economic parameters were most important?"
data[, 394] <- factor(data[, 394], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[394] <- "C42e_6"
# LimeSurvey Field type: A
data[, 395] <- as.character(data[, 395])
attributes(data)$variable.labels[395] <- "[Other] 3.4.2 You have evaluated economic enablers as very important for your initiative. Please specify which economic parameters were most important?"
names(data)[395] <- "C42e_other"
# LimeSurvey Field type: F
data[, 396] <- as.numeric(data[, 396])
attributes(data)$variable.labels[396] <- "[CAP Pillar I - Greening measure] 3.4.3 You have evaluated policy enablers as very important for your initiative. Please specify which policy instuments were most important?"
data[, 396] <- factor(data[, 396], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[396] <- "C43e_SQ001"
# LimeSurvey Field type: F
data[, 397] <- as.numeric(data[, 397])
attributes(data)$variable.labels[397] <- "[CAP Pillar II - Rural Development Program] 3.4.3 You have evaluated policy enablers as very important for your initiative. Please specify which policy instuments were most important?"
data[, 397] <- factor(data[, 397], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[397] <- "C43e_SQ002"
# LimeSurvey Field type: F
data[, 398] <- as.numeric(data[, 398])
attributes(data)$variable.labels[398] <- "[Nitrate Directive] 3.4.3 You have evaluated policy enablers as very important for your initiative. Please specify which policy instuments were most important?"
data[, 398] <- factor(data[, 398], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[398] <- "C43e_SQ003"
# LimeSurvey Field type: F
data[, 399] <- as.numeric(data[, 399])
attributes(data)$variable.labels[399] <- "[Water Framework Directive] 3.4.3 You have evaluated policy enablers as very important for your initiative. Please specify which policy instuments were most important?"
data[, 399] <- factor(data[, 399], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[399] <- "C43e_SQ004"
# LimeSurvey Field type: F
data[, 400] <- as.numeric(data[, 400])
attributes(data)$variable.labels[400] <- "[Sustainable Use of Pesticides Directive] 3.4.3 You have evaluated policy enablers as very important for your initiative. Please specify which policy instuments were most important?"
data[, 400] <- factor(data[, 400], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[400] <- "C43e_SQ005"
# LimeSurvey Field type: F
data[, 401] <- as.numeric(data[, 401])
attributes(data)$variable.labels[401] <- "[European Seed Law] 3.4.3 You have evaluated policy enablers as very important for your initiative. Please specify which policy instuments were most important?"
data[, 401] <- factor(data[, 401], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[401] <- "C43e_SQ006"
# LimeSurvey Field type: A
data[, 402] <- as.character(data[, 402])
attributes(data)$variable.labels[402] <- "[Other] 3.4.3 You have evaluated policy enablers as very important for your initiative. Please specify which policy instuments were most important?"
names(data)[402] <- "C43e_other"
# LimeSurvey Field type: F
data[, 403] <- as.numeric(data[, 403])
attributes(data)$variable.labels[403] <- "[Knowledge exchange with the actors of the initiative] 3.4.4 You have evaluated personal interactions as very important enablers for your initiative. Please specify which personal interactions were most important?"
data[, 403] <- factor(data[, 403], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[403] <- "C44e_1"
# LimeSurvey Field type: F
data[, 404] <- as.numeric(data[, 404])
attributes(data)$variable.labels[404] <- "[Public interaction/participation (e.g. field days)] 3.4.4 You have evaluated personal interactions as very important enablers for your initiative. Please specify which personal interactions were most important?"
data[, 404] <- factor(data[, 404], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[404] <- "C44e_2"
# LimeSurvey Field type: F
data[, 405] <- as.numeric(data[, 405])
attributes(data)$variable.labels[405] <- "[New contacts established through the initiative] 3.4.4 You have evaluated personal interactions as very important enablers for your initiative. Please specify which personal interactions were most important?"
data[, 405] <- factor(data[, 405], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[405] <- "C44e_3"
# LimeSurvey Field type: F
data[, 406] <- as.numeric(data[, 406])
attributes(data)$variable.labels[406] <- "[New cooperations established through the initiative] 3.4.4 You have evaluated personal interactions as very important enablers for your initiative. Please specify which personal interactions were most important?"
data[, 406] <- factor(data[, 406], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[406] <- "C44e_4"
# LimeSurvey Field type: F
data[, 407] <- as.numeric(data[, 407])
attributes(data)$variable.labels[407] <- "[New skills acquired through the initiative] 3.4.4 You have evaluated personal interactions as very important enablers for your initiative. Please specify which personal interactions were most important?"
data[, 407] <- factor(data[, 407], levels=c(1,0),labels=c("Yes", "Not selected"))
names(data)[407] <- "C44e_5"
# LimeSurvey Field type: A
data[, 408] <- as.character(data[, 408])
attributes(data)$variable.labels[408] <- "[Other] 3.4.4 You have evaluated personal interactions as very important enablers for your initiative. Please specify which personal interactions were most important?"
names(data)[408] <- "C44e_other"
# LimeSurvey Field type: A
data[, 409] <- as.character(data[, 409])
attributes(data)$variable.labels[409] <- "3.5 Would you like to share anything else about this diversification initiative?"
names(data)[409] <- "C5"
# LimeSurvey Field type: A
data[, 410] <- as.character(data[, 410])
attributes(data)$variable.labels[410] <- "3.6 Thank you for participating in this survey! If you wish to receive the results of this research, please provide your email address where we can send you the final report."
names(data)[410] <- "C6"
