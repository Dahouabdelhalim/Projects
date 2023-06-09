###Bambach MASS EXTINCTIONS###

#CALCULATES FINAL CMR COEFFICIENTS FROM INTERACTIVE TERMS IN MODELS FOR SUPPLEMENTARY FIGURES

coeffs <- read.delim("CMR_coeffs_bodysize_bambach.txt")

#Bivalves

load("Bivalvia.CMR_bodysize_model_bambach_results.Rdata")
G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[2,4] <- coeffs[1,4] + coeffs[2,4] 
coeffs[2,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))

#bony fish

load("bony fish.CMR_bodysize_model_bambach_results.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[86,88]

coeffs[4,2] <- coeffs[3,2] + coeffs[4,2] 
coeffs[4,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))

G.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[261,2]
G.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[263,2]
G.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[261,263]

coeffs[4,4] <- coeffs[3,4] + coeffs[4,4] 
coeffs[4,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Cephalopods

load("Cephalopoda.CMR_bodysize_model_bambach_results.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[86,88]

coeffs[6,2] <- coeffs[5,2] + coeffs[6,2] 
coeffs[6,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))

G.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[261,2]
G.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[263,2]
G.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[261,263]

coeffs[6,4] <- coeffs[5,4] + coeffs[6,4] 
coeffs[6,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Crinoids

load("Crinoidea.CMR_bodysize_model_bambach_results.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[86,88]

coeffs[8,2] <- coeffs[7,2] + coeffs[8,2] 
coeffs[8,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))

G.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[261,2]
G.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[263,2]
G.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[261,263]

coeffs[8,4] <- coeffs[7,4] + coeffs[8,4] 
coeffs[8,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Echinoids

load("Echinoidea.CMR_bodysize_model_bambach_results.Rdata")
G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[10,4] <- coeffs[9,4] + coeffs[10,4] 
coeffs[10,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Gastropods

load("Gastropoda.CMR_bodysize_model_bambach_results.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[86,88]

coeffs[12,2] <- coeffs[11,2] + coeffs[12,2] 
coeffs[12,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))

G.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[261,2]
G.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[263,2]
G.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[261,263]

coeffs[12,4] <- coeffs[11,4] + coeffs[12,4] 
coeffs[12,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Ostracods

load("Ostracoda.CMR_bodysize_model_bambach_results.Rdata")

P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[86,88]

coeffs[14,2] <- coeffs[13,2] + coeffs[14,2] 
coeffs[14,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))

G.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[261,2]
G.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[263,2]
G.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[261,263]

coeffs[14,4] <- coeffs[13,4] + coeffs[14,4] 
coeffs[14,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Rhynchonellates

load("Rhynchonellata.CMR_bodysize_model_bambach_results.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[86,88]

coeffs[16,2] <- coeffs[15,2] + coeffs[16,2] 
coeffs[16,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))

G.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[261,2]
G.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[263,2]
G.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[261,263]

coeffs[16,4] <- coeffs[15,4] + coeffs[16,4] 
coeffs[16,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Strophomenates

load("Strophomenata.CMR_bodysize_model_bambach_results.Rdata")

G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[18,4] <- coeffs[17,4] + coeffs[18,4] 
coeffs[18,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Trilobites

load("Trilobita.CMR_bodysize_model_bambach_results.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta.vcv[86,88]

coeffs[20,2] <- coeffs[19,2] + coeffs[20,2] 
coeffs[20,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))


write.table(coeffs, "CMR_coeffs_bodysize_bambach_final.txt", sep="\\t")
__________________________________________________________________________________________
###Big 3 MASS EXTINCTIONS###

coeffs <- read.delim("CMR_coeffs_bodysize_big3.txt")

#Bivalves

load("Bivalvia.CMR_bodysize_model_results_big3.Rdata")
G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[2,4] <- coeffs[1,4] + coeffs[2,4] 
coeffs[2,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Cephalopods

load("Cephalopoda.CMR_bodysize_model_results_big3.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta.vcv[86,88]

coeffs[6,2] <- coeffs[5,2] + coeffs[6,2] 
coeffs[6,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))


#Crinoids

load("Crinoidea.CMR_bodysize_model_results_big3.Rdata")

G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[8,4] <- coeffs[7,4] + coeffs[8,4] 
coeffs[8,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Gastropods

load("Gastropoda.CMR_bodysize_model_results_big3.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[86,88]

coeffs[12,2] <- coeffs[11,2] + coeffs[12,2] 
coeffs[12,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))

G.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[261,2]
G.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[263,2]
G.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[261,263]

coeffs[12,4] <- coeffs[11,4] + coeffs[12,4] 
coeffs[12,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Ostracods

load("Ostracoda.CMR_bodysize_model_results_big3.Rdata")

G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[14,4] <- coeffs[13,4] + coeffs[14,4] 
coeffs[14,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Rhynchonellates

load("Rhynchonellata.CMR_bodysize_model_results_big3.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[86,88]

coeffs[16,2] <- coeffs[15,2] + coeffs[16,2] 
coeffs[16,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))

G.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[261,2]
G.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[263,2]
G.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[261,263]

coeffs[16,4] <- coeffs[15,4] + coeffs[16,4] 
coeffs[16,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Strophomenates

load("Strophomenata.CMR_bodysize_model_results_big3.Rdata")

G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[18,4] <- coeffs[17,4] + coeffs[18,4] 
coeffs[18,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


write.table(coeffs, "CMR_coeffs_bodysize_big3_final.txt", sep="\\t")

__________________________________________________________________________________________

###FIRST 4 MASS EXTINCTIONS###

coeffs <- read.delim("CMR_coeffs_bodysize_first4.txt")

#Bivalves

load("Bivalvia.CMR_bodysize_model_bambach_results.Rdata")
G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[2,4] <- coeffs[1,4] + coeffs[2,4] 
coeffs[2,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Cephalopods

load("Cephalopoda.CMR_bodysize_model_bambach_results.Rdata")

P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[86,88]

coeffs[6,2] <- coeffs[5,2] + coeffs[6,2] 
coeffs[6,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))

G.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[261,2]
G.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta[263,2]
G.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta.vcv[261,263]

coeffs[6,4] <- coeffs[5,4] + coeffs[6,4] 
coeffs[6,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Crinoids

load("Crinoidea.CMR_bodysize_model_bambach_results.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta.vcv[86,88]

coeffs[8,2] <- coeffs[7,2] + coeffs[8,2] 
coeffs[8,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))


#Gastropods

load("Gastropoda.CMR_bodysize_model_bambach_results.Rdata")
G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[12,4] <- coeffs[11,4] + coeffs[12,4] 
coeffs[12,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Ostracods

load("Ostracoda.CMR_bodysize_model_bambach_results.Rdata")

G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[14,4] <- coeffs[13,4] + coeffs[14,4] 
coeffs[14,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Rhynchonellates

load("Rhynchonellata.CMR_bodysize_model_bambach_results.Rdata")

G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[16,4] <- coeffs[15,4] + coeffs[16,4] 
coeffs[16,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Strophomenates

load("Strophomenata.CMR_bodysize_model_bambach_results.Rdata")

G.B1_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[259,2]
G.B2_se <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta[261,2]
G.covB1_2 <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta.vcv[259,261]

coeffs[18,4] <- coeffs[17,4] + coeffs[18,4] 
coeffs[18,5] <- sqrt((G.B1_se^2 + G.B2_se^2) + (2*G.covB1_2))


#Trilobites

load("Trilobita.CMR_bodysize_model_bambach_results.Rdata")
P.B1_se <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta[86,2]
P.B2_se <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta[88,2]
P.covB1_2 <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta.vcv[86,88]

coeffs[20,2] <- coeffs[19,2] + coeffs[20,2] 
coeffs[20,3] <- sqrt((P.B1_se^2 + P.B2_se^2) + (2*P.covB1_2))


write.table(coeffs, "CMR_coeffs_bodysize_first4_final.txt", sep="\\t")
