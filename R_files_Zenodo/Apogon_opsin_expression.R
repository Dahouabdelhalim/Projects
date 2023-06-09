setwd("C:/Users/uqmluehr/Desktop/Paper_Opsins Diversity/03 Quantitative Opsin Expression")

# Load raw data (mapped reads)
opsin.reads <- read.csv("C:/Users/uqmluehr/Desktop/Paper_Opsins Diversity/03 Quantitative Opsin Expression/cardinalfish_retinal_transcriptome_raw_reads.csv", stringsAsFactors=FALSE) # Opsin expression data
# create total pools
opsin.reads$total.opsin <- with(opsin.reads, reads_sws2b + reads_sws2aa + reads_sws2ab + reads_rh2b1 + reads_rh2b2 + reads_rh2a1 + reads_rh2a2 + reads_lws + reads_rh1)
opsin.reads$total.cone <- with(opsin.reads, reads_sws2b + reads_sws2aa + reads_sws2ab + reads_rh2b1 + reads_rh2b2 + reads_rh2a1 + reads_rh2a2 + reads_lws)
opsin.reads$total.rod <- with(opsin.reads, reads_rh1)
opsin.reads$frac.cone <- with(opsin.reads, total.cone / total.opsin)
opsin.reads$frac.rod <- with(opsin.reads, total.rod / total.opsin)
# normalise reads by opsin gene length
opsin.reads$norm.sws2b <- with(opsin.reads, reads_sws2b / len_sws2b)
opsin.reads$norm.sws2aa <- with(opsin.reads, reads_sws2aa / len_sws2aa)
opsin.reads$norm.sws2ab <- with(opsin.reads, reads_sws2ab / len_sws2ab)
opsin.reads$norm.rh2b1 <- with(opsin.reads, reads_rh2b1 / len_rh2b1)
opsin.reads$norm.rh2b2 <- with(opsin.reads, reads_rh2b2 / len_rh2b2)
opsin.reads$norm.rh2a1 <- with(opsin.reads, reads_rh2a1 / len_rh2a1)
opsin.reads$norm.rh2a2 <- with(opsin.reads, reads_rh2a2 / len_rh2a2)
opsin.reads$norm.lws <- with(opsin.reads, reads_lws / len_lws)
opsin.reads$norm.rh1 <- with(opsin.reads, reads_rh1 / len_rh1)
# generate normalized total pools
opsin.reads$norm.total <- with(opsin.reads, norm.sws2b + norm.sws2aa + norm.sws2ab + norm.rh2b1 + norm.rh2b2 + norm.rh2a1 + norm.rh2a2 + norm.lws +norm.rh1)
opsin.reads$norm.totalcone <- with(opsin.reads, norm.sws2b + norm.sws2aa + norm.sws2ab + norm.rh2b1 + norm.rh2b2 + norm.rh2a1 + norm.rh2a2 + norm.lws)
opsin.reads$norm.single <- with(opsin.reads, norm.sws2b + norm.sws2aa + norm.sws2ab)
opsin.reads$norm.double <- with(opsin.reads, norm.rh2b1 + norm.rh2b2 + norm.rh2a1 + norm.rh2a2 + norm.lws)
# calculate fractions (%) of total single, total double and total opsin
opsin.reads$sws2b <- with(opsin.reads, (norm.sws2b / norm.single)*100)
opsin.reads$sws2aa <- with(opsin.reads, (norm.sws2aa / norm.single)*100)
opsin.reads$sws2ab <- with(opsin.reads, (norm.sws2ab / norm.single)*100)
opsin.reads$rh2b1 <- with(opsin.reads, (norm.rh2b1 / norm.double)*100)
opsin.reads$rh2b2 <- with(opsin.reads, (norm.rh2b2 / norm.double)*100)
opsin.reads$rh2a1 <- with(opsin.reads, (norm.rh2a1 / norm.double)*100)
opsin.reads$rh2a2 <- with(opsin.reads, (norm.rh2a2 / norm.double)*100)
opsin.reads$lws <- with(opsin.reads, (norm.lws / norm.double)*100)
opsin.reads$rh1 <- with(opsin.reads, (norm.rh1 / norm.total)*100)
opsin.reads$cone <- with(opsin.reads, (norm.totalcone / norm.total)*100)
opsin.reads$rh2b <- with(opsin.reads, rh2b1 + rh2b2)
opsin.reads$rh2a <- with(opsin.reads, rh2a1 + rh2a2)
################################################################################
write.csv(opsin.reads, "C:/Users/uqmluehr/Desktop/Paper_Opsins Diversity/03 Quantitative Opsin Expression/opsinreads.csv")

# aggregate by TRIBE using MEDIAN
mediansws2b <- aggregate(sws2b ~ tribe, opsin.reads, median, na.rm=TRUE) 
mediansws2aa <- aggregate(sws2aa ~ tribe, opsin.reads, median, na.rm=TRUE)
mediansws2ab <- aggregate(sws2ab ~ tribe, opsin.reads, median, na.rm=TRUE)
medianrh2b <- aggregate(rh2b ~ tribe, opsin.reads, median, na.rm=TRUE)
medianrh2b1 <- aggregate(rh2b1 ~ tribe, opsin.reads, median, na.rm=TRUE)
medianrh2b2 <- aggregate(rh2b2 ~ tribe, opsin.reads, median, na.rm=TRUE)
medianrh2a <- aggregate(rh2a ~ tribe, opsin.reads, median, na.rm=TRUE)
medianrh2a1 <- aggregate(rh2a1 ~ tribe, opsin.reads, median, na.rm=TRUE)
medianrh2a2 <- aggregate(rh2a2 ~ tribe, opsin.reads, median, na.rm=TRUE)
medianlws <- aggregate(lws ~ tribe, opsin.reads, median, na.rm=TRUE)
mediancone <- aggregate(cone ~ tribe, opsin.reads, median, na.rm=TRUE)
medianrh1 <- aggregate(rh1 ~ tribe, opsin.reads, median, na.rm=TRUE)

median <- mediansws2b
median$sws2aa <- mediansws2aa$sws2aa
median$sws2ab <- mediansws2ab$sws2ab
median$rh2b <- medianrh2b$rh2b 
median$rh2a <- medianrh2a$rh2a
median$lws <- medianlws$lws
median$cone <- mediancone$cone
median$rh1 <- medianrh1$rh1


mediansws2b
mediansws2aa 
mediansws2ab 
medianrh2b 
medianrh2a 
medianlws 
mediancone
medianrh1 


IQRsws2b <- aggregate(sws2b ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRsws2aa <- aggregate(sws2aa ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRsws2ab <- aggregate(sws2ab ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRrh2b <- aggregate(rh2b ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRrh2b1 <- aggregate(rh2b1 ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRrh2b2 <- aggregate(rh2b2 ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRrh2a <- aggregate(rh2a ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRrh2a1 <- aggregate(rh2a1 ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRrh2a2 <- aggregate(rh2a2 ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRlws <- aggregate(lws ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRcone <- aggregate(cone ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)
IQRrh1 <- aggregate(rh1 ~ tribe, opsin.reads, FUN=IQR, na.rm=TRUE)

IQRsws2b
IQRsws2aa
IQRsws2ab
IQRrh2b
IQRrh2a 
IQRlws
IQRcone 
IQRrh1


meansws2b <- aggregate(sws2b ~ tribe, opsin.reads, mean, na.rm=TRUE) 
meansws2aa <- aggregate(sws2aa ~ tribe, opsin.reads, mean, na.rm=TRUE)
meansws2ab <- aggregate(sws2ab ~ tribe, opsin.reads, mean, na.rm=TRUE)
meanrh2b <- aggregate(rh2b ~ tribe, opsin.reads, mean, na.rm=TRUE)
meanrh2b1 <- aggregate(rh2b1 ~ tribe, opsin.reads, mean, na.rm=TRUE)
meanrh2b2 <- aggregate(rh2b2 ~ tribe, opsin.reads, mean, na.rm=TRUE)
meanrh2a <- aggregate(rh2a ~ tribe, opsin.reads, mean, na.rm=TRUE)
meanrh2a1 <- aggregate(rh2a1 ~ tribe, opsin.reads, mean, na.rm=TRUE)
meanrh2a2 <- aggregate(rh2a2 ~ tribe, opsin.reads, mean, na.rm=TRUE)
meanlws <- aggregate(lws ~ tribe, opsin.reads, mean, na.rm=TRUE)
meancone <- aggregate(cone ~ tribe, opsin.reads, mean, na.rm=TRUE)
meanrh1 <- aggregate(rh1 ~ tribe, opsin.reads, mean, na.rm=TRUE)

# average by SPECIES using MEAN
meansws2b <- aggregate(sws2b ~ species, opsin.reads, mean, na.rm=TRUE) 
meansws2aa <- aggregate(sws2aa ~ species, opsin.reads, mean, na.rm=TRUE)
meansws2ab <- aggregate(sws2ab ~ species, opsin.reads, mean, na.rm=TRUE)
meanrh2b <- aggregate(rh2b ~ species, opsin.reads, mean, na.rm=TRUE)
meanrh2b1 <- aggregate(rh2b1 ~ species, opsin.reads, mean, na.rm=TRUE)
meanrh2b2 <- aggregate(rh2b2 ~ species, opsin.reads, mean, na.rm=TRUE)
meanrh2a <- aggregate(rh2a ~ species, opsin.reads, mean, na.rm=TRUE)
meanrh2a1 <- aggregate(rh2a1 ~ species, opsin.reads, mean, na.rm=TRUE)
meanrh2a2 <- aggregate(rh2a2 ~ species, opsin.reads, mean, na.rm=TRUE)
meanlws <- aggregate(lws ~ species, opsin.reads, mean, na.rm=TRUE)
meancone <- aggregate(cone ~ species, opsin.reads, mean, na.rm=TRUE)
meanrh1 <- aggregate(rh1 ~ species, opsin.reads, mean, na.rm=TRUE)
# calculate standard deviation of mean per species
sdsws2b <- aggregate(sws2b ~ species, opsin.reads, sd, na.rm=TRUE)
sdsws2aa <- aggregate(sws2aa ~ species, opsin.reads, sd, na.rm=TRUE)
sdsws2ab <- aggregate(sws2ab ~ species, opsin.reads, sd, na.rm=TRUE)
sdrh2b <- aggregate(rh2b ~ species, opsin.reads, sd, na.rm=TRUE)
sdrh2b1 <- aggregate(rh2b1 ~ species, opsin.reads, sd, na.rm=TRUE)
sdrh2b2 <- aggregate(rh2b2 ~ species, opsin.reads, sd, na.rm=TRUE)
sdrh2a <- aggregate(rh2a ~ species, opsin.reads, sd, na.rm=TRUE)
sdrh2a1 <- aggregate(rh2a1 ~ species, opsin.reads, sd, na.rm=TRUE)
sdrh2a2 <- aggregate(rh2a2 ~ species, opsin.reads, sd, na.rm=TRUE)
sdlws <- aggregate(lws ~ species, opsin.reads, sd, na.rm=TRUE)
sdcone <- aggregate(cone ~ species, opsin.reads, sd, na.rm=TRUE)
sdrh1 <- aggregate(rh1 ~ species, opsin.reads, sd, na.rm=TRUE)

mean <- meansws2b
mean$sdsws2b <- sdsws2b$sws2b
mean$sws2aa <- meansws2aa$sws2aa
mean$sdsws2aa <- sdsws2aa$sws2aa
mean$sws2ab <- meansws2ab$sws2ab
mean$sdsws2ab <- sdsws2ab$sws2ab
mean$rh2b <- meanrh2b$rh2b 
mean$sdrh2b <-sdrh2b$rh2b
mean$rh2b1 <- meanrh2b1$rh2b1
mean$sdrh2b1 <- sdrh2b1$rh2b1
mean$rh2b2 <- meanrh2b2$rh2b2
mean$sdrh2b2 <- sdrh2b2$rh2b2
mean$rh2a <- meanrh2a$rh2a
mean$sdrh2a <- sdrh2a$rh2a
mean$rh2a1 <- meanrh2a1$rh2a1
mean$sdrh2a1 <- sdrh2a1$rh2a1
mean$rh2a2 <- meanrh2a2$rh2a2
mean$sdrh2a2 <- sdrh2a2$rh2a2
mean$lws <- meanlws$lws
mean$sdlws <- sdlws$lws
mean$cone <- meancone$cone
mean$sdcone <- sdcone$cone
mean$rh1 <- meanrh1$rh1
mean$sdrh1 <- sdrh1$rh1



write.csv(mean, "C:/Users/uqmluehr/Desktop/Paper_Opsins Diversity/03 Quantitative Opsin Expression/cardinals_expression.csv")
write.csv(opsin.reads, "C:/Users/uqmluehr/Desktop/Paper_Opsins Diversity/03 Quantitative Opsin Expression/opsinreads.csv")



opsins <- read.csv("C:/Users/uqmluehr/Desktop/Paper_Opsins Diversity/03 Quantitative Opsin Expression/opsinreads.csv", stringsAsFactors=FALSE) # Opsin expression data

opsins$single <- with(opsins, sws2b + sws2aa + sws2ab) 
opsins$double <- with(opsins, rh2b + rh2a + lws)
opsins$normsws2b <- with(opsins, (sws2b / single)*100)
opsins$normsws2aa <- with(opsins, (sws2aa / single)*100)
opsins$normsws2ab <- with(opsins, (sws2ab / single)*100)
opsins$normrh2b <- with(opsins, (rh2b / double)*100)
opsins$normrh2b1 <- with(opsins, (rh2b1 / double)*100)
opsins$normrh2b2 <- with(opsins, (rh2b2 / double)*100)
opsins$normrh2a <- with(opsins, (rh2a / double)*100)
opsins$normrh2a1 <- with(opsins, (rh2a1 / double)*100)
opsins$normrh2a2 <- with(opsins, (rh2a2 / double)*100)
opsins$normlws <- with(opsins, (lws / double)*100)
