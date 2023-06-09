rm(list = ls())

#### PRELIMINARIES ############################################################################################# 

#*Uploads the needed libraries --------------------------------------------------------------------------------

require(ggplot2)

require(data.table)

require(plotly)

require(DT)

require(R2HTML)


#*Set the number of significant digits for the output --------------------------
sig_dig = 4

#*Set the working directory ---------------------------------------------------------------------------------

#get the default wd
default_wd <- getwd()

#Set the directory where all the saved outputs will be stored
setwd("~/Work_directory")

new_wd <- getwd()

#### DATA UPLOAD ###############################################################################################################
# - Uploads the .csv files containing the information related to the F2N, +ve datasets.
# - Separates the metadata from the main data
# - Picks only the "LL" and "NN" sets

LIV_F2N_p_ve_metadata <- read.csv("./Data/+ve, LIV, F2N.csv", stringsAsFactors=F)
LIV_F2N_p_ve <-read.csv("./Data/+ve, LIV, F2N.csv", stringsAsFactors=F , skip=10 ) 
LIV_F2N_p_ve <- LIV_F2N_p_ve[!is.na(LIV_F2N_p_ve$m.z),]
LIV_F2N_p_ve_LL <- LIV_F2N_p_ve[,grep(as.vector(LIV_F2N_p_ve_metadata[4, ]),pattern="LL")]
LIV_F2N_p_ve_NN <- LIV_F2N_p_ve[,grep(as.vector(LIV_F2N_p_ve_metadata[4, ]),pattern="NN")]
rownames(LIV_F2N_p_ve_LL) <- LIV_F2N_p_ve$Lipid.variable
rownames(LIV_F2N_p_ve_NN) <- LIV_F2N_p_ve$Lipid.variable

SER_F2N_p_ve_metadata <- read.csv("./Data/+ve, SER, F2N.csv", stringsAsFactors=F)
SER_F2N_p_ve <-read.csv("./Data/+ve, SER, F2N.csv", stringsAsFactors=F , skip=10 ) 
SER_F2N_p_ve <- SER_F2N_p_ve[!is.na(SER_F2N_p_ve$m.z),]
SER_F2N_p_ve_LL <- SER_F2N_p_ve[,grep(as.vector(SER_F2N_p_ve_metadata[4, ]),pattern="LL")]
SER_F2N_p_ve_NN <- SER_F2N_p_ve[,grep(as.vector(SER_F2N_p_ve_metadata[4, ]),pattern="NN")]
rownames(SER_F2N_p_ve_LL) <- SER_F2N_p_ve$Lipid.variable
rownames(SER_F2N_p_ve_NN) <- SER_F2N_p_ve$Lipid.variable

HEA_F2N_p_ve_metadata <- read.csv("./Data/+ve, HEA, F2N.csv", stringsAsFactors=F)
HEA_F2N_p_ve <-read.csv("./Data/+ve, HEA, F2N.csv", stringsAsFactors=F , skip=10 ) 
HEA_F2N_p_ve <- HEA_F2N_p_ve[!is.na(HEA_F2N_p_ve$m.z),]
HEA_F2N_p_ve_LL <- HEA_F2N_p_ve[,grep(as.vector(HEA_F2N_p_ve_metadata[4, ]),pattern="LL")]
HEA_F2N_p_ve_NN <- HEA_F2N_p_ve[,grep(as.vector(HEA_F2N_p_ve_metadata[4, ]),pattern="NN")]
rownames(HEA_F2N_p_ve_LL) <- HEA_F2N_p_ve$Lipid.variable
rownames(HEA_F2N_p_ve_NN) <- HEA_F2N_p_ve$Lipid.variable

CEB_F2N_p_ve_metadata <- read.csv("./Data/+ve, CEB, F2N.csv", stringsAsFactors=F)
CEB_F2N_p_ve <-read.csv("./Data/+ve, CEB, F2N.csv", stringsAsFactors=F , skip=10 ) 
CEB_F2N_p_ve <- CEB_F2N_p_ve[!is.na(CEB_F2N_p_ve$m.z),]
CEB_F2N_p_ve_LL <- CEB_F2N_p_ve[,grep(as.vector(CEB_F2N_p_ve_metadata[4, ]),pattern="LL")]
CEB_F2N_p_ve_NN <- CEB_F2N_p_ve[,grep(as.vector(CEB_F2N_p_ve_metadata[4, ]),pattern="NN")]
rownames(CEB_F2N_p_ve_LL) <- CEB_F2N_p_ve$Lipid.variable
rownames(CEB_F2N_p_ve_NN) <- CEB_F2N_p_ve$Lipid.variable

RiB_F2N_p_ve_metadata <- read.csv("./Data/+ve, RiB, F2N.csv", stringsAsFactors=F)
RiB_F2N_p_ve <-read.csv("./Data/+ve, RiB, F2N.csv", stringsAsFactors=F , skip=10 ) 
RiB_F2N_p_ve <- RiB_F2N_p_ve[!is.na(RiB_F2N_p_ve$m.z),]
RiB_F2N_p_ve_LL <- RiB_F2N_p_ve[,grep(as.vector(RiB_F2N_p_ve_metadata[4, ]),pattern="LL")]
RiB_F2N_p_ve_NN <- RiB_F2N_p_ve[,grep(as.vector(RiB_F2N_p_ve_metadata[4, ]),pattern="NN")]
rownames(RiB_F2N_p_ve_LL) <- RiB_F2N_p_ve$Lipid.variable
rownames(RiB_F2N_p_ve_NN) <- RiB_F2N_p_ve$Lipid.variable

##### DATA ANALYSIS ###########################################################################################################

# Sets a function that counts the number of zeros in a vector -----------------------------------------------------------
count_zeroes <- function(x){length(which(x==0))}

# Counts the zeroes in each row of the considered datasets -----------------------------------------------------------
LIV_F2N_p_ve_LL_zeroes <- apply(LIV_F2N_p_ve_LL[, c(3:ncol(LIV_F2N_p_ve_LL))], 1, count_zeroes)
LIV_F2N_p_ve_NN_zeroes <- apply(LIV_F2N_p_ve_NN[, c(3:ncol(LIV_F2N_p_ve_NN))], 1, count_zeroes)

SER_F2N_p_ve_LL_zeroes <- apply(SER_F2N_p_ve_LL[, c(3:ncol(SER_F2N_p_ve_LL))], 1, count_zeroes)
SER_F2N_p_ve_NN_zeroes <- apply(SER_F2N_p_ve_NN[, c(3:ncol(SER_F2N_p_ve_NN))], 1, count_zeroes)

HEA_F2N_p_ve_LL_zeroes <- apply(HEA_F2N_p_ve_LL[, c(3:ncol(HEA_F2N_p_ve_LL))], 1, count_zeroes)
HEA_F2N_p_ve_NN_zeroes <- apply(HEA_F2N_p_ve_NN[, c(3:ncol(HEA_F2N_p_ve_NN))], 1, count_zeroes)

CEB_F2N_p_ve_LL_zeroes <- apply(CEB_F2N_p_ve_LL[, c(3:ncol(CEB_F2N_p_ve_LL))], 1, count_zeroes)
CEB_F2N_p_ve_NN_zeroes <- apply(CEB_F2N_p_ve_NN[, c(3:ncol(CEB_F2N_p_ve_NN))], 1, count_zeroes)

RiB_F2N_p_ve_LL_zeroes <- apply(RiB_F2N_p_ve_LL[, c(3:ncol(RiB_F2N_p_ve_LL))], 1, count_zeroes)
RiB_F2N_p_ve_NN_zeroes <- apply(RiB_F2N_p_ve_NN[, c(3:ncol(RiB_F2N_p_ve_NN))], 1, count_zeroes)


# Eliminates all the rows tha have a number of zeroes greater than the 20% of the number of total entries in that row -----------------------------------------------------------
SER_F2N_p_ve_LL_nozeroes <- SER_F2N_p_ve_LL[-which(SER_F2N_p_ve_LL_zeroes > ncol(SER_F2N_p_ve_LL)*20/100) ,]
SER_F2N_p_ve_NN_nozeroes <- SER_F2N_p_ve_LL[-which(SER_F2N_p_ve_NN_zeroes > ncol(SER_F2N_p_ve_NN)*20/100) ,]

LIV_F2N_p_ve_LL_nozeroes <- LIV_F2N_p_ve_LL[-which(LIV_F2N_p_ve_LL_zeroes > ncol(LIV_F2N_p_ve_LL)*20/100) ,]
LIV_F2N_p_ve_NN_nozeroes <- LIV_F2N_p_ve_LL[-which(LIV_F2N_p_ve_NN_zeroes > ncol(LIV_F2N_p_ve_NN)*20/100) ,]

HEA_F2N_p_ve_LL_nozeroes <- HEA_F2N_p_ve_LL[-which(HEA_F2N_p_ve_LL_zeroes > ncol(HEA_F2N_p_ve_LL)*20/100) ,]
HEA_F2N_p_ve_NN_nozeroes <- HEA_F2N_p_ve_LL[-which(HEA_F2N_p_ve_NN_zeroes > ncol(HEA_F2N_p_ve_NN)*20/100) ,]

CEB_F2N_p_ve_LL_nozeroes <- CEB_F2N_p_ve_LL[-which(CEB_F2N_p_ve_LL_zeroes > ncol(CEB_F2N_p_ve_LL)*20/100) ,]
CEB_F2N_p_ve_NN_nozeroes <- CEB_F2N_p_ve_LL[-which(CEB_F2N_p_ve_NN_zeroes > ncol(CEB_F2N_p_ve_NN)*20/100) ,]

RiB_F2N_p_ve_LL_nozeroes <- RiB_F2N_p_ve_LL[-which(RiB_F2N_p_ve_LL_zeroes > ncol(RiB_F2N_p_ve_LL)*20/100) ,]
RiB_F2N_p_ve_NN_nozeroes <- RiB_F2N_p_ve_LL[-which(RiB_F2N_p_ve_NN_zeroes > ncol(RiB_F2N_p_ve_NN)*20/100) ,]

# Finds B-type_lipids ---------------------------------------------------------------------------------------


SER_LIV_F2N_p_ve_NN_common <- intersect(rownames(SER_F2N_p_ve_NN_nozeroes), rownames(LIV_F2N_p_ve_NN_nozeroes))
SER_LIV_F2N_p_ve_LL_common <- intersect(rownames(SER_F2N_p_ve_LL_nozeroes), rownames(LIV_F2N_p_ve_LL_nozeroes))
write.csv(SER_LIV_F2N_p_ve_NN_common,"./Results_F2N_p_ve/SER_LIV_F2N_p_ve_NN_common.csv")
write.csv(SER_LIV_F2N_p_ve_LL_common,"./Results_F2N_p_ve/SER_LIV_F2N_p_ve_LL_common.csv")

SER_HEA_F2N_p_ve_NN_common <- intersect(rownames(SER_F2N_p_ve_NN_nozeroes), rownames(HEA_F2N_p_ve_NN_nozeroes))
SER_HEA_F2N_p_ve_LL_common <- intersect(rownames(SER_F2N_p_ve_LL_nozeroes), rownames(HEA_F2N_p_ve_LL_nozeroes))
write.csv(SER_HEA_F2N_p_ve_NN_common,"./Results_F2N_p_ve/SER_HEA_F2N_p_ve_NN_common.csv")
write.csv(SER_HEA_F2N_p_ve_LL_common,"./Results_F2N_p_ve/SER_HEA_F2N_p_ve_LL_common.csv")

SER_CEB_F2N_p_ve_NN_common <- intersect(rownames(SER_F2N_p_ve_NN_nozeroes), rownames(CEB_F2N_p_ve_NN_nozeroes))
SER_CEB_F2N_p_ve_LL_common <- intersect(rownames(SER_F2N_p_ve_LL_nozeroes), rownames(CEB_F2N_p_ve_LL_nozeroes))
write.csv(SER_CEB_F2N_p_ve_NN_common,"./Results_F2N_p_ve/SER_CEB_F2N_p_ve_NN_common.csv")
write.csv(SER_CEB_F2N_p_ve_LL_common,"./Results_F2N_p_ve/SER_CEB_F2N_p_ve_LL_common.csv")

SER_RiB_F2N_p_ve_NN_common <- intersect(rownames(SER_F2N_p_ve_NN_nozeroes), rownames(RiB_F2N_p_ve_NN_nozeroes))
SER_RiB_F2N_p_ve_LL_common <- intersect(rownames(SER_F2N_p_ve_LL_nozeroes), rownames(RiB_F2N_p_ve_LL_nozeroes))
write.csv(SER_RiB_F2N_p_ve_NN_common,"./Results_F2N_p_ve/SER_RiB_F2N_p_ve_NN_common.csv")
write.csv(SER_RiB_F2N_p_ve_LL_common,"./Results_F2N_p_ve/SER_RiB_F2N_p_ve_LL_common.csv")


# Finds A-type_lipids ---------------------------------------------------------------------------------------

A_set_F2N_NN <- Reduce(intersect, list(rownames(SER_F2N_p_ve_NN_nozeroes),
                                        #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                        rownames(HEA_F2N_p_ve_NN_nozeroes),
                                        rownames(CEB_F2N_p_ve_NN_nozeroes),
                                        rownames(LIV_F2N_p_ve_NN_nozeroes),
                                        rownames(RiB_F2N_p_ve_NN_nozeroes)) 
)
write.csv(A_set_F2N_NN,"./Results_F2N_p_ve/A_set_F2N_NN.csv")


A_set_F2N_LL <- Reduce(intersect, list(rownames(SER_F2N_p_ve_LL_nozeroes),
                                       #rownames(ADI_F2N_p_ve_LL_nozeroes),
                                       rownames(HEA_F2N_p_ve_LL_nozeroes),
                                       rownames(CEB_F2N_p_ve_LL_nozeroes),
                                       rownames(LIV_F2N_p_ve_NN_nozeroes),
                                       rownames(RiB_F2N_p_ve_LL_nozeroes)) 
)
write.csv(A_set_F2N_LL,"./Results_F2N_p_ve/A_set_F2N_LL.csv")



# Finds U-type_lipids ---------------------------------------------------------------------------------------


union_not_LIV_F2N_p_ve_NN <- Reduce(union, list(rownames(SER_F2N_p_ve_NN_nozeroes),
                                                rownames(HEA_F2N_p_ve_NN_nozeroes),
                                                rownames(CEB_F2N_p_ve_NN_nozeroes),
                                                #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                                rownames(RiB_F2N_p_ve_NN_nozeroes))
)
LIV_F2N_p_ve_NN_unique <- setdiff(rownames(LIV_F2N_p_ve_NN_nozeroes),union_not_LIV_F2N_p_ve_NN)
write.csv(LIV_F2N_p_ve_NN_unique,"./Results_F2N_p_ve/LIV_F2N_p_ve_NN_unique.csv")

#
#
union_not_LIV_F2N_p_ve_LL <- Reduce(union, list(rownames(SER_F2N_p_ve_LL_nozeroes),
                                                rownames(HEA_F2N_p_ve_LL_nozeroes),
                                                rownames(CEB_F2N_p_ve_LL_nozeroes),
                                                #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                                rownames(RiB_F2N_p_ve_LL_nozeroes))
)
LIV_F2N_p_ve_LL_unique <- setdiff(rownames(LIV_F2N_p_ve_LL_nozeroes),union_not_LIV_F2N_p_ve_LL)
write.csv(LIV_F2N_p_ve_LL_unique,"./Results_F2N_p_ve/LIV_F2N_p_ve_LL_unique.csv")



union_not_HEA_F2N_p_ve_NN <- Reduce(union, list(rownames(SER_F2N_p_ve_NN_nozeroes),
                                                rownames(LIV_F2N_p_ve_NN_nozeroes),
                                                rownames(CEB_F2N_p_ve_NN_nozeroes),
                                                #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                                rownames(RiB_F2N_p_ve_NN_nozeroes))
)
HEA_F2N_p_ve_NN_unique <- setdiff(rownames(HEA_F2N_p_ve_NN_nozeroes),union_not_HEA_F2N_p_ve_NN)
write.csv(HEA_F2N_p_ve_NN_unique,"./Results_F2N_p_ve/HEA_F2N_p_ve_NN_unique.csv")

#
#
union_not_HEA_F2N_p_ve_LL <- Reduce(union, list(rownames(SER_F2N_p_ve_LL_nozeroes),
                                                rownames(LIV_F2N_p_ve_LL_nozeroes),
                                                rownames(CEB_F2N_p_ve_LL_nozeroes),
                                                #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                                rownames(RiB_F2N_p_ve_LL_nozeroes))
)
HEA_F2N_p_ve_LL_unique <- setdiff(rownames(HEA_F2N_p_ve_LL_nozeroes),union_not_HEA_F2N_p_ve_LL)
write.csv(HEA_F2N_p_ve_LL_unique,"./Results_F2N_p_ve/HEA_F2N_p_ve_LL_unique.csv")



union_not_CEB_F2N_p_ve_NN <- Reduce(union, list(rownames(SER_F2N_p_ve_NN_nozeroes),
                                                rownames(HEA_F2N_p_ve_NN_nozeroes),
                                                rownames(LIV_F2N_p_ve_NN_nozeroes),
                                                #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                                rownames(RiB_F2N_p_ve_NN_nozeroes))
)
CEB_F2N_p_ve_NN_unique <- setdiff(rownames(CEB_F2N_p_ve_NN_nozeroes),union_not_CEB_F2N_p_ve_NN)
write.csv(CEB_F2N_p_ve_NN_unique,"./Results_F2N_p_ve/CEB_F2N_p_ve_NN_unique.csv")
#
#
union_not_CEB_F2N_p_ve_LL <- Reduce(union, list(rownames(SER_F2N_p_ve_LL_nozeroes),
                                                rownames(HEA_F2N_p_ve_LL_nozeroes),
                                                rownames(LIV_F2N_p_ve_LL_nozeroes),
                                                #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                                rownames(RiB_F2N_p_ve_LL_nozeroes))
)
CEB_F2N_p_ve_LL_unique <- setdiff(rownames(CEB_F2N_p_ve_LL_nozeroes),union_not_CEB_F2N_p_ve_LL)
write.csv(CEB_F2N_p_ve_LL_unique,"./Results_F2N_p_ve/CEB_F2N_p_ve_LL_unique.csv")



union_not_RiB_F2N_p_ve_NN <- Reduce(union, list(rownames(SER_F2N_p_ve_NN_nozeroes),
                                                rownames(HEA_F2N_p_ve_NN_nozeroes),
                                                rownames(LIV_F2N_p_ve_NN_nozeroes),
                                                #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                                rownames(CEB_F2N_p_ve_NN_nozeroes))
)
RiB_F2N_p_ve_NN_unique <- setdiff(rownames(RiB_F2N_p_ve_NN_nozeroes),union_not_RiB_F2N_p_ve_NN)
write.csv(RiB_F2N_p_ve_NN_unique,"./Results_F2N_p_ve/RiB_F2N_p_ve_NN_unique.csv")
#
#
union_not_RiB_F2N_p_ve_LL <- Reduce(union, list(rownames(SER_F2N_p_ve_LL_nozeroes),
                                                rownames(HEA_F2N_p_ve_LL_nozeroes),
                                                rownames(LIV_F2N_p_ve_LL_nozeroes),
                                                #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                                rownames(CEB_F2N_p_ve_LL_nozeroes))
)
RiB_F2N_p_ve_LL_unique <- setdiff(rownames(RiB_F2N_p_ve_LL_nozeroes),union_not_RiB_F2N_p_ve_LL)
write.csv(RiB_F2N_p_ve_LL_unique,"./Results_F2N_p_ve/RiB_F2N_p_ve_LL_unique.csv")


union_not_SER_F2N_p_ve_NN <- Reduce(union, list(rownames(RiB_F2N_p_ve_NN_nozeroes),
                                                rownames(HEA_F2N_p_ve_NN_nozeroes),
                                                rownames(LIV_F2N_p_ve_NN_nozeroes),
                                                #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                                rownames(CEB_F2N_p_ve_NN_nozeroes))
)
SER_F2N_p_ve_NN_unique <- setdiff(rownames(SER_F2N_p_ve_NN_nozeroes),union_not_SER_F2N_p_ve_NN)
write.csv(SER_F2N_p_ve_NN_unique,"./Results_F2N_p_ve/SER_F2N_p_ve_NN_unique.csv")
#
#
union_not_SER_F2N_p_ve_LL <- Reduce(union, list(rownames(RiB_F2N_p_ve_LL_nozeroes),
                                                rownames(HEA_F2N_p_ve_LL_nozeroes),
                                                rownames(LIV_F2N_p_ve_LL_nozeroes),
                                                #rownames(ADI_F2N_p_ve_NN_nozeroes),
                                                rownames(CEB_F2N_p_ve_LL_nozeroes))
)
SER_F2N_p_ve_LL_unique <- setdiff(rownames(SER_F2N_p_ve_LL_nozeroes),union_not_SER_F2N_p_ve_LL)
write.csv(SER_F2N_p_ve_LL_unique,"./Results_F2N_p_ve/SER_F2N_p_ve_LL_unique.csv")


########## COMPUTES t-Tests ##################
# Computes averages and t tests for each subset

pval<- vector(length=nrow(LIV_F2N_p_ve))
names<- vector(length=nrow(LIV_F2N_p_ve))
meandif <- vector(length=nrow(LIV_F2N_p_ve))

for(i in 1:nrow(LIV_F2N_p_ve)){
  
  names[i] <- rownames(LIV_F2N_p_ve_LL)[i]
  pval[i]<-t.test(LIV_F2N_p_ve_LL[i,], LIV_F2N_p_ve_NN[i,])$p.value
  meandif[i]<- mean(as.numeric(LIV_F2N_p_ve_LL[i,]), na.rm=T) - mean(as.numeric(LIV_F2N_p_ve_NN[i,]), na.rm=T) 
}

LIV_F2N_p_ve_meandiffs <- data.frame(pval_LIV_F2N_p_ve=pval, meandif_LIV_F2N_p_ve=meandif, row.names = names) 


pval<- vector(length=nrow(SER_F2N_p_ve))
names<- vector(length=nrow(SER_F2N_p_ve))
meandif <- vector(length=nrow(SER_F2N_p_ve))

for(i in 1:nrow(SER_F2N_p_ve)){
  
  names[i] <- rownames(SER_F2N_p_ve_LL)[i]
  pval[i]<-t.test(SER_F2N_p_ve_LL[i,], SER_F2N_p_ve_NN[i,])$p.value
  meandif[i]<- mean(as.numeric(SER_F2N_p_ve_LL[i,]), na.rm=T) - mean(as.numeric(SER_F2N_p_ve_NN[i,]), na.rm=T) 
}

SER_F2N_p_ve_meandiffs <- data.frame(pval_SER_F2N_p_ve=pval, meandif_SER_F2N_p_ve=meandif, row.names = names) 


pval<- vector(length=nrow(HEA_F2N_p_ve))
names<- vector(length=nrow(HEA_F2N_p_ve))
meandif <- vector(length=nrow(HEA_F2N_p_ve))

for(i in 1:nrow(HEA_F2N_p_ve)){
  
  names[i] <- rownames(HEA_F2N_p_ve_LL)[i]
  pval[i]<-t.test(HEA_F2N_p_ve_LL[i,], HEA_F2N_p_ve_NN[i,])$p.value
  meandif[i]<- mean(as.numeric(HEA_F2N_p_ve_LL[i,]), na.rm=T) - mean(as.numeric(HEA_F2N_p_ve_NN[i,]), na.rm=T) 
}

HEA_F2N_p_ve_meandiffs <- data.frame(pval_HEA_F2N_p_ve=pval, meandif_HEA_F2N_p_ve=meandif, row.names = names) 


pval<- vector(length=nrow(CEB_F2N_p_ve))
names<- vector(length=nrow(CEB_F2N_p_ve))
meandif <- vector(length=nrow(CEB_F2N_p_ve))

for(i in 1:nrow(CEB_F2N_p_ve)){
  
  names[i] <- rownames(CEB_F2N_p_ve_LL)[i]
  pval[i]<-t.test(CEB_F2N_p_ve_LL[i,], CEB_F2N_p_ve_NN[i,])$p.value
  meandif[i]<- mean(as.numeric(CEB_F2N_p_ve_LL[i,]), na.rm=T) - mean(as.numeric(CEB_F2N_p_ve_NN[i,]), na.rm=T) 
}

CEB_F2N_p_ve_meandiffs <- data.frame(pval_CEB_F2N_p_ve=pval, meandif_CEB_F2N_p_ve=meandif, row.names = names) 

pval<- vector(length=nrow(RiB_F2N_p_ve))
names<- vector(length=nrow(RiB_F2N_p_ve))
meandif <- vector(length=nrow(RiB_F2N_p_ve))

for(i in 1:nrow(RiB_F2N_p_ve)){
  
  names[i] <- rownames(RiB_F2N_p_ve_LL)[i]
  pval[i]<-t.test(RiB_F2N_p_ve_LL[i,], RiB_F2N_p_ve_NN[i,])$p.value
  meandif[i]<- mean(as.numeric(RiB_F2N_p_ve_LL[i,]), na.rm=T) - mean(as.numeric(RiB_F2N_p_ve_NN[i,]), na.rm=T) 
}

RiB_F2N_p_ve_meandiffs <- data.frame(pval_RiB_F2N_p_ve=pval, meandif_RiB_F2N_p_ve=meandif, row.names = names) 



mbind<-function(...){Reduce( function(x,y){cbind(x,y[match(row.names(x),row.names(y)),])}, list(...) )}

F2N_p_ve_meandiffs<-mbind(#ADI_F2N_p_ve_meandiffs,
  LIV_F2N_p_ve_meandiffs,
  SER_F2N_p_ve_meandiffs,
  HEA_F2N_p_ve_meandiffs,
  #BRA_F2N_p_ve_meandiffs,
  CEB_F2N_p_ve_meandiffs,
  RiB_F2N_p_ve_meandiffs
)

write.csv(F2N_p_ve_meandiffs, file = "F2N_p_ve_meandiffs.csv")





