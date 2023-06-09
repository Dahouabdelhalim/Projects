###PLOTS FOR SUPPLEMENTARY FIGURES 1â€“4###

totalOccs <- read.csv("pbdb.cleaned.classes.no.body.size.csv")  
sizeOccs <- read.csv("pbdb.cleaned.classes.csv")  

x <- table(unlist(sizeOccs$class), sizeOccs$max_ma) #table genus occurrences by interval
x <- x[, match(colnames(x), rev(colnames(x)))]

w <- rev(table(sizeOccs$max_ma)) #table genus occurrences by interval
x <- rbind(w, x)
rownames(x)[1] <- "All Genera"

y <- table(unlist(totalOccs$class), totalOccs$max_ma) #table genus occurrences by interval
y <- y[, match(colnames(y), rev(colnames(y)))] 

z <- rev(table(totalOccs$max_ma)) #table genus occurrences by interval
y <- rbind(z, y)
rownames(y)[1] <- "All Genera"

coverage <- x/y

barplot(coverage[1,], ylim=c(0,100), las=1, main="All Genera PBDB Occurrence Body Size Coverage", ylab="% Occurrence Coverage", xlab="Geologic Age (Ma)", col="gray95")

par(mfrow=c(5,2))

barplot(coverage[2,], ylim=c(0,100),las=1, main="Bivalvia", ylab="% Occurrence Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(coverage[3,], ylim=c(0,100), las=1, main="Bony fish", ylab="% Occurrence Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(coverage[4,], ylim=c(0,100), las=1, main="Cephalopoda", ylab="% Occurrence Coverage",xlab="Geologic Age (Ma)", col="gray95")

barplot(coverage[5,], ylim=c(0,100), las=1, main="Crinoidea", ylab="% Occurrence Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(coverage[6,], ylim=c(0,100), las=1, main="Echinoidea", ylab="% Occurrence Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(coverage[7,], ylim=c(0,100), las=1, main="Gastropoda", ylab="% Occurrence Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(coverage[8,], ylim=c(0,100), las=1, main="Ostracoda", ylab="% Occurrence Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(coverage[9,], ylim=c(0,100), las=1, main="Rhynchonellata", ylab="% Occurrence Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(coverage[10,], ylim=c(0,100),las=1, main="Strophomenata", ylab="% Occurrence Coverage", xlab="Geologic Age (Ma)", col="gray95" )

barplot(coverage[11,], ylim=c(0,100), las=1,main="Trilobita", ylab="% Occurrence Coverage", xlab="Geologic Age (Ma)", col="gray95")

dev.new()
par(mfrow=c(5,2))

barplot(x[2,],las=1, ylim=c(0,250), main="Bivalvia", ylab="# of Occurrences", xlab="Geologic Age (Ma)", col="gray95")

barplot(x[3,], ylim=c(0,250),las=1, main="Bony fish", ylab="# of Occurrences", xlab="Geologic Age (Ma)", col="gray95")

barplot(x[4,], ylim=c(0,250),las=1, main="Cephalopoda", ylab="# of Occurrences",xlab="Geologic Age (Ma)", col="gray95")

barplot(x[5,], ylim=c(0,250),las=1, main="Crinoidea", ylab="# of Occurrences", xlab="Geologic Age (Ma)", col="gray95")

barplot(x[6,], ylim=c(0,250),las=1, main="Echinoidea", ylab="# of Occurrences", xlab="Geologic Age (Ma)", col="gray95")

barplot(x[7,], ylim=c(0,250),las=1, main="Gastropoda", ylab="# of Occurrences", xlab="Geologic Age (Ma)", col="gray95")

barplot(x[8,], ylim=c(0,250), las=1, main="Ostracoda", ylab="# of Occurrences", xlab="Geologic Age (Ma)", col="gray95")

barplot(x[9,], ylim=c(0,250), las=1, main="Rhynchonellata", ylab="# of Occurrences", xlab="Geologic Age (Ma)", col="gray95")

barplot(x[10,], ylim=c(0,250), las=1, main="Strophomenata", ylab="# of Occurrences", xlab="Geologic Age (Ma)", col="gray95" )

barplot(x[11,], ylim=c(0,250), las=1,main="Trilobita", ylab="# of Occurrences", xlab="Geologic Age (Ma)", col="gray95")
------------------------------------------------------------------------------------------

totalOccs <- read.csv("pbdb.cleaned.classes.no.body.size.csv")
sizeOccs <- read.csv("pbdb.cleaned.classes.csv")  

gencoverage <- data.frame(matrix(nrow=11, ncol=length(unique(sizeOccs$max_ma))))
colnames(gencoverage) <- rev(sort(unique(sizeOccs$max_ma)))
rownames(gencoverage) <- c("All Genera", "Bivalvia", "bony fish", "Cephalopoda", "Crinoidea", "Echinoidea", "Gastropoda", "Ostracoda", "Rhynchonellata", "Strophomenata", "Trilobita")

x <- table(unlist(sizeOccs$genus, sizeOccs$class), sizeOccs$max_ma) #table genus occurrences by interval
x <- x[, match(colnames(x), rev(colnames(x)))]
x[x > 0] <- 1 # set counts > 1 equal to 1; presence-absence matrix

y <- table(unlist(totalOccs$genus), totalOccs$max_ma) #table genus occurrences by interval
y <- y[, match(colnames(y), rev(colnames(y)))]
x[x > 0] <- 1 # set counts > 1 equal to 1; presence-absence matrix
ages <- colnames(x) 

coverage <- (colSums(x)/colSums(y))*100

gencoverage[1,] <- coverage 

barplot(coverage, ylim=c(0,100), las=1, main="Coverage of Genera with Body Size in PBDB", ylab="% Genera Coverage", xlab="Geologic Age (Ma)", col="gray95")

x <- table(unlist(sizeOccs$genus, sizeOccs$class), sizeOccs$max_ma) #table genus occurrences by interval
x <- x[, match(colnames(x), rev(colnames(x)))]
x[x > 0] <- 1 # set counts > 1 equal to 1; presence-absence matrix

rawgenera <- colSums(x)

barplot(rawgenera, las=1, main="Genera with Body Size Data and PBDB Occurrences", ylab="Number of Genera", xlab="Geologic Age (Ma)", col="gray95")


library(plyr)
size.genera <- ddply(sizeOccs, c("class", "max_ma"), function(df)c(length(unique(df$genus)))) #find median and mean size by class
names(size.genera) <- c("class", "max_ma", "ngenera") #name output dataset

size.genera <- reshape(size.genera, idvar="class", timevar="max_ma", direction="wide")
size.genera <- size.genera[, match(colnames(size.genera), rev(colnames(size.genera)))]
size.genera <- size.genera[,c(88,1:87)]
rownames(size.genera) <- size.genera[,1]
size.genera[,1] <- NULL
colnames(size.genera) <- ages
size.genera[is.na(size.genera)] <- 0


tot.genera <- ddply(totalOccs, c("class", "max_ma"), function(df)c(length(unique(df$genus)))) #find median and mean size by class
names(tot.genera) <- c("class", "max_ma", "ngenera") #name output dataset

tot.genera <- reshape(tot.genera, idvar="class", timevar="max_ma", direction="wide")
tot.genera <- tot.genera[, match(colnames(tot.genera), rev(colnames(tot.genera)))]
tot.genera <- tot.genera[,c(88,1:87)]
rownames(tot.genera) <- tot.genera[,1]
tot.genera[,1] <- NULL
colnames(tot.genera) <- ages
tot.genera[is.na(tot.genera)] <- 0

gencoverage <- (size.genera/tot.genera)*100

par(mfrow=c(5,2))

barplot(as.matrix(gencoverage[1,])), ylim=c(0,100),las=1, main="Bivalvia", ylab="% Genera Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(gencoverage[2,]), ylim=c(0,100), las=1, main="Bony fish", ylab="% Genera Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(gencoverage[3,]), ylim=c(0,100), las=1, main="Cephalopoda", ylab="% Genera Coverage",xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(gencoverage[4,]), ylim=c(0,100), las=1, main="Crinoidea", ylab="% Genera Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(gencoverage[5,]), ylim=c(0,100), las=1, main="Echinoidea", ylab="% Genera Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(gencoverage[6,]), ylim=c(0,100), las=1, main="Gastropoda", ylab="% Genera Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(gencoverage[7,]), ylim=c(0,100), las=1, main="Ostracoda", ylab="% Genera Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(gencoverage[8,]), ylim=c(0,100), las=1, main="Rhynchonellata", ylab="% Genera Coverage", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(gencoverage[9,]), ylim=c(0,100),las=1, main="Strophomenata", ylab="% Genera Coverage", xlab="Geologic Age (Ma)", col="gray95" )

barplot(as.matrix(gencoverage[10,]), ylim=c(0,100), las=1,main="Trilobita", ylab="% Genera Coverage", xlab="Geologic Age (Ma)", col="gray95")

dev.new()
par(mfrow=c(5,2))

barplot(as.matrix(size.genera[1,]), ylim=c(0,200),las=1, main="Bivalvia", ylab="# of Genera", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(size.genera[2,]), ylim=c(0,200), las=1, main="Bony fish", ylab="# of Genera", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(size.genera[3,]), ylim=c(0,200), las=1, main="Cephalopoda", ylab="# of Genera",xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(size.genera[4,]), ylim=c(0,200), las=1, main="Crinoidea", ylab="# of Genera", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(size.genera[5,]), ylim=c(0,200), las=1, main="Echinoidea", ylab="# of Genera", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(size.genera[6,]), ylim=c(0,200), las=1, main="Gastropoda", ylab="# of Genera", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(size.genera[7,]), ylim=c(0,200), las=1, main="Ostracoda", ylab="# of Genera", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(size.genera[8,]), ylim=c(0,200), las=1, main="Rhynchonellata", ylab="# of Genera", xlab="Geologic Age (Ma)", col="gray95")

barplot(as.matrix(size.genera[9,]), ylim=c(0,200),las=1, main="Strophomenata", ylab="# of Genera", xlab="Geologic Age (Ma)", col="gray95" )

barplot(as.matrix(size.genera[10,]), ylim=c(0,200), las=1,main="Trilobita", ylab="# of Genera", xlab="Geologic Age (Ma)", col="gray95")