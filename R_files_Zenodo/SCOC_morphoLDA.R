## Sceloporus Morphological Analysis
## Written by Aryeh H. Miller, edited by Hayden R. Davis on 5 October 2021

library(MASS)
library(klaR)
library(ape)

setwd("~/<directory>")
scelops <- na.exclude(read.csv(file = "SCOC_morphological_data.csv", 
                               header = TRUE, fill = T))

#function to normalize merisistic data
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

norm_scelop <- as.data.frame(lapply(scelops[11:16], min_max_norm))

scelops_new <- cbind(scelops, norm_scelop)
scelops_new <- scelops_new[, c(1:10, 17:22)]
head(scelops_new)

#Log transform data
TlogSVL <- log(scelops_new$svl)
TlogHeadL <- log(scelops_new$hl)
TlogHeadW <- log(scelops_new$hw) 
TlogLFing <- log(scelops_new$lfl) 
TlogRFing <- log(scelops_new$rfl) 
TlogLLongToe <- log(scelops_new$lrtl)
TlogRLongToe <- log(scelops_new$lltl)

#Correct for SVL and generate residuals
TlogHeadLc <- residuals(lm(TlogHeadL ~ TlogSVL))
TlogHeadWc <- residuals(lm(TlogHeadW ~ TlogSVL))
TlogLFingc <- residuals(lm(TlogLFing ~ TlogSVL))
TlogRFingc <- residuals(lm(TlogRFing ~ TlogSVL))
TlogLLongTc <- residuals(lm(TlogLLongToe ~ TlogSVL))
TlogRLongTc <- residuals(lm(TlogRLongToe ~ TlogSVL))

#subset samples that will not be size corrected
ms <- scelops_new[, 15]
ds <- scelops_new[, 16]
lfp <- scelops_new[, 11]
rfp <- scelops_new[, 12]
ltl <- scelops_new[, 13]
rtl <- scelops_new[, 14]

head(scelops_new)

#Combine into new dataframe
datameasurements <-  droplevels(subset(scelops_new, select = c(id, sex, population,  svl,   
                      hl,   hw,  lfl,  rfl, lrtl, lltl, ltl, rtl, rfp, lfp, ms, ds)))

TTable <-  cbind(datameasurements,TlogSVL, TlogHeadL,TlogHeadW,
                 TlogLFing,TlogRFing, TlogLLongToe, TlogRLongToe, 
                 ltl, rtl, ms, ds, TlogHeadLc ,TlogHeadWc,TlogLFingc,TlogRFingc, 
                 TlogLLongTc, TlogRLongTc)

popcols <- scelops_new$population
levels(popcols) <- c("indianred4", "darkgreen", "steelblue3", "mediumorchid4", "orangered2")
cols <- as.vector(popcols)

# DFA by population
data_corr <- TTable[ , 24:ncol(TTable)]
pops_dfa <- TTable[, 3] #isolate groupings
z <- lda(TTable$population ~ ., data = data_corr)
pz <- predict(z)
lda_table <- table(pz$class, TTable$population)
write.table(lda_table, "ldaTable_meristic_corr.csv", sep = "\\t")
accuracy <- pz$class == TTable$population
pz$class[accuracy == F]
TTable$population[accuracy == F]
data.frame(original = TTable$population[accuracy == F], prediction = pz$class[accuracy == F])

scores <- data.frame(TTable$population, pz$x[])
scores
write.table(scores, file = "lda_scores_corr.csv", sep = "\\t")

##to make plot using ggplot
myCol <- c("indianred4", "darkgreen", "steelblue3", "mediumorchid4", "orangered2")

ggplot(scores, aes(x = LD1, y = LD2, color = TTable.population)) + theme_bw() + 
  theme(legend.key = element_blank()) + scale_color_manual(values = myCol) +
  scale_shape_manual(values=c(20, 20, 20)) + geom_point(aes(color = scores$TTable.population), size = 5) + 
  theme(legend.title=element_blank(), axis.text = element_text(size = 16), axis.title=element_text(size = 16),  
          legend.text = element_text(size = 20)) + guides(colour = guide_legend(override.aes = list(size = 6)))
