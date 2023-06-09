# # Note to user: This analysis requires packages "compositions" and "robCompositions".  For convenience, we use the package "xlsx",
# allowing users to directly call appendices provided online with the manuscript.  Saving *.xlsx appendices to text files
# and using read.table() or read.csv() is also an option in R-base without requiring installation of package "xlsx".

library("xlsx")
library("compositions")
library("robCompositions")

# specify working directory where appendices have been saved
setwd("E:/Smyser_FeralSwineAncestry")

FeralSwine_K17_Ancestry <- read.xlsx("S8_InvasiveFeralSwine_CompleteK17RefereceSet_Ancestry.xlsx", sheetName = "S8_FeralSwine_K17ReferenceSet")[,1:19]
head(FeralSwine_K17_Ancestry)

FeralSwine_Metadata <- read.xlsx("S1_InvasiveFeralSwine_Metadata.xlsx", sheetName = "S1_InvasiveFeralSwine_Metadata") 
head(FeralSwine_Metadata)

# pair ancestry values with metadata
FeSw_Ancestry_Metadata <- merge(x = FeralSwine_Metadata, y = FeralSwine_K17_Ancestry, by.x = "Individual.ID", by.y = "Individual.ID")
# reorder data frame 
FeSw_Ancestry_Metadata <- FeSw_Ancestry_Metadata[order(FeSw_Ancestry_Metadata$Order.x), ] 
# remove redundant Order column from the ancestry data frame
FeSw_Ancestry_Metadata <- subset(FeSw_Ancestry_Metadata, select = -(Order.y))

# view combined data frame
head(FeSw_Ancestry_Metadata)

# return the number of individuals from long-established versus newly-emergent populations
table(FeSw_Ancestry_Metadata$Population.History)

# separate data frame into a compositional matrix of ancestry and invasion history for construction of a model comparing ancestry between long-established and newly-emergent populations
ancestry <- acomp(FeSw_Ancestry_Metadata[,9:25])
invasion <- FeSw_Ancestry_Metadata$Population.History

# run the compositional model to evaluate whether there are differences in ancestry between feral swine in long-established and newly-emergent populations
model = lm(ilr(ancestry) ~ invasion)
anova(model)

# given that the overall model demonstrated significant differences in ancestry between long-established and newly-emergent populations, 
# we are interetested in whether those differences were attributable to ancestral contributions from European wild boar, Western heritage breeds,
# or commercial breeds.  We need to reconfigure the ancestry matrix so that there are 4 ancestry groups: European wild boar, Western heritage breeds,
# commercial breeds, and Other (summation of the remaining ancestry groups).

ancestry_4 <- FeSw_Ancestry_Metadata[, c("Individual.ID", "Population.History")]
ancestry_4$WildBoar <- FeSw_Ancestry_Metadata$Reference.Cluster.17 
ancestry_4$HeritagePig <- FeSw_Ancestry_Metadata$Reference.Cluster.16
ancestry_4$Commercial <- rowSums(FeSw_Ancestry_Metadata[,c(9, 10, 16, 17, 21)])
ancestry_4$Other <- rowSums(FeSw_Ancestry_Metadata[,c(11, 12, 13, 14, 15, 18, 19, 20, 22, 23)])

ancestry_4_comp <- acomp(ancestry_4[,3:6])
# view reduced ancestry composition matrix
ancestry_4_comp[1:6,]

invasion <- ancestry_4$Population.History

# Using the nomenclature presented in  Filzmoser et al. 2018 (section 10.7.1)
# Filzmoser, P., Hron, K., & Temple, M. (2018) Applied compositional data analysis - with worked examples in R. 
#     Cham, Switzerland: Springer. doi 10.1007/978-3-319-96422-5

Regression_models <- vector("list", ncol(ancestry_4_comp))

for (j in 1:ncol(ancestry_4_comp)){
  zj <- pivotCoord(ancestry_4_comp, pivotvar = j)
  # use only first coordinate
  res <- lm(zj[,1] ~ invasion)
  # result for the first coordinate
  Regression_models[[j]] <- summary(res)
}

Regression_models

