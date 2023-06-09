## SPIDER ANALYSIS SUPP: Used for analyses in Supporting Information of publication
# Written by: Jun Ying Lim and Susan Kennedy
# From publication:
# Kennedy, S. R., Lim, J. Y., Clavel, J., Krehenwinkel, H., and Gillespie, R. G. 2019. Spider webs, stable isotopes and molecular gut content analysis: Multiple lines of evidence support trophic niche differentiation in a community of Hawaiian spiders. Functional Ecology.

## DIRECTORIES ==============================
#Set "main.dir" to your own directory containing subfolders:
main.dir <- setwd() #Insert your file path in parentheses
data.dir <- file.path(main.dir, "raw_data")
res.dir <- file.path(main.dir, "results")

## PACKAGES ============================== 
options(stringsAsFactors = FALSE)
library(reshape2)
library(stringr)

## INPUT FILE ==============================
# Input raw data
isotope2014 <- read.csv(file.path(data.dir, "isotopes_2014.csv"))
isotopeData <- read.csv(file.path(data.dir, "isotopes_2014_2016_upper.csv"))

# Function to clean up Tukey test output
summaryTukeyHSD <- function(x){
  res <- as.data.frame(x$species)
  temp <- str_split_fixed(rownames(res), pattern = "-", n = 2)
  res$comp1 <- temp[,1]
  res$comp2 <- temp[,2]
  res$`p adj` <- round(res$`p adj`, digits = 5)
  res$diff <- round(res$diff, digits = 3)
  diff <- acast(res, comp1 ~ comp2, value.var = "diff")
  sig <-  acast(res, comp1 ~ comp2, value.var = "p adj")
  return(list("diff"=  diff, "sig" = sig))
}

# Test for isotopic diffs between web-builders and Spiny Legs at high elevation in 2014:
NClade <- lm (delta.15.N ~ category, data = isotope2014)
anova(NClade)

CClade <- lm (delta.13.C ~ category, data = isotope2014)
anova(CClade)


# Differences in isotopes among species, high elevation only ---------------
anovaNWeb <- aov(delta.15.N ~ species, data = subset(isotopeData, category == "web"))
summary.lm(anovaNWeb)
anovaCWeb <- aov(delta.13.C ~ species, data = subset(isotopeData, category == "web"))
summary.lm(anovaCWeb)

tukeyNWeb <- TukeyHSD(anovaNWeb)
tukeyNWeb_summary <- summaryTukeyHSD(tukeyNWeb)
write.csv(tukeyNWeb_summary$diff, file = file.path(res.dir, "tukeyNWebUpper_diff.csv"))
write.csv(tukeyNWeb_summary$sig, file = file.path(res.dir, "tukeyNWebUpper_sig.csv"))

tukeyCWeb <- TukeyHSD(anovaCWeb)
tukeyCWeb_summary <- summaryTukeyHSD(tukeyCWeb)
write.csv(tukeyCWeb_summary$diff, file = file.path(res.dir, "tukeyCWebUpper_diff.csv"))
write.csv(tukeyCWeb_summary$sig, file = file.path(res.dir, "tukeyCWebUpper_sig.csv"))

anovaNSpiny <- aov(delta.15.N ~ species, data = subset(isotopeData, category == "spiny"))
summary.lm(anovaNSpiny)
anovaCSpiny <- aov(delta.13.C ~ species, data = subset(isotopeData, category == "spiny"))
summary.lm(anovaCSpiny)

tukeyNSpiny <- TukeyHSD(anovaNSpiny)
tukeyNSpiny_summary <- summaryTukeyHSD(tukeyNSpiny)
write.csv(tukeyNSpiny_summary$diff, file = file.path(res.dir, "tukeyNSpinyUpper_diff.csv"))
write.csv(tukeyNSpiny_summary$sig, file = file.path(res.dir, "tukeyNSpinyUpper_sig.csv"))

tukeyCSpiny <- TukeyHSD(anovaCSpiny)
tukeyCSpiny_summary <- summaryTukeyHSD(tukeyCSpiny)
write.csv(tukeyCSpiny_summary$diff, file = file.path(res.dir, "tukeyCSpinyUpper_diff.csv"))
write.csv(tukeyCSpiny_summary$sig, file = file.path(res.dir, "tukeyCSpinyUpper_sig.csv"))

