## SPIDER ANALYSIS 2: To run after "spiders_calcOverlap.R"
# Written by: Jun Ying Lim and Susan Kennedy
# From publication:
# Kennedy, S. R., Lim, J. Y., Clavel, J., Krehenwinkel, H., and Gillespie, R. G. 2019. Spider webs, stable isotopes and molecular gut content analysis: Multiple lines of evidence support trophic niche differentiation in a community of Hawaiian spiders. Functional Ecology.

## DIRECTORIES ==============================
#Set "main.dir" to your own directory containing subfolders:
main.dir <- setwd() #Insert your file path in parentheses
data.dir <- file.path(main.dir, "raw_data")
fig.dir <- file.path(main.dir, "figures")
res.dir <- file.path(main.dir, "results")

## PACKAGES ============================== 
options(stringsAsFactors = FALSE)
library(ggplot2)
library(reshape2)
library(stringr)
library(cowplot)

## INPUT FILE ==============================
# Input raw data
webData <- read.csv(file.path(data.dir, "webs_2013_2014.csv"))
isotopeData <- read.csv(file.path(data.dir, "isotopes_2014_2016.csv"))
webData$species <- gsub(webData$species, pattern = "T. ", replacement = "")
webAnalysis <- read.csv(file.path(res.dir, "webAnalysis.csv"))

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

# Differences in webs ---------------
aovmod <- aov(PCA1 ~ species, data = webAnalysis)
summary.lm(aovmod)

aovmod2 <- aov(PCA2 ~ species, data = webAnalysis)
summary.lm(aovmod2)

tukeyPCA1 <- TukeyHSD(aovmod)
tukeyPCA1summary <- summaryTukeyHSD(tukeyPCA1)
write.csv(tukeyPCA1summary$diff, file = file.path(res.dir, "tukeyPCA1_diff.csv"))
write.csv(tukeyPCA1summary$sig, file = file.path(res.dir, "tukeyPCA1_sig.csv"))

tukeyPCA2 <- TukeyHSD(aovmod2)
tukeyPCA2_summary <- summaryTukeyHSD(tukeyPCA2)
write.csv(tukeyPCA2_summary$diff, file = file.path(res.dir, "tukeyPCA2_diff.csv"))
write.csv(tukeyPCA2_summary$sig, file = file.path(res.dir, "tukeyPCA2_sig.csv"))

# Differences in site choice ---------------
webSites <- melt(webData, id.vars = "species", measure.vars = "attveg1.group")
webSites_count <- table(webSites$species, webSites$value)
webSites_freq <- webSites_count / rowSums(webSites_count)

webSites_gg <- melt(webSites_freq)

webSites_count_melt <- melt(webSites_count)
chisq.test(t(webSites_count), simulate.p.value = TRUE)

fisher.test(webSites_count, simulate.p.value = TRUE)

angleMod <- aov(angle ~ species, data = webData)
summary.lm(angleMod)
heightMod <- aov(height ~ species, data = webData)
summary.lm(heightMod)

# Test for interspecific diffs in height & angle
heightTukey <- TukeyHSD(heightMod)
heightTukey_summary <- summaryTukeyHSD(heightTukey)
write.csv(heightTukey_summary$diff, file = file.path(res.dir, "tukeyHeight_diff.csv"))
write.csv(heightTukey_summary$sig, file = file.path(res.dir, "tukeyHeight_sig.csv"))

angleTukey <- TukeyHSD(angleMod)
angleTukey_summary <- summaryTukeyHSD(angleTukey)
write.csv(angleTukey_summary$diff, file = file.path(res.dir, "tukeyAngle_diff.csv"))
write.csv(angleTukey_summary$sig, file = file.path(res.dir, "tukeyAngle_sig.csv"))

# Differences in isotopes across species ---------------
anovaNWeb <- aov(delta.15.N ~ species, data = subset(isotopeData, category == "web"))
summary.lm(anovaNWeb) # delta 15N differed significantly across species
anovaCWeb <- aov(delta.13.C ~ species, data = subset(isotopeData, category == "web"))
summary.lm(anovaCWeb) # delta 13C differed significantly across species

tukeyNWeb <- TukeyHSD(anovaNWeb)
tukeyNWeb_summary <- summaryTukeyHSD(tukeyNWeb)
write.csv(tukeyNWeb_summary$diff, file = file.path(res.dir, "tukeyNWeb_diff.csv"))
write.csv(tukeyNWeb_summary$sig, file = file.path(res.dir, "tukeyNWeb_sig.csv"))

tukeyCWeb <- TukeyHSD(anovaCWeb)
tukeyCWeb_summary <- summaryTukeyHSD(tukeyCWeb)
write.csv(tukeyCWeb_summary$diff, file = file.path(res.dir, "tukeyCWeb_diff.csv"))
write.csv(tukeyCWeb_summary$sig, file = file.path(res.dir, "tukeyCWeb_sig.csv"))

anovaNSpiny <- aov(delta.15.N ~ species, data = subset(isotopeData, category == "spiny"))
summary.lm(anovaNSpiny)# delta 15N differed significantly across species
anovaCSpiny <- aov(delta.13.C ~ species, data = subset(isotopeData, category == "spiny"))
summary.lm(anovaCSpiny)# delta 13C differed significantly across species

tukeyNSpiny <- TukeyHSD(anovaNSpiny)
tukeyNSpiny_summary <- summaryTukeyHSD(tukeyNSpiny)
write.csv(tukeyNSpiny_summary$diff, file = file.path(res.dir, "tukeyNSpiny_diff.csv"))
write.csv(tukeyNSpiny_summary$sig, file = file.path(res.dir, "tukeyNSpiny_sig.csv"))

tukeyCSpiny <- TukeyHSD(anovaCSpiny)
tukeyCSpiny_summary <- summaryTukeyHSD(tukeyCSpiny)
write.csv(tukeyCSpiny_summary$diff, file = file.path(res.dir, "tukeyCSpiny_diff.csv"))
write.csv(tukeyCSpiny_summary$sig, file = file.path(res.dir, "tukeyCSpiny_sig.csv"))


## HYPERVOLUME OVERLAP CORRELATION TESTS AND PLOT ==============================  
webHypervolOverlapJ <- readRDS(file.path(res.dir, "webHypervolOverlapJ.rds"))

isotopeHypervolOverlapWebJ <- readRDS(file.path(res.dir, "isotopeHypervolOverlapWebJ.rds"))
isotopeHypervolOverlapSpinyJ <- readRDS(file.path(res.dir, "isotopeHypervolOverlapSpinyJ.rds"))

gutBetaWebES <- readRDS(file.path(res.dir, "gutBetaWebES.rds"))
gutBetaSpinyES <- readRDS(file.path(res.dir, "gutBetaSpinyES.rds"))

webHypervolOverlapJ <- as.dist(webHypervolOverlapJ)
isotopeHypervolOverlapWebJ <- as.dist(isotopeHypervolOverlapWebJ)
isotopeHypervolOverlapSpinyJ <- as.dist(isotopeHypervolOverlapSpinyJ)
gutBetaWebES <- as.dist(gutBetaWebES)
gutBetaSpinyES <- as.dist(gutBetaSpinyES)

overlapResultsWebJES <- data.frame(webOverlapJ = as.vector(webHypervolOverlapJ), isotopeOverlapWebJ = as.vector(isotopeHypervolOverlapWebJ), gutBetaWebES = as.vector(gutBetaWebES))

overlapResultsSpinyJES <- data.frame(isotopeOverlapSpinyJ = as.vector(isotopeHypervolOverlapSpinyJ), gutBetaSpinyES = as.vector(gutBetaSpinyES))

corValES <- cor.test(overlapResultsWebJES$webOverlapJ, overlapResultsWebJES$isotopeOverlapWebJ)
corStatPlot <- paste0("Pearson's corr. = ", round(corValES$estimate,2), "\\n", "p = ", round(corValES$p.value,2), "\\n")
webByIsotopePlotES <- ggplot(data = overlapResultsWebJES, aes(y = webOverlapJ, x = isotopeOverlapWebJ)) + geom_point() + annotate("text", x = Inf, y = -Inf, label = corStatPlot, hjust = 1, vjust = 0, parse = FALSE) + labs(y = "Web architecture overlap\\n(Jaccard)", x = "Isotopic overlap\\n(Jaccard)") + geom_smooth(method = "lm", se = FALSE)

corValES <- cor.test(overlapResultsWebJES$webOverlapJ, overlapResultsWebJES$gutBetaWebES)
corStatPlot <- paste0("Pearson's corr. = ", round(corValES$estimate,2), "\\n", "p = ", round(corValES$p.value,2), "\\n")
webByGutPlotES <- ggplot(data = overlapResultsWebJES, aes(y = webOverlapJ, x = gutBetaWebES)) + geom_point() + annotate("text", x = Inf, y = Inf, label = corStatPlot, hjust = 1, vjust = 1, parse = FALSE) + labs(y = "Web architecture overlap\\n(Jaccard)", x = "Gut content dissimilarity \\n (Bray-Curtis)") + geom_smooth(method = "lm", se = FALSE)

corValES <- cor.test(overlapResultsWebJES$isotopeOverlapWebJ, overlapResultsWebJES$gutBetaWebES)
corStatPlot <- paste0("Pearson's corr. = ", round(corValES$estimate,2), "\\n", "p = ", round(corValES$p.value,2), "\\n")
isoByGutWebPlotES <- ggplot(data = overlapResultsWebJES, aes(y = isotopeOverlapWebJ, x = gutBetaWebES)) + geom_point() + annotate("text", x = Inf, y = Inf, label = corStatPlot, hjust = 1, vjust = 1, parse = FALSE) + labs(y = "Isotopic overlap\\n(Jaccard)", x = "Gut content dissimilarity \\n (Bray-Curtis)") + geom_smooth(method = "lm", se = FALSE)

corValES <- cor.test(overlapResultsSpinyJES$isotopeOverlapSpinyJ, overlapResultsSpinyJES$gutBetaSpinyES)
corStatPlot <- paste0("Pearson's corr. = ", round(corValES$estimate,2), "\\n", "p = ", round(corValES$p.value,2), "\\n")
isoByGutSpinyPlotES <- ggplot(data = overlapResultsSpinyJES, aes(y = isotopeOverlapSpinyJ, x = gutBetaSpinyES)) + geom_point() + annotate("text", x = Inf, y = Inf, label = corStatPlot, hjust = 1, vjust = 1, parse = FALSE) + labs(y = "Isotopic overlap\\n(Jaccard)", x = "Gut content dissimilarity \\n (Bray-Curtis)") + geom_smooth(method = "lm", se = FALSE)

overlapPlotCombinedES <- plot_grid(plotlist = list(webByIsotopePlotES, webByGutPlotES, isoByGutWebPlotES, isoByGutSpinyPlotES), labels = c("a)","b)", "c)", "d)"), nrow = 2, ncol = 2)
ggsave(overlapPlotCombinedES, filename = file.path(fig.dir, "fig7_overlapcombinedES.pdf"), height = 7, width = 9)

