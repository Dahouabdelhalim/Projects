#Figure 1E
library(base2grob)
library(survival)
library(survminer)

ISPY <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/ISPY/ISPY_dataset.txt", row.names = 1, header = T)
survival <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/ISPY/ISPY_survival.txt", row.names = 1, header = T)
cx_results <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/ISPY/CIBERSORTx_output/CIBERSORTx_Adjusted.txt", row.names =1, header = T)
cx_results <- cx_results[,1:10]
lmo2 <- cx_results[, "LMO2posBasal"]
groups <- rep("Low", length(lmo2))
groups[lmo2>median(lmo2)] <- "High"
ispy <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/ISPY/ISPY_pheno.txt", row.names = 1, header = T)
er_status <- as.character(unlist((ispy[,1])))


surv.object <- Surv(survival$DRFS_time, event = survival$DRFS_events)

surv.fit <- survfit(surv.object~groups)
survtest <- survdiff(surv.object~groups)
pval <- 1 - pchisq(survtest$chisq, 1)
p <- signif(pval, 2)
fit.coxph <- coxph(surv.object ~ groups)

mypalette2 <-  c("darkred", "navyblue")

pdf("lmo2exp_survival_single.pdf", width = 5, height = 6, useDingbats = FALSE)
plot(surv.fit, col = mypalette2, mark.time = TRUE, lwd = 3, las = 1, cex.axis = 1)
title(main = "Basal LMO2+", line = 2)
title(main= bquote(italic(P)~value == .(pval)), line=1)
legend("bottomleft", title = "Clusters", legend = c("High", "Low"), lty = 1, lwd = 3, col = mypalette2, cex = 1.5)
dev.off()

surv.fit <- survfit(surv.object~groups+er_status)
survtest <- survdiff(surv.object~groups+er_status)
pval <- 1 - pchisq(survtest$chisq, 1)
p <- signif(pval, 2)
fit.coxph <- coxph(surv.object ~ groups + er_status)

mypalette2 <-  c("darkred", "orangered",  "navyblue", "dodgerblue")

pdf("lmo2exp_survival_double.pdf", width = 5, height = 6, useDingbats = FALSE)
plot(surv.fit, col = mypalette2, mark.time = TRUE, lwd = 3, las = 1, cex.axis = 1)
title(main = "Basal LMO2+", line = 2) 
title(main= bquote(italic(P)~value == .(pval)), line=1)
legend("bottomleft", title = "Clusters", legend = c("HighER-", "HighER+","LowER-", "LowER+"), lty = 1, lwd = 3, col = mypalette2, cex = 1.5)
dev.off()

p <- base2grob(survPlot)
