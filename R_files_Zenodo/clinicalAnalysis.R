library(survival)
library(survminer)

###ISPY###
#ispy_mat <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/ISPY/ISPY_dataset.txt", row.names = 1, header = T)
ispy_survival <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/ISPY/ISPY_survival.txt", row.names = 1, header = T)
ispy_pheno <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/ISPY/ISPY_pheno.txt", row.names = 1, header = T)
ispy_cx_results <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/ISPY/CIBERSORTx_output/CIBERSORTx_Adjusted.txt", row.names =1, header = T)
lmo2 <- ispy_cx_results[, "LMO2posBasal"]
ispy_pheno2 <- readRDS("/data1/ggulati/mango/breastdata/clinical/ISPY/ISPY_metadata_temp.rds")

##boxplot
coolBoxplot <- function(y = NULL,  x= NULL, ylab = "y-axis", ymin = NULL, ymax = NULL, cols = "purple"){
tog = ifelse(is.null(ymin), round(min(y) - 0.1, 1), ymin)
tog2 = ifelse(is.null(ymax), round(max(y), 1), ymax)
par(oma = c(2,3,4,2))
par(mar = c(7,6.75,1,1), xpd = NA)
boxplot(y~x, outline = F, xaxt = "n", yaxt = "n",
        staplelwd = 1,
        medlwd = 2.5,
        whisklty = 1,
        border = "black",
        col = adjustcolor(cols, alpha.f = 0.6),
        ylab = "",
        xlab = "",
        cex.lab = 1.75,
        las = 1,
        frame.plot = F,
        ylim = c(tog,tog2),
        cex.axis = 2)
mtext(ylab,
      cex = 2, side = 2, line = 5.5)
title(main = "",font.main = 1, cex.main = 2.25, line = 3)
axis(1, pos = tog, at = 1:length(levels(x)), labels = FALSE)
axis(2, pos = 0, at = signif(as.numeric(formatC(seq(tog, tog2, 0.05), format = "f"), 2)), las = 1, cex.axis = 2)
segments(-0.15,tog,length(levels(x))+0.5,tog)
text(seq_along(1:length(levels(x))), par("usr")[3] - 0.0125, labels = levels(x),
     srt = 45, adj = 1, xpd = TRUE, cex = 2)
}


#LMO2 abundance by stage
stage <- ispy_pheno2$clinical_ajcc_stage
stage[-which(stage %in% c("I", "IIA", "IIB", "IIIA", "IIIB", "IIIC"))] <- NA
stage <- factor(stage, levels = c("I", "IIA", "IIB", "IIIA", "IIIB", "IIIC"))
pdf("ISPY_Stage.pdf", width = 10, height = 6)
coolBoxplot(y = lmo2, x = stage, ymin = 0, ymax = 0.35, ylab = "Estimated abundance of \\nLMO2+ basal cells")
dev.off()

summary(aov(lmo2~stage))
#LMO2 abundance by Pam50
pam50 <- ispy_pheno2$pam50_class
pam50 <- factor(pam50, levels = c("Normal", "LumA", "LumB", "Her2", "Basal"))
pdf("ISPY_pam50.pdf", width = 8, height = 6)
coolBoxplot(y = lmo2, x = pam50, ymin = 0, ymax = 0.35, ylab = "Estimated abundance of \\nLMO2+ basal cells")
dev.off()

summary(aov(lmo2~pam50))
###Metabric Discovery###
#metabric_d_mat <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/Metabric_discovery/Metabric_discovery_data.txt", row.names = 1, header = T)
metabric_d_survival <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/Metabric_discovery/Metabric_discovery_survival.txt", row.names = 1, header = T)
metabric_d_pheno <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/Metabric_discovery/Metabric_discovery_pheno.txt", row.names = 1, header = T, sep = "\\t")
metabric_d_cx_results <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/Metabric_discovery/CIBERSORTx_output/CIBERSORTx_Adjusted.txt", row.names =1, header = T)
lmo2 <- metabric_d_cx_results[, "LMO2posBasal"]

#Grade
grade <- metabric_d_pheno$Grade
grade <- factor(grade)
pdf("Metabric_D_Grade.pdf", width = 6, height = 6)
coolBoxplot(y = lmo2, x = grade, ymin = 0, ymax = 0.15, ylab = "Estimated abundance of \\nLMO2+ basal cells")
dev.off()
summary(aov(lmo2~grade))

#Pam50
pam50 <- metabric_d_pheno$pam50
pam50 <- factor(pam50, levels = c("Normal", "LumA", "LumB", "Her2", "Basal"))
pdf("Metabric_D_PAM50.pdf", width = 8, height = 6)
coolBoxplot(y = lmo2, x = pam50, ymin = 0, ymax = 0.15, ylab = "Estimated abundance of \\nLMO2+ basal cells")
dev.off()
summary(aov(lmo2~pam50))

#Overall survival
groups <- rep("Low", length(lmo2))
groups[lmo2>median(lmo2)] <- "High"
er_status <- metabric_d_pheno$ER_IHC_status
er_status[-which(er_status %in% c("Negative", "Positive"))] <- NA
surv.object <- Surv(metabric_d_survival$OS_time, metabric_d_survival$OS_events)

#
surv.fit <- survfit(surv.object~groups)
survtest <- survdiff(surv.object~groups)
pval <- 1 - pchisq(survtest$chisq, 1)
p <- signif(pval, 2)
fit.coxph <- coxph(surv.object ~ groups)
mypalette2 <-  c("darkred", "navyblue")
pdf("metabric_d_survival_single.pdf", width = 5, height = 6, useDingbats = FALSE)
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

###Metabric Validation###
#metabric_v_mat <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/Metabric_validation/Metabric_validation_data.txt", row.names = 1, header = T)
metabric_v_survival <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/Metabric_validation/Metabric_validation_survival.txt", row.names = 1, header = T)
metabric_v_pheno <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/Metabric_validation/Metabric_validation_pheno.txt", row.names = 1, header = T, sep = "\\t")
metabric_v_cx_results <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/Metabric_validation/CIBERSORTx_output/CIBERSORTx_Adjusted.txt", row.names =1, header = T)
lmo2 <- metabric_v_cx_results[, "LMO2posBasal"]

#Gradesummary(aov(lmo2~pam50))

grade <- metabric_v_pheno$Grade
grade <- factor(grade)
pdf("Metabric_Grade.pdf", width = 6, height = 6)
coolBoxplot(y = lmo2, x = grade, ymin = 0, ymax = 0.15, ylab = "Estimated abundance of \\nLMO2+ basal cells")
dev.off()
summary(aov(lmo2~grade))

#Pam50
pam50 <- metabric_v_pheno$pam50
pam50 <- factor(pam50, levels = c("Normal", "LumA", "LumB", "Her2", "Basal"))
pdf("Metabric_PAM50.pdf", width = 8, height = 6)
coolBoxplot(y = lmo2, x = pam50, ymin = 0, ymax = 0.15, ylab = "Estimated abundance of \\nLMO2+ basal cells")
dev.off()
summary(aov(lmo2~pam50))

#Overall survival
groups <- rep("Low", length(lmo2))
groups[lmo2>median(lmo2)] <- "High"
er_status <- metabric_v_pheno$ER_IHC_status
er_status[-which(er_status %in% c("Negative", "Positive"))] <- NA
surv.object <- Surv(metabric_v_survival$OS_time, metabric_v_survival$OS_events)

#
surv.fit <- survfit(surv.object~groups)
survtest <- survdiff(surv.object~groups)
pval <- 1 - pchisq(survtest$chisq, 1)
p <- signif(pval, 2)
fit.coxph <- coxph(surv.object ~ groups)
mypalette2 <-  c("darkred", "navyblue")
pdf("metabric_v_survival_single.pdf", width = 5, height = 6, useDingbats = FALSE)
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


###TCGA###
#tcga_mat <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/TCGA/TCGA_data.txt", row.names = 1, header = T)
tcga_survival <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/TCGA/TCGA_survival.txt", row.names = 1, header = T)
tcga_pheno <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/TCGA/TCGA_pheno.txt", row.names = 1, header = T, sep = "\\t")
tcga_cx_results <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/TCGA/CIBERSORTx_output/CIBERSORTx_Adjusted.txt", row.names =1, header = T)
lmo2 <- tcga_cx_results[, "LMO2posBasal"]

#Stage
stage <- tcga_pheno$Stage
stage[which(stage %in% c("Stage Tis", "Stage I", "Stage II", "Stage III"))] <- NA
stage <- factor(stage)
pdf("TCGA_Stage.pdf", width = 10, height = 6)
coolBoxplot(y = lmo2, x = stage, ymin = 0, ymax = 0.15, ylab = "Estimated abundance of \\nLMO2+ basal cells")
dev.off()

summary(aov(lmo2~stage))

#Pam50
pam50 <- tcga_pheno$PAM50
pam50 <- factor(pam50, levels = c("Normal", "LumA", "LumB", "Her2", "Basal"))
pdf("TCGA_PAM50.pdf", width = 8, height = 6)
coolBoxplot(y = lmo2, x = pam50, ymin = 0, ymax = 0.15, ylab = "Estimated abundance of \\nLMO2+ basal cells")
dev.off()

summary(aov(lmo2~pam50))
#Overall survival
groups <- rep("Low", length(lmo2))
groups[lmo2>median(lmo2)] <- "High"
er_status <- tcga_pheno$ER_score
er_status[-which(er_status %in% c("Negative", "Positive"))] <- NA
surv.object <- Surv(tcga_survival$OS_time, tcga_survival$OS_events)

#
surv.fit <- survfit(surv.object~groups)
survtest <- survdiff(surv.object~groups)
pval <- 1 - pchisq(survtest$chisq, 1)
p <- signif(pval, 2)
fit.coxph <- coxph(surv.object ~ groups)
mypalette2 <-  c("darkred", "navyblue")
pdf("tcga_survival_single.pdf", width = 5, height = 6, useDingbats = FALSE)
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

