#plotting iNEXT results
setwd("D:/Research/Lessepsian/Data_analysis/Annihilation/")

require(readxl)
results <- read_excel("input_files/iNEXT_results.xlsx", na="NA"); #results$LD_ratio <- round(results$LD_ratio, 2)
results$LD_ratio <- results$LD_ratio*100 #to plot as percentages

#select spatial scale
results <- subset(results, results$Scale=="Station") #stations only

#prepare subsets
Intertidal_hard_native <- subset(results, results$Group=="Intertidal_hard" & results$Component=="Native")
Intertidal_hard_NIS <- subset(results, results$Group=="Intertidal_hard" & results$Component=="NIS")

Subtidal_hard_native <- subset(results, results$Group=="Subtidal_hard" & results$Component=="Native")
Subtidal_hard_NIS <- subset(results, results$Group=="Subtidal_hard" & results$Component=="NIS")

Subtidal_soft_native <- subset(results, results$Group=="Subtidal_soft" & results$Component=="Native")
Subtidal_soft_NIS <- subset(results, results$Group=="Subtidal_soft" & results$Component=="NIS")

Subtidal <- subset(results, results$Depth>0 & results$Depth<80)
Subtidal_native <- subset(Subtidal, Subtidal$Component=="Native")
Subtidal_NIS <- subset(Subtidal, Subtidal$Component=="NIS")

Deep_hard_native <- subset(results, results$Group=="Deep_hard" & results$Component=="Native")
Deep_hard_NIS <- subset(results, results$Group=="Deep_hard" & results$Component=="NIS")

Deep_soft_native <- subset(results, results$Group=="Deep_soft" & results$Component=="Native")
Deep_soft_NIS <- subset(results, results$Group=="Deep_soft" & results$Component=="NIS")

Deep_native <- subset(results, results$Depth==80 & results$Component=="Native")
Deep_NIS <- subset(results, results$Depth==80 & results$Component=="NIS")

#prepare lists to input into boxplots
to_plot_native <- list(Intertidal_hard_native$LD_ratio, 
                Subtidal_hard_native$LD_ratio, Subtidal_soft_native$LD_ratio, 
                Deep_native$LD_ratio)

to_plot_NIS <- list(Intertidal_hard_NIS$LD_ratio, 
                Subtidal_hard_NIS$LD_ratio, Subtidal_soft_NIS$LD_ratio, 
                Deep_NIS$LD_ratio)

to_plot <- list(Intertidal_hard_native$LD_ratio, Intertidal_hard_NIS$LD_ratio, 
                Subtidal_hard_native$LD_ratio, Subtidal_hard_NIS$LD_ratio, 
                Subtidal_soft_native$LD_ratio, Subtidal_soft_NIS$LD_ratio, 
                Deep_native$LD_ratio, NA)

to_plot_depths <- list(Intertidal_hard_native$LD_ratio, Intertidal_hard_NIS$LD_ratio, 
                       Subtidal_native$LD_ratio, Subtidal_NIS$LD_ratio, 
                       Deep_native$LD_ratio, Deep_NIS$LD_ratio)

#boxplot white background
postscript("plots/annihilation_LD_white.eps", width=7, height=7)
win.metafile(filename="plots/annihilation_LD_white.wmf", width=7, height=7)
boxplot(to_plot, 
               names=c("", "", "", "", "", "", "", ""), 
               border=c("blue", "red","blue", "red","blue", "red","blue", "red"), ylim=c(0,150)
               )
mtext(2, text="current / historical species richness [%]", line=2.5, col="black")
mtext(1, text=c("Rocky intertidal", "Rocky subtidal", "Soft subtidal", "Mesophotic"), at=c(1.5, 3.5, 5.5, 7.5), cex=1, line=2.5)
abline(h=c(50,100), lty=c(2,1), col="grey")
boxplot(to_plot, 
        names=c("Native", "NIS", "Native", "NIS", "Native", "NIS", "Native", "NIS"), 
        border=c("blue", "red","blue", "red","blue", "red","blue", "red"), ylim=c(0,150),
        col=rgb(1,1,1,alpha=0), 
        add=TRUE)
dev.off()

