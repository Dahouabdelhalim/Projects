#barplot sizes
setwd("D:/Research/Lessepsian/Data_analysis/Annihilation/")

###### preparing input file ######
require(readxl)
sizes <- read_excel("input_files/size_results_ISRAEL_TOT_2020.02.19.xlsx", na="NA")

sizes <- subset(sizes, sizes$n>=10) #removes species with less than 10 specimens
sizes <- subset(sizes, sizes$Component!="Unknown") #removes species whose native/NIS status is unknown

###### computes metric ######
sizes$plot <- round(-(sizes$Max_sizeLit - sizes$Max_sizeL)/sizes$Max_sizeLit, 3)
sizes$plotIL <- round(-(sizes$Max_sizeIL - sizes$Max_sizeL)/sizes$Max_sizeIL, 3)

#write.csv(sizes2, file="results/sizes_with_taxonomy_2019.12.27.csv")

sizes_int_hard <- subset(sizes, sizes$Group=="Intertidal_hard")
sizes_sh_soft <- subset(sizes, sizes$Group=="Shallow_soft")
sizes_sh_hard <- subset(sizes, sizes$Group=="Shallow_hard")
sizes_dp_hard <- subset(sizes, sizes$Group=="Deep_hard")

to_plot <- list(sizes_int_hard$plot[sizes_int_hard$Component=="Native"],
                sizes_int_hard$plot[sizes_int_hard$Component=="NIS"],
                sizes_sh_soft$plot[sizes_sh_soft$Component=="Native"], 
                sizes_sh_soft$plot[sizes_sh_soft$Component=="NIS"],
                sizes_sh_hard$plot[sizes_sh_hard$Component=="Native"],
                sizes_sh_hard$plot[sizes_sh_hard$Component=="NIS"],
                NA,NA) #no mesophotic size data

to_plot <- list(sizes_int_hard$plotIL[sizes_int_hard$Component=="Native"],
                sizes_int_hard$plotIL[sizes_int_hard$Component=="NIS"],
                sizes_sh_soft$plotIL[sizes_sh_soft$Component=="Native"], 
                sizes_sh_soft$plotIL[sizes_sh_soft$Component=="NIS"],
                sizes_sh_hard$plotIL[sizes_sh_hard$Component=="Native"],
                sizes_sh_hard$plotIL[sizes_sh_hard$Component=="NIS"])

save(to_plot, sizes_int_hard, sizes_sh_hard, sizes_sh_soft, sizes_dp_hard, file="results/to_plot_depth_20200219.Rdata")



###### version with n under each boxplot ######
win.metafile(filename="plots/annihilation_LDsize_20200424.wmf", width=6, height=7)
#par(bg="black", fg="white")
boxplot(to_plot, 
        names=c(paste("native, n=", length(sizes_int_hard$Max_sizeL[sizes_int_hard$Component=="Native"]), sep=""), 
                paste("NIS, n=", length(sizes_int_hard$Max_sizeL[sizes_int_hard$Component=="NIS"]), sep=""), 
                paste("native, n=", length(sizes_sh_hard$Max_sizeL[sizes_sh_hard$Component=="Native"]), sep=""), 
                paste("NIS, n=", length(sizes_sh_hard$Max_sizeL[sizes_sh_hard$Component=="NIS"]), sep=""),
                paste("native, n=", length(sizes_sh_soft$Max_sizeL[sizes_sh_soft$Component=="Native"]), sep=""), 
                paste("NIS, n=", length(sizes_sh_soft$Max_sizeL[sizes_sh_soft$Component=="NIS"]), sep=""),
                "Native","NIS"), 
        border=c("blue", "red"), ylim=c(-1,1)) #main="LD max size ratio")
abline(h=0, col="grey"); abline(h=c(-0.5,0.5), lty=3, col="grey")
boxplot(to_plot, 
        names=c(paste("native, n=", length(sizes_int_hard$Max_sizeL[sizes_int_hard$Component=="Native"]), sep=""), 
                paste("NIS, n=", length(sizes_int_hard$Max_sizeL[sizes_int_hard$Component=="NIS"]), sep=""), 
                paste("native, n=", length(sizes_sh_hard$Max_sizeL[sizes_sh_hard$Component=="Native"]), sep=""), 
                paste("NIS, n=", length(sizes_sh_hard$Max_sizeL[sizes_sh_hard$Component=="NIS"]), sep=""),
                paste("native, n=", length(sizes_sh_soft$Max_sizeL[sizes_sh_soft$Component=="Native"]), sep=""), 
                paste("NIS, n=", length(sizes_sh_soft$Max_sizeL[sizes_sh_soft$Component=="NIS"]), sep=""),
                "Native","NIS"), 
        border=c("blue", "red"), ylim=c(-1,1), #main="LD max size ratio", 
        #col.axis="white", col.main="white", 
        col=rgb(1,1,1,alpha=0),
        #boxlwd=2, medlwd=2, whisklwd=2, staplelwd=2, outlwd=2,
        add=TRUE)
mtext(1, text=c("Rocky intertidal", "Rocky subtidal", "Soft subtidal", "Mesophotic"), at=c(1.5, 3.5, 5.5, 7.5), cex=1, line=2.5)
dev.off()
###### Wilcoxon test ######
#H0 distributions are the same
#p<0.05 => rejection of H0
###lit
wilcox.test(plot ~ Component, data=sizes_int_hard, 
            subset=(Component %in% c("Native", "NIS")))

wilcox.test(plot ~ Component, data=sizes_sh_hard, 
            subset=(Component %in% c("Native", "NIS")))

wilcox.test(plot ~ Component, data=sizes_sh_soft, 
            subset=(Component %in% c("Native", "NIS")))

#litIL
wilcox.test(plotIL ~ Component, data=sizes_int_hard, 
            subset=(Component %in% c("Native", "NIS")))

wilcox.test(plotIL ~ Component, data=sizes_sh_hard, 
            subset=(Component %in% c("Native", "NIS")))

wilcox.test(plotIL ~ Component, data=sizes_sh_soft, 
            subset=(Component %in% c("Native", "NIS")))

