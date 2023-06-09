### rnrb_admixture_plots.R
### Libby Natola
### Started 7 December 2021
### Plot admixture results in a reproduceable non-excel doc way

### set working directory
setwd("~/Documents/UBC/Bioinformatics/rnrb/admixture")
library(ggplot2)
library(dplyr)


### read in data
nbc_tbl=read.table("nbc_max-miss30_filtered.snps_unrelated_filtered.birds_noeasties_refiltered.pruned.2.Q")
nbc_tbl
nbc_fam <- read.table("nbc_max-miss30_filtered.snps_unrelated_filtered.birds_noeasties_refiltered.pruned.fam")
nbc_info <- read.table("../sample_info_nbc.txt", header=TRUE)

### add ID column from the .fam file
nbc_admix_data <- cbind(nbc_fam$V2, nbc_tbl)

### rename columns
colnames(nbc_admix_data) <- c("ID", "p", "q")

### add info
nbc_admix_data <- left_join(nbc_admix_data, nbc_info)
#write out
write.csv(nbc_admix_data, "nbc_admix_data.csv")

# check admixture columns
# V1 = RNSA high, V2 = RBSA high, switch columns
nbc_col.order <- c("V2","V1")
nbc_tbl <- nbc_tbl[,nbc_col.order]

### order data by column value so they're easily visualized, but these have the bars going the opposite way of the other two (top of bar is red instead of blue)
nbc_tbl <- nbc_tbl[order(nbc_tbl$V1,nbc_tbl$V2),]
nbc_tbl

### visualize
barplot(t(as.matrix(nbc_tbl)), col=c( "#FF4343", "#4343FF"),
        xlab="Individual", ylab="Ancestry", border=NA, space=0, xaxt='n', main = "NBC ADMIXTURE")
box(col=rgb(253/255, 231/255, 37/255), lwd=5, which = "plot")

pdf(file = "nbc_tassel_admixture_plot_k2.pdf")
barplot(t(as.matrix(nbc_tbl)), col=c( "#FF4343", "#4343FF"),
        xlab="Individual", ylab="Ancestry", border=NA, space=0, xaxt='n', main = "NBC ADMIXTURE")
dev.off()


### read in data
manning_tbl=read.table("man_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.2.Q")
manning_tbl
man_fam <- read.table("man_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.fam")
man_info <- read.table("../sample_info_man.txt", header=TRUE)

### add ID column from the .fam file
man_admix_data <- cbind(man_fam$V2, manning_tbl)

### rename columns
colnames(man_admix_data) <- c("ID", "q", "p")

### add info
man_admix_data <- left_join(man_admix_data, man_info)
#write out
write.csv(man_admix_data, "man_admix_data.csv")


# check admixture columns
# V1 = RBSA high, V2 = RNSA high, keep columns

### order data by column value so they're easily visualized
manning_tbl <- manning_tbl[order(manning_tbl$V2,manning_tbl$V1),]
manning_tbl
mean(manning_tbl$V1)
### visualize
barplot(t(as.matrix(manning_tbl)), col=c("#FF4343", "#4343FF"),
        xlab="Individual", ylab="Ancestry", border=NA, space=0, xaxt='n', main = "Manning ADMIXTURE")

pdf(file = "manning_tassel_admixture_plot_k2.pdf")
barplot(t(as.matrix(manning_tbl)), col=c("#FF4343", "#4343FF"),
        xlab="Individual", ylab="Ancestry", border=NA, space=0, xaxt='n')
dev.off()



### read in data
caor_tbl=read.table("caor_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.2.Q")
caor_tbl
caor_fam <- read.table("caor_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.fam")
caor_info <- read.table("../sample_info_caor.txt", header=TRUE)

### add ID column from the .fam file
caor_admix_data <- cbind(caor_fam$V2, caor_tbl)

### rename columns
colnames(caor_admix_data) <- c("ID", "p", "q")

### add info
caor_admix_data <- left_join(caor_admix_data, caor_info)
#write out
write.csv(caor_admix_data, "caor_admix_data.csv")

# check admixture columns
# V1 = RNSA high, V2 = RBSA high, switch columns
caor_col.order <- c("V2","V1")
caor_tbl <- caor_tbl[,caor_col.order]

### order data by column value so they're easily visualized
caor_tbl <- caor_tbl[order(caor_tbl$V1,caor_tbl$V2),]
caor_tbl

### visualize
barplot(t(as.matrix(caor_tbl)), col=c("#FF4343", "#4343FF"),
        xlab="Individual", ylab="Ancestry", border=NA, space=0, xaxt='n', main = "CA/OR ADMIXTURE")

barplot(t(as.matrix(caor_tbl)), col=c("red", "blue"),
        border=NA, space=0, xaxt='n')

pdf(file = "caor_tassel_admixture_plot_k2.pdf")
barplot(t(as.matrix(caor_tbl)), col=c("#FF4343", "#4343FF"),
        xlab="Individual", ylab="Ancestry", border=NA, space=0, xaxt='n')
dev.off()


# plot all three together

pdf(file="admixture_three_plots.pdf")
par(mar=c(1,1,1,1), mfrow=c(3,1))
barplot(t(as.matrix(nbc_tbl)), col=c( "#FF4343", "#4343FF"),
        xlab="Individual", ylab="Ancestry", border=NA, space=0, xaxt='n')
barplot(t(as.matrix(manning_tbl)), col=c("#FF4343", "#4343FF"),
        xlab="Individual", ylab="Ancestry", border=NA, space=0, xaxt='n')
barplot(t(as.matrix(caor_tbl)), col=c("#FF4343", "#4343FF"),
        xlab="Individual", ylab="Ancestry", border=NA, space=0, xaxt='n')





###########

### get parental groups for FSTs
nbc_rb <- filter(nbc_admix_data, q > 0.95)
print(nbc_rb$ID)
# [1] "KE21S01:HWI-ST765:6:ss_plate_6:B04" "KE21S02:HWI-ST765:6:ss_plate_6:B05"
# [3] "KE22S01:HWI-ST765:6:ss_plate_6:B06" "KE22S02:HWI-ST765:6:ss_plate_6:B07"
# [5] "KE22S03:HWI-ST765:6:ss_plate_6:B08" "KE22S04:HWI-ST765:6:ss_plate_6:B09"
# [7] "KE24S01:HWI-ST765:6:ss_plate_6:A11" "KE25S01:HWI-ST765:6:ss_plate_6:A03"
# [9] "KE25S03:HWI-ST765:6:ss_plate_6:A07" "KE27S02:HWI-ST765:6:ss_plate_6:A05"
# [11] "KE27S04:HWI-ST765:6:ss_plate_6:A06" "KE29S01:HWI-ST765:6:ss_plate_6:A08"
# [13] "KE29S02:HWI-ST765:6:ss_plate_6:A09" "KE29S04:HWI-ST765:6:ss_plate_6:A10"
# [15] "KE29S05:HWI-ST765:6:ss_plate_6:A12" "KE29S06:HWI-ST765:6:ss_plate_6:B01"
# [17] "KE30S01:HWI-ST765:6:ss_plate_6:B02" "KE30S02:HWI-ST765:6:ss_plate_6:B03"

nbc_rn <- filter(nbc_admix_data, q < 0.05)
print(nbc_rn$ID)
# [1] "JF14S01:HWI-ST765:6:ss_plate_6:F06" "JF14S02:HWI-ST765:6:ss_plate_6:F07"
# [3] "JF15S01:HWI-ST765:6:ss_plate_6:F10" "JF16S03:HWI-ST765:6:ss_plate_6:G03"
# [5] "JF16S04:HWI-ST765:6:ss_plate_6:G04" "JF17S01:HWI-ST765:6:ss_plate_6:F08"
# [7] "JF17S02:HWI-ST765:6:ss_plate_6:F09" "KE18S01:HWI-ST765:6:ss_plate_6:E08"
# [9] "KE18S02:HWI-ST765:6:ss_plate_6:E09" "KE18S03:HWI-ST765:6:ss_plate_6:E10"
# [11] "KE18S04:HWI-ST765:6:ss_plate_6:E11" "KE19S01:HWI-ST765:6:ss_plate_6:E12"
# [13] "KE19S02:HWI-ST765:6:ss_plate_6:F01" "KE31S02:HWI-ST765:6:ss_plate_6:F03"
# [15] "KF03S01:HWI-ST765:6:ss_plate_6:G06" "KF04S01:HWI-ST765:6:ss_plate_6:G07"
# [17] "KF04S07:HWI-ST765:6:ss_plate_6:G09" "KF05S01:HWI-ST765:6:ss_plate_6:G08"
# [19] "KF11S01:HWI-ST765:6:ss_plate_6:G10" "KF11S02:HWI-ST765:6:ss_plate_6:G11"
# [21] "KF12S04:HWI-ST765:6:ss_plate_6:F05" "KF12S05:HWI-ST765:6:ss_plate_6:F02"
# [23] "KF12S06:HWI-ST765:6:ss_plate_6:C12" "KF14S01:HWI-ST765:6:ss_plate_6:D01"
# [25] "KF15S05:HWI-ST765:6:ss_plate_6:D05" "KF16S03:HWI-ST765:6:ss_plate_6:D07"
# [27] "KF19S02:HWI-ST765:6:ss_plate_6:C01" "KF19S03:HWI-ST765:6:ss_plate_6:B10"
# [29] "KF19S05:HWI-ST765:6:ss_plate_6:B11" "KF20S01:HWI-ST765:6:ss_plate_6:B12"
# [31] "KF20S02:HWI-ST765:6:ss_plate_6:C02" "KF20S03:HWI-ST765:6:ss_plate_6:C03"
# [33] "KF20S05:HWI-ST765:6:ss_plate_6:D12" "KF21S02:HWI-ST765:6:ss_plate_6:E03"
# [35] "KF21S03:HWI-ST765:6:ss_plate_6:E04" "KF21S04:HWI-ST765:6:ss_plate_6:E05"

nbc_hyb <- filter(nbc_admix_data, q > 0.05 & q < 0.95)
nbc_hyb



man_rb <- filter(man_admix_data, q > 0.95)
print(man_rb$ID)
# [1] "QF16L01:HTHMNDRXX:1:plate3:G10" "RD30L01:H5J22BBXY:6:plate1:F07"
# [3] "RE05L01:H5J22BBXY:6:plate1:H02" "RE11L01:H5J22BBXY:6:plate1:E02"
# [5] "RE11L02:H5J22BBXY:6:plate1:G05" "RE16L01:H5J22BBXY:6:plate1:B11"
# [7] "RE17L01:H5J22BBXY:6:plate1:A06" "RE21L01:HTHMNDRXX:1:plate3:A11"
# [9] "RE23L01:H5J22BBXY:6:plate1:H06" "RE27L01:H5J22BBXY:6:plate1:H03"
# [11] "RE30L01:H5J22BBXY:6:plate1:D03" "RF18L01:H5J22BBXY:6:plate1:G11"
# [13] "RF22L01:HTHMNDRXX:1:plate3:C03" "RG09L02:H5J22BBXY:6:plate1:H04"
# [15] "UF09L01:HTHMNDRXX:1:plate3:B01" "UF10L01:HTHMNDRXX:1:plate3:H07"
# [17] "UF10L02:HTHMNDRXX:1:plate3:G06"

man_rn <- filter(man_admix_data, q < 0.05)
print(man_rn$ID)
# [1] "QF15L01:H5J22BBXY:6:plate1:A11" "QF17L01:H5J22BBXY:6:plate1:B04"
# [3] "QF19L01:H5J22BBXY:6:plate1:A01" "QF19L02:H5J22BBXY:6:plate1:A04"
# [5] "RD25L01:H5J22BBXY:6:plate1:A02" "RE02L01:H5J22BBXY:6:plate1:A12"
# [7] "RE08L01:H5J22BBXY:6:plate1:D08" "RF06L01:H5J22BBXY:6:plate1:C07"
# [9] "RF10L01:H5J22BBXY:6:plate1:B01" "RF19L01:H5J22BBXY:6:plate1:E04"
# [11] "RF21L01:H5J22BBXY:6:plate1:D01" "RF21L02:H5J22BBXY:6:plate1:C12"
# [13] "RF21L03:H5J22BBXY:6:plate1:C11" "RF21L04:H5J22BBXY:6:plate1:H07"
# [15] "RF24L01:H5J22BBXY:6:plate1:E09" "RG03L01:H5J22BBXY:6:plate1:G02"
# [17] "RG03L02:HTHMNDRXX:1:plate3:C05" "RG06L01:H5J22BBXY:6:plate1:G03"
# [19] "RG09L01:H5J22BBXY:6:plate1:B12" "UF08L01:HTHMNDRXX:1:plate3:G12"
# [21] "UF11L01:HTHMNDRXX:1:plate3:F06"

man_hyb <- filter(man_admix_data, q > 0.05 & q < 0.95)
man_hyb



caor_rb <- filter(caor_admix_data, q > 0.95)
print(caor_rb$ID)
# [1] "64453:C4Y4MACXX:7:250354918" "64581:C4Y4MACXX:7:250354939"
# [3] "B_205:C55VHACXX:1:250382565" "B_206:C55VHACXX:1:250382566"
# [5] "B_214:C55VHACXX:1:250382574" "B_217:C55VHACXX:1:250382577"
# [7] "B_223:C55VHACXX:1:250382583" "B_227:C55VHACXX:1:250382587"
# [9] "B_229:C55VHACXX:1:250382589" "B_233:C55VHACXX:1:250382593"
# [11] "B_236:C55VHACXX:1:250382596" "B_241:C55VHACXX:1:250382601"
# [13] "B_242:C55VHACXX:1:250382602" "B_266:C55VHACXX:1:250382626"
# [15] "B_271:C55VHACXX:1:250382631" "B_275:C55VHACXX:1:250382635"
# [17] "B_276:C55VHACXX:1:250382636" "B_277:C55VHACXX:1:250382637"
# [19] "B_278:C55VHACXX:1:250382638" "B_279:C55VHACXX:1:250382639"
# [21] "B_280:C55VHACXX:1:250382640" "B_281:C55VHACXX:1:250382641"
# [23] "B_283:C55VHACXX:1:250382643" "B_284:C55VHACXX:1:250382644"
# [25] "B_285:C55VHACXX:1:250382645" "B_307:C4Y4MACXX:7:250354886"
# [27] "B_308:C4Y4MACXX:7:250354887" "B_312:C4Y4MACXX:7:250354891"
# [29] "B_313:C4Y4MACXX:7:250354892" "B_314:C4Y4MACXX:7:250354893"
# [31] "B_315:C4Y4MACXX:7:250354894" "B_316:C4Y4MACXX:7:250354895"
# [33] "B_317:C4Y4MACXX:7:250354896" "B_320:C4Y4MACXX:7:250354899"
# [35] "B_321:C4Y4MACXX:7:250354900" "B_322:C4Y4MACXX:7:250354901"
# [37] "B_323:C4Y4MACXX:7:250354902" "B_325:C4Y4MACXX:7:250354904"
# [39] "B_326:C4Y4MACXX:7:250354905" "B_328:C4Y4MACXX:7:250354907"
# [41] "B_329:C4Y4MACXX:7:250354908" "B_332:C4Y4MACXX:7:250354911"


caor_rn <- filter(caor_admix_data, q < 0.05)
print(caor_rn$ID)
# [1] "B_243:C55VHACXX:1:250382603" "B_244:C55VHACXX:1:250382604"
# [3] "B_245:C55VHACXX:1:250382605" "B_246:C55VHACXX:1:250382606"
# [5] "B_248:C55VHACXX:1:250382608" "B_249:C55VHACXX:1:250382609"
# [7] "B_251:C55VHACXX:1:250382611" "B_252:C55VHACXX:1:250382612"
# [9] "B_253:C55VHACXX:1:250382613" "B_255:C55VHACXX:1:250382615"
# [11] "B_256:C55VHACXX:1:250382616" "B_258:C55VHACXX:1:250382618"
# [13] "B_260:C55VHACXX:1:250382620" "B_261:C55VHACXX:1:250382621"
# [15] "B_262:C55VHACXX:1:250382622" "B_263:C55VHACXX:1:250382623"
# [17] "B_268:C55VHACXX:1:250382628" "B_336:C4Y4MACXX:7:250354915"

caor_hyb <- filter(caor_admix_data, q > 0.05 & q < 0.95)
caor_hyb

### want to make a chi-square of hyb:parentals. first make dataframe, then run chisq
nbc_hyb2parental <- c((nrow(nbc_hyb)),((nrow(nbc_rb)+nrow(nbc_rn))))
man_hyb2parental <- c((nrow(man_hyb)),((nrow(man_rb)+nrow(man_rn))))
caor_hyb2parental <- c((nrow(caor_hyb)),((nrow(caor_rb)+nrow(caor_rn))))
hyb_prop_x2_df <- data.frame(nbc_hyb2parental, man_hyb2parental, caor_hyb2parental)
colnames(hyb_prop_x2_df) <- c("hybrid","parental")
chisq.test(hyb_prop_x2_df)

### okay yes significant effect of transect on population's hybrid proportion, let's look pairwise
nbc_man_x2_df <- data.frame(nbc_hyb2parental, man_hyb2parental)
chisq.test(nbc_man_x2_df)
nbc_caor_x2_df <- data.frame(nbc_hyb2parental, caor_hyb2parental)
chisq.test(nbc_caor_x2_df)
man_caor_X2_df <- data.frame(man_hyb2parental, caor_hyb2parental)
chisq.test(man_caor_X2_df)
