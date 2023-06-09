
# ******************************************************************************
# ******************************************************************************
# ======================== general =============================================
#
# A load data
# B correlation among variables
# C effect size dimorphism
#
# ======================== notes =============================================
#
# ======================== A load data ==============================================
# data file to load
name.t <- "all.traits.clean.size.corrected.outlier.rm.csv"

# load data files from data.clean
d <- read.csv(paste(cd.path,name.t, sep = ""), stringsAsFactors = FALSE)
str(d)


# *** get sample sizes for each character per island per sex ***
# all samples included (max sample size)
# make name for the sample size files
sample.size.max <- paste(re.path,'2a.sample.size.max.csv',sep = "")
str(d) # 591 obs
# save overview to file
write.csv(table(d$island, d$sex), sample.size.max, row.names = TRUE)

# remove all rows with missing values (min sample size)
# make name for the sample size files
sample.size.min <- paste(re.path,'2a.sample.size.min.csv', sep = "")
d.min <- na.omit(d) # exclude all rows with an NA
str(d.min) # 561 obs
# save overview to file
write.csv(table(d.min$island, d.min$sex), sample.size.min, row.names = TRUE)

# ======================== B correlation among variables ======================

# get covariance matrix to check correlation among variables
corr.d <- cor(d[,5:15],use = "pairwise.complete.obs")
corr.d.res <- cor(d[,c(5,16:ncol(d))],use = "pairwise.complete.obs")

# plot correlations
pdf(paste(fi.path,'2a.VarCovar.pdf', sep = ""), width = 5, height = 6)
corrplot(corr.d, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
dev.off()
# residuals
pdf(paste(fi.path,'2a.VarCovarVarCovarResidual.pdf', sep = ""), width = 5, height = 6)
corrplot(corr.d.res, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
dev.off()

# save var-covar table
write.csv(corr.d, paste(re.path,'2a.VarCovar.csv', sep = ""), row.names = TRUE)
# save var-covar table
write.csv(corr.d.res, paste(re.path,'2a.VarCovarResidual.csv', sep = ""), 
          row.names = TRUE)

# ======================== C  effect size dimorphism  =================
# estimate the effect size of the dimorphism for one and two species islands

mean.sexes <-aggregate(d[,5:15],list(d$island, d$sex, d$presence), 
                       mean, na.rm=TRUE)
colnames(mean.sexes)[1:3] <- c("island", "sex", "presence")
head(mean.sexes)
mean.f <- mean.sexes[mean.sexes$sex == "f",]
mean.m <- mean.sexes[mean.sexes$sex == "m",]

# manually check if islands align for males and females
cbind(mean.f$island, mean.m$island)

effect.size <- cbind.data.frame(mean.f[, c(1,3)], mean.m[,4:ncol(mean.m)] - 
                                  mean.f[,4:ncol(mean.f)] )
# mean effect size
mean.effect.size <- aggregate(effect.size[,3:ncol(effect.size)], 
                              list(effect.size$presence), mean)
mean.effect.size.1 <- cbind.data.frame(rep("mean", 2), mean.effect.size)

# sd of effect size
sd.effect.size <- aggregate(effect.size[,3:ncol(effect.size)], 
                            list(effect.size$presence), sd)
sd.effect.size.1 <- cbind.data.frame(rep("sd", 2), sd.effect.size)

# save both combined
colnames(mean.effect.size.1)[1:2] <- colnames(sd.effect.size.1)[1:2] <- 
  c("descr.stats", "presence")

descr.effect.size <- rbind.data.frame(mean.effect.size.1, sd.effect.size.1)

write.csv(descr.effect.size, paste(re.path,'2a.effect.size.dimorph.csv', sep = ""), row.names = TRUE)

# ******************************************************************************
# ******************************************************************************










