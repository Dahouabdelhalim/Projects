library(Biobase)
library(MSnbase)
library(limma)

f <- "female_list.csv"
e <- c(2:11)
x <- readMSnSet2(f, e, fnames=1)

Sample<-gsub("(w*).(d)", "2.1", colnames(exprs(x)))
# for pattern "generic word [.] digit", replace with ".2 or .1", for exprs(x)

tmp <- data.frame(do.call(rbind, strsplit(Sample, "[.]")))
#string split at "."

names(tmp) <- c("Sample","Replicate")

sampleNames(x)<-paste(tmp$Sample, tmp$Replicate, sep = "")
rownames(tmp) <- paste(tmp$Sample, tmp$Replicate, sep = "")
pData(x) <- tmp

x <- log(x,2)
x <- normalise(x,"diff.median")
boxplot(exprs(x))
write.csv(exprs(x),"norm_female_proteins.csv")

#########
#sim mated to virgin comparison
#need to add [.] before 1 or 2 in column headers
f <- "norm_female_proteins.csv"
e <- c(8:11)
x <- readMSnSet2(f, e, fnames=1)

Sample<-gsub("(w*).(d)", "2.1", colnames(exprs(x)))
tmp <- data.frame(do.call(rbind, strsplit(Sample, "[.]")))
names(tmp) <- c("Sample","Replicate")

sampleNames(x)<-paste(tmp$Sample, tmp$Replicate, sep = "")
rownames(tmp) <- paste(tmp$Sample, tmp$Replicate, sep = "")
pData(x) <- tmp

boxplot(exprs(x))

mat.design <- model.matrix(~ 0 + pData(x)$Sample)
colnames(mat.design) <- c('SS','SV')

contrast.matrix <- makeContrasts(SS-SV,
                                 levels = mat.design)

fit <- lmFit(exprs(x),design=mat.design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
tt <- topTable(fit3, adjust.method = "BH", n = Inf,coef=1)
write.csv(tt,"sssv_female.csv")

#mau mated to virgin comparison
f <- "norm_female_proteins.csv"
e <- c(2:5)
x <- readMSnSet2(f, e, fnames=1)

Sample<-gsub("(w*).(d)", "2.1", colnames(exprs(x)))
tmp <- data.frame(do.call(rbind, strsplit(Sample, "[.]")))
names(tmp) <- c("Sample","Replicate")

sampleNames(x)<-paste(tmp$Sample, tmp$Replicate, sep = "")
rownames(tmp) <- paste(tmp$Sample, tmp$Replicate, sep = "")
pData(x) <- tmp

boxplot(exprs(x))

mat.design <- model.matrix(~ 0 + pData(x)$Sample)
colnames(mat.design) <- c('MM','MV')

contrast.matrix <- makeContrasts(MM-MV,
                                 levels = mat.design)

fit <- lmFit(exprs(x),design=mat.design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
tt <- topTable(fit3, adjust.method = "BH", n = Inf,coef=1)
write.csv(tt,"mmmv_female.csv")

#########
#virgins
#rearrange columns of norm_female_proteins to: SV, MV, SM, SS, MM
#new arrangement is "norm_female_proteins_reordered"

f <- "norm_female_proteins_reordered.csv"
e <- c(2:5)
x <- readMSnSet2(f, e, fnames=1)

Sample<-gsub("(w*).(d)", "2.1", colnames(exprs(x)))
tmp <- data.frame(do.call(rbind, strsplit(Sample, "[.]")))
names(tmp) <- c("Sample","Replicate")

sampleNames(x)<-paste(tmp$Sample, tmp$Replicate, sep = "")
rownames(tmp) <- paste(tmp$Sample, tmp$Replicate, sep = "")
pData(x) <- tmp

boxplot(exprs(x))

mat.design <- model.matrix(~ 0 + pData(x)$Sample)
colnames(mat.design) <- c('SV','MV')

contrast.matrix <- makeContrasts(SV-MV,
                                 levels = mat.design)

fit <- lmFit(exprs(x),design=mat.design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
tt <- topTable(fit3, adjust.method = "BH", n = Inf,coef=1)
write.csv(tt,"svmv_female.csv")

#########
#mated

f <- "norm_female_proteins_reordered.csv"
e <- c(8:11)
x <- readMSnSet2(f, e, fnames=1)

Sample<-gsub("(w*).(d)", "2.1", colnames(exprs(x)))
tmp <- data.frame(do.call(rbind, strsplit(Sample, "[.]")))
names(tmp) <- c("Sample","Replicate")

sampleNames(x)<-paste(tmp$Sample, tmp$Replicate, sep = "")
rownames(tmp) <- paste(tmp$Sample, tmp$Replicate, sep = "")
pData(x) <- tmp

boxplot(exprs(x))

mat.design <- model.matrix(~ 0 + pData(x)$Sample)
colnames(mat.design) <- c('SS','MM')

contrast.matrix <- makeContrasts(SS-MM,
                                 levels = mat.design)

fit <- lmFit(exprs(x),design=mat.design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
tt <- topTable(fit3, adjust.method = "BH", n = Inf,coef=1)
write.csv(tt,"ssmm_female.csv")

#########
#het vs ss
#rearrange columns of norm_female_proteins to: SV, MV, SS, SM, MM
#new arrangement is "norm_female_proteins_reordered2"

f <- "norm_female_proteins_reordered2.csv"
e <- c(6:9)
x <- readMSnSet2(f, e, fnames=1)

Sample<-gsub("(w*).(d)", "2.1", colnames(exprs(x)))
tmp <- data.frame(do.call(rbind, strsplit(Sample, "[.]")))
names(tmp) <- c("Sample","Replicate")

sampleNames(x)<-paste(tmp$Sample, tmp$Replicate, sep = "")
rownames(tmp) <- paste(tmp$Sample, tmp$Replicate, sep = "")
pData(x) <- tmp

boxplot(exprs(x))

mat.design <- model.matrix(~ 0 + pData(x)$Sample)
colnames(mat.design) <- c('SS','SM')

contrast.matrix <- makeContrasts(SS-SM,
                                 levels = mat.design)

fit <- lmFit(exprs(x),design=mat.design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
tt <- topTable(fit3, adjust.method = "BH", n = Inf,coef=1)
write.csv(tt,"sssm_female.csv")

#########
#het vs mm

f <- "norm_female_proteins_reordered2.csv"
e <- c(8:11)
x <- readMSnSet2(f, e, fnames=1)

Sample<-gsub("(w*).(d)", "2.1", colnames(exprs(x)))
tmp <- data.frame(do.call(rbind, strsplit(Sample, "[.]")))
names(tmp) <- c("Sample","Replicate")

sampleNames(x)<-paste(tmp$Sample, tmp$Replicate, sep = "")
rownames(tmp) <- paste(tmp$Sample, tmp$Replicate, sep = "")
pData(x) <- tmp

boxplot(exprs(x))

mat.design <- model.matrix(~ 0 + pData(x)$Sample)
colnames(mat.design) <- c('MM','SM')

contrast.matrix <- makeContrasts(MM-SM,
                                 levels = mat.design)

fit <- lmFit(exprs(x),design=mat.design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
tt <- topTable(fit3, adjust.method = "BH", n = Inf,coef=1)
write.csv(tt,"mmsm_female.csv")

#########
#het vs sim virgin
#rearrange columns of norm_female_proteins to: SV, SM, MV, SS, MM
#new arrangement is "norm_female_proteins_reordered3"

f <- "norm_female_proteins_reordered3.csv"
e <- c(2:5)
x <- readMSnSet2(f, e, fnames=1)

Sample<-gsub("(w*).(d)", "2.1", colnames(exprs(x)))
tmp <- data.frame(do.call(rbind, strsplit(Sample, "[.]")))
names(tmp) <- c("Sample","Replicate")

sampleNames(x)<-paste(tmp$Sample, tmp$Replicate, sep = "")
rownames(tmp) <- paste(tmp$Sample, tmp$Replicate, sep = "")
pData(x) <- tmp

boxplot(exprs(x))

mat.design <- model.matrix(~ 0 + pData(x)$Sample)
colnames(mat.design) <- c('SM','SV')

contrast.matrix <- makeContrasts(SM-SV,
                                 levels = mat.design)

fit <- lmFit(exprs(x),design=mat.design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
tt <- topTable(fit3, adjust.method = "BH", n = Inf,coef=1)
write.csv(tt,"smsv_female.csv")