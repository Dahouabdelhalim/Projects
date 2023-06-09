library(tidyverse)

dat <- read.table("submission_files/qpcr_data.txt", header=T)
meta <- read.table("resubmission_files/sample_metadata_v2.txt", header=T)

## merge with metadata
dat.met <- inner_join(dat, meta, by="SampleID")

## multiply all Amount_SYBR_Copies * 80ul (the volume eluted per hindgut)
## this is Amount_SYBR_Copies per reaction, and 1 ul gDNA was added to each reaction
## this gives you a per-gut, rather than per-reaction, value
dat.met$Amount_SYBR_Copies_gut <- dat.met$Amount_SYBR_Copies * 80

## note need to adjust for Y-33-MG for which only ~1/3 was extracted (see qubit, p. 130)
## i.e. not 1 gut / 80 ul for this sample but rather 1/3 gut / 80 ul
dat.met$Amount_SYBR_Copies[dat.met$SampleID=="Y-33_MG"] <- 3 * dat.met$Amount_SYBR_Copies[dat.met$SampleID=="Y-33_MG"]

dat_TJH <- filter(dat.met, Project == "TJH_BAM")

dat_TJH_mean <- dat_TJH %>% filter(Caste != "NA") %>%     #remove extr blanks, qpcr-NTC
  group_by(SampleID, Age_days, Colony, Species, Caste, BeeID, GutType) %>%
  dplyr::summarize(mean_copies_gut = mean(Amount_SYBR_Copies_gut))  #average technical replicates

dat_TJH_mean_noQ <- filter(dat_TJH_mean, Caste != "queen")




######### stats/modeling

## code for one set of points (one colony, one gut type)
temp <- dat_TJH_mean_noQ %>% filter(Colony=="B", GutType=="MG")

# use SSlogis() function to find initial parameter estimates
# to be fed to the nls() function, which finds the best-fit logistic equation

time <- temp$Age_days
population <- log10(temp$mean_copies_gut)
plot(time, population, las=1, pch=16)

# Asym = K, xmid = x value at inflection point, scal = scaling parameter for x
logisticModelSS <- nls(population~SSlogis(time, Asym, xmid, scal))
summary(logisticModelSS)
coef(logisticModelSS)

## repeat for all colony X gut type combo's, and plot in ggplot
dat_TJH_mean_noQ$log10_mean_copies_gut <- log10(dat_TJH_mean_noQ$mean_copies_gut)

modellist <- list()
colonies <- unique(dat_TJH_mean_noQ$Colony)
guttypes <- unique(dat_TJH_mean_noQ$GutType)
counter <- 0

df_len <- length(colonies) * length(guttypes)
col_gut_factors <- data.frame(numeric(df_len),numeric(df_len),numeric(df_len)); names(col_gut_factors) <- c("modelnum","colony","guttype")

for (i in 1:length(colonies)) {
  tmp <- filter(dat_TJH_mean_noQ, Colony==colonies[i])
  
  for (j in 1:length(guttypes)) {
    tmp2 <- filter(tmp, GutType==guttypes[j])
    logisticModelSS <- nls(data=tmp2, log10_mean_copies_gut ~ SSlogis(Age_days, Asym, xmid, scal))
    
    counter <- counter + 1
    print(counter)
    modellist[[counter]] <- logisticModelSS
    
    col_gut_factors$modelnum[counter] <- counter
    col_gut_factors$colony[counter] <- colonies[i]
    col_gut_factors$guttype[counter] <- guttypes[j]
    
  }
}


#add colors based on colony
col_gut_factors$color[col_gut_factors$colony=="B"] <- "#9ecae1" 
col_gut_factors$color[col_gut_factors$colony=="Y"] <- "#ffeda0"
col_gut_factors$color[col_gut_factors$colony=="W"] <- "#f0f0f0"


## facet by gut type, show all 3 colonies together

#facet labels for gut type
guttype.lab <- c(
  HG = "Hindgut",
  MG = "Midgut"
)

qpcrplot <- ggplot(dat_TJH_mean_noQ, aes(x=Age_days, y=log10_mean_copies_gut,fill=Colony)) +
  geom_point(size=4, alpha=0.35, shape=21) +
  theme_bw() +
  ylab("log10(16S rRNA gene copies per gut)") + xlab("Age (days)") +
  scale_fill_manual(values=c("Y"="#ffeda0",
                              "B"="#9ecae1",
                              "W"="#f0f0f0")) +
  facet_wrap(~GutType,nrow=2, ncol=1, scales="fixed",
             labeller = labeller(GutType = guttype.lab)) +
  theme(axis.text = element_text(size=10)) +
  theme(strip.text = element_text(size=10))


grp1 <- filter(dat_TJH_mean_noQ, Colony==col_gut_factors$colony[1], GutType==col_gut_factors$guttype[1])
grp2 <- filter(dat_TJH_mean_noQ, Colony==col_gut_factors$colony[2], GutType==col_gut_factors$guttype[2])
grp3 <- filter(dat_TJH_mean_noQ, Colony==col_gut_factors$colony[3], GutType==col_gut_factors$guttype[3])
grp4 <- filter(dat_TJH_mean_noQ, Colony==col_gut_factors$colony[4], GutType==col_gut_factors$guttype[4])
grp5 <- filter(dat_TJH_mean_noQ, Colony==col_gut_factors$colony[5], GutType==col_gut_factors$guttype[5])
grp6 <- filter(dat_TJH_mean_noQ, Colony==col_gut_factors$colony[6], GutType==col_gut_factors$guttype[6])

qpcrplot_w_logcurves <- qpcrplot +
  geom_line(data=grp1,
            aes(y=predict(modellist[[1]])),
            color=col_gut_factors$color[1],
            size=2) +
  geom_line(data=grp2,
            aes(y=predict(modellist[[2]])),
            color=col_gut_factors$color[2],
            size=2) +
  geom_line(data=grp3,
            aes(y=predict(modellist[[3]])),
            color=col_gut_factors$color[3],
            size=2) +
  geom_line(data=grp4,
            aes(y=predict(modellist[[4]])),
            color=col_gut_factors$color[4],
            size=2) +
  geom_line(data=grp5,
            aes(y=predict(modellist[[5]])),
            color=col_gut_factors$color[5],
            size=2) +
  geom_line(data=grp6,
            aes(y=predict(modellist[[6]])),
            color=col_gut_factors$color[6],
            size=2)

