# run stats on glial processes, starting with exact tests of goodness-of-fit for proportions of processes across the IPL
library(corrplot);
library("dplyr")
library(rstatix)
 
base <-  "C:/..."   #change to appropriate file location
fname <- paste(base, "Figure 1-source data 1.csv",sep="") 

# load/generate contingency table of process counts in each sublayer for each cell at all ages
fig1D_summaryData <- read.csv(fname)[,21:40]

#pool cells from same experiment, p8/9
e1pooled <- colSums(fig1D_summaryData[1:3,1:5]);
e2pooled <- colSums(fig1D_summaryData[4:10,1:5]);
e3pooled <- colSums(fig1D_summaryData[11:13,1:5]);
e4pooled <- colSums(fig1D_summaryData[14:19,1:5]);
e5pooled <- colSums(fig1D_summaryData[20:26,1:5]); 
e6pooled <- colSums(fig1D_summaryData[27:32,1:5]);

#pool experiments into one large table
allp8p9 <- rbind(e1pooled, e2pooled, e3pooled, e4pooled,e5pooled,e6pooled);

#pool cells for P10/11
e1pooled <- colSums(fig1D_summaryData[1:3,6:10]);
e2pooled <- colSums(fig1D_summaryData[4:6,6:10]);
e3pooled <- colSums(fig1D_summaryData[7:9,6:10]);
e4pooled <- colSums(fig1D_summaryData[10:18,6:10]);
e5pooled <- colSums(fig1D_summaryData[19:24,6:10]); 

allp10p11 <- rbind(e1pooled, e2pooled, e3pooled, e4pooled,e5pooled);

#now P12-17
e1pooled <- colSums(fig1D_summaryData[1:4,11:15]);
e2pooled <- colSums(fig1D_summaryData[5:9,11:15]);
e3pooled <- colSums(fig1D_summaryData[10:12,11:15]);
e4pooled <- colSums(fig1D_summaryData[13:18,11:15]);
e5pooled <- colSums(fig1D_summaryData[19:24,11:15]);  
e6pooled <- colSums(fig1D_summaryData[25:30,11:15]); 

allp12p17 <- rbind(e1pooled, e2pooled, e3pooled, e4pooled,e5pooled,e6pooled);

#now P20+
e1pooled <- colSums(fig1D_summaryData[1:6,16:20]);
e2pooled <- colSums(fig1D_summaryData[7:10,16:20]);
e3pooled <- colSums(fig1D_summaryData[11:14,16:20]);
e4pooled <- colSums(fig1D_summaryData[15:16,16:20]);
e5pooled <- colSums(fig1D_summaryData[17:18,16:20]);  

allp20 <- rbind(e1pooled, e2pooled, e3pooled, e4pooled,e5pooled);

#goodness-of-fit tests for process counts across IPL at each age (Fig 1D)
prob <- c(0.2,0.2,0.2,0.2,0.2);
p8p9_gof <- chisq.test(colSums(allp8p9),p=prob); #chi2 gof test for young age
p10p11_gof <- chisq.test(colSums(allp10p11),p=prob); #chi2 gof test for p10/11 age
p12p17_gof <- chisq.test(colSums(allp12p17),p=prob); #chi2 gof test for p12-17 age
p20_gof <- chisq.test(colSums(allp20),p=prob); #chi2 gof test for p12-17 age
p_corr <- p.adjust(c(p8p9_gof$p.value,p10p11_gof$p.value,p12p17_gof$p.value,p20_gof$p.value), "BH")

#now test for independence between ages (Fig 1D)
p8_9pool = colSums(allp8p9);
p10_11pool = colSums(allp10p11);
p12_17pool = colSums(allp12p17);
p20pool = colSums(allp20);
allAges = rbind(p8_9pool,p10_11pool,p12_17pool,p20pool);

age_chi2 <- chisq.test(allAges);
corrplot(age_chi2$residuals, is.cor = FALSE);
contrib <- 100*age_chi2$residuals^2/age_chi2$statistic;
round(contrib, 3);
corrplot(contrib, is.cor = FALSE)

library(chisq.posthoc.test)
chisq.posthoc.test(allAges, method = "BH", round = 6);

#directly compare proportion in each sublayer to same sublayer across different ages using chisq test (Fig 1D)
s1all <- cbind(allAges[,1],rowSums(allAges[,2:5]));
s2all<- cbind(allAges[,2],rowSums(allAges[,c(1,3:5)]));
s3all<-cbind(allAges[,3],rowSums(allAges[,c(1:2,4:5)]));
s4all<-cbind(allAges[,4],rowSums(allAges[,c(1:3,5)]));
s5all<-cbind(allAges[,5],rowSums(allAges[,c(1:4)]));

s1xsq<- chisq.test(s1all)
s2xsq<- chisq.test(s2all)
s3xsq<- chisq.test(s3all)
s4xsq<- chisq.test(s4all)
s5xsq<- chisq.test(s5all)
p_corr <- p.adjust(c(s1xsq$p.value,s2xsq$p.value,s3xsq$p.value,s4xsq$p.value,s5xsq$p.value), "BH")

# load contingency table of process counts in each sublayer for each motility category  at all ages
fname <- paste(base, "Figure 1-source data 1.csv",sep="") 
fig1E_summaryData <- read.csv(fname)[,42:61]

View(fig1E_summaryData);

#do chi2 test for independence between sublayer at each age 
chisq.test(fig1E_summaryData[1:7,])
chisq.posthoc.test(fig1E_summaryData[1:7,], method = "BH", round = 6);

#overall tests for independent proportions of categories in each sublayer for each age group (Fig. 1E)
allCats1 <- chisq.test(fig1E_summaryData[1:7,1:5]);
allCats2 <- chisq.test(fig1E_summaryData[1:7,6:10]);
allCats3 <- chisq.test(fig1E_summaryData[1:7,11:15]);
allCats4 <- chisq.test(fig1E_summaryData[1:7,16:20]);

# extensions, age1 
e1 <- chisq.test(rbind(fig1E_summaryData[1,1:5], colSums(fig1E_summaryData[2:7,1:5])));
# extensions, age2
e2 <- chisq.test(rbind(fig1E_summaryData[1,6:10], colSums(fig1E_summaryData[2:7,6:10])));
# extensions, age3
e3 <- chisq.test(rbind(fig1E_summaryData[1,11:15], colSums(fig1E_summaryData[2:7,11:15])));
# extensions, age4
e4 <- chisq.test(rbind(fig1E_summaryData[1,16:20], colSums(fig1E_summaryData[2:7,16:20])));

# new process, age1
np1 <- chisq.test(rbind(fig1E_summaryData[2,1:5], colSums(fig1E_summaryData[c(1,3:7),1:5])));
# new process, age2
np2 <- chisq.test(rbind(fig1E_summaryData[2,6:10], colSums(fig1E_summaryData[c(1,3:7),6:10])));
# new process, age3
np3 <- chisq.test(rbind(fig1E_summaryData[2,11:15], colSums(fig1E_summaryData[c(1,3:7),11:15])));
# new process, age4
np4 <- chisq.test(rbind(fig1E_summaryData[2,16:20], colSums(fig1E_summaryData[c(1,3:7),16:20])));

# retraction, age1
r1 <- chisq.test(rbind(fig1E_summaryData[3,1:5], colSums(fig1E_summaryData[c(1:2,4:7),1:5])));
# retraction, age2
r2 <- chisq.test(rbind(fig1E_summaryData[3,6:10], colSums(fig1E_summaryData[c(1:2,4:7),6:10])));
# retraction, age3
r3 <- chisq.test(rbind(fig1E_summaryData[3,11:15], colSums(fig1E_summaryData[c(1:2,4:7),11:15])));
# retraction, age4
r4 <- chisq.test(rbind(fig1E_summaryData[3,16:20], colSums(fig1E_summaryData[c(1:2,4:7),16:20])));

# lost process, age1
l1 <- chisq.test(rbind(fig1E_summaryData[4,1:5], colSums(fig1E_summaryData[c(1:3,5:7),1:5])));
# lost process, age2
l2 <- chisq.test(rbind(fig1E_summaryData[4,6:10], colSums(fig1E_summaryData[c(1:3,5:7),6:10])));
# lost process, age3
l3 <- chisq.test(rbind(fig1E_summaryData[4,11:15], colSums(fig1E_summaryData[c(1:3,5:7),11:15])));
# lost process, age4
l4 <- chisq.test(rbind(fig1E_summaryData[4,16:20], colSums(fig1E_summaryData[c(1:3,5:7),16:20])));

# e>r, age1
er1 <- chisq.test(rbind(fig1E_summaryData[5,1:5], colSums(fig1E_summaryData[c(1:4,6:7),1:5])));
# e>r, age2
er2 <- chisq.test(rbind(fig1E_summaryData[5,6:10], colSums(fig1E_summaryData[c(1:4,6:7),6:10])));
# e>r, age3
er3 <- chisq.test(rbind(fig1E_summaryData[5,11:15], colSums(fig1E_summaryData[c(1:4,6:7),11:15])));
# e>r, age4
er4 <- chisq.test(rbind(fig1E_summaryData[5,16:20], colSums(fig1E_summaryData[c(1:4,6:7),16:20])));

# r>e, age1
re1 <- chisq.test(rbind(fig1E_summaryData[6,1:5], colSums(fig1E_summaryData[c(1:5,7),1:5])));
# r>e, age2
re2 <- chisq.test(rbind(fig1E_summaryData[6,6:10], colSums(fig1E_summaryData[c(1:5,7),6:10])));
# r>e, age3
re3 <- chisq.test(rbind(fig1E_summaryData[6,11:15], colSums(fig1E_summaryData[c(1:5,7),11:15])));
# r>e, age4
re4 <- chisq.test(rbind(fig1E_summaryData[6,16:20], colSums(fig1E_summaryData[c(1:5,7),16:20])));

# stable, age1
s1 <- chisq.test(rbind(fig1E_summaryData[7,1:5], colSums(fig1E_summaryData[c(1:6),1:5])));
#stable, age2
s2 <- chisq.test(rbind(fig1E_summaryData[7,6:10], colSums(fig1E_summaryData[c(1:6),6:10])));
# stable, age3
s3 <- chisq.test(rbind(fig1E_summaryData[7,11:15], colSums(fig1E_summaryData[c(1:6),11:15])));
# stable, age4
s4 <- chisq.test(rbind(fig1E_summaryData[7,16:20], colSums(fig1E_summaryData[c(1:6),16:20])));

#compare proportion stable processes across ages (Fig 1F)
#run kruskal wallis on this data, with pairwise wilcoxon post-hoc tests
stableProps <- read.csv(fname)[,62:63]

colnames(stableProps)[1] <- 'stable_props'
colnames(stableProps)[2] <- 'group'
stableProps$group <- ordered(stableProps$group, levels = c("P8-9", "P10-11", "P12-13","P14", "P16", "P17", "P20","P23-24","P40+"))
test <- kruskal.test(stable_props ~ group, data = stableProps)
pairTest <- pairwise.wilcox.test(stableProps$stable_props, stableProps$group, p.adj = 'fdr')
stableProps %>%
  group_by(group) %>%
  get_summary_stats(stable_props, type = "mean_se")

# now do the same process for cytoskeleton manipulations (Fig. 1G)
fname <- paste(base, "Figure 1-source data 2.csv",sep="") 
cytoSkel <- read.csv(fname)[,2:5]

#compare overall proportions of motility categories using chisq tests for independent proportions
acsf_cytoD <- chisq.test(cytoSkel[1:7,1:2]) #ACSF vs. cytochalasin-D

acsf_cytoD_ext <- chisq.test(rbind(cytoSkel[1,1:2], colSums(cytoSkel[2:7,1:2])));
acsf_cytoD_np <- chisq.test(rbind(cytoSkel[2,1:2], colSums(cytoSkel[c(1,3:7),1:2])));
acsf_cytoD_retr <- chisq.test(rbind(cytoSkel[3,1:2], colSums(cytoSkel[c(1:2,4:7),1:2])));
acsf_cytoD_lp <- chisq.test(rbind(cytoSkel[4,1:2], colSums(cytoSkel[c(1:3,5:7),1:2])));
acsf_cytoD_er <- chisq.test(rbind(cytoSkel[5,1:2], colSums(cytoSkel[c(1:4,6:7),1:2])));
acsf_cytoD_re <- chisq.test(rbind(cytoSkel[6,1:2], colSums(cytoSkel[c(1:5,7),1:2])));
acsf_cytoD_stab <- chisq.test(rbind(cytoSkel[7,1:2], colSums(cytoSkel[c(1:6),1:2])));
p_corr <- p.adjust(c(acsf_cytoD_ext$p.value,acsf_cytoD_np$p.value,acsf_cytoD_retr$p.value,acsf_cytoD_lp$p.value,acsf_cytoD_er$p.value,acsf_cytoD_re$p.value,acsf_cytoD_stab$p.value), "BH")


acsf_noco <- chisq.test(cytoSkel[1:7,c(1,3)]) #ACSF vs. nocodazole

acsf_noco_ext <- chisq.test(rbind(cytoSkel[1,c(1,3)], colSums(cytoSkel[2:7,c(1,3)])));
acsf_noco_np <- chisq.test(rbind(cytoSkel[2,c(1,3)], colSums(cytoSkel[c(1,3:7),c(1,3)])));
acsf_noco_retr <- chisq.test(rbind(cytoSkel[3,c(1,3)], colSums(cytoSkel[c(1:2,4:7),c(1,3)])));
acsf_noco_lp <- chisq.test(rbind(cytoSkel[4,c(1,3)], colSums(cytoSkel[c(1:3,5:7),c(1,3)])));
acsf_noco_er <- chisq.test(rbind(cytoSkel[5,c(1,3)], colSums(cytoSkel[c(1:4,6:7),c(1,3)])));
acsf_noco_re <- chisq.test(rbind(cytoSkel[6,c(1,3)], colSums(cytoSkel[c(1:5,7),c(1,3)])));
acsf_noco_stab <- chisq.test(rbind(cytoSkel[7,c(1,3)], colSums(cytoSkel[c(1:6),c(1,3)])));
p_corr <- p.adjust(c(acsf_noco_ext$p.value,acsf_noco_np$p.value,acsf_noco_retr$p.value,acsf_noco_lp$p.value,acsf_noco_er$p.value,acsf_noco_re$p.value,acsf_noco_stab$p.value), "BH")

acsf_cyto_noco <- chisq.test(cytoSkel[1:7,c(1,4)]) #ACSF vs. cytochalasin-D + nocodazole

acsf_cyto_noco_ext <- chisq.test(rbind(cytoSkel[1,c(1,4)], colSums(cytoSkel[2:7,c(1,4)])));
acsf_cyto_noco_np <- chisq.test(rbind(cytoSkel[2,c(1,4)], colSums(cytoSkel[c(1,3:7),c(1,4)])));
acsf_cyto_noco_retr <- chisq.test(rbind(cytoSkel[3,c(1,4)], colSums(cytoSkel[c(1:2,4:7),c(1,4)])));
acsf_cyto_noco_lp <- chisq.test(rbind(cytoSkel[4,c(1,4)], colSums(cytoSkel[c(1:3,5:7),c(1,4)])));
acsf_cyto_noco_er <- chisq.test(rbind(cytoSkel[5,c(1,4)], colSums(cytoSkel[c(1:4,6:7),c(1,4)])));
acsf_cyto_noco_re <- chisq.test(rbind(cytoSkel[6,c(1,4)], colSums(cytoSkel[c(1:5,7),c(1,4)])));
acsf_cyto_noco_stab <- chisq.test(rbind(cytoSkel[7,c(1,4)], colSums(cytoSkel[c(1:6),c(1,4)])));
p_corr <- p.adjust(c(acsf_cyto_noco_ext$p.value,acsf_cyto_noco_np$p.value,acsf_cyto_noco_retr$p.value,acsf_cyto_noco_lp$p.value,acsf_cyto_noco_er$p.value,acsf_cyto_noco_re$p.value,acsf_cyto_noco_stab$p.value), "BH")

acsf_cytoD_v_noco <- chisq.test(cytoSkel[1:7,c(2,3)]) #cytochalasin-D vs nocodazole

acsf_cytoD_v_noco_ext <- chisq.test(rbind(cytoSkel[1,c(2,3)], colSums(cytoSkel[2:7,c(2,3)])));
acsf_cytoD_v_noco_np <- chisq.test(rbind(cytoSkel[2,c(2,3)], colSums(cytoSkel[c(1,3:7),c(2,3)])));
acsf_cytoD_v_noco_retr <- chisq.test(rbind(cytoSkel[3,c(2,3)], colSums(cytoSkel[c(1:2,4:7),c(2,3)])));
acsf_cytoD_v_noco_lp <- chisq.test(rbind(cytoSkel[4,c(2,3)], colSums(cytoSkel[c(1:3,5:7),c(2,3)])));
acsf_cytoD_v_noco_er <- chisq.test(rbind(cytoSkel[5,c(2,3)], colSums(cytoSkel[c(1:4,6:7),c(2,3)])));
acsf_cytoD_v_noco_re <- chisq.test(rbind(cytoSkel[6,c(2,3)], colSums(cytoSkel[c(1:5,7),c(2,3)])));
acsf_cytoD_v_noco_stab <- chisq.test(rbind(cytoSkel[7,c(2,3)], colSums(cytoSkel[c(1:6),c(2,3)])));

acsf_cytoD_v_both <- fisher.test(cytoSkel[1:7,c(2,4)]) #cytochalasin-D vs. cytochalasin-D + nocodazole

acsf_cytoD_v_both_ext <- fisher.test(rbind(cytoSkel[1,c(2,4)], colSums(cytoSkel[2:7,c(2,4)])));
acsf_cytoD_v_both_np <- fisher.test(rbind(cytoSkel[2,c(2,4)], colSums(cytoSkel[c(1,3:7),c(2,4)])));
acsf_cytoD_v_both_retr <- fisher.test(rbind(cytoSkel[3,c(2,4)], colSums(cytoSkel[c(1:2,4:7),c(2,4)])));
acsf_cytoD_v_both_lp <- fisher.test(rbind(cytoSkel[4,c(2,4)], colSums(cytoSkel[c(1:3,5:7),c(2,4)])));
acsf_cytoD_v_both_er <- fisher.test(rbind(cytoSkel[5,c(2,4)], colSums(cytoSkel[c(1:4,6:7),c(2,4)])));
acsf_cytoD_v_both_re <- fisher.test(rbind(cytoSkel[6,c(2,4)], colSums(cytoSkel[c(1:5,7),c(2,4)])));
acsf_cytoD_v_both_stab <- fisher.test(rbind(cytoSkel[7,c(2,4)], colSums(cytoSkel[c(1:6),c(2,4)])));

acsf_noco_v_both <- fisher.test(cytoSkel[1:7,c(3,4)], workspace = 1000000) #nocodazole vs. cytochalasin-D + nocodazole

acsf_noco_v_both_ext <- fisher.test(rbind(cytoSkel[1,c(3,4)], colSums(cytoSkel[2:7,c(3,4)])));
acsf_noco_v_both_np <- fisher.test(rbind(cytoSkel[2,c(3,4)], colSums(cytoSkel[c(1,3:7),c(3,4)])));
acsf_noco_v_both_retr <- fisher.test(rbind(cytoSkel[3,c(3,4)], colSums(cytoSkel[c(1:2,4:7),c(3,4)])));
acsf_noco_v_both_lp <- fisher.test(rbind(cytoSkel[4,c(3,4)], colSums(cytoSkel[c(1:3,5:7),c(3,4)])));
acsf_noco_v_both_er <- fisher.test(rbind(cytoSkel[5,c(3,4)], colSums(cytoSkel[c(1:4,6:7),c(3,4)])));
acsf_noco_v_both_re <- fisher.test(rbind(cytoSkel[6,c(3,4)], colSums(cytoSkel[c(1:5,7),c(3,4)])));
acsf_noco_v_both_stab <- fisher.test(rbind(cytoSkel[7,c(3,4)], colSums(cytoSkel[c(1:6),c(3,4)])));

#do pairwise Wilcoxon test for difference in proportion stable processes during cytoskeletal manipulations (Fig 1G, right)
fname <- paste(base, "Figure 1-source data 2.csv",sep="") 
cyto_props <- read.csv(fname)[1:12,15:18]
colnames(cyto_props) <- c("acsf","cyto","noco","cyt_noc")

t1 <- wilcox.test(cyto_props$acsf, cyto_props$cyto)
t2 <- wilcox.test(cyto_props$acsf, cyto_props$noco)
t3 <- wilcox.test(cyto_props$acsf, cyto_props$cyt_noc)
t4 <- wilcox.test(cyto_props$cyto, cyto_props$noco)
t5 <- wilcox.test(cyto_props$cyt_noc, cyto_props$noco)
t6 <- wilcox.test(cyto_props$cyto, cyto_props$cyt_noc)
p_corr <- p.adjust(c(t1$p.value,t2$p.value,t3$p.value,t4$p.value,t5$p.value,t6$p.value), "BH")

data_long <- gather(cyto_props, factor_key=TRUE)
data_long%>% group_by(key)%>%
  summarise(median= median(value, na.rm = TRUE), q1 = quantile(value, prob = 0.25,na.rm = TRUE),q3 = quantile(value, prob = 0.75,na.rm = TRUE))

# now do the same process for EGF  (Fig. 1H)
fname <- paste(base, "Figure 1-source data 2.csv",sep="") 
egf <- read.csv(fname)[,19:22]

#compare overall proportions of motility categories, p8/9
acsf_EGF_p8 <- chisq.test(egf[,1:2])

acsf_EGF_p8_ext <- chisq.test(rbind(egf[1,1:2], colSums(egf[2:7,1:2])));
acsf_EGF_p8_np <- chisq.test(rbind(egf[2,1:2], colSums(egf[c(1,3:7),1:2])));
acsf_EGF_p8_retr <- chisq.test(rbind(egf[3,1:2], colSums(egf[c(1:2,4:7),1:2])));
acsf_EGF_p8_lp <- chisq.test(rbind(egf[4,1:2], colSums(egf[c(1:3,5:7),1:2])));
acsf_EGF_p8_er <- chisq.test(rbind(egf[5,1:2], colSums(egf[c(1:4,6:7),1:2])));
acsf_EGF_p8_re <- chisq.test(rbind(egf[6,1:2], colSums(egf[c(1:5,7),1:2])));
acsf_EGF_p8_stab <- chisq.test(rbind(egf[7,1:2], colSums(egf[c(1:6),1:2])));
p_corr <- p.adjust(c(acsf_EGF_p8_ext$p.value,acsf_EGF_p8_np$p.value,acsf_EGF_p8_retr$p.value,acsf_EGF_p8_lp$p.value,acsf_EGF_p8_er$p.value,acsf_EGF_p8_re$p.value,acsf_EGF_p8_stab$p.value), "BH")

#compare overall proportions of motility categories, p17
#acsf_EGF_p17 <- chisq.test(egf[,3:4])
acsf_EGF_p17 <- fisher.test(egf[,3:4])

#acsf_EGF_p17_ext <- chisq.test(rbind(egf[1,3:4], colSums(egf[2:7,3:4])));
 acsf_EGF_p17_ext <- fisher.test(rbind(egf[1,3:4], colSums(egf[2:7,3:4])));
#acsf_EGF_p17_np <- chisq.test(rbind(egf[2,3:4], colSums(egf[c(1,3:7),3:4])));
 acsf_EGF_p17_np <- fisher.test(rbind(egf[2,3:4], colSums(egf[c(1,3:7),3:4])));
#acsf_EGF_p17_retr <- chisq.test(rbind(egf[3,3:4], colSums(egf[c(1:2,4:7),3:4])));
 acsf_EGF_p17_retr <- fisher.test(rbind(egf[3,3:4], colSums(egf[c(1:2,4:7),3:4])));
#acsf_EGF_p17_lp <- chisq.test(rbind(egf[4,3:4], colSums(egf[c(1:3,5:7),3:4])));
 acsf_EGF_p17_lp <- fisher.test(rbind(egf[4,3:4], colSums(egf[c(1:3,5:7),3:4])));
#acsf_EGF_p17_er <- chisq.test(rbind(egf[5,3:4], colSums(egf[c(1:4,6:7),3:4])));
 acsf_EGF_p17_er <- fisher.test(rbind(egf[5,3:4], colSums(egf[c(1:4,6:7),3:4])));
#acsf_EGF_p17_re <- chisq.test(rbind(egf[6,3:4], colSums(egf[c(1:5,7),3:4])));
 acsf_EGF_p17_re <- fisher.test(rbind(egf[6,3:4], colSums(egf[c(1:5,7),3:4])));
#acsf_EGF_p17_stab <- chisq.test(rbind(egf[7,3:4], colSums(egf[c(1:6),3:4])));
 acsf_EGF_p17_stab <- fisher.test(rbind(egf[7,3:4], colSums(egf[c(1:6),3:4])));
p_corr <- p.adjust(c(acsf_EGF_p17_ext$p.value,acsf_EGF_p17_np$p.value,acsf_EGF_p17_retr$p.value,acsf_EGF_p17_lp$p.value,acsf_EGF_p17_er$p.value,acsf_EGF_p17_re$p.value,acsf_EGF_p17_stab$p.value), "BH")

#do pairwise Wilcoxon test for EGF (Fig 1H, right)
fname <- paste(base, "Figure 1-source data 2.csv",sep="") 
egf_props <- read.csv(fname)[,32:35]
colnames(egf_props) <- c("ctrl_p9","drug_p9","ctrl_p17","drug_p17")

p1 <- wilcox.test(egf_props$ctrl_p9, egf_props$drug_p9,paired = TRUE)
p2 <- wilcox.test(egf_props$ctrl_p17, egf_props$drug_p17,paired = TRUE)
p_corr <- p.adjust(c(p1$p.value,p2$p.value), "BH")

data_long <- gather(egf_props, factor_key=TRUE)
data_long%>% group_by(key)%>%
  summarise(median= median(value, na.rm = TRUE), q1 = quantile(value, prob = 0.25,na.rm = TRUE),q3 = quantile(value, prob = 0.75,na.rm = TRUE))

