#loading the library for multiple comparison in fishers exact test
# Running fishers exact test on preening data for Y45 and G123 

library(RVAideMemoire)

#Y45 data
dat = matrix(c(0,30,6,0,35,0,0,1,1,13,27,109,7,60,43,129, 4, 2, 3, 22), ncol = 5, byrow = FALSE)
colnames(dat) = c("CtrBrd","RdThrBrd","BlThrMirr","RdThrMirr", "RdHdMirr")
rownames(dat) = c("sess_1","sess_2", "sess_3", "sess_4")
dat = as.table(dat)

#Running the mutiple comparison test
fisher.multcomp(dat, p.method = 'bonferroni')

#G123 DATA
dat = matrix(c(6,2,0,1,7,0,1,6,0,0,5,25,0,1,29,59,0,0,0,12), ncol = 5, byrow = FALSE)
colnames(dat) = c("CtrBrd","RdThrBrd","BlThrMirr","RdThrMirr", "RdHdMirr")
rownames(dat) = c("sess_1","sess_2", "sess_3", "sess_4")
dat = as.table(dat)

#Running the mutiple comparison test
fisher.multcomp(dat, p.method = "bonferroni")

