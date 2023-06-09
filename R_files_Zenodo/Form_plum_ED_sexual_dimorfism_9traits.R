library(tidyr)
library(dplyr)


Form.plum.mean.ut<-read.csv("Form_plum_mean_ut_9traits.csv",row.names = 1)
head(Form.plum.mean.ut)

#log transformation
Form.plum.mean.ut.l<-log(Form.plum.mean.ut)
head(Form.plum.mean.ut.l)

#calculate euclidean distance
plum_ed<-dist(Form.plum.mean.ut.l, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

plum_ed_mat<-as.matrix(plum_ed)
head(plum_ed_mat)
write.table(plum_ed_mat,"Form_plum_ed_9traits.csv",sep = ",",row.names = F)
