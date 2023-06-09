library(tidyr)
library(dplyr)


Form.vocal.mean.ut<-read.csv("Form_vocal_mean_ut_9traits.csv",row.names = 1)
head(Form.vocal.mean.ut)

#log transformation
Form.vocal.mean.ut.l<-log(Form.vocal.mean.ut)
head(Form.vocal.mean.ut.l)

#calculate euclidean distance
vocal_ed<-dist(Form.vocal.mean.ut.l, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

vocal_ed_mat<-as.matrix(vocal_ed)
head(vocal_ed_mat)

write.table(vocal_ed_mat,"Form_vocal_ed_9traits.csv",sep = ",",row.names = F)
