#load package
library(corrplot)

####Males
voc.m<-read.csv("Form_vocal_mean_m_ok.csv",row.names = 1)
head(voc.m)
str(voc.m)

colnames(voc.m)
colnames(voc.m) <- c("sex","Ncou","Ntyp","MinF","MaxF","AggE","AvgE","Ldur","PeakF","Ndiv","Nrat","LBW","Lmr")
head(voc.m)

#remove sex and entropy columns
voc.mcc<-cor(voc.m[-c(1,6,7)])
voc.mcc

#Plot
col3 <- colorRampPalette(c("yellow1","red")) 
corrplot.mixed(voc.mcc,lower.col = col3(100),upper.col=col3(100),tl.col = "black",order = "alphabet",addgrid.col = "black",upper="ellipse")



#####################################################################################################################################################################################################################################################################


####Females
voc.f<-read.csv("Form_vocal_mean_f_ok.csv",row.names = 1)
head(voc.f)
str(voc.f)

colnames(voc.f)
colnames(voc.f) <- c("sex","Ncou","Ntyp","MinF","MaxF","AggE","AvgE","Ldur","PeakF","Ndiv","Nrat","LBW","Lmr")
head(voc.f)


#remove sex and entropy columns
voc.fcc<-cor(voc.f[-c(1,6,7)])


#Plot
col3 <- colorRampPalette(c("yellow1","red")) 
corrplot.mixed(voc.fcc,lower.col = col3(100),upper.col=col3(100),tl.col = "black",order = "alphabet",addgrid.col = "black",upper="ellipse")


