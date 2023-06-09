#load package
library(corrplot)


####Males
plum.m<-read.csv("Form_plum_mean.2.0_m.csv",row.names = 1)
head(plum.m)
str(plum.m)

colnames(plum.m)
colnames(plum.m) <- c("sex","Rdor","Gdor","Bdor","Ldor","MEdor","SEdor","CTdor","Rven","Gven","Bven","Lven","MEven","SEven","CTven","Rwin","Gwin","Bwin","Lwin","MEwin","SEwin","CTwin")
head(plum.m)

#remove sex column
plum.mc<-cor(plum.m[-1])
plum.mc

#Plot
col3 <- colorRampPalette(c("yellow1","red")) 
corrplot.mixed(plum.mc,lower.col = col3(100),upper.col=col3(100),tl.col = "black",upper="ellipse")


############################################################################################################################################################################################
##############################################################################################

####Females
plum.f<-read.csv("Form_plum_mean.2.0_f.csv",row.names = 1)
head(plum.f)
str(plum.f)

colnames(plum.f)
colnames(plum.f) <- c("sex","Rdor","Gdor","Bdor","Ldor","MEdor","SEdor","CTdor","Rven","Gven","Bven","Lven","MEven","SEven","CTven","Rwin","Gwin","Bwin","Lwin","MEwin","SEwin","CTwin")
head(plum.f)

#remove sex
plum.fc<-cor(plum.f[-1])
plum.fc


#Plot
col3 <- colorRampPalette(c("yellow1","red"))
corrplot.mixed(plum.fc,lower.col = col3(100),upper.col=col3(100),tl.col = "black",upper="ellipse")


