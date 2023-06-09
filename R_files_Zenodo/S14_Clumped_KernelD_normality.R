# Evaluating differences between chalky and foliated microstructures in oysters through Kernal density estimation.
setwd("E:/Dropbox/Research/Manuscripts/[Review] GCA - Oyster microstructures/presentations")
dat<-read.csv("Oyster_microstructures_data.csv",stringsAsFactors = FALSE)
D47<-dat[,which(colnames(dat)=="D47_final")]

D47chalky<-D47[which(dat$microstructure=="chalky")]
D47chalky_SW<-shapiro.test(D47chalky) # Shapiro-Wilk normality test
D47chalky_SWp<-D47chalky_SW$p.value # If >0.05 data is normal

D47foliated<-D47[which(dat$microstructure=="foliated")]
D47foliated_SW<-shapiro.test(D47foliated) # Shapiro-Wilk normality test
D47foliated_SWp<-D47foliated_SW$p.value # If >0.05 data is normal

# Create and plot kernel density distributions
D47chalky_d<-density(D47chalky,bw=0.04,from=0.5,to=0.9)
D47foliated_d<-density(D47foliated,bw=0.04,from=0.5,to=0.9)
x11();plot(D47chalky_d, main = "D47 of microstructures", col = "grey")
lines(D47foliated_d,col="black")

# Calculate statistical difference
Ttest <- t.test(D47chalky,
                D47foliated,
                alternative = "two-sided",
                paired = FALSE,
                var.equal = FALSE,
                conf.level = 0.95
)
KS<-ks.test(D47chalky, D47foliated) # Kolmogorov-Smirnov Test
KSp<-KS$p.value # If <0.05 the samples were sampled from the same distribution

# Create plot
pdf("Kernel_distributions_D47.pdf")
plot(parameterchalky_d, main = "D47 of microstructures", col = "grey")
lines(parameterfoliated_d,col="black")
points(parameterfoliated,rep(0,length(parameterfoliated)),col="black", pch = 3)
points(parameterchalky,rep(0,length(parameterchalky)),col="grey", pch = 3)
dev.off()