#### SCRIPT FROM PERRARD & LOOPE, 2015. Patriline differences reveal genetic influence on forewing size and shape in a yellowjacket wasp (Hymenoptera: Vespidae: Vespula flavopilosa Jacobson, 1978)
#### In order to run this script, the Rmorph library is required (available by request to Michel Baylac: baylac@mnhn.fr. If unsuccessful, you can also contact the first author: perrard@mnhn.fr)
#### The following function is also required:
# regallom: allometry analysis based on multivariate regression.
regallom <- function(mat, size, groups, logCS= TRUE) {
	n <- nrow(mat) ; p <- ncol(mat)
	if (length(size) != n) stop(cat("Rmorph(regallom): length of size (", length(size),") should be equal to the rows number of mat (", n, ") exiting...\\n"))
	if (logCS) lCS <- log(size)
	else lCS <- size
	MATS<-cbind(lCS,mat)
	matc <- cmeangr(MATS, groups)$gcmat[,2:(p+1)]
	lCSc<- cmeangr(MATS, groups)$gcmat[,1]
	lCSm<-mean(lCS)
	CAC <- lm(matc~lCSc)
	slopes<-CAC$coef[2,]
	mmat <- colMeans(mat)
	intp <- mmat - slopes*lCSm
	pred <- matrix(0,n,p)
	for (i in 1:n) {
		for (j in 1:p) {
			pred[i,j] <- slopes[j]*lCS[i]+intp[j]
		}
	}
	RSCg <- mat - pred	
	RSC <- mvreg(mat, lCS)$res
	pRS <- Rpca(RSC)
	pRSg <- Rpca(RSCg)
	return(list(logCS= logCS, lCS= lCS, intp=intp, slopes= slopes, fit=pred, RSC=RSC, RSCg=RSCg, pRS= pRS,pRSg = pRSg))
}	
	
#### ANALYSIS OF DATA FROM PERRARD & LOOPE, 2015, PLOS ONE.
source("Rmorph.R")
read.table("Landmark-Data-Perrard-Loope-2015.txt",h=T)->Data #File available on Dryad.
FWall<-Data[,6:43]
Wdata<-Data[,1:5]
liensFW<-c(1,2,3,4,5,12,13,19,18,17,16,15,14,7,8,9,10,2,NA,3,9,NA,4,6,7,NA,6,13,NA,8,11,15,NA,7,14,18,NA,1,5)

###################################
###### Test error measurement #####
###################################
# Selection of wing shapes prepared twice and measured twice.
FWtest<-FWall[-which(Wdata$Montage==0),]
Wdatatest<-Wdata[-which(Wdata$Montage==0),]

# Superimposition of the shapes.
gpa(FWtest,links=liensFW)->gpaFWtest


# Computation of the sum of squares.
TestF<-Wdatatest$Patriline
TestI<-Wdatatest$Indiv
TestM<-Wdatatest$Montage
TestD<-Wdatatest$Digitization
TestM2<-paste(TestI,TestM,sep="")
Test<-gpaFWtest$scores[,1:34]

mod<-lm(Test~TestF+TestI+TestM2)
 effects<-mod$effects
 SSfat<-sum(diag(crossprod(effects[which(mod$assign==1),,drop=F])))
 SSind<-sum(diag(crossprod(effects[which(mod$assign==2),,drop=F])))
 SSmont<-sum(diag(crossprod(effects[which(mod$assign==3),,drop=F])))
 SSres<-sum(diag(var(mod$residuals)*57))

 SSfat
 SSind
 SSmont
 SSres
 
 # Computation of the degrees of freedom:
 P<-38 #(=2 x number of landmarks)
 n<-39
 r<-m<-2
 dffat<-(P-4)*3
 dfind<-((n-1)*(P-4)*4)
 dfmont<-((m-1)*n*(P-4)*4)
 dfres<-((r-1)*m*n*(P-4)*4)

 dffat
 dfind
 dfmont
 dfres
 
 # Computation of the mean squares:
MSfat<-SSfat/dffat 
MSind<-SSind/dfind 
MSmont<-SSmont/dfmont 
MSres<-SSres/dfres 

MSfat
MSind
MSmont
MSres

##########################################################################
########## STUDY OF THE WING VARIATION BETWEEN PATRILINES ################
##########################################################################

FW<-FWall[c(1:181,seq(182,259,2)),] # Selection of a single measurement for every specimen.
FWdata<-Wdata[c(1:181,seq(182,259,2)),] # Selection of a single measurement for every specimen.
for (i in 1:dim(FWdata)[2]) FWdata[,i]<-as.factor(as.vector(FWdata[,i]))

# Superimposition of all the wings at once
gpa(FW,links=liensFW)->gpaFW

#Selection of balanced samples in the main patrilines of Colonies 1 and 2.
which(FWdata$Patriline=="F1_1")->FW11
sample(FW11,34)->fw11
# In the study, the following random sample was selected:
#fw11<-c(158,219,163,177,197,127,199,108,138,135,181,114,179,125,209,210,116,120,174,205,196,140,171,107,214,152,113,198,175,162,124,157,193,109) 
which(FWdata$Patriline=="F1_2")->fw12

which(FWdata$Patriline=="F2_1")->FW21
sample(FW21,34)->fw21
# In the study, the following random sample was selected:
#fw21<-c(6,9,12,14,16,21,23,25,28,29,30,42,43,45,47,48,50,55,56,63,65,68,70,72,73,77,81,84,85,88,90,94,95,101)
which(FWdata$Patriline=="F2_2")->FW22
sample(FW22,34)->fw22
# In the study, the following random sample was selected:
#fw22<-c(8,11,13,15,17,19,20,24,26,27,36,39,49,51,53,54,57,59,67,71,74,76,78,80,82,83,86,89,91,96,97,98,99,100) 

############################
##### Analysis of sizes #####
############################
csize<-log(gpaFW$size)

summary(aov(csize~FWdata$Colony*FWdata$Patriline)) # Comparison of wing sizes between colonies and patrilines

# Size variation inter-colony
ms2<-mean(csize[which(FWdata$Colony==1)])
ms1<-mean(csize[which(FWdata$Colony==2)])
var(c(ms1,ms2)) 

# Size variation inter-patrilines
meangr(csize,group=FWdata$Patriline)->msFath
var(msFath) #inter-patrilines

# Size variation intra-patrilines
sizIP<-csize
for (i in 1:8) sizIP[which(FWdata$Patriline==levels(FWdata$Patriline)[i])]<-csize[which(FWdata$Patriline==levels(FWdata$Patriline)[i])]-msFath[which(levels(FWdata$Patriline)[i]==rownames(msFath))]
var(sizIP) ##intra-patriline

# Comparison between main patrilines of Colony 1
SizefwH1<-log(gpaFW$size[c(fw11,fw12)])
t.test(SizefwH1[1:34],SizefwH1[35:68])

# Comparison between main patrilines of Colony 2
SizefwH2<-log(gpaFW$size[c(fw21,fw22)])
t.test(SizefwH2[1:34],SizefwH2[35:68])

#############################
##### Analysis of shapes #####
#############################

# Variation in the shape-space for all wing shapes
plotg(gpaFW$scores[,1:2],group=FWdata$Patriline)

summary(manova(gpaFW$scores[,1:34]~csize*FWdata$Colony*FWdata$Patriline))# influence of size, colony and patrilines on the wing shapes

#Reassignment with cross-validation using canonical variate analyses:
library(MASS)
for (i in 2:34) assign(paste("lda",i,"PC",sep=""),lda(gpaFW$scores[,1:i],FWdata$Patriline,CV=TRUE))
table(FWdata$Patriline,lda34PC$class) #the number in "lda34PC$class" indicate the number of Principal Components used for the CVA. It goes from 2 to 34.
100*sum(diag(table(FWdata$Patriline,lda34PC$class)))/220 # Percentage of correct reassignment of specimens to their patriline based the wing shape variables.
100*sum(diag(table(FWdata$Patriline,lda2PC$class)))/220
100*sum(diag(table(FWdata$Patriline,lda19PC$class)))/220

# difference in wing shapes between the main patrilines of Colony 1
library(Hotelling)
ShapefwH1<-gpaFW$scores[c(fw11,fw12),1:34]
print(hotelling.test(ShapefwH1[1:34,],ShapefwH1[35:68,]))
 
# difference in wing shapes between the main patrilines of Colony 2
ShapefwH2<-gpaFW$scores[c(fw21,fw22),1:34]
print( hotelling.test(ShapefwH2[1:34,],ShapefwH2[35:68,]))
 
rallFW2<-regallom(mat=ShapefwH2,logCS=F,size=SizefwH2,groups=c(rep(1,34),rep(2,34)))
print(hotelling.test(rallFW2$RSCg[1:34,],rallFW2$RSCg[35:68,]))
 
 # Comparison of shape changes using angles between shape-change vectors in the shape-space
 
Colo<-c(rep(1,68),rep(2,68))
Pat<-c(rep(1,34),rep(2,34),rep(3,34),rep(4,34))
 ShapefwH<-rbind(ShapefwH1,ShapefwH2)
 cmeangr(ShapefwH,group=Colo)$gmeans->COLO
 a<-(COLO[1,]-COLO[2,])
 cmeangr(ShapefwH,group=Pat)$gmeans->PATRI
 b<-(PATRI[1,]-PATRI[2,])
 c<-(PATRI[3,]-PATRI[4,])
 A<-a/sqrt(t(a)%*%a) #Vector "normed"
 B<-b/sqrt(t(b)%*%b) 
 C<-c/sqrt(t(c)%*%c)
 theta1 <- acos( sum(A*B) / ( sqrt(sum(A * A)) * sqrt(sum(B * B)) ) )
 theta2 <- acos( sum(A*C) / ( sqrt(sum(A * A)) * sqrt(sum(C * C)) ) )
 theta3 <- acos( sum(B*C) / ( sqrt(sum(B * B)) * sqrt(sum(C * C)) ) )
 AngleCP1<-180*theta1/pi #Angle between vectors between colonies and between main patrilines of colony 1 
 AngleCP2<-180*theta2/pi #Angle between vectors between colonies and between main patrilines of colony 2
 AngleP1P2<-180*theta3/pi #Angle between vectors between main patrilines of colony 1 and between main patrilines of colony 2

 AngleCP1
 AngleCP2
 AngleP1P2
 
ANG<-1:10000
for (i in 1:10000) {
ab<-sample(1:136,4,replace=FALSE)
 a<-(ShapefwH[ab[1],]-ShapefwH[ab[2],])
 b<-(ShapefwH[ab[3],]-ShapefwH[ab[4],])
 A<-a/sqrt(t(a)%*%a) #Vector "normed"
 B<-b/sqrt(t(b)%*%b) #Vector "normed"
 theta1 <- acos( sum(A*B) / ( sqrt(sum(A * A)) * sqrt(sum(B * B)) ) )
 ANG[i]<-180*theta1/pi
 }
 ANG<-ANG[which(ANG>1)]
 length(which(ANG<AngleCP1))/length(ANG)
 length(which(ANG<AngleCP2))/length(ANG)
 length(which(ANG<AngleP1P2))/length(ANG) 
 plot(density(ANG))
 abline(v=c(AngleCP1,AngleCP2,AngleP1P2))
 
 ##### Visualisations of shape changes
 library(geomorph)
 fwH<-c(fw11,fw12,fw21,fw22)
 cmeangr(gpaFW$res[fwH,],group=Colo)$gmeans->COLO
 cmeangr(gpaFW$res[fwH,],group=Pat)$gmeans->PATRI 
  plotRefToTarget(mat2col(colMeans(COLO[1:2,])),mat2col(COLO[1,]),mag=10)
x11()
 plotRefToTarget(mat2col(colMeans(COLO[1:2,])),mat2col(COLO[2,]),mag=10)
x11()
 plotRefToTarget(mat2col(colMeans(PATRI[1:2,])),mat2col(PATRI[1,]),mag=10)
x11()
 plotRefToTarget(mat2col(colMeans(PATRI[1:2,])),mat2col(PATRI[2,]),mag=10)
x11()
 plotRefToTarget(mat2col(colMeans(PATRI[3:4,])),mat2col(PATRI[3,]),mag=10)
x11()
 plotRefToTarget(mat2col(colMeans(PATRI[3:4,])),mat2col(PATRI[4,]),mag=10)
