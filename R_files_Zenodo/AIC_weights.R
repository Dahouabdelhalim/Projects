### Model averaging theta estimates
### Lagomarsino et al 2017

require(stats)
require(plotrix)
setwd("~/Desktop/For new OU analyses")
dat <- read.table("Camp_summary_poll22Feb2017.txt",row.names=NULL)
dat <- dat[,-1]
dat$trait <- as.character(dat$trait)
trait <- unique(dat$trait)    

##plotting AICc scores for different models and percent aerenchyma
#dat <- dat[which(dat$scenario%in%c(1)),]
#dat <- dat[which(dat$trait%in%trait),]
#dat[which(dat$scenario==5)[c(10:1350)],]$aic_c <- dat[which(dat$scenario==5)[c(10:1350)],]$aic_c-#dat[which(dat$scenario==1),]$aic_c
#quartz()
#x <- dat[which(dat$scenario==2),'aic_c']
#y <- dat[which(dat$scenario==5),'aic_c']
#plot(y~x,cex=.2,ylim=c(8,45),xlim=c(8,45),xlab="OU1 AICc",ylab="OU3 AICc")
#abline(a=0,b=1,col="red",lty=2,lwd=2)
#scenario <- c(2:6)
#model <- c("OU1","OU2","OU2","OU3","OU3")
#dat <- dat[which(dat$model=="scenario"),]

###function to calculate aikake weights

akaike_weights <- function(dat){
	
	aic_diff <- dat - min(dat)
	w <- exp(-0.5*(aic_diff))/sum(exp(-0.5*(aic_diff)))
	return(w)
	
}



###calcultating Aikake weights for all models

dat_weight <- as.data.frame(matrix(NA,nrow=length(trait),ncol=5))
colnames(dat_weight) <- c("trait","BM","OU1","OU2","OU3")
dat_weight$trait <- trait

theta_est <- as.data.frame(matrix(NA,nrow=0,ncol=10))
colnames(theta_est) <- c("trait","theta1","theta1_025","theta1_975","theta2","theta2_025","theta2_975","theta3","theta3_025","theta4_975")

mod_count <- as.data.frame(matrix(0,nrow=length(trait),ncol=5))
colnames(mod_count) <- c("trait","BM","OU1","OU2","OU3")
mod_count$trait <- trait


###for ace reconstruction model; run for pollination data

for(i in 1:length(trait)){
	
	temp <- dat[which(dat$trait==trait[i]),]
	temp <- unique(temp[which(temp$scenario%in%c(1,2,3,4)),])
	weights <- as.matrix(aggregate(temp$aic_c,by=list(temp$tree_num),FUN=akaike_weights))
	dat_weight[i,2:5] <-  apply(weights[,2:dim(weights)[2]],2,mean)
	
	
}

###counting the number of trees supporting each model -- model count

for(i in 1:length(trait)){
		temp <- dat[which(dat$trait==trait[i]),] 
		temp <- unique(temp[which(temp$scenario%in%c(1,2,3,4)),])
		weights <- as.matrix(aggregate(temp$aic_c,by=list(temp$tree_num),FUN=akaike_weights))
		dat_weight[i,2:5] <- apply(weights[,2:dim(weights)[2]],2,mean)
	for(j in 1:100){
		
		mod_count[i,(which(weights[j,c(2:5)]==max(weights[j,c(2:5)]))+1)] <- mod_count[i,(which(weights[j,c(2:5)]==max(weights[j,c(2:5)]))+1)] + 1
		
	}
	
}


###calculating model-weighted parameter estimates for theta, as well as .025 and .975 quantiles for the distribtuon of values across 1000 trees

fin_dat <- as.data.frame(matrix(NA,nrow=0,ncol=8))
colnames(fin_dat) <- c("trait","tree_num","theta1","theta1_se","theta2","theta2_se","theta3","theta3_se")

for(i in 1:length(trait)){
	
	for(j in 1:100){
		
		print(paste(i,j,sep="_"))
	    temp <- dat[which(dat$tree_num==j&dat$trait==trait[i]),]
	    temp <- unique(temp[which(temp$scenario%in%c(2,3,4)),])
	    aic_diff <- temp$aic_c - min(temp$aic_c)
		w <- exp(-0.5*(aic_diff))/sum(exp(-0.5*(aic_diff)))
		temp[1,c("theta2","theta3")] <- temp[1,]$theta1
		temp[1,c("theta2_se","theta3_se")] <- temp[1,]$theta1_se
		temp[2,c("theta3")] <- temp[2,c("theta1")]
		temp[2,c("theta3_se")] <- temp[2,c("theta1_se")]
		theta1 <- weighted.mean(temp$theta1,w)
		theta2 <- weighted.mean(temp$theta2,w)
		theta3 <- weighted.mean(temp$theta3,w)
		temp$theta1_se2 <- sqrt(((theta1 - temp$theta1)^2) + (temp$theta1_se^2))
		temp$theta2_se2 <- sqrt(((theta2 - temp$theta2)^2) + (temp$theta2_se^2))
		temp$theta3_se2 <- sqrt(((theta3 - temp$theta3)^2) + (temp$theta3_se^2))
		temp$theta1_se_w <- temp$theta1_se2 * w
		temp$theta2_se_w <- temp$theta2_se2 * w
		temp$theta3_se_w <- temp$theta3_se2 * w
		theta1_se <- sum(temp$theta1_se_w)
		theta2_se <- sum(temp$theta2_se_w)
		theta3_se <- sum(temp$theta3_se_w)
		fin_dat <- rbind(fin_dat,c(i,j,theta1,theta1_se,theta2,theta2_se,theta3,theta3_se))
		
	}
	
	 colnames(fin_dat) <- c("trait","tree_num","theta1","theta1_se","theta2","theta2_se","theta3","theta3_se")
  	 fin_dat[which(fin_dat$trait==i),]$trait <- as.character(trait[i])
		
}

##calculating weighted sums of SE across all 1000 trees
CI <- rep(qt(.95,15),3)
 se_final <- NULL
 
for(i in 2:length(trait)){
	
	temp <- fin_dat[which(fin_dat$trait==trait[i]),]
	sumSE <- c(sqrt(sum((temp$theta1_se)^2)/100),sqrt(sum((temp$theta2_se)^2)/100),sqrt(sum((temp$theta3_se)^2)/100))
	point_SE <- c(std.error(temp$theta1),std.error(temp$theta2),std.error(temp$theta3))
	mean_point <- c(mean(temp$theta1),mean(temp$theta2),mean(temp$theta3))
	ltr <- as.data.frame(cbind(rep(trait[i],3),mean_point,sumSE,point_SE,CI,c("theta1","theta2","theta3")))
	se_final <- rbind(se_final,ltr)
} 
 
 write.csv(se_final,file="Theta_summary.csv")
 
###plotting the mean point estimate and 95% CI based on weighted sum of errors across all 1000 trees
PerAir_J <- ggplot(ltr, aes(x=theta, y=mean_point)) +
    xlab("") +
    ylab("") +
    geom_errorbar(width=.1, aes(ymin=mean_point+sumSE, ymax=mean_point-sumSE),colour="red") +
    geom_errorbar(width=.1, aes(ymin=mean_point+point_SE, ymax=mean_point-point_SE)) +
    geom_point(shape=21, size=3, fill="black") +
    ggtitle("Percent Aerenchyma") +
    theme_bw()
   
    	
    
      
###  fin_dat[,c("theta1","theta2","theta3")] <- (10^fin_dat[,c("theta1","theta2","theta3")])

###summarizing theta data

trait_mean <- aggregate(fin_dat[,c(3,5,7)],by=list(fin_dat$trait),FUN=mean)
trait_std_err1 <- aggregate(fin_dat[,c(3,5,7)],by=list(fin_dat$trait),FUN=std.error)
trait_std_err2 <- aggregate(fin_dat[,c(4,6,8)],by=list(fin_dat$trait),FUN=std.error)
trait_025 <- aggregate(fin_dat[,c(3,5,7)],by=list(fin_dat$trait),FUN=quantile,probs=.025,na.rm=TRUE)
trait_975 <- aggregate(fin_dat[,c(3,5,7)],by=list(fin_dat$trait),FUN=quantile,probs=.975,na.rm=TRUE)

###boxplot figure for AIC weighted theta values
colnames(fin_dat) <- c("trait","tree_num","theta1","theta2","theta3")
trait2 <- unique(fin_dat$trait)
quartz(width=6,height=3)
layout(matrix(1:2,nrow=1,byrow=TRUE))
#par(mar=c(.75,.75,.75,.75))
#par(mfrow=c(10,2),mar=c(5.1,4.1,4.1,2.1),oma=c(0,0,0,0))

###reading in original trait values to add to the theta estimates on the box plot

dat_mean <- read.csv("~/Desktop/VP_revision_analyses/PCA_all/VP_PCA_traits.csv")

for(i in 1:length(trait)){
	
	temp <- fin_dat[which(fin_dat$trait==trait[i]),]
	colnames(temp)[3:5] <- c("T","A/T","VP")
	#temp[,c("nonVP","gen","VP")] <- (10^(temp[,c("nonVP","gen","VP")]))
	par(mar=c(1.0,2.5,1.0,1.0))
	boxplot(temp[,c("T","A/T","VP")],outline=FALSE,notch=TRUE,range=.95,cex.axis=1.5,xaxt='n')
	 
}

fin_weight <- as.data.frame(matrix(NA,nrow=length(trait),ncol=10))
colnames(fin_weight) <- c("trait","BM","OU1","OU2.1_lagrange","OU2.2_lagrange","OU3.1_lagrange","OU3.2_lagrange","OU2_ace","OU3.1_ace","OU3.2_ace")

###calculating mean AICc values to report in table  (along with model weights)
fin_weight[,'trait'] <- trait
 
for(i in 1:length(trait)){

		temp <- dat[which(dat$trait==trait[i]),]
		paste(trait[i])
		mean(mean(temp$aic_c))
		std.error(temp$aic_c)
		#fin_weight[,'trait'] <- trait
		#fin_weight[i,'BM'] <- mean(temp$aic_c)
		
}