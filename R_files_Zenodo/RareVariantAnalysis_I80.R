# This script contains functions for counting how many rare alleles are shared on same side vs. diff side of highway
# Created by Elizabeth Mandeville, modified by M.E.F. LaCava and C.A. Buerkle

####Count private alleles ####
# Count # private rare alleles vs. shared across populations 
count.private <- function(g, mac, group){
  heterozygote.inds<-apply(g[apply(g, 1, function(x) sum(abs(x-1)< 0.05, na.rm=T)) == mac,], 1, function(x){which(abs(x - 1) <0.05)})
  #count.diff<-sum(apply(heterozygote.inds, 2, function(x) sum(abs(diff(group[x])))))
  #again, diff function doesn't work properly: group=0,0,1 results in 1 but group=0,1,0 results in 2 -> over-counts loci
  count.diff<-sum(apply(heterozygote.inds, 2, function(x) var(group[x])!=0)) #var of group=0 when group is all same, so !=0 when diff
  return(data.frame(count.same=ncol(heterozygote.inds) - count.diff, 
                    count.diff=count.diff, 
                    prop.private=1-(count.diff/ncol(heterozygote.inds))))
}

####Test on fake data ####
#make fake genotype data - mac 2-5, 2 loci in each mac (1 same side, 1 diff side, except for mac=6 both diff side)
fakedata <- matrix(data=0,nrow=10,ncol=10)  # loci are in rows, individuals are in columns
fakedata[1,1:2] <- 1 #mac=2, same side   
fakedata[2,1:3] <- 1 #mac=3, same side
fakedata[3,1:4] <- 1 #mac=4, same side
fakedata[4,1:5] <- 1 #mac=5, same side
fakedata[5,1:6] <- 1 #mac=6, diff side
fakedata[6,5:6] <- 1 #mac=2, diff side
fakedata[7,5:7] <- 1 #mac=3, diff side
fakedata[8,5:8] <- 1 #mac=4, diff side
fakedata[9,5:9] <- 1 #mac=5, diff side
fakedata[10,5:10] <- 1 #mac=6, diff side
fakedata

#make fake side data (ind 1-5 on one side, ind 6-10 on other side)
fakeside <- c(rep(0,5), rep(1,5))

#apply to fake data
results <- data.frame(mac=2:6, count.same=rep(0,5), count.diff=rep(0,5), prop.private=rep(0,5))
for (i in 2:6) {
  priv <- count.private(fakedata,i,fakeside)
  results$count.same[results$mac==i] <- priv$count.same
  results$count.diff[results$mac==i] <- priv$count.diff
  results$prop.private[results$mac==i] <- priv$prop.private
}
results$nloci<-results$count.same + results$count.diff
results

#### PH data ####
##Import my data - subset samples around I-80 corridor
genest.i80 <- read.table("pntest_mean_i80_173.sorted.txt", header=F)
db <- read.csv("ph_i80.csv", header=T)
#db$south is the "group" data - 1=south of i80, 0=north of i80

##Loop count.private over mac values
results <- data.frame(mac=2:10, count.same=rep(0,9), count.diff=rep(0,9), prop.private=rep(0,9))
for (i in 2:10) {
  priv <- count.private(genest.i80,i,db$south)
  results$count.same[results$mac==i] <- priv$count.same
  results$count.diff[results$mac==i] <- priv$count.diff
  results$prop.private[results$mac==i] <- priv$prop.private
}
results$nloci<-results$count.same + results$count.diff
results
results.private <- results

plot(results.private$mac, results.private$prop.private, ylim=c(0, 1),
     xlab="Rare allele count", ylab="Proportion private alleles",
     main="Proportion of loci with private alleles (same side of I-80)")


#### Permute to create null model ####
#Import null model permutations (since it takes a while to run)
#nullmodel <- read.csv("nullI80_ind173_mac2-10_100reps.csv",header=T)
##Permute group data to generate null model
# genest.i80 <- read.table("pntest_mean_i80_173.sorted.txt", header=F)
# db <- read.csv("ph_i80.csv", header=T)
# reps <- 100
# mac <- 2:10
# null <- data.frame(matrix(nrow=reps*length(mac),ncol=5))
# names(null) <- c("iteration","mac","count.same","count.diff","prop.private")
# null$iteration <- rep(1:reps,each=length(mac))
# null$mac <- rep(mac,length.out=nrow(null))
# null[,3:5] <- 0
# for (j in 1:reps) {
#   group <- sample(db$south)
#   for (i in mac) {
#       priv <- count.private(genest.i80,i,group)
#       null$count.same[null$iteration==j & null$mac==i] <- priv$count.same
#       null$count.diff[null$iteration==j & null$mac==i] <- priv$count.diff
#       null$prop.private[null$iteration==j & null$mac==i] <- priv$prop.private
#     }
#   }
# null
#write.csv(null,"nullI80_ind173_mac2-10_10reps.csv",row.names=F,quote=F)

# tmp<-numeric(500)
# for(i in 1:500){
#   tmp[i]<-count.private(genest.i80, 2, sample(db$south))$prop.private
# }

##Plot data with permuted null model
#calc quantiles
# for (i in  2:10) {
#   results.private$lower[results.private$mac==i] <- quantile(nullmodel$prop.private[nullmodel$mac==i],probs=0.025,names=F)
#   results.private$upper[results.private$mac==i] <- quantile(nullmodel$prop.private[nullmodel$mac==i],probs=0.975,names=F)
# }
# results.private
# #plot
# plot(results.private$mac, results.private$prop.private, type="l", ylim=c(0, 1),
#      xlab="Rare allele count", ylab="Proportion private alleles",
#      main="Proportion of private alleles (same side of I-80)")
# lines(2:10, 2/(2:10 + 1),col=rgb(1,0,0)) # cab: I think this is the correct expectation and the graylines are the uncertainty around it
# lines(2:10, qbinom(0.975, size=results.private$nloci, prob=2/(2:10 + 1))/results.private$nloci, col=rgb(1,0,0,alpha=0.3))
# lines(2:10, qbinom(0.025, size=results.private$nloci, prob=2/(2:10 + 1))/results.private$nloci, col=rgb(1,0,0,alpha=0.3))
# par(new=T)
# plot(results.private$mac, results.private$lower, type="l", ylim=c(0, 1), col="darkgray",lty=2,xlab=NA,ylab=NA)
# par(new=T)
# plot(results.private$mac, results.private$upper, type="l", ylim=c(0, 1), col="darkgray", lty=2,xlab=NA,ylab=NA)
# par(new=F)



#### Analytical null model for expectation ####
dbinom(2:10, 2:10, prob=68/173) + dbinom(0, 2:10, prob=68/173)
## the variance around the expectation will be driven by the number of loci that fit each mac criterion
results.private$lower <- qbinom(0.025, results.private$nloci, prob=dbinom(2:10, 2:10, prob=68/173) + dbinom(0, 2:10, prob=68/173))/results.private$nloci
results.private$upper <- qbinom(0.975, results.private$nloci, prob=dbinom(2:10, 2:10, prob=68/173) + dbinom(0, 2:10, prob=68/173))/results.private$nloci

#Plot data with analytical null model shaded behind
plot(results.private$mac, results.private$prop.private, type="l", ylim=c(0, 1),
     xlab="Rare allele count", ylab="Proportion private alleles")
polygon(c(results.private$mac,rev(results.private$mac)),c(results.private$lower,rev(results.private$upper)),
        col="gray84",border=NA)
lines(results.private$mac,results.private$prop.private)
#legend("topright",legend=c("Proportion","Null"),bty="n",col=c("black","gray84"),
#       lty=c(1,0),lwd=c(1,0),pch=c(NA,22),pt.bg=c("black","gray84"),pt.cex=2)


#### Plot for manuscript ####
pdf(file="rarealleles_i80.pdf",width=5,height=5,useDingbats=FALSE)
par(mar=c(5,5,1,1), xpd=T) #default: c(bottom, left, top, right) c(5, 4, 4, 2)
plot(results.private$mac, results.private$prop.private, type="l", ylim=c(0, 1), bty="l", lab=c(9,5,7),
     xlab="Rare variant count", ylab="Proportion private variants",cex.axis=1.2, cex.lab=1.5)
polygon(c(results.private$mac,rev(results.private$mac)),c(results.private$lower,rev(results.private$upper)),
        col="gray84",border=NA)
lines(results.private$mac,results.private$prop.private,lwd=1.5)
#legend("topright",legend=c("Proportion","Null"),bty="n",col=c("black","gray84"),
#       lty=c(1,0),lwd=c(1,0),pch=c(NA,22),pt.bg=c("black","gray84"),cex=1.5, pt.cex=2)
dev.off()
