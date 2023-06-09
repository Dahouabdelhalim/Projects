library(PBD);
library(ape);
#trop
trop.sp.vec <- NULL;#to store the number of sp per sim
sis.dist <- NULL# to store distancs between sister taxa
for(i in 1:100){
	temp.trop.pbd.sim <- pbd_sim(par = c(1.13, 0.15, 1.13, 0.3, 0.3), age = 6);#parameters for trop birds
	temp.trimmed.tree <- temp.trop.pbd.sim$stree_random;
	trop.sp.vec <- c(trop.sp.vec, length(temp.trimmed.tree$tip.label));
	
	#calculate sister divergence
	list.sp <- NULL;
	temp.sis.dist<-NULL;
	for(i in 1:length(temp.trimmed.tree$tip.label)){
		if(sum(temp.trimmed.tree$tip.label[i]==list.sp)==0){
		for(j in 1:length(temp.trimmed.tree$tip.label)){
			if(i != j){
				if(is.monophyletic(temp.trimmed.tree, c(temp.trimmed.tree$tip.label[i],temp.trimmed.tree$tip.label[j])) == TRUE){
					list.sp <- c(list.sp, temp.trimmed.tree$tip.label[j])
					dist <- cophenetic(temp.trimmed.tree)[temp.trimmed.tree$tip.label[i],temp.trimmed.tree$tip.label[j]]
					temp.sis.dist<-c(temp.sis.dist,dist)
				}
			}
		}
			
		}	
		}
		sis.dist<-c(sis.dist,temp.sis.dist)
	
}
trop.sis.dist <- sis.dist;
trop.sp <- trop.sp.vec;
trop.sis <- trop.sis.dist;

#temp
temp.sp.vec <- NULL;#to store the number of sp per sim
sis.dist <- NULL;
for(i in 1:100){
	temp.temp.pbd.sim <- pbd_sim(par = c(1.16, 0.5, 1.16, 0.6, 0.6), age = 6);#parameters for trop birds
	temp.trimmed.tree <- temp.temp.pbd.sim$stree_random;
	temp.sp.vec <- c(temp.sp.vec, length(temp.trimmed.tree$tip.label));
	#calculate sister divergence
	list.sp <- NULL;
	temp.sis.dist<-NULL;
	for(i in 1:length(temp.trimmed.tree$tip.label)){
		if(sum(temp.trimmed.tree$tip.label[i]==list.sp)==0){
		for(j in 1:length(temp.trimmed.tree$tip.label)){
			if(i != j){
				if(is.monophyletic(temp.trimmed.tree, c(temp.trimmed.tree$tip.label[i],temp.trimmed.tree$tip.label[j])) == TRUE){
					list.sp <- c(list.sp, temp.trimmed.tree$tip.label[j])
					dist <- cophenetic(temp.trimmed.tree)[temp.trimmed.tree$tip.label[i],temp.trimmed.tree$tip.label[j]]
					temp.sis.dist<-c(temp.sis.dist,dist)
				}
			}
		}
			
		}	
		}
		sis.dist<-c(sis.dist,temp.sis.dist)
}
temp1.sis.dist <- sis.dist;

#temp-low converting high splitting
temp2.sp.vec <- NULL;#to store the number of sp per sim
sis.dist <- NULL;
for(i in 1:100){
	temp.temp.pbd.sim <- pbd_sim(par = c(1.3, 0.15, 1.3, 0.6, 0.6), age = 6);#parameters for trop birds
	temp.trimmed.tree <- temp.temp.pbd.sim$stree_random;
	temp2.sp.vec <- c(temp2.sp.vec, length(temp.trimmed.tree$tip.label));
	#calculate sister divergence
	list.sp <- NULL;
	temp.sis.dist<-NULL;
	for(i in 1:length(temp.trimmed.tree$tip.label)){
		if(sum(temp.trimmed.tree$tip.label[i]==list.sp)==0){
		for(j in 1:length(temp.trimmed.tree$tip.label)){
			if(i != j){
				if(is.monophyletic(temp.trimmed.tree, c(temp.trimmed.tree$tip.label[i],temp.trimmed.tree$tip.label[j])) == TRUE){
					list.sp <- c(list.sp, temp.trimmed.tree$tip.label[j])
					dist <- cophenetic(temp.trimmed.tree)[temp.trimmed.tree$tip.label[i],temp.trimmed.tree$tip.label[j]]
					temp.sis.dist<-c(temp.sis.dist,dist)
				}
			}
		}
			
		}	
		}
		sis.dist<-c(sis.dist,temp.sis.dist)
}
temp2.sis.dist <- sis.dist;





par(mfrow = c(1,2));
boxplot(trop.sp.vec, temp.sp.vec, temp2.sp.vec, names = c('tropical', 'temperate', 'temperate\\nlow converting'), ylab = 'Number of Species',cex.lab=1.5);
boxplot(trop.sis.dist, temp1.sis.dist, temp2.sis.dist, names = c('tropical', 'temperate', 'temperate\\nlow converting'), ylab = 'Divergence Times between Sister Taxa', cex.lab = 1.5);