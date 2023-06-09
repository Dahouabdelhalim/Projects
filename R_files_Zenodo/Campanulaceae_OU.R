###Campanulaceae OU floral trait analyses
### Lagomarsino et al 2017


 require(ape);require(picante);require(car);require(plotrix);require(nlme);require(geiger);require(abind);require(OUwie);require(doMC)
 
 setwd("~/Desktop/For new OU analyses")
 
 ###read in trait data
 ltr1 <- as.data.frame(read.csv("FloralData_Nov2016.csv"))
 ltr1$Species <- as.character(ltr1$Species)
 ltr1$Species <- gsub(" ","",ltr1$Species)
 trait_list <- colnames(ltr1)[2:length(colnames(ltr1))]
 
 ###reading in discrete traits used for reconstructions; should match reconstructed states in nexus trees
 ltr2=read.csv('qualtraits.csv')
 spp_wdata <- match(ltr1$Species,ltr2$Species) #matching datasets and only including those for which there is trait data; be sure to comment on this in methods - that recontructions were conducted on the whole tree, but continuous trait data was only available for a subset of the species 
ltr <- cbind(ltr1,ltr2[spp_wdata,2:4])
 
 ###reading in the phylogeny for OU analysis; start with fruit reconstructions
 phy2 <- read.nexus('fruit_recon.nex')

 tip.drop <- function(phy2){
      
        for(i in 1:length(phy2)){
        	
        	tips2drop <- which((!phy2[[1]]$tip.label%in%ltr$species_list))
        	phy2[[i]] <- drop.tip(phy2[[i]],tips2drop)
        	
        	}
 			
 		return(phy2)
    }

 ###function to get model fits for traits and set of trees

OU_fit <- function(trait,phy,model=c("BM","OU1","OU2","OU3"),tree_num=c(1:100),scenario_num=c(1:11)){

###setting up temp data frame for the scenario, tree, etc.

	temp_dat <- dat
	temp_dat[1,'scenario'] <- scenario_num
	temp_dat[1,'tree_num'] <- tree_num
	temp_dat[1,'model'] <- model
	
	
###looping through all traits to retrieve likelihood values, AIC.c's and thetas for the models and trees	

####Brownian single state model

	if(model=="BM"){
		dat_temp <- list(trait[,c('Species','Regime_2',names(trait[2]))])
			for(i in 3:17){
				dat_temp1<- list(trait[,c('Species','Regime_2',names(trait[i]))])
				dat_temp <- c(dat_temp,dat_temp1)
				}
					model_OUwie <- lapply(dat_temp,OUwie,phy=phy,model="BM1")
					for(j in 1:length(model_OUwie)){
						temp_dat[1,'trait']<-trait_list[j]
						temp_dat[1,'aic_c']<-model_OUwie[[j]]$AICc
						temp_dat[1,'lik']<-model_OUwie[[j]]$loglik
						write.table(temp_dat,file="Camp_summary_fruit.txt",append=TRUE,col.names=FALSE)
						}
			
			}else{}
							
	
	
###OU single state model	

	if(model=="OU1"){
		dat_temp <- list(trait[,c('Species','Regime_2',names(trait[2]))])
			for(i in 3:17){
				dat_temp1<- list(trait[,c('Species','Regime_2',names(trait[i]))])
				dat_temp <- c(dat_temp,dat_temp1)
				}
					model_OUwie <- lapply(dat_temp,OUwie,phy=phy,model="OU1",ub=10000000)
					for(j in 1:length(model_OUwie)){
						temp_dat[1,'trait']<-trait_list[j]
						temp_dat[1,'aic_c']<-model_OUwie[[j]]$AICc
						temp_dat[1,'lik']<-model_OUwie[[j]]$loglik
						temp_dat[1,'theta1']<-model_OUwie[[j]]$theta[1,1]
						temp_dat[1,'theta1_se'] <- model_OUwie[[j]]$theta[1,2]
						
						write.table(temp_dat,file="Camp_summary_fruit.txt",append=TRUE,col.names=FALSE)
						}
			
			}else{}
							
	
	
### OU 2 state model

	if(model=="OU2"){
		dat_temp <- list(trait[,c('Species','Regime_2',names(trait[2]))])
			for(i in 3:17){
				dat_temp1<- list(trait[,c('Species','Regime_2',names(trait[i]))])
				dat_temp <- c(dat_temp,dat_temp1)
				}
					model_OUwie <- lapply(dat_temp,OUwie,phy=phy,model="OUM")
					for(j in 1:length(model_OUwie)){
						temp_dat[1,'trait']<-trait_list[j]
						temp_dat[1,'aic_c']<-model_OUwie[[j]]$AICc
						temp_dat[1,'lik']<-model_OUwie[[j]]$loglik
						temp_dat[1,'theta1']<-model_OUwie[[j]]$theta[1,1]
						temp_dat[1,'theta1_se'] <- model_OUwie[[j]]$theta[1,2]
						temp_dat[1,'theta2']<-model_OUwie[[j]]$theta[2,1]
						temp_dat[1,'theta2_se'] <- model_OUwie[[j]]$theta[2,2]
						write.table(temp_dat,file="Camp_summary_fruit.txt",append=TRUE,col.names=FALSE)
						}
			
			}else{}
							


### OU 3 state model

	if(model=="OU3"){
		dat_temp <- list(trait[,c('Species','Regime_3',names(trait[5]))])
			for(i in 3:17){
				dat_temp1<- list(trait[,c('Species','Regime_3',names(trait[i]))])
				dat_temp <- c(dat_temp,dat_temp1)
				}
					model_OUwie <- lapply(dat_temp,OUwie,phy=phy,model="OUM")
					for(j in 1:length(model_OUwie)){
						temp_dat[1,'trait']<-trait_list[j]
						temp_dat[1,'aic_c']<-model_OUwie[[j]]$AICc
						temp_dat[1,'lik']<-model_OUwie[[j]]$loglik
						temp_dat[1,'theta1']<-model_OUwie[[j]]$theta[1,1]
						temp_dat[1,'theta1_se'] <- model_OUwie[[j]]$theta[1,2]
						temp_dat[1,'theta2']<-model_OUwie[[j]]$theta[2,1]
						temp_dat[1,'theta2_se'] <- model_OUwie[[j]]$theta[2,2]
						temp_dat[1,'theta3']<-model_OUwie[[j]]$theta[3,1]
						temp_dat[1,'theta3_se'] <- model_OUwie[[j]]$theta[3,2]
						write.table(temp_dat,file="Camp_summary_fruit.txt",append=TRUE,col.names=FALSE)
						}
			
			}else{}

							
	}

###create a vector for storing/writing a file; start with fruit reconstructions

dat <- as.data.frame(matrix(nrow=0,ncol=12))
names(dat) <- c("scenario","tree_num","model","trait","aic_c","lik","theta1","theta1_se","theta2","theta2_se","theta3","theta3_se")
write.table(dat,file="Camp_summary_fruit.txt")
dat <- as.data.frame(matrix(nrow=1,ncol=12))
names(dat) <- c("scenario","tree_num","model","trait","aic_c","lik","theta1","theta1_se","theta2","theta2_se","theta3","theta3_se")


scenarios<-c(1:3)
model <- c("BM","OU1","OU2")

###modifying trait data for specific analysis
colnames(ltr)[which(colnames(ltr)=='Fruit')] <- 'Regime_2'

for(i in 1:length(scenarios)){
	
	#string1 <- paste(scenarios[i],"fruit_tree.tre",sep="_")
	#phy2<- read.tree(string1)
	
		
	for(j in 1:100){
		
			phy_temp <- phy2[[j]]
			print(paste(scenarios[i],j,sep="_"))
			tips2drop <- which((!phy_temp$tip.label%in%ltr$Species))
        	phy_temp <- drop.tip(phy_temp,tips2drop)
			OU_fit(ltr,phy_temp,model[i],j,scenarios[i])
							
		}
	}
	
		

