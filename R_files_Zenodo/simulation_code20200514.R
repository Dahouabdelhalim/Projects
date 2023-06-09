	library(gstat)
	library(sp)		
	library(vegan)
	library(ade4)
	library(adehabitatHS)
	library(entropart)
	library(parallel)
#generating 2 environmental variables
grid.gen<-function(lat.length,long.length)
   {
	 # Spherical model was used
	 long<-seq(10,(long.length-1)*20+10,by=20)
	 lat<-seq(10,(lat.length-1)*20+10,by=20)
	 xy <- expand.grid(lat,long)
     colnames(xy) <- c('x','y')
     # define the gstat object (spatial model)
     geo.model1<- gstat(formula=z~1, locations=~x+y, dummy=T, beta=0,
            model=vgm(psill=120,model='Sph',range=90,nugget=10))
     # time one simulations based on the gstat object
     geo.model2<- gstat(formula=z~1, locations=~x+y, dummy=T, beta=0,
            model=vgm(psill=90,model='Sph',range=60,nugget=5))

     sim.res1 <- predict(geo.model1, newdata=xy, nsim=1)
     sim.res2 <- predict(geo.model2, newdata=xy, nsim=1)

     GridEnv1<-matrix(sim.res1$sim1,lat.length,long.length,byrow=T)
	 GridEnv2<-matrix(sim.res2$sim1,lat.length,long.length,byrow=T);
	 #GridEnv1<-100*(GridEnv1-min(GridEnv1))/(max(GridEnv1)-min(GridEnv1))
	 GridEnv1<-(GridEnv1-min(GridEnv1))/(max(GridEnv1)-min(GridEnv1))
	 GridEnv2<-(GridEnv2-min(GridEnv2))/(max(GridEnv2)-min(GridEnv2))
	
	 return(list(GridEnv1=GridEnv1,GridEnv2=GridEnv2))
 	}
  
  ##environmental landscape
	long=20
	lat=20
	habitat<-grid.gen(lat.length=lat,long.length=long)
  #environmental variable 1
	GridEnvValues1<-as.vector(habitat$GridEnv1)
  #environmental variable 2
	GridEnvValues2<-as.vector(habitat$GridEnv2)
	env=data.frame(Env1=GridEnvValues1,Env2=GridEnvValues2)

#calculating observed beta-diversity and beta-diviation using four indices
beta.deviation = function(obs.data,simdata){
     #effective numbers of species (ENS)/Hill numbers
     multiple.beta=function(obs.data)
       {
		 #Partitions the diversity of a metacommunity into alpha and beta components
         mc=MetaCommunity(as.data.frame(t(obs.data)))#weights are equal by default
         dp=DivPart(q=1,mc,Biased=FALSE)#bias correction choice:"Best" is "ChaoWangJost".
         multiple.beta=dp$TotalBetaDiversity
         return(multiple.beta)
        }
  
     #beta diversity-sums of squared deviation from column means
     hellinger.betadiv = function(obs.data)
       {
	     hellinger.betadiv = sum(scale(decostand(obs.data, "hellinger"),center=TRUE, 
         scale=FALSE)^2)/(nrow(obs.data)-1)
		 return(hellinger.betadiv)
		}
  
     #calculating pairwise Bray-Curtise/percentage difference dissimilarity
     percentage.diff<- function(obs.data)
	   { 
         D = vegdist(obs.data, "bray")
	     D = sqrt(D)#sqrt.D=T
	     SStotal <- sum(D^2)/dim(obs.data)[1]     # eq. 8
	     BDtotal <- SStotal/(dim(obs.data)[1]-1)   # eq. 3
	     return(BDtotal)
	    }
  
     #calculating pairwise Jaccard-Chao dissimilarity
     jaccard_chao = function(obs.data)
       {
         D = vegdist(obs.data,method="chao")
	     SStotal <- sum(D^2)/dim(obs.data)[1]     # eq. 8
	     BDtotal <- SStotal/(dim(obs.data)[1]-1)   # eq. 3
	     return(BDtotal)
	    }
  
     #calculating observed, deviation and standardized beta diversity
     #multiple-site beta diversity using Hill numbers
     multiple.obs = multiple.beta(obs.data)
	 #parallel computation
	 mc <- getOption("mc.cores", 6)
     multiple.null = unlist(mclapply(simdata$perm, multiple.beta,mc.cores=mc))
     multiple.deviation = multiple.obs-mean(multiple.null)
     multiple.sesdev = multiple.deviation/sd(multiple.null)

     #calculating observed, deviation and standardized 
     #sums of squared deviation 
     hellinger.obs = hellinger.betadiv(obs.data)
     hellinger.null = unlist(mclapply(simdata$perm, hellinger.betadiv,mc.cores=mc))
     hellinger.deviation = hellinger.obs-mean(hellinger.null)
     hellinger.sesdev = hellinger.deviation/sd(hellinger.null)
  
     #calculating observed, deviation and standardized 
     #pairwise total beta diversity of Bray-Curtise/percentage difference dissimilarity 
     percentage.obs = percentage.diff(obs.data)
     percentage.null = unlist(mclapply(simdata$perm, percentage.diff,mc.cores=mc))
     percentage.deviation = percentage.obs-mean(percentage.null)
     percentage.sesdev = percentage.deviation/sd(percentage.null)

     #calculating observed, deviation and standardized 
     #pairwise total beta diversity of Jaccard-Chao dissimilarity   
     chao.obs=jaccard_chao(obs.data)
     chao.null = unlist(mclapply(simdata$perm, jaccard_chao,mc.cores=mc))
     chao.deviation = chao.obs-mean(chao.null)
     chao.sesdev = chao.deviation/sd(chao.null)
  
     return(list(multiple.obs = multiple.obs, multiple.deviation = multiple.deviation, multiple.sesdev = multiple.sesdev,
     hellinger.obs = hellinger.obs, hellinger.deviation = hellinger.deviation, hellinger.sesdev = hellinger.sesdev,
     percentage.obs = percentage.obs, percentage.deviation = percentage.deviation, percentage.sesdev = percentage.sesdev,
     chao.obs = chao.obs, chao.deviation = chao.deviation, chao.sesdev = chao.sesdev))
    }	

#calculate normalized divergence when q=1	
Diff.qiu1=function(abun){
     norm.beta=NULL;
	 abun=as.matrix(abun)
     N=nrow(abun); 
     D=numeric(N); 
     gabun=colSums(abun);gI=which(gabun>0);gT=sum(abun)
  	   for(i in 1:N){
    	 T=sum(abun[i,]);I=which(abun[i,]>0);
         D[i]=exp(-sum((abun[i,][I]/T)*log(abun[i,][I]/T)));     
    	}
  	 gD=exp(-sum((gabun[gI]/gT)*log(gabun[gI]/gT)));
     aI=which(abun>0); 
     aD=exp(-sum(abun[aI]/gT*log(abun[aI]/gT)))/N; 
     bD=gD/aD;
     CqN=1-log(bD)/log(N);
     UqN=CqN;
	 norm.beta=1-CqN;
	 return(norm.beta)
    }

#calculate normalized divergence when q is not 1
Diff.qiu2=function(abun,q){
	 abun=as.matrix(abun)
     N=nrow(abun); 
     D=numeric(N);
     gabun=colSums(abun);gI=which(gabun>0);gT=sum(abun) 
  	   for(i in 1:N) {
  		 T=sum(abun[i,]);I=which(abun[i,]>0);
      	 D[i]=sum((abun[i,][I]/T)^q)^(1/(1-q));  
   	    } 
  	 gD=sum((gabun[gI]/gT)^q)^(1/(1-q));    
     aI=which(abun>0); 
     aD=sum((abun[aI]/gT)^q)^(1/(1-q))/N
     bD=gD/aD;
     CqN=1-(bD^(1-q)-1)/(N^(1-q)-1);
     UqN=1-(bD^(q-1)-1)/(N^(q-1)-1); 
     norm.sor=1-CqN
	 norm.jac=1-UqN
	 return(data.frame(norm.sor,norm.jac))
    }  

#calculate observed value, raw and standardized beta-null deviation of normalized divergence
norm.beta=function(data,q,data.null){
     
	 if(q==1){
	     norm.obs=Diff.qiu1(data)
	     mc <- getOption("mc.cores", 6)
		 norm.null = unlist(mclapply(data.null$perm, Diff.qiu1,mc.cores=mc))
	     #norm.null = unlist(lapply(data.null$perm, Diff.qiu1))
	     norm.deviation=norm.obs-mean(norm.null)
	     norm.sesdev=norm.deviation/sd(norm.null)
	     return(data.frame(norm.obs,norm.deviation,norm.sesdev))
	    }
	 else{
		 norm.obs=Diff.qiu2(data,q)
		 norsor.obs=norm.obs$norm.sor
		 norjac.obs=norm.obs$norm.jac
		 mc <- getOption("mc.cores", 6)
		 norm.null=do.call(rbind,mclapply(data.null$perm, Diff.qiu2,q,mc.cores=mc))
		 #norm.null=do.call(rbind,lapply(data.null$perm, Diff.qiu2,q))
		 norsor.deviation=norsor.obs-mean(norm.null$norm.sor)
		 norjac.deviation=norjac.obs-mean(norm.null$norm.jac)
	     norsor.sesdev=norsor.deviation/sd(norm.null$norm.sor)
		 norjac.sesdev=norjac.deviation/sd(norm.null$norm.jac)
	     return(data.frame(norsor.obs,norsor.deviation,norsor.sesdev,norjac.obs,norjac.deviation,norjac.sesdev))
        }
	}

norm.deviation=function(data,data.null)
   {
     #q=1,horn
	 norm1=norm.beta(data,q=1,data.null)
	 #q=0,sorenson and jaccard
	 norm0=norm.beta(data,q=0,data.null)
	 colnames(norm0)=c("norsor.obs0","norsor.deviation0","norsor.sesdev0","norjac.obs0",
	  "norjac.deviation0","norjac.sesdev0")
	 #q=2,horn-morisita and regional non-overlap
	 norm2=norm.beta(data,q=2,data.null)
	 colnames(norm2)=c("norsor.obs2","norsor.deviation2","norsor.sesdev2","norjac.obs2",
	 "norjac.deviation2","norjac.sesdev2")
	 all_beta=cbind(norm0,norm1,norm2)
	 return(all_beta)
    }

#generating simulated community data
simulate.comm=function(S, ind, varm_para, shape)
   {
     varm=NULL
     m=NULL
     #survival rate of each species
     lambdaE<-function(Env1,Env2,m1,m2,varm) 
       {
		 (1/(sqrt(2*pi*(varm^2)))) * exp(- (((Env1-m1)^2)/(2*(varm^2))+((Env2-m2)^2)/(2*(varm^2))) )
		}
  
     #defining niche breadth 
     varm = varm_para
     #defining the gradient of niche position using beta distribution
     m = rbeta(S, shape1 = shape, shape2 = shape)
  
     lmat<- matrix(0, nrow = long*lat, ncol = S)
     # Calculating The fundamental Niche
     dimnames(lmat)<-list(1:nrow(lmat), paste('sp',1:S, sep=''))
	   for(i in 1:S) lmat[,i]<- 
	     #for two niche axes in this paper
	     lambdaE(Env1=GridEnvValues1, m1=m[i], Env2=GridEnvValues2, m2=m[i], varm=varm)
	     lmat[which(lmat<1e-301)]<-1e-301
	     cell<-1:(long*lat)
	     sim.comm<-t(sapply(cell, function (i) table(sample(as.factor(1:S), ind, prob=lmat[i,]/sum(lmat[i,]), replace=T))))
	      if(sum(apply(sim.comm, MARGIN=2, FUN=sum)>0)) sim.comm[, apply(sim.comm, MARGIN=2, FUN=sum)>0]->sim.comm
	     dimnames(sim.comm)<-list(1:nrow(sim.comm), paste('sp', 1:ncol(sim.comm), sep=''))
		 return(data.frame(sim.comm))
    }
 
#calculating beta indices based on simulated community data
simulation=function(S,ind,perms,replicates,gridsize,plotdim)
   {
	 stems=NULL
	 gamma=NULL
     breadth=NULL
	 margin=NULL
	 sim.comm=sim.comm2=NULL
	 sim.beta.deviation=NULL
	 multiple.obs=multiple.deviation=multiple.sesdev=NULL
     hellinger.obs=hellinger.deviation=hellinger.sesdev=NULL
     percentage.obs=percentage.deviation=percentage.sesdev=NULL
     chao.obs=chao.deviation=chao.sesdev=NULL
	 norsor.obs0=norsor.deviation0= norsor.sesdev0=norjac.obs0=norjac.deviation0=norjac.sesdev0=norm.obs=
	 norm.deviation=norm.sesdev=norsor.obs2=norsor.deviation2= norsor.sesdev2=norjac.obs2=norjac.deviation2=
	 norjac.sesdev2=NULL
	 shape=c(2,6,10)
	 reps=replicates
      for(l in 1:length(ind))
	   {
         for(j in 1:length(S))
		   { 
			 varm_para=c(0.1,0.3,0.5)
	          for(z in 1:length(varm_para))
			   {
	             for(v in 1:length(shape))
				   {
		             for(k in 1:reps)
			           {
			             #simulating community data based on species number, individual number, niche breadth and niche position
		                 sim.comm=simulate.comm(S[j], ind[l], varm_para[z], shape[v])
			             #simulating community data using null model with fixed rows and columns for 499 times
	                     sim.comm2 = permatfull(sim.comm, fixedmar="both", shuffle="ind", times=perms)
        	            
						 #calculating observed value, beta-null deviation and standardized deviation of beta diversity using 4 kinds of indices
			             sim.beta.deviation=beta.deviation(obs.data=sim.comm,simdata=sim.comm2)
			             multiple.obs<-c(multiple.obs, sim.beta.deviation$multiple.obs)
			             multiple.deviation<-c(multiple.deviation, sim.beta.deviation$multiple.deviation)
			             multiple.sesdev<-c(multiple.sesdev, sim.beta.deviation$multiple.sesdev)
                         hellinger.obs<-c(hellinger.obs, sim.beta.deviation$hellinger.obs)
			             hellinger.deviation<-c(hellinger.deviation, sim.beta.deviation$hellinger.deviation)
			             hellinger.sesdev<-c(hellinger.sesdev, sim.beta.deviation$hellinger.sesdev)
			             percentage.obs<-c(percentage.obs, sim.beta.deviation$percentage.obs)
			             percentage.deviation<-c(percentage.deviation, sim.beta.deviation$percentage.deviation)
			             percentage.sesdev<-c(percentage.sesdev, sim.beta.deviation$percentage.sesdev)
			             chao.obs<-c(chao.obs, sim.beta.deviation$chao.obs)
			             chao.deviation<-c(chao.deviation, sim.beta.deviation$chao.deviation)
			             chao.sesdev<-c(chao.sesdev, sim.beta.deviation$chao.sesdev)
						
			             #calculating observed value, beta-null deviation and standardized deviation of beta diversity using the normalized divergence index
			             sim.norm.deviation=norm.deviation(data=sim.comm,data.null=sim.comm2)
			             norsor.obs0<-c(norsor.obs0, sim.norm.deviation$norsor.obs0)
			             norsor.deviation0<-c(norsor.deviation0, sim.norm.deviation$norsor.deviation0)
			             norsor.sesdev0<-c(norsor.sesdev0, sim.norm.deviation$norsor.sesdev0)
                         norjac.obs0<-c(norjac.obs0, sim.norm.deviation$norjac.obs0)
			             norjac.deviation0<-c(norjac.deviation0, sim.norm.deviation$norjac.deviation0)
			             norjac.sesdev0<-c(norjac.sesdev0, sim.norm.deviation$norjac.sesdev0)
			             norm.obs<-c(norm.obs, sim.norm.deviation$norm.obs)
			             norm.deviation<-c(norm.deviation, sim.norm.deviation$norm.deviation)
			             norm.sesdev<-c(norm.sesdev, sim.norm.deviation$norm.sesdev)
			             norsor.obs2<-c(norsor.obs2, sim.norm.deviation$norsor.obs2)
			             norsor.deviation2<-c(norsor.deviation2, sim.norm.deviation$norsor.deviation2)
			             norsor.sesdev2<-c(norsor.sesdev2, sim.norm.deviation$norsor.sesdev2)
                         norjac.obs2<-c(norjac.obs2, sim.norm.deviation$norjac.obs2)
			             norjac.deviation2<-c(norjac.deviation2, sim.norm.deviation$norjac.deviation2)
			             norjac.sesdev2<-c(norjac.sesdev2, sim.norm.deviation$norjac.sesdev2)
			 
	                     gamma<-c(gamma, S[j])
                         stems<-c(stems, ind[l])
			             breadth = c(breadth, varm_para[z])
	                     margin=c(margin, shape[v]) 
		                }
                    }
                }
            }
        }
  
     return(data.frame(gamma, stems, breadth, margin, multiple.obs, multiple.deviation, multiple.sesdev,
	 hellinger.obs,hellinger.deviation, hellinger.sesdev, percentage.obs, percentage.deviation, 
	 percentage.sesdev, chao.obs, chao.deviation, chao.sesdev,norsor.obs0,norsor.deviation0, 
	 norsor.sesdev0, norjac.obs0,norjac.deviation0,norjac.sesdev0,norm.obs,norm.deviation,norm.sesdev,
	 norsor.obs2,norsor.deviation2, norsor.sesdev2,norjac.obs2,norjac.deviation2, norjac.sesdev2))
    }
#calculating beta-diversity metrics with simulating metacommunities 
sim_data=simulation(S=c(50,100,150,200,300,400),ind=c(50,100,150,200,250,300),perms=999,replicates=50,gridsize=20,plotdim=c(400,400))
save.image("sim_data.RData")	
	

	
	


