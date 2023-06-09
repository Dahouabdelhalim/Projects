mkdata<-function(theta=NA,ploidy=NA,fis=NA,fst=NA,depth=NA,N=200,snps=10000,reps=10){
    sout<-vector("list",reps)    
    for(x in 1:reps){

        ## sample mean allele frequncies
        pi<-rbeta(snps,shape1=theta,shape2=theta)
        ## generate two gene pools based on Fst
        if(fst==0){
            p1<-pi
            p2<-pi
        }
        else{
            S<-(1-fst)/fst ## check this
            p1<-rbeta(snps,shape1=pi*S,shape2=(1-pi)*S)
            p2<-rbeta(snps,shape1=pi*S,shape2=(1-pi)*S)
        }    
        n1<-N/2
        ## sample ploidies based on ploidy vector
        pd<-sample(c(2,3,4),N,replace=TRUE,prob=ploidy)
            
        ## make individuals
        G<-matrix(NA,nrow=N,ncol=snps)
        for(j in 1:N){
            ## grap appropriate allele frequencies
            if(j<=n1){p<-p1}
            else{p<-p2}
            if(fis==0){ ## can sample all at once when fis==0
                g<-rbinom(n=snps,size=pd[j],prob=p)
            }
            else{
                g<-rep(NA,snps)
                for(i in 1:snps){
                    ## sample inbreeding state, multiple for trips and tetras
                    f<-sample(c(0,1),(pd[j]-1),replace=TRUE,prob=c((1-fis),fis))
                    ## diploids
                    if(pd[j]==2){
                        if(f==0){ ## not inbred
                            g[i]<-rbinom(n=1,size=2,prob=p[i])
                        }
                        else{
                            g[i]<-2 * rbinom(n=1,size=1,prob=p[i])
                        }
                    }      
                    ## triploids
                    if(pd[j]==3){
                        if(f[1]==0){ ## not inbred
                            g[i]<-rbinom(n=1,size=3,prob=p[i])
                        }
                        else if (f[1]==1 & f[2]==1){## doubly inbred
                            g[i]<-3 * rbinom(n=1,size=1,prob=p[i])
                        }
                        else{ ## only unreduced gamete inbred
                            g[i]<-2 * rbinom(n=1,size=1,prob=p[i]) + rbinom(n=1,size=1,prob=p[i])
                        }
                    }  
                    ## tetraploids
                    if(pd[j]==4){
                        if(f[1]==0 & f[2]==0){ ## not inbred
                            g[i]<-rbinom(n=1,size=4,prob=p[i])
                        }
                        else if (sum(f)==3){## triply inbred
                            g[i]<-4 * rbinom(n=1,size=1,prob=p[i])
                        }
                        else if (f[1]==1 & f[2]==1){ ## doubly inbred
                            g[i]<-2 * rbinom(n=1,size=1,prob=p[i]) + 2* rbinom(n=1,size=1,prob=p[i])
                        }
                        else { # singly inbred
                            g[i]<-2 * rbinom(n=1,size=1,prob=p[i]) + rbinom(n=1,size=2,prob=p[i])
                        }    
                    }    
                }
            }
            G[j,]<-g          
        }
        ## finished individual loop, now sample sequence
        a1<-matrix(NA,nrow=N,ncol=snps)
        a2<-a1
        for(j in 1:N){
            ## sample number of reads per SNP
            nrds<-rpois(n=snps,lambda=depth)
            ## generate read data
            rds<-rbinom(n=snps,size=nrds,prob=G[j,]/pd[j])
            hets<-which(rds >= 1 & rds < nrds)
            a1[j,hets]<-rds[hets]
            a2[j,hets]<-nrds[hets]-rds[hets]
        }     
        ## return results
        sout[[x]]<-list(a1,a2,pd)
    }   
    sout     
}

############### runs simulaitons ###############
## base conditions, high diversity, no structure or inbreeding, even sampling of ploidy and high coverage
outDipTrip<-mkdata(theta=0.8,ploidy=c(0.5,0.5,0.0),fis=0,fst=0,depth=8)
outDipTetra<-mkdata(theta=0.8,ploidy=c(0.5,0.0,0.5),fis=0,fst=0,depth=8)
outDipTripTetra<-mkdata(theta=0.8,ploidy=c(0.33,0.33,0.33),fis=0,fst=0,depth=8)

## lower diversity
outDipTripLdiv<-mkdata(theta=0.1,ploidy=c(0.5,0.5,0.0),fis=0,fst=0,depth=8)
outDipTetraLdiv<-mkdata(theta=0.1,ploidy=c(0.5,0.0,0.5),fis=0,fst=0,depth=8)
outDipTripTetraLdiv<-mkdata(theta=0.1,ploidy=c(0.33,0.33,0.33),fis=0,fst=0,depth=8)

## population stucture
outDipTripStr<-mkdata(theta=0.8,ploidy=c(0.5,0.5,0.0),fis=0,fst=0.3,depth=8)
outDipTetraStr<-mkdata(theta=0.8,ploidy=c(0.5,0.0,0.5),fis=0,fst=0.3,depth=8)
outDipTripTetraStr<-mkdata(theta=0.8,ploidy=c(0.33,0.33,0.33),fis=0,fst=0.3,depth=8)

## inbreeding
outDipTripIn<-mkdata(theta=0.8,ploidy=c(0.5,0.5,0.0),fis=0.2,fst=0,depth=8)
outDipTetraIn<-mkdata(theta=0.8,ploidy=c(0.5,0.0,0.5),fis=0.2,fst=0,depth=8)
outDipTripTetraIn<-mkdata(theta=0.8,ploidy=c(0.33,0.33,0.33),fis=0.2,fst=0,depth=8)

## uneven sampling
outDipTripExd<-mkdata(theta=0.8,ploidy=c(0.9,0.1,0.0),fis=0,fst=0,depth=8)
outDipTetraExd<-mkdata(theta=0.8,ploidy=c(0.9,0.0,0.1),fis=0,fst=0,depth=8)
outDipTripExt<-mkdata(theta=0.8,ploidy=c(0.1,0.9,0.0),fis=0,fst=0,depth=8)
outDipTetraExt<-mkdata(theta=0.8,ploidy=c(0.1,0.0,0.9),fis=0,fst=0,depth=8)

## low coverage
outDipTripLcov<-mkdata(theta=0.8,ploidy=c(0.5,0.5,0.0),fis=0,fst=0,depth=2)
outDipTetraLcov<-mkdata(theta=0.8,ploidy=c(0.5,0.0,0.5),fis=0,fst=0,depth=2)
outDipTripTetraLcov<-mkdata(theta=0.8,ploidy=c(0.33,0.33,0.33),fis=0,fst=0,depth=2)

## save workspace
save(list=ls(),file="sims.Rdata")

## example analysis
source("ploidy.R")

## version for allotetraploids, here Fst is between species
mkdata2<-function(theta=NA,ploidy=NA,propallo=1,fst=NA,depth=NA,N=200,snps=10000,reps=10){
    sout<-vector("list",reps)    
    for(x in 1:reps){

        ## sample mean allele frequncies
        pi<-rbeta(snps,shape1=theta,shape2=theta)
        ## generate two gene pools based on Fst
        if(fst==0){
            p1<-pi
            p2<-pi
        }
        else{
            S<-(1-fst)/fst ## check this
            p1<-rbeta(snps,shape1=pi*S,shape2=(1-pi)*S)
            p2<-rbeta(snps,shape1=pi*S,shape2=(1-pi)*S)
	
        }
	px<-list(p1,p2)    
        ## sample ploidies based on ploidy vector
        pd<-sample(c(2,3,4),N,replace=TRUE,prob=ploidy)
	## sample whether tetraploids are allotetraploids (sample ignored for non-tetraploids)
	allo<-sample(c(0,1),N,replace=TRUE,prob=c(1-propallo,propallo))            
        ## make individuals
        G<-matrix(NA,nrow=N,ncol=snps)
        pops<-rep(NA,N)
        for(j in 1:N){
        	g<-rep(NA,snps)
		pop<-sample(c(1,2),1,prob=c(0.5,0.5)) ## population for non-allos
                for(i in 1:snps){
                    ## diploids
		    if(pd[j]==2){
			g[i]<-rbinom(n=1,size=2,prob=px[[pop]][i])
                    }      
                    ## triploids
                    if(pd[j]==3){
                        g[i]<-rbinom(n=1,size=3,prob=px[[pop]][i])
                    }  
                    ## tetraploids
                    if(pd[j]==4){
			if(allo[j]==0){ ## autotetraploid
                            g[i]<-rbinom(n=1,size=4,prob=px[[pop]][i])
			}
			else{ ## allotet
			    g[i]<-rbinom(n=1,size=2,prob=px[[1]][i])+rbinom(n=1,size=2,prob=px[[2]][i])
			}
                    }    
                }
	        pops[j]<-pop	
                G[j,]<-g          
        }
        ## finished individual loop, now sample sequence
        a1<-matrix(NA,nrow=N,ncol=snps)
        a2<-a1
        for(j in 1:N){
            ## sample number of reads per SNP
            nrds<-rpois(n=snps,lambda=depth)
            ## generate read data
            rds<-rbinom(n=snps,size=nrds,prob=G[j,]/pd[j])
            hets<-which(rds >= 1 & rds < nrds)
            a1[j,hets]<-rds[hets]
            a2[j,hets]<-nrds[hets]-rds[hets]
        }     
        ## return results
        sout[[x]]<-list(a1,a2,pd,allo,pops)
    }   
    sout     
}

############### runs simulaitons ###############
## base conditions, high diversity, no structure or inbreeding, even sampling of ploidy and high coverage
outAlloTetra<-mkdata2(theta=0.8,ploidy=c(0.5,0.0,0.5),propallo=1,fst=0.8,depth=8)
outAlloAutoTetra<-mkdata2(theta=0.8,ploidy=c(0.5,0.0,0.5),propallo=0.5,fst=0.8,depth=8)

## save workspace
save(list=ls(),file="simsAlloTet.Rdata")

