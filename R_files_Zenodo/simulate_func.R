## Simulation code by E.G. Mandeville, 2018
## This is a function for individual-based simulations of genomic outcomes of hybridization. Ancestry of individuals is tracked using the same statistics generated from empirical data using entropy. For each individual, we track q, proportion of ancestry, and Q, interspecific ancestry.

simulate.hyb <- function(nind.start = 100, prop.sp1=0.5, n.generation=10, makeplot=TRUE, printoutput=TRUE, imm.sp1=0, imm.sp2=0){

    ## Number of immigrants in each generation. Default is 0.

    nind.imm <- round(nind.start * imm.sp1) + round(nind.start * imm.sp2)
    nind.imm.sp1 <- round(nind.start * imm.sp1)
    nind.imm.sp2 <- round(nind.start * imm.sp2)

    ## Set up matrix for keeping track of q, or proportion of ancestry. Note that in this simulation, q is calculated from Q values for each individual
    qfromQ.allgen <-matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    qfromQ.allgen[1,] <- c(rep(0, round(prop.sp1*nind.start)), rep(1, (nind.start-round(prop.sp1*nind.start))))

    ## Set up empty matrices for keeping track of Q in all 4 categories
    Q11.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    Q11.allgen[1,] <- c(rep(1,round(prop.sp1*nind.start, digits=0)), rep(0,(nind.start-round(prop.sp1*nind.start, digits=0))))

    Q12.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    Q12.allgen[1,] <- c(rep(0, nind.start))

    Q21.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    Q21.allgen[1,] <- c(rep(0, nind.start))

    Q22.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    Q22.allgen[1,] <- c(rep(0,round(prop.sp1*nind.start, digits=0)), rep(1,(nind.start-round(prop.sp1*nind.start, digits=0))))

    Qinter.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    Qinter.allgen[1,] <- rep(0, nind.start)

    if(makeplot == TRUE){
        pdf(paste("qQ_",n.generation,"gen_",prop.sp1,"sp1_",nind.start,"ind_", nind.imm.sp1, "immsp1_", nind.imm.sp2, "immsp2.pdf",sep=""))
##qfromQ <- Q11.allgen+((Q12.allgen+Q21.allgen)/2)
        par(mfrow=c(ceiling(n.generation/4),4))
    }

    for(i in 2:(n.generation+1)){


    ## sample a vector of integers 1:nind, then use that vector to index q.allgen and Q.allgen so we can get matching q and Q
        ## sampling WITH REPLACEMENT, i.e. the fact that an individual has already produced one offspring does not make it ineligible to produce another
        ## Keep population constant even in the face of immigration
        randvec1 <- sample(1:(nind.start+nind.imm), nind.start, replace=T)
        randvec2 <- sample(1:(nind.start+nind.imm), nind.start, replace=T)

        ## sample from prev. generation for parents of new offspring, assume random mating
        ## use random vectors to index both q and Q
        ## Updated Nov. 5, 2018 to include immigration

        parent1.Q11 <- c(Q11.allgen[i-1,], rep(1,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec1]
        parent1.Q12 <- c(Q12.allgen[i-1,], rep(0,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec1]
        parent1.Q21 <- c(Q21.allgen[i-1,], rep(0,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec1]
        parent1.Q22 <- c(Q22.allgen[i-1,], rep(0,nind.imm.sp1), rep(1,nind.imm.sp2))[randvec1]

        parent2.Q11 <- c(Q11.allgen[i-1,], rep(1,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec2]
        parent2.Q12 <- c(Q12.allgen[i-1,], rep(0,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec2]
        parent2.Q21 <- c(Q21.allgen[i-1,], rep(0,nind.imm.sp1), rep(0,nind.imm.sp2))[randvec2]
        parent2.Q22 <- c(Q22.allgen[i-1,], rep(0,nind.imm.sp1), rep(1,nind.imm.sp2))[randvec2]


        Q11.allgen[i,] <- (parent1.Q11+parent1.Q12)*(parent2.Q11+parent2.Q12)
        Q12.allgen[i,] <- (((parent1.Q11+parent1.Q12)*(parent2.Q21+parent2.Q22))+((parent1.Q21+parent1.Q22)*(parent2.Q11+parent2.Q12)))/2 ## Need to do this to retain symmetry of Q12,Q21
        Q21.allgen[i,] <- Q12.allgen[i,]
        Q22.allgen[i,] <- (parent1.Q21+parent1.Q22)*(parent2.Q21+parent2.Q22)

        qfromQ.allgen[i,] <- Q11.allgen[i,]+((Q12.allgen[i,]+Q21.allgen[i,])/2)

        Qinter.allgen[i,] <- Q12.allgen[i,]+Q21.allgen[i,]

        if(makeplot == TRUE){
            plot(qfromQ.allgen[i,], Qinter.allgen[i,], main=paste("Gen.", i-1, sep=" "), xlab="q", ylab="Q", xlim=c(0,1), ylim=c(0,1), type="n", axes=F)
            axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
            axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
            arrows(0,0,0.5,1, length=0, col="gray90")
            arrows(0.5,1,1,0, length=0, col="gray90")
            points(qfromQ.allgen[i,], Qinter.allgen[i,], col="gray45")
        }
    }
    if(makeplot == TRUE){
        dev.off()

    }

    if(printoutput == TRUE){
        return(list(qfromQ.allgen, Qinter.allgen))
    }
}

