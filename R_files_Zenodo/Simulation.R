#By Yi Zou
#Last modified 2019-Oct

#Loading library and "ESS" function
Pkg <- c("vegan", "CommEcol", "plyr")
if (length(setdiff(Pkg , rownames(installed.packages()))) > 0) 
    {install.packages(setdiff(Pkg , rownames(installed.packages())))}

library(vegan) #Version 2.5.3 in the simulation
library(CommEcol) #Version 1.6.5
library(plyr) #Version 1.8.4
source("ESS.R")

#1. Data generation
#Generate the control dataset, ranked from the most "rare" to the most "dominant" species as "sp1" to "sp100"
set.seed(1)
x <- rlnorm(100,meanlog=6.5,sdlog=1)#Log distribution of the abundance
x1 <- round(sort(x),0)
names(x1) <- paste0("sp",1:length(x))

#Generate dataset sharing 25% of the rare species with x
set.seed(2)
x3 <- round(sort(rlnorm(100,meanlog=6.5,sdlog=1)),0)
set.seed(2)
namex3_1_25 <- sample(paste0("sp",1:25))
names(x3) <- c(namex3_1_25,paste0("sp",226:300))

#Generate dataset, sharing 25% of the dominant species with x
set.seed(4)
x4 <- round(sort(rlnorm(100,meanlog=6.5,sdlog=1)),0)
set.seed(4)
namex4_1_25 <- sample(paste0("sp",76:100))
names(x4) <- c(paste0("sp",326:400),namex4_1_25)

#Generate dataset, sharing 50% of the rare species
set.seed(5)
x5 <- round(sort(rlnorm(100,meanlog=6.5,sdlog=1)),0)
sum(x5)
set.seed(5)
namex5_1_25 <- sample(paste0("sp",1:50))
names(x5) <- c(namex5_1_25,paste0("sp",451:500))

#Generate dataset, sharing 50% of the dominant species
set.seed(6)
x6 <- round(sort(rlnorm(100,meanlog=6.5,sdlog=1)),0)
set.seed(6)
namex6_1_25 <- sample(paste0("sp",51:100))
names(x6) <- c(paste0("sp",551:600),namex6_1_25)

#Generate dataset, sharing 75% of the rare species
set.seed(7)
x7 <- round(sort(rlnorm(100,meanlog=6.5,sdlog=1)),0)
set.seed(7)
namex7_1_25 <- sample(paste0("sp",1:75))
names(x7) <- c(namex7_1_25,paste0("sp",676:700))

#Generate dataset, sharing 75% of the dominant species
set.seed(8)
x8 <- round(sort(rlnorm(100,meanlog=6.5,sdlog=1)),0)
set.seed(8)
namex8_1_25 <- sample(paste0("sp",26:100))
names(x8) <- c(paste0("sp",776:800),namex8_1_25)

#Making a list of the dataset
Lsp <- list(x1,x3,x4,x5,x6,x7,x8)

#Combine as data.frame
Ct <- t(as.matrix(x1))
R25 <- t(as.matrix(x3))
D25 <- t(as.matrix(x4))
R50 <- t(as.matrix(x5))
D50 <- t(as.matrix(x6))
R75 <- t(as.matrix(x7))
D75 <- t(as.matrix(x8))

Group <- c("Control","R25","R50","R75","D25","D50","D75")
Df <- data.frame(Group,rbind.fill.matrix(Ct,R25,R50,R75,D25,D50,D75))
Df[is.na(Df)] <- 0

#2 Calculate CNESS values for different sample size parameters ‘m’ from 1 to 100000 with exponential increase
Treat <- unique(Df$Group)
Result <- data.frame() #To store the results
for (j in 2:length(Treat))
{ 
    Treat_j=Treat[j]
    subdf <- rbind(Df[1,],Df[Df$Group==Treat_j,]) 
    for (e in seq(from=0,to=5,by=0.1))#The smoother curve will require longer calculation time.
    {
        m=round(10^e,0)
        Dst <- ESS(subdf[,-1],m=m)
        Out <- data.frame(X=Treat[j],m=m,Dst=Dst[1])
        Result <- rbind(Out, Result)
    }
}


#3. Calculate dissimilarity indices with the change of sampling size 
#Note, calculation time for the following simulation can easily exceed 5 hours (Processor: 2.6 GHz, Intel Core i7), while the project can be split between processor cores to accelerate the process

Lsp1 <- Lsp[[1]]
sp1 <- names(Lsp1)
Pool1 <- rep(sp1,Lsp1)
Sum1 <- sum(Lsp1)

#3.1 strategy1: equal sampling coverage
Out <- vector() 
for (j in c(2:7)) 
{
    Lspj <- Lsp[[j]]
    spj <- names(Lspj )
    Poolj <- rep(spj,Lspj)
    Sumj <- sum(Lspj)
    
    for (N in c(0.01,0.1,1,10,100))
    {   
        
        Nitt=1000
        Mean1 <- vector() 
        Mean2 <- vector()
        Mean3 <- vector()
        Mean4 <- vector()
        MeanBC <- vector()
        MeanEuc <- vector()
        MeanChao <- vector() 
        for (k in 1:Nitt) 
        {
            N1 <- floor( Sum1*N/100) #Equal sampling strategy between control (S1) and treatment (S2)
            Nj <- floor( Sumj*N/100)
            Xsample <- sample(Pool1,N1,replace=F)#Sample from X1 pool
            S1 <- t(as.matrix(table(Xsample)))
            Xsample2 <-  sample(Poolj,Nj,replace=F)#Sample from X2 pool
            S2 <- t(as.matrix(table(Xsample2))) 
            S1S2 <-  rbind.fill.matrix(S1,S2)
            S1S2 [is.na(S1S2)] <- 0
            
            for (m in c(1,10,100,1000))
            {
                if (m>Nj){next}
                if(m==1){Mean1[k] <- ESS(S1S2,m=m)}
                if(m==10){Mean2[k] <- ESS(S1S2,m=m)}
                if(m==100){Mean3[k] <- ESS(S1S2,m=m)}
                if(m==1000){Mean4[k] <- ESS(S1S2,m=m)}
            }
            MeanBC[k]<- vegdist (S1S2,method="bray",binary=F)# Bray-Curtis index
            MeanEuc[k]<- vegdist (decostand(S1S2, method = "total"),method="euclidean")#Proportion based Euclidian distance 
            MeanChao[k] <- dis.chao(S1S2,index="sorensen")#Chao-Sorense index
        }
        Distmean1 <- mean(Mean1)
        Bounds1 <- quantile(Mean1,probs = c(0.025, 0.975))
        Lower1 <- Bounds1[1]; Upper1 <- Bounds1[2]
        CV1 <- sd(Mean1)/Distmean1
        Df1 <-  c(N=N,X=j,m=1,Dst=Distmean1,Lo=Lower1,Up=Upper1,CV=CV1)
        
        Distmean2 <- mean(Mean2)
        Bounds2 <- quantile(Mean2,probs = c(0.025, 0.975))
        Lower2 <- Bounds2[1]; Upper2 <- Bounds2[2]
        CV2 <- sd(Mean2)/Distmean2
        Df2 <-  c(N=N,X=j,m=10,Dst=Distmean2,Lo2=Lower2,Up=Upper2,CV=CV2)
        
        Distmean3 <- mean(Mean3)
        Bounds3 <- quantile(Mean3,probs = c(0.025, 0.975))
        Lower3 <- Bounds3[1]; Upper3 <- Bounds3[2]
        CV3 <- sd(Mean3)/Distmean3
        Df3 <-  c(N=N,X=j,m=100,Dst=Distmean3,Lo=Lower3,Up=Upper3,CV=CV3)
        
        Distmean4 <- mean(Mean4)
        Bounds4 <- quantile(Mean4,probs = c(0.025, 0.975))
        Lower4 <- Bounds4[1]; Upper4 <- Bounds4[2]
        CV4 <- sd(Mean4)/Distmean4
        Df4 <-  c(N=N,X=j,m=1000,Dst=Distmean4,Lo=Lower4,Up=Upper4,CV=CV4)
        
        DistBC <- mean(MeanBC)
        BoudnsBC <- quantile(MeanBC,probs = c(0.025, 0.975))
        LowerBC <- BoudnsBC[1]; UpperBC <- BoudnsBC[2]
        CV5 <- sd(MeanBC)/DistBC
        Df5 <-  c(N=N,X=j,m="BC",Dst=DistBC,Lo=LowerBC,Up=UpperBC,CV=CV5)
        
        DistEuc <- mean(MeanEuc)
        BoudnsEuc <- quantile(MeanEuc,probs = c(0.025, 0.975))
        LowerEuc <- BoudnsEuc[1]; UpperEuc <- BoudnsEuc[2]
        CV6 <- sd(MeanEuc)/DistEuc
        Df6 <-  c(N=N,X=j,m="Euc",Dst=DistEuc,Lo=LowerEuc,Up=UpperEuc,CV=CV6)
        
        DistChao <- mean(MeanChao)
        BoudnsChao <- quantile(MeanChao,probs = c(0.025, 0.975))
        LowerChao <- BoudnsChao[1]; UpperChao <- BoudnsChao[2]
        CV7 <- sd(MeanChao)/DistChao
        Df7 <-  c(N=N,X=j,m="SoAbu",Dst=DistChao,Lo=LowerChao,Up=UpperChao,CV=CV7)
        
        Df <- rbind(Df1,Df2,Df3,Df4,Df5,Df6,Df7)
        Out <- rbind(Out,Df)
    }
}

rownames(Out) <- 1:nrow(Out)


#3.2 strategy2: unequal sampling coverage
Out2 <- vector() 
for (j in c(2:7)) 
{
    Lspj <- Lsp[[j]]
    spj <- names(Lspj )
    Poolj <- rep(spj,Lspj)
    Sumj <- sum(Lspj)
    
    for (N in c(0.01,0.1,1,10,100))#Setting different N according sampling completeness
    {   
        Nitt=1000
        Mean1 <- vector() 
        Mean2 <- vector()
        Mean3 <- vector()
        Mean4 <- vector()
        MeanBC <- vector()
        MeanEuc <- vector()
        MeanChao <- vector() 
        for (k in 1:Nitt) 
        {
            N1 <- 1000  #Unequal sampling strategy, the control dataset contains 1000 individuals
            Nj <- floor( Sumj*N/100)
            Xsample <- sample(Pool1,N1,replace=F)#Sample from X1 pool,
            S1 <- t(as.matrix(table(Xsample)))
            Xsample2 <-  sample(Poolj,Nj,replace=F)#Sample from X2 pool
            S2 <- t(as.matrix(table(Xsample2))) 
            S1S2 <-  rbind.fill.matrix(S1,S2)
            S1S2 [is.na(S1S2)] <- 0
            
            for (m in c(1,10,100,1000))
            {
                if (m>Nj){next}
                if(m==1){Mean1[k] <- ESS(S1S2,m=m)}
                if(m==10){Mean2[k] <- ESS(S1S2,m=m)}
                if(m==100){Mean3[k] <- ESS(S1S2,m=m)}
                if(m==1000){Mean4[k] <- ESS(S1S2,m=m)}
            }
            MeanBC[k]<- vegdist (S1S2,method="bray",binary=F)# Bray-Curtis index
            MeanEuc[k]<- vegdist (decostand(S1S2, method = "total"),method="euclidean")#Proportion based Euclidian distance
            MeanChao[k] <- dis.chao(S1S2,index="sorensen")#Chao-Sorense index
        }
        Distmean1 <- mean(Mean1)
        Bounds1 <- quantile(Mean1,probs = c(0.025, 0.975))
        Lower1 <- Bounds1[1]; Upper1 <- Bounds1[2]
        CV1 <- sd(Mean1)/Distmean1
        Df1 <-  c(N=N,X=j,m=1,Dst=Distmean1,Lo=Lower1,Up=Upper1,CV=CV1)
        
        Distmean2 <- mean(Mean2)
        Bounds2 <- quantile(Mean2,probs = c(0.025, 0.975))
        Lower2 <- Bounds2[1]; Upper2 <- Bounds2[2]
        CV2 <- sd(Mean2)/Distmean2
        Df2 <-  c(N=N,X=j,m=10,Dst=Distmean2,Lo2=Lower2,Up=Upper2,CV=CV2)
        
        Distmean3 <- mean(Mean3)
        Bounds3 <- quantile(Mean3,probs = c(0.025, 0.975))
        Lower3 <- Bounds3[1]; Upper3 <- Bounds3[2]
        CV3 <- sd(Mean3)/Distmean3
        Df3 <-  c(N=N,X=j,m=100,Dst=Distmean3,Lo=Lower3,Up=Upper3,CV=CV3)
        
        Distmean4 <- mean(Mean4)
        Bounds4 <- quantile(Mean4,probs = c(0.025, 0.975))
        Lower4 <- Bounds4[1]; Upper4 <- Bounds4[2]
        CV4 <- sd(Mean4)/Distmean4
        Df4 <-  c(N=N,X=j,m=1000,Dst=Distmean4,Lo=Lower4,Up=Upper4,CV=CV4)
        
        DistBC <- mean(MeanBC)
        BoudnsBC <- quantile(MeanBC,probs = c(0.025, 0.975))
        LowerBC <- BoudnsBC[1]; UpperBC <- BoudnsBC[2]
        CV5 <- sd(MeanBC)/DistBC
        Df5 <-  c(N=N,X=j,m="BC",Dst=DistBC,Lo=LowerBC,Up=UpperBC,CV=CV5)
        
        DistEuc <- mean(MeanEuc)
        BoudnsEuc <- quantile(MeanEuc,probs = c(0.025, 0.975))
        LowerEuc <- BoudnsEuc[1]; UpperEuc <- BoudnsEuc[2]
        CV6 <- sd(MeanEuc)/DistEuc
        Df6 <-  c(N=N,X=j,m="Euc",Dst=DistEuc,Lo=LowerEuc,Up=UpperEuc,CV=CV6)
        
        DistChao <- mean(MeanChao)
        BoudnsChao <- quantile(MeanChao,probs = c(0.025, 0.975))
        LowerChao <- BoudnsChao[1]; UpperChao <- BoudnsChao[2]
        CV7 <- sd(MeanChao)/DistChao
        Df7 <-  c(N=N,X=j,m="SoAbu",Dst=DistChao,Lo=LowerChao,Up=UpperChao,CV=CV7)
        
        Df <- rbind(Df1,Df2,Df3,Df4,Df5,Df6,Df7)
        Out2 <- rbind(Out2,Df)
    }
}
rownames(Out2) <- 1:nrow(Out2)

