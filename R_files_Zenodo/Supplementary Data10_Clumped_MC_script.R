library(xlsx)
library(TTR)

# Load and order data
setwd("<path>")
dat<-read.csv("<name_data_file.csv>",stringsAsFactors = FALSE)

# Subsample data and create running matrices
N = 1000 # Set number of iterations for Monte Carlo simulations. 100 is enough for most data, creates >10.000 successful simulations
p = 0.05 # Set threshold for significant p-value
d18O<-dat[,c(6,2)] # Isolate d18O data
D47<-dat[,c(6,4)] # Isolate D47 data
Popt<-vector()
SD_d18O<-0.05 # Set standard deviation of error on d18O measurements
SD_D47<-0.044 # Set standard deviation of error on D47 measurements

# MONTE CARLO SIMULATION

for(i in 1:N){
    # Progress
    cat(paste("Modeling Iteration: ",i),"\\r")
    flush.console()

    d18O_mod<-rnorm(length(d18O[,2]),d18O[,2],SD_d18O) # Randomly resample d18O data using measurement uncertainty
    D47_mod<-rnorm(length(D47[,2]),D47[,2],SD_D47) # Randomly resample D47 data using measurement uncertainty
    
    # Combine d18O and D47
    X<-cbind(d18O_mod,D47_mod)
    # Sort by d18O
    X1<-X[order(X[,1]),]
    # Inverse sort by d18O
    X2<-X[order(X[,1],decreasing=TRUE),]

    # Expanding window seasonality and T-test
    # Keep record of summer and winter values for successful sample windows
    win<-seq(1,length(d18O_mod),1) # Create vector of sample size windows
    dsum<-runMean(X1[,1],1,cumulative=TRUE) # Calculate average summer d18O for progressively larger sample size windows
    dwin<-runMean(X2[,1],1,cumulative=TRUE) # Calculate average winter d18O for progressively larger sample size windows
    dsumsd<-runSD(X1[,1],1,cumulative=TRUE) # Calculate standard deviation within summer d18O values for progressively larger sample size windows
    dwinsd<-runSD(X2[,1],1,cumulative=TRUE) # Calculate standard deviation within winter d18O values for progressively larger sample size windows
    Dsum<-runMean(X1[,2],1,cumulative=TRUE) # Calculate average summer D47 for progressively larger sample size windows
    Dwin<-runMean(X2[,2],1,cumulative=TRUE) # Calculate average winter D47 for progressively larger sample size windows
    Dsumsd<-runSD(X1[,2],1,cumulative=TRUE) # Calculate standard deviation within summer D47 values for progressively larger sample size windows
    Dwinsd<-runSD(X2[,2],1,cumulative=TRUE) # Calculate standard deviation within winter D47 values for progressively larger sample size windows
    SDpool<-sqrt((Dsumsd^2+Dwinsd^2)/2) # Calculate pooled standard deviation for each window
    T<-(Dsum-Dwin)/(SDpool*sqrt(2/win)) # Calculate two-sample T-value for each window (equal sample size, equal variance)
    Pval<-pt(T,win-1) # Calculate p-value for each window
    res<-cbind(rep(i,length(which(Pval<p))),
        win[which(Pval<p)],
        Pval[which(Pval<p)],
        dsum[which(Pval<p)],
        dsumsd[which(Pval<p)]/sqrt(win[which(Pval<p)]),
        dwin[which(Pval<p)],
        dwinsd[which(Pval<p)]/sqrt(win[which(Pval<p)]),
        Dsum[which(Pval<p)],
        Dsumsd[which(Pval<p)]/sqrt(win[which(Pval<p)]),
        Dwin[which(Pval<p)],
        Dwinsd[which(Pval<p)]/sqrt(win[which(Pval<p)])
        ) # Combine results for export
    Popt<-rbind(Popt,res)
}

# POST PROCESSING

# Add temperature calculations of optimal runs (following Kele et al., 2015 recalculated by Bernasconi et al., 2018)
Popt<-cbind(Popt,sqrt((0.0449*10^6)/(Popt[,8]-0.167))-273.15)
Popt<-cbind(Popt,(sqrt((0.0449*10^6)/(Popt[,8]-Popt[,9]-0.167))-273.15)-(sqrt((0.0449*10^6)/(Popt[,8]-0.167))-273.15))
Popt<-cbind(Popt,sqrt((0.0449*10^6)/(Popt[,10]-0.167))-273.15)
Popt<-cbind(Popt,(sqrt((0.0449*10^6)/(Popt[,10]-Popt[,11]-0.167))-273.15)-(sqrt((0.0449*10^6)/(Popt[,10]-0.167))-273.15))

# Add seawater d18O calculations of optimal runs (following Kim and O'Neil, 1997)
Popt<-cbind(Popt,((Popt[,4]-(exp(((18.03*10^3)/(Popt[,12]+273)-32.42)/1000)-1)*1000)*1.03092+30.92))
Popt<-cbind(Popt,(((Popt[,4]+Popt[,5])-(exp(((18.03*10^3)/((Popt[,12]+Popt[,13])+273)-32.42)/1000)-1)*1000)*1.03092+30.92)-((Popt[,4]-(exp(((18.03*10^3)/(Popt[,12]+273)-32.42)/1000)-1)*1.03092+30.92))
Popt<-cbind(Popt,((Popt[,6]-(exp(((18.03*10^3)/(Popt[,14]+273)-32.42)/1000)-1)*1000)*1.03092+30.92))
Popt<-cbind(Popt,(((Popt[,6]+Popt[,7])-(exp(((18.03*10^3)/((Popt[,14]+Popt[,15])+273)-32.42)/1000)-1)*1000)*1.03092+30.92)-((Popt[,6]-(exp(((18.03*10^3)/(Popt[,14]+273)-32.42)/1000)-1)*1.03092+30.92))

# Add slopes and intercepts for D47-d18O conversion
Popt<-cbind(Popt,(Popt[,8]-Popt[,10])/(Popt[,4]-Popt[,6]))
Popt<-cbind(Popt,((Popt[,8]+Popt[,10])/2)-Popt[,20]*((Popt[,4]+Popt[,6])/2))

colnames(Popt)<-c("run#","optimal sample size","p-value","dOsum","dOsumse","dOwin","dOwinse","Dsum","Dsumse","Dwin","Dwinse","Tsum","Tsumse","Twin","Twinse","dswsum","dswsumse","dswwin","dswwinse","D47-d18O_slope","D47-d18O_int")

# Export results of optimized sample sizes
write.csv(Popt,"Popt.csv")


# Calculate clustering statistics
wtav<-c(sum(Popt[,4]*Popt[,2])/sum(Popt[,2]),
    sum(Popt[,6]*Popt[,2])/sum(Popt[,2]),
    sum(Popt[,8]*Popt[,2])/sum(Popt[,2]),
    sum(Popt[,10]*Popt[,2])/sum(Popt[,2]),
    sum(Popt[,12]*Popt[,2])/sum(Popt[,2]),
    sum(Popt[,14]*Popt[,2])/sum(Popt[,2]),
    sum(Popt[,16]*Popt[,2])/sum(Popt[,2]),
    sum(Popt[,18]*Popt[,2])/sum(Popt[,2])
    ) # Calculate weighed averages of clusters for summers and winters
wtsd<-c(sqrt((sum(Popt[,2]*(Popt[,5]*sqrt(Popt[,2]))^2)+sum(Popt[,2]*(Popt[,4]-wtav[1])^2))/sum(Popt[,2])),
    sqrt((sum(Popt[,2]*(Popt[,7]*sqrt(Popt[,2]))^2)+sum(Popt[,2]*(Popt[,6]-wtav[2])^2))/sum(Popt[,2])),
    sqrt((sum(Popt[,2]*(Popt[,9]*sqrt(Popt[,2]))^2)+sum(Popt[,2]*(Popt[,8]-wtav[3])^2))/sum(Popt[,2])),
    sqrt((sum(Popt[,2]*(Popt[,11]*sqrt(Popt[,2]))^2)+sum(Popt[,2]*(Popt[,10]-wtav[4])^2))/sum(Popt[,2])),
    sqrt((sum(Popt[,2]*(Popt[,13]*sqrt(Popt[,2]))^2)+sum(Popt[,2]*(Popt[,12]-wtav[5])^2))/sum(Popt[,2])),
    sqrt((sum(Popt[,2]*(Popt[,15]*sqrt(Popt[,2]))^2)+sum(Popt[,2]*(Popt[,14]-wtav[6])^2))/sum(Popt[,2])),
    sqrt((sum(Popt[,2]*(Popt[,17]*sqrt(Popt[,2]))^2)+sum(Popt[,2]*(Popt[,16]-wtav[7])^2))/sum(Popt[,2])),
    sqrt((sum(Popt[,2]*(Popt[,19]*sqrt(Popt[,2]))^2)+sum(Popt[,2]*(Popt[,18]-wtav[8])^2))/sum(Popt[,2]))
    ) # Calculate weighed standard deviations of clusters for summers and winters
clustersd<-c(sd(Popt[,4]),
    sd(Popt[,6]),
    sd(Popt[,8]),
    sd(Popt[,10]),
    sd(Popt[,12]),
    sd(Popt[,14]),
    sd(Popt[,16]),
    sd(Popt[,18])
    ) # Calculate weighed standard deviations of clusters for summers and winters
CL<-1.96*wtsd/sqrt(mean(Popt[,2]))
clusterCL<-1.96*clustersd/sqrt(mean(Popt[,2]))
cluster_stats<-cbind(wtav,wtsd,CL,clustersd,clusterCL)
colnames(cluster_stats)<-c("weighed average","weighed standard deviation","95% CL","Cluster standard deviation","Cluster standard error")
rownames(cluster_stats)<-c("dsum","dwin","Dsum","Dwin","Tsum","Twin","dwsum","dwwin")

write.csv(cluster_stats,"cluster_stats.csv")


# Open d18O data for D47 modelling
VUB18<-read.csv("VUB18.csv",stringsAsFactors = FALSE)
comb<-cbind(c(dat[,7],VUB18[,1]),c(dat[,2],VUB18[,2])) # Combine d18O results

# Use Popt matrix to model D47, T and d18Osw from d18O data
combD<-as.data.frame(cbind(comb,outer(comb[,2],Popt[,20])+rep(Popt[,21],each=nrow(comb))))
combT<-as.data.frame(cbind(comb,sqrt((0.0449*10^6)/(combD[,3:length(combD[1,])]-0.167))-273.15))
combd<-as.data.frame(cbind(comb,(matrix(combT[,2],nrow=length(comb[,2]),ncol=length(combT[1,])-2)-(exp(((18.03*10^3)/(combT[,3:length(combT[1,])]+273)-32.42)/1000)-1)*1000)*1.03092+30.92))

# Stack results for export
tot<-cbind(rep(comb[,1],ncol(combD)-2),
    rep(comb[,2],ncol(combD)-2),
    stack(combD[,3:ncol(combD)])[,1],
    stack(combT[,3:ncol(combT)])[,1],
    stack(combd[,3:ncol(combd)])[,1]
    )

# Add corresponding p-values
tot<-cbind(tot,rep(Popt[,3],each=nrow(combD)))
# Add modulo of age
tot<-cbind(tot[,1]%%1,tot)
colnames(tot)<-c("mod(year)","Age_yr","d18Occ","D47_mod","T_mod","d18Osw_mod","p_value")
# Export
write.csv(tot,"modcombined.csv")

# Use Popt matrix to model monthly D47, T and d18Osw from d18O data
# Calculate stats per month
means<-vector()
medians<-vector()
sds<-vector()
Nsam<-vector()
for(i in 1:12){
    means<-rbind(means,apply(tot[which(tot[,1]>=((i-1)/12) & tot[,1]<(i/12)),],2,mean))
    medians<-rbind(medians,apply(tot[which(tot[,1]>=((i-1)/12) & tot[,1]<(i/12)),],2,median))
    sds<-rbind(sds,apply(tot[which(tot[,1]>=((i-1)/12) & tot[,1]<(i/12)),],2,sd))
    Nsam<-append(Nsam,length(comb[which((comb[,1]%%1)>=((i-1)/12) & (comb[,1]%%1)<(i/12)),2]))
}
# Combine and format
monthly<-cbind(means,medians,sds,sds/sqrt(Nsam),1.96*sds/sqrt(Nsam),Nsam)
colnames(monthly)<-c(paste("Mean_",colnames(tot)),
    paste("Median_",colnames(tot)),
    paste("StDev_",colnames(tot)),
    paste("StErr_",colnames(tot)),
    paste("95% CL_",colnames(tot)),
    "N"
    )
# Export monthly values
write.csv(monthly,"monthly.csv")

# Calculate three month summer and winter averages

# Find highest and lowest D47 three month period
Dlist<-cbind(c(12,1:12,1),c(monthly[12,4],monthly[,4],monthly[1,4])) # List all monthly temperatures + months and add last before and first after to allow moving average
if(sum(is.na(Dlist[,2]))>0){Dlist<-Dlist[-as.numeric(which(is.na(Dlist[,2]))),]} # Remove NA's of underrepresented months
D3M<-cbind(Dlist[-c(12,13),1],Dlist[-c(1,13),1],Dlist[-c(1,2),1],(Dlist[-c(1,2),2]+Dlist[-c(1,13),2]+Dlist[-c(12,13),2])/3) # Calculate three month moving temperature averages and list by first, middle and last month

C3M<-D3M[which.max(D3M[,4]),1:3] # Extract coldest (highest D47) three months
W3M<-D3M[which.min(D3M[,4]),1:3] # Extract warmest (lowest D47) three months

# Calculate statistics of warmest and coldest 3 month period
wtav<-c(colSums(monthly[W3M,3:6]*monthly[W3M,36])/sum(monthly[W3M,36]),
    colSums(monthly[C3M,3:6]*monthly[C3M,36])/sum(monthly[C3M,36])
    ) # Weighted averages
wtsd<-c(sqrt((colSums(monthly[W3M,36]*monthly[W3M,17:20]^2)+colSums(monthly[W3M,36]*(monthly[W3M,3:6]-t(matrix(rep(wtav[1:4],3),nrow=4)))^2))/sum(monthly[W3M,36])),
    sqrt((colSums(monthly[C3M,36]*monthly[C3M,17:20]^2)+colSums(monthly[C3M,36]*(monthly[C3M,3:6]-t(matrix(rep(wtav[5:8],3),nrow=4)))^2))/sum(monthly[C3M,36]))
    ) #Weighted standard deviation
wtse<-wtsd/c(sqrt(sum(monthly[W3M,36])),sqrt(sum(monthly[C3M,36])))
CL<-1.96*wtse
W3Mstats<-c(paste("#",W3M[1],"-#",W3M[3],sep=""),wtav[1:4],wtsd[1:4],wtse[1:4],CL[1:4])
C3Mstats<-c(paste("#",C3M[1],"-#",C3M[3],sep=""),wtav[5:8],wtsd[5:8],wtse[5:8],CL[5:8])
three_month_stats<-rbind(W3Mstats,C3Mstats)
colnames(three_month_stats)<-c("Months",
    paste("weighed average_",colnames(tot)[3:6]),
    paste("weiged stdev_",colnames(tot)[3:6]),
    paste("weiged sterr_",colnames(tot)[3:6]),
    paste("95% CL_",colnames(tot)[3:6])
    )
rownames(three_month_stats)<-c("Summer","Winter")

write.csv(three_month_stats,"Three_month_stats.csv")

# Compile stat summary
statsummary1<-cbind(1:12,monthly[,36],monthly[,12],monthly[,33],monthly[,13],monthly[,34])
statsummary2<-cbind(c("Summer","Winter"),c(mean(Popt[,2]),mean(Popt[,2])),cluster_stats[5:6,1],cluster_stats[5:6,4],cluster_stats[7:8,1],cluster_stats[7:8,4])
statsummary3<-cbind(c("Summer","Winter"),c(sum(monthly[W3M,36]),sum(monthly[C3M,36])),three_month_stats[,4],three_month_stats[,16],three_month_stats[,5],three_month_stats[,17],three_month_stats[,1])
statsummary<-rbind(c("Month","N_sample","Median T","95% CL","Median d18Osw","95% CL",""),
    cbind(statsummary1,rep("",nrow(statsummary1))),
    rep("",7),
    c("Clusters","Mean_N_cluster","weighed average T","SD cluster","weighed average d18Osw","SD cluster",""),
    cbind(statsummary2,rep("",nrow(statsummary2))),
    rep("",7),
    c("Three month average","N_sam","weighed average T","95% CL","weighed average d18Osw","95% CL","Months"),
    statsummary3
)

write.csv(statsummary,"stat_summary.csv")