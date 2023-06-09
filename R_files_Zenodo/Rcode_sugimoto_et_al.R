##################################
##R code for Sugimoto et al. Hierarchical trait-based model reveals positive and negative effects of land abandonment on butterfly communities across climatic regions in Japan
##May 14,2021
##################################

#R.4.0.2 / rstan Version 2.21.1

library(tidyverse)
library(lubridate)
library(inline)
library(Rcpp)
library(RcppParallel)
library(rstan)
packageVersion("rstan")
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

#reading dataset

dataset<-read.csv("dataset.csv")	#dataset_maskRL.csv: a dataset with red list species masked

Y<-dataset$Y    #response variable
speciesID<-dataset$speciesID    #group ID of species
nsamp<-length(Y)    #number of samples
ns<-speciesID%>%unique%>%length #number of species
regionID<-dataset$regionID  #group ID of regions
nregion<-regionID%>%unique%>%length     #number of regions
settlementID<-dataset$settlementID      #group ID of settlements
nsettlement<-settlementID%>%unique%>%length #number of settlements

Xtemp<-scale(dataset$Xtemp)     #mean annual temperature
Xtemp.mu<-attr(Xtemp,"scaled:center")
Xtemp.sigma<-attr(Xtemp,"scaled:scale")
Xtemp<-c(Xtemp)
Xaban_bin<-dataset$Xaban.binary #abandonment (binary)
Xaban_year<-scale(dataset$Xaban.year) #years since abandonment
Xaban_year.mu<-attr(Xaban_year,"scaled:center")
Xaban_year.sd<-attr(Xaban_year,"scaled:scale")
Xaban_year<-c(Xaban_year)

Xmonth<-scale(dataset$Xmonth)
Xmonth.mu<-attr(Xmonth,"scaled:center")
Xmonth.sd<-attr(Xmonth,"scaled:scale")
Xmonth<-c(Xmonth)

Xdryfield<-dataset$Xdryfield
Xbuiltup<-dataset$Xbuiltup




###baseline model（binary variable of abandonment）
stanmodelcode <- "
data {
	int<lower=0> ns;
    int<lower=0> nsamp;
    int<lower=0> nsettlement;
	int<lower=0> Y[nsamp];
	int<lower=0> regionID[nsamp];
	int<lower=0> nregion;
    int<lower=0> speciesID[nsamp];
    int<lower=0> settlementID[nsamp];
	vector[nsamp] Xtemp;
	vector[nsamp] Xaban_bin;
	vector[nsamp] Xmonth;
	vector[nsamp] Xdryfield;
	vector[nsamp] Xbuiltup;
} 

parameters {
	real<lower=0,upper=pi()/2> sigmaa_unif;
	real<lower=0,upper=pi()/2> sigmas_unif;
	real<lower=0,upper=1> a0_unif[5];
	matrix[5,ns] beta_stdn;
	vector<lower=0,upper=1>[ns] beta0_unif;
	vector[nregion] ranef_a;
    vector[nsettlement] ranef_s;
	cholesky_factor_corr[2] L_Omega;
    vector<lower=0,upper=pi()/2>[5] sigmab_unif;
} 

transformed parameters {
	matrix[2,2] rho;
	matrix[ns,5] mu_beta;
    vector<lower=0>[5] sigmab;
    vector[ns] beta3;
    vector[ns] beta4;
    vector[ns] beta5;
    vector[nsamp] logitp;
    real<lower=0> sigmaa;
    real<lower=0> sigmas;
	matrix[ns,2] beta;
	vector[ns] beta0;
	real a0[5];

    for(k in 1:5) sigmab[k]=5*tan(sigmab_unif[k]);
    sigmaa = 5*tan(sigmaa_unif);
    sigmas = 5*tan(sigmas_unif);    
    rho = L_Omega * L_Omega';
    for(k in 1:5) a0[k] = 10 * inv_Phi(a0_unif[k]);
    for(i in 1:ns) beta0[i] = 10*inv_Phi(beta0_unif[i]);
	for(i in 1:ns){
		for(j in 1:5){
			mu_beta[i,j] = a0[j];
		}
	}
    beta = mu_beta[,1:2] + (diag_pre_multiply(sigmab[1:2],L_Omega)*beta_stdn[1:2,])';
    beta3 = mu_beta[,3] + sigmab[3]*beta_stdn[3,]';
    beta4 = mu_beta[,4] + sigmab[4]*beta_stdn[4,]';
    beta5 = mu_beta[,5] + sigmab[5]*beta_stdn[5,]';
    for(i in 1:nsamp) logitp[i] = beta0[speciesID[i]]+beta[speciesID[i],1]*Xtemp[i]+beta[speciesID[i],2]*Xaban_bin[i]+beta3[speciesID[i]]*Xmonth[i]+beta4[speciesID[i]]*Xdryfield[i]+beta5[speciesID[i]]*Xbuiltup[i]+ranef_s[settlementID[i]]*sigmas+ranef_a[regionID[i]]*sigmaa;
 
}

model {
	to_vector(beta_stdn)~std_normal();
	Y~bernoulli_logit(logitp);
	L_Omega~lkj_corr_cholesky(1.0); 
	ranef_a~std_normal();
	ranef_s~std_normal();

} 

generated quantities{
    vector[nsamp] log_lik;
    for(i in 1:nsamp) log_lik[i] = bernoulli_logit_lpmf(Y[i] | logitp[i]);
}
"

data<-list(Y=Y,speciesID=speciesID,nsamp=nsamp,ns=ns,regionID=regionID,nregion=nregion,settlementID=settlementID,nsettlement=nsettlement,Xtemp=Xtemp,Xaban_bin=Xaban_bin,Xmonth=Xmonth,Xdryfield=Xdryfield,Xbuiltup=Xbuiltup)
pars<-c("a0","beta0","beta","beta3","beta4","beta5","sigmaa","sigmas","sigmab","ranef_a","ranef_s","rho","log_lik")
res.baseline<-stan(model_code=stanmodelcode,model_name="abandonmentbutterfly",data=data,pars=pars,iter=2000,thin=1,chains = 3,control=list(adapt_delta=0.99,max_treedepth = 15))


###"year" model（years since abandonment）

stanmodelcode <- "
data {
	int<lower=0> ns;
    int<lower=0> nsamp;
    int<lower=0> nsettlement;
	int<lower=0> Y[nsamp];
	int<lower=0> regionID[nsamp];
	int<lower=0> nregion;
    int<lower=0> speciesID[nsamp];
    int<lower=0> settlementID[nsamp];
	vector[nsamp] Xtemp;
	vector[nsamp] Xaban_year;
	vector[nsamp] Xmonth;
	vector[nsamp] Xdryfield;
	vector[nsamp] Xbuiltup;
} 

parameters {
	real<lower=0,upper=pi()/2> sigmaa_unif;
	real<lower=0,upper=pi()/2> sigmas_unif;
	real<lower=0,upper=1> a0_unif[5];
	matrix[5,ns] beta_stdn;
	vector<lower=0,upper=1>[ns] beta0_unif;
	vector[nregion] ranef_a;
    vector[nsettlement] ranef_s;
	cholesky_factor_corr[2] L_Omega;
    vector<lower=0,upper=pi()/2>[5] sigmab_unif;
} 

transformed parameters {
	matrix[2,2] rho;
	matrix[ns,5] mu_beta;
    vector<lower=0>[5] sigmab;
    vector[ns] beta3;
    vector[ns] beta4;
    vector[ns] beta5;
    vector[nsamp] logitp;
    real<lower=0> sigmaa;
    real<lower=0> sigmas;
	matrix[ns,2] beta;
	vector[ns] beta0;
	real a0[5];

    for(k in 1:5) sigmab[k]=5*tan(sigmab_unif[k]);
    sigmaa = 5*tan(sigmaa_unif);
    sigmas = 5*tan(sigmas_unif);    
    rho = L_Omega * L_Omega';
    for(k in 1:5) a0[k] = 10 * inv_Phi(a0_unif[k]);
    for(i in 1:ns) beta0[i] = 10*inv_Phi(beta0_unif[i]);
	for(i in 1:ns){
		for(j in 1:5){
			mu_beta[i,j] = a0[j];
		}
	}
    beta = mu_beta[,1:2] + (diag_pre_multiply(sigmab[1:2],L_Omega)*beta_stdn[1:2,])';
    beta3 = mu_beta[,3] + sigmab[3]*beta_stdn[3,]';
    beta4 = mu_beta[,4] + sigmab[4]*beta_stdn[4,]';
    beta5 = mu_beta[,5] + sigmab[5]*beta_stdn[5,]';
    for(i in 1:nsamp) logitp[i] = beta0[speciesID[i]]+beta[speciesID[i],1]*Xtemp[i]+beta[speciesID[i],2]*Xaban_year[i]+beta3[speciesID[i]]*Xmonth[i]+beta4[speciesID[i]]*Xdryfield[i]+beta5[speciesID[i]]*Xbuiltup[i]+ranef_s[settlementID[i]]*sigmas+ranef_a[regionID[i]]*sigmaa;
 
}

model {
	to_vector(beta_stdn)~std_normal();
	Y~bernoulli_logit(logitp);
	L_Omega~lkj_corr_cholesky(1.0); 
	ranef_a~std_normal();
	ranef_s~std_normal();

} 

generated quantities{
    vector[nsamp] log_lik;
    for(i in 1:nsamp) log_lik[i] = bernoulli_logit_lpmf(Y[i] | logitp[i]);
}
"



data<-list(Y=Y,speciesID=speciesID,nsamp=nsamp,ns=ns,regionID=regionID,nregion=nregion,settlementID=settlementID,nsettlement=nsettlement,Xtemp=Xtemp,Xaban_year=Xaban_year,Xmonth=Xmonth,Xdryfield=Xdryfield,Xbuiltup=Xbuiltup)
pars<-c("a0","beta0","beta","sigmaa","sigmas","sigmab","ranef_a","ranef_s","rho","log_lik")
res.year<-stan(model_code=stanmodelcode,model_name="abandonmentbutterfly",data=data,pars=pars,iter=2000,thin=1,chains = 3,control=list(adapt_delta=0.99,max_treedepth = 15))

###WBIC calculation （baseline model）

stanmodelcode <- "
data {
	int<lower=0> ns;
    int<lower=0> nsamp;
    int<lower=0> nsettlement;
	int<lower=0> Y[nsamp];
	int<lower=0> regionID[nsamp];
	int<lower=0> nregion;
    int<lower=0> speciesID[nsamp];
    int<lower=0> settlementID[nsamp];
	vector[nsamp] Xtemp;
	vector[nsamp] Xaban_bin;
	vector[nsamp] Xmonth;
	vector[nsamp] Xdryfield;
	vector[nsamp] Xbuiltup;
} 

parameters {
	real<lower=0,upper=pi()/2> sigmaa_unif;
	real<lower=0,upper=pi()/2> sigmas_unif;
	real<lower=0,upper=1> a0_unif[5];
	matrix[5,ns] beta_stdn;
	vector<lower=0,upper=1>[ns] beta0_unif;
	vector[nregion] ranef_a;
    vector[nsettlement] ranef_s;
	cholesky_factor_corr[2] L_Omega;
    vector<lower=0,upper=pi()/2>[5] sigmab_unif;
} 

transformed parameters {
	matrix[2,2] rho;
	matrix[ns,5] mu_beta;
    vector<lower=0>[5] sigmab;
    vector[ns] beta3;
    vector[ns] beta4;
    vector[ns] beta5;
    vector[nsamp] logitp;
    real<lower=0> sigmaa;
    real<lower=0> sigmas;
	matrix[ns,2] beta;
	vector[ns] beta0;
	real a0[5];

    for(k in 1:5) sigmab[k]=5*tan(sigmab_unif[k]);
    sigmaa = 5*tan(sigmaa_unif);
    sigmas = 5*tan(sigmas_unif);    
    rho = L_Omega * L_Omega';
    for(k in 1:5) a0[k] = 10 * inv_Phi(a0_unif[k]);
    for(i in 1:ns) beta0[i] = 10*inv_Phi(beta0_unif[i]);
	for(i in 1:ns){
		for(j in 1:5){
			mu_beta[i,j] = a0[j];
		}
	}
    beta = mu_beta[,1:2] + (diag_pre_multiply(sigmab[1:2],L_Omega)*beta_stdn[1:2,])';
    beta3 = mu_beta[,3] + sigmab[3]*beta_stdn[3,]';
    beta4 = mu_beta[,4] + sigmab[4]*beta_stdn[4,]';
    beta5 = mu_beta[,5] + sigmab[5]*beta_stdn[5,]';
    for(i in 1:nsamp) logitp[i] = beta0[speciesID[i]]+beta[speciesID[i],1]*Xtemp[i]+beta[speciesID[i],2]*Xaban_bin[i]+beta3[speciesID[i]]*Xmonth[i]+beta4[speciesID[i]]*Xdryfield[i]+beta5[speciesID[i]]*Xbuiltup[i]+ranef_s[settlementID[i]]*sigmas+ranef_a[regionID[i]]*sigmaa;
 
}

model {
	to_vector(beta_stdn)~std_normal();
	for(i in 1:nsamp) target += 1/log(nsamp)*bernoulli_logit_lpmf(Y[i] | logitp[i]);
	L_Omega~lkj_corr_cholesky(1.0); 
	ranef_a~std_normal();
	ranef_s~std_normal();

} 

generated quantities{
    vector[nsamp] log_lik;
    for(i in 1:nsamp) log_lik[i] = bernoulli_logit_lpmf(Y[i] | logitp[i]);
}
"

data<-list(Y=Y,speciesID=speciesID,nsamp=nsamp,ns=ns,regionID=regionID,nregion=nregion,settlementID=settlementID,nsettlement=nsettlement,Xtemp=Xtemp,Xaban_bin=Xaban_bin,Xmonth=Xmonth,Xdryfield=Xdryfield,Xbuiltup=Xbuiltup)
pars<-c("a0","beta0","beta","sigmaa","sigmas","sigmab","ranef_a","ranef_s","rho","log_lik")
res.baselineBIC<-stan(model_code=stanmodelcode,model_name="abandonmentbutterfly",data=data,pars=pars,iter=2000,thin=1,chains = 3,control=list(adapt_delta=0.99,max_treedepth = 15))


###WBIC calculation（year model）
stanmodelcode <- "
data {
	int<lower=0> ns;
    int<lower=0> nsamp;
    int<lower=0> nsettlement;
	int<lower=0> Y[nsamp];
	int<lower=0> regionID[nsamp];
	int<lower=0> nregion;
    int<lower=0> speciesID[nsamp];
    int<lower=0> settlementID[nsamp];
	vector[nsamp] Xtemp;
	vector[nsamp] Xaban_year;
	vector[nsamp] Xmonth;
	vector[nsamp] Xdryfield;
	vector[nsamp] Xbuiltup;
} 

parameters {
	real<lower=0,upper=pi()/2> sigmaa_unif;
	real<lower=0,upper=pi()/2> sigmas_unif;
	real<lower=0,upper=1> a0_unif[5];
	matrix[5,ns] beta_stdn;
	vector<lower=0,upper=1>[ns] beta0_unif;
	vector[nregion] ranef_a;
    vector[nsettlement] ranef_s;
	cholesky_factor_corr[2] L_Omega;
    vector<lower=0,upper=pi()/2>[5] sigmab_unif;
} 

transformed parameters {
	matrix[2,2] rho;
	matrix[ns,5] mu_beta;
    vector<lower=0>[5] sigmab;
    vector[ns] beta3;
    vector[ns] beta4;
    vector[ns] beta5;
    vector[nsamp] logitp;
    real<lower=0> sigmaa;
    real<lower=0> sigmas;
	matrix[ns,2] beta;
	vector[ns] beta0;
	real a0[5];

    for(k in 1:5) sigmab[k]=5*tan(sigmab_unif[k]);
    sigmaa = 5*tan(sigmaa_unif);
    sigmas = 5*tan(sigmas_unif);    
    rho = L_Omega * L_Omega';
    for(k in 1:5) a0[k] = 10 * inv_Phi(a0_unif[k]);
    for(i in 1:ns) beta0[i] = 10*inv_Phi(beta0_unif[i]);
	for(i in 1:ns){
		for(j in 1:5){
			mu_beta[i,j] = a0[j];
		}
	}
    beta = mu_beta[,1:2] + (diag_pre_multiply(sigmab[1:2],L_Omega)*beta_stdn[1:2,])';
    beta3 = mu_beta[,3] + sigmab[3]*beta_stdn[3,]';
    beta4 = mu_beta[,4] + sigmab[4]*beta_stdn[4,]';
    beta5 = mu_beta[,5] + sigmab[5]*beta_stdn[5,]';
    for(i in 1:nsamp) logitp[i] = beta0[speciesID[i]]+beta[speciesID[i],1]*Xtemp[i]+beta[speciesID[i],2]*Xaban_year[i]+beta3[speciesID[i]]*Xmonth[i]+beta4[speciesID[i]]*Xdryfield[i]+beta5[speciesID[i]]*Xbuiltup[i]+ranef_s[settlementID[i]]*sigmas+ranef_a[regionID[i]]*sigmaa;
 
}

model {
	to_vector(beta_stdn)~std_normal();
	for(i in 1:nsamp) target += 1/log(nsamp)*bernoulli_logit_lpmf(Y[i] | logitp[i]);
	L_Omega~lkj_corr_cholesky(1.0); 
	ranef_a~std_normal();
	ranef_s~std_normal();

} 

generated quantities{
    vector[nsamp] log_lik;
    for(i in 1:nsamp) log_lik[i] = bernoulli_logit_lpmf(Y[i] | logitp[i]);
}
"


data<-list(Y=Y,speciesID=speciesID,nsamp=nsamp,ns=ns,regionID=regionID,nregion=nregion,settlementID=settlementID,nsettlement=nsettlement,Xtemp=Xtemp,Xaban_year=Xaban_year,Xmonth=Xmonth,Xdryfield=Xdryfield,Xbuiltup=Xbuiltup)
pars<-c("a0","beta0","beta","sigmaa","sigmas","sigmab","ranef_a","ranef_s","rho","log_lik")
res.yearBIC<-stan(model_code=stanmodelcode,model_name="abandonmentbutterfly",data=data,pars=pars,iter=2000,thin=1,chains = 3,control=list(adapt_delta=0.99,max_treedepth = 15))

resWBIC<-numeric(2)
resWBIC[1]<-(-1)*mean(rowSums(rstan::extract(res.baselineBIC)$log_lik))
resWBIC[2]<-(-1)*mean(rowSums(rstan::extract(res.yearBIC)$log_lik))
exp(diff(resWBIC))


##habitat trait model

stanmodelcode <- "
data {
	int<lower=0> ns;
    int<lower=0> nsamp;
    int<lower=0> nsettlement;
	int<lower=0> Y[nsamp];
	int<lower=0> regionID[nsamp];
	int<lower=0> nregion;
    int<lower=0> speciesID[nsamp];
    int<lower=0> settlementID[nsamp];
	vector[nsamp] Xtemp;
	vector[nsamp] Xaban_bin;
	vector[nsamp] Xmonth;
	vector[nsamp] Xdryfield;
	vector[nsamp] Xbuiltup;
    vector[ns] habitat; 
} 

parameters {
	real<lower=0,upper=pi()/2> sigmaa_unif;
	real<lower=0,upper=pi()/2> sigmas_unif;
	real<lower=0,upper=1> a0_unif[5];
    real<lower=0,upper=1> ah_unif[2];
	matrix[5,ns] beta_stdn;
	vector<lower=0,upper=1>[ns] beta0_unif;
	vector[nregion] ranef_a;
    vector[nsettlement] ranef_s;
	cholesky_factor_corr[2] L_Omega;
    vector<lower=0,upper=pi()/2>[5] sigmab_unif;
} 

transformed parameters {
	matrix[2,2] rho;
	matrix[ns,5] mu_beta;
    vector<lower=0>[5] sigmab;
    vector[ns] beta3;
    vector[ns] beta4;
    vector[ns] beta5;
    vector[nsamp] logitp;
    real<lower=0> sigmaa;
    real<lower=0> sigmas;
	matrix[ns,2] beta;
	vector[ns] beta0;
	real a0[5];
    real ah[2];

    for(k in 1:5) sigmab[k]=5*tan(sigmab_unif[k]);
    sigmaa = 5*tan(sigmaa_unif);
    sigmas = 5*tan(sigmas_unif);    
    rho = L_Omega * L_Omega';
    for(k in 1:5) a0[k] = 10 * inv_Phi(a0_unif[k]);
    for(k in 1:2) ah[k] = 10 * inv_Phi(ah_unif[k]);
    for(i in 1:ns) beta0[i] = 10*inv_Phi(beta0_unif[i]);
	for(i in 1:ns){
		for(j in 1:2){
			mu_beta[i,j] = a0[j]+ah[j]*habitat[i];
		}
		for(j in 3:5){
			mu_beta[i,j] = a0[j];
		}
        
	}
    beta = mu_beta[,1:2] + (diag_pre_multiply(sigmab[1:2],L_Omega)*beta_stdn[1:2,])';
    beta3 = mu_beta[,3] + sigmab[3]*beta_stdn[3,]';
    beta4 = mu_beta[,4] + sigmab[4]*beta_stdn[4,]';
    beta5 = mu_beta[,5] + sigmab[5]*beta_stdn[5,]';
    for(i in 1:nsamp) logitp[i] = beta0[speciesID[i]]+beta[speciesID[i],1]*Xtemp[i]+beta[speciesID[i],2]*Xaban_bin[i]+beta3[speciesID[i]]*Xmonth[i]+beta4[speciesID[i]]*Xdryfield[i]+beta5[speciesID[i]]*Xbuiltup[i]+ranef_s[settlementID[i]]*sigmas+ranef_a[regionID[i]]*sigmaa;
 
}

model {
	to_vector(beta_stdn)~std_normal();
	Y~bernoulli_logit(logitp);
	L_Omega~lkj_corr_cholesky(1.0); 
	ranef_a~std_normal();
	ranef_s~std_normal();

} 

generated quantities{
    vector[nsamp] log_lik;
    for(i in 1:nsamp) log_lik[i] = bernoulli_logit_lpmf(Y[i] | logitp[i]);
}
"


habitattable<-read.csv("habitat.type.csv")
habitattable<-habitattable[order(habitattable$ID),-1]
habitat.list<-colnames(habitattable)	
nhabitat<-length(habitat.list)					#生息環境リストの長さ
pars<-c("a0","ah","beta0","beta","sigmaa","sigmas","sigmab","ranef_a","ranef_s","rho")

res.habitat<-vector("list",nhabitat)
names(res.habitat)<-habitat.list

for(i in 1:nhabitat){
    habitat<-habitattable[,habitat.list[i]]
    data<-list(Y=Y,speciesID=speciesID,nsamp=nsamp,ns=ns,regionID=regionID,nregion=nregion,settlementID=settlementID,nsettlement=nsettlement,Xtemp=Xtemp,Xaban_bin=Xaban_bin,Xmonth=Xmonth,Xdryfield=Xdryfield,Xbuiltup=Xbuiltup,habitat=habitat)
    res.habitat[[i]]<-stan(model_code=stanmodelcode,model_name="abandonmentbutterfly",data=data,pars=pars,iter=2000,thin=1,chains = 3,control=list(adapt_delta=0.99,max_treedepth = 15))
}


###Ridge regression

stanmodelcode <- "
data {
	int<lower=0> ns;
    int<lower=0> nsamp;
    int<lower=0> nsettlement;
	int<lower=0> Y[nsamp];
	int<lower=0> regionID[nsamp];
	int<lower=0> nregion;
    int<lower=0> speciesID[nsamp];
    int<lower=0> settlementID[nsamp];
 	int<lower=0> nhabitat;
   	vector[nsamp] Xtemp;
	vector[nsamp] Xaban_bin;
	vector[nsamp] Xmonth;
	vector[nsamp] Xdryfield;
	vector[nsamp] Xbuiltup;
    matrix[ns,nhabitat] habitat; 
} 

parameters {
	real<lower=0,upper=pi()/2> sigmaa_unif;
	real<lower=0,upper=pi()/2> sigmas_unif;
    real<lower=0,upper=pi()/2> sigmaah_unif[2];
	real<lower=0,upper=1> a0_unif[5];
    real<lower=0,upper=1> ah_unif[2];
	matrix[5,ns] beta_stdn;
	vector<lower=0,upper=1>[ns] beta0_unif;
	vector[nregion] ranef_a;
    vector[nsettlement] ranef_s;
	cholesky_factor_corr[2] L_Omega;
    vector<lower=0,upper=pi()/2>[5] sigmab_unif;
    matrix[2,nhabitat] ah_stdn;
} 

transformed parameters {
	matrix[2,2] rho;
	matrix[ns,5] mu_beta;
    vector<lower=0>[5] sigmab;
    vector[ns] beta3;
    vector[ns] beta4;
    vector[ns] beta5;
    vector[nsamp] logitp;
    real<lower=0> sigmaa;
    real<lower=0> sigmas;
    real<lower=0> sigmaah[2];
	matrix[ns,2] beta;
	vector[ns] beta0;
	vector[5] a0;
    matrix[2,nhabitat] ah;

    for(k in 1:5) sigmab[k]=5*tan(sigmab_unif[k]);
    sigmaa = 5*tan(sigmaa_unif);
    sigmas = 5*tan(sigmas_unif);
    for(j in 1:2) sigmaah[j] = 5*tan(sigmaah_unif[j]);    
    rho = L_Omega * L_Omega';
    for(k in 1:5) a0[k] = 10 * inv_Phi(a0_unif[k]);
    for(i in 1:ns) beta0[i] = 10*inv_Phi(beta0_unif[i]);
    for(j in 1:2){
        for(m in 1:nhabitat){
            ah[j,m] = ah_stdn[j,m]*sigmaah[j];
        }
    }
	for(i in 1:ns){
		for(j in 1:2){
			mu_beta[i,j] = a0[j]+ah[j,]*habitat[i,]';
		}
		for(j in 3:5){
			mu_beta[i,j] = a0[j];
		}
        
	}
    beta = mu_beta[,1:2] + (diag_pre_multiply(sigmab[1:2],L_Omega)*beta_stdn[1:2,])';
    beta3 = mu_beta[,3] + sigmab[3]*beta_stdn[3,]';
    beta4 = mu_beta[,4] + sigmab[4]*beta_stdn[4,]';
    beta5 = mu_beta[,5] + sigmab[5]*beta_stdn[5,]';
    for(i in 1:nsamp) logitp[i] = beta0[speciesID[i]]+beta[speciesID[i],1]*Xtemp[i]+beta[speciesID[i],2]*Xaban_bin[i]+beta3[speciesID[i]]*Xmonth[i]+beta4[speciesID[i]]*Xdryfield[i]+beta5[speciesID[i]]*Xbuiltup[i]+ranef_s[settlementID[i]]*sigmas+ranef_a[regionID[i]]*sigmaa;
 
}

model {
    to_vector(ah_stdn)~ std_normal();
	to_vector(beta_stdn)~std_normal();
	Y~bernoulli_logit(logitp);
	L_Omega~lkj_corr_cholesky(1.0); 
	ranef_a~std_normal();
	ranef_s~std_normal();

} 

generated quantities{
    vector[nsamp] log_lik;
    for(i in 1:nsamp) log_lik[i] = bernoulli_logit_lpmf(Y[i] | logitp[i]);
}
"

pars<-c("a0","ah","beta0","beta","sigmaa","sigmas","sigmab","sigmaah","ranef_a","ranef_s","rho")
data<-list(Y=Y,speciesID=speciesID,nsamp=nsamp,ns=ns,regionID=regionID,nregion=nregion,settlementID=settlementID,nsettlement=nsettlement,Xtemp=Xtemp,Xaban_bin=Xaban_bin,Xmonth=Xmonth,Xdryfield=Xdryfield,Xbuiltup=Xbuiltup,habitat=habitattable,nhabitat=nhabitat)
res.ridge<-stan(model_code=stanmodelcode,model_name="abandonmentbutterfly",data=data,pars=pars,iter=2000,thin=1,chains = 3,control=list(adapt_delta=0.99,max_treedepth = 15))


####Harvest the results

library(export)
library(cowplot)
library(lemon)

spnames<-read.csv("spnames.csv",fileEncoding="UTF-8")

#Summary of estimates from the baseline model
summbaseline<-summary(res.baseline,pars=c("a0","beta0","beta","beta3","beta4","beta5","sigmaa","sigmas","sigmab","rho"))[[1]]

coefaban<-summbaseline[grep("^beta\\\\[.*2\\\\]$",rownames(summbaseline)),c("mean","2.5%","97.5%")]
colnames(coefaban)<-sub("%","",paste("aban_",colnames(coefaban),sep=""))
signifaban<-apply(coefaban[,c(2,3)],1,prod)>0

coeftemp<-summbaseline[grep("^beta\\\\[.*1\\\\]$",rownames(summbaseline)),c("mean","2.5%","97.5%")]
colnames(coeftemp)<-sub("%","",paste("temp_",colnames(coeftemp),sep=""))
signiftemp<-apply(coeftemp[,c(2,3)],1,prod)>0

coefmonth<-summbaseline[grep("^beta3",rownames(summbaseline)),c("mean","2.5%","97.5%")]
colnames(coefmonth)<-sub("%","",paste("month_",colnames(coefmonth),sep=""))
signifmonth<-apply(coefmonth[,c(2,3)],1,prod)>0

coefdrycrop<-summbaseline[grep("^beta4",rownames(summbaseline)),c("mean","2.5%","97.5%")]
colnames(coefdrycrop)<-sub("%","",paste("drycrop_",colnames(coefdrycrop),sep=""))
signifdrycrop<-apply(coefdrycrop[,c(2,3)],1,prod)>0

coefbuiltup<-summbaseline[grep("^beta5",rownames(summbaseline)),c("mean","2.5%","97.5%")]
colnames(coefbuiltup)<-sub("%","",paste("builtup_",colnames(coefbuiltup),sep=""))
signifbuiltup<-apply(coefbuiltup[,c(2,3)],1,prod)>0

intercept<-summbaseline[grep("^beta0",rownames(summbaseline)),c("mean","2.5%","97.5%")]
colnames(intercept)<-sub("%","",paste("intercept_",colnames(intercept),sep=""))


coeftable0<-cbind(spnames[,c("ID","family","subfamily","sciname","namejp")]%>%dplyr::arrange(ID),
			coefaban,signifaban,coeftemp,signiftemp,
			coefdrycrop,signifdrycrop,coefbuiltup,signifbuiltup,
			coefmonth,signifmonth,
			intercept)


coeftable<-coeftable0%>%as_tibble%>%dplyr::arrange(sciname)%>%dplyr::arrange(subfamily)%>%dplyr::arrange(family)%>%select(-ID,-subfamily)
sum((coeftable$aban_mean<0)&coeftable$signifaban)
sum((coeftable$aban_mean>0)&coeftable$signifaban)

coeftable%>%arrange(aban_mean)%>%select(sciname)
coeftable%>%arrange(desc(aban_mean))%>%select(sciname)

coeftableout<-coeftable%>%mutate_if(is.double,signif,digits=3)%>%mutate_if(is.double,formatC,digits=3,format="g",flag="#")

#95% equiprobable circle
mubeta<-summbaseline[grep("^a0\\\\[",rownames(summbaseline),"mean")][1:2]
sigmab<-diag(summbaseline[grep("^sigmab",rownames(summbaseline),"mean")][1:2])
rho<-matrix(summbaseline[grep("^rho",rownames(summbaseline),"mean")],2,2)
covmat<-sigmab%*%rho%*%sigmab

equiprob<-function(phi,p,mubeta,sigma,rho){ 
	r=sqrt(-2*(1-rho^2)*log(1-p)/(1-2*rho*sin(phi)*cos(phi)))
	x<-mubeta[1]+sigma[1]*r*cos(phi)
	y<-mubeta[2]+sigma[2]*r*sin(phi)
	return(data.frame(phi,x,y))
}

eq<-equiprob(seq(0,2*pi,0.01),0.95,mubeta,diag(sigmab),rho[1,2])

#principal component
eigenvec<-eigen(covmat)$vector[,1,drop=T]
slope<-sqrt(eigenvec[2]/eigenvec[1])
intercept<-mubeta[2]-mubeta[1]*slope

#Fig.2
coefplot<-ggplot(coeftable,aes(y=aban_mean,x=temp_mean))+
	geom_hline(yintercept=0)+
	geom_vline(xintercept=0)+
	geom_path(eq,mapping=aes(y=y,x=x),linetype="dashed")+
	geom_abline(slope=slope,intercept=intercept,linetype="dashed")+
	geom_errorbar(aes(ymin=aban_2.5,ymax=aban_97.5),col=gray(0.5,alpha=0.5))+
	geom_errorbarh(aes(xmin=temp_2.5,xmax=temp_97.5),col=gray(0.5,alpha=0.5))+
	geom_point()+
	theme_cowplot()+
	xlab("Coefficient of temperature")+
	ylab("Coefficient of abandonment")

coefplot

###Table S5

table.habitat<-res.habitat%>%lapply(summary)%>%lapply("[[",1)%>%sapply("[","ah[2]",)%>%t
signif.habitat<-apply(table.habitat[,c(4,8)],1,prod)>0
habitat.type<-c("forest","farmland","garden","built-up area","forest edge","grassland","riverside","open forest","alpine area","wetland","rocky area")
table.habitat<-data.frame(habitat.type=habitat.type,mean=table.habitat[,c(1)],signif.habitat,table.habitat[,c(4,8)])
table.habitat<-table.habitat%>%mutate_if(is.double,signif,digits=3)%>%mutate_if(is.double,formatC,digits=3,format="g",flag="#")



###Fig.3
sample.habitat<-res.habitat%>%lapply(extract,pars=c("a0[2]","ah[2]"))%>%lapply(as.data.frame)
sample.habitat2<-sample.habitat%>%lapply(apply,1,cumsum)
habitat.type2<-data.frame(names(sample.habitat2),habitat.type)
colnames(habitat.type2)<-c("habitat_type","habitat_type2")
habitat.mean<-sample.habitat2%>%
	sapply(apply,1,mean)%>%
	data.frame(treatment=c(0,1))%>%
	tidyr::pivot_longer(-treatment,names_to="habitat_type",values_to="mean")%>%
	left_join(habitat.type2,by="habitat_type")

habitat.lci<-sample.habitat2%>%
	sapply(apply,1,quantile,probs=0.025)%>%
	data.frame(treatment=c(0,1))%>%
	tidyr::pivot_longer(-treatment,names_to="habitat_type",values_to="ci2.5")

habitat.uci<-sample.habitat2%>%
	sapply(apply,1,quantile,probs=0.975)%>%
	data.frame(treatment=c(0,1))%>%
	tidyr::pivot_longer(-treatment,names_to="habitat_type",values_to="ci97.5")

habitat.mean<-cbind(habitat.mean,ci2.5=habitat.lci$ci2.5,ci97.5=habitat.uci$ci97.5)


betaname<-names(res.habitat[[1]])[grep("^beta\\\\[.*\\\\,2\\\\]$",names(res.habitat[[1]]))]
sample.species<-res.habitat%>%lapply(extract,pars=betaname)%>%lapply(as.data.frame)

species.mean<-sample.species%>%
		lapply(apply,2,mean)%>%
		as.data.frame
rownames(species.mean)<-colnames(Y)
species.mean<-species.mean%>%
		rownames_to_column("name")%>%
		tidyr::pivot_longer(-name,names_to="habitat_type",values_to="mean")%>%
		left_join(habitat.type2,by="habitat_type")


habitat_long<-habitattable%>%
		rownames_to_column("name")%>%
		tidyr::pivot_longer(-name,names_to="habitat_type",values_to="treatment")

species.mean<-cbind(species.mean,treatment=habitat_long$treatment)

habitatplot<-ggplot(species.mean%>%dplyr::filter(is.element(habitat_type,names(signif.habitat)[signif.habitat])),aes(y=mean,x=as.ordered(treatment)))+
	geom_hline(yintercept=0)+
	geom_point(alpha=0.5,col=gray(0.5))+
	geom_errorbar(data=habitat.mean%>%dplyr::filter(is.element(habitat_type,names(signif.habitat)[signif.habitat])),aes(x=treatment+1+0.1,ymin=ci2.5,ymax=ci97.5),width = 0.2,size=1)+
	geom_point(data=habitat.mean%>%dplyr::filter(is.element(habitat_type,names(signif.habitat)[signif.habitat])),aes(x=treatment+1+0.1,y=mean),cex=2)+
	facet_rep_wrap(~habitat_type2,repeat.tick.labels = TRUE,strip.position = "bottom")+
	theme_cowplot()+
	ylab("Coefficients of abandonment")+
	xlab("")+
	theme(strip.placement = "outside",strip.background = element_blank(),strip.text.x = element_text(size = 14))
	

habitatplot

###CG and CL
habitat.proj<-read.csv("habitat.proj.csv")
bpaprojlong<-read.csv("presence.proj.csv")	#presence.proj_maskRL.csv: a dataset with red list species masked

spnames.proj<-read.csv("spnames.proj.csv",fileEncoding="UTF-8")

bpaproj<-bpaprojlong%>%pivot_wider(names_from=projID,values_from=presence)
bpaprojmat<-bpaproj%>%select(-code,-lat,-long)%>%as.matrix(dimnames=NULL)

bpaprojmatRL<-bpaprojmat[,spnames.proj$RLrank!=""]


MCsamp<-extract(res.ridge,pars=c("a0","ah"))
a0mean<-mean(MCsamp$a0[,2])
ahmean<-apply(MCsamp$ah[,2,],2,mean)
names(ahmean)<-colnames(habitattable)
ahmean<-ahmean[colnames(habitat.proj)[-1]]
nsamp<-nrow(MCsamp$a0)



ahproj<-a0mean*rep(1,nrow(habitat.proj))+as.matrix(habitat.proj[,-1])%*%ahmean
ahproj.pos<-pmax(ahproj,0)
ahproj.neg<-pmin(ahproj,0)

CG.mean<-c(bpaprojmat%*%ahproj.pos)
CL.mean<-c(bpaprojmat%*%ahproj.neg)*(-1)

CG.RL.mean<-c(bpaprojmatRL%*%ahproj.pos[spnames.proj$RLrank!=""])
CL.RL.mean<-c(bpaprojmatRL%*%ahproj.neg[spnames.proj$RLrank!=""])*(-1)

#standard deviation
a0samp<-MCsamp$a0[,2]
ahsamp<-MCsamp$ah[,2,]

colnames(ahsamp)<-colnames(habitattable)
ahsamp<-ahsamp[,colnames(habitat.proj)[-1]]

ahproj<-outer(a0samp,rep(1,nrow(habitat.proj)))+ahsamp%*%t(as.matrix(habitat.proj[,-1]))
ahproj.pos<-pmax(ahproj,0)
ahproj.neg<-pmin(ahproj,0)

CG.sd1<-CG.sd2<-CL.sd1<-CL.sd2<-CG.RL.sd1<-CG.RL.sd2<-CL.RL.sd1<-CL.RL.sd2<-numeric(nrow(bpaprojmat))
for(i in 1:nsamp){
    cat(i,"\\n")
    CG.sd1<-CG.sd1+c(bpaprojmat%*%ahproj.pos[i,])^2
    CG.sd2<-CG.sd2+c(bpaprojmat%*%ahproj.pos[i,])

    CL.sd1<-CL.sd1+c(bpaprojmat%*%ahproj.neg[i,])^2
    CL.sd2<-CL.sd2+c(bpaprojmat%*%ahproj.neg[i,])

    CG.RL.sd1<-CG.RL.sd1+c(bpaprojmatRL%*%ahproj.pos[i,spnames.proj$RLrank!=""])^2
    CG.RL.sd2<-CG.RL.sd2+c(bpaprojmatRL%*%ahproj.pos[i,spnames.proj$RLrank!=""])

    CL.RL.sd1<-CL.RL.sd1+c(bpaprojmatRL%*%ahproj.neg[i,spnames.proj$RLrank!=""])^2
    CL.RL.sd2<-CL.RL.sd2+c(bpaprojmatRL%*%ahproj.neg[i,spnames.proj$RLrank!=""])
}


CG.sd<-sqrt(CG.sd1/(nsamp-1)-(CG.sd2/(nsamp-1))^2)
CL.sd<-sqrt(CL.sd1/(nsamp-1)-(CL.sd2/(nsamp-1))^2)

CG.RL.sd<-sqrt(CG.RL.sd1/(nsamp-1)-(CG.RL.sd2/(nsamp-1))^2)
CL.RL.sd<-sqrt(CL.RL.sd1/(nsamp-1)-(CL.RL.sd2/(nsamp-1))^2)


projtable<-tibble(bpaproj[,1:3],CG.mean,CG.sd,CL.mean,CL.sd,CG.RL.mean,CG.RL.sd,CL.RL.mean,CL.RL.sd)



#################Fig.4 hexiagonal bin
library(sf)
library(scales)

mesh3elev<-st_read(dsn = getwd(), layer = "mesh3Relev",options = "ENCODING=UTF-8")
colnames(mesh3elev)[1]<-"code"
mesh3elev<-mesh3elev%>%mutate(code=as.numeric(code))
mesh3<-mesh3elev%>%right_join(projtable,by="code")

ggplot(mesh3,aes(x=elev_mean,y=CL.mean))+
	geom_hex()+
	scale_fill_gradient(low="gray",high=muted("blue"),trans="log10")+
	theme_cowplot()+
	xlab("Elevation (m)")+
	ylim(c(0,25))+
	ylab("CL")

ggplot(mesh3,aes(x=elev_mean,y=CG.mean))+
	geom_hex()+
	scale_fill_gradient(low="gray",high=muted("red"),trans="log10")+
	theme_cowplot()+
	xlab("Elevation (m)")+
	ylim(c(0,25))+
	ylab("CG")


ggplot(mesh3,aes(x=elev_mean,y=CL.RL.mean))+
	geom_hex()+
	scale_fill_gradient(low="gray",high=muted("blue"),trans="log10")+
	theme_cowplot()+
	xlab("Elevation (m)")+
	ylim(c(0,7))+
	ylab("CL")

ggplot(mesh3,aes(x=elev_mean,y=CG.RL.mean))+
	geom_hex()+
	scale_fill_gradient(low="gray",high=muted("red"),trans="log10")+
	theme_cowplot()+
	xlab("Elevation (m)")+
	ylim(c(0,7))+
	ylab("CG")


