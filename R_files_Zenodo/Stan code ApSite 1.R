library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("~/Application Set site 1.R")

Nplots <- 16
mu_prior=rep(0,3)

# Application site 1, batch 1
lm11 <- c("Npars", "dr", "hr", "hx", "As11dr", "As11hr", "As11hx", "Ndata1", "Ndata11", "Nplots","subj","mu_prior")


Aps11_init<-list(
  list( a=c(1.2, 0.6, 3, 0.6, 0.4, 60, 4),
        sigma_dr=0.05,
        mu_dr=dr,
        sigma_u=c(0.01,0.1,1),
        u=matrix(0.01, ncol=3, nrow=16),
        b=c(1.2, 0.6, 3.0, 0.6, 0.4, 60, 4),
        sigma_As11dr=0.2,
        mu_As11dr=As11dr ),
  list( a=c(1.1, 0.5, 2, 0.7, 0.3, 70, 4.5),
        sigma_dr=0.05,
        mu_dr=dr,
        sigma_u=c(0.01,0.1,1),
        u=matrix(0.01, ncol=3, nrow=16),
        b=c(1.1, 0.5, 2.0, 0.7, 0.3, 70, 4.5),
        sigma_As11dr=0.2,
        mu_As11dr=As11dr )
)


As11mod <- "
data {
int<lower=1> Ndata1;                       // Number of data points in the fit set
real hx[Ndata1];
real hr[Ndata1];
real dr[Ndata1];
int<lower=0> Npars;                        // Number of parameters in the model (Npars=7)
int<lower=1> Nplots;
int<lower=1, upper=Nplots> subj[Ndata1];   // Subject id
vector[3] mu_prior;
int<lower=0> Ndata11;                      // Number of data points in the calibration set
real As11hx[Ndata11];
real As11hr[Ndata11];
real As11dr[Ndata11];
}

parameters {
vector<lower=0>[Npars] a;
real<lower=0> sigma_dr;
vector<lower=0>[3] sigma_u;
cholesky_factor_corr[3] L_u;
vector[3] z_u[Nplots];
vector<lower=0>[Npars] b;
real<lower=0> sigma_As11dr;
}

transformed parameters {
vector[3] u[Nplots];
vector[Npars] beta;
vector[Nplots] v1;
vector[Nplots] v2;
vector[Nplots] v3;

{
  matrix[3,3] Sigma_u;
  Sigma_u <- diag_pre_multiply(sigma_u, L_u); // Random effects
  for(j in 1:Nplots)
  u[j] <- Sigma_u*z_u[j];
}
for (n in 1:Nplots){
v1[n] <- u[n, 1];
v2[n] <- u[n, 2];
v3[n] <- u[n, 3];
}

beta<-a;
beta[4] <- a[4] + mean(v1);
beta[5] <- a[5] + mean(v2);
beta[6] <- a[6] + mean(v3);    

}

model {
real mu_dr[Ndata1];
real mu_As11dr[Ndata11];
L_u ~ lkj_corr_cholesky(1.0);

b ~ normal(mean(beta), sd(beta));
sigma_dr ~ cauchy(0, 2.5);
sigma_As11dr ~ cauchy(0, 2.5);

for(j in 1:Nplots) z_u[j] ~ normal(0,1);

for (i in 1:Ndata1){
mu_dr[i] <- (a[1] - (hx[i] * a[2]) / (hx[i] + a[3]))
* ((1 - (a[4]+u[subj[i],1]) * hr[i]) * (1 + (a[5]+u[subj[i],2]) * exp(-(a[6]+u[subj[i],3]) * hr[i])) - (1 - (a[4]+u[subj[i],1])) * pow(hr[i], a[7]));
}
for(k in 1:Ndata11){
mu_As11dr[k] <- (b[1] - (As11hx[k] * b[2]) / (As11hx[k] + b[3]))
* ((1 - b[4] * As11hr[k]) * (1 + b[5] * exp(-b[6] * As11hr[k])) - (1 - b[4]) * pow(As11hr[k], b[7]));
}

dr ~ normal(mu_dr, sigma_dr);
As11dr ~ normal(mu_As11dr, sigma_As11dr);
}
"


As11stan <- stan(model_code=As11mod, data=lm11, iter=2000, chains=2, init=Aps11_init, control=list(max_treedepth=12))