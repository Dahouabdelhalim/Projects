# Autor: Koyo Nakamura
# Date: 3/9/2020

library("rstan")
library("brms")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("Nakamura_Watanabe_2020.Rdata")
#################
# Experiment 1
#################

# Validation of sexual dimorphism dimension (Best-fitted model)

model1 <- brm(Rating ~ Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces +
                (Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces|Subject_ID) +
                (Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces|Face_ID),
              data=dat1,family="gaussian", chains=4, iter=13000, warmup=3000, thin=1, seed=1, control = list(adapt_delta = 0.99, max_treedepth = 15))

print(summary(model1))


# Correlation between perceived masculinity/femininity and attractiveness

cor_data <- list(N=200,male_data=cor_data_male[,c(2,3)],female_data=cor_data_female[,c(2,3)])
par = c("rho1","rho2","delta_r")

cor_anal <- '
data{
	int<lower=0>  N; 
	vector[2] male_data[N];
	vector[2] female_data[N];
}
parameters{
	vector[2] mu1;
	vector<lower=0>[2] sigma1;
	real<lower=-1,upper=1> rho1;
	vector[2] mu2;
	vector<lower=0>[2] sigma2;
	real<lower=-1,upper=1> rho2;
}
transformed parameters{
	vector<lower=0>[2] sig21;
	matrix[2,2] Sigma1;
	vector<lower=0>[2] sig22;
	matrix[2,2] Sigma2;

	sig21[1] = pow(sigma1[1],2);
	sig21[2] = pow(sigma1[2],2);
	Sigma1[1,1] = sig21[1];
	Sigma1[2,2] = sig21[2];
	Sigma1[1,2] = sigma1[1]*sigma1[2]*rho1;
	Sigma1[2,1] = sigma1[1]*sigma1[2]*rho1;

	sig22[1] = pow(sigma2[1],2);
	sig22[2] = pow(sigma2[2],2);
	Sigma2[1,1] = sig22[1];
	Sigma2[2,2] = sig22[2];
	Sigma2[1,2] = sigma2[1]*sigma2[2]*rho2;
	Sigma2[2,1] = sigma2[1]*sigma2[2]*rho2;
}
model{
	for(i in 1:N){
		male_data[i] ~ multi_normal(mu1,Sigma1);
		female_data[i] ~ multi_normal(mu2,Sigma2);
	}
}
generated quantities{
	real delta_r;
	delta_r = rho2 - rho1;
}
';

fit = stan(model_code = cor_anal, data = cor_data,chains=4, iter=13000, warmup=3000,seed=1,cores=4, control = list(adapt_delta = 0.99, max_treedepth = 15))
print(fit, pars=par,digits_summary=2)


#################
# Experiment 2
#################

# Attractiveness orthogonal to sexual dimorphism
# Attractiveness rating scores (Best-fitted model)

model2 <- brm(Rating ~ Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces +
                (Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces|Subject_ID) +
                (Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces|Face_ID),
              data=dat2,family="gaussian", chains=4, iter=13000, warmup=3000, thin=1, seed=1, control = list(adapt_delta = 0.99, max_treedepth = 15))

print(summary(model2))


# Attractiveness orthogonal to sexual dimorphism
# Sexual dimorphism rating scores (Best-fitted model)

model3 <- brm(Rating ~ Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces +
                (Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces|Subject_ID) +
                (Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces|Face_ID),
              data=dat3,family="gaussian", chains=4, iter=13000, warmup=3000, thin=1, seed=1, control = list(adapt_delta = 0.99, max_treedepth = 15))

print(summary(model3))


# Sexual dimorphism orthogonal to attractiveness
# Attractiveness rating scores (Best-fitted model)

model4 <- brm(Rating ~ Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces +
                (Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces|Subject_ID) +
                (Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces|Face_ID),
              data=dat4,family="gaussian", chains=4, iter=13000, warmup=3000, thin=1, seed=1, control = list(adapt_delta = 0.99, max_treedepth = 15))

print(summary(model4))


# Sexual dimorphism orthogonal to attractiveness
# Sexual dimorphism rating scores (Best-fitted model)

model5 <- brm(Rating ~ Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces +
                (Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces|Subject_ID) +
                (Face_Exaggeration + Sex_of_Faces + Face_Exaggeration:Sex_of_Faces + I(Face_Exaggeration^2) + I(Face_Exaggeration^2):Sex_of_Faces|Face_ID),
              data=dat5,family="gaussian", chains=4, iter=13000, warmup=3000, thin=1, seed=1, control = list(adapt_delta = 0.99, max_treedepth = 15))

print(summary(model5))


# Comparison of the effect of face exaggeration on perceived attractiveness and sexual dimorphism

bs_model2 = brms::posterior_samples(model2, "^b")
bs_model3 = brms::posterior_samples(model3, "^b")
bs_model4 = brms::posterior_samples(model4, "^b")
bs_model5 = brms::posterior_samples(model5, "^b")


# β attractiveness orthogonal to sexual dimorphism - β sexual dimorphism orthogonal to attractiveness
b_dif1 = bs_model2$b_Face_Exaggeration - bs_model4$b_Face_Exaggeration
b_dif2 = bs_model2$b_IFace_ExaggerationE2 - bs_model4$b_IFace_ExaggerationE2

print(c(EAP=mean(b_dif1),quantile(b_dif1, c(0.025 ,0.975))))
print(c(EAP=mean(b_dif2),quantile(b_dif2, c(0.025 ,0.975))))

# β sexual dimorphism orthogonal to attractiveness - β attractiveness orthogonal to sexual dimorphism
b_dif3 = bs_model5$b_Face_Exaggeration - bs_model3$b_Face_Exaggeration
b_dif4 = bs_model5$b_IFace_ExaggerationE2 - bs_model3$b_IFace_ExaggerationE2

print(c(EAP=mean(b_dif3),quantile(b_dif3, c(0.025 ,0.975))))
print(c(EAP=mean(b_dif4),quantile(b_dif4, c(0.025 ,0.975))))

