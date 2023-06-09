# Functions
logit <- function(p) log(p/(1-p))


# Libraries
library(lme4)
library(TMB)

# Compiling c++ code
compile( "eira_adaptdyn.cpp" )
dyn.load(dynlib("eira_adaptdyn") )

# DATA
DS <- read.csv("Discharge.txt", sep = "\\t")
discharge <- log(DS$DC_Jun_Sept[DS$Year >= (1940-8) & !is.na(DS$Year)]) 
n <- length(discharge)

d <- read.csv("Atlantic salmon of the river Eira 1925-2016.csv", sep = ";")
lm_coef <- coef(summary(lm(log(d$Mass) ~ -1 + factor(Year.of.catch), data = d)))
x<-tapply(log(d$Mass), d$Year.of.catch, mean)
d_yr1 <- data.frame(year = c(as.numeric(names(x))),
                   n = c(tapply(d$Mass, d$Year.of.catch, length)),
                   ln_mass = lm_coef[,1],
                   ln_mass_SE = lm_coef[,2]
)
d_yr1 <- d_yr1[order(d_yr1$year),]
d_yr <- subset(d_yr1, year>=1940)

d$tot.age <- d$Smolt.age+d$Sea.age
n_age <- tapply(d$tot.age, d$tot.age, length)[1:5]
n_age[5] <- sum(tapply(d$tot.age, d$tot.age, length)[5:9])
p_age <- n_age/sum(n_age)


afreq <- data.frame(year = c(1987, 2016),
           Avgll3 = c(137, 127),
           Tvgll3 = c(238, 264),
           Asix6 = c(195, 111),
           Tsix6 = c(238, 288))


# Genetic effects (Table S3)
b <- c(0.43, 0.27, 0.66, -0.25, 0.33, 0.28) 
V <- matrix(c( 0.25e-2,  0.01e-2,        0,        0, -0.11e-2,        0,
               0.01e-2,  0.59e-2,        0,        0, -0.04e-2,        0,
                     0,        0,  0.61e-2,        0,        0, -0.24e-2,
                     0,        0,        0,  1.22e-2,        0, -0.20e-2,
              -0.11e-2, -0.04e-2,        0,        0,  0.32e-2,        0,
                     0,        0, -0.24e-2, -0.20e-2,        0,  1.30e-02),
            ncol = 6, nrow = 6) 

x <- match(1940:2016, d_yr$year)
Data = list(n_yr= n,
            p_alder = as.vector(p_age),
            discharge = discharge, 
            z_observed = c(rep(NA, 8), d_yr$ln_mass[x]), # first values correspond to the 8 years before the time series
            Vp = 0.31,  # VarCorr(lmer(log(Mass)~1+(1|Year.of.catch), data = d, subset = Year.of.catch < 1954))
            SE = c(rep(NA, 8), d_yr$ln_mass_SE[x]),
            b = b, 
            V = V, 
            A_vgll3 = c(rep(NA, 8), afreq$Avgll3[match(1940:2016, afreq$year)]),
            T_vgll3 = c(rep(NA, 8), afreq$Tvgll3[match(1940:2016, afreq$year)]),
            A_six6 = c(rep(NA, 8), afreq$Asix6[match(1940:2016, afreq$year)]),
            T_six6 = c(rep(NA, 8), afreq$Tsix6[match(1940:2016, afreq$year)]),
            Prior_q = c(-2,1),
            Prior_logit_h2 = c(0,1)
)

#plot(Data$discharge, Data$z)
Data[["z_observed"]][is.na(Data[["z_observed"]])] <- 999 # this is NA
Data[["SE"]][is.na(Data[["SE"]])] <- 999 
Data[["A_vgll3"]][is.na(Data[["A_vgll3"]])] <- 999 
Data[["T_vgll3"]][is.na(Data[["T_vgll3"]])] <- 999 
Data[["A_six6"]][is.na(Data[["A_six6"]])] <- 999 
Data[["T_six6"]][is.na(Data[["T_six6"]])] <- 999 

# PARAMETERS  
Parameters = list("z_origin" = rep(mean(log(d$Mass)[d$Year.of.catch%in%c(1940:1953)]), 8), 
                  "logit_p1_0" = logit(0.72),
                  "logit_p2_0" = logit(0.95),#rep(logit(0.95), 8),
                  "q"=0, 
                  "Opt0"=9.3, 
                  "b0"=0,
                  "b1"=0,
                  "logit_h2"=0, 
                  "logit_Ps_before"= 1.4, 
                  "lnSigma_e" = 0,
                  "e" = rep(0, n)
)  

Random = c("z_origin", "logit_p1_0", "logit_p2_0", "q", "logit_h2", "e") 
Map = list(b1 = factor(NA)) # setting the quadratic for the optimum relationship to zero 
          
# Fitting model
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, map=Map,
                DLL = "eira_adaptdyn")
Opt = nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr)
SD = sdreport(Obj)

# Results
ranef <- summary(SD, "random")
report <- summary(SD, "report")

# Model output
dd <- data.frame(year = 1925:2016,
                 z_observed =c(d_yr1$ln_mass[1:2], rep(NA, 5+8),d_yr$ln_mass[x]),
                 z = c(rep(NA, 7+8), report[rownames(report)%in%"z",1][-c(1:8)]),
                 z_SE = c(rep(NA, 7+8), report[rownames(report)%in%"z",2][-c(1:8)]),
                 z_star = c(rep(NA, 7+8), report[rownames(report)%in%"z_star",1][-c(1:8)]),
                 opt = c(rep(NA, 7), report[rownames(report)%in%"Opt",1]),
                 opt_SE = c(rep(NA, 7), report[rownames(report)%in%"Opt",2]),
                 p1 = c(rep(NA, 7+8), report[rownames(report)%in%"p1",1][-c(1:8)]),
                 p1_SE = c(rep(NA, 7+8), report[rownames(report)%in%"p1",2][-c(1:8)]),
                 p1_delta = c(rep(NA, 7+8), report[rownames(report)%in%"delta_p1",1][-c(1:8)]),
                 p2 = c(rep(NA, 7+8), report[rownames(report)%in%"p2",1][-c(1:8)]),
                 p2_SE = c(rep(NA, 7+8), report[rownames(report)%in%"p2",2][-c(1:8)]),
                 p2_delta = c(rep(NA, 7+8), report[rownames(report)%in%"delta_p2",1][-c(1:8)]),
                 Va = c(rep(NA, 7+8), report[rownames(report)%in%"Va",1][-c(1:8)]),
                 Va_SE = c(rep(NA, 7+8), report[rownames(report)%in%"Va",2][-c(1:8)]),
                 Va_vgll3 = c(rep(NA, 7+8), report[rownames(report)%in%"Va_vgll3",1][-c(1:8)]),
                 Va_vgll3_SE = c(rep(NA, 7+8), report[rownames(report)%in%"Va_vgll3",2][-c(1:8)]),
                 Va_six6 = c(rep(NA, 7+8), report[rownames(report)%in%"Va_six6",1][-c(1:8)]),
                 Va_six6_SE = c(rep(NA, 7+8), report[rownames(report)%in%"Va_six6",1][-c(1:8)]),
                 Beta = c(rep(NA, 7+8), report[rownames(report)%in%"Beta",1][-c(1:8)]),
                 Beta_SE = c(rep(NA, 7+8), report[rownames(report)%in%"Beta",2][-c(1:8)]),
                 Discharge = c(rep(NA, 7),discharge),
                 seldiff = c(rep(NA, 7+8), report[rownames(report)%in%"seldiff",1][-c(1:8)]),
                 seldiff_SE = c(rep(NA, 7+8), report[rownames(report)%in%"seldiff",2][-c(1:8)]),
                 Load = c(rep(NA, 7), report[rownames(report)%in%"Load",1]),
                 Load_SE = c(rep(NA, 7), report[rownames(report)%in%"Load",2])
)
dd$year <- as.numeric(dd$year)
#dd[,c(2,3,5,6)] <- dd[,c(2,3,5,6)]-log(1000) #converting to kg

dd <- subset(dd, dd$year %in% 1940:2016)

write.csv(dd, file = "model_output.csv")
