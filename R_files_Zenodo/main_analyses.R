################################################################################################################################################################
# This code is for the main analyses of the model
# from 'Hierarchical variation in phenotypic flexibility across timescales and associated survival selection shape the dynamics of partial seasonal migration' 
# by Paul Acker, Francis Daunt, Sarah Wanless, Sarah J. Burthe, Mark A. Newell, Michael P. Harris, Carrie Gunn, Robert Swann, Ana Payo-Payo, and Jane M. Reid.
################################################################################################################################################################


### Load the capture-recapture data
data  <- read.csv('dataset.csv',header=T,colClasses = "character")
head(data)


### Critical values that will be used as input for analyses of the Stan model
# Number of occasions per complete year
N_occ <- 5
# Number of years
N_year <- 11
# Total number of occasions across the study period (51): 5 occasions in the first 10 years and 1 occasion in the last (11th) year.
N_time <- N_occ*(N_year-1)+1 # this is also the dimension of the capture-recapture history matrix
# Number of observation events
N_event <- 6 # this is also the maximum value that cells of the capture-recapture history matrix can take
#  Number of unique unique couples of capture-recapture histories (1718)
N_crh <- nrow(data)

### Re-structure into the input format for the Stan model
## Matrix of unique capture-recapture histories (i.e. capture-resighting sequences)
## (rows: focal history, columns: focal capture-resighting occasion (time), value: observation
crh <- matrix(nrow=N_crh,ncol=N_time)
for(h in 1:N_crh) {
  crh[h,] <- as.numeric(strsplit(data[h,'capture_recapture_history'],'')[[1]])
}
head(crh)
## Matrix of frequencies of the unique capture-recapture histories in each sex 
## (rows: focal history, columns: sex (first: females, second: males), value: count of individuals that have the focal history)
freq <- matrix(as.numeric(unlist(data[,c('count_females','count_males')])),nrow=N_crh,ncol=2)
head(freq)

### For information, let's retrieve the data size mentioned in the main text:
# Total number of individals
print(sum(freq))
# Total number of females
print(sum(freq[,1]))
# Total number of males
print(sum(freq[,2]))

### Bundle these elements into a list (i.e. the input format for Stan)
dataset <- list(N_year  = N_year,
                N_occ   = N_occ,
                N_time  = N_time,
                N_event = N_event,
                N_crh   = N_crh,
                crh     = crh,
                freq    = freq)

### Run the analysis of the model
# This will require the package rstan
library(rstan)
rstan_options(auto_write = TRUE)
# set nb of cores
options(mc.cores = 5)
# Number of chains to be run
nc <- 5
# Number of warmup iterations
nw <- 1000
# Total number of iterations (including warmup)
ni <- 2000
# Thinning rate
nt <- 1 # 1 means that there is no thinning: we store samples from all iterations
# Parameters to be monitored
parameters <- c("pi","kappa","yphi","phi","r","em","d","im","swi","p_r","p_m")
# Call Stan
out <- stan(file = 'model_code.stan',
            data = dataset, 
            pars = parameters,
            chains = nc, warmup = nw, iter = ni, thin = nt)
# Note that this run is expected to take several days (the exact duration will vary depending on available computing power)

### Extract the full posterior of each parameters
## Store these full posteriors in arrays named after the parameter 
## Note that 1st dimension will be the iteration (all chains merged and iterations permuted)
## Other indexes are the same as in the Stan code
for(p in parameters) {
  assign(p,extract(out,p,permuted=T)[[1]])
}

### Then calculate derived quantities
# Note that we derived the quantities in R after having run the analyses of the Stan model
# However, a more elegant way would have been to derive these quantities directly in the Stan program
# This would have required appropriate code in the "generated quantities" block of the Stan code

## Because we do these calculations in R, we need to recalculate all transformed data calculated in the Stan code
# Critical values
N_trans <- N_time-1
L_occ <- 1
N_mig <- N_event-1
N_omig <- N_event-2
N_sex <- 2
N_area <- N_event;
N_omig <- N_mig-1;
N_osite <- N_omig+1;
N_C2ML <- 2+N_mig;
N_C3R <- N_C2ML+1;
N_C3M1 <- N_C3R+1;
N_C3ML <- N_C3R+N_mig;
N_C3RX <- N_C3ML+1;
N_C3MX1 <- N_C3RX+1;
N_C3MXL <- N_C3RX+N_mig;
N_state <- N_C3MXL;
N_dead <- N_state+1;
N_fullyr <- N_year-1;
N_tactic <- N_tactic <- 3;
# Time frame
occ <- vector()
year <- vector()
time <- array(dim=c(N_occ,N_year))
for (y in 1:N_year) {
  year[(N_occ*(y-1)+1):ifelse(y<N_year,(N_occ*y),N_time)] = rep(y,ifelse(y<N_year,N_occ,L_occ))
  for (o in 1:ifelse(y<N_year,N_occ,L_occ)) {
    time[o,y] = N_occ*(y-1)+o
    occ[time[o,y]] = o
  }
}
# Time counters for the parameter arrays
N_time2 <- N_time3 <- N_time4 <- 0
t2 <- t3 <- t4 <- vector()
for (t in 1:N_trans) {
  if(occ[t]==1) {
    N_time2 = N_time2 <- N_time2 + 1
    t2[time[occ[t],year[t]]] = N_time2 
  }
  else {
    N_time4 = N_time4 + 1
    t4[time[occ[t],year[t]]] = N_time4
    if(occ[t]<N_occ) {
      N_time2 = N_time2 <- N_time2 + 1
      t2[time[occ[t],year[t]]] = N_time2 
      N_time3 = N_time3 <- N_time3 + 1
      t3[time[occ[t],year[t]]] = N_time3
    }
  }
}
# Function to recompute the state-transition matrix ps
computeps <- function(pi,kappa,p_r,p_m,phi,em,d,im,swi) {
  psi_r <- array(NA,dim=c(N_sex,N_time2,N_area))
  psi_m <- array(NA,dim=c(N_sex,N_mig,N_time3,N_area))
  ps <- array(0,dim=c(N_sex,N_dead,N_trans,N_dead))
  for (s in 1:N_sex) {
    for (t in 1:N_time2) {
      psi_r[s,t,1] = 1-em[s,t];
      psi_r[s,t,2:N_area] = em[s,t]*d[t,];
    }
    for (t in 1:N_time3) {
      for (k in 1:N_mig) {
        psi_m[s,k,t,1] = im[s,k,t];
        for (l in 2:N_area) {
          psi_m[s,k,t,l] = (1-im[s,k,t])*swi/(N_mig-1);
        }
        psi_m[s,k,t,k+1] = (1-im[s,k,t])*(1-swi);
      }
    }
    for (t in 1:N_trans) {
      ps[s,1,t,N_dead] = 1.0-phi[s,1,t];
      ps[s,N_dead,t,N_dead] = 1.0;
      if (occ[t]<N_occ) {
        ps[s,1,t,1] = phi[s,1,t];
        ps[s,N_C3R,t,N_dead] = 1.0-phi[s,3,t];
        if (occ[t]<(N_occ-1)) {
          ps[s,N_C3R,t,N_C3R:N_C3ML] = phi[s,3,t]*psi_r[s,t2[t],];
        }
      }
      if (occ[t]==1) {
        ps[s,2,t,N_dead] = 1-phi[s,2,t];
        ps[s,2,t,3:N_C2ML] = phi[s,2,t]*d[t2[t],];
      }
      else {
        for (k in 1:N_mig) {
          ps[s,2+k,t,N_dead] = 1-phi[s,2,t];
          ps[s,N_C3R+k,t,N_dead] = 1-phi[s,3,t];
        }
        if (occ[t]<N_occ) {
          if (occ[t]==2) {
            for (k in 1:N_mig) {
              ps[s,N_C3R+k,t,N_C3RX:N_C3MXL] = phi[s,3,t]*psi_m[s,k,t3[t],];
            }
          }
          for (k in 1:N_mig) {
            ps[s,2+k,t,3:N_C2ML] = rep(phi[s,2,t],N_mig)*swi/(N_mig-1);
            ps[s,2+k,t,2+k] = phi[s,2,t]*(1-swi);
          }
        }
        if (occ[t]>2) {
          ps[s,N_C3RX,t,N_dead] = 1.0-phi[s,3,t];
          if (occ[t]<N_occ) {
            for (k in 1:N_mig) {
              ps[s,N_C3R+k,t,N_C3RX] = phi[s,3,t]*psi_m[s,k,t3[t],1];
              ps[s,N_C3R+k,t,N_C3M1:N_C3ML] = phi[s,3,t]*psi_m[s,k,t3[t],2:N_area];
            }
            ps[s,N_C3RX,t,N_C3RX] = phi[s,3,t];
            for (k in 1:N_mig) {
              ps[s,N_C3RX+k,t,N_dead] = 1.0-phi[s,3,t];
            }
            if (occ[t]<(N_occ-1)) {
              for (k in 1:N_mig) {
                ps[s,N_C3RX+k,t,N_C3RX:N_C3MXL] = phi[s,3,t]*psi_m[s,k,t3[t],];
              }
            }
            else {
              ps[s,N_C3R,t,N_C3M1:N_C3ML] = phi[s,3,t]*d[t2[t],];
              for (k in 1:N_mig) {
                ps[s,N_C3RX+k,t,N_C3RX] = phi[s,3,t];
              }
            }
          }
          else {
            ps[s,1,t,1:2] = kappa[s,1,year[t],1:2]*phi[s,1,t];
            ps[s,1,t,N_C3R] = kappa[s,1,year[t],3]*phi[s,1,t];
            ps[s,N_C3RX,t,1:2] = kappa[s,3,year[t],1:2]*phi[s,3,t];
            ps[s,N_C3RX,t,N_C3R] = kappa[s,3,year[t],3]*phi[s,3,t];
            for (k in 1:N_mig) {
              ps[s,2+k,t,1:2] = kappa[s,2,year[t],1:2]*phi[s,2,t];
              ps[s,2+k,t,N_C3R] = kappa[s,2,year[t],3]*phi[s,2,t];
              ps[s,N_C3R+k,t,1:2] = kappa[s,3,year[t],1:2]*phi[s,3,t];
              ps[s,N_C3R+k,t,N_C3R] = kappa[s,3,year[t],3]*phi[s,3,t];
            }
          }
        }
      }
    }
  }
  return(ps)
}
# Use this function and do recompute ps
N_iter <- nc*(ni-nw) # total number of iterations: number of chains times monitored iterations (see run command above)
ps <- array(dim=c(N_iter,N_sex,N_dead,N_trans,N_dead))
for (it in 1:N_iter) {
  ps[it,,,,] <- computeps(pi[it,,,],
                          kappa[it,,,,],
                          p_r[it,,],
                          p_m[it,,,],
                          phi[it,,,],
                          em[it,,],
                          d[it,,],
                          im[it,,,],
                          swi[it])
}

## Other quantities that will be useful for calculations below
# Annual survival in each year (instead of year categories)
yphi_eachyear <- array(dim=c(N_iter,2,3,10))
yphi_eachyear[,,,1] <- yphi[,,,5]
for (y in c(2,3,6,7,8,10)) {
  yphi_eachyear[,,,y] <- yphi[,,,1]
}
yphi_eachyear[,,,4] <- yphi[,,,2]
yphi_eachyear[,,,5] <- yphi[,,,3]
yphi_eachyear[,,,9] <- yphi[,,,4]
# Number of new releases in each sex per year
get_F_time <- function(x) { min(which(c(x)==1)) }
F_time <- apply(crh,1,get_F_time)
new_rel <- array(0,dim=c(N_sex,N_year))
cnt <- 0
for(t in seq(1,N_time,N_occ)) {
  cnt <- cnt+1
  for (s in 1:N_sex) {
    new_rel[s,cnt] <- sum(freq[which(F_time==t),s])
  }
}
# Total number of individuals
N_ind <- sum(new_rel)
# Cumulative number of individuals released
cumrel <- array(dim=c(2,N_time)) 
for (s in 1:2) {
  cumrel[s,] <- rep(cumsum(new_rel[s,]),each=5)[1:N_time]
}
# Cumulative number of individuals released but without current year's new releases
cumrelbis <- array(dim=c(2,N_time)) 
for (s in 1:2) {
  cumrelbis[s,] <- cumrel[s,] - rep(new_rel[s,],each=5)[1:N_time]
}

## Now calculate key derived quantities: frequencies F_1, F_2 and F_3 (see main text and OSM)
# Function to calculate the proportion in each state given release in each year
getproportions <- function(pi,ps) {
  gamma <- array(0,dim=c(N_sex,N_fullyr,N_time,N_dead));
  for (s in 1:2) {
    for (yi in 1:N_fullyr) {
      stt <- seq(1,N_trans,5)[yi]
      for(ki in 1:3) {
        gamma[s,yi,stt,c(1,2,N_C3R)] = pi[s,yi,];
        for (t in (stt+1):N_time) {
          for (k in 1:N_dead) {
            acc <- rep(0,length=N_dead)
            for (j in 1:N_dead) {
              acc[j] <- gamma[s,yi,t-1,j] * ps[s,j,t-1,k];
            }
            gamma[s,yi,t,k] <- sum(acc)
          }
        }
      }
    }
  }
  return(gamma)
}
# Do this computation (Note: it takes several minutes)
gamma <- array(dim=c(N_iter,N_sex,N_fullyr,N_time,N_dead))
for (i in 1:N_iter) {
  gamma[i,,,,] <- getproportions(pi[i,,,],ps[i,,,,])
}
# Collapse all proportions per tactic (and keep the dead state)
proportions <- array(dim=c(N_iter,N_sex,N_fullyr,N_time,N_tactic+1))
for (i in 1:N_iter) {
  for (s in 1:N_sex) {
    for (yi in 1:N_fullyr) {
      proportions[i,s,yi,,1] <- gamma[i,s,yi,,1]
      proportions[i,s,yi,,2] <- apply(gamma[i,s,yi,,2:N_C2ML],1,sum)
      proportions[i,s,yi,,3] <- apply(gamma[i,s,yi,,N_C3R:N_state],1,sum)
      proportions[i,s,yi,,4] <- gamma[i,s,yi,,N_dead]
    }
  }
}
# Relative frequency of each tactic at the start of each year F_1
cohort_counts <- array(dim=dim(proportions))
year_counts <- array(dim=c(N_iter,N_sex,N_time,N_tactic+1))
F_1 <- array(dim=c(N_iter,N_sex,N_year,N_tactic))
for (i in 1:N_iter) {
  for (s in 1:N_sex) {
    for (y in 1:N_fullyr) {
      cohort_counts[i,s,y,,] <- proportions[i,s,y,,]*new_rel[s,y]
    }
    year_counts[i,s,,] <- apply(cohort_counts[i,s,,,],c(2,3),sum)/cumrel[s,]
    F_1[i,s,,] <- year_counts[i,s,seq(1,N_time,N_occ),1:3]/(1-year_counts[i,s,seq(1,N_time,N_occ),4])
  }
}
# Relative frequency of each tactic right after survival selection in each year F_2
F_2 <- array(dim=c(N_iter,N_sex,N_fullyr,N_tactic))
for (i in 1:N_iter) {
  for (s in 1:N_sex) {
    for (y in 1:(N_fullyr)) {
      for (k in 1:3) {
        F_2[i,s,y,k] <- F_1[i,s,y,k]*yphi_eachyear[i,s,k,y]
      }
      F_2[i,s,y,] <- F_2[i,s,y,]/sum(F_2[i,s,y,])
    }
  }
}
# Relative frequency of tactic after tactic switching (before entry of new individuals in the dataset)
# (we can use the above calculation of cohort counts but exclude new recruits, see below)
cohort_countsbis <- cohort_counts
year_countsbis <- array(dim=c(N_iter,N_sex,N_time,N_tactic+1))
F_3 <- array(dim=c(N_iter,N_sex,N_year,N_tactic))
for (i in 1:N_iter) {
  for (s in 1:N_sex) {
    for (y in 1:(N_fullyr)) {
      cohort_countsbis[i,s,y,which(year==y),] <- 0
    }
    year_countsbis[i,s,,] <- apply(cohort_countsbis[i,s,,,],c(2,3),sum)/cumrelbis[s,]
    F_3[i,s,,] <- year_countsbis[i,s,seq(1,N_time,5),1:3]/(1-year_countsbis[i,s,seq(1,N_time,5),4])
  }
}

## Finally, keep only the years for which results were shown in the paper (i.e. exclude the first year 2009-10)
F_1 <- F_1[,,2:10,]
F_2 <- F_2[,,2:10,]
F_3 <- F_3[,,3:10,]
dimnames(F_1) <- dimnames(F_2) <- dimnames(F_3) <- list(iterations=NULL,NULL,NULL,NULL)

## Now calculate the key quantities E_1, E_2, and E_3 (see main text and ESM)
E_1 <- array(dim=c(N_iter,N_sex,N_fullyr-1,N_tactic))
E_2 <- E_3 <- array(dim=c(N_iter,N_sex,N_fullyr-2,N_tactic))
for(i in 1:N_iter) {
  for (s in 1:N_sex) {
    for (y in 1:(N_fullyr-1)) {
      E_1[i,s,y,] <- F_2[i,s,y,] - F_1[i,s,y,]
      if (y<(N_fullyr-1)) {
        E_2[i,s,y,] <- F_3[i,s,y,] - F_2[i,s,y,]
        E_3[i,s,y,] <- F_1[i,s,y+1,] - F_3[i,s,y,]
      }
    }
  }
}
dimnames(E_1) <- dimnames(E_2) <- dimnames(E_3) <- list(iterations=NULL,NULL,NULL,NULL)


### Now save all posterior samples
deriv <- c('F_1','F_2','F_3','E_1','E_2','E_3')
save(list=c(parameters,deriv),file='posterior_samples.Rdata')

### And finally, calculate the posterior summaries and save them
library(posterior)
## Posterior summaries of the model parameters
## Store the summary for each model parameter in an element of the list "posterior_summaries" named after the parameter
for(p in parameters) {
  assign(p,extract(out,p,permuted=F))
  cilims <- c(0.0005,0.005,0.025,0.25,0.5,0.75,0.975,0.995,0.9995)
  summ <- as.data.frame(summarise_draws(get(p),'mean','sd',~quantile(.x, probs = cilims),'mcse_mean','ess_basic','rhat'))
  rownames(summ) <- summ[,1]
  summ <- summ[,-1]
  colnames(summ) <- c('mean','sd',paste0(cilims*100,'%'),'mcse_mean','n_eff','Rhat')
  assign(p,summ)
}

## Posterior summaries of derived parameters, F and E
summarize <- function(x) {
  output <- c(mean(x,na.rm=T),sd(x,na.rm=T),quantile(x,probs=cilims,na.rm=T))
  names(output) <- c('mean','sd',paste0(cilims*100,'%'))
  return(output)
}
for(p in deriv) {
  param <- get(p)
  dims <- dim(param)
  pnames <- vector()
  summ <- matrix(nrow=0,ncol=11)
  for (d1 in 1:dims[2]) {
    for (d2 in 1:dims[3]) { 
      for (d3 in 1:dims[4]) { 
        pnames <- c(pnames,paste0(p,'[',d1,',',d2,',',d3,']'))
        summ <- rbind(summ,summarize(param[,d1,d2,d3]))
      }
    }
  }
  rownames(summ) <- pnames
  summ <- data.frame(summ)
  colnames(summ) <- c('mean','sd',paste0(cilims*100,'%'))
  assign(p,summ)
}

## Save the posterior summaries
save(list=c(parameters,deriv),file='posterior_summaries.Rdata')
