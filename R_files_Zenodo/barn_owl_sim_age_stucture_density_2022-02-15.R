# What: Simulations of age structured population dynamics ----

# Will perform simulations to test how reproductive value weighting performes under density dependence and density-dependent selection. The characteristics of reproductive values should be conserved when density affects all age classes similarly.

# Density-independent dynamics
# n_t+1 = L*n_t, where L = l + epsilon_t, epsilon_t = demograpic and environmental variance

# Density-dependent dynamics
# n_t+1 = M(n_t)n_t (Liu and Cohen 1987), where M(n_t) = D(n_t)(L_max+epsilon_t), D(n_t) contain all density dependent effects on fecundities and survivals and epsilon_t = demograpic and environmental variance.

#



# Load: Libraries ----

library(ggplot2)
library(cowplot);theme_set(theme_cowplot())
library(reshape2)
library(grid)
library(gridExtra)
library(dplyr)

#



# Create: Define projection matrices ----

# Use the example from Engen et al 2014

# Leslie matrix
l_t1_14 <- matrix(c(0, 0, 2, 3, 2,
                    0.3, 0, 0, 0, 0,
                    0, 0.6, 0, 0, 0,
                    0, 0, 0.9, 0, 0,
                    0, 0, 0, 0.9, 0), ncol = 5, byrow=T)

# Get konstant for lambda == 1
konstant_lt1_14 <- eulerlotka(lambda = 1, leslie = l_t1_14, census = "pre", konstant = 1)$konstant

# Lambda == 1
l_t1_14 <- rbind(l_t1_14[1, ]*konstant_lt1_14, l_t1_14[-1,])
eigenl2(l_t1_14)

# Lambda_max
# Lambda == max (scaling konstant to get lambda ~ 1.50)
l_tmax_14 <- rbind(l_t1_14[1, ]/exp(-0.00010*12000),
                l_t1_14[2:5,]/exp(-0.00001*12000)) # -0.0005
eigenl2(l_tmax_14)

# Reproductive values (at K)
v_14 <- eigenl2(l_t1_14)$v

# Stable age distribution (at K)
u_14 <- eigenl2(l_t1_14)$u

#



# 1. Run deterministic: Simulate N and phenotype (z) without density dependence (Engen et al 2014) ----

# Population sizes
n <- c(2000, 2000, 2000, 2000, 2000)
time <- 40
nt <- data.frame(time = 1:time, a1 = rep(NA, time), a2 = rep(NA, time), a3 = rep(NA, time), a4 = rep(NA, time), a5 = rep(NA, time), n = NA, nv = NA)

# Phenotypes
z <- c(1200, 1200, 1300, 1400, 1500)
zt <- data.frame(time = 1:time, a1 = rep(NA, time), a2 = rep(NA, time), a3 = rep(NA, time), a4 = rep(NA, time), a5 = rep(NA, time), z = NA, zv = NA)

# Starting point
# N
nt[1, 2:6] <- n
nt$n[1] <- sum(nt[1, 2:6])
nt$nv[1] <- sum(nt[1, 2:6]*v_14)
# z
zt[1, 2:6] <- z
zt$z[1] <- weighted.mean(zt[1, 2:6], w = nt[1, 2:6])
zt$zv[1] <- weighted.mean(zt[1, 2:6], w = nt[1, 2:6]*v_14)

# Iterate
for (i in 2:time){
  d_var <- rnorm(25, 0, 0) # sd = 0.1
  e_var <- rnorm(1, 0, 0) # sd = 0.05
  e <- matrix(c(1,1,1,1,1,
                0.8,0,0,0,0,
                0,0.8,0,0,0,
                0,0,0.8,0,0,
                0,0,0,0.8,0.8),
              byrow = TRUE,
              ncol = 5) * (e_var+(d_var/sum(nt[i-1, 2:6])))
  lte <- l_t1_14+e
  nt[i, 2:6] <- lte%*%t(nt[i-1, 2:6])
  # Rowsums N
  nt$n[i] <- sum(nt[i, 2:6])
  nt$nv[i] <- sum(nt[i, 2:6]*v_14)
  # Phenotypes z
  zt[i, 2] <- weighted.mean(zt[i-1, 2:6], w = nt[i-1, 2:6]*lte[1, ])
  zt[i, 3:5] <- zt[i-1, 2:4]
  zt[i, 6] <- weighted.mean(zt[i-1, 5:6], w = nt[i-1, 5:6]*lte[5, 4:5])
  # Rowsums z
  zt$z[i] <- weighted.mean(zt[i, 2:6], w = nt[i, 2:6])
  zt$zv[i] <- weighted.mean(zt[i, 2:6], w = nt[i, 2:6]*v_14)
}

# Plot population size
ggplot(data = subset(nt, time > 0), mapping = aes(x = time, y = n)) +
  geom_line(mapping = aes(y = nv), linetype = 2) +
  geom_line()

# Plot phenotype
ggplot(data = subset(zt, time > 0), mapping = aes(x = time, y = z)) +
  geom_line(mapping = aes(y = zv), linetype = 2) +
  geom_line()

# Keep simulation
nt_1 <- nt
zt_1 <- zt


#



# 2. Run deterministic: Simulate N and phenotype (z) with density dependence (Engen et al 2014) ----

# Population sizes
n <- c(2000, 2000, 2000, 2000, 2000)
time <- 40
nt <- data.frame(time = 1:time, a1 = rep(NA, time), a2 = rep(NA, time), a3 = rep(NA, time), a4 = rep(NA, time), a5 = rep(NA, time), n = NA, nv = NA)

# Phenotypes
z <- c(1200, 1200, 1300, 1400, 1500)
zt <- data.frame(time = 1:time, a1 = rep(NA, time), a2 = rep(NA, time), a3 = rep(NA, time), a4 = rep(NA, time), a5 = rep(NA, time), z = NA, zv = NA)

# List for v's
v_list <- data.frame(time = 1:time, n = NA, v1 = NA, v2 = NA, v3 = NA, v4 = NA, v5 = NA)
v_list$n[1] <- sum(n)
v_list[1, 3:7] <- v_14

# Starting point
# N
nt[1, 2:6] <- n
nt$n[1] <- sum(nt[1, 2:6])
nt$nv[1] <- sum(nt[1, 2:6]*v_14)
# z
zt[1, 2:6] <- z
zt$z[1] <- weighted.mean(zt[1, 2:6], w = nt[1, 2:6])
zt$zv[1] <- weighted.mean(zt[1, 2:6], w = nt[1, 2:6]*v_14)

# Iterate
for (i in 2:time){
  d_var <- rnorm(25, 0, 0.0) # sd = 0.1
  e_var <- rnorm(1, 0, 0.0) # sd = 0.05
  e <- matrix(c(1,1,1,1,1,
                1,0,0,0,0,
                0,1,0,0,0,
                0,0,1,0,0,
                0,0,0,1,1), byrow = TRUE, ncol = 5) * (e_var+d_var/sum(nt[i-1, 2:6]))
  #lte <- l_tmax_14+e
  # Density functions for fecundity (Fn) and survival (En)
  Fn <- exp(-0.00010*sum(nt[i-1, 2:6]))
  En <- exp(-0.00001*sum(nt[i-1, 2:6]))
  # D(N)
  dn <- matrix(c(Fn,Fn,Fn,Fn,Fn,
                 En,0,0,0,0,
                 0,En,0,0,0,
                 0,0,En,0,0,
                 0,0,0,En,En), byrow = TRUE, ncol = 5)
  # N_t+1 = N_t+D(N)(M-1)N_t (Jensen 1995)
  lte <- dn*l_tmax_14+e
  nt[i, 2:6] <- (lte) %*% t(nt[i-1, 2:6])
  # Estimate v
  v_list[i, 2] <- sum(nt[i, 2:6])
  v_list[i, 3:7] <- eigenl2(lte)$v
  # Rowsums N
  nt$n[i] <- sum(nt[i, 2:6])
  nt$nv[i] <- sum(nt[i, 2:6]*v_14)
  # Phenotypes z
  zt[i, 2] <- weighted.mean(zt[i-1, 2:6], w = nt[i-1, 2:6]*lte[1, ])
  zt[i, 3:5] <- zt[i-1, 2:4]
  zt[i, 6] <- weighted.mean(zt[i-1, 5:6], w = nt[i-1, 5:6]*lte[5, 4:5])
  # Rowsums z
  zt$z[i] <- weighted.mean(zt[i, 2:6], w = nt[i, 2:6])
  zt$zv[i] <- weighted.mean(zt[i, 2:6], w = nt[i, 2:6]*v_14)
}

# To long format v_list
v_list <- melt(v_list, id.vars= c("time", "n"), variable.name = "age", value.name = "v")
head(v_list)
# Cumulative mean
v_list <- v_list %>%
  group_by(age) %>%
  mutate(vc = cummean(v))

# Temporal mean
v_list_mean <- v_list %>%
  group_by(age) %>%
  summarise(vm = mean(v),
            vse= se(v)) %>%
  left_join(., data.frame(age = factor(c("v1", "v2", "v3", "v4", "v5")), v = v_14))

# Plot reproductive value means
ggplot(data = v_list_mean, mapping = aes(x = age, y = vm-v, ymin = vm-v-vse*1.96, ymax=vm-v+vse*1.96)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = 3) +
  theme_cowplot()

# Plot reproductive values over time
ggplot(data = v_list, mapping = aes(x = time, y = v, linetype = age)) +
  geom_hline(data = data.frame(time = 1, age = c("v1", "v2", "v3", "v4", "v5"), v = v_14),
             mapping = aes(yintercept = v, linetype = age), colour = "grey70") + 
  geom_line()

# Plot population size
ggplot(data = subset(nt, time > 0), mapping = aes(x = time, y = n)) +
  geom_line(mapping = aes(y = nv), linetype = 2) +
  geom_line()

# Plot phenotype
ggplot(data = subset(zt, time > 0), mapping = aes(x = time, y = z)) +
  geom_line(mapping = aes(y = zv), linetype = 2) +
  geom_line()

# Keep simulation
nt_2 <- nt
zt_2 <- zt

#



# 3. Run deterministic: Simulate N and phenotype (z) with a single episode of fecundity selection without density dependence (Engen et al 2014) ----

# Define L_1j = mean fecundity in age class j, hence l_1j = L_1j + beta*(z-z_mean) is an expression for selection which does not change the mean fitness for weak selection (l_1j == L_1j) (test: round(mean(0.3 + 0.5*(rnorm(1000000, mean = 10, sd = 2)-10)), 2)). To simulate selection one only need to assume a phenotypic variance and define a selection differential. Response = selection differential * heriability, R = S*h^2. Hence, one can calculate the offspring generation using z_offspring = z_parents + response.

# The same arguments hold true for viablity selection except that this will change the mean phenotype of the age class under selection and thus change the mean fitness for this age class. Hence, for viability selection in age class L_21, we have l_21 = L_21 + beta*(z-z_mean) != L_21. The mean fitness after selection can be calculated by l_21 = L_21 + beta*(z_new_mean-z_mean).

# Population sizes
n <- c(2000, 2000, 2000, 2000, 2000)
time <- 40
nt <- data.frame(time = 1:time, a1 = rep(NA, time), a2 = rep(NA, time), a3 = rep(NA, time), a4 = rep(NA, time), a5 = rep(NA, time), n = NA, nv = NA)

# Phenotypes
z <- c(1200, 1200, 1300, 1400, 1500)
zt <- data.frame(time = 1:time, a1 = rep(NA, time), a2 = rep(NA, time), a3 = rep(NA, time), a4 = rep(NA, time), a5 = rep(NA, time), z = NA, zv = NA)

# Starting point
# N
nt[1, 2:6] <- n
nt$n[1] <- sum(nt[1, 2:6])
nt$nv[1] <- sum(nt[1, 2:6]*v_14)
# z
zt[1, 2:6] <- z
zt$z[1] <- weighted.mean(zt[1, 2:6], w = nt[1, 2:6])
zt$zv[1] <- weighted.mean(zt[1, 2:6], w = nt[1, 2:6]*v_14)

# Iterate
for (i in 2:time){
  d_var <- rnorm(25, 0, 0) # sd = 0.1
  e_var <- rnorm(1, 0, 0) # sd = 0.05
  e <- matrix(c(1,1,1,1,1,
                0.8,0,0,0,0,
                0,0.8,0,0,0,
                0,0,0.8,0,0,
                0,0,0,0.8,0.8),
              byrow = TRUE,
              ncol = 5) * (e_var+(d_var/sum(nt[i-1, 2:6])))
  lte <- l_t1_14+e
  nt[i, 2:6] <- lte%*%t(nt[i-1, 2:6])
  # Rowsums N
  nt$n[i] <- sum(nt[i, 2:6])
  nt$nv[i] <- sum(nt[i, 2:6]*v_14)
  # Phenotypes z
  # Define a single case of fecundity selection in the first age class (element L_11): beta=selection gradient=S/Var(P)=0.1, Var(P)=phenotypic variance=1000, S=selection differential=beta*Var(P)=0.1*1000=100, h^2=heritability=0.5, R=response=S*h^2=100*0.5=50
  zt[i, 2] <- weighted.mean(zt[i-1, 2:6]+ifelse(i==20, c(50,0,0,0,0), 0), w = nt[i-1, 2:6]*lte[1, ])
  zt[i, 3:5] <- zt[i-1, 2:4]
  zt[i, 6] <- weighted.mean(zt[i-1, 5:6], w = nt[i-1, 5:6]*lte[5, 4:5])
  # Rowsums z
  zt$z[i] <- weighted.mean(zt[i, 2:6], w = nt[i, 2:6])
  zt$zv[i] <- weighted.mean(zt[i, 2:6], w = nt[i, 2:6]*v_14)
}

# Plot population size
ggplot(data = subset(nt, time > 0), mapping = aes(x = time, y = n)) +
  geom_line(mapping = aes(y = nv), linetype = 2) +
  geom_line()

# Plot phenotype
ggplot(data = subset(zt, time > 0), mapping = aes(x = time, y = z)) +
  geom_line(mapping = aes(y = zv), linetype = 2) +
  geom_line()

# Keep simulation
nt_3 <- nt
zt_3 <- zt

#



# 4. Run deterministic: Simulate N and phenotype (z) with a single episode of fecundity selection with density dependence (Engen et al 2014) ----

# Define L_1j = mean fecundity in age class j, hence l_1j = L_1j + beta*(z-z_mean) is an expression for selection which does not change the mean fitness for weak selection (l_1j == L_1j) (test: round(mean(0.3 + 0.5*(rnorm(1000000, mean = 10, sd = 2)-10)), 2)). To simulate selection one only need to assume a phenotypic variance and define a selection differential. Response = selection differential * heriability, R = S*h^2. Hence, one can calculate the offspring generation using z_offspring = z_parents + response.

# The same arguments hold true for viablity selection except that this will change the mean phenotype of the age class under selection and thus change the mean fitness for this age class. Hence, for viability selection in age class L_21, we have l_21 = L_21 + beta*(z-z_mean) != L_21. The mean fitness after selection can be calculated by l_21 = L_21 + beta*(z_new_mean-z_mean).

# Population sizes
n <- c(2000, 2000, 2000, 2000, 2000)
time <- 40
nt <- data.frame(time = 1:time, a1 = rep(NA, time), a2 = rep(NA, time), a3 = rep(NA, time), a4 = rep(NA, time), a5 = rep(NA, time), n = NA, nv = NA)

# Phenotypes
z <- c(1200, 1200, 1300, 1400, 1500)
zt <- data.frame(time = 1:time, a1 = rep(NA, time), a2 = rep(NA, time), a3 = rep(NA, time), a4 = rep(NA, time), a5 = rep(NA, time), z = NA, zv = NA)

# Starting point
# N
nt[1, 2:6] <- n
nt$n[1] <- sum(nt[1, 2:6])
nt$nv[1] <- sum(nt[1, 2:6]*v_14)
# z
zt[1, 2:6] <- z
zt$z[1] <- weighted.mean(zt[1, 2:6], w = nt[1, 2:6])
zt$zv[1] <- weighted.mean(zt[1, 2:6], w = nt[1, 2:6]*v_14)

# Iterate
for (i in 2:time){
  d_var <- rnorm(25, 0, 0) # sd = 0.1
  e_var <- rnorm(1, 0, 0) # sd = 0.05
  e <- matrix(c(1,1,1,1,1,
                1,0,0,0,0,
                0,1,0,0,0,
                0,0,1,0,0,
                0,0,0,1,1), byrow = TRUE, ncol = 5) * (e_var+d_var/sum(nt[i-1, 2:6]))
  #lte <- l_tmax_14+e
  # Density functions for fecundity (Fn) and survival (En)
  Fn <- exp(-0.00010*sum(nt[i-1, 2:6]))
  En <- exp(-0.00001*sum(nt[i-1, 2:6]))
  # D(N)
  dn <- matrix(c(Fn,Fn,Fn,Fn,Fn,
                 En,0,0,0,0,
                 0,En,0,0,0,
                 0,0,En,0,0,
                 0,0,0,En,En), byrow = TRUE, ncol = 5)
  # N_t+1 = N_t+D(N)(M-1)N_t (Jensen 1995)
  lte <- dn*l_tmax_14+e
  nt[i, 2:6] <- (lte) %*% t(nt[i-1, 2:6])
  # Rowsums N
  nt$n[i] <- sum(nt[i, 2:6])
  nt$nv[i] <- sum(nt[i, 2:6]*v_14)
  # Phenotypes z
  # Define a single case of fecundity selection in the first age class (element L_11): beta=selection gradient=S/Var(P)=0.1, Var(P)=phenotypic variance=1000, S=selection differential=beta*Var(P)=0.1*1000=100, h^2=heritability=0.5, R=response=S*h^2=100*0.5=50
  zt[i, 2] <- weighted.mean(zt[i-1, 2:6]+ifelse(i==20, c(50,0,0,0,0), 0), w = nt[i-1, 2:6]*lte[1, ])
  zt[i, 3:5] <- zt[i-1, 2:4]
  zt[i, 6] <- weighted.mean(zt[i-1, 5:6], w = nt[i-1, 5:6]*lte[5, 4:5])
  # Rowsums z
  zt$z[i] <- weighted.mean(zt[i, 2:6], w = nt[i, 2:6])
  zt$zv[i] <- weighted.mean(zt[i, 2:6], w = nt[i, 2:6]*v_14)
}

# Plot population size
ggplot(data = subset(nt, time > 0), mapping = aes(x = time, y = n)) +
  geom_line(mapping = aes(y = nv), linetype = 2) +
  geom_line()

# Plot phenotype
ggplot(data = subset(zt, time > 0), mapping = aes(x = time, y = z)) +
  geom_line(mapping = aes(y = zv), linetype = 2) +
  geom_line()

# Keep simulation
nt_4 <- nt
zt_4 <- zt

#



# 5. Run deterministic: Simulate N and phenotype (z) with a single episode of density dependent fecundity selection with density dependence (Engen et al 2014) ----

# Define L_1j = mean fecundity in age class j, hence l_1j = L_1j * exp(beta_1*(z-z_mean) + beta_2*(n-n_mean) + beta_3*(z-z_mean)*(n-n_mean)) is an expression for selection which does not change the mean fitness for weak selection only change due to density changes in the second parameter. Hence, to simulate selection one only need to assume a phenotypic variance (Var(P)) and define a selection differential. Response = selection differential * heriability, R = S*h^2. Hence, one can calculate the offspring generation using z_offspring = z_parents + response. With density dependent selection selection differential=S=(beta_1+beta3*(n-n_mean))*Var(P)

# The same arguments hold true for viablity selection except that this will change the mean phenotype of the age class under selection and thus change the mean fitness for this age class. Hence, for viability selection in age class L_21, we have l_21 = L_21 + beta*(z-z_mean) != L_21. The mean fitness after selection can be calculated by l_21 = L_21 + beta*(z_new_mean-z_mean).

# Population sizes
n <- c(2000, 2000, 2000, 2000, 2000)
time <- 40
nt <- data.frame(time = 1:time, a1 = rep(NA, time), a2 = rep(NA, time), a3 = rep(NA, time), a4 = rep(NA, time), a5 = rep(NA, time), n = NA, nv = NA)

# Phenotypes
z <- c(1200, 1200, 1300, 1400, 1500)
zt <- data.frame(time = 1:time, a1 = rep(NA, time), a2 = rep(NA, time), a3 = rep(NA, time), a4 = rep(NA, time), a5 = rep(NA, time), z = NA, zv = NA)

# Starting point
# N
nt[1, 2:6] <- n
nt$n[1] <- sum(nt[1, 2:6])
nt$nv[1] <- sum(nt[1, 2:6]*v_14)
# z
zt[1, 2:6] <- z
zt$z[1] <- weighted.mean(zt[1, 2:6], w = nt[1, 2:6])
zt$zv[1] <- weighted.mean(zt[1, 2:6], w = nt[1, 2:6]*v_14)

# Iterate
for (i in 2:time){
  d_var <- rnorm(25, 0, 0) # sd = 0.1
  e_var <- rnorm(1, 0, 0) # sd = 0.05
  e <- matrix(c(1,1,1,1,1,
                1,0,0,0,0,
                0,1,0,0,0,
                0,0,1,0,0,
                0,0,0,1,1), byrow = TRUE, ncol = 5) * (e_var+d_var/sum(nt[i-1, 2:6]))
  #lte <- l_tmax_14+e
  # Density functions for fecundity (Fn) and survival (En)
  Fn <- exp(-0.00010*sum(nt[i-1, 2:6]))
  En <- exp(-0.00001*sum(nt[i-1, 2:6]))
  # D(N)
  dn <- matrix(c(Fn,Fn,Fn,Fn,Fn,
                 En,0,0,0,0,
                 0,En,0,0,0,
                 0,0,En,0,0,
                 0,0,0,En,En), byrow = TRUE, ncol = 5)
  # N_t+1 = N_t+D(N)(M-1)N_t (Jensen 1995)
  lte <- dn*l_tmax_14+e
  nt[i, 2:6] <- (lte) %*% t(nt[i-1, 2:6])
  # Rowsums N
  nt$n[i] <- sum(nt[i, 2:6])
  nt$nv[i] <- sum(nt[i, 2:6]*v_14)
  # Phenotypes z
  # Define a single case of density dependent fecundity selection in the first age class (element L_11): beta_1 + beta_3*(n-n_mean)=-0.5+0.0007*(n-12000), selection gradient=S/Var(P), Var(P)=phenotypic variance=1000, S=selection differential=beta*Var(P)=-0.5+0.0007*(n-12000)*1000, h^2=heritability=0.5, R=response=S*h^2=-0.5+0.0007*(n-12000)*1000*0.5
  zt[i, 2] <- weighted.mean(zt[i-1, 2:6]++ifelse(i==10, c(-0.5+0.001*(nt$n[i]-12000)*1000*0.5,0,0,0,0), 0), w = nt[i-1, 2:6]*lte[1, ])
  zt[i, 3:5] <- zt[i-1, 2:4]
  zt[i, 6] <- weighted.mean(zt[i-1, 5:6], w = nt[i-1, 5:6]*lte[5, 4:5])
  # Rowsums z
  zt$z[i] <- weighted.mean(zt[i, 2:6], w = nt[i, 2:6])
  zt$zv[i] <- weighted.mean(zt[i, 2:6], w = nt[i, 2:6]*v_14)
}

# Plot population size
ggplot(data = subset(nt, time > 0), mapping = aes(x = time, y = n)) +
  geom_line(mapping = aes(y = nv), linetype = 2) +
  geom_line()

# Plot phenotype
ggplot(data = subset(zt, time > 0), mapping = aes(x = time, y = z)) +
  geom_line(mapping = aes(y = zv), linetype = 2) +
  geom_line()

# Keep simulation
nt_5 <- nt
zt_5 <- zt

#



# Explore reproductive values for different N's (Engen et al 2014) ----

# At K
eigenl2(l_t1_14)

# Ns
n <- 10000:14000

# List for v's
v_list <- data.frame(N = n, v1 = rep(NA, length(n)), v2 = NA, v3 = NA, v4 = NA, v5 = NA)

# Loop
for (i in 1:length(n)){
 
  # Select N
  N <- (n)[i]
  
  # Max L
  l_tmax_14 <- rbind(l_t1_14[1, ]/exp(-0.00010*12000),
                  l_t1_14[2:5,]/exp(-0.00001*12000))
  
  # Density-dependence
  Fn <- exp(-0.00010*N)
  En <- exp(-0.00001*N)
  # D(N)
  dn <- matrix(c(Fn,Fn,Fn,Fn,Fn,
                 En,0,0,0,0,
                 0,En,0,0,0,
                 0,0,En,0,0,
                 0,0,0,En,En), byrow = TRUE, ncol = 5)
  # Estimate v
  v_list[i, 2:6] <- eigenl2(dn*l_tmax_14)$v
   
}

# To long format
v_list <- melt(v_list, id.vars= "N", variable.name = "v")
head(v_list)

# Plot
ggplot(data = v_list, mapping = aes(x=N, y=value, linetype = v)) +
  geom_line() +
  geom_vline(xintercept = 12000, colour = "red") +
  ylab("Reproductive value")

#



# Figure: Collect the figures into one figure with five panels ----

# Plot the five panels
p1 <- ggplot(data = zt_1, mapping = aes(x = time, y = z)) +
  geom_line(mapping = aes(y = zv), linetype = 2, size = 0.25) +
  geom_line(size = 0.25) +
  geom_text(mapping = aes(label = lab), data = data.frame(time = 1, z = 1385, lab = "a"), fontface = "bold", size = 10/(14/5)) +
  ggtitle("Density-independent") +
  scale_y_continuous(breaks = c(1375, 1350, 1325, 1300), limits = c(1280, 1390)) +
  xlab("") +
  ylab("") +
  theme_cowplot(font_size = 10, line_size = 0.25) +
  theme(plot.margin = unit(c(0.5,0,0,0), "cm"))
p2 <- ggplot(data = zt_2, mapping = aes(x = time, y = z)) +
  geom_line(mapping = aes(y = zv), linetype = 2, size = 0.25) +
  geom_line(size = 0.25) +
  geom_text(mapping = aes(label = lab), data = data.frame(time = 1, z = 1385, lab = "b"), fontface = "bold", size = 10/(14/5)) +
  ggtitle("Density-dependent") +
  scale_y_continuous(breaks = c(1375, 1350, 1325, 1300), limits = c(1280, 1390)) +
  xlab("") +
  ylab("") +
  theme_cowplot(font_size = 10, line_size = 0.25) +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "cm"))
p3 <- ggplot(data = zt_3, mapping = aes(x = time, y = z)) +
  geom_line(mapping = aes(y = zv), linetype = 2, size = 0.25) +
  geom_line(size = 0.25) +
  geom_text(mapping = aes(label = lab), data = data.frame(time = 1, z = 1385, lab = "c"), fontface = "bold", size = 10/(14/5)) +
  scale_y_continuous(breaks = c(1375, 1350, 1325, 1300), limits = c(1280, 1390)) +
  xlab("") +
  ylab("") +
  theme_cowplot(font_size = 10, line_size = 0.25) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
p4 <- ggplot(data = zt_4, mapping = aes(x = time, y = z)) +
  geom_line(mapping = aes(y = zv), linetype = 2, size = 0.25) +
  geom_line(size = 0.25) +
  geom_text(mapping = aes(label = lab), data = data.frame(time = 1, z = 1385, lab = "d"), fontface = "bold", size = 10/(14/5)) +
  scale_y_continuous(breaks = c(1375, 1350, 1325, 1300), limits = c(1280, 1390)) +
  xlab("") +
  ylab("") +
  theme_cowplot(font_size = 10, line_size = 0.25) +
  theme(plot.margin = unit(c(0,0.5,0,0), "cm"))
p5 <- ggplot(data = zt_5, mapping = aes(x = time, y = z)) +
  geom_blank() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))
p6 <- ggplot(data = zt_5, mapping = aes(x = time, y = z)) +
  geom_line(mapping = aes(y = zv), linetype = 2, size = 0.25) +
  geom_line(size = 0.25) +
  geom_text(mapping = aes(label = lab), data = data.frame(time = 1, z = 1385, lab = "e"), fontface = "bold", size = 10/(14/5)) +
  scale_y_continuous(breaks = c(1375, 1350, 1325, 1300), limits = c(1280, 1390)) +
  xlab("") +
  ylab("") +
  theme_cowplot(font_size = 10, line_size = 0.25) +
  theme(plot.margin = unit(c(0,0.5,0,0), "cm"))

# Collect into one figure
p_1to5 <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)

# Add common x and y labels
p_1to5 <- grid.arrange(arrangeGrob(p_1to5,
                         left = textGrob("Mean phenotype", gp=gpar(fontsize=10), rot=90),
                         bottom = textGrob("Time",gp=gpar(fontsize=10))))

# Save to file
ggsave(filename = "output/pheno_sim.pdf", plot = p_1to5, width = 150, height = 150, units = "mm")

#



# Clean ----
rm(nt, zt, l_t1, konstant_lt1, l_tmax, v, l_t1_14, konstant_lt1_14, l_tmax_14, v_14, v_list, n, N, dn, d_var, e_var, e)