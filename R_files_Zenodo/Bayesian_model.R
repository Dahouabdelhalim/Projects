source ("Library.R")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



#####melanogaster####

d <- read.csv("normmax_imaging_4C_melanogaster.csv")
d$x <- d$IPI
d$y <- d$normmax

observation <- list(N = nrow(d), x = d$x, y = d$y, Ninds = length(unique(d$ID)), id = d$ID)

params = c(a0 = 0.1, b0 = 0.01, c0 = 0.1, d0 = 0.1, SD_a = 0.005, SD_b = 0.005, SD_c = 0.005, SD_d = 0.005)
init0 = c(as.list(params), shape = 50)
init0$a = rep(0.1, observation$Ninds)
init0$b = rep(0.01, observation$Ninds)
init0$c = rep(0.1, observation$Ninds)
init0$d = rep(0.1, observation$Ninds)
init = setNames(rep(list(init0), 4), LETTERS[1:4])

expo = function(x, y0, lambda) y0 * exp(- lambda * x)

genfunc = function(params) {
  function(x) {
    expo(x, params["a0"] , params["b0"]) - expo(x, params["c0"], params["d0"])
  }
}

mod = rstan::stan_model(file = "Bayesian_model.stan")
fit = rstan::sampling(mod, data = observation, init = init, iter = 10000, seed = 1)
Rhat = summary(fit)$summary[,"Rhat"]
coef = broom.mixed::tidyMCMC(fit) %>% {setNames(.$estimate, .$term)} %>% print()

write.csv(coef, "coef_melanogaster.csv")
write.csv(Rhat, "Rhat_melanogaster.csv")


#####simulans####

d <- read.csv("normmax_imaging_4C_simulans.csv")
d$x <- d$IPI
d$y <- d$normmax

observation <- list(N = nrow(d), x = d$x, y = d$y, Ninds = length(unique(d$ID)), id = d$ID)

params = c(a0 = 0.1, b0 = 0.01, c0 = 0.1, d0 = 0.1, SD_a = 0.005, SD_b = 0.005, SD_c = 0.005, SD_d = 0.005)
init0 = c(as.list(params), shape = 50)
init0$a = rep(0.1, observation$Ninds)
init0$b = rep(0.01, observation$Ninds)
init0$c = rep(0.1, observation$Ninds)
init0$d = rep(0.1, observation$Ninds)
init = setNames(rep(list(init0), 4), LETTERS[1:4])

expo = function(x, y0, lambda) y0 * exp(- lambda * x)

genfunc = function(params) {
  function(x) {
    expo(x, params["a0"] , params["b0"]) - expo(x, params["c0"], params["d0"])
  }
}

mod = rstan::stan_model(file = "Bayesian_model.stan")
fit = rstan::sampling(mod, data = observation, init = init, iter = 10000, seed = 1)
Rhat = summary(fit)$summary[,"Rhat"]
coef = broom.mixed::tidyMCMC(fit) %>% {setNames(.$estimate, .$term)} %>% print()

write.csv(coef, "coef_simulans.csv")
write.csv(Rhat, "Rhat_simulans.csv")


