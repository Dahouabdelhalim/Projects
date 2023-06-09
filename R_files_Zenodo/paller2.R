d <- read.csv("data/Paller2.csv")
d$host <- 1  # number of tadpoles in each container
d$Meta <- d$Metacerc
d$Number <- d$Cerc

# fit models ----------------------------
n_models <- length(models())
model_fits <- vector(mode = "list", length = n_models)
names(model_fits) <- models()
for (i in 1:n_models) {
  model_fits[[i]] <- fit_model(d, model = models()[i], use_ode = TRUE)
}


# Frequency dependent 2 is problematic, handle manually -------------------
# simple line search
pvec <- seq(0.11, .12, length.out = 100)
nll <- rep(NA, 20)
for (i in seq_along(pvec)) {
  nll[i] <- NLL.freq2_ode(pvec[i], Np = d$Number, H = 1, 
                          v = d$Volume, t = d$time, metacerc = d$Metacerc)
}

par <- c(beta = .118)
optim_res <- optim(par, fn = NLL.freq2_ode, Np = d$Number, H = d$host, 
                   v = d$Volume, t = d$time, metacerc = d$Meta, method = "Brent", 
                   lower = 0, upper = .12)
optim_res$model <- fd2_ode
names(optim_res$par) <- "beta"
class(optim_res) <- c("optim_res", "list")

model_fits[["Frequency dependent 2"]] <- optim_res

# make aic table -------------------------------
header <- "Paller: varying parasite number"
aic_tab <- aic_table(model_fits, d)
save_aic_table(aic_tab, header, filename = "tables/paller2_table.md")
(best_models <- rownames(subset(aic_tab, dAICc < 2)))

# make parameter table ----------------------------------------------------
make_par_table(model_fits, name = "tables/parameter_estimates/paller2.csv")

