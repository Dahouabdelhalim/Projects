d <- read.csv("data/Karvonen2.csv")
d$host <- 1  # number of tadpoles in each container
d$Number <- d$Cerc
d$Meta <- d$Metacerc

# fit models ----------------------------
n_models <- length(models())
model_fits <- vector(mode = "list", length = n_models)
names(model_fits) <- models()
for (i in 1:n_models) {
  model_fits[[i]] <- fit_model(d, model = models()[i], use_ode = TRUE)
}

# make aic table -------------------------------
header <- "Karvonen: varying parasite number"
aic_tab <- aic_table(model_fits, d)
save_aic_table(aic_tab, header, filename = "tables/karvonen2_table.md")

# visualize best model fits ----------------------
(best_models <- rownames(subset(aic_tab, dAICc < 2)))

# make parameter table ----------------------------------------------------
make_par_table(model_fits, name = "tables/parameter_estimates/karvonen2.csv")
