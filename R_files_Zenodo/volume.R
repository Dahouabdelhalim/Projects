d <- read.csv("data/Volume.csv")
d$time <- 30
d$host <- 1  # number of tadpoles in each container
d$Number <- d$Cerc_Number

# fit models ----------------------------
n_models <- length(models())
model_fits <- vector(mode = "list", length = n_models)
names(model_fits) <- models()
for (i in 1:n_models) {
  model_fits[[i]] <- fit_model(d, model = models()[i])
}

# make aic table -------------------------------
header <- "Table 2a: varying parasite number (constant density)"
aic_tab <- aic_table(model_fits, d)
save_aic_table(aic_tab, header, filename = "tables/volume_table.md")
(best_models <- rownames(subset(aic_tab, dAICc < 2)))

# make parameter table ----------------------------------------------------
make_par_table(model_fits, name = "tables/parameter_estimates/volume.csv")
