# experiment varying parasite density
d <- read.csv("data/Density.csv")
d$host <- 1  # number of tadpoles in each container
d$time <- 30  # time elapsed
d$Volume <- 1
d$Number <- d$Cerc_Density

# fit models ----------------------------
n_models <- length(models())
model_fits <- vector(mode = "list", length = n_models)
names(model_fits) <- models()
for (i in 1:n_models) {
    model_fits[[i]] <- fit_model(d, model = models()[i])
}

# make aic table -------------------------------
header <- "Table 2d: varying parasite density"
aic_tab <- aic_table(model_fits, d)
save_aic_table(aic_tab, header, filename = "tables/paras_density_table.md")
(best_models <- rownames(subset(aic_tab, dAICc < 2)))

# make parameter table ----------------------------------------------------
make_par_table(model_fits, name = "tables/parameter_estimates/paras_density.csv")
