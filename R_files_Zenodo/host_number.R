# experiment varying density of hosts
d <- read.csv("data/Host_Density_Combined.csv")
d$Meta <- d$Meta_Number
d$host <- d$Host_Number  # number of tadpoles in each container
d$time <- 30  # time elapsed
d$Volume <- 2.1
d$Number <- d$Cerc_Number

# fit models ----------------------------
n_models <- length(models())
model_fits <- vector(mode = "list", length = n_models)
names(model_fits) <- models()
for (i in 1:n_models) {
  model_fits[[i]] <- fit_model(d, model = models()[i])
}

# make aic table -------------------------------
header <- "Table 2b: varying host density"
aic_tab <- aic_table(model_fits, d)
save_aic_table(aic_tab, header, filename = "tables/host_den_table.md")

# print best model fits -----------------------------
(best_models <- rownames(subset(aic_tab, dAICc < 2)))

# make parameter table ----------------------------------------------------
make_par_table(model_fits, name = "tables/parameter_estimates/host_den.csv")
