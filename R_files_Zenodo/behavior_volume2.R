# experiment varying # parasites and volume to maintain
# constant  vcx
d <- read.csv("data/Behavior Volume2.csv")
d$host <- 1  # number of tadpoles in each container
d$time <- 30  # time elapsed

# fit models ----------------------------
n_models <- length(models())
model_fits <- vector(mode = "list", length = n_models)
names(model_fits) <- models()
for (i in 1:n_models) {
    model_fits[[i]] <- fit_model(d, model = models()[i])
}

convergence_codes <- unlist(lapply(model_fits, function(x) x$convergence))
if (any(convergence_codes != 0)) {
  stop("Model did not converge")
}

# make aic table -------------------------------
header <- "Table 2f: varying parasite number, with tadpoles anaesthesized"
aic_tab <- aic_table(model_fits, d)
save_aic_table(aic_tab, header, filename = "tables/behavior_table.md")

# visualize best model fits ----------------------
(best_models <- rownames(subset(aic_tab, dAICc < 2)))


# make parameter table ----------------------------------------------------
make_par_table(model_fits, name = "tables/parameter_estimates/behavior.csv")
