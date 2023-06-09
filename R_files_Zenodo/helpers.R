library(knitr)
library(dplyr)

# computes AICc from optim output
aic_c <- function(optim_result, n) {
    k <- length(optim_result$par)
    nll <- optim_result$value
    AIC <- 2 * k + 2 * nll
    AIC + 2 * k * (k + 1)/(n - k - 1)
}

models <- function() {
    # lists all models in the model set
    c("Constant 1", "Constant 2", "Density dependent", 
        "Density independent 1", "Density independent 2", "Power C", 
        "Power H", "Power CH", "Negative binomial 1", "Negative binomial 2")
}

aic_table <- function(optim_l, d) {
    aic_list <- lapply(optim_l, aic_c, n = nrow(d))
    aic_vals <- sort(unlist(aic_list))
    daic <- aic_vals - min(aic_vals)
    aic_tab <- data.frame(AICc = aic_vals, dAICc = round(daic, 3))
    stopifnot(nrow(aic_tab) == length(optim_l))
    aic_tab
}

save_aic_table <- function(aic_tab, header, filename) {
    sink(filename)
    cat(paste("###", header))
    print(kable(aic_tab, format = "markdown"))
    cat("\\n")
    sink()
}

plot_best <- function(fits, best_models, d) {
  plot(jitter(d$Number / d$Volume), d$Meta, 
       xlab = "Density of Cercariae", 
       ylab = "Number of metacercariae")
  xvals <- seq(min(d$Number / d$Volume), max(d$Number / d$Volume), 
               length.out = 100)
  H <- unique(d$host)
  v <- unique(d$Volume)
  t <- unique(d$time)
  if (length(H) > 1) {
    warning("> 1 unique number of hosts in the data. Predictions are valid only 
              if none of the models use H as a parameter.") 
  }
  if (length(v) > 1) {
    warning("> 1 unique number of volumes in the data. Predictions are valid only 
              if none of the models use Volume as a parameter.") 
  }
  if (length(t) > 1) {
    warning("> 1 unique number of times in the data. Predictions are valid only 
              if none of the models use t as a parameter.") 
  }
  for (i in seq_along(best_models)) {
    model <- fits[[best_models[i]]]$model
    estimates <- fits[[best_models[i]]]$par
    mu_y <- model(p = estimates, Np = xvals, H = H[1], v = v[1], t = t[1])
    lines(xvals, mu_y, col = i, lty = i, lwd = 2)
  }
  legend('topleft', 
         col = seq_along(best_models), 
         lty = seq_along(best_models), 
         legend = best_models, bty = 'n', lwd = 2)
}

firstup <- function(x) {
  # capitalize first letter of a string
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

make_par_table <- function(model_fits, name) {
  lapply(model_fits, function(x) x[["par"]]) %>%
    lapply(function(x) data.frame(MLE = x, 
                                  Parameter = names(x), 
                                  stringsAsFactors = FALSE)) %>%
    bind_rows(.id = "Model") %>%
    mutate(Experiment = gsub("^[^:]+:*", "", header), 
           Experiment = trimws(Experiment), 
           Experiment = firstup(Experiment)) %>%
    select(Experiment, Model, Parameter, MLE) %>%
    write.csv(file = name, row.names = FALSE)
}
