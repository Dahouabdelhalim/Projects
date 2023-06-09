fit_model <- function(data, model, use_ode = FALSE) {
  if (!(model %in% models())) {
    stop(paste("Model", model, "is not in the model set."))
}
  if (model == "Constant 1") {
    par <- c(beta = 1e-06)
    optim_res <- optim(par, fn = NLL.const1, Np = data$Number, 
                       t = data$time, metacerc = data$Meta, method = "Brent", 
                       lower = 0, upper = 1)
    if (is.infinite(optim_res$value)) {
      print("Infinite NLL. Choosing better boundaries.")
      par <- c(beta = 0.01)
      optim_res <- optim(par, fn = NLL.const1, Np = data$Number, 
                         t = data$time, metacerc = data$Meta, method = "Brent", 
                         lower = 0, upper = .1)
    }
    optim_res$model <- const1
    names(optim_res$par) <- "beta"
  } else if (model == "Constant 2") {
    par <- c(beta = 1e-06)
    # following ensures we don't transmit more cercariae than we
    # have
    maxval <- min(data$Number/(data$host * data$time))
    optim_res <- optim(par, fn = NLL.const2, Np = data$Number, 
                       H = data$host, t = data$time, metacerc = data$Meta, 
                       method = "Brent", lower = 0, upper = maxval)
    optim_res$model <- const2
    names(optim_res$par) <- "beta"
  } else if (model == "Density dependent") {
    par <- c(beta = 1e-05)
    optim_res <- optim(par, fn = NLL.den2, Np = data$Number, 
                       H = data$host, v = data$Volume, t = data$time, metacerc = data$Meta, 
                       method = "Brent", lower = 0, upper = 0.1)
    optim_res$model <- densdep2
    names(optim_res$par) <- "beta"
  } else if (model == "Density independent 1") {
    par <- c(beta = 1e-03)
    optim_res <- optim(par, fn = NLL.freq1, Np = data$Number, 
                       H = data$host, t = data$time, metacerc = data$Meta, 
                       method = "Brent", lower = 0, upper = .8)
    if (is.infinite(optim_res$value)) {
      print("Infinite NLL. Choosing better boundaries.")
      par <- c(beta = 0.01)
      optim_res <- optim(par, fn = NLL.freq1, Np = data$Number, 
                         H = data$host, t = data$time, metacerc = data$Meta, 
                         method = "Brent", lower = 0, upper = .1)
    }
    optim_res$model <- freqdep1
    names(optim_res$par) <- "beta"
  } else if (model == "Density independent 2") {
    if (use_ode){
      par <- c(beta = 1)
      optim_res <- optim(par, fn = NLL.freq2_ode, Np = data$Number, H = data$host, 
                       v = data$Volume, t = data$time, metacerc = data$Meta, method = "Brent", 
                       lower = 0, upper = 1.5)
      optim_res$model <- fd2_ode
    } else {
      par <- c(beta = 1e-06)
      optim_res <- optim(par, fn = NLL.freq2, Np = data$Number, 
                         H = data$host, v = data$Volume, t = data$time, 
                         metacerc = data$Meta, 
                         method = "Brent", lower = 0, upper = 0.2)
      optim_res$model <- freqdep2
    }
    names(optim_res$par) <- "beta"
  } else if (model == "Power C") {
    par <- c(beta = 1, q = 2)
    optim_res <- optim(par, fn = NLL.powC, Np = data$Number, 
                       H = data$host, t = data$time, metacerc = data$Meta)
    optim_res$model <- powerC
  } else if (model == "Power H") {
    par <- c(beta = 1e-04, p = 1)
    optim_res <- optim(par, fn = NLL.powH, Np = data$Number, 
                       H = data$host, t = data$time, metacerc = data$Meta, 
                       control = list(maxit = 1e+06))
    optim_res$model <- powerH
  } else if (model == "Power CH") {
    par <- c(beta = 1, q = 2, p = 2)
    optim_res <- optim(par, fn = NLL.powCH, Np = data$Number, 
                       H = data$host, t = data$time, metacerc = data$Meta)
    optim_res$model <- powerCH
  } else if (model == "Negative binomial 1") {
    par <- c(beta = .5, k = .5)
    optim_res <- optim(par, fn = NLL.NB1, Np = data$Number, 
                       t = data$time, metacerc = data$Meta)
    optim_res$model <- NB1
  } else if (model == "Negative binomial 2") {
    par <- c(beta = .5, k = .5)
    optim_res <- optim(par, fn = NLL.NB2, Np = data$Number, 
                       H = data$host, t = data$time, metacerc = data$Meta)
    optim_res$model <- NB2
  }
  class(optim_res) <- c("optim_res", "list")
  optim_res
}