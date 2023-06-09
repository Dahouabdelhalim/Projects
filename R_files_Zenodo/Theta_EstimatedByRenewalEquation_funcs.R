# Author: Chris Wymant, Jan-March 2020

library(ggplot2)
library(tidyverse)
library(rriskDistributions)
library(gridExtra)
library(cubature)

# Abbreviations:
# serint = serial interval. Deprecated. It now in fact refers to generation time.
# incper = incubation period

# Fast and safe redefinition of integrate.
integrate2 <- function(...) {
  result <- tryCatch({integrate(...)$value},
                     error = function(e) cubintegrate(...)$integral)
  result
}

# Get the best fit shape & rate parameters for desribing the lower, central and 
# upper values for each input data parameter as coming from a gamma or lognormal
# distribution.
# If testing, calculate the 2.5th and 97.5th percentiles of the distribution
# after fitting (ideally they would be exactly the values that informed the
# fit, but we are overconstraining the distribution) and the mean too.
fit.pdf.to.param.cis <- function(params, use.lnorm, test.param.pdfs) {
  if (use.lnorm) {
    params <- params %>%
      mutate(lnorm.params = pmap(list(lower, central, upper),
                                 function(x, y, z) {
                                   get.lnorm.par(c(0.025, 0.5, 0.975),
                                                 c(x, y, z),
                                                 show.output = F, plot = F)}))
    params <- unnest(params, lnorm.params) %>%
      group_by(name) %>%
      mutate(col=seq_along(name)) %>%
      spread(key=col, value=lnorm.params) %>%
      rename(meanlog = `1`) %>%
      rename(sdlog = `2`)
    if (test.param.pdfs) {
      params <- params %>%
        mutate(lnorm.mean = pmap(list(meanlog, sdlog), function(m, s) {
          integrate2(function(x) {x * dlnorm(x, m, s)}, lower = 0, upper = Inf)}))
      params <- params %>%
        mutate(lnorm.lower = pmap(list(lower, meanlog, sdlog), function(x, y, z) plnorm(x, y, z)))
      params <- params %>%
        mutate(lnorm.upper = pmap(list(upper, meanlog, sdlog), function(x, y, z) plnorm(x, y, z)))
    }
  } else {
    params <- params %>%
      mutate(gamma.params = pmap(list(lower, central, upper),
                                 function(x, y, z) {
                                   get.gamma.par(c(0.025, 0.5, 0.975),
                                                 c(x, y, z),
                                                 show.output = F, plot = F)}))
    # That assigned two values to the same column. Split into two columns.
    params <- unnest(params, gamma.params) %>%
      group_by(name) %>%
      mutate(col=seq_along(name)) %>%
      spread(key=col, value=gamma.params) %>%
      rename(shape = `1`) %>%
      rename(rate = `2`)
    if (test.param.pdfs) {
      params$gamma.mean <- params$shape / params$rate
      params <- params %>%
        mutate(gamma.lower = pmap(list(lower, shape, rate), function(x, y, z) pgamma(x, y, z)))
      params <- params %>%
        mutate(gamma.upper = pmap(list(upper, shape, rate), function(x, y, z) pgamma(x, y, z)))
      
    }
  }
  params
}


# General functions
serint <- function(tau, serint.shape, serint.scale) {
  dweibull(x = tau, shape = serint.shape, scale = serint.scale)
  # HACK for different functional forms:
  #dgamma(x = tau, shape = serint.shape, rate = serint.scale)
  #dlnorm(x = tau, sdlog = serint.shape, meanlog = serint.scale)
}
prob.asymp <- function(tau, incper.meanlog, incper.sdlog) {
  1 - plnorm(tau, meanlog = incper.meanlog, sdlog = incper.sdlog)
}
prob.symp <- function(tau, incper.meanlog, incper.sdlog) {
  plnorm(tau, meanlog = incper.meanlog, sdlog = incper.sdlog)
}

# Functions specific to this model:

p.el <- function(l, env.decay.rate, env.constant.duration,
                 env.infectiousness.type) {
  if (env.infectiousness.type == "constant") {return(l < env.constant.duration)}
  else if (env.infectiousness.type == "exp.decay") {return(exp(-env.decay.rate * l))}
  else stop(paste0("Unsupported value '", env.infectiousness.type,
                   "' given as env.infectiousness.type argument to p.el function."))
}

model.gen.beta.s.eff.times.1minPa.div.by.R0 <- function(tau, frac.Rs, frac.Rp, P.a, xp,
                                                        incper.meanlog, incper.sdlog,
                                                        serint.shape,
                                                        serint.scale) {
  s.of.tau <- prob.symp(tau = tau, incper.meanlog = incper.meanlog,
                        incper.sdlog = incper.sdlog) 
  func <- (s.of.tau * (frac.Rs + frac.Rp) / (xp * (1 - s.of.tau) + s.of.tau) ) *
    serint(tau = tau, serint.shape = serint.shape,
           serint.scale = serint.scale) 
  func
}

model.gen.beta.p.times.a.1minPa.div.by.R0 <- function(tau, frac.Rs, frac.Rp, P.a, xp,
                                                      incper.meanlog, incper.sdlog,
                                                      serint.shape,
                                                      serint.scale) {
  #if (tau == 0) return(0)
  s.of.tau <- prob.symp(tau = tau, incper.meanlog = incper.meanlog,
                        incper.sdlog = incper.sdlog)
  func <- (xp * (1 - s.of.tau) / s.of.tau) *
    model.gen.beta.s.eff.times.1minPa.div.by.R0(tau = tau, frac.Rs = frac.Rs,
                                                frac.Rp = frac.Rp,
                                                P.a = P.a, xp = xp,
                                                incper.meanlog = incper.meanlog,
                                                incper.sdlog = incper.sdlog,
                                                serint.shape = serint.shape,
                                                serint.scale = serint.scale)
  func
}




model.gen.beta.s.div.by.RSorP <- function(tau, incper.meanlog, incper.sdlog,
                                          serint.shape, serint.scale,
                                          P.a, xp) {
  s.of.tau <- prob.symp(tau = tau, incper.meanlog = incper.meanlog,
                        incper.sdlog = incper.sdlog)  
  serint(tau = tau, serint.shape = serint.shape,
         serint.scale = serint.scale) /
    ((1 - P.a) * (s.of.tau + xp * (1 - s.of.tau)))
  #model.gen.f(tau = tau, incper.meanlog = incper.meanlog,
  #            incper.sdlog = incper.sdlog, P.a = P.a, xp = xp)
}

model.gen.beta.sym.tot <- function(tau, incper.meanlog, incper.sdlog,
                                   serint.shape, serint.scale,
                                   P.a, xp, RSorP) {
  (1 - P.a) * RSorP * model.gen.beta.s.div.by.RSorP(tau = tau,
                                                    incper.meanlog = incper.meanlog,
                                                    incper.sdlog = incper.sdlog,
                                                    serint.shape = serint.shape,
                                                    serint.scale = serint.scale,
                                                    P.a = P.a,
                                                    xp = xp) *
    prob.symp(tau = tau, incper.meanlog = incper.meanlog, incper.sdlog = incper.sdlog)
}

model.gen.beta.presym.tot <- function(tau, incper.meanlog, incper.sdlog,
                                      serint.shape, serint.scale,
                                      P.a, xp, RSorP) {
  xp * (1 - P.a) * RSorP * model.gen.beta.s.div.by.RSorP(tau = tau,
                                                         incper.meanlog = incper.meanlog,
                                                         incper.sdlog = incper.sdlog,
                                                         serint.shape = serint.shape,
                                                         serint.scale = serint.scale,
                                                         P.a = P.a,
                                                         xp = xp) *
    prob.asymp(tau = tau, incper.meanlog = incper.meanlog, incper.sdlog = incper.sdlog)
}


model.gen.beta.env.div.by.E.RSorP <-
  Vectorize(function(tau, incper.meanlog, incper.sdlog, serint.shape,
                     serint.scale, P.a, xp, env.decay.rate,
                     env.constant.duration, env.infectiousness.type) {
    integrate2(function(l) {
      model.gen.beta.s.div.by.RSorP(tau = tau - l,
                                    incper.meanlog = incper.meanlog,
                                    incper.sdlog = incper.sdlog,
                                    serint.shape = serint.shape,
                                    serint.scale = serint.scale,
                                    P.a = P.a,
                                    xp = xp) *
        p.el(l = l, env.decay.rate = env.decay.rate,
             env.constant.duration = env.constant.duration,
             env.infectiousness.type = env.infectiousness.type)},
      lower = 0, upper = tau)
  }, vectorize.args = "tau")

model.gen.f <- function(tau, incper.meanlog, incper.sdlog, P.a, xp) {
  s.of.tau <- prob.symp(tau = tau, incper.meanlog = incper.meanlog,
                        incper.sdlog = incper.sdlog)  
  (1 - P.a) * (s.of.tau + xp * (1 - s.of.tau))
}

model.gen.full.beta.div.by.RSorP <-
  function(tau, serint.shape, serint.scale, P.a, xa, incper.meanlog,
           incper.sdlog, xp, env.scale.constant, env.decay.rate,
           env.constant.duration, env.infectiousness.type) {
    serint(tau = tau, serint.shape = serint.shape,
           serint.scale = serint.scale) *
      (1 + P.a * xa / model.gen.f(tau = tau, incper.meanlog = incper.meanlog,
                                  incper.sdlog = incper.sdlog, P.a = P.a, xp = xp)) +
      env.scale.constant * model.gen.beta.env.div.by.E.RSorP(
        tau = tau, serint.shape = serint.shape,
        serint.scale = serint.scale, incper.meanlog = incper.meanlog,
        incper.sdlog = incper.sdlog, P.a = P.a, xp = xp, env.decay.rate = env.decay.rate,
        env.constant.duration = env.constant.duration,
        env.infectiousness.type = env.infectiousness.type)
  }


model.gen.solve <- function(frac.Re, P.a, doubling.time, xp, xa,
                            incper.meanlog, incper.sdlog, serint.scale,
                            serint.shape, theta.obs, env.decay.rate,
                            env.constant.duration, env.infectiousness.type) {
  
  r <- log(2) / doubling.time # units of per day
  
  integral.of.model.gen.beta.s.div.by.RSorP <-
    integrate2(model.gen.beta.s.div.by.RSorP, lower = 0, upper = Inf,
              incper.meanlog = incper.meanlog, incper.sdlog = incper.sdlog,
              serint.shape = serint.shape, serint.scale = serint.scale,
              P.a = P.a, xp = xp)
  
  env.scale.constant <- (1 + P.a * xa * integral.of.model.gen.beta.s.div.by.RSorP) / 
       (((1 / frac.Re) - 1) * integrate2(model.gen.beta.env.div.by.E.RSorP,
                 lower = 0, upper = Inf,
                 serint.shape = serint.shape,
                 serint.scale = serint.scale,
                 incper.meanlog = incper.meanlog,
                 incper.sdlog = incper.sdlog,
                 P.a = P.a,
                 xp = xp,
                 env.decay.rate = env.decay.rate,
                 env.constant.duration = env.constant.duration,
                 env.infectiousness.type = env.infectiousness.type))
  
  
  RSorP <- 1 / integrate2(function(tau) {
    model.gen.full.beta.div.by.RSorP(tau = tau, serint.shape = serint.shape,
                                     serint.scale = serint.scale,
                                     P.a = P.a, xa = xa, incper.meanlog = incper.meanlog,
                                     incper.sdlog = incper.sdlog, xp = xp,
                                     env.scale.constant = env.scale.constant,
                                     env.decay.rate = env.decay.rate,
                                     env.constant.duration = env.constant.duration,
                                     env.infectiousness.type = env.infectiousness.type) *
      exp(-r * tau)}, lower = 0, upper = Inf)
  
  R0 <- RSorP * integrate2(model.gen.full.beta.div.by.RSorP, lower = 0, upper = Inf,
                          serint.shape = serint.shape,
                          serint.scale = serint.scale,
                          P.a = P.a, xa = xa, incper.meanlog = incper.meanlog,
                          incper.sdlog = incper.sdlog, xp = xp,
                          env.scale.constant = env.scale.constant,
                          env.decay.rate = env.decay.rate,
                          env.constant.duration = env.constant.duration,
                          env.infectiousness.type = env.infectiousness.type)
  
  should.equal.one <- integrate2(function(tau) {exp(-r * tau) * RSorP *
      model.gen.full.beta.div.by.RSorP(tau = tau, serint.shape = serint.shape,
                                       serint.scale = serint.scale,
                                       P.a = P.a, xa = xa, incper.meanlog = incper.meanlog,
                                       incper.sdlog = incper.sdlog, xp = xp,
                                       env.scale.constant = env.scale.constant,
                                       env.decay.rate = env.decay.rate,
                                       env.constant.duration = env.constant.duration,
                                       env.infectiousness.type = env.infectiousness.type)},
      lower = 0, upper = Inf)
  
  RS <- integrate2(model.gen.beta.sym.tot, lower = 0, upper = Inf,
                  incper.meanlog = incper.meanlog,
                  incper.sdlog = incper.sdlog,
                  serint.shape = serint.shape,
                  serint.scale = serint.scale,
                  P.a = P.a, xp = xp, RSorP = RSorP)
  RP <- integrate2(model.gen.beta.presym.tot, lower = 0, upper = Inf,
                  incper.meanlog = incper.meanlog,
                  incper.sdlog = incper.sdlog,
                  serint.shape = serint.shape,
                  serint.scale = serint.scale,
                  P.a = P.a, xp = xp, RSorP = RSorP)
  RA <- RSorP * P.a * xa * integral.of.model.gen.beta.s.div.by.RSorP
  RE <- env.scale.constant * RSorP * integrate2(function(tau) {model.gen.beta.env.div.by.E.RSorP(
    tau = tau, serint.shape = serint.shape,
    serint.scale = serint.scale, incper.meanlog = incper.meanlog,
    incper.sdlog = incper.sdlog, P.a = P.a, xp = xp, env.decay.rate = env.decay.rate,
    env.constant.duration = env.constant.duration,
    env.infectiousness.type = env.infectiousness.type)},
    lower = 0, upper = Inf)
  
  diff.RSorP <- abs(RSorP - RS - RP)
  diff.R0 <- abs(R0 - RS - RP - RA - RE)
  diff.frac.E <- abs(frac.Re - RE / R0)
  diff.unit.integral <- abs(1 - should.equal.one)
  verbose <- FALSE
  if (verbose) {
    cat("Diff between RSorP and RS + RP:", diff.RSorP, "\\n")
    cat("Diff between R0 and RS + RP + RE + RA:", diff.R0, "\\n")
    cat("Diff between frac.E and RE / R0:", diff.frac.E, "(", frac.Re, "c.f.", RE / R0, "\\n")
    cat("Diff between integral of beta(tau) exp(-r tau) and 1:", diff.unit.integral, "\\n")
    cat("R0 = RS + RP + RA + RE\\n")
    cat(R0, "=", RS, "+", RP, "+", RA, "+", RE, "\\n")
  }
  tolerance <- 1e-3
  stopifnot(diff.RSorP < tolerance) #, "RSorP = RS + RP")
  stopifnot(diff.R0 < tolerance) #, "R0 = RS + RP + RE + RA")
  stopifnot(diff.frac.E < tolerance) #, "frac.E = RE / R0")
  stopifnot(diff.unit.integral < tolerance) #, "integral of beta(tau) exp(-r tau) = 1")
  

  
  theta.obs.predicted <- 1 - integrate2(function(tau) {
    model.gen.beta.sym.tot(tau = tau,
                           incper.meanlog = incper.meanlog,
                           incper.sdlog = incper.sdlog,
                           serint.shape = serint.shape,
                           serint.scale = serint.scale,
                           P.a = P.a, xp = xp, RSorP = RSorP) * exp(-r * tau)},
    lower = 0, upper = Inf)

  list(env.scale.constant = env.scale.constant, R0 = R0, RSorP = RSorP, RA = RA,
       RS = RS, RP = RP, RE = RE, theta.obs.predicted = theta.obs.predicted)
  
}
