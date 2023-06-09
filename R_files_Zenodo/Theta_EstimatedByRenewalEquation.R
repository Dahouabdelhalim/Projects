#####################################
# This code solves for the different 
# functional forms contributing to from specified
# beta(tau) in our model of infectiousness.

# (Ferretti, Wymant et al, Science 2020)
#####################################

library(ggplot2)
library(tidyverse)
library(rriskDistributions)
library(gridExtra)
library(RColorBrewer)
library(R.utils)

rm(list=ls())

# Abbreviations:
# serint = serial interval. Deprecated. It now in fact refers to generation time.
# incper = incubation period

# TODO: set this to the directory you want to work in. Make sure
# Theta_EstimatedByRenewalEquation_funcs.R is in there.
ncov.dir <- "" 

setwd(ncov.dir)
source("Theta_EstimatedByRenewalEquation_funcs.R")

# A tibble for parameters - their central values, and lower and upper CIs (for
# when we assign CIs to our calculations.)
data.params <- tibble(name = character(0),
                      central = numeric(0),
                      lower   = numeric(0),
                      upper   = numeric(0))
################################################################################
# INPUT DATA
# Add each parameter and 95% CIs as a new row to the data.params data frame.

# Doubling time. Units = days. Source: in-house fit to data from Li et al NEJM
# and China CDC.
data.params <- bind_rows(data.params,
                         list(name = "doubling.time", central = 4.951051, 
                              lower = 4.244315, upper = 6.362081)
)

# Incubation period parameters, from
# https://github.com/HopkinsIDD/ncov_incubation#parameter-estimates
# All have units of days
data.params <- bind_rows(data.params,
                         list(name = "incper.meanlog", central = 1.644,
                              lower = 1.495, upper = 1.798)
)
data.params <- bind_rows(data.params,
                         list(name = "incper.sdlog", central = 0.363,
                              lower = 0.201, upper = 0.521)
)

# All have units of days
data.params <- bind_rows(data.params,
                         list(name = "serint.shape", central = 2.826, # Weibull shape = 2.826, Gamma shape = 5.683, lnorm "shape" (sdlog) = 0.44
                              lower = 1.75, upper = 4.7)
) 
data.params <- bind_rows(data.params,
                         list(name = "serint.scale",  central = 5.665, 
                              lower = 4.7, upper = 6.9) 
) 


# Parameter not varied: Theta observed. Source: Li et al NEJM, merging 
# periods two and three, removing from the denominator people exposed to the market.
theta.obs <- (141 + 59) / (164 + 76)

# Params govered by priors, 
xp <- 1;
P.a <- 0.40; P.a.alpha <- 1.5; P.a.beta <- 1.75
xa <- 0.1; xa.alpha <- 1.5; xa.beta <- 5.5
frac.Re <- 0.1; frac.Re.alpha <- 1.5; frac.Re.beta <- 5.5
env.infectiousness.type <- "constant" # "exp.decay" or "constant": binary choice of the shape of intensinty of transmission from the environment after it's contaminated:
env.constant.duration <- 3; env.constant.duration.shape <- 4; env.constant.duration.rate <- 1
env.decay.rate <- log(10) # # prior params not defined, this is no longer used. log(10) means it decays by a factor 10 per day - that's the most intuitive way to frame it probably.

# Plotting params:
tau.test <- seq(from=0, to=13, by=0.05) # x axis 
colors <- brewer.pal(4, "Paired")
names(colors) <- c("pre-symptomatic", "symptomatic", "environmental", "asymptomatic")
base.font.size <- 27
plot.height <- 9
plot.width <- 15
minor_breaks <- seq(0, 13, 2)
breaks <- seq(0, 13, 2)

################################################################################

# Get the median of the serint distribution for file names. 
serint.median <- qweibull(p = 0.5, scale =
                            data.params[data.params$name == "serint.scale", ]$central,
                          shape = data.params[data.params$name == "serint.shape", ]$central)
#serint.median <- qlnorm(p = 0.5, meanlog = data.params[data.params$name == "serint.scale", ]$central, sdlog = data.params[data.params$name == "serint.shape", ]$central)

#plot(tau.test, serint(tau = tau.test,
#                      serint.shape = data.params[data.params$name == "serint.shape", ]$central,
#                      serint.scale = data.params[data.params$name == "serint.scale", ]$central))

# Get the best fit shape & rate parameters for desribing the lower, central and 
# upper values for each input data parameter as coming from a gamma or lognormal
# distribution.
# If testing, calculate the 2.5th and 97.5th percentiles of the distribution
# after fitting (ideally they would be exactly the values that informed the
# fit, but we are overconstraining the distribution) and the mean too.
use.lnorm <- T
test.param.pdfs <- F
data.params <- fit.pdf.to.param.cis(data.params, use.lnorm, test.param.pdfs)

# From now enforce that only lognorm is used, for convenience.
stopifnot(use.lnorm)

plot.priors <- T
if (plot.priors) {
  p4 <- ggplot(data = tibble(x = seq(0, 15, by = 0.1),
                             y = dgamma(x, shape = env.constant.duration.shape,
                                        rate = env.constant.duration.rate))) +
    labs(x = "Duration of infectiousness of\\ncontaminated environment (days)",
         y = "Prior probability density") +
    coord_cartesian(ylim = c(0, 0.25), expand = F) +
    theme_classic() +
    geom_line(aes(x = x, y = y), col = "blue")
  p2 <- ggplot(data = tibble(x = seq(0, 1, by = 0.01),
                             y = dbeta(x, shape1 = xa.alpha, shape2 = xa.beta))) +
    labs(x = "Relative infectiousness of\\nasymptomatic indivudals",
         y = "Prior probability density") +
    coord_cartesian(ylim = c(0, 3.2), expand = F) +
    theme_classic() +
    geom_line(aes(x = x, y = y), col = "blue")
  p3 <- ggplot(data = tibble(x = seq(0, 1, by = 0.01),
                             y = dbeta(x, shape1 = frac.Re.alpha, shape2 = frac.Re.beta))) +
    labs(x = "Fraction of all transmission\\nthat's environmentally mediated",
         y = "Prior probability density") +
    coord_cartesian(ylim = c(0, 3.2), expand = F) +
    theme_classic() +
    geom_line(aes(x = x, y = y), col = "blue")
  p1 <- ggplot(data = tibble(x = seq(0, 1, by = 0.01),
                             y = dbeta(x, shape1 = P.a.alpha, shape2 = P.a.beta))) +
    labs(x = "Fraction of infected individuals\\nwho are asymptomatic",
         y = "Prior probability density") +
    coord_cartesian(ylim = c(0, 1.5), expand = F) +
    theme_classic() +
    geom_line(aes(x = x, y = y), col = "blue")
  library(gridExtra)
  ptotal <- grid.arrange(p1, p2, p3, p4, ncol=2)
  ggsave(paste("priors.pdf",sep=""), ptotal, height = 7, width = 7) 
}

# Explore the CIs of the data params and the priors of the other params?
uncertainty.analysis <- F
if (uncertainty.analysis) {
  
  # Simulate parameter sets, or analyse previously simulated ones?
  # If simulate.dont.analyse, run this script from the command line with
  # for i in $(seq 1 100); do timeout 90 Rscript Theta_EstimatedByRenewalEquation.R; done
  # Or a sequence of different length than 100. Can run that command many
  # times in parallel - each will generate a set of files with random
  # parameters. For a small minority of
  # parameter sets, the numerical integration gets stuck forever, and cannot
  # be interrupted from within R (withTimeout, setTimeLimit can't properly
  # interrupt native code). Hence the need to kill from the command line
  # and try again.
  simulate.dont.analyse <- F
  if (simulate.dont.analyse) {
    
    # Populate a df with one row per random draw and one col per parameter
    num.draws <- 50 # Run this many times to get many sets of 50.
    results <- as_tibble(matrix(ncol = 0, nrow = num.draws))
    # Draw data param columns from previously fitted lognormals: 
    for (param.num in seq(from=1, to=nrow(data.params))) {
      results[[data.params$name[[param.num]]]] <-
        rlnorm(num.draws, meanlog = data.params$meanlog[[param.num]],
               sdlog = data.params$sdlog[[param.num]])
    }
    # Draw unconstrained params from their priors:
    results$env.constant.duration <- rgamma(n = num.draws,
                                            shape = env.constant.duration.shape,
                                            rate = env.constant.duration.rate)
    results$xa <- rbeta(n = num.draws, shape1 = xa.alpha, shape2 = xa.beta)
    results$frac.Re <- rbeta(n = num.draws, shape1 = frac.Re.alpha, shape2 = frac.Re.beta)
    results$P.a <- rbeta(n = num.draws, shape1 = P.a.alpha, shape2 = P.a.beta)
    ## Double check drawn values look right:
    #hist(results$serint.shape, breaks = 50, freq = F)
    #x <- seq(1, 10, by = 0.01)
    #lines(x, dlnorm(x, meanlog = data.params[data.params$name == "serint.shape", ]$meanlog, sdlog = data.params[data.params$name == "serint.shape", ]$sdlog))
    
    # Our results df currently has one col per input param. Initialise one empty col per output param.
    input.param.names <- names(results)
    output.param.names <- c("R0", "theta.obs.predicted", "RA", "RS", "RP", "RE")
    for (name in output.param.names) results[[name]] <- NA
    
    # Now for each draw, take that random set of input params and calculate the results
    for (draw in seq(num.draws)) {
      if (draw %% 10 == 0) cat("Now considering draw", draw,
                               "using these args:\\nfrac.Re = ",frac.Re,
                               "P.a =", P.a, "doubling.time = ",doubling.time,
                               "xp = ", xp, "xa = ", xa, "incper.meanlog = ", incper.meanlog,
                               "incper.sdlog = ", incper.sdlog,
                               "serint.scale = ", serint.scale,
                               "serint.shape = ", serint.shape, 
                               "theta.obs = ", theta.obs, "env.decay.rate = ", env.decay.rate,
                               "env.constant.duration = ", env.constant.duration, 
                               "env.infectiousness.type = ", env.infectiousness.type, "\\n")
      cat("Now considering draw", draw, "\\n")
      
      # For convenience, define variables e.g. called foo instead of
      # pararams[data.params$name == "foo", ]$central
      for (param.name in input.param.names) assign(param.name, results[draw,][[param.name]])
      
      dummy <- tryCatch({dummy <-
        withTimeout({model.gen.solve(frac.Re = frac.Re, P.a = P.a, doubling.time = doubling.time,
                                     xp = xp, xa = xa, incper.meanlog = incper.meanlog,
                                     incper.sdlog = incper.sdlog,
                                     serint.scale = serint.scale,
                                     serint.shape = serint.shape, 
                                     theta.obs = theta.obs, env.decay.rate = env.decay.rate,
                                     env.constant.duration = env.constant.duration, 
                                     env.infectiousness.type = env.infectiousness.type)}, timeout = 3)},
        error = function(e) {
          dummy <- rep(NA, length(output.param.names))
          names(dummy) <- output.param.names
          dummy
        })
      results[draw, output.param.names] <- dummy[output.param.names]
      
    }
    if (anyNA(results$R0)) {
      results <- results[seq(min(which(is.na(results$R0))) - 1), ]
    }
    write_csv(results, paste0("GeneralModelUncertaintyAnalysis_run_", Sys.time()))
    quit("no")
    
    # This next scope if we are analysing previously simulated parameter sets:
  } else {
    
    # Read in all the parameter output files and merge.
    files <- list.files(path = "./",
                        full.names = TRUE,
                        pattern = "GeneralModelUncertaintyAnalysis_*")
    results <- NULL
    for (file in files) {
      this.df <- read_csv(file)
      if (is.null(results)) {
        results <- this.df
      } else {
        results <- bind_rows(results, this.df)
      }
    }
    
    # Got one NA, unexpected, weird
    results <- results[complete.cases(results), ]
    
    # Be certain no parameter set was duplicated
    results <- results[!duplicated(results),]
    
    # Keep only 10k results for niceness in reporting.
    results <- results[1:min(nrow(results), 10000),] # min() inserted just for you, user!
    
    # Sanity check that it looks right
    hist(results$serint.shape, breaks = 50, freq = F)
    x <- seq(1, 10, by = 0.01)
    lines(x, dlnorm(x, meanlog = data.params[data.params$name == "serint.shape", ]$meanlog, sdlog = data.params[data.params$name == "serint.shape", ]$sdlog))
  }
  
  
  # Calculate the two separate contributions to R0 for each parameter draw
  #results <- results %>% mutate(Rp = pmap_dbl(list(xp, R0, serint.shape,
  #                                                 serint.scale),
  #                                            rp.model.presymp))
  #results <- results %>% mutate(Rs = pmap_dbl(list(xp, R0, serint.shape,
  #                                                 serint.scale),
  #                                            rs.model.presymp))
  
  # Calculate derived things from the reported model output
  results$frac.Rs <- results$RS / results$R0
  results$frac.Rp <- results$RP / results$R0
  results$frac.Ra <- results$RA / results$R0
  results$theta <- 1 - results$RS / results$R0
  results$Rp.div.by.RsRp <- results$RP / (results$RS + results$RP)
  results
  Rp.confints <- quantile(results$RP, c(0.025, 0.975, 0.5))
  Rp.confints
  Rs.confints <- quantile(results$RS, c(0.025, 0.975, 0.5))
  Rs.confints
  Ra.confints <- quantile(results$RA, c(0.025, 0.975, 0.5))
  Ra.confints
  Re.confints <- quantile(results$RE, c(0.025, 0.975, 0.5))
  Re.confints
  R0.confints <- quantile(results$R0, c(0.025, 0.975, 0.5))
  R0.confints
  frac.Rp.confints <- quantile(results$frac.Rp, c(0.025, 0.975, 0.5))
  frac.Rp.confints
  frac.Rs.confints <- quantile(results$frac.Rs, c(0.025, 0.975, 0.5))
  frac.Rs.confints
  frac.Ra.confints <- quantile(results$frac.Ra, c(0.025, 0.975, 0.5))
  frac.Ra.confints
  frac.Re.confints <- quantile(results$frac.Re, c(0.025, 0.975, 0.5))
  frac.Re.confints
  Rp.div.by.RsRp.confints <-  quantile(results$Rp.div.by.RsRp, c(0.025, 0.975, 0.5))
  Rp.div.by.RsRp.confints
  theta.obs.predicted.confints <- quantile(results$theta.obs.predicted, c(0.025, 0.975, 0.5))
  theta.obs.predicted.confints
  
  # Plot heat maps of R0 contribution correlations
  p1 <- ggplot(results, aes(x=RS, y=RP) ) +
    geom_bin2d(bins = 40, aes(fill = ..density..)) +
    scale_fill_continuous(type = "viridis") +
    theme_classic() + 
    coord_cartesian(expand = F) + 
    labs(x = expression("R"["S"]),
         y = expression("R"["P"]))
  p2 <- ggplot(results, aes(x=RS, y=RE) ) +
    geom_bin2d(bins = 40, aes(fill = ..density..)) +
    scale_fill_continuous(type = "viridis") +
    theme_classic() + 
    coord_cartesian(expand = F) + 
    labs(x = expression("R"["S"]),
         y = expression("R"["E"]))
  p3 <- ggplot(results, aes(x=RS, y=RA) ) +
    geom_bin2d(bins = 40, aes(fill = ..density..)) +
    scale_fill_continuous(type = "viridis") +
    theme_classic() + 
    coord_cartesian(expand = F) + 
    labs(x = expression("R"["S"]),
         y = expression("R"["A"]))
  p4 <- ggplot(results, aes(x=RS, y=R0) ) +
    geom_bin2d(bins = 40, aes(fill = ..density..)) +
    scale_fill_continuous(type = "viridis") +
    theme_classic() + 
    coord_cartesian(expand = F) + 
    labs(x = expression("R"["S"]),
         y = expression("R"["0"]))
  ptotal <- grid.arrange(p1, p2, p3, p4, ncol=2)
  ggsave("RS_RPetc_correlation.pdf", ptotal, height = 7, width = 8)
  
  #xmax <- max(max(results$RA), max(results$RS), max(results$RE), max(results$RP))
  
  # Set an overflow bin (manually reported in paper)
  xmax <- 1.8
  results$RS <- pmin(results$RS, xmax)
  results$RA <- pmin(results$RA, xmax)
  results$RP <- pmin(results$RP, xmax)
  results$RE <- pmin(results$RE, xmax)
  
  # Plot distribution of contributions to R0
  ggplot(results) + geom_histogram(aes(RP, y=..density..), color = colors[["pre-symptomatic"]], fill = colors[["pre-symptomatic"]]) + labs(x = expression("R"["P"]), y = "posterior density") + theme_classic(base_size = 25) + coord_cartesian(xlim=c(0, xmax), expand = F)
  
  p <- grid.arrange(nrow = 2, #ggplot(results) + geom_histogram(aes(theta)),
                    ggplot(results) + geom_histogram(aes(RP, y=..density..), color = colors[["pre-symptomatic"]], fill = colors[["pre-symptomatic"]]) + labs(x = expression("R"["P"]), y = "posterior density") + theme_classic(base_size = 25) + coord_cartesian(xlim=c(0, xmax), expand = F),
                    ggplot(results) + geom_histogram(aes(RS, y=..density..), color = colors[["symptomatic"]], fill = colors[["symptomatic"]]) + labs(x = expression("R"["S"]), y = "posterior density") + theme_classic(base_size = 25) + coord_cartesian(xlim=c(0, xmax),expand = F) ,
                    ggplot(results) + geom_histogram(aes(RE, y=..density..), color = colors[["environmental"]], fill = colors[["environmental"]]) + labs(x = expression("R"["E"]), y = "posterior density") + theme_classic(base_size = 25) + coord_cartesian(xlim=c(0, xmax),expand = F) ,
                    ggplot(results) + geom_histogram(aes(RA, y=..density..), color = colors[["asymptomatic"]], fill = colors[["asymptomatic"]]) + labs(x = expression("R"["A"]), y = "posterior density") + theme_classic(base_size = 25) + coord_cartesian(xlim=c(0, xmax),expand = F) ,
                    ggplot(results) + geom_histogram(aes(R0, y=..density..)) + labs(x = expression("R"["0"]), y = "posterior density") + theme_classic(base_size = 25) + coord_cartesian(expand = F) ,
                    ggplot(results) + geom_histogram(aes(frac.Rp, y=..density..), color = colors[["pre-symptomatic"]], fill = colors[["pre-symptomatic"]]) + labs(x = expression("R"["P"]*"/R"["0"]), y = "posterior density") + theme_classic(base_size = 25) + coord_cartesian(xlim=c(0, 1), expand = F) ,
                    ggplot(results) + geom_histogram(aes(frac.Rs, y=..density..), color = colors[["symptomatic"]], fill = colors[["symptomatic"]]) + labs(x = expression("R"["S"]*"/R"["0"]), y = "posterior density") + theme_classic(base_size = 25) + coord_cartesian(xlim=c(0, 1), expand = F) ,
                    ggplot(results) + geom_histogram(aes(frac.Re, y=..density..), color = colors[["environmental"]], fill = colors[["environmental"]]) + labs(x = expression("R"["E"]*"/R"["0"]), y = "posterior density") + theme_classic(base_size = 25) + coord_cartesian(xlim=c(0, 1), expand = F),
                    ggplot(results) + geom_histogram(aes(frac.Ra, y=..density..), color = colors[["asymptomatic"]], fill = colors[["asymptomatic"]]) + labs(x = expression("R"["A"]*"/R"["0"]), y = "posterior density") + theme_classic(base_size = 25) + coord_cartesian(xlim=c(0, 1), expand = F))
  ggsave(p, filename = "GeneralModel_CIs.pdf", height = 10, width = 25)
  
  # Plot RP/(RP+RS)
  p2 <- ggplot(results) + geom_histogram(aes(Rp.div.by.RsRp, y=..density..)) +
    labs(x = expression("R"["P"]*"/(R"["P"]*"+R"["S"]*")"), y = "posterior density") +
    theme_classic(base_size = 25) + coord_cartesian(expand = F) +
    xlim(c(0, 0.9))
  ggsave(p2, filename = "GeneralModel_RPdivbyRP+RS.pdf", height = 6, width = 7)
  
  # Next scope if we're not doing uncertainty analysis, just the central estimate
  # Make sure you don't run this code after running the scope above (i.e
  # respect the if-else set up) or at least one parameter is mis-set.
} else {
  
  # For convenience, define variables e.g. called foo instead of
  # pararams[data.params$name == "foo", ]$central
  for (param.num in seq(from=1, to=nrow(data.params))) {
    assign(data.params$name[[param.num]], data.params$central[[param.num]])
  }
  
  # Solve!
  dummy <- model.gen.solve(frac.Re = frac.Re, P.a = P.a, doubling.time = doubling.time,
                           xp = xp, xa = xa, incper.meanlog = incper.meanlog,
                           incper.sdlog = incper.sdlog,
                           serint.scale = serint.scale,
                           serint.shape = serint.shape, 
                           theta.obs = theta.obs, env.decay.rate = env.decay.rate,
                           env.constant.duration = env.constant.duration, 
                           env.infectiousness.type = env.infectiousness.type)
  env.scale.constant <- dummy$env.scale.constant
  R0 <- dummy$R0
  RSorP <- dummy$RSorP
  theta.obs.predicted <- dummy$theta.obs.predicted # the theta you would observe distorted by exponentially growing dynamics

  RA <- dummy$RA # contribution to R0 from asymps
  RS <- dummy$RS # contribution to R0 from symps
  RP <- dummy$RP # contribution to R0 from pre-symps
  RE <- dummy$RE # contribution to R0 from environment
  theta.true <- 1 - RS / R0 # the theta that's not distorted by dynamics
  
  RP / R0
  RS / R0
  RE / R0
  RA / R0
  RP / (RP + RS)
  
  # Plot the different contributions to beta(tau):
  
  df.beta.p <- data.frame(tau = tau.test, label = "pre-symptomatic", beta = vapply(
    tau.test, model.gen.beta.presym.tot, numeric(1), incper.meanlog = incper.meanlog,
    incper.sdlog = incper.sdlog, serint.shape = serint.shape,
    serint.scale = serint.scale, P.a = P.a, xp = xp, RSorP = RSorP))
  
  df.beta.s <- data.frame(tau = tau.test, label = "symptomatic", beta = vapply(
    tau.test, model.gen.beta.sym.tot, numeric(1), incper.meanlog = incper.meanlog,
    incper.sdlog = incper.sdlog, serint.shape = serint.shape,
    serint.scale = serint.scale, P.a = P.a, xp = xp, RSorP = RSorP))
  
  df.beta.a <- data.frame(tau = tau.test, label = "asymptomatic", beta = vapply(
    tau.test, function(tau) {RSorP * P.a * xa *
        model.gen.beta.s.div.by.RSorP(tau = tau,
                                      incper.meanlog = incper.meanlog, 
                                      incper.sdlog = incper.sdlog,
                                      serint.shape = serint.shape,
                                      serint.scale = serint.scale,
                                      P.a = P.a, xp = xp)} , numeric(1)))
  
  df.beta.e <- data.frame(tau = tau.test, label = "environmental", beta = vapply(
    tau.test, model.gen.beta.env.div.by.E.RSorP, numeric(1),
    serint.shape = serint.shape,
    serint.scale = serint.scale, 
    incper.meanlog = incper.meanlog,
    incper.sdlog = incper.sdlog, P.a = P.a,
    xp = xp, env.decay.rate = env.decay.rate,
    env.constant.duration = env.constant.duration,
    env.infectiousness.type = env.infectiousness.type))
  df.beta.e$beta <- df.beta.e$beta * env.scale.constant * RSorP
  
  df.plot <- rbind(df.beta.a, df.beta.p, df.beta.s, df.beta.e)
  df.plot$label <- factor(df.plot$label, levels = c("pre-symptomatic",
                                                    "symptomatic",
                                                    "environmental",
                                                    "asymptomatic"))
  
  # Go from long to wide format for ease of finding the maximum stacked value
  df.plot.wide <- spread(df.plot, label, beta)
  df.plot.wide$max <- df.plot.wide$`pre-symptomatic` + df.plot.wide$symptomatic +
    df.plot.wide$environmental + df.plot.wide$asymptomatic
  ymax <- max(df.plot.wide$max * 1.05)
  
  p <- ggplot(df.plot, aes(x=tau, y=beta)) + #, color=label)) +
    theme_bw(base_size = base.font.size) +
    geom_area(aes(fill=label)) +
    labs(x = expression(paste(tau, " (days)")),
         y = expression(paste(beta, "(", tau, ")  (transmissions per day)")),
         fill = bquote(paste('R'['0'] * ' = ' * .(format(round(R0, 1), nsmall = 1)) * ":" )),
         title = NULL,  #bquote(atop('Assumptions: median generation time = ' *
         subtitle = NULL) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    coord_cartesian(xlim = c(0, max(tau.test)), ylim = c(0, ymax), expand = F) +
    scale_fill_manual(values = colors[c("pre-symptomatic", "symptomatic", "environmental", "asymptomatic")],
                      labels = c(
                        bquote(paste('R'['p'] * ' = ' * .(format(round(RP, 1), nsmall = 1)) * " from pre-symptomatic")),
                        bquote(paste('R'['s'] * ' = ' * .(format(round(RS, 1), nsmall = 1)) * " from symptomatic")),
                        bquote(paste('R'['e'] * ' = ' * .(format(round(RE, 1), nsmall = 1)) * " from environmental")),
                        bquote(paste('R'['a'] * ' = ' * .(format(round(RA, 1), nsmall = 1)) * " from asymptomatic")))) +
    scale_x_continuous(minor_breaks = minor_breaks, breaks = breaks)
  p
  plot.name <- paste0("ModelGeneral_MedianSerint_", format(round(serint.median, 1), nsmall = 1), "_xp_", xp, ".pdf")
  ggsave(plot.name, p, height = plot.height, width = plot.width)
}
