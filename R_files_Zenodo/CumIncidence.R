# similar to cuminc function but allows more control over time points
"CumIncidence" <- function(ftime, fstatus, group, t, strata, rho = 0, 
                           cencode = 0, subset, na.action = na.omit, level,
                           xlab = "Time", ylab = "Probability", 
                           col, lty, lwd, digits = 4)
{
  # check for the required package
  if(!require("cmprsk"))
  { stop("Package `cmprsk' is required and must be installed.\\n 
           See help(install.packages) or write the following command at prompt
           and then follow the instructions:\\n
           > install.packages(\\"cmprsk\\")") } 
  # collect data
  mf  <- match.call(expand.dots = FALSE)
  mf[[1]] <- as.name("list")
  mf$t <- mf$digits <- mf$col <- mf$lty <- mf$lwd <- mf$level <- 
    mf$xlab <- mf$ylab <- NULL
  mf <- eval(mf, parent.frame())
  g <- max(1, length(unique(mf$group)))
  s <- length(unique(mf$fstatus))
  if(missing(t)) 
  { time <- pretty(c(0, max(mf$ftime)), 6)
  ttime <- time <- time[time < max(mf$ftime)] }
  else { ttime <- time <- t }
  # fit model and estimates at time points
  fit   <- do.call("cuminc", mf)
  tfit <- timepoints(fit, time)
  # print result
  cat("\\n+", paste(rep("-", 67), collapse=""), "+", sep ="")
  cat("\\n| Cumulative incidence function estimates from competing risks data |")
  cat("\\n+", paste(rep("-", 67), collapse=""), "+\\n", sep ="")
  tests <- NULL
  if(g > 1)
  { 
    tests <- data.frame(fit$Tests[,c(1,3,2)], check.names = FALSE)
    colnames(tests) <- c("Statistic", "df", "p-value")
    tests$`p-value` <- format.pval(tests$`p-value`)
    cat("Test equality across groups:\\n")
    print(tests, digits = digits) 
  }
  cat("\\nEstimates at time points:\\n")
  print(tfit$est, digits = digits)
  cat("\\nStandard errors:\\n")
  print(sqrt(tfit$var), digits = digits)
  #
  if(missing(level))
  { # plot cumulative incidence functions
    if(missing(t))
    { time <- sort(unique(c(ftime, time)))
    x <- timepoints(fit, time) }
    else x <- tfit
    col <- if(missing(col)) rep(1:(s-1), rep(g,(s-1))) else col
    lty <- if(missing(lty)) rep(1:g, s-1) else lty
    lwd <- if(missing(lwd)) rep(1, g*(s-1)) else lwd      
    matplot(time, base::t(x$est), type="s", ylim = c(0,1), 
            xlab = xlab, ylab = ylab, xaxs="i", yaxs="i", 
            col = col, lty = lty, lwd = lwd)
    legend("topleft", legend =  rownames(x$est), x.intersp = 2, 
           bty = "n", xjust = 1, col = col, lty = lty, lwd = lwd)
    out <- list(test = tests, est = tfit$est, se = sqrt(tfit$var))
  }
  else
  { if(level < 0 | level > 1) 
    error("level must be a value in the range [0,1]")
    # compute pointwise confidence intervals
    oldpar <- par(ask=TRUE)
    on.exit(par(oldpar))
    if(missing(t))
    { time <- sort(unique(c(ftime, time)))
    x <- timepoints(fit, time) }
    else x <- tfit
    z <- qnorm(1-(1-level)/2)
    lower <- x$est ^ exp(-z*sqrt(x$var)/(x$est*log(x$est)))
    upper <- x$est ^ exp(z*sqrt(x$var)/(x$est*log(x$est)))
    col <- if(missing(col)) rep(1:(s-1), rep(g,(s-1))) 
    else             rep(col, g*(s-1))
    lwd <- if(missing(lwd)) rep(1, g*(s-1)) 
    else             rep(lwd, g*(s-1))      
    # plot pointwise confidence intervals
    for(j in 1:nrow(x$est))
    { matplot(time, cbind(x$est[j,], lower[j,], upper[j,]), type="s", 
              xlab = xlab, ylab = ylab, xaxs="i", yaxs="i", 
              ylim = c(0,1), col = col[j], lwd = lwd[j], lty = c(1,3,3))
      legend("topleft", legend =  rownames(x$est)[j], bty = "n", xjust = 1) }
    # print pointwise confidence intervals
    i <- match(ttime, time)
    ci <- array(NA, c(2, length(i), nrow(lower)))
    ci[1,,] <- base::t(lower[,i])
    ci[2,,] <- base::t(upper[,i])
    dimnames(ci) <- list(c("lower", "upper"), ttime, rownames(lower))
    cat(paste("\\n", level*100, "% pointwise confidence intervals:\\n\\n", sep=""))
    print(ci, digits = digits)
    out <- list(test = tests, est = x$est, se = sqrt(tfit$var), ci = ci)
  }
  # return results
  invisible(out)
}