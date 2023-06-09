# Functions used by "VanBuskirk_Smith_2021_Rscript.R"
#
#
set.up.interval <- function(rnge, fine.grid) {
  # This function sets up a sequence that includes 0 in the middle.
  # It is used when creating a "newdata" data frame for predicting from a fitted model.
  interval <- (rnge[2] - rnge[1]) / fine.grid
  if (sum(sign(rnge)) == 0) {
    s1   <- seq(from = 0, to = (rnge[2] + interval), by = interval)
    s2   <- seq(from = -interval, to = (rnge[1] - interval), by = -interval)
    s    <- c(s1, s2)
    s    <- s[order(s)]
    s[1] <- rnge[1]
    s[length(s)] <- rnge[2]
    } else {
    s <- seq(from = rnge[1], to = rnge[2], by = interval)
    }
  return(s)
  }
#
#
make.dates <- function(x) {
  # This function returns the day and month for the numerical date system
  # used in the Isle Royale brood dataset (in that datset, 1 = 1 May).
  if(x < 32) output <- list(x, "May")
  if(x > 31) output <- list(x-31, "Jun")
  if(x > 61) output <- list(x-61, "Jul")
  if(x > 92) output <- list(x-92, "Aug")
  return(output)
  }
#
#
Chevin.2015 <- function(ps, meanz.t, allz.t, outliers = FALSE, threshold = 2) {
  # This function executes Eqn 4 in Chevin et al. (2015, Evolution).
  # Eqn 3 must have been fit already in brms.
  # For each iteration, it calculates the following parameters:
  #  * omega is the width of the fitness function.
  #  * theta is the phenotypic optimum for each year.
  #  * Wmax is the height if the fitness function for each year.
  #  * total.load.t is the proportional decline in fitness in year t due to the deviation of
  #       the observed trait from theta; this is calculated for each brood and then averaged.
  #  * lag.load.t is the proportional load caused by deviation of population mean from theta.
  #  * beta is the directional selection coefficient for each year.
  #  * fluct.load is a single measure of load for each iteration, due to fluctuations in
  #       theta (de Villemereuil et al. 2020, PNAS).
  #
  # Arguments:
  #   ps        = posterior samples from the brms model.
  #   meanz.t   = data frame with years and mean trait value in each year (15 rows, 2 columns).
  #   allz.t    = data frame with z values for individual broods. Trait should be in column 2.
  #   outliers  = TRUE/FALSE indicating whether "remove.outliers" should be run.
  #   threshold = SD units above/below the median to remove.
  #
  N                   <- length(ps[,1])
  years               <- c(1983,1985:1998)
  trait               <- unlist(strsplit(names(ps)[2], "_"))[2]
  two.beta.zz         <- 2 * ps[,3]
  omega               <- suppressWarnings(sqrt(-1 / two.beta.zz ))
  omega[is.na(omega)] <- 200             # if it's not stabilizing selection, assign a broad and flat fitness function.
  theta.t             <- beta.t <- Wmax.t <- lag.load.t <- total.load.t <- as.data.frame(matrix(0, nrow = N, ncol = length(years)))
  names(theta.t)      <- paste("theta.", years, sep = "")
  names(Wmax.t)       <- paste("Wmax.", years, sep = "")
  names(total.load.t) <- paste("load.all.", years, sep = "")
  names(lag.load.t)   <- paste("load.", years, sep = "")
  names(beta.t)       <- paste("beta.", years, sep = "")
  for (i in 1:length(years)) {
    allz        <- subset(allz.t, allz.t$year == years[i])
    which.col   <- paste("r_year[", years[i], ",", trait, "]", sep = "")
    theta.t[,i] <- (ps[,2] + ps[ , which.col]) / two.beta.zz
    Wmax.t[,i]  <- exp( ps[ , 1] - (ps[,2] + ps[ , which.col])^2 / (2 * two.beta.zz) )
    temp <- data.frame(matrix(NA, nrow = 4000, ncol = length(allz$year)))
    for (j in 1:length(allz$year)) {
      temp[,j] <- 1 - exp(-(allz[j,2] - theta.t[,i])^2 / (2 * omega^2))             # Fractional load for this specific brood, from Eqn 1 of Chevinet al. (2015, Evolution).
      }
    total.load.t[,i] <- apply(temp, 1, mean)
    lag.load.t[,i]   <- 1 - exp(-(meanz.t$x[i] - theta.t[,i])^2 / (2 * omega^2))    # Proportional lag load observed in year t, from Eqn 1 of Chevin et al. (2015, Evolution).
    beta.t[,i]       <- (theta.t[,i] - meanz.t$x[i]) / (omega^2 + 1)                # This from Lande (1976, Evolution). Not sure exactly where, but maybe about Eqn 14.
    }
  # Calculate the evolutionary load due to fluctuating theta, as per de Villemereuil et al. (2020, PNAS).
  # This is part of the lag load so it should be smaller than the lag load.
  sigma2.theta.t <- apply(theta.t, 1, var)
  fluct.load     <- sigma2.theta.t / (2 * (omega^2 + 1))
  # Remove outliers if you chose to do so.
  if (outliers == TRUE) {
    omega      <- remove.outliers(omega, units = threshold, remove.na = "na")
    theta      <- as.data.frame(apply(theta.t, 2, remove.outliers, units = threshold, remove.na = "na"))
    Wmax       <- as.data.frame(apply(Wmax.t, 2, remove.outliers, units = threshold, remove.na = "na"))
    total.load <- as.data.frame(apply(total.load.t, 2, remove.outliers, units = threshold, remove.na = "na"))
    lag.load   <- as.data.frame(apply(lag.load.t, 2, remove.outliers, units = threshold, remove.na = "na"))
    beta       <- as.data.frame(apply(beta.t, 2, remove.outliers, units = threshold, remove.na = "na"))
    fluct.load <- remove.outliers(fluct.load, units = threshold, remove.na = "na")
    }
  return(list(years = years, omega = omega, theta = theta, Wmax = Wmax, beta = beta, total.load = total.load, lag.load = lag.load, fluctuation.load = fluct.load))
  }    # end function Chevin.2015
#
#
draw.contour.labels <- function(xlimits, ylimits, xfraction, yfraction, x.corners, y.corners, contour.labels, label.bg = "white", text.size, draw.box = TRUE, label.line.width = 0.8) {
  # Dedicated internal function for drawing labels on contour intervals
  # Arguments:
  #  xlimits           vector with lower and upper limits of the x axis.
  #  ylimits           limits of the y axis.
  #  xfraction         width of the label as a fraction of the full x axis (usually 0.08-0.10).
  #  yfraction         height of the label as a fraction of the full y axis.
  #  x.corners         vector with x-axis coordinates of the lower left corner of each label.
  #  y.corners         y-axis coordinates of the lower left corner of each label.
  #  contour.labels    character vector of labels; must be same length as x.corners and y.corners.
  #  label.bg          Fill color of the contour labels. NA if you want no rectangle drawn.
  #  text.size         size of text to use.
  #  draw-box          TRUE of you want a line drawn around the label.
  #  label.line.width  thickness of the line around the box.
  xint <- xfraction * diff(xlimits)
  yint <- yfraction * diff(ylimits)
  for (i in 1:length(contour.labels)) {
    if (draw.box) { rect(x.corners[i], y.corners[i], x.corners[i]+xint, y.corners[i]+yint, col = label.bg, lwd = label.line.width)
      } else {
      rect(x.corners[i], y.corners[i], x.corners[i]+xint, y.corners[i]+yint, col = label.bg, border = NA)
      }
    text(contour.labels[i], x = x.corners[i]+(0.5*xint), y = y.corners[i]+(0.5*yint), adj = 0.5, cex = text.size)
    }
  }
#
#
samp.size = function(x) {
  # Calculate sample sizes (excluding NAs) of selected columns in a data frame.
  # This function is usually called by "data.frame.sample.size".
  n = length(x) - sum(is.na(x))
  nas = sum(is.na(x))
  out = c(n, nas)
  names(out) = c("", "NAs")
  out
  }
#
#
NA.to.zero <- function(x) {
  # Convert NA's to 0 (zero's).
  # Args:
  #  x = a number or vector of numbers
  x[is.na(x)] <- 0
  x
  }
#
#
rename <- function(dat, oldnames, newnames) {
     # Rename one or more columns, and return the entire data frame.
     # Examples of calling this function:
     #    rename(data, "old1", "new1")
     #    rename(data, c("old1", "old2"), c("new1", "new2"))
     #
  if( length(oldnames) != length(newnames) ) {
    message("Error: oldnames and newnames must have the same length")
    flush.console()
    break
    }
  for(i in 1:length(oldnames)) {
    colnames(dat)[which(names(dat) == oldnames[i])] <- newnames[i]
    }
  return(dat)
  }
#
#
format.p.val <- function(x, ndigits = 5) {
  # Format values to a specified number of digits.
  # x is a number or vector.
  ff <- format(x, scientific = FALSE)
  suppressWarnings(format(round((as.numeric(ff) * (10^ndigits)))/(10^ndigits), scientific = FALSE))
  }
#
#
filled.contour3 <- function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...)  {
  # Modification by Ian Taylor from the filled.contour function to remove the key and facilitate overplotting with contour()
  # Further modified by Carey McGilliard and Bridget Ferris to allow multiple plots on one page.
  # Came from here on Sept 2010: http://wiki.cbr.washington.edu/qerm/sites/qerm/images/b/bb/Example_4-panel_v1a.R
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
 # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
 # on.exit(par(par.orig))
 # w <- (3 + mar.orig[2]) * par("csi") * 2.54
 # par(las = las)
 # mar <- mar.orig
 plot.new()
 # par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
       col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}
#
#
critical.r <- function( n, alpha = 0.05 ) {
  # Calculate critical values of correlation coefficients.
  # I got this from: https://www.researchgate.net/post/What_is_the_formula_to_calculate_the_critical_value_of_correlation#:~:text=The critical value of the "Pearson Correlation Coefficient",t values are from the Critical t Table%29
  # Arguments:
  #   n      number of pairs
  #   alpha  critical alpha value
  # Call the function like this: critical.r(n = 102)
  df         <- n - 2
  critical.t <- qt( alpha/2, df, lower.tail = F )
  critical.r <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
  return( critical.r )
  }
#
#
remove.outliers <- function(x, med.mean = "med", units = 2, remove.na = "remove") {
   # Remove outliers that are greater than a certain number of SD units away from the mean or median.
   # Arguments:
   #  x         = a vector of observations from which outliers should be removed.
   #  med.mean  = "mean" or "med", indicating whether outliers should be measured from the median or mean.
   #  units     = threshold distancefrom the mean/med used to gauge outliers (SD units).
   #  remove.na = "remove" or "na", indicating whether outliers should be removed from x or replaced with NA.
   #
   std.dev     <- sd(x, na.rm = TRUE)
   threshold.x <- median(x, na.rm = TRUE)
   if(med.mean == "mean") threshold.x <- mean(x, na.rm = TRUE)
   scaled.dev <- abs(x - threshold.x) / std.dev
   if(remove.na == "remove") x <- x[ !(scaled.dev > units) ]    # units SD away from the med/mean is considered an outlier and removed.
   if(remove.na == "na") x[ (scaled.dev > units) ] <- NA        # units SD away from the med/mean is considered an outlier and replaced with NA.
   x
   }
#
#
aggregate.weighted.mean <- function(x, w, by, dat) {
  # Calculate weighted means of one or more variables by one or more groups.
  #
  # Requires packages data.table.
  # Arguments:
  #  x   = Names of one or more variables for which weighted means are desired.
  #  w   = A single variable name containing weights.
  #  by  = Names of one or more groups within which weighted means should be calculated.
  #  dat = data frame with all the variables.
  #
  require(data.table)
  n.var <- length(x)
  for (i in 1:n.var) {
    d        <- dat[,c(x[i], w, by)]
    names(d) <- c("x", "w", by)
    d        <- d[complete.cases(d),]
    setDT(d)
    d[, wt.mean.x := weighted.mean(x, w), by = by]
    d[, wt.mean.n := length(x), by = by]
    d   <- unique(d, by = c(by, "wt.mean.x", "wt.mean.n"))
    d   <- rename(d, c("wt.mean.x", "wt.mean.n"), c(paste("wtMean.", x[i], sep = ""), paste("N.", x[i], sep = "")))
    d$x <- d$w <- NULL
    if (i == 1) { res <- d
      } else {
      res <- merge(res, d)
      }
    }
  return(res)
  }
#
#
tryCatch.W.E <- function(expr) {
  # Run code, catch errors, and keep going. Works better than try().
  # I got this from demo(error.catching).
  # After running this, the following will be true if there was an error: "if(is.null(m.D.obs$value$message)==FALSE)"
  W <- NULL
  w.handler <- function(w) { # warning handler
    W <<- w
    invokeRestart("muffleWarning")
    }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler), warning = W)
  }
#
#
image.scale <- function(z, zlim, col = heat.colors(12), breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...) {
  if(!missing(breaks)) {         # This function adds a scale to an image plot.
                                 # downloaded from http://menugget.blogspot.de/2011/08/adding-scale-to-image-plot.html on 160202.
    if(length(breaks) != (length(col)+1)){stop("must have one more break than color")}
    }
  if(missing(breaks) & !missing(zlim)) {
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
    }
  if(missing(breaks) & missing(zlim)) {
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3) #adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
    }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)) { poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i]) }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz)  { YLIM <- c(0,1); XLIM <- range(breaks) }
  if(!horiz) { YLIM <- range(breaks); XLIM <- c(0,1) }
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)) {
    if(horiz) { polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA) }
    if(!horiz){ polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA) }
    }
  }
#


