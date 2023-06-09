### Functions used in 'main_TI_RSOS_2018.r' 
# Tsuyoshi Ito
# August 20th, 2018


### read csv files ------------------------------------------------------
# This is a slight modification of 'read.csv.folder' function in 'Morpho' package.
read.csv.folder2 <-
  function(folder,
           x = 1:28,
           y = 4:6,
           rownames = NULL,
           header = FALSE,
           dec = ".",
           sep = ",",
           pattern = "csv",
           addSpec = NULL,
           back = TRUE,
           skip = 1)
  {
    if (substr(folder, start = nchar(folder), stop = nchar(folder)) != "/")
    {
      folder <- paste(folder, "/", sep = "")
    }
    
    file.ext <- paste(".", pattern, sep = "")
    name <- list.files(folder, pattern = file.ext)
    xlen <- length(x)
    ylen <- length(y)
    NA.list <- NULL
    
    ln <- length(name)
    arr <- array(NA, dim = c(xlen, ylen, ln))
    if (is.factor(x))
    {
      x <- as.character(x)
    }
    if (is.character(x))
      # check if selection contains variable names
      for (i in 1:ln)
      {
        data <-
          read.table(
            paste(folder, name[i], sep = ""),
            header = header,
            dec = dec,
            sep = sep,
            skip = skip
          )
        dat <- NULL
        count <- 1
        if (is.null(rownames))
        {
          stop("please specify column containing Landmark names!")
        }
        rn <- data[, rownames]
        for (j in 1:length(x))
        {
          check <- which(rn == x[j])
          
          if (length(check) == 0)
          {
            warning(paste("dataset", i, "misses entry for Landmark", j))
            data[9999, y] <- rep(NA, ylen)
            dat[count] <- 9999
            
          }
          if (length(check) > 1)
          {
            warning(
              paste(
                "dataset",
                i,
                "contains landmark #",
                x[j],
                "with the same name - first match was used."
              )
            )
            dat[count] <- check[1]
          }
          else
          {
            empty <- which(rn == x[j])
            if (length(empty) != 0)
            {
              dat[count] <- which(rn == x[j])
            }
          }
          count <- count + 1
        }
        arr[, , i] <- as.matrix(data[dat, y])
        if (i == 1)
          rown <- x
        
      }
    
    else
    {
      for (i in 1:ln)
      {
        data <-
          read.table(
            paste(folder, name[i], sep = ""),
            header = header,
            dec = dec,
            sep = sep,
            skip = skip
          )
        arr[, , i] <- as.matrix(data[x, y])
        if (i == 1)
          if (is.null(rownames))
          {
            rown <- c(1:xlen)
          }
        else
        {
          rown <- data[x, rownames]
        }
      }
    }
    
    
    nas0 <-
      which(is.na(arr))	# check for NAs and store information about missing Landmark and individual
    nas1 <- as.integer(nas0 / (xlen * ylen)) + 1
    nas <- nas1[-(which(duplicated(nas1)))]
    
    if (length(nas) > 0)
    {
      NA.list <- list()
      for (i in 1:length(nas))
      {
        nas2 <- nas0[which(nas1 == nas[i])] %% (xlen * ylen)
        nas2 <- nas2 %% xlen
        nas2 <- nas2[-which(duplicated(nas2))]
        if (0 %in% nas2)
        {
          nas2[which(nas2 == 0)] <- xlen
        }
        if (length(nas2) > 0)
        {
          NA.list[[as.character(nas[i])]] <- sort(nas2)
        }
        else
        {
          NA.list[[as.character(nas[i])]] <- NULL
          nas <- nas[-i]
        }
      }
    }
    else
    {
      nas <- NULL
    }
    
    dim3 <- NULL
    if (back)
    {
      dim3 <- paste(sub(file.ext, "", name), addSpec, sep = "")
    }
    else
    {
      dim3 <- paste(addSpec, sub(file.ext, "", name), sep = "")
    }
    if (ylen == 2)
    {
      dimnames(arr) <- list(rown, c("X", "Y"), dim3)
    }
    else
      
    {
      dimnames(arr) <- list(rown, c("X", "Y", "Z"), dim3)
    }
    
    return(list(
      arr = arr,
      NAs = nas,
      NA.list = NA.list
    ))
    
  }





### Smirnovâ€“Grubbs test  -----------------------------------------------
# This is a slight modification of http://www.geocities.jp/ancientfishtree/NewFiles/sg.R.tar.gz
SG <- function(x, confidence.interval = 0.05, dell.value = NULL)
{
  while (TRUE)
  {
    branch.length <- range(x)
    x <- x[!is.na(x)]
    n <- length(x)
    t <- abs(range(x) - mean(x)) / sd(x)
    p <- n * pt(sqrt((n - 2) / ((n - 1) ^ 2 / t ^ 2 / n - 1)), n - 2, lower.tail =
                  FALSE)
    result <- list(parameter = c(df = n - 2))
    result1 <-
      structure(c(
        result,
        list(
          branch.length = branch.length[1],
          statistic = c(t = t[1]),
          p.value = p[1]
        )
      ), class = "htest")
    result2 <-
      structure(c(
        result,
        list(
          branch.length = branch.length[2],
          statistic = c(t = t[2]),
          p.value = p[2]
        )
      ), class = "htest")
    # Maximum
    if (result2$p.value < confidence.interval)
    {
      #          cat('res2 branch.length ',result2$branch.length,'\\n')
      #          cat('res2 p ',result2$p.value,'\\n')
      dell.value <- c(dell.value, result2$branch.length)
      #          cat('dell.value ',dell.value,'\\n')
      x <- x[x != result2$branch.length]
      next
    }
    # Minimum
    if (result1$p.value < confidence.interval)
    {
      #          cat('res1 branch.length ',result1$branch.length,'\\n')
      #          cat('res1 p ',result1$p.value,'\\n')
      dell.value <- c(dell.value, result1$branch.length)
      #          cat('dell.value ',dell.value,'\\n')
      x <- x[x != result1$branch.length]
      next
    }
    else{
      return(list(x = x, dellVal = dell.value))
      break
    }
  }
}




### matrix correlation between Procrustes adn PC Euclidian distances --------
library(shapes)
r.matrix <- function(arr = NULL, pcs = NULL) {
  n <- dim(arr)[3]
  symM <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i > j) {
        symM[i, j] <- procdist(arr[, , i], arr[, , j])
      } else if (i < j) {
        symM[i, j] <- procdist(arr[, , i], arr[, , j])
      }
    }
  }
  
  n <- nrow(pcs)
  l <- ncol(pcs)
  r <- NULL
  for (k in 1:l) {
    pcM <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i > j) {
          pcM[i, j] <- dist(rbind(pcs[i, 1:k],
                                  pcs[j, 1:k]))
        } else if (i < j) {
          pcM[i, j] <- dist(rbind(pcs[i, 1:k],
                                  pcs[j, 1:k]))
        }
      }
    }
    r <- c(r, cor(c(pcM), c(symM)))
  }
  return (r)
}




### rbind without column names ------------------------------------------
rbind2 <- function(..., cnames = NULL) {
  df <- list(...)
  x <- df[[1]]
  colnames(x) <- cnames
  for (i in 2:length(df)) {
    y <- df[[i]]
    colnames(y) <- cnames
    x <- rbind(x, y)
  }
  return (x)
}





### vector correlation --------------------------------------------------
# This is written with reference to https://doi.org/10.5061/dryad.r43k1.2
ang.vec <- function(x, y) {
  vec.angle <- NULL
  vec.cor <- NULL
  if (is.null(nrow(x))) {
    xi <- x
    yi <- y
    cpro <- (t(xi) %*% yi) / (sqrt(t(xi) %*% xi) * sqrt(t(yi) %*% yi))
    vec.angle <- acos(cpro) * (180 / pi)
    vec.cor <- cor(xi, yi)
  } else{
    for (i in 1:nrow(x)) {
      xi <- c(as.matrix(x[i, ]))
      yi <- c(as.matrix(y[i, ]))
      cpro <- (t(xi) %*% yi) / (sqrt(t(xi) %*% xi) * sqrt(t(yi) %*% yi))
      vec.angle <- c(vec.angle, acos(cpro) * (180 / pi))
      vec.cor <- c(vec.cor, cor(xi, yi))
    }
  }
  return(list(vec.cor = vec.cor, vec.angle = vec.angle))
}

ang.vec2 <- function(x, y) {
  x2 <- apply(x, 2, mean)
  y2 <- apply(y, 2, mean)
  vec.angle <- NULL
  vec.cor <- NULL
  for (i in 1:nrow(x)) {
    xi <- c(as.matrix(x[i, ]))
    cpro <- (t(xi) %*% y2) / (sqrt(t(xi) %*% xi) * sqrt(t(y2) %*% y2))
    vec.angle <- c(vec.angle, acos(cpro) * (180 / pi))
    vec.cor <- c(vec.cor, cor(xi, y2))
  }
  for (i in 1:nrow(y)) {
    yi <- c(as.matrix(y[i, ]))
    cpro <- (t(x2) %*% yi) / (sqrt(t(x2) %*% x2) * sqrt(t(yi) %*% yi))
    vec.angle <- c(vec.angle, acos(cpro) * (180 / pi))
    vec.cor <- c(vec.cor, cor(x2, yi))
  }
  return(list(vec.cor = vec.cor, vec.angle = vec.angle))
}




### PCA for geometric morphometric --------------------------------------
# This is written with reference to 'plotTangentSpace' function in 'geomorph' package
geom_pca <- function (A) {
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  
  d <- prcomp(x)$sdev ^ 2
  cd <- cumsum(d) / sum(d)
  cd <- length(which(cd < 1))
  if (length(cd) < length(d))
    cd <- cd + 1
  tol <- max(c(d[cd] / d[1], 0.005))
  
  pc.res <-
    prcomp(
      x,
      center = TRUE,
      scale = FALSE,
      retx = TRUE,
      tol = tol
    )
  pcdata <- pc.res$x
  
  shapes <- shape.names <- NULL
  for (i in 1:ncol(pcdata)) {
    pcaxis.min <-
      -3 * sd(pcdata[, i])
    pcaxis.max <- 3 * sd(pcdata[, i])
    pc.min <- pc.max <- rep(0, dim(pcdata)[2])
    pc.min[i] <- pcaxis.min
    pc.max[i] <- pcaxis.max
    pc.min <-
      as.matrix(pc.min %*% (t(pc.res$rotation))) + as.vector(t(ref))
    pc.max <-
      as.matrix(pc.max %*% (t(pc.res$rotation))) + as.vector(t(ref))
    shapes <- rbind(shapes, pc.min, pc.max)
    shape.names <-
      c(shape.names,
        paste("PC", i, "min", sep = ""),
        paste("PC", i, "max", sep = ""))
  }
  shapes <- arrayspecs(shapes, p, k)
  shapes <- lapply(seq(dim(shapes)[3]), function(x)
    shapes[, , x])
  names(shapes) <- shape.names
  
  out <-
    list(
      pc.summary = summary(pc.res),
      pc.scores = pcdata,
      pc.shapes = shapes,
      sdev = pc.res$sdev,
      rotation = pc.res$rotation
    )
  out
}




### shape score ---------------------------------------------------------
# This is written with reference to https://doi.org/10.5061/dryad.r43k1.2
shape_vec <- function(A, B) {
  library(geomorph)
  
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  
  d <- prcomp(x)$sdev ^ 2
  cd <- cumsum(d) / sum(d)
  cd <- length(which(cd < 1))
  if (length(cd) < length(d))
    cd <- cd + 1
  tol <- max(c(d[cd] / d[1], 0.005))
  
  pc.res <-
    prcomp(
      x,
      center = TRUE,
      scale = FALSE,
      retx = TRUE,
      tol = tol
    )
  pcdata <- pc.res$x
  
  Y <- pcdata[, 1:10]
  
  y_sd <- apply(Y, 2, sd)
  y_mean <- apply(Y, 2, mean)
  
  shape.score <- as.matrix(Y) %*% B %*% ((t(B) %*% B) ^ -0.5)
  
  # unscaling
  pc_vec_min <- c(-B * 10 * y_sd + y_mean , rep(0, 30))
  pc_vec_max <- c(B * 10 * y_sd + y_mean , rep(0, 30))
  
  pc.min <-
    as.matrix(pc_vec_min %*% (t(pc.res$rotation))) + as.vector(t(ref))
  pc.max <-
    as.matrix(pc_vec_max %*% (t(pc.res$rotation))) + as.vector(t(ref))
  
  shapes <- rbind(pc.min, pc.max)
  shapes <- arrayspecs(shapes, 28, 3)
  
  out <- list(shapes = shapes, shape.score = shape.score)
  out
}




### summary  ------------------------------------------------------------
# This is written with reference to https://github.com/marschmi/project02/blob/master/RCookbook_summarySE.R
summaryCI <-
  function(data = NULL,
           measurevar,
           groupvars = NULL,
           na.rm = FALSE,
           conf.interval = .95,
           .drop = TRUE) {
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm = FALSE) {
      if (na.rm)
        sum(!is.na(x))
      else
        length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(
      data,
      groupvars,
      .drop = .drop,
      .fun = function(xx, col) {
        c(
          N    = length2(xx[[col]], na.rm = na.rm),
          mean = mean   (xx[[col]], na.rm = na.rm),
          sd   = sd     (xx[[col]], na.rm = na.rm),
          ci_upper = quantile (xx[[col]], probs = conf.interval /
                                 2 + .5, na.rm = na.rm)[[1]],
          ci_lower = quantile (
            xx[[col]],
            probs = (1 - conf.interval) / 2,
            na.rm = na.rm
          )[[1]]
        )
      },
      measurevar
    )
    
    datac$se <-
      datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    return(datac)
  }
