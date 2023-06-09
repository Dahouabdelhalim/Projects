#######################################################################
# Creating functions and define parameters assumption
#######################################################################

alpha <- c(0.9718)
theta <- c(0.1)
sigmap <<- c(10, 1.000001, 0.1)

flexshares <- read.csv(file="inputs/flexshares.csv", header=TRUE, sep=",")

numberOfDays <- function(date) {
  m <- format(date, format="%m")
  while (format(date, format="%m") == m) {
    date <- date + 1
  }
  return(as.integer(format(date - 1, format="%d")))
}

set.scenario <- function(scen) {
  ifelse(scen==1, {#optimistic
    shareflex <<- c(0.67, 0.05, 0.28) #sigma level - highly flex, midflex, inflex
    shareother <<- c(0.15, 0.05, 0.8)},
    ifelse(scen==2, 
           {#2nd scenario=baseline/moderate
             shareflex <<- c(0.33, 0.05, 0.62)
             shareother <<- c(0.075, 0.05, 0.875)},
           ifelse(scen==3,
                  {#3rd scenario=pessimistic
                    shareflex <<- c(0.15, 0.05, 0.8)
                    shareother <<- c(0, 0.05, 0.95)},
                  print("Scenario is not available"))))
}

betaflexf <- function(xb, sigmath, month, scen) { #sigmath=1:highlyflex, 2:midflex, 3:inflex
  set.scenario(scen)
  flexload <- flexshares[,month]*shareflex[sigmath]
  otherload <- (1-flexshares[,month])*shareother[sigmath]
  loadshare <- xb*(flexload+otherload)
  betas <- loadshare/sum(loadshare)
  return(betas)
}
print("loaded beta function")


demandbygroupf <- function(xb, pd, pb, sigmath, month, scen) { 
  #pb = price at baseline; pd=new price; Mb=budget at baseline; scen=which scenario
  #sigmath:for each group highlyflex=1; midflex=2; inflex=3
  #calculate load shares that are flexible
  set.scenario(scen)
  flexload <- flexshares[,month]*shareflex[sigmath]
  otherload <- (1-flexshares[,month])*shareother[sigmath]
  loadshare<- xb*(flexload+otherload)
  #set parameters
  p             <- pd/pb
  Mb            <- sum(pb*loadshare)/(1-alpha)
  beta          <- betaflexf(xb, sigmath, month, scen)
  sigma         <- sigmap[sigmath]
  #expressions to construct function
  exp1 <- (1-theta)/(1-sigma)
  exp2 <- (sigma-theta)/(1-sigma)
  exp3 <- sum(beta*(p^(1-sigma)))
  e <- ((alpha+(1-alpha)*(exp3^(exp1)))^(-1))
  x <- (Mb* e *((1-alpha)*(exp3)^exp2)*beta*(p^(-1*sigma)))
  return(x/pb)
}

demandbygroupf.ces <- function(xb, pd, pb, sigmath, month, scen) { 
  #pb = price at baseline; pd=new price; Mb=budget at baseline; scen=which scenario
  #sigmath:for each group highlyflex=1; midflex=2; inflex=3
  #calculate load shares that are flexible
  set.scenario(scen)
  flexload <- flexshares[,month]*shareflex[sigmath]
  otherload <- (1-flexshares[,month])*shareother[sigmath]
  loadshare<- xb*(flexload+otherload)
  #set parameters
  p             <- pd/pb
  Mb            <- sum(pb*loadshare)
  beta          <- betaflexf(xb, sigmath, month, scen)
  sigma         <- sigmap[sigmath]
  #expressions to construct function
  e <- (sum(beta*(p^(1-sigma))))^(1/(1-sigma))
  v <- Mb/e
  x <- loadshare*(v/Mb)*((e/p)^sigma)
  return(x)
}

demandtotf <- function(xb, pd, pb, month, scen){
  tot <- 0
  for (i in (1:3)) {
    ifelse(i==1, tot <- demandbygroupf(xb, pd, pb, i,month,scen),
           tot <- tot + demandbygroupf(xb, pd, pb, i,month,scen))
  }
  return(tot)
}
print("loaded demand function")

#willingness to pay function total
wtpbygroupf <- function(xd, xb, pd, pb, sigmath, month, scen) { 
  #pb = price at baseline; pd=new price; Mb=budget at baseline; scen=which scenario
  #sigmath:for each group highlyflex=1; midflex=2; inflex=3
  #calculate load shares that are flexible
  set.scenario(scen)
  flexload <- flexshares[,month]*shareflex[sigmath]
  otherload <- (1-flexshares[,month])*shareother[sigmath]
  loadshare<- xb*(flexload+otherload)
  #set parameters
  p             <- pd/pb
  Mb            <- sum(pb*loadshare)/(1-alpha)
  beta          <- betaflexf(xb, sigmath, month, scen)
  sigma         <- sigmap[sigmath]
  #expressions to construct function
  exp1 <- (1-theta)/(1-sigma)
  exp2 <- 1/(1-theta)
  exp3 <- sum(beta*(p^(1-sigma)))
  e <- ((alpha+(1-alpha)*(exp3^(exp1)))^exp2)
  cs  <- Mb/e
  cost <- sum(pd*xd)
  #print(paste("beta:",beta,"e:",e,"cs:",cs,"cost:",cost))
  return(cs+cost)
}

csbygroupf <- function(xb, pd, pb, sigmath, month, scen) { 
  #pb = price at baseline; pd=new price; Mb=budget at baseline; scen=which scenario
  #sigmath:for each group highlyflex=1; midflex=2; inflex=3
  #calculate load shares that are flexible
  set.scenario(scen)
  flexload <- flexshares[,month]*shareflex[sigmath]
  otherload <- (1-flexshares[,month])*shareother[sigmath]
  loadshare<- xb*(flexload+otherload)
  #set parameters
  p             <- pd/pb
  Mb            <- sum(pb*loadshare)/(1-alpha)
  beta          <- betaflexf(xb, sigmath, month, scen)
  sigma         <- sigmap[sigmath]
  #expressions to construct function
  exp1 <- (1-theta)/(1-sigma)
  exp2 <- 1/(1-theta)
  exp3 <- sum(beta*(p^(1-sigma)))
  e <- ((alpha+(1-alpha)*(exp3^(exp1)))^exp2)
  cs  <- Mb/e
  return(cs)
}

csbygroupf.ces <- function(xb, pd, pb, sigmath, month, scen) { 
  #pb = price at baseline; pd=new price; Mb=budget at baseline; scen=which scenario
  #sigmath:for each group highlyflex=1; midflex=2; inflex=3
  #calculate load shares that are flexible
  set.scenario(scen)
  flexload <- flexshares[,month]*shareflex[sigmath]
  otherload <- (1-flexshares[,month])*shareother[sigmath]
  loadshare<- xb*(flexload+otherload)
  #set parameters
  p             <- pd/pb
  Mb            <- sum(pb*loadshare)/(1-alpha)
  beta          <- betaflexf(xb, sigmath, month, scen)
  sigma         <- sigmap[sigmath]
  #expressions to construct function
  e <- (sum(beta*(p^(1-sigma))))^(1/(1-sigma))
  v <- Mb/e
  return(v)
}


wtpbygroupf.ces <- function(xd, xb, pd, pb, sigmath, month, scen) { 
  #pb = price at baseline; pd=new price; Mb=budget at baseline; scen=which scenario
  #sigmath:for each group highlyflex=1; midflex=2; inflex=3
  #calculate load shares that are flexible
  set.scenario(scen)
  flexload <- flexshares[,month]*shareflex[sigmath]
  otherload <- (1-flexshares[,month])*shareother[sigmath]
  loadshare<- xb*(flexload+otherload)
  #set parameters
  p             <- pd/pb
  Mb            <- sum(pb*loadshare)/(1-alpha)
  beta          <- betaflexf(xb, sigmath, month, scen)
  sigma         <- sigmap[sigmath]
  #expressions to construct function
  e <- (sum(beta*(p^(1-sigma))))^(1/(1-sigma))
  v <- Mb/e
  cost <- sum(pd*xd)
  return(v+cost)
}


baseloadbygf <- function(xb, sigmath, month, scen) { 
  #pb = price at baseline; pd=new price; Mb=budget at baseline; scen=which scenario
  #sigmath:for each group highlyflex=1; midflex=2; inflex=3
  #calculate load shares that are flexible
  set.scenario(scen)
  flexload <- flexshares[,month]*shareflex[sigmath]
  otherload <- (1-flexshares[,month])*shareother[sigmath]
  loadshare<- xb*(flexload+otherload)
  return(loadshare)
}


integrand <- function(t, p0, p1) {
  p <- t * p1 + (1-t) * p0
  x <- x(p)
  return (sum(x * (p1-p0)))
}

wtp <- function(p) {
  p0 = 180+0*p # arbitrary baseline prices
  p1 = p
  wtp = -integrate(Vectorize(integrand, vectorize.args=c('t')), 0, 1, p0, p1)$value
  return (wtp)
}

q.dot.p <- function(q, p) {
  # calculate dot products for each of the components of the q and p lists
  # (energy, up reserves and down reserves), then add them together. 
  # This acts like one dot product across all energy and reserve products.
  return (Reduce("+", lapply(lmerge("*", q, p), sum)))
}

lmerge <- function (func, a, b) {
  # combine elements of lists a and b, using function func
  # works like mapply, but retains list structure
  list_len <- length(a)
  if (list_len != length(b))
    stop("merge.lists was called with two lists of different lengths")
  the_func <- match.fun(func)
  out <- vector("list", list_len)
  for (i in 1:list_len)
    out[[i]] <- the_func(a[[i]], b[[i]])
  return (out)
}

q.dot.p.for.t <- function(t, p0, p1, b.loads, b.prices, sigmath, month, scenario) {
  # calculate demand * (p1-p0) at a price equal to t*p1 + (1-t)*p0
  # This can be integrated for t running from 0 to 1 to calculate
  # the line integral of q.dp from p0 to p1.

    # make a list of prices for energy and reserve products
  p <- lmerge("+", lapply(p1, "*", t), lapply(p0, "*", 1-t))
  # get demand at these prices
  q <- demandbytype(b.loads, b.prices, p[[1]], p[[2]], p[[3]], sigmath, month, scenario)
  #print(paste(sum(q[[1]]),"sigma=",sigmath))
  # calculate p1-p0 for the line integral
  p_diff <- lmerge("-", p1, p0)
  # calculate dot product of q and p_diff
  return (q.dot.p(q, p_diff))
}

demandbytype <- function(b.loads, b.prices, prices, pup, pdown, sigmath, month, scenario) {
  # calculate the highest and lowest loads that could occur in each hour, 
  # holding prices constant in other hours
  demandflexmin <- vector('numeric', length(prices))
  demandflexmax <- vector('numeric', length(prices))
  testprices <- prices
  testprices[testprices < 1] <- 1
  for (i in (1:length(prices))) {
    # could find the "realistic" lower limit to demandflex as follows, but probably not necessary:
    # testprices[i] <- max(100 * prices[i], 100 * b.prices[i]) # must be large compared to current or baseline prices
    # test <- demandfflex(b[,timeseries,location], testprices, b.prices, Mb, month, scenario)/b.prices
    # demandflexmin[i] = test[i]
    demandflexmin[i] <- 0
    old_price <- testprices[i]
    testprices[i] <- 1
    test <- demandbygroupf(b.loads, testprices, b.prices, 1, month, scenario)
    demandflexmax[i] = test[i]
    testprices[i] <- old_price
  }
  # adjust prices and budget to reflect net price of electricity consumption
  # note: qup = -(qflex-qflexmin), qdown = -(qflexmax-qflex), 
  # so net cost to customer will be
  # qother*prices + qflex*prices + qup*pup + qdown*pdown
  # = qother*prices + qflex*prices - (qflex-qflexmin)*pup - (qflexmax-qflex)*pdown
  # = qother*prices + qflex * (prices - pup + pdown) + qflexmin*pup - qflexmax*pdown
  # so we adjust budget and prices accordingly.
  netprices <- (prices - pup + pdown)
  
  # set floors on prices
  netprices[netprices < 1] <- 1
  prices[prices < 1] <- 1
  
  # find demand levels
  ifelse(sigmath==0, {
  demandflex = demandbygroupf(b.loads, netprices, b.prices, 1, month, scenario)
  demandother = (demandbygroupf(b.loads, prices, b.prices, 2, month, scenario) + 
                   demandbygroupf(b.loads, prices, b.prices, 3, month, scenario))
  demand = demandflex + demandother
  demandup = -(demandflex - demandflexmin)
  demanddown = -(demandflexmax - demandflex)
  }, ifelse(sigmath==1, {
    demand = demandbygroupf(b.loads, netprices, b.prices, 1, month, scenario)
    demandup = -(demand - demandflexmin)
    demanddown = -(demandflexmax - demand)
  }, ifelse(sigmath==10, {#flexible demand without considering
    demand = demandbygroupf(b.loads, prices, b.prices, 1, month, scenario)
    demandup = 0
    demanddown = 0
    }, ifelse(sigmath==2, {
      demand = demandbygroupf(b.loads, prices, b.prices, 2, month, scenario)
      demandup = 0
      demanddown = 0}, {
        demand = demandbygroupf(b.loads, prices, b.prices, 3, month, scenario)
        demandup = 0
        demanddown = 0}
  ))))
  # note: reserve quantities are negative to indicate sales from the demand side
  return (list(demand, demandup, demanddown))

}

integralcs <- function(p1, p0, q0, sigmath, month, scenario) {
  # calculate the line integral of q dot p from p0 to p1, with the specified demand parameters
  iqp <- integrate(
    Vectorize(q.dot.p.for.t, vectorize.args=c('t')), 0, 1,
    p0, p1, q0[[1]], p0[[1]], sigmath, month, scenario)$value
  # get the quantities bid for energy and up and down reserves (assuming p0 and q0 are baselines)
  q1 <- demandbytype(q0[[1]], p0[[1]], p1[[1]], p1[[2]], p1[[3]], sigmath, month, scenario)
  # calculate wtp for this bid
  wtp <- q.dot.p(q1, p1) - q.dot.p(q0, p0) - iqp
  cs <- wtp - q.dot.p(q1, p1) #CS contains q.dot.p(q0, p0)
  return(list(wtp, cs, q.dot.p(q1, p1), q1))
}

