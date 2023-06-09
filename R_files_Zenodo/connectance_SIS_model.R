# Timestamp: 2020/03/30.

require(gridExtra) # For grid.arrange().
library(epiR) # For epi.prev().
library(extraDistr) # For rtpois().
library(hexbin) # For hexbinplot().
library(latticeExtra) # For panel.lines().
library(pagenum) # For pagenum().
library(plotwidgets) # For hsl2col().
library(rootSolve) # For uniroot().
library(vecsets) # For vsetdiff().
library(viridis) # For viridis_pal().









# ---------------
# Initialisation.
# ---------------
num.bee.species <- 13
num.flower.species <- 7

# Typical scales of parameters.
frac.day.bee.spend.foraging <- 1/12
single.bee.visitation.rate.in.per.hr <- 5 * 60
prob.infected.bee.poops.on.flower <- 0.01
prob.contaminated.flower.infects.bee <- 0.01
bee.recovery.rate.in.per.hr <- 1/(3 * 24)
flower.recovery.rate.in.per.hr <- 1/3
num.of.bees.per.species <- 50
num.of.flowers.per.species <- 2000

alpha.scale <- (single.bee.visitation.rate.in.per.hr
                * frac.day.bee.spend.foraging
                * prob.infected.bee.poops.on.flower
                * num.of.bees.per.species
                / num.of.flowers.per.species)

beta.scale  <- (single.bee.visitation.rate.in.per.hr
                * frac.day.bee.spend.foraging
                * prob.contaminated.flower.infects.bee)

gamma.scale <- bee.recovery.rate.in.per.hr

zeta.scale  <- flower.recovery.rate.in.per.hr

# How many order of magnitudes do we allow the parameters to vary?
ord.mag.range <- 1

# Number of networks to generate.
num.replicate <- 100000

# Split visits equally among all the flower species visited by the bee species?
equal.split <- F

# We model the bee degree distribution using a zero-truncated binomial.
# This is the range of proportions we will use for the binomial.
# We choose this range because the fitted proportions from the data varies between 0.13 and 0.31 for different sites.
prop.min <- 0.05
prop.max <- 0.4

# Define a function to perform the Gale-Ryser test.
# This is based on the version in Krause (1996).
# All degree vectors must already be sorted in descending order.
gale.ryser <- function(p, q){
  pstar <- apply( matrix(rep(as.vector(rbind(rep(1, length(p)), rep(0, length(p)))), 
                             as.vector(rbind(p, max(p) - p))),
                             nrow=length(p), byrow=T), 2, sum )
  
  # The initial degree distribution is designed such that length(pstar) = max(p) <= length(q).
  # However, the above might not be true for the remaining bee and flower degree distribution should we select the wrong flower species for a bee species.
  # So we cannot assume that the maximum integer to try is the length of q.
  m <- max(length(pstar), length(q))
  return ( all( cumsum( c( q, rep(0, m-length(q)) ) ) <= cumsum( c( pstar, rep(0, m-length(pstar)) ) ) ) )
}

# --------------------------------------------------------------------------------

# How R_0, maximum growth rate lambda and steady-state prevalence vary with connectance.

matrix.dim <- num.bee.species + num.flower.species

# Create dummy objects to store results.
result.df <- data.frame(num.realised.pair = rep(NA, num.replicate),
                        R0                = rep(NA, num.replicate),
                        lambda            = rep(NA, num.replicate),
                        mean.prevalence   = rep(NA, num.replicate))
ss.bee.prevalence.matrix <- c()
ss.flower.prevalence.matrix <- c()

set.seed(cat(Sys.time(),"\\n"))

start.time <- Sys.time()
for(replicate in 1:num.replicate){
  
  prop <- runif(1, prop.min, prop.max)
  
  flower.prop <- uniroot.all( function(x) {num.bee.species*x/(1-(1-x)^num.bee.species) - num.bee.species*num.flower.species*prop/(1 - (1-prop)^num.flower.species)/num.flower.species}, c(0,1) )

  # Create empty T and Sigma matrices to start.
  T.matrix <- matrix(0, nrow=matrix.dim, ncol=matrix.dim)
  Sigma.matrix <- matrix(0, nrow=matrix.dim, ncol=matrix.dim)
  
  # Randomise the bee and flower recovery rates. Insert them into the Sigma matrix.
  gamma.vec  <- 10^( runif( num.bee.species,
                            log10(gamma.scale) - ord.mag.range/2,
                            log10(gamma.scale) + ord.mag.range/2 ) )
  zeta.vec   <- 10^( runif( num.flower.species,
                            log10(zeta.scale) - ord.mag.range/2,
                            log10(zeta.scale) + ord.mag.range/2 ) )
  Sigma.matrix[1:num.bee.species, 1:num.bee.species] <- -diag(gamma.vec)
  Sigma.matrix[(num.bee.species+1):(num.bee.species+num.flower.species), (num.bee.species+1):(num.bee.species+num.flower.species)] <- -diag(zeta.vec)
  
  # Create vector for the random number of flower species each bee species visits.
  # We require this to be sorted in descending order for the Gale-Ryser test later on.
  bee.deg.vec <- sort( rtbinom(num.bee.species, num.flower.species, prop, a=0), decreasing=T )
  
  # Create vector for the random number of bee species visiting each flower species.
  
  flower.deg.vec <- rep(0, num.flower.species)
  gale.ryser.result <- F
  while(!gale.ryser.result){
    # Create flower.deg.vec with same total degrees as bee.deg.vec.
    while(sum(flower.deg.vec)!=sum(bee.deg.vec)){
      flower.deg.vec <- rtbinom(num.flower.species, num.bee.species, flower.prop, a=0)
    }
    # We require this to be sorted in descending order for the Gale-Ryser test later on.
    flower.deg.vec <- sort( flower.deg.vec, decreasing=T )
    # Perform Gale-Ryser test.
    gale.ryser.result <- gale.ryser(bee.deg.vec, flower.deg.vec)
  }

  # Flower degree pool.
  flower.deg.pool <- rep(1:num.flower.species, flower.deg.vec)

  # Loop over bee species.
  for(bee in 1:num.bee.species){
    
    # The specific flower species the bee species visits.
    # We want to make sure the remaining degree distributions still satisfy the Gale-Ryser test.
    gale.ryser.result <- F
    while(!gale.ryser.result){
      flower.visited.vec <- sample(x = unique(flower.deg.pool), size = bee.deg.vec[bee], replace = F)
      remaining.flower.deg.pool <- vsetdiff(flower.deg.pool, flower.visited.vec)
      # There is no need to perform the Gale-Ryser test for the last bee.
      if(bee < num.bee.species){
        gale.ryser.result <- gale.ryser(bee.deg.vec[(bee+1):num.bee.species],
                                        sort(as.vector(table(remaining.flower.deg.pool)), decreasing=T))
      } else {
        gale.ryser.result <- T
      }
    }
    flower.deg.pool <- remaining.flower.deg.pool
    
    # Splitting visits equally or randomly among flower species.
    if(equal.split){
      visit.fraction.vec <- rep(1/bee.deg.vec[bee], bee.deg.vec[bee])
    } else{
      visit.fraction.vec <- runif(bee.deg.vec[bee])
      visit.fraction.vec <- visit.fraction.vec/sum(visit.fraction.vec)
    }
    
    # Transmission parameters.
    alpha.vec <- 10^( runif( bee.deg.vec[bee],
                             log10(alpha.scale) - ord.mag.range/2,
                             log10(alpha.scale) + ord.mag.range/2 ) ) * visit.fraction.vec
    beta.vec  <- 10^( runif( bee.deg.vec[bee],
                             log10(beta.scale) - ord.mag.range/2,
                             log10(beta.scale) + ord.mag.range/2 ) ) * visit.fraction.vec
    T.matrix[bee, num.bee.species + flower.visited.vec] <- beta.vec
    T.matrix[num.bee.species + flower.visited.vec, bee] <- alpha.vec
  }
  
  # Calculating R0.
  # Note that solve() gives the inverse matrix.
  next.gen.matrix <- -T.matrix %*% solve(Sigma.matrix)
  second.gen.matrix <- next.gen.matrix %*% next.gen.matrix
  second.gen.matrix.flower <- second.gen.matrix[(num.bee.species+1):(num.bee.species+num.flower.species),
                                                (num.bee.species+1):(num.bee.species+num.flower.species)]
  # I verified that eigen() can return complex eigenvalues.
  second.gen.eigenval.vec <- eigen(second.gen.matrix.flower, symmetric=F, only.values=T)$values
  real.second.gen.eigenval.vec <- Re(second.gen.eigenval.vec[Im(second.gen.eigenval.vec)==0])
  R0 <- max(real.second.gen.eigenval.vec)
  
  # Calculating max rate.
  time.evol.matrix <- T.matrix + Sigma.matrix
  time.evol.eigenval.vec <- eigen(time.evol.matrix, symmetric=F, only.values=T)$values
  real.time.evol.eigenval.vec <- Re(time.evol.eigenval.vec[Im(time.evol.eigenval.vec)==0])
  lambda <- max(real.time.evol.eigenval.vec)
  
  # Calculating steady-state prevalence.
  steady.state.prevalence <- multiroot( function(x) {as.vector(T.matrix %*% matrix(x, nrow=matrix.dim)) * (1-x) + as.vector(Sigma.matrix %*% matrix(x, nrow=matrix.dim))},
                                        rep(1, matrix.dim), positive = T, maxiter = 1000)$root
  
  ss.bee.prevalence.matrix    <- rbind(ss.bee.prevalence.matrix,
                                       c(sum(bee.deg.vec), sort(steady.state.prevalence[1:num.bee.species])))
  ss.flower.prevalence.matrix <- rbind(ss.flower.prevalence.matrix,
                                       c(sum(bee.deg.vec), sort(steady.state.prevalence[(num.bee.species+1):(num.bee.species+num.flower.species)])))
  result.df[replicate, ] <- c(sum(bee.deg.vec), R0, lambda, mean(steady.state.prevalence[1:num.bee.species]))
}
print(Sys.time() - start.time)

# Restrict to range of connectance whose bins have sufficient counts in them (30).
bin.count.vec <- table(result.df$num.realised.pair/num.bee.species/num.flower.species)
x.min <- max(min(as.numeric(names(which(bin.count.vec>=30)))), 0.15)
x.max <- min(max(as.numeric(names(which(bin.count.vec>=30)))), 0.5)
x.lim <- c(x.min, x.max)
result.df <- result.df[result.df$num.realised.pair/num.bee.species/num.flower.species >= x.min &
                       result.df$num.realised.pair/num.bee.species/num.flower.species <= x.max, ]
ss.bee.prevalence.matrix <- ss.bee.prevalence.matrix[ss.bee.prevalence.matrix[, 1]/num.bee.species/num.flower.species >= x.min &
                                                     ss.bee.prevalence.matrix[, 1]/num.bee.species/num.flower.species <= x.max  , ]
ss.flower.prevalence.matrix <- ss.flower.prevalence.matrix[ss.flower.prevalence.matrix[, 1]/num.flower.species/num.flower.species >= x.min &
                                                           ss.flower.prevalence.matrix[, 1]/num.flower.species/num.flower.species <= x.max  , ]


# Generate plots.
median.R0.vec <- tapply(log10(result.df$R0), result.df$num.realised.pair/num.bee.species/num.flower.species, median)
plot1 <- hexbinplot(log10(R0) ~ num.realised.pair/num.bee.species/num.flower.species, data = result.df,
                    xlim=range(result.df$num.realised.pair/num.bee.species/num.flower.species),
                    ylim=c(-0.25, 1.25),
                    xbins=28, aspect="1",
                    xlab=list("Connectance", cex=1.5),
                    ylab=list(expression("Log reproduction number"~(log[10]*R[0])), cex=1.5),
                    scales=list(cex=1.5),
                    colorcut = seq(0,1,length.out=17),
                    # colramp = function(n) hsl2col(rbind(rep(179.4495413,17), rep(0.6228571,17), seq(0.9,0.1,length.out=17))),
                    colramp = function(n) colorRampPalette(c("white", "#218E8DFF"))(17),
                    panel=function(x, y, ...){
                      panel.hexbinplot(x,y,...)
                      panel.lines(loess.smooth(as.numeric(names(median.R0.vec)),
                                               median.R0.vec), col="black", lwd=3)
                    })


median.lambda.vec <- tapply(result.df$lambda, result.df$num.realised.pair/num.bee.species/num.flower.species, median)
plot2 <- hexbinplot(lambda ~ num.realised.pair/num.bee.species/num.flower.species, data = result.df,
           xlim=range(result.df$num.realised.pair/num.bee.species/num.flower.species),
           ylim=c(-0.01, 0.075),
           xbins=28, aspect="1",
           xlab=list("Connectance", cex=1.5),
           ylab=list(expression("Linearized growth rate"~(lambda)), cex=1.5),
           scales=list(cex=1.5),
           colorcut = seq(0,1,length.out=17),
           # colramp = function(n) hsl2col(rbind(rep(179.4495413,17), rep(0.6228571,17), seq(0.9,0.1,length.out=17))),
           colramp = function(n) colorRampPalette(c("white", "#218E8DFF"))(17),
           panel=function(x, y, ...){
             panel.hexbinplot(x,y,...)
             panel.lines(loess.smooth(as.numeric(names(median.lambda.vec)),
                                      median.lambda.vec), col="black", lwd=3)
           })

# Laura's colour.
col.line.vec <- rev(viridis_pal()(num.bee.species))
col.vec <- hsl2col(col2hsl(col.line.vec) + matrix(c(rep(0, 2*num.bee.species), rep(0.3, num.bee.species)), nrow=3, byrow=T))
  
plot3 <- xyplot(NA ~ NA,
                xlim=range(ss.bee.prevalence.matrix[,1]/num.bee.species/num.flower.species),
                ylim=range(ss.bee.prevalence.matrix[,2:(num.bee.species+1)]),
                aspect="1",
                xlab=list("Connectance", cex=1.5),
                ylab=list("Steady-state prevalence", cex=1.5),
                scales=list(cex=1.5),
                pch=20, col=col.vec, cex=0.1,
                panel=function(...)
                {
                  panel.xyplot(...)
                  for(i in 1:5000){
                    panel.points(rep(ss.bee.prevalence.matrix[i,1]+runif(1,-0.5,0.5), num.bee.species)/num.bee.species/num.flower.species,
                                 ss.bee.prevalence.matrix[i,-1],
                                 pch=20, col=col.vec, cex=0.1)
                  }
                  for(bee in 1:num.bee.species){
                    median.prevalence.vec <- tapply(ss.bee.prevalence.matrix[,bee+1], ss.bee.prevalence.matrix[,1]/num.bee.species/num.flower.species, median)
                    panel.lines(loess.smooth(as.numeric(names(median.prevalence.vec)), median.prevalence.vec),
                                col=col.line.vec[bee], lwd=3)
                  }
                })

# Height is decreased until it is now the the limiting dimension of both plots.
# As a result, both plot squares now have the same dimension.
dev.new(width=12, height=5.2, noRStudioGD = T)
grid.arrange(plot1, plot2, ncol=2)

pagenum(num="", text="(a)", x=.03, y=0.97, just=c("center", "center"), col="black", cex=1.5)
pagenum(num="", text="(b)", x=.53, y=0.97, just=c("center", "center"), col="black", cex=1.5)

dev.new(width=6, height=5.2, noRStudioGD = T)
plot3