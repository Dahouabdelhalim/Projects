#required packages ####

require(lhs)
require(asnipe)
require(compiler)
require(truncnorm)

# convenience functions ####

#convenience function to calculate p-values
p_val <- function(observed, random){

  full <- c(observed,random) #full distribution, including observed value
  ls <- mean(full <= observed) #probability of less than or equal to observed
  gr <- mean(full >= observed) #probability of greater than or equal to observed
  p.twotailed <- min(c(ls,gr))*2 #minimum times two to get two-tailed
  return(p.twotailed)

}

#suppress printed output of functions
quiet <- function(x){
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# data generation functions ####

#function to generate GBI with variations in gregariousness
generate_gbi <- cmpfun(function(n.ids, n_groups, maxim.size, sd_pref){

  mean_gs <- mean(1:maxim.size)
  val <- rtruncnorm(n = n.ids, mean = mean_gs, sd = sd_pref, a = 0, b = maxim.size)
  group_size <- sample(maxim.size,n_groups,rep=T) #sizes of each group

  gbi <- matrix(0, ncol = n.ids, nrow = n_groups) #set up the gbi

  for(i in 1:n_groups){
    p_group <- 1/((group_size[i] - val)^2) #get probability of group membership
    in_group <- sample(n.ids, group_size[i], prob = p_group) #sample
    gbi[i,in_group] <- 1 #record which individuals were in group
  }

  return(gbi)

})

#function to generate SP with variations in association probability
generate_sp <- cmpfun(function(n.ids, assoc.shape1, assoc.shape2, n_samp, sight.prob){
  assoc_probs <- matrix(0,nrow = n.ids, ncol = n.ids)
  assoc_probs[lower.tri(assoc_probs)] <- rbeta(n = sum(lower.tri(assoc_probs)), shape1 = assoc.shape1, shape2 = assoc.shape2)
  assoc_probs[upper.tri(assoc_probs)] <- t(assoc_probs)[upper.tri(assoc_probs)]
  sp <- array(0, dim = c(n_samp,n.ids,n.ids))
  sight <- matrix(0, nrow = n_samp, ncol = n.ids)
  for(i in 1:n_samp){
    sight[i,] <- rbinom(size = 1, n = n.ids, prob = sight.prob)
    sp.i <- sp[i,,]
    sp.i[lower.tri(sp.i)] <- rbinom(n = sum(lower.tri(sp.i)), size = 1, prob = assoc_probs[lower.tri(assoc_probs)])
    sp.i[upper.tri(sp.i)] <- t(sp.i)[upper.tri(sp.i)]
    sp.i[sight[i,]==0,] <- 0
    sp.i[,sight[i,]==0] <- 0
    sp[i,,] <- sp.i
  }
  return(list(sp,sight))
})

# permutation functions ####

#implementations of the "trial swap" method proposed by Miklos & Podani (2004), suggested for SNA by Kraus et al. (2009), and implemented in SOCPROG
trial_swap_sp <- cmpfun(function(sp, trials){
  n <- dim(sp)[2] #number of individuals
  t <- dim(sp)[1] #number of sampling periods
  for(i in 1:trials){
    s <- sample(t,1) # choose a random sampling period
    inds <- sample(n,4) #choose 4 unique individuals
    c <- inds[1:2] #columns
    r <- inds[3:4] #rows
    trial.matrix <- sp[s,c,r] #trial swap submatrix
    if(all(rowSums(trial.matrix) == 1) & all(colSums(trial.matrix) == 1)){ #check if the swap is possible
      sp[s,c,r] <- ifelse(sp[s,c,r] == 1, 0, 1) #if it is, do the swap in both upper and lower triangles
      sp[s,r,c] <- ifelse(sp[s,r,c] == 1, 0, 1)
    }
  }
  return(sp)
})
trial_swap_gbi <- cmpfun(function(gbi, trials){
  n <- ncol(gbi) #number of individuals
  t <- nrow(gbi) #number of groups
  for(i in 1:trials){
    c <- sample(n,2) #random columns
    r <- sample(t,2) #random rows
    trial.matrix <- gbi[r,c] #get the trial swap submatrix
    if(all(rowSums(trial.matrix) == 1) & all(colSums(trial.matrix) == 1)){ #is a swap possible
      gbi[r,c] <- ifelse(gbi[r,c] == 1, 0, 1) #swap
    }
  }
  return(gbi) #return the swapped matrix
})

# simulation functions ####

#Functions to run a single simulation for each scenario, given the dataset of LHS samples (without transformation) and a row index
run_edge_sim <- cmpfun(function(data, index){

  mu <- qunif(data[index,1],0.01,0.5) #mean association index
  phi <- qunif(data[index,2],1,10) #association index precision
  N <- floor(qunif(data[index,3],20,101)) #population size
  t <- floor(qunif(data[index,4],20,201)) #number of sampling periods
  pobs <- qunif(data[index,5],0.1,1) #observation probability

  alpha <- mu*phi
  beta <- (1-mu)*phi

  sim_data <- generate_sp(n.ids = N, assoc.shape1 = alpha, assoc.shape2 = beta, n_samp = t, sight.prob = pobs) #generate SP and sightings matrices

  #subset to individuals that were observed on at least one day
  sim_data[[1]] <- sim_data[[1]][,colSums(sim_data[[2]]) > 0, colSums(sim_data[[2]]) > 0]
  sim_data[[2]] <- sim_data[[2]][,colSums(sim_data[[2]]) > 0]

  x <- apply(sim_data[[1]], c(2,3), sum) #number of associations per dyad
  d <- apply(sim_data[[2]], 2, function(z){
    apply(sim_data[[2]], 2, function(x){
      sum(1 - ((1-z)*(1-x))) #this will be 1 if either present, 0 otherwise
    })
  }) #number of sampling periods per dyad

  sri <- x/d #SRI

  N <- ncol(sri) #re-define N based on remaining individuals

  trait <- runif(N,0,1) #randomly assign trait
  trait.diff <- sapply(trait,function(z)abs(trait-z)) # differences in trait values

  coef.obs <- coef(lm(sri[lower.tri(sri)] ~ trait.diff[lower.tri(trait.diff)]))[2] #fit linear model
  t.obs <- summary(lm(sri[lower.tri(sri)] ~ trait.diff[lower.tri(trait.diff)]))$coef[2,3]
  p.obs <-  summary(lm(sri[lower.tri(sri)] ~ trait.diff[lower.tri(trait.diff)]))$coef[2,4]
  
  sp.p <- sim_data[[1]] #copy the SP array

  coef.p <- t.p <- NA #vectors to hold permuted results

  coef.np <- t.np <- NA
  
  pb <- txtProgressBar(min = 1, max = 10000, style = 3)
  
  for(j in 1:10000){
    setTxtProgressBar(pb,j)
    np <- sample(N)
    sri.np <- sri[np,np]
    lm.np <- lm(sri.np[lower.tri(sri.np)] ~ trait.diff[lower.tri(trait.diff)])
    coef.np[j] <- coef(lm.np)[2]
    t.np[j] <- summary(lm.np)$coef[2,3]
    
    sp.p <- trial_swap_sp(sp.p,1000) #do the swaps
    sri.p <- apply(sp.p,c(2,3),sum)/d
    t.p[j] <- summary(lm(sri.p[lower.tri(sri.p)] ~ trait.diff[lower.tri(trait.diff)]))$coef[2,3]
    coef.p[j] <- coef(lm(sri.p[lower.tri(sri.p)] ~ trait.diff[lower.tri(trait.diff)]))[2] #permuted linear model
  }

  p.coef <- p_val(coef.obs,coef.p)
  p.t <- p_val(t.obs,t.p)
  p.t.np <- p_val(t.obs,t.np)
  p.coef.np <- p_val(coef.obs,coef.np)
  sd_sri <- sd(sri[lower.tri(sri)])
  mean_sri <- mean(sri[lower.tri(sri)])
  pop_size <- ncol(sim_data[[2]])
  assoc_per <- mean(colSums(apply(sim_data[[1]], c(2,3), sum)))
  days_per <- mean(colSums(sim_data[[2]]))
  
  return(c(p.coef,p.t,p.obs,p.t.np,p.coef.np,sd_sri,days_per,ncol(sri)))
  
})
run_node_sim <- cmpfun(function(data, index){

  N <- floor(qunif(data[index,1], 20, 101)) #pop size
  groups <- floor(qunif(data[index,2], 20, 501)) #number of observed groupings
  max_size <- floor(qunif(data[index,3], 5, 10)) #maximum group size
  sd_pref <- qunif(data[index,4], 0.01, 2) #SD of size preferences

  gbi <- generate_gbi(N,groups,max_size,sd_pref) #generate to GBI
  gbi <- gbi[,colSums(gbi) > 0] #subset to individuals that were observed at all
  sri <- quiet(get_network(gbi)) #get the SRI network
  N <- ncol(sri) #redefine N

  trait <- runif(N,0,1) #randomly generate trait
  coef.obs <- coef(lm(colSums(sri) ~ trait))[2]
  t.obs <- summary(lm(colSums(sri) ~ trait))$coef[2,3]
  p.obs <- summary(lm(colSums(sri) ~ trait))$coef[2,4]
  
  gbi.p <- gbi
  coef.p <-  t.p <- NA
  coef.np <- t.np <- NA
  
  pb <- txtProgressBar(min = 1, max = 10000, style = 3)
  
  for(j in 1:10000){
    
    setTxtProgressBar(pb,j)
    
    lm.np <- lm(sample(colSums(sri)) ~ trait)
    coef.np[j] <- coef(lm.np)[2]
    t.np[j] <- summary(lm.np)$coef[2,3]
    
    gbi.p <- trial_swap_gbi(gbi.p,1000)
    sri.p <- quiet(get_network(gbi.p))
    coef.p[j] <- coef(lm(colSums(sri.p) ~ trait))[2]
    t.p[j] <- summary(lm(colSums(sri.p) ~ trait))$coef[2,3]
    
  }

  p.coef <- p_val(coef.obs,coef.p)
  p.t <- p_val(t.obs,t.p)
  p.coef.np <- p_val(coef.obs,coef.np)
  p.t.np <- p_val(t.obs,t.np)
  sd_str <- sd(colSums(sri))
  days_per <- mean(colSums(gbi))
  mean_str <- mean(colSums(sri))
  
  return(c(p.coef,p.t,p.obs,p.t.np,p.coef.np,sd_str,days_per,ncol(sri)))

})

# Run simulations

node_lhs <- randomLHS(k = 4, n = 200)
edge_lhs <- randomLHS(k = 5, n = 200)

# p.ds.coef and p.ds.t are the datastream p-values using the coefficient and t-value
# same naming scheme for np, which are from node permutations.
# p.orig is the original p-value from the linear model

edge_data <- t(sapply(1:2, function(z) run_edge_sim(edge_lhs, z) ))
colnames(edge_data) <- c("p.ds.coef","p.ds.t","p.orig","p.np.t","p.np.coef","sd_sri","days_per","N")
edge_data <- as.data.frame(edge_data)

node_data <- t(sapply(1:2, function(z) run_node_sim(node_lhs, z) ))
colnames(node_data) <- c("p.ds.coef","p.ds.t","p.orig","p.np.t","p.np.coef","sd_str","days_per","N")
node_data <- as.data.frame(node_data)