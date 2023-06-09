library("stats", lib.loc="C:/Program Files/R/R-3.3.1/library")
library(nnet)
library(metafor)
library(MASS)

ss <- 750*5
si = 5000

res = matrix(0, nrow = si, ncol = 79)

for (z in 1:si){
  
  set.seed(z)
  
  se = diag(0.5,5)
  l = mvrnorm(n = ss, mu = c(-0.5,-0.25,0,0.25,0.5), Sigma = se)
  
  ## Simulate "study" by a multinomial logit model
  s = matrix(nrow = 5, ncol = ss)
  beta = rbind(c(0.25,    1, 0.55, -0.25, 0.75, -0.25, -0.5),
               c(0.25, 0.35, 0.45, -0.35, 0.65,  0.25,  0.5),
               c(0.25,    1, 0.35, -0.45, 0.55,  0.25,  0.5),
               c(0.25, 0.55, 0.25, -0.55, 0.45, -0.25, -0.5))
  L = cbind(1,l, l[,1]*l[,2])
  enum.lp = rep(1,ss) + exp(L%*%beta[1,]) + exp(L%*%beta[2,]) + exp(L%*%beta[3,]) + exp(L%*%beta[4,])
  lp = cbind(1/enum.lp,
             exp(L%*%beta[1,])/enum.lp,
             exp(L%*%beta[2,])/enum.lp,
             exp(L%*%beta[3,])/enum.lp,
             exp(L%*%beta[4,])/enum.lp)
  
  study = rep(0,ss)
  for (i in 1:ss){
    s[,i] = rmultinom(1, size =1, prob = lp[i,]) ##
    study[i] = which(s[,i]==1)  
  } 
  
  new = data.frame(l = l, study = study, l.12 = l[,1]*l[,2])
  
  new$x = rbinom(ss,1,0.5)
  
  lp = -0.25 - 0.5*new$x + 0.5*new$l.1 - 0.25*new$l.2 - 0.5*new$l.3 + 0.25*new$l.4 - 0.5*new$l.5 - 1.5*new$x*new$l.1
  
  new$y = rbinom(ss,1,exp(lp)/(1 + exp(lp)))
  
  ########################################################
  
  # Design matrix: X
  X <- cbind(1,new$x,new$l.1, new$l.2, new$l.3, new$l.4, new$l.5,new$x*new$l.1)
  
  # K logistic regression models
  a = c("x + l.1 + l.2 + l.3 + l.4 + l.5 + x:l.1",
        "x + l.1 + l.2 + l.3 + l.4 + l.5 + x:l.1",
        "x + l.1 + l.2 + l.3 + l.4 + l.5 + x:l.1",
        "x + l.1 + l.2 + l.3 + l.4 + l.5 + x:l.1",
        "x + l.1 + l.2 + l.3 + l.4 + l.5 + x:l.1")
  
  K = length(a)
  
  
  # fits the models
  lmod = lapply(1:K, 
                function(k) 
                  glm(as.formula(paste("y~", a[k])), family = "binomial", 
                      data = new, subset=study==k))
  # Then keep the coefficients
  beta <- lapply(1:K, function (k) coef(lmod[[k]]))
  
  # and the indices corresponding variables in the design matrix
  idx = list()
  idx[[1]] = 1:8
  idx[[2]] = 1:8
  idx[[3]] = 1:8
  idx[[4]] = 1:8
  idx[[5]] = 1:8
  
  
  new_pred1 = new
  new_pred1$x = 1
  pre1 = lapply(1:K, 
                function(k) 
                  predict(lmod[[k]], newdata = new_pred1, type="response"))
  
  new_pred0 = new
  new_pred0$x = 0
  pre0 = lapply(1:K, 
                function(k) 
                  predict(lmod[[k]], newdata = new_pred0, type="response"))
  
  # vector of all RRjk's for one j
  method2 = function(j) sapply(1:K, function(k) sum(pre1[[k]][which(new$study == j)])/
                                 sum(pre0[[k]][which(new$study == j)]))
  # all RRjk
  RR = as.vector(sapply(1:K,method2)) #RR11, RR12, RR21, RR22
  
  
  # y: vector of responses, X: design matrix, s: trials, betas: list of all betas, idxs: list of all indices, rrs: vector of RRs
  func <- function(y,X,s,betas,idxs,rrs){
    # Scores for betas
    B <- matrix(NA,nrow=length(unlist(betas)),ncol=length(y))
    i <- 0
    for(k in 1:length(betas))
    {
      lp <- t(betas[[k]]) %*% t(X[,idxs[[k]]]) #lp = matrix (1 x ss)
      B[i+(1:length(betas[[k]])),] <- matrix(unlist(lapply(idxs[[k]], function(u) X[,u]*as.numeric(s==k)*(y-exp(lp)/(1+exp(lp))))), 
                                             ncol = length(s), nrow = length(beta[[k]]),byrow = T)	
      i <- max(which(!is.na(B[,1])))
      #print(i)
    }
    
    # Score for RRs
    R <- matrix(NA,nrow=length(betas)^2,ncol=length(y))
    
    Xt <- X
    Xt[,2] <- 1 
    Xt[,8] <- X[,3] # to make X:L = L
    
    Xnt <- X
    Xnt[,2] <- 0 
    Xnt[,8] <- 0 # to make X:L = 0 
    
    for(k in 1:length(betas))
    {
      lp1 <- t(betas[[k]]) %*% t(Xt[,idxs[[k]]])
      lp0 <- t(betas[[k]]) %*% t(Xnt[,idxs[[k]]])
      # Calculate RR.k's
      for(j in 1:length(betas))
        R[k+(j-1)*length(betas),] <- as.numeric(s == j) * (exp(lp1)/(1+exp(lp1)) - rrs[k+(j-1)*length(betas)]*exp(lp0)/(1+exp(lp0)))
    }
    return(rbind(B,R))	
  }
  
  phi = func(y = new$y, X = X, s = new$study, betas = beta, idxs = idx, rrs = RR)
  
  # matrix B: 12*12
  n.para = length(phi[,1])
  b= matrix (0, nrow = n.para, ncol = n.para)
  for (i in 1:ss) b = b + phi[,i] %*% t(phi[,i])
  b = b/ss
  
  # matrix A:
  h = 1e-5
  
  # Part of betas
  a1 = list()
  for (i in 1:length(beta)){
    beta. = beta
    identity = diag(length(beta.[[i]]))
    dev1 = unlist(lapply(1:length(beta[[i]]), function(u){
      beta.[[i]] = complex(real = beta[[i]], imaginary = h*identity[u,])
      out = rowSums(Im(func(y = new$y, X = X, s = new$study, betas = beta., idxs = idx, rrs = RR))/h)
      return(out)}))
    a1[[i]] = -matrix(dev1, nrow = length(unlist(beta))+length(RR), ncol = length(beta[[i]]))/ss
  }
  
  
  # Part of RRs
  identity = diag(length(RR))
  dev2 = unlist(lapply(1:length(RR), function(u){
    rr. = complex(real = RR, imaginary = h*identity[u,])
    out = rowSums(Im(func(y = new$y, X = X, s = new$study, betas = beta, idxs = idx, rrs = rr.))/h)
    return(out)}))
  a2 = -matrix(dev2, nrow = length(unlist(beta))+length(RR), ncol = length(RR))/ss
  
  # Paste the 2 part together
  a = a1[[1]]
  for (u in 2:K){
    a = cbind(a,a1[[u]])
  }
  a = cbind(a,a2)
  
  # sandwich estimator of variance
  var.para= solve(a) %*% b %*% t(solve(a))
  
  # logRR and SE of logRR
  logRR = log(RR)
  n.rr = length(RR) # number of RRjk's
  n.para = length(RR) + length(unlist(beta))
  up = n.para - n.rr + 1
  var.logRR = diag(diag(1/RR)%*%var.para[up:n.para, up:n.para]%*%diag(1/RR)/ss)
  
  # Covariance matrix of logRR
  v = diag(1/RR) %*% var.para[up:n.para, up:n.para] %*% diag(1/RR)/ss
  
  ############################
  
  m_logRR = matrix(logRR, ncol = K, nrow = K, byrow = T)
  m_var = matrix(var.logRR, ncol = K, nrow = K, byrow = T)
  
  
  meta_j = function(j){
    m = rma(yi = m_logRR[j,], vi = m_var[j,], method = 'DL')
    
    l = matrix(0,ncol = K*K, nrow = K - 1)
    l[,1 + K*(j-1)] = 1
    nr = K - 1
    for (u in 1:nr) l[u, 1 + u + K*(j - 1)] = -1
    stat = t(l%*%logRR) %*% solve(l %*% v %*% t(l)) %*% (l%*%logRR)
    p.value = 1 - pchisq(stat, df = nr)
    
    return(c(as.numeric(m$b), p.value, summary(m)$'tau2', summary(m)$'I2'))
  }
  
  
  meta_k = function(k){
    l = matrix(0,ncol = K*K, nrow = K - 1)
    l[,k] = 1
    nr = K - 1
    for (u in 1:nr) l[u, k + K*u] = -1
    stat = t(l%*%logRR) %*% solve(l %*% v %*% t(l)) %*% (l%*%logRR)
    p.value = 1 - pchisq(stat, df = nr)
    return(p.value)
  }
  
  # result MA: j = 1
  res[z,1:4] = meta_j(1)
  
  # result MA: j = 2
  res[z,5:8] = meta_j(2)
  
  # result MA: j = 3
  res[z,9:12] = meta_j(3)
  
  # result MA: j = 4
  res[z,13:16] = meta_j(4)
  
  # result MA: j = 5
  res[z,17:20] = meta_j(5)
  
  # result MA: k = 1
  res[z,21] = meta_k(1)
  # result MA: k = 2
  res[z,22] = meta_k(2)
  # result MA: k = 3
  res[z,23] = meta_k(3)  
  # result MA: k = 4
  res[z,24] = meta_k(4)
  # result MA: k = 5
  res[z,25] = meta_k(5)
  
  
  # result MA: conventional
  meta.normal = rma(yi = diag(m_logRR), vi = diag(m_var), method = 'DL')
  res[z,26] = as.numeric(meta.normal$b)
  
  lt = matrix(0,ncol = K*K, nrow = K - 1)
  lt[,1] = 1
  lt[1,7] = -1; lt[2,13] = -1; lt[3,19] = -1; lt[4,25] = -1
  nr = K - 1
  stat = t(lt%*%logRR) %*% solve(lt %*% v %*% t(lt)) %*% (lt%*%logRR)
  res[z,27] = 1 - pchisq(stat, df = nr)
  
  res[z,28:29] = c(summary(meta.normal)$'tau2', summary(meta.normal)$'I2')
  
  # Save the 25 estimates
  res[z,30:54] = RR
  # Save the 25 variances
  res[z,55:79] = var.logRR
}  

# analysis

## boxplot of the 4 estimates
boxplot(res[,30:54])

# Test: BCM
sapply(c(2,6,10,14,18), function(i) length(res[,i][which(res[,i] <= 0.05)])/si)

# Test: CM
sapply(21:25, function(i) length(res[,i][which(res[,i] <= 0.05)])/si)

# Test: conventional
length(res[,27][which(res[,27] <= 0.05)])/si

# Tau2: j = 1 to 5; standard
summary(res[,c(3,7,11,15,19,28)])

#I2: j = 1 to 5; standard
summary(res[,c(4,8,12,16,20,29)])

# Estimate: logRR(j.)
sapply(c(1,5,9,13,17), function(i) mean(res[,i]))
sapply(c(1,5,9,13,17), function(i) var(res[,i]))

write.csv(res,"F:/PhD/PhD_article 1_5 trial setting/DS/DSRR_S1.csv")