##Continuous SI model for identifying equilibrium of two viral types

library(deSolve)

## Create an SIR function
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- (-beta * S * I) + (gamma * I) + (gamma1 * I1) + (-beta1 * S * I1)
    dI <-  (beta * S * I) - (gamma * I) - (delta * I)
    dI1<- (beta1 * S * I1) - (gamma1 * I1) + (delta * I)
    
    return(list(c(dS, dI,dI1)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001
init       <- c(S = 0.8, I = 0.1,I1 = 0.1)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 3, gamma = 1,delta = 0.000, beta1 = 3, gamma1 = 1)
## Time frame
times      <- seq(0, 1000, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data

##Simulations for different combinations of transmission rates, virulence rates and ratios of the two
tran_vir_types_ratio = na.omit(as.data.frame(matrix(ncol = 6)))
library(deSolve)
count = 1
for (i in c(seq(0.1,10,0.1))){
  for (j in c(seq(0.1,10,0.1))){
    for (k in c(seq(1.5,20,0.5))){
      ## Create an SIR function
      sir <- function(time, state, parameters) {
        
        with(as.list(c(state, parameters)), {
          
          dS <- (-beta * S * I) + (gamma * I) + (gamma1 * I1) + (-beta1 * S * I1)
          dI <-  (beta * S * I) - (gamma * I) - (delta * I)
          dI1<- (beta1 * S * I1) - (gamma1 * I1) + (delta * I)
          
          return(list(c(dS, dI,dI1)))
        })
      }
      
      ### Set parameters
      ## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
      init       <- c(S = 0.8, I = 0.1,I1 = 0.1)
      ## beta: infection parameter; gamma: recovery parameter
      parameters <- c(beta = i, gamma = i/k,delta = 0.000, beta1 = j, gamma1 = j/k)
      ## Time frame
      times      <- seq(0, 500, by = 1)
      
      ## Solve using ode (General Solver for Ordinary Differential Equations)
      out <- ode(y = init, times = times, func = sir, parms = parameters)
      ## change to data frame
      out <- as.data.frame(out)
      ## Delete time variable
      out$time <- NULL
      ## Show data
      tran_vir_types_ratio[count,] = c(k,i,j,out[500,])
      count = count + 1
    }
  }
}



HU_parameters <- c(beta = 0.004, gamma = 0.003,delta = 0.000, beta1 = 0.002, gamma1 = 0.0015)
SR_parameters <- c(beta = 0.005454545, gamma = 0.003,delta = 0.000, beta1 = 0.0027272727, gamma1 = 0.0015)
CH_parameters <- c(beta = 0.012, gamma = 0.003,delta = 0.000, beta1 = 0.036, gamma1 = 0.009)
PR_parameters <- c(beta = 0.015, gamma = 0.003,delta = 0.000, beta1 = 0.05, gamma1 = 0.01)


library(deSolve)
## Create an SIR function
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- (-beta * S * I) + (gamma * I) + (gamma1 * I1) + (-beta1 * S * I1)
    dI <-  (beta * S * I) - (gamma * I) - (delta * I)
    dI1<- (beta1 * S * I1) - (gamma1 * I1) + (delta * I)
    return(list(c(dS, dI,dI1)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 0.8, I = 0.1,I1 = 0.1)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 0.1, gamma = 0.2,delta = 0.000, beta1 = 0.05, gamma1 = 0.1)
## Time frame
times      <- seq(0, 1000, by = 1)
## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL

out_params = as.data.frame(matrix(ncol = 5))
init <- c(S = 0.8, I = 0.1,I1 = 0.1)
times      <- seq(0, 1000, by = 1)
for (i in c(0.1,0.2, 0.3, 0.4, 0.5,0.66, 0.75,1,1.5,2,2.5,3,4,5,6,7,8,9,10,12,15,20,25,30,35,40,45,50,75,100)){
  trans = 0.1
  vir = trans/i
  parameters <- c(beta = trans, gamma = vir,delta = 0.000, beta1 = trans, gamma1 = vir)
  out <- ode(y = init, times = times, func = sir, parms = parameters)
  out <- as.data.frame(out)
  out$time <- NULL
  out_freq = c(trans,vir,as.vector(t(out[1001,])))
  out_params = rbind(out_params,out_freq)
}
ggplot(out_params,aes(x=V1/V2,y=V4+V5)) + geom_line(size=2) + scale_x_log10() + theme_bw() + theme(text=element_text(family = "sans",size = 18)) + labs(x="Transmission/Virulence",y="Stable Infection Frequency")

ggplot(out_params,aes(x=V2/V1,y=V4+V5)) + geom_line(size=2) + scale_x_log10() + theme_bw() + theme(text=element_text(family = "sans",size = 18)) + labs(x="Virulence/Transmission",y="Stable Infection Frequency")
ggplot(out_params,aes(x=V2/V1,y=V4+V5)) + geom_line(size=2) + scale_y_log10() + scale_x_log10() + theme_bw() + theme(text=element_text(family = "sans",size = 18)) + labs(x="Virulence/Transmission",y="Stable Infection Frequency")


out_balance = as.data.frame(matrix(ncol = 5))
init <- c(S = 0.8, I = 0.1,I1 = 0.1)
times      <- seq(0, 1000, by = 1)
for (i in seq(0.02,2,0.02)){
  for (j in seq(0.02,2,0.02)){
    trans1=i
    vir1=i/5
    trans2=j
    vir2=j/5
    parameters <- c(beta = trans1, gamma = vir1,delta = 0.000, beta1 = trans2, gamma1 = vir2)
    out <- ode(y = init, times = times, func = sir, parms = parameters)
    out <- as.data.frame(out)
    out$time <- NULL
    out_freq = c(i,j,as.vector(t(out[1001,])))
    out_balance = rbind(out_balance,out_freq)
  }
}
vir_balance_plot = ggplot(out_balance,aes(x=V1,y=V2,fill=log(V5/V4))) + theme_bw() + theme(legend.position = "bottom",text=element_text(size=18)) + geom_tile() + scale_fill_gradient2(low="#377EB8",mid="white",high="#E41A1C",midpoint = 0)
ggplot(out_balance,aes(x=V1,y=V2,fill=V4/V5)) + theme_bw() + theme(legend.position = "bottom",text=element_text(size=18)) + geom_tile() + scale_fill_gradient2(low="#377EB8",mid="white",high="#E41A1C",midpoint = 1,limits =c(0.147,6.8461))

out_params_extra = as.data.frame(matrix(ncol = 5))
init <- c(S = 0.8, I = 0.1,I1 = 0.1)
times      <- seq(0, 1000, by = 1)
for (j in c(0.0001,0.0005,0.001,0.0025,0.005,0.01)){
  for (i in c(0.1,0.2, 0.3, 0.4, 0.5,0.66, 0.75,1,1.5,2,2.5,3,4,5,6,7,8,9,10,12,15,20,25,30,35,40,45,50,75,100)){
    trans = j
    vir = trans/i
    parameters <- c(beta = trans, gamma = vir,delta = 0.000, beta1 = trans, gamma1 = vir)
    out <- ode(y = init, times = times, func = sir, parms = parameters)
    out <- as.data.frame(out)
    out$time <- NULL
    out_freq = c(trans,vir,as.vector(t(out[1001,])))
    out_params_extra = rbind(out_params_extra,out_freq)
  }
}
out_params_extra = na.omit(out_params_extra)

vir_params_plot = ggplot(out_params_extra,aes(x=V2/V1,y=V4+V5,col=as.factor(V1),linetype=as.factor(V1))) + geom_line(size=2) + scale_x_log10() + theme_bw() + theme(legend.position = "bottom",text=element_text(family = "sans",size = 18)) + labs(x="Virulence/Transmission",y="Stable Infection Frequency",col="Transmission Rate") + scale_color_grey()
multiplot(vir_params_plot,vir_balance_plot)




out_balance_rob = as.data.frame(matrix(ncol = 5))
init <- c(S = 0.8, I = 0.1,I1 = 0.1)
times      <- seq(0, 1000, by = 1)
for (i in seq(0.02,2,0.02)){
  for (j in seq(0.02,2,0.02)){
    trans1=i
    vir1=i/5
    trans2=j
    vir2=j/5
    parameters <- c(beta = trans1, gamma = vir1,delta = 0.000, beta1 = trans2, gamma1 = vir2)
    out <- ode(y = init, times = times, func = sir, parms = parameters)
    out <- as.data.frame(out)
    out$time <- NULL
    out_freq = c(i,j,as.vector(t(out[1001,])))
    out_balance_rob = rbind(out_balance_rob,out_freq)
  }
}

library(deSolve)
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- (-beta * S * I) + (gamma * I) + (gamma1 * I1) + (-beta1 * S * I1)
    dI <- (beta * S * I) - (gamma * I)
    dI1<- (beta1 * S * I1) - (gamma1 * I1)
    return(list(c(dS, dI,dI1)))
  })
}


out_params_robZ = as.data.frame(matrix(ncol = 5))
init <- c(S = 0.899, I = 0.1,I1 = 0.001)
times      <- seq(0, 1000, by = 1)
for (j in seq(0.0001,0.1,0.0001)){
  vir1 = 5e-04
  trans1 = vir1*4
  vir2 = j
  trans2 = vir2*4
  parameters <- c(beta = trans1, gamma = vir1, beta1 = trans2, gamma1 = vir2)
  out <- ode(y = init, times = times, func = sir, parms = parameters)
  out <- as.data.frame(out)
  out$time <- NULL
  out_freq = c(trans2,vir2,as.vector(t(out[1000,])))
  out_params_robZ = rbind(out_params_robZ,out_freq)
}
out_params_robZ = na.omit(out_params_robZ)
colnames(out_params_robZ) = c("Transmission","Virulence","S","I1","I2")
out_params_robZ2 = melt(data = out_params_robZ,id.vars = "Virulence",measure.vars = c("S","I1","I2"))
out_params_robZ3 = droplevels(out_params_robZ2[out_params_robZ2$Virulence != 0,])

out_params_robX2 = melt(data = out_params_robX,id.vars = "Transmission",measure.vars = c("S","I1","I2"))





## Discrete simulation for wait time to fixation of high titer type

recomb = function(type1,type2,samp_val){
  r1 = type1
  if (samp_val == 1){
    r1[,c(5,6,7,8)] =  type2[,c(5,6,7,8)]
  }
  else if (samp_val == 2){
    r1[,c(6,7,8)] =  type2[,c(6,7,8)]
  }
  else if (samp_val == 3){
    r1[,c(7,8)] =  type2[,c(7,8)]
  }
  else if (samp_val == 4){
    r1[,c(8)] =  type2[,c(8)]
  }
  else if (samp_val == 5){
    r1[,c(4,5,6,7)] =  type2[,c(4,5,6,7)]
  }
  else if (samp_val == 6){
    r1[,c(4,5,6)] =  type2[,c(4,5,6)]
  }
  else if (samp_val == 7){
    r1[,c(4,5)] =  type2[,c(4,5)]
  }
  else if (samp_val == 8){
    r1[,c(4)] =  type2[,c(4)]
  }
  else if (samp_val == 9){
    r1[,c(5)] =  type2[,c(5)]
  }
  else if (samp_val == 10){
    r1[,c(6)] =  type2[,c(6)]
  }
  else if (samp_val == 11){
    r1[,c(7)] =  type2[,c(7)]
  }
  return(r1)
}

new_mut = function(titre,mut1,mut2,mut3,mut4,mut5){
  num_muts = sum(c(mut1,mut2,mut3,mut4,mut5) == 1)
  new_titre = titre*((1.5^num_muts))
  m1 = mut1
  m2 = mut2
  m3 = mut3
  m4 = mut4
  m5 = mut5
  if (m1 == 0){m1 = as.character(sample(c(0,1),1,prob = c((1-new_titre*0.00001),(new_titre*0.00001))))}
  if (m2 == 0){m2 = as.character(sample(c(0,1),1,prob = c((1-new_titre*0.00001),(new_titre*0.00001))))}
  if (m3 == 0){m3 = as.character(sample(c(0,1),1,prob = c((1-new_titre*0.00001),(new_titre*0.00001))))}
  if (m4 == 0){m4 = as.character(sample(c(0,1),1,prob = c((1-new_titre*0.00001),(new_titre*0.00001))))}
  if (m5 == 0){m5 = as.character(sample(c(0,1),1,prob = c((1-new_titre*0.00001),(new_titre*0.00001))))}
  if (m1 == 1){m1 = as.character(sample(c(1,0),1,prob = c((1-new_titre*0.00001),(new_titre*0.00001))))}
  if (m2 == 1){m2 = as.character(sample(c(1,0),1,prob = c((1-new_titre*0.00001),(new_titre*0.00001))))}
  if (m3 == 1){m3 = as.character(sample(c(1,0),1,prob = c((1-new_titre*0.00001),(new_titre*0.00001))))}
  if (m4 == 1){m4 = as.character(sample(c(1,0),1,prob = c((1-new_titre*0.00001),(new_titre*0.00001))))}
  if (m5 == 1){m5 = as.character(sample(c(1,0),1,prob = c((1-new_titre*0.00001),(new_titre*0.00001))))}
  return(c(new_titre,m1,m2,m3,m4,m5))
}

infsim_reps = na.omit(as.data.frame(matrix(ncol=5)))
infsim_types = na.omit(as.data.frame(matrix(ncol=7)))
colnames(infsim_types) = c("count","mut1","mut2","mut3","mut4","mut5","freq")

for (j in c(1:100)){
  mut_gain = c()
  mut_titre = c()
  new_m1 = 0
  new_m2 = 0
  new_m3 = 0
  new_m4 = 0
  new_m5 = 0
  n=1000
  inf_rate=0.49
  death_rate = 0.11
  mut = 0.000001
  u0 <- data.frame(species = rep("A",n),state  = c(rep("S",990),rep("I",10)),titre = 0,mut1 = 0,mut2 = 0,mut3 = 0,mut4 = 0,mut5 = 0)
  u0$titre[u0$state == "I"] = rpois(length(u0$titre[u0$state == "I"]),lambda = 10)
  u0[,c(3,4,5,6,7,8)] = as.data.frame(t(mapply(new_mut,u0$titre,u0$mut1,u0$mut2,u0$mut3,u0$mut4,u0$mut5)))
  u0$titre = as.numeric(as.character(u0$titre))
  u0$mut1 = as.numeric(as.character(u0$mut1))
  u0$mut2 = as.numeric(as.character(u0$mut2))
  u0$mut3 = as.numeric(as.character(u0$mut3))
  u0$mut4 = as.numeric(as.character(u0$mut4))
  u0$mut5 = as.numeric(as.character(u0$mut5))
  count = 1
  un_1 = u0
  for (i in c(1:1000)){
    count = count + 1
    x = ddply(un_1, .(state,mut1, mut2,mut3,mut4,mut5), summarize,mean_titre = mean(titre))
    un_2 = na.omit(as.data.frame(matrix(ncol=8)))
    colnames(un_2) = colnames(u0)
    for (v in c(1:length(x$mean_titre))){
      infs = un_1[un_1$state == "I" & un_1$mut1 == x$mut1[v] & un_1$mut2 == x$mut2[v] & un_1$mut3 == x$mut3[v] & un_1$mut4 == x$mut4[v] & un_1$mut5 == x$mut5[v],]
      num_inf = ((sqrt(x$mean_titre[v]) * inf_rate * (length(un_1$state[un_1$state == "S"])/n) * (length(infs$state)/n))*n)
      num_dead = ((sqrt(x$mean_titre[v]) * death_rate * (length(infs$state)/n))*n)
      diffs = rpois(1,num_inf) - rpois(1,num_dead)
      xv = c("species" = "A","state" = as.character(x$state[v]),"titre" = 0,"mut1" = x$mut1[v],"mut2" = x$mut2[v],"mut3" = x$mut3[v],"mut4" = x$mut4[v],"mut5" = x$mut5[v])
      if (diffs > 0){
        un_2 = rbind(un_2,t(matrix(xv,ncol=diffs,nrow = 8)))
      }
    }
    cc = c("species" = "A","state" = "S","titre" = 0,"mut1" = 0,"mut2" = 0,"mut3" = 0,"mut4" = 0,"mut5" = 0)
    if (dim(un_2)[1] < n){
      un_2 = rbind(un_2,t(matrix(cc,ncol=n-dim(un_2)[1],nrow = 8)))
    }
    un_1 = un_2
    colnames(un_2) = colnames(u0)
    colnames(un_1) = colnames(u0)
    un_1$titre = as.numeric(as.character(un_1$titre))
    un_1$mut1 = as.numeric(as.character(un_1$mut1))
    un_1$mut2 = as.numeric(as.character(un_1$mut2))
    un_1$mut3 = as.numeric(as.character(un_1$mut3))
    un_1$mut4 = as.numeric(as.character(un_1$mut4))
    un_1$mut5 = as.numeric(as.character(un_1$mut5))
    un_1$titre[un_1$state == "I"] = rpois(length(un_1$titre[un_1$state == "I"]),lambda = 10)
    un_1[,c(3,4,5,6,7,8)] = as.data.frame(t(mapply(new_mut,un_1$titre,un_1$mut1,un_1$mut2,un_1$mut3,un_1$mut4,un_1$mut5)))
    un_1$titre = as.numeric(as.character(un_1$titre))
    un_1$mut1 = as.numeric(as.character(un_1$mut1))
    un_1$mut2 = as.numeric(as.character(un_1$mut2))
    un_1$mut3 = as.numeric(as.character(un_1$mut3))
    un_1$mut4 = as.numeric(as.character(un_1$mut4))
    un_1$mut5 = as.numeric(as.character(un_1$mut5))
    num_recs = rpois(1,nrow(un_1[un_1$state == "I",])/10000)
    if (num_recs > 0){
      for (k in c(1:num_recs)){
        r1 = un_1[un_1$state == "I",][sample(nrow(un_1[un_1$state == "I",]), 1), ]
        r2 = un_1[un_1$state == "I",][sample(nrow(un_1[un_1$state == "I",]), 1), ]
        samp_val = sample(c(1,2,3,4,5,6,7,8,9,10,11),1)
        out_r1 = recomb(r1,r2,samp_val)
        un_1[as.numeric(rownames(out_r1)),] = out_r1
      }
    }
    print(c(count,length(un_1$state[un_1$state == "I"]),length(un_1$mut1[un_1$mut1 == 1]),length(un_1$mut1[un_1$mut2 == 1]),length(un_1$mut1[un_1$mut3 == 1]),length(un_1$mut1[un_1$mut4 == 1]),length(un_1$mut1[un_1$mut5 == 1])))
    if (length(un_1$mut1[un_1$mut1 == 1]) > 10 & new_m1 == 0){
      mut_gain = c(mut_gain,as.numeric(as.character(i)))
      new_m1 = 1
    }
    if (length(un_1$mut2[un_1$mut2 == 1]) > 10 & new_m2 == 0){
      mut_gain = c(mut_gain,as.numeric(as.character(i)))
      new_m2 = 1
    }
    if (length(un_1$mut3[un_1$mut3 == 1]) > 10 & new_m3 == 0){
      mut_gain = c(mut_gain,as.numeric(as.character(i)))
      new_m3 = 1
    }
    if (length(un_1$mut2[un_1$mut4 == 1]) > 10 & new_m4 == 0){
      mut_gain = c(mut_gain,as.numeric(as.character(i)))
      new_m4 = 1
    }
    if (length(un_1$mut3[un_1$mut5 == 1]) > 10 & new_m5 == 0){
      mut_gain = c(mut_gain,as.numeric(as.character(i)))
      new_m5 = 1
    }
  }
  if (length(un_1$state[un_1$state == "I"]) > 0){
    if (length(mut_gain) < 5){mut_gain = c(mut_gain,rep(1000,5-length(mut_gain)))}
    print(mut_gain)
    infsim_reps[j,] = mut_gain
    kk = as.data.frame(c(j,count(un_1[un_1$state == "I",], c("mut1","mut2","mut3","mut4","mut5"))))
    colnames(kk) = colnames(infsim_types)
    infsim_types = rbind(infsim_types,kk)
  }
}

