## First, derive the predator and prey encounter functions
## estimates of mean probability of encounter per 3 h sampling interval.
ProbEncPred.per.t = c(0.002642006, 0.018348104, 0.007905097, 0.018348104, 0.010526218, 0.002642006, 0.000000000, 0.000000000)
ProbEncPrey.per.t = c(0.01574770, 0.08359932, 0.13997485, 0.19714144, 0.11223447, 0.10279025, 0.05404053, 0.01574770)

## fit 3rd order polynomials
ToDs = c(63,189,315,441,567,693,819,945) ### these are chosen to represent the midpoint of each 3 h block, beginning with the 3h block from 0600-0900, with midpoint at sixty-three 85.8 second time steps into the day, and each subsequent midpoint 3h (126 time steps) later

fit.pred = lm(ProbEncPred.per.t~ poly(ToDs,3,raw=T))
fit.prey = lm(ProbEncPrey.per.t~ poly(ToDs,3,raw=T))

### create new data frame for predicting encounter probabilities on single time step level
newdata.TOD = data.frame(ToDs = seq(1, 1007, 1))  ## new data frame
# predict for the predators
predict.fit.pred = predict(fit.pred, newdata = newdata.TOD)  ### predictions from new data frame
plot(predict.fit.pred~newdata.TOD$ToDs)  ### plotting, check for over-wigglyness/negative values
predict.fit.pred[which(predict.fit.pred<0)] = 0  ### get rid of negatives
plot(predict.fit.pred~newdata.TOD$ToDs)  ## check
predict.fit.pred[960:1007] = 0   ### get rid of artifactual wigglyness
plot(predict.fit.pred~newdata.TOD$ToDs)  ### check

which(predict.fit.pred == max(predict.fit.pred)) ## what time step does predation peak?

# now for the prey
predict.fit.prey = predict(fit.prey, newdata = newdata.TOD)
plot(predict.fit.prey~newdata.TOD$ToDs)
predict.fit.prey[which(predict.fit.prey<0)] = 0
plot(predict.fit.prey~newdata.TOD$ToDs)











### create a function with the following inputs
### pred = predator encounter probability, function of t (P(t) in MS), 
### foreff = probability that attack on prey by focal individual succeeds, (f in MS)
### m = probability that a predator encounter results in predation of focal individual (m in MS)
### c = basal metabolic rate in terms of proportion of energy in typical prey consumed per t (c1 in MS)
### value = the energetic value of a prey item, = 1 (energy is in units of prey items)
### metcostforage = scalar by which basal metabolic rate increases with foraging activity (s in MS)
### maxt = how many time steps per day (determined by total time period / handling time)
### maxc = maximum energetic state
### prey = prey encounter probability, function of t (R(t) in MS)
### c0 = inflection point of sigmoidal fitness function
### k = steepness of the fitness function near c0

foragedielsearch = function(pred, foreff, m, c, value, metcostforage, maxt, maxc, prey, c0, k){       
  
  ##### begin by creating the fitness function fx that can be used to initialize reward matrix
  Cs=0:maxc
  fx= 1/(1+(exp(-k*(Cs-c0))))  ### sigmoidal fitness function
  #plot(fx~Cs)
  
  ### optional plotting to show encounter functions and model predictions as a 2-panel plot.
  
  #quartz(width=4, height=6.5)
  #par(mfcol=c(2,1), mar=c(5,4,1.5,2))    ### set up a plot window to hold the encounter functions and predictions as panels
  
  #plot(pred*pforage~c(1:maxt), ylim=c(0,0.1), ylab="Probability of Event", xlab="Time of Day (hours)", type="n", xaxt="n")
 # lines(pred*pforage~c(1:maxt), col="red")    ### lines for predator encounter probability
 # lines(prey*foreff~c(1:maxt), col="blue")    ### lines for prey encounter probability
 # legend("topright", legend = c("Predation", "Prey Capture"), lty=1, col=c("red", "blue"), cex=0.8)
#  axis(1, at = c(0,127,254,381,508,635,762,889,1016), labels=c("0600", "0900", "1200", "1500", "1800", "2100", "0000", "0300", "0600"))
  
  Reward = matrix(0,nrow = maxc+1,ncol=maxt+1)   ### create empty reward matrix (holds expected fitness values for each combination of energetic reserves (rows) and time step t (column), assuming optimal behavior from that t forward) ###
  Reward[,maxt+1]=(fx)     ### set last column to fitness value (fx) for each level of reserves at T
  Result=matrix(0,nrow = maxc+1,ncol=maxt)     ### create empty result matrix, hold optimal forage decision (0,1) for each combination of energetic reserves and t
  
  
  for (t in maxt:1){          ### work backward from the end of the time interval.
    
    ## Calculate prob of predation event due to foraging  at t ##
    
    dforage = m * pred[t]
    
    ### Calculate probability of not dying and gaining and not dying and losing for each choice###
    
    p.forage.enc.up=(1-dforage)*prey[t]*foreff ## forages, is not predated, finds prey, captures
    p.forage.enc.lose=(1-dforage)*prey[t]*(1-foreff) ## forages, is not predated, finds prey, fails
    
    p.forage.noenc = (1-dforage)*(1-prey[t])  ## forages, is not predated, finds no prey.
    
    p.rest.lose=(1)      ### rests, has no chance of predation, but must probability of losing energy is 1.
    
    ### Use above probabilities to calculate exp. terminal reward for each choice for each possible level of energetic reserves ####
    
    for(i in 2:maxc){ ### loop over values of energetic reserves. start with i = 2, because i = 1 corresponds to energy reserves = 0 --> dead individual, no decision to make.
      
      ### calculate expected state for each possibility ###
      state.rest = i - c   ### individual rests, loses energy
      state.yes.enc.success = i + value - (metcostforage * c) ## individual gains energy, pays foraging cost
      state.yes.enc.failure = i - (metcostforage * c)  ### individual fails to capture, pays foraging cost
      state.yes.noenc = i - (metcostforage * c)  ### individual does not encounter prey, pays cost
      
      ### Need to find Expected Future Fitness (Reward) for each of the options using interpolation
      
      ### first, find the next lowest discrete value of i for each decision
      li.state.rest = floor(state.rest)
      li.state.yes.enc.success = floor(state.yes.enc.success)
      li.state.yes.enc.failure = floor(state.yes.enc.failure)
      li.state.yes.noenc = floor(state.yes.noenc)
      
      #### next, use interpolation to find expected reward for each state
      
      ### first, for not foraging
      if (state.rest == li.state.rest){
        Fit.rest = Reward[state.rest,t+1]
      }
      else{
        Fit.rest = ((state.rest-li.state.rest)*Reward[li.state.rest+1,t+1] + (((li.state.rest+1)-state.rest)*Reward[li.state.rest,t+1])) / ((li.state.rest+1) - li.state.rest)
      }
      
      ### for foraging, encounter, and success
      
      if (state.yes.enc.success == li.state.yes.enc.success){
        Fit.yes.enc.success = Reward[state.yes.enc.success,t+1]
      }
      else{
        Fit.yes.enc.success = ((state.yes.enc.success-li.state.yes.enc.success)*Reward[li.state.yes.enc.success+1,t+1] + (((li.state.yes.enc.success+1)-state.yes.enc.success)*Reward[li.state.yes.enc.success,t+1])) / ((li.state.yes.enc.success+1) - li.state.yes.enc.success)
      }
      
      ### for foraging, encounter, no success
      
      if (state.yes.enc.failure == li.state.yes.enc.failure){
        Fit.yes.enc.failure = Reward[state.yes.enc.failure,t+1]
      }
      else{
        Fit.yes.enc.failure = ((state.yes.enc.failure-li.state.yes.enc.failure)*Reward[li.state.yes.enc.failure+1,t+1] + (((li.state.yes.enc.failure+1)-state.yes.enc.failure)*Reward[li.state.yes.enc.failure,t+1])) / ((li.state.yes.enc.failure+1) - li.state.yes.enc.failure)
      }
      
      ### for foraging, no encounter
      
      if (state.yes.noenc == li.state.yes.noenc){
        Fit.yes.noenc = Reward[state.yes.enc.failure,t+1]
      }
      else{
        Fit.yes.noenc = ((state.yes.noenc-li.state.yes.noenc)*Reward[li.state.yes.noenc+1,t+1] + (((li.state.yes.noenc+1)-state.yes.noenc)*Reward[li.state.yes.noenc,t+1])) / ((li.state.yes.noenc+1) - li.state.yes.noenc)
      }
      
      
      #### now, set fitness for i and t (Reward) according to the likelihood of prey encounter and the fitness payoffs of whichever decision maximizes future fitness if prey is encountered.
      
      #Reward[i,t] = ifelse((p.rest.lose * Fit.rest)> ((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc)), (p.rest.lose * Fit.rest), ((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc)))
      
      ### can add some rounding to reduce artifacts
      Reward[i,t] = ifelse((p.rest.lose * Fit.rest)> ((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc))+ 0.0002, (p.rest.lose * Fit.rest), ((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc)))
      
      ### now, record what the optimal decision was in the Result matrix (1 = forage, 0 = rest)
      
      #Result[i,t] = ifelse(((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc)) > (p.rest.lose * Fit.rest), 1, 0)
      
      Result[i,t] = ifelse(((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc))+0.0002 > (p.rest.lose * Fit.rest), 1, 0) 
      
    }    ### this ends the loop over i, still in the t-loop
    
    
    ### individuals in top condition cannot gain energy reserves, but can lose them ###
    
    ### need new state change possibilities
    
    state.rest = i - c
    state.yes.enc.success = i + value - (metcostforage * c)
    state.yes.enc.failure = i - (metcostforage * c)
    state.yes.noenc = i - (metcostforage * c)
    
    state.rest = (maxc+1) - c    ### loses c from max condition
    state.yes.enc.success = (maxc+1) ### success covers costs of foraging, stays in max condition
    state.yes.enc.failure = ((maxc+1) - (metcostforage * c))
    state.yes.noenc = ((maxc+1) - (metcostforage * c)) ### loses cost of foraging from max condition
    
    
    ### Need to find Expected Future Fitness (Reward) for each of the options using interpolation
    
    ### first, find the next lowest discrete value of i for each decision
    li.state.rest = floor(state.rest)
    li.state.yes.enc.success = floor(state.yes.enc.success)
    li.state.yes.enc.failure = floor(state.yes.enc.failure)
    li.state.yes.noenc = floor(state.yes.noenc)
    
    #### next, use interpolation to find expected reward for each state
    
    ### first, for not foraging
    if (state.rest == li.state.rest){
      Fit.rest = Reward[state.rest,t+1]
    }
    else{
      Fit.rest = ((state.rest-li.state.rest)*Reward[li.state.rest+1,t+1] + (((li.state.rest+1)-state.rest)*Reward[li.state.rest,t+1])) / ((li.state.rest+1) - li.state.rest)
    }
    
    ### for foraging, encounter, and success
    
    if (state.yes.enc.success == li.state.yes.enc.success){
      Fit.yes.enc.success = Reward[state.yes.enc.success,t+1]
    }
    else{
      Fit.yes.enc.success = ((state.yes.enc.success-li.state.yes.enc.success)*Reward[li.state.yes.enc.success+1,t+1] + (((li.state.yes.enc.success+1)-state.yes.enc.success)*Reward[li.state.yes.enc.success,t+1])) / ((li.state.yes.enc.success+1) - li.state.yes.enc.success)
    }
    
    ### for foraging, encounter, no success
    
    if (state.yes.enc.failure == li.state.yes.enc.failure){
      Fit.yes.enc.failure = Reward[state.yes.enc.failure,t+1]
    }
    else{
      Fit.yes.enc.failure = ((state.yes.enc.failure-li.state.yes.enc.failure)*Reward[li.state.yes.enc.failure+1,t+1] + (((li.state.yes.enc.failure+1)-state.yes.enc.failure)*Reward[li.state.yes.enc.failure,t+1])) / ((li.state.yes.enc.failure+1) - li.state.yes.enc.failure)
    }
    
    ### for foraging, no encounter
    
    if (state.yes.noenc == li.state.yes.noenc){
      Fit.yes.noenc = Reward[state.yes.enc.failure,t+1]
    }
    else{
      Fit.yes.noenc = ((state.yes.noenc-li.state.yes.noenc)*Reward[li.state.yes.noenc+1,t+1] + (((li.state.yes.noenc+1)-state.yes.noenc)*Reward[li.state.yes.noenc,t+1])) / ((li.state.yes.noenc+1) - li.state.yes.noenc)
    }
    
    #### now, set fitness for i and t (Reward) according to the likelihood of prey encounter and the fitness payoffs of whichever decision maximizes future fitness if prey is encountered.
    
    #Reward[i,t] = ifelse((p.rest.lose * Fit.rest)> ((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc)), (p.rest.lose * Fit.rest), ((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc)))
    
    #rounding
    Reward[i,t] = ifelse((p.rest.lose * Fit.rest)> ((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc))+ 0.0002, (p.rest.lose * Fit.rest), ((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc)))
    
    ### now, record what the optimal decision was in the Result matrix (1 = forage, 0 = rest)
    
      #Result[i,t] = ifelse(((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc))> (p.rest.lose * Fit.rest), 1, 0) 
    
      ## rounding
     Result[i,t] = ifelse(((p.forage.enc.up * Fit.yes.enc.success) + (p.forage.enc.lose * Fit.yes.enc.failure) + (p.forage.noenc * Fit.yes.noenc))+0.0002 > (p.rest.lose * Fit.rest), 1, 0) 
    
    
  }   ### this ends the t-loop.
  
  
  boundary=vector("numeric", length=length(maxt))
  time=1:maxt
  for(t in 1:maxt){
    boundary[t]=min(which(Result[2:(maxc+1),t]==0))
  }
  plot(boundary~time, ylim=c(0,100), ylab="Energetic State", xlab="Time of Day (hours)", type="n", xaxt="n")
  lines(boundary~time)
  axis(1, at = c(0,126,252,378,504,630,756,882,1007), labels=c("0600", "0900", "1200", "1500", "1800", "2100", "0000", "0300", "0600"))
 # text(time[316],c((boundary[which.min(boundary)]+20),(boundary[which.min(boundary)]-20)), labels=c("Rest","Forage"), cex=0.8)
  text(time[316],c((60),(20)), labels=c("Rest","Forage"), cex=0.8)
  
  
  #Plot results
  #quartz()
  #image(x=1:maxt,y=0:maxc,t(Result),xlab="Time of day",ylab="Condition",tck=0)
  
}

### foragedielsearch = function(pred, foreff, m, c, value, metcostforage, maxt, maxc, N, prey, c0, k)

## start using the full-scale encounter functions from the malaise traps
foragedielsearch(predict.fit.pred, 0.5, 0.5, 0.000521, 1, 4, 1007, 100, predict.fit.prey, 50, 0.125)

## can then try some reasonable estimates of what hits the web, like 50% of the trap estimate
foragedielsearch(predict.fit.pred*0.5, 0.5, 0.5, 0.000521, 1, 4, 1007, 100, predict.fit.prey*0.5, 50, 0.125)

### checking effects of variation in predator and/or prey abundance
foragedielsearch (rep.int((mean(predict.fit.pred)),1007)*0.5, 0.5, 0.5, 0.000521, 1, 4, 1007, 100, predict.fit.prey*0.5, 50, 0.125)
foragedielsearch (predict.fit.pred*0.5, 0.5, 0.5, 0.000521, 1, 4, 1007, 100, rep.int((mean(predict.fit.prey))*0.5,1007), 50, 0.125)
foragedielsearch (rep.int((mean(predict.fit.pred)),1007), 0.5, 0.5, 0.000521, 1, 4, 1007, 100, rep.int((mean(predict.fit.prey)),1007), 50, 0.125)
