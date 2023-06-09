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

foragedielround2 = function(pred, foreff, m, c, value, metcostforage, maxt, maxc, prey, c0, k){                     
  
  ##### begin by creating the fitness function fx that can be used to initialize reward matrix
  Cs=0:maxc   ### individuals can have reserves of 0, which implies death by starvation.
  fx= 1/(1+(exp(-k*(Cs-c0))))  ### sigmoidal fitness function
  #plot(fx~Cs)
  
  ### optional plotting to show encounter functions and model predictions as a 2-panel plot.
  
  #quartz(width=4, height=6.5)
  #par(mfcol=c(2,1), mar=c(5,4,1.5,2))    ### set up a plot window to hold the encounter functions and predictions as panels
  
  #plot(pred*pforage~c(1:maxt), ylim=c(0,0.1), ylab="Probability of Event", xlab="Time of Day (hours)", type="n", xaxt="n")
 # lines(pred*pforage~c(1:maxt), col="red")    ### lines for predator encounter probability
 # lines(prey*foreff~c(1:maxt), col="blue")    ### lines for prey encounter probability
 # legend("topright", legend = c("Predation", "Prey Capture"), lty=1, col=c("red", "blue"), cex=0.8)
  ### axis for time-of-day values is hard-coded here... easy enough to set up a sequence to generate the ticks (i.e., at = ) for whatever time intervals (i.e., labels = ) you like though.
#  axis(1, at = c(0,127,254,381,508,635,762,889,1016), labels=c("0600", "0900", "1200", "1500", "1800", "2100", "0000", "0300", "0600"))
  
  Reward = matrix(0,nrow = maxc+1,ncol=maxt+1)   ### create empty reward matrix (holds expected fitness values for each combination of energetic reserves (rows) and time step t (column), assuming optimal behavior from that t forward)
  Reward[,maxt+1]=(fx)      ### set last column to fitness value (fx) for each level of reserves at T 
  Result=matrix(0,nrow = maxc+1,ncol=maxt)      ### create empty result matrix, hold optimal forage decision (0,1) for each combination of energetic reserves and t
  
  
  for (t in maxt:1){     ### work backward from the end of the time interval.
    
    ## Calculate prob of predation event due to foraging  at t ##
    
    dforage = m * pred[t]
    
    ### Calculate probability of not dying and gaining and not dying and losing for each choice###
    ### Note, can exclude scenarios where predation occurs, because predation corresponds to 
    ### expected fitness of 0.
    
    p.forage.up=(1-dforage)*foreff   ### individual attacks, is not predated, and captures a prey 
    p.forage.lose=(1-dforage)*(1-foreff)  ### individual attacks, is not predated, and fails to capture a prey
    
    p.rest.lose=(1)      ### Individual rests, has no chance of predation, probability of gaining energy is 0, probability of losing is 1.
    
    ### Use above probabilities to calculate exp. terminal reward for each choice for each possible level of energetic reserves ####
    
    for(i in 2:maxc){      ### loop over values of energetic reserves. start with i = 2, because i = 1 corresponds to energy reserves = 0 --> dead individual, no decision to make.
      
      ### calculate expected state for each possibility ###
      state.no.enc = i - c   ### no prey are encountered, lose energy
      state.yes.rest = i - c   ### prey are encountered but not attacked, lose energy
      state.yes.attack.success = i + value - (metcostforage * c)   ### prey encountered, attacked, and captured, pay added metabolic cost of foraging
      state.yes.attack.failure = i - (metcostforage * c)     ### prey encountered, attacked, but failed to capture, still pay added metabolic cost of foraging
      
      ### Need to find Expected Future Fitness (Reward) for each of the possibilities with interpolation
      
      ### first, find the next lowest discrete value of i for each decision
      li.state.no.enc = floor(state.no.enc)
      li.state.yes.rest = floor(state.yes.rest)
      li.state.yes.attack.success = floor(state.yes.attack.success)
      li.state.yes.attack.failure = floor(state.yes.attack.failure)
      
      #### next, use interpolation to find expected reward for each state
      
      ### first, for no encounter
      if (state.no.enc == li.state.no.enc){  ### is the future state if no enc. a discrete value?
        Fit.no = Reward[state.no.enc,t+1]  ## then expected fitness from this possible outcome is equal to that of an individual in the state.no.enc level of i in the next time step, t+1.
      }
      else{  ### if not, interpolate using the future state and the two nearest discrete values.
        Fit.no = ((state.no.enc-li.state.no.enc)*Reward[li.state.no.enc+1,t+1] + (((li.state.no.enc+1)-state.no.enc)*Reward[li.state.no.enc,t+1])) / ((li.state.no.enc+1) - li.state.no.enc)
      }
      
      ### for encounter and resting
      
      if (state.yes.rest == li.state.yes.rest){
        Fit.yes.rest = Reward[state.yes.rest,t+1]
      }
      else{
        Fit.yes.rest = ((state.yes.rest-li.state.yes.rest)*Reward[li.state.yes.rest+1,t+1] + (((li.state.yes.rest+1)-state.yes.rest)*Reward[li.state.yes.rest,t+1])) / ((li.state.yes.rest+1) - li.state.yes.rest)
      }
      
      ### for encounter and attack and success
      
      if (state.yes.attack.success == li.state.yes.attack.success){
        Fit.yes.attack.success = Reward[state.yes.attack.success,t+1]
      }
      else{
        Fit.yes.attack.success = ((state.yes.attack.success-li.state.yes.attack.success)*Reward[li.state.yes.attack.success+1,t+1] + (((li.state.yes.attack.success+1)-state.yes.attack.success)*Reward[li.state.yes.attack.success,t+1])) / ((li.state.yes.attack.success+1) - li.state.yes.attack.success)
      }
      
      ### for encounter and attack and failure
      
      if (state.yes.attack.failure == li.state.yes.attack.failure){
        Fit.yes.attack.failure = Reward[state.yes.attack.failure,t+1]
      }
      else{
        Fit.yes.attack.failure = ((state.yes.attack.failure-li.state.yes.attack.failure)*Reward[li.state.yes.attack.failure+1,t+1] + (((li.state.yes.attack.failure+1)-state.yes.attack.failure)*Reward[li.state.yes.attack.failure,t+1])) / ((li.state.yes.attack.failure+1) - li.state.yes.attack.failure)
      }
      
      
      #### set expected fitness for i and t (Reward[i,t]) according to the likelihood of prey encounter and the fitness payoffs of whichever decision maximizes future fitness if prey is encountered.
      
      #Reward[i,t] = ((1 - prey[t]) * p.rest.lose * Fit.no) + (prey[t]) * ifelse ((p.rest.lose * Fit.yes.rest) > ((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure)), (p.rest.lose * Fit.yes.rest), ((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure)))
      
      ### some rounding can be included here to prevent artifacts that arise under certain parameter combinations where differences in expected fitness between options become very small toward the end of the time period.
      
    Reward[i,t] = ((1 - prey[t]) * p.rest.lose * Fit.no) + (prey[t]) * ifelse ((p.rest.lose * Fit.yes.rest) + 0.0002> ((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure)), (p.rest.lose * Fit.yes.rest), ((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure)))
      
      ### now, record what the optimal decision was in the Result matrix (1 = forage, 0 = rest)
      
      #Result[i,t] = ifelse(((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure)) > (p.rest.lose * Fit.yes.rest), 1, 0) 
      
      ### again, some rounding can be included. Be sure that both the reward and result matrix calculations are rounding equivalently.
      Result[i,t] = ifelse(((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure)) +0.0002 > (p.rest.lose * Fit.yes.rest), 1, 0) 
      
    }    ### this ends the loop over i, still in the t-loop
    
    
    ### individuals in top condition cannot gain energy reserves, but can lose them ###
    
    ### need new state change possibilities
    
    state.no.enc = (maxc+1) - c    ### loses c from max condition
    state.yes.rest = (maxc+1) - c   ### loses c from max condition
    state.yes.attack.success = maxc+1   ### success covers costs of foraging, stays in max condition
    state.yes.attack.failure = ((maxc+1) - (metcostforage * c)) ### loses cost of foraging from max condition
    
    
    ### Need to find Expected Future Fitness (Reward) for each of the options using interpolation
    
    ### first, find the next lowest discrete value of i for each decision
    li.state.no.enc = floor(state.no.enc)
    li.state.yes.rest = floor(state.yes.rest)
    li.state.yes.attack.success = floor(state.yes.attack.success)
    li.state.yes.attack.failure = floor(state.yes.attack.failure)
    
    
    #### next, use interpolation to find expected reward for each state
    
    ### first, for no encounter
    if (state.no.enc == li.state.no.enc){
      Fit.no = Reward[state.no.enc,t+1]
    }
    else{
      Fit.no = ((state.no.enc-li.state.no.enc)*Reward[li.state.no.enc+1,t+1] + (((li.state.no.enc+1)-state.no.enc)*Reward[li.state.no.enc,t+1])) / ((li.state.no.enc+1) - li.state.no.enc)
    }
    
    ### for encounter and resting
    
    if (state.yes.rest == li.state.yes.rest){
      Fit.yes.rest = Reward[state.yes.rest,t+1]
    }
    else{
      Fit.yes.rest = ((state.yes.rest-li.state.yes.rest)*Reward[li.state.yes.rest+1,t+1] + (((li.state.yes.rest+1)-state.yes.rest)*Reward[li.state.yes.rest,t+1])) / ((li.state.yes.rest+1) - li.state.yes.rest)
    }
    
    ### for encounter and attack and success
    
    if (state.yes.attack.success == li.state.yes.attack.success){
      Fit.yes.attack.success = Reward[state.yes.attack.success,t+1]
    }
    else{
      Fit.yes.attack.success = ((state.yes.attack.success-li.state.yes.attack.success)*Reward[li.state.yes.attack.success+1,t+1] + (((li.state.yes.attack.success+1)-state.yes.attack.success)*Reward[li.state.yes.attack.success,t+1])) / ((li.state.yes.attack.success+1) - li.state.yes.attack.success)
    }
    
    ### for encounter and attack and failure
    
    if (state.yes.attack.failure == li.state.yes.attack.failure){
      Fit.yes.attack.failure = Reward[state.yes.attack.failure,t+1]
    }
    else{
      Fit.yes.attack.failure = ((state.yes.attack.failure-li.state.yes.attack.failure)*Reward[li.state.yes.attack.failure+1,t+1] + (((li.state.yes.attack.failure+1)-state.yes.attack.failure)*Reward[li.state.yes.attack.failure,t+1])) / ((li.state.yes.attack.failure+1) - li.state.yes.attack.failure)
    }
    
    
    #### now, set fitness for i and t (Reward) according to the likelihood of prey encounter and the fitness payoffs of whichever decision maximizes future fitness if prey is encountered.
    
    #Reward[i,t] = ((1 - prey[t]) * p.rest.lose * Fit.no) + (prey[t]) * ifelse ((p.rest.lose * Fit.yes.rest) > ((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure)), (p.rest.lose * Fit.yes.rest), ((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure)))
    
    # can round again if needed
    Reward[i,t] = ((1 - prey[t]) * p.rest.lose * Fit.no) + (prey[t]) * ifelse ((p.rest.lose * Fit.yes.rest) + 0.0002> ((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure)), (p.rest.lose * Fit.yes.rest), ((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure)))
    
    ### now, record what the optimal decision was in the Result matrix (1 = forage, 0 = rest)
    
    #Result[i,t] = ifelse(((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure))> (p.rest.lose * Fit.yes.rest), 1, 0)
    
    Result[i,t] = ifelse(((p.forage.up * Fit.yes.attack.success) + (p.forage.lose * Fit.yes.attack.failure))+ 0.0002> (p.rest.lose * Fit.yes.rest), 1, 0)
    
    
  }   ### this ends the t-loop.
  
  
  boundary=vector("numeric", length=length(maxt))
  time=1:maxt
  for(t in 1:maxt){
    boundary[t]=min(which(Result[2:(maxc+1),t]==0))
  }
  plot(boundary~time, ylim=c(0,100), ylab="Energetic State", xlab="Time of Day (hours)", type="n", xaxt="n")
  lines(boundary~time)
  axis(1, at = c(0,126,252,378,504,630,756,882,1007), labels=c("0600", "0900", "1200", "1500", "1800", "2100", "0000", "0300", "0600"))
  text(time[which.min(boundary)],c((boundary[which.min(boundary)]+20),(boundary[which.min(boundary)]-20)), labels=c("Rest","Forage"), cex=0.8)
  
  
}  




debug(foragedielround2)
### foragediel = function(pred, foreff, m, c, value, metcostforage, maxt, maxc, prey, c0, k)

## start using the full-scale encounter functions from the malaise traps
foragedielround2(predict.fit.pred, 0.5, 0.5, 0.000521, 1, 4, 1007, 100, predict.fit.prey, 50, 0.125) 

## can then consider more reasonable estimates of what hits the web, like %50 of those estimates
### just hard-coding some scalars alongside the functions as input
### these represent the values at which each function was

foragedielround2(predict.fit.pred*0.5, 0.5, 0.5, 0.000521, 1, 4, 1007, 100, predict.fit.pred*0.5, 50, 0.125) 

### quick and dirty way to visualize effects of variation predator or prey abundances
foragedielround2 (rep.int((mean(predict.fit.pred))*0.5,1007), 0.5, 0.5, 0.000521, 1, 4, 1007, 100, predict.fit.pred*0.5, 50, 0.125)
foragedielround2 (predict.fit.pred*0.5, 0.5, 0.5, 0.000521, 1, 4, 1007, 100, rep.int((mean(predict.fit.pred))*0.5,1007), 50, 0.125)
foragedielround2 (rep.int((mean(predict.fit.pred))*0.5,1016), 0.5, 0.5, 0.000521, 1, 4, 1007, 100, rep.int((mean(predict.fit.pred))*0.5,1016), 50, 0.125)

