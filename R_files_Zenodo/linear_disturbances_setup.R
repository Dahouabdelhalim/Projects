disturbance_regime_trial=list()
for (l in 1:num_trials){

disturbance_regime=matrix(1,nrow=max_time,ncol = total_patches)
for (q in 1:max_time){
disturbance_regime_temp = rbinom(n = (total_patches),size=1,prob = 1-disturbance_prob_param)



# This causes a disturbance in patches nearby the original disturbance location
  for (index in which(disturbance_regime_temp==0)){
    patches_that_are_hit = (index-disturbance_size_param):(index+disturbance_size_param)
    patches_that_are_hit=patches_that_are_hit[patches_that_are_hit>0]
    patches_that_are_hit=patches_that_are_hit[patches_that_are_hit<=length(disturbance_regime_temp)]
  disturbance_regime_temp[patches_that_are_hit]=0
  }


disturbance_regime_temp[disturbance_regime_temp==0]=1-disturbance_magnitude

disturbance_regime[q,]=disturbance_regime_temp
}


# Insert optional code here to modify where disturbances can occur
if (biased_disturbance_side == TRUE){disturbance_regime[,32:39]=1}
if (biased_disturbance_rotating == TRUE){disturbance_regime[,seq(32,47,2)]=1}
#disturbance_regime[,32:41]=1
#disturbance_regime[,seq(30,48,2)]=1
#colSums(disturbance_regime==0)
#disturbance_regime[,1:10]=1
#disturbance_regime[,30:42]=1

disturbance_regime_trial[[l]]=disturbance_regime

}


#probability_of_disturbance = 100*median(colSums(disturbance_regime==1-disturbance_magnitude)/(200))
