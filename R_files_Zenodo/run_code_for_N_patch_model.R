


  
# This line specifies which reserve combinations someone wants to run. It is best to leave out reserves on the ends as they would have high mortality sending recruits outside the reserve
configs_to_try = t(combn((extra_patches+3):(extra_patches+num_patches-2), num_reserves))
## Optional line to stick one reserve in the middle
#configs_to_try = configs_to_try[1:(which(configs_to_try[,1]==(extra_patches+3+1))[1]-1),]
configs_to_try = configs_to_try[configs_to_try[,round(num_reserves/2)]==round(mean(range(configs_to_try))),]
###

result=matrix(0,nrow=nrow(configs_to_try),ncol=(num_reserves+1))

for (index1 in 1:nrow(configs_to_try)){
  x1=rep(1,total_patches);x1[as.numeric(configs_to_try[index1,])]=0
   result[index1,]=c(which(x1==0),evalFunc(x1))
 
 #result=matrix(0,nrow=5,ncol=(num_reserves+1))
# for (index1 in 1:5){
#   x1=rep(1,total_patches);x1[c(15,15+index1)]=0
#   result[index1,]=c(which(x1==0),evalFunc(x1))
#   print(result)
 }
#}
  ranked_configs = result[order(result[,num_reserves+1]),1:(num_reserves+1)]
  num_near_top = which(ranked_configs[,(num_reserves+1)] <= ranked_configs[1,(num_reserves+1)]-0.001*ranked_configs[1,(num_reserves+1)])
  
  
  # Check to ensure optimal spacing works
    #dude=ranked_configs[num_near_top,1:num_reserves]
    #spacing <- apply(dude,1,FUN=diff) -1
    #print(min(colMeans(spacing)))
  # End chunk 

  #print(ranked_configs[1:5,])
    if (is.na(ranked_configs[1,1+num_reserves]) | ranked_configs[1,1+num_reserves] == 0){
      ranked_configs = matrix(NA,nrow=nrow(ranked_configs),ncol=ncol(ranked_configs))
    }
    optimal_reserves= ranked_configs[1,1:num_reserves]
#}
   # optimal_spacing[value] = diff(optimal_reserves)#mean(mean(apply(optimal_reserves,1,FUN=diff)))
   # mean_objective_value[value]=-ranked_configs[1,1+num_reserves]#mean(-ranked_configs[num_near_top,(num_reserves+1)])/(max_time*num_trials)
    
    #optimal_spacing[value] = diff(optimal_reserves[1,])
    #print(optimal_spacing[value])
    #print(length(num_near_top))
  
    #test=rep(1,num_patches)
    #test[optimal_reserves[value,]]=0
    #print(evalFunc(test))
  
