library(dagitty)
library(rethinking)

dag1 <- dagitty("dag{
  
  socialintegration -> Survival
  
  rank -> socialintegration
  rank -> Survival
  
  age -> Survival
  age -> socialintegration
  age -> rank
  age -> mother
  age -> daughters
  
  mother -> Survival
  
  daughters -> Survival
  
  groupsize -> socialintegration
  groupsize -> Survival
  
  
}" )


coordinates( dag1 ) <- list( x=c(Survival=0.5, socialintegration=0.5, rank = 0.25, age = 0.25, mother = 0.5, daughters = 0.25, groupsize = 0.75 ) ,
                             y=c(Survival=0.6, socialintegration=0.2, rank = 0.2, age = 0.6, mother = 1, daughters = 1, groupsize = 0.6 ) )


drawdag( dag1 )

# predictors included in the final models
adjustmentSets(dag1,exposure="socialintegration" , outcome="Survival" )


