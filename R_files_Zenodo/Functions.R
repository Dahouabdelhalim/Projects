Directional.Preference <- function(Data, Label,PRINT)
  # Function to test for directional preference
  # Both relative preference and difference are calculated
{
  # Create relative differences based on focal male (=left male first half, right male second half)
  Left.Male          <- Data[,1]
  Right.Male         <- Data[,2]
  # Calculate directional preference using relative preference
  Rel.Diff           <- Left.Male/(Left.Male + Right.Male )
  linreg             <- lm(Data[,3]~Rel.Diff)                 # linear regression
  lin.Result         <- summary(linreg)                       # Summary stats    
  P1.Directional     <- lin.Result$coefficients[8]            # Store coefficient
  y                  <- predict(linreg, type="response")      # Calculate predicted line
  plot(Rel.Diff,Data[,3],ylim = c(-.2,1.2),xlab="Directional Preference Rel Prop",ylab="Probability",main=Label)
  points(Rel.Diff,y, col="red")
  # Calculate directional selection using difference 
  Diff               <- Left.Male - Right.Male
  linreg             <- lm(Data[,3]~Diff)                     # linear regression
  lin.Result2        <- summary(linreg)                       # Summary stats    
  P2.Directional     <- lin.Result2$coefficients[8]           # Store coefficient
  y                  <- predict(linreg, type="response")      # Calculate predicted line
  plot(Diff,Data[,3],ylim = c(-.2,1.2),xlab="Directional Preference Difference",ylab="Probability",main=Label)
  points(Diff,y, col="red")
  if(PRINT=="YES")
  {
    print("********Directional Preference Analysis********")    
    print(c("Directional selection using relative proportion"))
    print(lin.Result)    # Print summary stats
    print(c("Directional selection using difference"))
    print(lin.Result2)  #Print summary stats
  }  
  return(c(P1.Directional, P2.Directional, lin.Result$coefficients[2]))
}
############################################################################################### 
Directional.Preference01 <- function(Data)
{
  Left.Male          <- Data[,1]
  Right.Male         <- Data[,2]
  # Calculate directional preference using relative preference
  Rel.Diff           <- Left.Male/(Left.Male + Right.Male )
  logreg             <- glm(Data[,3]~Rel.Diff,family=binomial())  # logistic regression
  log.Result         <- summary(logreg)                       # Summary stats    
  P1.Directional     <- log.Result$coefficients[8]
  y                  <- predict(logreg, type="response")
  plot(Rel.Diff,Data[,3],ylim = c(-.2,1.2),xlab="Rel Prop",ylab="State",main="Directional Selection")
  points(Rel.Diff,y, col="red")
  # Calculate directional selection using difference 
  Diff               <- Left.Male - Right.Male
  logreg             <- glm(Data[,3]~Diff,family=binomial())  # logistic regression
  log.Result         <- summary(logreg)                       # Summary stats    
  P2.Directional     <- log.Result$coefficients[8]
  y                  <- predict(logreg, type="response")
  plot(Diff,Data[,3],ylim = c(-.2,1.2),xlab="Difference",ylab="State",main="Directional Selection")
  points(Diff,y, col="red")    
  return(c(P1.Directional, P2.Directional,log.Result$coefficients[2]))
}
###############################################################################################  
Stabilizing.Preference <- function(X.trial,Data, N.X, Min.Trait, Max.Trait, Label, N.pairs, Pred, PRINT)
{
  Results <- matrix(0,N.X,2)                        # Matrix to hold absolute trait value, P.log, P.lin
  # Iterate through the proposed preference values
  for( Ith.X in 1:N.X){Results[Ith.X,]<-X.Estimate(X.trial[Ith.X], Data, N.pairs, Pred)} # end of the Ith.X loop
  ############## Find the best fitting values based on the highest number of correct predictions
  Best       <- order(Results[,2], na.last=TRUE, decreasing=TRUE) #FALSE) # Find the largest correlation (note sign changed)
  Y          <- Results[Best[1],1]                                         # Female Preference at this value
  # Fit linear regression at best fit to ensure the slope is negative & significant
  Left.male  <- abs(Data[,1]-Y[1])
  Right.male <- abs(Data[,2]-Y[1])
  Pred       <- Left.male/(Left.male+Right.male) #ifelse(Left.male<Right.male,1,0)
  linreg     <- lm(Data[,3]~Pred)
  lin.Result <- summary(linreg)                 # Summary stats
  y          <- predict(linreg, type="response")
  # Export data
  Directional.Preference.metric <- Data[,1]/(Data[,1]+Data[,2])
  Export.Data                   <- matrix(0,N.X,7)  # Set matrix to take output
  n                             <- nrow(Data)
  Export.Data[1:n,1:5]          <- cbind(Data[1:n,1:3],Pred[1:n],Directional.Preference.metric[1:n])
  Export.Data[,6:7]             <- Results
  col.Names                     <- c("left","right","Observed Preference", "Predicted Preference",
                                     "Directional Preference1","Proposed X", "r2")
  colnames(Export.Data)         <- col.Names
  write.table(x=Export.Data, file = "Metric.Table.txt", row.names=FALSE,append = FALSE ) 
  # Plot distribution of preference metric to compare ED vs AD models
  hist(Pred, main="Distribution of Preference metric") 
  if(PRINT=="YES")
  {      
    print(c("Stabilizing selection at best X"))
    print(lin.Result)
  }
  Slope      <- lin.Result$coefficients[2] 
  P.Slope    <- lin.Result$coefficients[8] 
  ####################################### Plot results ############################################  
  n       <- nrow(Data)
  Tvalue  <- qt(0.025,n-2)
  r2.crit <- Tvalue^2/(n-2)
  Ymax    <- max(r2.crit, Results[,2])
  Ymax    <- Ymax*1.1
  Ymin    <- min(Results[,2])
  Ymin    <- Ymin-0.2*abs(Ymin)
  plot(Results[,1],  (Results[,2])  ,xlim=c(Min.Trait, Max.Trait),ylim=c(Ymin,Ymax), xlab="Proposed Female Preference",ylab="Signed r^2 (best fit when +)",main=Label)
  lines(Results[,1], (Results[,2]) )
  lines(c(Min.Trait, Max.Trait), c(r2.crit, r2.crit))
  return(c(Y, Slope, P.Slope)) #
} # end of function    
############################# END OF FUNCTION ###################################  
BOOTSTRAP <- function(Data,indices) #,Data, N.X, Min.Trait, Max.Trait, Label, N.pairs, Pred)
{
  #   *********************Stabilizing Preference*********************
  Data1 <- Data[indices,]
  Results <- matrix(0,N.X,2)                        # Matrix to hold absolute trait value, P.log, P.lin
  # Iterate through the absolute preference values
  for( Ith.X in 1:N.X){Results[Ith.X,]<-X.Estimate(X.trial[Ith.X], Data1, N.pairs, Pred)} # end of the Ith.X loop
  ############## Find the best fitting values based on the signed r-squared value
  Best       <- order(Results[,2], na.last=TRUE, decreasing=TRUE)   # Find the largest correct number predicted
  Y          <- Results[Best[1],1]                                  # Female Preference at this value
  return(c(Y, 0,0)) # return the last to to conform with previous function
} # end of function
############################# END OF FUNCTION ###################################  
###############################################################################################  
Stabilizing.Preference01 <- function(X.trial,Data, N.X, Min.Trait, Max.Trait, Label, N.pairs, Pred)
{
  Results <- matrix(0,N.X,2)                        # Matrix to hold absolute trait value, P.log, P.lin
  # Iterate through the absolute preference values
  for( Ith.X in 1:N.X){Results[Ith.X,]<-X.Estimate01(X.trial[Ith.X], Data, N.pairs, Pred)} # end of the Ith.X loop
  ############## Find the best fitting values based on the highest number of correct predictions
  Best       <- order(Results[,2], na.last=TRUE, decreasing=TRUE)  #FALSE) # Find the largest loglikelihood
  Y          <- Results[Best[1],1]                                         # Female Preference at this value
  #Y          <- c(Y,Results[Best[1],2])                                    # Number correct
  # Fit logistic regression at best fit to ensure the slope is negative
  Diff       <- abs(Data[,1]-Y[1]) - abs(Data[,2]-Y[1]) 
  logreg     <- glm(Data[,3]~Diff,family=binomial())     # logistic regressiont logistic to best fit
  log.Result <- summary(logreg)                          # Summary stats
  Slope      <- log.Result$coefficients[2]
  Prob       <- log.Result$coefficients[8]
  ####################################### Plot results ############################################  
  plot(Results[,1],  (Results[,2])  ,xlim=c(Min.Trait, Max.Trait), xlab="Proposed female Preference",ylab="Correct N",main=Label)
  lines(Results[,1], (Results[,2]) )
  #write.table(x=Results, file = "Results.txt", row.names=FALSE,append = FALSE )
  return(c(Y, Slope, Prob)) #
} # end of function
BOOTSTRAP <- function(Data,indices) #,Data, N.X, Min.Trait, Max.Trait, Label, N.pairs, Pred)
{
  #   *********************Stabilizing Preference*********************
  Data1 <- Data[indices,]
  Results <- matrix(0,N.X,2)                        # Matrix to hold absolute trait value, P.log, P.lin
  # Iterate through the absolute preference values
  for( Ith.X in 1:N.X){Results[Ith.X,]<-X.Estimate(X.trial[Ith.X], Data1, N.pairs, Pred)} # end of the Ith.X loop
  ############## Find the best fitting values based on the signed r-squared value
  Best       <- order(Results[,2], na.last=TRUE, decreasing=TRUE)   # Find the largest correct number predicted
  Y          <- Results[Best[1],1]                                  # Female Preference at this value
  return(c(Y, 0,0)) # return the last to to conform with previous function
} # end of function
############################# END OF FUNCTION ################################### 
X.Estimate <- function(X, Data,N.pairs,Pred)
  # Function to calculate the regression of the relative preference on the estimated relative male trait value for a given X      
{
  Left.male  <- abs(Data[,1]-X)
  Right.male <- abs(Data[,2]-X)
  Pred       <- Left.male/(Left.male+Right.male) #ifelse(Left.male<Right.male,1,0)
  linreg     <- lm(Data[,3]~Pred)
  linreg     <- summary(linreg)
  R          <- linreg$r.squared
  Slope      <- linreg$coefficients[2]
  R          <- -sign(Slope)*R
  return(c(X, R))
}
######################################END OF FUNCTION#############################################
X.Estimate01 <- function(X, Data,N.pairs,Pred)
{
  Left.male  <- abs(Data[,1]-X)
  Right.male <- abs(Data[,2]-X)
  Pred       <- ifelse(Left.male<Right.male,1,0)
  N          <- N.pairs - sum(abs(Data[,3]- Pred))
  return(c(X, N))
}
######################################END OF FUNCTION#############################################