#' Rank reliability estimate
#'
#' Estimates the  reliability of rank using dominance data
#' @param rank.df A dataframe containing the identity in the first coloumn, with their rank (or rank score) in the second column. Note: additional rank estimates can be added to subsequent columns, e.g., id, elo_ranking, david_score_ranks,  ... etc.
#' @param testing.df A dataframe of dominance interactions: winner, loser, datetime
#' @examples 
#' @import dplyr
#' @export
reliability_check<-function(rank.df,testing.df){
  
  rank.df$ID<- as.factor(rank.df$ID)
  rank.df$Winner <- factor(rank.df$ID, levels=levels(testing.df$winner))
  rank.df$Loser <- factor(rank.df$ID, levels=levels(testing.df$loser))
  rank.df <- rank.df %>% select( ncol(rank.df), ncol(rank.df)-1,everything())
  
  
  is.match.true<- vector()
  result.reliability<- vector(length = ncol(rank.df))
  total <- nrow(testing.df)
  global.reliability.df<- data.frame(Method=as.character("method"), Reliability= as.numeric(100), low.ci.boot=as.numeric(0.01),high.ci.boot=as.numeric(0.99))
  
  for (j in 4:ncol(rank.df)){
    
    for(i in 1:total){
      
      winner_ID_i <- testing.df$winner[i]
      loser_ID_i <- testing.df$loser[i]
      match_winner <- which(unlist(rank.df$Winner) == winner_ID_i)
      match_loser<- which (unlist(rank.df$Loser)== loser_ID_i)
      
      if (testing.df$draw[i] == "TRUE"  ){
        
        if (length(rank.df[match_loser,j])==0 | length(rank.df[match_winner,j])==0){
          is.match.true[i] <- NA
        } else{
          
          if (as.numeric(as.character(rank.df[match_winner,j]))== as.numeric(as.character(rank.df[match_loser,j]))){
            is.match.true[i] <- 1
            
          } else {
            is.match.true[i] <- 0
            
          }
        }
      }
      if (testing.df$draw[i] == "FALSE"){
        
        if (length(rank.df[match_loser,j])==0 | length(rank.df[match_winner,j])==0 ){
          is.match.true[i] <- NA
          
        } else if (as.numeric(as.character(rank.df[match_winner,j])) <=as.numeric(as.character(rank.df[match_loser,j]))){
          is.match.true[i] <- 1
          
        } else if(as.numeric(as.character(rank.df[match_winner,j])) >=as.numeric(as.character(rank.df[match_loser,j]))){
          is.match.true[i] <- 0
          
        }
      }
    }
    
    # Store match in df
    #testing.df<-add_column(testing.df, match=is.match.true)
    testing.df$match<-is.match.true
    names(testing.df)[ncol(testing.df)] <- paste("Match", colnames(rank.df)[j], sep = "_")
    
    # Calculate global % of reliability and the CIs
    global.reliability<- sum(na.omit(testing.df[ncol(testing.df)]))/nrow(na.omit(testing.df[ncol(testing.df)]))
    
    boot.nb <- 1000
    sample.size <- nrow(testing.df)
    is.match.true<-vector()
    
    #where to temporarily store bootstrapped samples from the original data
    boot.relia <- vector()
    
    for (n in 1:boot.nb){
      
      #random sample with replacement from the observed interactions
      random.rows <- sample(1:sample.size, sample.size,replace=T)
      Boot.List <- testing.df[random.rows,]
      #calculate and store the measure calculated from the bootstrapped sample
      boot.relia[length(boot.relia)+1] <- sum(na.omit(Boot.List[ncol(Boot.List)]))/nrow(na.omit(Boot.List[ncol(Boot.List)]))
      
    }
    CI.result.boot<-quantile(boot.relia, probs = c(0.0275,0.975),na.rm=T)
    
    # Store % and CI
    temporary.df<-data.frame(Method=colnames(testing.df)[ncol(testing.df)],Reliability=global.reliability,low.ci.boot=CI.result.boot[1],high.ci.boot=CI.result.boot[2])
    global.reliability.df<- rbind(global.reliability.df, temporary.df)
    
  } 
  global.reliability.df<-global.reliability.df[-1,]
  rownames(global.reliability.df)<-c()
  testing.clean<-testing.df
  
  return(list(global.reliability.df,testing.clean))
}




#' Create a testing dataset
#'
#' This function creates a testing dataset that can then be used by the function 'reliability_check'.
#' @param testing.raw A dataframe containing the observed dominance data with the first and second coloumns containing the individual, a third column with the outcomes, and a fourth with a date (as ymd).
#' @param ties Wether to remove or keep ties/draws in the testing dataset.
#' @examples 
#' @importFrom lubridate ymd
#' @export
set_testing<-function (testing.raw, ties="remove"){
  
  total<- nrow(testing.raw)
  
  for ( i in 1:total){
    
    from_ID_i <- testing.raw$from[i]
    to_ID_i <- testing.raw$to[i]
    #winner[i]<- 0
    #loser[i]<-0
    
    if (ties=="remove"){
      if (testing.raw$result[i]== "1"){
        testing.raw$winner[i]<- as.character(from_ID_i)
        testing.raw$loser[i]<- as.character(to_ID_i)
        testing.raw$draw[i]<- FALSE
        
      }
      
      if(testing.raw$result[i]=="2" ){
        
        testing.raw$winner[i]<- as.character(to_ID_i)
        testing.raw$loser[i]<- as.character(from_ID_i)
        testing.raw$draw[i]<- FALSE
        
      }
      
      if(testing.raw$result[i]=="3"){
        ## Note that result =3 are the draws, here it doesnt matter who is the winner/loser in the TESTING dataset, we are just assesing if, when there is a draw, both individuals have the same rank... thats why I arbitraly attribue the agg column as the winners.
        testing.raw$winner[i]<- NA
        testing.raw$loser[i]<-NA
        testing.raw$testing.raw$draw[i]<- NA
        
      }
    }
    if (ties=="keep"){
      
      if (testing.raw$result[i]== "1"){
        testing.raw$winner[i]<- as.character(from_ID_i)
        testing.raw$loser[i]<- as.character(to_ID_i)
        testing.raw$draw[i]<- FALSE
        
      }
      
      if(testing.raw$result[i]=="2" ){
        
        testing.raw$winner[i]<- as.character(to_ID_i)
        testing.raw$loser[i]<- as.character(from_ID_i)
        testing.raw$draw[i]<- FALSE
        
      }
      if(testing.raw$result[i]=="3"){
        ## Note that result =3 are the draws, here it doesnt matter who is the winner/loser in the TESTING dataset, we are just assesing if, when there is a draw, both individuals have the same rank... thats why I arbitraly attribue the agg column as the winners.
        
        testing.raw$winner[i]<- as.character(to_ID_i)
        testing.raw$loser[i]<- as.character(from_ID_i)
        testing.raw$draw[i]<- TRUE
        
      }
      
    }
    
  }
  testing.data<- testing.raw[,c("date","winner","loser","result","draw")]
  testing.data$date<- lubridate::ymd(testing.data$date)
  testing.data$winner<- as.factor(testing.data$winner)
  testing.data$loser<- as.factor(testing.data$loser)
  
  return(testing.data)
}


#' A plotting function for output from reliability_check()
#'
#' This function creates a plot using the outputs from reliability_check()
#' @param global.reliability A dataframe containing the estimates of overall reliability from the reliability_check() function.
#' @param df.reliability A dataframe containing the dyad level estimates of reliability from the reliability_check() function.
#' @examples 
#' @import ggplot2 reshape2 scales
#' @export
reliability_plot <- function (global.reliability,df.reliability,method="loess" ){
  
  ## plot bootstrap and CIs 
  ci.plot<-ggplot(global.reliability, aes(x=Method, y=Reliability, color=Method)) + 
    geom_pointrange(aes(ymin=low.ci.boot, ymax=high.ci.boot))
  
  ## Jitter plot
  reliability.melt<- melt(df.reliability,id.vars = c("date","winner","loser", "result", "draw"))
  colnames(reliability.melt)[6:7]<-c("Method", "Match")
  jitter.plot<-ggplot(reliability.melt, aes(date, Match, colour = Method)) + geom_jitter() + stat_smooth(method = method)+ scale_x_date(labels = date_format("%b %Y"), date_breaks = "1 month")
  
  return(list(ci.plot,jitter.plot))   
}

