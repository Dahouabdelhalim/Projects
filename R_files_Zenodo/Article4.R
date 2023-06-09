# Read data file
dat1<-read.csv("house1.csv",header = TRUE)
dat2<-read.csv("all_trends_3.csv",header = TRUE)
dat3<-read.csv("量表.csv",header = TRUE)
price <- dat1[1:131,2];
data <- dat2[51:181,];


plot(dat3$V.Mean.Sum,dat3$A.Mean.Sum)
plot(dat3$D.Mean.Sum,dat3$A.Mean.Sum)
plot(dat3$V.Mean.Sum,dat3$D.Mean.Sum)
cor.test(dat3$V.Mean.Sum,dat3$D.Mean.Sum)

data_all <- dat3[ , c(3, 6, 9)]
data_final <- data.frame(data_all,dat3$Word)
data_final
#rownames(data_all) <-  dat3$Word 

head(data_final[,2:4])
km <- kmeans(data_all[,1:2],3)
plot(data_all[c("data_final$V.Mean.Sum", "data_final$A.Mean.Sum")], col = km$cluster, pch = as.integer(data_final$V.Mean.Sum))
####################
### main program ###
####################
sum_r = vector()
name_r = vector()
for (word_seq in 2:ncol(data)) {
  cor <- cor.test(dat1[ , 2], data[ , word_seq])
  if (cor$p.value > 1e-20)  print(paste("P value of", colnames(data)[word_seq], "> 0.001"))
  else{
    print(paste("P value of", colnames(data)[word_seq], "< 0.001"));
    
    list_r <- main_function(price, data, word_seq)
    name_r <- c(name_r, main_function(price, data, word_seq)[[1]])
    sum_r <- c(sum_r, main_function(price, data, word_seq)[[2]])
    print(paste("Sum r =", sum_r[length(sum_r)]))
  }
}


names(sum_r) <- name_r;

write.csv(sum_r,file = "Results_001.csv")
getwd()
dat2<-read.csv("Results_001.csv",header = TRUE)
plot(dat2$D.Mean.Sum,dat2$profit)
plot(dat2$A.Mean.Sum,dat2$profit)
plot(dat2$V.Mean.Sum,dat2$profit)
cor.test(dat2$A.Mean.Sum,dat2$profit)
cor.test(dat2$V.Mean.Sum,dat2$profit)
cor.test(dat2$D.Mean.Sum,dat2$profit)
wilcox.test(dat2$A.Mean.Sum,dat2$profit)
wilcox.test(dat2$D.Mean.Sum,dat2$profit)
wilcox.test(dat2$V.Mean.Sum,dat2$profit)

####################
###     end      ###
####################

##main function
main_function <- function(price, data, keyword_seq){
  keyword = data[ , keyword_seq]
  r<-rep(1,nrow(data))
  deltat <- 6 # Half a Year
  for(i in 1:nrow(data)){
    if (i>1) {r[i] <- r[i-1]}
    if (i>deltat){ 
      if(i<nrow(data)){
        now <- keyword[i]
        previous <- 0  
        for (t in 1:deltat){
          previous <- (previous+keyword[i-t])
        }
        previous <- (previous/deltat)
        value<-(now-previous)
        index_now <- price[i]
        index_next <- price[i+1]
        index_r <- (index_next/index_now)
        if (value>=0) {r[i]<-r[i-1]/index_r}
        if (value<0)  {r[i]<-r[i-1]*index_r}
      }
    }  
  }
  
  return(list(name = colnames(data)[keyword_seq], sum = sum (r)))
  
  
}
####end

#Col names
#colnames(data)<-c("Time","Search","Price","Unknown","Sales","Real-Estate","Short_Sale")
keyword<-data$Garden
keyword
cor.test(data$house.price,data$Price)
cor.test(data$Unemployment,data$Price)
cor.test(data$Cash,data$Price)
cor.test(data$Wedding,data$Price)
cor.test(data$Loss,data$Price)
cor.test(data$Luck,data$Price)
cor.test(data$Debt,data$Price)
cor.test(data$Divorce,data$Price)
cor.test(data$Credit,data$Price)
cor.test(data$Salary,data$Price)
cor.test(data$Tourism,data$Price)
cor.test(data$Money,data$Price)
cor.test(data$Town,data$Price)
cor.test(data$Unfair,data$Price)
cor.test(data$House,data$Price)
cor.test(data$Profit,data$Price)
cor.test(data$Garden,data$Price)
#cor.test(data$Search,data$Listing)
#cor.test(data$Short.Sale,data$Listing)
#cor.test(data$Debt,data$Listing)
#cor.test(data$Marriage,data$Listing)

#cor.test(data$Bank.Rate,data$Listing)
#cor.test(data$Car,data$Listing)
#cor.test(data$Home.inspection,data$Listing)
#cor.test(data$Home.staging,data$Listing)
#cor.test(data$Credit,data$Listing)
#cor.test(data$Lottery,data$Listing)
#cor.test(data$Job.Interview,data$Listing)
# Init trading account

r<-rep(1,nrow(data))
deltat <- 6 # Half a Year
for(i in 1:nrow(data)){
  if (i>1) {r[i] <- r[i-1]}
  if (i>deltat){ 
    if(i<nrow(data)){
      now <- keyword[i]
      previous <- 0  
      for (t in 1:deltat){
        previous <- (previous+keyword[i-t])
      }
      previous <- (previous/deltat)
      value<-(now-previous)
      index_now <- data$Price[i]
      index_next <- data$Price[i+1]
      index_r <- (index_next/index_now)
      if (value>=0) {r[i]<-r[i-1]/index_r}
      if (value<0)  {r[i]<-r[i-1]*index_r}
    }
  }  
}

sum (r)

plot(100*(r-1),type="l",col="blue",
     xlab="Time, t [Months]", ylab="Profit and Loss [%]")

plot(x=data$Time,y=data$Price)
# buy & hold
invest <- rep(0,nrow(data)) 

 for (i in 1:nrow(data)){
   if (i>1) {invest[i] <- invest[i-1]}
    if(i<nrow(data)){
     index_now <- data$Price[i]
     index_next <- data$Price[i+1]
     invest[i] <- ((index_next-data$Price[1])/data$Price[1])
                    }
                        }
plot(100*(invest),type="l",col="red",
     xlab="Time, t [Months]", ylab="Profit and Loss [%]")

(data$Listing[109]-data$Listing[1])/data$Listing[1]

# Print result
r
invest
t  <- 1:131
n1 <- 100*invest
n2 <- 100*(r-1) 
n1
n2
sum(n1)
sum(n2)
data1 <- data.frame(t,n1,n2)

library(ggplot2)
p = ggplot() + 
     geom_line(data = data1, aes(x = t, y = n1), color = "blue") +
    geom_line(data = data1, aes(x = t, y = n2), color = "red") +
    xlab('Time, t [Months]') +
    ylab('Profit and Loss [%]')
p
sum(r)
sum(invest)

r
dat2<-read.csv("results.csv",header = TRUE)
plot(dat2$A.Mean.Sum,dat2$Profit)
cor.test(dat2$Profit,dat2$A.Mean.Sum)
 
data_01<-read.csv("Article_04.csv",header = TRUE)
data_01
plot(data_01$A.Mean.Sum,data_01$Profit)
plot(data_01$V.Mean.Sum,data_01$Profit)
plot(data_01$D.Mean.Sum,data_01$Profit)
cor.test(data_01$A.Mean.Sum,data_01$Profit)
cor.test(data_01$V.Mean.Sum,data_01$Profit)
cor.test(data_01$D.Mean.Sum,data_01$Profit)
