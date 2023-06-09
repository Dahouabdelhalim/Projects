##Load the data

##Rename data - All Adults (p=0.116)
data <- AllAdultsRelatedness

##Calculate the empirical value of W
wilco <- wilcox.test(relate ~ category,
                     data = data, paired = FALSE, alternative = "two.sided")
W_empirical <- wilco$statistic

##Calculate the distribution of Ws based on 10,000 permutations of the ranks
distrib <- numeric(10000) #empty vector for storing the Ws
data_test <- data #create a dummy data set
proportion <- as.vector(table(data$category)[2]/
                          (table(data$category)[1] + table(data$category)[2]))
  #calculate the original proportion of Within vs Between in category

for (i in 1:10000) {
  for (j in 1:length(data_test$category)) {
    data_test$category[j] <- if (runif(1) < proportion) {
      "Within"
    } else {
      "Between"
    }
  } #assign randomly a rank respecting original proportions
  
  wilco_test <- wilcox.test(relate ~ category, #compute the Wilcoxon rank sum test
                            data = data_test, paired = FALSE, alternative = "two.sided")
  distrib[i] <- wilco_test$statistic #store the W value in the vector
}


# ##Compute the proportion of Ws from the permutated test that are equal or larger than the
# ##empirical value
# sum(distrib >= W_empirical) / length(distrib) #can be interpreted as the p-value
# #the probability of committing a type I error (Manly 1997)

##Test if the empirical value is significantly part of this distribution
if (W_empirical != mean(distrib)) {
  if (W_empirical > mean(distrib)) {
    sum(distrib >= W_empirical) / length(distrib)
    #Because we are on the upper part of the distribution
  } else {
    sum(distrib <= W_empirical) / length(distrib)
    #Because we are on the lower part of the distribution
  }
} else {
  print("The empirical value is equal to the mean of the distribution!")
}
#The output value can be interpreted as the p-value
#the probability of committing a type I error (Manly 1997)

#citation()


##########################################################################################################################

##Rename data - All Offspring (p<0.001, though functionally 0)
data <- OffspringRelatedness

##Calculate the empirical value of W
wilco <- wilcox.test(relate ~ category,
                     data = data, paired = FALSE, alternative = "two.sided")
W_empirical <- wilco$statistic

##Calculate the distribution of Ws based on 10,000 permutations of the ranks
distrib <- numeric(10000) #empty vector for storing the Ws
data_test <- data #create a dummy data set
proportion <- as.vector(table(data$category)[2]/
                          (table(data$category)[1] + table(data$category)[2]))
#calculate the original proportion of Within vs Between in category

for (i in 1:10000) {
  for (j in 1:length(data_test$category)) {
    data_test$category[j] <- if (runif(1) < proportion) {
      "Within"
    } else {
      "Between"
    }
  } #assign randomly a rank respecting original proportions
  
  wilco_test <- wilcox.test(relate ~ category, #compute the Wilcoxon rank sum test
                            data = data_test, paired = FALSE, alternative = "two.sided")
  distrib[i] <- wilco_test$statistic #store the W value in the vector
}

# ##Compute the proportion of Ws from the permutated test that are equal or larger than the
# ##empirical value
# sum(distrib >= W_empirical) / length(distrib) #can be interpreted as the p-value
# #the probability of committing a type I error (Manly 1997)

##Test if the empirical value is significantly part of this distribution
if (W_empirical != mean(distrib)) {
  if (W_empirical > mean(distrib)) {
    sum(distrib >= W_empirical) / length(distrib)
    #Because we are on the upper part of the distribution
  } else {
    sum(distrib <= W_empirical) / length(distrib)
    #Because we are on the lower part of the distribution
  }
} else {
  print("The empirical value is equal to the mean of the distribution!")
}
#The output value can be interpreted as the p-value
#the probability of committing a type I error (Manly 1997)

#citation()


##########################################################################################################################

##Rename data - Just Adult Females (p=0.44)
data <- FemaleRelatedness

##Calculate the empirical value of W
wilco <- wilcox.test(relate ~ category,
                     data = data, paired = FALSE, alternative = "two.sided")
W_empirical <- wilco$statistic

##Calculate the distribution of Ws based on 10,000 permutations of the ranks
distrib <- numeric(10000) #empty vector for storing the Ws
data_test <- data #create a dummy data set
proportion <- as.vector(table(data$category)[2]/
                          (table(data$category)[1] + table(data$category)[2]))
#calculate the original proportion of Within vs Between in category

for (i in 1:10000) {
  for (j in 1:length(data_test$category)) {
    data_test$category[j] <- if (runif(1) < proportion) {
      "Within"
    } else {
      "Between"
    }
  } #assign randomly a rank respecting original proportions
  
  wilco_test <- wilcox.test(relate ~ category, #compute the Wilcoxon rank sum test
                            data = data_test, paired = FALSE, alternative = "two.sided")
  distrib[i] <- wilco_test$statistic #store the W value in the vector
}

# ##Compute the proportion of Ws from the permutated test that are equal or larger than the
# ##empirical value
# sum(distrib >= W_empirical) / length(distrib) #can be interpreted as the p-value
# #the probability of committing a type I error (Manly 1997)

##Test if the empirical value is significantly part of this distribution
if (W_empirical != mean(distrib)) {
  if (W_empirical > mean(distrib)) {
    sum(distrib >= W_empirical) / length(distrib)
    #Because we are on the upper part of the distribution
  } else {
    sum(distrib <= W_empirical) / length(distrib)
    #Because we are on the lower part of the distribution
  }
} else {
  print("The empirical value is equal to the mean of the distribution!")
}
#The output value can be interpreted as the p-value
#the probability of committing a type I error (Manly 1997)

#citation()

##########################################################################################################################

##Rename data - Just Adult Males (p=0.37)
data <- MaleRelatedness

##Calculate the empirical value of W
wilco <- wilcox.test(relate ~ category,
                     data = data, paired = FALSE, alternative = "two.sided")
W_empirical <- wilco$statistic

##Calculate the distribution of Ws based on 10,000 permutations of the ranks
distrib <- numeric(10000) #empty vector for storing the Ws
data_test <- data #create a dummy data set
proportion <- as.vector(table(data$category)[2]/
                          (table(data$category)[1] + table(data$category)[2]))
#calculate the original proportion of Within vs Between in category

for (i in 1:10000) {
  for (j in 1:length(data_test$category)) {
    data_test$category[j] <- if (runif(1) < proportion) {
      "Within"
    } else {
      "Between"
    }
  } #assign randomly a rank respecting original proportions
  
  wilco_test <- wilcox.test(relate ~ category, #compute the Wilcoxon rank sum test
                            data = data_test, paired = FALSE, alternative = "two.sided")
  distrib[i] <- wilco_test$statistic #store the W value in the vector
}

# ##Compute the proportion of Ws from the permutated test that are equal or larger than the
# ##empirical value
# sum(distrib >= W_empirical) / length(distrib) #can be interpreted as the p-value
# #the probability of committing a type I error (Manly 1997)

##Test if the empirical value is significantly part of this distribution
if (W_empirical != mean(distrib)) {
  if (W_empirical > mean(distrib)) {
    sum(distrib >= W_empirical) / length(distrib)
    #Because we are on the upper part of the distribution
  } else {
    sum(distrib <= W_empirical) / length(distrib)
    #Because we are on the lower part of the distribution
  }
} else {
  print("The empirical value is equal to the mean of the distribution!")
}
#The output value can be interpreted as the p-value
#the probability of committing a type I error (Manly 1997)

#citation()