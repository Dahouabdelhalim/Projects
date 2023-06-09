# get data

setwd("C:/Users/mleib1/Desktop/nvm")

d <- read.table("nvm_adam_1995_2017.dsv", sep = "|", header = T)

######## looking at the ones with status = sold (coded as "5")

d <- d[d$STATUS == 5, ]

### coding for year of transaction

d$DATUM_AFMELDING <- as.factor(d$DATUM_AFMELDING)

levels(d$DATUM_AFMELDING)[as.numeric(d$DATUM_AFMELDING)]

year <- as.numeric(unlist(strsplit(levels(d$DATUM_AFMELDING)[as.numeric(d$DATUM_AFMELDING)], split = "-"))[seq(3, 3*148500, 3)])
length(year)

table(year)
d$yearoftran <- year

# code distrcite 

d$district <- substr(d$BUURT, 1, 1)


###### calculate selling minus asking price

d$sale.minus.ask <- d$TRANSACTIEPRIJS - d$LAATSTVRKOOPPR

# code whether sell is above ask (1) or no (0)

d$above <- 0
d$above[d$sale.minus.ask > 0] <- 1

#### restricting to year = 2017

d17 <- d[d$yearoftran == 17, ]

### from http://statline.cbs.nl/Statweb/publication/?VW=T&DM=SLEN&PA=83913ENG&D1=0%2c3&D2=17&D3=0-3%2c5-8%2c10-13%2c15-18%2c20-23%2c25-28%2c30-33%2c35-38%2c40-43%2c45-48%2c50-53%2c55-58%2c60-63%2c65-68%2c70-73%2c75-78%2c80-83%2c85-88%2c90-93%2c95-98%2c100-103%2c105-108%2c110-113%2c115-117&HD=181128-1454&LA=EN&HDR=T%2cG1&STB=G2
# in 2017 there were 2740 + 2703 + 3002 + 2841 = 11286 houses sold
# the data (d17 = all the houses in 2017) contains 8278 listings = 73.34% of all houses sold in Amsterdam. 
## since the exact number of sold housescan change based on different estimations, in the paper we state the data set covers ~70% of the real-estate transactions

## removing 3 values where the asking price is -1 

d17 <- d17[d17$LAATSTVRKOOPPR > 0,  ]


### coding for age of property

sum(d17$BOUWJAAR == -1)
sum(d17$BOUWJAAR != -1)
table(d17$BOUWJAAR)

d17$BOUWJAAR[d17$BOUWJAAR == -1] <- NA 
is.numeric(d17$BOUWJAAR)

d17$ageofhouse <- 2017 - d17$BOUWJAAR 

sum(d17$BWPER == 0) 

# BWPER - contraction year period 


###### creating data sets for sold above, at, and below asking price  

d17sellers<- d17[d17$LAATSTVRKOOPPR-d17$TRANSACTIEPRIJS < 0,  ]   # sold above asking price

d17sameprice <- d17[d17$LAATSTVRKOOPPR-d17$TRANSACTIEPRIJS == 0,  ] # sold at asking price

d17buyers<- d17[d17$LAATSTVRKOOPPR-d17$TRANSACTIEPRIJS > 0,  ] # solf below asking price

###### which % sold above, at, below asking price

length(d17sellers$ID)/length(d17$ID)
length(d17sameprice$ID)/length(d17$ID)
length(d17buyers$ID)/length(d17$ID)

##### in 2017, out of 8278 listings, 5879 (71.01%) were ones where the selling price is higher than the asking price
##### 10.50 sold at the asking price (together it's 81% = like what the article we use in the presentation says)
## 18.43 sold below

#summary(d17sellers$LAATSTVRKOOPPR)
#mean(d17sellers$LAATSTVRKOOPPR)
#sd(d17sellers$LAATSTVRKOOPPR)

#summary(d17$LAATSTVRKOOPPR)
#mean(d17$LAATSTVRKOOPPR)
#sd(d17$LAATSTVRKOOPPR)


####### code for precision level such that Nzeros counts the number of zeros at the end of the asking price 
d17sellers$Nzero <- (d17sellers$LAATSTVRKOOPPR %% 10 == 0) +  
  (d17sellers$LAATSTVRKOOPPR %% 100 == 0) + 
  (d17sellers$LAATSTVRKOOPPR %% 1000 == 0) +
  (d17sellers$LAATSTVRKOOPPR %% 10000 == 0) + 
  (d17sellers$LAATSTVRKOOPPR %% 100000 == 0) + 
  (d17sellers$LAATSTVRKOOPPR %% 1000000 == 0) + 
  (d17sellers$LAATSTVRKOOPPR %% 10000000 == 0)


######## code number of zeros also for all the 2017 data

d17$Nzero <- (d17$LAATSTVRKOOPPR %% 10 == 0) +  
  (d17$LAATSTVRKOOPPR %% 100 == 0) + 
  (d17$LAATSTVRKOOPPR %% 1000 == 0) +
  (d17$LAATSTVRKOOPPR %% 10000 == 0) + 
  (d17$LAATSTVRKOOPPR %% 100000 == 0) + 
  (d17$LAATSTVRKOOPPR %% 1000000 == 0) + 
  (d17$LAATSTVRKOOPPR %% 10000000 == 0)


############### code for the number of digits the asking price has 

d17sellers$Ndigits <- nchar(d17sellers$LAATSTVRKOOPPR, type = "bytes")

d17$Ndigits <- nchar(d17$LAATSTVRKOOPPR, type = "bytes")

#summary(d17$LAATSTVRKOOPPR[d17$Ndigits == 5])
#summary(d17$LAATSTVRKOOPPR[d17$Ndigits == 6])
#summary(d17$LAATSTVRKOOPPR[d17$Ndigits == 7])

### n of listing with 5 , 6, and 7 digits among listings sold above asking price

sum(d17sellers$Ndigits == 5)
sum(d17sellers$Ndigits == 7)
sum(d17sellers$Ndigits == 6)

#prop of listings with 6 digits

sum(d17sellers$Ndigits == 6)/length(d17sellers$ID)


### n of listing with 6, digits among listings sold below asking price

sum(d17$sale.minus.ask < 0 & d17$Ndigits ==6)

#prop of listings with 6 digits

sum(d17$sale.minus.ask < 0 & d17$Ndigits ==6)/
  sum(d17$sale.minus.ask < 0)



#sum(d17$Ndigits == 5 & d17$sale.minus.ask != 0)
#sum(d17$Ndigits == 7 & d17$sale.minus.ask != 0)
#sum(d17$Ndigits == 6 & d17$sale.minus.ask != 0)
#sum(d17$Ndigits == 6 & d17$sale.minus.ask != 0)/length(d17$ID[d17$sale.minus.ask != 0])

#mean(d17$TRANSACTIEPRIJS[d17$Ndigits ==6 & d17$sale.minus.ask != 0])
#summary(d17$TRANSACTIEPRIJS[d17$Ndigits ==6 & d17$sale.minus.ask != 0])
#sd(d17$TRANSACTIEPRIJS[d17$Ndigits ==6 & d17$sale.minus.ask != 0])




#mean(d17$TRANSACTIEPRIJS[d17$Ndigits ==6 & d17$sale.minus.ask != 0])
#summary(d17$TRANSACTIEPRIJS[d17$Ndigits ==6 & d17$sale.minus.ask != 0])
#sd(d17$TRANSACTIEPRIJS[d17$Ndigits ==6 & d17$sale.minus.ask != 0])


# ----------------------------------------------------------------------------------------
## creating new data set with 6 digits at the asking price
# ----------------------------------------------------------------------------------------


## for listings sold above the asking price

d17sellersnew <- d17sellers[d17sellers$Ndigits == 6, ]

#summary(d17sellersnew$LAATSTVRKOOPPR)
#sd(d17sellersnew$LAATSTVRKOOPPR)

## for all listings 

d17new <- d17[d17$Ndigits == 6 & d17$sale.minus.ask != 0, ]

#summary(d17new$LAATSTVRKOOPPR)
#sd(d17new$LAATSTVRKOOPPR)

### proportions of precision level  

(sum(d17sellersnew$Nzero == 0) + sum(d17sellersnew$Nzero == 1))/length(d17sellersnew$Nzero)  # precise to the 1s and 10s
sum(d17sellersnew$Nzero == 2)/length(d17sellersnew$Nzero) # precise to the 100s
sum(d17sellersnew$Nzero == 3)/length(d17sellersnew$Nzero) # precise to the 1000s
sum(d17sellersnew$Nzero == 4)/length(d17sellersnew$Nzero) # precise to the 10000s
sum(d17sellersnew$Nzero == 5)/length(d17sellersnew$Nzero) # precise to the 100000s


#sample(d17sellersnew$LAATSTVRKOOPPR[d17sellersnew$Nzero == 0], 1)
#sample(d17sellersnew$LAATSTVRKOOPPR[d17sellersnew$Nzero == 1], 1)
#sample(d17sellersnew$LAATSTVRKOOPPR[d17sellersnew$Nzero == 2], 1)
#sample(d17sellersnew$LAATSTVRKOOPPR[d17sellersnew$Nzero == 3], 1)
#sample(d17sellersnew$LAATSTVRKOOPPR[d17sellersnew$Nzero == 4], 1)
#sample(d17sellersnew$LAATSTVRKOOPPR[d17sellersnew$Nzero == 5], 1)


############## pie chart ####################


pdf("Amsterdam_2017_preicsion.pdf")

slices <- c(10.99,20.98, 62.99, 5.02)
lbls <- c("Round", "Precise (1000s)", "Precise (100s)", "Precise (10s and 1s)")
pie(slices, labels = c("10.99%", "20.98%", "62.99", "5.02%"), cex = 1.2, main= "N = 5713", col = c("deepskyblue2", "gray83" ,"gray70", "gray57"), init.angle = 90)
dev.off()


# ----------------------------------------------------------------------------------------
#  main analyses (the analyses reported in the  paper)  

## predicting the log of selling price for:          
# model 1: number of zeros at the end of the asking price 
# model 2: number of zeros + control variables
# model 3: number of zeros + control variables + log of asking price 
# model 4: adding sold below the asking price, 

#  in models 2, 3, and 4 we  control for the following 
#(1) asking price [LAATSTVRKOOPPR] (log, in models 3 and 4)
# (2) number of days the property was on the market [TOM] (models 2-4)
# (3) district in which the property was in [BUURT] (models 2-4)
# (4) type of property [SOORTWONING] (models 2-4)
# (5) size of the property [WOON] (models 2-4)
# (6) Construction year period [BWPER] (instead of age of the property [ageofhouse]) because there are less missing values (models 2-4)
# (7) number of rooms in the property[NKAMERS] (models 2-4)
# (8) the real-estate agent that dealt with the property [MAKELAARSKANTOOR].  (models 2-4)


# ----------------------------------------------------------------------------------------
# for prices with 6 digits in 2017 (the data reported in the paper)


# calculate a log transform of the transaction price

d17sellersnew$logtransaction <- log(d17sellersnew$TRANSACTIEPRIJS)


### model 1

summary(lm(logtransaction ~ Nzero, data = d17sellersnew))

### model 2

summary(lm(logtransaction ~ Nzero + TOM + BWPER + WOON + NKAMERS +as.factor(SOORTWONING) + as.factor(district) + as.factor(MAKELAARSKANTOOR), data = d17sellersnew))

### model 3

summary(lm(logtransaction ~ Nzero + log(LAATSTVRKOOPPR) + TOM + BWPER + WOON + NKAMERS +as.factor(SOORTWONING) + as.factor(district) + as.factor(MAKELAARSKANTOOR), data = d17sellersnew))

### for model 4: 

#### adding the data for the houses that are sold below the asking price as well, and the interaction

options(max.print=9999)

########### log transform the transaction price in the dataset with all listings (with 6 digits)

d17new$logtransaction <- log(d17new$TRANSACTIEPRIJS)

## code for the selling price above the asking price (0) or below it (1) [in the analysis we exclude those where the asking price = selling price]

d17new$buyersmarket[d17new$LAATSTVRKOOPPR < d17new$TRANSACTIEPRIJS] <- 0
d17new$buyersmarket[d17new$LAATSTVRKOOPPR > d17new$TRANSACTIEPRIJS ] <- 1

table(d17new$buyersmarket, useNA = "always")

summary(lm(logtransaction ~ Nzero + buyersmarket + Nzero*buyersmarket + log(LAATSTVRKOOPPR)+ TOM + BWPER + WOON + NKAMERS +as.factor(SOORTWONING) + as.factor(district) + as.factor(MAKELAARSKANTOOR), data = d17new[!is.na(d17new$buyersmarket), ]))

## differences between prices of those who sold above and below 

summary(lm(logtransaction ~ buyersmarket, data = d17new))

#summary(lm(logtransaction ~ buyersmarket, data = d17new[!is.na(d17new$buyersmarket), ]))

aggregate(d17new$LAATSTVRKOOPPR, by = list(d17new$buyersmarket), FUN = mean)
aggregate(d17new$LAATSTVRKOOPPR, by = list(d17new$buyersmarket), FUN = sd)
aggregate(d17new$LAATSTVRKOOPPR, by = list(d17new$buyersmarket), FUN = length)

1287/(1287+5713) #proportion of properties sold below the asking price


# ----------------------------------------------------------------------------------------
# for som
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# regression analysis reported in SOM (Table S1) 
# ----------------------------------------------------------------------------------------


# calculate level of roundness (NUMBER OF ZEROS AT THE END OF ASKING PRICE DEVIDED BY N OF DIGITS
# and log transform the transaction price

# for dataset only of those sold above asking price

d17sellers$roundnessLevel <- d17sellers$Nzero/d17sellers$Ndigits

d17sellers$logtransaction <- log(d17sellers$TRANSACTIEPRIJS)

# for all dataset
d17$roundnessLevel <- d17$Nzero/d17$Ndigits

d17$logtransaction <- log(d17$TRANSACTIEPRIJS)

### model 1

summary(lm(logtransaction ~ roundnessLevel, data = d17[d17$sale.minus.ask > 0, ]))

### model 2

summary(lm(logtransaction ~ roundnessLevel + TOM + BWPER + WOON + NKAMERS +as.factor(SOORTWONING) + as.factor(district) + as.factor(MAKELAARSKANTOOR), data = d17[d17$sale.minus.ask > 0, ]))

### model 3

summary(lm(logtransaction ~ roundnessLevel + log(LAATSTVRKOOPPR)+ TOM + BWPER + WOON + NKAMERS +as.factor(SOORTWONING) + as.factor(district) + as.factor(MAKELAARSKANTOOR), data = d17[d17$sale.minus.ask > 0, ]))

### model 4 

## code for the selling price above the asking price (0) or below it (1) [in the analysis we exclude those where the asking price = selling price]

d17$buyersmarket[d17$LAATSTVRKOOPPR < d17$TRANSACTIEPRIJS] <- 0
d17$buyersmarket[d17$LAATSTVRKOOPPR > d17$TRANSACTIEPRIJS ] <- 1
table(d17$buyersmarket, useNA = "always")


summary(lm(logtransaction ~ roundnessLevel + buyersmarket + roundnessLevel*buyersmarket + log(LAATSTVRKOOPPR)+ TOM + BWPER + WOON + NKAMERS +as.factor(SOORTWONING) + as.factor(district) + as.factor(MAKELAARSKANTOOR), data = d17[!is.na(d17$buyersmarket), ]))


# ----------------------------------------------------------------------------------------
# summary statistics reported in SOM (Table S2) 
# ----------------------------------------------------------------------------------------


############### summary statistics

#asking price

length(d17$LAATSTVRKOOPPR) 
mean(d17$LAATSTVRKOOPPR)
sd(d17$LAATSTVRKOOPPR)

#selling price

length(d17$TRANSACTIEPRIJS) 
mean(d17$TRANSACTIEPRIJS)
sd(d17$TRANSACTIEPRIJS)

#n zeros at the end of the asking price

length(d17$Nzero) 
mean(d17$Nzero)
sd(d17$Nzero)

#n digits of the asking price

length(d17$Ndigits) 
mean(d17$Ndigits)
sd(d17$Ndigits)

# days on the market

length(d17$TOM) - sum(is.na(d17$TOM))
mean(d17$TOM, na.rm = T)
sd(d17$TOM, na.rm = T)

# size of property

length(d17$WOON) 
mean(d17$WOON, na.rm = T)
sd(d17$WOON, na.rm = T)

#n rooms

length(d17$NKAMERS) 
mean(d17$NKAMERS, na.rm = T)
sd(d17$NKAMERS, na.rm = T)

