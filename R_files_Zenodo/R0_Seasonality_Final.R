##### Marta Shocket, Stanford University
##### Updated April 2018

##### Purposes: 1) Calculate R0 as a function of temperature for monthly mean temperatures from example cities across Australia
##              2) Calculate R0 as a function of temperature for nationwide, population-weighted monthly mean temperature and 
##                 compare to data of average monthly RRV cases
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Calculate R0(T) seasonality for 15 largest urban areas in Australia
##           3) Calculate R0(T) seasonality nationwide (weighted by population)
##           4) Manuscript Figures 5 and 6


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("")

# Get R0 function
R0.quantiles <- read.csv("R0Quantiles.csv")
R0.med <- R0.quantiles$R0_full_inf_median

# Get data frame of monthly avg temps for all cities
city.temps.df <- read.csv("RRVCitiesmean_Transposed_Rounded.csv")

# Get RRV human case data
# Data from: https://cameronwebb.wordpress.com/2014/04/09/why-is-mosquito-borne-disease-risk-greater-in-autumn/
# "Average monthly notifications of Ross River virus across Australia between 1992-2013 (National Notifiable Diseases Surveillance System)"
rrv.cases <- read.csv("MonthlyAvgRRVCases.csv")

# Temp list
Temp.xs <- seq(5, 45, 0.2)
N.Temp.xs <-length(Temp.xs)
Temp.xs.ch <- as.character(Temp.xs) # which() isn't working for numbers, so convert to characters


#######
##### 2. Calculate R0(T) seasonality for 15 largest urban areas in Australia
#######

# Create output dataframe for R0(T) calcs
city.R0s <- data.frame(matrix(nrow = 12, ncol = 15))
colnames(city.R0s) <- c("Adelaide", "Brisbane", "Cairns", "Canberra", "Darwin", "Geelong", "Gold Coast", "Hobart",
                        "Melbourne", "Newcastle", "Perth", "Sunshine Coast", "Sydney", "Townsville", "Wollongong")

# Calculate R0 for each month in each city
for(k in 1:ncol(city.R0s)){ # Loop through cities
  
  city.temps <- city.temps.df[, 1+k] # Get temps for a city
  
  for(i in 1:length(city.temps)) { # Loop through months
    j <- which(Temp.xs.ch == as.character(city.temps[i])) # get temp and turn it into correct index
    city.R0s[i, k] <- R0.med[j] # get R0 for index
  }
  
}


#######
##### 3. Calculate R0(T) seasonality nationwide (weighted by population)
#######

# Proportion of population in each metro area, order matching the city.R0s
city.pop.perc <- c(0.0716, 0.1277, 0.0081, 0.0235, 0.0079, 0.0104, 0.0350, 0.0121, 0.2556, 0.0236, 0.1094, 0.0172, 0.2721, 0.0097, 0.0160)
sum(city.pop.perc) # check to see this equals ~1

# Create output dataframe for pop-weighted R0(T) calcs
pop.weighted.R0s <- data.frame(matrix(nrow = 12, ncol = 15))
colnames(pop.weighted.R0s) <- c("Adelaide", "Brisbane", "Cairns", "Canberra", "Darwin", "Geelong", "Gold Coast", "Hobart",
                                "Melbourne", "Newcastle", "Perth", "Sunshine Coast", "Sydney", "Townsville", "Wollongong")

# Calculate population-weighted seasonal R0
for(i in 1:nrow(city.R0s)){ # Loop through months
    
  pop.weighted.R0s[i, ] <- city.R0s[i, ] * city.pop.perc
  
}

# Sum all cities together
pop.weighted.R0.sums <- rowSums(pop.weighted.R0s)


#######
##### Manuscript Figures 5 and 6
#######

##### Sort data frames and vectors so they go from July-June instead of Jan-Dec
city.temps.df.sorted <- rbind(city.temps.df[7:12,], city.temps.df[1:6,])
city.R0s.sorted <- rbind(city.R0s[7:12,], city.R0s[1:6,])
pop.weighted.R0.sums.sorted <- c(pop.weighted.R0.sums[7:12], pop.weighted.R0.sums[1:6])
rrv.cases.sorted <- rbind(rrv.cases[7:12,], rrv.cases[1:6,])
  
###### Figure 5
par(mfrow = c(2,1), mar = c(0.5, 4.25, 1, 1), oma = c(3.25, 0, 0, 0), las = 1, mgp = c(2.6, 1, 0))

plot(city.temps.df$Brisbane ~ seq(1,12,1), xaxt = "n", ylim = c(0,35), type = "l", lwd = "2", col = "white",
     ylab = expression(paste("Monthly mean temperature (",degree,"C)")), xlab = "", cex.axis = 0.9, cex.lab = 1.15)
axis(side = 1, at = seq(1,12,1), labels = c(rep("",12)), cex.axis = 0.9)
cord.x.outer <- c(1, 1, 12, 12)
cord.y.outer <- c(15.8, 33, 33, 15.8)
polygon(cord.x.outer, cord.y.outer, col = "gray93", border = NA)
cord.x.med <- c(1, 1, 12, 12)
cord.y.med <- c(17, 31.4, 31.4, 17)
polygon(cord.x.med, cord.y.med, col = "gray85", border = NA)
cord.x.inner <- c(1, 1, 12, 12)
cord.y.inner <- c(18.4, 30.4, 30.4, 18.4)
polygon(cord.x.inner, cord.y.inner, col = "gray78", border = NA)
segments(1, 26.4, 12, 26.4, lwd = "2", col = "gray50", lty = 2)
lines(city.temps.df.sorted$Darwin ~ seq(1,12,1), lwd = "2", col = "darkred")
lines(city.temps.df.sorted$Cairns ~ seq(1,12,1), lwd = "2", col = "red")
lines(city.temps.df.sorted$Brisbane ~ seq(1,12,1), lwd = "2", col = "sienna2")
lines(city.temps.df.sorted$Perth ~ seq(1,12,1), lwd = "2", col = "orange")
lines(city.temps.df.sorted$Sydney ~ seq(1,12,1), lwd = "2", col = "cyan4")
lines(city.temps.df.sorted$Melbourne ~ seq(1,12,1), lwd = "2", col = "blue")
lines(city.temps.df.sorted$Hobart ~ seq(1,12,1), lwd = "2", col = "darkblue")

legend("topleft", legend = "A", bty = "n", adj = c(1.5, 0))
legend(x = 3.5, y = 11, lwd = 2, col = c("darkred", "red", "sienna2", "orange"),
       legend = c("Darwin", "Cairns", "Brisbane", "Perth"), cex = 0.9, bty = "n")
legend(x = 7.15, y = 11, lwd = 2, col = c("cyan4", "blue", "darkblue"),
       legend = c("Sydney", "Melbourne", "Hobart"), cex = 0.9, bty = "n")

plot(city.R0s.sorted$Brisbane ~ seq(1,12,1), xaxt = "n", ylim = c(0,1), type = "l", lwd = "2", col = "darkorange",
     ylab = expression(paste("Relative ",italic(R)[0],"(",italic(T),")")), xlab = "", cex.axis = 0.9, cex.lab = 1.15)
axis(side = 1, at = seq(1,12,1), labels = c("J", "A", "S", "O", "N", "D", "J", "F", "M", "A", "M", "J"), cex.axis = 0.9)
lines(city.R0s.sorted$Darwin ~ seq(1,12,1), lwd = "2", col = "darkred")
lines(city.R0s.sorted$Cairns ~ seq(1,12,1), lwd = "2", col = "red")
lines(city.R0s.sorted$Brisbane ~ seq(1,12,1), lwd = "2", col = "sienna2")
lines(city.R0s.sorted$Perth ~ seq(1,12,1), lwd = "2", col = "orange")
lines(city.R0s.sorted$Sydney ~ seq(1,12,1), lwd = "2", col = "cyan4")
lines(city.R0s.sorted$Melbourne ~ seq(1,12,1), lwd = "2", col = "blue")
lines(city.R0s.sorted$Hobart ~ seq(1,12,1), lwd = "2", col = "darkblue")

legend("topleft", legend = "B", bty = "n", adj =  c(1.5, 0))
mtext(text = "Month of year (July - June)", side = 1, line = 2.5, cex = 1.15)

###### Figure 6

# Create vectors of modified x-values and y-values so line plot and bar plot line up properly
x.values = c(0.5, seq(1, 12, 1), 12.5)
pop.weighted.R0.sums.sorted.ends <- c(pop.weighted.R0.sums.sorted[1], pop.weighted.R0.sums.sorted, pop.weighted.R0.sums.sorted[12])

# Smaller version - for manuscript
par(mfrow = c(1,1), mar = c(4.5, 4.5, 1, 5.2), oma = c(0, 0, 0, 0), las = 1, mgp = c(3, 1, 0))
barplot(height = rrv.cases.sorted$RRV.Cases, ylim = c(0, 1000), ylab = "Mean monthly RRV cases", col = "darkgray",
        main = "", cex.lab = 1.25, cex.axis = 0.95, names.arg = c("J", "A", "S", "O", "N", "D", "J", "F", "M", "A", "M", "J"))
par(new = TRUE)
plot(pop.weighted.R0.sums.sorted ~ seq(1.5, 12.5, 1), type = "l", lwd = 2, ylim = c(0,0.56), xlim = c(1.1,12.9), main = "", axes = FALSE, bty = "n",
     ylab = "", xlab = "", xaxt = "n", cex.lab = 1.25, cex.axis = 0.95)
axis(side = 4, cex.lab = 0.9)
text(16.1, 0.25, expression(paste("Relative ",italic(R)[0],"(",italic(T),")")), xpd = NA, srt = 270, cex = 1.25)
mtext(text = "Month of year (July - June)", side = 1, line = 2.5, cex = 1.2)