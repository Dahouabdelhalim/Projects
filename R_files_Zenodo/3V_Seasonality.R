## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Updated August 2020
##
## Purpose: Plot month-of-onset case data for WNV, SLEV, and EEEV along with model predicted R0(T)
##          
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Add vector percentages to WNV case data 
##           3) Calculate R0 for 3 viruses
##           4) Plot Figure 9


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("~/WNV Case Data")

# Read in WNV case data
wnv.data <- read.csv("countiesnew_monthlytemps.csv")

# Read in WNV+ vector proportion data for different states
vec.data <- read.csv("WNVPosMosqbyStateYrSpp_fromPaull.csv")


##########
###### 2. Add vector percentages to WNV case data 
##########

####################### Calculate average vector % calculations across years
state.list <- unique(vec.data$State)
vec.data.out <- data.frame(State = state.list, PIPPercAvg = numeric(length(state.list)), QUIPercAvg = numeric(length(state.list)), TARPercAvg = numeric(length(state.list)))

for(i in 1:length(state.list)){
  vec.data.sub <- subset(vec.data, State == state.list[i])
  vec.data.out$PIPPercAvg[i] <- mean(vec.data.sub$PIPPerc)
  vec.data.out$QUIPercAvg[i] <- mean(vec.data.sub$QUIPerc)
  vec.data.out$TARPercAvg[i] <- mean(vec.data.sub$TARPerc)
}

#### Add in missing data for DE, MD, and TN
# DE and ME's neighboring states (NH, MD) are all 100% PIP
# TN is bordered by KY in the north (100% PIP), and MS, AL, and GA on the south (all 99-100% QUI), so we'll split 50/50 between them
vec.data.out.2 <- data.frame(State = c("DE", "ME", "TN"), PIPPercAvg = c(1, 1, 0.5), QUIPercAvg = c(0, 0, 0.5), TARPercAvg = c(0, 0, 0))
vec.data.out <- rbind(vec.data.out, vec.data.out.2)

write.csv(vec.data.out, file = "WNVPosMosqStateAvgs.csv")

####################### Add to WNV county data set
# Add in empty percent vec columns
wnv.data$PIPPercAvg <- numeric(nrow(wnv.data))
wnv.data$QUIPercAvg <- numeric(nrow(wnv.data))
wnv.data$TARPercAvg <- numeric(nrow(wnv.data))

# Add data from correct state for each row
for(i in 1:nrow(wnv.data)){
  state <- wnv.data$Abbr[i]
  mozzie.perc <- subset(vec.data.out, State == as.character(state))
  
  wnv.data$PIPPercAvg[i] <- mozzie.perc$PIPPercAvg
  wnv.data$QUIPercAvg[i] <- mozzie.perc$QUIPercAvg
  wnv.data$TARPercAvg[i] <- mozzie.perc$TARPercAvg
}


##########
###### 3. Calculate R0 for 3 viruses
##########

# Round monthly temperatures
wnv.data$tmp01 <- round(wnv.data$tmp01, 1)
wnv.data$tmp02 <- round(wnv.data$tmp02, 1)
wnv.data$tmp03 <- round(wnv.data$tmp03, 1)
wnv.data$tmp04 <- round(wnv.data$tmp04, 1)
wnv.data$tmp05 <- round(wnv.data$tmp05, 1)
wnv.data$tmp06 <- round(wnv.data$tmp06, 1)
wnv.data$tmp07 <- round(wnv.data$tmp07, 1)
wnv.data$tmp08 <- round(wnv.data$tmp08, 1)
wnv.data$tmp09 <- round(wnv.data$tmp09, 1)
wnv.data$tmp10 <- round(wnv.data$tmp10, 1)
wnv.data$tmp11 <- round(wnv.data$tmp11, 1)
wnv.data$tmp12 <- round(wnv.data$tmp12, 1)

# Create empty columns for monthly R0s for WNV, SLEV, and EEEV
wnv.data$WNV.R0.01 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.02 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.03 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.04 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.05 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.06 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.07 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.08 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.09 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.10 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.11 <- numeric(nrow(wnv.data))
wnv.data$WNV.R0.12 <- numeric(nrow(wnv.data))

wnv.data$SLEV.R0.01 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.02 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.03 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.04 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.05 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.06 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.07 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.08 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.09 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.10 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.11 <- numeric(nrow(wnv.data))
wnv.data$SLEV.R0.12 <- numeric(nrow(wnv.data))

wnv.data$EEEV.R0.01 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.02 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.03 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.04 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.05 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.06 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.07 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.08 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.09 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.10 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.11 <- numeric(nrow(wnv.data))
wnv.data$EEEV.R0.12 <- numeric(nrow(wnv.data))

# Get and scale median R0 values from R0 calculations (file = "6V_R0Calculation_inf.R")
WNVPIP.R0.median.scaled <- R0.CpipWNV.out$median / max(R0.CpipWNV.out$median)
WNVQUI.R0.median.scaled <- R0.CquiWNV.out$median / max(R0.CquiWNV.out$median)
WNVTAR.R0.median.scaled <- R0.CtarWNV.out$median / max(R0.CtarWNV.out$median)

SLEV.R0.median.scaled <- R0.CtarSLEV.out$median / max(R0.CtarSLEV.out$median)
EEEV.R0.median.scaled <- R0.AtriEEEV.out$median / max(R0.AtriEEEV.out$median)

# Make sure we have the right column index for temperature $tmp
wnv.data <- wnv.data[,-1]
colnames(wnv.data)[53]

######### Calculate R0s for for avg temp for each month
# Loop through months
for(k in 1:12){
  
  # Loop through counties
  for(i in 1:nrow(wnv.data)){
    
    # Get temperature and translate it into an index
    temp <- wnv.data[i, 53 + k]
    if(temp < 2) temp <- 2
    # temp.index <- which(Temp.xs == temp) # This works most of the time, but intermittently fails because of the way R stores floating point numbers
    temp.index <- as.integer((temp - 0.9) * 10) # This formula only works for when Temp.xs = seq(1, x, 0.1)
    
    # Get R0 for each vector + virus
    WNVPIP.R0 <- WNVPIP.R0.median.scaled[temp.index]
    WNVQUI.R0 <- WNVQUI.R0.median.scaled[temp.index]
    WNVTAR.R0 <- WNVTAR.R0.median.scaled[temp.index]
    
    # Weight WNV R0s based on the species in each county
    wnv.data[i, 68 + k] <- WNVPIP.R0 * wnv.data$PIPPercAvg[i] + WNVQUI.R0 * wnv.data$QUIPercAvg[i] + WNVTAR.R0 * wnv.data$TARPercAvg[i]
    wnv.data[i, 80 + k] <- SLEV.R0.median.scaled[temp.index]
    wnv.data[i, 92 + k] <- EEEV.R0.median.scaled[temp.index]
    
  }
  
}


# Subset to only counties with observed WNV cases 
wnv.data.WNVsub <- wnv.data[wnv.data$Timeadjsum>0,]

# Subset to states with observed EEEV
wnv.data.EEEVsub <- subset(wnv.data, Abbr == "AL" | Abbr == "AR" | Abbr == "CT" | Abbr == "FL" | Abbr == "GA" | Abbr == "LA" |
                             Abbr == "ME" |Abbr == "MD" | Abbr == "MA" | Abbr == "MI" | Abbr == "MO" | Abbr == "MT" | Abbr == "NH" |
                             Abbr == "NJ" | Abbr == "NY" |Abbr == "NC" | Abbr == "TX" | Abbr == "RI" | Abbr == "VT" | Abbr == "VA")

# Subset to states with observed SLEV
wnv.data.SLEVsub <- subset(wnv.data, Abbr == "AL" | Abbr == "AZ" | Abbr == "AR" | Abbr == "CA" | Abbr == "FL" | Abbr == "IL" | 
                             Abbr == "IN" | Abbr == "LA" | Abbr == "MI" | Abbr == "MS" | Abbr == "MO" | Abbr == "NV" |
                             Abbr == "NC" | Abbr == "TX" | Abbr == "UT" | Abbr == "WA")

# Calculate county population weights
wnv.data$pop.weight <- wnv.data$pop / sum(wnv.data$pop)
wnv.data.WNVsub$pop.weight <- wnv.data.WNVsub$pop / sum(wnv.data.WNVsub$pop)
wnv.data.SLEVsub$pop.weight <- wnv.data.SLEVsub$pop / sum(wnv.data.SLEVsub$pop)
wnv.data.EEEVsub$pop.weight <- wnv.data.EEEVsub$pop / sum(wnv.data.EEEVsub$pop)

seasonalR0.WNV <- numeric(12)
seasonalR0.SLEV <- numeric(12)
seasonalR0.EEEV <- numeric(12)

seasonalR0.WNV.sub <- numeric(12)
seasonalR0.SLEV.sub <- numeric(12)
seasonalR0.EEEV.sub <- numeric(12)

seasonalR0.WNV[1] <- sum(wnv.data$WNV.R0.01*wnv.data$pop.weight)
seasonalR0.WNV[2] <- sum(wnv.data$WNV.R0.02*wnv.data$pop.weight)
seasonalR0.WNV[3] <- sum(wnv.data$WNV.R0.03*wnv.data$pop.weight)
seasonalR0.WNV[4] <- sum(wnv.data$WNV.R0.04*wnv.data$pop.weight)
seasonalR0.WNV[5] <- sum(wnv.data$WNV.R0.05*wnv.data$pop.weight)
seasonalR0.WNV[6] <- sum(wnv.data$WNV.R0.06*wnv.data$pop.weight)
seasonalR0.WNV[7] <- sum(wnv.data$WNV.R0.07*wnv.data$pop.weight)
seasonalR0.WNV[8] <- sum(wnv.data$WNV.R0.08*wnv.data$pop.weight)
seasonalR0.WNV[9] <- sum(wnv.data$WNV.R0.09*wnv.data$pop.weight)
seasonalR0.WNV[10] <- sum(wnv.data$WNV.R0.10*wnv.data$pop.weight)
seasonalR0.WNV[11] <- sum(wnv.data$WNV.R0.11*wnv.data$pop.weight)
seasonalR0.WNV[12] <- sum(wnv.data$WNV.R0.12*wnv.data$pop.weight)

seasonalR0.WNV.sub[1] <- sum(wnv.data.WNVsub$WNV.R0.01*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[2] <- sum(wnv.data.WNVsub$WNV.R0.02*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[3] <- sum(wnv.data.WNVsub$WNV.R0.03*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[4] <- sum(wnv.data.WNVsub$WNV.R0.04*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[5] <- sum(wnv.data.WNVsub$WNV.R0.05*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[6] <- sum(wnv.data.WNVsub$WNV.R0.06*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[7] <- sum(wnv.data.WNVsub$WNV.R0.07*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[8] <- sum(wnv.data.WNVsub$WNV.R0.08*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[9] <- sum(wnv.data.WNVsub$WNV.R0.09*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[10] <- sum(wnv.data.WNVsub$WNV.R0.10*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[11] <- sum(wnv.data.WNVsub$WNV.R0.11*wnv.data.WNVsub$pop.weight)
seasonalR0.WNV.sub[12] <- sum(wnv.data.WNVsub$WNV.R0.12*wnv.data.WNVsub$pop.weight)

seasonalR0.SLEV[1] <- sum(wnv.data$SLEV.R0.01*wnv.data$pop.weight)
seasonalR0.SLEV[2] <- sum(wnv.data$SLEV.R0.02*wnv.data$pop.weight)
seasonalR0.SLEV[3] <- sum(wnv.data$SLEV.R0.03*wnv.data$pop.weight)
seasonalR0.SLEV[4] <- sum(wnv.data$SLEV.R0.04*wnv.data$pop.weight)
seasonalR0.SLEV[5] <- sum(wnv.data$SLEV.R0.05*wnv.data$pop.weight)
seasonalR0.SLEV[6] <- sum(wnv.data$SLEV.R0.06*wnv.data$pop.weight)
seasonalR0.SLEV[7] <- sum(wnv.data$SLEV.R0.07*wnv.data$pop.weight)
seasonalR0.SLEV[8] <- sum(wnv.data$SLEV.R0.08*wnv.data$pop.weight)
seasonalR0.SLEV[9] <- sum(wnv.data$SLEV.R0.09*wnv.data$pop.weight)
seasonalR0.SLEV[10] <- sum(wnv.data$SLEV.R0.10*wnv.data$pop.weight)
seasonalR0.SLEV[11] <- sum(wnv.data$SLEV.R0.11*wnv.data$pop.weight)
seasonalR0.SLEV[12] <- sum(wnv.data$SLEV.R0.12*wnv.data$pop.weight)

seasonalR0.SLEV.sub[1] <- sum(wnv.data.SLEVsub$SLEV.R0.01*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[2] <- sum(wnv.data.SLEVsub$SLEV.R0.02*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[3] <- sum(wnv.data.SLEVsub$SLEV.R0.03*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[4] <- sum(wnv.data.SLEVsub$SLEV.R0.04*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[5] <- sum(wnv.data.SLEVsub$SLEV.R0.05*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[6] <- sum(wnv.data.SLEVsub$SLEV.R0.06*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[7] <- sum(wnv.data.SLEVsub$SLEV.R0.07*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[8] <- sum(wnv.data.SLEVsub$SLEV.R0.08*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[9] <- sum(wnv.data.SLEVsub$SLEV.R0.09*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[10] <- sum(wnv.data.SLEVsub$SLEV.R0.10*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[11] <- sum(wnv.data.SLEVsub$SLEV.R0.11*wnv.data.SLEVsub$pop.weight)
seasonalR0.SLEV.sub[12] <- sum(wnv.data.SLEVsub$SLEV.R0.12*wnv.data.SLEVsub$pop.weight)

seasonalR0.EEEV[1] <- sum(wnv.data$EEEV.R0.01*wnv.data$pop.weight)
seasonalR0.EEEV[2] <- sum(wnv.data$EEEV.R0.02*wnv.data$pop.weight)
seasonalR0.EEEV[3] <- sum(wnv.data$EEEV.R0.03*wnv.data$pop.weight)
seasonalR0.EEEV[4] <- sum(wnv.data$EEEV.R0.04*wnv.data$pop.weight)
seasonalR0.EEEV[5] <- sum(wnv.data$EEEV.R0.05*wnv.data$pop.weight)
seasonalR0.EEEV[6] <- sum(wnv.data$EEEV.R0.06*wnv.data$pop.weight)
seasonalR0.EEEV[7] <- sum(wnv.data$EEEV.R0.07*wnv.data$pop.weight)
seasonalR0.EEEV[8] <- sum(wnv.data$EEEV.R0.08*wnv.data$pop.weight)
seasonalR0.EEEV[9] <- sum(wnv.data$EEEV.R0.09*wnv.data$pop.weight)
seasonalR0.EEEV[10] <- sum(wnv.data$EEEV.R0.10*wnv.data$pop.weight)
seasonalR0.EEEV[11] <- sum(wnv.data$EEEV.R0.11*wnv.data$pop.weight)
seasonalR0.EEEV[12] <- sum(wnv.data$EEEV.R0.12*wnv.data$pop.weight)

seasonalR0.EEEV.sub[1] <- sum(wnv.data.EEEVsub$EEEV.R0.01*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[2] <- sum(wnv.data.EEEVsub$EEEV.R0.02*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[3] <- sum(wnv.data.EEEVsub$EEEV.R0.03*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[4] <- sum(wnv.data.EEEVsub$EEEV.R0.04*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[5] <- sum(wnv.data.EEEVsub$EEEV.R0.05*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[6] <- sum(wnv.data.EEEVsub$EEEV.R0.06*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[7] <- sum(wnv.data.EEEVsub$EEEV.R0.07*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[8] <- sum(wnv.data.EEEVsub$EEEV.R0.08*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[9] <- sum(wnv.data.EEEVsub$EEEV.R0.09*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[10] <- sum(wnv.data.EEEVsub$EEEV.R0.10*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[11] <- sum(wnv.data.EEEVsub$EEEV.R0.11*wnv.data.EEEVsub$pop.weight)
seasonalR0.EEEV.sub[12] <- sum(wnv.data.EEEVsub$EEEV.R0.12*wnv.data.EEEVsub$pop.weight)

Seasonal.R0s <- data.frame(month = seq(1,12,1), seasonalR0.WNV.all = seasonalR0.WNV, seasonalR0.WNV.sub = seasonalR0.WNV.sub,
                           seasonalR0.SLEV = seasonalR0.SLEV, seasonalR0.SLEV.sub = seasonalR0.SLEV.sub,
                           seasonalR0.EEEV = seasonalR0.EEEV, seasonalR0.EEEV.sub = seasonalR0.EEEV.sub)

# Save output
write.csv(wnv.data, file = "countiesnew_monthlytemps.csv")
write.csv(Seasonal.R0s, file = "SeasonalR0s.csv")


##########
###### 4. Plot Figure 9
##########

# Set working directory
setwd("~/WNV Case Data")

# Get data
month.data <- read.csv("MonthOfOnset.csv")
Seasonal.R0s <- read.csv("SeasonalR0s.csv")

#### Plot
par(mfrow = c(1,1), mar = c(5, 5, 1, 5), las = 1)
plot(WNV/100 ~ Month, data = month.data, type = "l", xlab = "Month of Onset", ylab = expression(paste("Cases (SLEV & EEEV 10"^1,", WNV 10"^-2,")")), lwd = 2, xaxt = "n", cex.lab = 1.3)
axis(side = 1, at = c(2,4,6,8,10,12), labels = c("Feb", "April", "June", "Aug", "Oct", "Dec"))
lines(SLEV ~ Month, data = month.data, col = "gray40", lwd = 2)
lines(EEEV ~ Month, data = month.data, col = "gray70", lwd = 2)
text(x = 2.6, y = 192, labels = "Cases:")
legend(x = 1, y = 190, legend = c("WNV", "SLEV", "EEEV"), lwd = 2, col = c("black", "gray40", "gray70"), bty = "n")
text(x = 2.6, y = 140, labels = expression(paste(italic(R)[0],"(",italic(T),"):")))
legend(x = 1, y = 140, legend = c("WNV", "SLEV", "EEEV"), lwd = 2, col = c("black", "gray40", "gray70"), bty = "n", lty = 2)

par(new = TRUE)
plot(seasonalR0.WNV.sub ~ month, data = Seasonal.R0s, type = "l", col = "black", lwd = 2, lty = 2, yaxt = "n", xaxt = "n", ylab = "", xlab = "")
lines(seasonalR0.SLEV.sub ~ month, data = Seasonal.R0s, col = "gray40", lwd = 2, lty = 2)
lines(seasonalR0.EEEV.sub ~ month, data = Seasonal.R0s, col = "gray70", lwd = 2, lty = 2)

axis(side = 4)
text(14.3, 0.4, expression(paste("Relative ",italic(R)[0],"(",italic(T),")")), xpd = NA, srt = 270, cex = 1.25)