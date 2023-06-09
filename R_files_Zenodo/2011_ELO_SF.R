# ANO 2011 ELO RANKING

setwd("C:/Users/Selene/Documents/SANTA CRUZ/2011 ELO/R FILES")

#dataSheet <- read.csv('2011 ELO DATA ALL.csv', header = T)
dataSheet <- read.csv('2011 ELO DATA WO BS.csv', header = T)
#dataSheet <- read.csv('2011 ELO CALC DEC 21 DATA.csv', header = T)

names <- union(unique(dataSheet[,3]), unique(dataSheet[,4]))
n <- length(names); t <- length(dataSheet[,1])
dim <- n*t
elo.mat <- matrix(rep(0, dim + n), nrow = t + 1, ncol = n, 
	dimnames = list(c(1:(t+1)),names))
animals <- data.frame(rbind (rep(0,5), dataSheet), elo.mat)
animals[1,c(1:(n+5))] <- rep(1000,n+5)


for(i in 2:(t+1)) {
	numW <- which(names == dataSheet$WINN[i-1]) + 5
	rankW <- animals[i-1,numW]
	numL <- which(names == dataSheet$LOSE[i-1]) + 5
	rankL <- animals[i-1,numL]
	
	diffR <- rankL - rankW
	P <- 1/(1+10^(diffR/400))	
		
	animals[i,numW] <- rankW + (1 - P)*100
	animals[i,numL] <- rankL + (P-1)*100
			
	others <- c(6:(n+5))[-c(numW - 5,numL - 5)]	
	for(j in others) animals[i,j] <- animals[i-1,j]
}

#write.csv(animals[-1,], '2011 ELO RESULTS ALL.csv')
write.csv(animals[-1,], '2011 ELO RESULTS WO BS.csv')
#write.csv(animals[-1,], 'Elo_results_Dec_21_final.csv')

