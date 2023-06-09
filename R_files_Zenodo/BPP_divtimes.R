#Read in data
timing_posterior <- read.table("subsampled_divtime.txt", header = TRUE, sep = "\\t")
#Add columns for generation time and mutation rate samples. These arguments can be changed
new <- cbind(timing_posterior,g=rgamma(10000, 16, 16),u=rgamma(10000, 5.975, 10.864))
#Add column for per year mutation rate
r <- numeric()
for (i in 1:length(new$g)) {
  r <- c(r,new$u[i]/new$g[i])
}
new <- cbind(new,r)

t_laredoensisA <- numeric()
for (i in 1:length(timing_posterior$laredoensisA_tau_6H)) {
  t_laredoensisA <- c(t_laredoensisA,(new$laredoensisA_tau_6H[i]/new$r[i])*(10^8))
}
t_laredoensisB <- numeric()
for (i in 1:length(timing_posterior$laredoensisB_tau_6H)) {
  t_laredoensisB <- c(t_laredoensisB,(new$laredoensisB_tau_6H[i]/new$r[i])*(10^8))
}

new <- cbind(new,t_laredoensisA,t_laredoensisB)
new$t_laredoensisA <- as.numeric(as.character(new$t_laredoensisA))
new$t_laredoensisB <- as.numeric(as.character(new$t_laredoensisB))

divtime_table = data.frame(matrix(ncol = 2, nrow = 0))
names(divtime_table) <- c("divergence_time","species")

for (i in new$t_laredoensisA) {
  if( (i > quantile(new$t_laredoensisA,0.025)) & (i < quantile(new$t_laredoensisA,0.975)) )
    divtime_table[nrow(divtime_table)+1,] <- c(i,"laredoensisA")
}
for (i in new$t_laredoensisB) {
  if( (i > quantile(new$t_laredoensisB,0.025)) & (i < quantile(new$t_laredoensisB,0.975)) )
    divtime_table[nrow(divtime_table)+1,] <- c(i,"laredoensisB")
}

write.table(divtime_table, "calibrated_divtimes.txt", sep = "\\t")
