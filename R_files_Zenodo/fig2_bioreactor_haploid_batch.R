# Figure 2
library(ggplot2)
library(reshape2)

#next generation function
p_prime  = function(p, s, mu) {
  return( p*(1-s)*(1-mu) / (1-p*s))
}

#set up data frame
df3 <- data.frame(time=rep((1:100),4), s=rep(c(0.2,0.4),each=200), mu=rep(c(1e-4,1e-8),each=100))
df3$p <- 0
df3[df3$time==1,]$p <-1-df3[df3$time==1,]$mu
threshold = 0.8

for(t in 1:nrow(df3)) {
  if(df3[t,]$time != 1) {
    pp <- p_prime(df3[t-1,]$p, df3[t-1,]$s, df3[t-1,]$mu)
    if(pp < threshold) { pp <- 1-df3[t,]$mu }
    df3[t,]$p <- pp
  }
}

df3$s <- factor(df3$s, levels=c(0.4,0.2))
levels(df3$s) <- c("s~'='~0.4","s~'='~0.2")

df3$mu <- factor(df3$mu, levels=c(1e-4,1e-8))
levels(df3$mu) <- c("mu~'='~1e-4","mu~'='~1e-8")

fig2 <- ggplot(data=df3, aes(x=time, y=p)) + theme_bw() +
  geom_line() + 
  facet_grid(s ~ mu, labeller = label_parsed) + 
  xlab("Generation (t)") +
  ylab( expression(paste('Frequency (',italic("p")['t'],')'))) +
  ylim(0,1)
fig2

ggsave("fig2.pdf", fig2)  