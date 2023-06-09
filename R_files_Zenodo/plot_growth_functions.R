


# Plot G function
#G_1 = \\frac{r_1 N_1(t)^\\omega}{1 + N_1(t)^\\omega/K_1}

for (omega in 1.4){
curve((r*x^omega)/(1 + (x^omega)/K) ,from=0,to=2,ylab='N(t+1)',xlab='N(t)',
      ylim=c(0,2),xlim=c(0,2),las=1)
abline(coef = c(0,1),lty=2)
mtext(text = paste('omega = ',omega,sep=''),side = 3,line = 1)
while((as.numeric(Sys.time()) - as.numeric(date_time))<2){}
}

for (omega in seq(1,2,0.1)){
r=2
K=1
#omega=1.5
#curve((r*x^omega)/(1 + (x^omega)/K) ,from=0,to=2,ylab='N(t+1)',xlab='N(t)',
#      ylim=c(0,2),xlim=c(0,2),las=1)
curve(-1+(r*x^omega)/(1 + (x^omega)/K)/x ,from=-3,to=2,ylab='Growth rate',xlab='N(t)',
     ylim=c(-1,1),xlim=c(0,2),las=1)
abline(h=0)
#abline(coef = c(0,1),lty=2)
mtext(text = paste('omega = ',omega,sep=''),side = 3,line = 1)

bv_function <- function(x,r,omega) x - r*(x^omega)/(1+(x^omega))
#allee_threshold <- uniroot(bv_function,r=r,omega=omega,interval=c(0.0001,0.9))$root
#abline(v=allee_threshold,col='blue')
print(100*allee_threshold/2)
date_time<-Sys.time()
while((as.numeric(Sys.time()) - as.numeric(date_time))<1.1){}
}
