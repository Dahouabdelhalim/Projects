##########################################
#R code for generating NEW Figs. 1-2
#Extinction debt in local habitats: 
#quantifying the roles of random drift, immigration, and emigration
###################
library(Rcpp)
#

##############################
#C version of Sim.Point()
sourceCpp(code='
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector Sim_Point_C(double v, double u, int p0,int N,double tp,double dt)
{
//Rcpp::NumericVector i0=Rcpp::rep_each(0.0,N+1);
Rcpp::NumericVector i0(N+1);
for(int i=0;i<N+1;i++)
{
i0[i]=0;
}//for
i0[p0]=1;

Rcpp::NumericVector m=i0;
Rcpp::NumericVector w=i0;
double t=0.0;
while(t < tp) //gap is important! t<tp is wrong!
{
for(int s=0;s<N+1;s++)
{
if(s==0)
{
m[s]=(w[s+1]-v*w[s])*dt+w[s];
}else if(s==1)
{
m[s]=(v*w[s-1]+2.0*w[s+1]-(2-u)*w[s])*dt+w[s];
}else if(s>1 && s<N)
{
m[s]=((s-1)*(1-u)*w[s-1]+(s+1)*w[s+1]-s*(2-u)*w[s])*dt+w[s];
}else if(s==N)
{
m[s]=((s-1)*(1-u)*w[s-1]-s*w[s])*dt+w[s];
}//if
}//for
w=m;
t=t+dt;
}//while

//final output
if(t+dt>tp)
{
return((w+m)/2);
}else if(t+dt==tp)
{
return(m);
}//if

}
')




##################################
#simulate point value at specific time point
Sim.Point<-function(v,u,p0,N,tp,dt)
{
#dt=min(.002,tp/10000)
#
i0=m=rep(0,1,N+1)
i0[p0+1]=1
w=i0
t=0
while(t<tp)
{
for(s in 1:(N+1)) #all the states including 0
{
if(s==1) #p0
{
#p[s,t]=(p[s+1,t-1]-v*p[s,t-1])*dt+p[s,t-1]
m[s]=(w[s+1]-v*w[s])*dt+w[s]
}else if(s==2) #p1
{
#p[s,t]=(v*p[s-1,t-1]+2*p[s+1,t-1]-(2-u)*p[s,t-1])*dt+p[s,t-1]
m[s]=(v*w[s-1]+2*w[s+1]-(2-u)*w[s])*dt+w[s]
}else if(s>2 & s<N+1) #p(s-1)
{
#p[s,t]=((s-2)*(1-u)*p[s-1,t-1]+s*p[s+1,t-1]-(s-1)*(2-u)*p[s,t-1])*dt+p[s,t-1]
m[s]=((s-2)*(1-u)*w[s-1]+s*w[s+1]-(s-1)*(2-u)*w[s])*dt+w[s]
}else if(s==N+1) #last state
{
#p[s,t]=((s-2)*(1-u)*p[s-1,t-1]-(s-1)*p[s,t-1])*dt+p[s,t-1]
m[s]=((s-2)*(1-u)*w[s-1]-(s-1)*w[s])*dt+w[s] #2019-9-7: correct version!!!
}
}#s
#########
w=m #update w
t=t+dt #update time
}#while
if(t+dt>tp)
{
return((w+m)/2) #mean
}else if(t+dt==tp)
{
return(m) #the exact point
}
}#
#



#########################
#Rosindell's analytical solution
Rosindell<-function(u,p0,T,dt)
{
tt=(0:(T+1))*dt
if(u==0)
{
p=(tt/(1+tt))^p0
}else
{
p=((1-exp(-u*tt))/(1-(1-u)*exp(-u*tt)))^p0
}
return(p)
}#
#




#metacommunity Fisher's logseries
phi.meta<-function(n,A,omega)
{
aA=1/(log(1+A/omega))
#
if(n==0)
{
return(0)
}else
{
return(aA/n*(1-exp(-1/aA))^n)
}
}#
#
#


#local community Fisher's logseries
phi.local<-function(n,A,a,omega)
{
aA=1/(log(1+A/omega))
aa=1/(log(1+a/omega))
#
if(n==0)
{
return(1-aA/aa)
}else
{
return(aA/n*(1-exp(-1/aa))^n)
}
}#
#


########################################
#Fisher's logseries extinction debt
#mu is mutation rate, omega is parameter for Fisher's logseries
EDa<-function(t,Stot,a,A,mu,omega,abundance.upperbound=1000)
{
v=0
if(mu==0)
{
for(j in 1:abundance.upperbound)
{
v=v+(t/(1+t))^j*Stot*phi.local(n=j,A=A,a=a,omega=omega)
}#j
}else
{
for(j in 1:abundance.upperbound)
{
v=v+((1-exp(-mu*t))/(1-(1-mu)*exp(-mu*t)))^j*Stot*phi.local(n=j,A=A,a=a,omega=omega)
}#j
}
return(v)
}#end
#



########################################
#cumulative extinction for the case when immigration=0
#Fisher's logseries extinction debt with instant habitat destruction
#tao is the time point when the habitat b is destructed
#b is the size of destructed habitat
EDb<-function(t,tao,Stot,a,b,A,mu,omega,abundance.upperbound=100)
{
if(t<tao)
{
v=0
if(mu==0)
{
for(j in 1:abundance.upperbound)
{
v=v+(t/(1+t))^j*Stot*phi.local(n=j,A=A,a=a,omega=omega)
}#j
}else
{
for(j in 1:abundance.upperbound)
{
v=v+((1-exp(-mu*t))/(1-(1-mu)*exp(-mu*t)))^j*Stot*phi.local(n=j,A=A,a=a,omega=omega)
}
}#else
########
}else if(t>=tao)
{
#sps.num=Stot*(1-phi.local(n=0,A=A,a=A-b,omega=omega))*phi.local(n=j,A=A-b,a=a,omega=omega)
v=0
if(mu==0)
{
for(j in 1:abundance.upperbound)
{
sps.num=Stot*(1-phi.local(n=0,A=A,a=A-b,omega=omega))*phi.local(n=j,A=A,a=a,omega=omega)
#sps.num=Stot*(1-phi.local(n=0,A=A,a=A-b,omega=omega))*phi.local(n=j,A=A-b,a=a,omega=omega)

v=v+(t/(1+t))^j*sps.num
}#j
}else
{
for(j in 1:abundance.upperbound)
{
sps.num=Stot*(1-phi.local(n=0,A=A,a=A-b,omega=omega))*phi.local(n=j,A=A,a=a,omega=omega)
#sps.num=Stot*(1-phi.local(n=0,A=A,a=A-b,omega=omega))*phi.local(n=j,A=A-b,a=a,omega=omega)

v=v+((1-exp(-mu*t))/(1-(1-mu)*exp(-mu*t)))^j*sps.num
}#j
}
}#
#
#imminent loss of species
ims=Stot*phi.local(n=0,A=A,a=A-b,omega=omega)
#
return(list(extinction.debt=v,imminent.extinction=ims))
}#end
#


#####################################
#cumulative extinction using simulation method
ED<-function(t,tao,Stot,a,b,A,mv,mu,omega,abundance.upperbound=100,dt=.002)
{
if(t<tao)
{
v=0
for(j in 1:abundance.upperbound)
{
p=Sim_Point_C(v=mv,u=mu,p0=j,N=2*abundance.upperbound,tp=t,dt=dt)
v=v+p[1]*Stot*phi.local(n=j,A=A,a=a,omega=omega)
}#j
########
}else if(t>=tao)
{
#sps.num=Stot*(1-phi.local(n=0,A=A,a=A-b,omega=omega))*phi.local(n=j,A=A-b,a=a,omega=omega)
v=0
for(j in 1:abundance.upperbound)
{
p=Sim_Point_C(v=mv,u=mu,p0=j,N=2*abundance.upperbound,tp=t,dt=dt)
#
sps.num=Stot*(1-phi.local(n=0,A=A,a=A-b,omega=omega))*phi.local(n=j,A=A,a=a,omega=omega)
v=v+p[1]*sps.num
}#j
}#else
#
#imminent loss of species
ims=Stot*phi.local(n=0,A=A,a=A-b,omega=omega)
#
return(list(extinction.debt=v,imminent.extinction=ims))
}#end
#

#####################################
#CSF logseries
Fisher<-function(n,a,A,omega)
{
aA=1/log(1+A/omega)
aa=1/log(1+a/omega)
#
if(n==0)
{
return(1-aA/aa)
}else
{
return(aA/n*(1-exp(-1/aa))^n)
}
}#
#


##########################################
#species richness dynamic
SR1<-function(t,tao,Stot,a,b,A,v,mu,omega,abundance.upperbound=1000)
{
if(t<tao)
{
v=0
for(j in 1:abundance.upperbound)
{
p=MatSol.Point(v,u=mu,p0=j,N=abundance.upperbound+1,tp=t)
v=v+(1-p[1,1])*Stot*phi.local(n=j,A=A,a=a,omega=omega)
}#j

########
}else if(t>=tao)
{
#sps.num=Stot*(1-phi.local(n=0,A=A,a=A-b,omega=omega))*phi.local(n=j,A=A-b,a=a,omega=omega)
v=0
for(j in 1:abundance.upperbound)
{
sps.num=Stot*(1-phi.local(n=0,A=A,a=A-b,omega=omega))*phi.local(n=j,A=A,a=a,omega=omega)
p=MatSol.Point(v,u=mu,p0=j,N=abundance.upperbound+1,tp=t)
v=v+(1-p[1,1])*sps.num
}#j

}#else
#
#imminent loss of species
ims=Stot*phi.local(n=0,A=A,a=A-b,omega=omega)
#
return(list(species.richness=v,imminent.extinction=ims))
}#end
#













##############################################################################
##############################################################################
##############################################################################
#####################empirical analyses#######################################
##############################################################################
#Fig. 1: impact of destructed habitat size
#Note: the computation can be very slow
out=out1=out2=out0=vector()
for(i in 0:1000)
{
out0[i]=ED(t=i,tao=100,Stot=1000,a=1,b=0,A=100,mv=.02,mu=.01,omega=.001)$extinction.debt
out[i]=ED(t=i,tao=100,Stot=1000,a=1,b=10,A=100,mv=.02,mu=.01,omega=.001)$extinction.debt
out1[i]=ED(t=i,tao=100,Stot=1000,a=1,b=50,A=100,mv=.02,mu=.01,omega=.001)$extinction.debt
out2[i]=ED(t=i,tao=100,Stot=1000,a=1,b=80,A=100,mv=.02,mu=.01,omega=.001)$extinction.debt
}#i
#
windows()
plot(1:1000,out0,xlab="time",ylab="total species loss",type="l")
lines(1:1000,out,col=2)
lines(1:1000,out1,col=3)
lines(1:1000,out2,col=4)
legend(600,300,legend=c(0,10,50,80),title="Destructed area size",col=c(1,2,3,4),lwd=1,bty="n",border="green")
lines(rep(100,1,500),1:500,lwd=3,lty=3,col="grey")
#
#


############################################
#Fig. 2: impact of original SAD
#Note: the computation can be very slow
out=out1=out2=out0=vector()
for(i in 0:1000)
{
out0[i]=ED(t=i,tao=100,Stot=1000,a=1,b=0,A=100,mv=.02,mu=.001,omega=.001)$extinction.debt
out[i]=ED(t=i,tao=100,Stot=1000,a=1,b=50,A=100,mv=.02,mu=.001,omega=.001)$extinction.debt
out1[i]=ED(t=i,tao=100,Stot=1000,a=1,b=0,A=100,mv=.02,mu=.001,omega=.1)$extinction.debt
out2[i]=ED(t=i,tao=100,Stot=1000,a=1,b=50,A=100,mv=.02,mu=.001,omega=.1)$extinction.debt
}#i
#
windows()
plot(1:1000,out0,xlab="time",ylab="total species loss",type="l")
lines(1:1000,out,col=2)
lines(1:1000,out1,col=3)
lines(1:1000,out2,col=4)
legend(120,300,legend=c("destructed area size=0 & omega=0.001",
"destructed area size=50 & omega=0.001","destructed area size=0 & omega=0.1",
"destructed area size=50 & omega=0.1"),col=c(1,2,3,4),lwd=1,bty="n",border="green")
lines(rep(100,1,500),1:500,lwd=3,lty=3,col="grey")
#




