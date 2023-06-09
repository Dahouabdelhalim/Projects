#############################################################################################
#Please see supplementary file S2 of the main paper for a detailed tutorial using this code
#Please direct inquiries to Jordan S. Martin (jsm.primatology@gmail.com)
#############################################################################################

#prepare workspace

#load packages
library(pacman)
p_load(LaplacesDemon,igraph,brms,lavaan,lvnet,HDInterval,plyr,glasso,qgraph, viridis,Deriv,
       robustbase,ggplot2,cowplot,ggjoy,png,reshape2,psych,xlsx,rstan,Matrix,bayesplot, ggstance)

#parallel processing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


###############################################################################################

#set workspace

#load data
data<-read.csv("Martin et al MEE 2018 marmoset data.csv")

#add small constant to duration measures with zeroes
data$Groominitiate <- data$Groominitiate + 0.001
data$Contact.sit <- data$Contact.sit + 0.001

#trait names
trait<-c("groom","proximity","contact","active","scent","gnaw")

#df for repeatability
rep.est<-data.frame(trait=trait, mRint=NA, sdRint=NA,
                                 mRadj=NA, sdRadj=NA,
                                 mRadjshort=NA, sdRadjshort=NA)


#Derive observation-level variance (OLV)

#logit derivation
D(expression(log(u)-log(1-u)),"u")

#derivation
deriv<-"1/u + 1/(1 - u)"

#var x
var_x<-"u*(1 - u)/(1 + phi)"

#OLV
OLV<-paste0(var_x, "*(", deriv,")^2")
Simplify(expression(u*(1 - u)/(1 + phi)*(1/u + 1/(1 - u))^2) )

#function
OLV.Beta<- function(u, phi){
  
  u * (1 - u) * (1/(1 - u) + 1/u)^2/(1 + phi)
  
}

#Poisson
OLV.Poisson<- function(u, OLRE){
  
  OLRE + log(1 + 1/u) 

}

###univariate models###

#####groom############################
m.groom<-           brm(formula=
                        Groominitiate  ~ Month + 
                        (1|Subject) + (1|Group) + (1|I_series),
                         
                         data=data, family=Beta(), 
                         prior=c(prior(normal(0,10),  class= Intercept),
                                 prior(normal(0,10),  class= b),
                                 prior(cauchy(0,3),  class= sd)),
                       
                         warmup=3000,iter=5500,chains=2,seed=9,
                         control=list(adapt_delta=0.99))  

summary(m.groom, prob=0.90)

#extract posterior
groom.post<-posterior_samples(m.groom)

#R
g.int<-groom.post$b_Intercept
g.phi<-groom.post$phi
g.id<-(groom.post$sd_Subject__Intercept)^2
g.group<-(groom.post$sd_Group__Intercept)^2
g.series<- groom.post$sd_I_series__Intercept^2

#sample mean
g.u<-mean(data$Groominitiate)

#OLV
g.OLV <- OLV.Beta(u=g.u, phi = g.phi)


###R_int
g.R_int<- g.id/(g.id + g.series)

median(g.R_int)
rep.est[trait=="groom","mRint"]<-median(g.R_int)

sd(g.R_int)
rep.est[trait=="groom","sdRint"]<-sd(g.R_int)

###R_adj
g.R_adj <- g.id/(g.id+g.series+g.group+g.OLV)

median(g.R_adj)
rep.est[trait=="groom","mRadj"]<-median(g.R_adj)

sd(g.R_adj)
rep.est[trait=="groom","sdRadj"]<-sd(g.R_adj)


###R_adjshort
g.R_adjshort <- (g.id+g.series)/(g.id+g.series+g.group+g.OLV)

median(g.R_adjshort)
rep.est[trait=="groom","mRadjshort"]<-median(g.R_adjshort)

sd(g.R_adj)
rep.est[trait=="groom","sdRadjshort"]<-sd(g.R_adjshort)



#####prox#############################
m.prox<-           brm(formula=
                          Proximity  ~ Month + 
                          (1|Subject) + (1|Group) + (1|I_series),
                        
                        data=data, family=Beta(), 
                        prior=c(prior(normal(0,10), class= b),
                                prior(normal(0,10),  class= Intercept),
                                prior(cauchy(0,3),  class= sd)),
                        
                        warmup=2000,iter=5000,chains=4,seed=9,
                        control=list(adapt_delta=0.99))  


summary(m.prox, prob=0.90)

#extract posterior
prox.post<-posterior_samples(m.prox)

#R
p.int<-prox.post$b_Intercept
p.phi<-prox.post$phi
p.id<-(prox.post$sd_Subject__Intercept)^2
p.group<-(prox.post$sd_Group__Intercept)^2
p.series<- prox.post$sd_I_series__Intercept^2

p.u<-mean(data$Proximity)
  
#OLV
p.OLV <- OLV.Beta(u=p.u, phi = p.phi)

###R_int
p.R_int<- p.id/(p.id + p.series)

median(p.R_int)
rep.est[trait=="proximity","mRint"]<-median(p.R_int)

sd(p.R_int)
rep.est[trait=="proximity","sdRint"]<-sd(p.R_int)

###R_adj
p.R_adj <- p.id/(p.id+p.series+p.group+p.OLV)

median(p.R_adj)
rep.est[trait=="proximity","mRadj"]<-median(p.R_adj)

sd(p.R_adj)
rep.est[trait=="proximity","sdRadj"]<-sd(p.R_adj)

###R_adjshort
p.R_adjshort <- (p.id+p.series)/(p.id+p.series+p.group+p.OLV)

median(p.R_adjshort)
rep.est[trait=="proximity","mRadjshort"]<-median(p.R_adjshort)

sd(p.R_adjshort)
rep.est[trait=="proximity","sdRadjshort"]<-sd(p.R_adjshort)

#####contact###########################
m.cts<-           brm(formula=
                         Contact.sit  ~ Month +
                         (1|Subject) + (1|Group) + (1|I_series),
                       
                       data=data, family=Beta(), 
                       prior=c(prior(normal(0,10),  class= b),
                               prior(normal(0,10),  class= Intercept),
                               prior(cauchy(0,3),  class= sd)),
                       
                       warmup=2000,iter=5000,chains=4,seed=9,
                       control=list(adapt_delta=0.99))  

#extract posterior
cts.post<-posterior_samples(m.cts)

#R
cs.int<-cts.post$b_Intercept
cs.phi<-cts.post$phi
cs.id<-(cts.post$sd_Subject__Intercept)^2
cs.group<-(cts.post$sd_Group__Intercept)^2
cs.series<- cts.post$sd_I_series__Intercept^2

cs.u <- mean(data$Contact.sit)

#OLV
cs.OLV <- OLV.Beta(u=cs.u, phi = cs.phi)

###R_int
cs.R_int<- cs.id/(cs.id + cs.series)

median(cs.R_int)
rep.est[trait=="contact","mRint"]<-median(cs.R_int)

sd(cs.R_int)
rep.est[trait=="contact","sdRint"]<-sd(cs.R_int)

###R_adj
cs.R_adj <- cs.id/(cs.id+cs.series+cs.group+cs.OLV)

median(cs.R_adj)
rep.est[trait=="contact","mRadj"]<-median(cs.R_adj)

sd(cs.R_adj)
rep.est[trait=="contact","sdRadj"]<-sd(cs.R_adj)


###R_adjshort
cs.R_adjshort <- (cs.id+cs.series)/(cs.id+cs.series+cs.group+cs.OLV)

median(cs.R_adjshort)
rep.est[trait=="contact","mRadjshort"]<-median(cs.R_adjshort)

sd(cs.R_adjshort)
rep.est[trait=="contact","sdRadjshort"]<-sd(cs.R_adjshort)


#####active###########################
m.active<-           brm(formula=
                          Active  ~ Month  + 
                          (1|Subject) + (1|Group) + (1|OLRE) + (1|I_series),
                        
                        data=data, family=poisson(), 
                        prior=c(prior(normal(0,10),  class= b),
                                prior(normal(0,10),  class= Intercept),
                                prior(cauchy(0,3),  class= sd)),
                        
                        warmup=2500,iter=5500,chains=4,seed=9,
                        control=list(adapt_delta=0.99)) 


summary(m.active, prob=0.90)

#extract posterior
active.post<-posterior_samples(m.active)

#R
act.int<-active.post$b_Intercept
act.id<-(active.post$sd_Subject__Intercept)^2
act.group<-(active.post$sd_Group__Intercept)^2
act.OLRE<-(active.post$sd_OLRE__Intercept)^2
act.series<- active.post$sd_I_series__Intercept^2

act.u <- mean(data$Active)

#OLV
act.OLV <- OLV.Poisson(u=act.u, OLRE = act.OLRE)

###R_int
act.R_int<- act.id/(act.id + act.series)

median(act.R_int)
rep.est[trait=="active","mRint"]<-median(act.R_int)

sd(act.R_int)
rep.est[trait=="active","sdRint"]<-sd(act.R_int)


###R_adj
act.R_adj <- act.id/(act.id+act.series+act.group+act.OLV)

median(act.R_adj)
rep.est[trait=="active","mRadj"]<-median(act.R_adj)

sd(act.R_adj)
rep.est[trait=="active","sdRadj"]<-sd(act.R_adj)

###R_adjshort
act.R_adjshort <- (act.id+act.series)/(act.id+act.series+act.group+act.OLV)

median(act.R_adjshort)
rep.est[trait=="active","mRadjshort"]<-median(act.R_adjshort)

sd(act.R_adj)
rep.est[trait=="active","sdRadjshort"]<-sd(act.R_adjshort)


#####gnawing###########################
m.gnaw<-           brm(formula=
                           Gnawing  ~ Month + 
                           (1|Subject) + (1|Group) + (1|I_series) + (1|OLRE),
                         
                         data=data, family=poisson(), 
                         prior=c(prior(normal(0,10),  class= b),
                                 prior(normal(0,10),  class= Intercept),
                                 prior(cauchy(0,3),  class= sd)),
                         
                         warmup=2000,iter=5000,chains=4,seed=9,
                         control=list(adapt_delta=0.99)) 


summary(m.gnaw, prob=0.90)

#extract posterior
gnaw.post<-posterior_samples(m.gnaw)

#R
gnw.int<-gnaw.post$b_Intercept
gnw.id<-(gnaw.post$sd_Subject__Intercept)^2
gnw.group<-(gnaw.post$sd_Group__Intercept)^2
gnw.OLRE<-(gnaw.post$sd_OLRE__Intercept)^2
gnw.series<- gnaw.post$sd_I_series__Intercept^2

gnw.u<-mean(data$Gnawing)

#OLV
gnw.OLV <- OLV.Poisson(u=gnw.u, OLRE = gnw.OLRE)

###R_int
gnw.R_int<- gnw.id/(gnw.id + gnw.series)

median(gnw.R_int)
rep.est[trait=="gnaw","mRint"]<-median(gnw.R_int)

sd(gnw.R_int)
rep.est[trait=="gnaw","sdRint"]<-sd(gnw.R_int)

###R_adj
gnw.R_adj <- gnw.id/(gnw.id+gnw.series+gnw.group+gnw.OLV)

median(gnw.R_adj)
rep.est[trait=="gnaw","mRadj"]<-median(gnw.R_adj)

sd(gnw.R_adj)
rep.est[trait=="gnaw","sdRadj"]<-sd(gnw.R_adj)

###R_adjshort
gnw.R_adjshort <- (gnw.id+gnw.series)/(gnw.id+gnw.series+gnw.group+gnw.OLV)

median(gnw.R_adjshort)
rep.est[trait=="gnaw","mRadjshort"]<-median(gnw.R_adjshort)

sd(gnw.R_adjshort)
rep.est[trait=="gnaw","sdRadjshort"]<-sd(gnw.R_adjshort)


#####scentmark###########################
m.scent<-           brm(formula=
                           Scentmark  ~ Month + 
                           (1|Subject) + (1|Group) + (1|I_series) + (1|OLRE),
                         
                         data=data, family=poisson(), 
                         prior=c(prior(normal(0,10),  class= b),
                                 prior(normal(0,10),  class= Intercept),
                                 prior(cauchy(0,3),  class= sd)),
                         
                         warmup=2000,iter=5000,chains=4,seed=9,
                         control=list(adapt_delta=0.99)) 


summary(m.scent, prob=0.90)

#extract posterior
scent.post<-posterior_samples(m.scent)

#R
scm.int<-scent.post$b_Intercept
scm.phi<-scent.post$phi
scm.id<-(scent.post$sd_Subject__Intercept)^2
scm.group<-(scent.post$sd_Group__Intercept)^2
scm.OLRE<-(scent.post$sd_OLRE__Intercept)^2
scm.series<- scent.post$sd_I_series__Intercept^2

scm.u<-mean(data$Scentmark)

#OLV
scm.OLV <- OLV.Poisson(u=scm.u, OLRE = scm.OLRE)


###R_int
scm.R_int<- scm.id/(scm.id + scm.series)

median(scm.R_int)
rep.est[trait=="scent","mRint"]<-median(scm.R_int)

sd(scm.R_int)
rep.est[trait=="scent","sdRint"]<-sd(scm.R_int)

###R_adj
scm.R_adj <- scm.id/(scm.id+scm.series+scm.group+scm.OLV)

median(scm.R_adj)
rep.est[trait=="scent","mRadj"]<-median(scm.R_adj)

sd(scm.R_adj)
rep.est[trait=="scent","sdRadj"]<-sd(scm.R_adj)

###R_adj
scm.R_adjshort <- (scm.id+scm.series)/(scm.id+scm.series+scm.group+scm.OLV)

median(scm.R_adjshort)
rep.est[trait=="scent","mRadjshort"]<-median(scm.R_adjshort)

sd(scm.R_adjshort)
rep.est[trait=="scent","sdRadjshort"]<-sd(scm.R_adjshort)


#overview of repeatability
rep.est

######multivariate model###############################################################

#build each univariate model
mv.agrm<-        bf(formula=
                     Groominitiate ~ Month + 
                     (1|id|Subject) + (1|Group)) + Beta()

mv.prox<-         bf(formula=
                       Proximity ~ Month +
                       (1|id|Subject) + (1|Group)) + Beta()

mv.cont<-         bf(formula=
                       Contact.sit ~ Month +
                       (1|id|Subject) + (1|Group)) + Beta()

mv.actv<-          bf(formula=
                       Active ~ Month +
                       (1|id|Subject) + (1|Group) + (1|OLRE)) + poisson()

mv.scnt<-          bf(formula=
                       Scentmark ~ Month +
                       (1|id|Subject) + (1|Group)+ (1|OLRE)) + poisson()

mv.gnaw<-          bf(formula=
                       Gnawing ~ Month +
                       (1|id|Subject) + (1|Group) + (1|OLRE)) + poisson()

#full model
mv.marmoset<-           brm(formula=
                            mv.agrm + mv.prox + mv.cont + 
                            mv.actv + mv.scnt + mv.gnaw + set_rescor(FALSE),
                         
                         data = data,
                         prior=c(prior(normal(0,10),  class = b),
                                 prior(normal(0,10),  class = Intercept),
                                 prior(lkj(1), class = cor),
                                 prior(cauchy(0,3),  class = sd, resp="Groominitiate"),
                                 prior(cauchy(0,3),  class = sd, resp="Proximity"),
                                 prior(cauchy(0,3),  class = sd, resp="Contactsit"),
                                 prior(cauchy(0,3),  class = sd, resp="Active"),
                                 prior(cauchy(0,3),  class = sd, resp="Scentmark"),
                                 prior(cauchy(0,3),  class = sd, resp="Gnawing")),
                         
                         warmup=3000,iter=5500,chains=4,seed=9,thin=10,
                         control=list(adapt_delta=0.99))  

#model summary
summary(mv.marmoset, prob = 0.9)

#extract posterior
mv.post<-posterior_samples(mv.marmoset)

#autocorrelation
#see str(mv.post) for parameter names
mcmc_acf(mv.post, pars = c(colnames(mv.post[40:54])))

#trait names
trait<-c("groom","proximity","contact","active","scent","gnaw")

#correlations of interest
mv.cor<-mv.post[40:54]

#lower triangle elements
#we add a column of 1's to each for the diagonal
mv.cor_grm<-cbind(rep(1,nrow(mv.cor)),mv.cor[,c(1,2,4,7,11)])
mv.cor_prox<-cbind(rep(1,nrow(mv.cor)),mv.cor[,c(3,5,8,12)])
mv.cor_cnt<-cbind(rep(1,nrow(mv.cor)),mv.cor[,c(6,9,13)])
mv.cor_act<-cbind(rep(1,nrow(mv.cor)),mv.cor[,c(10,14)])
mv.cor_sct<-cbind(rep(1,nrow(mv.cor)),mv.cor[,c(15)])

#trait names
trait<-c("groom","proximity","contact","active","scent","gnaw")

#empty correlation matrix
med.cor<-matrix(0,6,6,dimnames=list(trait,trait))

#fill in median elements
med.cor[1:6,1]<-apply(mv.cor_grm,2,median)
med.cor[2:6,2]<-apply(mv.cor_prox,2,median)
med.cor[3:6,3]<-apply(mv.cor_cnt,2,median)
med.cor[4:6,4]<-apply(mv.cor_act,2,median)
med.cor[5:6,5]<-apply(mv.cor_sct,2,median)
med.cor[6,6]<-1

#create symmetric matrix
med.cor<-forceSymmetric(med.cor,uplo="L")
med.cor<-as.matrix(med.cor)

#round for clarity
round(med.cor,digits=2)

#########model generation################################################################
#parallel analysis
fa.parallel(med.cor, n.obs = 24, fa = "fa", fm="gls", n.iter = 1000, error.bars = TRUE)
#efa
efa<-
  fa(med.cor, n.obs=24, nfactors=2, fm="gls", rotate="geominQ")
efa

#ega
ggm<- as.matrix(EBICglasso(med.cor,n=24,gamma=0.25,refit=TRUE))

round(ggm,digits=2)

#igraph
ggm.g<- graph.adjacency(abs(ggm), mode = "undirected", weighted = TRUE, diag = FALSE)

#community detection
cluster_walktrap(ggm.g)
cluster_louvain(ggm.g)


###construct models##########################################################################
#null model
null.mod <- '
Contact =~ contact
Groom =~ groom
Prox =~ proximity
Act =~ active
Scent =~ scent
Gnaw =~ gnaw'

#initial hypothesis
init.mod<-'
Sociability =~ contact + groom + proximity
Arousal =~ active + scent + gnaw'

#oblique model
obq.mod<-'
Sociability =~ contact + groom + proximity
Arousal =~ active + scent + gnaw
Sociability ~~ Arousal'

#CFA w/ residual zero-order correlation
cfa.mod<-'
Sociability =~ contact + groom + proximity
Arousal =~ active + scent + gnaw
groom ~~ active'

#EFA (0.3 threshold)
thres.mod<-'
Sociability =~ contact + groom + proximity + active
Arousal =~ active + scent + gnaw'

#ESEM
esem.mod<-'
Sociability =~ contact + groom + proximity + active + scent + -0.01*gnaw
Arousal =~ active + scent + gnaw + contact + 0.05*groom + proximity'

#convert to `lvnet` syntax
null.mod<-lav2lvnet(null.mod,data=med.cor,lavaanifyOps=list(auto.fix.first=TRUE,std.lv=TRUE))
init.mod<-lav2lvnet(init.mod,data=med.cor,lavaanifyOps=list(auto.fix.first=FALSE,std.lv=TRUE))
obq.mod<-lav2lvnet(obq.mod,data=med.cor,lavaanifyOps=list(auto.fix.first=FALSE,std.lv=TRUE))
cfa.mod<-lav2lvnet(cfa.mod,data=med.cor,lavaanifyOps=list(auto.fix.first=FALSE,std.lv=TRUE))
thres.mod<-lav2lvnet(thres.mod,data=med.cor,lavaanifyOps=list(auto.fix.first=FALSE,std.lv=TRUE))
esem.mod<-lav2lvnet(esem.mod,data=med.cor,lavaanifyOps=list(auto.fix.first=FALSE,std.lv=TRUE))

##fix theta matrix for cfa
diag(cfa.mod$theta)<-NA #set diagonal to zero

#create GNM
#duplicate model
gnm.mod<-init.mod

#add omega_theta matrix
gnm.mod$omega_theta<-matrix(0,6,6,dimnames=list(trait,trait))
gnm.mod$omega_theta[c("groom","active"),c("groom","active")]<-NA #ALG ~ ACT
diag(gnm.mod$omega_theta)<-0

#pure GGM
#add in omega_theta matrix
omega<-matrix(0,6,6,dimnames=list(trait,trait))
omega[c("groom","contact"),c("groom","contact")]<-NA #ALG ~ CNT
omega[c("groom","active"),c("groom","active")]<-NA #ALG ~ ACT
omega[c("proximity","contact"),c("proximity","contact")]<-NA #PRX ~ CNT
omega[c("active","gnaw"),c("active","gnaw")]<-NA #ACT ~ GNW 
omega[c("scent","gnaw"),c("scent","gnaw")]<-NA #SCM ~ GNW 

diag(omega)<-NA 


#####Model comparison##################################################################

#estimate model to extract fit measure labels
temp<-lvnet(data=med.cor,lambda=null.mod$lambda,ebicTuning=0.25,
            psi=null.mod$psi,sampleSize=24)

#29 fit measure labels
fit<-names(temp$fitMeasures[1:29])

#create list of posterior matrices
cor.post <-
  alply(mv.cor,1, function(x){
    
    post.cor <- matrix(0,6,6,dimnames=list(trait,trait))
    
    post.cor[1:6,1] <- as.numeric(c(1,x[c(1,2,4,7,11)]))
    post.cor[2:6,2] <- as.numeric(c(1,x[c(3,5,8,12)]))
    post.cor[3:6,3] <- as.numeric(c(1,x[c(6,9,13)]))
    post.cor[4:6,4] <- as.numeric(c(1,x[c(10,14)]))
    post.cor[5:6,5] <- as.numeric(c(1,x[15]))
    post.cor[6,6]<-1
    
    post.cor<- as.matrix(
      forceSymmetric(post.cor,uplo="L") ) 
    
    return(post.cor) 
  }) 


###null##################################################################################
null.F1<-matrix(NA,1000,6,dimnames=list(1:1000,trait))
null.F2<-matrix(NA,1000,6,dimnames=list(1:1000,trait))

#hold fit measures
null.fit<-matrix(NA,1000,29,dimnames=list(1:1000,fit))

#progress bar
pb <- txtProgressBar(min = 1, max = 1000, style = 3)

#estimate model across posterior
for(i in 1:1000){
  
  #correlation matrix sample
  cor <- cor.post[[i]]
  
  tryCatch({
    
    #run model
    mod.est <-lvnet(data = cor,lambda = null.mod$lambda,
                    psi = null.mod$psi,
                    sampleSize=24, ebicTuning = 0.25)
    
    #collect factor loadings
    null.F1[i,]<-mod.est$matrices$lambda[,1]
    null.F2[i,]<-mod.est$matrices$lambda[,2]
    
    #fit measures
    null.fit[i,]<-unlist(mod.est$fitMeasures[1:29])
    
    setTxtProgressBar(pb, i)},
    
    warning=function(w){} ) }


##init################################################################################
init.F1<-matrix(NA,1000,6,dimnames=list(1:1000,trait))
init.F2<-matrix(NA,1000,6,dimnames=list(1:1000,trait))

#hold fit measures
init.fit<-matrix(NA,1000,29,dimnames=list(1:1000,fit))

#progress bar
pb <- txtProgressBar(min = 1, max = 1000, style = 3)

#estimate model across posterior
for(i in 1:1000){
  
  #correlation matrix sample
  cor <- cor.post[[i]]
  
  tryCatch({
    
    #run model
    mod.est <-lvnet(data = cor,lambda = init.mod$lambda,
                    psi = init.mod$psi,
                    sampleSize=24, ebicTuning = 0.25)
    
    #collect factor loadings
    init.F1[i,]<-mod.est$matrices$lambda[,1]
    init.F2[i,]<-mod.est$matrices$lambda[,2]
    
    #fit measures
    init.fit[i,]<-unlist(mod.est$fitMeasures[1:29])
    
    setTxtProgressBar(pb, i)},
    
    warning=function(w){} ) }


##oblique#################################################################################
obq.F1<-matrix(NA,1000,6,dimnames=list(1:1000,trait))
obq.F2<-matrix(NA,1000,6,dimnames=list(1:1000,trait))

#hold fit measures
obq.fit<-matrix(NA,1000,29,dimnames=list(1:1000,fit))

#progress bar
pb <- txtProgressBar(min = 1, max = 1000, style = 3)

#estimate model across posterior
for(i in 1:1000){
  
  #correlation matrix sample
  cor <- cor.post[[i]]
  
  tryCatch({
    
    #run model
    mod.est <-lvnet(data = cor,lambda = obq.mod$lambda,
                    psi = obq.mod$psi,
                    sampleSize=24, ebicTuning = 0.25)
    
    #collect factor loadings
    obq.F1[i,]<-mod.est$matrices$lambda[,1]
    obq.F2[i,]<-mod.est$matrices$lambda[,2]
    
    #fit measures
    obq.fit[i,]<-unlist(mod.est$fitMeasures[1:29])
    
    setTxtProgressBar(pb, i)},
    
    warning=function(w){} ) }


###threshold model######################################################################
thres.F1<-matrix(NA,1000,6,dimnames=list(1:1000,trait))
thres.F2<-matrix(NA,1000,6,dimnames=list(1:1000,trait))

#hold fit measures
thres.fit<-matrix(NA,1000,29,dimnames=list(1:1000,fit))

#progress bar
pb <- txtProgressBar(min = 1, max = 1000, style = 3)

#estimate model across posterior
for(i in 1:1000){
  
  #correlation matrix sample
  cor <- cor.post[[i]]
  
  tryCatch({
    
    #run model
    mod.est <-lvnet(data = cor,lambda = thres.mod$lambda,
                    psi = thres.mod$psi,
                    sampleSize=24, ebicTuning = 0.25)
    
    #collect factor loadings
    thres.F1[i,]<-mod.est$matrices$lambda[,1]
    thres.F2[i,]<-mod.est$matrices$lambda[,2]
    
    #fit measures
    thres.fit[i,]<-unlist(mod.est$fitMeasures[1:29])
    
    setTxtProgressBar(pb, i)},
    
    warning=function(w){} ) }


##CFA#####################################################################################
#factor loadings
cfa.F1<-matrix(NA,1000,6,dimnames=list(1:1000,trait))
cfa.F2<-matrix(NA,1000,6,dimnames=list(1:1000,trait))

#fit measures
cfa.fit<-matrix(NA,1000,29,dimnames=list(1:1000,fit))

#residual correlation
cfa.resid<-list()


#progress bar
pb <- txtProgressBar(min = 1, max = 1000, style = 3)

#estimate model across posterior
for(i in 1:1000){
  
  #correlation matrix sample
  cor <- cor.post[[i]]
  
  tryCatch({
    
    #run model
    #theta matrix estimates residual correlation
    mod.est <-lvnet(data = cor,lambda = cfa.mod$lambda,
                    psi = cfa.mod$psi,theta = cfa.mod$theta,
                    sampleSize=24, ebicTuning = 0.25)
    
    #collect factor loadings
    cfa.F1[i,]<-mod.est$matrices$lambda[,1]
    cfa.F2[i,]<-mod.est$matrices$lambda[,2]
    
    #fit measures
    cfa.fit[i,]<-unlist(mod.est$fitMeasures[1:29])
    
    #residual correlation
    cfa.resid[[i]]<-matrix(mod.est$matrices$theta,6,6,dim=list(trait,trait))
    
    
    setTxtProgressBar(pb, i)},
    
    warning=function(w){} ) }


##ESEM######################################################################################

#grab starting values
esem.start<-
matrix(matrix(efa$loadings,1),6,2)

esem.F1<-matrix(NA,1000,6,dimnames=list(1:1000,trait))
esem.F2<-matrix(NA,1000,6,dimnames=list(1:1000,trait))

#hold fit measures
esem.fit<-matrix(NA,1000,29,dimnames=list(1:1000,fit))

#progress bar
pb <- txtProgressBar(min = 1, max = 1000, style = 3)

#estimate model across posterior
for(i in 1:1000){
  
  #correlation matrix sample
  cor <- cor.post[[i]]
  
  tryCatch({
    
    #run model
    mod.est <-lvnet(data = cor,lambda = esem.mod$lambda,
                    psi = esem.mod$psi,startValues = list(esem.start), 
                    sampleSize=24, ebicTuning = 0.25)
    
    #collect factor loadings
    esem.F1[i,]<-mod.est$matrices$lambda[,1]
    esem.F2[i,]<-mod.est$matrices$lambda[,2]
    
    #fit measures
    esem.fit[i,]<-unlist(mod.est$fitMeasures[1:29])

    setTxtProgressBar(pb, i)},
    
    warning=function(w){} ) }


##GGM###################################################################################
#capture parameter estimates
ggm.est<-list()

#capture EBIC values
ggm.fit<-matrix(NA,1000,1,dimnames=list(1:1000,"EBIC"))

#progress bar
pb <- txtProgressBar(min = 1, max = 1000, style = 3)

#estimate model across posterior
for(i in 1:1000){
  
  covMat<-cor.post[[i]]
  
  #estimate precision
  #zero entries are fixed by the zeroes in 'omega'
  glassoRes <- suppressWarnings(
    glasso::glasso(covMat, 0, zero = which(omega == 0, arr.ind = TRUE)))
  
  #precision
  invSigma <- glassoRes$wi
  
  #GGM
  p.cor <- matrix(qgraph::wi2net(invSigma), 6, 6, dimnames=list(trait,trait))
  ggm.est[[i]] <- p.cor
  
  #EBIC
  ggm.fit[i,]<- ggmFit(p.cor, covMat, sampleSize=24, ebicTuning=0.25)$fitMeasures$ebic
  
  setTxtProgressBar(pb, i)      
  
}


##GNM############################################################################################

#hold factor loadings
gnm.F1<-matrix(NA,1000,6,dimnames=list(1:1000,trait))
gnm.F2<-matrix(NA,1000,6,dimnames=list(1:1000,trait))

#hold fit measures
gnm.fit<-matrix(NA,1000,29,dimnames=list(1:1000,fit))

#residual GGM posterior
gnm.resid<-list()

#progress bar
pb <- txtProgressBar(min = 1, max = 1000, style = 3)

#estimate model across posterior
for(i in 1:1000){
  
  #correlation matrix sample
  cor <- cor.post[[i]]
  
  tryCatch({
    
    #run model
    mod.est <-lvnet(data = cor,lambda = gnm.mod$lambda,
                    psi = gnm.mod$psi, omega_theta = gnm.mod$omega_theta,
                    sampleSize=24, ebicTuning=0.25)
    
    #collect factor loadings
    gnm.F1[i,]<-mod.est$matrices$lambda[,1]
    gnm.F2[i,]<-mod.est$matrices$lambda[,2]
    
    #fit measures
    gnm.fit[i,]<-unlist(mod.est$fitMeasures[1:29])
    
    #residual GGM
    gnm.resid[[i]]<-matrix(mod.est$matrices$omega_theta,6,6,dim=list(trait,trait))
    
    setTxtProgressBar(pb, i)},
    
    warning=function(w){} ) }


###model comparison############################################################################

#M2
init.loadings<-cbind(init.F1[,1:3],init.F2[,4:6]) #combine factor loadings

#collect posterior rows containing Heywood cases
init.heywood<-
  rownames(init.loadings[apply(init.loadings, MARGIN = 1, function(x) any(x >= 1 | x <= -1)), ])

#proportion of admissible solution
1000 - length(init.heywood)

#M3
obq.loadings<-cbind(obq.F1[,1:3],obq.F2[,4:6]) 
obq.heywood<-
  rownames(obq.loadings[apply(obq.loadings, MARGIN = 1, function(x) any(x >= 1 | x <= -1)), ])
1000 - length(obq.heywood)

#M4
cfa.loadings<-cbind(cfa.F1[,1:3],cfa.F2[,4:6]) 
cfa.heywood<-
  rownames(cfa.loadings[apply(cfa.loadings, MARGIN = 1, function(x) any(x >= 1 | x <= -1)), ])
1000 - length(cfa.heywood)

#M5
thres.loadings<-cbind(thres.F1[,1:4],thres.F2[,4:6]) #include activy cross-loading
thres.heywood<-
  rownames(thres.loadings[apply(thres.loadings, MARGIN = 1, function(x) any(x >= 1 | x <= -1)), ])
1000 - length(thres.heywood)

#M6
esem.loadings<-cbind(esem.F1[,1:6],esem.F2[,1:6]) #fully estimated factor matrix
esem.heywood<-
  rownames(esem.loadings[apply(esem.loadings, MARGIN = 1, function(x) any(x >= 1 | x <= -1)), ])
1000 - length(esem.heywood)

#M8
gnm.loadings<-cbind(gnm.F1[,1:3],gnm.F2[,4:6]) 
gnm.heywood<-
  rownames(gnm.loadings[apply(gnm.loadings, MARGIN = 1, function(x) any(x >= 1 | x <= -1)), ])
1000 - length(gnm.heywood)

#ebic vectors
null.fit<-data.frame(null.fit[,"ebic"])
init.fit<-data.frame(init.fit[,"ebic"])
obq.fit<-data.frame(obq.fit[,"ebic"])
cfa.fit<-data.frame(cfa.fit[,"ebic"])
thres.fit<-data.frame(thres.fit[,"ebic"])
esem.fit<-data.frame(esem.fit[,"ebic"])
ggm.fit<-ggm.fit
gnm.fit<-data.frame(gnm.fit[,"ebic"])

#remove heywood cases
#M1
null.ebic <- median(null.fit[,])

#M2
init.ebic<- median(init.fit[setdiff(1:1000,init.heywood),])

#M3
obq.ebic<-  median(obq.fit[setdiff(1:1000,obq.heywood),])

#M4
cfa.ebic<-  median(cfa.fit[setdiff(1:1000,cfa.heywood),])

#M5
thres.ebic<-median(thres.fit[setdiff(1:1000,thres.heywood),])

#M6
esem.ebic<- median(esem.fit[setdiff(1:1000,esem.heywood),])

#M7
ggm.ebic <- median(ggm.fit[,])

#M8
gnm.ebic<-  median(gnm.fit[setdiff(1:1000,gnm.heywood),])


#minimum ebic
c(null.ebic, init.ebic, obq.ebic, cfa.ebic, thres.ebic, esem.ebic, ggm.ebic, gnm.ebic)

#gnm (M8) is EBIC_min
c(null.ebic, init.ebic, obq.ebic, cfa.ebic, thres.ebic, esem.ebic, ggm.ebic) - gnm.ebic


###model results######################################################################
#M7--sparse GGM
#median
apply(simplify2array(ggm.est), c(1,2), median)

#p>0
apply(simplify2array(ggm.est), c(1,2), function(x) sum(x>0)/1000 )

#M8--GNM
#remove Heywood cases
gnm.adm.loadings <- gnm.loadings[setdiff(1:1000,gnm.heywood),]
gnm.adm.resid <- (simplify2array(gnm.resid))[,,setdiff(1:1000,gnm.heywood)]

#median
apply(gnm.adm.loadings,2, median)
apply(gnm.adm.resid, c(1,2), median)

#p>0
apply(gnm.adm.loadings, 2, function(x) sum(x>0)/(1000-length(gnm.heywood)))
sum(gnm.adm.resid[1,4,] > 0)/(1000-length(gnm.heywood)) #ALG ~~ ACT
