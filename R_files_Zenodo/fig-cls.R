library(praznik)
library(rFerns)
library(randomForest)
library(kernlab)
library(e1071)
library(parallel)
library(ggplot2)

RNGkind("L'Ecuyer-CMRG")

auroc<-function(score,cls)
 mean(rank(score)[cls]-1:sum(cls))/sum(!cls)

readRDS('morphine.RDS')->Q
Q[,1:90]->X
kTransform(X)->Xk
kTransform(Q$yUSV)->Yk

cvMasks<-function(n)
 sapply(1:n,function(i){
  rep(TRUE,n)->ans
  ans[i]<-FALSE
  ans
 })

cvKtMasks<-function(n){
 cvMasks(n)->p
 p[!p]<-NA
 !is.na(as.matrix(kTransform(data.frame(p))))
}

ktMods<-list(
 kt_rforest_prob=function(X,Y,Xp,...) predict(randomForest(X,Y,...),Xp,type="prob")[,">"],
 kt_rforest_cat=function(X,Y,Xp,...) predict(randomForest(X,Y,...),Xp),
 kt_nb_cat=function(X,Y,Xp,...) predict(naiveBayes(X,Y,...),Xp),
 kt_nb_prob=function(X,Y,Xp,...) predict(naiveBayes(X,Y,...),Xp,type="raw")[,">"],
 kt_svm_prob=function(X,Y,Xp,...) predict(ksvm(Y~.,data=X,prob.model=TRUE),Xp,type="prob")[,">"],
 kt_svm_cat=function(X,Y,Xp,...) predict(ksvm(Y~.,data=X),Xp)
)

cvkt<-function(x,y,model=ktMods$rfe,...){
 cvKtMasks(nrow(x))->M
 M<-M[,sample.int(ncol(M))]
 X<-kTransform(x)
 Y<-kTransform(if(!is.factor(y)) y else y!=levels(y)[1])

 Yp<-mclapply(1:ncol(M),function(fold){
  cM<-M[,fold]
  model(X[cM,],Y[cM],X[!cM,],...)->ans
  data.frame(idx=which(!cM),p=ans)
 },mc.cores=detectCores(),mc.preschedule=FALSE)
 do.call(rbind,Yp)->Yp
 #Order by idx with random duplicate breaking
 Yp[order(Yp$idx,runif(nrow(Yp))),]->Yp
 Yp$p[!duplicated(Yp$idx)]->Yp
 kInverse(Yp)->yp
 if(!is.factor(y)) 
  c(cor=cor(rank(y),rank(yp)))
 else
  c(auroc=auroc(rank(yp),y==levels(y)[2]))
}

rgMods<-list(
 rg_rforest=function(X,Y,Xp,...) 
  if(!is.factor(Y))
   predict(randomForest(X,Y,...),Xp) else 
   predict(randomForest(X,Y,...),Xp,type="prob")[,levels(Y)[2]],
 rg_svm=function(X,Y,Xp,...) 
  if(!is.factor(Y))
   predict(ksvm(Y~.,data=X),Xp) else
   predict(ksvm(Y~.,data=X,prob.model=TRUE),Xp,type="prob")[,levels(Y)[2]]
)

cv<-function(x,y,model=function(X,Y,Xp,...) predict(randomForest(X,Y,...),Xp),...){
 cvMasks(nrow(x))->m
 yp<-mclapply(1:ncol(m),function(fold){
  cm<-m[,fold]
  model(x[cm,],y[cm],x[!cm,],...)->ans
  data.frame(idx=which(!cm),p=ans)
 },mc.cores=detectCores(),mc.preschedule=FALSE)

 do.call(rbind,yp)->yp
 stopifnot(all(!duplicated(yp$idx)))
 yp$p[order(yp$idx)]->yp

 if(!is.factor(y))
  c(cor=cor(rank(y),rank(yp)))
 else
  c(auroc=auroc(rank(yp),y==levels(y)[2]))
}

doOne<-function(seed=17,yn="yUSV"){
 set.seed(seed)
 lapply(names(ktMods),function(ktmod){
  message(ktmod)
  mod<-ktMods[[ktmod]]
  cvkt(Q[,1:90],Q[,yn],mod)->sco
  data.frame(seed=seed,mod=ktmod,kt=TRUE,score=sco,yn=yn)
 })->kt
 lapply(names(rgMods),function(rgmod){
  message(rgmod)
  mod<-rgMods[[rgmod]]
  cv(Q[,1:90],Q[,yn],mod)->sco
  data.frame(seed=seed,mod=rgmod,kt=FALSE,score=sco,yn=yn)
 })->rg
 do.call(rbind,c(kt,rg))
}

#This is a lengthy calculation, henceforth we use file cache to 
# make the process idempotent. Remove files with *cache* in name
# to force full regeneration.

calculateAll<-function(){
 for(seed in 1:30){
  for(yn in c("yMorph","yUSV","yWithdrawal")){
   fN<-sprintf("fig-cls-cache-part-%s-%d.RDS",yn,seed)
   if(!file.exists(fN)){
    try(doOne(seed,yn))->Z
    saveRDS(Z,fN)
   }else{
    message("Found partial cache ",fN,", skipping regeneration")
   }
  }
 }
}

collectAll<-function()
 do.call(rbind,lapply(list.files(patt='fig-cls-cache-part-.*RDS$'),readRDS))

fetchAll<-function(){
 if(file.exists('fig-cls-cache.RDS')) return(readRDS('fig-cls-cache.RDS'))
 message("No cache found, recreating (this will take same time)")
 calculateAll()
 collectAll()->Ans
 saveRDS(Ans,'fig-cls-cache.RDS')
 unlink(list.files(patt='^fig-cls-cache-part-.*RDS$'))
 return(Ans)
}

#Plotting

enhance<-function(Q){
 grepl("prob$",Q$mod)->Q$kts
 ifelse(Q$kt,ifelse(Q$kts,"kts","ktv"),ifelse(grepl("rank",as.character(Q$mod)),"rank","dir"))->Q$kind
 sapply(strsplit(as.character(Q$mod),split='_'),'[',2)->Q$model

 Q$Model<-factor(as.character(Q$model),c("nb","rforest","svm"),c("Naive Bayes","Random Forest","SVM"))
 Q$Prediction<-factor(as.character(Q$kind),c("dir","ktv","kts"),c("Direct","KT class","KT score"))
 Q$Case<-factor(as.character(Q$yn),c("yUSV","yMorph","yWithdrawal"),c("USV [Spearman cc.]","Morphine [AUROC]","Withdrawal [AUROC]"))
 Q$Value<-Q$score

 Q
}

plotAll<-function(Q)
 ggplot(Q,aes(y=Value,fill=Model,x=Prediction))+geom_boxplot()+facet_wrap(~Case,scales="free_y")+theme(legend.position="bottom")

makeFigure<-function(Q)
 ggsave("fig-cls.pdf",device=cairo_pdf,width=8,height=5,plot=plotAll(enhance(Q)))

if(!interactive()) makeFigure(fetchAll())
