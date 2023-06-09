OS='W'
##############################################################################
##################### Analisis communities in time ###########################
##############################################################################
iter=5000;cores=1
library(ineq);library(parallel);library(textstem);library(tm);library(parallel)
CUs=c('Deforestation', 'Ecological Networks','Invasive Species', 'Trade')
x=lapply(c(1:4),function(File){CU=c('Deforestation', 'Ecological Networks','Invasive Species', 'Trade')[File]
if(OS=='M'){setwd(paste("/Users/felberarroyave/OneDrive - University of California Merced/PhD/research/Neglected/analysis",CU,sep="/"))}else{setwd(paste("C:/Files/OneDrive - University of California Merced/PhD/research/Neglected/analysis",CU,sep="/"))}
as.matrix(read.csv('BaseShort.csv'))})
if(OS=='M'){setwd(paste("/Users/felberarroyave/OneDrive - University of California Merced/PhD/research/Neglected/analysis",sep="/"))}else{setwd(paste("C:/Files/OneDrive - University of California Merced/PhD/research/Neglected/analysis",sep="/"))}
entrp=function(p){p=na.omit(p)/sum(na.omit(p));sum(-p*log(p))/sum(-rep(1/length(p),length(p))*log(rep(1/length(p),length(p))))}
#contry=as.matrix(read.csv('countries.csv'))
nams=c('VarCom','VarSc','VarCon','VarAuth','VarKW','VarCoAu','VarCoKW','VarSO','VarNS')
for(z in 1:length(x)){x[[z]]=x[[z]][!is.na(x[[z]][,'PY']),]};x1=x
###################
for (p in 9){
if(nams[p]=='VarCom'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){paste(z,x[[z]][,'PY'],x[[z]][,'UT'],x[[z]][,'COM'],sep='#')})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})}
#
if(nams[p]=='VarSc'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){paste(rep(z,length(unlist(strsplit(x[[z]][i,'SC'],'; ')))),rep(x[[z]][i,'PY'],length(unlist(strsplit(x[[z]][i,'SC'],'; ')))),rep(x[[z]][i,'UT'],length(unlist(strsplit(x[[z]][i,'SC'],'; ')))),unlist(strsplit(x[[z]][i,'SC'],'; ')),sep='#')}))})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})}
#
if(nams[p]=='VarSO'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){paste(z,x[[z]][,'PY'],x[[z]][,'UT'],x[[z]][,'SO'],sep='#')})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})}
#
if(nams[p]=='VarCon'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){paste(rep(z,length(unlist(strsplit(x[[z]][i,'CU'],'#')))),rep(x[[z]][i,'PY'],length(unlist(strsplit(x[[z]][i,'CU'],'#')))),rep(x[[z]][i,'UT'],length(unlist(strsplit(x[[z]][i,'CU'],'#')))),unlist(strsplit(x[[z]][i,'CU'],'#')),sep='#')}))})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})}
#
if(nams[p]=='VarNS'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){paste(rep(z,length(unlist(strsplit(x[[z]][i,'NS'],'#')))),rep(x[[z]][i,'PY'],length(unlist(strsplit(x[[z]][i,'NS'],'#')))),rep(x[[z]][i,'UT'],length(unlist(strsplit(x[[z]][i,'NS'],'#')))),unlist(strsplit(x[[z]][i,'NS'],'#')),sep='#')}))})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})}
#
if(nams[p]=='VarAuth'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){paste(rep(z,length(unlist(strsplit(x[[z]][i,'AU'],'##')))),rep(x[[z]][i,'PY'],length(unlist(strsplit(x[[z]][i,'AU'],'##')))),rep(x[[z]][i,'UT'],length(unlist(strsplit(x[[z]][i,'AU'],'##')))),
unlist(lapply(unlist(strsplit(x[[z]][i,'AU'],'##')),function(k){unlist(strsplit(k,','))[1]})),sep='#')}))})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y1=y1[y1[,4]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})
for(z in 1:length(y)){y[[z]]=y[[z]][rowSums(ifelse(y[[z]][,]!=0,1,0))>2,]}}
#
if(nams[p]=='VarKW'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){paste(rep(z,length(unlist(strsplit(gsub('-',' ',unlist(strsplit(x[[z]][i,'ID'],'; '))),' ')))),rep(x[[z]][i,'PY'],length(unlist(strsplit(gsub('-',' ',unlist(strsplit(x[[z]][i,'ID'],'; '))),' ')))),rep(x[[z]][i,'UT'],length(unlist(strsplit(gsub('-',' ',unlist(strsplit(x[[z]][i,'ID'],'; '))),' ')))),
unlist(lapply(unlist(strsplit(gsub('-',' ',unlist(strsplit(x[[z]][i,'ID'],'; '))),' ')),function(k){unlist(strsplit(k,','))[1]})),sep='#')}))})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y1=y1[y1[,4]!="NA",];y1[,4]=unlist(lapply(tolower(y1[,4]),function(z){stem_words(z)}))
y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})
for(z in 1:length(y)){y[[z]]=y[[z]][rowSums(ifelse(y[[z]][,]!=0,1,0))>2,]}}
#
if(nams[p]=='VarCoAu'){y1=unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){nam=na.omit(sort(unlist(lapply(unlist(strsplit(x[[z]][i,'AU'],'##')),function(k){unlist(strsplit(k,', '))[1]}))))
if(length(nam)>1){paste(z,x[[z]][i,'PY'],x[[z]][i,'UT'],unlist(lapply(1:(length(nam)-1),function(k){paste(nam[k],nam[(k+1):length(nam)],sep='_')})),sep='#')}}))}))
y1=matrix(unlist(strsplit(y1,'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})
for(z in 1:length(y)){y[[z]]=y[[z]][rowSums(ifelse(y[[z]][,]!=0,1,0))>2,]}}
#
if(nams[p]=='VarCoKW'){y1=unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){
nam=stem_words(tolower(sort(unlist(strsplit(unlist(strsplit(x[[z]][i,'ID'],'; ')),' ')))))
if(length(nam)>1){paste(z,x[[z]][i,'PY'],x[[z]][i,'UT'],unlist(lapply(1:(length(nam)-1),function(k){paste(nam[k],nam[(k+1):length(nam)],sep='_')})),sep='#')}}))}))
y1=matrix(unlist(strsplit(y1,'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})
for(z in 1:length(y)){y[[z]]=y[[z]][rowSums(ifelse(y[[z]][,]!=0,1,0))>2,]}}
##
##
##
for(z in 1:length(y)){y[[z]][y[[z]]==0]=NA}
m=unlist(lapply(1:length(y),function(z){paste(z,colnames(y[[z]]),apply(y[[z]],2,entrp),apply(y[[z]],2,ineq),sep='#')}))
m=matrix(unlist(strsplit(m,'#')),byrow=T,ncol=4)
colnames(m)=c('Level','Year','Entrop','GINI')
write.csv(m,paste(nams[p],'.csv',sep=''),row.names=F)
#
for(z in 1:length(y)){y[[z]][is.na(y[[z]])]=0
y[[z]]=t(apply(y[[z]],1,cumsum))
y[[z]][y[[z]]==0]=NA}
m=unlist(lapply(1:length(y),function(z){paste(z,colnames(y[[z]]),apply(y[[z]],2,entrp),apply(y[[z]],2,ineq),sep='#')}))
m=matrix(unlist(strsplit(m,'#')),byrow=T,ncol=4)
colnames(m)=c('Level','Year','Entrop','GINI')
write.csv(m,paste(nams[p],'_C','.csv',sep=''),row.names=F)}

##########################################################################
##########################################################################
###########. NULL MODELS.  ###############################################
##########################################################################
##########################################################################

for (p in 9){
m=unlist(mclapply(1:iter,function(k){
for(z in 1:length(x)){x[[z]][,'PY']=sample(x1[[z]][,'PY'])}
if(nams[p]=='VarCom'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){paste(z,x[[z]][,'PY'],x[[z]][,'UT'],x[[z]][,'COM'],sep='#')})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})}
#
if(nams[p]=='VarSc'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){paste(rep(z,length(unlist(strsplit(x[[z]][i,'SC'],'; ')))),rep(x[[z]][i,'PY'],length(unlist(strsplit(x[[z]][i,'SC'],'; ')))),rep(x[[z]][i,'UT'],length(unlist(strsplit(x[[z]][i,'SC'],'; ')))),unlist(strsplit(x[[z]][i,'SC'],'; ')),sep='#')}))})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})}
#
if(nams[p]=='VarSO'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){paste(z,x[[z]][,'PY'],x[[z]][,'UT'],x[[z]][,'SO'],sep='#')})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})}
#
if(nams[p]=='VarCon'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){paste(rep(z,length(unlist(strsplit(x[[z]][i,'CU'],'#')))),rep(x[[z]][i,'PY'],length(unlist(strsplit(x[[z]][i,'CU'],'#')))),rep(x[[z]][i,'UT'],length(unlist(strsplit(x[[z]][i,'CU'],'#')))),unlist(strsplit(x[[z]][i,'CU'],'#')),sep='#')}))})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})}
#
if(nams[p]=='VarNS'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){paste(rep(z,length(unlist(strsplit(x[[z]][i,'NS'],'#')))),rep(x[[z]][i,'PY'],length(unlist(strsplit(x[[z]][i,'NS'],'#')))),rep(x[[z]][i,'UT'],length(unlist(strsplit(x[[z]][i,'NS'],'#')))),unlist(strsplit(x[[z]][i,'NS'],'#')),sep='#')}))})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})}
#
if(nams[p]=='VarAuth'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){paste(rep(z,length(unlist(strsplit(x[[z]][i,'AU'],'##')))),rep(x[[z]][i,'PY'],length(unlist(strsplit(x[[z]][i,'AU'],'##')))),rep(x[[z]][i,'UT'],length(unlist(strsplit(x[[z]][i,'AU'],'##')))),
unlist(lapply(unlist(strsplit(x[[z]][i,'AU'],'##')),function(k){unlist(strsplit(k,','))[1]})),sep='#')}))})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y1=y1[y1[,4]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})
for(z in 1:length(y)){y[[z]]=y[[z]][rowSums(ifelse(y[[z]][,]!=0,1,0))>2,]}}
#
if(nams[p]=='VarKW'){y1=matrix(unlist(strsplit(unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){paste(rep(z,length(unlist(strsplit(gsub('-',' ',unlist(strsplit(x[[z]][i,'ID'],'; '))),' ')))),rep(x[[z]][i,'PY'],length(unlist(strsplit(gsub('-',' ',unlist(strsplit(x[[z]][i,'ID'],'; '))),' ')))),rep(x[[z]][i,'UT'],length(unlist(strsplit(gsub('-',' ',unlist(strsplit(x[[z]][i,'ID'],'; '))),' ')))),
unlist(lapply(unlist(strsplit(gsub('-',' ',unlist(strsplit(x[[z]][i,'ID'],'; '))),' ')),function(k){unlist(strsplit(k,','))[1]})),sep='#')}))})),'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y1=y1[y1[,4]!="NA",];y1[,4]=unlist(lapply(tolower(y1[,4]),function(z){stem_words(z)}))
y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})
for(z in 1:length(y)){y[[z]]=y[[z]][rowSums(ifelse(y[[z]][,]!=0,1,0))>2,]}}
#
if(nams[p]=='VarCoAu'){y1=unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){nam=na.omit(sort(unlist(lapply(unlist(strsplit(x[[z]][i,'AU'],'##')),function(k){unlist(strsplit(k,', '))[1]}))))
if(length(nam)>1){paste(z,x[[z]][i,'PY'],x[[z]][i,'UT'],unlist(lapply(1:(length(nam)-1),function(k){paste(nam[k],nam[(k+1):length(nam)],sep='_')})),sep='#')}}))}))
y1=matrix(unlist(strsplit(y1,'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})
for(z in 1:length(y)){y[[z]]=y[[z]][rowSums(ifelse(y[[z]][,]!=0,1,0))>2,]}}
#
if(nams[p]=='VarCoKW'){y1=unlist(lapply(1:length(x),function(z){unlist(lapply(1:dim(x[[z]])[1],function(i){
nam=stem_words(tolower(sort(unlist(strsplit(unlist(strsplit(x[[z]][i,'ID'],'; ')),' ')))))
if(length(nam)>1){paste(z,x[[z]][i,'PY'],x[[z]][i,'UT'],unlist(lapply(1:(length(nam)-1),function(k){paste(nam[k],nam[(k+1):length(nam)],sep='_')})),sep='#')}}))}))
y1=matrix(unlist(strsplit(y1,'#')),byrow=T,ncol=4)
y1=y1[y1[,2]!="NA",];y=lapply(1:length(x),function(z){table(y1[y1[,1]==z,4],y1[y1[,1]==z,2])})
for(z in 1:length(y)){y[[z]]=y[[z]][rowSums(ifelse(y[[z]][,]!=0,1,0))>2,]}}
##
##
##
for(z in 1:length(y)){y[[z]][y[[z]]==0]=NA}
m=unlist(lapply(1:length(y),function(z){paste('N',k,z,colnames(y[[z]]),apply(y[[z]],2,entrp),apply(y[[z]],2,ineq),sep='#')}))
for(z in 1:length(y)){y[[z]][is.na(y[[z]])]=0
y[[z]]=t(apply(y[[z]],1,cumsum))
y[[z]][y[[z]]==0]=NA}
c(m,unlist(lapply(1:length(y),function(z){paste('C',k,z,colnames(y[[z]]),apply(y[[z]],2,entrp),apply(y[[z]],2,ineq),sep='#')})))},mc.cores=cores))
m=matrix(unlist(strsplit(m,'#')),byrow=T,ncol=6)
colnames(m)=c('Type','Iter','Level','Year','Entrop','GINI')
write.csv(m,paste(nams[p],'_N','.csv',sep=''),row.names=F)}





((20+20+20+16.75+20+15+20)/7)/20
((15+15+10+15+15+10+15)/7)/15