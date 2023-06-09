library(caret)
library(Matrix)
library(xgboost)
#read in data
read.csv("train.csv")->train
read.csv("test.csv")->test
#read.csv("Ftrain.csv")->train
#prep data
train.ID<-train$Sample
train.ID<-as.character(train.ID)
train$Sample<-NULL
test.ID<-test$Sample
test.ID<-as.character(test.ID)
test$Sample<-NULL


train$Seq <- as.integer(train$Seq)
train$Target<-train$Seq-1
Target<-train$Target
train$Seq<-NULL
test$Seq<-NULL


test$Target<-0
test$test<-1
train$test<-0
all.data<-rbind(train,test)

dummies <- dummyVars(Target ~ ., data = all.data)
 all.data<-predict(dummies, newdata = all.data)
all.data<- data.frame(all.data)
test<-subset(all.data,test==1)
train<-subset(all.data,test==0)

train$test<-NULL
test$test<-NULL

train[train== -9]<- NA
test[test== -9]<- NA



dtrain<- xgb.DMatrix(data=as.matrix(train),label=Target,missing=NA)
dtest<- xgb.DMatrix(data=as.matrix(test),missing=NA)
param <- list(  objective           = "binary:logistic", 
				        gamma                 =0.01,
                booster             = "gbtree",
                eval_metric         = "logloss",
                eta                 = 0.01,
                max_depth           = 4,
                subsample           = 1,
                colsample_bytree    = 0.3
)

CV<-xgb.cv(params=param,nrounds=1000,nfold=53,missing=NA,prediction=T,data=dtrain,label=Target)
XGBm<-xgb.train( params=param,nrounds=1000,missing=NA,data=dtrain,label=Target)   
Prob<-predict(XGBm,newdata=dtest)
Prob_cv<-CV$pred
mCV<-cbind(train.ID,Prob_cv)
mMod<-cbind(test.ID,Prob)
Model<-rbind(mCV,mMod)

#importance_matrix <- xgb.importance(names(train), model = XGBm,label=Target, data=dtrain)
importance_matrix <- xgb.importance(names(train), model = XGBm)
 xgb.plot.importance(importance_matrix)
 Model<-data.frame(Model)

Model$Prob_cv<-as.numeric(levels(Model$Prob_cv))[Model$Prob_cv]

names(Model)<-c("Sample","pS.trenchii")
Model$pS.glynnii<-1-Model$pS.trenchii
Model<-Model[!duplicated(Model$Sample),]
mModel<-melt(Model)
names(mModel)<-c("Sample","Species","Probability")
read.csv("PlotOrder.csv")->order


mModel$Sample <- factor(mModel$Sample, levels = mModel$Sample[match(order$Sample,mModel$Sample)])
mModel[!mModel$Sample=="NA",]->mModel
ggplot(mModel ,aes(x=Sample,y=Probability,fill=Species))+geom_bar(stat="identity",position="stack")+ scale_fill_brewer(palette="Spectral")+theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour="black"))

read.csv("train.csv")->train
read.csv("test.csv")->test
rbind(train,test)->data
data<-data[match(order[,1],data$Sample),]
data$Sample<-factor(data$Sample,levels=data$Sample)
 ITS<-data[,c('Sample','ITS_D1','ITS_D4','ITS_D6')]
ITS<-transform(ITS,ITS=paste0(ITS_D1,ITS_D4,ITS_D6))[,c('Sample','ITS')]

ggplot(ITS,aes(x=Sample,fill=as.factor(ITS)))+geom_bar()+theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour="black"))
ggplot(data,aes(x=Sample,y=HetroLoci))+geom_bar(stat="identity")+theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour="black"))
ggplot(data,aes(x=Sample,fill=as.factor(Host)))+geom_bar()+theme(axis.text.x = element_text(size=5, angle = 90, hjust = 1,colour="black"))
