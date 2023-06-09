### Multivariate model yield predindviion 2016 with ndvi
### Running stepwise regression model with ndvi
rm(list=ls())
cat('\\f')
setwd("~/Documents/BHEARD_documents/Dissertation_research/Data_Analysis")
require(pastecs) 
require("lme4")
require(ggplot2)
library(leaps)
library(caret) 
library(dplyr)
BLUE2016=read.csv("BLUE_2016.csv")
ndvi_data=BLUE2016[,c(1:2,11:19,28)]
stepwise.ndvi.2016 <- lm(GRYLD ~ NDVI_20160121+NDVI_20160130+NDVI_20160203+NDVI_20160207+NDVI_20160223+NDVI_20160228+NDVI_20160303+NDVI_20160310+NDVI_20160315
                         ,data=ndvi_data)  ## writing the model with ndvi
stepwise.ndvi.2016 <- step(stepwise.ndvi.2016, direndviion = "both") ## running the stepwise regressing model
formula(stepwise.ndvi.2016) ## getting the selendvied variables from the stepwise regression model

### Predindvied yield from stepwise regression model with ndvi
step.predyield.ndvi.2016=NULL 
for(i in 1:10){
        step.prediction.ndvi.2016=ndvi_data[ndvi_data$trial==i,] #get just one trial to predindvi
        step.training.ndvi.2016=ndvi_data[ndvi_data$trial!=i,] #get all other trials to train on
        fit.ndvi.2016 = lm(formula = formula(stepwise.ndvi.2016), data =step.training.ndvi.2016) ## direndvily import variables into the model from stepwise regression 
        step.predvalue.ndvi.2016=predict(object=fit.ndvi.2016, newdata=step.prediction.ndvi.2016) ## Predindvied value of individual trial
        predicted.ndvi.trial2016=data.frame(entry=step.prediction.ndvi.2016$plot,step.predvalue.ndvi.2016) ## Predindvied value of individual trial in data frame
        step.predyield.ndvi.2016=rbind(step.predyield.ndvi.2016, predicted.ndvi.trial2016) ## Combinig all trial's predindvied value into one frame
}
step.predyield.ndvi.2016 = data.frame(trial=rep(1:10,each=60),step.predyield.ndvi.2016)
step.predyield.ndvi.2016 <- merge(step.predyield.ndvi.2016, ndvi_data[, c(1:2,12)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extrandviing predindviion accuracy from stepwise regression model with ndvi
cor.ndvi.step2016 <- NULL #start dataframe to hold correlation values
for(i in 1:10){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ndvi.2016[step.predyield.ndvi.2016$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ndvi.2016, trial_data$GRYLD) #get correlations
        cor.ndvi.step2016 <- c(cor.ndvi.step2016, trial_cor) #write results outside of the for loop
}

### Predindviion from shrinkage models
predindvi.ndvi.las2016=c()  #Intialize to capture predindvied values from LASSO model
predindvi.ndvi.rdg2016=c() #Intialize to capture predindvied values from Ridge model
predindvi.ndvi.enet2016=c() #Intialize to capture predindvied values from ElasticNet model
cor.ndvi.las2016=c() #Intialize to capture correlatio from LASSO model
cor.ndvi.rdg2016=c() #Intialize to capture correlatio from Ridge model
cor.ndvi.enet2016=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:10){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ndvi.data2016 = ndvi_data[ndvi_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ndvi.data2016 = ndvi_data[ndvi_data$trial!=trial,] #sets up predindviion populations
        rndvirl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ndvi.model2016 = train(GRYLD ~ ., data = train.ndvi.data2016[,-c(1,2)], #Looks like we are fitting a lasso model.
                                     method = "lasso",
                                     trControl = rndvirl,
                                     preProc = c("center", "scale"))
        ridge.ndvi.model2016 = train(GRYLD ~ ., data = train.ndvi.data2016[,-c(1,2)], #looks like we are fitting a ridge model.
                                     method = "ridge",
                                     trControl = rndvirl,
                                     preProc = c("center", "scale"))
        enet.ndvi.model2016 = train(GRYLD ~ ., data = train.ndvi.data2016[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                    method = "enet",
                                    trControl = rndvirl,
                                    preProc = c("center", "scale"))
        #model coefficients
        coef.ndvi.lasso2016 = lasso.ndvi.model2016$finalModel #extrandviing some coefficients from LASSO model
        coef.ndvi.ridge2016 = ridge.ndvi.model2016$finalModel #extrandviing some coefficients from Ridge model
        coef.ndvi.enet2016 = enet.ndvi.model2016$finalModel #extrandviing some coefficients from ElasticNet model
        #predindviion
        pred.ndvi.lasso2016 = predict(lasso.ndvi.model2016,test.ndvi.data2016[,-c(1,2)]) #making a predindvi statement.  This looks promising
        pred.ndvi.ridge2016 = predict(ridge.ndvi.model2016,test.ndvi.data2016[,-c(1,2)]) #appears the same with other models
        pred.ndvi.elasticnet2016 = predict(enet.ndvi.model2016,test.ndvi.data2016[,-c(1,2)])
        #extrandviing correlations
        cor.ndvi.l2016=cor(ndvi_data[ndvi_data$trial==trial,12],pred.ndvi.lasso2016) # Correlaiton between predindvied value from LASSO regression and BLUE
        cor.ndvi.las2016=c(cor.ndvi.las2016,cor.ndvi.l2016) #writing correlation out to cor.ndvi.las2016
        cor.ndvi.r2016=cor(ndvi_data[ndvi_data$trial==trial,12],pred.ndvi.ridge2016) # Correlaiton between predindvied value from Ridge regression and BLUE
        cor.ndvi.rdg2016=c(cor.ndvi.rdg2016,cor.ndvi.r2016) #writing correlation out to cor.ndvi.rdg2016
        cor.ndvi.e2016=cor(ndvi_data[ndvi_data$trial==trial,12],pred.ndvi.elasticnet2016) # Correlaiton between predindvied value from ElasticNet regression and BLUE
        cor.ndvi.enet2016=c(cor.ndvi.enet2016,cor.ndvi.e2016) #writing correlation out to cor.ndvi.enet2016
        #extrandviing predindvied values
        predindvi.ndvi.las2016 = c(predindvi.ndvi.las2016,pred.ndvi.lasso2016) # Combinig all predindvied values from LASSO
        predindvi.ndvi.rdg2016 = c(predindvi.ndvi.rdg2016,pred.ndvi.ridge2016) # Combinig all predindvied values from Ridge
        predindvi.ndvi.enet2016 = c(predindvi.ndvi.enet2016,pred.ndvi.elasticnet2016) # Combinig all predindvied values from ElasticNet
} #end for loop 

#print predindviion accuracies
cor.ndvi.step2016 #individual trial correlations for step regression
predaccuracy.ndvi.step2016=mean(cor.ndvi.step2016) #average correlation from stepwise model
cor.ndvi.las2016 #individual trial correlations for LASSO regression
predaccuracy.ndvi.lasso2016=mean(cor.ndvi.las2016) #average correlation from LASSO model
cor.ndvi.rdg2016 #individual trial correlations for Ridge regression
predaccuracy.ndvi.ridge2016=mean(cor.ndvi.rdg2016) #average correlation from Ridge model
cor.ndvi.enet2016 #individual trial correlations for ElasticNet regression
predaccuracy.ndvi.enet2016=mean(cor.ndvi.enet2016) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ndvi
plot(step.predyield.ndvi.2016$GRYLD, step.predyield.ndvi.2016$step.predvalue.ndvi.2016, xlab = "Observed Yield",ylab = "Predindvied Yield", main="Stepwise regression with ndvi 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.step2016 <-lm(step.predyield.ndvi.2016$step.predvalue.ndvi.2016 ~ step.predyield.ndvi.2016$GRYLD) #fit the model for observed values
abline(fit.ndvi.step2016, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.step2016 = vector('expression',2)
rp.ndvi.step2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.ndvi.step2016)$adj.r.squared,dig=2)))[2] # extrandviing R2 value
rp.ndvi.step2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.ndvi.step2016, digits = 2)))[2] # r (predindviion accuracy) value
legend("bottomright", bty="n",legend=rp.ndvi.step2016)
R2.ndvi.step2016=summary(fit.ndvi.step2016)$adj.r.squared
R2.ndvi.step2016=round(R2.ndvi.step2016,digits=2)
r.ndvi.step2016=round(predaccuracy.ndvi.step2016,digits=2)
models="stepwise"
traits="ndvi"
pred.accuracy.ndvi.step2016=data.frame(cbind(traits,models,r.ndvi.step2016,R2.ndvi.step2016))
colnames(pred.accuracy.ndvi.step2016)[3:4]=c("r","R2")

### Plotting LASSO regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.las2016, xlab = "Observed Yield",ylab = "Predindvied Yield", main="LASSO regression with ndvi 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.las2016 <-lm(predindvi.ndvi.las2016 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.las2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.las2016 = vector('expression',2)
rp.ndvi.las2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(summary(fit.ndvi.las2016)$adj.r.squared,dig=2)))[2]
rp.ndvi.las2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(predaccuracy.ndvi.lasso2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.las2016)
R2.ndvi.las2016=summary(fit.ndvi.las2016)$adj.r.squared
R2.ndvi.las2016=round(R2.ndvi.las2016,digits=2)
r.ndvi.las2016=round(predaccuracy.ndvi.lasso2016,digits=2)
models="LASSO"
traits="ndvi"
pred.accuracy.ndvi.las2016=data.frame(cbind(traits,models,r.ndvi.las2016,R2.ndvi.las2016))
colnames(pred.accuracy.ndvi.las2016)[3:4]=c("r","R2")

### Plotting ridge regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.rdg2016, xlab = "Observed Yield",ylab = "Predindvied Yield", main="Ridge regression with ndvi 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.rdg2016 <-lm(predindvi.ndvi.rdg2016 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.rdg2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.rdg2016 = vector('expression',2)
rp.ndvi.rdg2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(summary(fit.ndvi.rdg2016)$adj.r.squared,dig=2)))[2]
rp.ndvi.rdg2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(predaccuracy.ndvi.ridge2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.rdg2016)
R2.ndvi.rdg2016=summary(fit.ndvi.rdg2016)$adj.r.squared
R2.ndvi.rdg2016=round(R2.ndvi.rdg2016,digits=2)
r.ndvi.rdg2016=round(predaccuracy.ndvi.ridge2016,digits=2)
models="ridge"
traits="ndvi"
pred.accuracy.ndvi.rdg2016=data.frame(cbind(traits,models,r.ndvi.rdg2016,R2.ndvi.rdg2016))
colnames(pred.accuracy.ndvi.rdg2016)[3:4]=c("r","R2")

### Plotting enet regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.enet2016, xlab = "Observed Yield",ylab = "Predindvied Yield", main="ElasticNet regression with ndvi 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.enet2016 <-lm(predindvi.ndvi.enet2016 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.enet2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.enet2016 = vector('expression',2)
rp.ndvi.enet2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.ndvi.enet2016)$adj.r.squared,dig=2)))[2]
rp.ndvi.enet2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.ndvi.enet2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.enet2016)
R2.ndvi.enet2016=summary(fit.ndvi.enet2016)$adj.r.squared
R2.ndvi.enet2016=round(R2.ndvi.enet2016,digits=2)
r.ndvi.enet2016=round(predaccuracy.ndvi.enet2016,digits=2)
models="elasticnet"
traits="ndvi"
pred.accuracy.ndvi.enet2016=data.frame(cbind(traits,models,r.ndvi.enet2016,R2.ndvi.enet2016))
colnames(pred.accuracy.ndvi.enet2016)[3:4]=c("r","R2")

multivar.ndvi.prediction.accuracy2016=data.frame(rbind(pred.accuracy.ndvi.step2016,pred.accuracy.ndvi.las2016,pred.accuracy.ndvi.rdg2016,pred.accuracy.ndvi.enet2016))
# write.csv(multivar.ndvi.prediction.accuracy2016,file="Multivariate_predaccuracy_ndvi_2016.csv",row.names = FALSE,quote = FALSE)



### Multivariate model yield predindviion 2017 with ndvi
### Running stepwise regression model with ndvi
BLUE2017=read.csv("BLUE_2017.csv")
ndvi_data=BLUE2017[,c(1:2,17:30,39)]
stepwise.ndvi.2017 <- lm(GRYLD ~ NDVI_20170103+NDVI_20170108+NDVI_20170114+NDVI_20170120+NDVI_20170125+NDVI_20170131+NDVI_20170205+NDVI_20170210+NDVI_20170215+NDVI_20170220+NDVI_20170225+NDVI_20170302+NDVI_20170307+NDVI_20170313
                         ,data=ndvi_data)  ## writing the model with ndvi
stepwise.ndvi.2017 <- step(stepwise.ndvi.2017, direndviion = "both") ## running the stepwise regressing model
formula(stepwise.ndvi.2017) ## getting the selendvied variables from the stepwise regression model

### Predindvied yield from stepwise regression model with ndvi
step.predyield.ndvi.2017=NULL 
for(i in 1:11){
        step.prediction.ndvi.2017=ndvi_data[ndvi_data$trial==i,] #get just one trial to predindvi
        step.training.ndvi.2017=ndvi_data[ndvi_data$trial!=i,] #get all other trials to train on
        fit.ndvi.2017 = lm(formula = formula(stepwise.ndvi.2017), data =step.training.ndvi.2017) ## direndvily import variables into the model from stepwise regression 
        step.predvalue.ndvi.2017=predict(object=fit.ndvi.2017, newdata=step.prediction.ndvi.2017) ## Predindvied value of individual trial
        predicted.ndvi.trial2017=data.frame(entry=step.prediction.ndvi.2017$plot,step.predvalue.ndvi.2017) ## Predindvied value of individual trial in data frame
        step.predyield.ndvi.2017=rbind(step.predyield.ndvi.2017, predicted.ndvi.trial2017) ## Combinig all trial's predindvied value into one frame
}
step.predyield.ndvi.2017 = data.frame(trial=rep(1:11,each=60),step.predyield.ndvi.2017)
step.predyield.ndvi.2017 <- merge(step.predyield.ndvi.2017, ndvi_data[, c(1:2,17)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extrandviing predindviion accuracy from stepwise regression model with ndvi
cor.ndvi.step2017 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ndvi.2017[step.predyield.ndvi.2017$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ndvi.2017, trial_data$GRYLD) #get correlations
        cor.ndvi.step2017 <- c(cor.ndvi.step2017, trial_cor) #write results outside of the for loop
}

### Predindviion from shrinkage models
predindvi.ndvi.las2017=c()  #Intialize to capture predindvied values from LASSO model
predindvi.ndvi.rdg2017=c() #Intialize to capture predindvied values from Ridge model
predindvi.ndvi.enet2017=c() #Intialize to capture predindvied values from ElasticNet model
cor.ndvi.las2017=c() #Intialize to capture correlatio from LASSO model
cor.ndvi.rdg2017=c() #Intialize to capture correlatio from Ridge model
cor.ndvi.enet2017=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ndvi.data2017 = ndvi_data[ndvi_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ndvi.data2017 = ndvi_data[ndvi_data$trial!=trial,] #sets up predindviion populations
        rndvirl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ndvi.model2017 = train(GRYLD ~ ., data = train.ndvi.data2017[,-c(1,2)], #Looks like we are fitting a lasso model.
                                     method = "lasso",
                                     trControl = rndvirl,
                                     preProc = c("center", "scale"))
        ridge.ndvi.model2017 = train(GRYLD ~ ., data = train.ndvi.data2017[,-c(1,2)], #looks like we are fitting a ridge model.
                                     method = "ridge",
                                     trControl = rndvirl,
                                     preProc = c("center", "scale"))
        enet.ndvi.model2017 = train(GRYLD ~ ., data = train.ndvi.data2017[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                    method = "enet",
                                    trControl = rndvirl,
                                    preProc = c("center", "scale"))
        #model coefficients
        coef.ndvi.lasso2017 = lasso.ndvi.model2017$finalModel #extrandviing some coefficients from LASSO model
        coef.ndvi.ridge2017 = ridge.ndvi.model2017$finalModel #extrandviing some coefficients from Ridge model
        coef.ndvi.enet2017 = enet.ndvi.model2017$finalModel #extrandviing some coefficients from ElasticNet model
        #predindviion
        pred.ndvi.lasso2017 = predict(lasso.ndvi.model2017,test.ndvi.data2017[,-c(1,2)]) #making a predindvi statement.  This looks promising
        pred.ndvi.ridge2017 = predict(ridge.ndvi.model2017,test.ndvi.data2017[,-c(1,2)]) #appears the same with other models
        pred.ndvi.elasticnet2017 = predict(enet.ndvi.model2017,test.ndvi.data2017[,-c(1,2)])
        #extrandviing correlations
        cor.ndvi.l2017=cor(ndvi_data[ndvi_data$trial==trial,17],pred.ndvi.lasso2017) # Correlaiton between predindvied value from LASSO regression and BLUE
        cor.ndvi.las2017=c(cor.ndvi.las2017,cor.ndvi.l2017) #writing correlation out to cor.ndvi.las2017
        cor.ndvi.r2017=cor(ndvi_data[ndvi_data$trial==trial,17],pred.ndvi.ridge2017) # Correlaiton between predindvied value from Ridge regression and BLUE
        cor.ndvi.rdg2017=c(cor.ndvi.rdg2017,cor.ndvi.r2017) #writing correlation out to cor.ndvi.rdg2017
        cor.ndvi.e2017=cor(ndvi_data[ndvi_data$trial==trial,17],pred.ndvi.elasticnet2017) # Correlaiton between predindvied value from ElasticNet regression and BLUE
        cor.ndvi.enet2017=c(cor.ndvi.enet2017,cor.ndvi.e2017) #writing correlation out to cor.ndvi.enet2017
        #extrandviing predindvied values
        predindvi.ndvi.las2017 = c(predindvi.ndvi.las2017,pred.ndvi.lasso2017) # Combinig all predindvied values from LASSO
        predindvi.ndvi.rdg2017 = c(predindvi.ndvi.rdg2017,pred.ndvi.ridge2017) # Combinig all predindvied values from Ridge
        predindvi.ndvi.enet2017 = c(predindvi.ndvi.enet2017,pred.ndvi.elasticnet2017) # Combinig all predindvied values from ElasticNet
} #end for loop 

#print predindviion accuracies
cor.ndvi.step2017 #individual trial correlations for step regression
predaccuracy.ndvi.step2017=mean(cor.ndvi.step2017) #average correlation from stepwise model
cor.ndvi.las2017 #individual trial correlations for LASSO regression
predaccuracy.ndvi.lasso2017=mean(cor.ndvi.las2017) #average correlation from LASSO model
cor.ndvi.rdg2017 #individual trial correlations for Ridge regression
predaccuracy.ndvi.ridge2017=mean(cor.ndvi.rdg2017) #average correlation from Ridge model
cor.ndvi.enet2017 #individual trial correlations for ElasticNet regression
predaccuracy.ndvi.enet2017=mean(cor.ndvi.enet2017) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ndvi
plot(step.predyield.ndvi.2017$GRYLD, step.predyield.ndvi.2017$step.predvalue.ndvi.2017, xlab = "Observed Yield",ylab = "Predindvied Yield", main="Stepwise regression with ndvi 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.step2017 <-lm(step.predyield.ndvi.2017$step.predvalue.ndvi.2017 ~ step.predyield.ndvi.2017$GRYLD) #fit the model for observed values
abline(fit.ndvi.step2017, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.step2017 = vector('expression',2)
rp.ndvi.step2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.ndvi.step2017)$adj.r.squared,dig=2)))[2] # extrandviing R2 value
rp.ndvi.step2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.ndvi.step2017, digits = 2)))[2] # r (predindviion accuracy) value
legend("bottomright", bty="n",legend=rp.ndvi.step2017)
R2.ndvi.step2017=summary(fit.ndvi.step2017)$adj.r.squared
R2.ndvi.step2017=round(R2.ndvi.step2017,digits=2)
r.ndvi.step2017=round(predaccuracy.ndvi.step2017,digits=2)
models="stepwise"
traits="ndvi"
pred.accuracy.ndvi.step2017=data.frame(cbind(traits,models,r.ndvi.step2017,R2.ndvi.step2017))
colnames(pred.accuracy.ndvi.step2017)[3:4]=c("r","R2")

### Plotting LASSO regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.las2017, xlab = "Observed Yield",ylab = "Predindvied Yield", main="LASSO regression with ndvi 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.las2017 <-lm(predindvi.ndvi.las2017 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.las2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.las2017 = vector('expression',2)
rp.ndvi.las2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(summary(fit.ndvi.las2017)$adj.r.squared,dig=2)))[2]
rp.ndvi.las2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(predaccuracy.ndvi.lasso2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.las2017)
R2.ndvi.las2017=summary(fit.ndvi.las2017)$adj.r.squared
R2.ndvi.las2017=round(R2.ndvi.las2017,digits=2)
r.ndvi.las2017=round(predaccuracy.ndvi.lasso2017,digits=2)
models="LASSO"
traits="ndvi"
pred.accuracy.ndvi.las2017=data.frame(cbind(traits,models,r.ndvi.las2017,R2.ndvi.las2017))
colnames(pred.accuracy.ndvi.las2017)[3:4]=c("r","R2")

### Plotting ridge regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.rdg2017, xlab = "Observed Yield",ylab = "Predindvied Yield", main="Ridge regression with ndvi 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.rdg2017 <-lm(predindvi.ndvi.rdg2017 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.rdg2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.rdg2017 = vector('expression',2)
rp.ndvi.rdg2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(summary(fit.ndvi.rdg2017)$adj.r.squared,dig=2)))[2]
rp.ndvi.rdg2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(predaccuracy.ndvi.ridge2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.rdg2017)
R2.ndvi.rdg2017=summary(fit.ndvi.rdg2017)$adj.r.squared
R2.ndvi.rdg2017=round(R2.ndvi.rdg2017,digits=2)
r.ndvi.rdg2017=round(predaccuracy.ndvi.ridge2017,digits=2)
models="ridge"
traits="ndvi"
pred.accuracy.ndvi.rdg2017=data.frame(cbind(traits,models,r.ndvi.rdg2017,R2.ndvi.rdg2017))
colnames(pred.accuracy.ndvi.rdg2017)[3:4]=c("r","R2")

### Plotting enet regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.enet2017, xlab = "Observed Yield",ylab = "Predindvied Yield", main="ElasticNet regression with ndvi 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.enet2017 <-lm(predindvi.ndvi.enet2017 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.enet2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.enet2017 = vector('expression',2)
rp.ndvi.enet2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.ndvi.enet2017)$adj.r.squared,dig=2)))[2]
rp.ndvi.enet2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.ndvi.enet2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.enet2017)
R2.ndvi.enet2017=summary(fit.ndvi.enet2017)$adj.r.squared
R2.ndvi.enet2017=round(R2.ndvi.enet2017,digits=2)
r.ndvi.enet2017=round(predaccuracy.ndvi.enet2017,digits=2)
models="elasticnet"
traits="ndvi"
pred.accuracy.ndvi.enet2017=data.frame(cbind(traits,models,r.ndvi.enet2017,R2.ndvi.enet2017))
colnames(pred.accuracy.ndvi.enet2017)[3:4]=c("r","R2")

multivar.ndvi.prediction.accuracy2017=data.frame(rbind(pred.accuracy.ndvi.step2017,pred.accuracy.ndvi.las2017,pred.accuracy.ndvi.rdg2017,pred.accuracy.ndvi.enet2017))
# write.csv(multivar.ndvi.prediction.accuracy2017,file="Multivariate_predaccuracy_ndvi_2017.csv",row.names = FALSE,quote = FALSE)



### Multivariate model yield predindviion 2018 with ndvi
### Running stepwise regression model with ndvi
BLUE2018=read.csv("BLUE_2018.csv")
ndvi_data=BLUE2018[,c(1:2,15:26,35)]
stepwise.ndvi.2018 <- lm(GRYLD ~ NDVI_20180126+NDVI_20180131+NDVI_20180204+NDVI_20180210+NDVI_20180214+NDVI_20180219+NDVI_20180226+NDVI_20180301+NDVI_20180305+NDVI_20180310+NDVI_20180315+NDVI_20180320
                         ,data=ndvi_data)  ## writing the model with ndvi
stepwise.ndvi.2018 <- step(stepwise.ndvi.2018, direndviion = "both") ## running the stepwise regressing model
formula(stepwise.ndvi.2018) ## getting the selendvied variables from the stepwise regression model

### Predindvied yield from stepwise regression model with ndvi
step.predyield.ndvi.2018=NULL 
for(i in 1:11){
        step.prediction.ndvi.2018=ndvi_data[ndvi_data$trial==i,] #get just one trial to predindvi
        step.training.ndvi.2018=ndvi_data[ndvi_data$trial!=i,] #get all other trials to train on
        fit.ndvi.2018 = lm(formula = formula(stepwise.ndvi.2018), data =step.training.ndvi.2018) ## direndvily import variables into the model from stepwise regression 
        step.predvalue.ndvi.2018=predict(object=fit.ndvi.2018, newdata=step.prediction.ndvi.2018) ## Predindvied value of individual trial
        predicted.ndvi.trial2018=data.frame(entry=step.prediction.ndvi.2018$plot,step.predvalue.ndvi.2018) ## Predindvied value of individual trial in data frame
        step.predyield.ndvi.2018=rbind(step.predyield.ndvi.2018, predicted.ndvi.trial2018) ## Combinig all trial's predindvied value into one frame
}
step.predyield.ndvi.2018 = data.frame(trial=rep(1:11,each=60),step.predyield.ndvi.2018)
step.predyield.ndvi.2018 <- merge(step.predyield.ndvi.2018, ndvi_data[, c(1:2,15)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extrandviing predindviion accuracy from stepwise regression model with ndvi
cor.ndvi.step2018 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ndvi.2018[step.predyield.ndvi.2018$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ndvi.2018, trial_data$GRYLD) #get correlations
        cor.ndvi.step2018 <- c(cor.ndvi.step2018, trial_cor) #write results outside of the for loop
}

### Predindviion from shrinkage models
predindvi.ndvi.las2018=c()  #Intialize to capture predindvied values from LASSO model
predindvi.ndvi.rdg2018=c() #Intialize to capture predindvied values from Ridge model
predindvi.ndvi.enet2018=c() #Intialize to capture predindvied values from ElasticNet model
cor.ndvi.las2018=c() #Intialize to capture correlatio from LASSO model
cor.ndvi.rdg2018=c() #Intialize to capture correlatio from Ridge model
cor.ndvi.enet2018=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ndvi.data2018 = ndvi_data[ndvi_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ndvi.data2018 = ndvi_data[ndvi_data$trial!=trial,] #sets up predindviion populations
        rndvirl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ndvi.model2018 = train(GRYLD ~ ., data = train.ndvi.data2018[,-c(1,2)], #Looks like we are fitting a lasso model.
                                     method = "lasso",
                                     trControl = rndvirl,
                                     preProc = c("center", "scale"))
        ridge.ndvi.model2018 = train(GRYLD ~ ., data = train.ndvi.data2018[,-c(1,2)], #looks like we are fitting a ridge model.
                                     method = "ridge",
                                     trControl = rndvirl,
                                     preProc = c("center", "scale"))
        enet.ndvi.model2018 = train(GRYLD ~ ., data = train.ndvi.data2018[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                    method = "enet",
                                    trControl = rndvirl,
                                    preProc = c("center", "scale"))
        #model coefficients
        coef.ndvi.lasso2018 = lasso.ndvi.model2018$finalModel #extrandviing some coefficients from LASSO model
        coef.ndvi.ridge2018 = ridge.ndvi.model2018$finalModel #extrandviing some coefficients from Ridge model
        coef.ndvi.enet2018 = enet.ndvi.model2018$finalModel #extrandviing some coefficients from ElasticNet model
        #predindviion
        pred.ndvi.lasso2018 = predict(lasso.ndvi.model2018,test.ndvi.data2018[,-c(1,2)]) #making a predindvi statement.  This looks promising
        pred.ndvi.ridge2018 = predict(ridge.ndvi.model2018,test.ndvi.data2018[,-c(1,2)]) #appears the same with other models
        pred.ndvi.elasticnet2018 = predict(enet.ndvi.model2018,test.ndvi.data2018[,-c(1,2)])
        #extrandviing correlations
        cor.ndvi.l2018=cor(ndvi_data[ndvi_data$trial==trial,15],pred.ndvi.lasso2018) # Correlaiton between predindvied value from LASSO regression and BLUE
        cor.ndvi.las2018=c(cor.ndvi.las2018,cor.ndvi.l2018) #writing correlation out to cor.ndvi.las2018
        cor.ndvi.r2018=cor(ndvi_data[ndvi_data$trial==trial,15],pred.ndvi.ridge2018) # Correlaiton between predindvied value from Ridge regression and BLUE
        cor.ndvi.rdg2018=c(cor.ndvi.rdg2018,cor.ndvi.r2018) #writing correlation out to cor.ndvi.rdg2018
        cor.ndvi.e2018=cor(ndvi_data[ndvi_data$trial==trial,15],pred.ndvi.elasticnet2018) # Correlaiton between predindvied value from ElasticNet regression and BLUE
        cor.ndvi.enet2018=c(cor.ndvi.enet2018,cor.ndvi.e2018) #writing correlation out to cor.ndvi.enet2018
        #extrandviing predindvied values
        predindvi.ndvi.las2018 = c(predindvi.ndvi.las2018,pred.ndvi.lasso2018) # Combinig all predindvied values from LASSO
        predindvi.ndvi.rdg2018 = c(predindvi.ndvi.rdg2018,pred.ndvi.ridge2018) # Combinig all predindvied values from Ridge
        predindvi.ndvi.enet2018 = c(predindvi.ndvi.enet2018,pred.ndvi.elasticnet2018) # Combinig all predindvied values from ElasticNet
} #end for loop 

#print predindviion accuracies
cor.ndvi.step2018 #individual trial correlations for step regression
predaccuracy.ndvi.step2018=mean(cor.ndvi.step2018) #average correlation from stepwise model
cor.ndvi.las2018 #individual trial correlations for LASSO regression
predaccuracy.ndvi.lasso2018=mean(cor.ndvi.las2018) #average correlation from LASSO model
cor.ndvi.rdg2018 #individual trial correlations for Ridge regression
predaccuracy.ndvi.ridge2018=mean(cor.ndvi.rdg2018) #average correlation from Ridge model
cor.ndvi.enet2018 #individual trial correlations for ElasticNet regression
predaccuracy.ndvi.enet2018=mean(cor.ndvi.enet2018) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ndvi
plot(step.predyield.ndvi.2018$GRYLD, step.predyield.ndvi.2018$step.predvalue.ndvi.2018, xlab = "Observed Yield",ylab = "Predindvied Yield", main="Stepwise regression with ndvi 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.step2018 <-lm(step.predyield.ndvi.2018$step.predvalue.ndvi.2018 ~ step.predyield.ndvi.2018$GRYLD) #fit the model for observed values
abline(fit.ndvi.step2018, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.step2018 = vector('expression',2)
rp.ndvi.step2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.ndvi.step2018)$adj.r.squared,dig=2)))[2] # extrandviing R2 value
rp.ndvi.step2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.ndvi.step2018, digits = 2)))[2] # r (predindviion accuracy) value
legend("bottomright", bty="n",legend=rp.ndvi.step2018)
R2.ndvi.step2018=summary(fit.ndvi.step2018)$adj.r.squared
R2.ndvi.step2018=round(R2.ndvi.step2018,digits=2)
r.ndvi.step2018=round(predaccuracy.ndvi.step2018,digits=2)
models="stepwise"
traits="ndvi"
pred.accuracy.ndvi.step2018=data.frame(cbind(traits,models,r.ndvi.step2018,R2.ndvi.step2018))
colnames(pred.accuracy.ndvi.step2018)[3:4]=c("r","R2")

### Plotting LASSO regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.las2018, xlab = "Observed Yield",ylab = "Predindvied Yield", main="LASSO regression with ndvi 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.las2018 <-lm(predindvi.ndvi.las2018 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.las2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.las2018 = vector('expression',2)
rp.ndvi.las2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(summary(fit.ndvi.las2018)$adj.r.squared,dig=2)))[2]
rp.ndvi.las2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(predaccuracy.ndvi.lasso2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.las2018)
R2.ndvi.las2018=summary(fit.ndvi.las2018)$adj.r.squared
R2.ndvi.las2018=round(R2.ndvi.las2018,digits=2)
r.ndvi.las2018=round(predaccuracy.ndvi.lasso2018,digits=2)
models="LASSO"
traits="ndvi"
pred.accuracy.ndvi.las2018=data.frame(cbind(traits,models,r.ndvi.las2018,R2.ndvi.las2018))
colnames(pred.accuracy.ndvi.las2018)[3:4]=c("r","R2")

### Plotting ridge regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.rdg2018, xlab = "Observed Yield",ylab = "Predindvied Yield", main="Ridge regression with ndvi 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.rdg2018 <-lm(predindvi.ndvi.rdg2018 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.rdg2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.rdg2018 = vector('expression',2)
rp.ndvi.rdg2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(summary(fit.ndvi.rdg2018)$adj.r.squared,dig=2)))[2]
rp.ndvi.rdg2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(predaccuracy.ndvi.ridge2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.rdg2018)
R2.ndvi.rdg2018=summary(fit.ndvi.rdg2018)$adj.r.squared
R2.ndvi.rdg2018=round(R2.ndvi.rdg2018,digits=2)
r.ndvi.rdg2018=round(predaccuracy.ndvi.ridge2018,digits=2)
models="ridge"
traits="ndvi"
pred.accuracy.ndvi.rdg2018=data.frame(cbind(traits,models,r.ndvi.rdg2018,R2.ndvi.rdg2018))
colnames(pred.accuracy.ndvi.rdg2018)[3:4]=c("r","R2")

### Plotting enet regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.enet2018, xlab = "Observed Yield",ylab = "Predindvied Yield", main="ElasticNet regression with ndvi 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.enet2018 <-lm(predindvi.ndvi.enet2018 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.enet2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.enet2018 = vector('expression',2)
rp.ndvi.enet2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.ndvi.enet2018)$adj.r.squared,dig=2)))[2]
rp.ndvi.enet2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.ndvi.enet2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.enet2018)
R2.ndvi.enet2018=summary(fit.ndvi.enet2018)$adj.r.squared
R2.ndvi.enet2018=round(R2.ndvi.enet2018,digits=2)
r.ndvi.enet2018=round(predaccuracy.ndvi.enet2018,digits=2)
models="elasticnet"
traits="ndvi"
pred.accuracy.ndvi.enet2018=data.frame(cbind(traits,models,r.ndvi.enet2018,R2.ndvi.enet2018))
colnames(pred.accuracy.ndvi.enet2018)[3:4]=c("r","R2")

multivar.ndvi.prediction.accuracy2018=data.frame(rbind(pred.accuracy.ndvi.step2018,pred.accuracy.ndvi.las2018,pred.accuracy.ndvi.rdg2018,pred.accuracy.ndvi.enet2018))
# write.csv(multivar.ndvi.prediction.accuracy2018,file="Multivariate_predaccuracy_ndvi_2018.csv",row.names = FALSE,quote = FALSE)



### Multivariate model yield predindviion 2019 with ndvi
### Running stepwise regression model with ndvi
BLUE2019=read.csv("BLUE_2019.csv")
ndvi_data=BLUE2019[,c(1:2,16:28,37)]
stepwise.ndvi.2019 <- lm(GRYLD ~ NDVI_20190121+NDVI_20190127+NDVI_20190131+NDVI_20190205+NDVI_20190211+NDVI_20190218+NDVI_20190222+NDVI_20190228+NDVI_20190305+NDVI_20190311+NDVI_20190315+NDVI_20190320+NDVI_20190325
                         ,data=ndvi_data)  ## writing the model with ndvi
stepwise.ndvi.2019 <- step(stepwise.ndvi.2019, direndviion = "both") ## running the stepwise regressing model
formula(stepwise.ndvi.2019) ## getting the selendvied variables from the stepwise regression model

### Predindvied yield from stepwise regression model with ndvi
step.predyield.ndvi.2019=NULL 
for(i in 1:10){
        step.prediction.ndvi.2019=ndvi_data[ndvi_data$trial==i,] #get just one trial to predindvi
        step.training.ndvi.2019=ndvi_data[ndvi_data$trial!=i,] #get all other trials to train on
        fit.ndvi.2019 = lm(formula = formula(stepwise.ndvi.2019), data =step.training.ndvi.2019) ## direndvily import variables into the model from stepwise regression 
        step.predvalue.ndvi.2019=predict(object=fit.ndvi.2019, newdata=step.prediction.ndvi.2019) ## Predindvied value of individual trial
        predicted.ndvi.trial2019=data.frame(entry=step.prediction.ndvi.2019$plot,step.predvalue.ndvi.2019) ## Predindvied value of individual trial in data frame
        step.predyield.ndvi.2019=rbind(step.predyield.ndvi.2019, predicted.ndvi.trial2019) ## Combinig all trial's predindvied value into one frame
}
step.predyield.ndvi.2019 = data.frame(trial=rep(1:10,each=60),step.predyield.ndvi.2019)
step.predyield.ndvi.2019 <- merge(step.predyield.ndvi.2019, ndvi_data[, c(1:2,16)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extrandviing predindviion accuracy from stepwise regression model with ndvi
cor.ndvi.step2019 <- NULL #start dataframe to hold correlation values
for(i in 1:10){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ndvi.2019[step.predyield.ndvi.2019$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ndvi.2019, trial_data$GRYLD) #get correlations
        cor.ndvi.step2019 <- c(cor.ndvi.step2019, trial_cor) #write results outside of the for loop
}

### Predindviion from shrinkage models
predindvi.ndvi.las2019=c()  #Intialize to capture predindvied values from LASSO model
predindvi.ndvi.rdg2019=c() #Intialize to capture predindvied values from Ridge model
predindvi.ndvi.enet2019=c() #Intialize to capture predindvied values from ElasticNet model
cor.ndvi.las2019=c() #Intialize to capture correlatio from LASSO model
cor.ndvi.rdg2019=c() #Intialize to capture correlatio from Ridge model
cor.ndvi.enet2019=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:10){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ndvi.data2019 = ndvi_data[ndvi_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ndvi.data2019 = ndvi_data[ndvi_data$trial!=trial,] #sets up predindviion populations
        rndvirl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ndvi.model2019 = train(GRYLD ~ ., data = train.ndvi.data2019[,-c(1,2)], #Looks like we are fitting a lasso model.
                                     method = "lasso",
                                     trControl = rndvirl,
                                     preProc = c("center", "scale"))
        ridge.ndvi.model2019 = train(GRYLD ~ ., data = train.ndvi.data2019[,-c(1,2)], #looks like we are fitting a ridge model.
                                     method = "ridge",
                                     trControl = rndvirl,
                                     preProc = c("center", "scale"))
        enet.ndvi.model2019 = train(GRYLD ~ ., data = train.ndvi.data2019[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                    method = "enet",
                                    trControl = rndvirl,
                                    preProc = c("center", "scale"))
        #model coefficients
        coef.ndvi.lasso2019 = lasso.ndvi.model2019$finalModel #extrandviing some coefficients from LASSO model
        coef.ndvi.ridge2019 = ridge.ndvi.model2019$finalModel #extrandviing some coefficients from Ridge model
        coef.ndvi.enet2019 = enet.ndvi.model2019$finalModel #extrandviing some coefficients from ElasticNet model
        #predindviion
        pred.ndvi.lasso2019 = predict(lasso.ndvi.model2019,test.ndvi.data2019[,-c(1,2)]) #making a predindvi statement.  This looks promising
        pred.ndvi.ridge2019 = predict(ridge.ndvi.model2019,test.ndvi.data2019[,-c(1,2)]) #appears the same with other models
        pred.ndvi.elasticnet2019 = predict(enet.ndvi.model2019,test.ndvi.data2019[,-c(1,2)])
        #extrandviing correlations
        cor.ndvi.l2019=cor(ndvi_data[ndvi_data$trial==trial,16],pred.ndvi.lasso2019) # Correlaiton between predindvied value from LASSO regression and BLUE
        cor.ndvi.las2019=c(cor.ndvi.las2019,cor.ndvi.l2019) #writing correlation out to cor.ndvi.las2019
        cor.ndvi.r2019=cor(ndvi_data[ndvi_data$trial==trial,16],pred.ndvi.ridge2019) # Correlaiton between predindvied value from Ridge regression and BLUE
        cor.ndvi.rdg2019=c(cor.ndvi.rdg2019,cor.ndvi.r2019) #writing correlation out to cor.ndvi.rdg2019
        cor.ndvi.e2019=cor(ndvi_data[ndvi_data$trial==trial,16],pred.ndvi.elasticnet2019) # Correlaiton between predindvied value from ElasticNet regression and BLUE
        cor.ndvi.enet2019=c(cor.ndvi.enet2019,cor.ndvi.e2019) #writing correlation out to cor.ndvi.enet2019
        #extrandviing predindvied values
        predindvi.ndvi.las2019 = c(predindvi.ndvi.las2019,pred.ndvi.lasso2019) # Combinig all predindvied values from LASSO
        predindvi.ndvi.rdg2019 = c(predindvi.ndvi.rdg2019,pred.ndvi.ridge2019) # Combinig all predindvied values from Ridge
        predindvi.ndvi.enet2019 = c(predindvi.ndvi.enet2019,pred.ndvi.elasticnet2019) # Combinig all predindvied values from ElasticNet
} #end for loop 

#print predindviion accuracies
cor.ndvi.step2019 #individual trial correlations for step regression
predaccuracy.ndvi.step2019=mean(cor.ndvi.step2019) #average correlation from stepwise model
cor.ndvi.las2019 #individual trial correlations for LASSO regression
predaccuracy.ndvi.lasso2019=mean(cor.ndvi.las2019) #average correlation from LASSO model
cor.ndvi.rdg2019 #individual trial correlations for Ridge regression
predaccuracy.ndvi.ridge2019=mean(cor.ndvi.rdg2019) #average correlation from Ridge model
cor.ndvi.enet2019 #individual trial correlations for ElasticNet regression
predaccuracy.ndvi.enet2019=mean(cor.ndvi.enet2019) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ndvi
plot(step.predyield.ndvi.2019$GRYLD, step.predyield.ndvi.2019$step.predvalue.ndvi.2019, xlab = "Observed Yield",ylab = "Predindvied Yield", main="Stepwise regression with ndvi 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.step2019 <-lm(step.predyield.ndvi.2019$step.predvalue.ndvi.2019 ~ step.predyield.ndvi.2019$GRYLD) #fit the model for observed values
abline(fit.ndvi.step2019, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.step2019 = vector('expression',2)
rp.ndvi.step2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.ndvi.step2019)$adj.r.squared,dig=2)))[2] # extrandviing R2 value
rp.ndvi.step2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.ndvi.step2019, digits = 2)))[2] # r (predindviion accuracy) value
legend("bottomright", bty="n",legend=rp.ndvi.step2019)
R2.ndvi.step2019=summary(fit.ndvi.step2019)$adj.r.squared
R2.ndvi.step2019=round(R2.ndvi.step2019,digits=2)
r.ndvi.step2019=round(predaccuracy.ndvi.step2019,digits=2)
models="stepwise"
traits="ndvi"
pred.accuracy.ndvi.step2019=data.frame(cbind(traits,models,r.ndvi.step2019,R2.ndvi.step2019))
colnames(pred.accuracy.ndvi.step2019)[3:4]=c("r","R2")

### Plotting LASSO regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.las2019, xlab = "Observed Yield",ylab = "Predindvied Yield", main="LASSO regression with ndvi 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.las2019 <-lm(predindvi.ndvi.las2019 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.las2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.las2019 = vector('expression',2)
rp.ndvi.las2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(summary(fit.ndvi.las2019)$adj.r.squared,dig=2)))[2]
rp.ndvi.las2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(predaccuracy.ndvi.lasso2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.las2019)
R2.ndvi.las2019=summary(fit.ndvi.las2019)$adj.r.squared
R2.ndvi.las2019=round(R2.ndvi.las2019,digits=2)
r.ndvi.las2019=round(predaccuracy.ndvi.lasso2019,digits=2)
models="LASSO"
traits="ndvi"
pred.accuracy.ndvi.las2019=data.frame(cbind(traits,models,r.ndvi.las2019,R2.ndvi.las2019))
colnames(pred.accuracy.ndvi.las2019)[3:4]=c("r","R2")

### Plotting ridge regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.rdg2019, xlab = "Observed Yield",ylab = "Predindvied Yield", main="Ridge regression with ndvi 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.rdg2019 <-lm(predindvi.ndvi.rdg2019 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.rdg2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.rdg2019 = vector('expression',2)
rp.ndvi.rdg2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(summary(fit.ndvi.rdg2019)$adj.r.squared,dig=2)))[2]
rp.ndvi.rdg2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(predaccuracy.ndvi.ridge2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.rdg2019)
R2.ndvi.rdg2019=summary(fit.ndvi.rdg2019)$adj.r.squared
R2.ndvi.rdg2019=round(R2.ndvi.rdg2019,digits=2)
r.ndvi.rdg2019=round(predaccuracy.ndvi.ridge2019,digits=2)
models="ridge"
traits="ndvi"
pred.accuracy.ndvi.rdg2019=data.frame(cbind(traits,models,r.ndvi.rdg2019,R2.ndvi.rdg2019))
colnames(pred.accuracy.ndvi.rdg2019)[3:4]=c("r","R2")

### Plotting enet regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.enet2019, xlab = "Observed Yield",ylab = "Predindvied Yield", main="ElasticNet regression with ndvi 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.enet2019 <-lm(predindvi.ndvi.enet2019 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.enet2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.enet2019 = vector('expression',2)
rp.ndvi.enet2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.ndvi.enet2019)$adj.r.squared,dig=2)))[2]
rp.ndvi.enet2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.ndvi.enet2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.enet2019)
R2.ndvi.enet2019=summary(fit.ndvi.enet2019)$adj.r.squared
R2.ndvi.enet2019=round(R2.ndvi.enet2019,digits=2)
r.ndvi.enet2019=round(predaccuracy.ndvi.enet2019,digits=2)
models="elasticnet"
traits="ndvi"
pred.accuracy.ndvi.enet2019=data.frame(cbind(traits,models,r.ndvi.enet2019,R2.ndvi.enet2019))
colnames(pred.accuracy.ndvi.enet2019)[3:4]=c("r","R2")

multivar.ndvi.prediction.accuracy2019=data.frame(rbind(pred.accuracy.ndvi.step2019,pred.accuracy.ndvi.las2019,pred.accuracy.ndvi.rdg2019,pred.accuracy.ndvi.enet2019))
# write.csv(multivar.ndvi.prediction.accuracy2019,file="Multivariate_predaccuracy_ndvi_2019.csv",row.names = FALSE,quote = FALSE)



### Multivariate model yield predindviion 2020 with ndvi
### Running stepwise regression model with ndvi
BLUE2020=read.csv("BLUE_2020.csv")
ndvi_data=BLUE2020[,c(1:2,18:32,44)]
stepwise.ndvi.2020 <- lm(GRYLD ~ NDVI_20200112+NDVI_20200116+NDVI_20200121+NDVI_20200126+NDVI_20200130+NDVI_20200205+NDVI_20200210+NDVI_20200215+NDVI_20200220+NDVI_20200226+NDVI_20200302+NDVI_20200308+NDVI_20200313+NDVI_20200318+NDVI_20200323
                         ,data=ndvi_data)  ## writing the model with ndvi
stepwise.ndvi.2020 <- step(stepwise.ndvi.2020, direndviion = "both") ## running the stepwise regressing model
formula(stepwise.ndvi.2020) ## getting the selendvied variables from the stepwise regression model

### Predindvied yield from stepwise regression model with ndvi
step.predyield.ndvi.2020=NULL 
for(i in 1:11){
        step.prediction.ndvi.2020=ndvi_data[ndvi_data$trial==i,] #get just one trial to predindvi
        step.training.ndvi.2020=ndvi_data[ndvi_data$trial!=i,] #get all other trials to train on
        fit.ndvi.2020 = lm(formula = formula(stepwise.ndvi.2020), data =step.training.ndvi.2020) ## direndvily import variables into the model from stepwise regression 
        step.predvalue.ndvi.2020=predict(object=fit.ndvi.2020, newdata=step.prediction.ndvi.2020) ## Predindvied value of individual trial
        predicted.ndvi.trial2020=data.frame(entry=step.prediction.ndvi.2020$plot,step.predvalue.ndvi.2020) ## Predindvied value of individual trial in data frame
        step.predyield.ndvi.2020=rbind(step.predyield.ndvi.2020, predicted.ndvi.trial2020) ## Combinig all trial's predindvied value into one frame
}
step.predyield.ndvi.2020 = data.frame(trial=rep(1:11,each=60),step.predyield.ndvi.2020)
step.predyield.ndvi.2020 <- merge(step.predyield.ndvi.2020, ndvi_data[, c(1:2,18)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extrandviing predindviion accuracy from stepwise regression model with ndvi
cor.ndvi.step2020 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ndvi.2020[step.predyield.ndvi.2020$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ndvi.2020, trial_data$GRYLD) #get correlations
        cor.ndvi.step2020 <- c(cor.ndvi.step2020, trial_cor) #write results outside of the for loop
}

### Predindviion from shrinkage models
predindvi.ndvi.las2020=c()  #Intialize to capture predindvied values from LASSO model
predindvi.ndvi.rdg2020=c() #Intialize to capture predindvied values from Ridge model
predindvi.ndvi.enet2020=c() #Intialize to capture predindvied values from ElasticNet model
cor.ndvi.las2020=c() #Intialize to capture correlatio from LASSO model
cor.ndvi.rdg2020=c() #Intialize to capture correlatio from Ridge model
cor.ndvi.enet2020=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ndvi.data2020 = ndvi_data[ndvi_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ndvi.data2020 = ndvi_data[ndvi_data$trial!=trial,] #sets up predindviion populations
        rndvirl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ndvi.model2020 = train(GRYLD ~ ., data = train.ndvi.data2020[,-c(1,2)], #Looks like we are fitting a lasso model.
                                     method = "lasso",
                                     trControl = rndvirl,
                                     preProc = c("center", "scale"))
        ridge.ndvi.model2020 = train(GRYLD ~ ., data = train.ndvi.data2020[,-c(1,2)], #looks like we are fitting a ridge model.
                                     method = "ridge",
                                     trControl = rndvirl,
                                     preProc = c("center", "scale"))
        enet.ndvi.model2020 = train(GRYLD ~ ., data = train.ndvi.data2020[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                    method = "enet",
                                    trControl = rndvirl,
                                    preProc = c("center", "scale"))
        #model coefficients
        coef.ndvi.lasso2020 = lasso.ndvi.model2020$finalModel #extrandviing some coefficients from LASSO model
        coef.ndvi.ridge2020 = ridge.ndvi.model2020$finalModel #extrandviing some coefficients from Ridge model
        coef.ndvi.enet2020 = enet.ndvi.model2020$finalModel #extrandviing some coefficients from ElasticNet model
        #predindviion
        pred.ndvi.lasso2020 = predict(lasso.ndvi.model2020,test.ndvi.data2020[,-c(1,2)]) #making a predindvi statement.  This looks promising
        pred.ndvi.ridge2020 = predict(ridge.ndvi.model2020,test.ndvi.data2020[,-c(1,2)]) #appears the same with other models
        pred.ndvi.elasticnet2020 = predict(enet.ndvi.model2020,test.ndvi.data2020[,-c(1,2)])
        #extrandviing correlations
        cor.ndvi.l2020=cor(ndvi_data[ndvi_data$trial==trial,18],pred.ndvi.lasso2020) # Correlaiton between predindvied value from LASSO regression and BLUE
        cor.ndvi.las2020=c(cor.ndvi.las2020,cor.ndvi.l2020) #writing correlation out to cor.ndvi.las2020
        cor.ndvi.r2020=cor(ndvi_data[ndvi_data$trial==trial,18],pred.ndvi.ridge2020) # Correlaiton between predindvied value from Ridge regression and BLUE
        cor.ndvi.rdg2020=c(cor.ndvi.rdg2020,cor.ndvi.r2020) #writing correlation out to cor.ndvi.rdg2020
        cor.ndvi.e2020=cor(ndvi_data[ndvi_data$trial==trial,18],pred.ndvi.elasticnet2020) # Correlaiton between predindvied value from ElasticNet regression and BLUE
        cor.ndvi.enet2020=c(cor.ndvi.enet2020,cor.ndvi.e2020) #writing correlation out to cor.ndvi.enet2020
        #extrandviing predindvied values
        predindvi.ndvi.las2020 = c(predindvi.ndvi.las2020,pred.ndvi.lasso2020) # Combinig all predindvied values from LASSO
        predindvi.ndvi.rdg2020 = c(predindvi.ndvi.rdg2020,pred.ndvi.ridge2020) # Combinig all predindvied values from Ridge
        predindvi.ndvi.enet2020 = c(predindvi.ndvi.enet2020,pred.ndvi.elasticnet2020) # Combinig all predindvied values from ElasticNet
} #end for loop 

#print predindviion accuracies
cor.ndvi.step2020 #individual trial correlations for step regression
predaccuracy.ndvi.step2020=mean(cor.ndvi.step2020) #average correlation from stepwise model
cor.ndvi.las2020 #individual trial correlations for LASSO regression
predaccuracy.ndvi.lasso2020=mean(cor.ndvi.las2020) #average correlation from LASSO model
cor.ndvi.rdg2020 #individual trial correlations for Ridge regression
predaccuracy.ndvi.ridge2020=mean(cor.ndvi.rdg2020) #average correlation from Ridge model
cor.ndvi.enet2020 #individual trial correlations for ElasticNet regression
predaccuracy.ndvi.enet2020=mean(cor.ndvi.enet2020) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ndvi
plot(step.predyield.ndvi.2020$GRYLD, step.predyield.ndvi.2020$step.predvalue.ndvi.2020, xlab = "Observed Yield",ylab = "Predindvied Yield", main="Stepwise regression with ndvi 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.step2020 <-lm(step.predyield.ndvi.2020$step.predvalue.ndvi.2020 ~ step.predyield.ndvi.2020$GRYLD) #fit the model for observed values
abline(fit.ndvi.step2020, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.step2020 = vector('expression',2)
rp.ndvi.step2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.ndvi.step2020)$adj.r.squared,dig=2)))[2] # extrandviing R2 value
rp.ndvi.step2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.ndvi.step2020, digits = 2)))[2] # r (predindviion accuracy) value
legend("bottomright", bty="n",legend=rp.ndvi.step2020)
R2.ndvi.step2020=summary(fit.ndvi.step2020)$adj.r.squared
R2.ndvi.step2020=round(R2.ndvi.step2020,digits=2)
r.ndvi.step2020=round(predaccuracy.ndvi.step2020,digits=2)
models="stepwise"
traits="ndvi"
pred.accuracy.ndvi.step2020=data.frame(cbind(traits,models,r.ndvi.step2020,R2.ndvi.step2020))
colnames(pred.accuracy.ndvi.step2020)[3:4]=c("r","R2")

### Plotting LASSO regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.las2020, xlab = "Observed Yield",ylab = "Predindvied Yield", main="LASSO regression with ndvi 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.las2020 <-lm(predindvi.ndvi.las2020 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.las2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.las2020 = vector('expression',2)
rp.ndvi.las2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(summary(fit.ndvi.las2020)$adj.r.squared,dig=2)))[2]
rp.ndvi.las2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(predaccuracy.ndvi.lasso2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.las2020)
R2.ndvi.las2020=summary(fit.ndvi.las2020)$adj.r.squared
R2.ndvi.las2020=round(R2.ndvi.las2020,digits=2)
r.ndvi.las2020=round(predaccuracy.ndvi.lasso2020,digits=2)
models="LASSO"
traits="ndvi"
pred.accuracy.ndvi.las2020=data.frame(cbind(traits,models,r.ndvi.las2020,R2.ndvi.las2020))
colnames(pred.accuracy.ndvi.las2020)[3:4]=c("r","R2")

### Plotting ridge regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.rdg2020, xlab = "Observed Yield",ylab = "Predindvied Yield", main="Ridge regression with ndvi 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.rdg2020 <-lm(predindvi.ndvi.rdg2020 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.rdg2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.rdg2020 = vector('expression',2)
rp.ndvi.rdg2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                list(MYVALUE = format(summary(fit.ndvi.rdg2020)$adj.r.squared,dig=2)))[2]
rp.ndvi.rdg2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                list(MYOTHERVALUE = format(predaccuracy.ndvi.ridge2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.rdg2020)
R2.ndvi.rdg2020=summary(fit.ndvi.rdg2020)$adj.r.squared
R2.ndvi.rdg2020=round(R2.ndvi.rdg2020,digits=2)
r.ndvi.rdg2020=round(predaccuracy.ndvi.ridge2020,digits=2)
models="ridge"
traits="ndvi"
pred.accuracy.ndvi.rdg2020=data.frame(cbind(traits,models,r.ndvi.rdg2020,R2.ndvi.rdg2020))
colnames(pred.accuracy.ndvi.rdg2020)[3:4]=c("r","R2")

### Plotting enet regression model with ndvi
plot(ndvi_data$GRYLD, predindvi.ndvi.enet2020, xlab = "Observed Yield",ylab = "Predindvied Yield", main="ElasticNet regression with ndvi 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi.enet2020 <-lm(predindvi.ndvi.enet2020 ~ ndvi_data$GRYLD) #fit the model for observed values
abline(fit.ndvi.enet2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi.enet2020 = vector('expression',2)
rp.ndvi.enet2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.ndvi.enet2020)$adj.r.squared,dig=2)))[2]
rp.ndvi.enet2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.ndvi.enet2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi.enet2020)
R2.ndvi.enet2020=summary(fit.ndvi.enet2020)$adj.r.squared
R2.ndvi.enet2020=round(R2.ndvi.enet2020,digits=2)
r.ndvi.enet2020=round(predaccuracy.ndvi.enet2020,digits=2)
models="elasticnet"
traits="ndvi"
pred.accuracy.ndvi.enet2020=data.frame(cbind(traits,models,r.ndvi.enet2020,R2.ndvi.enet2020))
colnames(pred.accuracy.ndvi.enet2020)[3:4]=c("r","R2")

multivar.ndvi.prediction.accuracy2020=data.frame(rbind(pred.accuracy.ndvi.step2020,pred.accuracy.ndvi.las2020,pred.accuracy.ndvi.rdg2020,pred.accuracy.ndvi.enet2020))
# write.csv(multivar.ndvi.prediction.accuracy2020,file="Multivariate_predaccuracy_ndvi_2020.csv",row.names = FALSE,quote = FALSE)

multivar.ndvi.allyear=cbind(multivar.ndvi.prediction.accuracy2016,multivar.ndvi.prediction.accuracy2017[,3:4],multivar.ndvi.prediction.accuracy2018[,3:4],multivar.ndvi.prediction.accuracy2019[,3:4],multivar.ndvi.prediction.accuracy2020[,3:4])
write.csv(multivar.ndvi.allyear,file="multivar.yieldpred.ndvi.allYears.csv",row.names=FALSE,quote=FALSE)
