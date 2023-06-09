### Multivariate model yield prediction 2016 with ct
### Running stepwise regression model with ct
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
CT_data=BLUE2016[,c(1:10,28)]
stepwise.ct.2016 <- lm(GRYLD ~ CT_20160123+CT_20160204+CT_20160212+CT_20160223+CT_20160228+CT_20160302+CT_20160309+CT_20160315
                       ,data=CT_data)  ## writing the model with ct
stepwise.ct.2016 <- step(stepwise.ct.2016, direction = "both") ## running the stepwise regressing model
formula(stepwise.ct.2016) ## getting the selected variables from the stepwise regression model

### Predicted yield from stepwise regression model with ct
step.predyield.ct.2016=NULL 
for(i in 1:10){
        step.prediction.ct.2016=CT_data[CT_data$trial==i,] #get just one trial to predict
        step.training.ct.2016=CT_data[CT_data$trial!=i,] #get all other trials to train on
        fit.ct.2016 = lm(formula = formula(stepwise.ct.2016), data =step.training.ct.2016) ## directly import variables into the model from stepwise regression 
        step.predvalue.ct.2016=predict(object=fit.ct.2016, newdata=step.prediction.ct.2016) ## Predicted value of individual trial
        predicted.ct.trial2016=data.frame(entry=step.prediction.ct.2016$plot,step.predvalue.ct.2016) ## Predicted value of individual trial in data frame
        step.predyield.ct.2016=rbind(step.predyield.ct.2016, predicted.ct.trial2016) ## Combinig all trial's predicted value into one frame
}
step.predyield.ct.2016 = data.frame(trial=rep(1:10,each=60),step.predyield.ct.2016)
step.predyield.ct.2016 <- merge(step.predyield.ct.2016, CT_data[, c(1:2,11)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extracting r from stepwise regression model with ct
cor.ct.step2016 <- NULL #start dataframe to hold correlation values
for(i in 1:10){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ct.2016[step.predyield.ct.2016$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ct.2016, trial_data$GRYLD) #get correlations
        cor.ct.step2016 <- c(cor.ct.step2016, trial_cor) #write results outside of the for loop
}

### Prediction from shrinkage models
predict.ct.las2016=c()  #Intialize to capture predicted values from LASSO model
predict.ct.rdg2016=c() #Intialize to capture predicted values from Ridge model
predict.ct.enet2016=c() #Intialize to capture predicted values from ElasticNet model
cor.ct.las2016=c() #Intialize to capture correlatio from LASSO model
cor.ct.rdg2016=c() #Intialize to capture correlatio from Ridge model
cor.ct.enet2016=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:10){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ct.data2016 = CT_data[CT_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ct.data2016 = CT_data[CT_data$trial!=trial,] #sets up prediction populations
        rctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ct.model2016 = train(GRYLD ~ ., data = train.ct.data2016[,-c(1,2)], #Looks like we are fitting a lasso model.
                                   method = "lasso",
                                   trControl = rctrl,
                                   preProc = c("center", "scale"))
        ridge.ct.model2016 = train(GRYLD ~ ., data = train.ct.data2016[,-c(1,2)], #looks like we are fitting a ridge model.
                                   method = "ridge",
                                   trControl = rctrl,
                                   preProc = c("center", "scale"))
        enet.ct.model2016 = train(GRYLD ~ ., data = train.ct.data2016[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                  method = "enet",
                                  trControl = rctrl,
                                  preProc = c("center", "scale"))
        #model coefficients
        coef.ct.lasso2016 = lasso.ct.model2016$finalModel #extracting some coefficients from LASSO model
        coef.ct.ridge2016 = ridge.ct.model2016$finalModel #extracting some coefficients from Ridge model
        coef.ct.enet2016 = enet.ct.model2016$finalModel #extracting some coefficients from ElasticNet model
        #prediction
        pred.ct.lasso2016 = predict(lasso.ct.model2016,test.ct.data2016[,-c(1,2)]) #making a predict statement.  This looks promising
        pred.ct.ridge2016 = predict(ridge.ct.model2016,test.ct.data2016[,-c(1,2)]) #appears the same with other models
        pred.ct.elasticnet2016 = predict(enet.ct.model2016,test.ct.data2016[,-c(1,2)])
        #extracting correlations
        cor.ct.l2016=cor(CT_data[CT_data$trial==trial,11],pred.ct.lasso2016) # Correlaiton between predicted value from LASSO regression and BLUE
        cor.ct.las2016=c(cor.ct.las2016,cor.ct.l2016) #writing correlation out to cor.ct.las2016
        cor.ct.r2016=cor(CT_data[CT_data$trial==trial,11],pred.ct.ridge2016) # Correlaiton between predicted value from Ridge regression and BLUE
        cor.ct.rdg2016=c(cor.ct.rdg2016,cor.ct.r2016) #writing correlation out to cor.ct.rdg2016
        cor.ct.e2016=cor(CT_data[CT_data$trial==trial,11],pred.ct.elasticnet2016) # Correlaiton between predicted value from ElasticNet regression and BLUE
        cor.ct.enet2016=c(cor.ct.enet2016,cor.ct.e2016) #writing correlation out to cor.ct.enet2016
        #extracting predicted values
        predict.ct.las2016 = c(predict.ct.las2016,pred.ct.lasso2016) # Combinig all predicted values from LASSO
        predict.ct.rdg2016 = c(predict.ct.rdg2016,pred.ct.ridge2016) # Combinig all predicted values from Ridge
        predict.ct.enet2016 = c(predict.ct.enet2016,pred.ct.elasticnet2016) # Combinig all predicted values from ElasticNet
} #end for loop 

#print prediction accuracies
cor.ct.step2016 #individual trial correlations for step regression
predaccuracy.ct.step2016=mean(cor.ct.step2016) #average correlation from stepwise model
cor.ct.las2016 #individual trial correlations for LASSO regression
predaccuracy.ct.lasso2016=mean(cor.ct.las2016) #average correlation from LASSO model
cor.ct.rdg2016 #individual trial correlations for Ridge regression
predaccuracy.ct.ridge2016=mean(cor.ct.rdg2016) #average correlation from Ridge model
cor.ct.enet2016 #individual trial correlations for ElasticNet regression
predaccuracy.ct.enet2016=mean(cor.ct.enet2016) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ct
plot(step.predyield.ct.2016$GRYLD, step.predyield.ct.2016$step.predvalue.ct.2016, xlab = "Observed Yield",ylab = "Predicted Yield", main="Stepwise regression with CT 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.step2016 <-lm(step.predyield.ct.2016$step.predvalue.ct.2016 ~ step.predyield.ct.2016$GRYLD) #fit the model for observed values
abline(fit.ct.step2016, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.step2016 = vector('expression',2)
rp.ct.step2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                               list(MYVALUE = format(summary(fit.ct.step2016)$adj.r.squared,dig=2)))[2] # extracting R2 value
rp.ct.step2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(predaccuracy.ct.step2016, digits = 2)))[2] # r (r) value
legend("bottomright", bty="n",legend=rp.ct.step2016)
R2.ct.step2016=summary(fit.ct.step2016)$adj.r.squared
R2.ct.step2016=round(R2.ct.step2016,digits=2)
r.ct.step2016=round(predaccuracy.ct.step2016,digits=2)
models="stepwise"
traits="ct"
pred.accuracy.ct.step2016=data.frame(cbind(traits,models,r.ct.step2016,R2.ct.step2016))
colnames(pred.accuracy.ct.step2016)[3:4]=c("r","R2")

### Plotting LASSO regression model with ct
plot(CT_data$GRYLD, predict.ct.las2016, xlab = "Observed Yield",ylab = "Predicted Yield", main="LASSO regression with CT 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.las2016 <-lm(predict.ct.las2016 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.las2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.las2016 = vector('expression',2)
rp.ct.las2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(summary(fit.ct.las2016)$adj.r.squared,dig=2)))[2]
rp.ct.las2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                              list(MYOTHERVALUE = format(predaccuracy.ct.lasso2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.las2016)
R2.ct.las2016=summary(fit.ct.las2016)$adj.r.squared
R2.ct.las2016=round(R2.ct.las2016,digits=2)
r.ct.las2016=round(predaccuracy.ct.lasso2016,digits=2)
models="LASSO"
traits="ct"
pred.accuracy.ct.las2016=data.frame(cbind(traits,models,r.ct.las2016,R2.ct.las2016))
colnames(pred.accuracy.ct.las2016)[3:4]=c("r","R2")

### Plotting ridge regression model with ct
plot(CT_data$GRYLD, predict.ct.rdg2016, xlab = "Observed Yield",ylab = "Predicted Yield", main="Ridge regression with CT 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.rdg2016 <-lm(predict.ct.rdg2016 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.rdg2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.rdg2016 = vector('expression',2)
rp.ct.rdg2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(summary(fit.ct.rdg2016)$adj.r.squared,dig=2)))[2]
rp.ct.rdg2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                              list(MYOTHERVALUE = format(predaccuracy.ct.ridge2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.rdg2016)
R2.ct.rdg2016=summary(fit.ct.rdg2016)$adj.r.squared
R2.ct.rdg2016=round(R2.ct.rdg2016,digits=2)
r.ct.rdg2016=round(predaccuracy.ct.ridge2016,digits=2)
models="ridge"
traits="ct"
pred.accuracy.ct.rdg2016=data.frame(cbind(traits,models,r.ct.rdg2016,R2.ct.rdg2016))
colnames(pred.accuracy.ct.rdg2016)[3:4]=c("r","R2")

### Plotting enet regression model with ct
plot(CT_data$GRYLD, predict.ct.enet2016, xlab = "Observed Yield",ylab = "Predicted Yield", main="ElasticNet regression with CT 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.enet2016 <-lm(predict.ct.enet2016 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.enet2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.enet2016 = vector('expression',2)
rp.ct.enet2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                               list(MYVALUE = format(summary(fit.ct.enet2016)$adj.r.squared,dig=2)))[2]
rp.ct.enet2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(predaccuracy.ct.enet2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.enet2016)
R2.ct.enet2016=summary(fit.ct.enet2016)$adj.r.squared
R2.ct.enet2016=round(R2.ct.enet2016,digits=2)
r.ct.enet2016=round(predaccuracy.ct.enet2016,digits=2)
models="elasticnet"
traits="ct"
pred.accuracy.ct.enet2016=data.frame(cbind(traits,models,r.ct.enet2016,R2.ct.enet2016))
colnames(pred.accuracy.ct.enet2016)[3:4]=c("r","R2")

multivar.ct.prediction.accuracy2016=data.frame(rbind(pred.accuracy.ct.step2016,pred.accuracy.ct.las2016,pred.accuracy.ct.rdg2016,pred.accuracy.ct.enet2016))
# write.csv(multivar.ct.prediction.accuracy2016,file="Multivariate_predaccuracy_ct_2016.csv",row.names = FALSE,quote = FALSE)


### Multivariate model yield prediction 2017 with ct
### Running stepwise regression model with ct
BLUE2017=read.csv("BLUE_2017.csv")
CT_data=BLUE2017[,c(1:16,39)]
stepwise.ct.2017 <- lm(GRYLD ~ CT_20170104+CT_20170109+CT_20170114+CT_20170120+CT_20170125+CT_20170131+CT_20170205+CT_20170210+CT_20170215+CT_20170221+CT_20170225+CT_20170302+CT_20170307+CT_20170313
                       ,data=CT_data)  ## writing the model with ct
stepwise.ct.2017 <- step(stepwise.ct.2017, direction = "both") ## running the stepwise regressing model
formula(stepwise.ct.2017) ## getting the selected variables from the stepwise regression model

### Predicted yield from stepwise regression model with ct
step.predyield.ct.2017=NULL 
for(i in 1:11){
        step.prediction.ct.2017=CT_data[CT_data$trial==i,] #get just one trial to predict
        step.training.ct.2017=CT_data[CT_data$trial!=i,] #get all other trials to train on
        fit.ct.2017 = lm(formula = formula(stepwise.ct.2017), data =step.training.ct.2017) ## directly import variables into the model from stepwise regression 
        step.predvalue.ct.2017=predict(object=fit.ct.2017, newdata=step.prediction.ct.2017) ## Predicted value of individual trial
        predicted.ct.trial2017=data.frame(entry=step.prediction.ct.2017$plot,step.predvalue.ct.2017) ## Predicted value of individual trial in data frame
        step.predyield.ct.2017=rbind(step.predyield.ct.2017, predicted.ct.trial2017) ## Combinig all trial's predicted value into one frame
}
step.predyield.ct.2017 = data.frame(trial=rep(1:11,each=60),step.predyield.ct.2017)
step.predyield.ct.2017 <- merge(step.predyield.ct.2017, CT_data[, c(1:2,17)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extracting r from stepwise regression model with ct
cor.ct.step2017 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ct.2017[step.predyield.ct.2017$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ct.2017, trial_data$GRYLD) #get correlations
        cor.ct.step2017 <- c(cor.ct.step2017, trial_cor) #write results outside of the for loop
}

### Prediction from shrinkage models
predict.ct.las2017=c()  #Intialize to capture predicted values from LASSO model
predict.ct.rdg2017=c() #Intialize to capture predicted values from Ridge model
predict.ct.enet2017=c() #Intialize to capture predicted values from ElasticNet model
cor.ct.las2017=c() #Intialize to capture correlatio from LASSO model
cor.ct.rdg2017=c() #Intialize to capture correlatio from Ridge model
cor.ct.enet2017=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ct.data2017 = CT_data[CT_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ct.data2017 = CT_data[CT_data$trial!=trial,] #sets up prediction populations
        rctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ct.model2017 = train(GRYLD ~ ., data = train.ct.data2017[,-c(1,2)], #Looks like we are fitting a lasso model.
                                   method = "lasso",
                                   trControl = rctrl,
                                   preProc = c("center", "scale"))
        ridge.ct.model2017 = train(GRYLD ~ ., data = train.ct.data2017[,-c(1,2)], #looks like we are fitting a ridge model.
                                   method = "ridge",
                                   trControl = rctrl,
                                   preProc = c("center", "scale"))
        enet.ct.model2017 = train(GRYLD ~ ., data = train.ct.data2017[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                  method = "enet",
                                  trControl = rctrl,
                                  preProc = c("center", "scale"))
        #model coefficients
        coef.ct.lasso2017 = lasso.ct.model2017$finalModel #extracting some coefficients from LASSO model
        coef.ct.ridge2017 = ridge.ct.model2017$finalModel #extracting some coefficients from Ridge model
        coef.ct.enet2017 = enet.ct.model2017$finalModel #extracting some coefficients from ElasticNet model
        #prediction
        pred.ct.lasso2017 = predict(lasso.ct.model2017,test.ct.data2017[,-c(1,2)]) #making a predict statement.  This looks promising
        pred.ct.ridge2017 = predict(ridge.ct.model2017,test.ct.data2017[,-c(1,2)]) #appears the same with other models
        pred.ct.elasticnet2017 = predict(enet.ct.model2017,test.ct.data2017[,-c(1,2)])
        #extracting correlations
        cor.ct.l2017=cor(CT_data[CT_data$trial==trial,17],pred.ct.lasso2017) # Correlaiton between predicted value from LASSO regression and BLUE
        cor.ct.las2017=c(cor.ct.las2017,cor.ct.l2017) #writing correlation out to cor.ct.las2017
        cor.ct.r2017=cor(CT_data[CT_data$trial==trial,17],pred.ct.ridge2017) # Correlaiton between predicted value from Ridge regression and BLUE
        cor.ct.rdg2017=c(cor.ct.rdg2017,cor.ct.r2017) #writing correlation out to cor.ct.rdg2017
        cor.ct.e2017=cor(CT_data[CT_data$trial==trial,17],pred.ct.elasticnet2017) # Correlaiton between predicted value from ElasticNet regression and BLUE
        cor.ct.enet2017=c(cor.ct.enet2017,cor.ct.e2017) #writing correlation out to cor.ct.enet2017
        #extracting predicted values
        predict.ct.las2017 = c(predict.ct.las2017,pred.ct.lasso2017) # Combinig all predicted values from LASSO
        predict.ct.rdg2017 = c(predict.ct.rdg2017,pred.ct.ridge2017) # Combinig all predicted values from Ridge
        predict.ct.enet2017 = c(predict.ct.enet2017,pred.ct.elasticnet2017) # Combinig all predicted values from ElasticNet
} #end for loop 

#print prediction accuracies
cor.ct.step2017 #individual trial correlations for step regression
predaccuracy.ct.step2017=mean(cor.ct.step2017) #average correlation from stepwise model
cor.ct.las2017 #individual trial correlations for LASSO regression
predaccuracy.ct.lasso2017=mean(cor.ct.las2017) #average correlation from LASSO model
cor.ct.rdg2017 #individual trial correlations for Ridge regression
predaccuracy.ct.ridge2017=mean(cor.ct.rdg2017) #average correlation from Ridge model
cor.ct.enet2017 #individual trial correlations for ElasticNet regression
predaccuracy.ct.enet2017=mean(cor.ct.enet2017) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ct
plot(step.predyield.ct.2017$GRYLD, step.predyield.ct.2017$step.predvalue.ct.2017, xlab = "Observed Yield",ylab = "Predicted Yield", main="Stepwise regression with CT 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.step2017 <-lm(step.predyield.ct.2017$step.predvalue.ct.2017 ~ step.predyield.ct.2017$GRYLD) #fit the model for observed values
abline(fit.ct.step2017, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.step2017 = vector('expression',2)
rp.ct.step2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                               list(MYVALUE = format(summary(fit.ct.step2017)$adj.r.squared,dig=2)))[2] # extracting R2 value
rp.ct.step2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(predaccuracy.ct.step2017, digits = 2)))[2] # r (r) value
legend("bottomright", bty="n",legend=rp.ct.step2017)
R2.ct.step2017=summary(fit.ct.step2017)$adj.r.squared
R2.ct.step2017=round(R2.ct.step2017,digits=2)
r.ct.step2017=round(predaccuracy.ct.step2017,digits=2)
models="stepwise"
traits="ct"
pred.accuracy.ct.step2017=data.frame(cbind(traits,models,r.ct.step2017,R2.ct.step2017))
colnames(pred.accuracy.ct.step2017)[3:4]=c("r","R2")

### Plotting LASSO regression model with ct
plot(CT_data$GRYLD, predict.ct.las2017, xlab = "Observed Yield",ylab = "Predicted Yield", main="LASSO regression with CT 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.las2017 <-lm(predict.ct.las2017 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.las2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.las2017 = vector('expression',2)
rp.ct.las2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(summary(fit.ct.las2017)$adj.r.squared,dig=2)))[2]
rp.ct.las2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                              list(MYOTHERVALUE = format(predaccuracy.ct.lasso2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.las2017)
R2.ct.las2017=summary(fit.ct.las2017)$adj.r.squared
R2.ct.las2017=round(R2.ct.las2017,digits=2)
r.ct.las2017=round(predaccuracy.ct.lasso2017,digits=2)
models="LASSO"
traits="ct"
pred.accuracy.ct.las2017=data.frame(cbind(traits,models,r.ct.las2017,R2.ct.las2017))
colnames(pred.accuracy.ct.las2017)[3:4]=c("r","R2")

### Plotting ridge regression model with ct
plot(CT_data$GRYLD, predict.ct.rdg2017, xlab = "Observed Yield",ylab = "Predicted Yield", main="Ridge regression with CT 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.rdg2017 <-lm(predict.ct.rdg2017 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.rdg2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.rdg2017 = vector('expression',2)
rp.ct.rdg2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(summary(fit.ct.rdg2017)$adj.r.squared,dig=2)))[2]
rp.ct.rdg2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                              list(MYOTHERVALUE = format(predaccuracy.ct.ridge2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.rdg2017)
R2.ct.rdg2017=summary(fit.ct.rdg2017)$adj.r.squared
R2.ct.rdg2017=round(R2.ct.rdg2017,digits=2)
r.ct.rdg2017=round(predaccuracy.ct.ridge2017,digits=2)
models="ridge"
traits="ct"
pred.accuracy.ct.rdg2017=data.frame(cbind(traits,models,r.ct.rdg2017,R2.ct.rdg2017))
colnames(pred.accuracy.ct.rdg2017)[3:4]=c("r","R2")

### Plotting enet regression model with ct
plot(CT_data$GRYLD, predict.ct.enet2017, xlab = "Observed Yield",ylab = "Predicted Yield", main="ElasticNet regression with CT 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.enet2017 <-lm(predict.ct.enet2017 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.enet2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.enet2017 = vector('expression',2)
rp.ct.enet2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                               list(MYVALUE = format(summary(fit.ct.enet2017)$adj.r.squared,dig=2)))[2]
rp.ct.enet2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(predaccuracy.ct.enet2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.enet2017)
R2.ct.enet2017=summary(fit.ct.enet2017)$adj.r.squared
R2.ct.enet2017=round(R2.ct.enet2017,digits=2)
r.ct.enet2017=round(predaccuracy.ct.enet2017,digits=2)
models="elasticnet"
traits="ct"
pred.accuracy.ct.enet2017=data.frame(cbind(traits,models,r.ct.enet2017,R2.ct.enet2017))
colnames(pred.accuracy.ct.enet2017)[3:4]=c("r","R2")

multivar.ct.prediction.accuracy2017=data.frame(rbind(pred.accuracy.ct.step2017,pred.accuracy.ct.las2017,pred.accuracy.ct.rdg2017,pred.accuracy.ct.enet2017))
# write.csv(multivar.ct.prediction.accuracy2017,file="Multivariate_predaccuracy_ct_2017.csv",row.names = FALSE,quote = FALSE)




### Multivariate model yield prediction 2018 with ct
### Running stepwise regression model with ct
BLUE2018=read.csv("BLUE_2018.csv")
CT_data=BLUE2018[,c(1:14,35)]
stepwise.ct.2018 <- lm(GRYLD ~ CT_20180126+CT_20180131+CT_20180205+CT_20180210+CT_20180214+CT_20180219+CT_20180225+CT_20180301+CT_20180305+CT_20180310+CT_20180315+CT_20180320
                       ,data=CT_data)  ## writing the model with ct
stepwise.ct.2018 <- step(stepwise.ct.2018, direction = "both") ## running the stepwise regressing model
formula(stepwise.ct.2018) ## getting the selected variables from the stepwise regression model

### Predicted yield from stepwise regression model with ct
step.predyield.ct.2018=NULL 
for(i in 1:11){
        step.prediction.ct.2018=CT_data[CT_data$trial==i,] #get just one trial to predict
        step.training.ct.2018=CT_data[CT_data$trial!=i,] #get all other trials to train on
        fit.ct.2018 = lm(formula = formula(stepwise.ct.2018), data =step.training.ct.2018) ## directly import variables into the model from stepwise regression 
        step.predvalue.ct.2018=predict(object=fit.ct.2018, newdata=step.prediction.ct.2018) ## Predicted value of individual trial
        predicted.ct.trial2018=data.frame(entry=step.prediction.ct.2018$plot,step.predvalue.ct.2018) ## Predicted value of individual trial in data frame
        step.predyield.ct.2018=rbind(step.predyield.ct.2018, predicted.ct.trial2018) ## Combinig all trial's predicted value into one frame
}
step.predyield.ct.2018 = data.frame(trial=rep(1:11,each=60),step.predyield.ct.2018)
step.predyield.ct.2018 <- merge(step.predyield.ct.2018, CT_data[, c(1:2,15)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extracting r from stepwise regression model with ct
cor.ct.step2018 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ct.2018[step.predyield.ct.2018$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ct.2018, trial_data$GRYLD) #get correlations
        cor.ct.step2018 <- c(cor.ct.step2018, trial_cor) #write results outside of the for loop
}

### Prediction from shrinkage models
predict.ct.las2018=c()  #Intialize to capture predicted values from LASSO model
predict.ct.rdg2018=c() #Intialize to capture predicted values from Ridge model
predict.ct.enet2018=c() #Intialize to capture predicted values from ElasticNet model
cor.ct.las2018=c() #Intialize to capture correlatio from LASSO model
cor.ct.rdg2018=c() #Intialize to capture correlatio from Ridge model
cor.ct.enet2018=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ct.data2018 = CT_data[CT_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ct.data2018 = CT_data[CT_data$trial!=trial,] #sets up prediction populations
        rctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ct.model2018 = train(GRYLD ~ ., data = train.ct.data2018[,-c(1,2)], #Looks like we are fitting a lasso model.
                                   method = "lasso",
                                   trControl = rctrl,
                                   preProc = c("center", "scale"))
        ridge.ct.model2018 = train(GRYLD ~ ., data = train.ct.data2018[,-c(1,2)], #looks like we are fitting a ridge model.
                                   method = "ridge",
                                   trControl = rctrl,
                                   preProc = c("center", "scale"))
        enet.ct.model2018 = train(GRYLD ~ ., data = train.ct.data2018[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                  method = "enet",
                                  trControl = rctrl,
                                  preProc = c("center", "scale"))
        #model coefficients
        coef.ct.lasso2018 = lasso.ct.model2018$finalModel #extracting some coefficients from LASSO model
        coef.ct.ridge2018 = ridge.ct.model2018$finalModel #extracting some coefficients from Ridge model
        coef.ct.enet2018 = enet.ct.model2018$finalModel #extracting some coefficients from ElasticNet model
        #prediction
        pred.ct.lasso2018 = predict(lasso.ct.model2018,test.ct.data2018[,-c(1,2)]) #making a predict statement.  This looks promising
        pred.ct.ridge2018 = predict(ridge.ct.model2018,test.ct.data2018[,-c(1,2)]) #appears the same with other models
        pred.ct.elasticnet2018 = predict(enet.ct.model2018,test.ct.data2018[,-c(1,2)])
        #extracting correlations
        cor.ct.l2018=cor(CT_data[CT_data$trial==trial,15],pred.ct.lasso2018) # Correlaiton between predicted value from LASSO regression and BLUE
        cor.ct.las2018=c(cor.ct.las2018,cor.ct.l2018) #writing correlation out to cor.ct.las2018
        cor.ct.r2018=cor(CT_data[CT_data$trial==trial,15],pred.ct.ridge2018) # Correlaiton between predicted value from Ridge regression and BLUE
        cor.ct.rdg2018=c(cor.ct.rdg2018,cor.ct.r2018) #writing correlation out to cor.ct.rdg2018
        cor.ct.e2018=cor(CT_data[CT_data$trial==trial,15],pred.ct.elasticnet2018) # Correlaiton between predicted value from ElasticNet regression and BLUE
        cor.ct.enet2018=c(cor.ct.enet2018,cor.ct.e2018) #writing correlation out to cor.ct.enet2018
        #extracting predicted values
        predict.ct.las2018 = c(predict.ct.las2018,pred.ct.lasso2018) # Combinig all predicted values from LASSO
        predict.ct.rdg2018 = c(predict.ct.rdg2018,pred.ct.ridge2018) # Combinig all predicted values from Ridge
        predict.ct.enet2018 = c(predict.ct.enet2018,pred.ct.elasticnet2018) # Combinig all predicted values from ElasticNet
} #end for loop 

#print prediction accuracies
cor.ct.step2018 #individual trial correlations for step regression
predaccuracy.ct.step2018=mean(cor.ct.step2018) #average correlation from stepwise model
cor.ct.las2018 #individual trial correlations for LASSO regression
predaccuracy.ct.lasso2018=mean(cor.ct.las2018) #average correlation from LASSO model
cor.ct.rdg2018 #individual trial correlations for Ridge regression
predaccuracy.ct.ridge2018=mean(cor.ct.rdg2018) #average correlation from Ridge model
cor.ct.enet2018 #individual trial correlations for ElasticNet regression
predaccuracy.ct.enet2018=mean(cor.ct.enet2018) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ct
plot(step.predyield.ct.2018$GRYLD, step.predyield.ct.2018$step.predvalue.ct.2018, xlab = "Observed Yield",ylab = "Predicted Yield", main="Stepwise regression with CT 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.step2018 <-lm(step.predyield.ct.2018$step.predvalue.ct.2018 ~ step.predyield.ct.2018$GRYLD) #fit the model for observed values
abline(fit.ct.step2018, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.step2018 = vector('expression',2)
rp.ct.step2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                               list(MYVALUE = format(summary(fit.ct.step2018)$adj.r.squared,dig=2)))[2] # extracting R2 value
rp.ct.step2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(predaccuracy.ct.step2018, digits = 2)))[2] # r (r) value
legend("bottomright", bty="n",legend=rp.ct.step2018)
R2.ct.step2018=summary(fit.ct.step2018)$adj.r.squared
R2.ct.step2018=round(R2.ct.step2018,digits=2)
r.ct.step2018=round(predaccuracy.ct.step2018,digits=2)
models="stepwise"
traits="ct"
pred.accuracy.ct.step2018=data.frame(cbind(traits,models,r.ct.step2018,R2.ct.step2018))
colnames(pred.accuracy.ct.step2018)[3:4]=c("r","R2")

### Plotting LASSO regression model with ct
plot(CT_data$GRYLD, predict.ct.las2018, xlab = "Observed Yield",ylab = "Predicted Yield", main="LASSO regression with CT 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.las2018 <-lm(predict.ct.las2018 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.las2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.las2018 = vector('expression',2)
rp.ct.las2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(summary(fit.ct.las2018)$adj.r.squared,dig=2)))[2]
rp.ct.las2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                              list(MYOTHERVALUE = format(predaccuracy.ct.lasso2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.las2018)
R2.ct.las2018=summary(fit.ct.las2018)$adj.r.squared
R2.ct.las2018=round(R2.ct.las2018,digits=2)
r.ct.las2018=round(predaccuracy.ct.lasso2018,digits=2)
models="LASSO"
traits="ct"
pred.accuracy.ct.las2018=data.frame(cbind(traits,models,r.ct.las2018,R2.ct.las2018))
colnames(pred.accuracy.ct.las2018)[3:4]=c("r","R2")

### Plotting ridge regression model with ct
plot(CT_data$GRYLD, predict.ct.rdg2018, xlab = "Observed Yield",ylab = "Predicted Yield", main="Ridge regression with CT 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.rdg2018 <-lm(predict.ct.rdg2018 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.rdg2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.rdg2018 = vector('expression',2)
rp.ct.rdg2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(summary(fit.ct.rdg2018)$adj.r.squared,dig=2)))[2]
rp.ct.rdg2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                              list(MYOTHERVALUE = format(predaccuracy.ct.ridge2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.rdg2018)
R2.ct.rdg2018=summary(fit.ct.rdg2018)$adj.r.squared
R2.ct.rdg2018=round(R2.ct.rdg2018,digits=2)
r.ct.rdg2018=round(predaccuracy.ct.ridge2018,digits=2)
models="ridge"
traits="ct"
pred.accuracy.ct.rdg2018=data.frame(cbind(traits,models,r.ct.rdg2018,R2.ct.rdg2018))
colnames(pred.accuracy.ct.rdg2018)[3:4]=c("r","R2")

### Plotting enet regression model with ct
plot(CT_data$GRYLD, predict.ct.enet2018, xlab = "Observed Yield",ylab = "Predicted Yield", main="ElasticNet regression with CT 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.enet2018 <-lm(predict.ct.enet2018 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.enet2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.enet2018 = vector('expression',2)
rp.ct.enet2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                               list(MYVALUE = format(summary(fit.ct.enet2018)$adj.r.squared,dig=2)))[2]
rp.ct.enet2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(predaccuracy.ct.enet2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.enet2018)
R2.ct.enet2018=summary(fit.ct.enet2018)$adj.r.squared
R2.ct.enet2018=round(R2.ct.enet2018,digits=2)
r.ct.enet2018=round(predaccuracy.ct.enet2018,digits=2)
models="elasticnet"
traits="ct"
pred.accuracy.ct.enet2018=data.frame(cbind(traits,models,r.ct.enet2018,R2.ct.enet2018))
colnames(pred.accuracy.ct.enet2018)[3:4]=c("r","R2")

multivar.ct.prediction.accuracy2018=data.frame(rbind(pred.accuracy.ct.step2018,pred.accuracy.ct.las2018,pred.accuracy.ct.rdg2018,pred.accuracy.ct.enet2018))
# write.csv(multivar.ct.prediction.accuracy2018,file="Multivariate_predaccuracy_ct_2018.csv",row.names = FALSE,quote = FALSE)



### Multivariate model yield prediction 2019 with ct
### Running stepwise regression model with ct
BLUE2019=read.csv("BLUE_2019.csv")
CT_data=BLUE2019[,c(1:15,37)]
stepwise.ct.2019 <- lm(GRYLD ~ CT_20190123+CT_20190127+CT_20190131+CT_20190205+CT_20190211+CT_20190218+CT_20190223+CT_20190301+CT_20190305+CT_20190311+CT_20190316+CT_20190320+CT_20190325
                       ,data=CT_data)  ## writing the model with ct
stepwise.ct.2019 <- step(stepwise.ct.2019, direction = "both") ## running the stepwise regressing model
formula(stepwise.ct.2019) ## getting the selected variables from the stepwise regression model

### Predicted yield from stepwise regression model with ct
step.predyield.ct.2019=NULL 
for(i in 1:10){
        step.prediction.ct.2019=CT_data[CT_data$trial==i,] #get just one trial to predict
        step.training.ct.2019=CT_data[CT_data$trial!=i,] #get all other trials to train on
        fit.ct.2019 = lm(formula = formula(stepwise.ct.2019), data =step.training.ct.2019) ## directly import variables into the model from stepwise regression 
        step.predvalue.ct.2019=predict(object=fit.ct.2019, newdata=step.prediction.ct.2019) ## Predicted value of individual trial
        predicted.ct.trial2019=data.frame(entry=step.prediction.ct.2019$plot,step.predvalue.ct.2019) ## Predicted value of individual trial in data frame
        step.predyield.ct.2019=rbind(step.predyield.ct.2019, predicted.ct.trial2019) ## Combinig all trial's predicted value into one frame
}
step.predyield.ct.2019 = data.frame(trial=rep(1:10,each=60),step.predyield.ct.2019)
step.predyield.ct.2019 <- merge(step.predyield.ct.2019, CT_data[, c(1:2,16)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extracting r from stepwise regression model with ct
cor.ct.step2019 <- NULL #start dataframe to hold correlation values
for(i in 1:10){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ct.2019[step.predyield.ct.2019$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ct.2019, trial_data$GRYLD) #get correlations
        cor.ct.step2019 <- c(cor.ct.step2019, trial_cor) #write results outside of the for loop
}

### Prediction from shrinkage models
predict.ct.las2019=c()  #Intialize to capture predicted values from LASSO model
predict.ct.rdg2019=c() #Intialize to capture predicted values from Ridge model
predict.ct.enet2019=c() #Intialize to capture predicted values from ElasticNet model
cor.ct.las2019=c() #Intialize to capture correlatio from LASSO model
cor.ct.rdg2019=c() #Intialize to capture correlatio from Ridge model
cor.ct.enet2019=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:10){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ct.data2019 = CT_data[CT_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ct.data2019 = CT_data[CT_data$trial!=trial,] #sets up prediction populations
        rctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ct.model2019 = train(GRYLD ~ ., data = train.ct.data2019[,-c(1,2)], #Looks like we are fitting a lasso model.
                                   method = "lasso",
                                   trControl = rctrl,
                                   preProc = c("center", "scale"))
        ridge.ct.model2019 = train(GRYLD ~ ., data = train.ct.data2019[,-c(1,2)], #looks like we are fitting a ridge model.
                                   method = "ridge",
                                   trControl = rctrl,
                                   preProc = c("center", "scale"))
        enet.ct.model2019 = train(GRYLD ~ ., data = train.ct.data2019[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                  method = "enet",
                                  trControl = rctrl,
                                  preProc = c("center", "scale"))
        #model coefficients
        coef.ct.lasso2019 = lasso.ct.model2019$finalModel #extracting some coefficients from LASSO model
        coef.ct.ridge2019 = ridge.ct.model2019$finalModel #extracting some coefficients from Ridge model
        coef.ct.enet2019 = enet.ct.model2019$finalModel #extracting some coefficients from ElasticNet model
        #prediction
        pred.ct.lasso2019 = predict(lasso.ct.model2019,test.ct.data2019[,-c(1,2)]) #making a predict statement.  This looks promising
        pred.ct.ridge2019 = predict(ridge.ct.model2019,test.ct.data2019[,-c(1,2)]) #appears the same with other models
        pred.ct.elasticnet2019 = predict(enet.ct.model2019,test.ct.data2019[,-c(1,2)])
        #extracting correlations
        cor.ct.l2019=cor(CT_data[CT_data$trial==trial,16],pred.ct.lasso2019) # Correlaiton between predicted value from LASSO regression and BLUE
        cor.ct.las2019=c(cor.ct.las2019,cor.ct.l2019) #writing correlation out to cor.ct.las2019
        cor.ct.r2019=cor(CT_data[CT_data$trial==trial,16],pred.ct.ridge2019) # Correlaiton between predicted value from Ridge regression and BLUE
        cor.ct.rdg2019=c(cor.ct.rdg2019,cor.ct.r2019) #writing correlation out to cor.ct.rdg2019
        cor.ct.e2019=cor(CT_data[CT_data$trial==trial,16],pred.ct.elasticnet2019) # Correlaiton between predicted value from ElasticNet regression and BLUE
        cor.ct.enet2019=c(cor.ct.enet2019,cor.ct.e2019) #writing correlation out to cor.ct.enet2019
        #extracting predicted values
        predict.ct.las2019 = c(predict.ct.las2019,pred.ct.lasso2019) # Combinig all predicted values from LASSO
        predict.ct.rdg2019 = c(predict.ct.rdg2019,pred.ct.ridge2019) # Combinig all predicted values from Ridge
        predict.ct.enet2019 = c(predict.ct.enet2019,pred.ct.elasticnet2019) # Combinig all predicted values from ElasticNet
} #end for loop 

#print prediction accuracies
cor.ct.step2019 #individual trial correlations for step regression
predaccuracy.ct.step2019=mean(cor.ct.step2019) #average correlation from stepwise model
cor.ct.las2019 #individual trial correlations for LASSO regression
predaccuracy.ct.lasso2019=mean(cor.ct.las2019) #average correlation from LASSO model
cor.ct.rdg2019 #individual trial correlations for Ridge regression
predaccuracy.ct.ridge2019=mean(cor.ct.rdg2019) #average correlation from Ridge model
cor.ct.enet2019 #individual trial correlations for ElasticNet regression
predaccuracy.ct.enet2019=mean(cor.ct.enet2019) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ct
plot(step.predyield.ct.2019$GRYLD, step.predyield.ct.2019$step.predvalue.ct.2019, xlab = "Observed Yield",ylab = "Predicted Yield", main="Stepwise regression with CT 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.step2019 <-lm(step.predyield.ct.2019$step.predvalue.ct.2019 ~ step.predyield.ct.2019$GRYLD) #fit the model for observed values
abline(fit.ct.step2019, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.step2019 = vector('expression',2)
rp.ct.step2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                               list(MYVALUE = format(summary(fit.ct.step2019)$adj.r.squared,dig=2)))[2] # extracting R2 value
rp.ct.step2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(predaccuracy.ct.step2019, digits = 2)))[2] # r (r) value
legend("bottomright", bty="n",legend=rp.ct.step2019)
R2.ct.step2019=summary(fit.ct.step2019)$adj.r.squared
R2.ct.step2019=round(R2.ct.step2019,digits=2)
r.ct.step2019=round(predaccuracy.ct.step2019,digits=2)
models="stepwise"
traits="ct"
pred.accuracy.ct.step2019=data.frame(cbind(traits,models,r.ct.step2019,R2.ct.step2019))
colnames(pred.accuracy.ct.step2019)[3:4]=c("r","R2")

### Plotting LASSO regression model with ct
plot(CT_data$GRYLD, predict.ct.las2019, xlab = "Observed Yield",ylab = "Predicted Yield", main="LASSO regression with CT 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.las2019 <-lm(predict.ct.las2019 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.las2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.las2019 = vector('expression',2)
rp.ct.las2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(summary(fit.ct.las2019)$adj.r.squared,dig=2)))[2]
rp.ct.las2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                              list(MYOTHERVALUE = format(predaccuracy.ct.lasso2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.las2019)
R2.ct.las2019=summary(fit.ct.las2019)$adj.r.squared
R2.ct.las2019=round(R2.ct.las2019,digits=2)
r.ct.las2019=round(predaccuracy.ct.lasso2019,digits=2)
models="LASSO"
traits="ct"
pred.accuracy.ct.las2019=data.frame(cbind(traits,models,r.ct.las2019,R2.ct.las2019))
colnames(pred.accuracy.ct.las2019)[3:4]=c("r","R2")

### Plotting ridge regression model with ct
plot(CT_data$GRYLD, predict.ct.rdg2019, xlab = "Observed Yield",ylab = "Predicted Yield", main="Ridge regression with CT 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.rdg2019 <-lm(predict.ct.rdg2019 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.rdg2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.rdg2019 = vector('expression',2)
rp.ct.rdg2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(summary(fit.ct.rdg2019)$adj.r.squared,dig=2)))[2]
rp.ct.rdg2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                              list(MYOTHERVALUE = format(predaccuracy.ct.ridge2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.rdg2019)
R2.ct.rdg2019=summary(fit.ct.rdg2019)$adj.r.squared
R2.ct.rdg2019=round(R2.ct.rdg2019,digits=2)
r.ct.rdg2019=round(predaccuracy.ct.ridge2019,digits=2)
models="ridge"
traits="ct"
pred.accuracy.ct.rdg2019=data.frame(cbind(traits,models,r.ct.rdg2019,R2.ct.rdg2019))
colnames(pred.accuracy.ct.rdg2019)[3:4]=c("r","R2")

### Plotting enet regression model with ct
plot(CT_data$GRYLD, predict.ct.enet2019, xlab = "Observed Yield",ylab = "Predicted Yield", main="ElasticNet regression with CT 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.enet2019 <-lm(predict.ct.enet2019 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.enet2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.enet2019 = vector('expression',2)
rp.ct.enet2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                               list(MYVALUE = format(summary(fit.ct.enet2019)$adj.r.squared,dig=2)))[2]
rp.ct.enet2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(predaccuracy.ct.enet2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.enet2019)
R2.ct.enet2019=summary(fit.ct.enet2019)$adj.r.squared
R2.ct.enet2019=round(R2.ct.enet2019,digits=2)
r.ct.enet2019=round(predaccuracy.ct.enet2019,digits=2)
models="elasticnet"
traits="ct"
pred.accuracy.ct.enet2019=data.frame(cbind(traits,models,r.ct.enet2019,R2.ct.enet2019))
colnames(pred.accuracy.ct.enet2019)[3:4]=c("r","R2")

multivar.ct.prediction.accuracy2019=data.frame(rbind(pred.accuracy.ct.step2019,pred.accuracy.ct.las2019,pred.accuracy.ct.rdg2019,pred.accuracy.ct.enet2019))
# write.csv(multivar.ct.prediction.accuracy2019,file="Multivariate_predaccuracy_ct_2019.csv",row.names = FALSE,quote = FALSE)




### Multivariate model yield prediction 2020 with ct
### Running stepwise regression model with ct
BLUE2020=read.csv("BLUE_2020.csv")
CT_data=BLUE2020[,c(1:17,44)]
stepwise.ct.2020 <- lm(GRYLD ~ CT_20200112+CT_20200116+CT_20200121+CT_20200126+CT_20200130+CT_20200205+CT_20200210+CT_20200215+CT_20200220+CT_20200226+CT_20200302+CT_20200308+CT_20200313+CT_20200318+CT_20200323
                       ,data=CT_data)  ## writing the model with ct
stepwise.ct.2020 <- step(stepwise.ct.2020, direction = "both") ## running the stepwise regressing model
formula(stepwise.ct.2020) ## getting the selected variables from the stepwise regression model

### Predicted yield from stepwise regression model with ct
step.predyield.ct.2020=NULL 
for(i in 1:11){
        step.prediction.ct.2020=CT_data[CT_data$trial==i,] #get just one trial to predict
        step.training.ct.2020=CT_data[CT_data$trial!=i,] #get all other trials to train on
        fit.ct.2020 = lm(formula = formula(stepwise.ct.2020), data =step.training.ct.2020) ## directly import variables into the model from stepwise regression 
        step.predvalue.ct.2020=predict(object=fit.ct.2020, newdata=step.prediction.ct.2020) ## Predicted value of individual trial
        predicted.ct.trial2020=data.frame(entry=step.prediction.ct.2020$plot,step.predvalue.ct.2020) ## Predicted value of individual trial in data frame
        step.predyield.ct.2020=rbind(step.predyield.ct.2020, predicted.ct.trial2020) ## Combinig all trial's predicted value into one frame
}
step.predyield.ct.2020 = data.frame(trial=rep(1:11,each=60),step.predyield.ct.2020)
step.predyield.ct.2020 <- merge(step.predyield.ct.2020, CT_data[, c(1:2,18)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extracting r from stepwise regression model with ct
cor.ct.step2020 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ct.2020[step.predyield.ct.2020$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ct.2020, trial_data$GRYLD) #get correlations
        cor.ct.step2020 <- c(cor.ct.step2020, trial_cor) #write results outside of the for loop
}

### Prediction from shrinkage models
predict.ct.las2020=c()  #Intialize to capture predicted values from LASSO model
predict.ct.rdg2020=c() #Intialize to capture predicted values from Ridge model
predict.ct.enet2020=c() #Intialize to capture predicted values from ElasticNet model
cor.ct.las2020=c() #Intialize to capture correlatio from LASSO model
cor.ct.rdg2020=c() #Intialize to capture correlatio from Ridge model
cor.ct.enet2020=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ct.data2020 = CT_data[CT_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ct.data2020 = CT_data[CT_data$trial!=trial,] #sets up prediction populations
        rctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ct.model2020 = train(GRYLD ~ ., data = train.ct.data2020[,-c(1,2)], #Looks like we are fitting a lasso model.
                                   method = "lasso",
                                   trControl = rctrl,
                                   preProc = c("center", "scale"))
        ridge.ct.model2020 = train(GRYLD ~ ., data = train.ct.data2020[,-c(1,2)], #looks like we are fitting a ridge model.
                                   method = "ridge",
                                   trControl = rctrl,
                                   preProc = c("center", "scale"))
        enet.ct.model2020 = train(GRYLD ~ ., data = train.ct.data2020[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                  method = "enet",
                                  trControl = rctrl,
                                  preProc = c("center", "scale"))
        #model coefficients
        coef.ct.lasso2020 = lasso.ct.model2020$finalModel #extracting some coefficients from LASSO model
        coef.ct.ridge2020 = ridge.ct.model2020$finalModel #extracting some coefficients from Ridge model
        coef.ct.enet2020 = enet.ct.model2020$finalModel #extracting some coefficients from ElasticNet model
        #prediction
        pred.ct.lasso2020 = predict(lasso.ct.model2020,test.ct.data2020[,-c(1,2)]) #making a predict statement.  This looks promising
        pred.ct.ridge2020 = predict(ridge.ct.model2020,test.ct.data2020[,-c(1,2)]) #appears the same with other models
        pred.ct.elasticnet2020 = predict(enet.ct.model2020,test.ct.data2020[,-c(1,2)])
        #extracting correlations
        cor.ct.l2020=cor(CT_data[CT_data$trial==trial,18],pred.ct.lasso2020) # Correlaiton between predicted value from LASSO regression and BLUE
        cor.ct.las2020=c(cor.ct.las2020,cor.ct.l2020) #writing correlation out to cor.ct.las2020
        cor.ct.r2020=cor(CT_data[CT_data$trial==trial,18],pred.ct.ridge2020) # Correlaiton between predicted value from Ridge regression and BLUE
        cor.ct.rdg2020=c(cor.ct.rdg2020,cor.ct.r2020) #writing correlation out to cor.ct.rdg2020
        cor.ct.e2020=cor(CT_data[CT_data$trial==trial,18],pred.ct.elasticnet2020) # Correlaiton between predicted value from ElasticNet regression and BLUE
        cor.ct.enet2020=c(cor.ct.enet2020,cor.ct.e2020) #writing correlation out to cor.ct.enet2020
        #extracting predicted values
        predict.ct.las2020 = c(predict.ct.las2020,pred.ct.lasso2020) # Combinig all predicted values from LASSO
        predict.ct.rdg2020 = c(predict.ct.rdg2020,pred.ct.ridge2020) # Combinig all predicted values from Ridge
        predict.ct.enet2020 = c(predict.ct.enet2020,pred.ct.elasticnet2020) # Combinig all predicted values from ElasticNet
} #end for loop 

#print prediction accuracies
cor.ct.step2020 #individual trial correlations for step regression
predaccuracy.ct.step2020=mean(cor.ct.step2020) #average correlation from stepwise model
cor.ct.las2020 #individual trial correlations for LASSO regression
predaccuracy.ct.lasso2020=mean(cor.ct.las2020) #average correlation from LASSO model
cor.ct.rdg2020 #individual trial correlations for Ridge regression
predaccuracy.ct.ridge2020=mean(cor.ct.rdg2020) #average correlation from Ridge model
cor.ct.enet2020 #individual trial correlations for ElasticNet regression
predaccuracy.ct.enet2020=mean(cor.ct.enet2020) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ct
plot(step.predyield.ct.2020$GRYLD, step.predyield.ct.2020$step.predvalue.ct.2020, xlab = "Observed Yield",ylab = "Predicted Yield", main="Stepwise regression with CT 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.step2020 <-lm(step.predyield.ct.2020$step.predvalue.ct.2020 ~ step.predyield.ct.2020$GRYLD) #fit the model for observed values
abline(fit.ct.step2020, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.step2020 = vector('expression',2)
rp.ct.step2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                               list(MYVALUE = format(summary(fit.ct.step2020)$adj.r.squared,dig=2)))[2] # extracting R2 value
rp.ct.step2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(predaccuracy.ct.step2020, digits = 2)))[2] # r (r) value
legend("bottomright", bty="n",legend=rp.ct.step2020)
R2.ct.step2020=summary(fit.ct.step2020)$adj.r.squared
R2.ct.step2020=round(R2.ct.step2020,digits=2)
r.ct.step2020=round(predaccuracy.ct.step2020,digits=2)
models="stepwise"
traits="ct"
pred.accuracy.ct.step2020=data.frame(cbind(traits,models,r.ct.step2020,R2.ct.step2020))
colnames(pred.accuracy.ct.step2020)[3:4]=c("r","R2")

### Plotting LASSO regression model with ct
plot(CT_data$GRYLD, predict.ct.las2020, xlab = "Observed Yield",ylab = "Predicted Yield", main="LASSO regression with CT 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.las2020 <-lm(predict.ct.las2020 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.las2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.las2020 = vector('expression',2)
rp.ct.las2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(summary(fit.ct.las2020)$adj.r.squared,dig=2)))[2]
rp.ct.las2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                              list(MYOTHERVALUE = format(predaccuracy.ct.lasso2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.las2020)
R2.ct.las2020=summary(fit.ct.las2020)$adj.r.squared
R2.ct.las2020=round(R2.ct.las2020,digits=2)
r.ct.las2020=round(predaccuracy.ct.lasso2020,digits=2)
models="LASSO"
traits="ct"
pred.accuracy.ct.las2020=data.frame(cbind(traits,models,r.ct.las2020,R2.ct.las2020))
colnames(pred.accuracy.ct.las2020)[3:4]=c("r","R2")

### Plotting ridge regression model with ct
plot(CT_data$GRYLD, predict.ct.rdg2020, xlab = "Observed Yield",ylab = "Predicted Yield", main="Ridge regression with CT 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.rdg2020 <-lm(predict.ct.rdg2020 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.rdg2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.rdg2020 = vector('expression',2)
rp.ct.rdg2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                              list(MYVALUE = format(summary(fit.ct.rdg2020)$adj.r.squared,dig=2)))[2]
rp.ct.rdg2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                              list(MYOTHERVALUE = format(predaccuracy.ct.ridge2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.rdg2020)
R2.ct.rdg2020=summary(fit.ct.rdg2020)$adj.r.squared
R2.ct.rdg2020=round(R2.ct.rdg2020,digits=2)
r.ct.rdg2020=round(predaccuracy.ct.ridge2020,digits=2)
models="ridge"
traits="ct"
pred.accuracy.ct.rdg2020=data.frame(cbind(traits,models,r.ct.rdg2020,R2.ct.rdg2020))
colnames(pred.accuracy.ct.rdg2020)[3:4]=c("r","R2")

### Plotting enet regression model with ct
plot(CT_data$GRYLD, predict.ct.enet2020, xlab = "Observed Yield",ylab = "Predicted Yield", main="ElasticNet regression with CT 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ct.enet2020 <-lm(predict.ct.enet2020 ~ CT_data$GRYLD) #fit the model for observed values
abline(fit.ct.enet2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ct.enet2020 = vector('expression',2)
rp.ct.enet2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                               list(MYVALUE = format(summary(fit.ct.enet2020)$adj.r.squared,dig=2)))[2]
rp.ct.enet2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                               list(MYOTHERVALUE = format(predaccuracy.ct.enet2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ct.enet2020)
R2.ct.enet2020=summary(fit.ct.enet2020)$adj.r.squared
R2.ct.enet2020=round(R2.ct.enet2020,digits=2)
r.ct.enet2020=round(predaccuracy.ct.enet2020,digits=2)
models="elasticnet"
traits="ct"
pred.accuracy.ct.enet2020=data.frame(cbind(traits,models,r.ct.enet2020,R2.ct.enet2020))
colnames(pred.accuracy.ct.enet2020)[3:4]=c("r","R2")

multivar.ct.prediction.accuracy2020=data.frame(rbind(pred.accuracy.ct.step2020,pred.accuracy.ct.las2020,pred.accuracy.ct.rdg2020,pred.accuracy.ct.enet2020))
# write.csv(multivar.ct.prediction.accuracy2020,file="Multivariate_predaccuracy_ct_2020.csv",row.names = FALSE,quote = FALSE)

multivar.yieldpred.ct.allYear=cbind(multivar.ct.prediction.accuracy2016,multivar.ct.prediction.accuracy2017[,3:4],multivar.ct.prediction.accuracy2018[,3:4],multivar.ct.prediction.accuracy2019[,3:4],multivar.ct.prediction.accuracy2020[,3:4])
write.csv(multivar.yieldpred.ct.allYear,file="Multivar.yieldpred.ct.allYear.csv",row.names=FALSE,quote=FALSE)

