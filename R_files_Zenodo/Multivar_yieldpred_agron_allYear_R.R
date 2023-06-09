### Multivariate model yield prediction 2016 with agron
### Running stepwise regression model with agron
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
agron_data=BLUE2016[,c(1:2,20:28)]
stepwise.agron.2016 <- lm(GRYLD ~ DTHD+DAYSMT+PH+SN+SPKLNG+SPLN+GRNSPK+TGW,data=agron_data)  ## writing the model with agron
stepwise.agron.2016 <- step(stepwise.agron.2016, direagronion = "both") ## running the stepwise regressing model
formula(stepwise.agron.2016) ## getting the seleagroned variables from the stepwise regression model

### Prediagroned yield from stepwise regression model with agron
step.predyield.agron.2016=NULL 
for(i in 1:10){
        step.prediction.agron.2016=agron_data[agron_data$trial==i,] #get just one trial to prediagron
        step.training.agron.2016=agron_data[agron_data$trial!=i,] #get all other trials to train on
        fit.agron.2016 = lm(formula = formula(stepwise.agron.2016), data =step.training.agron.2016) ## direagronly import variables into the model from stepwise regression 
        step.predvalue.agron.2016=predict(object=fit.agron.2016, newdata=step.prediction.agron.2016) ## Prediagroned value of individual trial
        predicted.agron.trial2016=data.frame(entry=step.prediction.agron.2016$plot,step.predvalue.agron.2016) ## Prediagroned value of individual trial in data frame
        step.predyield.agron.2016=rbind(step.predyield.agron.2016, predicted.agron.trial2016) ## Combinig all trial's prediagroned value into one frame
}
step.predyield.agron.2016 = data.frame(trial=rep(1:10,each=60),step.predyield.agron.2016)
step.predyield.agron.2016 <- merge(step.predyield.agron.2016, agron_data[, c(1:2,11)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extraagroning prediagronion accuracy from stepwise regression model with agron
cor.agron.step2016 <- NULL #start dataframe to hold correlation values
for(i in 1:10){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.agron.2016[step.predyield.agron.2016$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.agron.2016, trial_data$GRYLD) #get correlations
        cor.agron.step2016 <- c(cor.agron.step2016, trial_cor) #write results outside of the for loop
}

### Prediagronion from shrinkage models
prediagron.agron.las2016=c()  #Intialize to capture prediagroned values from LASSO model
prediagron.agron.rdg2016=c() #Intialize to capture prediagroned values from Ridge model
prediagron.agron.enet2016=c() #Intialize to capture prediagroned values from ElasticNet model
cor.agron.las2016=c() #Intialize to capture correlatio from LASSO model
cor.agron.rdg2016=c() #Intialize to capture correlatio from Ridge model
cor.agron.enet2016=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:10){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.agron.data2016 = agron_data[agron_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.agron.data2016 = agron_data[agron_data$trial!=trial,] #sets up prediagronion populations
        ragronrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.agron.model2016 = train(GRYLD ~ ., data = train.agron.data2016[,-c(1,2)], #Looks like we are fitting a lasso model.
                                      method = "lasso",
                                      trControl = ragronrl,
                                      preProc = c("center", "scale"))
        ridge.agron.model2016 = train(GRYLD ~ ., data = train.agron.data2016[,-c(1,2)], #looks like we are fitting a ridge model.
                                      method = "ridge",
                                      trControl = ragronrl,
                                      preProc = c("center", "scale"))
        enet.agron.model2016 = train(GRYLD ~ ., data = train.agron.data2016[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                     method = "enet",
                                     trControl = ragronrl,
                                     preProc = c("center", "scale"))
        #model coefficients
        coef.agron.lasso2016 = lasso.agron.model2016$finalModel #extraagroning some coefficients from LASSO model
        coef.agron.ridge2016 = ridge.agron.model2016$finalModel #extraagroning some coefficients from Ridge model
        coef.agron.enet2016 = enet.agron.model2016$finalModel #extraagroning some coefficients from ElasticNet model
        #prediagronion
        pred.agron.lasso2016 = predict(lasso.agron.model2016,test.agron.data2016[,-c(1,2)]) #making a prediagron statement.  This looks promising
        pred.agron.ridge2016 = predict(ridge.agron.model2016,test.agron.data2016[,-c(1,2)]) #appears the same with other models
        pred.agron.elasticnet2016 = predict(enet.agron.model2016,test.agron.data2016[,-c(1,2)])
        #extraagroning correlations
        cor.agron.l2016=cor(agron_data[agron_data$trial==trial,11],pred.agron.lasso2016) # Correlaiton between prediagroned value from LASSO regression and BLUE
        cor.agron.las2016=c(cor.agron.las2016,cor.agron.l2016) #writing correlation out to cor.agron.las2016
        cor.agron.r2016=cor(agron_data[agron_data$trial==trial,11],pred.agron.ridge2016) # Correlaiton between prediagroned value from Ridge regression and BLUE
        cor.agron.rdg2016=c(cor.agron.rdg2016,cor.agron.r2016) #writing correlation out to cor.agron.rdg2016
        cor.agron.e2016=cor(agron_data[agron_data$trial==trial,11],pred.agron.elasticnet2016) # Correlaiton between prediagroned value from ElasticNet regression and BLUE
        cor.agron.enet2016=c(cor.agron.enet2016,cor.agron.e2016) #writing correlation out to cor.agron.enet2016
        #extraagroning prediagroned values
        prediagron.agron.las2016 = c(prediagron.agron.las2016,pred.agron.lasso2016) # Combinig all prediagroned values from LASSO
        prediagron.agron.rdg2016 = c(prediagron.agron.rdg2016,pred.agron.ridge2016) # Combinig all prediagroned values from Ridge
        prediagron.agron.enet2016 = c(prediagron.agron.enet2016,pred.agron.elasticnet2016) # Combinig all prediagroned values from ElasticNet
} #end for loop 

#print prediagronion accuracies
cor.agron.step2016 #individual trial correlations for step regression
predaccuracy.agron.step2016=mean(cor.agron.step2016) #average correlation from stepwise model
cor.agron.las2016 #individual trial correlations for LASSO regression
predaccuracy.agron.lasso2016=mean(cor.agron.las2016) #average correlation from LASSO model
cor.agron.rdg2016 #individual trial correlations for Ridge regression
predaccuracy.agron.ridge2016=mean(cor.agron.rdg2016) #average correlation from Ridge model
cor.agron.enet2016 #individual trial correlations for ElasticNet regression
predaccuracy.agron.enet2016=mean(cor.agron.enet2016) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with agron
plot(step.predyield.agron.2016$GRYLD, step.predyield.agron.2016$step.predvalue.agron.2016, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="Stepwise regression with agron 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.step2016 <-lm(step.predyield.agron.2016$step.predvalue.agron.2016 ~ step.predyield.agron.2016$GRYLD) #fit the model for observed values
abline(fit.agron.step2016, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.step2016 = vector('expression',2)
rp.agron.step2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(summary(fit.agron.step2016)$adj.r.squared,dig=2)))[2] # extraagroning R2 value
rp.agron.step2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(predaccuracy.agron.step2016, digits = 2)))[2] # r (prediagronion accuracy) value
legend("bottomright", bty="n",legend=rp.agron.step2016)
R2.agron.step2016=summary(fit.agron.step2016)$adj.r.squared
R2.agron.step2016=round(R2.agron.step2016,digits=2)
r.agron.step2016=round(predaccuracy.agron.step2016,digits=2)
models="stepwise"
traits="agron"
pred.accuracy.agron.step2016=data.frame(cbind(traits,models,r.agron.step2016,R2.agron.step2016))
colnames(pred.accuracy.agron.step2016)[3:4]=c("prediction accuracy","R2")

### Plotting LASSO regression model with agron
plot(agron_data$GRYLD, prediagron.agron.las2016, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="LASSO regression with agron 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.las2016 <-lm(prediagron.agron.las2016 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.las2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.las2016 = vector('expression',2)
rp.agron.las2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.agron.las2016)$adj.r.squared,dig=2)))[2]
rp.agron.las2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.agron.lasso2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.las2016)
R2.agron.las2016=summary(fit.agron.las2016)$adj.r.squared
R2.agron.las2016=round(R2.agron.las2016,digits=2)
r.agron.las2016=round(predaccuracy.agron.lasso2016,digits=2)
models="LASSO"
traits="agron"
pred.accuracy.agron.las2016=data.frame(cbind(traits,models,r.agron.las2016,R2.agron.las2016))
colnames(pred.accuracy.agron.las2016)[3:4]=c("prediction accuracy","R2")

### Plotting ridge regression model with agron
plot(agron_data$GRYLD, prediagron.agron.rdg2016, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="Ridge regression with agron 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.rdg2016 <-lm(prediagron.agron.rdg2016 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.rdg2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.rdg2016 = vector('expression',2)
rp.agron.rdg2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.agron.rdg2016)$adj.r.squared,dig=2)))[2]
rp.agron.rdg2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.agron.ridge2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.rdg2016)
R2.agron.rdg2016=summary(fit.agron.rdg2016)$adj.r.squared
R2.agron.rdg2016=round(R2.agron.rdg2016,digits=2)
r.agron.rdg2016=round(predaccuracy.agron.ridge2016,digits=2)
models="ridge"
traits="agron"
pred.accuracy.agron.rdg2016=data.frame(cbind(traits,models,r.agron.rdg2016,R2.agron.rdg2016))
colnames(pred.accuracy.agron.rdg2016)[3:4]=c("prediction accuracy","R2")

### Plotting enet regression model with agron
plot(agron_data$GRYLD, prediagron.agron.enet2016, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="ElasticNet regression with agron 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.enet2016 <-lm(prediagron.agron.enet2016 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.enet2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.enet2016 = vector('expression',2)
rp.agron.enet2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(summary(fit.agron.enet2016)$adj.r.squared,dig=2)))[2]
rp.agron.enet2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(predaccuracy.agron.enet2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.enet2016)
R2.agron.enet2016=summary(fit.agron.enet2016)$adj.r.squared
R2.agron.enet2016=round(R2.agron.enet2016,digits=2)
r.agron.enet2016=round(predaccuracy.agron.enet2016,digits=2)
models="elasticnet"
traits="agron"
pred.accuracy.agron.enet2016=data.frame(cbind(traits,models,r.agron.enet2016,R2.agron.enet2016))
colnames(pred.accuracy.agron.enet2016)[3:4]=c("prediction accuracy","R2")

multivar.agron.prediction.accuracy2016=data.frame(rbind(pred.accuracy.agron.step2016,pred.accuracy.agron.las2016,pred.accuracy.agron.rdg2016,pred.accuracy.agron.enet2016))
# write.csv(multivar.agron.prediction.accuracy2016,file="Multivariate_predaccuracy_agron_2016.csv",row.names = FALSE,quote = FALSE)




### Multivariate model yield prediction 2017 with agron
### Running stepwise regression model with agron
BLUE2017=read.csv("BLUE_2017.csv")
agron_data=BLUE2017[,c(1:2,31:39)]
stepwise.agron.2017 <- lm(GRYLD ~ DTHD+DAYSMT+PH+SN+SPKLNG+SPLN+GRNSPK+TGW,data=agron_data)  ## writing the model with agron
stepwise.agron.2017 <- step(stepwise.agron.2017, direagronion = "both") ## running the stepwise regressing model
formula(stepwise.agron.2017) ## getting the seleagroned variables from the stepwise regression model

### Prediagroned yield from stepwise regression model with agron
step.predyield.agron.2017=NULL 
for(i in 1:11){
        step.prediction.agron.2017=agron_data[agron_data$trial==i,] #get just one trial to prediagron
        step.training.agron.2017=agron_data[agron_data$trial!=i,] #get all other trials to train on
        fit.agron.2017 = lm(formula = formula(stepwise.agron.2017), data =step.training.agron.2017) ## direagronly import variables into the model from stepwise regression 
        step.predvalue.agron.2017=predict(object=fit.agron.2017, newdata=step.prediction.agron.2017) ## Prediagroned value of individual trial
        predicted.agron.trial2017=data.frame(entry=step.prediction.agron.2017$plot,step.predvalue.agron.2017) ## Prediagroned value of individual trial in data frame
        step.predyield.agron.2017=rbind(step.predyield.agron.2017, predicted.agron.trial2017) ## Combinig all trial's prediagroned value into one frame
}
step.predyield.agron.2017 = data.frame(trial=rep(1:11,each=60),step.predyield.agron.2017)
step.predyield.agron.2017 <- merge(step.predyield.agron.2017, agron_data[, c(1:2,11)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extraagroning prediagronion accuracy from stepwise regression model with agron
cor.agron.step2017 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.agron.2017[step.predyield.agron.2017$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.agron.2017, trial_data$GRYLD) #get correlations
        cor.agron.step2017 <- c(cor.agron.step2017, trial_cor) #write results outside of the for loop
}

### Prediagronion from shrinkage models
prediagron.agron.las2017=c()  #Intialize to capture prediagroned values from LASSO model
prediagron.agron.rdg2017=c() #Intialize to capture prediagroned values from Ridge model
prediagron.agron.enet2017=c() #Intialize to capture prediagroned values from ElasticNet model
cor.agron.las2017=c() #Intialize to capture correlatio from LASSO model
cor.agron.rdg2017=c() #Intialize to capture correlatio from Ridge model
cor.agron.enet2017=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.agron.data2017 = agron_data[agron_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.agron.data2017 = agron_data[agron_data$trial!=trial,] #sets up prediagronion populations
        ragronrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.agron.model2017 = train(GRYLD ~ ., data = train.agron.data2017[,-c(1,2)], #Looks like we are fitting a lasso model.
                                      method = "lasso",
                                      trControl = ragronrl,
                                      preProc = c("center", "scale"))
        ridge.agron.model2017 = train(GRYLD ~ ., data = train.agron.data2017[,-c(1,2)], #looks like we are fitting a ridge model.
                                      method = "ridge",
                                      trControl = ragronrl,
                                      preProc = c("center", "scale"))
        enet.agron.model2017 = train(GRYLD ~ ., data = train.agron.data2017[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                     method = "enet",
                                     trControl = ragronrl,
                                     preProc = c("center", "scale"))
        #model coefficients
        coef.agron.lasso2017 = lasso.agron.model2017$finalModel #extraagroning some coefficients from LASSO model
        coef.agron.ridge2017 = ridge.agron.model2017$finalModel #extraagroning some coefficients from Ridge model
        coef.agron.enet2017 = enet.agron.model2017$finalModel #extraagroning some coefficients from ElasticNet model
        #prediagronion
        pred.agron.lasso2017 = predict(lasso.agron.model2017,test.agron.data2017[,-c(1,2)]) #making a prediagron statement.  This looks promising
        pred.agron.ridge2017 = predict(ridge.agron.model2017,test.agron.data2017[,-c(1,2)]) #appears the same with other models
        pred.agron.elasticnet2017 = predict(enet.agron.model2017,test.agron.data2017[,-c(1,2)])
        #extraagroning correlations
        cor.agron.l2017=cor(agron_data[agron_data$trial==trial,11],pred.agron.lasso2017) # Correlaiton between prediagroned value from LASSO regression and BLUE
        cor.agron.las2017=c(cor.agron.las2017,cor.agron.l2017) #writing correlation out to cor.agron.las2017
        cor.agron.r2017=cor(agron_data[agron_data$trial==trial,11],pred.agron.ridge2017) # Correlaiton between prediagroned value from Ridge regression and BLUE
        cor.agron.rdg2017=c(cor.agron.rdg2017,cor.agron.r2017) #writing correlation out to cor.agron.rdg2017
        cor.agron.e2017=cor(agron_data[agron_data$trial==trial,11],pred.agron.elasticnet2017) # Correlaiton between prediagroned value from ElasticNet regression and BLUE
        cor.agron.enet2017=c(cor.agron.enet2017,cor.agron.e2017) #writing correlation out to cor.agron.enet2017
        #extraagroning prediagroned values
        prediagron.agron.las2017 = c(prediagron.agron.las2017,pred.agron.lasso2017) # Combinig all prediagroned values from LASSO
        prediagron.agron.rdg2017 = c(prediagron.agron.rdg2017,pred.agron.ridge2017) # Combinig all prediagroned values from Ridge
        prediagron.agron.enet2017 = c(prediagron.agron.enet2017,pred.agron.elasticnet2017) # Combinig all prediagroned values from ElasticNet
} #end for loop 

#print prediagronion accuracies
cor.agron.step2017 #individual trial correlations for step regression
predaccuracy.agron.step2017=mean(cor.agron.step2017) #average correlation from stepwise model
cor.agron.las2017 #individual trial correlations for LASSO regression
predaccuracy.agron.lasso2017=mean(cor.agron.las2017) #average correlation from LASSO model
cor.agron.rdg2017 #individual trial correlations for Ridge regression
predaccuracy.agron.ridge2017=mean(cor.agron.rdg2017) #average correlation from Ridge model
cor.agron.enet2017 #individual trial correlations for ElasticNet regression
predaccuracy.agron.enet2017=mean(cor.agron.enet2017) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with agron
plot(step.predyield.agron.2017$GRYLD, step.predyield.agron.2017$step.predvalue.agron.2017, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="Stepwise regression with agron 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.step2017 <-lm(step.predyield.agron.2017$step.predvalue.agron.2017 ~ step.predyield.agron.2017$GRYLD) #fit the model for observed values
abline(fit.agron.step2017, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.step2017 = vector('expression',2)
rp.agron.step2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(summary(fit.agron.step2017)$adj.r.squared,dig=2)))[2] # extraagroning R2 value
rp.agron.step2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(predaccuracy.agron.step2017, digits = 2)))[2] # r (prediagronion accuracy) value
legend("bottomright", bty="n",legend=rp.agron.step2017)
R2.agron.step2017=summary(fit.agron.step2017)$adj.r.squared
R2.agron.step2017=round(R2.agron.step2017,digits=2)
r.agron.step2017=round(predaccuracy.agron.step2017,digits=2)
models="stepwise"
traits="agron"
pred.accuracy.agron.step2017=data.frame(cbind(traits,models,r.agron.step2017,R2.agron.step2017))
colnames(pred.accuracy.agron.step2017)[3:4]=c("prediction accuracy","R2")

### Plotting LASSO regression model with agron
plot(agron_data$GRYLD, prediagron.agron.las2017, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="LASSO regression with agron 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.las2017 <-lm(prediagron.agron.las2017 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.las2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.las2017 = vector('expression',2)
rp.agron.las2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.agron.las2017)$adj.r.squared,dig=2)))[2]
rp.agron.las2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.agron.lasso2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.las2017)
R2.agron.las2017=summary(fit.agron.las2017)$adj.r.squared
R2.agron.las2017=round(R2.agron.las2017,digits=2)
r.agron.las2017=round(predaccuracy.agron.lasso2017,digits=2)
models="LASSO"
traits="agron"
pred.accuracy.agron.las2017=data.frame(cbind(traits,models,r.agron.las2017,R2.agron.las2017))
colnames(pred.accuracy.agron.las2017)[3:4]=c("prediction accuracy","R2")

### Plotting ridge regression model with agron
plot(agron_data$GRYLD, prediagron.agron.rdg2017, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="Ridge regression with agron 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.rdg2017 <-lm(prediagron.agron.rdg2017 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.rdg2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.rdg2017 = vector('expression',2)
rp.agron.rdg2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.agron.rdg2017)$adj.r.squared,dig=2)))[2]
rp.agron.rdg2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.agron.ridge2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.rdg2017)
R2.agron.rdg2017=summary(fit.agron.rdg2017)$adj.r.squared
R2.agron.rdg2017=round(R2.agron.rdg2017,digits=2)
r.agron.rdg2017=round(predaccuracy.agron.ridge2017,digits=2)
models="ridge"
traits="agron"
pred.accuracy.agron.rdg2017=data.frame(cbind(traits,models,r.agron.rdg2017,R2.agron.rdg2017))
colnames(pred.accuracy.agron.rdg2017)[3:4]=c("prediction accuracy","R2")

### Plotting enet regression model with agron
plot(agron_data$GRYLD, prediagron.agron.enet2017, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="ElasticNet regression with agron 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.enet2017 <-lm(prediagron.agron.enet2017 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.enet2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.enet2017 = vector('expression',2)
rp.agron.enet2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(summary(fit.agron.enet2017)$adj.r.squared,dig=2)))[2]
rp.agron.enet2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(predaccuracy.agron.enet2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.enet2017)
R2.agron.enet2017=summary(fit.agron.enet2017)$adj.r.squared
R2.agron.enet2017=round(R2.agron.enet2017,digits=2)
r.agron.enet2017=round(predaccuracy.agron.enet2017,digits=2)
models="elasticnet"
traits="agron"
pred.accuracy.agron.enet2017=data.frame(cbind(traits,models,r.agron.enet2017,R2.agron.enet2017))
colnames(pred.accuracy.agron.enet2017)[3:4]=c("prediction accuracy","R2")

multivar.agron.prediction.accuracy2017=data.frame(rbind(pred.accuracy.agron.step2017,pred.accuracy.agron.las2017,pred.accuracy.agron.rdg2017,pred.accuracy.agron.enet2017))
# write.csv(multivar.agron.prediction.accuracy2017,file="Multivariate_predaccuracy_agron_2017.csv",row.names = FALSE,quote = FALSE)




### Multivariate model yield prediction 2018 with agron
### Running stepwise regression model with agron
BLUE2018=read.csv("BLUE_2018.csv")
agron_data=BLUE2018[,c(1:2,27:35)]
stepwise.agron.2018 <- lm(GRYLD ~ DTHD+DAYSMT+PH+SN+SPKLNG+SPLN+GRNSPK+TGW,data=agron_data)  ## writing the model with agron
stepwise.agron.2018 <- step(stepwise.agron.2018, direagronion = "both") ## running the stepwise regressing model
formula(stepwise.agron.2018) ## getting the seleagroned variables from the stepwise regression model

### Prediagroned yield from stepwise regression model with agron
step.predyield.agron.2018=NULL 
for(i in 1:11){
        step.prediction.agron.2018=agron_data[agron_data$trial==i,] #get just one trial to prediagron
        step.training.agron.2018=agron_data[agron_data$trial!=i,] #get all other trials to train on
        fit.agron.2018 = lm(formula = formula(stepwise.agron.2018), data =step.training.agron.2018) ## direagronly import variables into the model from stepwise regression 
        step.predvalue.agron.2018=predict(object=fit.agron.2018, newdata=step.prediction.agron.2018) ## Prediagroned value of individual trial
        predicted.agron.trial2018=data.frame(entry=step.prediction.agron.2018$plot,step.predvalue.agron.2018) ## Prediagroned value of individual trial in data frame
        step.predyield.agron.2018=rbind(step.predyield.agron.2018, predicted.agron.trial2018) ## Combinig all trial's prediagroned value into one frame
}
step.predyield.agron.2018 = data.frame(trial=rep(1:11,each=60),step.predyield.agron.2018)
step.predyield.agron.2018 <- merge(step.predyield.agron.2018, agron_data[, c(1:2,11)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extraagroning prediagronion accuracy from stepwise regression model with agron
cor.agron.step2018 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.agron.2018[step.predyield.agron.2018$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.agron.2018, trial_data$GRYLD) #get correlations
        cor.agron.step2018 <- c(cor.agron.step2018, trial_cor) #write results outside of the for loop
}

### Prediagronion from shrinkage models
prediagron.agron.las2018=c()  #Intialize to capture prediagroned values from LASSO model
prediagron.agron.rdg2018=c() #Intialize to capture prediagroned values from Ridge model
prediagron.agron.enet2018=c() #Intialize to capture prediagroned values from ElasticNet model
cor.agron.las2018=c() #Intialize to capture correlatio from LASSO model
cor.agron.rdg2018=c() #Intialize to capture correlatio from Ridge model
cor.agron.enet2018=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.agron.data2018 = agron_data[agron_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.agron.data2018 = agron_data[agron_data$trial!=trial,] #sets up prediagronion populations
        ragronrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.agron.model2018 = train(GRYLD ~ ., data = train.agron.data2018[,-c(1,2)], #Looks like we are fitting a lasso model.
                                      method = "lasso",
                                      trControl = ragronrl,
                                      preProc = c("center", "scale"))
        ridge.agron.model2018 = train(GRYLD ~ ., data = train.agron.data2018[,-c(1,2)], #looks like we are fitting a ridge model.
                                      method = "ridge",
                                      trControl = ragronrl,
                                      preProc = c("center", "scale"))
        enet.agron.model2018 = train(GRYLD ~ ., data = train.agron.data2018[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                     method = "enet",
                                     trControl = ragronrl,
                                     preProc = c("center", "scale"))
        #model coefficients
        coef.agron.lasso2018 = lasso.agron.model2018$finalModel #extraagroning some coefficients from LASSO model
        coef.agron.ridge2018 = ridge.agron.model2018$finalModel #extraagroning some coefficients from Ridge model
        coef.agron.enet2018 = enet.agron.model2018$finalModel #extraagroning some coefficients from ElasticNet model
        #prediagronion
        pred.agron.lasso2018 = predict(lasso.agron.model2018,test.agron.data2018[,-c(1,2)]) #making a prediagron statement.  This looks promising
        pred.agron.ridge2018 = predict(ridge.agron.model2018,test.agron.data2018[,-c(1,2)]) #appears the same with other models
        pred.agron.elasticnet2018 = predict(enet.agron.model2018,test.agron.data2018[,-c(1,2)])
        #extraagroning correlations
        cor.agron.l2018=cor(agron_data[agron_data$trial==trial,11],pred.agron.lasso2018) # Correlaiton between prediagroned value from LASSO regression and BLUE
        cor.agron.las2018=c(cor.agron.las2018,cor.agron.l2018) #writing correlation out to cor.agron.las2018
        cor.agron.r2018=cor(agron_data[agron_data$trial==trial,11],pred.agron.ridge2018) # Correlaiton between prediagroned value from Ridge regression and BLUE
        cor.agron.rdg2018=c(cor.agron.rdg2018,cor.agron.r2018) #writing correlation out to cor.agron.rdg2018
        cor.agron.e2018=cor(agron_data[agron_data$trial==trial,11],pred.agron.elasticnet2018) # Correlaiton between prediagroned value from ElasticNet regression and BLUE
        cor.agron.enet2018=c(cor.agron.enet2018,cor.agron.e2018) #writing correlation out to cor.agron.enet2018
        #extraagroning prediagroned values
        prediagron.agron.las2018 = c(prediagron.agron.las2018,pred.agron.lasso2018) # Combinig all prediagroned values from LASSO
        prediagron.agron.rdg2018 = c(prediagron.agron.rdg2018,pred.agron.ridge2018) # Combinig all prediagroned values from Ridge
        prediagron.agron.enet2018 = c(prediagron.agron.enet2018,pred.agron.elasticnet2018) # Combinig all prediagroned values from ElasticNet
} #end for loop 

#print prediagronion accuracies
cor.agron.step2018 #individual trial correlations for step regression
predaccuracy.agron.step2018=mean(cor.agron.step2018) #average correlation from stepwise model
cor.agron.las2018 #individual trial correlations for LASSO regression
predaccuracy.agron.lasso2018=mean(cor.agron.las2018) #average correlation from LASSO model
cor.agron.rdg2018 #individual trial correlations for Ridge regression
predaccuracy.agron.ridge2018=mean(cor.agron.rdg2018) #average correlation from Ridge model
cor.agron.enet2018 #individual trial correlations for ElasticNet regression
predaccuracy.agron.enet2018=mean(cor.agron.enet2018) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with agron
plot(step.predyield.agron.2018$GRYLD, step.predyield.agron.2018$step.predvalue.agron.2018, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="Stepwise regression with agron 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.step2018 <-lm(step.predyield.agron.2018$step.predvalue.agron.2018 ~ step.predyield.agron.2018$GRYLD) #fit the model for observed values
abline(fit.agron.step2018, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.step2018 = vector('expression',2)
rp.agron.step2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(summary(fit.agron.step2018)$adj.r.squared,dig=2)))[2] # extraagroning R2 value
rp.agron.step2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(predaccuracy.agron.step2018, digits = 2)))[2] # r (prediagronion accuracy) value
legend("bottomright", bty="n",legend=rp.agron.step2018)
R2.agron.step2018=summary(fit.agron.step2018)$adj.r.squared
R2.agron.step2018=round(R2.agron.step2018,digits=2)
r.agron.step2018=round(predaccuracy.agron.step2018,digits=2)
models="stepwise"
traits="agron"
pred.accuracy.agron.step2018=data.frame(cbind(traits,models,r.agron.step2018,R2.agron.step2018))
colnames(pred.accuracy.agron.step2018)[3:4]=c("prediction accuracy","R2")

### Plotting LASSO regression model with agron
plot(agron_data$GRYLD, prediagron.agron.las2018, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="LASSO regression with agron 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.las2018 <-lm(prediagron.agron.las2018 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.las2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.las2018 = vector('expression',2)
rp.agron.las2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.agron.las2018)$adj.r.squared,dig=2)))[2]
rp.agron.las2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.agron.lasso2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.las2018)
R2.agron.las2018=summary(fit.agron.las2018)$adj.r.squared
R2.agron.las2018=round(R2.agron.las2018,digits=2)
r.agron.las2018=round(predaccuracy.agron.lasso2018,digits=2)
models="LASSO"
traits="agron"
pred.accuracy.agron.las2018=data.frame(cbind(traits,models,r.agron.las2018,R2.agron.las2018))
colnames(pred.accuracy.agron.las2018)[3:4]=c("prediction accuracy","R2")

### Plotting ridge regression model with agron
plot(agron_data$GRYLD, prediagron.agron.rdg2018, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="Ridge regression with agron 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.rdg2018 <-lm(prediagron.agron.rdg2018 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.rdg2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.rdg2018 = vector('expression',2)
rp.agron.rdg2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.agron.rdg2018)$adj.r.squared,dig=2)))[2]
rp.agron.rdg2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.agron.ridge2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.rdg2018)
R2.agron.rdg2018=summary(fit.agron.rdg2018)$adj.r.squared
R2.agron.rdg2018=round(R2.agron.rdg2018,digits=2)
r.agron.rdg2018=round(predaccuracy.agron.ridge2018,digits=2)
models="ridge"
traits="agron"
pred.accuracy.agron.rdg2018=data.frame(cbind(traits,models,r.agron.rdg2018,R2.agron.rdg2018))
colnames(pred.accuracy.agron.rdg2018)[3:4]=c("prediction accuracy","R2")

### Plotting enet regression model with agron
plot(agron_data$GRYLD, prediagron.agron.enet2018, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="ElasticNet regression with agron 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.enet2018 <-lm(prediagron.agron.enet2018 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.enet2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.enet2018 = vector('expression',2)
rp.agron.enet2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(summary(fit.agron.enet2018)$adj.r.squared,dig=2)))[2]
rp.agron.enet2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(predaccuracy.agron.enet2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.enet2018)
R2.agron.enet2018=summary(fit.agron.enet2018)$adj.r.squared
R2.agron.enet2018=round(R2.agron.enet2018,digits=2)
r.agron.enet2018=round(predaccuracy.agron.enet2018,digits=2)
models="elasticnet"
traits="agron"
pred.accuracy.agron.enet2018=data.frame(cbind(traits,models,r.agron.enet2018,R2.agron.enet2018))
colnames(pred.accuracy.agron.enet2018)[3:4]=c("prediction accuracy","R2")

multivar.agron.prediction.accuracy2018=data.frame(rbind(pred.accuracy.agron.step2018,pred.accuracy.agron.las2018,pred.accuracy.agron.rdg2018,pred.accuracy.agron.enet2018))
# write.csv(multivar.agron.prediction.accuracy2018,file="Multivariate_predaccuracy_agron_2018.csv",row.names = FALSE,quote = FALSE)




### Multivariate model yield prediction 2019 with agron
### Running stepwise regression model with agron
BLUE2019=read.csv("BLUE_2019.csv")
agron_data=BLUE2019[,c(1:2,29:37)]
stepwise.agron.2019 <- lm(GRYLD ~ DTHD+DAYSMD+PH+SN+SPKLNG+SPLN+GRNSPK+TGW
                          ,data=agron_data)  ## writing the model with agron
stepwise.agron.2019 <- step(stepwise.agron.2019, direagronion = "both") ## running the stepwise regressing model
formula(stepwise.agron.2019) ## getting the seleagroned variables from the stepwise regression model

### Prediagroned yield from stepwise regression model with agron
step.predyield.agron.2019=NULL 
for(i in 1:10){
        step.prediction.agron.2019=agron_data[agron_data$trial==i,] #get just one trial to prediagron
        step.training.agron.2019=agron_data[agron_data$trial!=i,] #get all other trials to train on
        fit.agron.2019 = lm(formula = formula(stepwise.agron.2019), data =step.training.agron.2019) ## direagronly import variables into the model from stepwise regression 
        step.predvalue.agron.2019=predict(object=fit.agron.2019, newdata=step.prediction.agron.2019) ## Prediagroned value of individual trial
        predicted.agron.trial2019=data.frame(entry=step.prediction.agron.2019$plot,step.predvalue.agron.2019) ## Prediagroned value of individual trial in data frame
        step.predyield.agron.2019=rbind(step.predyield.agron.2019, predicted.agron.trial2019) ## Combinig all trial's prediagroned value into one frame
}
step.predyield.agron.2019 = data.frame(trial=rep(1:10,each=60),step.predyield.agron.2019)
step.predyield.agron.2019 <- merge(step.predyield.agron.2019, agron_data[, c(1:2,11)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extraagroning prediagronion accuracy from stepwise regression model with agron
cor.agron.step2019 <- NULL #start dataframe to hold correlation values
for(i in 1:10){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.agron.2019[step.predyield.agron.2019$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.agron.2019, trial_data$GRYLD) #get correlations
        cor.agron.step2019 <- c(cor.agron.step2019, trial_cor) #write results outside of the for loop
}

### Prediagronion from shrinkage models
prediagron.agron.las2019=c()  #Intialize to capture prediagroned values from LASSO model
prediagron.agron.rdg2019=c() #Intialize to capture prediagroned values from Ridge model
prediagron.agron.enet2019=c() #Intialize to capture prediagroned values from ElasticNet model
cor.agron.las2019=c() #Intialize to capture correlatio from LASSO model
cor.agron.rdg2019=c() #Intialize to capture correlatio from Ridge model
cor.agron.enet2019=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:10){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.agron.data2019 = agron_data[agron_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.agron.data2019 = agron_data[agron_data$trial!=trial,] #sets up prediagronion populations
        ragronrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.agron.model2019 = train(GRYLD ~ ., data = train.agron.data2019[,-c(1,2)], #Looks like we are fitting a lasso model.
                                      method = "lasso",
                                      trControl = ragronrl,
                                      preProc = c("center", "scale"))
        ridge.agron.model2019 = train(GRYLD ~ ., data = train.agron.data2019[,-c(1,2)], #looks like we are fitting a ridge model.
                                      method = "ridge",
                                      trControl = ragronrl,
                                      preProc = c("center", "scale"))
        enet.agron.model2019 = train(GRYLD ~ ., data = train.agron.data2019[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                     method = "enet",
                                     trControl = ragronrl,
                                     preProc = c("center", "scale"))
        #model coefficients
        coef.agron.lasso2019 = lasso.agron.model2019$finalModel #extraagroning some coefficients from LASSO model
        coef.agron.ridge2019 = ridge.agron.model2019$finalModel #extraagroning some coefficients from Ridge model
        coef.agron.enet2019 = enet.agron.model2019$finalModel #extraagroning some coefficients from ElasticNet model
        #prediagronion
        pred.agron.lasso2019 = predict(lasso.agron.model2019,test.agron.data2019[,-c(1,2)]) #making a prediagron statement.  This looks promising
        pred.agron.ridge2019 = predict(ridge.agron.model2019,test.agron.data2019[,-c(1,2)]) #appears the same with other models
        pred.agron.elasticnet2019 = predict(enet.agron.model2019,test.agron.data2019[,-c(1,2)])
        #extraagroning correlations
        cor.agron.l2019=cor(agron_data[agron_data$trial==trial,11],pred.agron.lasso2019) # Correlaiton between prediagroned value from LASSO regression and BLUE
        cor.agron.las2019=c(cor.agron.las2019,cor.agron.l2019) #writing correlation out to cor.agron.las2019
        cor.agron.r2019=cor(agron_data[agron_data$trial==trial,11],pred.agron.ridge2019) # Correlaiton between prediagroned value from Ridge regression and BLUE
        cor.agron.rdg2019=c(cor.agron.rdg2019,cor.agron.r2019) #writing correlation out to cor.agron.rdg2019
        cor.agron.e2019=cor(agron_data[agron_data$trial==trial,11],pred.agron.elasticnet2019) # Correlaiton between prediagroned value from ElasticNet regression and BLUE
        cor.agron.enet2019=c(cor.agron.enet2019,cor.agron.e2019) #writing correlation out to cor.agron.enet2019
        #extraagroning prediagroned values
        prediagron.agron.las2019 = c(prediagron.agron.las2019,pred.agron.lasso2019) # Combinig all prediagroned values from LASSO
        prediagron.agron.rdg2019 = c(prediagron.agron.rdg2019,pred.agron.ridge2019) # Combinig all prediagroned values from Ridge
        prediagron.agron.enet2019 = c(prediagron.agron.enet2019,pred.agron.elasticnet2019) # Combinig all prediagroned values from ElasticNet
} #end for loop 

#print prediagronion accuracies
cor.agron.step2019 #individual trial correlations for step regression
predaccuracy.agron.step2019=mean(cor.agron.step2019) #average correlation from stepwise model
cor.agron.las2019 #individual trial correlations for LASSO regression
predaccuracy.agron.lasso2019=mean(cor.agron.las2019) #average correlation from LASSO model
cor.agron.rdg2019 #individual trial correlations for Ridge regression
predaccuracy.agron.ridge2019=mean(cor.agron.rdg2019) #average correlation from Ridge model
cor.agron.enet2019 #individual trial correlations for ElasticNet regression
predaccuracy.agron.enet2019=mean(cor.agron.enet2019) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with agron
plot(step.predyield.agron.2019$GRYLD, step.predyield.agron.2019$step.predvalue.agron.2019, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="Stepwise regression with agron 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.step2019 <-lm(step.predyield.agron.2019$step.predvalue.agron.2019 ~ step.predyield.agron.2019$GRYLD) #fit the model for observed values
abline(fit.agron.step2019, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.step2019 = vector('expression',2)
rp.agron.step2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(summary(fit.agron.step2019)$adj.r.squared,dig=2)))[2] # extraagroning R2 value
rp.agron.step2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(predaccuracy.agron.step2019, digits = 2)))[2] # r (prediagronion accuracy) value
legend("bottomright", bty="n",legend=rp.agron.step2019)
R2.agron.step2019=summary(fit.agron.step2019)$adj.r.squared
R2.agron.step2019=round(R2.agron.step2019,digits=2)
r.agron.step2019=round(predaccuracy.agron.step2019,digits=2)
models="stepwise"
traits="agron"
pred.accuracy.agron.step2019=data.frame(cbind(traits,models,r.agron.step2019,R2.agron.step2019))
colnames(pred.accuracy.agron.step2019)[3:4]=c("prediction accuracy","R2")

### Plotting LASSO regression model with agron
plot(agron_data$GRYLD, prediagron.agron.las2019, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="LASSO regression with agron 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.las2019 <-lm(prediagron.agron.las2019 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.las2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.las2019 = vector('expression',2)
rp.agron.las2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.agron.las2019)$adj.r.squared,dig=2)))[2]
rp.agron.las2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.agron.lasso2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.las2019)
R2.agron.las2019=summary(fit.agron.las2019)$adj.r.squared
R2.agron.las2019=round(R2.agron.las2019,digits=2)
r.agron.las2019=round(predaccuracy.agron.lasso2019,digits=2)
models="LASSO"
traits="agron"
pred.accuracy.agron.las2019=data.frame(cbind(traits,models,r.agron.las2019,R2.agron.las2019))
colnames(pred.accuracy.agron.las2019)[3:4]=c("prediction accuracy","R2")

### Plotting ridge regression model with agron
plot(agron_data$GRYLD, prediagron.agron.rdg2019, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="Ridge regression with agron 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.rdg2019 <-lm(prediagron.agron.rdg2019 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.rdg2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.rdg2019 = vector('expression',2)
rp.agron.rdg2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.agron.rdg2019)$adj.r.squared,dig=2)))[2]
rp.agron.rdg2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.agron.ridge2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.rdg2019)
R2.agron.rdg2019=summary(fit.agron.rdg2019)$adj.r.squared
R2.agron.rdg2019=round(R2.agron.rdg2019,digits=2)
r.agron.rdg2019=round(predaccuracy.agron.ridge2019,digits=2)
models="ridge"
traits="agron"
pred.accuracy.agron.rdg2019=data.frame(cbind(traits,models,r.agron.rdg2019,R2.agron.rdg2019))
colnames(pred.accuracy.agron.rdg2019)[3:4]=c("prediction accuracy","R2")

### Plotting enet regression model with agron
plot(agron_data$GRYLD, prediagron.agron.enet2019, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="ElasticNet regression with agron 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.enet2019 <-lm(prediagron.agron.enet2019 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.enet2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.enet2019 = vector('expression',2)
rp.agron.enet2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(summary(fit.agron.enet2019)$adj.r.squared,dig=2)))[2]
rp.agron.enet2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(predaccuracy.agron.enet2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.enet2019)
R2.agron.enet2019=summary(fit.agron.enet2019)$adj.r.squared
R2.agron.enet2019=round(R2.agron.enet2019,digits=2)
r.agron.enet2019=round(predaccuracy.agron.enet2019,digits=2)
models="elasticnet"
traits="agron"
pred.accuracy.agron.enet2019=data.frame(cbind(traits,models,r.agron.enet2019,R2.agron.enet2019))
colnames(pred.accuracy.agron.enet2019)[3:4]=c("prediction accuracy","R2")

multivar.agron.prediction.accuracy2019=data.frame(rbind(pred.accuracy.agron.step2019,pred.accuracy.agron.las2019,pred.accuracy.agron.rdg2019,pred.accuracy.agron.enet2019))
write.csv(multivar.agron.prediction.accuracy2019,file="Multivariate_predaccuracy_agron_2019.csv",row.names = FALSE,quote = FALSE)

multivar.yield.prediction.2019=rbind(multivar.ct.prediction.accuracy2019,multivar.ndvi.prediction.accuracy2019,multivar.ndvi_ct.prediction.accuracy2019,multivar.agron.prediction.accuracy2019,multivar.alltraits.prediction.accuracy2019)
# write.csv(multivar.yield.prediction.2019,file="Multivariate Yield Prediction 2019.csv",row.names=FALSE,quote=FALSE)



### Multivariate model yield prediction 2020 with agron
### Running stepwise regression model with agron
BLUE2020=read.csv("BLUE_2020.csv")
agron_data=BLUE2020[,c(1:2,33:44)]
stepwise.agron.2020 <- lm(GRYLD ~ GrndCov_20200112+GrndCov_20200206+DLA_Feb26+DLA_Mar09+DTHD+DAYSMT+PH+SN+SPLN+GRNSPK+TGW
                          ,data=agron_data)  ## writing the model with agron
stepwise.agron.2020 <- step(stepwise.agron.2020, direagronion = "both") ## running the stepwise regressing model
formula(stepwise.agron.2020) ## getting the seleagroned variables from the stepwise regression model

### Prediagroned yield from stepwise regression model with agron
step.predyield.agron.2020=NULL 
for(i in 1:11){
        step.prediction.agron.2020=agron_data[agron_data$trial==i,] #get just one trial to prediagron
        step.training.agron.2020=agron_data[agron_data$trial!=i,] #get all other trials to train on
        fit.agron.2020 = lm(formula = formula(stepwise.agron.2020), data =step.training.agron.2020) ## direagronly import variables into the model from stepwise regression 
        step.predvalue.agron.2020=predict(object=fit.agron.2020, newdata=step.prediction.agron.2020) ## Prediagroned value of individual trial
        predicted.agron.trial2020=data.frame(entry=step.prediction.agron.2020$plot,step.predvalue.agron.2020) ## Prediagroned value of individual trial in data frame
        step.predyield.agron.2020=rbind(step.predyield.agron.2020, predicted.agron.trial2020) ## Combinig all trial's prediagroned value into one frame
}
step.predyield.agron.2020 = data.frame(trial=rep(1:11,each=60),step.predyield.agron.2020)
step.predyield.agron.2020 <- merge(step.predyield.agron.2020, agron_data[, c(1:2,14)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extraagroning prediagronion accuracy from stepwise regression model with agron
cor.agron.step2020 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.agron.2020[step.predyield.agron.2020$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.agron.2020, trial_data$GRYLD) #get correlations
        cor.agron.step2020 <- c(cor.agron.step2020, trial_cor) #write results outside of the for loop
}

### Prediagronion from shrinkage models
prediagron.agron.las2020=c()  #Intialize to capture prediagroned values from LASSO model
prediagron.agron.rdg2020=c() #Intialize to capture prediagroned values from Ridge model
prediagron.agron.enet2020=c() #Intialize to capture prediagroned values from ElasticNet model
cor.agron.las2020=c() #Intialize to capture correlatio from LASSO model
cor.agron.rdg2020=c() #Intialize to capture correlatio from Ridge model
cor.agron.enet2020=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.agron.data2020 = agron_data[agron_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.agron.data2020 = agron_data[agron_data$trial!=trial,] #sets up prediagronion populations
        ragronrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.agron.model2020 = train(GRYLD ~ ., data = train.agron.data2020[,-c(1,2)], #Looks like we are fitting a lasso model.
                                      method = "lasso",
                                      trControl = ragronrl,
                                      preProc = c("center", "scale"))
        ridge.agron.model2020 = train(GRYLD ~ ., data = train.agron.data2020[,-c(1,2)], #looks like we are fitting a ridge model.
                                      method = "ridge",
                                      trControl = ragronrl,
                                      preProc = c("center", "scale"))
        enet.agron.model2020 = train(GRYLD ~ ., data = train.agron.data2020[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                     method = "enet",
                                     trControl = ragronrl,
                                     preProc = c("center", "scale"))
        #model coefficients
        coef.agron.lasso2020 = lasso.agron.model2020$finalModel #extraagroning some coefficients from LASSO model
        coef.agron.ridge2020 = ridge.agron.model2020$finalModel #extraagroning some coefficients from Ridge model
        coef.agron.enet2020 = enet.agron.model2020$finalModel #extraagroning some coefficients from ElasticNet model
        #prediagronion
        pred.agron.lasso2020 = predict(lasso.agron.model2020,test.agron.data2020[,-c(1,2)]) #making a prediagron statement.  This looks promising
        pred.agron.ridge2020 = predict(ridge.agron.model2020,test.agron.data2020[,-c(1,2)]) #appears the same with other models
        pred.agron.elasticnet2020 = predict(enet.agron.model2020,test.agron.data2020[,-c(1,2)])
        #extraagroning correlations
        cor.agron.l2020=cor(agron_data[agron_data$trial==trial,14],pred.agron.lasso2020) # Correlaiton between prediagroned value from LASSO regression and BLUE
        cor.agron.las2020=c(cor.agron.las2020,cor.agron.l2020) #writing correlation out to cor.agron.las2020
        cor.agron.r2020=cor(agron_data[agron_data$trial==trial,14],pred.agron.ridge2020) # Correlaiton between prediagroned value from Ridge regression and BLUE
        cor.agron.rdg2020=c(cor.agron.rdg2020,cor.agron.r2020) #writing correlation out to cor.agron.rdg2020
        cor.agron.e2020=cor(agron_data[agron_data$trial==trial,14],pred.agron.elasticnet2020) # Correlaiton between prediagroned value from ElasticNet regression and BLUE
        cor.agron.enet2020=c(cor.agron.enet2020,cor.agron.e2020) #writing correlation out to cor.agron.enet2020
        #extraagroning prediagroned values
        prediagron.agron.las2020 = c(prediagron.agron.las2020,pred.agron.lasso2020) # Combinig all prediagroned values from LASSO
        prediagron.agron.rdg2020 = c(prediagron.agron.rdg2020,pred.agron.ridge2020) # Combinig all prediagroned values from Ridge
        prediagron.agron.enet2020 = c(prediagron.agron.enet2020,pred.agron.elasticnet2020) # Combinig all prediagroned values from ElasticNet
} #end for loop 

#print prediagronion accuracies
cor.agron.step2020 #individual trial correlations for step regression
predaccuracy.agron.step2020=mean(cor.agron.step2020) #average correlation from stepwise model
cor.agron.las2020 #individual trial correlations for LASSO regression
predaccuracy.agron.lasso2020=mean(cor.agron.las2020) #average correlation from LASSO model
cor.agron.rdg2020 #individual trial correlations for Ridge regression
predaccuracy.agron.ridge2020=mean(cor.agron.rdg2020) #average correlation from Ridge model
cor.agron.enet2020 #individual trial correlations for ElasticNet regression
predaccuracy.agron.enet2020=mean(cor.agron.enet2020) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with agron
plot(step.predyield.agron.2020$GRYLD, step.predyield.agron.2020$step.predvalue.agron.2020, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="Stepwise regression with agron 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.step2020 <-lm(step.predyield.agron.2020$step.predvalue.agron.2020 ~ step.predyield.agron.2020$GRYLD) #fit the model for observed values
abline(fit.agron.step2020, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.step2020 = vector('expression',2)
rp.agron.step2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(summary(fit.agron.step2020)$adj.r.squared,dig=2)))[2] # extraagroning R2 value
rp.agron.step2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(predaccuracy.agron.step2020, digits = 2)))[2] # r (prediagronion accuracy) value
legend("bottomright", bty="n",legend=rp.agron.step2020)
R2.agron.step2020=summary(fit.agron.step2020)$adj.r.squared
R2.agron.step2020=round(R2.agron.step2020,digits=2)
r.agron.step2020=round(predaccuracy.agron.step2020,digits=2)
models="stepwise"
traits="agron"
pred.accuracy.agron.step2020=data.frame(cbind(traits,models,r.agron.step2020,R2.agron.step2020))
colnames(pred.accuracy.agron.step2020)[3:4]=c("prediction accuracy","R2")

### Plotting LASSO regression model with agron
plot(agron_data$GRYLD, prediagron.agron.las2020, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="LASSO regression with agron 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.las2020 <-lm(prediagron.agron.las2020 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.las2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.las2020 = vector('expression',2)
rp.agron.las2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.agron.las2020)$adj.r.squared,dig=2)))[2]
rp.agron.las2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.agron.lasso2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.las2020)
R2.agron.las2020=summary(fit.agron.las2020)$adj.r.squared
R2.agron.las2020=round(R2.agron.las2020,digits=2)
r.agron.las2020=round(predaccuracy.agron.lasso2020,digits=2)
models="LASSO"
traits="agron"
pred.accuracy.agron.las2020=data.frame(cbind(traits,models,r.agron.las2020,R2.agron.las2020))
colnames(pred.accuracy.agron.las2020)[3:4]=c("prediction accuracy","R2")

### Plotting ridge regression model with agron
plot(agron_data$GRYLD, prediagron.agron.rdg2020, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="Ridge regression with agron 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.rdg2020 <-lm(prediagron.agron.rdg2020 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.rdg2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.rdg2020 = vector('expression',2)
rp.agron.rdg2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                 list(MYVALUE = format(summary(fit.agron.rdg2020)$adj.r.squared,dig=2)))[2]
rp.agron.rdg2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                 list(MYOTHERVALUE = format(predaccuracy.agron.ridge2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.rdg2020)
R2.agron.rdg2020=summary(fit.agron.rdg2020)$adj.r.squared
R2.agron.rdg2020=round(R2.agron.rdg2020,digits=2)
r.agron.rdg2020=round(predaccuracy.agron.ridge2020,digits=2)
models="ridge"
traits="agron"
pred.accuracy.agron.rdg2020=data.frame(cbind(traits,models,r.agron.rdg2020,R2.agron.rdg2020))
colnames(pred.accuracy.agron.rdg2020)[3:4]=c("prediction accuracy","R2")

### Plotting enet regression model with agron
plot(agron_data$GRYLD, prediagron.agron.enet2020, xlab = "Observed Yield",ylab = "Prediagroned Yield", main="ElasticNet regression with agron 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.agron.enet2020 <-lm(prediagron.agron.enet2020 ~ agron_data$GRYLD) #fit the model for observed values
abline(fit.agron.enet2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.agron.enet2020 = vector('expression',2)
rp.agron.enet2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                  list(MYVALUE = format(summary(fit.agron.enet2020)$adj.r.squared,dig=2)))[2]
rp.agron.enet2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                  list(MYOTHERVALUE = format(predaccuracy.agron.enet2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.agron.enet2020)
R2.agron.enet2020=summary(fit.agron.enet2020)$adj.r.squared
R2.agron.enet2020=round(R2.agron.enet2020,digits=2)
r.agron.enet2020=round(predaccuracy.agron.enet2020,digits=2)
models="elasticnet"
traits="agron"
pred.accuracy.agron.enet2020=data.frame(cbind(traits,models,r.agron.enet2020,R2.agron.enet2020))
colnames(pred.accuracy.agron.enet2020)[3:4]=c("prediction accuracy","R2")

multivar.agron.prediction.accuracy2020=data.frame(rbind(pred.accuracy.agron.step2020,pred.accuracy.agron.las2020,pred.accuracy.agron.rdg2020,pred.accuracy.agron.enet2020))
# write.csv(multivar.agron.prediction.accuracy2020,file="Multivariate_predaccuracy_agron_2020.csv",row.names = FALSE,quote = FALSE)

multivar.agron.allyear=cbind(multivar.agron.prediction.accuracy2016,multivar.agron.prediction.accuracy2017[,3:4],multivar.agron.prediction.accuracy2018[,3:4],multivar.agron.prediction.accuracy2019[,3:4],multivar.agron.prediction.accuracy2020[,3:4])
write.csv(multivar.agron.allyear,file="multivar.yieldpred.agron.allYears.csv",row.names=FALSE,quote=FALSE)
