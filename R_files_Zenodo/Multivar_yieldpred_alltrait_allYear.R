### Multivariate model yield prediction 2016 with all traits together
### Running stepwise regression model with all traits
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
stepwise.alltraits.2016 <- lm(GRYLD ~ CT_20160123+CT_20160204+CT_20160212+CT_20160223+CT_20160228+CT_20160302+CT_20160309+CT_20160315+
                                      NDVI_20160121+NDVI_20160130+NDVI_20160203+NDVI_20160207+NDVI_20160223+NDVI_20160228+NDVI_20160303+NDVI_20160310+NDVI_20160315+
                                      DTHD+DAYSMT+PH+SN+SPKLNG+SPLN+GRNSPK+TGW,data=BLUE2016)  ## writing the model with all traits together
stepwise.alltraits.2016 <- step(stepwise.alltraits.2016, direction = "both") ## running the stepwise regressing model
formula(stepwise.alltraits.2016) ## getting the selected variables from the stepwise regression model

### Predicted yield from stepwise regression model with all traits together
step.predyield.alltraits.2016=NULL 
for(i in 1:10){
        step.prediction.alltraits.2016=BLUE2016[BLUE2016$trial==i,] #get just one trial to predict
        step.training.alltraits.2016=BLUE2016[BLUE2016$trial!=i,] #get all other trials to train on
        fit.alltraits.2016 = lm(formula = formula(stepwise.alltraits.2016), data =step.training.alltraits.2016) ## directly import variables into the model from stepwise regression 
        step.predvalue.alltraits.2016=predict(object=fit.alltraits.2016, newdata=step.prediction.alltraits.2016) ## Predicted value of individual trial
        predicted.alltraits.trial2016=data.frame(entry=step.prediction.alltraits.2016$plot,step.predvalue.alltraits.2016) ## Predicted value of individual trial in data frame
        step.predyield.alltraits.2016=rbind(step.predyield.alltraits.2016, predicted.alltraits.trial2016) ## Combinig all trial's predicted value into one frame
}
step.predyield.alltraits.2016 = data.frame(trial=rep(1:10,each=60),step.predyield.alltraits.2016)
step.predyield.alltraits.2016 <- merge(step.predyield.alltraits.2016, BLUE2016[, c(1:2,28)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extracting prediction accuracy from stepwise regression model with all traits together
cor.alltraits.step2016 <- NULL #start dataframe to hold correlation values
for(i in 1:10){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.alltraits.2016[step.predyield.alltraits.2016$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.alltraits.2016, trial_data$GRYLD) #get correlations
        cor.alltraits.step2016 <- c(cor.alltraits.step2016, trial_cor) #write results outside of the for loop
}

### Prediction from shrinkage models
predict.alltraits.las2016=c()  #Intialize to capture predicted values from LASSO model
predict.alltraits.rdg2016=c() #Intialize to capture predicted values from Ridge model
predict.alltraits.enet2016=c() #Intialize to capture predicted values from ElasticNet model
cor.alltraits.las2016=c() #Intialize to capture correlatio from LASSO model
cor.alltraits.rdg2016=c() #Intialize to capture correlatio from Ridge model
cor.alltraits.enet2016=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:10){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.alltraits.data2016 = BLUE2016[BLUE2016$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.alltraits.data2016 = BLUE2016[BLUE2016$trial!=trial,] #sets up prediction populations
        rctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.alltraits.model2016 = train(GRYLD ~ ., data = train.alltraits.data2016[,-c(1,2)], #Looks like we are fitting a lasso model.
                                          method = "lasso",
                                          trControl = rctrl,
                                          preProc = c("center", "scale"))
        ridge.alltraits.model2016 = train(GRYLD ~ ., data = train.alltraits.data2016[,-c(1,2)], #looks like we are fitting a ridge model.
                                          method = "ridge",
                                          trControl = rctrl,
                                          preProc = c("center", "scale"))
        enet.alltraits.model2016 = train(GRYLD ~ ., data = train.alltraits.data2016[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                         method = "enet",
                                         trControl = rctrl,
                                         preProc = c("center", "scale"))
        #model coefficients
        coef.alltraits.lasso2016 = lasso.alltraits.model2016$finalModel #extracting some coefficients from LASSO model
        coef.alltraits.ridge2016 = ridge.alltraits.model2016$finalModel #extracting some coefficients from Ridge model
        coef.alltraits.enet2016 = enet.alltraits.model2016$finalModel #extracting some coefficients from ElasticNet model
        #prediction
        pred.alltraits.lasso2016 = predict(lasso.alltraits.model2016,test.alltraits.data2016[,-c(1,2)]) #making a predict statement.  This looks promising
        pred.alltraits.ridge2016 = predict(ridge.alltraits.model2016,test.alltraits.data2016[,-c(1,2)]) #appears the same with other models
        pred.alltraits.elasticnet2016 = predict(enet.alltraits.model2016,test.alltraits.data2016[,-c(1,2)])
        #extracting correlations
        cor.alltraits.l2016=cor(BLUE2016[BLUE2016$trial==trial,28],pred.alltraits.lasso2016) # Correlaiton between predicted value from LASSO regression and BLUE
        cor.alltraits.las2016=c(cor.alltraits.las2016,cor.alltraits.l2016) #writing correlation out to cor.alltraits.las2016
        cor.alltraits.r2016=cor(BLUE2016[BLUE2016$trial==trial,28],pred.alltraits.ridge2016) # Correlaiton between predicted value from Ridge regression and BLUE
        cor.alltraits.rdg2016=c(cor.alltraits.rdg2016,cor.alltraits.r2016) #writing correlation out to cor.alltraits.rdg2016
        cor.alltraits.e2016=cor(BLUE2016[BLUE2016$trial==trial,28],pred.alltraits.elasticnet2016) # Correlaiton between predicted value from ElasticNet regression and BLUE
        cor.alltraits.enet2016=c(cor.alltraits.enet2016,cor.alltraits.e2016) #writing correlation out to cor.alltraits.enet2016
        #extracting predicted values
        predict.alltraits.las2016 = c(predict.alltraits.las2016,pred.alltraits.lasso2016) # Combinig all predicted values from LASSO
        predict.alltraits.rdg2016 = c(predict.alltraits.rdg2016,pred.alltraits.ridge2016) # Combinig all predicted values from Ridge
        predict.alltraits.enet2016 = c(predict.alltraits.enet2016,pred.alltraits.elasticnet2016) # Combinig all predicted values from ElasticNet
} #end for loop 

#print prediction accuracies
cor.alltraits.step2016 #individual trial correlations for step regression
predaccuracy.alltraits.step2016=mean(cor.alltraits.step2016) #average correlation from stepwise model
cor.alltraits.las2016 #individual trial correlations for LASSO regression
predaccuracy.alltraits.lasso2016=mean(cor.alltraits.las2016) #average correlation from LASSO model
cor.alltraits.rdg2016 #individual trial correlations for Ridge regression
predaccuracy.alltraits.ridge2016=mean(cor.alltraits.rdg2016) #average correlation from Ridge model
cor.alltraits.enet2016 #individual trial correlations for ElasticNet regression
predaccuracy.alltraits.enet2016=mean(cor.alltraits.enet2016) #average correlation from ElasticNet model

### Plotting different models
### Plotting stepwise regression model with all traits together
plot(step.predyield.alltraits.2016$GRYLD, step.predyield.alltraits.2016$step.predvalue.alltraits.2016, xlab = "Observed Yield",ylab = "Predicted Yield", main="Stepwise regression with all traits 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.step2016 <-lm(step.predyield.alltraits.2016$step.predvalue.alltraits.2016 ~ step.predyield.alltraits.2016$GRYLD) #fit the model for observed values
abline(fit.alltraits.step2016, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.step2016 = vector('expression',2)
rp.alltraits.step2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                      list(MYVALUE = format(summary(fit.alltraits.step2016)$adj.r.squared,dig=2)))[2] # extracting R2 value
rp.alltraits.step2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                      list(MYOTHERVALUE = format(predaccuracy.alltraits.step2016, digits = 2)))[2] # r (prediction accuracy) value
legend("bottomright", bty="n",legend=rp.alltraits.step2016)
R2.alltraits.step2016=summary(fit.alltraits.step2016)$adj.r.squared
R2.alltraits.step2016=round(R2.alltraits.step2016,digits=2)
r.alltraits.step2016=round(predaccuracy.alltraits.step2016,digits=2)
models="stepwise"
traits="alltraits"
pred.accuracy.alltraits.step2016=data.frame(cbind(traits,models,r.alltraits.step2016,R2.alltraits.step2016))
colnames(pred.accuracy.alltraits.step2016)[3:4]=c("r","R2")
# dev.off()
### Plotting LASSO regression model with all traits together
plot(BLUE2016$GRYLD, predict.alltraits.las2016, xlab = "Observed Yield",ylab = "Predicted Yield", main="LASSO regression with all traits 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.las2016 <-lm(predict.alltraits.las2016 ~ BLUE2016$GRYLD) #fit the model for observed values
abline(fit.alltraits.las2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.las2016 = vector('expression',2)
rp.alltraits.las2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(summary(fit.alltraits.las2016)$adj.r.squared,dig=2)))[2]
rp.alltraits.las2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                     list(MYOTHERVALUE = format(predaccuracy.alltraits.lasso2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.las2016)
R2.alltraits.las2016=summary(fit.alltraits.las2016)$adj.r.squared
R2.alltraits.las2016=round(R2.alltraits.las2016,digits=2)
r.alltraits.las2016=round(predaccuracy.alltraits.lasso2016,digits=2)
models="LASSO"
traits="alltraits"
pred.accuracy.alltraits.las2016=data.frame(cbind(traits,models,r.alltraits.las2016,R2.alltraits.las2016))
colnames(pred.accuracy.alltraits.las2016)[3:4]=c("r","R2")

### Plotting ridge regression model with all traits together
plot(BLUE2016$GRYLD, predict.alltraits.rdg2016, xlab = "Observed Yield",ylab = "Predicted Yield", main="Ridge regression with all traits 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.rdg2016 <-lm(predict.alltraits.rdg2016 ~ BLUE2016$GRYLD) #fit the model for observed values
abline(fit.alltraits.rdg2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.rdg2016 = vector('expression',2)
rp.alltraits.rdg2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(summary(fit.alltraits.rdg2016)$adj.r.squared,dig=2)))[2]
rp.alltraits.rdg2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                     list(MYOTHERVALUE = format(predaccuracy.alltraits.ridge2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.rdg2016)
R2.alltraits.rdg2016=summary(fit.alltraits.rdg2016)$adj.r.squared
R2.alltraits.rdg2016=round(R2.alltraits.rdg2016,digits=2)
r.alltraits.rdg2016=round(predaccuracy.alltraits.ridge2016,digits=2)
models="ridge"
traits="alltraits"
pred.accuracy.alltraits.rdg2016=data.frame(cbind(traits,models,r.alltraits.rdg2016,R2.alltraits.rdg2016))
colnames(pred.accuracy.alltraits.rdg2016)[3:4]=c("r","R2")

### Plotting enet regression model with all traits together
plot(BLUE2016$GRYLD, predict.alltraits.enet2016, xlab = "Observed Yield",ylab = "Predicted Yield", main="ElasticNet regression with all traits 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.enet2016 <-lm(predict.alltraits.enet2016 ~ BLUE2016$GRYLD) #fit the model for observed values
abline(fit.alltraits.enet2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.enet2016 = vector('expression',2)
rp.alltraits.enet2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                      list(MYVALUE = format(summary(fit.alltraits.enet2016)$adj.r.squared,dig=2)))[2]
rp.alltraits.enet2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                      list(MYOTHERVALUE = format(predaccuracy.alltraits.enet2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.enet2016)
R2.alltraits.enet2016=summary(fit.alltraits.enet2016)$adj.r.squared
R2.alltraits.enet2016=round(R2.alltraits.enet2016,digits=2)
r.alltraits.enet2016=round(predaccuracy.alltraits.enet2016,digits=2)
models="elasticnet"
traits="alltraits"
pred.accuracy.alltraits.enet2016=data.frame(cbind(traits,models,r.alltraits.enet2016,R2.alltraits.enet2016))
colnames(pred.accuracy.alltraits.enet2016)[3:4]=c("r","R2")

multivar.alltraits.prediction.accuracy2016=data.frame(rbind(pred.accuracy.alltraits.step2016,pred.accuracy.alltraits.las2016,pred.accuracy.alltraits.rdg2016,pred.accuracy.alltraits.enet2016))
# write.csv(multivar.alltraits.prediction.accuracy2016,file="Multivariate_predaccuracy_alltraits_2016.csv",row.names = FALSE,quote = FALSE)



### Multivariate model yield prediction 2017 with all traits together
### Running stepwise regression model with all traits
BLUE2017=read.csv("BLUE_2017.csv")
stepwise.alltraits.2017 <- lm(GRYLD ~ CT_20170104+CT_20170109+CT_20170114+CT_20170120+CT_20170125+CT_20170131+CT_20170205+CT_20170210+CT_20170215+CT_20170221+CT_20170225+CT_20170302+CT_20170307+CT_20170313+
                                      NDVI_20170103+NDVI_20170108+NDVI_20170114+NDVI_20170120+NDVI_20170125+NDVI_20170131+NDVI_20170205+NDVI_20170210+NDVI_20170215+NDVI_20170220+NDVI_20170225+NDVI_20170302+NDVI_20170307+NDVI_20170313+
                                      DTHD+DAYSMT+PH+SN+SPKLNG+SPLN+GRNSPK+TGW,data=BLUE2017)  ## writing the model with all traits together
stepwise.alltraits.2017 <- step(stepwise.alltraits.2017, direction = "both") ## running the stepwise regressing model
formula(stepwise.alltraits.2017) ## getting the selected variables from the stepwise regression model

### Predicted yield from stepwise regression model with all traits together
step.predyield.alltraits.2017=NULL 
for(i in 1:11){
        step.prediction.alltraits.2017=BLUE2017[BLUE2017$trial==i,] #get just one trial to predict
        step.training.alltraits.2017=BLUE2017[BLUE2017$trial!=i,] #get all other trials to train on
        fit.alltraits.2017 = lm(formula = formula(stepwise.alltraits.2017), data =step.training.alltraits.2017) ## directly import variables into the model from stepwise regression 
        step.predvalue.alltraits.2017=predict(object=fit.alltraits.2017, newdata=step.prediction.alltraits.2017) ## Predicted value of individual trial
        predicted.alltraits.trial2017=data.frame(entry=step.prediction.alltraits.2017$plot,step.predvalue.alltraits.2017) ## Predicted value of individual trial in data frame
        step.predyield.alltraits.2017=rbind(step.predyield.alltraits.2017, predicted.alltraits.trial2017) ## Combinig all trial's predicted value into one frame
}
step.predyield.alltraits.2017 = data.frame(trial=rep(1:11,each=60),step.predyield.alltraits.2017)
step.predyield.alltraits.2017 <- merge(step.predyield.alltraits.2017, BLUE2017[, c(1:2,39)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extracting prediction accuracy from stepwise regression model with all traits together
cor.alltraits.step2017 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 11 trials for 10 fold cross validation
        trial_data <- step.predyield.alltraits.2017[step.predyield.alltraits.2017$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.alltraits.2017, trial_data$GRYLD) #get correlations
        cor.alltraits.step2017 <- c(cor.alltraits.step2017, trial_cor) #write results outside of the for loop
}

### Prediction from shrinkage models
predict.alltraits.las2017=c()  #Intialize to capture predicted values from LASSO model
predict.alltraits.rdg2017=c() #Intialize to capture predicted values from Ridge model
predict.alltraits.enet2017=c() #Intialize to capture predicted values from ElasticNet model
cor.alltraits.las2017=c() #Intialize to capture correlatio from LASSO model
cor.alltraits.rdg2017=c() #Intialize to capture correlatio from Ridge model
cor.alltraits.enet2017=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.alltraits.data2017 = BLUE2017[BLUE2017$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.alltraits.data2017 = BLUE2017[BLUE2017$trial!=trial,] #sets up prediction populations
        rctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.alltraits.model2017 = train(GRYLD ~ ., data = train.alltraits.data2017[,-c(1,2)], #Looks like we are fitting a lasso model.
                                          method = "lasso",
                                          trControl = rctrl,
                                          preProc = c("center", "scale"))
        ridge.alltraits.model2017 = train(GRYLD ~ ., data = train.alltraits.data2017[,-c(1,2)], #looks like we are fitting a ridge model.
                                          method = "ridge",
                                          trControl = rctrl,
                                          preProc = c("center", "scale"))
        enet.alltraits.model2017 = train(GRYLD ~ ., data = train.alltraits.data2017[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                         method = "enet",
                                         trControl = rctrl,
                                         preProc = c("center", "scale"))
        #model coefficients
        coef.alltraits.lasso2017 = lasso.alltraits.model2017$finalModel #extracting some coefficients from LASSO model
        coef.alltraits.ridge2017 = ridge.alltraits.model2017$finalModel #extracting some coefficients from Ridge model
        coef.alltraits.enet2017 = enet.alltraits.model2017$finalModel #extracting some coefficients from ElasticNet model
        #prediction
        pred.alltraits.lasso2017 = predict(lasso.alltraits.model2017,test.alltraits.data2017[,-c(1,2)]) #making a predict statement.  This looks promising
        pred.alltraits.ridge2017 = predict(ridge.alltraits.model2017,test.alltraits.data2017[,-c(1,2)]) #appears the same with other models
        pred.alltraits.elasticnet2017 = predict(enet.alltraits.model2017,test.alltraits.data2017[,-c(1,2)])
        #extracting correlations
        cor.alltraits.l2017=cor(BLUE2017[BLUE2017$trial==trial,39],pred.alltraits.lasso2017) # Correlaiton between predicted value from LASSO regression and BLUE
        cor.alltraits.las2017=c(cor.alltraits.las2017,cor.alltraits.l2017) #writing correlation out to cor.alltraits.las2017
        cor.alltraits.r2017=cor(BLUE2017[BLUE2017$trial==trial,39],pred.alltraits.ridge2017) # Correlaiton between predicted value from Ridge regression and BLUE
        cor.alltraits.rdg2017=c(cor.alltraits.rdg2017,cor.alltraits.r2017) #writing correlation out to cor.alltraits.rdg2017
        cor.alltraits.e2017=cor(BLUE2017[BLUE2017$trial==trial,39],pred.alltraits.elasticnet2017) # Correlaiton between predicted value from ElasticNet regression and BLUE
        cor.alltraits.enet2017=c(cor.alltraits.enet2017,cor.alltraits.e2017) #writing correlation out to cor.alltraits.enet2017
        #extracting predicted values
        predict.alltraits.las2017 = c(predict.alltraits.las2017,pred.alltraits.lasso2017) # Combinig all predicted values from LASSO
        predict.alltraits.rdg2017 = c(predict.alltraits.rdg2017,pred.alltraits.ridge2017) # Combinig all predicted values from Ridge
        predict.alltraits.enet2017 = c(predict.alltraits.enet2017,pred.alltraits.elasticnet2017) # Combinig all predicted values from ElasticNet
} #end for loop 

#print prediction accuracies
cor.alltraits.step2017 #individual trial correlations for step regression
predaccuracy.alltraits.step2017=mean(cor.alltraits.step2017) #average correlation from stepwise model
cor.alltraits.las2017 #individual trial correlations for LASSO regression
predaccuracy.alltraits.lasso2017=mean(cor.alltraits.las2017) #average correlation from LASSO model
cor.alltraits.rdg2017 #individual trial correlations for Ridge regression
predaccuracy.alltraits.ridge2017=mean(cor.alltraits.rdg2017) #average correlation from Ridge model
cor.alltraits.enet2017 #individual trial correlations for ElasticNet regression
predaccuracy.alltraits.enet2017=mean(cor.alltraits.enet2017) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with all traits together
plot(step.predyield.alltraits.2017$GRYLD, step.predyield.alltraits.2017$step.predvalue.alltraits.2017, xlab = "Observed Yield",ylab = "Predicted Yield", main="Stepwise regression with all traits 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.step2017 <-lm(step.predyield.alltraits.2017$step.predvalue.alltraits.2017 ~ step.predyield.alltraits.2017$GRYLD) #fit the model for observed values
abline(fit.alltraits.step2017, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.step2017 = vector('expression',2)
rp.alltraits.step2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                      list(MYVALUE = format(summary(fit.alltraits.step2017)$adj.r.squared,dig=2)))[2] # extracting R2 value
rp.alltraits.step2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                      list(MYOTHERVALUE = format(predaccuracy.alltraits.step2017, digits = 2)))[2] # r (prediction accuracy) value
legend("bottomright", bty="n",legend=rp.alltraits.step2017)
R2.alltraits.step2017=summary(fit.alltraits.step2017)$adj.r.squared
R2.alltraits.step2017=round(R2.alltraits.step2017,digits=2)
r.alltraits.step2017=round(predaccuracy.alltraits.step2017,digits=2)
models="stepwise"
traits="alltraits"
pred.accuracy.alltraits.step2017=data.frame(cbind(traits,models,r.alltraits.step2017,R2.alltraits.step2017))
colnames(pred.accuracy.alltraits.step2017)[3:4]=c("r","R2")

### Plotting LASSO regression model with all traits together
plot(BLUE2017$GRYLD, predict.alltraits.las2017, xlab = "Observed Yield",ylab = "Predicted Yield", main="LASSO regression with all traits 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.las2017 <-lm(predict.alltraits.las2017 ~ BLUE2017$GRYLD) #fit the model for observed values
abline(fit.alltraits.las2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.las2017 = vector('expression',2)
rp.alltraits.las2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(summary(fit.alltraits.las2017)$adj.r.squared,dig=2)))[2]
rp.alltraits.las2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                     list(MYOTHERVALUE = format(predaccuracy.alltraits.lasso2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.las2017)
R2.alltraits.las2017=summary(fit.alltraits.las2017)$adj.r.squared
R2.alltraits.las2017=round(R2.alltraits.las2017,digits=2)
r.alltraits.las2017=round(predaccuracy.alltraits.lasso2017,digits=2)
models="LASSO"
traits="alltraits"
pred.accuracy.alltraits.las2017=data.frame(cbind(traits,models,r.alltraits.las2017,R2.alltraits.las2017))
colnames(pred.accuracy.alltraits.las2017)[3:4]=c("r","R2")

### Plotting ridge regression model with all traits together
plot(BLUE2017$GRYLD, predict.alltraits.rdg2017, xlab = "Observed Yield",ylab = "Predicted Yield", main="Ridge regression with all traits 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.rdg2017 <-lm(predict.alltraits.rdg2017 ~ BLUE2017$GRYLD) #fit the model for observed values
abline(fit.alltraits.rdg2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.rdg2017 = vector('expression',2)
rp.alltraits.rdg2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(summary(fit.alltraits.rdg2017)$adj.r.squared,dig=2)))[2]
rp.alltraits.rdg2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                     list(MYOTHERVALUE = format(predaccuracy.alltraits.ridge2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.rdg2017)
R2.alltraits.rdg2017=summary(fit.alltraits.rdg2017)$adj.r.squared
R2.alltraits.rdg2017=round(R2.alltraits.rdg2017,digits=2)
r.alltraits.rdg2017=round(predaccuracy.alltraits.ridge2017,digits=2)
models="ridge"
traits="alltraits"
pred.accuracy.alltraits.rdg2017=data.frame(cbind(traits,models,r.alltraits.rdg2017,R2.alltraits.rdg2017))
colnames(pred.accuracy.alltraits.rdg2017)[3:4]=c("r","R2")

### Plotting enet regression model with all traits together
plot(BLUE2017$GRYLD, predict.alltraits.enet2017, xlab = "Observed Yield",ylab = "Predicted Yield", main="ElasticNet regression with all traits 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.enet2017 <-lm(predict.alltraits.enet2017 ~ BLUE2017$GRYLD) #fit the model for observed values
abline(fit.alltraits.enet2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.enet2017 = vector('expression',2)
rp.alltraits.enet2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                      list(MYVALUE = format(summary(fit.alltraits.enet2017)$adj.r.squared,dig=2)))[2]
rp.alltraits.enet2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                      list(MYOTHERVALUE = format(predaccuracy.alltraits.enet2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.enet2017)
R2.alltraits.enet2017=summary(fit.alltraits.enet2017)$adj.r.squared
R2.alltraits.enet2017=round(R2.alltraits.enet2017,digits=2)
r.alltraits.enet2017=round(predaccuracy.alltraits.enet2017,digits=2)
models="elasticnet"
traits="alltraits"
pred.accuracy.alltraits.enet2017=data.frame(cbind(traits,models,r.alltraits.enet2017,R2.alltraits.enet2017))
colnames(pred.accuracy.alltraits.enet2017)[3:4]=c("r","R2")

multivar.alltraits.prediction.accuracy2017=data.frame(rbind(pred.accuracy.alltraits.step2017,pred.accuracy.alltraits.las2017,pred.accuracy.alltraits.rdg2017,pred.accuracy.alltraits.enet2017))
# write.csv(multivar.alltraits.prediction.accuracy2017,file="Multivariate_predaccuracy_alltraits_2017.csv",row.names = FALSE,quote = FALSE)


### Multivariate model yield prediction 2018 with all traits together
### Running stepwise regression model with all traits
BLUE2018=read.csv("BLUE_2018.csv")
stepwise.alltraits.2018 <- lm(GRYLD ~ CT_20180126+CT_20180131+CT_20180205+CT_20180210+CT_20180214+CT_20180219+CT_20180225+CT_20180301+CT_20180305+CT_20180310+CT_20180315+CT_20180320+
                                      NDVI_20180126+NDVI_20180131+NDVI_20180204+NDVI_20180210+NDVI_20180214+NDVI_20180219+NDVI_20180226+NDVI_20180301+NDVI_20180305+NDVI_20180310+NDVI_20180315+NDVI_20180320+
                                      DTHD+DAYSMT+PH+SN+SPKLNG+SPLN+GRNSPK+TGW,data=BLUE2018)  ## writing the model with all traits together
stepwise.alltraits.2018 <- step(stepwise.alltraits.2018, direction = "both") ## running the stepwise regressing model
formula(stepwise.alltraits.2018) ## getting the selected variables from the stepwise regression model

### Predicted yield from stepwise regression model with all traits together
step.predyield.alltraits.2018=NULL 
for(i in 1:11){
        step.prediction.alltraits.2018=BLUE2018[BLUE2018$trial==i,] #get just one trial to predict
        step.training.alltraits.2018=BLUE2018[BLUE2018$trial!=i,] #get all other trials to train on
        fit.alltraits.2018 = lm(formula = formula(stepwise.alltraits.2018), data =step.training.alltraits.2018) ## directly import variables into the model from stepwise regression 
        step.predvalue.alltraits.2018=predict(object=fit.alltraits.2018, newdata=step.prediction.alltraits.2018) ## Predicted value of individual trial
        predicted.alltraits.trial2018=data.frame(entry=step.prediction.alltraits.2018$plot,step.predvalue.alltraits.2018) ## Predicted value of individual trial in data frame
        step.predyield.alltraits.2018=rbind(step.predyield.alltraits.2018, predicted.alltraits.trial2018) ## Combinig all trial's predicted value into one frame
}
step.predyield.alltraits.2018 = data.frame(trial=rep(1:11,each=60),step.predyield.alltraits.2018)
step.predyield.alltraits.2018 <- merge(step.predyield.alltraits.2018, BLUE2018[, c(1:2,35)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extracting prediction accuracy from stepwise regression model with all traits together
cor.alltraits.step2018 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.alltraits.2018[step.predyield.alltraits.2018$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.alltraits.2018, trial_data$GRYLD) #get correlations
        cor.alltraits.step2018 <- c(cor.alltraits.step2018, trial_cor) #write results outside of the for loop
}

### Prediction from shrinkage models
predict.alltraits.las2018=c()  #Intialize to capture predicted values from LASSO model
predict.alltraits.rdg2018=c() #Intialize to capture predicted values from Ridge model
predict.alltraits.enet2018=c() #Intialize to capture predicted values from ElasticNet model
cor.alltraits.las2018=c() #Intialize to capture correlatio from LASSO model
cor.alltraits.rdg2018=c() #Intialize to capture correlatio from Ridge model
cor.alltraits.enet2018=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.alltraits.data2018 = BLUE2018[BLUE2018$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.alltraits.data2018 = BLUE2018[BLUE2018$trial!=trial,] #sets up prediction populations
        rctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.alltraits.model2018 = train(GRYLD ~ ., data = train.alltraits.data2018[,-c(1,2)], #Looks like we are fitting a lasso model.
                                          method = "lasso",
                                          trControl = rctrl,
                                          preProc = c("center", "scale"))
        ridge.alltraits.model2018 = train(GRYLD ~ ., data = train.alltraits.data2018[,-c(1,2)], #looks like we are fitting a ridge model.
                                          method = "ridge",
                                          trControl = rctrl,
                                          preProc = c("center", "scale"))
        enet.alltraits.model2018 = train(GRYLD ~ ., data = train.alltraits.data2018[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                         method = "enet",
                                         trControl = rctrl,
                                         preProc = c("center", "scale"))
        #model coefficients
        coef.alltraits.lasso2018 = lasso.alltraits.model2018$finalModel #extracting some coefficients from LASSO model
        coef.alltraits.ridge2018 = ridge.alltraits.model2018$finalModel #extracting some coefficients from Ridge model
        coef.alltraits.enet2018 = enet.alltraits.model2018$finalModel #extracting some coefficients from ElasticNet model
        #prediction
        pred.alltraits.lasso2018 = predict(lasso.alltraits.model2018,test.alltraits.data2018[,-c(1,2)]) #making a predict statement.  This looks promising
        pred.alltraits.ridge2018 = predict(ridge.alltraits.model2018,test.alltraits.data2018[,-c(1,2)]) #appears the same with other models
        pred.alltraits.elasticnet2018 = predict(enet.alltraits.model2018,test.alltraits.data2018[,-c(1,2)])
        #extracting correlations
        cor.alltraits.l2018=cor(BLUE2018[BLUE2018$trial==trial,35],pred.alltraits.lasso2018) # Correlaiton between predicted value from LASSO regression and BLUE
        cor.alltraits.las2018=c(cor.alltraits.las2018,cor.alltraits.l2018) #writing correlation out to cor.alltraits.las2018
        cor.alltraits.r2018=cor(BLUE2018[BLUE2018$trial==trial,35],pred.alltraits.ridge2018) # Correlaiton between predicted value from Ridge regression and BLUE
        cor.alltraits.rdg2018=c(cor.alltraits.rdg2018,cor.alltraits.r2018) #writing correlation out to cor.alltraits.rdg2018
        cor.alltraits.e2018=cor(BLUE2018[BLUE2018$trial==trial,35],pred.alltraits.elasticnet2018) # Correlaiton between predicted value from ElasticNet regression and BLUE
        cor.alltraits.enet2018=c(cor.alltraits.enet2018,cor.alltraits.e2018) #writing correlation out to cor.alltraits.enet2018
        #extracting predicted values
        predict.alltraits.las2018 = c(predict.alltraits.las2018,pred.alltraits.lasso2018) # Combinig all predicted values from LASSO
        predict.alltraits.rdg2018 = c(predict.alltraits.rdg2018,pred.alltraits.ridge2018) # Combinig all predicted values from Ridge
        predict.alltraits.enet2018 = c(predict.alltraits.enet2018,pred.alltraits.elasticnet2018) # Combinig all predicted values from ElasticNet
} #end for loop 

#print prediction accuracies
cor.alltraits.step2018 #individual trial correlations for step regression
predaccuracy.alltraits.step2018=mean(cor.alltraits.step2018) #average correlation from stepwise model
cor.alltraits.las2018 #individual trial correlations for LASSO regression
predaccuracy.alltraits.lasso2018=mean(cor.alltraits.las2018) #average correlation from LASSO model
cor.alltraits.rdg2018 #individual trial correlations for Ridge regression
predaccuracy.alltraits.ridge2018=mean(cor.alltraits.rdg2018) #average correlation from Ridge model
cor.alltraits.enet2018 #individual trial correlations for ElasticNet regression
predaccuracy.alltraits.enet2018=mean(cor.alltraits.enet2018) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with all traits together
plot(step.predyield.alltraits.2018$GRYLD, step.predyield.alltraits.2018$step.predvalue.alltraits.2018, xlab = "Observed Yield",ylab = "Predicted Yield", main="Stepwise regression with all traits 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.step2018 <-lm(step.predyield.alltraits.2018$step.predvalue.alltraits.2018 ~ step.predyield.alltraits.2018$GRYLD) #fit the model for observed values
abline(fit.alltraits.step2018, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.step2018 = vector('expression',2)
rp.alltraits.step2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                      list(MYVALUE = format(summary(fit.alltraits.step2018)$adj.r.squared,dig=2)))[2] # extracting R2 value
rp.alltraits.step2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                      list(MYOTHERVALUE = format(predaccuracy.alltraits.step2018, digits = 2)))[2] # r (prediction accuracy) value
legend("bottomright", bty="n",legend=rp.alltraits.step2018)
R2.alltraits.step2018=summary(fit.alltraits.step2018)$adj.r.squared
R2.alltraits.step2018=round(R2.alltraits.step2018,digits=2)
r.alltraits.step2018=round(predaccuracy.alltraits.step2018,digits=2)
models="stepwise"
traits="alltraits"
pred.accuracy.alltraits.step2018=data.frame(cbind(traits,models,r.alltraits.step2018,R2.alltraits.step2018))
colnames(pred.accuracy.alltraits.step2018)[3:4]=c("r","R2")

### Plotting LASSO regression model with all traits together
plot(BLUE2018$GRYLD, predict.alltraits.las2018, xlab = "Observed Yield",ylab = "Predicted Yield", main="LASSO regression with all traits 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.las2018 <-lm(predict.alltraits.las2018 ~ BLUE2018$GRYLD) #fit the model for observed values
abline(fit.alltraits.las2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.las2018 = vector('expression',2)
rp.alltraits.las2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(summary(fit.alltraits.las2018)$adj.r.squared,dig=2)))[2]
rp.alltraits.las2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                     list(MYOTHERVALUE = format(predaccuracy.alltraits.lasso2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.las2018)
R2.alltraits.las2018=summary(fit.alltraits.las2018)$adj.r.squared
R2.alltraits.las2018=round(R2.alltraits.las2018,digits=2)
r.alltraits.las2018=round(predaccuracy.alltraits.lasso2018,digits=2)
models="LASSO"
traits="alltraits"
pred.accuracy.alltraits.las2018=data.frame(cbind(traits,models,r.alltraits.las2018,R2.alltraits.las2018))
colnames(pred.accuracy.alltraits.las2018)[3:4]=c("r","R2")

### Plotting ridge regression model with all traits together
plot(BLUE2018$GRYLD, predict.alltraits.rdg2018, xlab = "Observed Yield",ylab = "Predicted Yield", main="Ridge regression with all traits 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.rdg2018 <-lm(predict.alltraits.rdg2018 ~ BLUE2018$GRYLD) #fit the model for observed values
abline(fit.alltraits.rdg2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.rdg2018 = vector('expression',2)
rp.alltraits.rdg2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(summary(fit.alltraits.rdg2018)$adj.r.squared,dig=2)))[2]
rp.alltraits.rdg2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                     list(MYOTHERVALUE = format(predaccuracy.alltraits.ridge2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.rdg2018)
R2.alltraits.rdg2018=summary(fit.alltraits.rdg2018)$adj.r.squared
R2.alltraits.rdg2018=round(R2.alltraits.rdg2018,digits=2)
r.alltraits.rdg2018=round(predaccuracy.alltraits.ridge2018,digits=2)
models="ridge"
traits="alltraits"
pred.accuracy.alltraits.rdg2018=data.frame(cbind(traits,models,r.alltraits.rdg2018,R2.alltraits.rdg2018))
colnames(pred.accuracy.alltraits.rdg2018)[3:4]=c("r","R2")

### Plotting enet regression model with all traits together
plot(BLUE2018$GRYLD, predict.alltraits.enet2018, xlab = "Observed Yield",ylab = "Predicted Yield", main="ElasticNet regression with all traits 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.enet2018 <-lm(predict.alltraits.enet2018 ~ BLUE2018$GRYLD) #fit the model for observed values
abline(fit.alltraits.enet2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.enet2018 = vector('expression',2)
rp.alltraits.enet2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                      list(MYVALUE = format(summary(fit.alltraits.enet2018)$adj.r.squared,dig=2)))[2]
rp.alltraits.enet2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                      list(MYOTHERVALUE = format(predaccuracy.alltraits.enet2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.enet2018)
R2.alltraits.enet2018=summary(fit.alltraits.enet2018)$adj.r.squared
R2.alltraits.enet2018=round(R2.alltraits.enet2018,digits=2)
r.alltraits.enet2018=round(predaccuracy.alltraits.enet2018,digits=2)
models="elasticnet"
traits="alltraits"
pred.accuracy.alltraits.enet2018=data.frame(cbind(traits,models,r.alltraits.enet2018,R2.alltraits.enet2018))
colnames(pred.accuracy.alltraits.enet2018)[3:4]=c("r","R2")

multivar.alltraits.prediction.accuracy2018=data.frame(rbind(pred.accuracy.alltraits.step2018,pred.accuracy.alltraits.las2018,pred.accuracy.alltraits.rdg2018,pred.accuracy.alltraits.enet2018))
# write.csv(multivar.alltraits.prediction.accuracy2018,file="Multivariate_predaccuracy_alltraits_2018.csv",row.names = FALSE,quote = FALSE)


### Multivariate model yield prediction 2019 with all traits together
### Running stepwise regression model with all traits
BLUE2019=read.csv("BLUE_2019.csv")
stepwise.alltraits.2019 <- lm(GRYLD ~ CT_20190123+CT_20190127+CT_20190131+CT_20190205+CT_20190211+CT_20190218+CT_20190223+CT_20190301+CT_20190305+CT_20190311+CT_20190316+CT_20190320+CT_20190325+
                                      NDVI_20190121+NDVI_20190127+NDVI_20190131+NDVI_20190205+NDVI_20190211+NDVI_20190218+NDVI_20190222+NDVI_20190228+NDVI_20190305+NDVI_20190311+NDVI_20190315+NDVI_20190320+NDVI_20190325+
                                      DTHD+DAYSMD+PH+SN+SPKLNG+SPLN+GRNSPK+TGW,data=BLUE2019)  ## writing the model with all traits together
stepwise.alltraits.2019 <- step(stepwise.alltraits.2019, direction = "both") ## running the stepwise regressing model
formula(stepwise.alltraits.2019) ## getting the selected variables from the stepwise regression model

### Predicted yield from stepwise regression model with all traits together
step.predyield.alltraits.2019=NULL 
for(i in 1:10){
        step.prediction.alltraits.2019=BLUE2019[BLUE2019$trial==i,] #get just one trial to predict
        step.training.alltraits.2019=BLUE2019[BLUE2019$trial!=i,] #get all other trials to train on
        fit.alltraits.2019 = lm(formula = formula(stepwise.alltraits.2019), data =step.training.alltraits.2019) ## directly import variables into the model from stepwise regression 
        step.predvalue.alltraits.2019=predict(object=fit.alltraits.2019, newdata=step.prediction.alltraits.2019) ## Predicted value of individual trial
        predicted.alltraits.trial2019=data.frame(entry=step.prediction.alltraits.2019$plot,step.predvalue.alltraits.2019) ## Predicted value of individual trial in data frame
        step.predyield.alltraits.2019=rbind(step.predyield.alltraits.2019, predicted.alltraits.trial2019) ## Combinig all trial's predicted value into one frame
}
step.predyield.alltraits.2019 = data.frame(trial=rep(1:10,each=60),step.predyield.alltraits.2019)
step.predyield.alltraits.2019 <- merge(step.predyield.alltraits.2019, BLUE2019[, c(1:2,37)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extracting prediction accuracy from stepwise regression model with all traits together
cor.alltraits.step2019 <- NULL #start dataframe to hold correlation values
for(i in 1:10){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.alltraits.2019[step.predyield.alltraits.2019$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.alltraits.2019, trial_data$GRYLD) #get correlations
        cor.alltraits.step2019 <- c(cor.alltraits.step2019, trial_cor) #write results outside of the for loop
}

### Prediction from shrinkage models
predict.alltraits.las2019=c()  #Intialize to capture predicted values from LASSO model
predict.alltraits.rdg2019=c() #Intialize to capture predicted values from Ridge model
predict.alltraits.enet2019=c() #Intialize to capture predicted values from ElasticNet model
cor.alltraits.las2019=c() #Intialize to capture correlatio from LASSO model
cor.alltraits.rdg2019=c() #Intialize to capture correlatio from Ridge model
cor.alltraits.enet2019=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:10){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.alltraits.data2019 = BLUE2019[BLUE2019$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.alltraits.data2019 = BLUE2019[BLUE2019$trial!=trial,] #sets up prediction populations
        rctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.alltraits.model2019 = train(GRYLD ~ ., data = train.alltraits.data2019[,-c(1,2)], #Looks like we are fitting a lasso model.
                                          method = "lasso",
                                          trControl = rctrl,
                                          preProc = c("center", "scale"))
        ridge.alltraits.model2019 = train(GRYLD ~ ., data = train.alltraits.data2019[,-c(1,2)], #looks like we are fitting a ridge model.
                                          method = "ridge",
                                          trControl = rctrl,
                                          preProc = c("center", "scale"))
        enet.alltraits.model2019 = train(GRYLD ~ ., data = train.alltraits.data2019[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                         method = "enet",
                                         trControl = rctrl,
                                         preProc = c("center", "scale"))
        #model coefficients
        coef.alltraits.lasso2019 = lasso.alltraits.model2019$finalModel #extracting some coefficients from LASSO model
        coef.alltraits.ridge2019 = ridge.alltraits.model2019$finalModel #extracting some coefficients from Ridge model
        coef.alltraits.enet2019 = enet.alltraits.model2019$finalModel #extracting some coefficients from ElasticNet model
        #prediction
        pred.alltraits.lasso2019 = predict(lasso.alltraits.model2019,test.alltraits.data2019[,-c(1,2)]) #making a predict statement.  This looks promising
        pred.alltraits.ridge2019 = predict(ridge.alltraits.model2019,test.alltraits.data2019[,-c(1,2)]) #appears the same with other models
        pred.alltraits.elasticnet2019 = predict(enet.alltraits.model2019,test.alltraits.data2019[,-c(1,2)])
        #extracting correlations
        cor.alltraits.l2019=cor(BLUE2019[BLUE2019$trial==trial,37],pred.alltraits.lasso2019) # Correlaiton between predicted value from LASSO regression and BLUE
        cor.alltraits.las2019=c(cor.alltraits.las2019,cor.alltraits.l2019) #writing correlation out to cor.alltraits.las2019
        cor.alltraits.r2019=cor(BLUE2019[BLUE2019$trial==trial,37],pred.alltraits.ridge2019) # Correlaiton between predicted value from Ridge regression and BLUE
        cor.alltraits.rdg2019=c(cor.alltraits.rdg2019,cor.alltraits.r2019) #writing correlation out to cor.alltraits.rdg2019
        cor.alltraits.e2019=cor(BLUE2019[BLUE2019$trial==trial,37],pred.alltraits.elasticnet2019) # Correlaiton between predicted value from ElasticNet regression and BLUE
        cor.alltraits.enet2019=c(cor.alltraits.enet2019,cor.alltraits.e2019) #writing correlation out to cor.alltraits.enet2019
        #extracting predicted values
        predict.alltraits.las2019 = c(predict.alltraits.las2019,pred.alltraits.lasso2019) # Combinig all predicted values from LASSO
        predict.alltraits.rdg2019 = c(predict.alltraits.rdg2019,pred.alltraits.ridge2019) # Combinig all predicted values from Ridge
        predict.alltraits.enet2019 = c(predict.alltraits.enet2019,pred.alltraits.elasticnet2019) # Combinig all predicted values from ElasticNet
} #end for loop 

#print prediction accuracies
cor.alltraits.step2019 #individual trial correlations for step regression
predaccuracy.alltraits.step2019=mean(cor.alltraits.step2019) #average correlation from stepwise model
cor.alltraits.las2019 #individual trial correlations for LASSO regression
predaccuracy.alltraits.lasso2019=mean(cor.alltraits.las2019) #average correlation from LASSO model
cor.alltraits.rdg2019 #individual trial correlations for Ridge regression
predaccuracy.alltraits.ridge2019=mean(cor.alltraits.rdg2019) #average correlation from Ridge model
cor.alltraits.enet2019 #individual trial correlations for ElasticNet regression
predaccuracy.alltraits.enet2019=mean(cor.alltraits.enet2019) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with all traits together
plot(step.predyield.alltraits.2019$GRYLD, step.predyield.alltraits.2019$step.predvalue.alltraits.2019, xlab = "Observed Yield",ylab = "Predicted Yield", main="Stepwise regression with all traits 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.step2019 <-lm(step.predyield.alltraits.2019$step.predvalue.alltraits.2019 ~ step.predyield.alltraits.2019$GRYLD) #fit the model for observed values
abline(fit.alltraits.step2019, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.step2019 = vector('expression',2)
rp.alltraits.step2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                      list(MYVALUE = format(summary(fit.alltraits.step2019)$adj.r.squared,dig=2)))[2] # extracting R2 value
rp.alltraits.step2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                      list(MYOTHERVALUE = format(predaccuracy.alltraits.step2019, digits = 2)))[2] # r (prediction accuracy) value
legend("bottomright", bty="n",legend=rp.alltraits.step2019)
R2.alltraits.step2019=summary(fit.alltraits.step2019)$adj.r.squared
R2.alltraits.step2019=round(R2.alltraits.step2019,digits=2)
r.alltraits.step2019=round(predaccuracy.alltraits.step2019,digits=2)
models="stepwise"
traits="alltraits"
pred.accuracy.alltraits.step2019=data.frame(cbind(traits,models,r.alltraits.step2019,R2.alltraits.step2019))
colnames(pred.accuracy.alltraits.step2019)[3:4]=c("r","R2")

### Plotting LASSO regression model with all traits together
plot(BLUE2019$GRYLD, predict.alltraits.las2019, xlab = "Observed Yield",ylab = "Predicted Yield", main="LASSO regression with all traits 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.las2019 <-lm(predict.alltraits.las2019 ~ BLUE2019$GRYLD) #fit the model for observed values
abline(fit.alltraits.las2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.las2019 = vector('expression',2)
rp.alltraits.las2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(summary(fit.alltraits.las2019)$adj.r.squared,dig=2)))[2]
rp.alltraits.las2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                     list(MYOTHERVALUE = format(predaccuracy.alltraits.lasso2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.las2019)
R2.alltraits.las2019=summary(fit.alltraits.las2019)$adj.r.squared
R2.alltraits.las2019=round(R2.alltraits.las2019,digits=2)
r.alltraits.las2019=round(predaccuracy.alltraits.lasso2019,digits=2)
models="LASSO"
traits="alltraits"
pred.accuracy.alltraits.las2019=data.frame(cbind(traits,models,r.alltraits.las2019,R2.alltraits.las2019))
colnames(pred.accuracy.alltraits.las2019)[3:4]=c("r","R2")

### Plotting ridge regression model with all traits together
plot(BLUE2019$GRYLD, predict.alltraits.rdg2019, xlab = "Observed Yield",ylab = "Predicted Yield", main="Ridge regression with all traits 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.rdg2019 <-lm(predict.alltraits.rdg2019 ~ BLUE2019$GRYLD) #fit the model for observed values
abline(fit.alltraits.rdg2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.rdg2019 = vector('expression',2)
rp.alltraits.rdg2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(summary(fit.alltraits.rdg2019)$adj.r.squared,dig=2)))[2]
rp.alltraits.rdg2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                     list(MYOTHERVALUE = format(predaccuracy.alltraits.ridge2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.rdg2019)
R2.alltraits.rdg2019=summary(fit.alltraits.rdg2019)$adj.r.squared
R2.alltraits.rdg2019=round(R2.alltraits.rdg2019,digits=2)
r.alltraits.rdg2019=round(predaccuracy.alltraits.ridge2019,digits=2)
models="ridge"
traits="alltraits"
pred.accuracy.alltraits.rdg2019=data.frame(cbind(traits,models,r.alltraits.rdg2019,R2.alltraits.rdg2019))
colnames(pred.accuracy.alltraits.rdg2019)[3:4]=c("r","R2")

### Plotting enet regression model with all traits together
plot(BLUE2019$GRYLD, predict.alltraits.enet2019, xlab = "Observed Yield",ylab = "Predicted Yield", main="ElasticNet regression with all traits 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.enet2019 <-lm(predict.alltraits.enet2019 ~ BLUE2019$GRYLD) #fit the model for observed values
abline(fit.alltraits.enet2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.enet2019 = vector('expression',2)
rp.alltraits.enet2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                      list(MYVALUE = format(summary(fit.alltraits.enet2019)$adj.r.squared,dig=2)))[2]
rp.alltraits.enet2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                      list(MYOTHERVALUE = format(predaccuracy.alltraits.enet2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.enet2019)
R2.alltraits.enet2019=summary(fit.alltraits.enet2019)$adj.r.squared
R2.alltraits.enet2019=round(R2.alltraits.enet2019,digits=2)
r.alltraits.enet2019=round(predaccuracy.alltraits.enet2019,digits=2)
models="elasticnet"
traits="alltraits"
pred.accuracy.alltraits.enet2019=data.frame(cbind(traits,models,r.alltraits.enet2019,R2.alltraits.enet2019))
colnames(pred.accuracy.alltraits.enet2019)[3:4]=c("r","R2")

multivar.alltraits.prediction.accuracy2019=data.frame(rbind(pred.accuracy.alltraits.step2019,pred.accuracy.alltraits.las2019,pred.accuracy.alltraits.rdg2019,pred.accuracy.alltraits.enet2019))
# write.csv(multivar.alltraits.prediction.accuracy2019,file="Multivariate_predaccuracy_alltraits_2019.csv",row.names = FALSE,quote = FALSE)


### Multivariate model yield prediction 2020 with all traits together
### Running stepwise regression model with all traits
BLUE2020=read.csv("BLUE_2020.csv")
stepwise.alltraits.2020 <- lm(GRYLD ~ CT_20200112+CT_20200116+CT_20200121+CT_20200126+CT_20200130+CT_20200205+CT_20200210+CT_20200215+CT_20200220+CT_20200226+CT_20200302+CT_20200308+CT_20200313+CT_20200318+CT_20200323
                              +NDVI_20200112+NDVI_20200116+NDVI_20200121+NDVI_20200126+NDVI_20200130+NDVI_20200205+NDVI_20200210+NDVI_20200215+NDVI_20200220+NDVI_20200226+NDVI_20200302+NDVI_20200308+NDVI_20200313+NDVI_20200318+NDVI_20200323
                              +GrndCov_20200112+GrndCov_20200206+DLA_Feb26+DLA_Mar09+DTHD+DAYSMT+PH+SN+SPLN+GRNSPK+TGW,data=BLUE2020)  ## writing the model with all traits together
stepwise.alltraits.2020 <- step(stepwise.alltraits.2020, direction = "both") ## running the stepwise regressing model
formula(stepwise.alltraits.2020) ## getting the selected variables from the stepwise regression model

### Predicted yield from stepwise regression model with all traits together
step.predyield.alltraits.2020=NULL 
for(i in 1:11){
        step.prediction.alltraits.2020=BLUE2020[BLUE2020$trial==i,] #get just one trial to predict
        step.training.alltraits.2020=BLUE2020[BLUE2020$trial!=i,] #get all other trials to train on
        fit.alltraits.2020 = lm(formula = formula(stepwise.alltraits.2020), data =step.training.alltraits.2020) ## directly import variables into the model from stepwise regression 
        step.predvalue.alltraits.2020=predict(object=fit.alltraits.2020, newdata=step.prediction.alltraits.2020) ## Predicted value of individual trial
        predicted.alltraits.trial2020=data.frame(entry=step.prediction.alltraits.2020$plot,step.predvalue.alltraits.2020) ## Predicted value of individual trial in data frame
        step.predyield.alltraits.2020=rbind(step.predyield.alltraits.2020, predicted.alltraits.trial2020) ## Combinig all trial's predicted value into one frame
}
step.predyield.alltraits.2020 = data.frame(trial=rep(1:11,each=60),step.predyield.alltraits.2020)
step.predyield.alltraits.2020 <- merge(step.predyield.alltraits.2020, BLUE2020[, c(1:2,44)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extracting prediction accuracy from stepwise regression model with all traits together
cor.alltraits.step2020 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.alltraits.2020[step.predyield.alltraits.2020$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.alltraits.2020, trial_data$GRYLD) #get correlations
        cor.alltraits.step2020 <- c(cor.alltraits.step2020, trial_cor) #write results outside of the for loop
}

### Prediction from shrinkage models
predict.alltraits.las2020=c()  #Intialize to capture predicted values from LASSO model
predict.alltraits.rdg2020=c() #Intialize to capture predicted values from Ridge model
predict.alltraits.enet2020=c() #Intialize to capture predicted values from ElasticNet model
cor.alltraits.las2020=c() #Intialize to capture correlatio from LASSO model
cor.alltraits.rdg2020=c() #Intialize to capture correlatio from Ridge model
cor.alltraits.enet2020=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.alltraits.data2020 = BLUE2020[BLUE2020$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.alltraits.data2020 = BLUE2020[BLUE2020$trial!=trial,] #sets up prediction populations
        rctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.alltraits.model2020 = train(GRYLD ~ ., data = train.alltraits.data2020[,-c(1,2)], #Looks like we are fitting a lasso model.
                                          method = "lasso",
                                          trControl = rctrl,
                                          preProc = c("center", "scale"))
        ridge.alltraits.model2020 = train(GRYLD ~ ., data = train.alltraits.data2020[,-c(1,2)], #looks like we are fitting a ridge model.
                                          method = "ridge",
                                          trControl = rctrl,
                                          preProc = c("center", "scale"))
        enet.alltraits.model2020 = train(GRYLD ~ ., data = train.alltraits.data2020[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                         method = "enet",
                                         trControl = rctrl,
                                         preProc = c("center", "scale"))
        #model coefficients
        coef.alltraits.lasso2020 = lasso.alltraits.model2020$finalModel #extracting some coefficients from LASSO model
        coef.alltraits.ridge2020 = ridge.alltraits.model2020$finalModel #extracting some coefficients from Ridge model
        coef.alltraits.enet2020 = enet.alltraits.model2020$finalModel #extracting some coefficients from ElasticNet model
        #prediction
        pred.alltraits.lasso2020 = predict(lasso.alltraits.model2020,test.alltraits.data2020[,-c(1,2)]) #making a predict statement.  This looks promising
        pred.alltraits.ridge2020 = predict(ridge.alltraits.model2020,test.alltraits.data2020[,-c(1,2)]) #appears the same with other models
        pred.alltraits.elasticnet2020 = predict(enet.alltraits.model2020,test.alltraits.data2020[,-c(1,2)])
        #extracting correlations
        cor.alltraits.l2020=cor(BLUE2020[BLUE2020$trial==trial,44],pred.alltraits.lasso2020) # Correlaiton between predicted value from LASSO regression and BLUE
        cor.alltraits.las2020=c(cor.alltraits.las2020,cor.alltraits.l2020) #writing correlation out to cor.alltraits.las2020
        cor.alltraits.r2020=cor(BLUE2020[BLUE2020$trial==trial,44],pred.alltraits.ridge2020) # Correlaiton between predicted value from Ridge regression and BLUE
        cor.alltraits.rdg2020=c(cor.alltraits.rdg2020,cor.alltraits.r2020) #writing correlation out to cor.alltraits.rdg2020
        cor.alltraits.e2020=cor(BLUE2020[BLUE2020$trial==trial,44],pred.alltraits.elasticnet2020) # Correlaiton between predicted value from ElasticNet regression and BLUE
        cor.alltraits.enet2020=c(cor.alltraits.enet2020,cor.alltraits.e2020) #writing correlation out to cor.alltraits.enet2020
        #extracting predicted values
        predict.alltraits.las2020 = c(predict.alltraits.las2020,pred.alltraits.lasso2020) # Combinig all predicted values from LASSO
        predict.alltraits.rdg2020 = c(predict.alltraits.rdg2020,pred.alltraits.ridge2020) # Combinig all predicted values from Ridge
        predict.alltraits.enet2020 = c(predict.alltraits.enet2020,pred.alltraits.elasticnet2020) # Combinig all predicted values from ElasticNet
} #end for loop 

#print prediction accuracies
cor.alltraits.step2020 #individual trial correlations for step regression
predaccuracy.alltraits.step2020=mean(cor.alltraits.step2020) #average correlation from stepwise model
cor.alltraits.las2020 #individual trial correlations for LASSO regression
predaccuracy.alltraits.lasso2020=mean(cor.alltraits.las2020) #average correlation from LASSO model
cor.alltraits.rdg2020 #individual trial correlations for Ridge regression
predaccuracy.alltraits.ridge2020=mean(cor.alltraits.rdg2020) #average correlation from Ridge model
cor.alltraits.enet2020 #individual trial correlations for ElasticNet regression
predaccuracy.alltraits.enet2020=mean(cor.alltraits.enet2020) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with all traits together
plot(step.predyield.alltraits.2020$GRYLD, step.predyield.alltraits.2020$step.predvalue.alltraits.2020, xlab = "Observed Yield",ylab = "Predicted Yield", main="Stepwise regression with all traits 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.step2020 <-lm(step.predyield.alltraits.2020$step.predvalue.alltraits.2020 ~ step.predyield.alltraits.2020$GRYLD) #fit the model for observed values
abline(fit.alltraits.step2020, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.step2020 = vector('expression',2)
rp.alltraits.step2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                      list(MYVALUE = format(summary(fit.alltraits.step2020)$adj.r.squared,dig=2)))[2] # extracting R2 value
rp.alltraits.step2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                      list(MYOTHERVALUE = format(predaccuracy.alltraits.step2020, digits = 2)))[2] # r (prediction accuracy) value
legend("bottomright", bty="n",legend=rp.alltraits.step2020)
R2.alltraits.step2020=summary(fit.alltraits.step2020)$adj.r.squared
R2.alltraits.step2020=round(R2.alltraits.step2020,digits=2)
r.alltraits.step2020=round(predaccuracy.alltraits.step2020,digits=2)
models="stepwise"
traits="alltraits"
pred.accuracy.alltraits.step2020=data.frame(cbind(traits,models,r.alltraits.step2020,R2.alltraits.step2020))
colnames(pred.accuracy.alltraits.step2020)[3:4]=c("r","R2")

### Plotting LASSO regression model with all traits together
plot(BLUE2020$GRYLD, predict.alltraits.las2020, xlab = "Observed Yield",ylab = "Predicted Yield", main="LASSO regression with all traits 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.las2020 <-lm(predict.alltraits.las2020 ~ BLUE2020$GRYLD) #fit the model for observed values
abline(fit.alltraits.las2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.las2020 = vector('expression',2)
rp.alltraits.las2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(summary(fit.alltraits.las2020)$adj.r.squared,dig=2)))[2]
rp.alltraits.las2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                     list(MYOTHERVALUE = format(predaccuracy.alltraits.lasso2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.las2020)
R2.alltraits.las2020=summary(fit.alltraits.las2020)$adj.r.squared
R2.alltraits.las2020=round(R2.alltraits.las2020,digits=2)
r.alltraits.las2020=round(predaccuracy.alltraits.lasso2020,digits=2)
models="LASSO"
traits="alltraits"
pred.accuracy.alltraits.las2020=data.frame(cbind(traits,models,r.alltraits.las2020,R2.alltraits.las2020))
colnames(pred.accuracy.alltraits.las2020)[3:4]=c("r","R2")

### Plotting ridge regression model with all traits together
plot(BLUE2020$GRYLD, predict.alltraits.rdg2020, xlab = "Observed Yield",ylab = "Predicted Yield", main="Ridge regression with all traits 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.rdg2020 <-lm(predict.alltraits.rdg2020 ~ BLUE2020$GRYLD) #fit the model for observed values
abline(fit.alltraits.rdg2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.rdg2020 = vector('expression',2)
rp.alltraits.rdg2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                     list(MYVALUE = format(summary(fit.alltraits.rdg2020)$adj.r.squared,dig=2)))[2]
rp.alltraits.rdg2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                     list(MYOTHERVALUE = format(predaccuracy.alltraits.ridge2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.rdg2020)
R2.alltraits.rdg2020=summary(fit.alltraits.rdg2020)$adj.r.squared
R2.alltraits.rdg2020=round(R2.alltraits.rdg2020,digits=2)
r.alltraits.rdg2020=round(predaccuracy.alltraits.ridge2020,digits=2)
models="ridge"
traits="alltraits"
pred.accuracy.alltraits.rdg2020=data.frame(cbind(traits,models,r.alltraits.rdg2020,R2.alltraits.rdg2020))
colnames(pred.accuracy.alltraits.rdg2020)[3:4]=c("r","R2")

### Plotting enet regression model with all traits together
plot(BLUE2020$GRYLD, predict.alltraits.enet2020, xlab = "Observed Yield",ylab = "Predicted Yield", main="ElasticNet regression with all traits 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.alltraits.enet2020 <-lm(predict.alltraits.enet2020 ~ BLUE2020$GRYLD) #fit the model for observed values
abline(fit.alltraits.enet2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.alltraits.enet2020 = vector('expression',2)
rp.alltraits.enet2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                      list(MYVALUE = format(summary(fit.alltraits.enet2020)$adj.r.squared,dig=2)))[2]
rp.alltraits.enet2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                      list(MYOTHERVALUE = format(predaccuracy.alltraits.enet2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.alltraits.enet2020)
R2.alltraits.enet2020=summary(fit.alltraits.enet2020)$adj.r.squared
R2.alltraits.enet2020=round(R2.alltraits.enet2020,digits=2)
r.alltraits.enet2020=round(predaccuracy.alltraits.enet2020,digits=2)
models="elasticnet"
traits="alltraits"
pred.accuracy.alltraits.enet2020=data.frame(cbind(traits,models,r.alltraits.enet2020,R2.alltraits.enet2020))
colnames(pred.accuracy.alltraits.enet2020)[3:4]=c("r","R2")

multivar.alltraits.prediction.accuracy2020=data.frame(rbind(pred.accuracy.alltraits.step2020,pred.accuracy.alltraits.las2020,pred.accuracy.alltraits.rdg2020,pred.accuracy.alltraits.enet2020))
# write.csv(multivar.alltraits.prediction.accuracy2020,file="Multivariate_predaccuracy_alltraits_2020.csv",row.names = FALSE,quote = FALSE)

multivar.alltrait.allYear=cbind(multivar.alltraits.prediction.accuracy2016,multivar.alltraits.prediction.accuracy2017[,3:4],multivar.alltraits.prediction.accuracy2018[,3:4],multivar.alltraits.prediction.accuracy2019[,3:4],multivar.alltraits.prediction.accuracy2020[,3:4])
write.csv(multivar.alltrait.allYear,file="multivar.yielpred.alltrait.allyear.csv",row.names=FALSE,quote=FALSE)
# multivar.alltrait.allYear$prediction.accuracy=as.numeric(as.character(multivar.alltrait.allYear$prediction.accuracy))
# multivar.alltrait.allYear$models=as.character(multivar.alltrait.allYear$models)
# barplot(multivar.alltrait.allYear$prediction.accuracy,names.arg = multivar.alltrait.allYear$models,las=2,ylim=c(0,1))

