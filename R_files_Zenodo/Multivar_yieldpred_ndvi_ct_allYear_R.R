### Multivariate model yield predindvi_ction 2016 with ndvi_ct
### Running stepwise regression model with ndvi_ct
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
ndvi_ct_data=BLUE2016[,c(1:19,28)]
stepwise.ndvi_ct.2016 <- lm(GRYLD ~ CT_20160123+CT_20160204+CT_20160212+CT_20160223+CT_20160228+CT_20160302+CT_20160309+CT_20160315+NDVI_20160121+NDVI_20160130+NDVI_20160203+NDVI_20160207+NDVI_20160223+NDVI_20160228+NDVI_20160303+NDVI_20160310+NDVI_20160315
                            ,data=ndvi_ct_data)  ## writing the model with ndvi_ct
stepwise.ndvi_ct.2016 <- step(stepwise.ndvi_ct.2016, direndvi_ction = "both") ## running the stepwise regressing model
formula(stepwise.ndvi_ct.2016) ## getting the selendvi_cted variables from the stepwise regression model

### Predindvi_cted yield from stepwise regression model with ndvi_ct
step.predyield.ndvi_ct.2016=NULL 
for(i in 1:10){
        step.prediction.ndvi_ct.2016=ndvi_ct_data[ndvi_ct_data$trial==i,] #get just one trial to predindvi_ct
        step.training.ndvi_ct.2016=ndvi_ct_data[ndvi_ct_data$trial!=i,] #get all other trials to train on
        fit.ndvi_ct.2016 = lm(formula = formula(stepwise.ndvi_ct.2016), data =step.training.ndvi_ct.2016) ## direndvi_ctly import variables into the model from stepwise regression 
        step.predvalue.ndvi_ct.2016=predict(object=fit.ndvi_ct.2016, newdata=step.prediction.ndvi_ct.2016) ## Predindvi_cted value of individual trial
        predicted.ndvi_ct.trial2016=data.frame(entry=step.prediction.ndvi_ct.2016$plot,step.predvalue.ndvi_ct.2016) ## Predindvi_cted value of individual trial in data frame
        step.predyield.ndvi_ct.2016=rbind(step.predyield.ndvi_ct.2016, predicted.ndvi_ct.trial2016) ## Combinig all trial's predindvi_cted value into one frame
}
step.predyield.ndvi_ct.2016 = data.frame(trial=rep(1:10,each=60),step.predyield.ndvi_ct.2016)
step.predyield.ndvi_ct.2016 <- merge(step.predyield.ndvi_ct.2016, ndvi_ct_data[, c(1:2,20)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extrandvi_cting predindvi_ction accuracy from stepwise regression model with ndvi_ct
cor.ndvi_ct.step2016 <- NULL #start dataframe to hold correlation values
for(i in 1:10){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ndvi_ct.2016[step.predyield.ndvi_ct.2016$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ndvi_ct.2016, trial_data$GRYLD) #get correlations
        cor.ndvi_ct.step2016 <- c(cor.ndvi_ct.step2016, trial_cor) #write results outside of the for loop
}

### Predindvi_ction from shrinkage models
predindvi_ct.ndvi_ct.las2016=c()  #Intialize to capture predindvi_cted values from LASSO model
predindvi_ct.ndvi_ct.rdg2016=c() #Intialize to capture predindvi_cted values from Ridge model
predindvi_ct.ndvi_ct.enet2016=c() #Intialize to capture predindvi_cted values from ElasticNet model
cor.ndvi_ct.las2016=c() #Intialize to capture correlatio from LASSO model
cor.ndvi_ct.rdg2016=c() #Intialize to capture correlatio from Ridge model
cor.ndvi_ct.enet2016=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:10){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ndvi_ct.data2016 = ndvi_ct_data[ndvi_ct_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ndvi_ct.data2016 = ndvi_ct_data[ndvi_ct_data$trial!=trial,] #sets up predindvi_ction populations
        rndvi_ctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ndvi_ct.model2016 = train(GRYLD ~ ., data = train.ndvi_ct.data2016[,-c(1,2)], #Looks like we are fitting a lasso model.
                                        method = "lasso",
                                        trControl = rndvi_ctrl,
                                        preProc = c("center", "scale"))
        ridge.ndvi_ct.model2016 = train(GRYLD ~ ., data = train.ndvi_ct.data2016[,-c(1,2)], #looks like we are fitting a ridge model.
                                        method = "ridge",
                                        trControl = rndvi_ctrl,
                                        preProc = c("center", "scale"))
        enet.ndvi_ct.model2016 = train(GRYLD ~ ., data = train.ndvi_ct.data2016[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                       method = "enet",
                                       trControl = rndvi_ctrl,
                                       preProc = c("center", "scale"))
        #model coefficients
        coef.ndvi_ct.lasso2016 = lasso.ndvi_ct.model2016$finalModel #extrandvi_cting some coefficients from LASSO model
        coef.ndvi_ct.ridge2016 = ridge.ndvi_ct.model2016$finalModel #extrandvi_cting some coefficients from Ridge model
        coef.ndvi_ct.enet2016 = enet.ndvi_ct.model2016$finalModel #extrandvi_cting some coefficients from ElasticNet model
        #predindvi_ction
        pred.ndvi_ct.lasso2016 = predict(lasso.ndvi_ct.model2016,test.ndvi_ct.data2016[,-c(1,2)]) #making a predindvi_ct statement.  This looks promising
        pred.ndvi_ct.ridge2016 = predict(ridge.ndvi_ct.model2016,test.ndvi_ct.data2016[,-c(1,2)]) #appears the same with other models
        pred.ndvi_ct.elasticnet2016 = predict(enet.ndvi_ct.model2016,test.ndvi_ct.data2016[,-c(1,2)])
        #extrandvi_cting correlations
        cor.ndvi_ct.l2016=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,20],pred.ndvi_ct.lasso2016) # Correlaiton between predindvi_cted value from LASSO regression and BLUE
        cor.ndvi_ct.las2016=c(cor.ndvi_ct.las2016,cor.ndvi_ct.l2016) #writing correlation out to cor.ndvi_ct.las2016
        cor.ndvi_ct.r2016=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,20],pred.ndvi_ct.ridge2016) # Correlaiton between predindvi_cted value from Ridge regression and BLUE
        cor.ndvi_ct.rdg2016=c(cor.ndvi_ct.rdg2016,cor.ndvi_ct.r2016) #writing correlation out to cor.ndvi_ct.rdg2016
        cor.ndvi_ct.e2016=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,20],pred.ndvi_ct.elasticnet2016) # Correlaiton between predindvi_cted value from ElasticNet regression and BLUE
        cor.ndvi_ct.enet2016=c(cor.ndvi_ct.enet2016,cor.ndvi_ct.e2016) #writing correlation out to cor.ndvi_ct.enet2016
        #extrandvi_cting predindvi_cted values
        predindvi_ct.ndvi_ct.las2016 = c(predindvi_ct.ndvi_ct.las2016,pred.ndvi_ct.lasso2016) # Combinig all predindvi_cted values from LASSO
        predindvi_ct.ndvi_ct.rdg2016 = c(predindvi_ct.ndvi_ct.rdg2016,pred.ndvi_ct.ridge2016) # Combinig all predindvi_cted values from Ridge
        predindvi_ct.ndvi_ct.enet2016 = c(predindvi_ct.ndvi_ct.enet2016,pred.ndvi_ct.elasticnet2016) # Combinig all predindvi_cted values from ElasticNet
} #end for loop 

#print predindvi_ction accuracies
cor.ndvi_ct.step2016 #individual trial correlations for step regression
predaccuracy.ndvi_ct.step2016=mean(cor.ndvi_ct.step2016) #average correlation from stepwise model
cor.ndvi_ct.las2016 #individual trial correlations for LASSO regression
predaccuracy.ndvi_ct.lasso2016=mean(cor.ndvi_ct.las2016) #average correlation from LASSO model
cor.ndvi_ct.rdg2016 #individual trial correlations for Ridge regression
predaccuracy.ndvi_ct.ridge2016=mean(cor.ndvi_ct.rdg2016) #average correlation from Ridge model
cor.ndvi_ct.enet2016 #individual trial correlations for ElasticNet regression
predaccuracy.ndvi_ct.enet2016=mean(cor.ndvi_ct.enet2016) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ndvi_ct
plot(step.predyield.ndvi_ct.2016$GRYLD, step.predyield.ndvi_ct.2016$step.predvalue.ndvi_ct.2016, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="Stepwise regression with ndvi_ct 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.step2016 <-lm(step.predyield.ndvi_ct.2016$step.predvalue.ndvi_ct.2016 ~ step.predyield.ndvi_ct.2016$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.step2016, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.step2016 = vector('expression',2)
rp.ndvi_ct.step2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                    list(MYVALUE = format(summary(fit.ndvi_ct.step2016)$adj.r.squared,dig=2)))[2] # extrandvi_cting R2 value
rp.ndvi_ct.step2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.step2016, digits = 2)))[2] # r (predindvi_ction accuracy) value
legend("bottomright", bty="n",legend=rp.ndvi_ct.step2016)
R2.ndvi_ct.step2016=summary(fit.ndvi_ct.step2016)$adj.r.squared
R2.ndvi_ct.step2016=round(R2.ndvi_ct.step2016,digits=2)
r.ndvi_ct.step2016=round(predaccuracy.ndvi_ct.step2016,digits=2)
models="stepwise"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.step2016=data.frame(cbind(traits,models,r.ndvi_ct.step2016,R2.ndvi_ct.step2016))
colnames(pred.accuracy.ndvi_ct.step2016)[3:4]=c("r","R2")

### Plotting LASSO regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.las2016, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="LASSO regression with ndvi_ct 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.las2016 <-lm(predindvi_ct.ndvi_ct.las2016 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.las2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.las2016 = vector('expression',2)
rp.ndvi_ct.las2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(summary(fit.ndvi_ct.las2016)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.las2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                   list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.lasso2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.las2016)
R2.ndvi_ct.las2016=summary(fit.ndvi_ct.las2016)$adj.r.squared
R2.ndvi_ct.las2016=round(R2.ndvi_ct.las2016,digits=2)
r.ndvi_ct.las2016=round(predaccuracy.ndvi_ct.lasso2016,digits=2)
models="LASSO"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.las2016=data.frame(cbind(traits,models,r.ndvi_ct.las2016,R2.ndvi_ct.las2016))
colnames(pred.accuracy.ndvi_ct.las2016)[3:4]=c("r","R2")

### Plotting ridge regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.rdg2016, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="Ridge regression with ndvi_ct 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.rdg2016 <-lm(predindvi_ct.ndvi_ct.rdg2016 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.rdg2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.rdg2016 = vector('expression',2)
rp.ndvi_ct.rdg2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(summary(fit.ndvi_ct.rdg2016)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.rdg2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                   list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.ridge2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.rdg2016)
R2.ndvi_ct.rdg2016=summary(fit.ndvi_ct.rdg2016)$adj.r.squared
R2.ndvi_ct.rdg2016=round(R2.ndvi_ct.rdg2016,digits=2)
r.ndvi_ct.rdg2016=round(predaccuracy.ndvi_ct.ridge2016,digits=2)
models="ridge"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.rdg2016=data.frame(cbind(traits,models,r.ndvi_ct.rdg2016,R2.ndvi_ct.rdg2016))
colnames(pred.accuracy.ndvi_ct.rdg2016)[3:4]=c("r","R2")

### Plotting enet regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.enet2016, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="ElasticNet regression with ndvi_ct 2016",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.enet2016 <-lm(predindvi_ct.ndvi_ct.enet2016 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.enet2016, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.enet2016 = vector('expression',2)
rp.ndvi_ct.enet2016[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                    list(MYVALUE = format(summary(fit.ndvi_ct.enet2016)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.enet2016[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.enet2016, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.enet2016)
R2.ndvi_ct.enet2016=summary(fit.ndvi_ct.enet2016)$adj.r.squared
R2.ndvi_ct.enet2016=round(R2.ndvi_ct.enet2016,digits=2)
r.ndvi_ct.enet2016=round(predaccuracy.ndvi_ct.enet2016,digits=2)
models="elasticnet"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.enet2016=data.frame(cbind(traits,models,r.ndvi_ct.enet2016,R2.ndvi_ct.enet2016))
colnames(pred.accuracy.ndvi_ct.enet2016)[3:4]=c("r","R2")

multivar.ndvi_ct.prediction.accuracy2016=data.frame(rbind(pred.accuracy.ndvi_ct.step2016,pred.accuracy.ndvi_ct.las2016,pred.accuracy.ndvi_ct.rdg2016,pred.accuracy.ndvi_ct.enet2016))
# write.csv(multivar.ndvi_ct.prediction.accuracy2016,file="Multivariate_predaccuracy_ndvi_ct_2016.csv",row.names = FALSE,quote = FALSE)


### Multivariate model yield predindvi_ction 2017 with ndvi_ct
### Running stepwise regression model with ndvi_ct
BLUE2017=read.csv("BLUE_2017.csv")
ndvi_ct_data=BLUE2017[,c(1:30,39)]
stepwise.ndvi_ct.2017 <- lm(GRYLD ~ CT_20170104+CT_20170109+CT_20170114+CT_20170120+CT_20170125+CT_20170131+CT_20170205+CT_20170210+CT_20170215+CT_20170221+CT_20170225+CT_20170302+CT_20170307+CT_20170313+
                                    NDVI_20170103+NDVI_20170108+NDVI_20170114+NDVI_20170120+NDVI_20170125+NDVI_20170131+NDVI_20170205+NDVI_20170210+NDVI_20170215+NDVI_20170220+NDVI_20170225+NDVI_20170302+NDVI_20170307+NDVI_20170313
                            ,data=ndvi_ct_data)  ## writing the model with ndvi_ct
stepwise.ndvi_ct.2017 <- step(stepwise.ndvi_ct.2017, direndvi_ction = "both") ## running the stepwise regressing model
formula(stepwise.ndvi_ct.2017) ## getting the selendvi_cted variables from the stepwise regression model

### Predindvi_cted yield from stepwise regression model with ndvi_ct
step.predyield.ndvi_ct.2017=NULL 
for(i in 1:11){
        step.prediction.ndvi_ct.2017=ndvi_ct_data[ndvi_ct_data$trial==i,] #get just one trial to predindvi_ct
        step.training.ndvi_ct.2017=ndvi_ct_data[ndvi_ct_data$trial!=i,] #get all other trials to train on
        fit.ndvi_ct.2017 = lm(formula = formula(stepwise.ndvi_ct.2017), data =step.training.ndvi_ct.2017) ## direndvi_ctly import variables into the model from stepwise regression 
        step.predvalue.ndvi_ct.2017=predict(object=fit.ndvi_ct.2017, newdata=step.prediction.ndvi_ct.2017) ## Predindvi_cted value of individual trial
        predicted.ndvi_ct.trial2017=data.frame(entry=step.prediction.ndvi_ct.2017$plot,step.predvalue.ndvi_ct.2017) ## Predindvi_cted value of individual trial in data frame
        step.predyield.ndvi_ct.2017=rbind(step.predyield.ndvi_ct.2017, predicted.ndvi_ct.trial2017) ## Combinig all trial's predindvi_cted value into one frame
}
step.predyield.ndvi_ct.2017 = data.frame(trial=rep(1:11,each=60),step.predyield.ndvi_ct.2017)
step.predyield.ndvi_ct.2017 <- merge(step.predyield.ndvi_ct.2017, ndvi_ct_data[, c(1:2,31)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extrandvi_cting predindvi_ction accuracy from stepwise regression model with ndvi_ct
cor.ndvi_ct.step2017 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ndvi_ct.2017[step.predyield.ndvi_ct.2017$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ndvi_ct.2017, trial_data$GRYLD) #get correlations
        cor.ndvi_ct.step2017 <- c(cor.ndvi_ct.step2017, trial_cor) #write results outside of the for loop
}

### Predindvi_ction from shrinkage models
predindvi_ct.ndvi_ct.las2017=c()  #Intialize to capture predindvi_cted values from LASSO model
predindvi_ct.ndvi_ct.rdg2017=c() #Intialize to capture predindvi_cted values from Ridge model
predindvi_ct.ndvi_ct.enet2017=c() #Intialize to capture predindvi_cted values from ElasticNet model
cor.ndvi_ct.las2017=c() #Intialize to capture correlatio from LASSO model
cor.ndvi_ct.rdg2017=c() #Intialize to capture correlatio from Ridge model
cor.ndvi_ct.enet2017=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ndvi_ct.data2017 = ndvi_ct_data[ndvi_ct_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ndvi_ct.data2017 = ndvi_ct_data[ndvi_ct_data$trial!=trial,] #sets up predindvi_ction populations
        rndvi_ctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ndvi_ct.model2017 = train(GRYLD ~ ., data = train.ndvi_ct.data2017[,-c(1,2)], #Looks like we are fitting a lasso model.
                                        method = "lasso",
                                        trControl = rndvi_ctrl,
                                        preProc = c("center", "scale"))
        ridge.ndvi_ct.model2017 = train(GRYLD ~ ., data = train.ndvi_ct.data2017[,-c(1,2)], #looks like we are fitting a ridge model.
                                        method = "ridge",
                                        trControl = rndvi_ctrl,
                                        preProc = c("center", "scale"))
        enet.ndvi_ct.model2017 = train(GRYLD ~ ., data = train.ndvi_ct.data2017[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                       method = "enet",
                                       trControl = rndvi_ctrl,
                                       preProc = c("center", "scale"))
        #model coefficients
        coef.ndvi_ct.lasso2017 = lasso.ndvi_ct.model2017$finalModel #extrandvi_cting some coefficients from LASSO model
        coef.ndvi_ct.ridge2017 = ridge.ndvi_ct.model2017$finalModel #extrandvi_cting some coefficients from Ridge model
        coef.ndvi_ct.enet2017 = enet.ndvi_ct.model2017$finalModel #extrandvi_cting some coefficients from ElasticNet model
        #predindvi_ction
        pred.ndvi_ct.lasso2017 = predict(lasso.ndvi_ct.model2017,test.ndvi_ct.data2017[,-c(1,2)]) #making a predindvi_ct statement.  This looks promising
        pred.ndvi_ct.ridge2017 = predict(ridge.ndvi_ct.model2017,test.ndvi_ct.data2017[,-c(1,2)]) #appears the same with other models
        pred.ndvi_ct.elasticnet2017 = predict(enet.ndvi_ct.model2017,test.ndvi_ct.data2017[,-c(1,2)])
        #extrandvi_cting correlations
        cor.ndvi_ct.l2017=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,31],pred.ndvi_ct.lasso2017) # Correlaiton between predindvi_cted value from LASSO regression and BLUE
        cor.ndvi_ct.las2017=c(cor.ndvi_ct.las2017,cor.ndvi_ct.l2017) #writing correlation out to cor.ndvi_ct.las2017
        cor.ndvi_ct.r2017=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,31],pred.ndvi_ct.ridge2017) # Correlaiton between predindvi_cted value from Ridge regression and BLUE
        cor.ndvi_ct.rdg2017=c(cor.ndvi_ct.rdg2017,cor.ndvi_ct.r2017) #writing correlation out to cor.ndvi_ct.rdg2017
        cor.ndvi_ct.e2017=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,31],pred.ndvi_ct.elasticnet2017) # Correlaiton between predindvi_cted value from ElasticNet regression and BLUE
        cor.ndvi_ct.enet2017=c(cor.ndvi_ct.enet2017,cor.ndvi_ct.e2017) #writing correlation out to cor.ndvi_ct.enet2017
        #extrandvi_cting predindvi_cted values
        predindvi_ct.ndvi_ct.las2017 = c(predindvi_ct.ndvi_ct.las2017,pred.ndvi_ct.lasso2017) # Combinig all predindvi_cted values from LASSO
        predindvi_ct.ndvi_ct.rdg2017 = c(predindvi_ct.ndvi_ct.rdg2017,pred.ndvi_ct.ridge2017) # Combinig all predindvi_cted values from Ridge
        predindvi_ct.ndvi_ct.enet2017 = c(predindvi_ct.ndvi_ct.enet2017,pred.ndvi_ct.elasticnet2017) # Combinig all predindvi_cted values from ElasticNet
} #end for loop 

#print predindvi_ction accuracies
cor.ndvi_ct.step2017 #individual trial correlations for step regression
predaccuracy.ndvi_ct.step2017=mean(cor.ndvi_ct.step2017) #average correlation from stepwise model
cor.ndvi_ct.las2017 #individual trial correlations for LASSO regression
predaccuracy.ndvi_ct.lasso2017=mean(cor.ndvi_ct.las2017) #average correlation from LASSO model
cor.ndvi_ct.rdg2017 #individual trial correlations for Ridge regression
predaccuracy.ndvi_ct.ridge2017=mean(cor.ndvi_ct.rdg2017) #average correlation from Ridge model
cor.ndvi_ct.enet2017 #individual trial correlations for ElasticNet regression
predaccuracy.ndvi_ct.enet2017=mean(cor.ndvi_ct.enet2017) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ndvi_ct
plot(step.predyield.ndvi_ct.2017$GRYLD, step.predyield.ndvi_ct.2017$step.predvalue.ndvi_ct.2017, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="Stepwise regression with ndvi_ct 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.step2017 <-lm(step.predyield.ndvi_ct.2017$step.predvalue.ndvi_ct.2017 ~ step.predyield.ndvi_ct.2017$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.step2017, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.step2017 = vector('expression',2)
rp.ndvi_ct.step2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                    list(MYVALUE = format(summary(fit.ndvi_ct.step2017)$adj.r.squared,dig=2)))[2] # extrandvi_cting R2 value
rp.ndvi_ct.step2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.step2017, digits = 2)))[2] # r (predindvi_ction accuracy) value
legend("bottomright", bty="n",legend=rp.ndvi_ct.step2017)
R2.ndvi_ct.step2017=summary(fit.ndvi_ct.step2017)$adj.r.squared
R2.ndvi_ct.step2017=round(R2.ndvi_ct.step2017,digits=2)
r.ndvi_ct.step2017=round(predaccuracy.ndvi_ct.step2017,digits=2)
models="stepwise"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.step2017=data.frame(cbind(traits,models,r.ndvi_ct.step2017,R2.ndvi_ct.step2017))
colnames(pred.accuracy.ndvi_ct.step2017)[3:4]=c("r","R2")

### Plotting LASSO regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.las2017, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="LASSO regression with ndvi_ct 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.las2017 <-lm(predindvi_ct.ndvi_ct.las2017 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.las2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.las2017 = vector('expression',2)
rp.ndvi_ct.las2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(summary(fit.ndvi_ct.las2017)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.las2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                   list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.lasso2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.las2017)
R2.ndvi_ct.las2017=summary(fit.ndvi_ct.las2017)$adj.r.squared
R2.ndvi_ct.las2017=round(R2.ndvi_ct.las2017,digits=2)
r.ndvi_ct.las2017=round(predaccuracy.ndvi_ct.lasso2017,digits=2)
models="LASSO"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.las2017=data.frame(cbind(traits,models,r.ndvi_ct.las2017,R2.ndvi_ct.las2017))
colnames(pred.accuracy.ndvi_ct.las2017)[3:4]=c("r","R2")

### Plotting ridge regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.rdg2017, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="Ridge regression with ndvi_ct 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.rdg2017 <-lm(predindvi_ct.ndvi_ct.rdg2017 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.rdg2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.rdg2017 = vector('expression',2)
rp.ndvi_ct.rdg2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(summary(fit.ndvi_ct.rdg2017)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.rdg2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                   list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.ridge2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.rdg2017)
R2.ndvi_ct.rdg2017=summary(fit.ndvi_ct.rdg2017)$adj.r.squared
R2.ndvi_ct.rdg2017=round(R2.ndvi_ct.rdg2017,digits=2)
r.ndvi_ct.rdg2017=round(predaccuracy.ndvi_ct.ridge2017,digits=2)
models="ridge"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.rdg2017=data.frame(cbind(traits,models,r.ndvi_ct.rdg2017,R2.ndvi_ct.rdg2017))
colnames(pred.accuracy.ndvi_ct.rdg2017)[3:4]=c("r","R2")

### Plotting enet regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.enet2017, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="ElasticNet regression with ndvi_ct 2017",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.enet2017 <-lm(predindvi_ct.ndvi_ct.enet2017 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.enet2017, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.enet2017 = vector('expression',2)
rp.ndvi_ct.enet2017[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                    list(MYVALUE = format(summary(fit.ndvi_ct.enet2017)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.enet2017[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.enet2017, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.enet2017)
R2.ndvi_ct.enet2017=summary(fit.ndvi_ct.enet2017)$adj.r.squared
R2.ndvi_ct.enet2017=round(R2.ndvi_ct.enet2017,digits=2)
r.ndvi_ct.enet2017=round(predaccuracy.ndvi_ct.enet2017,digits=2)
models="elasticnet"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.enet2017=data.frame(cbind(traits,models,r.ndvi_ct.enet2017,R2.ndvi_ct.enet2017))
colnames(pred.accuracy.ndvi_ct.enet2017)[3:4]=c("r","R2")

multivar.ndvi_ct.prediction.accuracy2017=data.frame(rbind(pred.accuracy.ndvi_ct.step2017,pred.accuracy.ndvi_ct.las2017,pred.accuracy.ndvi_ct.rdg2017,pred.accuracy.ndvi_ct.enet2017))
# write.csv(multivar.ndvi_ct.prediction.accuracy2017,file="Multivariate_predaccuracy_ndvi_ct_2017.csv",row.names = FALSE,quote = FALSE)


### Multivariate model yield predindvi_ction 2018 with ndvi_ct
### Running stepwise regression model with ndvi_ct
BLUE2018=read.csv("BLUE_2018.csv")
ndvi_ct_data=BLUE2018[,c(1:26,35)]
stepwise.ndvi_ct.2018 <- lm(GRYLD ~ CT_20180126+CT_20180131+CT_20180205+CT_20180210+CT_20180214+CT_20180219+CT_20180225+CT_20180301+CT_20180305+CT_20180310+CT_20180315+CT_20180320+
                                    NDVI_20180126+NDVI_20180131+NDVI_20180204+NDVI_20180210+NDVI_20180214+NDVI_20180219+NDVI_20180226+NDVI_20180301+NDVI_20180305+NDVI_20180310+NDVI_20180315+NDVI_20180320
                            ,data=ndvi_ct_data)  ## writing the model with ndvi_ct
stepwise.ndvi_ct.2018 <- step(stepwise.ndvi_ct.2018, direndvi_ction = "both") ## running the stepwise regressing model
formula(stepwise.ndvi_ct.2018) ## getting the selendvi_cted variables from the stepwise regression model

### Predindvi_cted yield from stepwise regression model with ndvi_ct
step.predyield.ndvi_ct.2018=NULL 
for(i in 1:11){
        step.prediction.ndvi_ct.2018=ndvi_ct_data[ndvi_ct_data$trial==i,] #get just one trial to predindvi_ct
        step.training.ndvi_ct.2018=ndvi_ct_data[ndvi_ct_data$trial!=i,] #get all other trials to train on
        fit.ndvi_ct.2018 = lm(formula = formula(stepwise.ndvi_ct.2018), data =step.training.ndvi_ct.2018) ## direndvi_ctly import variables into the model from stepwise regression 
        step.predvalue.ndvi_ct.2018=predict(object=fit.ndvi_ct.2018, newdata=step.prediction.ndvi_ct.2018) ## Predindvi_cted value of individual trial
        predicted.ndvi_ct.trial2018=data.frame(entry=step.prediction.ndvi_ct.2018$plot,step.predvalue.ndvi_ct.2018) ## Predindvi_cted value of individual trial in data frame
        step.predyield.ndvi_ct.2018=rbind(step.predyield.ndvi_ct.2018, predicted.ndvi_ct.trial2018) ## Combinig all trial's predindvi_cted value into one frame
}
step.predyield.ndvi_ct.2018 = data.frame(trial=rep(1:11,each=60),step.predyield.ndvi_ct.2018)
step.predyield.ndvi_ct.2018 <- merge(step.predyield.ndvi_ct.2018, ndvi_ct_data[, c(1:2,27)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extrandvi_cting predindvi_ction accuracy from stepwise regression model with ndvi_ct
cor.ndvi_ct.step2018 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ndvi_ct.2018[step.predyield.ndvi_ct.2018$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ndvi_ct.2018, trial_data$GRYLD) #get correlations
        cor.ndvi_ct.step2018 <- c(cor.ndvi_ct.step2018, trial_cor) #write results outside of the for loop
}

### Predindvi_ction from shrinkage models
predindvi_ct.ndvi_ct.las2018=c()  #Intialize to capture predindvi_cted values from LASSO model
predindvi_ct.ndvi_ct.rdg2018=c() #Intialize to capture predindvi_cted values from Ridge model
predindvi_ct.ndvi_ct.enet2018=c() #Intialize to capture predindvi_cted values from ElasticNet model
cor.ndvi_ct.las2018=c() #Intialize to capture correlatio from LASSO model
cor.ndvi_ct.rdg2018=c() #Intialize to capture correlatio from Ridge model
cor.ndvi_ct.enet2018=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ndvi_ct.data2018 = ndvi_ct_data[ndvi_ct_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ndvi_ct.data2018 = ndvi_ct_data[ndvi_ct_data$trial!=trial,] #sets up predindvi_ction populations
        rndvi_ctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ndvi_ct.model2018 = train(GRYLD ~ ., data = train.ndvi_ct.data2018[,-c(1,2)], #Looks like we are fitting a lasso model.
                                        method = "lasso",
                                        trControl = rndvi_ctrl,
                                        preProc = c("center", "scale"))
        ridge.ndvi_ct.model2018 = train(GRYLD ~ ., data = train.ndvi_ct.data2018[,-c(1,2)], #looks like we are fitting a ridge model.
                                        method = "ridge",
                                        trControl = rndvi_ctrl,
                                        preProc = c("center", "scale"))
        enet.ndvi_ct.model2018 = train(GRYLD ~ ., data = train.ndvi_ct.data2018[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                       method = "enet",
                                       trControl = rndvi_ctrl,
                                       preProc = c("center", "scale"))
        #model coefficients
        coef.ndvi_ct.lasso2018 = lasso.ndvi_ct.model2018$finalModel #extrandvi_cting some coefficients from LASSO model
        coef.ndvi_ct.ridge2018 = ridge.ndvi_ct.model2018$finalModel #extrandvi_cting some coefficients from Ridge model
        coef.ndvi_ct.enet2018 = enet.ndvi_ct.model2018$finalModel #extrandvi_cting some coefficients from ElasticNet model
        #predindvi_ction
        pred.ndvi_ct.lasso2018 = predict(lasso.ndvi_ct.model2018,test.ndvi_ct.data2018[,-c(1,2)]) #making a predindvi_ct statement.  This looks promising
        pred.ndvi_ct.ridge2018 = predict(ridge.ndvi_ct.model2018,test.ndvi_ct.data2018[,-c(1,2)]) #appears the same with other models
        pred.ndvi_ct.elasticnet2018 = predict(enet.ndvi_ct.model2018,test.ndvi_ct.data2018[,-c(1,2)])
        #extrandvi_cting correlations
        cor.ndvi_ct.l2018=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,27],pred.ndvi_ct.lasso2018) # Correlaiton between predindvi_cted value from LASSO regression and BLUE
        cor.ndvi_ct.las2018=c(cor.ndvi_ct.las2018,cor.ndvi_ct.l2018) #writing correlation out to cor.ndvi_ct.las2018
        cor.ndvi_ct.r2018=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,27],pred.ndvi_ct.ridge2018) # Correlaiton between predindvi_cted value from Ridge regression and BLUE
        cor.ndvi_ct.rdg2018=c(cor.ndvi_ct.rdg2018,cor.ndvi_ct.r2018) #writing correlation out to cor.ndvi_ct.rdg2018
        cor.ndvi_ct.e2018=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,27],pred.ndvi_ct.elasticnet2018) # Correlaiton between predindvi_cted value from ElasticNet regression and BLUE
        cor.ndvi_ct.enet2018=c(cor.ndvi_ct.enet2018,cor.ndvi_ct.e2018) #writing correlation out to cor.ndvi_ct.enet2018
        #extrandvi_cting predindvi_cted values
        predindvi_ct.ndvi_ct.las2018 = c(predindvi_ct.ndvi_ct.las2018,pred.ndvi_ct.lasso2018) # Combinig all predindvi_cted values from LASSO
        predindvi_ct.ndvi_ct.rdg2018 = c(predindvi_ct.ndvi_ct.rdg2018,pred.ndvi_ct.ridge2018) # Combinig all predindvi_cted values from Ridge
        predindvi_ct.ndvi_ct.enet2018 = c(predindvi_ct.ndvi_ct.enet2018,pred.ndvi_ct.elasticnet2018) # Combinig all predindvi_cted values from ElasticNet
} #end for loop 

#print predindvi_ction accuracies
cor.ndvi_ct.step2018 #individual trial correlations for step regression
predaccuracy.ndvi_ct.step2018=mean(cor.ndvi_ct.step2018) #average correlation from stepwise model
cor.ndvi_ct.las2018 #individual trial correlations for LASSO regression
predaccuracy.ndvi_ct.lasso2018=mean(cor.ndvi_ct.las2018) #average correlation from LASSO model
cor.ndvi_ct.rdg2018 #individual trial correlations for Ridge regression
predaccuracy.ndvi_ct.ridge2018=mean(cor.ndvi_ct.rdg2018) #average correlation from Ridge model
cor.ndvi_ct.enet2018 #individual trial correlations for ElasticNet regression
predaccuracy.ndvi_ct.enet2018=mean(cor.ndvi_ct.enet2018) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ndvi_ct
plot(step.predyield.ndvi_ct.2018$GRYLD, step.predyield.ndvi_ct.2018$step.predvalue.ndvi_ct.2018, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="Stepwise regression with ndvi_ct 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.step2018 <-lm(step.predyield.ndvi_ct.2018$step.predvalue.ndvi_ct.2018 ~ step.predyield.ndvi_ct.2018$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.step2018, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.step2018 = vector('expression',2)
rp.ndvi_ct.step2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                    list(MYVALUE = format(summary(fit.ndvi_ct.step2018)$adj.r.squared,dig=2)))[2] # extrandvi_cting R2 value
rp.ndvi_ct.step2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.step2018, digits = 2)))[2] # r (predindvi_ction accuracy) value
legend("bottomright", bty="n",legend=rp.ndvi_ct.step2018)
R2.ndvi_ct.step2018=summary(fit.ndvi_ct.step2018)$adj.r.squared
R2.ndvi_ct.step2018=round(R2.ndvi_ct.step2018,digits=2)
r.ndvi_ct.step2018=round(predaccuracy.ndvi_ct.step2018,digits=2)
models="stepwise"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.step2018=data.frame(cbind(traits,models,r.ndvi_ct.step2018,R2.ndvi_ct.step2018))
colnames(pred.accuracy.ndvi_ct.step2018)[3:4]=c("r","R2")

### Plotting LASSO regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.las2018, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="LASSO regression with ndvi_ct 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.las2018 <-lm(predindvi_ct.ndvi_ct.las2018 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.las2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.las2018 = vector('expression',2)
rp.ndvi_ct.las2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(summary(fit.ndvi_ct.las2018)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.las2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                   list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.lasso2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.las2018)
R2.ndvi_ct.las2018=summary(fit.ndvi_ct.las2018)$adj.r.squared
R2.ndvi_ct.las2018=round(R2.ndvi_ct.las2018,digits=2)
r.ndvi_ct.las2018=round(predaccuracy.ndvi_ct.lasso2018,digits=2)
models="LASSO"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.las2018=data.frame(cbind(traits,models,r.ndvi_ct.las2018,R2.ndvi_ct.las2018))
colnames(pred.accuracy.ndvi_ct.las2018)[3:4]=c("r","R2")

### Plotting ridge regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.rdg2018, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="Ridge regression with ndvi_ct 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.rdg2018 <-lm(predindvi_ct.ndvi_ct.rdg2018 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.rdg2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.rdg2018 = vector('expression',2)
rp.ndvi_ct.rdg2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(summary(fit.ndvi_ct.rdg2018)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.rdg2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                   list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.ridge2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.rdg2018)
R2.ndvi_ct.rdg2018=summary(fit.ndvi_ct.rdg2018)$adj.r.squared
R2.ndvi_ct.rdg2018=round(R2.ndvi_ct.rdg2018,digits=2)
r.ndvi_ct.rdg2018=round(predaccuracy.ndvi_ct.ridge2018,digits=2)
models="ridge"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.rdg2018=data.frame(cbind(traits,models,r.ndvi_ct.rdg2018,R2.ndvi_ct.rdg2018))
colnames(pred.accuracy.ndvi_ct.rdg2018)[3:4]=c("r","R2")

### Plotting enet regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.enet2018, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="ElasticNet regression with ndvi_ct 2018",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.enet2018 <-lm(predindvi_ct.ndvi_ct.enet2018 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.enet2018, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.enet2018 = vector('expression',2)
rp.ndvi_ct.enet2018[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                    list(MYVALUE = format(summary(fit.ndvi_ct.enet2018)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.enet2018[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.enet2018, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.enet2018)
R2.ndvi_ct.enet2018=summary(fit.ndvi_ct.enet2018)$adj.r.squared
R2.ndvi_ct.enet2018=round(R2.ndvi_ct.enet2018,digits=2)
r.ndvi_ct.enet2018=round(predaccuracy.ndvi_ct.enet2018,digits=2)
models="elasticnet"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.enet2018=data.frame(cbind(traits,models,r.ndvi_ct.enet2018,R2.ndvi_ct.enet2018))
colnames(pred.accuracy.ndvi_ct.enet2018)[3:4]=c("r","R2")

multivar.ndvi_ct.prediction.accuracy2018=data.frame(rbind(pred.accuracy.ndvi_ct.step2018,pred.accuracy.ndvi_ct.las2018,pred.accuracy.ndvi_ct.rdg2018,pred.accuracy.ndvi_ct.enet2018))
# write.csv(multivar.ndvi_ct.prediction.accuracy2018,file="Multivariate_predaccuracy_ndvi_ct_2018.csv",row.names = FALSE,quote = FALSE)


### Multivariate model yield predindvi_ction 2019 with ndvi_ct
### Running stepwise regression model with ndvi_ct
BLUE2019=read.csv("BLUE_2019.csv")
ndvi_ct_data=BLUE2019[,c(1:28,37)]
stepwise.ndvi_ct.2019 <- lm(GRYLD ~ CT_20190123+CT_20190127+CT_20190131+CT_20190205+CT_20190211+CT_20190218+CT_20190223+CT_20190301+CT_20190305+CT_20190311+CT_20190316+CT_20190320+CT_20190325+
                                    NDVI_20190121+NDVI_20190127+NDVI_20190131+NDVI_20190205+NDVI_20190211+NDVI_20190218+NDVI_20190222+NDVI_20190228+NDVI_20190305+NDVI_20190311+NDVI_20190315+NDVI_20190320+NDVI_20190325
                            ,data=ndvi_ct_data)  ## writing the model with ndvi_ct
stepwise.ndvi_ct.2019 <- step(stepwise.ndvi_ct.2019, direndvi_ction = "both") ## running the stepwise regressing model
formula(stepwise.ndvi_ct.2019) ## getting the selendvi_cted variables from the stepwise regression model

### Predindvi_cted yield from stepwise regression model with ndvi_ct
step.predyield.ndvi_ct.2019=NULL 
for(i in 1:10){
        step.prediction.ndvi_ct.2019=ndvi_ct_data[ndvi_ct_data$trial==i,] #get just one trial to predindvi_ct
        step.training.ndvi_ct.2019=ndvi_ct_data[ndvi_ct_data$trial!=i,] #get all other trials to train on
        fit.ndvi_ct.2019 = lm(formula = formula(stepwise.ndvi_ct.2019), data =step.training.ndvi_ct.2019) ## direndvi_ctly import variables into the model from stepwise regression 
        step.predvalue.ndvi_ct.2019=predict(object=fit.ndvi_ct.2019, newdata=step.prediction.ndvi_ct.2019) ## Predindvi_cted value of individual trial
        predicted.ndvi_ct.trial2019=data.frame(entry=step.prediction.ndvi_ct.2019$plot,step.predvalue.ndvi_ct.2019) ## Predindvi_cted value of individual trial in data frame
        step.predyield.ndvi_ct.2019=rbind(step.predyield.ndvi_ct.2019, predicted.ndvi_ct.trial2019) ## Combinig all trial's predindvi_cted value into one frame
}
step.predyield.ndvi_ct.2019 = data.frame(trial=rep(1:10,each=60),step.predyield.ndvi_ct.2019)
step.predyield.ndvi_ct.2019 <- merge(step.predyield.ndvi_ct.2019, ndvi_ct_data[, c(1:2,29)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extrandvi_cting predindvi_ction accuracy from stepwise regression model with ndvi_ct
cor.ndvi_ct.step2019 <- NULL #start dataframe to hold correlation values
for(i in 1:10){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ndvi_ct.2019[step.predyield.ndvi_ct.2019$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ndvi_ct.2019, trial_data$GRYLD) #get correlations
        cor.ndvi_ct.step2019 <- c(cor.ndvi_ct.step2019, trial_cor) #write results outside of the for loop
}

### Predindvi_ction from shrinkage models
predindvi_ct.ndvi_ct.las2019=c()  #Intialize to capture predindvi_cted values from LASSO model
predindvi_ct.ndvi_ct.rdg2019=c() #Intialize to capture predindvi_cted values from Ridge model
predindvi_ct.ndvi_ct.enet2019=c() #Intialize to capture predindvi_cted values from ElasticNet model
cor.ndvi_ct.las2019=c() #Intialize to capture correlatio from LASSO model
cor.ndvi_ct.rdg2019=c() #Intialize to capture correlatio from Ridge model
cor.ndvi_ct.enet2019=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:10){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ndvi_ct.data2019 = ndvi_ct_data[ndvi_ct_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ndvi_ct.data2019 = ndvi_ct_data[ndvi_ct_data$trial!=trial,] #sets up predindvi_ction populations
        rndvi_ctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ndvi_ct.model2019 = train(GRYLD ~ ., data = train.ndvi_ct.data2019[,-c(1,2)], #Looks like we are fitting a lasso model.
                                        method = "lasso",
                                        trControl = rndvi_ctrl,
                                        preProc = c("center", "scale"))
        ridge.ndvi_ct.model2019 = train(GRYLD ~ ., data = train.ndvi_ct.data2019[,-c(1,2)], #looks like we are fitting a ridge model.
                                        method = "ridge",
                                        trControl = rndvi_ctrl,
                                        preProc = c("center", "scale"))
        enet.ndvi_ct.model2019 = train(GRYLD ~ ., data = train.ndvi_ct.data2019[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                       method = "enet",
                                       trControl = rndvi_ctrl,
                                       preProc = c("center", "scale"))
        #model coefficients
        coef.ndvi_ct.lasso2019 = lasso.ndvi_ct.model2019$finalModel #extrandvi_cting some coefficients from LASSO model
        coef.ndvi_ct.ridge2019 = ridge.ndvi_ct.model2019$finalModel #extrandvi_cting some coefficients from Ridge model
        coef.ndvi_ct.enet2019 = enet.ndvi_ct.model2019$finalModel #extrandvi_cting some coefficients from ElasticNet model
        #predindvi_ction
        pred.ndvi_ct.lasso2019 = predict(lasso.ndvi_ct.model2019,test.ndvi_ct.data2019[,-c(1,2)]) #making a predindvi_ct statement.  This looks promising
        pred.ndvi_ct.ridge2019 = predict(ridge.ndvi_ct.model2019,test.ndvi_ct.data2019[,-c(1,2)]) #appears the same with other models
        pred.ndvi_ct.elasticnet2019 = predict(enet.ndvi_ct.model2019,test.ndvi_ct.data2019[,-c(1,2)])
        #extrandvi_cting correlations
        cor.ndvi_ct.l2019=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,29],pred.ndvi_ct.lasso2019) # Correlaiton between predindvi_cted value from LASSO regression and BLUE
        cor.ndvi_ct.las2019=c(cor.ndvi_ct.las2019,cor.ndvi_ct.l2019) #writing correlation out to cor.ndvi_ct.las2019
        cor.ndvi_ct.r2019=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,29],pred.ndvi_ct.ridge2019) # Correlaiton between predindvi_cted value from Ridge regression and BLUE
        cor.ndvi_ct.rdg2019=c(cor.ndvi_ct.rdg2019,cor.ndvi_ct.r2019) #writing correlation out to cor.ndvi_ct.rdg2019
        cor.ndvi_ct.e2019=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,29],pred.ndvi_ct.elasticnet2019) # Correlaiton between predindvi_cted value from ElasticNet regression and BLUE
        cor.ndvi_ct.enet2019=c(cor.ndvi_ct.enet2019,cor.ndvi_ct.e2019) #writing correlation out to cor.ndvi_ct.enet2019
        #extrandvi_cting predindvi_cted values
        predindvi_ct.ndvi_ct.las2019 = c(predindvi_ct.ndvi_ct.las2019,pred.ndvi_ct.lasso2019) # Combinig all predindvi_cted values from LASSO
        predindvi_ct.ndvi_ct.rdg2019 = c(predindvi_ct.ndvi_ct.rdg2019,pred.ndvi_ct.ridge2019) # Combinig all predindvi_cted values from Ridge
        predindvi_ct.ndvi_ct.enet2019 = c(predindvi_ct.ndvi_ct.enet2019,pred.ndvi_ct.elasticnet2019) # Combinig all predindvi_cted values from ElasticNet
} #end for loop 

#print predindvi_ction accuracies
cor.ndvi_ct.step2019 #individual trial correlations for step regression
predaccuracy.ndvi_ct.step2019=mean(cor.ndvi_ct.step2019) #average correlation from stepwise model
cor.ndvi_ct.las2019 #individual trial correlations for LASSO regression
predaccuracy.ndvi_ct.lasso2019=mean(cor.ndvi_ct.las2019) #average correlation from LASSO model
cor.ndvi_ct.rdg2019 #individual trial correlations for Ridge regression
predaccuracy.ndvi_ct.ridge2019=mean(cor.ndvi_ct.rdg2019) #average correlation from Ridge model
cor.ndvi_ct.enet2019 #individual trial correlations for ElasticNet regression
predaccuracy.ndvi_ct.enet2019=mean(cor.ndvi_ct.enet2019) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ndvi_ct
plot(step.predyield.ndvi_ct.2019$GRYLD, step.predyield.ndvi_ct.2019$step.predvalue.ndvi_ct.2019, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="Stepwise regression with ndvi_ct 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.step2019 <-lm(step.predyield.ndvi_ct.2019$step.predvalue.ndvi_ct.2019 ~ step.predyield.ndvi_ct.2019$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.step2019, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.step2019 = vector('expression',2)
rp.ndvi_ct.step2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                    list(MYVALUE = format(summary(fit.ndvi_ct.step2019)$adj.r.squared,dig=2)))[2] # extrandvi_cting R2 value
rp.ndvi_ct.step2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.step2019, digits = 2)))[2] # r (predindvi_ction accuracy) value
legend("bottomright", bty="n",legend=rp.ndvi_ct.step2019)
R2.ndvi_ct.step2019=summary(fit.ndvi_ct.step2019)$adj.r.squared
R2.ndvi_ct.step2019=round(R2.ndvi_ct.step2019,digits=2)
r.ndvi_ct.step2019=round(predaccuracy.ndvi_ct.step2019,digits=2)
models="stepwise"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.step2019=data.frame(cbind(traits,models,r.ndvi_ct.step2019,R2.ndvi_ct.step2019))
colnames(pred.accuracy.ndvi_ct.step2019)[3:4]=c("r","R2")

### Plotting LASSO regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.las2019, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="LASSO regression with ndvi_ct 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.las2019 <-lm(predindvi_ct.ndvi_ct.las2019 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.las2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.las2019 = vector('expression',2)
rp.ndvi_ct.las2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(summary(fit.ndvi_ct.las2019)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.las2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                   list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.lasso2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.las2019)
R2.ndvi_ct.las2019=summary(fit.ndvi_ct.las2019)$adj.r.squared
R2.ndvi_ct.las2019=round(R2.ndvi_ct.las2019,digits=2)
r.ndvi_ct.las2019=round(predaccuracy.ndvi_ct.lasso2019,digits=2)
models="LASSO"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.las2019=data.frame(cbind(traits,models,r.ndvi_ct.las2019,R2.ndvi_ct.las2019))
colnames(pred.accuracy.ndvi_ct.las2019)[3:4]=c("r","R2")

### Plotting ridge regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.rdg2019, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="Ridge regression with ndvi_ct 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.rdg2019 <-lm(predindvi_ct.ndvi_ct.rdg2019 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.rdg2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.rdg2019 = vector('expression',2)
rp.ndvi_ct.rdg2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(summary(fit.ndvi_ct.rdg2019)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.rdg2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                   list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.ridge2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.rdg2019)
R2.ndvi_ct.rdg2019=summary(fit.ndvi_ct.rdg2019)$adj.r.squared
R2.ndvi_ct.rdg2019=round(R2.ndvi_ct.rdg2019,digits=2)
r.ndvi_ct.rdg2019=round(predaccuracy.ndvi_ct.ridge2019,digits=2)
models="ridge"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.rdg2019=data.frame(cbind(traits,models,r.ndvi_ct.rdg2019,R2.ndvi_ct.rdg2019))
colnames(pred.accuracy.ndvi_ct.rdg2019)[3:4]=c("r","R2")

### Plotting enet regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.enet2019, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="ElasticNet regression with ndvi_ct 2019",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.enet2019 <-lm(predindvi_ct.ndvi_ct.enet2019 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.enet2019, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.enet2019 = vector('expression',2)
rp.ndvi_ct.enet2019[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                    list(MYVALUE = format(summary(fit.ndvi_ct.enet2019)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.enet2019[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.enet2019, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.enet2019)
R2.ndvi_ct.enet2019=summary(fit.ndvi_ct.enet2019)$adj.r.squared
R2.ndvi_ct.enet2019=round(R2.ndvi_ct.enet2019,digits=2)
r.ndvi_ct.enet2019=round(predaccuracy.ndvi_ct.enet2019,digits=2)
models="elasticnet"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.enet2019=data.frame(cbind(traits,models,r.ndvi_ct.enet2019,R2.ndvi_ct.enet2019))
colnames(pred.accuracy.ndvi_ct.enet2019)[3:4]=c("r","R2")

multivar.ndvi_ct.prediction.accuracy2019=data.frame(rbind(pred.accuracy.ndvi_ct.step2019,pred.accuracy.ndvi_ct.las2019,pred.accuracy.ndvi_ct.rdg2019,pred.accuracy.ndvi_ct.enet2019))
# write.csv(multivar.ndvi_ct.prediction.accuracy2019,file="Multivariate_predaccuracy_ndvi_ct_2019.csv",row.names = FALSE,quote = FALSE)


### Multivariate model yield predindvi_ction 2020 with ndvi_ct
### Running stepwise regression model with ndvi_ct
BLUE2020=read.csv("BLUE_2020.csv")
ndvi_ct_data=BLUE2020[,c(1:32,44)]
stepwise.ndvi_ct.2020 <- lm(GRYLD ~ CT_20200112+CT_20200116+CT_20200121+CT_20200126+CT_20200130+CT_20200205+CT_20200210+CT_20200215+CT_20200220+CT_20200226+CT_20200302+CT_20200308+CT_20200313+CT_20200318+CT_20200323+
                                    NDVI_20200112+NDVI_20200116+NDVI_20200121+NDVI_20200126+NDVI_20200130+NDVI_20200205+NDVI_20200210+NDVI_20200215+NDVI_20200220+NDVI_20200226+NDVI_20200302+NDVI_20200308+NDVI_20200313+NDVI_20200318+NDVI_20200323
                            ,data=ndvi_ct_data)  ## writing the model with ndvi_ct
stepwise.ndvi_ct.2020 <- step(stepwise.ndvi_ct.2020, direndvi_ction = "both") ## running the stepwise regressing model
formula(stepwise.ndvi_ct.2020) ## getting the selendvi_cted variables from the stepwise regression model

### Predindvi_cted yield from stepwise regression model with ndvi_ct
step.predyield.ndvi_ct.2020=NULL 
for(i in 1:11){
        step.prediction.ndvi_ct.2020=ndvi_ct_data[ndvi_ct_data$trial==i,] #get just one trial to predindvi_ct
        step.training.ndvi_ct.2020=ndvi_ct_data[ndvi_ct_data$trial!=i,] #get all other trials to train on
        fit.ndvi_ct.2020 = lm(formula = formula(stepwise.ndvi_ct.2020), data =step.training.ndvi_ct.2020) ## direndvi_ctly import variables into the model from stepwise regression 
        step.predvalue.ndvi_ct.2020=predict(object=fit.ndvi_ct.2020, newdata=step.prediction.ndvi_ct.2020) ## Predindvi_cted value of individual trial
        predicted.ndvi_ct.trial2020=data.frame(entry=step.prediction.ndvi_ct.2020$plot,step.predvalue.ndvi_ct.2020) ## Predindvi_cted value of individual trial in data frame
        step.predyield.ndvi_ct.2020=rbind(step.predyield.ndvi_ct.2020, predicted.ndvi_ct.trial2020) ## Combinig all trial's predindvi_cted value into one frame
}
step.predyield.ndvi_ct.2020 = data.frame(trial=rep(1:11,each=60),step.predyield.ndvi_ct.2020)
step.predyield.ndvi_ct.2020 <- merge(step.predyield.ndvi_ct.2020, ndvi_ct_data[, c(1:2,33)], by.x = c('trial', 'entry'), by.y = c('trial', 'plot')) # added observed values

### extrandvi_cting predindvi_ction accuracy from stepwise regression model with ndvi_ct
cor.ndvi_ct.step2020 <- NULL #start dataframe to hold correlation values
for(i in 1:11){ #loop over the 10 trials for 10 fold cross validation
        trial_data <- step.predyield.ndvi_ct.2020[step.predyield.ndvi_ct.2020$trial == i, ] #get data
        trial_cor <- cor(trial_data$step.predvalue.ndvi_ct.2020, trial_data$GRYLD) #get correlations
        cor.ndvi_ct.step2020 <- c(cor.ndvi_ct.step2020, trial_cor) #write results outside of the for loop
}

### Predindvi_ction from shrinkage models
predindvi_ct.ndvi_ct.las2020=c()  #Intialize to capture predindvi_cted values from LASSO model
predindvi_ct.ndvi_ct.rdg2020=c() #Intialize to capture predindvi_cted values from Ridge model
predindvi_ct.ndvi_ct.enet2020=c() #Intialize to capture predindvi_cted values from ElasticNet model
cor.ndvi_ct.las2020=c() #Intialize to capture correlatio from LASSO model
cor.ndvi_ct.rdg2020=c() #Intialize to capture correlatio from Ridge model
cor.ndvi_ct.enet2020=c() #Intialize to capture correlatio from ElasticNet model
for (trial in 1:11){ #this starts a for loop over the values 1:10, must be looping over the trials
        test.ndvi_ct.data2020 = ndvi_ct_data[ndvi_ct_data$trial==trial,] #this appears to be more subsetting for model preparation (cross fold validation) test population
        train.ndvi_ct.data2020 = ndvi_ct_data[ndvi_ct_data$trial!=trial,] #sets up predindvi_ction populations
        rndvi_ctrl <- trainControl(method = "CV")
        #different models: lasso, ridge, elastic net
        lasso.ndvi_ct.model2020 = train(GRYLD ~ ., data = train.ndvi_ct.data2020[,-c(1,2)], #Looks like we are fitting a lasso model.
                                        method = "lasso",
                                        trControl = rndvi_ctrl,
                                        preProc = c("center", "scale"))
        ridge.ndvi_ct.model2020 = train(GRYLD ~ ., data = train.ndvi_ct.data2020[,-c(1,2)], #looks like we are fitting a ridge model.
                                        method = "ridge",
                                        trControl = rndvi_ctrl,
                                        preProc = c("center", "scale"))
        enet.ndvi_ct.model2020 = train(GRYLD ~ ., data = train.ndvi_ct.data2020[,-c(1,2)], #looks like we are fitting a elasticNet model.
                                       method = "enet",
                                       trControl = rndvi_ctrl,
                                       preProc = c("center", "scale"))
        #model coefficients
        coef.ndvi_ct.lasso2020 = lasso.ndvi_ct.model2020$finalModel #extrandvi_cting some coefficients from LASSO model
        coef.ndvi_ct.ridge2020 = ridge.ndvi_ct.model2020$finalModel #extrandvi_cting some coefficients from Ridge model
        coef.ndvi_ct.enet2020 = enet.ndvi_ct.model2020$finalModel #extrandvi_cting some coefficients from ElasticNet model
        #predindvi_ction
        pred.ndvi_ct.lasso2020 = predict(lasso.ndvi_ct.model2020,test.ndvi_ct.data2020[,-c(1,2)]) #making a predindvi_ct statement.  This looks promising
        pred.ndvi_ct.ridge2020 = predict(ridge.ndvi_ct.model2020,test.ndvi_ct.data2020[,-c(1,2)]) #appears the same with other models
        pred.ndvi_ct.elasticnet2020 = predict(enet.ndvi_ct.model2020,test.ndvi_ct.data2020[,-c(1,2)])
        #extrandvi_cting correlations
        cor.ndvi_ct.l2020=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,33],pred.ndvi_ct.lasso2020) # Correlaiton between predindvi_cted value from LASSO regression and BLUE
        cor.ndvi_ct.las2020=c(cor.ndvi_ct.las2020,cor.ndvi_ct.l2020) #writing correlation out to cor.ndvi_ct.las2020
        cor.ndvi_ct.r2020=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,33],pred.ndvi_ct.ridge2020) # Correlaiton between predindvi_cted value from Ridge regression and BLUE
        cor.ndvi_ct.rdg2020=c(cor.ndvi_ct.rdg2020,cor.ndvi_ct.r2020) #writing correlation out to cor.ndvi_ct.rdg2020
        cor.ndvi_ct.e2020=cor(ndvi_ct_data[ndvi_ct_data$trial==trial,33],pred.ndvi_ct.elasticnet2020) # Correlaiton between predindvi_cted value from ElasticNet regression and BLUE
        cor.ndvi_ct.enet2020=c(cor.ndvi_ct.enet2020,cor.ndvi_ct.e2020) #writing correlation out to cor.ndvi_ct.enet2020
        #extrandvi_cting predindvi_cted values
        predindvi_ct.ndvi_ct.las2020 = c(predindvi_ct.ndvi_ct.las2020,pred.ndvi_ct.lasso2020) # Combinig all predindvi_cted values from LASSO
        predindvi_ct.ndvi_ct.rdg2020 = c(predindvi_ct.ndvi_ct.rdg2020,pred.ndvi_ct.ridge2020) # Combinig all predindvi_cted values from Ridge
        predindvi_ct.ndvi_ct.enet2020 = c(predindvi_ct.ndvi_ct.enet2020,pred.ndvi_ct.elasticnet2020) # Combinig all predindvi_cted values from ElasticNet
} #end for loop 

#print predindvi_ction accuracies
cor.ndvi_ct.step2020 #individual trial correlations for step regression
predaccuracy.ndvi_ct.step2020=mean(cor.ndvi_ct.step2020) #average correlation from stepwise model
cor.ndvi_ct.las2020 #individual trial correlations for LASSO regression
predaccuracy.ndvi_ct.lasso2020=mean(cor.ndvi_ct.las2020) #average correlation from LASSO model
cor.ndvi_ct.rdg2020 #individual trial correlations for Ridge regression
predaccuracy.ndvi_ct.ridge2020=mean(cor.ndvi_ct.rdg2020) #average correlation from Ridge model
cor.ndvi_ct.enet2020 #individual trial correlations for ElasticNet regression
predaccuracy.ndvi_ct.enet2020=mean(cor.ndvi_ct.enet2020) #average correlation from ElasticNet model

### Plotting different models

### Plotting stepwise regression model with ndvi_ct
plot(step.predyield.ndvi_ct.2020$GRYLD, step.predyield.ndvi_ct.2020$step.predvalue.ndvi_ct.2020, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="Stepwise regression with ndvi_ct 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.step2020 <-lm(step.predyield.ndvi_ct.2020$step.predvalue.ndvi_ct.2020 ~ step.predyield.ndvi_ct.2020$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.step2020, col = 'blue') #add trendline to graph (abline)
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.step2020 = vector('expression',2)
rp.ndvi_ct.step2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                    list(MYVALUE = format(summary(fit.ndvi_ct.step2020)$adj.r.squared,dig=2)))[2] # extrandvi_cting R2 value
rp.ndvi_ct.step2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.step2020, digits = 2)))[2] # r (predindvi_ction accuracy) value
legend("bottomright", bty="n",legend=rp.ndvi_ct.step2020)
R2.ndvi_ct.step2020=summary(fit.ndvi_ct.step2020)$adj.r.squared
R2.ndvi_ct.step2020=round(R2.ndvi_ct.step2020,digits=2)
r.ndvi_ct.step2020=round(predaccuracy.ndvi_ct.step2020,digits=2)
models="stepwise"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.step2020=data.frame(cbind(traits,models,r.ndvi_ct.step2020,R2.ndvi_ct.step2020))
colnames(pred.accuracy.ndvi_ct.step2020)[3:4]=c("r","R2")

### Plotting LASSO regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.las2020, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="LASSO regression with ndvi_ct 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.las2020 <-lm(predindvi_ct.ndvi_ct.las2020 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.las2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.las2020 = vector('expression',2)
rp.ndvi_ct.las2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(summary(fit.ndvi_ct.las2020)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.las2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                   list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.lasso2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.las2020)
R2.ndvi_ct.las2020=summary(fit.ndvi_ct.las2020)$adj.r.squared
R2.ndvi_ct.las2020=round(R2.ndvi_ct.las2020,digits=2)
r.ndvi_ct.las2020=round(predaccuracy.ndvi_ct.lasso2020,digits=2)
models="LASSO"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.las2020=data.frame(cbind(traits,models,r.ndvi_ct.las2020,R2.ndvi_ct.las2020))
colnames(pred.accuracy.ndvi_ct.las2020)[3:4]=c("r","R2")

### Plotting ridge regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.rdg2020, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="Ridge regression with ndvi_ct 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.rdg2020 <-lm(predindvi_ct.ndvi_ct.rdg2020 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.rdg2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.rdg2020 = vector('expression',2)
rp.ndvi_ct.rdg2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                   list(MYVALUE = format(summary(fit.ndvi_ct.rdg2020)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.rdg2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                   list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.ridge2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.rdg2020)
R2.ndvi_ct.rdg2020=summary(fit.ndvi_ct.rdg2020)$adj.r.squared
R2.ndvi_ct.rdg2020=round(R2.ndvi_ct.rdg2020,digits=2)
r.ndvi_ct.rdg2020=round(predaccuracy.ndvi_ct.ridge2020,digits=2)
models="ridge"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.rdg2020=data.frame(cbind(traits,models,r.ndvi_ct.rdg2020,R2.ndvi_ct.rdg2020))
colnames(pred.accuracy.ndvi_ct.rdg2020)[3:4]=c("r","R2")

### Plotting enet regression model with ndvi_ct
plot(ndvi_ct_data$GRYLD, predindvi_ct.ndvi_ct.enet2020, xlab = "Observed Yield",ylab = "Predindvi_cted Yield", main="ElasticNet regression with ndvi_ct 2020",col.main="red",cex.main=1,font.main=2,pch=20,col="red4",cex.lab=1,cex.axis=1, font.lab=2,font.axis=2,col.lab="blue",col.axis="red")
fit.ndvi_ct.enet2020 <-lm(predindvi_ct.ndvi_ct.enet2020 ~ ndvi_ct_data$GRYLD) #fit the model for observed values
abline(fit.ndvi_ct.enet2020, col = 'blue') #add trendline to graph
#code from https://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/
rp.ndvi_ct.enet2020 = vector('expression',2)
rp.ndvi_ct.enet2020[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                                    list(MYVALUE = format(summary(fit.ndvi_ct.enet2020)$adj.r.squared,dig=2)))[2]
rp.ndvi_ct.enet2020[2] = substitute(expression(italic(r) == MYOTHERVALUE), 
                                    list(MYOTHERVALUE = format(predaccuracy.ndvi_ct.enet2020, digits = 2)))[2]
legend("bottomright", bty="n",legend=rp.ndvi_ct.enet2020)
R2.ndvi_ct.enet2020=summary(fit.ndvi_ct.enet2020)$adj.r.squared
R2.ndvi_ct.enet2020=round(R2.ndvi_ct.enet2020,digits=2)
r.ndvi_ct.enet2020=round(predaccuracy.ndvi_ct.enet2020,digits=2)
models="elasticnet"
traits="ndvi_ct"
pred.accuracy.ndvi_ct.enet2020=data.frame(cbind(traits,models,r.ndvi_ct.enet2020,R2.ndvi_ct.enet2020))
colnames(pred.accuracy.ndvi_ct.enet2020)[3:4]=c("r","R2")

multivar.ndvi_ct.prediction.accuracy2020=data.frame(rbind(pred.accuracy.ndvi_ct.step2020,pred.accuracy.ndvi_ct.las2020,pred.accuracy.ndvi_ct.rdg2020,pred.accuracy.ndvi_ct.enet2020))
# write.csv(multivar.ndvi_ct.prediction.accuracy2020,file="Multivariate_predaccuracy_ndvi_ct_2020.csv",row.names = FALSE,quote = FALSE)

multivar.ndvi_ct.allyear=cbind(multivar.ndvi_ct.prediction.accuracy2016,multivar.ndvi_ct.prediction.accuracy2017[,3:4],multivar.ndvi_ct.prediction.accuracy2018[,3:4],multivar.ndvi_ct.prediction.accuracy2019[,3:4],multivar.ndvi_ct.prediction.accuracy2020[,3:4])
write.csv(multivar.ndvi_ct.allyear,file="multivar.yieldpred.ndvi_ct.allYear.csv",row.names=FALSE,quote=FALSE)
