setwd("~/Dropbox/Termites/HumidityPaper/Dryad")
rm(list=ls())
library('Matrix')
library('lme4')
library('sjPlot')

All.data<-read.table('OnMoundData_Shield_Unshielded_Experiments.csv',header = TRUE,sep = ',',stringsAsFactors=TRUE)
Wet.data<-All.data[All.data$Treatment=='Wet',];
Dry.data<-All.data[All.data$Treatment=='Dry',];

model.null   <-lmer(Weight ~            (1| Colony),data= All.data,REML=FALSE);
model.Shield <-lmer(Weight ~ Shielded + (1| Colony),data= All.data,REML=FALSE);
summary(model.Shield)
anova(model.null,model.Shield)
tab_model(model.Shield,show.p = FALSE)

##########
rm(list=ls())
library('Matrix')
library('lme4')
library('sjPlot')

All.data<-read.table('HumidityControlledEnviroment.csv',header = TRUE,sep = ',',stringsAsFactors=TRUE)
M.data<-All.data[All.data$Species=='michaelseni',];

model.null     <-lmer(weight ~            (1| Colony),data= M.data,REML=FALSE);
model.Humidity <-lmer(weight ~ Humidity + (1| Colony),data= M.data,REML=FALSE);
anova(model.null,model.Humidity)
summary(model.Humidity)
tab_model(model.Humidity,show.p = FALSE)

######
rm(list=ls())
All.data<-read.table('HumidityControlledEnviroment.csv',header = TRUE,sep = ',',stringsAsFactors=TRUE)
M.data<-All.data[All.data$Species=='michaelseni',];

model.null     <-lmer(Area ~            (1| Colony),data= M.data,REML=FALSE);
model.Humidity <-lmer(Area ~ Humidity + (1| Colony),data= M.data,REML=FALSE);
anova(model.null,model.Humidity)
summary(model.Humidity)
tab_model(model.Humidity,show.p = FALSE)

######

rm(list=ls())
library('Matrix')
library('ggplot2')
library('lme4')

All.data<-read.table('HumidityControlledEnviroment.csv',header = TRUE,sep = ',',stringsAsFactors=TRUE)

M.data<-All.data[All.data$Species=='michaelseni',];

model.null               <-lmer(weight ~                       (1| Colony),data= M.data,REML=FALSE);
model.Humidity          <-lmer(weight ~  Humidity +          (1| Colony),data= M.data,REML=FALSE);
anova(model.null,model.Humidity)
tab_model(model.Humidity,show.p = FALSE)

# Create a blank dataset with the Humiditys we want
Pred.Hum <- expand.grid(Humidity = seq(33, 92, 1), weight = 0)  

# Create matrix of relevant effect sizes
MM <- model.matrix(terms(model.Humidity), Pred.Hum)  

# Calculate weight based on the relevant effect sizes
Pred.Hum$weight <- MM %*% fixef(model.Humidity)  
pvar.Pred.Hum <- diag(MM %*% tcrossprod(vcov(model.Humidity), MM))

# Add standard errors to the dataframe
Pred.Hum <- data.frame(Pred.Hum, 
                       plo = Pred.Hum$weight - 1.96*sqrt(pvar.Pred.Hum),
                       phi = Pred.Hum$weight + 1.96*sqrt(pvar.Pred.Hum)) 



(mm_basic <- ggplot(M.data, aes(Humidity, weight)) +
    geom_ribbon(data = Pred.Hum, mapping = aes(x = Humidity, ymin = plo, ymax = phi),
                fill = "rosybrown1", alpha = 0.4) +
    geom_line(data = Pred.Hum, mapping = aes(x = Humidity),size=2) +
    geom_point(data = M.data, aes(colour = Colony,size=Species),size=5) +
    theme_bw()+ 
    labs(title = "Weight as a function of humidity", 
         x = "Humidity (%RH)", y = "Weight (g)") +
    scale_y_continuous( limits =c(0,11))+
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 12), 
          axis.title = element_text(size = 12), 
          plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.text = element_text(size = 11),
          legend.title = element_text(face = "bold"),
          legend.position = "bottom", 
          legend.box.background = element_rect(color = "grey", size = 0.3)))




