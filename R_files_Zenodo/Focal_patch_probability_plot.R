
#=======FOCAL PATCH SELECTIVITY========#=======HYPOTHESIS 5

setwd("") #Smith_focalpatch_data.csv


#install.packages("tidyverse")
#install.packages("magrittr")
library(tidyverse)
library(magrittr)


set.seed(1234)

 
mydat <- data.frame(read.csv("Smith_focalpatch_data.csv"))

# logistic regression model: 
mod1 <- glm(mydat$focal_patch~mydat$GI_avg, data=mydat, family=binomial(link="logit"))

# new predictor values to use for prediction: 
plotdat <- data.frame(GI=(mydat$GI_avg))

# df with predictions, lower and upper limits of CIs: 
preddat <- predict(mod1,
                   type = "link",
                   newdata=plotdat,
                   se.fit=TRUE) %>% 
  as.data.frame() %>% 
  mutate(GI = (mydat$GI_avg), 
         
         # model object mod1 has a component called linkinv that 
         # is a function that inverts the link function of the GLM:
         lower = mod1$family$linkinv(fit - 1.96*se.fit), 
         point.estimate = mod1$family$linkinv(fit), 
         upper = mod1$family$linkinv(fit + 1.96*se.fit)) 



# plotting with ggplot: 
gonad.plot<-preddat %>% ggplot(aes(x = mydat$GI_avg, 
                       y = point.estimate)) + 
  geom_line(colour = "darkgreen", size = 1) + 
  geom_ribbon(aes(ymin = lower,
                  ymax = upper), fill = "green",
              alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits= c(0,20)) +
  #scale_y_continuous(limits = c(0,1)) OLD CODE
  scale_y_continuous(expand = c(0,0), limits= c(0,1))


print(gonad.plot + theme_bw(base_size = 20) + theme(panel.grid.minor = element_blank()) + geom_hline(yintercept=0.5, linetype="dashed", color = "red") + labs(y="probability (CI)", x= "localized mean gonad index (%)"))


