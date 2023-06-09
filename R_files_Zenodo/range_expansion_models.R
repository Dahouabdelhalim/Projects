## Range Expansion models
# this script determines model fit of three competing models of range expansion described in Table 1 and Table S4

#set working directory where the the output of thetagenerator_2020.R and the haversine distances between populations and range expansion origins as depicted in Figure 5
df <- read.csv("~/allmodels_diversity.csv")

## Model 1: Nested Origins
# thetaW
m1_W <- lm(df$thetaw ~ df$coldist)
summary(m1_W) # Adjusted R-squared: 0.6112 P-value: 3.906e-06
AIC(m1_W) # -292.6228
BIC(m1_W) # -289.0886
# thetapi
m1_pi <- lm(df$thetapi ~ df$coldist)
summary(m1_pi) # Adjusted R-squared: 0.4944 P-value: 7.639e-05
AIC(m1_pi) # -288.9176
BIC(m1_pi) # -285.3835 
# tajimaD
m1_D <- lm(df$tajimaD ~ df$coldist)
summary(m1_D) # Adjusted R-squared: 0.5451 P-Value: 2.297e-05
AIC(m1_D) # -58.5981
BIC(m1_D) # -55.06394

## Model 2: Independent Origins
# thetaW
m2_W <- lm(df$thetaw ~ df$origin2)
summary(m2_W) # Adjusted R-squared: 0.3108 P-value: 0.002748
AIC(m2_W) # -278.8854
BIC(m2_W) # -275.3512
# thetapi
m2_pi <- lm(df$thetapi ~ df$origin2)
summary(m2_pi) # Adjusted R-squared: 0.275  P-value: 0.005007
AIC(m2_pi) # -280.2661
BIC(m2_pi) #  -276.7319
# tajimaD
m2_D <- lm(df$tajimaD ~ df$origin2)
summary(m2_D) # Adjusted R-squared: 0.4504 P-Value: 0.0001985
AIC(m2_D) # -54.05752
BIC(m2_D) # -50.52335

## Model 3: Single Origin
# thetaW
m3_W <- lm(df$thetaw ~ df$origin1)
summary(m3_W) # Adjusted R-squared: 0.397 P-value: 0.0005772
AIC(m3_W) # -282.0941
BIC(m3_W) # -278.56
# thetapi
m3_pi <- lm(df$thetapi ~ df$origin1)
summary(m3_pi) # Adjusted R-squared: 0.3418 P-value: 0.0016
AIC(m3_pi) # -282.5873
BIC(m3_pi) # -279.0532
# tajimaD
m3_D <- lm(df$tajimaD ~ df$origin1)
summary(m3_D) # Adjusted R-squared: 0.2966 P-Value: 0.003498
AIC(m3_D) # -48.13504
BIC(m3_D) # -44.60088

library(AICcmodavg)
W_list <- list(m1_W, m2_W, m3_W)
pi_list <- list(m1_pi, m2_pi, m3_pi)
D_list <- list (m1_D, m2_D, m3_D)

aictab(W_list, second.ord = FALSE)
aictab(pi_list, second.ord = FALSE)
aictab(D_list, second.ord = FALSE)

bictab(W_list)
bictab(pi_list)
bic(D_list)