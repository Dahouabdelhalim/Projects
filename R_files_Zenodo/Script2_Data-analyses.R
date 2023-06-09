
#### Script 2 - Data Analysis    

# Import data from existing script (Relative Path)
source("Script1_Data-preparation.R")

# Load packages
library(car)
library(arm)
library(coefplot)
library(ggplot2)

# Mean and median of impact factor
# Mean:
regression_data$if. <- as.numeric(regression_data$if.)
mean(regression_data$if., na.rm = TRUE)
# Median:
median(regression_data$if., na.rm = TRUE)

# Remove journals from countries that 'Animal Protection Index' does not score (n = 8)

regression_data <- regression_data[!(regression_data$country=="Taiwan"),]
regression_data <- regression_data[!(regression_data$country=="Slovakia"),]
regression_data <- regression_data[!(regression_data$country=="Singapore"),]
regression_data <- regression_data[!(regression_data$country=="Hungary"),]
regression_data <- regression_data[!(regression_data$country=="Finland"),]
regression_data <- regression_data[!(regression_data$country=="Czech Republic"),]
regression_data <- regression_data[!(regression_data$country=="Belgium"),]
regression_data <- regression_data[!(regression_data$country=="Bulgaria"),]

## MODEL 1

# Run glm binomial model to determine whether we can predict fulfillment of the criterion "Does the journal have any statement realted to animal care" 

regression_data$final_any_ac_state_YN <- as.numeric(regression_data$final_any_ac_state_YN)
regression_data$welfare_law <- as.numeric(regression_data$welfare_law)
regression_data$if. <- as.numeric(regression_data$if.)
regression_data$conserv_YN <- as.numeric(regression_data$conserv_YN)

# Fit glm (the variables in this model have not been centered or scaled)
model1 <- glm(final_any_ac_state_YN ~ if. + oa_YN + welfare_law + conserv_YN + welfare_law*oa_YN, binomial, data = regression_data)

summary(model1)
plot(model1)

# Variance and normality 
par(mfrow = c(2,2))
plot(model1)

# variance inflation factors
vif(model1)

## MODEL 1b - Centered and Scaled 
# Center and scale impact factor, and center binary variables
regression_data_rescale <- mutate(regression_data, 
                        impactfact.center = rescale(regression_data$if.),
                        welfarelaw.center = rescale(regression_data$welfare_law),
                        openaccess.center = rescale(regression_data$oa_YN),
                        conservation.center = rescale(regression_data$conserv_YN)
)

# Fit glm with independent variables centeres and scaled
model_rescale <- glm(final_any_ac_state_YN ~ impactfact.center + openaccess.center + welfarelaw.center + conservation.center + welfarelaw.center*openaccess.center, binomial, data = regression_data_rescale)

saveRDS(model_rescale, file = "model_rescale.rds")

summary(model_rescale)

# Variance and normality 
x11()
par(mfrow = c(2,2))
plot(model_rescale)

# Variance inflation factors
vif(model_rescale)

# Visually explore coefficient plot (Script 4 contains code to produce coefficient plot in manuscript)

coefplot:::buildModelCI(model_rescale)
cp <- coefplot(model_rescale, decreasing = FALSE, innerCI = 1, outerCI = 2, sort = "magnitude", color = "black", lwdInner = 1.5, lwdOuter = 1, pointSize = 4.5, newNames = c(welfarelaw.center = "Animal Welfare Legislation", impactfact.center = "Impact Factor", conservation.center = "Conservation Status", openaccess.center = "Open Access Status"), title = "")
cp <- cp + theme_classic() 
cp <- cp + geom_point(size = 4, pch = 21, fill = "white")
print(cp)
ggsave("coef_plot.png", width = 5, height = 3)

## HISTOGRAM- FREQUENCY OF JOURNAL CRITERIA FULFILLMENT

# Remove the nested or co-dependent criteria because they are not independent (violate GLM assumptions). Remove 'final_ac_num' and 'final_ac_inst' because their fulfillment depends on whether "final_ac_approve" was fulfilled.

regression_data$final_ac_num <- NULL
regression_data$final_ac_inst <- NULL

# Remove rows (journals) that score 0 for the criteria, "Does the journal have any statement related to animal care" because we are only including journals that DO have a statement related to animal care.

regression_data <- regression_data[!(regression_data$final_any_ac_state_YN=="0"),]

# Create a "response score" column. Response score is the number of criteria each journal fulfills.
regression_data<-mutate(regression_data, 
                        final_ac_approve_present=ifelse(as.numeric(final_ac_approve)>0, 1, 0),
                        final_field_YN_present=ifelse(as.numeric(final_field_YN)>0, 1, 0),
                        final_three_r_present=ifelse(as.numeric(final_three_r)>0, 1, 0),
                        condition_present=ifelse(as.numeric(condition)>0, 1, 0))

# Add up number of criteria fulfilled (response score)
response_score <- rowSums(regression_data [ ,13:16])
response_score <- as.numeric(response_score)

freq <- table(regression_data$response_score)

regression_data <- mutate(regression_data, response_score)

x11()
hist <- ggplot(regression_data, aes(x = regression_data$response_score))+
  geom_histogram(binwidth = 0.2,   colour = "black", fill = "lightgrey") + theme_classic() + geom_bar(position = "dodge") + xlab("Number of Criteria Fulfilled") + ylab("Frequency")
print(hist)

### MODEL 2

# Run poisson model to determine whether we can predict the number of criteria journals fulfill.
model2 <- glm(response_score ~ as.numeric(if.)+ oa_YN+ welfare_law+ conserv_YN+ welfare_law*oa_YN, poisson, data = regression_data)
summary(model2)

#variance inflation factors 
vif(model2)

#Variance and normality 
par(mfrow = c(2,2))
plot(model2)

