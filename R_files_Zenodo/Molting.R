# Read file or use file.choose()
data <- read.csv("Molting updated for R.csv")
str(data)

# Run generalized linear regression (binomial)
m <- glm(infection_status ~ time_interval_between_exposure_and_molting_h * age_at_exposure_d, family=binomial(link=logit), data=data)
summary(m)

# Calculate various pseudo R^2
library(DescTools)
PseudoR2(m, which="all")
