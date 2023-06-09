### Read in MA data
data <- read.csv("MAData.csv")

### Load required library
library(lme4)

### Analysis corresponding to Table 2
## session as a proxy for time from release
m0 <- glmer(oviposit~mating+hoppers+day+ants+X.N+(1|indiv),
            data=data, family=binomial, subset=host==1,
            control=glmerControl(optimizer="bobyqa"))

### Results presented in Table 2
summary(m0)
