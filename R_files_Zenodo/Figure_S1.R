## FIGURE S1A ##
# G-test goodness-of-fit
library(DescTools)

# Plant+Lin vs Plant+Hexane
observed=c(2,19)    
expected=c(0.5, 0.5)
GTest(x=observed,p=expected,correct="none")

# Left vs right
observed=c(10,11)    
expected=c(0.5, 0.5)
GTest(x=observed,p=expected,correct="none")


## FIGURE S1B##
library(MASS) 

# Plant+Lin vs Plant+Hexane, file: mirabilis_lin.csv
wilcox.test(mirabilis_lin$plant.Linalool, mirabilis_lin$plant.Hexane, paired=TRUE)

# Left vs right, file: mirabilis_side.csv
wilcox.test(mirabilis_side$left, mirabilis_side$right, paired=TRUE)
