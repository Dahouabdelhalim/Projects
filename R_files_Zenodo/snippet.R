library(randomForest)
load('Curated_Pull_Request_Data.csv')
# Process Raw Data
merged = data$merged
data$merged = NULL
data$contain_issue_fix = factor(data$contain_issue_fix)
data$contain_test_code = factor(data$contain_test_code)
data$user_accepted_repo = factor(data$user_accepted_repo)
data$dependency = factor(data$dependency)

# Log Transforming Skewed Continuous variables
v = !sapply(data, is.factor)
v[c(2,4)] = FALSE
data[,v] = data.frame(sapply(data[,v], function(x) log(as.numeric(x) +1) ))

# Break data into Training and Test sets: an example
indtrain <- sample(nrow(data), 0.7*nrow(data), replace = FALSE)
Train <- data[indtrain,]
Test <- data[-indtrain,]

yTrain = factor(merged[indtrain])
yTest = factor(merged[-indtrain])
yAll = factor(merged)

# Train Random Forest Model
fit = randomForest(Train, yTrain, importance = T)

