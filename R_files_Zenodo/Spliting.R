#50%
setwd("C:/Users/l756n005/Desktop/MODELLING/Splitting/P_ayeaye/apos_TOWN")

all <- read.csv("Payeaye_occ_inside.csv")

all <- unique(all)

all$check <- paste(all[,2], all[,3], sep = "_")
train <- all[sample(nrow(all), round((length(all[,1])/2))), ]
test <- all[!all[,4] %in% train[,4],]

all$check <- NULL
train$check <- NULL
test$check <- NULL

write.csv(all,"Payeaye_joint.csv",row.names = FALSE)
write.csv(train,"Payeaye_train.csv",row.names = FALSE)
write.csv(test,"Payeaye_test.csv",row.names = FALSE)



#75%
setwd("C:/Users/l756n005/Desktop/MODELLING/Splitting/P_ayeaye/apos_TOWN")

all <- read.csv("Payeaye_occ_inside.csv")

all <- unique(all)

all$check <- paste(all[,2], all[,3], sep = "_")
train <- all[sample(nrow(all), round((length(all[,1])/4 *3))), ]
test <- all[!all[,4] %in% train[,4],]

all$check <- NULL
train$check <- NULL
test$check <- NULL

write.csv(all,"Payeaye_joint.csv",row.names = FALSE)
write.csv(train,"Payeaye_train.csv",row.names = FALSE)
write.csv(test,"Payeaye_test.csv",row.names = FALSE)
