###Pollen-based biome reconstruction in R
#read data
setwd('...') #set path

pollen <- read.csv("pollen_data.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

variable <- read.csv("pollen_variable.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
colnames(variable)[1] <- "taxa"

bio <- read.csv("pft-biome.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
colnames(bio)[1] <- "biome"

pft <- read.csv("pollen-pft.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
colnames(pft)[1] <- "taxa"

####calibrate pollen percentage using the threhold and weight
pollen.no <- nrow(variable)
pollennames <- variable$taxa

#### threhold and weight
pollen1 <- pollen[ ,names(pollen) %in% pollennames]

for(k in 1:pollen.no){
  a <- as.vector(unlist(pollen[pollennames[k]])) # subtract pollen percentage 
  
  a[(a < variable[k,'threhold'])] <- 0  # set the pollen percentage below the threshold to 0 
  
  pollen1[ ,pollennames[k]] <- a * variable[k,"weight"] # multiply the percentage value by weight
}

pollen[ ,names(pollen) %in% pollennames] <- pollen1

####calculate the biome score sample by sample
###set the result matrix
sample.no <- nrow(pollen)
info.no <- ncol(pollen)-pollen.no
biome.no <- nrow(bio)
biomenames <- bio$biome
res <- matrix(0, nrow=sample.no, ncol=info.no+biome.no,byrow=TRUE)
res <- data.frame(res)
res[,1:info.no] <- pollen[1:info.no]
colnames(res) <- c(colnames(pollen)[1:info.no],biomenames)
####calculation
for(i in 1:sample.no){for(j in 1:biome.no){biome1 <- bio$biome[j]
                                           pft1 <- colnames(bio)[bio[j,]==1]
                                           pollen2 <- pollennames[rowSums(pft[,pft1])>0]
                                           res[i,biome1] <- sum(sqrt(pollen[i,pollen2])) 
                                           }
                      print(i)
                      }


for(i in 1:sample.no){res[i,"Best"]<- biomenames[which.max(res[i,biomenames])]
                      res[i,"Best.num"]<- which.max(res[i,biomenames])}
write.csv(res, file="Biome_result.csv")

##END


