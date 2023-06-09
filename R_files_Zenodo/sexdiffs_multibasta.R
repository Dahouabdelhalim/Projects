require(BaSTA)

############# ANALYSES #############

#####
furRegNC <- read.delim("furDrySex.txt")
attach(furRegNC)
inputMat<-(furRegNC)
newData <- DataCheck(inputMat, studyStart = 1, studyEnd = 66, autofix = rep(1, 7), silent = FALSE)

multiNF6 <- multibasta(object = furRegNC, studyStart = 1, studyEnd = 66, minAge = 6, models = c("EX","GO", "LO", "WE"), shapes = c("simple", "Makeham","bathtub"), niter = 150000, burnin = 15001, thinning = 50, nsim = 4, parallel = TRUE, ncpus = 4, updateJumps = TRUE)
summary(multiNF6)
coef(multiNF6,showAll=TRUE)
detach(furRegNC)

furRegNC <- read.delim("furWetSex.txt")
attach(furRegNC)
inputMat<-(furRegNC)
newData <- DataCheck(inputMat, studyStart = 1, studyEnd = 66, autofix = rep(1, 7), silent = FALSE)

multiNF6 <- multibasta(object = furRegNC, studyStart = 1, studyEnd = 66, minAge = 6, models = c("EX","GO", "LO", "WE"), shapes = c("simple", "Makeham","bathtub"), niter = 150000, burnin = 15001, thinning = 50, nsim = 4, parallel = TRUE, ncpus = 4, updateJumps = TRUE)
summary(multiNF6)
coef(multiNF6,showAll=TRUE)

#######
kadRegNC <- read.delim("kadDrySex.txt")
attach(kadRegNC)
inputMat<-(kadRegNC)
newData <- DataCheck(inputMat, studyStart = 1, studyEnd = 66, autofix = rep(1, 7), silent = FALSE)

multiNK6 <- multibasta(object = kadRegNC, studyStart = 1, studyEnd = 66, minAge = 6, models = c("EX","GO", "LO", "WE"), shapes = c("simple", "Makeham","bathtub"),  niter = 150000, burnin = 15001, thinning = 50, nsim = 4, parallel = TRUE, ncpus = 4, updateJumps = TRUE)
summary(multiNK6)
coef(multiNK6,showAll=TRUE)
detach(kadRegNC)


kadRegNC <- read.delim("kadWetSex.txt")
attach(kadRegNC)
inputMat<-(kadRegNC)
newData <- DataCheck(inputMat, studyStart = 1, studyEnd = 66, autofix = rep(1, 7), silent = FALSE)

multiNK6 <- multibasta(object = kadRegNC, studyStart = 1, studyEnd = 66, minAge = 6, models = c("EX","GO", "LO", "WE"), shapes = c("simple", "Makeham","bathtub"),  niter = 150000, burnin = 15001, thinning = 50, nsim = 4, parallel = TRUE, ncpus = 4, updateJumps = TRUE)
summary(multiNK6)
coef(multiNK6,showAll=TRUE)

######
ortRegNC <- read.delim("ortDrySex.txt")
attach(ortRegNC)
inputMat<-(ortRegNC)
newData <- DataCheck(inputMat, studyStart = 1, studyEnd = 88, autofix = rep(1, 7), silent = FALSE)

multiNO <- multibasta(object = ortRegNC, studyStart = 1, studyEnd = 88, minAge = 6, models = c("EX","GO", "LO","WE"), shapes = c("simple", "Makeham","bathtub"), niter = 150000, burnin = 15001, thinning = 50, nsim = 4, parallel = TRUE, ncpus = 4, updateJumps = TRUE)
summary(multiNO)
coef(multiNO,showAll=TRUE)

detach(ortDrySex)

ortRegNC <- read.delim("ortWetSex.txt")
attach(ortRegNC)
inputMat<-(ortRegNC)
newData <- DataCheck(inputMat, studyStart = 1, studyEnd = 88, autofix = rep(1, 7), silent = FALSE)

multiNO <- multibasta(object = ortRegNC, studyStart = 1, studyEnd = 88, minAge = 6, models = c("EX","GO", "LO","WE"), shapes = c("simple", "Makeham","bathtub"), niter = 150000, burnin = 15001, thinning = 50, nsim = 4, parallel = TRUE, ncpus = 4, updateJumps = TRUE)
summary(multiNO)
coef(multiNO,showAll=TRUE)

########
pieRegNC <- read.delim("pieDrySex.txt")
attach(pieRegNC)
inputMat<-(pieRegNC)
newData <- DataCheck(inputMat, studyStart = 1, studyEnd = 95, autofix = rep(1, 7), silent = FALSE)

multiNP <- multibasta(object = pieRegNC, studyStart = 1, studyEnd = 95, minAge = 8, models = c("EX","GO", "LO", "WE"), shapes = c("simple", "Makeham","bathtub"), niter = 150000, burnin = 15001, thinning = 50, nsim = 4, parallel = TRUE, ncpus = 4, updateJumps = TRUE)
summary(multiNP)
coef(multiNP,showAll=TRUE)

detach(pieRegNC)

pieRegNC <- read.delim("pieWetSex.txt")
attach(pieRegNC)
inputMat<-(pieRegNC)
newData <- DataCheck(inputMat, studyStart = 1, studyEnd = 95, autofix = rep(1, 7), silent = FALSE)

multiNP <- multibasta(object = pieRegNC, studyStart = 1, studyEnd = 95, minAge = 8, models = c("EX","GO", "LO", "WE"), shapes = c("simple", "Makeham","bathtub"), niter = 150000, burnin = 15001, thinning = 50, nsim = 4, parallel = TRUE, ncpus = 4, updateJumps = TRUE)
summary(multiNP)
coef(multiNP,showAll=TRUE)
