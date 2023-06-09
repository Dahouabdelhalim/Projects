# SET DIRECTORIES (START) ####################################################################
setwd("")
localDir = "."
dataDir = file.path(localDir, "data")
ModelDir = file.path(localDir, "models")
library(Hmsc)
# SET DIRECTORIES AND LOAD LIBRARIES (END) ####################################################################

# LOAD AND PRE-PROCESS DATA (START) ####################################################################
load(file.path(dataDir, "HMSC data"))
# The R-file HMSC data contains the following data matrices:
# Y is the samples times OTUs matrix of the sequence data 
# XData is the samples time covariates dataframe of environmental covariates
# TrData is the OTUs times trais dataframe of the OTU classifications

ny = dim(Y)[1]

prev = colSums(1*(Y>0))
selSP = prev>=ny*0.05
ns = sum(selSP)
Y=Y[,selSP]
TrData = TrData[selSP,]

Yabu = Y
Yabu[Y==0] = NA
Yabu = log(Yabu)
for (i in 1:ns){
  Yabu[,i] = Yabu[,i] - mean(Yabu[,i],na.rm=TRUE)
  Yabu[,i] = Yabu[,i] / sqrt(var(Yabu[,i],na.rm=TRUE))
}
# LOAD AND PRE-PROCESS DATA (END) ####################################################################


# SET UP THE MODEL (START) ###################################################################
# STUDY DESIGN
studyDesign = data.frame(ID = as.factor(as.character(XData$ID..)), transect = XData$transect)
rL.transect = HmscRandomLevel(units = levels(studyDesign$transect))
rL.ID = HmscRandomLevel(units = studyDesign$ID)

# REGRESSION MODEL FOR ENVIRONMENTAL COVARIATES.


# REGRESSION MODEL FOR TRAITS
TrFormula = ~class

# CONSTRUCT AND FIT MODELS (START)  ###################################################################
# MODELTYPE = 1: PRESENCE-ABSENCE DATA; MODELTYPE = 2: ABUNDANCE DATA
# MODEL = 1: NULL; uS; cS
samples = 1000
nChains = 4
for (thin in c(1,10,100,1000)){
  for (modeltype in 1:2){
    for (model in 1:2){
      if(model==1){
        XFormula = ~ sequencing.depth
      }
      if(model==2){
        XFormula = ~ Species + altitude + altitude2 + pH + humidity + active.layer + vegetation + sequencing.depth
      }
      for (dataset in 1:3){
        print(paste0("thin = ",as.character(thin),", modeltype = ",c("pa","abu")[modeltype],", model = ",as.character(model), ", dataset = ",as.character(dataset)))
        if(dataset==1){
          selSP = TrData$class == "Mycorrhizal"
        }
        if(dataset==2){
          selSP = TrData$class == "Endophyte"
        }
        if(dataset==3){
          selSP = !((TrData$class == "Mycorrhizal") | (TrData$class == "Endophyte"))
        }
        if(dataset==4){
          selSP = 1:ns
        }
        set.seed(1)
        m = Hmsc(Y=if(modeltype==1){1*(Y[,selSP]>0)} else {Yabu[,selSP]},
                 XData = XData,  XFormula = XFormula,
                 distr=if(modeltype==1){"probit"} else {"normal"},
                 studyDesign=studyDesign, ranLevels=if(model==1){list(ID=rL.ID)} else {list(ID=rL.ID,transect=rL.transect)})

        ptm = proc.time()
        m = sampleMcmc(m, samples = samples, verbose=samples, thin=thin, adaptNf=rep(ceiling(0.4*samples*thin),c(1,2)[model]), transient = ceiling(0.5*samples*thin),
                       nChains = nChains, nParallel = nChains)
        computational.time =  proc.time() - ptm
        
        print(computeWAIC(m))

        filename = file.path(ModelDir, paste("dataset_",as.character(dataset),"_",c("pa","abundance")[modeltype],
                                             "_model",as.character(model),
                                             "_thin_", as.character(thin),"_samples_", as.character(samples),
                                             "_chains_",as.character(nChains),
                                             ".Rdata",sep = ""))
        save(m,file=filename,computational.time)
      }
    }
  }
}
# CONSTRUCT AND FIT MODELS (END)  ###################################################################
