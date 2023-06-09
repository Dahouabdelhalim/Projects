# SET DIRECTORIES AND LOAD LIBRARIES (START) ##################################################################
setwd("")
localDir = "."
dataDir = file.path(localDir, "data")
ModelDir = file.path(localDir, "models")
MixingDir = file.path(localDir, "mixing")
library(Hmsc)
# SET DIRECTORIES AND LOAD LIBRARIES (END) ####################################################################

# READ MODELS (START) ####################################################################

samples = 1000
nChains = 4

thin_max=1000

for (modeltype in 1:2){
   for (dataset in 1:4){
      for (model in 1:2){
        cat(modeltype,dataset,model)
         cont=TRUE
         thin=thin_max
         while(cont){
            filename = file.path(ModelDir, paste("dataset_",as.character(dataset),"_",c("pa","abundance")[modeltype],
                                                 "_model",as.character(model),
                                                 "_thin_", as.character(thin),"_samples_", as.character(samples),
                                                 "_chains_",as.character(nChains),
                                                 ".Rdata",sep = ""))
            if(file.exists(filename)){
               load(filename)

               mpost = convertToCodaObject(m)

               es.beta = effectiveSize(mpost$Beta)
               ge.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf

               es.gamma = effectiveSize(mpost$Gamma)
               ge.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf

               es.V = effectiveSize(mpost$V)
               ge.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf

               mpost$temp = mpost$Omega[[1]]
               for(i in 1:length(mpost$temp)){
                  mpost$temp[[i]] = mpost$temp[[i]][,1:1000]
               }
               es.omega = effectiveSize(mpost$temp)
               ge.omega = gelman.diag(mpost$temp,multivariate=FALSE)$psrf


               mixing = list(es.beta=es.beta, ge.beta=ge.beta,
                             es.gamma=es.gamma, ge.gamma=ge.gamma,
                             es.V=es.V, ge.V=ge.V,
                             es.omega=es.omega, ge.omega=ge.omega)

               filename = file.path(MixingDir, paste("dataset_",as.character(dataset),"_",c("pa","abundance")[modeltype],
                                                     "_model",as.character(model),".Rdata",sep = ""))

               save(file=filename, mixing)

               cont=FALSE
            } else {
               thin=thin/10}
            if (thin<1){
               print("not found")
               cont=FALSE
            }
         }
      }
   }
}


