# SET DIRECTORIES AND LOAD LIBRARIES (START) ##################################################################
setwd("")
localDir = "."
dataDir = file.path(localDir, "data")
ModelDir = file.path(localDir, "models")
MixingDir = file.path(localDir, "mixing")
library(Hmsc)
# SET DIRECTORIES AND LOAD LIBRARIES (END) ####################################################################

# PLOT MIXING STATISTICS (START) ####################################################################
for (zzz in 1:2){
  tiff(filename=c(paste0("panels/effective_sample_size.tif"),paste0("panels/psrf.tif"))[zzz], res=300, unit="cm", height=17, width=17)
  
  par(mfrow=c(2,2))
  for (modeltype in 1:2){
    for (model in 1:2){
      dataset = 4
      
      filename = file.path(MixingDir, paste("dataset_",as.character(dataset),"_",c("pa","abundance")[modeltype],
                                            "_model",as.character(model),".Rdata",sep = ""))
      load(filename)
      resa=list()
      for (param in 1:2){
        if(param == 1){
          if(zzz==1){
            resa[[param]]=as.vector(mixing$es.beta)
          } else {
            resa[[param]]=as.vector(mixing$ge.beta)
          }
        } else {
          if(zzz==1){
            resa[[param]]=as.vector(mixing$es.omega)
          } else {
            resa[[param]]=as.vector(mixing$ge.omega)
          } 
        }
        
      }
      boxplot(resa, main = paste0(c("P-A","Abu")[[modeltype]]," (model ",as.character(model),")"),names=c("beta","omega"),ylim=if(zzz==1){c(0,6000)} else {c(0,2)})
    }
  } 
  dev.off()
}