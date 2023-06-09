# SET DIRECTORIES AND LOAD LIBRARIES (START) ##################################################################
setwd("")
localDir = "."
dataDir = file.path(localDir, "data")
ModelDir = file.path(localDir, "models")
library(Hmsc)
library(corrplot)
# SET DIRECTORIES AND LOAD LIBRARIES (END) ####################################################################


samples = 1000
nChains = 4

load(file.path(ModelDir,"TrData.RData"))
sp1 = TrData$class == "Mycorrhizal"
sp2 = TrData$class == "Endophyte"
sp3 = !((TrData$class == "Mycorrhizal") | (TrData$class == "Endophyte"))
n1=sum(sp1)
n2=sum(sp2)
n3=sum(sp3)

thin_max=1000
for (modeltype in 2:2){
   for (dataset in 4:4){
      for (model in c(2,1)){
         thin=thin_max
         cont=TRUE
         while(cont){
            filename = file.path(ModelDir, paste("dataset_",as.character(dataset),"_",c("pa","abundance")[modeltype],
                                                 "_model",as.character(model),
                                                 "_thin_", as.character(thin),"_samples_", as.character(samples),
                                                 "_chains_",as.character(nChains),
                                                 ".Rdata",sep = ""))
            if(file.exists(filename)){
               load(filename)
               cont=FALSE
            } else {
               thin=thin/10}
            if (thin<1){
               print("not found")
               cont=FALSE
            }
         }
         

         OmegaCor = computeAssociations(m)
         supportLevel = 0.95
         toPlot = ((OmegaCor[[1]]$support > supportLevel) + (OmegaCor[[1]]$support < (1-supportLevel)) > 0) * OmegaCor[[1]]$mean
 
         pos = (sum((OmegaCor[[1]]$support > supportLevel))-m$ns)/(m$ns*(m$ns-1))
         neg = mean((OmegaCor[[1]]$support < (1-supportLevel)))

         cat("modeltype = ",modeltype,"dataset = ",dataset, "model = ", model,"pos = ",pos, "neg = ",neg,"\\n")
         if(dataset==4){
            
            if(TRUE){
               predY = computePredictedValues(m, expected=TRUE)
               MF = evaluateModelFit(hM=m, predY=predY)
               if(modeltype==1){
                  ta = cbind(MF$AUC,MF$TjurR2)
                  colnames(ta) = c("AUC","TjurR2")
               } else {
                  ta = cbind(MF$R2)
                  colnames(ta) = c("R2")
               }
               filnam = paste0("panels/MF_",c("pa","abu")[modeltype],"_coc_",c("raw_","res_")[model])
               row.names(ta) = m$spNames
               write.csv(ta,file = paste0(filnam,".csv"))
            }
            
            plotOrder1 = which(sp1)[corrMatOrder(toPlot[sp1,sp1], order = "AOE")]
            plotOrder2 = which(sp2)[corrMatOrder(toPlot[sp2,sp2], order = "AOE")]
            plotOrder3 = which(sp3)[corrMatOrder(toPlot[sp3,sp3], order = "AOE")]

            plotOrder = c(plotOrder1,plotOrder2,plotOrder3)
            
            filnam = paste0("panels/",c("pa","abu")[modeltype],"_coc_",c("raw_","res_")[model],"order_groups")
            tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
            corrplot(toPlot[plotOrder,plotOrder], method = "color",  col = c("blue","white","red"),tl.cex=0.2)
            dev.off()
            namelist=as.data.frame(m$spNames[plotOrder])
            write.csv(namelist,file = paste0(filnam,".csv"))
            
            if(model==2){
              specialOrder = plotOrder
            }
            
            filnam = paste0("panels/",c("pa","abu")[modeltype],"_coc_",c("raw_","res_")[model],"order_groups_residual")
            tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
            corrplot(toPlot[specialOrder,specialOrder], method = "color",  col = c("blue","white","red"),tl.cex=0.2)
            dev.off()
            namelist=as.data.frame(m$spNames[specialOrder])
            write.csv(namelist,file = paste0(filnam,".csv"))
            
            plotOrder = corrMatOrder(toPlot)
            filnam = paste0("panels/",c("pa","abu")[modeltype],"_coc_",c("raw_","res_")[model],"order_global")
            tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
            corrplot(toPlot[plotOrder,plotOrder], method = "color",  col = c("blue","white","red"),tl.cex=0.2)
            dev.off()
            namelist=as.data.frame(m$spNames[plotOrder])
            write.csv(namelist,file = paste0(filnam,".csv"))

            if (model==2){
               groupnames = c("plant","altitude","fixed: transect","sequencing depth")
               group = c(1,1,1,1,1,2,2,3,3,3,3,4)
               VP = computeVariancePartitioning(m,group = group,groupnames = groupnames)
               filnam = paste0("panels/",c("pa","abu")[modeltype],"_VP_",c("raw_","res_")[model])
               tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
               plotVariancePartitioning(m,VP)
               dev.off()

               VP = computeVariancePartitioning(m)
               filnam = paste0("panels/",c("pa","abu")[modeltype],"_VP_not_grouped_",c("raw_","res_")[model])
               tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
               plotVariancePartitioning(m,VP)
               dev.off()
            }


            p11=(sum(OmegaCor[[1]]$support[sp1,sp1] > supportLevel)-n1)/(n1*(n1-1))
            p22=(sum(OmegaCor[[1]]$support[sp2,sp2] > supportLevel)-n2)/(n2*(n2-1))
            p33=(sum(OmegaCor[[1]]$support[sp3,sp3] > supportLevel)-n3)/(n3*(n3-1))
            p12=(sum(OmegaCor[[1]]$support[sp1,sp2] > supportLevel))/(n1*n2)
            p13=(sum(OmegaCor[[1]]$support[sp1,sp3] > supportLevel))/(n1*n3)
            p23=(sum(OmegaCor[[1]]$support[sp2,sp3] > supportLevel))/(n2*n3)
            res = c(p11,p22,p33,p12,p13,p23)
            cat("p11,p22,p33,p12,p13,p23 = ",res,"\\n")
            n11 = mean((OmegaCor[[1]]$support[sp1,sp1] < (1-supportLevel)))
            n22 = mean((OmegaCor[[1]]$support[sp2,sp2] < (1-supportLevel)))
            n33 = mean((OmegaCor[[1]]$support[sp3,sp3] < (1-supportLevel)))
            n12 = mean((OmegaCor[[1]]$support[sp1,sp2] < (1-supportLevel)))
            n13 = mean((OmegaCor[[1]]$support[sp1,sp3] < (1-supportLevel)))
            n23 = mean((OmegaCor[[1]]$support[sp2,sp3] < (1-supportLevel)))
            res = c(n11,n22,n33,n12,n13,n23)
            cat("n11,n22,n33,n12,n13,n23 = ",res,"\\n")
         }
      }
   }
}


