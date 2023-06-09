### RMark Extinction and Origination Models ###
### NOTE: This process can take 96 hours or more###

load(file="Logvol_encounter_history_Big3.Rdata")
 
library(RMark)

#create external dataframes to hold coefficients and model fit table for each class
holder <- data.frame(matrix(ncol=8, nrow=1))
names(holder) <- c("class", "Phi.coef", "Phi.se", "Gamma.coef","Gamma.se", "p.coef", "p.se", "regime")
AIC.class <- data.frame(matrix(ncol=7, nrow=4))
names(AIC.class) <- c("class", "model", "npar", "AICc", "DeltaAICc", "weight", "deviance")

#create internal dataframes to hold coefficients each class
BI1 <- data.frame(matrix(ncol=8, nrow=1))
BI2 <- data.frame(matrix(ncol=8, nrow=1))
BI3 <- data.frame(matrix(ncol=8, nrow=1))
BI4 <- data.frame(matrix(ncol=8, nrow=1))
MRI1 <- data.frame(matrix(ncol=8, nrow=1))
MRI2 <- data.frame(matrix(ncol=8, nrow=1))
names(BI1) <- c("class", "Phi.coef", "Phi.se", "Gamma.coef","Gamma.se", "p.coef", "p.se", "regime")
names(BI2) <- c("class", "Phi.coef", "Phi.se", "Gamma.coef","Gamma.se", "p.coef", "p.se", "regime")
names(BI3)  <- c("class", "Phi.coef", "Phi.se", "Gamma.coef","Gamma.se", "p.coef", "p.se", "regime")
names(BI4)  <- c("class", "Phi.coef", "Phi.se", "Gamma.coef","Gamma.se", "p.coef", "p.se", "regime")
names(MRI1)  <- c("class", "Phi.coef", "Phi.se", "Gamma.coef","Gamma.se", "p.coef", "p.se", "regime")
names(MRI2)  <- c("class", "Phi.coef", "Phi.se", "Gamma.coef","Gamma.se", "p.coef", "p.se", "regime")

#loop CMR models by class
classes <- as.character(sort(unique(nd4$class)))

for (i in 1:length(classes)) {
	nd7 <- subset(nd4, class==classes[i])
	nd7.prad <- process.data(nd7, model="Pradsen")
    print(classes[i])
        
    tryCatch({     
	bodysize.bias <-function() {

   	   	#extinction models 
       	Phi.ME <- list(formula=~time + logvol * mass.ext)
       	Phi.BE <- list(formula=~time + logvol)
   
   		#sampling model 
       	p.logvol <- list(formula=~time + logvol)
   
       	#origination models
       	Gamma.ME <- list(formula=~time + logvol * rec.mass)
       	Gamma.BE <- list(formula=~time + logvol)
       	
   
       	#create models to run by combining the above formulas 
       	model.list <- create.model.list("Pradsen")
   
   		#run all the models
   		design.parameters=list(Phi=list(formula=~mass.ext), Gamma=list(formula=~rec.mass))
   		results <- mark.wrapper(model.list, data=nd7.prad, design.parameters=design.parameters, output=FALSE, silent=FALSE)
   		return(results)
		}

		bodysize.models <- bodysize.bias() #results saved to object 'bodysize.models'
		save(bodysize.models, file=paste(classes[i],"CMR_bodysize_model_big3_results","Rdata", sep="."))
	
		#extract CMR coefficients for each class
		BI1[1,1] <- classes[i]
    	BI1[1,2:3] <- bodysize.models$Phi.BE.p.logvol.Gamma.BE$results$beta["Phi:logvol",1:2]
    	BI1[1,4:5] <- bodysize.models$Phi.BE.p.logvol.Gamma.BE$results$beta["Gamma:logvol",1:2]
    	BI1[1,6:7] <- bodysize.models$Phi.BE.p.logvol.Gamma.BE$results$beta["p:logvol",1:2]
    	BI1[1,8] <- "background-combined"
    	
    	BI2[1,1] <- classes[i]
    	BI2[1,2:3] <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta["Phi:logvol",1:2]
    	BI2[1,4:5] <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta["Gamma:logvol",1:2]
    	BI2[1,6:7] <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta["p:logvol",1:2]
    	BI2[1,8] <- "background-MEcombined"
    	
    	BI3[1,1] <- classes[i]
    	BI3[1,2:3] <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta["Phi:logvol",1:2]
    	BI3[1,4:5] <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta["Gamma:logvol",1:2]
    	BI3[1,6:7] <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta["p:logvol",1:2]
    	BI3[1,8] <- "background-BEME"
    	
    	BI4[1,1] <- classes[i]
    	BI4[1,2:3] <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta["Phi:logvol",1:2]
    	BI4[1,4:5] <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta["Gamma:logvol",1:2]
    	BI4[1,6:7] <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta["p:logvol",1:2]
    	BI4[1,8] <- "background-MEBE"
    
    	MRI1[1,1] <- classes[i]
    	MRI1[1,2:3] <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta["Phi:logvol:mass.ext",1:2]
    	MRI1[1,4:5] <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta["Gamma:logvol:rec.mass",1:2]
    	MRI1[1,6:7] <- bodysize.models$Phi.ME.p.logvol.Gamma.ME$results$beta["p:logvol",1:2]
    	MRI1[1,8] <- "mass.extinction/recovery-combined"
    	
    	MRI2[1,1] <- classes[i]
    	MRI2[1,2:3] <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta["Phi:logvol:mass.ext",1:2]
    	MRI2[1,4:5] <- bodysize.models$Phi.BE.p.logvol.Gamma.ME$results$beta["Gamma:logvol:rec.mass",1:2]
    	MRI2[1,6:7] <- bodysize.models$Phi.ME.p.logvol.Gamma.BE$results$beta["p:logvol",1:2]
    	MRI2[1,8] <- "mass.extinction/recovery-mixed"
 

		coefficients.classes <- rbind(BI1, BI2, BI3, BI4, MRI1, MRI2)
		holder <- rbind(holder, coefficients.classes)
		
		#extract model selection table for each class
		holder2 <-  data.frame(matrix(ncol=7, nrow=nrow(bodysize.models$model.table)))
		names(holder2) <- c("class", "model", "npar", "AICc", "DeltaAICc", "weight", "deviance")
    	holder2[,1] <- classes[i]
    	holder2[,2:7] <- bodysize.models$model.table[4:9]
    	
    	AIC.class <- rbind(holder2, AIC.class)
    	}, error=function(e){cat("Model ERROR :",classes[i], "\\n")})
}

coefficients.classes <- na.omit(holder)
coefficients.classes

write.table(coefficients.classes, file="CMR_coeffs_bodysize_big3.txt", sep="\\t")

AIC.class <- na.omit(AIC.class)
AIC.class

write.table(AIC.class, file="AIC_table_bodysize_big3.txt", sep="\\t")

