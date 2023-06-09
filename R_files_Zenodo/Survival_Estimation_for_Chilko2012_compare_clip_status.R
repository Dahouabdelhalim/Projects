#//******************************************************
#//	Programmer: Kintama Research Services (kintama.com)
#//	Project Name: Ken Jeffries et al. Immune response genes and pathogen presence predict migration survival in wild salmon smolts. Accepted by Molecular Ecology. 
#//	Program Date:  Aug 2013
#//	Version:1
#//	Comments:  RMark Survival analysis of 2012 V7 tagged Chilko Lake sockeye salmon submitted to Dryad (http://datadryad.org/)
#//	
#//	
#//******************************************************/


#Note the capture history (CH) sequence here is as follows:
  #CH= "RELEASE","HENRYS_BRIDGE","SIWASH_BRIDGE","DERBY","SOUTHARM","NSOG","QCS"

###############################
#Packages
###############################
library(RMark)


#################################
#Import datafile and process for analysis
#################################
#Note: removed 8 fish that were clipped but that did not have molecular data in
  #order to match the survival results with the publication.

#Import capture history file
Chilko2012V7<-import.chdata("CH_RMark2013-08-13_V7_NO_FAR_or_MISS.txt", header=TRUE,
                            field.types=c("f"))

summary(Chilko2012V7) 

#Process the data and make the design data list (ddl)
Chilko2012V7.process=process.data(Chilko2012V7,model="CJS", groups=c("treatment_type"))

Chilko2012V7.ddl=make.design.data(Chilko2012V7.process, remove.unused=TRUE)



###############################################
#Estimate survival and test the effect of gill clipping
###############################################
#Note:  c-hat estimate from phi(time*treatment_type) p(time*treatment_type) model,
  #estimated in Program Mark using 10 points and 100 replicates (1200 simulations).

#Note: Fix the detection efficiency (p) of the QCS detection site to 0.67
  #which is the measured value for V7-tagged fish 2005-2007 as per 
  #Welch, DW, Melnychuk, MC, Payne, JC, Rechisky, EL, Porter, AD, Jackson, GD, Ward, 
  #BR, Vincent, SP, Semmens, J (2011) In situ measurement of coastal ocean movements 
  #and survival of juvenile Pacific salmon. Proceedings of the National Academy of 
  #Sciences 108:8708-8713.

#Survival function
Chilko2012V7.models=function()
{
  #---Phi---
  Phi.time=list(formula=~-1+time)
  Phi.time.treatment=list(formula=~-1+time:treatment_type)
  
  #---p---
  pQCS.indices=as.numeric(row.names(Chilko2012V7.ddl$p[Chilko2012V7.ddl$p$time==7,]))
  pQCS.values=rep(0.67,length(pQCS.indices))
  
  p.time=list(formula=~-1+time, fixed=list(index=pQCS.indices, value=pQCS.values))
  
  cml=create.model.list("CJS")
  results=mark.wrapper(cml,data=Chilko2012V7.process,ddl=Chilko2012V7.ddl,
                       adjust=TRUE, realvcv=TRUE, profile.int=T,chat=1.43)
  return(results)
}


#Run models
Chilko2012V7.results=Chilko2012V7.models()

#Return results
Chilko2012V7.results

#Adjust parameter count  
Chilko2012V7.results[[1]]<-adjust.parameter.count(Chilko2012V7.results[[1]],11)
Chilko2012V7.results[[2]]<-adjust.parameter.count(Chilko2012V7.results[[2]],16)

#Recalculate table
Chilko2012V7.results$model.table=model.table(Chilko2012V7.results)

#Return results
Chilko2012V7.results



###############################################
##Output the segment survivals to csv
###############################################
order<-c(1,2)

for(i in 1:length(order))
{
model<-order[i]
#The confidence intervals were corrected for c-hat using the chat argument in the
#mark.wrapper function.
write.table(Chilko2012V7.results[[model]]$model.name, "real_V7.csv", sep=",",
            append=TRUE, quote=TRUE)
write.table(Chilko2012V7.results[[model]]$results$real, "real_V7.csv", sep=",",
            append=TRUE, quote=FALSE)

#To obtain standard error estimates corrected for c-hat, used the compute.real function.
write.table(compute.real(Chilko2012V7.results[[model]], data=Chilko2012V7, se=TRUE),
            "real_V7.csv", sep=",", append=TRUE, quote=FALSE)
}

#Finally, output the vcv matrix for use calculating the error in the cumulative
#survival estimates via the delta method.
write.table(compute.real(Chilko2012V7.results[[2]], data=Chilko2012V7, vcv=TRUE),
            "realVCmatrix.csv", sep=",", quote=FALSE)