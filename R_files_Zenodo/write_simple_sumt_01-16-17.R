# write simple sumt

# The standard MrBayes output format for summary trees is difficult to read into and depict in R with support values at nodes. 
# This script creates a NEXUS file with a MrBayes block that would load the output of existing MrBayes analyses, 
# then run the 'sumt' command with the correct formatting for use in R, thus quickly generating all desired summary trees in the correct format.


#change to dir for MB files
setwd("C:/dave/workspace/mrbayes")

#get files
files<-list.files()    
#get tree files
treeFiles<-files[grep(x=files,pattern="run..t")]  
#get the unique identifiers
inputs<-sort(unique(sub(".run..t","",treeFiles)))
#get the nexus files
inputs<-sort(unique(sub(".tree.","",inputs)))

# this should catch partitioned unlinked analyses 
# (with more than one tree) as separate instances
# for calculating consensus trees as the rest of the nexus file is read...

#now start constructing sumt lines
bayesBlock<-"#NEXUS"
#which are run I
isI<-grep("^.{3}I", inputs)
for(i in (1:length(inputs))[-isI]){
	origNex<-scan(inputs[i],what="character",sep="\\n",blank.lines.skip=FALSE)
	cropNex<-origNex[2:grep(origNex,pattern="mcmcp")]
	#remove log line
	cropNex<-cropNex[-grep(cropNex,pattern="log start")]
	sumtLine<-paste0("sumt filename=",inputs[i]," outputname=",
		paste0(inputs[i],".simple")," conformat=simple;")
	bayesBlock<-c(bayesBlock,cropNex,sumtLine,"","end;","")
	}

#now run I (simulations)
	#find the runI nexus file
nexFiles<-files[grep(x=files,pattern=".nex$")]  
runINexus<-nexFiles[grep(x=nexFiles,pattern="runI")]
origNex<-scan(runINexus,what="character",sep="\\n",
	blank.lines.skip=FALSE)
cropNex<-origNex[2:grep(origNex,pattern="mcmcp")[1]]
#remove log line
cropNex<-cropNex[-grep(cropNex,pattern="log start")]
#now generate lines for all inputs
for(i in (1:length(inputs))[isI]){
	sumtLine<-paste0("sumt filename=",inputs[i]," outputname=",
		paste0(inputs[i],".simple")," conformat=simple;")
	bayesBlock<-c(bayesBlock,cropNex,sumtLine,"","end;","")
	}

write(bayesBlock,file="simple_sumt.nex")


