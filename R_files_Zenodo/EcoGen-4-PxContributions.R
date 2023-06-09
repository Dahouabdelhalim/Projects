###################################################################################################
########      Assignment of allopolyploid haplotypes to their corresponding subgenomes      #######
###################################################################################################


########################
###  Load libraries  ###
library(stringr)
library(sjmisc)


###########################################
###  Set working directories and paths  ###


workDir = "PATH/to/workingDirectory/"
infile = "PATH/to/file/EcoGen-4-PxContributions.txt"
setwd(workDir)



#################################
###       Load datasets       ###


sp = c("Ca","Co","Ta","Um","Cr","Cy","Ge","Tr")

loc = c("3Piso","AcylcoA","Aldolase","BRCA1","CENP_E","Cyclin","DMC1","E2F","eif3k","elongTS","galac","Glusyn","GTP","HS","KinMot","lipase","mtPorin","PCNA_2L","PHY_C","RFC1","Ribprot","RING","Sec24","SMC4","sucP","TransFac","vesATPase","VHS_GAT","waxy","zinc")


# Table of individual genotypes
df = read.table(infile, header=TRUE, sep="\\t", colClasses="character")
colnames(df)[which(colnames(df)=="X3Piso")] = "3Piso"
df[,loc] = apply(df[,loc],2,function(x) str_pad(x,3,pad="0"))
df[,loc] = apply(df[,loc],2,function(x) sub(pattern="00-",replacement="-",x=x))


# Assign polyploid haplotypes to their respective subgenomes based on the table (taking into account gene trees branch lengths for derived/lost haplotypes)

res = NULL

for(j in loc){
  nhap = unique(df[,j])
  
  for(i in sp){
    s = df[which(df$sp==i),] # subset at sp i at loc j
    nind = nrow(s) #nb of ind within sp i (accounting for ploidy level)
    
    for(n in nhap){
      if(n!="-"){
        x = sum(s[,j]==n)/nind # % of hap n in sp 1 at loc j
        res = rbind(res,cbind(loc=j,sp=i,hap=n,val=x,contrib=NA))
      }else{
        x1 = nrow(s[which(s[,j]==n & s$AllopN=="A"),])/nind
        x2 = nrow(s[which(s[,j]==n & s$AllopN=="B"),])/nind
        res = rbind(res,cbind(loc=j,sp=i,hap=n,val=x1,contrib="na1"))
        res = rbind(res,cbind(loc=j,sp=i,hap=n,val=x2,contrib="na2"))
      }
    }
  }
}

res = data.frame(res,stringsAsFactors=F)


for(i in c("Cr","Cy","Ge","Tr")){
  
  if(i=="Cr"){       p1 = "Ta" ; p2 = "Co"
  }else if(i=="Cy"){ p1 = "Ta" ; p2 = "Ca"
  }else if(i=="Ge"){ p1 = "Um" ; p2 = "Co"
  }else if(i=="Tr"){ p1 = "Um" ; p2 = "Ca"
  }
  
  nind = which(res$sp==i & res$val>0) # rows of all haplotypes present in Px
  
  for(n in nind){
    
    lc = res[n,"loc"] ; hap = res[n,"hap"]
    nsubG = df[which(df$sp==i & df[,lc]==hap),"AllopN"]
    
    if(hap=="-"){
    }else{
      
      hap.fq1 = res[which(res$sp==p1 & res$hap==hap & res$loc==lc),"val"]
      hap.fq2 = res[which(res$sp==p2 & res$hap==hap & res$loc==lc),"val"]
      
      if(hap.fq1>0 & hap.fq2==0)       { res[n,"contrib"] = "p1"  # add: p1 hap
      }else if(hap.fq1==0 & hap.fq2>0) { res[n,"contrib"] = "p2"  # add: p2 hap
      }else if(hap.fq1>0 & hap.fq2>0)  { res[n,"contrib"] = "p12" # add: p1/p2 hap
      }else if(hap.fq1==0 & hap.fq2==0){
        
        if("A"%in%nsubG & !("B"%in%nsubG))          { res[n,"contrib"] = "n1"  # new: p1
        }else if(!("A"%in%nsubG) & "B"%in%nsubG)    { res[n,"contrib"] = "n2"  # new: p2
        }else if("A"%in%nsubG & "B"%in%nsubG)       { res[n,"contrib"] = "n12" # new: p1/p2
        }else if(!("A"%in%nsubG) & !("B"%in%nsubG)) { res[n,"contrib"] = "na"  # new: unknown
        }
      }
      
    }
    
  }
}

res.px = res[which(res$sp %in% c("Cr","Cy","Ge","Tr") & res$val>0),]
res.px$val = as.numeric(res.px$val)



# Table of counts

resf = matrix(nrow=4,ncol=8,dimnames=list(c("Cr","Cy","Ge","Tr"),c("p1","p2","p12","n1","n2","n12","na1","na2")))

for(i in c("Cr","Cy","Ge","Tr")){
  s.px = res.px[which(res.px$sp==i),]
  resf[i,"p1"] = sum(s.px[which(s.px$contrib=="p1"), "val"])
  resf[i,"p2"] = sum(s.px[which(s.px$contrib=="p2"), "val"])
  resf[i,"p12"] = sum(s.px[which(s.px$contrib=="p12"), "val"])
  resf[i,"n1"] = sum(s.px[which(s.px$contrib=="n1"), "val"])
  resf[i,"n2"] = sum(s.px[which(s.px$contrib=="n2"), "val"])
  resf[i,"n12"] = sum(s.px[which(s.px$contrib=="n12"), "val"])
  resf[i,"na1"] = sum(s.px[which(s.px$contrib=="na1"), "val"])
  resf[i,"na2"] = sum(s.px[which(s.px$contrib=="na2"), "val"])
}

resf = as.data.frame(resf)
resf ; rowSums(resf)

nresf = resf
nresf$p1 = nresf$p1 + nresf$p12/2
nresf$p2 = nresf$p2 + nresf$p12/2
nresf$n1 = nresf$n1 + nresf$n12/2
nresf$n2 = nresf$n2 + nresf$n12/2
nresf = nresf[,-which(names(nresf) %in% c("p12","n12"))]
nresf ; rowSums(nresf)



# Table of proportions

nresf.px = nresf
nresf.px[] = lapply(nresf, function(x) x/30)
nresf.px ; rowSums(nresf.px)













