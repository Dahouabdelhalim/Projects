library(dartR) 
library(vcfR)
library(poppr)
library(pegas)
library(ape)
library(adegenet)
library(ade4)
library(genepopedit)
library(geoR)
library(radiatior)
library(spdep)
library(adehabitat)
library(splancs)
library(raster)
library(poppr)
library(geoR)
library(varhandle)
library(logisticPCA)
library(geosphere)
library(caret)
library(fossil)
library(vegan)
library(BEDASSLE)
library(hierfstat)
library(broom)
library(corrplot) 
library(BAMMtools)

###get slopes for GBS data
dist_gen_geo<-function(locality,genetic)

{
setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S4_Gemsnake_pop_localities")

read.table( locality,header=T,fill=TRUE)->lat
lat[,c(1,4,3)]->latX
names(latX)<-c("Sample", "lat", "lon")


setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S2_VCF_Files_Gemsnakes")
read.vcfR(genetic  )->vmaf

##this will make the dataset only biallelic
vcfR2genlight(vmaf)->gen
vcfR2DNAbin(vmaf,consensus=TRUE,extract.haps = FALSE)->genX

gl.filter.allna(gen, by.pop = FALSE, verbose = NULL)->gen2

gl2gi(gen2)->geni

join_gen_lat<-function(obj,lat)
  
{ names(obj$tab[,1])->n2
  lat[match(n2, lat$Sample),]->m
  data.frame(m$lon,m$lat)->coor

  obj@other$xy<-coor
  data<-data.frame(indNames(obj),coor)
  names(data)<-c("Samples","long","lat")
  
  results<-list()
  results$obj<-obj
  results$loc<-data
  return(results)
}

join_gen_lat(geni,latX)->end

end$obj->obj
obj@other<-NULL

gi2gl(obj)->obj



earth.dist(data.frame(end$loc$long,end$loc$lat))->Gdist
Gdist[which(Gdist==0)]<-1
as.genpop(end$obj@tab)->genP
dist.genpop(genP)->genDist
mean(genDist)->meanGenDist
sd(genDist)->sdGenDist
mean(Gdist)->meanGeoDist
sd(Gdist)->sdGeoDist

plot(log(Gdist),genDist, pch=19,main= genetic,cex=1.2)
lm(as.numeric(genDist)~log(as.numeric(Gdist)))->reg
reg$coefficients[2]->genDistSlope
abline(reg,col="darkred",lwd=3)
cor.test(genDist,log(Gdist))$p.value->p.value_GenDistSlope

dist(end$obj@tab)->Eucldis
mean(Eucldis)->Mean_Eucldis
sd(Eucldis)->sdEucldis
mean(Gdist)->meanGeoDist
sd(Gdist)->sdGeoDist

plot(log(Gdist),Eucldis, pch=19,main= genetic,cex=1.2)
lm(as.numeric(Eucldis)~log(as.numeric(Gdist)))->regN
regN$coefficients[2]->EuclDistSlope
abline(regN,col="darkred",lwd=3)
cor.test(Eucldis,log(Gdist))$p.value->p.value_EuclSLope

pop(obj)<-obj@ind.names
gl2gi(obj)->obj2


##function from https://dfzljdn9uc3pi.cloudfront.net/2018/5089/1/Article_S3.pdf
get_pairwise_Fst <-function(genmat, grp) {
  grp = factor(grp)
  if (length(levels(grp)) > 1) {
    allele.counts = apply(genmat@tab, 2, function(x) {
      pops = grp
      counts = tapply(x, pops, function(y) {
        sum(y, na.rm = T)
      })
      return(counts)
    })
    sample.sizes = apply(genmat@tab, 2, function(x) {
      pops = grp
      sizes = tapply(x, pops, function(y) {
        2 * sum(!is.na(y))
      })
      return(sizes)
    })
    Fst = calculate.all.pairwise.Fst(allele.counts, sample.sizes)
    colnames(Fst) = rownames(Fst) = levels(grp)
    # print(Fst)
    return(tidy(as.dist(Fst, upper = FALSE)))
  } else {
    print("Only one pop")
    return(matrix())
  }
}

get_pairwise_Fst(obj2,obj2@pop)->fstD
fstD$distance->fstD
fstD/(1-fstD)->fstDx
plot(log(Gdist),fstDx,pch=19,main= genetic,cex=1.2)
lm(as.numeric(fstDx)~log(as.numeric(Gdist)))->reg2
reg2$coefficients[2]->one_minus_fst_slope
abline(reg2,col="darkred",lwd=3)
cor.test(fstDx,log(Gdist))$p.value->p.value_Fst_Slope
mean(fstD)->meanfstD
sd(fstD)->sdfst

length(gen@ind.names)->n
 nuc.div(genX)->pi
 #tajima.test(genX)$D->tajima

data.frame(n,one_minus_fst_slope,p.value_Fst_Slope, meanfstD,sdfst,genDistSlope,p.value_GenDistSlope,meanGenDist,sdGenDist,meanGeoDist,sdGeoDist,EuclDistSlope,Mean_Eucldis,sdEucldis, p.value_EuclSLope,pi)->data_fin

return(data_fin)
}

##make sure to properly set your working directies to the folders with locality and VCF data---see the two setwd() in the function above set to the two tables at:
###setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S4_Gemsnake_pop_localities")
###setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Datat/Data_S2_VCF_Files_Gemsnakes")
setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S4_Gemsnake_pop_localities")
dir()->l
setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S2_VCF_Files_Gemsnakes")
dir()->v

lapply(c(1:36), function(x) dist_gen_geo (l[x],v[x]) )->out

##format the output into a table
names(out)<-v
do.call("rbind", out)->GBS_Results



###get genetic distances and other metrics for mtDNA

mtDist<-function(locality, mtgenetic)
  
{
  
  setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S4_Gemsnake_pop_localities")
  read.table( locality,header=T,fill=TRUE)->lat
  lat[,c(1,4,3)]->latX
  names(latX)<-c("Sample", "lat", "lon")
  
  setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S3_mtDNA_Files_Gemsnakes")
  read.FASTA(mtgenetic)->mt
  
  names(mt)<- gsub("\\t","",names(mt))
  names(mt)<-sub("^[^_]*_[^_]*_", "\\\\1", names(mt))
  
  latX$Sample<-gsub("_.*","\\\\1",latX$Sample)
  
  latX[match(names(mt), latX$Sample),]->m
  
  dist.dna(mt,pairwise.deletion=TRUE, as.matrix=FALSE)->Gendistmt
  earth.dist(data.frame(m$lon,m$lat))->Geomtdist
  
  Geomtdist[which(Geomtdist==0)]<-1
  
  plot(log(Geomtdist),Gendistmt, pch=19,main= locality,cex=1.2)
  lm(Gendistmt~log(Geomtdist))->regNN
  regNN$coefficients[2]->mtDistSlope
  abline(regNN,col="darkred",lwd=3)
  cor.test(Gendistmt,log(Geomtdist))$p.value->p.value_mtSLope
  nuc.div(mt)->pimtDNA
  
  mean(Gendistmt,na.rm=TRUE)->meandistmtDNA
  
  sd(Gendistmt,na.rm=TRUE)->sddistmtDNA
  
  length(mt)->number
  
  
  data.frame(number, mtDistSlope, p.value_mtSLope,meandistmtDNA,sddistmtDNA,pimtDNA)->mtDNA_Final
  
  return(mtDNA_Final)}

##again for the mtDNA, set the working directories correctly
##setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S4_Gemsnake_pop_localities")
##setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S3_mtDNA_Files_Gemsnakes")
setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S4_Gemsnake_pop_localities")
dir()->L

setwd("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S3_mtDNA_Files_Gemsnakes")
dir()->mtg


lapply(c(1:length(L)), function(x) tryCatch(mtDist(L[x], mtg[x]),error=function(e)(NA)))->tabmt

names(tabmt)<-L
do.call("rbind", tabmt)->mtDNA_Results

data.frame(row.names(GBS_Results),GBS_Results,mtDNA_Results)->Pop_Structure_Results
names(Pop_Structure_Results)[1]<-"Pop_Name"


##use clads in Julia and get output here
###getting estimates from clads

isfar2 <- get(load('~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S8_Clads_output'))
data.frame(isfar2$tree$tip.label,isfar2$lambdatip_map)->clads_Results
names(clads_Results)<-c("Tree_name", "cladsNew"  )


##use BAMM and get output here
###getting estimates from BAMM

ed <- getEventData("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S11_Tree_for_sp_rates_estimates.tre" , "~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S9_BAMM_output.txt" )
plot.bammdata(ed, lwd=2, labels = T, cex = 0.5)->gemR
addBAMMshifts(gemR, cex=2)
addBAMMlegend(gemR)

data.frame(ed$tip.label, ed$meanTipLambda)->BAMM_Results
names(BAMM_Results)<-c("Tree_name", "BAMM"  )


###get DR estimates from Singhal here: https://github.com/singhal/brazil_IBD/blob/main/analyses/speciesSpecificDivRate.R
jetzDivRates <- function(tree) {
  
  spRate <- function(sp, tree) {
    #get branch lengths from root to tip
    edges <- vector()
    daughterNode <- match(sp, tree$tip.label)
    while (daughterNode != (length(tree$tip.label) + 1)) {
      parentNode <- tree$edge[which(tree$edge[,2] == daughterNode), 1]
      edges <- c(edges, tree$edge.length[which(tree$edge[,1] == parentNode & tree$edge[,2] == daughterNode)])
      daughterNode <- parentNode
    }
    
    res <- sum(sapply(1:length(edges), function(x) edges[x] * (1/(2 ^ (x-1)))))
    res <- res ^ (-1)
    
    return(res)
  }
  
  rates <- unlist(lapply(tree$tip.label, function(x) spRate(x, tree)))
  names(rates) <- tree$tip.label
  
  return(rates)
  
}

read.tree("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S11_Tree_for_sp_rates_estimates.tre")->tree

jetzDivRates(tree)->DR_Results
data.frame(names(DR_Results),DR_Results)->DR_Results
rownames(DR_Results)<-NULL
names(DR_Results)<-c("Tree_name", "DR"  )

merge(clads_Results, BAMM_Results,by="Tree_name")->M1
merge(M1,DR_Results, by="Tree_name")->Tip_States

##plot their rates

###Making some graphs

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

t_col("green",50)->gr
t_col("blue",50)->bl
t_col("red",50)->rd
# First distribution
hist(DR_Results$DR,breaks=20, xlim=c(0,.8), col=rd)

# Second with add=T to plot on top
hist(BAMM_Results$BAMM, breaks=20, xlim=c(0,.8), col ="black", add=T)

hist(clads_Results$cladsNew,breaks=20,col=bl,add=T)

###you can merge the Tip_States with Pop_Structure_Results using the conversion bridging Tree_name in the former with Pop_Name in the later from the following table
###or of course just use this table that has all of the population and speciation rate metrics in it already
read.table("~/Dropbox/Gemsnakes_div_st/Supplementary_Data/Data_S1_Final_Data.txt"   ,header=T)->table



###plot correlations among metrics
#w/out mtDNA
table[c(4:28)]->res

cor(res,use="pairwise.complete.obs")->res2

testRes = cor.mtest(res2, conf.level = 0.95)

corrplot(res2, method = "circle", shade.col = NA, order="original",type="upper")







