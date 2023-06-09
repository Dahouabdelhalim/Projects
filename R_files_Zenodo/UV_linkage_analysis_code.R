library(qtl)


## read in data 
fullset <- read.cross('csvr','.', 'LinkageMap.txt', genotypes=NULL, estimate.map=FALSE, crosstype="4way", convertXdata=TRUE)

 ### THERE IS A BUG meaning the map positions get imported incorrectly. Do this to fix: 
for(i in seq_along(fullset$geno)) {
   fullset$geno[[i]]$map <- rbind(fullset$geno[[i]]$map, fullset$geno[[i]]$map)
 }
 
fullset <- jittermap(fullset)
fullset <- calc.genoprob(fullset, step=0.1)

males <- fullset[fullset$pheno$sex==0]

out.males <- scanone(males, model='binary',pheno.col='uv')
out.perm <- scanone(males, model='binary', n.perm=1000,pheno.col='uv')
a <- summary(out.perm)[1]
plot(out.males)
abline(h=a,col='red')
