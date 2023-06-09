setwd("")
library(vegan)
da = read.csv("Families of cooccurring groups.csv",row.names = 1)
da = da[da$Family!="unk" & da$Family!="No hit",]
da = droplevels(da)

# TEST 1: ARE CERTAIN FAMILIES OVER- OR UNDERREPRESENTED IN THE ASSOCIATED GROUPS
nperm = 10000
ta = table(da$Family)
sel = which(ta>4)
fams = names(ta)[sel]
nsel = length(sel)
nn = rep(NA,nsel)
na = rep(NA,nsel)
neg = rep(NA,nsel)
pos = rep(NA,nsel)
for(i in 1:nsel){
   po = which(da$Family == fams[i])
   npo = length(po)
   test = sum(da$Group[po]=="No group")
   rtest = rep(NA,nperm)
   for(j in 1:nperm){
      rtest[j] = sum(sample(da$Group,npo)=="No group")
   }
   nn[i] = npo
   na[i] = npo-test
   neg[i] = mean(test<=rtest)
   pos[i] = mean(test>=rtest)
}
res = data.frame(nn,na,neg,pos)
rownames(res) = fams
colnames(res) = c("n","na","Pr(na<=null)","Pr(na>=null)")
write.csv(res,file="permutation test 1.csv")

# TEST 2: ARE GROUPS 1 & 2 PHYLOGENETICALLY SIMILAR OR DISSIMILAR?; ARE GROUPS 3 & 4 PHYLOGENETICALLY SIMILAR OR DISSIMILAR?
res = matrix(NA,2,4)
for(z in 1:2){
   if(z==1){
      t1 = "Group 1"
      t2 = "Group 2"
      
   } else {
      t1 = "Group 3"
      t2 = "Group 4"      
   }
   da1 = droplevels(da[da$Group == t1 | da$Group == t2,])
   fams = levels(da1$Family)
   nfams = length(fams)
   ta = matrix(NA,2,nfams)
   for(i in 1:nfams){
      ta[1,i] = sum(da1[da1$Group==t1,]$Family==fams[i])
      ta[2,i] = sum(da1[da1$Group==t2,]$Family==fams[i])
   }
   test = vegdist(ta,method="bray")
   rtest = rep(NA,nperm)
   for(j in 1:nperm){
      rgroup = sample(da1$Group)
      for(i in 1:nfams){
         ta[1,i] = sum(da1[rgroup==t1,]$Family==fams[i])
         ta[2,i] = sum(da1[rgroup==t2,]$Family==fams[i])
      }
      rtest[j] = vegdist(ta,method="bray")
   }
   res[z,] = c(test,mean(rtest),mean(test<=rtest),mean(test>=rtest))
}
colnames(res) = c("BC","E[BC]","Pr(BC<=null)","Pr(BC>=null")
rownames(res) = c("Group 1 vs Group 2","Group 3 vs Group 4")
write.csv(res,file="permutation test 2.csv")
