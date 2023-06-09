# Figure 6
# relationship between Divergence magnitude and direction
# when species diverge from their ancestor more closely aligned to CREA, they diverge farther




# the plot:
par(mfrow=c(1,3), mar = c(5,6,4,2))
plot(NA, xlim = range(div.crea.angle.sf),ylim = range(subfam.divergence),
     xlab = "", ylab = "Morphological Divergence",
     cex.lab = 3, main = "Ruminantia", cex.main = 3, cex.axis = 2)
box()
points(div.crea.angle.sf,subfam.divergence,pch=pts.pch,bg=pts.col,
     cex = rescale.numeric(R.df$csize,c(1,3)))
abline(lm(subfam.divergence~div.crea.angle.sf), col = "gray", lwd = 2)
# legend("topright", c("Tragulidae","Cervidae","Moschidae","Bovidae", "Random Angles"), col = 'black',
#        pt.lwd = 0.5, pt.cex = 3,pt.bg = c(family.colors[c(3,1,4,2)], "gray75"), pch = c(24,21,25,22,22),title.adj = -1, cex = 2)
legend("topright", c("Tragulidae","Cervidae","Moschidae","Bovidae"), col = 'black',
       pt.lwd = 0.5, pt.cex = 3,pt.bg = c(family.colors[c(3,1,4,2)]), pch = c(24,21,25,22),title.adj = -1, cex = 2)

# family-specific, colored by subfamily
par(mar=c(5,2,4,2))
plot(div.crea.angle.sf[which(rt[,1]=="Cervidae")],
     subfam.divergence[which(rt[,1]=="Cervidae")],
     pch=21,
     bg=as.factor(rt[which(rt[,1]=="Cervidae"),2]),
     main = "Cervidae",
     cex.main = 3,
     cex = rescale.numeric(cervid.size,c(1.5,3.5)),
     xlab = "Divergence~CREA Angle",
     cex.lab = 3,
     axes = F)
box() ; axis(1, cex.axis = 2)
legend("topright",
       legend = sort(as.character(unique(as.factor(rt[which(rt[,1]=="Cervidae"),2])))),
       fill = CBP,
       cex = 2)
for(i in sort(as.character(unique(as.factor(rt[which(rt[,1]=="Cervidae"),2]))))){
  COLS <- setNames(CBP,sort(as.character(unique(as.factor(rt[which(rt[,1]=="Cervidae"),2])))))
  if(length(which(rt[,2]==i)) > 1){
  X<-subfam.divergence[which(rt[,2]==i)] ; Y <- div.crea.angle.sf[which(rt[,2]==i)]
  f <- lm(X~Y)
  x <- Y
  y <- predict(f,data.frame(x))
  lines(x,y,col = COLS[i], lwd =2)}
}


par(mar=c(5,2,4,4))
plot(div.crea.angle.sf[which(rt[,1]=="Bovidae")],
     subfam.divergence[which(rt[,1]=="Bovidae")],
     pch=22,
     bg=as.factor(rt[which(rt[,1]=="Bovidae"),2]),
     main = "Bovidae",
     cex.main = 3,
     cex = rescale.numeric(cervid.size,c(1.5,3.5)),
     xlab = "", ylab = "",
     cex.lab = 3,
     axes = F)
box() ; axis(1, cex.axis = 2)
legend("topright", legend = sort(as.character(unique(as.factor(rt[which(rt[,1]=="Bovidae"),2])))), fill = CBP, cex = 2, bg = "transparent")
for(i in sort(as.character(unique(as.factor(rt[which(rt[,1]=="Bovidae"),2]))))){
  COLS <- setNames(CBP,sort(as.character(unique(as.factor(rt[which(rt[,1]=="Bovidae"),2])))))
  if(length(which(rt[,2]==i)) > 1){
  X<-subfam.divergence[which(rt[,2]==i)] ; Y <- div.crea.angle.sf[which(rt[,2]==i)]
  f <- lm(X~Y)
  x <- Y
  y <- predict(f,data.frame(x))
  lines(x,y,col = COLS[i], lwd =2)}
}



# there IS a positive relationship between alignment of divergence.traj~CREA and the magnitude of divergence
# in other words, species that diverge from their ancestor in the direction of CREA tend to diverge farther than those diverging elsewhere
anova(traj.crea.GLS <- gls(subfam.divergence~div.crea.angle.sf,correlation=corBrownian(phy=R.df$tree)))


# family specific
dacC <- div.crea.angle.sf[which(rt[,1]=="Cervidae")] ; sfdC <- subfam.divergence[which(rt[,1]=="Cervidae")]
anova(traj.crea.GLS.C <- gls(sfdC~dacC,correlation=corBrownian(phy=cervid.tree)))

dacB <- div.crea.angle.sf[which(rt[,1]=="Bovidae")] ; sfdB <- subfam.divergence[which(rt[,1]=="Bovidae")]
anova(traj.crea.GLS.C <- gls(sfdB~dacB,correlation=corBrownian(phy=bovid.tree)))


# subfamilies
big.subfams <- names(which(table(rt[,2]) > 9))

for(i in 1:length(big.subfams)){
  part <- which(rt[,2] == big.subfams[i])
  name <- names(part)
  angle <- div.crea.angle.sf[name]
  div <- subfam.divergence[name]
  tree <- keep.tip(R.df$tree, tip= name)
  GLS <- gls(div~angle, correlation = corBrownian(phy = tree))
  print(big.subfams[i])
  print(anova(GLS))
}

for(i in 1:length(big.subfams)){
  part <- which(rt[,2] == big.subfams[i])
  name <- names(part)
  angle <- div.crea.angle.sf[name]
  div <- subfam.divergence[name]
  tree <- keep.tip(R.df$tree, tip= name)
  LM <- lm(div~angle)
  print(big.subfams[i])
  print(anova(LM))
}

