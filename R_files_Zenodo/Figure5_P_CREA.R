# Figure 5
# the relationship between properties of P (as it relates to CREA) and morphological divergence




# Divergence~P plots
sf.div <- subfam.divergence[names(Vproj)]

# together:
par(mfrow=c(1,3),
    mar = c(6,6,2,0.1))
plot(NA,ylim = range(sf.div), xlim = range(Vproj),
     ylab = "Morphological Divergence", xlab = "Projected Variance of CREA", cex.lab = 2)
points(cbind(Vproj,subfam.divergence[names(Vproj)]), bg = pts.col[tree.20$tip.label],
       pch = pts.pch[tree.20$tip.label], cex = rescale.numeric(R.size[tree.20$tip.label], c(1,3)))
legend("topright", c("Tragulidae","Cervidae","Moschidae","Bovidae"), col = 'black',
       pt.lwd = 1, pt.cex = 2,pt.bg = c(family.colors[c(3,1,4,2)]), pch = c(24,21,25,22),title.adj = -1, cex = 1.5)
abline(lm(subfam.divergence[names(Vproj)]~Vproj), col ="gray", lty = 2)

par(mar=c(6,1.9,2,1.9))
plot(NA,ylim = range(sf.div), xlim = range(Vproj.div),
     ylab = "", xlab = "Projected Variance of Species' Divergence", cex.lab = 2, axes = FALSE)
box() ; axis(1)
points(cbind(Vproj.div,subfam.divergence[names(Vproj)]), bg = pts.col[tree.20$tip.label],
       pch = pts.pch[tree.20$tip.label], cex = rescale.numeric(R.size[tree.20$tip.label], c(1,3)))
abline(lm(subfam.divergence[names(Vproj)]~Vproj.div), col ="gray", lty = 1)

par(mar=c(6,0.1,2,6))
plot(NA,ylim = range(sf.div), xlim = range(Vrel),
     ylab = "", xlab = "Morphological Integration (Z Vrel)", cex.lab = 2, axes = FALSE)
box() ; axis(1)
points(cbind(Vrel,subfam.divergence[names(Vproj)]), bg = pts.col[tree.20$tip.label],
       pch = pts.pch[tree.20$tip.label], cex = rescale.numeric(R.size[tree.20$tip.label], c(1,3)))
abline(lm(subfam.divergence[names(Vproj)]~Vrel), col ="gray", lty =2)







# REGRESSIONS OF WITHIN-POPULATION METRICS AND DIVERGENCE
div.vproj.sf <- gls(sf.div~Vproj, correlation = corBrownian(phy=tree.20))
anova(div.vproj.sf)

div.VPROJDIV.sf <- gls(sf.div~Vproj.div, correlation = corBrownian(phy=tree.20))
anova(div.VPROJDIV.sf)

div.evd.sf <- gls(sf.div~Vrel, correlation = corBrownian(phy=tree.20))
anova(div.evd.sf)



# integration~divergence without Odocoileus
sf.div2 <- sf.div[-c(9,10)]
Vrel2 <- Vrel[-c(9,10)]
tree.20.2 <- keep.tip(tree.20, tip = tree.20$tip.label[-c(9,10)])
div.evd.sf2 <- gls(sf.div2~Vrel2, correlation = corBrownian(phy=tree.20.2))
anova(div.evd.sf2)

Vproj2 <- Vproj[-c(9,10)]
div.vproj.sf2 <- gls(sf.div2~Vproj2, correlation = corBrownian(phy=tree.20.2))
anova(div.vproj.sf2)



# just Bovidae
sf.div.bovids <- sf.div[which(rt2[,1] == "Bovidae")]
Vproj.bovids <- Vproj[which(rt2[,1] == "Bovidae")]
bovid.tree.intra <- keep.tip(rownames(rt2)[which(rt2[,1] == "Bovidae")],phy=tree.20)
div.vproj.sf.BOVIDAE <- gls(sf.div.bovids~Vproj.bovids, correlation = corBrownian(phy=bovid.tree.intra))
anova(div.vproj.sf.BOVIDAE)

Vproj.div.bovids <- Vproj.div[which(rt2[,1] == "Bovidae")]
div.vprojDIV.sf.BOVIDAE <- gls(sf.div.bovids~Vproj.div.bovids, correlation = corBrownian(phy=bovid.tree.intra))
anova(div.vprojDIV.sf.BOVIDAE)

evd.bovids <- Vrel[which(rt2[,1] == "Bovidae")]
div.evd.sf.BOVIDAE <- gls(sf.div.bovids~evd.bovids, correlation = corBrownian(phy=bovid.tree.intra))
anova(div.evd.sf.BOVIDAE)


# Cervidae
sf.div.cervids <- sf.div[which(rt2[,1] == "Cervidae")]
Vproj.cervids <- Vproj[which(rt2[,1] == "Cervidae")]
cervid.tree.intra <- keep.tip(rownames(rt2)[which(rt2[,1] == "Cervidae")],phy=tree.20)
div.vproj.sf.CERVIDAE <- gls(sf.div.cervids~Vproj.cervids, correlation = corBrownian(phy=cervid.tree.intra))
anova(div.vproj.sf.CERVIDAE)


VprojDIV.cervids <- Vproj.div[which(rt2[,1] == "Cervidae")]
div.vprojDIV.sf.CERVIDAE <- gls(sf.div.cervids~VprojDIV.cervids, correlation = corBrownian(phy=cervid.tree.intra))
anova(div.vprojDIV.sf.CERVIDAE)


evd.cervids <- Vrel[which(rt2[,1] == "Cervidae")]
div.evd.sf.CERVIDAE <- gls(sf.div.cervids~evd.cervids, correlation = corBrownian(phy=cervid.tree.intra))
anova(div.evd.sf.CERVIDAE)
