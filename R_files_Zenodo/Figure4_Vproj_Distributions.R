# Figure 4
# distribution of Vproj for different evolutionary axes


par(mfrow=c(1,3),
    mar = c(6,4,4,1))

plot(NA, xlab = "",ylab = "", main = "Ruminantia", xlim = c(0.01,0.15), ylim = c(-0.5,9.5), cex.main = 3, axes = F)
lines(c(0,.15),c(0,0))
meds <- apply(vproj.df,2,median)
for(i in 1:9){
  lines(c(meds[i],meds[i]),c(-0.2,0.2), lwd = 6)
  lines(c(meds[i],meds[i]),c(-0.2,0.2),
      col = c(palette("Okabe-Ito"),"gray50")[i], lwd = 3)}
#points(cbind(apply(vproj.df,2,median),0), pch = 21, bg = c(palette("Okabe-Ito"),"gray50"), cex = 3, lwd = 2)
vioplot(vproj.df, col = c(palette("Okabe-Ito"),"gray50"),
        horizontal = TRUE,  yaxt= 'n', add = T)
axis(2, cex.axis = 1.6, labels = c("Medians",names(vproj.df)), at = 0:9)
axis(1, cex.axis = 1.6, at = c(0.02,0.06,0.1,0.14))


par(mar = c(6,1,4,1))
plot(NA, xlab = "",ylab = "", sub = "Projected Variance", cex.sub = 3,main = "Cervidae", xlim = c(0.01,0.15), ylim = c(-0.5,9.5), cex.main = 3, axes = F)
lines(c(0,.15),c(0,0))
meds <- apply(vproj.df.cervidae,2,median)
for(i in 1:9){
  lines(c(meds[i],meds[i]),c(-0.2,0.2), lwd = 6)
  lines(c(meds[i],meds[i]),c(-0.2,0.2),
      col = c(palette("Okabe-Ito"),"gray50")[i], lwd = 3)}
vioplot(vproj.df.cervidae, col = c(palette("Okabe-Ito"),"gray50"),
        horizontal = TRUE,  yaxt= 'n', add = T)
axis(1, cex.axis = 1.6, at = c(0.02,0.06,0.1,0.14))

par(mar = c(6,1,4,1))
plot(NA, xlab = "", main = "Bovidae",  xlim = c(0.01,0.15), ylim = c(-0.5,9.5), cex.main = 3, axes = F)
lines(c(0,.15),c(0,0))
meds <- apply(vproj.df.bovidae,2,median)
for(i in 1:9){
  lines(c(meds[i],meds[i]),c(-0.2,0.2), lwd = 6)
  lines(c(meds[i],meds[i]),c(-0.2,0.2),
      col = c(palette("Okabe-Ito"),"gray50")[i], lwd = 3)}
vioplot(vproj.df.bovidae, col = c(palette("Okabe-Ito"),"gray50"),
        horizontal = TRUE, add = TRUE, yaxt = 'n')
axis(1, cex.axis = 1.6, at = c(0.02,0.06,0.1,0.14))





