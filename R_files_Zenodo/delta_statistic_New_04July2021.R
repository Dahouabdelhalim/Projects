library(ape)
source("code.R.R")

tree <- read.tree("C_tree.nwk")
plot(tree)


##slender
slender <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

delta_slender <- delta(slender,tree,0.1,0.0589,10000,10,100)


delta_slender

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(slender)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_slender <- sum(random_delta>delta_slender)/length(random_delta)
boxplot(random_delta)
abline(h=delta_slender,col="red")

p_value_slender


##appear
appear <- c(0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0)

delta_appear <- delta(appear,tree,0.1,0.0589,10000,10,100)


delta_appear

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(appear)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_appear <- sum(random_delta>delta_appear)/length(random_delta)
boxplot(random_delta)
abline(h=delta_appear,col="red")

p_value_appear


##rachis
rachis <- c(0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0)

delta_rachis <- delta(rachis,tree,0.1,0.0589,10000,10,100)


delta_rachis

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(rachis)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_rachis <- sum(random_delta>delta_rachis)/length(random_delta)
boxplot(random_delta)
abline(h=delta_rachis,col="red")

p_value_rachis


##bract
bract <- c(0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

delta_bract <- delta(bract,tree,0.1,0.0589,10000,10,100)


delta_bract

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(bract)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_bract <- sum(random_delta>delta_bract)/length(random_delta)
boxplot(random_delta)
abline(h=delta_bract,col="red")

p_value_bract


##density
density <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,0,0,1,0,1,1,0)

delta_density <- delta(density,tree,0.1,0.0589,10000,10,100)


delta_density

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(density)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_density <- sum(random_delta>delta_density)/length(random_delta)
boxplot(random_delta)
abline(h=delta_density,col="red")

p_value_density


##cincinnus
cincinnus <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0)

delta_cincinnus <- delta(cincinnus,tree,0.1,0.0589,10000,10,100)


delta_cincinnus

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(cincinnus)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_cincinnus <- sum(random_delta>delta_cincinnus)/length(random_delta)
boxplot(random_delta)
abline(h=delta_cincinnus,col="red")

p_value_cincinnus


##flowers
flowers <- c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,0,1,0,0,1,0,1,1,0,1,1,1,1,0,1,1,1,0,1,0,1,0,0,0,0,1,0,1,0,1,0,0,1,0,1,0,1,0,0)

delta_flowers <- delta(flowers,tree,0.1,0.0589,10000,10,100)


delta_flowers

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(flowers)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_flowers <- sum(random_delta>delta_flowers)/length(random_delta)
boxplot(random_delta)
abline(h=delta_flowers,col="red")

p_value_flowers


##calyx
calyx <- c(0,1,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0)

delta_calyx <- delta(calyx,tree,0.1,0.0589,10000,10,100)


delta_calyx

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(calyx)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_calyx <- sum(random_delta>delta_calyx)/length(random_delta)
boxplot(random_delta)
abline(h=delta_calyx,col="red")

p_value_calyx


##bending
bending <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

delta_bending <- delta(bending,tree,0.1,0.0589,10000,10,100)


delta_bending

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(bending)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_bending <- sum(random_delta>delta_bending)/length(random_delta)
boxplot(random_delta)
abline(h=delta_bending,col="red")

p_value_bending


##claw
claw <- c(0,0,1,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,1,0,1,0,0,1,0,0,0,0,0,1,1)

delta_claw <- delta(claw,tree,0.1,0.0589,10000,10,100)


delta_claw

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(claw)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_claw <- sum(random_delta>delta_claw)/length(random_delta)
boxplot(random_delta)
abline(h=delta_claw,col="red")

p_value_claw


##labellum
labellum <- c(1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0,1,0,0,1,1,1,1,0,1,1,0,1,0,1,0,0,1,1,1,1,0,1,1,0,0,1,0,0,0,0,0,0,0)

delta_labellum <- delta(labellum,tree,0.1,0.0589,10000,10,100)


delta_labellum

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(labellum)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_labellum <- sum(random_delta>delta_labellum)/length(random_delta)
boxplot(random_delta)
abline(h=delta_labellum,col="red")

p_value_labellum


##blotch
blotch <- c(0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,1,0,0,1,1,1,1,1,1,0,0,1,0,0,0,1,0,1,0,1,0,1,1,0,0,1,0,0,1,1,1,1,0)

delta_blotch <- delta(blotch,tree,0.1,0.0589,10000,10,100)


delta_blotch

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(blotch)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_blotch <- sum(random_delta>delta_blotch)/length(random_delta)
boxplot(random_delta)
abline(h=delta_blotch,col="red")

p_value_blotch


##stigma
stigma <- c(0,0,1,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,1,0,0,0,0,1,0,0,0,1,0,0,0)

delta_stigma <- delta(stigma,tree,0.1,0.0589,10000,10,100)


delta_stigma

random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(stigma)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value_stigma <- sum(random_delta>delta_stigma)/length(random_delta)
boxplot(random_delta)
abline(h=delta_stigma,col="red")

p_value_stigma
