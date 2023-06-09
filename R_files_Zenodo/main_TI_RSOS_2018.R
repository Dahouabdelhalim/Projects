### Main code for Kamaluddin et al. 2018. Royal Socety Open Science.
# Tsuyoshi Ito
# August 20th, 2018
# Slightly modified on (May 15th, 2019)

library(Morpho)
library(geomorph)
library(car)
library(raster)
library(brms)
library(mice)
library(ggplot2)
library(maps)
library(gridExtra)
library(grid)
library(reshape2)
library(coda)
library(tidyverse)
library(ggbeeswarm)
source("./functions_TI_RSOS_2018.r")

ggplot() + theme_set(theme_bw(base_size = 14, base_family = ""))
par.old <- par(no.readonly = T)


### import and prepare data ---------------------------------------------
# read landmark data
raw <- read.csv.folder2("./coordinates_of_landmarks/")

# estimate missing landmarks
raw$arr[raw$arr == 9999] <- NA
fixed.lm <- fixLMtps(raw$arr)

for (i in fixed.lm$check) {
  cat(c(dimnames(fixed.lm$out)[[3]][i], ":\\t", fixed.lm$checklist[[i]], "\\n"))
}

# find amd remove outliers
outliers <- find.outliers(fixed.lm$out)
clean.lm <-
  fixed.lm$out[, ,!outliers$dist.sort[, 2] %in% SG(outliers$dist.sort[, 2])$dellVal]

# check landmarks
# checkLM(raw$arr, path = "../../surfaces/")
# checkLM(clean.lm, text.lm = TRUE, path = "../../surfaces/")

# import data
data0pre <- read.csv("./data_TI_RSOS_2018.csv")
lm.data <- data.frame(ID = dimnames(clean.lm)[[3]])
data0 <- merge(lm.data, data0pre)

# generalized Procrustes analysis
gpa <-
  procSym(clean.lm, reflect = FALSE, pairedLM = cbind(c(5:16), c(17:28)))

# principal component analysis
pca <- geom_pca(gpa$Sym)

# matrix correlation between Procrustes distances and the Euclid distances of PC scores
rmtx <- r.matrix(gpa$Sym, pca$pc.scores)

# summary of PCA
PC_summary <- data.frame(
  PC = rep(1:40),
  Cor = rmtx,
  Var = pca$pc.summary$importance[2, 1:40],
  Cum = pca$pc.summary$importance[3, 1:40]
)

write.csv(PC_summary, './results/PC_summary.csv')

p1 <- ggplot(PC_summary)
p1 <-
  p1 + geom_bar(aes(x = PC, y = Cum),
                stat = "identity",
                position = "identity",
                fill = "gray")
p1 <-
  p1 + geom_bar(aes(x = PC, y = Var),
                stat = "identity",
                position = "identity",
                fill = "black")
p1 <- p1 + xlab("PC") + ylab("% variance")
p1 <- p1 + xlim(0, 41)

p2 <- ggplot(PC_summary)
p2 <- p2 + geom_line(aes(x = PC, y = Cor))
p2 <- p2 + geom_point(aes(x = PC, y = Cor))
p2 <- p2 + xlab("PC") + ylab("correlation coefficient")
p2 <- p2 + xlim(0, 41)

p1 <- p1 + labs(title = "(a)")
p2 <- p2 + labs(title = "(b)")

g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))
g1$widths <- g2$widths <- unit.pmax(g1$widths, g2$widths)

quartz(width = 3.54, height = 5.5) # 3.54 or 7.25
grid.arrange(g1, g2, layout_matrix = rbind(1, 2))
dev.copy2eps(file = "./results/pc.summary.eps")
dev.off()


# prepair dataset for statistical analyses
a.data0 <-
  cbind(data0, pca$pc.scores, data.frame(size = log(gpa$size)))
w.data0 <- a.data0[a.data0$environment == "Wild",]


plot(PC4 ~ size, col = as.numeric(environment), pch = as.numeric(age), a.data0)


### The maps of WorldClim and sampling sites ----------------------------
# bioclim
wc <- getData("worldclim", var = "bio", res = 10, path = './WorldClim')
wc <- wc[[c(1, 12)]]
names(wc) <- c("temperature", "precipitation")
wc.values <- extract(wc, w.data0[, c("longitude", "latitude")])
w.data1 <- cbind(w.data0, wc.values)

# for map
wc_2.5 <- getData("worldclim", var = "bio", res = 2.5, path = './WorldClim')
quartz(width = 7.25, height = 3)

par(
  mfrow = c(1, 2),
  mar = c(3, 3, 2, 5),
  mgp = c(1.7, 0.5, 0)
)

plot(
  wc_2.5[[1]],
  xlim = c(125, 145),
  ylim = c(30, 45),
  xlab = "longitude (°)",
  ylab = "latitude (°)"
)
mtext("(a)",
      font = 1,
      adj = 0,
      line = 0.2)
points(
  w.data1$longitude,
  w.data1$latitude,
  bg = "black",
  pch = 21,
  cex = 0.5
)

plot(
  wc_2.5[[3]],
  xlim = c(125, 145),
  ylim = c(30, 45),
  xlab = "longitude (°)",
  ylab = "latitude (°)"
)
mtext("(b)",
      font = 1,
      adj = 0,
      line = 0.2)
points(
  w.data1$longitude,
  w.data1$latitude,
  bg = "black",
  pch = 21,
  cex = 0.5
)

par(par.old)
dev.copy2eps(file = "./results/bioclim.eps")
dev.off()


# sampling sites
sampling_site <-
  unique(a.data0[, c("population", "environment", "latitude", "longitude")])
levels(sampling_site$environment) <- c("Captive", "Wild", "Both")
sampling_site[sampling_site$population == "Yakushima", "environment"] <-
  "Both"
m.data <- map_data("world2", "Japan")
p3 <- ggplot()
p3 <-
  p3 + geom_polygon(data = m.data,
                    aes(long, lat, group = subregion),
                    fill = "lightgray") + coord_fixed(1.3)
p3 <- p3 + xlim(128, 146) + ylim(30, 46)
p3 <-
  p3 + geom_point(
    data = data.frame(latitude = 35.38, longitude = 136.94),
    aes(y = latitude, x = longitude),
    size = 4,
    pch = 1
  )
p3 <-
  p3 + geom_point(
    data = sampling_site,
    aes(
      y = latitude,
      x = longitude,
      shape = environment,
      colour = environment
    ),
    size = 3
  )
p3 <-
  p3 + xlab("longitude (°)") + ylab("latitude (°)")

p3 <-
  p3 + geom_segment(
    mapping = aes(
      x = 140,
      y = 32.5,
      xend = 137.1,
      yend = 35.2
    ),
    arrow = arrow(),
    size = 0.5,
    color = "black"
  )
p3 <-
  p3 + annotate("text",
                x = 140,
                y = 32,
                label = "Primate Research Institute, Kyoto University")
p3 <- p3 + theme(legend.title = element_blank())

quartz(width = 7.25, height = 6.25) # 3.54 or 7.25
plot(p3)
dev.copy2eps(file = "./results/sampling.sites.eps", family = "sans")
dev.off()


### A Bayesian linear mixed models for captive–wild test --------------------------------------------------------------------
# to be numeric and scaling
a.data1 <-
  a.data0[, c(
    "ID",
    "population",
    "sex",
    "environment",
    "age",
    "latitude",
    "longitude",
    "captivity",
    "size",
    paste("PC", 1:40, sep = "")
  )]
a.data1[, c("sex", "environment", "age")] <-
  sapply(a.data1[, c("sex", "environment", "age")], as.numeric)
a.data1[, -c(1:2)] <- scale(a.data1[, -c(1:2)])

w.data2 <-
  w.data1[, c(
    "ID",
    "population",
    "sex",
    "environment",
    "age",
    "latitude",
    "longitude",
    "captivity",
    "size",
    "temperature",
    "precipitation",
    paste("PC", 1:40, sep = "")
  )]
w.data2[, c("sex", "environment", "age")] <-
  sapply(w.data2[, c("sex", "environment", "age")], as.numeric)
w.data2[, -c(1:2)] <- scale(w.data2[, -c(1:2)])


# captive–wild test for size
fit_a_size1 <-
  brm(
    size ~ sex + age + captivity + (1 |
                                      population),
    data = a.data1,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.9)
  )

fit_a_size2 <-
  brm(
    size ~ sex + age + (1 | population),
    data = a.data1,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.9)
  )

fit_a_size1_waic <- add_ic(fit_a_size1, ic = c("waic", "r2"))
fit_a_size2_waic <- add_ic(fit_a_size2, ic = c("waic", "r2"))

sink('./results/fit_a_size_summary.txt')
WAIC(fit_a_size1_waic, fit_a_size2_waic)
cat('\\n')
cat('fit_a_size1 R2', mean(fit_a_size1_waic$R2), '\\n')
cat('fit_a_size2 R2', mean(fit_a_size2_waic$R2), '\\n')
sink()  

# plot(fit_a_size1)
sink('./results/fit_a_size1.txt')
print (fit_a_size1) 
sink()  

# plot(fit_a_size2)
sink('./results/fit_a_size2.txt')
print (fit_a_size2) 
sink()  


# captive–wild test for shape
fit_a_shape1 <-
  brm(
    cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) ~ size + sex + age + captivity + (1 |
                                                                                                 population),
    data = a.data1,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.9)
  )

fit_a_shape2 <-
  brm(
    cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) ~ size + sex + age + (1 |
                                                                                     population),
    data = a.data1,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.9)
  )

fit_a_shape1_waic <- add_ic(fit_a_shape1, ic = c("waic", "r2"))
fit_a_shape2_waic <- add_ic(fit_a_shape2, ic = c("waic", "r2"))

sink('./results/fit_a_shape_summary.txt')
WAIC(fit_a_shape1_waic, fit_a_shape2_waic)
cat('\\n')
cat('fit_a_shape1 R2', mean(fit_a_shape1_waic$R2), '\\n')
cat('fit_a_shape2 R2', mean(fit_a_shape2_waic$R2), '\\n')
sink()  
# bayes_R2(fit_a_shape1_waic)
# bayes_R2(fit_a_shape2_waic)

# plot(fit_a_shape1)
sink('./results/fit_a_shape1.txt')
print (fit_a_shape1) 
sink()  

# plot(fit_a_shape2)
sink('./results/fit_a_shape2.txt')
print (fit_a_shape2) 
sink()  



# mcmc plot
mcmc_a_size <- NULL
for (i in 1:4) {
  mcmc_a_size <- rbind(mcmc_a_size, as.mcmc(fit_a_size1)[[i]])
}

mcmc_a_shape <- NULL
for (i in 1:4) {
  mcmc_a_shape <- rbind(mcmc_a_shape, as.mcmc(fit_a_shape1)[[i]])
}

mcmc_a_size_sum <-
  rbind2(
    cbind(data.frame(variable = "size"), rep(NA, 20000)),
    cbind(data.frame(variable = "sex"), mcmc_a_size[, 2]),
    cbind(data.frame(variable = "age"), mcmc_a_size[, 3]),
    cbind(data.frame(variable = "captivity"), mcmc_a_size[, 4]),
    cnames = c("explanatory", "size")
  )

mcmc_a_shape_sum <-
  rbind2(
    cbind(data.frame(variable = "size"), mcmc_a_shape[, c(11, 15, 19, 23, 27, 31, 35, 39, 43, 47)]),
    cbind(data.frame(variable = "sex"), mcmc_a_shape[, c(11, 15, 19, 23, 27, 31, 35, 39, 43, 47) + 1]),
    cbind(data.frame(variable = "age"), mcmc_a_shape[, c(11, 15, 19, 23, 27, 31, 35, 39, 43, 47) + 2]),
    cbind(data.frame(variable = "captivity"), mcmc_a_shape[, c(11, 15, 19, 23, 27, 31, 35, 39, 43, 47) + 3]),
    cnames = c("explanatory", paste("PC", 1:10, sep = ""))
  )

mcmc_a_sum <- cbind(mcmc_a_shape_sum, mcmc_a_size_sum)

p4 <-
  ggplot(
    data = summaryCI(
      melt(mcmc_a_sum),
      measurevar = "value",
      groupvars = c("explanatory", "variable"),
      na.rm = TRUE
    ),
    mapping = aes(x = explanatory, y = mean)
  )
p4 <-
  p4 + geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) + geom_point(position = position_dodge(0.2))
p4 <- p4 + geom_hline(yintercept = 0)
p4 <- p4 + facet_wrap(~ variable, ncol = 5)
p4 <-
  p4 + ylim(-0.85, 0.85) + coord_flip() + ylab("coefficient") + xlab("variables")

quartz(width = 7.25, height = 4)
plot (p4)
dev.copy2eps(file = "./results/stan_cap.eps")
dev.off()




### A Bayesian linear mixed models for ecogeographic variation --------------------------------------------------------------------
# to be numeric and scaling
w.data2 <-
  w.data1[, c(
    "ID",
    "population",
    "sex",
    "environment",
    "age",
    "latitude",
    "longitude",
    "captivity",
    "size",
    "temperature",
    "precipitation",
    paste("PC", 1:40, sep = "")
  )]
w.data2[, c("sex", "environment", "age")] <-
  sapply(w.data2[, c("sex", "environment", "age")], as.numeric)
w.data2[, -c(1:2)] <- scale(w.data2[, -c(1:2)])

# ecogeographic variation for size
fit_w_size1 <-
  brm(
    size ~ sex + age + temperature + precipitation + (1 |
                                                        population),
    data = w.data2,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.9)
  )

fit_w_size2 <-
  brm(
    size ~ sex + age + precipitation + (1 |
                                          population),
    data = w.data2,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.9)
  )

fit_w_size3 <-
  brm(
    size ~ sex + age + temperature + (1 |
                                        population),
    data = w.data2,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.9)
  )

fit_w_size4 <-
  brm(
    size ~ sex + age + (1 |
                          population),
    data = w.data2,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.9)
  )

fit_w_size1_waic <- add_ic(fit_w_size1, ic = c("waic", "r2"))
fit_w_size2_waic <- add_ic(fit_w_size2, ic = c("waic", "r2"))
fit_w_size3_waic <- add_ic(fit_w_size3, ic = c("waic", "r2"))
fit_w_size4_waic <- add_ic(fit_w_size4, ic = c("waic", "r2"))

sink('./results/fit_w_size_summary.txt')
WAIC(fit_w_size1_waic,
     fit_w_size2_waic,
     fit_w_size3_waic,
     fit_w_size4_waic)
cat('\\n')
cat('fit_w_size1 R2', mean(fit_w_size1_waic$R2), '\\n')
cat('fit_w_size2 R2', mean(fit_w_size2_waic$R2), '\\n')
cat('fit_w_size3 R2', mean(fit_w_size3_waic$R2), '\\n')
cat('fit_w_size4 R2', mean(fit_w_size4_waic$R2), '\\n')
sink()
# bayes_R2(fit_w_size1_waic)
# bayes_R2(fit_w_size2_waic)
# bayes_R2(fit_w_size3_waic)
# bayes_R2(fit_w_size4_waic)

# plot(fit_w_size1)
sink('./results/fit_w_size1.txt')
print (fit_w_size1) 
sink()  

# plot(fit_w_size2)
sink('./results/fit_w_size2.txt')
print (fit_w_size2) 
sink()  

# plot(fit_w_size3)
sink('./results/fit_w_size3.txt')
print (fit_w_size3) 
sink()  

# plot(fit_w_size4)
sink('./results/fit_w_size4.txt')
print (fit_w_size4) 
sink()  


# ecogeographic variation for shape
fit_w_shape1 <-
  brm(
    cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) ~ size + sex + age + temperature + precipitation + (1 |
                                                                                                                   population),
    data = w.data2,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.99)
  )

fit_w_shape2 <-
  brm(
    cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) ~ size + sex + age + precipitation + (1 |
                                                                                                     population),
    data = w.data2,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.99)
  )

fit_w_shape3 <-
  brm(
    cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) ~ size + sex + age + temperature + (1 |
                                                                                                   population),
    data = w.data2,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.99)
  )

fit_w_shape4 <-
  brm(
    cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) ~ size + sex + age + (1 |
                                                                                     population),
    data = w.data2,
    chains = 4,
    iter = 10000,
    cores = 4,
    control = list(adapt_delta = 0.99)
  )

fit_w_shape1_waic <- add_ic(fit_w_shape1, ic = c("waic", "r2"))
fit_w_shape2_waic <- add_ic(fit_w_shape2, ic = c("waic", "r2"))
fit_w_shape3_waic <- add_ic(fit_w_shape3, ic = c("waic", "r2"))
fit_w_shape4_waic <- add_ic(fit_w_shape4, ic = c("waic", "r2"))

sink('./results/fit_w_shape_summary.txt')
WAIC(fit_w_shape1_waic,
     fit_w_shape2_waic,
     fit_w_shape3_waic,
     fit_w_shape4_waic)
cat('\\n')
cat('fit_w_shape1 R2', mean(fit_w_shape1_waic$R2), '\\n')
cat('fit_w_shape2 R2', mean(fit_w_shape2_waic$R2), '\\n')
cat('fit_w_shape3 R2', mean(fit_w_shape3_waic$R2), '\\n')
cat('fit_w_shape4 R2', mean(fit_w_shape4_waic$R2), '\\n')
sink()
# bayes_R2(fit_w_shape1_waic)
# bayes_R2(fit_w_shape2_waic)
# bayes_R2(fit_w_shape3_waic)
# bayes_R2(fit_w_shape4_waic)

# plot(fit_w_shape1)
sink('./results/fit_w_shape1.txt')
print (fit_w_shape1) 
sink() 

# plot(fit_w_shape2)
sink('./results/fit_w_shape2.txt')
print (fit_w_shape2) 
sink() 

# plot(fit_w_shape3)
sink('./results/fit_w_shape3.txt')
print (fit_w_shape3) 
sink() 

# plot(fit_w_shape4)
sink('./results/fit_w_shape4.txt')
print (fit_w_shape4) 
sink() 



# mcmc plot
mcmc_w_size <- NULL
for (i in 1:4) {
  mcmc_w_size <- rbind(mcmc_w_size, as.mcmc(fit_w_size1)[[i]])
}

mcmc_w_shape <- NULL
for (i in 1:4) {
  mcmc_w_shape <- rbind(mcmc_w_shape, as.mcmc(fit_w_shape1)[[i]])
}

mcmc_w_size_sum <-
  rbind2(
    cbind(data.frame(variable = "size"), rep(NA, 20000)),
    cbind(data.frame(variable = "sex"), mcmc_w_size[, 2]),
    cbind(data.frame(variable = "age"), mcmc_w_size[, 3]),
    cbind(data.frame(variable = "temperature"), mcmc_w_size[, 4]),
    cbind(data.frame(variable = "precipitation"), mcmc_w_size[, 5]),
    cnames = c("explanatory", "size")
  )

mcmc_w_shape_sum <-
  rbind2(
    cbind(data.frame(variable = "size"), mcmc_w_shape[, c(11, 16, 21, 26, 31, 36, 41, 46, 51, 56)]),
    cbind(data.frame(variable = "sex"), mcmc_w_shape[, c(11, 16, 21, 26, 31, 36, 41, 46, 51, 56) + 1]),
    cbind(data.frame(variable = "age"), mcmc_w_shape[, c(11, 16, 21, 26, 31, 36, 41, 46, 51, 56) + 2]),
    cbind(data.frame(variable = "temperature"), mcmc_w_shape[, c(11, 16, 21, 26, 31, 36, 41, 46, 51, 56) + 3]),
    cbind(data.frame(variable = "precipitation"), mcmc_w_shape[, c(11, 16, 21, 26, 31, 36, 41, 46, 51, 56) + 4]),
    cnames = c("explanatory", paste("PC", 1:10, sep = ""))
  )

mcmc_w_sum <- cbind(mcmc_w_shape_sum, mcmc_w_size_sum)

p5 <-
  ggplot(
    data = summaryCI(
      melt(mcmc_w_sum),
      measurevar = "value",
      groupvars = c("explanatory", "variable"),
      na.rm = TRUE
    ),
    mapping = aes(x = explanatory, y = mean)
  )
p5 <-
  p5 + geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0) + geom_point(position = position_dodge(0.2))
p5 <- p5 + geom_hline(yintercept = 0)
p5 <- p5 + facet_wrap(~ variable, ncol = 5)
p5 <-
  p5 + ylim(-1, 1) + coord_flip() + ylab("coefficient") + xlab("variables")
p5 <- p5 + scale_y_continuous(breaks = seq(-0.8, 0.8, 0.8))
quartz(width = 7.25, height = 4)
plot (p5)
dev.copy2eps(file = "./results/stan_wild.eps")
dev.off()





### vector angles --------------------------------------------------------
vec.ang.cap_temp2 <-
  ang.vec2(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "captivity", -1],
           mcmc_w_shape_sum[mcmc_w_shape_sum$explanatory == "temperature", -1])

vec.ang.cap_prec2 <-
  ang.vec2(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "captivity", -1],
           mcmc_w_shape_sum[mcmc_w_shape_sum$explanatory == "precipitation", -1])

vec.ang.cap_size <-
  ang.vec(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "captivity", -1],
          mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "size", -1])

vec.ang.cap_sex <-
  ang.vec(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "captivity", -1],
          mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "sex", -1])

vec.ang.cap_age <-
  ang.vec(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "captivity", -1],
          mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "age", -1])

vec.ang.sex_age2 <-
  ang.vec2(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "sex", -1],
           mcmc_w_shape_sum[mcmc_w_shape_sum$explanatory == "age", -1])

vec.ang.sex_age <-
  ang.vec(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "sex", -1],
          mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "age", -1])

vec.ang.sex_size2 <-
  ang.vec2(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "sex", -1],
           mcmc_w_shape_sum[mcmc_w_shape_sum$explanatory == "size", -1])

vec.ang.sex_size <-
  ang.vec(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "sex", -1],
          mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "size", -1])

vec.ang.size_age2 <-
  ang.vec2(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "size", -1],
           mcmc_w_shape_sum[mcmc_w_shape_sum$explanatory == "age", -1])

vec.ang.size_age <-
  ang.vec(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "size", -1],
          mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "age", -1])

vec.angle_df <-
  rbind(
    data.frame(
      vec.angle = vec.ang.sex_age[[2]],
      vector_comparison = rep("sex-age", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.angle = vec.ang.sex_age2[[2]],
      vector_comparison = rep("sex-age", 40000),
      type = rep("between", 40000)
    ),
    data.frame(
      vec.angle = vec.ang.sex_size[[2]],
      vector_comparison = rep("sex-size", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.angle = vec.ang.sex_size2[[2]],
      vector_comparison = rep("sex-size", 40000),
      type = rep("between", 40000)
    ),
    data.frame(
      vec.angle = vec.ang.size_age[[2]],
      vector_comparison = rep("size-age", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.angle = vec.ang.size_age2[[2]],
      vector_comparison = rep("size-age", 40000),
      type = rep("between", 40000)
    ),
    data.frame(
      vec.angle = vec.ang.cap_size[[2]],
      vector_comparison = rep("cap-size", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.angle = vec.ang.cap_sex[[2]],
      vector_comparison = rep("cap-sex", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.angle = vec.ang.cap_age[[2]],
      vector_comparison = rep("cap-age", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.angle = vec.ang.cap_temp2[[2]],
      vector_comparison = rep("cap-temp", 40000),
      type = rep("between", 40000)
    ),
    data.frame(
      vec.angle = vec.ang.cap_prec2[[2]],
      vector_comparison = rep("cap-prec", 40000),
      type = rep("between", 40000)
    )
  )

vec.angle_df2 <-
  summaryCI(
    vec.angle_df,
    measurevar = "vec.angle",
    groupvars = c("vector_comparison", "type")
  )

p6 <-
  ggplot(vec.angle_df2, aes(x = vector_comparison, y = mean, colour = type))
p6 <- p6 + geom_hline(yintercept = 90)
p6 <-
  p6 + geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                     width = 0,
                     position = position_dodge(0.2)) +
  geom_point(position = position_dodge(0.2))
p6 <-
  p6 + scale_y_continuous(breaks = seq(30, 150, 30), limits = c(40, 140))
p6 <- p6 + labs(x = "vector comparison", y = "vector angle")
p6 <- p6 + theme(legend.position = 'none')

quartz(width = 7.25, height = 3.25)
p6
dev.copy2pdf(file = "./results/vec_angle.pdf")
dev.off()

write.csv(vec.angle_df2, './results/vec_angle.csv')

vec.cor_df <-
  rbind(
    data.frame(
      vec.cor = vec.ang.sex_age[[1]],
      vector_comparison = rep("sex-age", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.cor = vec.ang.sex_age2[[1]],
      vector_comparison = rep("sex-age", 40000),
      type = rep("between", 40000)
    ),
    data.frame(
      vec.cor = vec.ang.sex_size[[1]],
      vector_comparison = rep("sex-size", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.cor = vec.ang.sex_size2[[1]],
      vector_comparison = rep("sex-size", 40000),
      type = rep("between", 40000)
    ),
    data.frame(
      vec.cor = vec.ang.size_age[[1]],
      vector_comparison = rep("size-age", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.cor = vec.ang.size_age2[[1]],
      vector_comparison = rep("size-age", 40000),
      type = rep("between", 40000)
    ),
    data.frame(
      vec.cor = vec.ang.cap_size[[1]],
      vector_comparison = rep("cap-size", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.cor = vec.ang.cap_sex[[1]],
      vector_comparison = rep("cap-sex", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.cor = vec.ang.cap_age[[1]],
      vector_comparison = rep("cap-age", 20000),
      type = rep("within", 20000)
    ),
    data.frame(
      vec.cor = vec.ang.cap_temp2[[1]],
      vector_comparison = rep("cap-temp", 40000),
      type = rep("between", 40000)
    ),
    data.frame(
      vec.cor = vec.ang.cap_prec2[[1]],
      vector_comparison = rep("cap-prec", 40000),
      type = rep("between", 40000)
    )
  )


vec.cor_df2 <-
  summaryCI(
    vec.cor_df,
    measurevar = "vec.cor",
    groupvars = c("vector_comparison", "type")
  )

p10 <-
  ggplot(vec.cor_df2, aes(x = vector_comparison, y = mean, colour = type))
p10 <- p10 + geom_hline(yintercept = 0)
p10 <-
  p10 + geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                     width = 0,
                     position = position_dodge(0.2)) +
  geom_point(position = position_dodge(0.2))
p10 <-
  p10 + scale_y_continuous(breaks = seq(-0.75, 0.75, 0.25), limits = c(-0.8, 0.8))
p10 <- p10 + labs(x = "vector comparison", y = "vector correlation")
p10 <- p10 + theme(legend.position = 'none')

quartz(width = 7.25, height = 3.25)
p10
dev.copy2pdf(file = "./results/vec_cor.pdf")
dev.off()

# plot
vec.comp_df <- rbind(vec.angle_df2, vec.cor_df2) %>%
  mutate(method = c(rep("angle", 11), rep("correlation", 11))) %>%
  mutate(base = c(rep(90, 11), rep(0, 11)))

p11 <- ggplot(vec.comp_df, aes(x = vector_comparison, y = mean, colour = type))
p11 <- p11 + geom_hline(aes(yintercept = base))
p11 <-
  p11 + geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                      width = 0,
                      position = position_dodge(0.2)) +
  geom_point(position = position_dodge(0.2))
p11 <- p11 + theme(legend.position = 'none')
p11 <- p11 + facet_wrap(~ method, ncol = 1, scales = "free_y")
p11 <- p11 + labs(x = "vector comparison")
p11 <- p11 + theme(axis.title.y = element_blank())

quartz(width = 7.25, height = 6)
p11
dev.copy2pdf(file = "./results/vec_comp.pdf")
dev.off()

### shape score ----------------------------------------------------------------------
# shape score
shape_cap <- shape_vec(A = gpa$Sym ,
                       B = apply(mcmc_a_shape_sum[mcmc_a_shape_sum$explanatory == "captivity", -1], 2, mean))

# correlation between shape score and PC4
cor.test(shape_cap$shape.score, a.data0$PC4)

# boxplot for PC4
bxPC4_df <- a.data0[, c("PC4", "sex", "captivity")]
bxPC4_df[bxPC4_df$captivity == 0,]$captivity <- "0"
bxPC4_df[bxPC4_df$captivity == 1,]$captivity <- "1"
bxPC4_df[bxPC4_df$captivity == 2,]$captivity <- "2"
p7 <-
  ggplot(data = bxPC4_df, aes(y = PC4, x = captivity, fill = sex))
p7 <- p7 + geom_boxplot()
p7 <- p7 + labs(title = "(a)")

# marginal effect for PC4
me_fit_a_shape1 <-
  marginal_effects(
    fit_a_shape1,
    effects = "captivity:sex",
    int_conditions = list(sex = setNames(c(
      -0.8692996, 1.1438153
    ), c("F", "M"))),
    plot = F,
    resp = "PC4"
  )

me_df <-
  me_fit_a_shape1$`PC4.PC4_captivity:sex`[, c("sex", "captivity", "estimate__", "lower__", "upper__")]

me_df$captivity <-
  me_df$captivity * sd(a.data0$captivity) + mean(a.data0$captivity)
me_df[, c("estimate__", "lower__", "upper__")] <-
  me_df[, c("estimate__", "lower__", "upper__")] * sd(a.data0$PC4) + mean(a.data0$PC4)

p8 <- ggplot(data = me_df, aes(x = captivity))
p8 <-
  p8 + geom_ribbon(
    aes(ymin = lower__, ymax = upper__),
    data = me_df[me_df$sex == "M", ],
    alpha = 0.4,
    fill = "#00BFC4"
  )
p8 <-
  p8 + geom_ribbon(
    aes(ymin = lower__, ymax = upper__),
    data = me_df[me_df$sex == "F", ],
    alpha = 0.4,
    fill = "#F8766D"
  )
p8 <-
  p8 + geom_line(aes(y = estimate__),
                 data = me_df[me_df$sex == "M", ],
                 colour = "#00BFC4",
                 size = 1)
p8 <-
  p8 + geom_line(aes(y = estimate__),
                 data = me_df[me_df$sex == "F", ],
                 colour = "#F8766D",
                 size = 1)
p8 <-
  p8 + geom_point(aes(x = captivity, y = PC4, color = sex), data = a.data0)
p8 <- p8 + ylab("PC4")
p8 <-
  p8 + scale_x_continuous(breaks = seq(0, 2, by = 1), limits = c(0, 2))
p8 <- p8 + labs(title = "(b)")


# boxplot for shape score
temp_df <-
  cbind(data.frame(shape.score = shape_cap$shape.score), a.data0[, c("sex", "captivity")])
temp_df[temp_df$captivity == 0,]$captivity <- "0"
temp_df[temp_df$captivity == 1,]$captivity <- "1"
temp_df[temp_df$captivity == 2,]$captivity <- "2"
p9 <-
  ggplot(data = temp_df, aes(y = shape.score, x = captivity, fill = sex))
p9 <- p9 + geom_boxplot()
p9 <- p9 + ylab("shape score")


# plot
p7 <- p7 + theme(legend.position = 'none')
p8 <- p8 + theme(legend.position = 'none')
p9 <- p9 + theme(legend.position = 'none') + labs(title = "(c)")

quartz(width = 7.25, height = 3.54)
grid.arrange(p7, p8, p9, ncol = 3)
dev.copy2pdf(file = "./results/cap_effect_on_shape.pdf")
dev.off()



### visualization -------------------------------------------------------
# create mesh
ref_mesh <- file2mesh("./reference_shape/ref_mesh.ply")
ref_land <- readLandmarks.csv(file = './reference_shape/ref_land.csv', x = 1:28, y = 1:3, rownames = NULL, header = FALSE,
                              dec = ".", sep = ",")$LM

# shape change along shape score
mesh2ply(tps3d(ref_mesh, ref_land, shape_cap$shapes[, , 1]),
         "./results/shape_cap_n")
mesh2ply(tps3d(ref_mesh, ref_land, shape_cap$shapes[, , 2]),
         "./results/shape_cap_p")

# shape changes along PC axes
for (i in 1:10) {
  pci_n <-
    tps3d(ref_mesh, ref_land, pca$pc.shapes[[2 * i - 1]])
  pci_p <-
    tps3d(ref_mesh, ref_land, pca$pc.shapes[[2 * i]])
  mesh2ply(pci_p, paste("./results/PC", i, "_p", sep = ""))
  mesh2ply(pci_n, paste("./results/PC", i, "_n", sep = ""))
}

# mean shape
mean_mesh <- tps3d(ref_mesh, ref_land, gpa$mshape)
mesh2ply(mean_mesh, "./results/mean")






### basic presentations of the results ----------------------------------

temp_df <- a.data0 %>%
  gather(key = "variable", value = "value", paste("PC", 1:10, sep=""), "size") %>%
  mutate(population = fct_relevel(population, "Shimokita", "Kinkazan", "Hakusan", "Nagatoro", "Kamimatsu",
                                  "Matsukawa", "Takamori", "Takahama", "WakasaF", "Inuyama", "WakasaS",
                                  "Mitoya", "Kotsugu", "Yoshida", "Arashiyama", "Koga", "Hasumi", "Izu",
                                  "Minoo", "Shodoshima", "Miyajima", "Noguchi", "Yakushima")) %>%
  mutate(variable = fct_relevel(variable, paste("PC", 1:10, sep=""), "size"))
  

p10 <- ggplot(temp_df, aes(x = population, y = value, color = environment)) + 
  geom_quasirandom(dodge.width = 1, cex = 1) +
  facet_wrap(variable ~ ., ncol = 2, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'none')

quartz(width = 7.25, height = 8.75)
plot (p10)
dev.copy2eps(file = "./results/basic_prop.eps")
dev.off()


### allometry ---------------------------------------------------------------
procD.lm(shape ~ size, data = geomorph.data.frame(shape = clean.lm, size = a.data0$size))



# correlation test -----------------------------
df_for_cor <- cbind(a.data0[,c("ID", "sex", "age", "age_y", "generation", "size", "PC4")], data.frame(shape.score = shape_cap$shape.score))

summary(lm(PC4 ~ size + sex + age + age_y, data = df_for_cor))
summary(lm(shape.score ~ size + sex + age + age_y, data = df_for_cor))

summary(lm(PC4 ~ size + sex + age + age_y + generation, data = df_for_cor[df_for_cor$generation > 0,]))
summary(lm(shape.score ~ size + sex + age + age_y + generation, data = df_for_cor[df_for_cor$generation > 0,]))

cor.test(
  df_for_cor[df_for_cor$age == "full",]$age_y,
  df_for_cor[df_for_cor$age == "full",]$PC4
)

cor.test(
  df_for_cor[df_for_cor$age == "full",]$age_y,
  df_for_cor[df_for_cor$age == "full",]$shape.score
)



