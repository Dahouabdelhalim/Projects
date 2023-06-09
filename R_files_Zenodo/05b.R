##02c.phylosignal
labels.traits.continuous <- c("Seed mass", "Leaf dry matter content", "SLA",
                              "Leaf N content", "Leaf C content", "Leaf P content",
                              "Stem dry matter content", "Leaf circularity",
                              "Vegetative height", "Leaf length", "Leaf thickness", "Leaf width",
                              "Petiole length","Genome size")
###
labels.signal <- c(labels.productivity, labels.traits.continuous)

tr.prairie.phylosig <- tr.prairie.biomassPlot
tr.prairie.phylosig$node.label <- NULL
message('Doing Blomberg\\'s K')
prairie.phylosignal <- lapply(labels.signal, function(x) {
  phylosignal(structure(all.prairie.ordi[tr.prairie.phylosig$tip.label, x], names = tr.prairie.phylosig$tip.label),
              tr.prairie.phylosig)[1,]
}) %>%
  do.call('rbind', .)
# row.names(prairie.phylosignal) <- names(all.prairie.small)

message('Doing Pagel\\'s Lambda')
prairie.lambda <- list(estimated = fitContinuous(tr.prairie.phylosig,
                                                 all.prairie.ordi[, labels.signal],
                                                 model = 'lambda'),
                       zero = fitContinuous(geiger::rescale(tr.prairie.phylosig, 'lambda', 0),
                                            all.prairie.ordi[, labels.signal])
)
prairie.phylosignal <- cbind(lambda = sapply(prairie.lambda$estimated, function(x) x$opt$lambda),
                             L.ratio = 2 * (sapply(prairie.lambda$estimated, function(x) x$opt$lnL) -
                                              sapply(prairie.lambda$zero, function(x) x$opt$lnL)),
                             L.ratio.p = pchisq(2 * (sapply(prairie.lambda$estimated, function(x) x$opt$lnL) -
                                                       sapply(prairie.lambda$zero, function(x) x$opt$lnL)),
                                                df = 1,
                                                lower.tail = F),
                             prairie.phylosignal)
prairie.phylosignal.pretty <- prairie.phylosignal[order(prairie.phylosignal$lambda, decreasing = T), ]
prairie.phylosignal.pretty<- round(prairie.phylosignal.pretty, 5)
prairie.phylosignal.pretty[prairie.phylosignal.pretty == 0] <- "< 0.00001"
prairie.phylosignal.pretty[prairie.phylosignal.pretty == 1] <- "> 0.99999"
write.csv(prairie.phylosignal.pretty, file = '../OUT/TABLE.phylosignal.csv')

##########
#04.

#dat.traits.cor <- data.frame(
#  NDVI = cor(all.prairie.ordi$NDVI, all.prairie.ordi[, labels.traits.continuous])[1, ],
#  Biomass = cor(all.prairie.ordi$Biomass, all.prairie.ordi[, labels.traits.continuous])[1, ],
#  lambda = prairie.phylosignal[labels.traits.continuous, 'lambda']
#)

prairie.phylosignal <- read.csv(file = "../DATA/TABLE.phylosignal.csv")

dat.traits.cor <- data.frame(
    simper = cor(all.prairie.ordi$NDVI, all.prairie.ordi[, labels.traits.continuous])[1, ],
    coverChange = cor(all.prairie.ordi$Biomass, all.prairie.ordi[, labels.traits.continuous])[1, ],
    lambda = prairie.phylosignal[1:13, 1:2] #labels.traits.continuous as first 13 rows and second column=lambda
  )

