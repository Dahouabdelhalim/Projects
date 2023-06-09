devtools::install_github("https://github.com/revbayes/RevGadgets/tree/development")

Kt_mean <- readOBDP( start_time_trace_file="start_time_trace_N85_S200.txt", 
                     popSize_distribution_matrices_file="Kt_trace_N85_S200.txt", 
                     trees_trace_file="mcmc_OBDP_Cetaceans_wellMixedTrees.trees" )

p <- plotDiversityOBDP( Kt_mean,
                        xlab="Time (My)",
                        ylab="Number of lineages",
                        xticks_n_breaks=21,
                        col_Hidden="dodgerblue3",
                        col_LTT="gray25",
                        col_Total="forestgreen",
                        col_Hidden_interval="dodgerblue2",
                        col_Total_interval="darkolivegreen4",
                        palette_Hidden=c("transparent", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black"),
                        palette_Total=c("transparent", "green4", "forestgreen", "black"),
                        line_size=0.7,
                        interval_line_size=0.5,
                        show_Hidden=FALSE,
                        show_LTT=TRUE,
                        show_Total=TRUE,
                        show_intervals=TRUE,
                        show_densities=TRUE,
                        show_expectations=TRUE,
                        use_interpolate=TRUE )
p

library(deeptime)

q <- gggeo_scale(p, dat="periods", height=ggplot2::unit(1.3, "line"), abbrv=F, size=4.5, neg=T)
r <- gggeo_scale(q, dat="epochs", height=ggplot2::unit(1.1, "line"), abbrv=F, size=3.5, neg=T, skip=c("Paleocene", "Pliocene", "Pleistocene", "Holocene"))
s <- gggeo_scale(r, dat="stages", height=ggplot2::unit(1, "line"), abbrv=T, size=2.5, neg=T)
s
ggsave(filename="nbLineages_Cetacea_genera_Total.svg", plot=s, device="svg", width=12, height=6)
ggsave(filename="nbLineages_Cetacea_genera_Total.pdf", plot=s, device="pdf", width=12, height=6)
ggsave(filename="nbLineages_Cetacea_genera_Total.png", plot=s, device="png", width=12, height=6)

library("ggtree")

library("treeio")
tree <- read.beast(beast_file)
tree

tree <- read.mrbayes("mcmc_OBDP_Cetaceans.tre")
mit_rates <- read.table("mcmc_OBDP_Cetaceans_mit_rates.out", header=T, row.names=1)

ggtree(tree, aes(color=c(colMeans(mit_rates), rep(.015,16)))) +
  scale_color_continuous(low='darkgreen', high='red', name="Mit. Rates") +
  theme(legend.position="right")
