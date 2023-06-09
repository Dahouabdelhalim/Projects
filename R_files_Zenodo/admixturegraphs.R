library(admixturegraph)

laredo_fs <- read.csv("laredo_f.csv", header=TRUE)
plot(f4stats(laredo_fs))

#laredoensis single_origin outgroup
sin_leaves <- c("gularis", "laredoA", "laredoB", "sexlin", "tigris")
sin_inner_nodes <- c("a", "c", "b", "d", "e", "root")
sin_edges <- parent_edges(c(edge("gularis", "c"), edge("sexlin", "d"), edge("laredoA", "a"), edge("laredoB", "a"), edge("tigris", "root"), edge("a", "b"), edge("b", "c"), edge("b", "d"), edge("c", "e"), edge("d", "e"), admixture_edge("b", "c", "d", "z"), edge("e", "root")))
sin_graph <- agraph(sin_leaves, sin_inner_nodes, sin_edges)

#laredoensis multi_origin old A
mult1_leaves <- c("gularis", "laredoA", "laredoB", "sexlin", "tigris")
mult1_inner_nodes <- c("a", "b", "c", "d", "e", "f", "g", "root")
mult1_edges <- parent_edges(c(edge("gularis", "c"), edge("sexlin", "b"), edge("laredoA", "d"), edge("laredoB", "a"), edge("tigris", "root"), admixture_edge("a", "c", "b", "z"), edge("c", "e"), edge("b", "f"), admixture_edge("d", "e", "f", "y"), edge("e", "g"), edge("f", "g"), edge("g", "root")))
mult1_graph <- agraph(mult1_leaves, mult1_inner_nodes, mult1_edges)

#laredoensis multi_origin old B
mult2_leaves <- c("gularis", "laredoA", "laredoB", "sexlin", "tigris")
mult2_inner_nodes <- c("a", "b", "c", "d", "e", "f", "g", "root")
mult2_edges <- parent_edges(c(edge("gularis", "c"), edge("sexlin", "b"), edge("laredoB", "d"), edge("laredoA", "a"), edge("tigris", "root"), admixture_edge("a", "c", "b", "z"), edge("c", "e"), edge("b", "f"), admixture_edge("d", "e", "f", "y"), edge("e", "g"), edge("f", "g"), edge("g", "root")))
mult2_graph <- agraph(mult2_leaves, mult2_inner_nodes, mult2_edges)


laredo_sin_fit <- fit_graph(laredo_fs, sin_graph)
summary(laredo_sin_fit)

mcmc_sin <- make_mcmc_model(sin_graph, laredo_fs)
initial_sin <- rep(0.5, length(mcmc_sin$parameter_names))
sin_chain1 <- run_metropolis_hasting(mcmc_sin, initial_sin, iterations = 1000000, no_temperatures=4, verbose = FALSE,cores=4)
write.table(sin_chain1,file="lared_single_mcmc.txt", quote=FALSE, sep= "\\t")
