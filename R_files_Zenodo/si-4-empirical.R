# Empirical (mammal) analyses

#Headers
source("si-2-functions.R")
require(vegan)
require(pez)
require(ape)
require(picante)
require(OUwie)
require(phytools)
require(plotrix)
require(suppdata)

#Basic phylogenetic/clade method
make.clade.comm <- function(comm, tree, mean=TRUE){
    clade.mat <- clade.matrix(tree)$clade.matrix
    output <- counts <- matrix(0, nrow=nrow(comm), ncol=nrow(clade.mat))
    #Loop over species
    for(i in seq(ncol(comm))){
        output[ ,clade.mat[,i]==1] <- output[ ,clade.mat[,i]==1] + comm[,i]
        counts[ ,clade.mat[,i]==1] <- counts[ ,clade.mat[,i]==1] + 1
    }
    if(mean)
        output[output != 0] <- output[output != 0] / counts[output != 0]
    return(output)
}

#Load community data and make into site-species matrix
data <- read.csv(suppdata("E092-201", "MCDB_communities.csv", "esa_data_archives"))
data$Abundance <- as.numeric(as.character(data$Abundance))
data <- data[!is.na(data$Abundance) & data$Abundance > 0, ]
species <- read.csv(suppdata("E092-201", "MCDB_species.csv", "esa_data_archives"))
data$species <- paste(species$Genus, species$Species)[match(data$Species_ID, species$Species_ID)]
data$species <- gsub("  ", " ", data$species)
comm <- as.matrix(sample2matrix(data.frame(data$Site_ID, as.numeric(as.character(data$Abundance)), data$species)))

#Load phylogeny
tree <- read.nexus(suppdata("10.1111/j.1461-0248.2009.01307.x", 1))[[1]]
tree$tip.label <- gsub("_", " ", tree$tip.label)

#Load traits data (PanTHERIA)
t <- read.delim(suppdata("E090-184", "PanTHERIA_1-0_WR05_Aug2008.txt", "esa_archives"))
rownames(t) <- t$MSW05_Binomial
traits <- data.frame(body.mass = log10(t$X5.1_AdultBodyMass_g), row.names=rownames(t))
traits <- na.omit(traits)

#Merge datasets together
c.data <- comparative.comm(tree, comm, traits)

#Calculate variance et al.
var <- clade.var(c.data, 0.05)
null <- matrix(NA, ncol=9999, nrow=nrow(var))
t.data <- c.data
for(i in seq_len(ncol(null))){
    prog.bar(i, ncol(null))
    t.data$comm <- t.data$comm[,sample(ncol(t.data$comm))]
    null[,i] <- clade.var(t.data, 0.05)$variance
}

# This takes a while so you may want to save out the workspace
#save.image("si-4-empirical.RData")

# Grab data for a plot
var$rank <- numeric(nrow(var))
for(i in seq_len(nrow(var)))
    var$rank[i] <- rank(c(var$variance[i],null[i,]))[1]

# Phylo side-plot
node.tree.plot <- function(tree, ranks, foc.nodes, cols=c("red","blue"), base.node.cex=1.5, foc.node.cex=5, tip.cex=.5, n.colors=10000, ...){
    ranks <- round(ranks)
    cols <- colorRampPalette(cols)(n.colors)
    sizes <- ifelse(seq_along(ranks) %in% foc.nodes, foc.node.cex, base.node.cex)
    sizes <- ifelse(seq_along(sizes) <= length(tree$tip.label), tip.cex, sizes)
    
    t <- tree$edge[,1][tree$edge[,1] > length(tree$tip.label)]
    matching <- match(t, seq_along(ranks))
    
    to.return <- plot(tree, show.tip.label=FALSE, edge.color=cols[ranks[matching]], ...)
    nodelabels(pch=20, node=seq_along(ranks), cex=sizes, col=cols[ranks])
    invisible(to.return)
}

pdf("combined-plot.pdf")
# Setup
layout(matrix(1:3, nrow=1, ncol=3), width=c(.2,.6,.2))
# First phylo plot (using wrapper)
node.tree.plot(c.data$phy, var$rank, c(609,618), edge.width=1.5, tip.cex=0, base.node.cex=0, no.margin=TRUE, n.colors=10000)
gradient.rect(10, 50, 20, 300, col=colorRampPalette(c("red","blue"))(1001), gradient="up", border=NA)
text(rep(21,4), seq(50,300,length.out=6), c(seq(0,1,length.out=6)), adj=0)
text(rep(10,2), c(40,310), c(expression(beta~overdispersed),expression(beta~clustered)), adj=0, font=2)
# Central matrix plot
details <- plot(c.data$phy, edge.color=NA, show.tip.label=FALSE)
image(log10(c.data$comm[order(rowSums(c.data$comm>0),decreasing=TRUE),]), y=1:483, x=seq(0,147.1,length.out=939),add=TRUE, col=colorRampPalette(c("yellow","red"))(1000))
gradient.rect(60, -15, 140, -10, col=colorRampPalette(c("yellow","red"))(9999), border=NA)
text(seq(60,140,length.out=5), rep(-5, 2), 10^round(c(seq(-1,3.84,length.out=5))), adj=0.5)
text(10, -12.5, "Relative abundance", font=2, adj=0)
# Second phylo plot
details <- plotBranchbyTrait(c.data$phy, c.data$data$body.mass, palette=colorRampPalette(c("grey10","grey90")), show.tip.label=FALSE, direction="left", edge.width=1.5, legend=FALSE)
gradient.rect(130, 50, 140, 300, col=colorRampPalette(c("grey10","grey90"))(9999), gradient="up", border=NA)
text(rep(129,4), seq(50,300,length.out=5), round(c(10^seq(0.4,3.9,length.out=5))), adj=1)
text(rep(130,2), c(310,40), c("body mass","(g)"), adj=1, font=2)
nodelabels(node=c(609,618), pch=20, col="red", cex=5)
dev.off()
