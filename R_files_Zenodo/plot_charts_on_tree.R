# # # # plotting pie charts on phylogenetic tree (nodes), e.g. BioGeoBEARS output (exported per area probabilities)

library(ggtree)
library(ggimage)

# set working directory
setwd("/path/to/BioGeoBEARS/output/")

# load tree and probability of occupancy per area for each tip/node data
mytree <- read.tree("BioGeoBEARS_input_tree_pruned.newick")
nodes <- read.csv("area_prob_nodes.csv", header = TRUE)
tips <- read.csv("area_prob_tips.csv", header = TRUE)


# create base tree 
p <- ggtree(mytree, size=0.05, ladderize = TRUE)

# create pie chart infromation for nodes and bar chart information for tips
pies <- nodepie(nodes, cols=2:9, color=c(Pal="#1B9E77", Nea="#FF4500", Sea="#0000CD", Neo="#666666", Afr="#FFD700", Aus="#66A61E", Ant="#A6761D", Ind="#CD2990"))
bars <- nodebar(tips, cols=2:9, position='dodge', color=c(Pal="#1B9E77", Nea="#FF4500", Sea="#0000CD", Neo="#666666", Afr="#FFD700", Aus="#66A61E", Ant="#A6761D", Ind="#CD2990"))

# col2rgb(c(PAL="#1B9E77", NEA="#FF4500", SEA="#0000CD", NEO="#666666", AFR="#FFD700", AUS="#66A61E", ANT="#A6761D", IND="#CD2990"))



pdffn = "BioGeoBEARS_results.pdf"
pdf(pdffn, width=80, height=80)

p + 
  geom_inset(pies, width = .015, height = .015) +
  geom_inset(bars, width = .05, height = .0045, hjust=-6) + 
  geom_tiplab(offset = 5)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)



