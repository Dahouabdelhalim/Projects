## # # # # # # # # # # # # # # # # # # # # # # # #
##
## Project: CDKL5-NLS Phosphoproteomics
##
## Script name: Script_File_S4_Interactor_analysis_Phospho-TMT.R
##
## Purpose of script: Visualize interaction network
##
## Author: Florian Weiland
##
## Date Created: 2020-05-13
##
## # # # # # # # # # # # # # # # # # # # # # # # #

## Libraries ----

suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggnetwork, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(wesanderson, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggrepel, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(STRINGdb, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(igraph, warn.conflicts = FALSE))

#Set palette ----

pal <- wes_palette("Darjeeling1")

## Read in database ----

message("Loading local STRING database")

inter.df <- read.csv(
  "./Databases/9606.protein.links.v11.0.txt",
  sep = " "
)

message("Establishing connection to online STRING database")

string_db <- STRINGdb$new(
  version = "11",
  species = 9606,
  score_threshold = 150,
  input_directory = "./Databases/"
)

## Read in data and map ENSEMBL names ----
message("Loading TMT data")

data <- read.csv(
  "./Output_R/Phospho-TMT_results_clean.csv"
)

## Only proteins with down-regulated phospho-peptides

data <- data[which(data$logFC > 0), ]

message("Mapping gene names to ENSP")
data <- string_db$map(data, "Gene.names")

## Get interaction partners
uni.gene <- unique(data$STRING_id)

int.poi <- inter.df[inter.df$protein1 %in% uni.gene, ] 

## Only direct interactions
int.poi.direct <- int.poi[int.poi$protein2 %in% uni.gene, ]

## Network analysis ----
message("Building Network")

set.seed(8)

params <- list()
params$seed.coord <- matrix(rnorm(nrow(int.poi.direct) * 2), ncol = 2)

poi.net <- ggnetwork(int.poi.direct[, 1:2], layout = "kamadakawai", arrow.gap = 0, layout.par = params)

## Ammend with gene names

poi.net$Gene.names <- data[match(poi.net$vertex.names, data$STRING_id), "Gene.names"]

## Remove self-interactions

rem <- which(poi.net$x == poi.net$xend & poi.net$y == poi.net$yend)

poi.net <- poi.net[-rem, ]

## Get combined_score and confidence level for each interaction
message("Ammending data")

for (cur.row in 1:nrow(poi.net)) {
  
  int.row <- which(poi.net[cur.row, "xend"] == poi.net$x & poi.net[cur.row, "yend"] == poi.net$y)
  int.row <- int.row[1]
  
  poi.net[cur.row, "Interactor"] <- poi.net[int.row, "vertex.names"]
  
  df.row.1 <- grep(poi.net[cur.row, "vertex.names"], int.poi$protein1)
  df.row.2 <- grep(poi.net[cur.row, "Interactor"], int.poi$protein2)
  
  df.row <- df.row.1[which(df.row.1 %in% df.row.2 == TRUE)]
  
  poi.net[cur.row, "Score"] <- int.poi[df.row, "combined_score"]
  
  if (poi.net[cur.row, "Score"] < 400) {
    
    poi.net[cur.row, "Confidence"] <- "Low"
    
  }
  
  if (700 > poi.net[cur.row, "Score"] & poi.net[cur.row, "Score"] >= 400) {
    
    poi.net[cur.row, "Confidence"] <- "Medium"
    
  } 
  
  if (900 > poi.net[cur.row, "Score"] & poi.net[cur.row, "Score"] >= 700) {
    
    poi.net[cur.row, "Confidence"] <- "High"
    
  }
  
  if (poi.net[cur.row, "Score"] >= 900) {
    
    poi.net[cur.row, "Confidence"] <- "Very High"
    
  }
  
}

## Calculate p-value ----
message("Calculating Statistics")

poi <- data$STRING_id

stats <- string_db$get_ppi_enrichment(poi)

## Get subnetwork and plot ----
message("Subnetwork Analysis")

sub <- string_db$get_subnetwork(poi)

sub.s <- simplify(sub, remove.multiple = TRUE)

ensp.labels <- V(sub.s)$name

labels <- data[match(ensp.labels, data$STRING_id), "Gene.names"]

message("Subnetwork Analysis - Writing graph")

svg(
  "Output_R/Subnetwork.svg",
  bg = "transparent"
)

plot(
  sub.s,
  vertex.label = labels,
  vertex.label.color = "black",
  vertex.label.cex = 1.1,
  vertex.label.font = 2
)

dev.off()

## Plot data ----

message("Network Analysis - Writing graph")

poi.net$Confidence <- factor(poi.net$Confidence, levels = c("Low", "Medium", "High", "Very High"))

## Adjust Gene names to harmonize with manuscript

poi.net$Gene.names <- sub("TCEB3", "ELOA", poi.net$Gene.names)
poi.net$Gene.names <- sub("MPLKIP", "TTDN1", poi.net$Gene.names)
poi.net$Gene.names <- sub("YLPM1", "ZAP3", poi.net$Gene.names)
poi.net$Gene.names <- sub("MAPRE2", "EB2", poi.net$Gene.names)

label.data <- poi.net[!duplicated(poi.net$Gene.names), ]

## Plot direct interactions ----

ggplot(poi.net, aes(x, y, xend = xend, yend = yend)) +
  scale_x_reverse() + 
  geom_edges(
    aes(colour = Confidence),
    size = 2
  ) +
  scale_color_manual(
    name = "Confidence",
    breaks = c("Low", "Medium", "Very High"),
    values = c("Low" = pal[5], "Medium" = pal[4], "Very High" = "black")
  ) +
  geom_nodes(
    size = 8,
    color = pal[3]
  ) +
  geom_text_repel(
    data = label.data,
    aes(label = Gene.names),
    size = 12,
    fontface = "bold",
    color = "white",
    bg.color = "black",
    bg.r = 0.10,
    max.overlaps = 100,
    nudge_y = 0.04,
    min.segment.length = 5
  ) + 
  annotate(
    geom = "text",
    y = 0.47 - 0.06 *(0:2),
    x = 1,
    hjust = 0,
    label = c(
      paste0("Interactions: ", stats$edges),
      paste0("Expected interactions: ", stats$lambda),
      paste0("\\U1D45D-Value: ", round(stats$enrichment, 5))
    ),
    size = 10
  ) +  
  theme_bw() +
  theme(panel.border = element_blank(),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent",colour = "black"),
        legend.position = c(0.15, 0.18),
        legend.text = element_text(size = 32),
        legend.title = element_text(size = 36),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
  ) + 
  guides(linetype = guide_legend(override.aes = list(size = 2)))

ggsave(
  "./Output_R/Interaction_direct_String.svg",
  width = 308,
  height = 308,
  unit = "mm",
  dpi = 300,
  bg = "transparent"
)

message("Done")
