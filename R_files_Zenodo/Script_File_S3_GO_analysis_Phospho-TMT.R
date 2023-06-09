## # # # # # # # # # # # # # # # # # # # # # # # #
##
## Project: CDKL5-NLS Phosphoproteomics
##
## Script name: Script_File_S3_GO_analysis_Phospho-TMT.R
##
## Purpose of script: Test for significant enrichment of GO terms
##
## Author: Florian Weiland
##
## Date Created: 2020-05-13
##
## # # # # # # # # # # # # # # # # # # # # # # # #

## Libraries ----

suppressPackageStartupMessages(library(GO.db, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(extrafont, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE))

theme_set(theme_classic(base_family = "Arial"))

## Read in data
message("Loading and preparing data")

data.bg <- read.csv(
  "./Output_R/Phospho-TMT_results_complete.csv"
)

data.fg <- read.csv(
  "./Output_R/Phospho-TMT_results_clean.csv"
)

## Unique protein entries only

data.bg.clean <- data.bg[!duplicated(data.bg$Parent.Entry), ]
data.fg.clean <- data.fg[!duplicated(data.fg$Parent.Entry), ]

## Get GO terms

go.fg <- strsplit(data.fg.clean$Gene.ontology.IDs, "; ")
go.bg <- strsplit(data.bg.clean$Gene.ontology.IDs, "; ")

uni.go.fg <- unique(unlist(strsplit(data.fg.clean$Gene.ontology.IDs, "; ")))

## Fisher's exact test on unique GO terms in foreground against background
message("Statistical testing")

go.df <- data.frame(
  GO.ID = character(),
  GO.Term = character(),
  Uniprot.Acc = character(),
  Gene.Name = character(),
  Foreground.Hits = numeric(),
  Background.Hits = numeric(),
  Foreground.Size = numeric(),
  Background.Size = numeric(),
  p.Value = numeric()
)

pb.1 <- txtProgressBar(min = 1, max = length(uni.go.fg), style = 3)
counter = 0

for (cur.go in uni.go.fg) {
  
  counter <- counter + 1
  
  count.fg <- length(which(grepl(cur.go, go.fg) == TRUE))
  count.bg <- length(which(grepl(cur.go, go.bg) == TRUE))
  
  fg.in.GO     <- count.fg
  fg.not.in.GO <- length(go.fg) - count.fg
  
  bg.in.GO <- count.bg
  bg.not.in.GO <- length(go.bg) - count.bg
  
  con.table <- matrix(nrow = 2, ncol = 2)
  rownames(con.table) <- c("Contains.GO", "Does.not.contain.GO")
  colnames(con.table) <- c("Foreground", "Background")
  
  con.table["Contains.GO", "Foreground"]         <- fg.in.GO
  con.table["Does.not.contain.GO", "Foreground"] <- fg.not.in.GO
  
  con.table["Contains.GO", "Background"]         <- bg.in.GO
  con.table["Does.not.contain.GO", "Background"] <- bg.not.in.GO
  
  p.temp <- fisher.test(con.table, alternative = "greater")$p.value
  
  ## Need to catch obsolete terms
  
  if (is.null(GOTERM[[cur.go]]) != TRUE) {
  
  go.df.temp <- data.frame(
    GO.ID = cur.go,
    Go.Term = GOTERM[[cur.go]]@Term,
    Uniprot.Acc = paste(data.fg.clean[grep(cur.go, data.fg.clean$Gene.ontology.IDs), "Parent.Entry" ], collapse = "; "),
    Gene.Name = paste(data.fg.clean[grep(cur.go, data.fg.clean$Gene.ontology.IDs), "Gene.names" ], collapse = "; "),
    Foreground.Hits = fg.in.GO,
    Background.Hits = bg.in.GO,
    Foreground.Size = length(go.fg),
    Background.Size = length(go.bg),
    p.Value = p.temp
  )
  
  } else {
    
    go.df.temp <- data.frame(
      GO.ID = cur.go,
      Go.Term = "NULL",
      Uniprot.Acc = paste(data.fg.clean[grep(cur.go, data.fg.clean$Gene.ontology.IDs), "Parent.Entry" ], collapse = "; "),
      Gene.Name = paste(data.fg.clean[grep(cur.go, data.fg.clean$Gene.ontology.IDs), "Gene.names" ], collapse = "; "),
      Foreground.Hits = fg.in.GO,
      Background.Hits = bg.in.GO,
      Foreground.Size = length(go.fg),
      Background.Size = length(go.bg),
      p.Value = p.temp
    )
    
    
  }
  
  go.df <- rbind(go.df, go.df.temp)
  
  setTxtProgressBar(pb.1, counter)
  
}

go.sig <- go.df[which(go.df$p.Value <= 0.01 & go.df$Foreground.Hits >= 3), ]
go.sig <- go.sig[order(go.sig$p.Value), ]
go.sig$GO.ID <- factor(go.sig$GO.ID, levels = c(go.sig[order(go.sig$p.Value), "GO.ID"]))

go.sig$Enrichment.Factor <- (go.sig$Foreground.Hits / go.sig$Foreground.Size) /
  (go.sig$Background.Hits / go.sig$Background.Size)

message("\\nPlotting GO terms")

ggplot(go.sig, aes(x = Enrichment.Factor, y = GO.ID)) +
  geom_col(width = 0.5, fill = "grey50", colour = "black", size = 1) +
  geom_text(
    aes(x = 0.05, label = paste0(Go.Term, " (", Foreground.Hits, " proteins)")),
    hjust = 0,
    size = 4,
    colour = "white"
    ) +
  xlab("GO term enrichment factor") +
  theme(
    axis.text.x = element_text(size = 14, family = "Arial"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

ggsave(
  "./Output_R/GO_analysis.svg",
  height = 105,
  width = 148,
  unit = "mm"
)


message("Writing results")

write.csv(
  go.sig,
  "./Output_R/GO_Analysis.csv",
  row.names = FALSE
)

message("Analysis done \\n")
