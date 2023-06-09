## PERMANOVA OF GUT CONTENT DATA
# Written by: Jun Ying Lim and Susan Kennedy
# From publication:
# Kennedy, S. R., Lim, J. Y., Clavel, J., Krehenwinkel, H., and Gillespie, R. G. 2019. Spider webs, stable isotopes and molecular gut content analysis: Multiple lines of evidence support trophic niche differentiation in a community of Hawaiian spiders. Functional Ecology.

#Set "main.dir" to your own directory containing subfolders:
main.dir <- setwd() #Insert your file path in parentheses
data.dir <- file.path(main.dir, "raw_data")
fig.dir <- file.path(main.dir, "figures")
res.dir <- file.path(main.dir, "results")

#### loading packages ####

library(vegan)

#### load data ####

# load order data
# Note that this dataset is randomly subsampled to 5 indivs per web species, and 6 indivs per Spiny Leg species
ord <- read.csv(file.path(data.dir, "gutES.csv"), header =  TRUE, stringsAsFactors = FALSE)

env_ord <- ord[ , 1:4]
ord <- ord[ , -c(1:4)]

# Broad stroke: sig diffs by clade?
adonis(ord ~ env_ord$clade)

# Sig diffs by species?
adonis(ord ~ env_ord$species)

# Pairwise PERMANOVA in web-builders
acueur <- c("acuta", "eurychasma")
ord_acueur <- ord[which(env_ord$species %in% acueur), ]
env_acueur <- env_ord[which(env_ord$species %in% acueur), ]

adonis(ord_acueur ~ env_acueur$species)


acufil <- c("acuta", "filiciphilia")
ord_acufil <- ord[which(env_ord$species %in% acufil), ]
env_acufil <- env_ord[which(env_ord$species %in% acufil), ]

adonis(ord_acufil ~ env_acufil$species)


acuste <- c("acuta", "stelarobusta")
ord_acuste <- ord[which(env_ord$species %in% acuste), ]
env_acuste <- env_ord[which(env_ord$species %in% acuste), ]

adonis(ord_acuste ~ env_acuste$species)


acutri <- c("acuta", "trituberculata")
ord_acutri <- ord[which(env_ord$species %in% acutri), ]
env_acutri <- env_ord[which(env_ord$species %in% acutri), ]

adonis(ord_acutri ~ env_acutri$species)


eurfil <- c("eurychasma", "filiciphilia")
ord_eurfil <- ord[which(env_ord$species %in% eurfil), ]
env_eurfil <- env_ord[which(env_ord$species %in% eurfil), ]

adonis(ord_eurfil ~ env_eurfil$species)


eurste <- c("eurychasma", "stelarobusta")
ord_eurste <- ord[which(env_ord$species %in% eurste), ]
env_eurste <- env_ord[which(env_ord$species %in% eurste), ]

adonis(ord_eurste ~ env_eurste$species)


eurtri <- c("eurychasma", "trituberculata")
ord_eurtri <- ord[which(env_ord$species %in% eurtri), ]
env_eurtri <- env_ord[which(env_ord$species %in% eurtri), ]

adonis(ord_eurtri ~ env_eurtri$species)


filste <- c("filiciphilia", "stelarobusta")
ord_filste <- ord[which(env_ord$species %in% filste), ]
env_filste <- env_ord[which(env_ord$species %in% filste), ]

adonis(ord_filste ~ env_filste$species)


filtri <- c("filiciphilia", "trituberculata")
ord_filtri <- ord[which(env_ord$species %in% filtri), ]
env_filtri <- env_ord[which(env_ord$species %in% filtri), ]

adonis(ord_filtri ~ env_filtri$species)


stetri <- c("stelarobusta", "trituberculata")
ord_stetri <- ord[which(env_ord$species %in% stetri), ]
env_stetri <- env_ord[which(env_ord$species %in% stetri), ]

adonis(ord_stetri ~ env_stetri$species)


# Pairwise PERMANOVA in Spiny Legs

brekam <- c("brevignatha", "kamakou")
ord_brekam <- ord[which(env_ord$species %in% brekam), ]
env_brekam <- env_ord[which(env_ord$species %in% brekam), ]

adonis(ord_brekam ~ env_brekam$species)


brequa <- c("brevignatha", "quasimodo")
ord_brequa <- ord[which(env_ord$species %in% brequa), ]
env_brequa <- env_ord[which(env_ord$species %in% brequa), ]

adonis(ord_brequa ~ env_brequa$species)


brewai <- c("brevignatha", "waikamoi")
ord_brewai <- ord[which(env_ord$species %in% brewai), ]
env_brewai <- env_ord[which(env_ord$species %in% brewai), ]

adonis(ord_brewai ~ env_brewai$species)


kamqua <- c("kamakou", "quasimodo")
ord_kamqua <- ord[which(env_ord$species %in% kamqua), ]
env_kamqua <- env_ord[which(env_ord$species %in% kamqua), ]

adonis(ord_kamqua ~ env_kamqua$species)


kamwai <- c("kamakou", "waikamoi")
ord_kamwai <- ord[which(env_ord$species %in% kamwai), ]
env_kamwai <- env_ord[which(env_ord$species %in% kamwai), ]

adonis(ord_kamwai ~ env_kamwai$species)



quawai <- c("quasimodo", "waikamoi")
ord_quawai <- ord[which(env_ord$species %in% quawai), ]
env_quawai <- env_ord[which(env_ord$species %in% quawai), ]

adonis(ord_quawai ~ env_quawai$species)

#load OTU data

# Note that this dataset is subsampled to 5 indivs per web species, and 6 indivs per Spiny Leg species
otu <- read.csv(file.path(data.dir, "otuES.csv"), header =  TRUE, stringsAsFactors = FALSE)

env_otu <- otu[ , 1:4]
otu <- otu[ , -c(1:4)]

# Broad stroke: sig diffs by clade?
adonis(otu ~ env_otu$clade)

# Sig diffs by species?
adonis(otu ~ env_otu$species)

# Pairwise PERMANOVA in web-builders
acueur <- c("acuta", "eurychasma")
otu_acueur <- otu[which(env_otu$species %in% acueur), ]
env_acueur <- env_otu[which(env_otu$species %in% acueur), ]

adonis(otu_acueur ~ env_acueur$species)


acufil <- c("acuta", "filiciphilia")
otu_acufil <- otu[which(env_otu$species %in% acufil), ]
env_acufil <- env_otu[which(env_otu$species %in% acufil), ]

adonis(otu_acufil ~ env_acufil$species)


acuste <- c("acuta", "stelarobusta")
otu_acuste <- otu[which(env_otu$species %in% acuste), ]
env_acuste <- env_otu[which(env_otu$species %in% acuste), ]

adonis(otu_acuste ~ env_acuste$species)


acutri <- c("acuta", "trituberculata")
otu_acutri <- otu[which(env_otu$species %in% acutri), ]
env_acutri <- env_otu[which(env_otu$species %in% acutri), ]

adonis(otu_acutri ~ env_acutri$species)


eurfil <- c("eurychasma", "filiciphilia")
otu_eurfil <- otu[which(env_otu$species %in% eurfil), ]
env_eurfil <- env_otu[which(env_otu$species %in% eurfil), ]

adonis(otu_eurfil ~ env_eurfil$species)


eurste <- c("eurychasma", "stelarobusta")
otu_eurste <- otu[which(env_otu$species %in% eurste), ]
env_eurste <- env_otu[which(env_otu$species %in% eurste), ]

adonis(otu_eurste ~ env_eurste$species)


eurtri <- c("eurychasma", "trituberculata")
otu_eurtri <- otu[which(env_otu$species %in% eurtri), ]
env_eurtri <- env_otu[which(env_otu$species %in% eurtri), ]

adonis(otu_eurtri ~ env_eurtri$species)


filste <- c("filiciphilia", "stelarobusta")
otu_filste <- otu[which(env_otu$species %in% filste), ]
env_filste <- env_otu[which(env_otu$species %in% filste), ]

adonis(otu_filste ~ env_filste$species)


filtri <- c("filiciphilia", "trituberculata")
otu_filtri <- otu[which(env_otu$species %in% filtri), ]
env_filtri <- env_otu[which(env_otu$species %in% filtri), ]

adonis(otu_filtri ~ env_filtri$species)


stetri <- c("stelarobusta", "trituberculata")
otu_stetri <- otu[which(env_otu$species %in% stetri), ]
env_stetri <- env_otu[which(env_otu$species %in% stetri), ]

adonis(otu_stetri ~ env_stetri$species)


# Pairwise PERMANOVA in Spiny Legs

brekam <- c("brevignatha", "kamakou")
otu_brekam <- otu[which(env_otu$species %in% brekam), ]
env_brekam <- env_otu[which(env_otu$species %in% brekam), ]

adonis(otu_brekam ~ env_brekam$species)


brequa <- c("brevignatha", "quasimodo")
otu_brequa <- otu[which(env_otu$species %in% brequa), ]
env_brequa <- env_otu[which(env_otu$species %in% brequa), ]

adonis(otu_brequa ~ env_brequa$species)


brewai <- c("brevignatha", "waikamoi")
otu_brewai <- otu[which(env_otu$species %in% brewai), ]
env_brewai <- env_otu[which(env_otu$species %in% brewai), ]

adonis(otu_brewai ~ env_brewai$species)


kamqua <- c("kamakou", "quasimodo")
otu_kamqua <- otu[which(env_otu$species %in% kamqua), ]
env_kamqua <- env_otu[which(env_otu$species %in% kamqua), ]

adonis(otu_kamqua ~ env_kamqua$species)


kamwai <- c("kamakou", "waikamoi")
otu_kamwai <- otu[which(env_otu$species %in% kamwai), ]
env_kamwai <- env_otu[which(env_otu$species %in% kamwai), ]

adonis(otu_kamwai ~ env_kamwai$species)



quawai <- c("quasimodo", "waikamoi")
otu_quawai <- otu[which(env_otu$species %in% quawai), ]
env_quawai <- env_otu[which(env_otu$species %in% quawai), ]

adonis(otu_quawai ~ env_quawai$species)

