############################
#HEADERS####################
############################
require(caper)
require(fulltext)
require(testdat)
#simple cap fixer (from R help!...)
.simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)), sep = "", collapse = " ")
}

############################
#BASIC DATA#################
############################
#Load migration data
data <- read.delim("mammals_data.csv", as.is=TRUE)
data$Migration <- sanitize_text(data$Migration)
data$Bin_omial <- with(data, paste(Genus, Species, sep="_"))
data$tax.lookup <- with(data, paste(Class, Order, Suborder, Family))
data$order <- sapply(data$Order, .simpleCap)
data$uncertain <- grepl("unknown", data$Migration, ignore.case=TRUE) | grepl("unclear", data$Migration, ignore.case=TRUE)
data <- data[data$uncertain == FALSE,]

#Remove columns we don't need
data <- data[,c("Locomotion", "Type", "tax.lookup", "Bin_omial", "order", "uncertain")]
names(data)[1:2] <- c("movement", "migration.type")
data$movement[data$movement==""] <- "W"
data$movement[data$movement=="A"] <- "W"
data$migration <- ifelse(data$migration.type=="N" | data$migration.type=="U", 0, 1)

#Load phylogenies and make ultrametric
trees <- scan("faurby_2015_posterior.tre", what="character")
trees <- lapply(trees, function(x) read.tree(text=x))
for(i in seq_along(trees)){
    terminal <- which(trees[[i]]$edge[,2] %in% seq_along(trees[[i]]$tip.label))
    tips <- length(trees[[i]]$tip.label)
    trees[[i]]$edge.length[terminal[order(trees[[i]]$edge[terminal,2])]] <- trees[[i]]$edge.length[terminal[order(trees[[i]]$edge[terminal,2])]] - (dist.nodes(trees[[i]])[1:tips,tips+1] - median(dist.nodes(trees[[i]])[1:tips,tips+1]))
}

#Migration status
data$refuge <- as.numeric(grepl("R", data$migration.type))
data$breeding <- as.numeric(grepl("B", data$migration.type))
data$tracking <- as.numeric(grepl("T", data$migration.type))

#Make comparative.data objects and save
c.datas <- lapply(trees, function(x) comparative.data(x, data, Bin_omial))
saveRDS(c.datas, "basic_workspace.RDS")

############################
#REGRESSION DATA############
############################
#Load
pantheria <- read.delim(ft_get_si("E090-184", "PanTHERIA_1-0_WR05_Aug2008.txt", "esa_archives"), row.names=NULL, as.is=TRUE)
pantheria <- pantheria[,c("MSW05_Binomial", "X5.1_AdultBodyMass_g", "X12.1_HabitatBreadth", "X6.2_TrophicLevel", "X6.1_DietBreadth")]
pantheria[pantheria==-999] <- NA
names(pantheria) <- c("Binomial", "mass", "habitat", "trophic", "diet")
redlist <- read.csv("red_list_ordered.csv")
redlist$Binomial <- with(redlist, paste(Genus, Species, sep=" "))
redlist$threatened <- !as.character(redlist$Red.List.status) %in% c("LC","NT")
pantheria$redlist <- redlist$threatened[match(pantheria$Binomial, redlist$Binomial)]
pantheria$Binomial <- gsub(" ", "_", pantheria$Binomial)
data <- merge(data, pantheria, by.x="Bin_omial", by.y="Binomial")

#Make comparative.data objects and save
c.datas <- lapply(trees, function(x) comparative.data(x, data, Bin_omial))
saveRDS(c.datas, "regression_workspace.RDS")
