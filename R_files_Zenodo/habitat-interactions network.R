#load packages
library(reshape2)
library(plyr)
library(dplyr)
library(bipartite)
library(vegan)
library(glmmTMB)
library(emmeans)
library(igraph)
library(data.table)
library(ggplot2)
library(DHARMa)
library(multcompView)
library(betalink)
library(rpgm)
library(splitstackshape)
library(MuMIn)
library(ggpubr)

##################################
#Read in data----
##################################

network <- read.csv("data/Pollen networks.csv", header=TRUE)
land.use <- read.csv("data/Land_use_Pollinators.csv", header=TRUE)

##################################
#data wrangling for metrics----
##################################

#remove insects with no pollen
network$pollen.found[network$pollen.found == 0] <- NA
network <- na.omit(network)
network.melt <- melt(network, id.vars=c(1,3,8), measure.vars=c(14:29))
network.melt$value[network.melt$value == 0] <- NA
network.melt <- na.omit(network.melt)
colnames(network.melt)[3] <- "pollinator"
colnames(network.melt)[4] <- "plant"
#network.melt <-network.melt[network.melt$plant!="Poaceae",]#this removes poaceae...

#combine plants and pollinators
network.melt$int <- paste(network.melt$pollinator,"_",network.melt$plant)

binary.links <- network.melt %>%
  group_by(Site, Land.use) %>%
  summarise(links.int = n_distinct(int))

#interactions unique to only one site
xy <- dcast(network.melt, Site ~ int)
xy2 <- xy[,-1]
rownames(xy2) <- xy[,1]
xy2[xy2>0] <-1
xy3 = xy2[,colSums(xy2) < 2]
xy4 <- cbind(xy3, rowSums(xy3))
xy4 <- xy4[64]
colnames(xy4)[1] <- "number.unique.interactions"
unique.ints <- setDT(as.data.frame(xy4), keep.rownames = TRUE)[]
colnames(unique.ints)[1] <- "Site"

binary.links <- merge(binary.links, unique.ints)

binary.links$ratio <-  binary.links$number.unique.interactions/binary.links$links.int

#cast data by int
network.int <- dcast(network.melt, Site + Land.use ~ int, fun.aggregate = sum, na.rm =T, value.var="value", fill = 0)
network.plant.nodf <- dcast(network.melt, Land.use ~ plant, fun.aggregate = sum, na.rm =T, value.var="value", fill = 0)
network.pol.nodf <- dcast(network.melt, Land.use ~ pollinator, fun.aggregate = sum, na.rm =T, value.var="value", fill = 0)

#cast by landuse only
network.int.comb <- dcast(network.melt, Land.use ~ int, fun.aggregate = sum, na.rm =T, value.var="value", fill = 0)

######################################
#plot the network data----
######################################

network.int1 <- network.int[,-c(1:2)]
network.int2 <- as.matrix(network.int1)
patch <- unique(network.int[, c("Site", "Land.use")])
patch$color <- NA

#define landuse colour
patch[c(1:3,11,22), "color"] <- "gold"
patch[c(4,8:10,12,15), "color"] <- "blue"
patch[c(5:6,14,16,18,23), "color"] <- "red"
patch[c(7,13,17,19:21), "color"] <- "green"
rownames(network.int2) <- patch$Site

#prepare legend and make plot
legend <- unique(patch[,c("Land.use", "color")])
par(xpd = T) #allow plotting outside the plot
plotweb(network.int2, col.low = as.character(patch$color))
legend(x=0, y=0.25, as.character(legend$Land.use), pch=21,
       col="#777777", pt.bg=as.character(legend$color),
       pt.cex=1, cex=.6, bty="n", ncol=2)

#try visweb
visweb(network.int2, prednames = T, preynames = T, labsize = 0.6)

#prepare the data for igraph
links <- network.melt[,c("int", "Site", "value")]
colnames(links)[3] <- "weight"
node1 <-  unique(network.melt[,c("Site", "Land.use")])
colnames(node1) <- c("node", "attribute")
node1$type <- "habitat"

#define landuse colour
patch2 <- node1
patch2$color <- NA
patch2$rn <- 1:23
patch2[c(1,6,7,9,13,18), "color"] <- "gold"
patch2[c(2,3,5,8,12,17), "color"] <- "blue"
patch2[c(4,10,11,15,23), "color"] <- "red"
patch2[c(14,16,19:22), "color"] <- "green"

node2 <-  data.frame(node = unique(network.melt[,c("int")]),
                     attribute = NA,
                     type = "species")
nodes <- rbind(node1, node2)

#create igraph object
net <- graph_from_data_frame(d=links,
                             vertices=nodes, directed=F)
# Generate colors based habitat:
clrs <- data.frame(nod = V(net)$attribute,
                   cols = c(patch2$color, rep("pink", 99)))
V(net)$color <- as.character(clrs$cols)
# Compute node degrees (#links) and use that to set node size:
deg <- degree(net, mode="all")
V(net)$size <- deg*0.3
# Setting them to NA will render no labels:
V(net)$label <- as.character(nodes$node)
# Set edge width based on weight:
E(net)$width <- 1 #E(net)$weight/30
#change arrow size and edge color:
E(net)$arrow.size <- .2 #but note no arrows in Unidirected graphs like this
E(net)$edge.color <- "black"
#prepare colors
cl <- unique(clrs)
cl$nod <- as.character(cl$nod)
cl$nod[which(is.na(cl$nod))] <- "interactions"
plot(net, vertex.label = ifelse(degree(net) > 10, V(net)$label, NA))#plots nodes with 10+ links
legend(x=-1.5, y=-1.1, cl$nod, pch=21,
       col="#777777", pt.bg=as.character(cl$cols),
       pt.cex=2, cex=.8, bty="n", ncol=2)

###########################################
#nestedness analysis----
#are plants, pollinators and interactions in less intensive landuses nested within more intensive landuses?
###########################################

#create landuse "gradient' lowest intensity to highest intensity
gradient <- c("Forest", "Avocado", "Dairy", "Cropping")

#reorder networks based on gradient above
#interactions
network.int.comb <- network.int.comb %>%
  mutate(Land.use =  factor(Land.use, levels = gradient)) %>%
  arrange(Land.use)

#plants
network.plant.nodf <- network.plant.nodf %>%
  mutate(Land.use =  factor(Land.use, levels = gradient)) %>%
  arrange(Land.use)

#pollinators
network.pol.nodf <- network.pol.nodf %>%
  mutate(Land.use =  factor(Land.use, levels = gradient)) %>%
  arrange(Land.use)

#run null model
#Null model II from Bascompte et al. (2003). Creates random networks by probabilistically
#fixing row and column marginal totals. The expected number of links is same as observed number.
#Rodriguez-Girona and Santamaria (2006) showed that this null model has the best compromise between Type I and Type II errors
# Run this code once to create the null model function
null.model.II <- function(web){
  web <- as.matrix(web > 0) + 0
  # calculate the probability based on row marginals. Creates matrix same size as web, with row sums divided by number of columns (to get probability of a 1 in each cell of each row), repeated across all columns for each row.
  row.probs <- matrix(rowSums(web)/ncol(web),nrow(web),ncol(web))
  # calculate the probability based on column marginals (again, repeated for whole column). Transpose used instead of byrow=T
  col.probs <- t(matrix(colSums(web)/nrow(web),ncol(web),nrow(web)))
  # calculate the element by element mean of this probabilities
  mat.probs <- (row.probs + col.probs) / 2.0
  # generate a random matrix with 1s proportional to the above probabilities. rbinom(n, size, prob) n is number of observations, size is number of trials, prob is prob of success in each trial
  mat.null <- matrix(rbinom(nrow(web)*ncol(web),1,as.vector(mat.probs)),nrow(web),ncol(web))
  # return that matrix in all its glory
  return(mat.null)
}

#Begin permutation test (two tailed)
reps <- 9999#set number of permutations

#Create a list with spaces for each output matrix
nulls<-vector("list",reps)
for (i in 1:reps) {
  nulls[[i]]<-null.model.II(network.int.comb[,-1])
}

#call any individual matrix from that list using nulls[[x]], where x is the number of the matrix you want to call
null.NODF <- matrix(, nrow=reps, ncol=1)
for (i in 1:reps) {
  null.NODF[i, ] <- nestednodf(nulls[[i]], order = FALSE, weighted = FALSE, wbinary = TRUE)$statistic[3]
}

#Calculate NODF of observed web
int.nest <- as.vector(nestednodf(network.int.comb[,-1], order = FALSE, wbinary=TRUE, weighted=FALSE))$statistic[3]

NODF.int <- c(null.NODF)
#generate p-value
pval <- sum(NODF.int>int.nest) / length(NODF.int)
ifelse(pval > 0.5, 1-pval, pval)

sd <- apply(null.NODF, 2, sd)#calculate standard deviation
z_int <- ((int.nest - colMeans(null.NODF))/sd)#can use this plot plot z score

#now do the same for plants
#Create a list with spaces for each output matrix
nulls<-vector("list",reps)
for (i in 1:reps) {
  nulls[[i]]<-null.model.II(network.plant.nodf[,-1])
}

#call any individual matrix from that list using nulls[[x]], where x is the number of the matrix you want to call
null.NODF <- matrix(, nrow=reps, ncol=1)
for (i in 1:reps) {
  null.NODF[i, ] <- nestednodf(nulls[[i]], order = FALSE, weighted = FALSE, wbinary = TRUE)$statistic[3]
}

#Calculate NODF of observed web
plant.nest <- as.vector(nestednodf(network.plant.nodf[,-1], order = FALSE, weighted=FALSE, wbinary=TRUE))$statistic[3]

NODFs_plant <- c(null.NODF)
#generate p-value
pval_plant <- sum(NODFs_plant>plant.nest) / length(NODFs_plant)
ifelse(pval_plant > 0.5, 1-pval_plant, pval_plant)

sd <- apply(null.NODF, 2, sd)#calculate standard deviation
z_plant <- ((plant.nest - colMeans(null.NODF))/sd)#can use this plot plot z score

#now do the same for pollinators
#Create a list with spaces for each output matrix
nulls<-vector("list",reps)
for (i in 1:reps) {
  nulls[[i]]<-null.model.II(network.pol.nodf[,-1])
}

#call any individual matrix from that list using nulls[[x]], where x is the number of the matrix you want to call
null.NODF <- matrix(, nrow=reps, ncol=1)
for (i in 1:reps) {
  null.NODF[i, ] <- nestednodf(nulls[[i]], order = FALSE, weighted = FALSE, wbinary = TRUE)$statistic[3]
}

#Calculate NODF of observed web
pol.nest <- as.vector(nestednodf(network.pol.nodf[,-1], order = FALSE, weighted=FALSE, wbinary=TRUE))$statistic[3]

NODFs_pol <- c(null.NODF)
#generate p-value
pval_pol <- sum(NODFs_pol>pol.nest) / length(NODFs_pol)
ifelse(pval_pol > 0.5, 1-pval_pol, pval_pol)

sd <- apply(null.NODF, 2, sd)#calculate standard deviation
z_pol <- ((pol.nest - colMeans(null.NODF))/sd)#can use this plot plot z score

#combine all z scores in a single dataframe
lu <- as.data.frame(cbind(z_pol,z_int,z_plant))
lu$scale <- "Site-level land-use"
setnames(lu, old = c('z_pol','z_int','z_plant'), new = c('Pollinator','Interaction','Plant'))
div <- as.data.frame(cbind(z_pla.div,z_int.div,z_pol.div))
div$scale <- "Landscape-level habitat diversity"
setnames(div, old = c('z_pla.div','z_int.div','z_pol.div'), new = c('Plant','Interaction','Pollinator'))
fr <- as.data.frame(cbind(z_pol.for,z_pla.for,z_int.for))
fr$scale <- "Landscape-level forest cover"
setnames(fr, old = c('z_pol.for','z_pla.for','z_int.for'), new = c('Pollinator','Plant','Interaction'))

#combine the dataframes
zval.all <- rbind(lu,div,fr)
zval.all <- melt(zval.all)
zval.all$variable <- factor(zval.all$variable, levels=c("Pollinator", "Plant", "Interaction"))
zval.all$scale <- factor(zval.all$scale, levels=c("Site-level land-use", "Landscape-level forest cover", "Landscape-level habitat diversity"))

#plot NODF zscores
p <- ggplot()
p <- p + xlab(NULL) + ylab("NODF (Z score)")
p <- p + geom_point(data=zval.all, aes(x=variable, y=value, color=variable),
                     alpha=1, size=3.5)
p <- p + geom_hline(yintercept = 1.96, linetype="dashed")
p <- p + geom_hline(yintercept = -1.96, linetype="dashed")
p <- p + geom_hline(yintercept = 0)
p <- p +  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96, fill = "grey", alpha = .5, color = NA)
p <- p + facet_wrap(~scale)
p <- p + theme(panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black")) +
  theme(panel.border=element_rect(colour = "black", fill = "NA", size = 1)) +
  theme(axis.text.x=element_text(angle= 90, hjust = 0.5, vjust = 0.5, size =16),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =16),
        axis.title.y=element_text(size=24, vjust = 1),
        axis.title.x=element_text(size=24, vjust = 0.5),
        axis.text=element_text(colour = "black"))+
  theme(strip.background = element_rect(colour="NA", fill=NA),
        strip.text = element_text(size=12)) +
  theme(axis.ticks.length = unit(2, "mm"))
p <- p + theme(panel.spacing.x=unit(1, "lines"),panel.spacing.y=unit(1, "lines"))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p <- p + theme(legend.position="none")
p

#combine null and true NODF values for plotting
#not sure these density plots are necessary so have only plotted for sites
#this only works for most recent in the script (pollinators). need to rename to run for all three LU analyses
nest.all <- as.data.frame(cbind(int.nest,pol.nest,plant.nest))
nest.all <- melt(nest.all)
nest.all$variable <- revalue(nest.all$variable, c("int.nest"="NODF.int", "pol.nest"="NODFs_pol", "plant.nest"="NODFs_plant"))
colnames(nest.all)[2] <- "NODF_true"
null.all <- as.data.frame(cbind(NODFs_pol,NODFs_plant,NODF.int))
null.all <- melt(null.all)
null.all <- merge(null.all,nest.all, by="variable")
null.all <- null.all %>%
  group_by(variable) %>%
  mutate(sd = sd(value),
            error = qnorm(0.975)*sd,
            ci.low = mean(value) - error,
            ci.upp = mean(value) + error)

#plot density distribution
p <- ggplot()
p <- p + xlab("NODF") + ylab("Density")
p <- p + theme(text = element_text(size=18))
p <- p + geom_density(data=null.all, aes(value))
p <- p + facet_wrap(~variable, scales="free")
p <- p + geom_vline(data=null.all, aes(xintercept = NODF_true), linetype="dashed", colour="red")
p <- p + geom_vline(data=null.all, aes(xintercept = ci.low), linetype="dashed", colour="blue")
p <- p + geom_vline(data=null.all, aes(xintercept = ci.upp), linetype="dashed", colour="blue")
p <- p + theme(panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black")) +
  theme(panel.border=element_rect(colour = "black", fill = "NA", size = 1)) +
  theme(axis.text.x=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =18),
        axis.title.x=element_text(size=30, vjust = 1),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =18),
        axis.title.y=element_text(size=30, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(strip.background = element_rect(colour="NA", fill=NA),
        strip.text = element_text(size=20))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p

###########################################
#modularity----
###########################################
#c : among-module connectivity
#z : within-module connectivity

#Compute the modules
res <- computeModules(network.int2)

#plot the modules
plotModuleWeb(res, displayAlabels = T)#not really useful

#null models as per Dormann & Strauss for claculating critical c & z values
nullnet <- nullmodel(network.int2, N=100, method=3)#using vaznull
modules.nulls <- sapply(nullnet, computeModules)
cz.nulls <- sapply(modules.nulls, czvalues, weighted=TRUE,level="lower")
cz.nulls.int <- sapply(modules.nulls, czvalues, weighted=TRUE,level="higher")

#calculate c critical
c.nulls <- as.data.frame(unlist(cz.nulls[1,]))
colnames(c.nulls)[1] <- "cval"
quantile(c.nulls$cval,probs=c(0.975))

c.nulls.int <- as.data.frame(unlist(cz.nulls.int[1,]))
colnames(c.nulls.int)[1] <- "cval"
quantile(c.nulls.int$cval,probs=c(0.975))

#calculate z critical
z.nulls <- as.data.frame(unlist(cz.nulls[2,]))
colnames(z.nulls)[1] <- "zval"
z.nulls <- na.omit(z.nulls)
quantile(z.nulls$zval,probs=c(0.975))

z.nulls.int <- as.data.frame(unlist(cz.nulls.int[2,]))
colnames(z.nulls.int)[1] <- "zval"
z.nulls.int <- na.omit(z.nulls.int)
quantile(z.nulls.int$zval,probs=c(0.975))

#is the network more modular than expected?
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
(z <- (res@likelihood - mean(like.nulls))/sd(like.nulls))

#calculate actual c and z values for each node (for interactions and sites)
cz <- as.data.frame(czvalues(res, weighted = TRUE, level = "lower"))
cz_interactions <- as.data.frame(czvalues(res, weighted = TRUE, level = "higher"))
setDT(cz, keep.rownames = TRUE)[]
setDT(cz_interactions, keep.rownames = TRUE)[]
colnames(cz)[1] <- "Site"
colnames(cz_interactions)[1] <- "Interaction"

#############################
#calculate site strength----
#############################

#transpose the matrix to calculate strength for the lower level of the network
patch$strength <- bipartite::strength(t(network.int2), type="Bascompte")

#need to install dev version of glmmTMB to link with emmeans
#dev_mode(on=T)
#install_github("glmmTMB/glmmTMB/glmmTMB@effects")
#dev_mode(on=F)
library(glmmTMB)

#run model for strength
strength.mod <- glmmTMB(strength ~ Land.use + (1|Site),
                                family=gaussian(link=log),
                                data=patch)

#plot residuals
plot(residuals(stength.mod, "pearson"))

#calculate lsmeans
stength.ls <- emmeans(strength.mod, pairwise ~ Land.use, level = .95, adjust = "fdr", int.adjust="mvt", type = "response")

#generate letters for groups
stength.CLD <- CLD(stength.ls, Letters = letters, level = .95, adjust = "fdr")
stength.CLD <- stength.CLD[,c(1,7)]
max.stength <- patch %>% group_by(Land.use) %>% summarise(max = max(strength+1))
stength.CLD <- merge(max.stength, stength.CLD)
stength.CLD <- merge(stength.ls$emmeans, stength.CLD)

#run model for number of links
links.mod <- glmmTMB(links.int ~ Land.use + (1|Site),
                        truncated_poisson(link = "log"),
                        data=binary.links)
plot(residuals(links.mod, "pearson"))

#calculate lsmeans
links.ls <- emmeans(links.mod, pairwise ~ Land.use, level = .95, adjust = "fdr", int.adjust="none", type = "response")

#generate letters for groups
links.CLD <- CLD(links.ls, Letters = letters, level = .95, adjust = "fdr")
links.CLD <- links.CLD[,c(1,7)]
max.links <- binary.links %>% group_by(Land.use) %>% summarise(max = max(links.int+3))
links.CLD <- merge(max.links, links.CLD)
links.CLD <- merge(links.ls$emmeans, links.CLD)
links.CLD$max[3] <- 25.16

uni.int.mod <- glmmTMB(number.unique.interactions ~ Land.use + (1|Site),
              family=poisson(link="log"),
              data=binary.links)

#calculate lsmeans
uni.ls <- emmeans(uni.int.mod, pairwise ~ Land.use, level = .95, adjust = "fdr", int.adjust="none", type = "response")

#generate letters for groups
uni.CLD <- CLD(uni.ls, Letters = letters, level = .95, adjust = "fdr")
uni.CLD <- uni.CLD[,c(1,7)]
max.uni <- binary.links %>% group_by(Land.use) %>% summarise(max = max(number.unique.interactions+1))
uni.CLD <- merge(max.uni, uni.CLD)
uni.CLD <- merge(uni.ls$emmeans, uni.CLD)

#run models comparing various landuse metrics to strength and number of links
#extract relevant columns
for.div <- land.use[c(1,3,5)]
for.div <- unique(for.div)

#merge wih number of links and weighted strength
binary.landscape <- merge(for.div,binary.links)
binary.landscape <- merge(binary.landscape,patch)

#run model for number of links
links.bl.mod <- glmmTMB(links.int ~ LU_Div+percent_Treecover+Land.use + (1|Site),
                        family=truncated_poisson(link="log"),
                        data=binary.landscape)
plot(residuals(links.bl.mod, "pearson"))

links.bl.dredge <- dredge(links.bl.mod)

#write file to git folder
write.csv(links.bl.dredge, "/Users/macuser/Library/Mobile Documents/com~apple~CloudDocs/H_drive_DT/Mareeba_Pollen_Network/no.links_model_selection.csv")

#run model for number distinct interactions
dist.mod <- glmmTMB(number.unique.interactions ~ LU_Div+percent_Treecover+Land.use + (1|Site),
                        family=poisson(link="log"),
                        data=binary.landscape)
plot(residuals(dist.mod, "pearson"))

dist.dredge <- dredge(dist.mod)

#write file to git folder
write.csv(dist.dredge, "/Users/macuser/Library/Mobile Documents/com~apple~CloudDocs/H_drive_DT/Mareeba_Pollen_Network/no.distinct_interactions_model_selection.csv")

#run model for strength
stength.bl.mod <- glmmTMB(strength ~ LU_Div+percent_Treecover+Land.use + (1|Site),
                          family=Gamma(link="log"),
                          data=binary.landscape)
plot(residuals(stength.bl.mod, "pearson"))

stength.bl.dredge <- dredge(stength.bl.mod)

#merge modularity and landuse type dfs for ploting
patch_plot <- merge(cz,patch)
#write file to git folder
write.csv(stength.bl.dredge, "/Users/macuser/Library/Mobile Documents/com~apple~CloudDocs/H_drive_DT/Mareeba_Pollen_Network/strength_model_selection.csv")

#reorder landuse for plotting
new.order <- c("Forest","Avocado","Dairy", "Cropping")
stength.CLD$Land.use <- factor(as.character(stength.CLD$Land.use), levels=new.order)
stength.CLD <- stength.CLD[order(stength.CLD$Land.use),]
patch_plot$Land.use <- factor(as.character(patch_plot$Land.use), levels=new.order)
patch_plot <- patch_plot[order(patch_plot$Land.use),]
links.CLD$Land.use <- factor(as.character(links.CLD$Land.use), levels=new.order)
links.CLD <- links.CLD[order(links.CLD$Land.use),]
binary.links$Land.use <- factor(as.character(binary.links$Land.use), levels=new.order)
binary.links <- binary.links[order(binary.links$Land.use),]
uni.CLD$Land.use <- factor(as.character(uni.CLD$Land.use), levels=new.order)
uni.CLD <- uni.CLD[order(uni.CLD$Land.use),]
binary.links$Land.use <- factor(as.character(binary.links$Land.use), levels=new.order)
binary.links <- binary.links[order(binary.links$Land.use),]

#plot strength
p <- ggplot()
p <- p + xlab(NULL) + ylab("Strength (sum site dependencies)")
p <- p + geom_point(data=stength.CLD, aes(x=Land.use, y=response, color=Land.use, fill=Land.use),
            alpha=1, size=5)
p <- p + geom_errorbar(data=stength.CLD, aes(x=Land.use, ymin=lower.CL, ymax=upper.CL, color=Land.use),
       width = 0)
p <- p + geom_point(data=patch_plot, aes(x=Land.use, y=strength, color=Land.use, fill=Land.use),
            size=2.5, shape=1)
p <- p + geom_text(data = stength.CLD, aes(x = Land.use, y=max, label=.group), size =5)
p <- p + geom_text(data = patch_plot, aes(x=Land.use, y=strength, label=Site), size =3, hjust=-0.25, vjust=0)
p <- p + theme(axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text=element_text(colour = "black"))+
  theme(axis.ticks.length = unit(2, "mm"),
        axis.ticks = element_line(colour = 'black', size = 0.4))+
  theme(axis.text.x=element_text(angle= 90, hjust = 1, vjust = 0.5, size =14),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =14),
        axis.title.y=element_text(size=18, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(strip.background = element_rect(colour="NA", fill=NA))
p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.4))
p <- p + theme(panel.spacing.x=unit(1, "lines"),panel.spacing.y=unit(1, "lines"))
p <- p + theme(axis.title.y=element_text(margin=margin(0,0,0,0)))
p <- p + scale_color_brewer(palette="Set1")
p <- p + scale_fill_brewer(palette="Set1")
p1 <- p + theme(legend.position="none")
p1

#plot number of interaction links
p <- ggplot()
p <- p + xlab(NULL) + ylab("Plant-pollinator interaction richness")
p <- p + geom_point(data=links.CLD, aes(x=Land.use, y=rate, color=Land.use, fill=Land.use),
                    alpha=1, size=5)
p <- p + geom_errorbar(data=links.CLD, aes(x=Land.use, ymin=lower.CL, ymax=upper.CL, color=Land.use),
                       width = 0)
p <- p + geom_point(data=binary.links, aes(x=Land.use, y=links.int, color=Land.use, fill=Land.use),
                    size=2.5, shape=1)
p <- p + geom_text(data = links.CLD, aes(x = Land.use, y=max, label=.group), size =5)
p <- p + geom_text(data = binary.links, aes(x=Land.use, y=links.int, label=Site), size =3, hjust=-0.25, vjust=0)
p <- p + theme(axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text=element_text(colour = "black"))+
  theme(axis.ticks.length = unit(2, "mm"),
        axis.ticks = element_line(colour = 'black', size = 0.4))+
  theme(axis.text.x=element_text(angle= 90, hjust = 1, vjust = 0.5, size =14),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =14),
        axis.title.y=element_text(size=18, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(strip.background = element_rect(colour="NA", fill=NA))
p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.4))
p <- p + theme(panel.spacing.x=unit(1, "lines"),panel.spacing.y=unit(1, "lines"))
p <- p + theme(axis.title.y=element_text(margin=margin(0,0,0,0)))
p <- p + scale_color_brewer(palette="Set1")
p <- p + scale_fill_brewer(palette="Set1")
p2 <- p + theme(legend.position="none")
p2

#plot number of unique interactions
p <- ggplot()
p <- p + xlab(NULL) + ylab("Number of unique plant-pollinator interactions")
p <- p + geom_point(data=uni.CLD, aes(x=Land.use, y=rate, color=Land.use, fill=Land.use),
                    alpha=1, size=5)
p <- p + geom_errorbar(data=uni.CLD, aes(x=Land.use, ymin=lower.CL, ymax=upper.CL, color=Land.use),
                       width = 0)
p <- p + geom_point(data=binary.links, aes(x=Land.use, y=number.unique.interactions, color=Land.use, fill=Land.use),
                    size=2.5, shape=1)
p <- p + geom_text(data = uni.CLD, aes(x = Land.use, y=max, label=.group), size =5)
p <- p + geom_text(data = binary.links, aes(x=Land.use, y=number.unique.interactions, label=Site), size =3, hjust=-0.25, vjust=0)
p <- p + theme(axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text=element_text(colour = "black"))+
  theme(axis.ticks.length = unit(2, "mm"),
        axis.ticks = element_line(colour = 'black', size = 0.4))+
  theme(axis.text.x=element_text(angle= 90, hjust = 1, vjust = 0.5, size =14),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =14),
        axis.title.y=element_text(size=18, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(strip.background = element_rect(colour="NA", fill=NA))
p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.4))
p <- p + theme(panel.spacing.x=unit(1, "lines"),panel.spacing.y=unit(1, "lines"))
p <- p + theme(axis.title.y=element_text(margin=margin(0,0,0,0)))
p <- p + scale_color_brewer(palette="Set1")
p <- p + scale_fill_brewer(palette="Set1")
p3 <- p + theme(legend.position="none")
p3

#merge plots
figure1 <- ggarrange(p2, p3, p1, ncol = 3, nrow = 1)
figure1
ggsave(figure1,file="graphs/Figs for paper_V1/figure-1.pdf", device = "pdf",dpi=320,width=14,height=10,units = c("in"))
ggsave(figure1,file="graphs/Figs for paper_V1/figure-1.png", device = "png",dpi=320,width=14,height=10,units = c("in"))

#plot modularity for sites
p <- ggplot(patch_plot, aes(y=z, x=c))
p <- p + xlab("Among-module connectivity") + ylab("Within-module connectivity")
p <- p + theme(text = element_text(size=18))
p <- p + geom_point(size=2.5, aes(color=Land.use))
p <- p + geom_text(aes(label=Site),hjust=0, vjust=0)
p <- p + geom_hline(yintercept=2.189251, linetype="dashed")
p <- p + geom_vline(xintercept=0.6207565, linetype="dashed")
p <- p + theme(panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black")) +
  theme(panel.border=element_rect(colour = "black", fill = "NA", size = 1)) +
  theme(axis.text.x=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =16),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =16),
        axis.title.y=element_text(size=24, vjust = 1),
        axis.title.x=element_text(size=24, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(strip.background = element_rect(colour="NA", fill=NA),
        strip.text = element_text(size=20))
p <- p + theme(panel.spacing.x=unit(1, "lines"),panel.spacing.y=unit(1, "lines"))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p <- p + scale_color_brewer(palette="Set1")
p

#plot modularity for interactions
p <- ggplot(cz_interactions, aes(y=z, x=c))
p <- p + xlab("Among-module connectivity") + ylab("Within-module connectivity")
p <- p + theme(text = element_text(size=18))
p <- p + geom_point(size=2.5)
p <- p + geom_text(aes(label=Interaction),hjust=0, vjust=0)
p <- p + geom_hline(yintercept=3.666025 , linetype="dashed")
p <- p + geom_vline(xintercept=0.6846496, linetype="dashed")
p <- p + scale_x_continuous(limits = c(0.0001, 1.1))
p <- p + scale_y_continuous(limits = c(0.0001, 4.5))
p <- p + theme(panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black")) +
  theme(panel.border=element_rect(colour = "black", fill = "NA", size = 1)) +
  theme(axis.text.x=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =16),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =16),
        axis.title.y=element_text(size=24, vjust = 1),
        axis.title.x=element_text(size=24, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(strip.background = element_rect(colour="NA", fill=NA),
        strip.text = element_text(size=20))
p <- p + theme(panel.spacing.x=unit(1, "lines"),panel.spacing.y=unit(1, "lines"))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p <- p + scale_color_brewer(palette="Set1")
p

##################################################
#calculate interaction turnover across sites----
#################################################

int.matrices <- with(network.melt, tapply(value, list(plant, pollinator, Site),sum))  # Use tapply to create a set of bird-plant interaction matrices.
#use length to count total number of transects during which that interaction was observed, or sum to count total number of fruits eaten.
int.matrices[is.na(int.matrices)]<-0  # Replace NA values with zero

apply(int.matrices,1,sum) # See how many total observations are in each web

# Create a function to count the number of non-zero cells in the matrix (i.e. unique links)
nonzero<-function(x){
  return(length(which(x>0)))
}

apply(int.matrices,3,nonzero) # For each site, return the number of unique links

links<-apply(int.matrices,3,nonzero) # Put above result into an object. Note that 3 is the third dimension of the array, so don't change that number.
minsize<-names(links[links>1]) # Get list of sites with greater than a certain number of links (currently 1)
min.matrices<-int.matrices[,,minsize] # Retreive only those interaction matrices over the minimum size threshold

min.mat.list<-lapply(seq(dim(min.matrices)[3]), function(x) min.matrices[ , , x]) # Turn 3D array from tapply into a list
names(min.mat.list)<-dimnames(min.matrices)[[3]] # Rename each matrix within the list with its original site name
min.mat.prep <- prepare_networks(min.mat.list, directed = TRUE)  # takes list of networks as matrices, returns a list of igraph objects

beta.networks <-network_betadiversity(min.mat.prep, complete=F)
beta.networks <- melt(beta.networks)
beta.networks$comp <- paste(beta.networks$i,"_",beta.networks$j)
i <- patch[,c(1,2)]
colnames(i)[1] <- "i"
colnames(i)[2] <- "land.use.i"
j <- patch[,c(1,2)]
colnames(j)[1] <- "j"
colnames(j)[2] <- "land.use.j"
beta.networks <- merge(beta.networks,i, by="i")
beta.networks <- merge(beta.networks,j, by="j")
beta.networks$land.use.pair <- paste(beta.networks$land.use.i,"_",beta.networks$land.use.j)

int_turn <- filter(beta.networks, variable == "WN") %>% droplevels()

#remove duplicates
int_turn$land.use.pair[int_turn$land.use.pair == "Cropping _ Avocado"] <- "Avocado _ Cropping"
int_turn$land.use.pair[int_turn$land.use.pair == "Forest _ Avocado"] <- "Avocado _ Forest"
int_turn$land.use.pair[int_turn$land.use.pair == "Dairy _ Avocado"] <- "Avocado _ Dairy"
int_turn$land.use.pair[int_turn$land.use.pair == "Forest _ Cropping"] <- "Cropping _ Forest"
int_turn$land.use.pair[int_turn$land.use.pair == "Dairy _ Cropping"] <- "Cropping _ Dairy"
int_turn$land.use.pair[int_turn$land.use.pair == "Forest _ Dairy"] <- "Dairy _ Forest"

int_turn$land.use.pair <- factor(int_turn$land.use.pair,
                                 levels=c("Dairy _ Dairy",
                                          "Cropping _ Cropping",
                                          "Avocado _ Avocado",
                                          "Forest _ Forest",
                                          "Cropping _ Dairy",
                                          "Avocado _ Dairy",
                                          "Dairy _ Forest",
                                          "Avocado _ Cropping",
                                          "Cropping _ Forest",
                                          "Avocado _ Forest"))

#run model
turn.mod <- glmmTMB(value ~ land.use.pair + (1|comp),
                    family=gaussian,
                    data=int_turn)
plot(residuals(turn.mod, "pearson"))

#calculate difference between site landuse diversity
div.turn <- as.data.frame(for.div[,c(1,2)])
div.turn.dist=dist(div.turn)
div.turn.dist=as.matrix(div.turn.dist, labels=TRUE)
colnames(div.turn.dist) <- rownames(div.turn.dist) <- div.turn[['Site']]
div.turn.dist[upper.tri(div.turn.dist)] <- NA
div.turn.dist <- matrix(div.turn.dist, dimnames=list(t(outer(colnames(div.turn.dist), rownames(div.turn.dist), FUN=paste)), NULL))
div.turn.dist <- setDT(as.data.frame(div.turn.dist), keep.rownames = TRUE)[]
div.turn.dist <- cSplit(div.turn.dist, "rn", sep = " ", direction = "wide")
div.turn.dist <- div.turn.dist[!div.turn.dist$rn_1==div.turn.dist$rn_2,]
div.turn.dist <- na.omit(div.turn.dist)
div.turn.dist$comp <- paste(div.turn.dist$rn_1,"_",div.turn.dist$rn_2)
colnames(div.turn.dist)[1] <- "dif.div"

#calculate difference between site percent forest
for.turn <- as.data.frame(for.div[,c(1,3)])
for.turn.dist=dist(for.turn)
for.turn.dist=as.matrix(for.turn.dist, labels=TRUE)
colnames(for.turn.dist) <- rownames(for.turn.dist) <- div.turn[['Site']]
for.turn.dist[upper.tri(for.turn.dist)] <- NA
for.turn.dist <- matrix(for.turn.dist, dimnames=list(t(outer(colnames(for.turn.dist), rownames(for.turn.dist), FUN=paste)), NULL))
for.turn.dist <- setDT(as.data.frame(for.turn.dist), keep.rownames = TRUE)[]
for.turn.dist <- cSplit(for.turn.dist, "rn", sep = " ", direction = "wide")
for.turn.dist <- for.turn.dist[!for.turn.dist$rn_1==for.turn.dist$rn_2,]
for.turn.dist <- na.omit(for.turn.dist)
for.turn.dist$comp <- paste(for.turn.dist$rn_1,"_",for.turn.dist$rn_2)
colnames(for.turn.dist)[1] <- "dif.for"

#merge dataframes
int_turn <- merge(div.turn.dist,int_turn)
int_turn <- merge(int_turn,for.turn.dist)

turn.mod <- glmmTMB(value ~ land.use.pair + (1|comp),
                        family=gaussian(link="log"),
                        data=int_turn)

turn.mod.all <- glmmTMB(value ~ land.use.pair + dif.for + dif.div + (1|comp),
                    family=gaussian(link="log"),
                    data=int_turn)

dredge(turn.mod.all)

#run lsmeans on landuse model
turn.ls <- emmeans(turn.mod, pairwise ~ land.use.pair, level = .95, adjust = "fdr", int.adjust="none", type = "response")

#generate letters for groups
CLD <- CLD(turn.ls, Letters = letters, level = .95, adjust = "fdr", int.adjust="none")
CLD <- CLD[,c(1,7)]
max <- int_turn %>% group_by(land.use.pair) %>% summarise(max = max(value+0.03))
CLD <- merge(max, CLD)
CLD <- merge(turn.ls$emmeans, CLD)

p <- ggplot()
p <- p + xlab(NULL) + ylab("Plant-pollinator interaction turnover between land-uses")
p <- p + theme(text = element_text(size=18))
p <- p + geom_jitter(data=int_turn, aes(land.use.pair, value),
                     color="grey36", alpha=1, size=1.5, shape=1, position = position_jitter(width = 0.2))
p <- p + geom_point(data=CLD, aes(x=land.use.pair, y=emmean),
                    alpha=1, size=3.5)
p <- p + geom_errorbar(data=CLD, aes(x=land.use.pair, ymin=lower.CL, ymax=upper.CL),
                       width = 0)
p <- p + scale_y_continuous(breaks = seq(0.5, 1, by=0.1), limits=c(0.45,1.04))
p <- p + geom_text(data = CLD, aes(x = land.use.pair, y=max, label=.group), size=5)
p <- p + theme(axis.line.x = element_line(size=0, colour = "black"),
               axis.line.y = element_line(size=0, colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(size=.8, colour = "black", fill=NA),
               panel.background = element_blank()) +
  theme(axis.text.x=element_text(angle= 90, hjust = 1, vjust = 0.5, size =12),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =12),
        axis.title.y=element_text(size=15, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(axis.ticks.length = unit(2.5, "mm"))
p <- p + theme(panel.spacing.x=unit(1, "lines"),panel.spacing.y=unit(1, "lines"))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p <- p + theme(legend.position="none")
p

###############################
#multivariate analysis of site by interactions matrix----
###############################

#cast data by int
network.int.or <- dcast(network.melt, Site + Land.use ~ int, fun.aggregate = sum, na.rm =T, value.var="value", fill = 0)

#give ID to df for ordination
network.int.or$ID <- paste(network.int.or$Site,"_",network.int.or$Land.use)
ord <- network.int.or[,-102]
rownames(ord) <- network.int.or[,102]
ord <- ord[-c(1:2)]

#convert ordnation to binary
#ord[ord>0] <-1

#create vector with method levels (needs to be numeric)
fjT <- c(1,1,1,2,3,3,4,2,2,2,1,2,4,3,2,3,4,3,4,3,3,1,3,1,3,2,4) #method

##################################################
#run the ordination
##################################################

#run vegdist
vj <- vegdist(ord, method = 'bray')
mj <- metaMDS(vj, autotransform = FALSE, wascores = FALSE, k=2)

#make NMDS data frame with method column for plotting
NMDS = data.frame(MDS1 = mj$points[,1], MDS2 = mj$points[,2],site=network.int$Site, landuse=network.int$Land.use)

#create ellipse function (this is hidden in vegan)
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#create df with ellipse values
#drop site from NMDS
df_ell <- data.frame()
for(g in levels(NMDS$landuse)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$landuse==g,],
            veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
            ,group=g))
}

##################################################
#plot the ordination with ellipses
##################################################

p <- ggplot(data = NMDS, aes(MDS1, MDS2))
p <- p + geom_point(aes(colour = landuse), size=2)
p <- p + geom_text(aes(label = site, colour=landuse), vjust = -0.8)
p <- p + geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)
p <- p + theme(axis.line.x = element_line(size=0, colour = "black"),
               axis.line.y = element_line(size=0, colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(size=.8, colour = "black", fill=NA),
               panel.background = element_blank()) +
  theme(axis.text.x=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =12),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =12),
        axis.title.x=element_text(size=15, vjust = 1),
        axis.title.y=element_text(size=15, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(axis.ticks.length = unit(2.5, "mm"))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p <- p + scale_colour_brewer(palette = "Dark2")
p

##################################################
#run the permanova
##################################################

#create method ID dataframe
ID <- select(network.int, one_of(c("Site", "Land.use")))

#Run the permanova
adonis(ord ~ Land.use, data = ID, permutations = 9999, method = "bray")
#is a significant difference in composition of interactions between landuse types

#compute pairwise differences in interaction composition

#first create the function
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='fdr')
{
  library(vegan)
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
}

#run pairwise function on the matrix
pairwise.perm <- pairwise.adonis(network.int.or[,3:101],network.int.or$Land.use)
write.csv(pairwise.perm, "/Users/macuser/Library/Mobile Documents/com~apple~CloudDocs/H_drive_DT/Mareeba_Pollen_Network/model selection tables/pairwise.perm.csv")

#have a look at dispersion
## Calculate multivariate dispersions
mod <- betadisper(vj, ID$Land.use)
mod

## Perform test
anova(mod)

## Permutation test for F
permutest(mod, pairwise = TRUE)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))#all good!

#######################################
#calculate specialisation----
#######################################

web <- network.int.or[, -c(1,2,102)]
web <- web[,colSums(web) > 19]#remove species with no links at each site
obs <- unlist(specieslevel(web, index="PDI", level="higher"))
obs <- t(as.data.frame(obs))
nm <- nullmodel(web = web, N=999, method=3)#using Vazquez method - check this
null <- sapply(nm, specieslevel, index="PDI", level="higher")
null2 <- rbind.fill(lapply(null, as.data.frame))
null3 <- as.data.frame(matrix(null2[,1], nrow=999, ncol=19))

p1 <- sum(null3[,1]>obs[,1]) / length(null3[,1])
p2 <- sum(null3[,2]>obs[,2]) / length(null3[,2])
p3 <- sum(null3[,3]>obs[,3]) / length(null3[,3])
p4 <- sum(null3[,4]>obs[,4]) / length(null3[,4])
p5 <- sum(null3[,5]>obs[,5]) / length(null3[,5])
p6 <- sum(null3[,6]>obs[,6]) / length(null3[,6])
p7 <- sum(null3[,7]>obs[,7]) / length(null3[,7])
p8 <- sum(null3[,8]>obs[,8]) / length(null3[,8])
p9 <- sum(null3[,9]>obs[,9]) / length(null3[,9])
p10 <- sum(null3[,10]>obs[,10]) / length(null3[,10])
p11 <- sum(null3[,11]>obs[,11]) / length(null3[,11])
p11 <- ifelse(p11 > 0.5, 1-p11, p11)
p12 <- sum(null3[,12]>obs[,12]) / length(null3[,12])
p13 <- sum(null3[,13]>obs[,13]) / length(null3[,13])
p14 <- sum(null3[,14]>obs[,14]) / length(null3[,14])
p15 <- sum(null3[,15]>obs[,15]) / length(null3[,15])
p16 <- sum(null3[,16]>obs[,16]) / length(null3[,16])
p17 <- sum(null3[,17]>obs[,17]) / length(null3[,17])
p18 <- sum(null3[,18]>obs[,18]) / length(null3[,18])
p18 <- ifelse(p18 > 0.5, 1-p18, p18)
p19 <- sum(null3[,19]>obs[,19]) / length(null3[,19])

p.all <- as.data.frame(rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19))
colnames(p.all) <- "pvalue"
obs2 <- t(obs)
PDI <- cbind(obs2,p.all)
rownames(PDI) <- colnames(web)
colnames(PDI)[1] <- "PDI_value" #PDI  0 (perfect generalist) and 1 (perfect specialist)

#add column with which landuse type each interaction is specialised to
network.int.or$ID <- paste(network.int.or$Site, network.int.or$Land.use, sep = "_")
sub.int <- network.int.or[, -c(1,2,102)]
rownames(sub.int) <- network.int.or[,102]
sub.int <- sub.int[,colSums(sub.int) > 19]#remove species with no links at each site
setDT(sub.int, keep.rownames = TRUE)[]
mx <- melt(sub.int)
mx <- mx %>%
  group_by(variable) %>%
  filter(value == max(value))
mx$rn <- mx$rn %>% str_replace(".*_", "")
setDT(PDI, keep.rownames = TRUE)[]
colnames(PDI)[1] <- "variable"
#merge with origional specialisation df
PDI <- merge(PDI, mx, by="variable")
colnames(PDI)[1] <- "Interaction"
colnames(PDI)[4] <- "Land-use"
PDI <- PDI[,-5]

write.csv(PDI, "/Users/macuser/Library/Mobile Documents/com~apple~CloudDocs/H_drive_DT/Mareeba_Pollen_Network/model selection tables/specialisation.csv")

##################################################
#END
##################################################