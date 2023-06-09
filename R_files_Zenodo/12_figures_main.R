# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Simonet & McNally 2020           #
#   Figures produced from  models output   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#local_project_dir='/path/to/where/repo/is/cloned'
setwd(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/'))
source('./scripts/analysis/sourced_ggthemes.R')
source('./scripts/analysis/sourced_packages.R')
library("ape")

# GET DATA ----

d<- read.table('./output/tables/ANALYSIS_DATA_ASSEMBLED.txt', header=TRUE, stringsAsFactors = FALSE) %>% mutate(first = !duplicated(species_id)) %>% rename(species = species_id)

d_mean<- d %>%
  group_by(species) %>%
  mutate(mean_relative_abundance = mean(within_host_relative_abundance)) %>%
  select(-host, -within_host_relatedness, -within_host_relative_abundance) %>%
  subset(first == TRUE) %>%
  ungroup() %>%
  as.data.frame()

midas.tree<- read.tree('./data/species_info_files/midas_tree_renamed.newick')
phylogeny<- drop.tip(midas.tree, midas.tree$tip.label[which(!is.element(midas.tree$tip.label, d_mean$species))])
phylogeny<- chronopl(phylogeny, lambda = 0)
phylogeny<-makeNodeLabel(phylogeny)

# Made shorter species names for plots
plot_ids<- read.table("data/species_info_files/species_plot_names.txt", header=TRUE, sep ="\\t")

d_mean<- left_join(d_mean, plot_ids, by = 'species')


# FIGURE 1 ----

# Make figure 1C
# Make species as levels ordered by mean relatedness
#d_mean$species<- factor(d_mean$species, d_mean[order(d_mean$mean_relatedness),'species'])


d_mean$plot_names<- factor(d_mean$plot_names, d_mean[order(d_mean$mean_relatedness),'plot_names'])
d_mean$species<- factor(d_mean$species, d_mean[order(d_mean$mean_relatedness),'species'])

d_mean<- d_mean %>% arrange(mean_relatedness)


ylim = c(1, length(unique(d$species)))
xlim = c(-1, 1)

# Output figure
ppi <- 1200# 800

png('./output/figures/FIG_1C.png', width=(8/2.54)*ppi, height=(17.8/2.54)*ppi, res=ppi)

# Set the plot
par(mar = c(2, 3.8, 0, 1.5))
plot(y = d_mean$plot_id, x = d_mean$mean_relatedness, pch = 16, xlim = xlim, ylim = ylim, cex = 0.1,
     yaxt='n', ylab = '', xlab = '', bty = 'n', xaxt = 'n')

# Background shades and within host relatedness
for(i in 1:nrow(d_mean)){

# rectangle
rect(xleft = xlim[1], ybottom = i-0.5, xright = xlim[2], ytop = i+0.5,
     border = NA,
     col = ifelse(i %% 2 == 0, 'gray98', 'gray93'))
  
  
foo<- d[d$species == levels(d_mean$species)[i],] 

for(h in 1:nrow(foo)){
  lines(y = c(i+0.5,i-0.5), x = rep(foo$within_host_relatedness[h], 2),
        col = ifelse(foo$within_host_relatedness[h] >= quantile(foo$within_host_relatedness, 0.25) & 
                     foo$within_host_relatedness[h] <= quantile(foo$within_host_relatedness, 0.75),
                     'dodgerblue', 'gray66'),
        lwd = 0.5)
}

}

# Mean relatedness
points(y = d_mean$species, x = d_mean$mean_relatedness, pch = 16, xlim = c(-1.5, 2), cex = 0.6)

# Axes labels and ticks
axis(2, at = seq(1,101), labels = rep('', 101), las = 2, cex.axis = 0.4, font = 4, tck = -.02, line = -0.1)
axis(2, at = seq(1,101), labels = d_mean$plot_name, las = 2, cex.axis = 0.4, font = 4, line = -0.8, lwd = 0)
axis(4, at = seq(1,101), labels = paste0('n = ', d_mean$nb_host), las = 2, cex.axis = 0.4, padj = 0, lwd = 0, line = -1)
axis(1, line = -0.9, lwd = 1.1, cex.axis = 0.6, tck = -.04, padj = -2)
mtext('Mean relatedness', 1, font = 2, line = 0.3, cex = 0.8)

dev.off()


# Load images and assemble panel 1
library(magick)

fig1A <- image_read("./output/figures/Fig_1A.png")
fig1B <- image_read("./output/figures/FIG_1B.png")
fig1C <- image_read("./output/figures/FIG_1C.png")


width_1A<- round(((8.8)*100)/17.8)
width_1C<- 100 - width_1A
width_1B<- round(((8.8)*100)/17.8)
height_1A<- round(((5.5)*100)/17.8)
height_1B<-  100 - height_1A

m<- matrix(data =  c(
rep(c(rep(1, width_1A), rep(3, width_1C)), height_1A),
rep(c(rep(2, width_1B), rep(3, width_1C)), height_1B)),
nrow = 100,
byrow = TRUE)


pdf(file = paste0("./output/figures/FIGURE_1.pdf"), width=(17.8/2.54), height=(17.8/2.54))

layout(m)#;layout.show(3)

par(mar = c(1, 1.5, 1, 1.5)) # bottom, left, top, right
plot(NA, xlim=0:1, ylim=0:1, bty="n", axes=0, xaxs = 'i', yaxs='i')
rasterImage(fig1A, 0, 0, 1, 0.95)
text(x = 0.02, y = 0.97, 'A', font = 2, cex = 1.5) # FIG 1A

par(mar = c(0, 1.5, 1, 1.5))
plot(NA, xlim=0:1, ylim=0:1, bty="n", axes=0, xaxs = 'i', yaxs='i')
rasterImage(fig1B, 0.01, 0.1, 1, 1)
text(x = 0.02, y = 0.97, 'B', font = 2, cex = 1.5) # FIG 1B

par(mar = c(0, 1, 0, 0))
plot(NA, xlim=0:1, ylim=0:1, bty="n", axes=0, xaxs = 'i', yaxs='i')
rasterImage(fig1C, 0, 0, 0.96,1)
text(x = 0.03, y = 0.985, 'C', font = 2, cex = 1.5) # FIG 1C

dev.off()


# FIGURE 2 ----

library(ggtree)
library(phytools)
plot_character_2<- function(s, d, l, adj, col, nbtips, trans){
  a<- ((2*pi)/nbtips)
  alpha = s*a
  eps = (d/(tan(asin(d))))-1
  
  if(i <= nbtips/4){
    xt = cos(alpha)
    yt = sin(alpha)
    
    x2 = cos(asin(d) + alpha) * (1 + eps)
    y2 = sin(asin(d) + alpha) * (1 + eps)
    
    x3 = cos(alpha) * (1 + l)
    y3 = sin(alpha) * (1 + l)
    
    x4 = -xt + x2 + x3
    y4 = y2 - yt + y3
    
    x5 = xt + abs(xt - x2)
    y5 = yt - abs(y2 - yt)
    
    x6 = x3 + abs(x5 - xt)
    y6 = y3 - abs(y5 - yt)
    
    
    tx = abs(abs(x4) - abs(x2))
    ty = abs(abs(y4) - abs(y2))
    
    x2 = x2+adj
    x4 = x4+adj
    x6 = x6+adj
    x5 = x5+adj
    
    y2 = y2+adj
    y4 = y4+adj
    y6 = y6+adj
    y5 = y5+adj
    
    
    polygon(x = c(x2 + trans*tx,
                  x4 + trans*tx,
                  x6 + trans*tx,
                  x5 + trans*tx),
            
            y = c( y2 + trans*ty,
                   y4 + trans*ty,
                   y6 + trans*ty,
                   y5 + trans*ty),
            
            col = col, border = NA)
    
  }else if (i >= nbtips/4 && i < nbtips/2){
    xt = cos(alpha)
    yt = sin(alpha)
    
    x2 = cos(asin(d) + alpha) * (1 + eps)
    y2 = sin(asin(d) + alpha) * (1 + eps)
    
    x3 = cos(alpha) * (1 + l)
    y3 = sin(alpha) * (1 + l)
    
    x4 =  x3-abs(abs(xt) - abs(x2))
    y4 =  y3-abs(abs(yt) - abs(y2))
    
    x5 = xt + abs(abs(xt) - abs(x2))
    y5 = yt + abs(abs(y2) - abs(yt))
    
    x6 = x3 + abs(abs(x5) - abs(xt))
    y6 = y3 + abs(abs(y5) - abs(yt))
    
    
    tx = abs(abs(x6) - abs(x5))
    ty = abs(abs(y6) - abs(y5))
    
    x2 = x2-adj
    x4 = x4-adj
    x6 = x6-adj
    x5 = x5-adj
    
    y2 = y2+adj
    y4 = y4+adj
    y6 = y6+adj
    y5 = y5+adj
    
    
    polygon(x = c(x2 - trans*tx,
                  x4 - trans*tx,
                  x6 - trans*tx,
                  x5 - trans*tx),
            
            y = c( y2 + trans*ty,
                   y4 + trans*ty,
                   y6 + trans*ty,
                   y5 + trans*ty),
            
            col = col, border = NA)
    
    
  }else if (i >= nbtips/2 && i < nbtips*0.75){
    xt = cos(alpha)
    yt = sin(alpha)
    
    x2 = cos(asin(d) + alpha) * (1 + eps)
    y2 = sin(asin(d) + alpha) * (1 + eps)
    
    x3 = cos(alpha) * (1 + l)
    y3 = sin(alpha) * (1 + l)
    
    x4 =  x3+abs(abs(xt) - abs(x2))
    y4 =  y3-abs(abs(yt) - abs(y2))
    
    x5 = xt - abs(abs(xt) - abs(x2))
    y5 = yt + abs(abs(y2) - abs(yt))
    
    x6 = x3 - abs(abs(x5) - abs(xt))
    y6 = y3 + abs(abs(y5) - abs(yt))
    
    tx = abs(abs(x6) - abs(x5))
    ty = abs(abs(y6) - abs(y5))
    
    x2 = x2-adj
    x4 = x4-adj
    x6 = x6-adj
    x5 = x5-adj
    
    y2 = y2-adj
    y4 = y4-adj
    y6 = y6-adj
    y5 = y5-adj
    
    
    polygon(x = c(x2 - trans*tx,
                  x4 - trans*tx,
                  x6 - trans*tx,
                  x5 - trans*tx),
            
            y = c( y2 - trans*ty,
                   y4 - trans*ty,
                   y6 - trans*ty,
                   y5 - trans*ty),
            
            col = col, border = NA)
    
    
  }else if (i >= nbtips*0.75 && i <= nbtips){
    xt = cos(alpha)
    yt = sin(alpha)
    
    x2 = cos(asin(d) + alpha) * (1 + eps)
    y2 = sin(asin(d) + alpha) * (1 + eps)
    
    x3 = cos(alpha) * (1 + l)
    y3 = sin(alpha) * (1 + l)
    
    x4 =  x3+abs(abs(xt) - abs(x2))
    y4 =  y3+abs(abs(yt) - abs(y2))
    
    x5 = xt - abs(abs(xt) - abs(x2))
    y5 = yt - abs(abs(y2) - abs(yt))
    
    x6 = x3 - abs(abs(x5) - abs(xt))
    y6 = y3 - abs(abs(y5) - abs(yt))
    
    tx = abs(abs(x4) - abs(x2))
    ty = abs(abs(y4) - abs(y2))
    
    x2 = x2+adj
    x4 = x4+adj
    x6 = x6+adj
    x5 = x5+adj
    
    y2 = y2-adj
    y4 = y4-adj
    y6 = y6-adj
    y5 = y5-adj
    
    
    polygon(x = c(x2 + trans*tx,
                  x4 + trans*tx,
                  x6 + trans*tx,
                  x5 + trans*tx),
            
            y = c( y2 - trans*ty,
                   y4 - trans*ty,
                   y6 - trans*ty,
                   y5 - trans*ty),
            
            col = col, border = NA)
    
  }else{
    print('It seems you are attempting to plot more species than there is in your tree')
  }
  
  
  
  
}


# get the order of species in the same as tree from ggtree
g<- ggtree(phylogeny)
spl<- as.data.frame(g$data)
spl<- spl[spl$isTip==TRUE,]


# build a color map from it
names(spl)<- c(names(spl)[1:3], 'species', names(spl)[5:9])


wideall<- d %>% filter(first == TRUE)
spl2<- left_join(spl, wideall, by = 'species')


rvals<- spl2[,c('species', 'mean_relatedness')]
rvals<-setNames(rvals[,2],rvals[,1])
obj<-contMap(phylogeny,rvals,plot=FALSE)
mycols<- rep('black', 1001)
names(mycols)<- 0:1000
obj$cols[]<-mycols
obj<-setMap(obj,invert=TRUE)

colmin = 'snow'
colmax = 'navy'
bias = 0.01
trait = 'secretion_system_no4'

map2color2<- function(colmin, colmax, bias, trait){
  colortest3<- colorRampPalette(c(colmin, colmax), bias = bias)(length(unique(spl2[,which(colnames(spl2)==trait)][!is.na(spl2[,which(colnames(spl2)==trait)])])))
  test2<- spl2[,c(which(colnames(spl2)== 'species'), which(colnames(spl2)== trait))]
  
  mapping<- data.frame(values = sort(unique(test2[,2])), colors = colortest3)
  colnames(mapping)<- c(trait, 'colors')
  
  
  yo<- data.frame(value = NA, colors = 'black')
  colnames(yo)<- c(trait, 'colors')
  mapping<- rbind(mapping, yo)
  
  
  test2<- left_join(test2, mapping, trait)
  test2$colors<- as.character(test2$colors)
  
  return(list(test2$colors, mapping))
}

my.arc.cladelabel<- function(plotted_obj, phylogeny, node, arc.offset, lab.offset, text_color, arc.lwd, marknode, clade_label, cex, text.offset){
  
  g<- ggtree(phylogeny)
  gdata<- as.data.frame(g$data)
  tree = plotted_obj$tree
  
  objtree <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  
  # Add the arc
  arc.cladelabels(tree=obj$tree,' ', node, ln.offset = arc.offset, lwd = 2, mark.node=marknode)
  
  # Add name of clade
  
  h <- max(sqrt(objtree$xx^2 + objtree$yy^2))
  d <- getDescendants(tree, node)
  d <- sort(d[d <= Ntip(tree)])
  deg <- atan(objtree$yy[d]/objtree$xx[d]) * 180/pi
  ii <- intersect(which(objtree$yy[d] >= 0), which(objtree$xx[d] < 0))
  deg[ii] <- 180 + deg[ii]
  ii <- intersect(which(objtree$yy[d] < 0), which(objtree$xx[d] < 0))
  deg[ii] <- 180 + deg[ii]
  ii <- intersect(which(objtree$yy[d] < 0), which(objtree$xx[d] >=0))
  deg[ii] <- 360 + deg[ii]
  
  x0 <- lab.offset * cos(median(deg) * pi/180) * h
  y0 <- lab.offset * sin(median(deg) * pi/180) * h
  
  #text.offset=1+(nchar(clade_label)*0.015)
  
  text(x = text.offset*x0, y = text.offset*y0, clade_label, col = text_color, srt = atan(y0/x0)*(180/pi), cex = cex, font = 3, adj = ifelse(x0 < 0, 1, 0))
  
}

legendbar<- function(colormap, x, y, ytop, digits, label, cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust){
  
  # Continuous legend colorbar for relatedness
  colmap<- colormap[[2]][-(nrow(colormap[[2]])),]
  colmap$colors<- as.character(colmap$colors)
  # clbr2<-matrix(ncol=nrow(colmap),nrow=2)
  # clbr2[1,]<-seq(0,1,length.out=nrow(colmap))
  # clbr2[2,]<-seq(0,1,length.out=nrow(colmap))
  
  xb2 = x
  yb2 = y
  yb2top = ytop
  
  
  # image(z=clbr2,y=seq(yb2, yb2+yb2top,length.out=nrow(colmap)),x=c(xb2,xb2+0.05),col=colmap$colors,zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
  
  
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj), c(yb2, yb2+yb2top), lwd = 1.2)
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj+tick_length), c(yb2,yb2), lwd = 1.2)
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj+tick_length), c(yb2+0.25*yb2top,yb2+0.25*yb2top), lwd = 1.2)
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj+tick_length), c(yb2+0.5*yb2top,yb2+0.5*yb2top), lwd = 1.2)
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj+tick_length), c(yb2+0.75*yb2top,yb2+0.75*yb2top), lwd = 1.2)
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj+tick_length), c(yb2+yb2top,yb2+yb2top), lwd = 1.2)
  
  
  labs<- c(round(min(colmap[,1]),digits),
           round(quantile(colmap[,1], 0.25),digits),
           round(quantile(colmap[,1], 0.5),digits),
           round(quantile(colmap[,1], 0.75),digits),
           round(max(colmap[,1]),digits))
  
  text(y = c(yb2, yb2+0.25*yb2top, yb2+0.5*yb2top, yb2+0.75*yb2top, yb2+yb2top), x = xb2+ticks_labels_adjust, pos = 4, labels = labs, font = 1, cex = cex_ticks)
  
  
  text(y = yb2+0.5*yb2top, x = xb2-label_adjust, labels = label, font = font_label, cex = cex_label, srt = 90)
  
  
}

coltoplot1<- map2color2('snow', 'darkred', 0.01, 'mean_relatedness')
coltoplot2<- map2color2('snow', 'darkorchid4', 0.01, 'nb_extracellular')
coltoplot3<- map2color2('snow', 'navy', 0.01, 'secretion_system_no4')
coltoplot4<- map2color2('snow', 'darkorange', 0.01, 'siderophores')
coltoplot5<- map2color2('snow', 'yellow', 0.01, 'quorum_sensing')
coltoplot6<- map2color2('snow', 'forestgreen', 0.01, 'biofilm')
coltoplot7<- map2color2('snow', 'magenta', 0.01, 'ab_degradation')


spl2_mock<- spl2[c(1,2),]
spl2_mock[c(1,2),]<- NA
spl2_mock[,'nb_extracellular']<- c(0,1)
spl2_mock[,'mean_relatedness']<- c(0,1)

spl2<- rbind(spl2, spl2_mock)

coltoplot1<- map2color2('snow', 'darkred', 0.01, 'mean_relatedness')
coltoplot2<- map2color2('snow', 'darkorchid4', 0.01, 'nb_extracellular')

# Fig 2: tree and clades

#pdf(file = './output/figures/Fig2_tree_and_clades2.pdf', width = (11.3/2.54), height = (11.3/2.54))

ppi=800
png('./output/figures/FIG_2.tree.png', width=(11.3/2.54)*ppi, height=(11.3/2.54)*ppi, res=ppi)

par(mar = c(0,0,0,0),
    oma = c(0,0,0,0))

plot(obj,type="fan",ftype="off",lwd=c(2,6),outline=FALSE,
     xlim=c(-1,1),
     ylim = c(-2.5,1.8), legend=FALSE)


for(i in 1:101){
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot1[[1]][i] , 101, 0)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot2[[1]][i] , 101, 1.2)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot3[[1]][i] , 101, 2.4)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot4[[1]][i] , 101, 3.6)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot5[[1]][i] , 101, 4.8)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot6[[1]][i] , 101, 6.0)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot7[[1]][i] , 101, 7.5)
  
}

plotted_obj = obj
#phylogeny = mytree
arc.offset = 1.3
lab.offset = 1.31
text_color = 'black'
arc.lwd = 2
marknode = FALSE
cex = 0.4


nodes_and_genus<- read.table('./data/species_info_files/nodes_genus.txt', colClasses = c('numeric', 'character','numeric'),sep = '\\t', header=TRUE)


nodes_and_genus$text.offset<- 1.02#+(nchar(nodes_and_genus$genus)*0.015)


for(i in 1:nrow(nodes_and_genus)){
  
  my.arc.cladelabel(plotted_obj,
                    phylogeny,
                    node = nodes_and_genus$node[i],
                    arc.offset,
                    lab.offset, 'black',
                    arc.lwd, marknode,
                    clade_label = nodes_and_genus$genus[i],
                    cex, text.offset = nodes_and_genus$text.offset[i])
}

dev.off()



# Fig 2: legends
colors_legend<- data.frame(colors_1 = colorRampPalette(c('snow', 'darkred'), bias = 0.01)(101),
                           colors_2 = colorRampPalette(c('snow', 'darkorchid4'), bias = 0.01)(101),
                           colors_3 = colorRampPalette(c('snow', 'navy'), bias = 0.01)(101),
                           colors_4 = colorRampPalette(c('snow', 'darkorange'), bias = 0.01)(101),
                           colors_5 = colorRampPalette(c('snow', 'yellow'), bias = 0.01)(101),
                           colors_6 = colorRampPalette(c('snow', 'forestgreen'), bias = 0.01)(101),
                           colors_7 = colorRampPalette(c('snow', 'magenta'), bias = 0.01)(101)
)


y = -2.5
ytop = 0.65

colors_legend$colors_1<- as.character(colors_legend$colors_1)
colors_legend$colors_2<- as.character(colors_legend$colors_2)
colors_legend$colors_3<- as.character(colors_legend$colors_3)
colors_legend$colors_4<- as.character(colors_legend$colors_4)
colors_legend$colors_5<- as.character(colors_legend$colors_5)
colors_legend$colors_6<- as.character(colors_legend$colors_6)
colors_legend$colors_7<- as.character(colors_legend$colors_7)

clbr2<-matrix(ncol=nrow(colors_legend),nrow=2)
clbr2[1,]<-seq(0,1,length.out=nrow(colors_legend))
clbr2[2,]<-seq(0,1,length.out=nrow(colors_legend))


x = 0.5
y = 0
ytop = 1
digits = 1
cex_ticks = 0.45
cex_label = 0.5
vert_bar_adj = 0.15
tick_length = 0.03
label_adjust = 0.43
ticks_labels_adjust = -0.03
font_label = 2

png('./output/figures/FIG_2.bar1.png', width=(1/2.54)*ppi, height=(2.5/2.54)*ppi, res=ppi)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_1, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot1, x, y, ytop, digits, label = 'Relatedness', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()


digits = 0

png('./output/figures/FIG_2.bar2.png', width=(1/2.54)*ppi, height=(3/2.54)*ppi, res=ppi)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_2, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot2, x, y, ytop, digits, label = 'Secretome size', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()


png('./output/figures/FIG_2.bar3.png', width=(1/2.54)*ppi, height=(3/2.54)*ppi, res=ppi)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_3, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot3, x, y, ytop, digits, label = 'Secretion systems', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()


png('./output/figures/FIG_2.bar4.png', width=(1/2.54)*ppi, height=(3/2.54)*ppi, res=ppi)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_4, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot4, x, y, ytop, digits, label = 'Siderophores', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()


png('./output/figures/FIG_2.bar5.png', width=(1/2.54)*ppi, height=(3/2.54)*ppi, res=ppi)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_5, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot5, x, y, ytop, digits, label = 'Quorum sensing', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()



png('./output/figures/FIG_2.bar6.png', width=(1/2.54)*ppi, height=(3/2.54)*ppi, res=ppi)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_6, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot6, x, y, ytop, digits, label = 'Biofilm', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()

png('./output/figures/FIG_2.bar7.png', width=(1/2.54)*ppi, height=(3/2.54)*ppi, res=ppi)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_7, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot7, x, y, ytop, digits, label = 'Antibiotic degradation', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()


# Figure 2: ASSEMBLE

fig2 <- image_read("./output/figures/FIG_2.tree.png")

bar1 <- image_read("./output/figures/FIG_2.bar1.png")
bar2 <- image_read("./output/figures/FIG_2.bar2.png")
bar3 <- image_read("./output/figures/FIG_2.bar3.png")
bar4 <- image_read("./output/figures/FIG_2.bar4.png")
bar5 <- image_read("./output/figures/FIG_2.bar5.png")
bar6 <- image_read("./output/figures/FIG_2.bar6.png")
bar7 <- image_read("./output/figures/FIG_2.bar7.png")


height_bar<- round((3*100)/11.4)
height_tree<-  round(100 - height_bar)
width_bar<- round(100/7)

fig2_crop <- image_crop(fig2, "2900x2900+350+150")
plot(NA, xlim=0:1, ylim=0:1, bty="n", axes=0, xaxs = 'i', yaxs='i')
rasterImage(fig2_crop, 0.15, 0, 0.85,1)


m<- matrix(
  data = c(rep(rep(1,7), 5),
           rep(rep(2,7), 2)),
  byrow = TRUE,
  nrow = 7
)



pdf(file = paste0("./output/figures/FIGURE_2.pdf"), width=(11.4/2.54), height=(11.4/2.54))
layout(m)#; layout.show(2)

par(mar = c(0, 0, 0, 0))
plot(NA, xlim=0:1, ylim=0:1, bty="n", axes=0, xaxs = 'i', yaxs='i')
rasterImage(fig2_crop, 0.15, 0, 0.85,1)

par(mar = c(0, 0, 0, 0))
plot(NA, xlim=0:1, ylim=0:1, bty="n", axes=0, xaxs = 'i', yaxs='i')

eps = 0.02
rasterImage(bar1, 0.15+eps, 0.25, 0.23+eps, 1)
rasterImage(bar2, 0.25+eps, 0.25, 0.33+eps, 1)
rasterImage(bar3, 0.35+eps, 0.25, 0.43+eps, 1)
rasterImage(bar4, 0.45+eps, 0.25, 0.53+eps, 1)
rasterImage(bar5, 0.55+eps, 0.25, 0.63+eps, 1)
rasterImage(bar6, 0.65+eps, 0.25, 0.73+eps, 1)
rasterImage(bar7, 0.75+eps, 0.25, 0.83+eps, 1)


dev.off()



 # FIGURE 3 ----

load("./output/analyses/MODEL1_CHAIN_1.RData")

# Code wrapper to do one plot per GO cooperative trait
plotLM<- function(response, D, col, length.yseq, cex.axis = 0.6, cex = 0.6, tck = -0.05, line.x = -1, line.y = -0.7){


mod<-lm(I(D[,response]/total_cds) ~ mean_relatedness,data=D)

yseq<- seq(
  range(D[,response]/D$total_cds, na.rm = TRUE)[1],
  range(D[,response]/D$total_cds, na.rm = TRUE)[2], length.out = length.yseq)

yseq.labels<- formatC(yseq*1000, digits = 0, format = 'f')

#yseq.labels[yseq == 0]<- '0'


plot(I(D[,response]/total_cds)~ mean_relatedness,data=D,
     col=t_col(col),
     pch=16, bty="l", xlab="", ylab="", cex=cex, las = 1,
     xaxt='n', yaxt='n',
     xlim = c(0,1),
     ylim=range(yseq))

abline(coef(mod)[1],coef(mod)[2],lwd=1)

axis(1, cex.axis = cex.axis, font = 1, tck = tck, at = c(0, 0.25, 0.5, 0.75, 1),  labels = rep('', 5))
axis(1, cex.axis = cex.axis, font = 1, tck = tck, lwd = 0, line = line.x, at = c(0, 0.25, 0.5, 0.75, 1), labels = c('0',' ',"0.5",' ', "1"))
axis(2, cex.axis = cex.axis, font = 1, tck = tck, at = yseq, labels = rep('', length(yseq)))
axis(2, cex.axis = cex.axis, font = 1, tck = tck, lwd = 0, line = line.y, at = yseq, labels = yseq.labels
  )


}



# Append meta-analysisoutput to model results output to get a common table used for plotting
ma.1.df<- data.frame(cooperative_trait = 'Meta-analysis',
                     effect = as.numeric(MA.MODELS_1[which(MA.MODELS_1[,1] == 'mean_relatedness'),'estimate']),
                     hpd_lower = as.numeric(MA.MODELS_1[which(MA.MODELS_1[,1] == 'mean_relatedness'),'ci.lower']),
                     hpd_higher = as.numeric(MA.MODELS_1[which(MA.MODELS_1[,1] == 'mean_relatedness'),'ci.upper']))


dat.R <- dat.R %>% arrange(effect)


dat.R.relatedness<- 
  rbind(dat.R[dat.R$predictor_variable == 'mean_relatedness',c('cooperative_trait', 'effect', 'hpd_lower', 'hpd_higher')],
        ma.1.df)


pdf('./output/figures/FIGURE_3.pdf', width=(11.3/2.54), height=(6.5/2.54))

m<-  rbind(c(1,1,2,3), c(1,1,4,5),c(1,1,6,7))
layout(m)


par(mar = c(1, 4, 0.15, 1),
    oma=c(1,0,1,0))
plot(NA, xlim=c(-2.5,3), ylim=c(1,7), bty = 'L', xlab = 'Estimated effect\\n(phylogenetic mixed model posterior mean)', ylab = '', xaxt = 'n', yaxt = 'n')
axis(1, at = seq(-2,3,1), labels = rep('', 6),  font = 1, cex.axis = 0.6,  tck = -.02)
axis(1, at = seq(-2,3,1), font = 1, cex.axis = 0.6, tck = -.02, lwd = 0, line = -1)

axis(2, at = seq(1,7,1), labels = rep('', 7), tck = -.02)
axis(2, at = seq(1,7,1),
     labels = c('Quorum\\nsensing', 'Secretion\\nsystems', 'Secretome', 'Antibiotic\\n degradation', 'Biofilm', 'Siderophores', 'Overall\\neffect'), font = 2,
     las = 2, lwd = 0, line = -0.5,
     cex.axis = 0.6)

cols = c('yellow', 'navy', 'darkorchid4', 'magenta', 'forestgreen', 'darkorange', 'black')

# to add pvalues
p.secretome = round(summary(mods.R$secretome)$solutions['mean_relatedness','pMCMC'], 3)
p.biofilm = round(summary(mods.R$biofilm)$solutions['mean_relatedness','pMCMC'], 3)
p.sid = summary(mods.R$siderophores)$solutions['mean_relatedness','pMCMC']
p.qs = round(summary(mods.R$quorum_sensing)$solutions['mean_relatedness','pMCMC'], 3)
p.ab = round(summary(mods.R$ab_degradation)$solutions['mean_relatedness','pMCMC'], 3)
p.ssyst = round(summary(mods.R$secretion_system_no4)$solutions['mean_relatedness','pMCMC'], 3)
p.overall = format(as.numeric(MA.MODELS_1[1,'p.value']), digits = 2)


ps<- c(round(p.qs, 2), round(p.ssyst, 2), format(p.secretome, scientific = T, digits = 2), round(p.ab, 2), format(p.biofilm, scientific = T, digits = 2), format(p.sid, scientific = T, digits = 2), p.overall)


for(i in 1:nrow(dat.R.relatedness)){
  points(x = dat.R.relatedness$effect[i], y = i, pch = 16, cex = 0.8, col = cols[i]) 
  lines(x = dat.R.relatedness[i,c('hpd_lower', 'hpd_higher')], y = rep(i,2), col =  cols[i])
  lines(x = rep(dat.R.relatedness[i,c('hpd_lower')], 2), y = c(i+0.05, i-0.05), col =  cols[i])
  lines(x = rep(dat.R.relatedness[i,c('hpd_higher')], 2), y = c(i+0.05, i-0.05), col =  cols[i])

  
  ptext = ifelse(i<7, 'pval: ', 'pval: ')
  
  #text(paste0('(', ptext, ps[i], ')'), x = dat.R.relatedness[i,c('hpd_higher')], y  = i+0.15, cex = 0.4)
  
}


abline(v = 0, lty = 'dotted')
mtext(side = 3, 'A', font = 2, at = -4.5, cex = 0.5)



par(mar = c(1, 2.5, 0.15, 0.5))


plotLM('siderophores', D = d_mean, col = 'darkorange', length.yseq = 3)
mtext(side = 3, 'B', font = 2, at = -0.35, cex = 0.5)
plotLM('biofilm', D = d_mean, col = 'forestgreen', length.yseq = 3)
plotLM('ab_degradation', D = d_mean, col = 'magenta', length.yseq = 3)

mtext(side = 2, text = expression(bold(paste("Proportion of cooperative genes (10"^" -03", ")"))), cex = 0.4, line = 1.2)

D<-d_mean
mod<-lm(I(nb_extracellular/total_cds)~ mean_relatedness+gram_profile,data=D[D$gram_profile != 'gram0',])


yseq<- seq(range(D$nb_extracellular/D$total_cds, na.rm = TRUE)[1],
           range(D$nb_extracellular/D$total_cds, na.rm = TRUE)[2], length.out = 3)

yseq.labels<- formatC(yseq*1000, digits = 0, format = 'f')


plot(I(nb_extracellular[gram_profile=="n"]/total_cds[gram_profile=="n"])~ mean_relatedness[gram_profile=="n"],
     data=D,
     col=t_col('darkorchid4'),
     ylim=range(yseq), xlim = c(0,1), bty="l", xlab="", ylab="", cex = 0.6, las = 1,
     xaxt = 'n', yaxt='n')


points(I(nb_extracellular[gram_profile=="p"]/total_cds[gram_profile=="p"])~ mean_relatedness[gram_profile=="p"],
       data=D,
       col=t_col('darkorchid4'), pch=16, cex=0.6)

abline(coef(mod)[1],coef(mod)[2],lwd=1,lty=2,col=rgb(r=0,g=0,b=0))
abline(coef(mod)[1]+coef(mod)[3],coef(mod)[2],lwd=1,col=rgb(r=0,g=0,b=0))

cex.axis = 0.6
cex = 0.6
tck = -0.05
line.x = -1
line.y = -0.7

axis(1, cex.axis = cex.axis, font = 1, tck = tck, at = c(0, 0.25, 0.5, 0.75, 1),  labels = rep('', 5))
axis(1, cex.axis = cex.axis, font = 1, tck = tck, lwd = 0, line = line.x, at = c(0, 0.25, 0.5, 0.75, 1), labels = c('0',' ',"0.5",' ', "1"))
axis(2, cex.axis = cex.axis, font = 1, tck = tck, at = yseq, labels = rep('', length(yseq)))
axis(2, cex.axis = cex.axis, font = 1, tck = tck, lwd = 0, line = line.y, at = yseq, labels = yseq.labels
)


plotLM('secretion_system_no4', D = d_mean, col = 'navy', length.yseq = 3)
plotLM('quorum_sensing', D = d_mean, col = 'yellow', length.yseq = 3)

mtext(side=1, "Mean relatedness", outer = TRUE, line = -0.2, adj = 0.85, col ='black', cex = 0.4, font =2)
mtext(side=1, "Estimated effect", outer = TRUE, line = -0.2, adj = 0.25, col ='black', cex = 0.4, font =2)


dev.off()


 # FIGURE 4 ----

load("./output/analyses/MODEL3_CHAIN_1.RData")

intercept = mean(m3$Sol[,1])
b_sporulation = mean(m3$Sol[,'sporulation_score'])
b_abundance = mean(m3$Sol[,'within_host_relative_abundance'])

p_sporulation = round(summary(m3)$solutions['sporulation_score','pMCMC'], 2)
p_abundance = round(summary(m3)$solutions['within_host_relative_abundance','pMCMC'], 2)


pdf('./output/figures/FIGURE_4.pdf', width=(8.7/2.54), height=(8.7/2.54))

par(mfrow=c(2,2),mar=c(2,2,1,1))

# THEORY PLOT #1
m<-seq(0,0.05,length.out=101)
n<-50
r = 1/(n-(n-1)*(1-m)^2)
plot(NA,NA,xlim=c(min(m),max(m)),ylim=c(0,1),xlab="",ylab="",bty="l",
     cex.axis=0.4, xaxt = 'n', yaxt = 'n')

mtext(side = 2, 'Relatedness', font = 2, line = 1, cex = 0.5)
mtext(side = 1, 'Migration rate', font = 2, line = 1, cex = 0.5)
axis(1, at = seq(0,0.05, 0.01), labels = rep('', 6), tck = -0.02)
axis(1, at = seq(0,0.05, 0.01), lwd = 0, line  = -1, cex.axis =  0.5)
axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = rep('', 5), tck = -0.02)
axis(2, at = c(0, 0.25, 0.5, 0.75, 1), lwd = 0, line  = -0.7, cex.axis =  0.5)

lines(r~m,lwd=1)
m<-0.05
n<-seq(1,50,length.out=101)
r = 1/(n-(n-1)*(1-m)^2)


# THEORY PLOT #2

plot(NA,NA,xlim=c(min(n),max(n)),ylim=c(0,1),xlab="",ylab="",bty="l",
     cex.axis=0.4, xaxt = 'n', yaxt = 'n')

mtext(side = 2, 'Relatedness', font = 2, line = 1, cex = 0.5)
mtext(side = 1, 'Group size', font = 2, line = 1, cex = 0.5)
axis(1, at = seq(0,50, 10), labels = rep('', 6), tck = -0.02)
axis(1, at = seq(0,50, 10), lwd = 0, line  = -1, cex.axis =  0.5)
axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = rep('', 5), tck = -0.02)
axis(2, at = c(0, 0.25, 0.5, 0.75, 1), lwd = 0, line  = -0.7, cex.axis =  0.5)

lines(r~n,lwd=1)

# EMPIRICAL PLOT #1

plot(within_host_relatedness ~ sporulation_score, cex = 0.5,
     data=d,
     bty="l", col =  adjustcolor('forestgreen', 0.1),
     ylim = range(d$within_host_relatedness), cex.axis=0.4, xaxt = 'n', yaxt = 'n')
     
abline(intercept+mean(d_mean_109$mean_relative_abundance) * b_abundance,b_sporulation,lwd=1)

mtext(side = 2, 'Relatedness', font = 2, line = 1, cex = 0.5)
mtext(side = 1, 'Sporulation score', font = 2, line = 1, cex = 0.5)
axis(1, at = seq(0,0.5, 0.1), labels = rep('', 6), tck = -0.02)
axis(1, at = seq(0,0.5, 0.1), lwd = 0, line  = -1, cex.axis =  0.5)
axis(2, at = seq(-1,1,0.5), labels = rep('', 5), tck = -0.02)
axis(2, at = seq(-1,1,0.5), lwd = 0, line  = -0.7, cex.axis =  0.5)
text(paste0('(pMCMC = ', p_sporulation, ')'), x = 0.48, y = 0.43, cex = 0.4)


# EMPIRICAL PLOT #2

d<- d %>%
  group_by(species) %>%
  mutate(mean_relative_abundance = mean(within_host_relative_abundance)) %>%
  ungroup()

plot(within_host_relatedness ~ within_host_relative_abundance, cex = 0.5,
     ylim = range(d$within_host_relatedness),cex.axis=0.4, xaxt = 'n', yaxt = 'n',
     data=d,
     bty="l",
     col = adjustcolor('darkorange', 0.1))

abline(intercept+mean(d$sporulation_score)*b_sporulation,b_abundance,lwd=1)

mtext(side = 2, 'Relatedness', font = 2, line = 1, cex = 0.5)
mtext(side = 1, 'Within host relative abundance', font = 2, line = 1, cex = 0.5)
axis(1, at = seq(0,0.8, 0.1), labels = rep('', 9), tck = -0.02)
axis(1, at = seq(0,0.8, 0.1), lwd = 0, line  = -1, cex.axis =  0.5)

axis(2, at = seq(-1,1,0.5), labels = rep('', 5), tck = -0.02)
axis(2, at = seq(-1,1,0.5), lwd = 0, line  = -0.7, cex.axis =  0.5)

text(paste0('(pMCMC = ', p_abundance, ')'), x = 0.71, y = 0.69, cex = 0.4)


dev.off()












