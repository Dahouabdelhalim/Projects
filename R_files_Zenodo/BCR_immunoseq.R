# Kenneth B. Hoehn
# 8/8/2021
# Assign clonal clusters and build phylogenetic trees for Mlynarczyk et al. Science
# from mouse immunoseq data
# Figures 7G and S16A-G

# Run these commands first
# VERSION="1.17.1"
# wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/${VERSION}/ncbi-igblast-${VERSION}-x64-linux.tar.gz
# tar -zxf ncbi-igblast-${VERSION}-x64-linux.tar.gz
# 
# # Download reference databases and setup IGDATA directory
# fetch_igblastdb.sh -o igblast
# cp -r ncbi-igblast-${VERSION}/internal_data igblast
# cp -r ncbi-igblast-${VERSION}/optional_file igblast
# # Build IgBLAST database from IMGT reference sequences
# fetch_imgtdb.sh -o germlines/imgt

library(alakazam)
library(ggplot2)
library(dplyr)
library(shazam)
library(scoper)
library(ggtree)
library(dowser)
library(scoper)
library(ggpubr)
library(apTreeshape)
library(tidyr)

print(sessionInfo())

t=read.table("../Immunoseq_Bcl2.Btg1_tumors_RearrangementDetails.tsv",sep="\\t",head=TRUE)

glimpse(t)

# Filter sequences to have at least 5 templates (optional)
ft = filter(t, templates > 4)

# Convert data table to fasta file
n = paste0(">",1:nrow(ft),"|templates=",ft$templates,"|sample_name=",ft$sample_name,"\\n",ft$rearrangement)
writeLines(n,con="immunoseq.fasta")

# filter files
s = data.frame(Biostrings::readDNAStringSet("../germlines/mouse/vdj/imgt_mouse_IGHV.fasta"))
s$seq = s[,1]
s$id = rownames(s)
s = s[grepl("C57",s$id),]

writeLines(paste0(">",s$id,"\\n",s$seq),con="germlines/imgt/mouse/vdj/imgt_mouse_IGHV.fasta")

#Now run these commands in command line within Docker image to align sequences to 
# VDJ references. Could also set these up to run locally:
# https://changeo.readthedocs.io/en/stable/examples/igblast.html

# # Download and extract IgBLAST#
# Make new blast database for mice
# imgt2igblast.sh -i germlines/imgt -o igblast
#
# #Run docker (need to pull image if not already there)
# sudo docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.2.0 bash
#
# #Align to VDJ sequences 
# AssignGenes.py igblast -s immunoseq.fasta -b igblast \\
#     --organism mouse --loci ig --format blast
# MakeDb.py igblast -i immunoseq_igblast.fmt7 -s immunoseq.fasta \\
#     -r germlines/imgt/mouse/vdj/imgt_mouse_IGH* \\
#     --extended --failed
# exit

ft = readChangeoDb("immunoseq_igblast_db-pass.tsv")

group_palette = c("Bcl2"="#000000", "Bcl2+Q36H"="#FF9200") 

ft = filter(ft, productive)

# check clonal thresholds
dist_cross = distToNearest(filter(ft,locus=="IGH"),
        sequenceColumn="junction", 
        vCallColumn="v_call", jCallColumn="j_call",
        model="ham", normalize="len", nproc=1,
        cross="sample_name")

pdf("intermediates/crossDistance.pdf",height=20,width=8)
ggplot(subset(dist_cross, !is.na(cross_dist_nearest)), 
             aes(x=cross_dist_nearest)) + 
    theme_bw() + 
    xlab("Cross-sample_id Hamming distance") + 
    ylab("Count") +
    geom_histogram(color="white", binwidth=0.02) +
    geom_vline(xintercept=0.12, color="firebrick", linetype=2) +
    facet_grid( ~ ., scales="free_y")
dev.off()

plots = list()
pdf("intermediates/dist_to_nearest.pdf",width=6,height=6)
for(sample_name in unique(dist_cross$sample_name)){
   print(sample_name)
   temp = filter(dist_cross,!!sample_name==sample_name)
   dist_ham <- distToNearest(filter(temp,locus=="IGH"), sequenceColumn="junction", 
                    vCallColumn="v_call", jCallColumn="j_call",
                    model="ham", normalize="len", nproc=2)
   output <- findThreshold(dist_ham$dist_nearest, method="density")
   threshold <- output@threshold
   g = ggplot(subset(dist_ham, !is.na(dist_nearest)),aes(x=dist_nearest,
   ,y = ..density..)) + 
     theme_bw() + 
     xlab("Hamming distance") + 
     ylab("Count") +
     scale_x_continuous(breaks=seq(0, 1, 0.1)) +
     geom_histogram(color="white", binwidth=0.02) +
     ggtitle(paste(sample_name,threshold))+
     geom_histogram(
       aes(x=cross_dist_nearest,y = -..density..),
       color="white", binwidth=0.02,fill="black")+
     xlim(0,max(filter(dist_cross,
       !is.na(cross_dist_nearest))$cross_dist_nearest))+
     geom_vline(xintercept=0.1,color="grey")
   if(!is.na(threshold)){
       g = g + geom_vline(xintercept=threshold, color="firebrick", linetype=2)
   }
   plots[[sample_name]] = g
   print(g)
}
dev.off()

# Identify clones
clones = ft %>%
    group_by(sample_name) %>%
    do(as.data.frame(
    hierarchicalClones(., 0.1, only_heavy = TRUE, 
      cdr3 = TRUE, nproc = 7,verbose = FALSE, log = NULL,
      summarize_clones = TRUE)))

clones = ungroup(clones)
clones$clone_id = paste0(clones$sample_name,"-",clones$clone_id)

# Reconstruct germlines
references = readIMGT(dir = "germlines/imgt/mouse/vdj")
germ = createGermlines(clones,references,nproc=1,organism="mouse")

# calculate SHM frequency in the V gene
comb_germline <- observedMutations(germ, sequenceColumn="sequence_alignment",
       germlineColumn="germline_alignment_d_mask",
       regionDefinition=IMGT_V,
       frequency=TRUE,
       combine=TRUE, 
       nproc=3)

writeChangeoDb(comb_germline,"intermediates/all_cloned_data.tsv")
comb_germline = readChangeoDb("intermediates/all_cloned_data.tsv")

comb_germline$clone_id = paste0(comb_germline$sample_name, comb_germline$clone_id)

# group sequences into cohorts
comb_germline$Cohort = ""
comb_germline$Cohort[grepl("tumor231d_BCL2_",comb_germline$sample_name)] = "Bcl2"
comb_germline$Cohort[grepl("tumor231d_BCL2\\\\+Q36H_",comb_germline$sample_name)] = "Bcl2+Q36H"
comb_germline$Cohort[grepl("tumorendpoint_BCL2_",comb_germline$sample_name)] = "tumorendpoint_BCL2"
comb_germline$Cohort[grepl("tumorendpoint_BCL2\\\\+Q36H_",comb_germline$sample_name)] = "tumorendpoint_BCL2+Q36H"

sample2mouse = 
c("tumor231d_BCL2_284"="mouse 1",
"tumor231d_BCL2_285"="mouse 2",
"tumor231d_BCL2_286"="mouse 3",
"tumor231d_BCL2_287"="mouse 4",
"tumor231d_BCL2_288"="mouse 5",
"tumor231d_BCL2+Q36H_289"="mouse 10",
"tumor231d_BCL2+Q36H_290"="mouse 7",
"tumor231d_BCL2+Q36H_291"="mouse 8",
"tumor231d_BCL2+Q36H_292"="mouse 9",
"tumor231d_BCL2+Q36H_293"="mouse 6")

my_comparisons = list(c("Bcl2", "Bcl2+Q36H"))

data = filter(comb_germline, Cohort %in% c("Bcl2","Bcl2+Q36H"))
data$mouse = sample2mouse[data$sample_name]
data$mouse = factor(data$mouse, levels = paste("mouse",1:10))

boxsize=0.15
bracketsize=0.15
axissize=0.15
panelsize=0.25
sigsize=2
pointsize=0.1
fontsize=9

max_shm = max(data$mu_freq)
# plot SHM by mouse
pdf("results/shm_big.pdf",width=3.5,height=2.5,useDingbats=FALSE)
print(
  ggplot(data,aes(x=mouse,y=mu_freq))+
  geom_boxplot(outlier.shape=NA)+theme_bw()+
  geom_jitter(width=0.01,alpha=0.1,height=0,
    aes(size=templates, color=Cohort))+
  theme(axis.text.x = element_text(angle = 90, 
    vjust = 0.5, hjust=1))+
  xlab("") + ylab("V gene SHM/site")+
  labs(size="Templates")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(color="Cohort")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=group_palette)+
  theme(axis.text.x = element_text(angle = 90, 
  vjust = 0.5, hjust=1))+
  ylim(0,max_shm)
)
dev.off()

means = data %>%
  group_by(Cohort,mouse) %>%
  summarize(mean_shm = mean(mu_freq),
    weighted_mean_shm = weighted.mean(x=mu_freq,
      w=templates))

max_mean_shm = max(means$mean_shm)

pdf("results/shm_small.pdf",width=2.5,height=2.5,useDingbats=FALSE)
print(
  ggboxplot(filter(means),
    x = "Cohort", y = "mean_shm", outlier.shape=NA,size=boxsize) +
  geom_point(size= 1.5, aes(color=Cohort))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("V gene SHM/site")+
  xlab("")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(color="Cohort")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=group_palette)+
  theme(axis.text.x = element_text(angle = 90, 
  vjust = 0.5, hjust=1))+
  ylim(0, max_mean_shm*1.07)
  )
dev.off()

# sequences are unique
data %>% group_by(mouse) %>% summarize(n=n(), u=n_distinct(sequence),diff=n-u)

# do diversity analysis with unique sequences
sample_curve <- alphaDiversity(data,
  group="sample_name", clone="clone_id",
  min_q=1.9, max_q=2.1, step_q=0.1,min_n=30,
  ci=0.95, nboot=1000, uniform=TRUE)

aq2 = ungroup(filter(sample_curve@diversity,
  q==2))

# group sequences into cohorts
aq2$Cohort = ""
aq2$Cohort[grepl("tumor231d_BCL2_",aq2$sample_name)] = "Bcl2"
aq2$Cohort[grepl("tumor231d_BCL2\\\\+Q36H_",aq2$sample_name)] = "Bcl2+Q36H"
aq2$Cohort[grepl("tumorendpoint_BCL2_",aq2$sample_name)] = "tumorendpoint_BCL2"
aq2$Cohort[grepl("tumorendpoint_BCL2\\\\+Q36H_",aq2$sample_name)] = "tumorendpoint_BCL2+Q36H"
aq2$mouse = factor(sample2mouse[aq2$sample_name], levels = paste("mouse",1:10))

max_d = max(aq2$d_upper)

# plot SHM by mouse
pdf("results/diversity_big.pdf",width=3.5,height=2.5,useDingbats=FALSE)
print(
  ggplot(aq2,aes(x=mouse,y=d,ymin=d_lower,ymax=d_upper))+
  geom_point(size=2,aes(color=Cohort))+
  geom_errorbar(width=0.3)+
  theme(axis.text.x = element_text(angle = 90, 
    vjust = 0.5, hjust=1))+
  xlab("") + ylab("Simpson's diversity")+
  labs(size="Templates")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(color="Cohort")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=group_palette)+
  theme(axis.text.x = element_text(angle = 90, 
  vjust = 0.5, hjust=1))+
  ylim(0, max_d*1.02)
)
dev.off()

max_mean_d = max(aq2$d)
pdf("results/diversity_small.pdf",width=2.5,height=2.5,useDingbats=FALSE)
print(
  ggboxplot(filter(aq2),
    x = "Cohort", y = "d", outlier.shape=NA,size=boxsize) +
  geom_point(size= 1.5, aes(color=Cohort))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("Simpson's diversity")+
  xlab("")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(color="Cohort")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=group_palette)+
  theme(axis.text.x = element_text(angle = 90, 
  vjust = 0.5, hjust=1))+
  coord_cartesian(expand=TRUE)+
  ylim(0, max_mean_d*1.07)
  )
dev.off()

# remove CDR3s for SHM and tree building
mgerm = maskSequences(data, mask_codons=FALSE, mask_cdr3=TRUE, nproc=3)

# Format clones for tree building
clones = formatClones(mgerm, columns=c("mouse","Cohort","sample_name"),
  num_fields="templates", seq="sequence_masked")

# tally number of templates within each clone
clones$templates = unlist(lapply(clones$data,
  function(x)sum(x@data$templates)))

trees = tibble()
for(sample in unique(clones$sample_name)){
  # build IgPhyML trees for clones with at least 100 templates
  print(sample)
  tt = getTrees(filter(clones, templates > 100 &
    sample_name == sample), nproc=6, build="igphyml", omega="ce",
  exec="~/Dropbox/Projects/IgPhyML_development/igphyml/src/igphyml")
  trees = bind_rows(trees, tt)
}
saveRDS(trees, "intermediates/igphyml_trees.rds")
trees = readRDS("intermediates/igphyml_trees.rds")

# Make branches equal to total SHM, rather than SHM/site
trees = scaleBranches(trees)

# order by total template count
trees = trees[order(trees$templates,decreasing=TRUE),]
max = max(unlist(lapply(trees$data,function(x)x@data$templates)))
max_div = max(unlist(lapply(trees$trees,function(x){
    max(cophenetic(x)[,"Germline"])
  })))

for(sample in unique(trees$sample_name)){
  pt = filter(trees,sample_name==sample)
  p = plotTrees(pt,scale=5)
  p = lapply(p,function(x){
    x + geom_tippoint(aes(size=templates)) +
    scale_size(limits=c(1,max)) + 
    xlim(0, max_div)
    })
  treesToPDF(p,file=paste0("results/sample_trees/",sample,"_igphyml_trees.pdf"))
}

# plot the top tree in each mouse
toptrees = trees %>%
  group_by(Cohort,mouse) %>%
  slice_max(templates,n=1)

toptrees = toptrees %>%
  arrange(toptrees$Cohort,toptrees$mouse)

max = max(unlist(lapply(toptrees$data,function(x)x@data$templates)))
max_div = max(unlist(lapply(toptrees$trees,function(x){
    max(cophenetic(x)[,"Germline"])
  })))
max_temps = max(unlist(lapply(toptrees$data,function(x){
    max(x@data$templates)
  })))
min_temps = min(unlist(lapply(toptrees$data,function(x){
    min(x@data$templates)
  })))
 
p = plotTrees(toptrees,scale=5)
p = lapply(1:length(p),function(x){
    p[[x]] + geom_tippoint(aes(size=templates),
      pch=21, fill=group_palette[toptrees$Cohort[x]]) +
    scale_size(limits=c(min_temps, max_temps),
      range=c(0.2,2)) +
    xlim(0, max_div) +
    ggtitle(toptrees$mouse[x])
    })

pnoscale = lapply(p, function(x){
  x + theme(legend.position="none")
})

pdf("results/toptrees.pdf", width=5, height=4, useDingbats=FALSE)
gridExtra::grid.arrange(grobs=pnoscale,ncol=5)
gridExtra::grid.arrange(grobs=p,ncol=5)
dev.off()

# Get IgPhyML parameter values
params = bind_rows(toptrees$parameters)
params$mouse = toptrees$mouse
params$Cohort = toptrees$Cohort

gr = params %>%
  select(-clone,-lhood,-nsite) %>%
  gather("key","value", -mouse,-Cohort,-nseq, -tree_length)

gr$parameter = toupper(lapply(strsplit(gr$key, split="_"),function(x)x[1]))
gr$parameter = factor(gr$parameter,
  levels=c("OMEGA","KAPPA","WRC","GYW","WA","TW","GRS","SYC"))

pdf("results/params.pdf",width=3.5,height=2.5)
  ggboxplot(filter(gr,nseq > 0 & parameter != "OMEGA"),
    x = "parameter", y = "value", outlier.shape=NA, size=boxsize) +
  geom_jitter(size= 1, height=0, width=0.1, aes(color=Cohort))+
  ylab("Estimate")+
  xlab("HLP19 model parameter")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(color="Cohort")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=group_palette)+
  theme(axis.text.x = element_text(angle = 90, 
  vjust = 0.5, hjust=1))+
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-1,NA)

ggboxplot(filter(gr,nseq > 10 & parameter != "OMEGA"),
    x = "parameter", y = "value", outlier.shape=NA, size=boxsize) +
  geom_jitter(size= 1, height=0, width=0.1, aes(color=Cohort))+
  ylab("Estimate")+
  xlab("HLP19 model parameter")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(color="Cohort")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=group_palette)+
  theme(axis.text.x = element_text(angle = 90, 
  vjust = 0.5, hjust=1))+
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-1,NA)
dev.off()

pdf("results/omega_big.pdf",width=3.5,height=2.5,useDingbats=FALSE)
print(
  ggplot(params,aes(x=mouse,y=omega_mle,ymin=omega_lci,ymax=omega_uci))+
  geom_point(size=2,aes(color=Cohort))+
  geom_errorbar(width=0.3)+
  theme(axis.text.x = element_text(angle = 90, 
    vjust = 0.5, hjust=1))+
  xlab("") + ylab("dN/dS")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(color="Cohort")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=group_palette)+
  theme(axis.text.x = element_text(angle = 90, 
  vjust = 0.5, hjust=1))
)
dev.off()

max_o = max(params$omega_mle)
min_o = min(params$omega_mle)
pdf("results/omega_small.pdf",width=2.5,height=2.5,useDingbats=FALSE)
print(
  ggboxplot(params,
    x = "Cohort", y = "omega_mle", outlier.shape=NA,size=boxsize) +
  geom_point(size= 1.5, aes(color=Cohort))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("dN/dS")+
  xlab("")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(color="Cohort")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=group_palette)+
  theme(axis.text.x = element_text(angle = 90, 
  vjust = 0.5, hjust=1))+
  coord_cartesian(expand=TRUE)+
  ylim(min_o, max_o*1.07)
  )
dev.off()
