library(phylowood)
library(ape)
library(tidyverse)
library(ggtree)
library(viridis)
library(phytools)
library(adephylo)
library(ggnewscale)
library(wesanderson)
library(forcats)
library(gridExtra)



###################################
#Part A Visualize the number of species per family on a cladogram of families
###################################

# load tree from the Harris 2016 paoper
tr <- read.tree("input/knb.1177.1_Harris2016.data")

## prune to eudicots
drp <- listTips(tr)

keep <- drp[[which(unlist(lapply(drp, function(k){"eupteleaceae" %in% names(k) & "salicaceae" %in% names(k) & !"ceratophyllaceae" %in% names(k) })))]]
keep <- names(keep)
tr <- keep.tip(tr, keep)

# get the number of species per family
fams <- phylowood_sp %>% 
  filter(geo_insularwoody) %>% 
  dplyr::select(family = tax_family, tax_name) %>% 
  distinct() %>% 
  mutate(family = tolower(family)) %>% 
  group_by(family) %>% 
  dplyr::summarize(iw_species = n()) %>% 
  ungroup() %>% 
  right_join(as_tibble(tr$tip.label), by = c("family" = "value")) %>% 
  mutate(iw_species = ifelse(is.na(iw_species), 0, iw_species)) %>% 
  arrange(factor(family, levels = tr$tip.label))%>%  
  mutate(iw_species_log = log(iw_species + 0.1)) %>% 
  mutate(iw_species_log = ifelse(is.infinite(iw_species_log), 0, iw_species_log))

# save for the phylogenetic signal analysis
save(fams, file = "output/number_of_IW_species_per_family.rda")
save(tr, file = "output/eudicot_phylogeny.rda")

dat <- fams$iw_species_log
names(dat) <- fams$family
dat <- data.frame(dat)

# plot the tree and identify superclades
ggtree(tr) +
  geom_tiplab(size=2)

ggsave(file = "output/family_Tree.pdf", height = 30)

#supersosids: Peridiscaceae & salicaceae

# superasterids: berberidopsidaceae & acanthaceae
# clade labels
d <- data.frame(node=c(getMRCA(phy = tr, tip = c("peridiscaceae", "salicaceae")),
                       getMRCA(phy = tr, tip = c("berberidopsidaceae", "acanthaceae"))),
                type=c("Super-rosids", "super-asterids"))


# ONE FAMILY IS MISSING HERE, add to figure labels
phylowood_sp[!tolower(phylowood_sp$tax_family) %in% fams$family, "tax_family"] %>%  distinct()
phylowood_sp[phylowood_sp$tax_family == "Kewaceae",]

# because it does not have insular woody representatives
sel <- rownames(dat)[which(dat$dat > 0)]

lab <- tr$tip.label
lab[!lab %in% sel] <- " "

tr$tip.label <- str_to_title(lab)
rownames(dat) <- str_to_title(rownames(dat))

# Number of shifts
fam_num <- phylowood_isl %>% 
  left_join(phylowood_gen %>%  dplyr::select(tax_genus, tax_family, tax_apgIV), by = "tax_genus") %>% 
  group_by(tax_family) %>% 
  summarize(shifts = sum(trt_minshift_nr)) %>% 
  filter(!is.na(tax_family)) %>% 
  bind_rows(data.frame(tax_family = rownames(dat)[!rownames(dat) %in% .$tax_family])) %>% 
  mutate(log_shifts = log(shifts + 0.1)) %>% 
  replace_na(list(log_shifts = 0))

shifts <- data.frame(fam_num$log_shifts,
                     row.names = fam_num$tax_family)


#plot
p11 <- ggtree(tr, layout="circular", size = 0.1) +
  geom_tiplab(offset=60, size=2)+ 
  geom_hilight(data=d, aes(node=node, fill=type), alpha = 0.3)+
  scale_fill_manual(values=c("steelblue", "darkgreen"))

p11 <- p11 +  new_scale_fill()
p11 <-  gheatmap(p11,
                 dat, 
                 width=0.15, 
                 low="blue",
                 high="red",
                 colnames = FALSE, 
                 offset = -10, 
                 font.size = 12)+
  scale_fill_gradient(name = "Number of\\nIW\\nspecies",
                      high= "black",
                      low = "grey90",
                    breaks = c(0, log(10), log(100), log(300)),
                    labels = c(0, 10, 100, 300),
                    limits = c(0,log(350)),
                    na.value = "grey95")
#  scale_fill_viridis(option="D", name="Number of\\nDW\\nspecies", direction =-1)+
  # scale_fill_viridis(option="D",
  #                    name = "Number of\\nIW\\nspecies",
  #                    direction =1,
  #                    breaks = c(0, log(10), log(100)),
  #                    labels = c(0, 10, 100),
  #                    na.value = "grey95")

p2 <- p11 +  new_scale_fill()
p3 <- gheatmap(p2, 
         shifts,
         offset=15,
         width=.15, 
         colnames = FALSE)+
  scale_fill_gradient(high= "black",
                      low = "grey90",
                       name="Number of\\ntransitions",
                       breaks = c(0, log(1), log(10), log(40), log(50)),
                       labels = c(0, 1, 10, 40, 50),
                       na.value = "grey95")

ggsave(plot = p3, "output/figure_cladogram_species_per_family.pdf", height = 7, width = 8)


# top 12 
shifts10 <- fam_num %>% 
  arrange(desc(shifts))%>% 
  slice(1:10) %>% 
  mutate(tax_family = tolower(tax_family))

species10 <- fams %>% 
  arrange(desc(iw_species))%>% 
  slice(1:10)

top10 <- full_join(shifts10, species10, by = c("tax_family" = "family"))

# chose families for plots four top shifts and two top species
li <- c("asteraceae", "amaranthaceae", "brassicaceae", "rubiaceae", "campanulaceae", "gesneriaceae")

bub_pos <- data.frame(x = c(1,1.1,1.2),
                      y = c(1,1.1,1))

# load the data from the preparation
load(file = "output/iw_perarchipelago_family.rda")
load(file = "output/iw_pergenus_family.rda")
load(file = "output/iw_perhabitat_family.rda")
load(file = "output/iw_permoisturepreference_family.rda")

# get the colour palette for the archipelagos
arch_num %>% 
  filter(tolower(tax_family) %in% li)
pal <- wes_palette("Zissou1", 53, type = "continuous")
pal2 <- wes_palette("Cavalcanti1", 3)[c(2,3,1)]
pal3 <- wes_palette("Darjeeling1", 5)[c(5,2,3)]
pal4 <- wes_palette("GrandBudapest2", 4)[3]

for(i in 1:length(li)){
  print(i)
  #Species proportion top archipelagos
  sub <- arch_num %>% filter(tolower(tax_family) == li[i]) %>% 
    bind_cols(bub_pos) %>% 
    mutate(lab = paste(geo_archipelago,
                       paste(round(frac,0),"%", sep = " "), sep = "\\n"))
  
  
  #Archipelago
 arch <-  ggplot(sub)+
    geom_point(aes(x= x, y = y, size = frac, color = frac))+
    geom_text(aes(x = x, y = y, 
                  color = frac, 
                  label = lab),
              nudge_y = 0.07)+
    scale_colour_gradientn(colors = pal, limits=c(8,53))+
    scale_size(range = c(min(sub$frac), max(sub$frac)))+
    scale_x_continuous(limits = c(0.95, 1.6))+
    scale_y_continuous(limits = c(0.95, 1.6))+
    ggtitle(paste(unique(sub$tax_family),
                  " (",
                  fam_num %>% filter(tolower(tax_family) == li[i]) %>% pull(shifts),
                  " shifts, ",
                  top10 %>% filter(tax_family == li[i]) %>% pull(iw_species),
                  " species)",
                  sep = ""))+
    theme_void()+
    theme(legend.position = "none",
          plot.title = element_text(vjust = -70))
  
  
  #genera
  sub <- gen_num%>% filter(tolower(tax_family) == li[i])
  
  
  gen <- ggplot(sub)+
    geom_bar(aes(y = fct_reorder(tax_genus,frac), x = frac), 
             stat = 'identity', 
             fill= pal4)+
    theme_bw()+
    xlab("Proportion of species")+
    theme(axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # habitat open/forest
  sub <- for_op%>% filter(tolower(tax_family) == li[i]) %>% 
    mutate(geo_habitat_broad = factor(geo_habitat_broad,
                                       levels = c(NA, "forest", "widespread", "open"),
                                      exclude = NULL))
  
  open <- ggplot(sub)+
    geom_bar(aes(y = tax_family, x = n, 
                 fill = geo_habitat_broad), 
             position="stack", stat="identity")+
   # scale_fill_manual(values = c( "darkgreen", "orange", "lightyellow", "grey"))+
    scale_fill_manual(values = pal2)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    theme_bw()+
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none")
  
  #moisture preference
  sub <- moist%>% filter(tolower(tax_family) == li[i]) %>% 
    mutate(geo_moisture_preference = factor(geo_moisture_preference,
                                      levels = c(NA, "humid", "mesic", "arid"),
                                      exclude = NULL))
  
  moi <- ggplot(sub)+
    geom_bar(aes(y = tax_family, x = n, 
                 fill = geo_moisture_preference), 
             position="stack", stat="identity")+
    scale_fill_manual(values = pal3)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    theme_bw()+
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none")
  
  com <- grid.arrange(
    grobs = list(arch ,gen, open, moi),
    heights = c(8, 1, 1),
    layout_matrix = rbind(c(1, 2),
                          c(3, 3),
                          c(4, 4)))
  
  ggsave(com, file = paste("output/family_figures/", li[i], ".pdf", sep = ""))
}
