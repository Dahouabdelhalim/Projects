library(ape)
#library(phylowood)
library(tidyverse)
library(phytools)

load(file = "output/species_list_islands_insularwoody.rda")
dat <- sp_list %>% 
  mutate(onisland = ifelse(onisland, "b", "a"))

# GIFT island species
load(file = "output/GIFT_island_checklist_including_incomplete.rda")

gift <- checklists_broad %>% 
  filter(native == 1) %>% 
  select(species) %>% 
  mutate(species = gsub(" " , "_", species)) %>% 
  distinct()

# Phylogenetic tree
tr <- read.tree("input/GBMB_v01_20201027.tre")

li <- unique(sp_list$clade_nr)

#derived woody species
load("input/dw.rda")

# a for loop across families
## define function to calculate node probabilities from Liam Revell's blog
foo<-function(x){
  y<-sapply(x$maps,function(x) names(x)[1])
  names(y)<-x$edge[,1]
  y<-y[as.character(length(x$tip)+1:x$Nnode)]
  return(y)
}

# a function to extract times per state in time bins
get_state_ages_make.simmap <- function(x, bin_size = 2, option = c("shifts", "time")){
  node_ages <- data.frame(depth = node.depth.edgelength(x), 
                          node_id = 1:nrow(data.frame(node.depth.edgelength(x))))
  
  ages <- data.frame(x$edge) %>% 
    left_join(node_ages, by = c("X1" = "node_id")) %>% 
    left_join(node_ages, by = c("X2" = "node_id")) %>% 
    mutate(depth.x = depth.x - max(depth.y))%>% 
    mutate(depth.y = depth.y - max(depth.y))
  
  res_shifts <-  data.frame()
  res_time <-  data.frame()
  
  for(z in 1:length(x$maps)){
    sub <- x$maps[[z]]
    
    if(length(sub) > 1){
      out_time <- data.frame(start = c(ages$depth.x[z], ages$depth.x[z] +  cumsum(sub)[-length(sub)]),
                             end = c(ages$depth.x[z] + cumsum(sub)),
                             state = names(sub))

      if("shifts" %in% option){
        out_shifts <- out_time[-1,]
        out_shifts <-out_shifts %>% 
          select(shift_age = start,
                 shift_type = state) %>% 
          mutate(shift_type = paste("to_", shift_type, sep = ""))
      }
      
    }else{
      out_time <- data.frame(start = ages$depth.x[z],
                        end = ages$depth.y[z],
                        state = names(sub))
      out_shifts <- data.frame(shift_age = ages$depth.x[z],
                               shift_type = "none")
    }
    
    if("time" %in% option){
      res_time <- rbind(res_time, out_time)
    }
    
    if("shifts" %in% option){
      res_shifts <- rbind(res_shifts, out_shifts) %>% 
        filter(shift_type != "none")
    }
  }
  
  bins <- seq((min(ages$depth.x) %/% bin_size)* bin_size,0, bin_size)
  
  if("time" %in% option){
    for(z in 1:(length(bins)-1)){
      bin <- ifelse(res_time$start < bins[z] & res_time$end < bins[z] & !res_time$start > bins[z], res_time$end - res_time$start, 0)
      bin <- ifelse(res_time$start < bins[z] & res_time$end > bins[z], (bins[z] - res_time$start) , bin)
      
      res_time <- cbind(res_time, bin)
    }
      
      names(res_time) <- c(names(res_time[1:3]), paste("midpoint", bins[-length(bins)] + bin_size/2, sep = ""))
      res_time <- as.data.frame(res_time)
      rownames(res_time) <- NULL
      
      #summarize per state and bin
      suma <- res_time %>% 
        pivot_longer(cols = starts_with("midpoint"),
                     names_to = "bin", 
                     values_to = "time") %>% 
        group_by(state,  bin) %>%
        summarize(time = sum(time)) %>% 
        pivot_wider(id_cols = bin, names_from = state, values_from = time) %>% 
        mutate(rel_a = a / (a + b)) %>% 
        mutate(rel_b = b / (a + b)) %>% 
        mutate(rel_a = ifelse(is.nan(rel_a), 0, rel_a)) %>% 
        mutate(rel_b = ifelse(is.nan(rel_b), 0, rel_b))
  }
  
  if("shifts" %in% option){
    suma_shifts <- res_shifts %>% 
      mutate(shifts_bin = cut(shift_age, breaks = bins, labels = (bins[-length(bins)]+ bin_size/2)))
  }

  
  if(length(option) == 1){
    if(option == "time"){
      return(suma)
    }
    
    if(option == "shifts"){
      return(suma_shifts)
    }
  }

  if(length(option) == 2){
    return(list(suma, suma_shifts))
  }
}

# Run the analyses
for(i in 1:length(li)){
  print(paste(i, "/", length(li), sep = ""))
  
  # subset list to one genus
  sub <- dat %>% filter(clade_nr %in% li[i]) %>%
    filter(!is.na(species)) %>% 
    distinct()
  
  # get MRCA of this family
  mrca <- ape::getMRCA(phy = tr, tip = sub$species)
  
  # in case there is only one tip for family, get the closest relatives
  # if the entire family is not in the tree, do nothing
  if(all(sub$dwood == "b") & length(sub$species) > 1){
    
    mrca <- phytools::getParent(tr, mrca)
    
  }else if(length(sub$species) == 1){
    
    sis <- getSisters(tr, node = sub$value)
    mrca <- phytools::getParent(tr, sis)
    mrca <- getParent(tr, mrca)
    
  }else if(length(sub$species) > 1){
    crown <- extract.clade(phy = tr, node = mrca)
    
    # delete names with .sp or .cf
    del <- crown$tip.label %>% 
      as_tibble() %>% 
      filter(grepl("sp\\\\.", value) | grepl("cf\\\\.", value) | grepl("\\\\.", value)| grepl("_x_", value)) %>% 
      unlist()
    
    crown <- drop.tip(crown, del)
    
    # remove the species removed from the phylogeny from sub
    sub <- sub %>% 
      filter(species %in% crown$tip.label)
    
    # clades with only 3 species dont run properly
    while(length(crown$tip.label) < 4){
      mrca <- getParent(tr, mrca)
      crown <- extract.clade(phy = tr, node = mrca)
    }
    
    # Derived woodiness
    #match the order of the traits to the tree lables
    trait <- sub %>%
      filter(species %in% crown$tip.label) %>% 
      select(dwood) %>%
      unlist()
    
    names(trait) <- sub$species
    
    # add potential species that are in the phylogeny but not in sub
    missing <- crown$tip.label[!crown$tip.label %in% names(trait)] %>%
      as_tibble() %>%
      mutate(trait = ifelse(value %in% gsub(" ", "_", dw$tax_name), "b", "a"))
    
    add <- unlist(missing$trait)
    names(add) <- missing$value
    
    trait <- c(trait, add)
    
    # match to the names in the phylogeny
    trait <- trait[match(names(trait), crown$tip.label)]
    
    ## reconstruct ancestral area, http://blog.phytools.org/2013/03/new-totally-rewritten-version-of.html
    ##  maybe also try BiSSE: https://www.zoology.ubc.ca/~fitzjohn/diversitree.docs/asr-bisse.html
    mtrees_wood <-make.simmap(tree = crown,
                              x = trait,
                              model="ARD",
                              nsim=100)
    

    
    #calculate node probabilities
    plo_wood <- summary(mtrees_wood)
    
    # if there is only one derived woody species the model sometimes does not converge properly
    k <- 1
    
    while(!is.matrix(plo_wood$count) & k<10){
      mtrees_wood<-make.simmap(tree = crown,
                               x = trait,
                               model="ARD",
                               nsim=100)
      
      #calculate node probabilities
      plo <- summary(mtrees_wood)
      k <- k+1
    }
    
    # return the time in each state per branch segment per timebin
    state_per_bin <- lapply(mtrees_wood, get_state_ages_make.simmap)

    # fro the time per state per bin
    time <- lapply(state_per_bin, "[[", 1)
    time <- bind_rows(time, .id = "replicate") %>% 
      group_by(bin) %>% 
      summarize(a = mean(a),
                b = mean(b),
                rel_a = mean(rel_a),
                rel_b = mean(rel_b)) %>% 
      ungroup() %>% 
      mutate(age = parse_number(bin)) %>% 
      mutate(clade_ID = li[i])
    
    save(time, file = file.path("output","ancestral_state_reconstruction", 
                                          paste(li[i], "_state_per_bin_woodiness_time", ".rda", sep = "")))
    
    #shifts per bin
    shifts <- lapply(state_per_bin, "[[", 2)
    shifts <- bind_rows(shifts, .id = "replicate") %>% 
      group_by(shifts_bin, shift_type, replicate) %>%
      count() %>% 
      ungroup() %>% 
      group_by(shifts_bin, shift_type) %>% 
      summarize(n = mean(n))
    
    save(time, file = file.path("output","ancestral_state_reconstruction", 
                                paste(li[i], "_state_per_bin_woodiness_shifts", ".rda", sep = "")))
    
    # island occurrence
    #match the order of the traits to the tree labels
    trait <- sub %>%
      select(onisland) %>%
      unlist()
    
    names(trait) <- sub$species
    
    # add potential species that are in the phylogeny but not in sub
    missing <- crown$tip.label[!crown$tip.label %in% names(trait)] %>%
      as_tibble() %>%
      mutate(trait = ifelse(value %in% gift$species, "b", "a"))
    
    add <- unlist(missing$trait)
    names(add) <- missing$value
    
    trait <- c(trait, add)
    
    
    # match to the names in the phylogeny
    trait <- trait[match(names(trait), crown$tip.label)]
    
    mtrees_island <-make.simmap(tree = crown,
                                x = trait,
                                model="ARD",
                                nsim=100)
    
    #calculate node probabilities
    plo_island <- summary(mtrees_island)
    
    # if there is only one derived woody species the model sometimes does not converge properly
    k <- 1
    
    while(!is.matrix(plo_island$count) & k<10){
      mtrees_island<-make.simmap(tree = crown,
                                 x = trait,
                                 model="ARD",
                                 nsim=100)
      
      #calculate node probabilities
      plo <- summary(mtrees_island)
      k <- k+1
    }
    
    # return the time in each s tate per branch segment per timebin
    state_per_bin <- lapply(mtrees_island, get_state_ages_make.simmap)
    state_per_bin <- bind_rows(state_per_bin, .id = "replicate") %>% 
      group_by(bin) %>% 
      summarize(a = mean(a),
                b = mean(b),
                rel_a = mean(rel_a),
                rel_b = mean(rel_b)) %>% 
      ungroup() %>% 
      mutate(age = parse_number(bin)) %>% 
      mutate(clade_ID = li[i])
    
    save(state_per_bin, file = file.path("output","ancestral_state_reconstruction", 
                                         paste(li[i], "_state_per_bin_island", ".rda", sep = "")))
    
    # create density maps
    obj_wood<-densityMap(mtrees_wood,states=c("a","b"),plot=FALSE)
    obj_island<-densityMap(mtrees_island,states=c("a","b"),plot=FALSE)
    obj_island$tree$tip.label <- rep("", length(obj_island$tree$tip.label))
    
    # plot to .pdf
    pdf(file.path("output","ancestral_state_reconstruction", paste(li[i], ".pdf", sep = "")),
        height = 45, width = 11)
    
    par(mfcol = c(1,2))
    plot(obj_wood,lwd=4,outline=TRUE,fsize=c(0.5),legend=50, direction = "rightwards")
    title("Derived woodiness", line = -1)
    plot(obj_island, lwd=4,outline=TRUE,legend=50, direction = "leftwards")
    title("Presence on island", line = -1)
    dev.off()
  }
}
