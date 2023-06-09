library(igraph)


# calculate modularity of simulated graphs

sim_graph_mod <- function(edges_df) {
  
  #library(dplyr)
  #library(tidyr)
  library(igraph)
  
  # get number of edges and nodes for simlation
  
  modularity <- list()
  
  all_nodes <- c(edges_df$author_1_ran_id, edges_df$author_2_ran_id)
  n_nodes = length(unique(all_nodes))
  #print(n_nodes)
  n_edges =  nrow(edges_df)
  #print(n_edges)
  weight_vector = edges_df$cosine_similarity
  #print(weight_vector)
  
  # generate random graphs with equal number of nodes, edges and edge weights
  
  for (i in 1:1000) {
    g <- sample_gnm(n_nodes, n_edges)
    E(g)$weight <- weight_vector
    
    fc <- cluster_fast_greedy(g, weights = E(g)$weight)
    mod <- modularity(fc)
    modularity[i] <- mod
    
  }
  
  hist(unlist(modularity))
  
  unlist(modularity)
  
  
}

# initiate dataframe to save data to
mydata <- data.frame('ind' = c(1:1000), stringasFactor = FALSE)

# call sim_mod function on files

path = 'P:/thesis/thesis_final_data/produced_data/co_authors_graph_data/2_temporal_co_authors_cos_sim/cos_sim_filt_0_1/'
file.names <- dir(path, pattern =".csv")

for(i in 1:length(file.names)){
  dec = '.'
  delim = ','  
  edges_df <- read.csv(paste(path, file.names[i], sep = ''), dec=dec, sep = delim, stringsAsFactors=FALSE, header = TRUE)
  
  str(edges_df)
  print(unique(edges_df$cosine_similarity))
  
  a <- sim_graph_mod(edges_df)
  
  time_period = paste('time_period_', as.character(i), sep='')
  mydata[time_period] = a
}

# save mod sim data to csv

path2 = 'P:/thesis/thesis_final_data/produced_data/co_authors_graph_data/5_coauthor_mod_sim_july_2020/'
write.csv(mydata, paste(path2, 'coautor_sim_mod_cos001.csv', sep = ''), row.names = FALSE)







# calculate modularity of (non simulated) co-authorship graphs

graph_mod <- function(edges_df) {
  
  #library(dplyr)
  #library(tidyr)
  library(igraph)
  
  # generate graph object
  
  # rename columns in edges df
  names(edges_df)[names(edges_df)=="author_1_ran_id"] <- "from"
  names(edges_df)[names(edges_df)=="author_2_ran_id"] <- "to"
  names(edges_df)[names(edges_df)=="cosine_similarity"] <- "weight"
  
  col_order <- c("from", "to", "weight")
  edges_df2 <- edges_df[, col_order]
  
  # generate graph
  g <- graph_from_data_frame(edges_df2, directed = FALSE)
  E(g)$weight <- edges_df2$weight
  
  # calculate modularity
  fc <- cluster_fast_greedy(g, weights = E(g)$weight)
  mod <- modularity(fc)
  
  mod
  
  
}

# call mod function

path2 = 'P:/thesis/thesis_final_data/produced_data/co_authors_graph_data/5_coauthor_mod_sim_july_2020/'
file.names <- dir(path, pattern =".csv")

mod_list_2 = list()

for(i in 1:length(file.names)){
  dec = '.'
  delim = ','  
  edges_df <- read.csv(paste(path, file.names[i], sep = ''), dec=dec, sep = delim, stringsAsFactors=FALSE, header = TRUE)
  
  #str(edges_df)
  #print(unique(edges_df$cosine_similarity))
  
  a <- graph_mod(edges_df)
  print(a)
  
  time_period = paste('time_period_', as.character(i), sep = '')
  mod_list_2[time_period] = a
}

# initiate dataframe to save data to
mod_list_2






#plot

library(ggplot2)

#http://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization

# prep data for ggplot

mod_value <- apply(mydata, 2, FUN=min) # this creates a named numeric vector, which I later turn into a dataframe
min_df <- data.frame(mod_value)
min_df['group'] = 'min mod sim'
min_df
df <- cbind(time_period = rownames(min_df), min_df)
rownames(df) <- 1:nrow(df)
df2 = df[-1,]
df3 = df2[-1,]
df3
df3['legend']  = 'simulated min max bounds'



mod_value <- apply(mydata, 2, FUN=max)
max_df <- data.frame(mod_value)
max_df['group'] = 'max mod sim'
max_df
df_2 <- cbind(time_period = rownames(max_df), max_df)
rownames(df_2) <- 1:nrow(df_2)
df_3 = df_2[-1,]
df_4 = df_3[-1,]
df_4
df_4['legend']  = 'simulated min max bounds'


# make dataframe from original modularities
timeper <- names(mod_list_2)
mod_list_3 <- as.numeric(mod_list_2)
mod_list_3
mod_df_4 <- data.frame("time_period" = timeper, "mod_value" = mod_list_3)
mod_df_4$ group = 'original'
mod_df_4
mod_df_4['legend']  = 'mod co-authorship network'


new <- rbind(df3, df_4, mod_df_4)
new


p<-ggplot(new, aes(x=time_period, y=mod_value, group=group)) +
  geom_line(aes(color=legend))+
  geom_point(aes(color=legend)) +
  labs(title="Modularity of co-authorship network",x="time period", y = "modularity")
p
