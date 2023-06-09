#####################################################################################
# calculate_network_metrics.R                                                       #
#   Parses soil networks stored as text files and calculates network metrics for    #
#     each main network (all nodes and links) and each disconnected subnetwork.     #
#####################################################################################
source("network_metrics_functions.R")

#################################################################
# PARSING FUNCTIONS                                             #
#################################################################
# Main function for parsing network files, calculating metrics for main net and subnets
parse_network_file <- function(net_file) {
  # ID of photo
  soil_net_id <- strsplit(net_file, split = "_")[[1]][2]
  
  # Read in data file - first determine how many columns each row could have
  max_cols <- max(count.fields(net_file, sep = " "))
  # Read in max columns for each row - delete the other ones later
  network_data <- read.table(net_file, sep = " ", col.names = paste0("V", seq_len(max_cols)), 
                             fill = TRUE)
  # Metadata rows have number of nodes and number of links - these split up the data file
  metadata_rows <- which(is.na(network_data$V4))
  # First section of the data file is network membership
  network_members <- as.data.frame(network_data[metadata_rows[1]:metadata_rows[2], ])
  # Second section of the data file is the lengths of the links between nodes - node ID is from
  #   order of the network membership list
  link_lengths <- as.data.frame(network_data[metadata_rows[2]:length(network_data$V1), ])
  
  # Clean up network membership list - keep only real columns
  network_members <- network_members[ , c("V2", "V4", "V6")]
  # Rename columns - assumes structure of data to be network ID, x coordinate, y coordinate
  names(network_members) <- c("network_id", "x", "y")
  # Delete metadata row
  network_members <- network_members[which(! is.na(network_members$x)), ]
  # Add node ID
  network_members$node_id <- seq_along(network_members$x)
  
  # Clean up link lengths list- keep only real columns
  link_lengths <- link_lengths[ , c("V2", "V4", "V6")]
  # Rename columns - assumes structure of data to be node_from ID, node_to ID, link length
  names(link_lengths) <- c("node_from", "node_to", "length")
  # Delete metadata row - will only have one column
  link_lengths <- link_lengths[which(! is.na(link_lengths$node_to)), ]
  # Remove duplicates, e.g link from 1,2 is the same as 2,1
  link_lengths <- link_lengths[!duplicated(t(apply(link_lengths, 1, sort))), ]
  
  # Create igraph network from soil network structure data
  net <- graph_from_data_frame(link_lengths, directed = FALSE)
  # Set link lengths as edge weights so weights don't have to be specified as 'length' in later algorithms
  E(net)$weight <- E(net)$length 
  
  # Calculate metrics for each subnet
  subnets_data <- do.call(rbind, lapply(unique(network_members$network_id), parse_subnet, 
                                              network_members, link_lengths, soil_net_id))

  # Calculate exponent and fit of power law of link length distribution
  pl_linklength <- powerlaw_linklength_distribution(link_lengths)
  
  # Calculate exponent and fit of power law of spatial entropy
  pl_sp_entropy <- get_spatial_entropy(net)
  
  # Calculate convex hull area and density (uses hull area)
  chull_area <- get_chull_area(network_members)
  network_density <- vcount(net)/chull_area
  
  # Run metrics and output CSV data for main network
  net_data <- data.frame(
    subnet_id = NA,
    net_id = soil_net_id,
    sw = is_smallworld(net),
    n_nodes = vcount(net),
    n_links = ecount(net),
    avg_degree = get_avg_node_degree(net),
    gamma_index = get_gamma_index(net),
    beta_index = get_beta_index(net),
    diameter = get_network_diameter(net),
    cost = get_network_cost(net),
    grc = global_reach_centrality(net),
    chull_area = chull_area,
    network_density = network_density,
    linklength_powerlaw_exp = pl_linklength$power_law_exponent,
    linklength_powerlaw_xmin = pl_linklength$xmin,
    linklength_powerlaw_ks_stat = pl_linklength$ks_stat,
    linklength_powerlaw_ks_p = pl_linklength$ks_p,
    sp_entropy_powerlaw_exp = pl_sp_entropy$power_law_exponent,
    sp_entropy_powerlaw_xmin = pl_sp_entropy$xmin,
    sp_entropy_powerlaw_ks_stat = pl_sp_entropy$ks_stat,
    sp_entropy_powerlaw_ks_p = pl_sp_entropy$ks_p
  )
  
  # Combine main network metrics and subnet metrics to return
  all_net_data <- rbind(net_data, subnets_data)
  
  # Update progress bar
  i <<- i + 1
  setTxtProgressBar(pb, i)
  
  return(all_net_data)
  
}

# Function called by parse_network_file to parse each of the subnets that comprise the main network
parse_subnet <- function(subnet_id, network_members, link_lengths, network_id) {
  # Determine members of the subnet
  subnet_members <- network_members[network_members$network_id==subnet_id, ]

  # Get only link_lengths rows where the 'to' or 'from' node is in the subnet_members list
  subnet_links <- link_lengths[link_lengths$node_from %in% subnet_members$node_id |
                                link_lengths$node_to %in% subnet_members$node_id, ]

  # This shouldn't happen, but in case there's a single-node network, break early
  if (length(subnet_links$node_from) < 1) {
    null_csv_data <- data.frame(
      subnet_id = subnet_id,
      net_id = network_id,
      sw = FALSE,
      n_nodes = 1,
      n_links = 0,
      avg_degree = NA,
      gamma_index = NA,
      beta_index = NA,
      diameter = NA,
      cost = NA,
      grc = NA,
      chull_area = NA,
      network_density = NA,
      linklength_powerlaw_exp = NA,
      linklength_powerlaw_xmin = NA,
      linklength_powerlaw_ks_stat = NA,
      linklength_powerlaw_ks_p = NA,
      sp_entropy_powerlaw_exp = NA,
      sp_entropy_powerlaw_xmin = NA,
      sp_entropy_powerlaw_ks_stat = NA,
      sp_entropy_powerlaw_ks_p = NA 
    )
    return(null_csv_data)
  }
  
  # Build igraph object
  subnet <- graph_from_data_frame(subnet_links, directed = FALSE)
  E(subnet)$weight <- E(subnet)$length

  # Calculate exponent and fit of power law of link length distribution
  pl_linklength <- powerlaw_linklength_distribution(subnet_links)
  
  # Calculate exponent and fit of power law of spatial entropy
  pl_sp_entropy <- get_spatial_entropy(subnet)
  
  # Calculate convex hull area and density (uses hull area)
  chull_area <- get_chull_area(subnet_members)
  network_density <- vcount(subnet)/chull_area

  # Run metrics and output CSV data
  csv_data <- data.frame(
    subnet_id = subnet_id,
    net_id = network_id,
    sw = is_smallworld(subnet),
    n_nodes = vcount(subnet),
    n_links = ecount(subnet),
    avg_degree = get_avg_node_degree(subnet),
    gamma_index = get_gamma_index(subnet),
    beta_index = get_beta_index(subnet),
    diameter = get_network_diameter(subnet),
    cost = get_network_cost(subnet),
    grc = global_reach_centrality(subnet),
    chull_area = chull_area,
    network_density = network_density,
    linklength_powerlaw_exp = pl_linklength$power_law_exponent,
    linklength_powerlaw_xmin = pl_linklength$xmin,
    linklength_powerlaw_ks_stat = pl_linklength$ks_stat,
    linklength_powerlaw_ks_p = pl_linklength$ks_p,
    sp_entropy_powerlaw_exp = pl_sp_entropy$power_law_exponent,
    sp_entropy_powerlaw_xmin = pl_sp_entropy$xmin,
    sp_entropy_powerlaw_ks_stat = pl_sp_entropy$ks_stat,
    sp_entropy_powerlaw_ks_p = pl_sp_entropy$ks_p
  )

  return(csv_data)
 
}

# Different versions of above functions, with calculations of mean and SD link length and not
#   entropy or link length power law statistics, small-worldness
parse_network_file2 <- function(net_file) {
  # ID of photo
  soil_net_id <- strsplit(net_file, split = "_")[[1]][2]
  
  # Read in data file - first determine how many columns each row could have
  max_cols <- max(count.fields(net_file, sep = " "))
  # Read in max columns for each row - delete the other ones later
  network_data <- read.table(net_file, sep = " ", col.names = paste0("V", seq_len(max_cols)), 
                             fill = TRUE)
  # Metadata rows have number of nodes and number of links - these split up the data file
  metadata_rows <- which(is.na(network_data$V4))
  # First section of the data file is network membership
  network_members <- as.data.frame(network_data[metadata_rows[1]:metadata_rows[2], ])
  # Second section of the data file is the lengths of the links between nodes - node ID is from
  #   order of the network membership list
  link_lengths <- as.data.frame(network_data[metadata_rows[2]:length(network_data$V1), ])
  
  # Clean up network membership list - keep only real columns
  network_members <- network_members[ , c("V2", "V4", "V6")]
  # Rename columns - assumes structure of data to be network ID, x coordinate, y coordinate
  names(network_members) <- c("network_id", "x", "y")
  # Delete metadata row
  network_members <- network_members[which(! is.na(network_members$x)), ]
  # Add node ID
  network_members$node_id <- seq_along(network_members$x)
  
  # Clean up link lengths list- keep only real columns
  link_lengths <- link_lengths[ , c("V2", "V4", "V6")]
  # Rename columns - assumes structure of data to be node_from ID, node_to ID, link length
  names(link_lengths) <- c("node_from", "node_to", "length")
  # Delete metadata row - will only have one column
  link_lengths <- link_lengths[which(! is.na(link_lengths$node_to)), ]
  # Remove duplicates, e.g link from 1,2 is the same as 2,1
  link_lengths <- link_lengths[!duplicated(t(apply(link_lengths, 1, sort))), ]
  
  # Create igraph network from soil network structure data
  net <- graph_from_data_frame(link_lengths, directed = FALSE)
  # Set link lengths as edge weights so weights don't have to be specified as 'length' in later algorithms
  E(net)$weight <- E(net)$length 
  
  # Calculate metrics for each subnet
  subnets_data <- do.call(rbind, lapply(unique(network_members$network_id), parse_subnet2, 
                                        network_members, link_lengths, soil_net_id))
  
  # Calculate convex hull area and density (uses hull area)
  chull_area <- get_chull_area(network_members)
  network_density <- vcount(net)/chull_area
  
  # Run metrics and output CSV data for main network
  net_data <- data.frame(
    subnet_id = NA,
    net_id = soil_net_id,
    n_nodes = vcount(net),
    n_links = ecount(net),
    avg_degree = get_avg_node_degree(net),
    avg_link_length = get_avg_link_length(link_lengths),
    sd_link_length = get_sd_link_length(link_lengths),
    gamma_index = get_gamma_index(net),
    beta_index = get_beta_index(net),
    diameter = get_network_diameter(net),
    cost = get_network_cost(net),
    grc = global_reach_centrality(net),
    chull_area = chull_area,
    network_density = network_density
  )
  
  # Combine main network metrics and subnet metrics to return
  all_net_data <- rbind(net_data, subnets_data)
  
  # Update progress bar
  i <<- i + 1
  setTxtProgressBar(pb, i)
  
  return(all_net_data)
  
}
parse_subnet2 <- function(subnet_id, network_members, link_lengths, network_id) {
  # Determine members of the subnet
  subnet_members <- network_members[network_members$network_id==subnet_id, ]
  
  # Get only link_lengths rows where the 'to' or 'from' node is in the subnet_members list
  subnet_links <- link_lengths[link_lengths$node_from %in% subnet_members$node_id |
                                 link_lengths$node_to %in% subnet_members$node_id, ]
  
  # This shouldn't happen, but in case there's a single-node network, break early
  if (length(subnet_links$node_from) < 1) {
    null_csv_data <- data.frame(
      subnet_id = subnet_id,
      net_id = network_id,
      n_nodes = 1,
      n_links = 0,
      avg_degree = NA,
      avg_link_length = NA,
      sd_link_length = NA,
      gamma_index = NA,
      beta_index = NA,
      diameter = NA,
      cost = NA,
      grc = NA,
      chull_area = NA,
      network_density = NA
    )
    return(null_csv_data)
  }
  
  # Build igraph object
  subnet <- graph_from_data_frame(subnet_links, directed = FALSE)
  E(subnet)$weight <- E(subnet)$length
  
  # Calculate convex hull area and density (uses hull area)
  chull_area <- get_chull_area(subnet_members)
  network_density <- vcount(subnet)/chull_area
  
  # Run metrics and output CSV data
  csv_data <- data.frame(
    subnet_id = subnet_id,
    net_id = network_id,
    n_nodes = vcount(subnet),
    n_links = ecount(subnet),
    avg_degree = get_avg_node_degree(subnet),
    avg_link_length = get_avg_link_length(subnet_links),
    sd_link_length = get_sd_link_length(subnet_links),
    gamma_index = get_gamma_index(subnet),
    beta_index = get_beta_index(subnet),
    diameter = get_network_diameter(subnet),
    cost = get_network_cost(subnet),
    grc = global_reach_centrality(subnet),
    chull_area = chull_area,
    network_density = network_density
  )
  
  return(csv_data)
  
}

# Make list of files to parse to create network structure and calculate metrics (change path as required)
network_files <- list.files(path = "data/all_networks/", pattern = "IMG_*", full.names = TRUE)

# Set up progress bar
i <<- 0
pb <- txtProgressBar(min = 0, max = length(network_files), style = 3)
setTxtProgressBar(pb, i)

# Main function loop to parse files
soil_networks_data <- do.call(rbind, lapply(network_files, parse_network_file))
soil_networks_data2 <- do.call(rbind, lapply(network_files, parse_network_file2))

# Store output
write.csv(soil_networks_data, "output/all_networks_data.csv", row.names = FALSE)
write.csv(soil_networks_data2, "output/all_networks_data_linklengths.csv", row.names = FALSE)

