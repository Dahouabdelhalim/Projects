########################
### Custom functions ###
########################

# Modify transitions with user provided values and recovery

modified_transition_custom <- function(fire_stack_early = NULL,
                                       fire_stack_late = NULL,
                                       early_surv_mult = 1,
                                       late_surv_mult = 1,
                                       early_tran_mult = 1,
                                       late_tran_mult = 1,
                                       early_recr_mult = 1,
                                       late_recr_mult = 1) {
  
  fun <- function (transition_array, landscape, timestep) {
    
    if(timestep == 1) {
      early_fire_bool_demo <<- FALSE
      late_fire_bool_demo <<- FALSE
    }
    
    values <- raster::getValues(landscape[[fire_stack_early]][[timestep]])
    values <- values[which(!is.na(values))]
    idx <- which(values == 1)
    
    if(length(idx) > 0) {
      early_fire_bool_demo <<- TRUE
      idx_early_demo <<- idx
      for (i in idx) {
        transition_array[1, 1, i] <- transition_array[1, 1, i] * early_surv_mult
        transition_array[2, 2, i] <- transition_array[2, 2, i] * early_surv_mult
        transition_array[2, 1, i] <- transition_array[2, 1, i] * early_tran_mult
        transition_array[1, 2, i] <- transition_array[1, 2, i] * early_recr_mult
      }
    }
    
    if(early_fire_bool_demo & timestep %in% seq(4, 124, 6)) {
      for (i in idx_early_demo) {
        transition_array[1, 1, i] <- transition_array[1, 1, i] * (early_surv_mult + (1 - early_surv_mult) / 4)
        transition_array[2, 2, i] <- transition_array[2, 2, i] * (early_surv_mult + (1 - early_surv_mult) / 4)
        transition_array[2, 1, i] <- transition_array[2, 1, i] * (early_tran_mult + (1 - early_tran_mult) / 4)
        transition_array[1, 2, i] <- transition_array[1, 2, i] * (early_recr_mult + (1 - early_recr_mult) / 4)
      }
    }
    
    if(early_fire_bool_demo & timestep %in% seq(5, 125, 6)) {
      for (i in idx_early_demo) {
        transition_array[1, 1, i] <- transition_array[1, 1, i] * (early_surv_mult + (1 - early_surv_mult) / 2)
        transition_array[2, 2, i] <- transition_array[2, 2, i] * (early_surv_mult + (1 - early_surv_mult) / 2)
        transition_array[2, 1, i] <- transition_array[2, 1, i] * (early_tran_mult + (1 - early_tran_mult) / 2)
        transition_array[1, 2, i] <- transition_array[1, 2, i] * (early_recr_mult + (1 - early_recr_mult) / 2)
      }
    }
    
    if(early_fire_bool_demo & timestep %in% seq(6, 126, 6)) {
      for (i in idx_early_demo) {
        transition_array[1, 1, i] <- transition_array[1, 1, i] * (early_surv_mult + 3 * (1 - early_surv_mult) / 4)
        transition_array[2, 2, i] <- transition_array[2, 2, i] * (early_surv_mult + 3 * (1 - early_surv_mult) / 4)
        transition_array[2, 1, i] <- transition_array[2, 1, i] * (early_tran_mult + 3 * (1 - early_tran_mult) / 4)
        transition_array[1, 2, i] <- transition_array[1, 2, i] * (early_recr_mult + 3 * (1 - early_recr_mult) / 4)
      }
    }
    
    values <- raster::getValues(landscape[[fire_stack_late]][[timestep]])
    values <- values[which(!is.na(values))]
    idx <- which(values == 1)
    
    if(length(idx) > 0) {
      late_fire_bool_demo <<- TRUE
      idx_late_demo <<- idx
      for (i in idx) {
        transition_array[1, 1, i] <- transition_array[1, 1, i] * late_surv_mult
        transition_array[2, 2, i] <- transition_array[2, 2, i] * late_surv_mult
        transition_array[2, 1, i] <- transition_array[2, 1, i] * late_tran_mult
        transition_array[1, 2, i] <- transition_array[1, 2, i] * late_recr_mult
      }
    }
    
    if(late_fire_bool_demo & timestep %in% seq(6, 126, 6)) {
      for (i in idx_late_demo) {
        transition_array[1, 1, i] <- transition_array[1, 1, i] * (late_surv_mult + (1 - late_surv_mult) / 2)
        transition_array[2, 2, i] <- transition_array[2, 2, i] * (late_surv_mult + (1 - late_surv_mult) / 2)
        transition_array[2, 1, i] <- transition_array[2, 1, i] * (late_tran_mult + (1 - late_tran_mult) / 2)
        transition_array[1, 2, i] <- transition_array[1, 2, i] * (late_recr_mult + (1 - late_recr_mult) / 2)
      }
    }
    
    if (timestep %in% seq(2, 122, 6)) {
      early_fire_bool_demo <<- FALSE
      late_fire_bool_demo <<- FALSE
      idx_early_demo <<- NULL
      idx_late_demo <<- NULL
    }
    
    # print(transition_array[ , , 28860])
    transition_array
    
  }
}

# Modify habitat with user provided values and recovery

habitat_suitability_mod <- function(fire_stack_early = NULL,
                                    fire_stack_late = NULL,
                                    reduction_early = 1,
                                    reduction_late = 1,
                                    noise = 0) {
  
  fun <- function (landscape, timestep) {
    
    if(timestep == 1) {
      early_fire_bool <<- FALSE
      late_fire_bool <<- FALSE
    }
    
    cell_idx <- which(raster::getValues(landscape[[fire_stack_early]][[timestep]]) == 1)
    
    if(length(cell_idx) > 0) {
      early_fire_bool <<- TRUE
      cell_idx_early <<- cell_idx
      reduction_early_draw <<- pmin(pmax(rnorm(1, reduction_early, noise), 0), 1)
      landscape$suitability[cell_idx] <- landscape$suitability[cell_idx] * reduction_early_draw
    }
    
    if(early_fire_bool & timestep %in% seq(4, 124, 6)) {
      landscape$suitability[cell_idx_early] <- pmin(landscape$suitability[cell_idx_early] + (1 - reduction_early_draw) / 4, 1)
    }
    
    if(early_fire_bool & timestep %in% seq(5, 125, 6)) {
      landscape$suitability[cell_idx_early] <- pmin(landscape$suitability[cell_idx_early] + (1 - reduction_early_draw) / 4, 1)
    }
    
    if(early_fire_bool & timestep %in% seq(6, 126, 6)) {
      landscape$suitability[cell_idx_early] <- pmin(landscape$suitability[cell_idx_early] + (1 - reduction_early_draw) / 4, 1)
    }
    
    if(early_fire_bool & timestep %in% seq(7, 121, 6)) {
      landscape$suitability[cell_idx_early] <- pmin(landscape$suitability[cell_idx_early] + (1 - reduction_early_draw) / 4, 1)
    }
    
    cell_idx <- which(raster::getValues(landscape[[fire_stack_late]][[timestep]]) == 1)
    
    if(length(cell_idx) > 0) {
      late_fire_bool <<- TRUE
      cell_idx_late <<- cell_idx
      reduction_late_draw <<- pmin(pmax(rnorm(1, reduction_late, noise), 0), 1)
      landscape$suitability[cell_idx] <- landscape$suitability[cell_idx] * reduction_late_draw
    }
    
    if(late_fire_bool & timestep %in% seq(6, 126, 6)) {
      landscape$suitability[cell_idx_late] <- pmin(landscape$suitability[cell_idx_late] + (1 - reduction_late_draw) / 2, 1)
    }
    
    if(late_fire_bool & timestep %in% seq(7, 121, 6)) {
      landscape$suitability[cell_idx_late] <- pmin(landscape$suitability[cell_idx_late] + (1 - reduction_late_draw) / 2, 1)
    }
    
    if (timestep %in% seq(2, 122, 6)) {
      early_fire_bool <<- FALSE
      late_fire_bool <<- FALSE
      cell_idx_early <<- NULL
      cell_idx_late <<- NULL
      reduction_late_draw <<- NULL
      reduction_early_draw <<- NULL
    }
    
    landscape
    
  }
}

# Relate carrying capacity to habitat suitability

k_function <- function(landscape, timestep) {
  
  suit <- landscape$suitability
  
  max_ind / (1 + exp(-(suit - 0.5) / 0.05))
}

# Extract populations

return_pop_data <- function(results_object){
  n_sims <- length(results_object)
  n_timesteps <- length(results_object[[1]])
  idx <- which(!is.na(raster::getValues(results_object[[1]][[1]]$population[[1]])))
  pop_mat <- matrix(NA, nrow = n_sims, ncol = n_timesteps)
  for (i in seq_len(n_sims)){
    for (j in seq_len(n_timesteps)) {
      pop_mat[i, j] <- sum(raster::extract(results_object[[i]][[j]]$population, idx))
    }
  }
  pop_mat
}

# Plot population trend

plot_pop <- function(pop_data, upper_limit, title = "") {
  max_pop <- max(pop_data)
  min_pop <- min(pop_data)
  mean_trend <- colSums(pop_data) / nrow(pop_data)
  matplot(t(pop_data),
          type = 'l',
          lwd = 0.5,
          lty = 1,
          ylim = c(0, 1.05 * upper_limit),
          col = 'grey',
          xlab = 'Year',
          ylab = 'Population size',
          xaxt = 'n')
  lines(1:126, mean_trend, lwd = 1.5)
  abline(h = round(mean(apply(pop_data, 1, function(x) min(x))), 0), lwd = 1, lty = 2)
  axis(1, at = c(0, 30, 60, 90, 120), labels = c(0, 5, 10, 15, 20))
  text(0.75, upper_limit, labels = title, pos = 4, cex = 0.8)
}

# Transform maps for 3D plotting

transform_polygons <- function(spatial_object, shift_matrix) {
  out_fort <- fortify(spatial_object, region = "id")
  out_fort <- left_join(out_fort, spatial_object@data, c("id" = "id"))
  xy <- as.matrix(out_fort[ , c("long", "lat")]) %*% shift_matrix
  out_fort$x <- xy[ , 1]; out_fort$y = xy[ , 2]
  out_fort
}

# Plot expected minimum populations

plot_ema <- function (pop_data_list, p = 0.05, all_points = FALSE) {
  
  # extract number of simulations
  n_species <- length(pop_data_list)
  species_names <- names(pop_data_list)
  n_sims <- length(pop_data_list[[1]])
  sim_names <- names(pop_data_list[[1]])
  
  # calculate z-score
  z <- qnorm(1 - p / 2)
  
  # populate table with ema mean and error values
  df <- foreach (i = seq_len(n_species), .combine = rbind) %do% {
    foreach (j = seq_len(n_sims), .combine = rbind) %do% {
      pops <- pop_data_list[[i]][[j]]
      min_total_pops <- apply(pops, 1, function(x) min(x))
      ema_mean <- mean(min_total_pops)
      ema_sd <- sd(min_total_pops)
      ema_lower <- ema_mean - z * (ema_sd / sqrt(length(min_total_pops)))
      ema_upper <- ema_mean + z * (ema_sd / sqrt(length(min_total_pops)))
      data.frame("name" = as.factor(paste0(species_names[[i]])),
                 "sim" = as.factor(paste0(sim_names[[j]])),
                 "ema_mean" = ema_mean,
                 "ema_sd" = ema_sd,
                 "ema_lower" = ema_lower,
                 "ema_upper" = ema_upper)
    }
  }
  
  max_limits <- aggregate(df$ema_upper, by = df["name"], function(x) ceiling(max(x)))
  
  # get the raw values also if all_points is set to TRUE
  if (all_points) {
    df_raw <- foreach (i = seq_len(n_species), .combine = rbind) %do% {
      foreach (j = seq_len(n_sims), .combine = rbind) %do% {
        pops <- pop_data_list[[i]][[j]]
        min_total_pops <- apply(pops, 1, function(x) min(x))
        data.frame("name" = as.factor(paste0(species_names[[i]])),
                   "sim" = as.factor(paste0(sim_names[[j]])),
                   "min_pops" = min_total_pops)
      }
    }
    max_limits_points <- aggregate(df_raw$min_pops, by = df_raw["name"], function(x) ceiling(max(x)))
  }
  
  p <- list()
  for (i in seq_len(n_species)) {
    p[[i]] <- ggplot(df[df$name == species_names[i], ], aes(x = sim, y = ema_mean, ymin = ema_lower, ymax = ema_upper)) +
      {if (all_points)geom_jitter(data = df_raw[df_raw$name == species_names[i], ], mapping = aes(y = min_pops), position = position_jitter(0.2), col = "lightgrey", size = 0.1)} +
      geom_errorbar(width = .2, position = position_dodge(0.05)) +
      geom_point(mapping = aes(shape = sim, fill = sim), color = "black", size = 1.7) +
      geom_hline(yintercept = 0, size = 1.2) +
      # facet_grid(.~name,
      #            switch = "x", # Moves the labels from the top to the bottom
      #            #labeller = label_both # Adds the labels to the year and X variables
      # ) +
      xlab(species_names[i]) +
      {if (i == 1)ylab(paste0("Expected minimum abundance"))} +
      {if (i != 1)ylab(paste0(""))} +
      scale_y_continuous(limits = c(0, max_limits[i, 2]), expand = c(0, 0)) +
      {if (all_points)scale_y_continuous(limits = c(0, max_limits_points[i, 2]), expand = c(0, 0))} +
      scale_fill_manual(values = c('black', 'white', 'black', 'white')) +
      scale_shape_manual(values = c(21, 21, 24, 24)) +
      theme(
        strip.background = element_blank(),
        #axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.line.y = element_line(color="black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position="none"
      )
  }
  
  do.call(grid.arrange, list("grobs" = p, "ncol" = 3, "padding" = 0))
}


## Combined plots
# plot_ema <- function (pop_data_list, interval = 95, all_points = FALSE) {
#   
#   # extract number of simulations
#   n_species <- length(pop_data_list)
#   species_names <- names(pop_data_list)
#   n_sims <- length(pop_data_list[[1]])
#   sim_names <- names(pop_data_list[[1]])
#   
#   # calculate interval range
#   interval_range <- c((100 - interval) / 2, 100 - (100 - interval) / 2) / 100
#   
#   # populate table with ema mean and error values
#   df <- foreach (i = seq_len(n_species), .combine = rbind) %do% {
#     foreach (j = seq_len(n_sims), .combine = rbind) %do% {
#       pops <- pop_data_list[[i]][[j]]
#       min_total_pops <- apply(pops, 1, function(x) min(x))
#       ema_mean <- mean(min_total_pops)
#       ema_lower <- stats::quantile(min_total_pops, interval_range)[1]
#       ema_upper <- stats::quantile(min_total_pops, interval_range)[2]
#       data.frame("name" = as.factor(paste0(species_names[[i]])),
#                  "sim" = as.factor(paste0(sim_names[[j]])),
#                  "ema_mean" = ema_mean,
#                  "ema_lower" = ema_lower,
#                  "ema_upper" = ema_upper)
#     }
#   }
#   
#   max_limit <- ceiling(max(df$ema_upper))
#   
#   # get the raw values also if all_points is set to TRUE
#   if (all_points) {
#     df_raw <- foreach (i = seq_len(n_species), .combine = rbind) %do% {
#       foreach (j = seq_len(n_sims), .combine = rbind) %do% {
#         pops <- pop_data_list[[i]][[j]]
#         min_total_pops <- apply(pops, 1, function(x) min(x))
#         data.frame("name" = as.factor(paste0(species_names[[i]])),
#                    "sim" = as.factor(paste0(sim_names[[j]])),
#                    "min_pops" = min_total_pops)
#       }
#     }
#     max_limit <- ceiling(max(df_raw[, -c(1:2)]))
#   }
#   
#   ggplot(df, aes(x = sim, y = ema_mean, ymin = ema_lower, ymax = ema_upper)) +
#     {if (all_points)geom_jitter(data = df_raw, mapping = aes(y = min_pops), position = position_jitter(0.2), col = "lightgrey", size = 0.1)} +
#     geom_errorbar(width = .2, position = position_dodge(0.05)) +
#     geom_point(mapping = aes(shape = sim, fill = sim), color = "black", size = 2.0) +
#     geom_hline(yintercept = 0, size = 1.2) +
#     facet_grid(.~name,
#                switch = "x", # Moves the labels from the top to the bottom
#                #labeller = label_both # Adds the labels to the year and X variables
#     ) + 
#     xlab("") +
#     ylab(paste0("Expected minimum abundance")) +
#     scale_y_continuous(limits = c(0, max_limit), expand = c(0, 0)) +
#     scale_fill_manual(values = c('black', 'white', 'black', 'white')) +
#     scale_shape_manual(values = c(21, 21, 24, 24)) +
#     theme(
#       strip.background = element_blank(),
#       axis.title.x = element_blank(),
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_rect(fill = "white"),
#       panel.border = element_blank(),
#       axis.line.y = element_line(color="black"),
#       legend.title = element_blank(),
#       legend.key = element_blank()
#     )
# }