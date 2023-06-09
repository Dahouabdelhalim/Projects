source ("Library.R")

Setting <- list(theme(axis.title.x = element_text(size = 50, vjust = -0.5)),
                 theme(axis.title.y = element_text(size = 50, vjust = 2)),
                 theme(axis.text.x = element_text(size = 40, colour = "black", vjust = 0.5)),
                 theme(axis.text.y = element_text(size = 40, colour = "black")))


#2A Accumulative distribution --------------------------------------------------

Map = readr::read_csv("Accumulative_distribution.csv") %>%
  dplyr::select(-1)

Type_name = c("Calocasiae", "Dmelanogaster")
FLY = c("C. alocasiae", "D. melanogaster")
gl = list()
gl_zoom = list()
for (i in 1:2) {

  data = Map %>%
    filter(Type == Type_name[i])
  
  g <- ggplot(data, aes(x = Difference_x, y = Difference_y)) +
    stat_bin2d(bins = 240) +
    scale_fill_gradient(low = "grey80", high = "black", limits = c(0, 30)) +
    coord_fixed() +
    theme_bw() +
    xlab("Distance (mm)") +
    ylab("Distance (mm)") +
    scale_x_continuous(limits = c(-120, 120), breaks = seq(-120, 120, by = 60)) +
    scale_y_continuous(limits = c(-120, 120), breaks = seq(-120, 120, by = 60)) +
    labs(fill = "Number of flies") +
    theme(legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 18),
          legend.direction = "horizontal", 
          legend.key.size = grid::unit(1.0, "cm"), 
          legend.position = c(0.76, 0.86)) +
    guides(fill = guide_colorbar(title.position = "top")) +
    Setting
  gl <- c(gl, list(g))
  
  g_zoom <- ggplot(data, aes(x = Difference_x, y = Difference_y)) +
    stat_bin2d(bins = 40) +
    scale_fill_gradient(low = "grey80", high = "black", limits = c(0, 30)) +
    coord_fixed() +
    theme_bw() +
    xlab("Distance (mm)") +
    ylab("Distance (mm)") +
    scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
    scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10)) +
    guides(fill = FALSE) +
    Setting
  gl_zoom <- c(gl_zoom, list(g_zoom))
}

A1 <- cowplot::plot_grid(gl[[1]], NULL, gl_zoom[[1]],
                         nrow = 1,
                         rel_widths = c(10, 0.5, 10)) +
  ggtitle("C. alocasiae") +
  theme(plot.title = element_text(face = "italic", size = 60, colour = "deepskyblue4", hjust = 0.5, vjust = 2))

A2 <- cowplot::plot_grid(gl[[2]], NULL, gl_zoom[[2]],
                         nrow = 1,
                         rel_widths = c(10, 0.5, 10)) +
  ggtitle("D. melanogaster") +
  theme(plot.title = element_text(face = "italic", size = 60, colour = "orange2", hjust = 0.5, vjust = 2))


Fig2A <- cowplot::plot_grid(A1, NULL, A2,
                            nrow = 1,
                            rel_widths = c(10, 1, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))

ggsave("Fig2/2A.png", Fig2A, width = 30, height = 8, dpi = 100)


#2B Inter individual distance --------------------------------------------------

Distance = readr::read_csv("Inter_individual_distance.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Calocasiae(IR)")

Distance_random = readr::read_csv("Inter_individual_distance_random.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Calocasiae(IR)") %>%
  rename(Movie_number = Trial_number)

Size_mean = readr::read_csv("Size_mean_sd.csv") %>%
  dplyr::select(-1)
size_Calo = Size_mean$C.alo_mean
size_Dmel = Size_mean$D.mel_mean
Size <- c(size_Calo, size_Dmel)

data_all = bind_rows(data.frame(Data = as.character("Original"), Distance),
                     data.frame(Data = as.character("Random"), Distance_random))

Type_name = c("Calocasiae", "Dmelanogaster")
Data_name = c("Original", "Random")
COLOR = c("deepskyblue4", "orange2")
FLY = c("C. alocasiae", "D. melanogaster")
gl = list()
for (i in 1:2) {
  for (j in 1:2) {
    
    data = data_all %>%
      filter(Type == Type_name[i], Data == Data_name[j])
    
    #Convert vertical axis from count to percent
    g_each <- ggplot(data, aes(x = Distance, group = factor(Movie_number))) + 
      geom_histogram(binwidth = 2, boundary = 0)
    g_each_info <- ggplot_build(g_each)
    data_each = as.data.frame(g_each_info$data) %>%
      group_by(group) %>%
      mutate(percent = (count/sum(count)) * 100) %>%
      ungroup() %>%
      dplyr::select(group, x, percent) %>%
      rename(Movie_number = group)
    
    g_mean <- ggplot(data, aes(x = Distance)) + 
      geom_histogram(binwidth = 2, boundary = 0)
    g_mean_info <- ggplot_build(g_mean)
    data_mean = as.data.frame(g_mean_info$data) %>%
      mutate(percent = (count/sum(count)) * 100) %>%
      dplyr::select(x, percent)
    
    
    Range = Size[i] * 2
    if(j == 2){
      k = "black"
    }else{
      k = COLOR[i]
    }
    
    g <- ggplot() +
      geom_line(data = data_each, aes(x = x, y = percent, group = factor(Movie_number)), size = 1.0, alpha = 0.3, color = k) +
      geom_point(data = data_each, aes(x = x, y = percent, group = factor(Movie_number)), size = 3, alpha = 0.4, color = k) +
      geom_line(data = data_mean, aes(x = x, y = percent), size = 2.0, color = k) +
      geom_point(data = data_mean, aes(x = x, y = percent), size = 4.5, color = k) +
      geom_vline(xintercept = Range, linetype = "dashed", size = 1) +
      theme_bw() +
      xlab("Distance (mm)") +
      scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 20)) +
      scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 10)) +
      Setting +
      theme(axis.title.y = element_blank())
    g
    gl <- c(gl, list(g))
  }
}

Fig2B <- cowplot::plot_grid(NULL, gl[[1]], NULL, gl[[3]],
                            NULL, NULL, NULL, NULL,
                            NULL, gl[[2]], NULL, gl[[4]],
                            nrow = 3,
                            rel_widths = c(2, 20, 4, 20),
                            rel_heights = c(8, -1.6, 8)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
Fig2B

ggsave("Fig2/2B.png", Fig2B, width = 30, height = 12, dpi = 100)


#Overall -----------------------------------------------------------------------
Fig2 <- cowplot::plot_grid(Fig2A, NULL, Fig2B,
                             nrow = 3, 
                             rel_heights = c(8, 2, 12))

ggsave("Fig2/Fig2.png", Fig2, width = 30, height = 22, dpi = 100)