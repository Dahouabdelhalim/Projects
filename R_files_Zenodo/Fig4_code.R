source ("Library.R")

Setting <- list(theme(axis.title.x = element_text(size = 50, vjust = -0.5)),
                 theme(axis.title.y = element_text(size = 50, vjust = 2)),
                 theme(axis.text.x = element_text(size = 40, colour = "black", vjust = 0.5)),
                 theme(axis.text.y = element_text(size = 40, colour = "black")))


#4AB Correlation between speed and distance ------------------------------------

frame_adjust = 30/128

Parameters = readr::read_csv("Average_params.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Calocasiae(IR)") %>%
  mutate(across(Speed, ~ .x * frame_adjust)) #convert to every second

Type_name = c("Calocasiae", "Dmelanogaster")
COLOR = c("deepskyblue4", "orange2")

gl <- list()
for (i in 1:length(Type_name)) {
  
  dis_max = 80
  if(i == 1){
    vel_max = 1.5
    com_max = 1
    xscale = 3
  }else{
    vel_max = 3
    com_max = 1.5
    xscale = 3
  }
  
  data = Parameters %>%
    filter(Type == Type_name[i]) %>%
    mutate(Movie_number = as.character(Movie_number))
  
  #All movies
  g_all <- ggplot(data, aes(x = Speed, y = Distance, col = Movie_number)) +
    geom_path(size = 1.2,
              alpha = 0.5) +
    theme_bw() +
    labs(x = "Speed (mm/sec)",
         y = "Distance (mm)") +
    scale_color_brewer(palette = "Paired") +
    scale_x_continuous(limits = c(0, vel_max),
                       breaks = seq(0, vel_max, by = vel_max/xscale)) +
    scale_y_continuous(limits = c(0, dis_max),
                       breaks = seq(0, dis_max, by = dis_max/2)) +
    guides(col = FALSE) +
    Setting
  
  
  #One movie
  if(i == 1){
    data_one = data %>%
      filter(Movie_number == 3)
  }else{
    data_one = data %>%
      filter(Movie_number == 1)
  }
  
  g_one <- ggplot(data_one, aes(x = Speed, y = Distance, col = Time)) +
    geom_path(size = 1, alpha = 0.5) +
    theme_bw() +
    labs(x = "Speed (mm/sec)",
         y = "Distance (mm)") + 
    scale_x_continuous(limits = c(0, vel_max),
                       breaks = seq(0, vel_max, by = vel_max/xscale)) +
    scale_y_continuous(limits = c(0, dis_max),
                       breaks = seq(0, dis_max, by = dis_max/2)) +
    scale_color_gradient(name = "Time (min)",
                         low = "green3",
                         high = "magenta3",
                         limits = c(0, 20),
                         breaks = seq(0, 20, by = 10)) +
    guides(col = guide_colorbar(title.position = "top")) +
    theme(legend.title = element_text(size = 25),
          legend.text = element_text(size = 23),
          legend.direction = "horizontal",
          legend.key.size = grid::unit(1, "cm"),
          legend.position = c(0.80, 0.13)) +
    Setting
  
  if(i == 1){
    Title <- list(ggtitle("C. alo"),
                  theme(plot.title = element_text(face = "italic", size = 60, hjust = 0.5, vjust = 1, colour = COLOR[i])))
  }else{
    Title <- list(ggtitle("D. mel"),
                  theme(plot.title = element_text(face = "italic", size = 60, hjust = 0.5, vjust = 1, colour = COLOR[i])))
  }
  g <- cowplot::plot_grid(g_one, NULL, g_all,
                          nrow = 1,
                          rel_widths = c(10, 1, 10))
  
  gl <- c(gl, list(g))
}


Fig4AB <- cowplot::plot_grid(gl[[1]], NULL, gl[[2]],
                             nrow = 3,
                             rel_heights = c(10, 1, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) 

ggsave("Fig4/4AB.png", Fig4AB, width = 20, height = 20, dpi = 100)


#4CD Nearby populations by state -----------------------------------------------

numNN = readr::read_csv("Nearest_neighbor_population.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Calocasiae(IR)", Definition == 4)

fit = readr::read_csv("Nearest_neighbor_population_lambda.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Calocasiae(IR)", Definition == 4)

pvalue = c(2.83*10^-7, 3.66*10^-4)
Type_name = c("Calocasiae", "Dmelanogaster")
gl = list()
for (i in 1:length(Type_name)) {
  
  data_active = numNN %>%
    filter(Type == Type_name[i]) %>%
    filter(Status == "active")
  
  data_rest = numNN %>%
    filter(Type == Type_name[i]) %>%
    filter(Status == "rest")
  
  COLOR = c("deeppink", "deepskyblue")
  g <- ggplot() + 
    geom_histogram(data = data_active,
                   aes(x = numNN, y = stat(count / sum(count))),
                   position = "identity", binwidth = 1, fill = COLOR[1], alpha = 0.3) +
    geom_histogram(data = data_rest,
                   aes(x = numNN, y = stat(count / sum(count))),
                   position = "identity", binwidth = 1, fill = COLOR[2], alpha = 0.3) +
    annotate("text", x = 13, y = 0.94, label = "Walking", colour = "black", size = 15, hjust = 0, vjust = 0) +
    annotate("segment", x = 11.8, xend = 11.8, y = 0.93, yend = 1.0, colour = COLOR[1], alpha = 0.3, size = 17) +
    annotate("text", x = 13, y = 0.84, label = "Stopping", colour = "black", size = 15, hjust = 0, vjust = 0) +
    annotate("segment", x = 11.8, xend = 11.8, y = 0.83, yend = 0.9, colour = COLOR[2], alpha = 0.3, size = 17) +
    theme_bw() +
    xlab("Number of neighbours") +
    ylab("Probability mass function") +
    scale_x_continuous(limits = c(-0.5, 20.5), breaks = seq(0, 20, by = 5)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    Setting
  g
  
  data_fit = fit %>%
    filter(Type == Type_name[i])
  
  pv <- paste("p", sprintf('%.2e', pvalue[i]), sep = " = ")
  
  LABEL = c("W", "S")
  gs  <- ggplot(data_fit, aes(x = factor(Status), y = Lambda)) + 
    geom_boxplot(aes(fill = Status), width =.5, outlier.colour = NA, size = 1.5)
  gs
  
  dat <- ggplot_build(gs)$data[[1]]
  
  p <- gs + geom_segment(data = dat, aes(x = xmin, xend = xmax,
                                            y = middle, yend = middle), colour = "white", size = 2) +
    geom_dotplot(binaxis = "y", binwidth = 0.1, stackdir = "center", fill = "black", dotsize = 1, alpha = 1) +
    scale_x_discrete(labels = LABEL) +
    scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) +
    theme_cowplot() +
    ylab("ƒÉ") +
    scale_fill_manual(values = COLOR) +
    guides(fill = FALSE) +
    Setting +
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(vjust = 0)) +
    geom_signif(
      y_position = 4,
      xmin = 1.0,
      xmax = 2.0,
      annotation = "",
      tip_length = 0,
      size = 2) +
    annotate("text", 
             x = 1.5,
             y = 4.5,
             label = pv,
             vjust = 0.5,
             size = 10)
  p
  
  gg <- g + cowplot::draw_plot(p, x = 5, y = 0, width = 16, height = 0.8)
  
  gl <- c(gl, list(gg))
}

Fig4CD <- cowplot::plot_grid(gl[[1]], NULL, gl[[2]], nrow = 3, align = "hv", rel_heights = c(10, 1, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) 
Fig4CD

ggsave("Fig4/4CD.png", Fig4CD, width = 10, height = 20, dpi = 100)


#Overall -----------------------------------------------------------------------

Fig4 <- cowplot::plot_grid(NULL, Fig4AB, NULL, Fig4CD, 
                               nrow = 1, 
                               rel_widths = c(1.5, 20, 2.5, 10))

ggsave("Fig4/Fig4.png", Fig4, width = 30, height = 20, dpi = 100)