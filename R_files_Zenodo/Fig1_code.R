source ("Library.R")

Setting <- list(theme(axis.title.x = element_text(size = 50, vjust = -0.5)),
                 theme(axis.title.y = element_text(size = 50, vjust = 2)),
                 theme(axis.text.x = element_text(size = 40, colour = "black", vjust = 0.5)),
                 theme(axis.text.y = element_text(size = 40, colour = "black")))


#1C-1 Photo after 20 minutes ---------------------------------------------------

photo_Calo <- readPNG("Calo_20min.png")
image_Calo <- ggdraw() + draw_image(photo_Calo, scale = 1)

photo_Dmel <- readPNG("Dmel_20min.png")
image_Dmel <- ggdraw() + draw_image(photo_Dmel, scale = 1)


#1C-2 Changes in group population ----------------------------------------------

Group = readr::read_csv("Group_judgment_correction.csv") %>%
  dplyr::select(-1) %>%
  group_by(Type, Movie_number, Time) %>%
  summarise(Group_Size = sum(Group_in)) %>%
  mutate(Min = Time/60) %>%
  ungroup()

interval = 10
Time <- c(seq(0, 1200, interval))

Group_10s = NULL #Cut every 10 seconds
for (t in Time) {
  df = Group %>%
    filter(Time == t)
  Group_10s = bind_rows(Group_10s, df)
}

#Median
Group_median = Group_10s %>%
  group_by(Type, Min) %>%
  summarise(Group_Size = median(Group_Size))

Type_name = c("Calocasiae", "Dmelanogaster")
COLOR = c("deepskyblue4", "orange2")
gl = list()
for (i in 1:2) {
  
  data = Group %>%
    filter(Type == Type_name[i])
  
  data_median = Group_median %>%
    filter(Type == Type_name[i])
  
  g <- ggplot() +
    geom_line(data = data, aes(x = Min, y = Group_Size, group = factor(Movie_number)), size = 1.0, alpha = 0.2) +
    geom_line(data = data_median, aes(x = Min, y = Group_Size), size = 2.0, colour = "black") +
    theme_bw() +
    xlab("Time (min)") +
    ylab("# Grouped flies") +
    scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5)) +
    scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5)) +
    Setting +
    theme(axis.title.y = element_text(hjust = 1))
  gl <- c(gl, list(g))
}


#Overall -----------------------------------------------------------------------

Fig1C <- cowplot::plot_grid(NULL, image_Calo, NULL, gl[[1]],
                            NULL, NULL, NULL, NULL,
                            NULL, image_Dmel, NULL, gl[[2]],
                            nrow = 3, 
                            rel_widths = c(8, 6, 1.5, 15),
                            rel_heights = c(10, 1, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) 

ggsave("Fig1/1C.png", Fig1C, width = 30, height = 12.5, dpi = 100)