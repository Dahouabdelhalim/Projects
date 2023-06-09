source ("Library.R")

Setting <- list(theme(axis.title.x = element_text(size = 50, vjust = -0.5)),
                 theme(axis.title.y = element_text(size = 50, vjust = 2)),
                 theme(axis.text.x = element_text(size = 40, colour = "black", vjust = 0.5)),
                 theme(axis.text.y = element_text(size = 40, colour = "black")))


#5A Change in speed and distance over time -------------------------------------

frame_adjust = 30/128

Parameters = readr::read_csv("Average_params.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Calocasiae(IR)") %>%
  mutate(across(Speed, ~ .x * frame_adjust)) #convert to every second

data = Parameters %>%
  filter(Type == "Calocasiae",
         Movie_number == 3) %>%
  dplyr::select(-c(Type, Movie_number)) %>%
  pivot_longer(-Time, names_to = "param", values_to = "val")

data_distance = data %>%
  filter(param == "Distance")

data_speed = data %>%
  filter(param == "Speed")

Time_point = c(0, 5.42, 8.08, 10.43)

g1 <- ggplot(data_distance, aes(x = Time, y = val, col = param)) +
  geom_vline(xintercept = Time_point, linetype = "dashed", size = 1.5, alpha = 0.5) +
  geom_line(size = 1.5) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 20),
                     breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(limits = c(0, 52),
                     breaks = seq(0, 50, by = 50)) +
  scale_color_manual(values = "green3", 
                     name = "") +
  labs(x = "Time (min)", y = "Distance\\n(mm)") +
  guides(col = FALSE) +
  Setting +
  theme(axis.text.x = element_text(colour = "white"),
        axis.title.x = element_text(colour = "white"),
        axis.text.y = element_blank(),
        axis.title.y = element_text(vjust = 8))

g2 <- ggplot(data_speed, aes(x = Time, y = val, col = param)) +
  geom_vline(xintercept = Time_point, linetype = "dashed", size = 1.5, alpha = 0.5) +
  geom_line(size = 1.5) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 20),
                     breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(limits = c(0, 0.8),
                     breaks = seq(0, 0.8, by = 0.8)) +
  scale_color_manual(values = "magenta3",
                     name = "") +
  labs(x = "Time (min)",
       y = "Speed  \\n(mm/sec)  ") +
  guides(col = FALSE) +
  Setting +
  theme(axis.title.y = element_text(vjust = 8),
        axis.text.y = element_blank())

Fig5A <- cowplot::plot_grid(g1, NULL, g2,
                          nrow = 3,
                          rel_heights = c(10, -2.9, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 6), "lines"))


#5B Photographs showing the passage of time ------------------------------------

photo1 <- readPNG("Calocasiae_f_c_0min0sec.png")
image1 <- ggdraw() + draw_image(photo1, scale = 1)

photo2 <- readPNG("Calocasiae_f_c_5min25sec.png")
image2 <- ggdraw() + draw_image(photo2, scale = 1)

photo3 <- readPNG("Calocasiae_f_c_8min05sec.png")
image3 <- ggdraw() + draw_image(photo3, scale = 1)

photo4 <- readPNG("Calocasiae_f_c_10min26sec.png")
image4 <- ggdraw() + draw_image(photo4, scale = 1)

Fig5B <- cowplot::plot_grid(image1, NULL, image2, NULL, image3, NULL, image4,
                             NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                             nrow = 2,
                             rel_widths = c(10, 0.5, 10, 0.5, 10, 0.5, 10),
                             rel_heights = c(10, 1)) +
  theme(plot.margin = unit(c(0, 1, 0, 0), "lines")) 


#Overall -----------------------------------------------------------------------

Fig5AB <- cowplot::plot_grid(Fig5A, NULL, Fig5B,
                            nrow = 1,
                            rel_widths = c(18, 0.5, 40))

ggsave("Fig5/Fig5.png", Fig5AB, width = 30, height = 8, dpi = 100)