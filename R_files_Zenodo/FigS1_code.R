source ("Library.R")

Setting <- list(theme(axis.title.x = element_text(size = 50, vjust = -0.5)), 
                 theme(axis.title.y = element_text(size = 50, vjust = 2)), 
                 theme(axis.text.x = element_text(size = 40, colour = "black", vjust = 0.5)), 
                 theme(axis.text.y = element_text(size = 40, colour = "black")))


#S1B Definition of group --------------------------------------------------------

sample1 = readr::read_csv("Group_judgment.csv") %>%
  dplyr::select(-1) %>%
  filter(Type == "Calocasiae", Movie_number == 2, Time == 929)

NAME = c("Fly_number", "X", "Y")
adjust_Calo = 127/823

sample = readr::read_csv("Calo_f_b.csv", col_names = F) %>%
  slice(28160)
s0 = NULL
for (i in 0:19) {
  s1 = sample %>%
    dplyr::select((1:3) + 6 * i) %>%
    rename_all(vars(NAME)) %>%
    mutate(x = X * adjust_Calo, y = Y * adjust_Calo) %>%
    dplyr::select(-c(X, Y))
  s0 = bind_rows(s0, s1)
}

Func_Dis = function(x1, y1, x2, y2)
{
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

size_Calo = readr::read_csv("Size_mean_sd.csv") %>%
  dplyr::select(C.alo_mean)
Threshold_Distance = size_Calo * 2
Def_Distance = function(Distance)
{
  if(Distance < Threshold_Distance){
    Distance <- 1
  }
  else{
    Distance <- 0
  }
}

Threshold_Population = 3
Def_Population = function(Population)
{
  if(Population < Threshold_Population){
    Population <- 0
  }
  else{
    Population <- 1
  }
}

s7 = NULL
for (i in 0:19) {
  s2 = s0 %>%
    filter(Fly_number == i)
  
  for (j in 0:19) {
    s3 = s0 %>%
      filter(Fly_number == j)
    s4 <- as.data.frame(mapply(Func_Dis, s2$x, s2$y, s3$x, s3$y))
    s5 <- as.data.frame(mapply(Def_Distance, s4[, 1]))
    colnames(s5) = j
    s2 = bind_cols(s2, s5)
  }
  s6 = s2 %>%
    gather(-c(Fly_number, x, y), key = other_fly, value = distance) %>%
    group_by(Fly_number, x, y) %>%
    summarise(Proximate_population = sum(distance)) %>%
    ungroup()
  s7 = bind_rows(s7, s6)
}
s8 <- as.data.frame(mapply(Def_Population, s7$Proximate_population))
colnames(s8) <- "Group_in"
s9 = bind_cols(s7, s8)

FigS1B <- ggplot(s9, aes(x = x, y = y)) +
  geom_point(aes(color = factor(Group_in)), size = 3) +
  coord_fixed() +
  theme_classic()
FigS1B

ggsave("FigS1/S1B.png", FigS1B, width = 10, height = 10, dpi = 100)


#S1C Accumulative distribution -------------------------------------------------

Map = readr::read_csv("Accumulative_distribution.csv") %>%
  dplyr::select(-1)

Map_random = readr::read_csv("Accumulative_distribution_random.csv") %>%
  dplyr::select(-1)

Type_name = c("Calocasiae", "Dmelanogaster", "Calocasiae(IR)")

gl = list()
gl_zoom = list()
for (i in 1:4) {
  
  if(i == 4){
    data = Map %>%
      filter(Type == "Calocasiae(IR)")
  }
  else{
    data = Map_random %>%
      filter(Type == Type_name[i])
  }
  
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
    Setting +
    theme(axis.title.y = element_text(vjust = 1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 18),
          legend.direction = "horizontal", 
          legend.key.size = grid::unit(1.0, "cm"), 
          legend.position = c(0.75, 0.86)) +
    guides(fill = guide_colorbar(title.position = "top"))
  g
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
    Setting +
    theme(axis.title.y = element_text(vjust = 1))
  gl_zoom <- c(gl_zoom, list(g_zoom))
  
}

FigS1C1 <- cowplot::plot_grid(gl[[1]], NULL, gl_zoom[[1]],
                              nrow = 1,
                              rel_widths = c(10, 0.5, 10)) +
  ggtitle(expression(paste(italic("C. alocasiae"), 
                           " [Random] ",
                           sep=""))) +
  theme(plot.title = element_text(size = 60, colour = "deepskyblue4", hjust = 0.5, vjust = 1))
FigS1C1

FigS1C2 <- cowplot::plot_grid(gl[[2]], NULL, gl_zoom[[2]],
                              nrow = 1,
                              rel_widths = c(10, 0.5, 10)) +
  ggtitle(expression(paste(italic("D. melanogaster"), 
                           " [Random] ",
                           sep=""))) +
  theme(plot.title = element_text(size = 60, colour = "orange2", hjust = 0.5, vjust = 1))
FigS1C2

FigS1C3 <- cowplot::plot_grid(gl[[3]], NULL, gl_zoom[[3]],
                              nrow = 1,
                              rel_widths = c(10, 0.5, 10)) +
  ggtitle(expression(paste(italic("C. alocasiae"), 
                           " (IR) [Random] ",
                           sep=""))) +
  theme(plot.title = element_text(size = 60, colour = "red4",hjust = 0.5, vjust = 1))
FigS1C3

FigS1C4 <- cowplot::plot_grid(gl[[4]], NULL, gl_zoom[[4]],
                              nrow = 1,
                              rel_widths = c(10, 0.5, 10)) +
  ggtitle(expression(paste(italic("C. alocasiae"), 
                           " (IR) ",
                           sep=""))) +
  theme(plot.title = element_text(size = 60, colour = "red4", hjust = 0.5, vjust = 1))
FigS1C4

FigS1C <- cowplot::plot_grid(FigS1C1, NULL, FigS1C2,
                             NULL, NULL, NULL,
                             FigS1C4, NULL, FigS1C3,
                             nrow = 3,
                             rel_widths = c(10, 1, 10),
                             rel_heights = c(10, 1, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) 

ggsave("FigS1/S1C.png", FigS1C, width = 30, height = 16, dpi = 100)


#S1D Inter individual distance -------------------------------------------------

Distance = readr::read_csv("Inter_individual_distance.csv") %>%
  dplyr::select(-1) %>%
  filter(Type == "Calocasiae(IR)")

Distance_random = readr::read_csv("Inter_individual_distance_random.csv") %>%
  dplyr::select(-1) %>%
  filter(Type == "Calocasiae(IR)") %>%
  rename(Movie_number = Trial_number)

size_IR = readr::read_csv("Size_mean_sd.csv") %>%
  dplyr::select(C.alo_mean)

data_all = bind_rows(data.frame(Data = as.character("Original"), Distance),
                     data.frame(Data = as.character("Random"), Distance_random))

Data_name = c("Original", "Random")
COLOR = c("red4")
gl = list()
for (j in 1:2) {
  
  data = data_all %>%
    filter(Data == Data_name[j])
  
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
  
  Range = as.numeric(size_IR * 2)
  if(j == 2){
    k = "black"
  }else{
    k = COLOR[1]
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

FigS1D <- cowplot::plot_grid(NULL, gl[[1]],
                             NULL, NULL, 
                             NULL, gl[[2]],
                             nrow = 3,
                             rel_widths = c(2, 20),
                             rel_heights = c(8, -1.4, 8)) +
  ggtitle(expression(paste(italic("C. alocasiae"), 
                           " (IR) ",
                           sep=""))) +
  theme(plot.title = element_text(size = 60, colour = "black", hjust = 0.5, vjust = 1)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
FigS1D

ggsave("FigS1/S1D.png", FigS1D, width = 10, height = 15, dpi = 100)


#S1E Speed ---------------------------------------------------------------------

Speed = readr::read_csv("Speed.csv") %>%
  dplyr::select(-1)

Type_name = c("Calocasiae", "Dmelanogaster", "Calocasiae(IR)")
Fly_name = c("C. alo", "D. mel")

COLOR = c("deepskyblue4", "orange2", "red4")
THRESHOLD = c(1.0, 1.0, 2.0)
gl = list()
for (i in 1:3) {
  
  data = Speed %>%
    filter(Type == Type_name[i])
  
  g <- ggplot(data, aes(x = Speed, y = stat(count / sum(count))*100)) +
    geom_histogram(binwidth = 0.5, fill = COLOR[[i]], boundary = 0.0) +
    geom_vline(xintercept = THRESHOLD[i], linetype = "dashed", size = 1.5) +
    theme_bw() +
    xlab("Speed (mm/s)") +
    ylab("Percentage (%)") +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 50)) +
    scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5)) +
    Setting
  
  if(i != 3){
    gg <- g +
      ggtitle(Fly_name[i]) +
      theme(plot.title = element_text(face = "italic", size = 60, hjust = 0.5, vjust = 1))
  }
  else{
    gg <- g +
      ggtitle(expression(paste(italic("C. alo"), 
                               " (IR) ",
                               sep=""))) +
      theme(plot.title = element_text(size = 60, hjust = 0.5, vjust = 1))
  }
  
  gl <- c(gl, list(gg))
  
}

FigS1E <- cowplot::plot_grid(gl[[1]], NULL, gl[[2]], NULL, gl[[3]],
                             nrow = 1,
                             rel_widths = c(10, 1.5, 10, 1.5, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) 
FigS1E

ggsave("FigS1/S1E.png", FigS1E, width = 20, height = 6, dpi = 100)


#S1F Probability of joining by group size --------------------------------------

Inout_size = readr::read_csv("Group_inout_size.csv") %>%
  dplyr::select(-1) %>%
  filter(Type == "Calocasiae(IR)") %>%
  group_by(Type, Group_size) %>%
  summarise(Ratio = mean(Ratio), Count = mean(Number_of_events_per_group_size))

event_num = 11
Name = c("C. alocasiae (IR)")
Type = c(rep(Name, 20 - event_num))
Group_size = c(11:19)
Ratio = c(rep(0, 20 - event_num))
Count = c(rep(0, 20 - event_num))
df = data.frame(Type, Group_size, Ratio, Count)

data = bind_rows(Inout_size, df) %>%
  mutate(Count = paste0("(", Count, ")")) %>%
  arrange(Type, Group_size)

FigS1F <- ggplot(data) + 
  geom_bar(aes(x = Group_size, y = Ratio), stat = "identity", width = 0.3, fill = "red4") +
  geom_text(aes(x = Group_size, y = 105, label = Count), size = 10) +
  theme_bw() +
  xlab("Number of individuals in a group") +
  ylab("Average percentage \\n joining a group (%)") +
  scale_x_continuous(limits = c(1.7, 19.3), breaks = seq(2, 19, by = 1)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, by = 25)) +
  ggtitle(expression(paste(italic("C. alocasiae"), 
                           " (IR) ",
                           sep=""))) +
  theme(plot.title = element_text(size = 60, hjust = 0.5, vjust = 1)) +
  Setting +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"),
        axis.title.y = element_text(lineheight = 1.2)) 
FigS1F

ggsave("FigS1/S1F.png", FigS1F, width = 20, height = 8, dpi = 100)


#Overall -----------------------------------------------------------------------

FigS1EF <- cowplot::plot_grid(FigS1E,
                              NULL,
                              FigS1F,
                              nrow = 3,
                              rel_heights = c(6, 1, 8))


FigS1DEF <- cowplot::plot_grid(FigS1D, NULL, FigS1EF,
                               nrow = 1,
                               rel_widths = c(10, 1.5, 20))


FigS1CDEF <- cowplot::plot_grid(FigS1C,
                               NULL,
                               FigS1DEF,
                               nrow = 3,
                               rel_heights = c(16, 2, 15))
FigS1CDEF

ggsave("FigS1/S1CDEF.png", FigS1CDEF, width = 30, height = 33, dpi = 100)