source ("Library.R")

Setting1 <- list(theme(axis.title.x = element_text(size = 50, vjust = -0.5)),
                 theme(axis.title.y = element_text(size = 50, vjust = 2)),
                 theme(axis.text.x = element_text(size = 40, colour = "black", vjust = 0.5)),
                 theme(axis.text.y = element_text(size = 40, colour = "black")))

Setting2 <- list(theme(axis.title.x = element_blank()),
                 theme(axis.title.y = element_text(size = 50, vjust = 2)),
                 theme(axis.text.x = element_text(size = 40, colour = "black", vjust = -0.5)),
                 theme(axis.text.y = element_text(size = 40, colour = "black")))


#6A Changes in group population ------------------------------------------------

Group = readr::read_csv("Group_judgment_correction.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Dmelanogaster") %>%
  group_by(Type, Movie_number, Time) %>%
  summarise(Group_Size = sum(Group_in)) %>%
  mutate(Min = Time/60) %>%
  ungroup()

interval = 10
Time <- c(seq(0, 1200, interval))

Group_10s = NULL
for (t in Time) {
  df = Group %>%
    filter(Time == t)
  Group_10s = bind_rows(Group_10s, df)
}

Group_median = Group_10s %>%
  group_by(Type, Min) %>%
  summarise(Group_Size = median(Group_Size))

Type_name = c("Calocasiae", "Calocasiae(IR)")
LABEL = c("White", "IR")
COLOR = c("deepskyblue4", "red4")
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
    ggtitle(LABEL[i]) +
    theme(plot.title = element_text(size = 60, hjust = 0.5, colour = COLOR[i])) +
    Setting1
  gl <- c(gl, list(g))
}

Fig6A <- cowplot::plot_grid(gl[[1]], NULL, gl[[2]], nrow = 3, align = "hv", rel_heights = c(10, 1, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) 
Fig6A

ggsave("Fig6/6A.png", Fig6A, width = 30, height = 15, dpi = 100)


#6B Probability of joining the group -------------------------------------------

data = readr::read_csv("Group_inout.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Dmelanogaster")

Sample_size = data %>%
  group_by(Type) %>%
  tally() %>%
  mutate(Count = paste0("(", n, ")"))

pv <- paste("p", 1.83*10^-36, sep = " = ")

LABEL = c("White", "IR")
COLOR = c("deepskyblue4", "red4")

g  <- ggplot(data, aes(x = factor(Type), y = Ratio)) + 
  geom_boxplot(aes(fill = Type), width =.5, outlier.colour = NA, size = 1.5) 

dat <- ggplot_build(g)$data[[1]]

Fig6B <- g + geom_segment(data = dat, aes(x = xmin, xend = xmax,
                                          y = middle, yend = middle), colour = "white", size = 2) +
  geom_dotplot(binaxis = "y", binwidth = 1, stackdir = "center", fill = "black", dotsize = 0.7, alpha = 1) +
  geom_text(data = Sample_size, aes(x = factor(Type), y = 105, label = Count), size = 10) +
  scale_fill_manual(values = COLOR) +
  theme_bw() +
  ylab("Percentage joining a group (%)") +
  scale_x_discrete(labels = LABEL) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, by = 25)) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 20)) +
  guides(fill = FALSE) +
  Setting2 +
  geom_signif(
    y_position = 113,
    xmin = 1.0,
    xmax = 2.0,
    annotation = "",
    tip_length = 0,
    size = 2) +
  annotate("text",
           x = 1.5,
           y = 118,
           label = pv,
           vjust = 0.5,
           size = 10)
Fig6B

ggsave("Fig6/6B.png", Fig6B, width = 20, height = 15, dpi = 100)


#6C Length of stay in the group ------------------------------------------------

data = readr::read_csv("Length_of_stay.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Dmelanogaster")

Sample_size = data %>%
  group_by(Type) %>%
  tally() %>%
  mutate(Count = paste0("(", n, ")"))

pv <- paste("p", 9.14*10^-37, sep = " = ")

tateziku = c(0, 1, 5, 20)
LABEL = c("White", "IR")
COLOR = c("deepskyblue4", "red4")
ylim = 50
g  <- ggplot(data, aes(x = factor(Type), y = Mean_duration)) + 
  geom_boxplot(aes(fill = Type), width =.5, outlier.colour = NA, size = 1.5) 

dat <- ggplot_build(g)$data[[1]]

Fig6C <- g + geom_segment(data = dat, aes(x = xmin, xend = xmax,
                                          y = middle, yend = middle), colour = "white", size = 2) +
  geom_dotplot(binaxis = "y", binwidth = 1/60, stackdir = "center", fill = "black", dotsize = 1, alpha = 1) +
  geom_text(data = Sample_size, aes(x = factor(Type), y = 25, label = Count), size = 10) +
  scale_y_log10(limits = c(0.16, ylim), breaks = tateziku) +
  scale_x_discrete(labels = LABEL) +
  theme_bw() +
  ylab("Average grouping duration (min) ") +
  scale_fill_manual(values = COLOR) +
  guides(fill = FALSE) +
  Setting2 +
  geom_signif(
    y_position = 1.57,
    xmin = 1.0,
    xmax = 2.0,
    annotation = "",
    tip_length = 0,
    size = 2) +
  annotate("text",
           x = 1.5,
           y = 45,
           label = pv,
           vjust = 0.5,
           size = 10)
Fig6C

ggsave("Fig6/6C.png", Fig6C, width = 10, height = 15, dpi = 100)


#6D Correlation between speed and distance -------------------------------------

frame_adjust = 30/128

Parameters = readr::read_csv("Average_params.csv") %>%
  dplyr::select(-1) %>%
  filter(Type == "Calocasiae(IR)") %>%
  mutate(across(Speed, ~ .x * frame_adjust)) #convert to every second

dis_max = 80
vel_max = 1.5
com_max = 1
xscale = 3

data = Parameters %>%
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
  Setting1

#One movie
data_one = data %>%
  filter(Movie_number == 3)

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
  Setting1


Fig6D <- cowplot::plot_grid(g_one, NULL, g_all,
                            nrow = 1,
                            rel_widths = c(10, 1, 10))

ggsave("Fig6/6D.png", Fig6D, width = 20, height = 10, dpi = 100)


#6E Nearby populations by state ------------------------------------------------

numNN = readr::read_csv("Nearest_neighbor_population.csv") %>%
  dplyr::select(-1) %>%
  filter(Type == "Calocasiae(IR)", Definition == 4)

fit = readr::read_csv("Nearest_neighbor_population_lambda.csv") %>%
  dplyr::select(-1) %>%
  filter(Type == "Calocasiae(IR)", Definition == 4)

data_active = numNN %>%
  filter(Status == "active")

data_rest = numNN %>%
  filter(Status == "rest")

g <- ggplot() + 
  geom_histogram(data = data_active,
                 aes(x = numNN, y = stat(count / sum(count))),
                 position = "identity", binwidth = 1, fill = "deeppink", alpha = 0.3) +
  geom_histogram(data = data_rest,
                 aes(x = numNN, y = stat(count / sum(count))),
                 position = "identity", binwidth = 1, fill = "deepskyblue", alpha = 0.3) +
  annotate("text", x = 13, y = 0.94, label = "Walking", colour = "black", size = 15, hjust = 0, vjust = 0) +
  annotate("segment", x = 11.8, xend = 11.8, y = 0.93, yend = 1.0, colour = "deeppink", alpha = 0.3, size = 17) +
  annotate("text", x = 13, y = 0.84, label = "Stopping", colour = "black", size = 15, hjust = 0, vjust = 0) +
  annotate("segment", x = 11.8, xend = 11.8, y = 0.83, yend = 0.9, colour = "deepskyblue", alpha = 0.3, size = 17) +
  theme_bw() +
  xlab("Number of neighbours") +
  ylab("Probability mass function") +
  scale_x_continuous(limits = c(-0.5, 20.5), breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  Setting1
g

data_fit = fit

pv <- paste("p", sprintf('%.2e', 1.06*10^-3), sep = " = ")

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
  Setting1 +
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

Fig6E <- g + cowplot::draw_plot(p, x = 5, y = 0, width = 16, height = 0.8)

ggsave("Fig6/6E.png", Fig6E, width = 10, height = 10, dpi = 100)


#Overall -----------------------------------------------------------------------

Fig6BC <- cowplot::plot_grid(Fig6B, NULL, Fig6C,
                             nrow = 1, 
                             rel_widths = c(20, 2, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
Fig6BC

Fig6DE<- cowplot::plot_grid(Fig6D, NULL, Fig6E,
                            nrow = 1,
                            rel_widths = c(20, 2.5, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) 

Fig6 <- cowplot::plot_grid(Fig6A, 
                              NULL,
                              Fig6BC,
                              NULL,
                              Fig6DE,
                              nrow = 5,
                              rel_heights = c(15, 2, 15, 2, 10))

ggsave("Fig6/Fig6.png", Fig6, width = 30, height = 45, dpi = 100)