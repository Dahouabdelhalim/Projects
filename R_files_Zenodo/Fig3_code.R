source ("Library.R")

Setting1 <- list(theme(axis.title.x = element_text(size = 50, vjust = -0.5)),
                 theme(axis.title.y = element_text(size = 50, vjust = 2)),
                 theme(axis.text.x = element_text(size = 40, colour = "black", vjust = 0.5)),
                 theme(axis.text.y = element_text(size = 40, colour = "black")))

Setting2 <- list(theme(axis.title.x = element_blank()),
                 theme(axis.title.y = element_text(size = 50, vjust = 2)),
                 theme(axis.text.x = element_text(size = 40, colour = "black", vjust = -0.5, face = "italic")),
                 theme(axis.text.y = element_text(size = 40, colour = "black")))


#3A Probability of joining the group -------------------------------------------

data = readr::read_csv("Group_inout.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Calocasiae(IR)")

Sample_size = data %>%
  group_by(Type) %>%
  tally() %>%
  mutate(Count = paste0("(", n, ")"))

pv <- paste("p", 1.64*10^-33, sep = " = ")

LABEL = c("C. alo", "D. mel")
COLOR = c("deepskyblue4", "orange2")

g  <- ggplot(data, aes(x = factor(Type), y = Ratio)) + 
  geom_boxplot(aes(fill = Type), width =.5, outlier.colour = NA, size = 1.5) 

dat <- ggplot_build(g)$data[[1]]

Fig3A <- g + geom_segment(data = dat, aes(x = xmin, xend = xmax,
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
Fig3A

ggsave("Fig3/3A.png", Fig3A, width = 20, height = 15, dpi = 100)


#3B Length of stay in the group ------------------------------------------------

data = readr::read_csv("Length_of_stay.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Calocasiae(IR)")

Sample_size = data %>%
  group_by(Type) %>%
  tally() %>%
  mutate(Count = paste0("(", n, ")"))

pv <- paste("p", 9.72*10^-224, sep = " = ")

tateziku = c(0, 1, 5, 20)
COLOR = c("deepskyblue4", "orange2")
LABEL = c("C. alo", "D. mel")
ylim = 50
g  <- ggplot(data, aes(x = factor(Type), y = Mean_duration)) + 
  geom_boxplot(aes(fill = Type), width =.5, outlier.colour = NA, size = 1.5) 

dat <- ggplot_build(g)$data[[1]]

Fig3B <- g + geom_segment(data = dat, aes(x = xmin, xend = xmax,
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
Fig3B

ggsave("Fig3/3B1.png", Fig3B, width = 10, height = 15, dpi = 100)


#3C Probability of joining by group size ---------------------------------------
Inout_size = readr::read_csv("Group_inout_size.csv") %>%
  dplyr::select(-1) %>%
  filter(Type != "Calocasiae(IR)") %>%
  group_by(Type, Group_size) %>%
  summarise(Ratio = mean(Ratio), Count = mean(Number_of_events_per_group_size))

event_num = 15
Name = c("Dmelanogaster")
Type = c(rep(Name, 20 - event_num))
Group_size = c(14:15, 17:19)
Ratio = c(rep(0, 20 - event_num))
Count = c(rep(0, 20 - event_num))
df = data.frame(Type, Group_size, Ratio, Count)

Inout_size_all = bind_rows(Inout_size, df) %>%
  mutate(Count = paste0("(", Count, ")")) %>%
  arrange(Type, Group_size)

Type_name = c("Calocasiae", "Dmelanogaster")
COLOR = c("deepskyblue4", "orange2")
LABEL = c("C. alocasiae", "D. melanogaster")
gl = list()
for (i in 1:2) {
  data = Inout_size_all %>%
    filter(Type == Type_name[i])
  
  g <- ggplot(data) + 
    geom_bar(aes(x = Group_size, y = Ratio), stat = "identity", width = 0.3, fill = COLOR[i]) +
    geom_text(aes(x = Group_size, y = 105, label = Count), size = 10) +
    theme_bw() +
    xlab("Number of individuals in a group") +
    ylab("Average percentage \\n joining a group (%)") +
    scale_x_continuous(limits = c(1.7, 19.3), breaks = seq(2, 19, by = 1)) +
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, by = 25)) +
    ggtitle(LABEL[i]) +
    Setting1 +
    theme(plot.title = element_text(face = "italic", size = 60, hjust = 0.5, vjust = 2),
          axis.title.x = element_text(vjust = 0),
          axis.title.y = element_text(lineheight = 1.2))
  
  gl <- c(gl, list(g))
  
}

Fig3C <- cowplot::plot_grid(gl[[1]],
                            NULL,
                            gl[[2]],
                            nrow = 3,
                            rel_heights = c(10, 1.5, 10))
Fig3C
ggsave("Fig3/3C.png", Fig3C, width = 30, height = 20, dpi = 100)


#Overall -----------------------------------------------------------------------

Fig3AB <- cowplot::plot_grid(Fig3A, NULL, Fig3B,
                             nrow = 1, 
                             rel_widths = c(20, 2, 10)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
Fig3AB

Fig3 <- cowplot::plot_grid(Fig3AB,
                              NULL,
                              Fig3C,
                              nrow = 3,
                              rel_heights = c(15, 2, 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
Fig3ABC

ggsave("Fig3/Fig3.png", Fig3, width = 30, height = 37, dpi = 100)