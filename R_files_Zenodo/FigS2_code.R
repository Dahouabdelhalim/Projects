source ("Library.R")

Setting <- list(theme(axis.title.x = element_text(size = 50, vjust = -0.5)), 
                 theme(axis.title.y = element_text(size = 50, vjust = 2)), 
                 theme(axis.text.x = element_text(size = 40, colour = "black", vjust = 0.5)), 
                 theme(axis.text.y = element_text(size = 40, colour = "black")))


#S2 Nearby populations by state ------------------------------------------------

numNN = readr::read_csv("Nearest_neighbor_population.csv") %>%
  dplyr::select(-1)

fit = readr::read_csv("Nearest_neighbor_population_lambda.csv") %>%
  dplyr::select(-1)

pvalue = c(2.83*10^-7, 
           2.39*10^-7, 
           4.05*10^-7, 
           3.66*10^-4, 
           1.47*10^-4, 
           2.71*10^-4, 
           1.06*10^-3, 
           8.73*10^-4, 
           6.71*10^-4)

Type_name = c("Calocasiae", "Dmelanogaster", "Calocasiae(IR)")
Def = c(4, 5, 6)
k = 0
for (i in 1:length(Type_name)) {
  gl = list()
  for(j in 1:length(Def)){
    
    data_active = numNN %>%
      filter(Type == Type_name[i], 
             Definition == Def[j], 
             Status == "active")
    
    data_rest = numNN %>%
      filter(Type == Type_name[i], 
             Definition == Def[j], 
             Status == "rest")
    
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
      ggtitle(paste("Def : ",
                    bquote(. (Def[j])),
                    "th body length",
                    sep = "")) +
      theme(plot.title = element_text(size = 50, hjust = 0.5, vjust = 1)) +
      Setting
    g
    
    data_fit = fit %>%
      filter(Type == Type_name[i],
             Definition == Def[j])
    
    k = k + 1
    pv <- paste("p", sprintf('%.2e', pvalue[k]), sep = " = ")
    
    COLOR = c("deeppink", "deepskyblue")
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
  if(i == 1){
    gl_Calo <- gl
  }
  if(i == 2){
    gl_Dmel <- gl
  }
  if(i == 3){
    gl_IR <- gl
  }
}

gl[[1]]

g_Calo <- cowplot::plot_grid(gl_Calo[[1]], NULL, gl_Calo[[2]], NULL, gl_Calo[[3]],
                             nrow = 1, 
                             rel_widths = c(10, 2, 10, 2, 10)) +
  ggtitle("C. alocasiae") +
  theme(plot.title = element_text(face = "italic", colour = "deepskyblue4", size = 60, hjust = 0.5, vjust = 3))


g_Dmel <- cowplot::plot_grid(gl_Dmel[[1]], NULL, gl_Dmel[[2]], NULL, gl_Dmel[[3]], 
                             nrow = 1, 
                             rel_widths = c(10, 2, 10, 2, 10)) +
  ggtitle("D. melanogaster") +
  theme(plot.title = element_text(face = "italic", colour = "orange2", size = 60, hjust = 0.5, vjust = 3))


g_IR <- cowplot::plot_grid(gl_IR[[1]], NULL, gl_IR[[2]], NULL, gl_IR[[3]], 
                           nrow = 1, 
                           rel_widths = c(10, 2, 10, 2, 10)) +
  ggtitle(expression(paste(italic("C. alocasiae"), 
                           " (IR) ",
                           sep=""))) +
  theme(plot.title = element_text(size = 60, colour = "red4", hjust = 0.5, vjust = 3))

FigS2 <- cowplot::plot_grid(g_Calo, NULL, g_Dmel, NULL, g_IR,
                            nrow = 5,
                            rel_heights = c(10, 2, 10, 2, 10)) +
  theme(plot.margin = unit(c(2, 1, 1, 1), "lines")) 

ggsave("FigS2/FigS2.png", FigS2, width = 30, height = 35, dpi = 100)