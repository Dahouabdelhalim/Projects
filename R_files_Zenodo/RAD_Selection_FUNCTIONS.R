# Plot theme --------------------------------------------------------------

theme_claudia <- function(addlegend = FALSE)
{
  theme_bw()+
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          panel.spacing = unit(0.8, "lines"),
          axis.title.x = element_text(size = 10, colour = "black", vjust = 0),
          axis.title.y = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(size = 8, colour = "black"),
          axis.text.y = element_text(size = 8, colour = "black"),
          axis.ticks = element_line(colour = "black", size = 0.5),
          axis.line = element_blank(),
          legend.title = element_blank(),
          legend.text = if(addlegend == TRUE) element_text(size = 9) else element_blank(),
          legend.background = if(addlegend == TRUE) element_rect(fill = "transparent") else element_blank(),
          legend.key.size = if(addlegend == TRUE) unit(0.7, "line") else element_blank())
}



# Plot colours ------------------------------------------------------------

yellow <- "#ffd32a"
dark_blue <- "#2980b9"
dark_green <- "#0A6B28"
purple <- "#833471"
fire_orange <- "#EE5A24"
sky_blue <- "#87CEEB"
grey <- "#C0C0C0"
light_green <- "#2ecc71"
pink <- "#ffc0cb"
mid_green <- "#008000"
navy <- "#123456"
pale <- "#eee8aa"
dark_red <- "#8b0000"

