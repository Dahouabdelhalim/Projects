
#>> Installing and loading required packages ----------------------------------------------------

install.packages("pacman")
library(pacman)
p_load(ggplot2, readxl, plyr, reshape2, viridis, effsize, ggsci)

#>> Loading and pre-processing the EV dataset ----------------------------------------------------

download.file("https://datadryad.org/stash/downloads/file_stream/763633", "EV.xlsx", mode = "wb")
EV <-  as.data.frame(read_excel("EV.xlsx"))
EV[,2:4] <- lapply(EV[,2:4], as.factor)
EV <- EV[,c(1:7,11:12,8:10,13:14)]

EV.2C <- melt(EV, id.vars=names(EV)[2:4], measure.vars = names(EV)[5:12],
              variable.name = "Gene", value.name = "Expr")
EV.2C$Type <- EV.2C$Gene
levels(EV.2C$Type) <- rep(c("Protein coding","lncRNA"), c(5,3))
levels(EV.2C$Gene)[8] <- "AGAP2-\\nAS1"

#>> Rendering Figure 2B ---------------------------------------------------------

SIZE = 10

ggplot(ddply(melt(EV, id.vars=names(EV)[2:4], measure.vars = names(EV)[5:12],
                  variable.name = "Gene", value.name = "Expr"),
                  .(Gene, X_study, cat2), plyr::summarize, Expr = median(Expr)), aes(x = X_study, y = cat2)) +
  
  geom_point(shape=16, aes(colour = Expr, size = Expr^2)) +
  
  labs(x = "",y = "", title = "") +
  facet_wrap(.~Gene, nrow = 1, scales = "fixed") +
  
  scale_x_discrete(labels = c("non-tumour","tumour"), expand = c(0.9,0)) +
  scale_color_viridis_c(begin = 0, end = 1, direction = -1, option = "B", alpha = 0.9, na.value = "white",
                        limits = c(0,18.7), breaks=c(0,11.6,18.7), labels=c(0,11.6,18.7)) +
  scale_size(name="", guide="legend", range = c(1,5.5),
             limits = c(0,18.7)^2, breaks=c(0,11.6,18.7)^2, labels=c(0,12,19))+
  guides(size = guide_legend(title.position="left",title.hjust = 0.5,order = 1, label.position = "bottom"),
         color = guide_colorbar(barwidth = 2.4, title = NULL, label = FALSE)) +
  
  theme(text=element_text(size=SIZE, color = "grey10", face = "plain"),
      axis.text.x=element_text(size=SIZE*0.65, color = "grey20", angle=45,hjust=1,vjust=1),
      axis.ticks=element_line(colour="grey25",size=0.05)) +
  theme(panel.grid = element_line(color="grey90",size=0.1))+
  theme(legend.key.size=unit(0.5,"lines"),legend.position = c(-0.17,1.05), legend.direction = "horizontal",
        legend.box = "vertical", legend.box.just = "left",
        legend.background=element_blank(),
        legend.key=element_blank(),legend.text=element_text(size=SIZE*0.7, face="plain"),
        legend.spacing = unit(0,"mm"), legend.margin = margin(0,0,0,0,"mm"))+
  theme(strip.background=element_rect(fill= "grey60"),
        strip.text.x=element_text(color="white",face="bold",size=SIZE*0.8,angle=90,vjust=0.5,hjust=0))+
  theme(panel.spacing = unit(0.1,"lines"),panel.background = element_rect(fill=NA),
        panel.border=element_rect(fill=NA,color=NA))

  ggsave(file = paste0(Sys.Date(),"_Figure_2B.pdf"), dpi=300, scale = 1,
         width=8, height=15, units="cm")

#>> Rendering Figure 2C and calculating the Hedges g effect sizes  -----------------------------------------
  
# Defining the split violin plot visualisation
  
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                             draw_group = function(self, data, ..., draw_quantiles = NULL) {
                               data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                               grp <- data[1, "group"]
                               newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                               newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                               newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                               
                               if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                 stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                           1))
                                 quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                 aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                 aesthetics$alpha <- rep(1, nrow(quantiles))
                                 both <- cbind(quantiles, aesthetics)
                                 quantile_grob <- GeomPath$draw_panel(both, ...)
                                 ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                               }
                               else {
                                 ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                               }
                             })
  
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                                draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                                show.legend = NA, inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
          position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
          params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# Plotting

SIZE = 9

ggplot(subset(EV.2C, cat2 == "Breast"), aes(x= Gene, y=Expr, colour = X_study, fill = X_study)) +
  theme_minimal(base_size = 14) +
  
  geom_split_violin(trim = FALSE, scale = "width", adjust = 0.5, size = 0.1,  width = 1) +
  stat_summary(fun = mean, geom = "point", size = 2, shape = 23, stroke = 0.5, alpha = 1,
               show.legend = FALSE, fill = "white",
               position = position_dodge(width = 0)) +
  
  facet_grid(.~Type, scales = "free_x", space = "free") +
  
  scale_fill_jco(name = "Breast tissue: ", alpha = 0.3, labels = c("non-tumour (n=178) ","tumour (n=1099)")) +
  scale_color_jco(guide = FALSE, alpha = 0.9) +
  scale_y_continuous(position = "left", expand = c(0,0.5), minor_breaks = NULL) +
  scale_x_discrete(position = "top")+
  
  guides(fill = guide_legend(override.aes = list(color = adjustcolor(get_palette("jco",2), alpha.f = 0.8),
                                                 fill = adjustcolor(get_palette("jco",2), alpha.f = 0.5)),
                             title.position = "left", title.vjust = 0.5)) +
  
  theme(strip.background=element_rect(fill= "grey90", color = "white", size = 1),
        strip.text.x=element_text(color="black",face="bold",size=SIZE,angle=0,vjust=0.5,hjust=0.5)) +
  theme(panel.spacing = unit(0,"lines"),panel.background = element_rect(fill=NA, color = NA),
        panel.border=element_rect(fill=NA,color=NA)) +
  theme(axis.ticks=element_line(colour="grey35",size=0.1), axis.ticks.length = unit(3,"pt"),
        panel.grid.major.x = element_line(color="grey80",size=0.2, linetype = 2), 
        text = element_text(colour="grey35", size = 12),
        panel.grid.major.y = element_blank(),
        legend.title = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(face = "plain", vjust = 0, hjust = 0.5, size = 10, color = "black"),
        legend.position = "bottom", legend.direction = "horizontal", legend.key.size=unit(1,"lines")) +
  labs(x = NULL, y = expression("log"[2]*" normalised RNA count"), title = "") +
  coord_cartesian(ylim = c(-3,20), clip = "off") +
  
# Hedges g effect size calculation and plot labeling  
  
  geom_label(data = setNames(cbind(ddply(subset(EV.2C, cat2 == "Breast"),.(Type, Gene),
                                         function(z) {cohen.d(Expr~X_study, data = z,
                                                              hedges.correction = TRUE, paired=FALSE)}[["estimate"]]*-1), value = -3.4),
                             c("Type", "Gene","lab","Expr")), aes(label = sprintf("%0.2f", round(lab, digits = 2)),
                                                                  x=Gene,y=Expr),
             colour = "grey50", inherit.aes = F, size = 3.5, label.size = unit(0,"lines"), label.padding = unit(2,"pt"),
             parse = T) +
  geom_text(data = subset(EV.2C, cat2 == "Breast" & Gene %in% c("EPCAM") )[1,],
            x = 0.2, y=-3.2, label = "Hedges g:", color = "grey 50", size = 3.5)

# Saving the plot

ggsave(file = paste0(Sys.Date(),"_Figure_2C.pdf"), dpi=300, scale = 1,
       width=16, height=10, units="cm")
