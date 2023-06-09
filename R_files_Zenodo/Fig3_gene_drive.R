library(ggplot2)
library(metR)
library(RColorBrewer)

results <- read.csv("fig3_DW.csv")
results$sigma[results$sigma < 0] <- NA
mybreaks <- seq(0,1,length.out=10)
results$mu <- as.factor(results$mu)
results$c <- as.factor(results$c)
levels(results$mu) <- c("mu~'='~'1e'^-7","mu~'='~'1e'^-4")
levels(results$c) <- c("c~'='~0.7","c~'='~0.95")

results <- na.omit(results)

# Define the number of colors you want
nb.cols <- 10
mycolors <- colorRampPalette(c("cornsilk2", brewer.pal(8, "GnBu")[3:8],"Darkblue"))(nb.cols)

# Come up w/ palette that doesn't require so much ink to print
DW <- ggplot(results, aes(h, s)) + 
      facet_grid(vars(c), vars(mu), labeller=label_parsed) +
      geom_contour_filled(aes(z=sigma)) + 
      geom_smooth(data=results[results$sigma == min(results$sigma[results$sigma>0], na.rm=T), ], aes(h, s), col="dimgrey", lty="dashed", se=F) +
      scale_fill_manual(name=expression(sigma),
                        values = mycolors,
                        guide = guide_coloursteps(show.limits=T, frame.colour = "black",
                                                  ticks.colour = "black",
                                                  barwidth=20)) +
      theme(legend.position = "bottom") +
      xlim(c(0,1)) + 
      ylim(c(0,1)) +
      ylab("Fitness cost of drive (s)") +
      xlab("Heterozygote effect (h)")

ggsave("Fig3_DW.pdf", DW)
                           
