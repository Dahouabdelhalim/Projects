
#### Script 4 - Coefficient Plot

# Load packages
library(dplyr)

# Import RDS file created in Script 2.
model_rescale<-readRDS("model_rescale.rds")

# Function to get CIs from models. Written by Sean Anderson, from Artelle et al 2016 Ecology of Conflict manuscript:  https://github.com/kartelle/Ecology-of-Conflict/blob/master/2-plot.R
get_cis <- function(x) {
  fe <- coef(x)
  ci <- confint(x, level = .95)
  ci2 <- confint(x, level = .5)
  d <- data.frame(effect = names(fe), fe = as.numeric(fe), l = ci[,1],
                  ll = ci2[,1], uu = ci2[,2], u = ci[,2], stringsAsFactors = FALSE)
  rownames(d) <- NULL
  d <- dplyr::filter(d, effect != "(Intercept)")
  d
}

model_rescale_cis<-get_cis(model_rescale)

# Extract label names from external CSV (set manually)
labs <- read.csv("model-rescale-coef-names.csv", stringsAsFactors = FALSE, strip.white = TRUE)

plot_model_rescale_coefs<-inner_join(labs, model_rescale_cis)%>%
  arrange(fe)%>%
  mutate(pos=1:length(label))

# Function to create individual panel 
make_panel <- function(dat, xlim=default_xlim ,ticks=default_ticks ,ylim=default_ylim,
                       add_importance = FALSE, label = "",
                       y_axis_label_offset=2.4) {
  
  plot(1, 1, xlim = xlim, ylim = ylim, xlab = "", ylab = "", type = "n",
       axes = FALSE)
  par(xpd = FALSE)
  abline(v = 0, col = "grey40", lty = "22")
  with(dat, segments(l, pos, u, pos, lwd = 0.5))
  with(dat, segments(ll, pos, uu, pos, lwd = 1.7))
  with(dat, points(fe, pos, pch = 21, bg = "white", col = "grey20"))
  axis(1, at = ticks, labels = ticks, col = "grey40")
  par(xpd = NA)
  with(dat, text(rep(xlim[1] * y_axis_label_offset, length(label)),
                 pos + 0.1, label, pos = 4, col = "grey40"))
  text(xlim[1] * y_axis_label_offset, ylim[2]+1, label, pos = 4, col = "grey10")
  
}

# Set x and ylim based on maxes of all data
default_xlim=c(-2.7, 4.5)
default_ylim = c(0.8, 5.4)
default_ticks = c(-2, -1, 0, 1, 2, 3, 4) 

tiff("coefficient plot.tif", width = 6, height = 1.6, units="in", res=300)
par(mfrow = c(1, 1))
par(las = 1, cex = 0.7, mar = c(2, 18, 2, 0), oma = c(0, 0,0,0),
    mgp = c(2, 0.25, 0), tck = -0.025, col.axis = "grey40")
make_panel(plot_model_rescale_coefs, y_axis_label_offset=3)
dev.off()
