usage_msg <- "
Usage: Rscript Figure_readlengths.R indatafile stylefile readlenmax outdir outprefix
"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
    errormsg <- sprintf('Error: Invalid number of arguments %d', length(args))
    print(errormsg)
    print(usage_msg)
    quit(save="no", status=1)
}
inputfile <- args[1]
stylefile <- args[2]
readlenmax <- args[3]
outdir <- args[4]
outprefix <- args[5]

imagefile_counts <- sprintf("%s/%s.png", outdir, outprefix)

library(RColorBrewer)
library(reshape2)
library(methods)
library(ggplot2)
library(grid)

source(stylefile)
plot_width <- 210
plot_height <- 180
plot_units <- "mm"
plot_resolution <- 200
line_width <- 0.6
legend_key_width <- 1.0
legend_key_height <- 1.0
legend_width <- 500
legend_width_units <- "mm"
vert_24h_line_width <- 0.25
text_size <- 13
font_family <- "Helvetica"
grid_major_colour <- "lightgrey"
grid_major_size <- 0.30
plot_margin_bottom <- 0
plot_title_font <- 10
std_point_size <- 4
subplotlabel_offset <- -0.06
linecolour <- "black"

style <- theme_bw(base_size=text_size, base_family=font_family)
style <- style + theme(panel.grid.major = element_line(colour=grid_major_colour,size=grid_major_size))
style <- style + theme(axis.ticks.x = element_blank())
style <- style + theme(axis.ticks.y = element_blank())
style <- style + theme(plot.margin=unit(c(0.0,0.0,0,0), "cm"))

Construct_Figure <- function()
{
    data <- read.table(inputfile, header=TRUE, sep='\\t')
    D <- data
    D$Experiment <- factor(D$Experiment, levels=experimentorder)
    D$readtype <- factor(D$readtype, levels=c("Template", "2D"))

  # Count histogram thing - only reads <= readlenmax bases
    pngpath <- imagefile_counts
    p <- ggplot(data=D[data$value <= readlenmax,], aes(x=value/1000, fill=Experiment, colour=Experiment))
    p <- p + geom_histogram(stat="bin", binwidth=1000/1000)
    p <- p + scale_colour_manual(values=exptpalette)
    p <- p + scale_fill_manual(values=exptpalette)
    p <- p + labs(x="Read length (K)", y="Frequency")
    p <- p + style
    p <- p + theme(legend.position="bottom")
    p <- p + guides(colour = guide_legend(override.aes = list(size=1.5)))
    p <- p + theme(plot.margin=unit(c(0,0,0,0), "cm"))
    p <- p + theme(plot.title = element_text(hjust = subplotlabel_offset, size=plot_title_font)) # -0.071
    p <- p + theme(axis.text.x = element_text(angle=0, hjust = 1))
    p <- p + scale_y_continuous(breaks = c(2000,4000,6000,8000,10000,12000))
    p <- p + facet_grid(readtype ~ Experiment)
    ggsave(pngpath, width=plot_width, height=plot_height, units=plot_units)
}

Construct_Figure()
