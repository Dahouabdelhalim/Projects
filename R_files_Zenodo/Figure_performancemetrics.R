usage_msg <- "
Usage: Rscript Figure_performancemetrics.R indatafile stylefile outdir outprefix
"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    errormsg <- sprintf('Error: Invalid number of arguments %d', length(args))
    print(errormsg)
    print(usage_msg)
    quit(save="no", status=1)
}
inputfile <- args[1]
stylefile <- args[2]
outdir <- args[3]
outprefix <- args[4]

imagefile <- sprintf("%s/%s.png", outdir, outprefix)

library(RColorBrewer)
library(reshape2)
library(methods)
library(ggplot2)
library(grid)

source(stylefile)
plot_width <- 210
plot_height <- 210
plot_units <- "mm"
plot_resolution <- 200
line_width <- 0.6
legend_key_width <- 1.0
legend_key_height <- 1.0
legend_width <- 500
legend_width_units <- "mm"
vert_24h_line_width <- 0.25
vert_24h_line_type <- "dashed"
vert_24h_line_colour <- medgrey
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
style <- style + theme(panel.grid.major = element_line(
    colour=grid_major_colour,size=grid_major_size))
style <- style + theme(axis.ticks.x = element_blank())
style <- style + theme(axis.ticks.y = element_blank())
style <- style + theme(plot.margin=unit(c(0.0,0.0,0,0), "cm"))

Construct_Figure <- function()
{
    data <- read.table(inputfile, header=TRUE, sep='\\t')
    D <- data[data$valuetype == "Mean",]
    D$Experiment <- factor(D$Experiment, levels=experimentorder)
    D$metric <- factor(D$metric, levels = 
        c("Length", "Q-score", "BQ", "GC", "GC (1D)", "Speed (1D)", "Count"))
    pngpath <- imagefile

    p <- ggplot(D, aes(x=time, y=value, colour=Experiment))
    p <- p + geom_line(data=D, size=line_width)
    p <- p + scale_colour_manual(values=exptpalette)
    p <- p + labs(x="Time (h)", y="")
    p <- p + style
    p <- p + theme(legend.position="bottom") +
        guides(colour = guide_legend(override.aes = list(size=3)))
    p <- p + theme(plot.margin=unit(c(0,0,0,0), "cm"))
    p <- p + theme(plot.title = element_text(
        hjust = subplotlabel_offset, size=plot_title_font)) # -0.071
    p <- p + theme(axis.text.x = element_text(angle=90, hjust = 1))
    p <- p + facet_grid(metric ~ Experiment, scales = "free")
    if (max(D$time) > 24) {
        p <- p + geom_vline(aes(xintercept=24),
            linetype=vert_24h_line_type, size=vert_24h_line_width,
            colour=vert_24h_line_colour)
    }
    ggsave(pngpath, width=plot_width, height=plot_height, units=plot_units)
    return(p)
}

Construct_Figure()
