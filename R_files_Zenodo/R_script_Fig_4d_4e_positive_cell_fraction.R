library(readxl) # To run this R-script, the R-library readxl needs to be installed and functional
library(multcomp) # library for multiple comparisons

# Change this to wherever you decide to store positive_cell_fraction.xlsx
root_folder="/path/to/excel_file"

setwd(root_folder)

cell_positive_data = as.data.frame(read_excel("positive_cell_fraction.xlsx"))

# Helper function for plotting (Thomas Braschler)

significance_labels<-function (x, levels = c(0.1, 0.05, 0.01, 0.001), codes = c(".",
"*", "**", "***"))
{
    s = x
    s[] = ""
    for (level_index in 1:length(levels)) {
        s[x <= levels[level_index]] = codes[level_index]
    }
    return(s)
}


barplot_with_errorbars<-function (height, sd_height = NULL, beside = FALSE, horiz = FALSE,
group_order = NULL, sig_codes = NULL, sig.cex = 1, ...)
{
    args = list(...)
    if (!is.null(group_order)) {
        if (is.matrix(height)) {
            height = height[, group_order]
            sd_height = sd_height[, group_order]
            if (!is.null(sig_codes)) {
                sig_codes = sig_codes[, group_order]
            }
        }
        else {
            height = height[group_order]
            sd_height = sd_height[group_order]
            if (!is.null(sig_codes)) {
                sig_codes = sig_codes[group_order]
            }
        }
        if (!is.null(args[["names.arg"]])) {
            args[["names.arg"]] = args[["names.arg"]][group_order]
        }
    }
    y = c(height)
    if (is.null(sd_height)) {
        sd_y = 0
    }
    else {
        sd_y = c(sd_height)
    }
    if (is.null(args$ylim)) {
        upper_y = max(c(y, max(y + sd_y * 2, na.rm = TRUE)),
        na.rm = TRUE)
        lower_y = 0
        limit_y = c(lower_y, upper_y)
    }
    else {
        limit_y = args$ylim
    }
    args$ylim = limit_y
    args[["height"]] = height
    args[["beside"]] = beside
    args[["horiz"]] = horiz
    theLength = NULL
    if (!is.null(args[["length"]])) {
        theLength = args[["length"]]
        args[["length"]] = NULL
    }
    w = c(do.call(barplot, args))
    if (!beside & is.matrix(height)) {
        xpoints = matrix(data = w, nrow = nrow(height), ncol = length(w),
        byrow = TRUE)
    }
    else {
        xpoints = w
    }
    ypoints_start = height
    if (!beside) {
        if (is.matrix(height)) {
            ypoints_start = apply(height, 2, cumsum)
        }
    }
    if (horiz) {
        errorbar_args_positive = list(x = ypoints_start[ypoints_start >=
        0], y = xpoints[ypoints_start >= 0], sd_y = sd_y[ypoints_start >=
        0], code = 2)
        errorbar_args_negative = list(x = ypoints_start[ypoints_start <
        0], y = xpoints[ypoints_start < 0], sd_y = sd_y[ypoints_start <
        0], code = 1)
        errorbar_args_positive[["horiz"]] = TRUE
        errorbar_args_negative[["horiz"]] = TRUE
    }
    else {
        errorbar_args_positive = list(x = xpoints[ypoints_start >=
        0], y = ypoints_start[ypoints_start >= 0], sd_y = sd_y[ypoints_start >=
        0], code = 2)
        errorbar_args_negative = list(x = xpoints[ypoints_start <
        0], y = ypoints_start[ypoints_start < 0], sd_y = sd_y[ypoints_start <
        0], code = 1)
        if (!is.null(args[["offset"]])) {
            errorbar_args_positive[["y"]] = errorbar_args_positive[["y"]] +
            args[["offset"]]
            errorbar_args_negative[["y"]] = errorbar_args_negative[["y"]] +
            args[["offset"]]
        }
    }
    if (!is.null(theLength)) {
        errorbar_args_positive[["length"]] = theLength
        errorbar_args_negative[["length"]] = theLength
    }
    do.call(errorbars, errorbar_args_positive)
    do.call(errorbars, errorbar_args_negative)
    if (horiz) {
        if (!is.null(sig_codes)) {
            text(y + sd_y + limit_y[2] * 0.03, c(xpoints), c(sig_codes))
        }
    }
    else {
        if (!is.null(sig_codes)) {
            if (beside) {
                text(c(xpoints), y + sd_y + limit_y[2] * 0.03,
                c(sig_codes), cex = sig.cex)
            }
            else {
                text(c(xpoints), ypoints_start + sd_y + limit_y[2] *
                0.03, c(sig_codes), cex = sig.cex)
            }
        }
    }
    invisible(w)
}

errorbars<-function (x, y, sd_y = 0, angle = 90, code = 3, horiz = FALSE,
...)
{
    ystart = y
    ystop = y
    if (code == 2 || code == 3) {
        ystop = y + sd_y
    }
    if (code == 1 || code == 3) {
        ystart = y - sd_y
    }
    if (horiz) {
        arrows(x0 = ystart, y0 = c(x), x1 = ystop, y1 = c(x),
        angle = angle, code = code, ...)
    }
    else {
        arrows(x0 = x, y0 = ystart, x1 = x, y1 = ystop, angle = angle,
        code = code, ...)
    }
}



cell_positive_data$Substrate = factor(cell_positive_data$Substrate,levels=c("TCP","PDMS","imprinted"))

nestin_data = cell_positive_data[cell_positive_data$Marker=="Nestin",]

tuj1_data = cell_positive_data[cell_positive_data$Marker=="Tuj1",]




aov_for_nestin_variance=aov(Fraction_positive~Substrate,nestin_data)
aov_for_tuj1_variance=aov(Fraction_positive~Substrate,tuj1_data)




summary(glht(aov_for_nestin_variance,linfct=mcp(Substrate="Dunnett")))

summary(glht(aov_for_tuj1_variance,linfct=mcp(Substrate="Dunnett")))



nestin_agg = aggregate(Fraction_positive~Substrate,nestin_data,FUN=mean)

tuj1_agg = aggregate(Fraction_positive~Substrate,tuj1_data,FUN=mean)

nestin_agg_sd = aggregate(Fraction_positive~Substrate,nestin_data,FUN=sd)

tuj1_agg_sd = aggregate(Fraction_positive~Substrate,tuj1_data,FUN=sd)


barplot_with_errorbars(nestin_agg$Fraction_positive,col=grey(c(1,0.625,0)),sd_height=nestin_agg_sd$Fraction_positive,sig_codes=c("",significance_labels(summary(glht(aov_for_nestin_variance,linfct=mcp(Substrate="Dunnett")))$test$pvalues)),sig.cex=3,main="Fraction Nestin positive",ylim=c(0,100))

dev.new()

barplot_with_errorbars(tuj1_agg$Fraction_positive,col=grey(c(1,0.625,0)),sd_height=tuj1_agg_sd$Fraction_positive,sig_codes=c("",significance_labels(summary(glht(aov_for_tuj1_variance,linfct=mcp(Substrate="Dunnett")))$test$pvalues)),sig.cex=3,main="Fraction tuj1 positive",ylim=c(0,100))

