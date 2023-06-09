library(readxl) # To run this R-script, the R-library readxl needs to be installed and functional
library(multcomp) # library for multiple comparisons

# Change this to wherever you decide to store qPCR_Comparison.xlsx
root_folder="/path/to/excel_file"

setwd(root_folder)

qPCR_data = as.data.frame(read_excel("qPCR_Comparison.xlsx"))


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



# Statistical analysis on the log transformed expression, this normalizes the variance (and also corresponds to the original data which is Ct values)

qPCR_data$log_expression = log(qPCR_data$Normalized_to_TCP_control)

# only ADSC to be analyzed here, and exclude house-keeping gene (this is just for normalization)
qPCR_data_ADSC = qPCR_data[qPCR_data$Cell == "ADSC" & qPCR_data$Gene != "GAPDH",]

aov_for_overall_variance=aov(log_expression~paste(Substrate,Gene),qPCR_data_ADSC[ qPCR_data_ADSC$Substrate != "TCP",])



# Estimate a common standard deviation from this
# From TukeyHSD (i.e. getAnywhere(TukeyHSD.aov))
 mm <- model.tables(aov_for_overall_variance, "means")

# Mean tables
tabs <- mm$tables

if (names(tabs)[1L] == "Grand mean") tabs <- tabs[-1L]

# number of replicates for the various conditions


# mean squared error of the residuals. This is to estimate the common variance
MSE <- sum(aov_for_overall_variance$residuals^2)/aov_for_overall_variance$df.residual
# Table with the means
tab <- tabs[["paste(Substrate, Gene)"]]

# Names of all the conditions being compared
nms=names(tab)
means <- as.vector(tab)
names(means)=nms
nn <- rep(mm$n[names(tabs)],length(means))
# Number of replicates
n <- nn
# We are only interested in the comparison within each gene, not from gene to gene
center = means
effective_n = n

# We are doing only m=length(means)/2 comparisons, while usually, the Tukey test would do
# n*(n-1)/2 comparisons
# Quadratic resolution formula (-b + sqrt(b^2-4*a*c))/2/a
# so we'd have n^2/2-n/2 - m = 0 and so n=(1/2+sqrt((1/2)^2+2*m))/1=1/2+sqrt((1/2)^2+2*m)
n_means_equivalent =1/2+sqrt((1/2)^2+length(means))
est <- center/sqrt((MSE/2) /effective_n)
pvals <- ptukey(abs(est), n_means_equivalent, aov_for_overall_variance$df.residual,
lower.tail = FALSE)

# Now we have pvalues and can plot

qPCR_data_ADSC$Substrate = factor(qPCR_data_ADSC$Substrate,levels=c("TCP","PDMS","imprinted"))

qPCR_data_ADSC_agg=aggregate(Normalized_to_TCP_control ~ Substrate+Gene,qPCR_data_ADSC,FUN=mean)
qPCR_data_ADSC_agg_sd=aggregate(Normalized_to_TCP_control ~ Substrate+Gene,qPCR_data_ADSC,FUN=sd)

for_barplot=matrix(nrow=length(unique(qPCR_data_ADSC_agg$Gene)),ncol=length(unique(qPCR_data_ADSC_agg$Substrate)))

rownames(for_barplot)=c("Nestin","Tuj1","MAP2","Col1")
colnames(for_barplot)=unique(qPCR_data_ADSC_agg$Substrate)

for_barplot_sd=for_barplot
for_barplot_pval = for_barplot

for(theSubstrate in colnames(for_barplot))
{
    for(theGene in rownames(for_barplot))
    {
        for_barplot[theGene,theSubstrate]=qPCR_data_ADSC_agg$Normalized_to_TCP_control[qPCR_data_ADSC_agg$Gene==theGene & qPCR_data_ADSC_agg$Substrate==theSubstrate]
        for_barplot_sd[theGene,theSubstrate]=qPCR_data_ADSC_agg_sd$Normalized_to_TCP_control[qPCR_data_ADSC_agg_sd$Gene==theGene & qPCR_data_ADSC_agg_sd$Substrate==theSubstrate]
        if(theSubstrate=="TCP")
        {
            for_barplot_pval[theGene,theSubstrate]=1 # By definition, this is normalized data
        } else
        {
            for_barplot_pval[theGene,theSubstrate]=pvals[paste(theSubstrate,theGene,sep=" ")]
        }
        
        
        
    }
}



barplot_with_errorbars(t(for_barplot),sd_height=t(for_barplot_sd),beside=TRUE,col=grey(c(1,0.5,0.2)),sig_codes=t(significance_labels(for_barplot_pval)),sig.cex=3,xlab="Gene",ylab="Gene expression normalized to GAPDH, relative to TCP")










# Statistical assessment on the significant genes whether PDMS also has an effect


genes_for_efficiency_assessment_PDMS = c("Nestin","Tuj1","MAP2")



qPCR_efficiency_PDMS = qPCR_data_ADSC[qPCR_data_ADSC$Gene %in% genes_for_efficiency_assessment_PDMS,]

qPCR_efficiency_PDMS$Substrate = factor(qPCR_efficiency_PDMS$Substrate,levels=c("TCP","PDMS","imprinted"))

a1<-aov(log_expression ~ Substrate +Gene, qPCR_efficiency_PDMS)

summary(glht(a1,linfct=mcp(Substrate="Dunnett")))







