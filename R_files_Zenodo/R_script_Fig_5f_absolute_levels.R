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


qPCR_data$log_expression = log(qPCR_data$Normalized_expression)

qPCR_data$condition=factor(rep("other",dim(qPCR_data)[1]),levels=c("negative","test","positive","other"))

qPCR_data$condition[(qPCR_data$Cell=="ADSC" & qPCR_data$Substrate=="imprinted")]="test"

qPCR_data$condition[(qPCR_data$Cell=="Rencell" & qPCR_data$Substrate=="TCP") ]="positive"

qPCR_data$condition[(qPCR_data$Cell=="ADSC" & qPCR_data$Substrate=="TCP") ]="negative"

genes_Rencell=unique(qPCR_data$Gene[qPCR_data$Cell=="Rencell"])
genes_ADSC = unique(qPCR_data$Gene[qPCR_data$Cell=="ADSC"])

common_genes = genes_Rencell[genes_Rencell %in%genes_ADSC]
common_genes = common_genes [common_genes != "GAPDH"] # Remove housekeeping

qPCR_common = qPCR_data[qPCR_data$Gene %in% common_genes,]

aov_contrast_neg_pos=aov(log_expression~paste(condition,Gene),qPCR_common[ qPCR_common$condition %in% c("negative","positive"),])

TukeyHSD(aov_contrast_neg_pos)

# Estimate a common standard deviation from this
# From TukeyHSD (i.e. getAnywhere(TukeyHSD.aov))
 mm <- model.tables(aov_contrast_neg_pos, "means")

# Mean tables
tabs <- mm$tables

if (names(tabs)[1L] == "Grand mean") tabs <- tabs[-1L]

# number of replicates for the various conditions
nn <- mm$n[names(tabs)]

# mean squared error of the residuals. This is to estimate the common variance
MSE <- sum(aov_contrast_neg_pos$residuals^2)/aov_contrast_neg_pos$df.residual
# Table with the means
tab <- tabs[["paste(condition, Gene)"]]
# Names of all the conditions being compared
nms=names(tab)
means <- as.vector(tab)
# Number of replicates
n <- nn[["paste(condition, Gene)"]]
# We are only interested in the comparison within each gene, not from gene to gene
center = means[(length(means)/2+1):length(means)]-means[1:(length(means)/2)]
effective_n = 1/(1/n[(length(n)/2+1):length(n)]+1/n[1:(length(n)/2)])

# We are doing only m=length(means)/2 comparisons, while usually, the Tukey test would do
# n*(n-1)/2 comparisons
# Quadratic resolution formula (-b + sqrt(b^2-4*a*c))/2/a
# so we'd have n^2/2-n/2 - m = 0 and so n=(1/2+sqrt((1/2)^2+2*m))/1=1/2+sqrt((1/2)^2+2*m)
n_means_equivalent =1/2+sqrt((1/2)^2+length(means))
est <- center/sqrt((MSE/2) /effective_n)
pvals <- ptukey(abs(est), n_means_equivalent, aov_contrast_neg_pos$df.residual,
lower.tail = FALSE)

genes_for_efficiency_assessment = common_genes[pvals<=0.05]


per_condition=aggregate (Normalized_expression ~ Cell + Substrate + Gene, qPCR_common,FUN=mean)

per_condition_sd=aggregate (Normalized_expression ~ Cell + Substrate + Gene, qPCR_common,FUN=sd)






qPCR_direct_comparison = qPCR_common[qPCR_common$Substrate=="TCP" | (qPCR_common$Cell=="ADSC" & qPCR_common$Substrate=="imprinted"),]

qPCR_direct_comparison$condition=factor(as.character(qPCR_direct_comparison$condition),levels=levels(qPCR_direct_comparison$condition)[1:(length(levels(qPCR_direct_comparison$condition))-1)])






positive_control = per_condition[per_condition$Substrate=="TCP" & per_condition$Cell=="Rencell",]
positive_control_sd = per_condition_sd[per_condition_sd$Substrate=="TCP" & per_condition_sd$Cell=="Rencell",]


negative_control = per_condition[per_condition$Substrate=="TCP" & per_condition$Cell=="ADSC",]
negative_control_sd = per_condition_sd[per_condition_sd$Substrate=="TCP" & per_condition_sd$Cell=="ADSC",]


test_condition = per_condition[per_condition$Substrate=="imprinted" & per_condition$Cell=="ADSC",]
test_condition_sd = per_condition_sd[per_condition_sd$Substrate=="imprinted" & per_condition_sd$Cell=="ADSC",]

for_barplot = matrix(ncol=3,nrow=dim(test_condition)[1])

colnames(for_barplot)=c("neg control","test","pos control")

rownames(for_barplot)=test_condition$Gene

for_barplot_sd = for_barplot


for_barplot[,"neg control"]=negative_control$Normalized_expression
for_barplot_sd[,"neg control"]=negative_control_sd$Normalized_expression

for_barplot[,"pos control"]=positive_control$Normalized_expression
for_barplot_sd[,"pos control"]=positive_control_sd$Normalized_expression

for_barplot[,"test"]=test_condition$Normalized_expression
for_barplot_sd[,"test"]=test_condition_sd$Normalized_expression

# We should do efficiency evaluation only on the genes that actually show a difference between positive and negative control
qPCR_statistics=qPCR_direct_comparison[qPCR_direct_comparison$Gene %in% genes_for_efficiency_assessment,]
qPCR_statistics$condition=factor(as.character(qPCR_statistics$condition),levels=c("test","negative","positive"))

a1<-aov(log_expression ~ condition +Gene, qPCR_statistics)

summary(glht(aov(log_expression ~ condition +Gene, qPCR_statistics),linfct=mcp(condition="Dunnett")))

for_barplot=for_barplot[c("Nestin","Tuj1","MAP2"),]
for_barplot_sd=for_barplot_sd[c("Nestin","Tuj1","MAP2"),]

barplot_with_errorbars(height=t(for_barplot),beside=TRUE,sd_height=t(for_barplot_sd),legend=TRUE,args.legend=list(x="topleft",cex=2))





