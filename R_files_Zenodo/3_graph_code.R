
# --------------------------------------------------------------------------------
# Install and load packages

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

packages <- c("corrplot", "tidyverse", "geomorph", "ggplot2", "ggpubr", "remotes", "ggConvexHull")
ipak(packages)


# --------------------------------------------------------------------------------
# plot correlation matrices for raw and size normalized data

setwd("path-to-directory") # edit file path for working directory

load("gls_models_and_data.Rdata")


responses <- c("zygomatic_width", "cranial_vault_width", "upper_jaw_width", "total_skull_length", "snout_length", "cranial_vault_height", "endocranial_volume")

pdf("cor_response.pdf", height = 4.5, width = 4.5)
corrplot(cor(dat[, responses]), method = "ellipse", type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.7, diag = FALSE, 
         cl.offset = 0.2, cl.ratio = 0.25, cl.align.text = "l", cl.cex = 0.7)
dev.off()

responses_norm <- paste0(responses, "_norm") 

pdf("cor_response_norm.pdf", height = 4.5, width = 4.5)
corrplot(cor(dat[, responses_norm]), method = "ellipse", type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.7, diag = FALSE,
         cl.offset = 0.2, cl.ratio = 0.25, cl.align.text = "l", cl.cex = 0.7)
dev.off()



# --------------------------------------------------------------------------------
# load fitted models and extract diagnostics

setwd("path-to-directory") # edit file path for working directory

load("gls_models_and_data.Rdata")

# create data frame with model diagnostics
diagnostics <- data.frame(
  model = c(rep("Centroid Size", 73), rep("Population Differences", 511), rep("Sexual Dimorphism", 511)),
  resid = c(resid(model_c, type = "p"), resid(model1, type = "p"), resid(model2, type = "p")),
  fitted = c(fitted(model_c), fitted(model1), fitted(model2))
)

str(diagnostics)


# plot model diagnostics

plot1 <- ggplot(diagnostics, aes(x = fitted, y = resid)) +
  geom_point(shape = 16, alpha = 0.15, color = "darkred") +
  geom_smooth(se = FALSE) +
  labs(x = "Fitted Values", y = "Standardized residuals") +
  theme_bw() +
  facet_wrap(~ model, scales = "free_x") +
  theme(axis.title = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"))
ggsave(plot1, file = "resid_fitted.pdf", height = 2, width = 5.5)


plot2 <- ggplot(diagnostics, aes(sample = resid)) +
  geom_qq(shape = 16, alpha = 0.1, color = "darkred") +
  geom_qq_line(color = "blue") + 
  labs(x = "Theoretical Gaussian Quantiles", y = "Sample Quantiles") +
  theme_bw() +
  coord_fixed() +
  facet_wrap(~ model) +
  theme(axis.title = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"))
ggsave(plot2, file = "qq_plot.pdf", height = 2, width = 6)



# --------------------------------------------------------------------------------
# load allometric data

setwd("path-to-directory") # edit file path for working directory

load("shape_objects.Rdata")


# plot allometric relationships

allometry_plot <- ggplot(allometry, aes(x = csize_ln, y = shape_score, 
                                        color = population, fill = population)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_brewer(palette = "Set1", name = "Population") +
  scale_fill_brewer(palette = "Set1", name = "Population") +
  facet_grid(shape_metric ~ allometry_type) +
  labs(x = "Centroid Size", y = "Shape Score") +
  theme_bw() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.spacing = unit(0, "cm"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.3, "lines"))

ggsave(allometry_plot, file = "allometry.pdf", height = 3.2, width = 5.5)

## ----------------------------------------------------------------------------
 
## Box and Jitter Plots for Centroid Size and Linear and Endocranial Measurements
#using data already loaded (gls_models_and_data.Rdata)

### Stats bars in these graphs come from plots in the model outputs
p_value_o <- data.frame(
  x = c(1, 1, 3, 3),
  y = c(10.6, 10.65, 10.65, 10.6)
)

p_value_t <- data.frame(
  x = c(2, 2, 3, 3),
  y = c(10.48, 10.53, 10.53, 10.48)
)

## centroid size
C.size <- ggplot(dat, aes(x = population, y = centroid_size, color = sex)) + 
  geom_boxplot(
    aes(color = sex), 
    position = position_dodge(0.7)) +
  geom_point(position = position_jitterdodge(),alpha=.3) +
  labs(x = "Fox Population", 
       y = "Centroid Size") +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2', name = "Sex", labels = c("Female", "Male")) +
  scale_x_discrete(breaks=c("Domesticated","Unselected","Wild"),
                   labels=c("Domesticated", "Unselected", "Wild")) +
  scale_y_continuous (breaks = seq(8.2,10.8, by=.2), limits = c(8.2,10.8)) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10, vjust = 0),
        axis.title.y = element_text(size = 10, vjust = 0)) +
  geom_line(data = p_value_o, 
            aes(x = x, y = y), inherit.aes = FALSE) +
  geom_line(data = p_value_t, 
            aes(x = x, y = y), inherit.aes = FALSE) +
  annotate("text", x = 2, y = 10.66, 
           label = "***", size = 4) + 
  annotate("text", x = 2.5, y = 10.54, 
           label = "***", size = 4)

## Snout Length
p_value_one <- data.frame(
  x = c(1, 1, 3, 3),
  y = c(8.33, 8.36, 8.36, 8.33)
)

p_value_two <- data.frame(
  x = c(2, 2, 3, 3),
  y = c(8.26, 8.29, 8.29, 8.26)
)

SL <- ggplot(dat, aes(x = population, y = snout_length_norm, color = sex)) + 
  geom_boxplot(
    aes(color = sex), 
    position = position_dodge(0.7)) +
  geom_point(position = position_jitterdodge(),alpha=.3) +
  labs(x = "Fox Population", 
       y = "Normalized Snout Length") +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2', name = "Sex", labels = c("Female", "Male")) +
  scale_x_discrete(breaks=c("Domesticated","Unselected","Wild"),
                   labels=c("Domesticated", "Unselected", "Wild")) +
  scale_y_continuous (breaks = seq(7.4,9, by=.2), limits = c(7.4,8.6)) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10, vjust = 0),
        axis.title.y = element_text(size = 10, vjust = 0)) +
  geom_line(data = p_value_one, 
            aes(x = x, y = y), inherit.aes = FALSE) +
  geom_line(data = p_value_two, 
            aes(x = x, y = y), inherit.aes = FALSE) +
  annotate("text", x = 2, y = 8.37, 
           label = "***", size = 4) + 
  annotate("text", x = 2.5, y = 8.30, 
           label = "***", size = 4)

## Cranial Vault Width 
p_value_three <- data.frame(
  x = c(1, 1, 3, 3),
  y = c(5.87, 5.90, 5.90, 5.87)
)

p_value_four <- data.frame(
  x = c(2, 2, 3, 3),
  y = c(5.79, 5.82, 5.82, 5.79)
)

CVW <- ggplot(dat, aes(x = population, y = cranial_vault_width_norm, color = sex)) + 
  geom_boxplot(
    aes(color = sex), 
    position = position_dodge(0.7)) +
  geom_point(position = position_jitterdodge(),alpha=.3) +
  labs(x = "Fox Population", 
       y = "Normalized Cranial Vault Width") +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2', name = "Sex", labels = c("Female", "Male")) +
  scale_x_discrete(breaks=c("Domesticated","Unselected","Wild"),
                   labels=c("Domesticated", "Unselected", "Wild")) +
  scale_y_continuous (breaks = seq(4.4,6, by=.2), limits = c(4.4,6.0)) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10, vjust = 0),
        axis.title.y = element_text(size = 10, vjust = 0)) +
  geom_line(data = p_value_three, 
            aes(x = x, y = y), inherit.aes = FALSE) +
  geom_line(data = p_value_four, 
            aes(x = x, y = y), inherit.aes = FALSE) +
  annotate("text", x = 2, y = 5.91, 
           label = "***", size = 4) + 
  annotate("text", x = 2.5, y = 5.83, 
           label = "***", size = 4)

## Endocranial Volume

p_value_five <- data.frame(
  x = c(1, 1, 3, 3),
  y = c(.455, .46, .46, .455)
  )

p_value_six <- data.frame(
    x = c(2, 2, 3, 3),
    y = c(.445, .450, .450, .445)
    )

ECV <- ggplot(dat, aes(x = population, y = endocranial_volume_cr_norm, color = sex)) + 
  geom_boxplot(
    aes(color = sex), 
    position = position_dodge(0.7)) +
  geom_point(position = position_jitterdodge(),alpha=.3) +
  labs(x = "Fox Population", 
       y = "Normalized Endocranial Volume") +
  theme_classic() +
  scale_color_brewer(palette = 'Dark2', name = "Sex", labels = c("F", "M")) +
  scale_x_discrete(breaks=c("Domesticated","Unselected","Wild"),
                   labels=c("Domesticated", "Unselected", "Wild")) +
  scale_y_continuous (breaks = seq(.3,.5, by=.02), limits = c(.3,.5)) +
  theme(legend.position = c(.93, .17),
            legend.box = "horizontal",
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 8),
            axis.title.x = element_text(size = 10, vjust = 0),
            axis.title.y = element_text(size = 10, vjust = 0)) +
  guides(color = guide_legend(override.aes = list(size = .5))) +
  geom_line(data = p_value_five, 
            aes(x = x, y = y), inherit.aes = FALSE) +
  geom_line(data = p_value_six, 
            aes(x = x, y = y), inherit.aes = FALSE) +
  annotate("text", x = 2, y = .461, 
               label = "***", size = 4) + 
  annotate("text", x = 2.5, y = .451, 
               label = "***", size = 4)

BoxJitter_plot <- ggarrange(C.size, CVW, SL, ECV, 
                            labels = c("A", "B", "C", "D"),
                            ncol = 2, nrow = 2)

ggsave(BoxJitter_plot, file = "BoxJit.pdf", height = 6.5, width = 6.5, units = "in")

#### --------------------------------------------------------------------------------
## Fox Sex Dimorphism Graph
#using data already loaded (gls_models_and_data.Rdata)

##Plot models
DomestSD <- subset(tab_sub, contrast == "M Domesticated - F Domesticated", select = c(estimate))
UnselSD <- subset(tab_sub, contrast == "M Unselected - F Unselected", select = c(estimate))
WildSD <- subset(tab_sub, contrast == "M Wild - F Wild", select = c(estimate))

SD <- cbind(DomestSD,UnselSD,WildSD)

DomestLCL <- subset(tab_sub, contrast == "M Domesticated - F Domesticated", select = c(lower.CL))
UnselLCL <- subset(tab_sub, contrast == "M Unselected - F Unselected", select = c(lower.CL))
WildLCL <- subset(tab_sub, contrast == "M Wild - F Wild", select = c(lower.CL))

LCL <- cbind(DomestLCL,UnselLCL,WildLCL)

DomestUCL <- subset(tab_sub, contrast == "M Domesticated - F Domesticated", select = c(upper.CL))
UnselUCL <- subset(tab_sub, contrast == "M Unselected - F Unselected", select = c(upper.CL))
WildUCL <- subset(tab_sub, contrast == "M Wild - F Wild", select = c(upper.CL))

UCL <- cbind(DomestUCL,UnselUCL,WildUCL)

LCI.values <- LCL

trait <- c(rep("CVH" , 3) , rep("CVW" , 3) , rep("ECV" , 3) ,
           rep("SL" , 3) , rep("TSL", 3), rep ("UJW", 3),
           rep("ZW", 3))
population <- rep(c("Domesticated" , "Unselected" , "Wild") , 7)

groups <- c("Domesticated", "Unselected", "Wild")

value <- unlist(c(SD[1,1:3], SD[2,1:3], SD[3,1:3], SD[4,1:3], SD[5,1:3], SD[6,1:3], SD[7,1:3]))  
value <- value
UCI <- unlist(c(UCL[1,1:3], UCL[2,1:3], UCL[3,1:3], UCL[4,1:3], UCL[5,1:3], UCL[6,1:3], UCL[7,1:3]))
UCI <- UCI
LCI <- unlist(c(LCL[1,1:3], LCL[2,1:3], LCL[3,1:3], LCL[4,1:3], LCL[5,1:3], LCL[6,1:3], LCL[7,1:3]))
LCI <- LCI
frame <- data.frame(trait, population, value, UCI, LCI)

##adding signficance codes after looking at model output
sig.lab <- c("*", "*", " "," "," "," ", "*", " ", " ", "*", "*", "*","*","*", 
             "*", "*", "*","*", "*", "*","*")
height <- c(rep.int(5,6), rep.int(2,3), rep.int(8,3), rep.int(15,3), rep.int(5,3), rep.int(8.25,3))


Dimorphism_graph <- ggplot(frame, aes(fill=population, y=value, x=trait)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width =.2, position = position_dodge(.9)) +
  scale_y_continuous (breaks = seq(0,15, by=1), limits = c(-.5,15)) +
  geom_text(aes(x = trait, y = height, label = sig.lab), 
            position = position_dodge(.9), size = 5) +
  scale_fill_brewer(palette = "BuGn", name = "Population") +
  xlab("Skull Measurement") +
  ylab("Average Difference Between Sexes Â± 95% Confidence Intervals") + 
  theme_classic() +
  theme(legend.position = "top",
        axis.title.x = element_text(size = 10, vjust = -2),
        axis.title.y = element_text(size = 10, vjust = 2))

ggsave(Dimorphism_graph, file = "Sex Dimorphism Graph.pdf", height = 5, width = 5)

##---------------------------------------------------------------------------------
##shape visualizations
# Shape change visualization
#create PCA
PCA <- gm.prcomp(A = FoxGPA$coords)

##PCA code

FoxPC1 <- PCA$x[,1]
FoxPC2 <- PCA$x[,2]

PC.dat <- data.frame(FoxPC1, FoxPC2, FoxDF$population, FoxDF$sex)

PCA_plot <- ggplot(PC.dat, aes(x = FoxPC1, y = FoxPC2, fill = FoxDF$population, color = FoxDF$population)) +
  geom_vline(xintercept = 0, color = "#f0f0f0") +
  geom_hline(yintercept = 0, color = "#f0f0f0") +
  geom_convexhull(alpha = 0.2) + 
  geom_point(aes(shape=FoxDF$sex, color = FoxDF.population), size = 3) +
  scale_shape_manual(values = c(16,4), name = "Sex", labels = c("Male", "Female")) +
  scale_color_brewer(palette = "Dark2", name = "Population", labels = c("Domesticated", "Unselected", "Wild")) +
  scale_fill_grey(guide = FALSE) +
  coord_fixed() +
  labs(x = "PC1 (31.9% of total variance explained)", 
       y = "PC2 (8.3% of total variance explained)") +
  theme_classic() +
  theme(legend.position = "top",
        axis.title.x = element_text(size = 12, vjust = -2),
        axis.title.y = element_text(size = 12, vjust = 2))

ggsave(PCA_plot, file = "Fox_PCA.pdf", height = 5.5, width = 8.5, units = "in")

## Wireframes links

Lat.wireframe.links <- Fox3D$wireframe
Dor.wireframe.links <- Lat.wireframe.links[-c(10,14,19:32),] 

## Grab skull warps for mean shape and the extreme ends of PC1 axes

PCA1.min <-PCA$shapes$shapes.comp1$min
Dor.PCA1.min <- PCA1.min[-c(19:29),-c(2)]*-1
Dor.PCA1.min <- mshape(Dor.PCA1.min)
Lat.PCA1.min <- PCA1.min[,-c(3)]*-1
Lat.PCA1.min <- mshape(Lat.PCA1.min)

PCA1.max <- PCA$shapes$shapes.comp1$max
Dor.PCA1.max <- PCA1.max[-c(19:29),-c(2)]*-1
Dor.PCA1.max <- mshape(Dor.PCA1.max)
Lat.PCA1.max <- PCA1.max[,-c(3)]*-1
Lat.PCA1.max <- mshape(Lat.PCA1.max)

#plot wireframes

#Dorsal View
plotRefToTarget(Dor.PCA1.min, Dor.PCA1.max, method = "points", mag = 2, gridPars=gridPar(tar.pt.bg = "red", tar.link.col="red"), links = Dor.wireframe.links)

#Lateral View
plotRefToTarget(Lat.PCA1.min, Lat.PCA1.max, method = "points", mag = 2, gridPars=gridPar(tar.pt.bg = "red", tar.link.col="red"),  links = Lat.wireframe.links)

#plot
png(file = "PC1_dorview.png", res = 120)
plotRefToTarget(Dor.PCA1.min, Dor.PCA1.max, method = "points", mag = 2, gridPars=gridPar(tar.pt.bg = "red", tar.link.col="red"), links = Dor.wireframe.links)
dev.off()

png(file = "PC1_latview.png", res = 100)
plotRefToTarget(Lat.PCA1.min, Lat.PCA1.max, method = "points", mag = 2, gridPars=gridPar(tar.pt.bg = "red", tar.link.col="red"),  links = Lat.wireframe.links)
legend("bottom", legend = c("PC1 min", "PC1 max"), horiz = TRUE, pch = c(16,16), cex = 1.2, col = c("gray", "red"), bty = "n")
dev.off()

##------------------------------------------------------------------------------------
## ProcD comparison plot

#grab table for plot
labels <- c("D-U", "D-W", "U-W")
mean_shape_table <- cbind(labels,unlist(mean_shape_diff$summary.table$d),
                          unlist(mean_shape_diff$summary.table$`UCL (95%)`))
colnames(mean_shape_table) <- c("Contrasts", "Dist.diff", "UCL")
mean_shape_table <- as.data.frame(mean_shape_table)
mean_shape_table$Dist.diff <- as.numeric(mean_shape_table$Dist.diff)
mean_shape_table$UCL <- as.numeric(mean_shape_table$UCL)

##Procrustes Distance Plot
sig.lab <- c("*", "***", "***") # as observed from pairwise comparison summary table
height <- c(rep.int(0.055, 3))

PDplot <- ggplot(mean_shape_table, aes(fill= Contrasts, y=Dist.diff, x=Contrasts)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin = UCL, ymax = 1), width = .9, position = position_dodge(.9)) +
  geom_text(aes(x = Contrasts, y = height, label = sig.lab), 
            position = position_dodge(.9), size = 5) + 
  scale_y_continuous(limits = c(0,0.06)) +
  xlab("Contrasts") +
  ylab("Estimate of difference in Procrustes Distance") +
  scale_fill_brewer(palette = "BuGn", name = "Contrasts") +
  theme_classic() +
  theme(legend.position = "top",
        axis.title.x = element_text(size = 10, vjust = -2),
        axis.title.y = element_text(size = 10, vjust = 2))

ggsave(PDplot, file = "ProcDplot.pdf", height = 5.5, width = 5.5)

