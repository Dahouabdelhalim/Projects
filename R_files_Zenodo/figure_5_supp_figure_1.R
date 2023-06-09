rm(list = ls())

pcks <- c("magrittr", "plyr", "ggplot2", "ggthemes", "knitr", "gridExtra", "broom", "lme4")
sapply(pcks, require, char = TRUE)


setwd(".../Data_and_Code")

abies <- read.csv("abies.csv")
rhodo <- read.csv("rhododendron.csv")

rhodo <- mutate(rhodo, DTreeLine = (elevation - AbiesTLelev))
abies <- mutate(abies, DTreeLine = (elevation - AbiesTLelev))

head(rhodo)
table(rhodo$tq)


# model Rhodo abundance
# =============================

head(rhodo)
nrow(rhodo)
rm(rhodo.data)
rhodo.data <- rhodo[,c("AbiesTL","sizeclass", "freq_live", "DTreeLine", "AbiesTLabsDist", "Slope","Canopy","Stumps","aspect1","transect")] %>% na.omit
head(rhodo.data)
nrow(rhodo.data)

rm(fit.small.above, fit.small.below, fit.medium.above, fit.medium.below, fit.large.above, fit.large.below, summary.small, summary.medium, summary.large)

fit.small.above <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(rhodo.data, sizeclass == "0-1m" & AbiesTL == "above treeline"), REML = FALSE)

fit.small.below <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(rhodo.data, sizeclass == "0-1m" & AbiesTL == "below treeline"), REML = FALSE)

fit.medium.above <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(rhodo.data, sizeclass == "1-2m" & AbiesTL == "above treeline"), REML = FALSE)

fit.medium.below <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(rhodo.data, sizeclass == "1-2m" & AbiesTL == "below treeline"), REML = FALSE)

fit.large.above <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(rhodo.data, sizeclass == "2m+" & AbiesTL == "above treeline"), REML = FALSE)

fit.large.below <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(rhodo.data, sizeclass == "2m+" & AbiesTL == "below treeline"), REML = FALSE)



prepdf <- function(df){
    df %>% subset(term != "(Intercept)") %>% subset(term != "sd_(Intercept).transect") %>% 
        mutate(low = estimate - 1.96*std.error, high = estimate + 1.96*std.error)
}

summary.small <- rbind(data.frame(where = "above treeline", tidy(fit.small.above)),
                       data.frame(where = "below treeline", tidy(fit.small.below))) %>% prepdf
summary.medium <- rbind(data.frame(where = "above treeline", tidy(fit.medium.above)),
                        data.frame(where = "below treeline", tidy(fit.medium.below))) %>% prepdf 
summary.large <- rbind(data.frame(where = "above treeline", tidy(fit.large.above)),
                       data.frame(where = "below treeline", tidy(fit.large.below))) %>% prepdf 

rhodo.summary <- rbind(data.frame(size = "0-1m", summary.small), 
                       data.frame(size = "1-2m", summary.medium), 
                       data.frame(size = "2m+", summary.large))
rhodo.summary

# function to change name
namechange <- function(mydf) {
    mydf$Term <- NA
    mydf$Term[mydf$term == "scale(AbiesTLabsDist)"] <- "Distance\\ntreeline"
    mydf$Term[mydf$term == "scale(Slope)"] <- "Slope"
    mydf$Term[mydf$term == "scale(Canopy)"] <- "Canopy"
    mydf$Term[mydf$term == "scale(Stumps)"] <- "Stumps"
    mydf$Term[mydf$term == "scale(aspect1)"] <- "Aspect"
    return(mydf)
}

rhodo.summary <- namechange(rhodo.summary); rhodo.summary

# merge all output in a df
bubble.df <- data.frame(response = "rhodo.ab", rhodo.summary); bubble.df


ggplot(rhodo.summary, 
       aes(Term, estimate, col = where)) + 
    geom_hline(yintercept = 0, col="lightgrey") + # zero line
    geom_errorbar(aes(ymin = low, ymax = high, lty = p.value > 0.05), width = 0.3) + # draw error bars
    ggtitle("Rhododendron density") + 
    geom_point() + 
    geom_line(aes(group = where)) + 
    facet_wrap(~size) + 
    theme_few() +  # come from package ggthemes
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_hue(l=40)

# function to put one decimal to y axis
scaleFUN <- function(x) sprintf("%.1f", x)


rh.abun.ggplot <- ggplot(rhodo.summary, 
                         aes(Term, estimate, col = where)) + 
    geom_hline(yintercept = 0, col="lightgrey") + # zero line
    geom_errorbar(aes(ymin = low, ymax = high, lty = p.value > 0.05), width = 0.3) + # draw error bars
    ggtitle("(a) Rhododendron density") +
    geom_point() + 
    geom_line(aes(group = where)) + 
    scale_y_continuous(labels=scaleFUN) +
    facet_wrap(~size) + 
    theme_few() +  # come from package ggthemes
    scale_colour_hue(l=40) +
    theme(axis.text=element_text(size=10), 
          axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
          axis.title = element_text(size = 10.5),
          strip.text.x = element_text(hjust = 0, size = 10, face = "bold"),
          plot.margin=unit(c(.35,0,.35,0), "cm"),
          plot.title = element_text(size = 11.5, face = "bold", hjust = 0),
          legend.position = c(0.825, 0.225),
          legend.title=element_text(size=9),
          legend.text=element_text(size=9),
          legend.key.height=unit(0.8,"line"), # space between two entries within a legend type
          legend.spacing.y = unit(-0.1, "cm"), # space between two types of legend entries
          axis.title.x=element_blank())

plot(rh.abun.ggplot)



# model Abies abundance
# =============================

head(abies)
nrow(abies)
rm(abies.data)
abies.data <- abies[,c("AbiesTL","sizeclass", "freq_live", "DTreeLine", "AbiesTLabsDist", "Slope","Canopy","Stumps","aspect1","transect")] %>% na.omit
head(abies.data)
nrow(abies.data)

rm(fit.small.above, fit.small.below, fit.medium.above, fit.medium.below, fit.large.above, fit.large.below, summary.small, summary.medium, summary.large)

fit.small.above <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(abies.data, sizeclass == "0-1m" & AbiesTL == "above treeline"), REML = FALSE)

fit.small.below <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(abies.data, sizeclass == "0-1m" & AbiesTL == "below treeline"), REML = FALSE)

# aspect removed to attain convergence
fit.medium.above <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps)) + (1| transect), family = poisson(), data = subset(abies.data, sizeclass == "1-2m" & AbiesTL == "above treeline"), REML = FALSE)

fit.medium.below <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(abies.data, sizeclass == "1-2m" & AbiesTL == "below treeline"), REML = FALSE)

fit.large.above <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(abies.data, sizeclass == "2m+" & AbiesTL == "above treeline"), REML = FALSE)

fit.large.below <- glmer(freq_live ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), family = poisson(), data = subset(abies.data, sizeclass == "2m+" & AbiesTL == "below treeline"), REML = FALSE)


summary.small <- rbind(data.frame(where = "above treeline", tidy(fit.small.above)),
                       data.frame(where = "below treeline", tidy(fit.small.below))) %>% prepdf
summary.medium <- rbind(data.frame(where = "above treeline", tidy(fit.medium.above)),
                        data.frame(where = "below treeline", tidy(fit.medium.below))) %>% prepdf 
summary.large <- rbind(data.frame(where = "above treeline", tidy(fit.large.above)),
                       data.frame(where = "below treeline", tidy(fit.large.below))) %>% prepdf 

abies.summary <- rbind(data.frame(size = "0-1m", summary.small), 
                       data.frame(size = "1-2m", summary.medium), 
                       data.frame(size = "2m+", summary.large))
abies.summary
abies.summary <- namechange(abies.summary); abies.summary


bubble.df <- rbind(bubble.df, data.frame(response = "abies.ab", abies.summary)); bubble.df


ggplot(abies.summary, 
       aes(Term, estimate, col = where)) + 
    geom_hline(yintercept = 0, col="lightgrey") + # zero line
    geom_errorbar(aes(ymin = low, ymax = high, lty = p.value > 0.05), width = 0.3) + # draw error bars
    ggtitle("Abies density") + 
    geom_point() + 
    geom_line(aes(group = where)) + 
    facet_wrap(~size) + 
    theme_few() +  # come from package ggthemes
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_hue(l=40)



ab.abun.ggplot <- ggplot(abies.summary, 
                         aes(Term, estimate, col = where)) + coord_cartesian(ylim = c(-5, 4)) +
    geom_hline(yintercept = 0, col="lightgrey") + # zero line
    geom_errorbar(aes(ymin = low, ymax = high, lty = p.value > 0.05), width = 0.3) + # draw error bars
    ggtitle("(b) Abies density") +
    geom_point() + 
    geom_line(aes(group = where)) + 
    scale_y_continuous(labels=scaleFUN) +
    facet_wrap(~size) + 
    theme_few() +  # come from package ggthemes
    scale_colour_hue(l=40) +
    theme(axis.text=element_text(size=10), 
          axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
          axis.title.y = element_text(colour = 'white'),
          strip.text.x = element_text(hjust = 0, size = 10, face = "bold"),
          plot.margin=unit(c(.35,0,.35,0), "cm"),
          plot.title = element_text(size = 11.5, face = "bold", hjust = 0),
          legend.position = "none",
          axis.title.x=element_blank())

plot(ab.abun.ggplot)







# model Rhodo mortality
# =============================

head(rhodo)
nrow(rhodo)
rm(rhodo.data)
rhodo.data <- rhodo[,c("AbiesTL","sizeclass", "freq_live", "DTreeLine", "AbiesTLabsDist", "Slope","Canopy","Stumps","aspect1","transect", "mortality")] %>% na.omit
head(rhodo.data)
nrow(rhodo.data)


# weight by live density so that 1 total individual that is dead in a plot with mortality of 1.0 does not have much influence

rm(fit.small.above, fit.small.below, fit.medium.above, fit.medium.below, fit.large.above, fit.large.below, summary.small, summary.medium, summary.large)

fit.small.above <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), 
                         family = binomial(), 
                         data = subset(rhodo.data, sizeclass == "0-1m" & AbiesTL == "above treeline"),
                         weights = freq_live)
fit.small.below <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), 
                         family = binomial(), 
                         data = subset(rhodo.data, sizeclass == "0-1m" & AbiesTL == "below treeline"),
                         weights = freq_live)

fit.medium.above <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), 
                          family = binomial(), 
                          data = subset(rhodo.data, sizeclass == "1-2m" & AbiesTL == "above treeline"),
                          weights = freq_live)
fit.medium.below <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), 
                          family = binomial(), 
                          data = subset(rhodo.data, sizeclass == "1-2m" & AbiesTL == "below treeline"),
                          weights = freq_live)

fit.large.above <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), 
                         family = binomial(), 
                         data = subset(rhodo.data, sizeclass == "2m+" & AbiesTL == "above treeline"),
                         weights = freq_live)
fit.large.below <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), 
                         family = binomial(), 
                         data = subset(rhodo.data, sizeclass == "2m+" & AbiesTL == "below treeline"),
                         weights = freq_live)



summary.small <- rbind(data.frame(where = "above treeline", tidy(fit.small.above)),
                       data.frame(where = "below treeline", tidy(fit.small.below))) %>% prepdf
summary.medium <- rbind(data.frame(where = "above treeline", tidy(fit.medium.above)),
                        data.frame(where = "below treeline", tidy(fit.medium.below))) %>% prepdf 
summary.large <- rbind(data.frame(where = "above treeline", tidy(fit.large.above)),
                       data.frame(where = "below treeline", tidy(fit.large.below))) %>% prepdf 

rm(rhodo.summary)
rhodo.summary <- rbind(data.frame(size = "0-1m", summary.small), 
                       data.frame(size = "1-2m", summary.medium), 
                       data.frame(size = "2m+", summary.large))
rhodo.summary
rhodo.summary <- namechange(rhodo.summary); rhodo.summary

bubble.df <- rbind(bubble.df, data.frame(response = "rhodo.mort", rhodo.summary)); bubble.df


ggplot(rhodo.summary, 
       aes(Term, estimate, col = where)) + 
    geom_hline(yintercept = 0, col="lightgrey") + # zero line
    geom_errorbar(aes(ymin = low, ymax = high, lty = p.value > 0.05), width = 0.3) + # draw error bars
    ggtitle("Rhododendron mortality") + 
    geom_point() + 
    geom_line(aes(group = where)) + 
    facet_wrap(~size) + 
    theme_few() +  # come from package ggthemes
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_hue(l=40)


rh.mort.ggplot <- ggplot(rhodo.summary, 
                         aes(Term, estimate, col = where)) + 
    geom_hline(yintercept = 0, col="lightgrey") + # zero line
    geom_errorbar(aes(ymin = low, ymax = high, lty = p.value > 0.05), width = 0.3) + # draw error bars
    ggtitle("(c) Rhododendron mortality") +
    geom_point() + 
    geom_line(aes(group = where)) + 
    scale_y_continuous(labels=scaleFUN) +
    facet_wrap(~size) + 
    theme_few() +  # come from package ggthemes
    scale_colour_hue(l=40) +
    theme(axis.text=element_text(size=10), 
          axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
          axis.title = element_text(size = 10.5),
          strip.text.x = element_text(hjust = 0, size = 10, face = "bold"),
          plot.margin=unit(c(.35,0,.35,0), "cm"),
          plot.title = element_text(size = 11.5, face = "bold", hjust = 0),
          legend.position = "none",
          axis.title.x=element_blank())

plot(rh.mort.ggplot)








# model Abies mortality
# =============================

head(abies)
nrow(abies)
rm(abies.data)
abies.data <- abies[,c("AbiesTL","sizeclass", "freq_dead", "freq_live", "DTreeLine", "AbiesTLabsDist", "Slope","Canopy","Stumps","aspect1","transect", "mortality")] %>% na.omit
head(abies.data)
nrow(abies.data)



rm(fit.small.above, fit.small.below, fit.medium.above, fit.medium.below, fit.large.above, fit.large.below, summary.small, summary.medium, summary.large)

fit.small.above <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(aspect1)) + (1| transect), 
                         family = binomial(), 
                         data = subset(abies.data, sizeclass == "0-1m" & AbiesTL == "above treeline"),
                         weights = freq_live)

fit.small.below <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), 
                         family = binomial(), 
                         data = subset(abies.data, sizeclass == "0-1m" & AbiesTL == "below treeline"),
                         weights = freq_live)

fit.medium.above <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), 
                          family = binomial(), 
                          data = subset(abies.data, sizeclass == "1-2m" & AbiesTL == "above treeline"),
                          weights = freq_live)
fit.medium.below <- glmer(mortality ~ (scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), 
                          family = binomial(), 
                          data = subset(abies.data, sizeclass == "1-2m" & AbiesTL == "below treeline"),
                          weights = freq_live)

fit.large.above <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) +  scale(aspect1)) + (1| transect), 
                         family = binomial(), 
                         data = subset(abies.data, sizeclass == "2m+" & AbiesTL == "above treeline"),
                         weights = freq_live)
fit.large.below <- glmer(mortality ~ (scale(AbiesTLabsDist) + scale(Slope) + scale(Canopy) + scale(Stumps) + scale(aspect1)) + (1| transect), 
                         family = binomial(), 
                         data = subset(abies.data, sizeclass == "2m+" & AbiesTL == "below treeline"),
                         weights = freq_live)


summary.small <- rbind(data.frame(where = "above treeline", tidy(fit.small.above)),
                       data.frame(where = "below treeline", tidy(fit.small.below))) %>% prepdf
summary.medium <- rbind(data.frame(where = "below treeline", tidy(fit.medium.below))) %>% prepdf 
# no fit.medium.above due to only one point
summary.large <- rbind(data.frame(where = "above treeline", tidy(fit.large.above)),
                       data.frame(where = "below treeline", tidy(fit.large.below))) %>% prepdf 


rm(abies.summary)
abies.summary <- rbind(data.frame(size = "0-1m", summary.small), 
                       data.frame(size = "1-2m", summary.medium), 
                       data.frame(size = "2m+", summary.large))
abies.summary
abies.summary <- namechange(abies.summary); abies.summary

bubble.df <- rbind(bubble.df, data.frame(response = "abies.mort", abies.summary)); bubble.df


ggplot(abies.summary, 
       aes(Term, estimate, col = where)) + 
    geom_hline(yintercept = 0, col="lightgrey") + # zero line
    geom_errorbar(aes(ymin = low, ymax = high, lty = p.value > 0.05), width = 0.3) + # draw error bars
    ggtitle("Abies mortality") + 
    geom_point() + 
    geom_line(aes(group = where)) + 
    facet_wrap(~size) + 
    theme_few() +  # come from package ggthemes
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_hue(l=40)


ab.mort.ggplot <- ggplot(abies.summary, 
                         aes(Term, estimate, col = where)) + coord_cartesian(ylim = c(-7, 7)) +
    geom_hline(yintercept = 0, col="lightgrey") + # zero line
    geom_errorbar(aes(ymin = low, ymax = high, lty = p.value > 0.05), width = 0.3) + # draw error bars
    ggtitle("(d) Abies mortality") +          
    geom_point() + 
    geom_line(aes(group = where)) + 
    scale_y_continuous(labels=scaleFUN) +
    facet_wrap(~size) + 
    theme_few() +  # come from package ggthemes
    scale_colour_hue(l=40) +
    theme(axis.text=element_text(size=10), 
          axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
          axis.title = element_text(size = 10.5),
          axis.title.y = element_text(colour = 'white'),
          strip.text.x = element_text(hjust = 0, size = 10, face = "bold"),
          plot.margin=unit(c(.35,0,.35,0), "cm"),
          plot.title = element_text(size = 11.5, face = "bold", hjust = 0),
          legend.position = "none",
          axis.title.x=element_blank())

plot(ab.mort.ggplot)







# ===================================

# FIGURE 5 -- bubble plot

# ===================================


bubble.df
dim(bubble.df)
table(bubble.df$response)

bubble.df$estimate[bubble.df$p.value > 0.05] <- 0; bubble.df
bubble.df <- subset(bubble.df, select = c("response", "size", "where", "Term", "estimate"))
bubble.df

maxef <- max(bubble.df$estimate, na.rm = T); maxef
maxef <- 1.15    # scaling factor to fit circle in the space
myrectcol <- 'gray80'
myeffsizetextcol <- 'white'
order.df <- data.frame(Term = c("Distance\\ntreeline", "Canopy", "Slope", "Aspect", "Stumps"), order = 1:5); order.df



require(plotrix)
require(grid)
require(grDevices)
require(reshape)



getwd()

pdf("Figure_5_bubble_plot.pdf", width=5.95, height=5.8)
par(mfrow=c(1,2))
par(oma = c(0,0,0,0)) # make room  for the overall x and y axis titles
par(mar = c(0,0,0,0)) # make the plots be closer together



    # Rcamp density above treeline
    # -----------------------------------
    
    rm(plotdata)
    plotdata <- subset(bubble.df, response == "rhodo.ab" & where == "above treeline"); plotdata
    plotdata <- plotdata[,c(4,2,5)]; plotdata
    plotdata <- cast(plotdata, Term ~ size); plotdata
    
    plotdata <- merge(plotdata, order.df, by = "Term", all = T); plotdata
    plotdata <- plotdata[order(plotdata$order),]; plotdata
    
    row.names(plotdata) <- plotdata$Term; plotdata
    plotdata <- plotdata[,-c(1,5)]; plotdata
    plotdata <- plotdata[,3:1]; plotdata
    plotdata[is.na(plotdata)] <- 0
    
    par(ps = 8.5, cex = 1, cex.main = 1)
    plot(c(-1, 8), c(-12,6.25), type = "n", asp=1, axes=FALSE, xlab="", ylab="")
    # plot(c(-1, 8), c(-10,5.5),  asp=1, xlab="", ylab="")
    
    
    for (m in c(0, 1.25, 2.5)) {
        
        for (i in 1:5) {
            
            x1 <- i * 1.25
            x2 <- x1 + 1
            rect(x1, m, x2, m+1, col = myrectcol, border="transparent") 
            
            colnum <- trunc(m) + 1; colnum
            
            # draw.circle later plots circle based on radius; calculate radius based on area = effect size
            r <- sqrt(abs(plotdata[i,colnum]) / 3.14159)/maxef
            mycirccol <- ifelse(plotdata[i,colnum] > 0, 'blue', 'red')
            draw.circle((x1+x2)/2, (m+m+1)/2, r, col=mycirccol, border='transparent')
            
            #plot effect size 
            if (abs(plotdata[i,colnum]) > 0) {
                text((x1+x2)/2, (m+m+1)/2, labels = round(plotdata[i,colnum],2), col=myeffsizetextcol, cex=0.5)
            }
            
            if(m == 0) {
                text((x1+x2)/2, 4, cex=0.95, srt=90, adj=c(0,0.4), labels = row.names(plotdata)[i])
            }
            
        }    # close row of block loop
        
        text(x=0.5, y=(m+m+1)/2, labels = colnames(plotdata)[colnum], font=1, cex=1)
        
    }	# close size class loop
    
    
    m <- 1.25
    text(x=-0.75, y=(m+m+1)/2, labels = "Density\\nabove treeline", font=1, cex=1, srt=90)
    
    text(x=0.7, y=4.75, labels = "Explanatory\\nvariables", font=1, cex=1, adj=1)
    arrows(0.8, 4.75, 1.2, 4.75, length = 0.05, angle = 30, code = 2, col = "gray", lty = 1, lwd = 2)
    
    text(x=0.7, y=4, labels = "Size classes", font=1, cex=1, adj=1)
    arrows(0.5, 3.7, 0.5, 3.3, length = 0.05, angle = 30, code = 2, col = "gray", lty = 1, lwd = 2)
    
    text(x=3, y=6.25, cex=1, adj=0, "R. campanulatum", font=4) 
    text(x=1.25, y=6.25, cex=1, adj=0, "(a)", font=2) 
    
    
    
    
    # Rcamp density below treeline
    # -----------------------------------
    
    rm(plotdata)
    plotdata <- subset(bubble.df, response == "rhodo.ab" & where == "below treeline"); plotdata
    plotdata <- plotdata[,c(4,2,5)]; plotdata
    plotdata <- cast(plotdata, Term ~ size); plotdata
    
    plotdata <- merge(plotdata, order.df, by = "Term", all = T); plotdata
    plotdata <- plotdata[order(plotdata$order),]; plotdata
    
    row.names(plotdata) <- plotdata$Term; plotdata
    plotdata <- plotdata[,-c(1,5)]; plotdata
    plotdata <- plotdata[,3:1]; plotdata 
    plotdata[is.na(plotdata)] <- 0; plotdata 
    
    # par(ps = 8.5, cex = 1, cex.main = 1)
    # plot(c(-1, 8), c(-10,5.5), type = "n", asp=1, axes=FALSE, xlab="", ylab="")
    
    for (m in c(0, 1.25, 2.5)-4) {
        
        for (i in 1:5) {
            
            x1 <- i * 1.25
            x2 <- x1 + 1
            rect(x1, m, x2, m+1, col = myrectcol, border="transparent") 
            
            colnum <- trunc(m+4) + 1; colnum
            
            # draw.circle later plots circle based on radius; calculate radius based on area = effect size
            r <- sqrt(abs(plotdata[i,colnum]) / 3.14159)/maxef
            mycirccol <- ifelse(plotdata[i,colnum] > 0, 'blue', 'red')
            draw.circle((x1+x2)/2, (m+m+1)/2, r, col=mycirccol, border='transparent')
            
            #plot effect size 
            if (abs(plotdata[i,colnum]) > 0) {
                text((x1+x2)/2, (m+m+1)/2, labels = round(plotdata[i,colnum],2), col=myeffsizetextcol, cex=0.5)
            }
            
            if(m == 0) {
                text((x1+x2)/2, 4, cex=0.95, srt=90, adj=c(0,0.25), labels = row.names(plotdata)[i])
            }
            
        }    # close row of block loop
        
        text(x=0.5, y=(m+m+1)/2, labels = colnames(plotdata)[colnum], font=1, cex=1)
        
    }	# close size class loop
    
    
    m <- 1.25-4
    text(x=-0.75, y=(m+m+1)/2, labels = "Density\\nbelow treeline", font=1, cex=1, srt=90)
    
    
    
    
    # Rcamp mortality above treeline
    # -----------------------------------
    
    rm(plotdata)
    plotdata <- subset(bubble.df, response == "rhodo.mort" & where == "above treeline"); plotdata
    plotdata <- plotdata[,c(4,2,5)]; plotdata
    plotdata <- cast(plotdata, Term ~ size); plotdata
    
    plotdata <- merge(plotdata, order.df, by = "Term", all = T); plotdata
    plotdata <- plotdata[order(plotdata$order),]; plotdata
    
    row.names(plotdata) <- plotdata$Term; plotdata
    plotdata <- plotdata[,-c(1,5)]; plotdata
    plotdata <- plotdata[,3:1]; plotdata 
    plotdata[is.na(plotdata)] <- 0; plotdata 
    
    
    for (m in c(0, 1.25, 2.5)-8) {
        
        for (i in 1:5) {
            
            x1 <- i * 1.25
            x2 <- x1 + 1
            rect(x1, m, x2, m+1, col = myrectcol, border="transparent") 
            
            colnum <- trunc(m+8) + 1; colnum
            
            # draw.circle later plots circle based on radius; calculate radius based on area = effect size
            r <- sqrt(abs(plotdata[i,colnum]) / 3.14159)/maxef
            mycirccol <- ifelse(plotdata[i,colnum] > 0, 'blue', 'red')
            draw.circle((x1+x2)/2, (m+m+1)/2, r, col=mycirccol, border='transparent')
            
            #plot effect size 
            if (abs(plotdata[i,colnum]) > 0) {
                text((x1+x2)/2, (m+m+1)/2, labels = round(plotdata[i,colnum],2), col=myeffsizetextcol, cex=0.5)
            }
            
            if(m == 0) {
                text((x1+x2)/2, 4, cex=0.95, srt=90, adj=c(0,0.25), labels = row.names(plotdata)[i])
            }
            
        }    # close row of block loop
        
        text(x=0.5, y=(m+m+1)/2, labels = colnames(plotdata)[colnum], font=1, cex=1)
        
    }	# close size class loop
    
    
    m <- 1.25-8
    text(x=-0.75, y=(m+m+1)/2, labels = "Mortality\\nabove treeline", font=1, cex=1, srt=90)
    
    
    
    # Rcamp mortality below treeline
    # -----------------------------------
    
    rm(plotdata)
    plotdata <- subset(bubble.df, response == "rhodo.mort" & where == "below treeline"); plotdata
    plotdata <- plotdata[,c(4,2,5)]; plotdata
    plotdata <- cast(plotdata, Term ~ size); plotdata
    
    plotdata <- merge(plotdata, order.df, by = "Term", all = T); plotdata
    plotdata <- plotdata[order(plotdata$order),]; plotdata
    
    row.names(plotdata) <- plotdata$Term; plotdata
    plotdata <- plotdata[,-c(1,5)]; plotdata
    plotdata <- plotdata[,3:1]; plotdata 
    plotdata[is.na(plotdata)] <- 0; plotdata 
    
    
    for (m in c(0, 1.25, 2.5)-12) {
        
        for (i in 1:5) {
            
            x1 <- i * 1.25
            x2 <- x1 + 1
            rect(x1, m, x2, m+1, col = myrectcol, border="transparent") 
            
            colnum <- trunc(m+12) + 1; colnum
            
            # draw.circle later plots circle based on radius; calculate radius based on area = effect size
            r <- sqrt(abs(plotdata[i,colnum]) / 3.14159)/maxef
            mycirccol <- ifelse(plotdata[i,colnum] > 0, 'blue', 'red')
            draw.circle((x1+x2)/2, (m+m+1)/2, r, col=mycirccol, border='transparent')
            
            #plot effect size 
            if (abs(plotdata[i,colnum]) > 0) {
                text((x1+x2)/2, (m+m+1)/2, labels = round(plotdata[i,colnum],2), col=myeffsizetextcol, cex=0.5)
            }
            
            if(m == 0) {
                text((x1+x2)/2, 4, cex=0.95, srt=90, adj=c(0,0.25), labels = row.names(plotdata)[i])
            }
            
        }    # close row of block loop
        
        text(x=0.5, y=(m+m+1)/2, labels = colnames(plotdata)[colnum], font=1, cex=1)
        
    }	# close size class loop
    
    
    m <- 1.25-12
    text(x=-0.75, y=(m+m+1)/2, labels = "Mortality\\nbelow treeline", font=1, cex=1, srt=90)
    
    
    

    
    # Abies density above treeline
    # -----------------------------------
    
    rm(plotdata)
    plotdata <- subset(bubble.df, response == "abies.ab" & where == "above treeline"); plotdata
    plotdata <- plotdata[,c(4,2,5)]; plotdata
    plotdata <- cast(plotdata, Term ~ size); plotdata
    
    plotdata <- merge(plotdata, order.df, by = "Term", all = T); plotdata
    plotdata <- plotdata[order(plotdata$order),]; plotdata
    
    row.names(plotdata) <- plotdata$Term; plotdata
    plotdata <- plotdata[,-c(1,5)]; plotdata
    plotdata <- plotdata[,3:1]; plotdata 
    plotdata[is.na(plotdata)] <- 0; plotdata 
    
    
    par(ps = 8.5, cex = 1, cex.main = 1)
    plot(c(4, 8), c(-12,6.25), type = "n", asp=1, axes=FALSE, xlab="", ylab="")
    # plot(c(-1, 8), c(-10,5.5),  asp=1, xlab="", ylab="")
    
    
    for (m in c(0, 1.25, 2.5)) {
        
        for (i in 1:5) {
            
            x1 <- i * 1.25
            x2 <- x1 + 1
            rect(x1, m, x2, m+1, col = myrectcol, border="transparent") 
            
            colnum <- trunc(m) + 1; colnum
            
            # draw.circle later plots circle based on radius; calculate radius based on area = effect size
            r <- sqrt(abs(plotdata[i,colnum]) / 3.14159)/maxef
            mycirccol <- ifelse(plotdata[i,colnum] > 0, 'blue', 'red')
            draw.circle((x1+x2)/2, (m+m+1)/2, r, col=mycirccol, border='transparent')
            
            #plot effect size 
            if (abs(plotdata[i,colnum]) > 0) {
                text((x1+x2)/2, (m+m+1)/2, labels = round(plotdata[i,colnum],2), col=myeffsizetextcol, cex=0.5)
            }
            
            if(m == 0) {
                text((x1+x2)/2, 4, cex=0.95, srt=90, adj=c(0,0.4), labels = row.names(plotdata)[i])
            }
            
        }    # close row of block loop
    }	# close size class loop
    
    
    text(x=3, y=6.25, cex=1, adj=0, "A. spectabilis", font=4) 
    text(x=1.25, y=6.25, cex=1, adj=0, "(b)", font=2) 
    
    
    
    
    # Abies density below treeline
    # -----------------------------------
    
    rm(plotdata)
    plotdata <- subset(bubble.df, response == "abies.ab" & where == "below treeline"); plotdata
    plotdata <- plotdata[,c(4,2,5)]; plotdata
    plotdata <- cast(plotdata, Term ~ size); plotdata
    
    plotdata <- merge(plotdata, order.df, by = "Term", all = T); plotdata
    plotdata <- plotdata[order(plotdata$order),]; plotdata
    
    row.names(plotdata) <- plotdata$Term; plotdata
    plotdata <- plotdata[,-c(1,5)]; plotdata
    plotdata <- plotdata[,3:1]; plotdata 
    plotdata[is.na(plotdata)] <- 0; plotdata 
    
    
    for (m in c(0, 1.25, 2.5)-4) {
        
        for (i in 1:5) {
            
            x1 <- i * 1.25
            x2 <- x1 + 1
            rect(x1, m, x2, m+1, col = myrectcol, border="transparent") 
            
            colnum <- trunc(m+4) + 1; colnum
            
            # draw.circle later plots circle based on radius; calculate radius based on area = effect size
            r <- sqrt(abs(plotdata[i,colnum]) / 3.14159)/maxef
            mycirccol <- ifelse(plotdata[i,colnum] > 0, 'blue', 'red')
            draw.circle((x1+x2)/2, (m+m+1)/2, r, col=mycirccol, border='transparent')
            
            #plot effect size 
            if (abs(plotdata[i,colnum]) > 0) {
                text((x1+x2)/2, (m+m+1)/2, labels = round(plotdata[i,colnum],2), col=myeffsizetextcol, cex=0.5)
            }
            
            if(m == 0) {
                text((x1+x2)/2, 4, cex=0.95, srt=90, adj=c(0,0.25), labels = row.names(plotdata)[i])
            }
            
        }    # close row of block loop
    }	# close size class loop
    
    
    
    
    
    # Abies mortality above treeline
    # -----------------------------------
    
    rm(plotdata)
    plotdata <- subset(bubble.df, response == "abies.mort" & where == "above treeline"); plotdata
    plotdata <- plotdata[,c(4,2,5)]; plotdata
    plotdata <- rbind(plotdata, plotdata[nrow(plotdata),]); plotdata
    # manually add and entirely missing sizeclass
    plotdata[nrow(plotdata),2] <- "1-2m"; plotdata
    plotdata[nrow(plotdata),3] <- 0; plotdata
    plotdata <- cast(plotdata, Term ~ size); plotdata
    
    plotdata <- merge(plotdata, order.df, by = "Term", all = T); plotdata
    plotdata <- plotdata[order(plotdata$order),]; plotdata
    
    row.names(plotdata) <- plotdata$Term; plotdata
    plotdata <- plotdata[,-c(1,5)]; plotdata
    plotdata <- plotdata[,3:1]; plotdata 
    plotdata[is.na(plotdata)] <- 0; plotdata 
    
    
    for (m in c(0, 1.25, 2.5)-8) {
        
        for (i in 1:5) {
            
            x1 <- i * 1.25
            x2 <- x1 + 1
            rect(x1, m, x2, m+1, col = myrectcol, border="transparent") 
            
            colnum <- trunc(m+8) + 1; colnum
            
            # draw.circle later plots circle based on radius; calculate radius based on area = effect size
            r <- sqrt(abs(plotdata[i,colnum]) / 3.14159)/maxef
            mycirccol <- ifelse(plotdata[i,colnum] > 0, 'blue', 'red')
            draw.circle((x1+x2)/2, (m+m+1)/2, r, col=mycirccol, border='transparent')
            
            #plot effect size 
            if (abs(plotdata[i,colnum]) > 0) {
                text((x1+x2)/2, (m+m+1)/2, labels = round(plotdata[i,colnum],2), col=myeffsizetextcol, cex=0.5)
            }
            
            if(m == 0) {
                text((x1+x2)/2, 4, cex=0.95, srt=90, adj=c(0,0.25), labels = row.names(plotdata)[i])
            }
            
        }    # close row of block loop
    }	# close size class loop
    
    
    
    
    # Abies mortality below treeline
    # -----------------------------------
    
    rm(plotdata)
    plotdata <- subset(bubble.df, response == "abies.mort" & where == "below treeline"); plotdata
    plotdata <- plotdata[,c(4,2,5)]; plotdata
    plotdata <- cast(plotdata, Term ~ size); plotdata
    
    plotdata <- merge(plotdata, order.df, by = "Term", all = T); plotdata
    plotdata <- plotdata[order(plotdata$order),]; plotdata
    
    row.names(plotdata) <- plotdata$Term; plotdata
    plotdata <- plotdata[,-c(1,5)]; plotdata
    plotdata <- plotdata[,3:1]; plotdata 
    plotdata[is.na(plotdata)] <- 0; plotdata 
    
    
    for (m in c(0, 1.25, 2.5)-12) {
        
        for (i in 1:5) {
            
            x1 <- i * 1.25
            x2 <- x1 + 1
            rect(x1, m, x2, m+1, col = myrectcol, border="transparent") 
            
            colnum <- trunc(m+12) + 1; colnum
            
            # draw.circle later plots circle based on radius; calculate radius based on area = effect size
            r <- sqrt(abs(plotdata[i,colnum]) / 3.14159)/maxef
            mycirccol <- ifelse(plotdata[i,colnum] > 0, 'blue', 'red')
            draw.circle((x1+x2)/2, (m+m+1)/2, r, col=mycirccol, border='transparent')
            
            #plot effect size 
            if (abs(plotdata[i,colnum]) > 0) {
                text((x1+x2)/2, (m+m+1)/2, labels = round(plotdata[i,colnum],2), col=myeffsizetextcol, cex=0.5)
            }
            
            if(m == 0) {
                text((x1+x2)/2, 4, cex=0.95, srt=90, adj=c(0,0.25), labels = row.names(plotdata)[i])
            }
            
        }    # close row of block loop
    }	# close size class loop
    

dev.off()

