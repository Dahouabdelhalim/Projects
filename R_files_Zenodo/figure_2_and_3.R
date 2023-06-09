rm(list = ls())

setwd(".../Data_and_Code")
abies <- read.csv("abies.csv")
rhodo <- read.csv("rhododendron.csv")

head(abies)
head(rhodo)

rhodo <- mutate(rhodo, DTreeLine = (elevation - AbiesTLelev))
abies <- mutate(abies, DTreeLine = (elevation - AbiesTLelev))

levels(rhodo$sizeclass) <- paste(c("(a)", "(b)", "(c)"), levels(rhodo$sizeclass))
levels(abies$sizeclass) <- paste(c("(d)", "(e)", "(f)"), levels(abies$sizeclass))

pcks <- c("magrittr", "plyr", "ggplot2", "ggthemes", "knitr", "gridExtra")
sapply(pcks, require, char = TRUE)



# FIGURE 2 -- Abundance plot
# ==========================

rm(r.df, a.df)
head(rhodo)
nrow(rhodo)
r.df <- rhodo[!is.na(rhodo$freq_live),]
nrow(r.df)

head(abies)
nrow(abies)
a.df <- abies[!is.na(abies$freq_live),]
nrow(a.df)


p1 <- ggplot(r.df, aes(DTreeLine, freq_live)) + facet_wrap(~sizeclass, nrow  = 3) + geom_point(size=0.75) + stat_smooth() + geom_vline(xintercept = 0, lty = 3,  col = 'gray') + xlab("Meters above tree line") + ylab(expression(Density~(n/{100~m^2}))) + ggtitle("Rhododendron campanulatum") + theme_few()

p2 <- ggplot(a.df, aes(DTreeLine, freq_live)) + facet_wrap(~sizeclass, nrow  = 3, scales="free_y") + geom_point(size=0.75, col = 2) + stat_smooth() + geom_vline(xintercept = 0, lty = 3, col = 'gray') + xlab("Meters above tree line") + ylab("") + ggtitle("Abies spectabilis") + theme_few()

# save plot
pdf("Figure_2_RhodoAbiesAbundance.pdf", width=5.4, height=6, family='Helvetica')

grid.arrange(p1 + scale_x_continuous(breaks = seq(-200,350,50)) + 
                 scale_y_continuous(breaks = seq(-25,100,25)) + 
                 theme(axis.text=element_text(size=8),
                       axis.text.x=element_text(angle=45, hjust=1),
                       axis.title = element_text(size = 8.5),
                       strip.text.x = element_text(hjust = 0, size = 8, face = "bold"),
                       plot.margin=unit(c(.5,0,.5,0), "cm"),
                       plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5)), 
             p2 + scale_x_continuous(breaks = seq(-200,350,50)) + 
                 theme(axis.text=element_text(size=8), 
                       axis.text.x=element_text(angle=45, hjust=1),
                       axis.title = element_text(size = 8.5),
                       strip.text.x = element_text(hjust = 0, size = 8, face = "bold"),
                       plot.margin=unit(c(.5,0,.5,0), "cm"),
                       plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5)),
             ncol = 2)

dev.off()



# FIGURE 3 -- Mortality plot
# ==========================
# avoid plotting zero as freq_dead when freq_live is also zero

rm(r.df, a.df)
head(rhodo)
nrow(rhodo)
r.df <- rhodo[!is.na(rhodo$mortality),]
nrow(r.df)

head(abies)
nrow(abies)
a.df <- abies[!is.na(abies$mortality),]
nrow(a.df)

p1 <- ggplot(r.df, aes(DTreeLine, mortality)) + facet_wrap(~sizeclass, nrow  = 3) + geom_point(size=0.75) + stat_smooth() + geom_vline(xintercept = 0, lty = 3,  col = 'gray') + xlab("Meters above tree line") + ylab(expression(Mortality)) + ggtitle("Rhododendron campanulatum") + theme_few()

p2 <- ggplot(a.df, aes(DTreeLine, mortality)) + facet_wrap(~sizeclass, nrow  = 3, scales="free_y") + geom_point(size=0.75, col = 2) + stat_smooth() + geom_vline(xintercept = 0, lty = 3, col = 'gray') + xlab("Meters above tree line") + ylab("") + ggtitle("Abies spectabilis") + theme_few()

# save plot
pdf("Figure_3_RhodoAbiesMortality.pdf", width=5.4, height=6, family='Helvetica')

grid.arrange(p1 + scale_x_continuous(breaks = seq(-200,350,50)) + 
                 theme(axis.text=element_text(size=8),
                       axis.text.x=element_text(angle=45, hjust=1),
                       axis.title = element_text(size = 8.5),
                       strip.text.x = element_text(hjust = 0, size = 8, face = "bold"),
                       plot.margin=unit(c(.5,0,.5,0), "cm"),
                       plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5)), 
             p2 + scale_x_continuous(breaks = seq(-200,350,50)) + 
                 theme(axis.text=element_text(size=8), 
                       axis.text.x=element_text(angle=45, hjust=1),
                       axis.title = element_text(size = 8.5),
                       strip.text.x = element_text(hjust = 0, size = 8, face = "bold"),
                       plot.margin=unit(c(.5,0,.5,0), "cm"),
                       plot.title = element_text(size = 9, face = "bold.italic", hjust = 0.5)),
             ncol = 2)

dev.off()
