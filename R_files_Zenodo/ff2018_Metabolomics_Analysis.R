
# 2018 Fed, 1-hour Fasted Hummingbird Metabolomics #

# Functions # 
sem.func <- function(x, y, z){
  sem.fed <- sd(x, na.rm=T)/sqrt(length(x))
  sem.fast <- sd(y, na.rm=T)/sqrt(length(y))
  return(c(rep(sem.fed, z), rep(sem.fast, z)))
}

# Libraries 
library("ggplot2")
library("ggpubr")
library("plyr")


### Fructose
fed.fru <- c(5.29,5.87,4.817,5.967,4.306,5.193, 5.943)
fast.fru <- c(0.9178, 0.03812, 0.0321, 0.0156, 0.06271, NA, NA)

fru.df <- data.frame(group = rep(c("Fed Fructose", "Fasted Fructose"), each = 7),
                     treatment = c(fed.fru, fast.fru),
                     sem = sem.func(fed.fru, fast.fru, 7))

fru.sig <- t.test(treatment ~ group, data = fru.df, var.equal = F, na.rm = T)
print(fru.sig)
# Welch Two Sample T-test: t = -17.253, df = 9.9044, p-value = 1.023e-08
# Mean Fasted Fructose: 0.213266 +/- 0.15 SEM
# Mean Fed Fructose: 5.340857 +/- 0.24 SEM

### Glucose
fed.glu <- c(33.11, 21, 34.08, 28.95, 29.68, 26.29, 37.16)
fast.glu <- c(31.62, 27.73, 28.65, 34.33, 26.01, NA, NA)

glu.df <- data.frame(group = rep(c("Fed Glucose", "Fasted Glucose"), each = 7),
                     treatment = c(fed.glu, fast.glu))

glu.sig <- t.test(treatment ~ group, data = glu.df, var.equal = F, na.rm = T)
print(glu.sig)
# Welch Two Sample T-test: t = -0.14743, df = 9.8799, p-value = 0.8858
# Mean Fasted Glucose: 29.7 +/- 1.25 SEM
# Mean Fed Glucose: 30.0  +/- 2.03 SEM

### Lactic Acid
fed.lac <- c(4468, 1588, 7828, 4115, 2082, 1704,8394)
fast.lac <- c(5346, 7397, 4528, 4906, 9567, NA, NA)
lac.df <- data.frame(group = rep(c("Fed Lactate", "Fasted Lactate"), each = 7),
                     treatment = c(fed.lac, fast.lac))

lac.sig <- t.test(treatment ~ group, data = lac.df, var.equal = F, na.rm = T)
print(lac.sig)
# Welch Two Sample t-test: t = 1.4257, df = 9.9428, p-value = 0.1846
# Mean Fasted Lactate: 6348.8 +/- 798.6 SEM, +/- 2112.9 SD
# Mean Fed Lactate: 4311.3 +/- 1072.2 SEM, +/- 2836.7 SD

### Pyruvic Acid
fed.pyr <- c(202.8, 100.5, 248.8, 186.3, 162.2, 189.7, 358.2)
fast.pyr <- c(211, 260.4, 212.3, 203.5, 224.2, NA, NA)
pyr.df <- data.frame(group = rep(c("Fed Pyruvate", "Fasted Pyruvate"), each = 7),
                     treatment = c(fed.pyr, fast.pyr))

pyr.sig <- t.test(treatment ~ group, data = pyr.df, var.equal = F, na.rm = T)
print(pyr.sig)
# Welch Two Sample t-test: t = 0.47989, df = 7.2661, p-value = 0.6454
# Mean Fasted Pyruvate: 222.3 +/- 8.6 SEM
# Mean Fed Pyruvate: 206.9 +/- 30.4 SEM


### Citric Acid 
fed.cit <- c(75.08, 7.315, 51.22, 45.66, 73.63, 27.09, 124.9)
fast.cit <- c(26.8, 22.03, 119.7, 134.9, 105.5, NA, NA)
cit.df <- data.frame(group = rep(c("Fed Citrate", "Fasted Citrate"), each = 7),
                     treatment = c(fed.cit, fast.cit))

cit.sig <- t.test(treatment ~ group, data = cit.df, var.equal = F, na.rm = T)
print(cit.sig)
boxplot(fed.cit, fast.cit)
# Welch Two Sample t-test: 
# Mean Fasted Citrate:  +/-  SEM, +/- SD
# Mean Fed Citrate:  +/-  SEM, +/- SD

### Malic Acid 
fed.mal <- c(30.75, 2.708, 39.26, 14.62, 21.27, 4.766, 41.22)
fast.mal <- c(16.52, 16.78, 20.99, 23.46, 37.58, NA, NA)
mal.df <- data.frame(group = rep(c("Fed Malate", "Fasted Malate"), each = 7),
                     treatment = c(fed.mal, fast.mal))

mal.sig <- t.test(treatment ~ group, data = mal.df, var.equal = F, na.rm = T)
print(mal.sig)
boxplot(fed.mal, fast.mal)
# Welch Two Sample t-test: 
# Mean Fasted malate:  +/-  SEM, +/- SD
# Mean Fed malate:  +/-  SEM, +/- SD


# ggBoxplots # 
# But they're actually barplots #

#scale_fill_manual(values = c("lightpink", "indianred2",
#                             "slategray1", "slategray3",
#                             "cornsilk", "burlywood1")) +   

fru.box <- ggplot(fru.df, (aes(x = group, y = treatment, fill = group))) + 
  geom_bar(stat = "summary", fun.y = "mean", color = "black") + 
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = .2) +
  labs(title = "Fructose") + 
  scale_x_discrete(name = "Treatment", 
                   labels = c("Fed", "Fasted"), 
                   limits = c("Fed Fructose", "Fasted Fructose")) + 
  scale_y_continuous(name = "Circulating Concentration (mM)", 
                     labels = scales::comma_format(accuracy=1)) + 
  scale_fill_manual(values = c("slategray1", "slategray3")) + 
  annotate("text", label = c("*"), x = c(1.5), y = 6.3) +
  geom_segment(aes(x = 1, xend = 2, y = 6.1, yend = 6.1)) + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.4, vjust = 0.75, size = 12),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = 0),
        panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = 0),
        axis.title.y = element_blank()) 

  

glu.box <- ggplot(glu.df, (aes(x = group, y = treatment, fill = group))) + 
  geom_bar(stat = "summary", fun.y = "mean", color = "black") + 
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = .2) +
 labs(title = "Glucose") +
  scale_x_discrete(name = "Treatment", 
                   labels = c("Fed", "Fasted"), 
                   limits = c("Fed Glucose", "Fasted Glucose")) + 
  scale_y_continuous(name = "Circulating Concentration (mM)", 
                     labels = scales::comma_format(accuracy=1)) + 
  scale_fill_manual(values = c("darkseagreen3", "darkseagreen4")) + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.4, vjust = 0.75, size = 12),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = 0),
        panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = 0)) 

lac.box <- ggplot(lac.df, (aes(x = group, y = treatment, fill = group))) + 
  geom_bar(stat = "summary", fun.y = "mean", color = "black") + 
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = .2) +
  labs(title = "Lactate") +
  scale_x_discrete(name = "Treatment", 
                   labels = c("Fed", "Fasted"), 
                   limits = c("Fed Lactate", "Fasted Lactate")) + 
  scale_y_continuous(name = "Concentration (mM)", 
                     labels = scales::comma_format(accuracy=1, scale=0.001)) + 
  scale_fill_manual(values = c("cornsilk", "burlywood1")) + 
  theme(legend.position = "none", 
        axis.text.y = element_text(angle = 0, hjust = 1),
        plot.title = element_text(hjust = 0.4, vjust = 0.75, size = 12),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = 0),
        panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = 0)) 

pyr.box <- ggplot(pyr.df, (aes(x = group, y = treatment, fill = group))) + 
  geom_bar(stat = "summary", fun.y = "mean", color = "black") + 
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = .2) +
 labs(title = "Pyruvate") +
  scale_x_discrete(name = "Treatment", 
                   labels = c("Fed", "Fasted"), 
                   limits = c("Fed Pyruvate", "Fasted Pyruvate")) + 
  scale_y_continuous(name = "Concentration (mM)", 
                     labels = scales::comma_format(accuracy=0.05, scale=0.001)) + 
  scale_fill_manual(values = c("lightpink", "indianred2")) + 
  theme(legend.position = "none", 
        axis.text.y = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.3, vjust = 0.75, size = 12),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = 0),
        panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = 0),
        axis.title.y = element_blank()) 



annotate_figure(ggarrange(glu.box, fru.box, lac.box, pyr.box, 
                                            ncol = 4, nrow = 1), 
                                  top = text_grob("Circulating Sugars and Metabolites", size = 15)) 


metabolite.fig <- ggarrange("", "", "", "","", glu.box, fru.box,"", lac.box, pyr.box,
          ncol = 5, nrow = 2, widths = c(1, 0.90, 0.05, 1, 1), heights = c(1, 10),
          labels = c("A) Glucose", "B) Fructose", "", "C) Lactate", "D) Pyruvate", "", "", "", ""),
          vjust = 2, hjust = c(-0.6, -0.25, 0, -0.6, -0.6),
          font.label = list(size = 12, family = "sans", face = "bold")) %>% 
  print()


metabolite.fig2 <- ggarrange(glu.box, fru.box,"", lac.box, pyr.box,
                            ncol = 5, nrow = 1, widths = c(1, 0.90, 0.05, 1, 1), heights = c(10),
                            labels = c("A)", "B)", "", "C)", "D)"),
#                            vjust = 1.64, hjust = c(-2.4, -0.9, 0, -2, -1.2),
                            vjust = 1.64, hjust = c(-2, -0.9, 0, -2, -1),
                            font.label = list(size = 13, family = "sans", face = "bold")) %>% 
  print()
# export to clipboard at 600x300


# ggsave(filename, plot = last_plot(), device = NULL, path = NULL,
#        scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
#        dpi = 300, limitsize = TRUE, ...)





#boxplots
dev.off()
dev.new()
par(mfrow = c(2, 3),
    mar = c(1, 1, 1, 1))
glu.box <- boxplot(fed.glu, fast.glu)
fru.box <- boxplot(fed.fru, fast.fru)
lac.box <- boxplot(fed.lac, fast.lac)
pyr.box <- boxplot(fed.pyr, fast.pyr)
cit.box <- boxplot(fed.cit, fast.cit)
mal.box <- boxplot(fed.mal, fast.mal)

#scatter plots
dev.off()
dev.new()
par(mfrow = c(2, 3),
    mar = c(1, 1, 1, 1))
glu.plot <- plot(c(fed.glu, fast.glu))
fru.plot <- plot(c(fed.fru, fast.fru))
lac.plot <- plot(c(fed.lac, fast.lac))
pyr.plot <- plot(c(fed.pyr, fast.pyr))
cit.plot <- plot(c(fed.cit, fast.cit))
mal.plot <- plot(c(fed.mal, fast.mal))




### Everything

fedfast.df <- setNames(data.frame (matrix (ncol = 2, nrow = 0)), c("group", "treatment"))
fedfast.df <- rbind(fedfast.df, glu.df, fru.df, lac.df, pyr.df, cit.df)
fedfast.df <- cbind(fedfast.df, condition = rep(c(rep("Fed", 7), rep("Fast", 7)), 5))

boxplot(fedfast.df)



fedfast.lm <- lm(treatment ~ group, data = fedfast.df)
fedfast.em <- emmeans(fedfast.lm, ~ group)
fedfast.mc <- contrast(fedfast.em, "tukey")

plot(fedfast.em)




fed.df <- data.frame(
  row.names = c("70", "85", "95", "96", "100", "101", "102"),
  fed.glu, fed.fru, fed.lac, fed.pyr, fed.cit)

fast.df <- data.frame(
  row.names = c("76", "86", "97", "98", "99", "x", "0"),
  fast.glu, fast.fru, fast.lac, fast.pyr, fast.cit)



dotchart(c(fed.glu, fast.glu, fed.fru, fast.fru, 
           fed.lac, fast.lac, fed.pyr, fast.pyr, 
           fed.cit, fast.cit))

ff.df <- data.frame(
  birdID = c("70", "85", "95", "96", "100", "101", "102",
                "76", "86", "97", "98", "99", "x", "0"),
  condition = c(rep("Fed", 7), rep("Fast", 7)),
  fed.glu, fed.fru, fed.lac, fed.pyr, fed.cit,
  fast.glu, fast.fru, fast.lac, fast.pyr, fast.cit)




library("ggplot2")
library("dplyr")
library("tidyverse")


ggplot(data = fed.df, 
       aes(x = c("fed.glu", "fast.glu"), y = c(fed.glu, fast.glu))) + 
  geom_bar(stat = 'identity')




fed.df %>%
ggplot(aes(colnames(fed.df[2:6]), fed.df[2:6])) + 
#  geom_boxplot()
  geom_line() + 
  geom_point() + 
  scale_y_continuous(sec.axis = sec_axis(~., name = ""))

ggplot() + 
  geom_point(position = "identity", data = fast.df, 
             aes(x = colnames(fast.df), 
                 y = fast.df))

ggplot() + 
  geom_point(position = "identity", data = fed.df, 
             aes(x = colnames(fed.df), 
                 y = c(fed.df[1:5])))
