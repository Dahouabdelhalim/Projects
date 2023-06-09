library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

relat <- figure2data
summary(relat)
levels(relat$Ind) <- gsub(" ", "\\n", levels(relat$Ind))

pd <- position_dodge(width = 0.4)

ggplot(relat, aes(x = Ind, y = Value, shape = Log)) +
  geom_errorbar(aes(ymax = Value + SD, ymin = Value - SD),
                    width = 0.2, position = pd) +
  geom_point(size = 5, position = pd) +
  scale_shape_manual(values = c(19, 21)) +
  labs(y = "Pairwise Relatedness\\n", x = "") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 25),
        legend.title=element_blank(),
        legend.text = element_text(size = 15),
        legend.position = c(.2, .8),
        axis.line.x = element_line(), axis.line.y = element_line()) +
  scale_x_discrete(limits = c("Adults", "Adult\\nMales", "Adult\\nFemales", "Larvae")) +
  scale_y_continuous(limits = c(-0.12, 0.20),
                     breaks = seq(-0.12, 0.20, 0.04),
                     labels = seq(-0.12, 0.20, 0.04))