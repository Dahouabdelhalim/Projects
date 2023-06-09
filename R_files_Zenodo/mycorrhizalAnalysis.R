library(tidyverse)
library(lme4)
library(wesanderson)
library(cowplot)
library(car)

sampleData <- read.csv(file = "sampleData.csv")
taxaTable <- read.csv(file = "taxaTable.csv")

filter(taxaTable, funcGroup == "mycorrhiza")$sapro %>% length()
filter(taxaTable, funcGroup == "saprotroph")$sapro %>% length()
filter(taxaTable, funcGroup == "pathogen")$sapro %>% length()

taxaTable <- filter(taxaTable, funcGroup == "mycorrhiza")

sum(taxaTable$ecto)/sum(taxaTable$ecto %>% length())
sum(taxaTable$amf)/sum(taxaTable$ecto %>% length())

sum(taxaTable$amf)

dataGrouped <- group_by(sampleData, ecosystemBurn) %>%
  summarize(ectoSD = sd(ecto),
            ecto = mean(ecto),
            amfSD = sd(amf),
            amf = mean(amf),
            ecosystem = first(ecosystem),
            burned = first(burned))

glmer(ecto ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% Anova(type = 3)
glmer(amf ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% Anova(type = 3)

glmer(ecto ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% lsmeans::lsmeans(pairwise~burned*ecosystem, adjust="tukey") %>%
  emmeans::CLD(Letters=letters)
dataGrouped$ectoGroups <- c("a", "b", "b", "b")

glmer(amf ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% lsmeans::lsmeans(pairwise~burned*ecosystem, adjust="tukey") %>%
  emmeans::CLD(Letters=letters)
dataGrouped$amfGroups <- c("a", "b", "bc", "c")

mycPlot <- plot_grid(
  ggplot(dataGrouped, aes(x = ecosystem, shape = ecosystem, y = ecto, color = burned, label = ectoGroups)) +
    geom_errorbar(aes(ymin=ecto-ectoSD, ymax=ecto+ectoSD),
                  width=0, 
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(x = ecosystem, y=ecto), size = 3, 
               position = position_dodge(width = 0.75)) +
    geom_text(color = "black", position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c(wes_palette("Rushmore1")[5], wes_palette("Rushmore1")[4])) +
    scale_x_discrete(labels = c("Evergreen forest", "Oak woodland")) +
    scale_shape_manual(values = c(17,16)) +
    theme_minimal() +
    xlab("") +
    ylab("Ectomycorrhizal richness") +
    ylim(5,35) +
    theme(legend.position = "NONE"),
  ggplot(dataGrouped, aes(x = ecosystem, shape = ecosystem, y = amf, color = burned, label = amfGroups)) +
    geom_errorbar(aes(ymin=amf-amfSD, ymax=amf+amfSD),
                  width=0, 
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(x = ecosystem, y=amf), size = 3, 
               position = position_dodge(width = 0.75)) +
    geom_text(color = "black", position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c(wes_palette("Rushmore1")[5], wes_palette("Rushmore1")[4])) +
    scale_x_discrete(labels = c("Evergreen forest", "Oak woodland")) +
    scale_shape_manual(values = c(17,16)) +
    theme_minimal() +
    xlab("") +
    ylab("Arbuscular mycorrhizal richness") +
    ylim(0, 17) +
    theme(legend.position = "NONE"),
  ncol = 2,
  labels = c("a","b")
)

mycPlot %>%
  save_plot(file = "mycPlot.svg", plot = .,
            base_width = 6,
            base_height = 3)
