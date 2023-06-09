library('tidyverse')
library('cowplot')
library('wesanderson')
library('lme4')
library('car')

sampleData <- read.csv(file = "sampleData.csv")

dataGrouped <- group_by(sampleData, ecosystemBurn) %>%
  summarize(mycorrhizaSD = sd(mycorrhiza),
            mycorrhiza = mean(mycorrhiza),
            richnessSD = sd(richness),
            richness = mean(richness),
            pathogenSD = sd(pathogen),
            pathogen = mean(pathogen),
            saproSD = sd(sapro),
            sapro = mean(sapro),
            mycorrhizaPropAbund = mean(mycorrhizaPropAbund),
            saprotrophPropAbund = mean(saprotrophPropAbund),
            pathogenPropAbund = mean(pathogenPropAbund),
            ecosystem = first(ecosystem),
            burned = first(burned))

glmer(richness ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% Anova()
glmer(richness ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% lsmeans::lsmeans(pairwise~burned*ecosystem, adjust="tukey") %>%
  emmeans::CLD(Letters=letters)
dataGrouped$richnessGroups <- c("a", "b", "ab", "b")

glmer(mycorrhiza ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% Anova(type = 3)
glmer(mycorrhiza ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% lsmeans::lsmeans(pairwise~burned*ecosystem, adjust="tukey") %>%
  emmeans::CLD(Letters=letters)
dataGrouped$mycorrhizaGroups <- c("a", "b", "b", "b")

glmer(sapro ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% Anova()
glmer(sapro ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% lsmeans::lsmeans(pairwise~burned*ecosystem, adjust="tukey") %>%
  emmeans::CLD(Letters=letters)
dataGrouped$saprotrophGroups <- c("a", "b", "ab", "b")

glmer(pathogen ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% Anova()
glmer(pathogen ~ burned*ecosystem + (1 | specSite), data = sampleData,
      family = "poisson") %>% lsmeans::lsmeans(pairwise~burned*ecosystem, adjust="tukey") %>%
  emmeans::CLD(Letters=letters)
dataGrouped$pathogenGroups <- c("a", "a", "a", "a")

richnessPlot <- plot_grid(
  ggplot(dataGrouped, aes(x = ecosystem, y = richness, color = burned, shape = ecosystem, label = richnessGroups)) +
    geom_errorbar(aes(ymin=richness-richnessSD, ymax=richness+richnessSD),
                  width=0, 
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(x = ecosystem, y=richness), size = 3, 
               position = position_dodge(width = 0.75)) +
    geom_text(color = "black", position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c(wes_palette("Rushmore1")[5], wes_palette("Rushmore1")[4])) +
    scale_x_discrete(labels = c("Evergreen forest", "Oak woodland")) +
    scale_shape_manual(values = c(17,16)) +
    theme_minimal() +
    xlab("") +
    ylab("Total richness") +
    ylim(100, 550) +
    theme(legend.position = "NONE"),
  ggplot(dataGrouped, aes(x = ecosystem, shape = ecosystem, y = mycorrhiza, color = burned, label = mycorrhizaGroups)) +
    geom_errorbar(aes(ymin=mycorrhiza-mycorrhizaSD, ymax=mycorrhiza+mycorrhizaSD),
                  width=0, 
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(x = ecosystem, y=mycorrhiza), size = 3, 
               position = position_dodge(width = 0.75)) +
    geom_text(color = "black", position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c(wes_palette("Rushmore1")[5], wes_palette("Rushmore1")[4])) +
    scale_x_discrete(labels = c("Evergreen forest", "Oak woodland")) +
    scale_shape_manual(values = c(17,16)) +
    theme_minimal() +
    xlab("") +
    ylab("Mycorrhiza") +
    ylim(5, 50) +
    theme(legend.position = "NONE"),
  ggplot(dataGrouped, aes(x = ecosystem, shape = ecosystem, y = sapro, color = burned, label = saprotrophGroups)) +
    geom_errorbar(aes(ymin=sapro-saproSD, ymax=sapro+saproSD),
                  width=0, 
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(x = ecosystem, y=sapro), size = 3, 
               position = position_dodge(width = 0.75)) +
    geom_text(color = "black", position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c(wes_palette("Rushmore1")[5], wes_palette("Rushmore1")[4])) +
    scale_x_discrete(labels = c("Evergreen forest", "Oak woodland")) +
    scale_shape_manual(values = c(17,16)) +
    theme_minimal() +
    xlab("") +
    ylab("Saprotrophs") +
    ylim(50, 200) +
    theme(legend.position = "NONE"),
  ggplot(dataGrouped, aes(x = ecosystem, shape = ecosystem, y = pathogen, color = burned, label = pathogenGroups)) +
    geom_errorbar(aes(ymin=pathogen-pathogenSD, ymax=pathogen+pathogenSD),
                  width=0, 
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(x = ecosystem, y=pathogen), size = 3, 
               position = position_dodge(width = 0.75)) +
    geom_text(color = "black", position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c(wes_palette("Rushmore1")[5], wes_palette("Rushmore1")[4])) +
    scale_x_discrete(labels = c("Evergreen forest", "Oak woodland")) +
    scale_shape_manual(values = c(17,16)) +
    theme_minimal() +
    xlab("") +
    ylab("Pathogens") +
    ylim(15, 60) +
    theme(legend.position = "NONE"),
  ncol = 2,
  labels = c("a","b","c","d")
)

richnessPlot %>%
  save_plot(file = "richnessPlot.svg", plot = .,
            base_width = 6,
            base_height = 5)

dataGrouped$otherPropAbund <- 1 - dataGrouped$mycorrhizaPropAbund - 
  dataGrouped$saprotrophPropAbund -
  dataGrouped$pathogenPropAbund

propAbundFrame <- dplyr::select(dataGrouped, ecosystemBurn, mycorrhizaPropAbund, saprotrophPropAbund,
                                pathogenPropAbund, otherPropAbund) %>% pivot_longer(-ecosystemBurn, names_to = "proportions") %>%
  left_join(dataGrouped, by = "ecosystemBurn")

propAbundFrame$ecosystemBurn <- propAbundFrame$ecosystemBurn %>% str_replace_all("DF", "Evergreen\\nforest,\\n") %>%
  str_replace_all("OW", "Oak\\nwoodland,\\n") %>%
  str_replace_all("B", "burned") %>%
  str_replace_all("U", "unburned")

propAbundFrame$proportions <- propAbundFrame$proportions %>% str_replace_all("mycorrhizaPropAbund", "Mycorrhiza") %>%
  str_replace_all("saprotrophPropAbund", "Saprotroph") %>%
  str_replace_all("pathogenPropAbund", "Pathogen") %>%
  str_replace_all("otherPropAbund", "Other")

lmer(mycorrhizaPropAbund ~ burned*ecosystem + (1 | specSite), data = sampleData) %>% 
  car::Anova(type = 3)
lmer(saprotrophPropAbund ~ burned*ecosystem + (1 | specSite), data = sampleData) %>% 
  car::Anova()
lmer(pathogenPropAbund ~ burned*ecosystem + (1 | specSite), data = sampleData) %>% 
  car::Anova()

propAbundPlot <- ggplot(propAbundFrame, aes(fill = proportions, y = value, x = ecosystemBurn)) +
  geom_bar(position = "fill", stat= "identity") +
  scale_fill_manual(values = c(wes_palette("Rushmore1")[3], 
                                "#bababa",
                                wes_palette("Rushmore1")[4], 
                                wes_palette("Rushmore1")[1])) +
  ylab("Proportion of sequences") +
  xlab("") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank())

ggsave(file = "Figure3.pdf",
       propAbundPlot,
       width = 5,
       height = 3)
