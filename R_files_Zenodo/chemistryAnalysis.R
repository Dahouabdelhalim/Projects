library('tidyverse')
library('cowplot')
library('wesanderson')

soilChem <- read.csv(file = "soilChemistry.csv")

dataGrouped <- filter(soilChem, is.na(pH) == FALSE) %>%
  group_by(ecosystemBurn) %>%
  summarize(cPercSD = sd(cPerc),
            cPerc = mean(cPerc),
            nPercSD = sd(nPerc),
            nPerc = mean(nPerc),
            pHSD = sd(pH),
            pH = mean(pH),
            burned.unburned = first(burned.unburned),
            ecosystem = first(ecosystem))

lme4::lmer(cPerc ~ burned.unburned*ecosystem + (1 | specSite), data = soilChem) %>% car::Anova()
lme4::lmer(nPerc ~ burned.unburned*ecosystem + (1 | specSite), data = soilChem) %>% car::Anova()
blme::blmer(pH ~ burned.unburned*ecosystem + (1 | specSite), data = soilChem,
            control = lme4::lmerControl(optimizer = "Nelder_Mead")) %>% car::Anova(type = 3)

lme4::lmer(cPerc ~ burned.unburned*ecosystem + (1 | specSite), data = soilChem) %>%
  lsmeans::lsmeans(pairwise~burned.unburned*ecosystem, adjust="tukey") %>%
  emmeans::CLD(Letters=letters)

lme4::lmer(nPerc ~ burned.unburned*ecosystem + (1 | specSite), data = soilChem) %>%
  lsmeans::lsmeans(pairwise~burned.unburned*ecosystem, adjust="tukey") %>%
  emmeans::CLD(Letters=letters)

blme::blmer(pH ~ burned.unburned*ecosystem + (1 | specSite), data = soilChem,
            control = lme4::lmerControl(optimizer = "Nelder_Mead")) %>%
  lsmeans::lsmeans(pairwise~burned.unburned*ecosystem, adjust="tukey") %>%
  emmeans::CLD(Letters=letters)

chemPlot <- plot_grid(
  ggplot(dataGrouped, aes(x = ecosystem, shape = ecosystem, y = cPerc, color = burned.unburned, label = c("a","a","a","a"))) +
    geom_errorbar(aes(ymin=cPerc-cPercSD, ymax=cPerc+cPercSD),
                  width=0, 
                  position = position_dodge(width = 0.75)) +
    geom_point(aes(x = ecosystem, y=cPerc), size = 3, 
               position = position_dodge(width = 0.75)) +
    geom_text(color = "black", position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c(wes_palette("Rushmore1")[5], wes_palette("Rushmore1")[4])) +
    scale_x_discrete(labels = c("Evergreen forest", "Oak woodland")) +
    scale_shape_manual(values = c(17,16)) +
    theme_minimal() +
    xlab("") +
    ylab("Soil C (%)") +
    ylim(2, 10) +
    theme(legend.position = "NONE"),
    ggplot(dataGrouped, aes(x = ecosystem, shape = ecosystem, y = nPerc, color = burned.unburned, label = c("a","a","a","a"))) +
      geom_errorbar(aes(ymin=nPerc-nPercSD, ymax=nPerc+nPercSD),
                    width=0, 
                    position = position_dodge(width = 0.75)) +
      geom_point(aes(x = ecosystem, y=nPerc), size = 3, 
                 position = position_dodge(width = 0.75)) +
      geom_text(color = "black", position = position_dodge(width = 0.75)) +
      scale_color_manual(values = c(wes_palette("Rushmore1")[5], wes_palette("Rushmore1")[4])) +
      scale_x_discrete(labels = c("Evergreen forest", "Oak woodland")) +
      scale_shape_manual(values = c(17,16)) +
      theme_minimal() +
      xlab("") +
      ylab("Soil N (%)") +
    ylim(0.1, 0.6) +
      theme(legend.position = "NONE"),
    ggplot(dataGrouped, aes(x = ecosystem, shape = ecosystem, y = pH, color = burned.unburned, label = c("b","a","a","a"))) +
      geom_errorbar(aes(ymin=pH-pHSD, ymax=pH+pHSD),
                    width=0, 
                    position = position_dodge(width = 0.75)) +
      geom_point(aes(x = ecosystem, y=pH), size = 3, 
                 position = position_dodge(width = 0.75)) +
      geom_text(color = "black", position = position_dodge(width = 0.75)) +
      scale_color_manual(values = c(wes_palette("Rushmore1")[5], wes_palette("Rushmore1")[4])) +
      scale_x_discrete(labels = c("Evergreen forest", "Oak woodland")) +
      scale_shape_manual(values = c(17,16)) +
      theme_minimal() +
      xlab("") +
      ylab("Soil pH") + 
    ylim(5, 7) +
      theme(legend.position = "NONE"),
  ncol = 3,
  labels = c("a","b","c")
)

chemPlot %>% save_plot(filename = "chemPlot.svg",
                       chemPlot, ncol = 3, nrow = 1,
                       base_height = 3, base_width = 3)
