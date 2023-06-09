# Packages ----------------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(likert)

# Prepare the data --------------------------------------------------------

# Identify the sheets in the Excel file
data_sheets <- excel_sheets('data/all_categories.xlsx')

# Extract the data from each
data_sheets %>%
  map(function(sheet){ # iterate through each sheet name
    assign(x = sheet,
           value = read_xlsx(path = 'data/all_categories.xlsx',
                             sheet = sheet),
           envir = .GlobalEnv)
  })

q_levels <- c('Never', 'Rarely', 'Does not apply to me',
              'Sometimes', 'Often')

# Q4 ----------------------------------------------------------------------

q4 <- as.data.frame(Q4) %>%
  mutate(role = if_else(role == "I am research data steward",
                        'Handle others data',
                        if_else(role == "I'm actively involved in research and manage my own data",
                                'Handle own data',
                                'Handle others data')))

for(i in seq_along(q4[1:12])) {
  q4[,i] <- factor(q4[ ,i], levels = q_levels)
}

# Question 4 overall ------------------------------------------------------

lgr <- likert(q4[1:12])

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white",
                            colors = c("#CC3333", "#FF6600", "#FFFF00",
                                       "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 50,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = paste0('n = ', nrow(q4)))

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q4_overall.png", lgr_plot,
#       width = 30, height = 40, units = "cm",
#       dpi = 300)


# Question 4 --------------------------------------------------------------

lgr <- likert(q4[1:12], grouping = q4$role)

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                    "#009900", "#006600"),
                plot.percent.low = TRUE, plot.percent.high = TRUE,
                plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                text.color = "black", centered = TRUE,
                include.center = TRUE, ordered = TRUE,
                 wrap = 200,
                legend = "Response", legend.position = "bottom", panel.arrange = "v",
                panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = 'Handle own data: n = 46\\nHandle others data: n = 16')

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q4.png", lgr_plot,
#        width = 30, height = 60, units = "cm",
#        dpi = 300)

# Q5 ----------------------------------------------------------------------

q5 <- as.data.frame(Q5) %>%
  mutate(role = if_else(role == "I am research data steward",
                        'Handle others data',
                        if_else(role == "I'm actively involved in research and manage my own data",
                                'Handle own data',
                                'Handle others data')))

for(i in seq_along(q5[1:5])) {
  q5[,i] <- factor(q5[ ,i], levels = q_levels)
}

# Question 5 overall ------------------------------------------------------

lgr <- likert(q5[1:5])

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 40,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = paste0('n = ', nrow(q5)))

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q5_overall.png", lgr_plot,
#        width = 30, height = 20, units = "cm",
#        dpi = 300)

# Question 5 --------------------------------------------------------------

lgr <- likert(q5[1:5], grouping = q5$role)

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 200,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = 'Handle own data: n = 46\\nHandle others data: n = 16')

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q5.png", lgr_plot,
#        width = 30, height = 30, units = "cm",
#        dpi = 300)

# Q6 ----------------------------------------------------------------------

q6 <- as.data.frame(Q6) %>%
  mutate(role = if_else(role == "I am research data steward",
                        'Handle others data',
                        if_else(role == "I'm actively involved in research and manage my own data",
                                'Handle own data',
                                'Handle others data')))

for(i in seq_along(q6[1:7])) {
  q6[,i] <- factor(q6[ ,i], levels = q_levels)
}

# Question 6 overall ------------------------------------------------------

lgr <- likert(q6[1:7])

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 40,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = paste0('n = ', nrow(q6)))

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q6_overall.png", lgr_plot,
#        width = 30, height = 25, units = "cm",
#        dpi = 300)

# Question 6 --------------------------------------------------------------

lgr <- likert(q6[1:7], grouping = q6$role)

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 200,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = 'Handle own data: n = 46\\nHandle others data: n = 16')

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q6.png", lgr_plot,
#        width = 30, height = 30, units = "cm",
#        dpi = 300)

# Q7 ----------------------------------------------------------------------

q7 <- as.data.frame(Q7) %>%
  mutate(role = if_else(role == "I am research data steward",
                        'Handle others data',
                        if_else(role == "I'm actively involved in research and manage my own data",
                                'Handle own data',
                                'Handle others data')))

for(i in seq_along(q7[1:4])) {
  q7[,i] <- factor(q7[ ,i], levels = q_levels)
}

# Question 7 overall ------------------------------------------------------

lgr <- likert(q7[1:4])

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 40,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = paste0('n = ', nrow(q7)))

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

ggsave("plots/q7_overall.png", lgr_plot,
       width = 30, height = 30, units = "cm",
       dpi = 300)

# Question 7 --------------------------------------------------------------

lgr <- likert(q7[1:4], grouping = q7$role)

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 200,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = 'Handle own data: n = 46\\nHandle others data: n = 16')

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q7b.png", lgr_plot,
#        width = 30, height = 20, units = "cm",
#        dpi = 300)

# Q7 ----------------------------------------------------------------------

q7 <- as.data.frame(Q7) %>%
  mutate(role = if_else(role == "I am research data steward",
                        'Handle others data',
                        if_else(role == "I'm actively involved in research and manage my own data",
                                'Handle own data',
                                'Handle others data')))

for(i in seq_along(q7[1:4])) {
  q7[,i] <- factor(q7[ ,i], levels = q_levels)
}

# Q8 ----------------------------------------------------------------------
# The levels change from this question on

q_levels <- c('Not important', 'Low importance', 'Neutral',
              'Important', 'Very important')

q8 <- as.data.frame(Q8) %>%
  mutate(role = if_else(role == "I am research data steward",
                        'Handle others data',
                        if_else(role == "I'm actively involved in research and manage my own data",
                                'Handle own data',
                                if_else(role == "I'm actively involved in research but don't directly manage data",
                                        "Don't directly manage data",
                                'Handle others data'))))

for(i in seq_along(q8[1:12])) {
  q8[,i] <- factor(q8[ ,i], levels = q_levels)
}

# Question 8 overall ------------------------------------------------------

lgr <- likert(q8[1:12])

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 40,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = paste0('n = ', nrow(q8)))

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q8_overall.png", lgr_plot,
#        width = 30, height = 40, units = "cm",
#        dpi = 300)

# Question 8 --------------------------------------------------------------

lgr <- likert(q8[1:12], grouping = q8$role)

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 200,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = "Handle own data: n = 46\\nHandle others data: n = 16\\nDon't directly manage data = 9")

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q8b.png", lgr_plot,
#        width = 30, height = 60, units = "cm",
#        dpi = 300)

# Q9 ----------------------------------------------------------------------

q9 <- as.data.frame(Q9) %>%
  mutate(role = if_else(role == "I am research data steward",
                        'Handle others data',
                        if_else(role == "I'm actively involved in research and manage my own data",
                                'Handle own data',
                                if_else(role == "I'm actively involved in research but don't directly manage data",
                                        "Don't directly manage data",
                                        'Handle others data'))))

for(i in seq_along(q9[1:5])) {
  q9[,i] <- factor(q9[ ,i], levels = q_levels)
}

# Question 9 overall ------------------------------------------------------

lgr <- likert(q9[1:5])

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 40,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = paste0('n = ', nrow(q9)))

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q9_overall.png", lgr_plot,
#        width = 30, height = 20, units = "cm",
#        dpi = 300)

# Question 9 --------------------------------------------------------------

lgr <- likert(q9[1:5], grouping = q9$role)

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 200,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = "Handle own data: n = 46\\nHandle others data: n = 16\\nDon't directly manage data = 9")

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q9b.png", lgr_plot,
#        width = 30, height = 30, units = "cm",
#        dpi = 300)

# Q10 ----------------------------------------------------------------------

q10 <- as.data.frame(Q10) %>%
  mutate(role = if_else(role == "I am research data steward",
                        'Handle others data',
                        if_else(role == "I'm actively involved in research and manage my own data",
                                'Handle own data',
                                if_else(role == "I'm actively involved in research but don't directly manage data",
                                        "Don't directly manage data",
                                        'Handle others data'))))

for(i in seq_along(q10[1:6])) {
  q10[,i] <- factor(q10[ ,i], levels = q_levels)
}

# Question 10 overall ------------------------------------------------------

lgr <- likert(q10[1:6])

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 40,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = paste0('n = ', nrow(q10)))

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q10_overall.png", lgr_plot,
#        width = 30, height = 20, units = "cm",
#        dpi = 300)

# Question 10 --------------------------------------------------------------

lgr <- likert(q10[1:6], grouping = q10$role)

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 200,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = "Handle own data: n = 46\\nHandle others data: n = 16\\nDon't directly manage data = 9")

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q10b.png", lgr_plot,
#        width = 30, height = 30, units = "cm",
#        dpi = 300)

# Q11 ----------------------------------------------------------------------

q11 <- as.data.frame(Q11) %>%
  mutate(role = if_else(role == "I am research data steward",
                        'Handle others data',
                        if_else(role == "I'm actively involved in research and manage my own data",
                                'Handle own data',
                                if_else(role == "I'm actively involved in research but don't directly manage data",
                                        "Don't directly manage data",
                                        'Handle others data'))))

for(i in seq_along(q11[1:4])) {
  q11[,i] <- factor(q11[ ,i], levels = q_levels)
}

# Question 11 overall ------------------------------------------------------

lgr <- likert(q11[1:4])

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 40,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = paste0('n = ', nrow(q11)))

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q11_overall.png", lgr_plot,
#        width = 30, height = 20, units = "cm",
#        dpi = 300)

# Question 11 --------------------------------------------------------------

lgr <- likert(q11[1:4], grouping = q11$role)

lgr_plot <- likert.bar.plot(lgr, low.color = "#CC3333", high.color = "#006600",
                            neutral.color = "#FFFF00", neutral.color.ramp = "white", colors = c("#CC3333", "#FF6600", "#FFFF00",
                                                                                                "#009900", "#006600"),
                            plot.percent.low = TRUE, plot.percent.high = TRUE,
                            plot.percent.neutral = TRUE, plot.percents = FALSE, text.size = 3,
                            text.color = "black", centered = TRUE,
                            include.center = TRUE, ordered = TRUE,
                            wrap = 200,
                            legend = "Response", legend.position = "bottom", panel.arrange = "v",
                            panel.strip.color = "#F0F0F0") +
  theme_bw() +
  theme(axis.line.y = element_line(colour = 'black', size = 0.8),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.line.x = element_line(colour = 'black', size = 0.8),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0, size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') +
  labs(caption = "Handle own data: n = 46\\nHandle others data: n = 16\\nDon't directly manage data = 9")

lgr_plot$layers[[2]]$geom_params$width = 0.3
lgr_plot$layers[[3]]$geom_params$width = 0.3

# ggsave("plots/q11b.png", lgr_plot,
#        width = 30, height = 20, units = "cm",
#        dpi = 300)
