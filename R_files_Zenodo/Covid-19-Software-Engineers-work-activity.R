# Data analysis Wave 1

library(lavaan)
library(semPlot)
library(car)
library(foreign)
#library(ggplot2)
#library(ggstatsplot)
library(psych)
library(effsize)
df <- read.spss("Wave1 and 2 activity project.sav", use.missings=TRUE, use.value.labels = F, to.data.frame = TRUE, check.names = FALSE) # 'use.missings' kicks all the missing values as defined in the SPSS data file out
summary(df$Age)
sd(df$Age, na.rm = T)
table(df$Gender) # 1: Men, 2: Women
summary(df$Duration__in_seconds_)
sd(df$Duration__in_seconds_)


# Internal consistencies and factor creation -----
# Well-being (adapted Satisfaction with Life Scale)
df$Well_being <- rowMeans(df[,21:25], na.rm = TRUE) - 14
alpha(df[,21:25])$total$raw_alpha

# Productivity
df$Productivity <- (df$Productivity_1/df$Productivity_4) * ((df$Productivity_1.1 + 100)/100) 

# Stress
df$Stress_2 <- 6 - df$Stress_2
df$Stress_3 <- 6 - df$Stress_3
df$Stress <- rowMeans(df[,78:81], na.rm = T)
alpha(df[,78:81])$total$raw_alpha

# Boredom
df$Boredom <- rowMeans(df[,49:56]) - 10
alpha(df[,49:56])$total$raw_alpha

# Need for relatedness
df$psychological_needs_2 <- 6 - df$psychological_needs_2
df$psychological_needs_4 <- 6 - df$psychological_needs_4
df$psychological_needs_6 <- 6 - df$psychological_needs_6
df$Relatedness <- rowMeans(df[,57:62], na.rm = T)
alpha(df[,57:62])$total$raw_alpha

# Need for competence
df$psychological_needs_8 <- 6 - df$psychological_needs_8
df$psychological_needs_10 <- 6 - df$psychological_needs_10
df$psychological_needs_12 <- 6 - df$psychological_needs_12
df$Competence <- rowMeans(df[,63:68], na.rm = T)
alpha(df[,63:68])$total$raw_alpha

# Need for autonomy
df$psychological_needs_14 <- 6 - df$psychological_needs_14
df$psychological_needs_16 <- 6 - df$psychological_needs_16
df$psychological_needs_18 <- 6 - df$psychological_needs_18
df$Autonomy <- rowMeans(df[,69:74], na.rm = T)
alpha(df[,69:74])$total$raw_alpha

# Communication with colleagues and line managers
df$Communication <- rowMeans(df[,75:77], na.rm = T)
alpha(df[,75:77])$total$raw_alpha

# Daily Routines
df$Daily_routines <- rowMeans(df[,82:84], na.rm = T)
alpha(df[,82:84])$total$raw_alpha

# Distractions at home
df$Distractions_at_home_1 <- df$Distractions_at_home_1 - 19
df$Distractions_at_home_2 <- 25 - df$Distractions_at_home_2
df$Distractions <- rowMeans(df[,85:86], na.rm = T)
alpha(df[,85:86])$total$raw_alpha



# wave2 ==========
# Internal consistencies and factor creation -----

# Well-beingW2 (adapted Satisfaction with Life Scale)
df$Well_beingW2 <- rowMeans(df[,127:131], na.rm = TRUE) - 14
alpha(df[,127:131])$total$raw_alpha
cor(df$Well_being, df$Well_beingW2, use = "p") # Test-retest reliability

# ProductivityW2
df$ProductivityW2 <- (df$Productivity_12/df$Productivity_42) * ((df$Productivity_1.12 + 100)/100)
cor(df$Productivity, df$ProductivityW2, use = "p")

# StressW2
df$Stress_22 <- 6 - df$Stress_22
df$Stress_32 <- 6 - df$Stress_32
df$StressW2 <- rowMeans(df[,188:191], na.rm = T)
alpha(df[,188:191])$total$raw_alpha
cor(df$Stress, df$StressW2, use = "p")

# BoredomW2
df$BoredomW2 <- rowMeans(df[,155:162]) - 10
alpha(df[,155:162])$total$raw_alpha
cor(df$Boredom, df$BoredomW2, use = "p")

# Need for relatednessW2
df$psychological_needs_22 <- 6 - df$psychological_needs_22
df$psychological_needs_42 <- 6 - df$psychological_needs_42
df$psychological_needs_62 <- 6 - df$psychological_needs_62
df$RelatednessW2 <- rowMeans(df[,167:172], na.rm = T)
alpha(df[,167:172])$total$raw_alpha
cor(df$Relatedness, df$RelatednessW2, use = "p")

# Need for competenceW2
df$psychological_needs_82 <- 6 - df$psychological_needs_82
df$psychological_needs_102 <- 6 - df$psychological_needs_102
df$psychological_needs_122 <- 6 - df$psychological_needs_122
df$CompetenceW2 <- rowMeans(df[,173:178], na.rm = T)
alpha(df[,173:178])$total$raw_alpha
cor(df$Competence, df$CompetenceW2, use = "p")

# Need for autonomyW2
df$psychological_needs_142 <- 6 - df$psychological_needs_142
df$psychological_needs_162 <- 6 - df$psychological_needs_162
df$psychological_needs_182 <- 6 - df$psychological_needs_182
df$AutonomyW2 <- rowMeans(df[,179:184], na.rm = T)
alpha(df[,179:184])$total$raw_alpha
cor(df$Autonomy, df$AutonomyW2, use = "p")

# Communication with colleagues and line managersW2
df$CommunicationW2 <- rowMeans(df[,185:187], na.rm = T)
alpha(df[,185:187])$total$raw_alpha
cor(df$Communication, df$CommunicationW2, use = "p")

# Daily RoutinesW2
df$Daily_routinesW2 <- rowMeans(df[,192:194], na.rm = T)
alpha(df[,192:194])$total$raw_alpha
cor(df$Daily_routines, df$Daily_routinesW2, use = "p")

# Distractions at homeW2
df$Distractions_at_home_12 <- df$Distractions_at_home_12 - 19
df$Distractions_at_home_22 <- 25 - df$Distractions_at_home_22
df$DistractionsW2 <- rowMeans(df[,195:196], na.rm = T)
alpha(df[,195:196])$total$raw_alpha
cor(df$Distractions, df$DistractionsW2, use = "p")






# Associations between task content and all other variables ####
# Normalise tasks for wave 1 (no need to do this for wave 2)
df[,30:44] <- (df[,30:44] / rowSums(df[,30:44], na.rm = T))*100

which( colnames(df)=="Q53_152")

cm_time <- cor(df[,c(30:44, 201:210)], use = "p")
# cm_time <- cor(df[,c(30:44, 201:210)], use = "p", method = "spearman")
round(cm_time[16:25,1:15], 2)
round(cm_time[1:15,1:15], 2) # are the tasks independent?
summary(as.numeric(round(cm_time[1:15,1:15], 2)))

cm_time2 <- cor(df[,c(136:150, 211:220)], use = "p", method = "spearman")
round(cm_time2[16:25,1:15], 2)
round(cm_time2[1:15,1:15], 2) # are the tasks independent?
summary(as.numeric(round(cm_time2[1:15,1:15], 2)))

diff <- df[,136:150] - df[,30:44]


 
### Pre-post comparisons of activities/tasks ####
# Loop
df1 <- df[,c(30:44, 136:150)]
ex <- 1:15
x <- matrix(data = NA, nrow = max(ex), ncol = 10)
for(i in ex){
  t <- t.test(df1[,i], df1[,i+15], paired = T)
  t.value <- round(as.numeric(t[1]),3)
  p.value <- round(as.numeric(t[3]),4)
  m1 <- round(mean(df1[,i], na.rm = T),3)
  m2 <- round(mean(df1[,i + 15], na.rm = T),3)
  sd1 <- round(sd(df1[,i], na.rm = T),3)
  sd2 <- round(sd(df1[,i + 15], na.rm = T),3)
  d <- round(as.numeric(effsize::cohen.d(df1[,i], df1[,i+15], paired = T, na.rm = T)[[3]]),3)
  Greater <- length(which(df1[,i] < df1[,i+15])) # How many participants score higher at t2 than t1?
  Lower <- length(which(df1[,i] > df1[,i+15]))
  Equal <- length(which(df1[,i] == df1[,i+15]))
  
  x[i,] <- c(m1, sd1, m2, sd2, t.value, p.value, d, Greater, Lower, Equal)
}
View(x[ex,])


Meyer <- c(17, 14, 8, 3, 5, 1, 15, 10, 4, 5, 2, 3, 2, 8, 3) # Meyer et al. results (typical workday). 
ex <- 1:15
x <- matrix(data = NA, nrow = max(ex), ncol = 7)
for(i in ex){
  t1 <- t.test(df1[,i], mu = Meyer[i])
  t2 <- t.test(df1[,i + 15], mu = Meyer[i])
  t.value1 <- round(as.numeric(t1[1]),3)
  t.value2 <- round(as.numeric(t2[1]),3)
  p.value1 <- round(as.numeric(t1[3]),4)
  p.value2 <- round(as.numeric(t2[3]),4)
  Mey <- Meyer[i]
  m1 <- round(mean(df1[,i], na.rm = T),3)
  m2 <- round(mean(df1[,i+15], na.rm = T),3)
  x[i,] <- c(Mey, m1, m2, t.value1, t.value2, p.value1, p.value2)
}  
View(x[ex,])  
t.test(df1[,7], mu = 15)

cor.test(Meyer, colMeans(df[,30:44], na.rm = T), use = "p")
cor.test(Meyer, colMeans(df[,136:150], na.rm = T), use = "p")
cor.test(colMeans(df[,30:44], na.rm = T), colMeans(df[,136:150], na.rm = T), use = "p")



###
### Graph ####
library("dplyr")
library("tidyr")
library("viridis")
library("ggthemes")
library(ggplot2)

#Week 1 & 2

df <- df %>% 
  rename(
    Q53_1_coding = Q53_1,
    Q53_1_coding2 = Q53_1_2,
    Q53_2_bugfixing = Q53_2,
    Q53_2_bugfixing2 = Q53_22,
    Q53_3_testing = Q53_3,
    Q53_3_testing2 = Q53_32,
    Q53_4_specification = Q53_4,
    Q53_4_specification2 = Q53_42,
    Q53_5_reviewingcode = Q53_5,
    Q53_5_reviewingcode2 = Q53_52,
    Q53_6_writedocumentation = Q53_6,
    Q53_6_writedocumentation2 = Q53_62,
    Q53_7_meetings = Q53_7,
    Q53_7_meetings2 = Q53_72,
    Q53_8_email = Q53_8,
    Q53_8_email2 = Q53_82,
    Q53_9_interruptions = Q53_9,
    Q53_9_interruptions2 = Q53_92,
    Q53_10_helping = Q53_10,
    Q53_10_helping2 = Q53_102,
    Q53_11_networking = Q53_11,
    Q53_11_networking2 = Q53_112,
    Q53_12_learning = Q53_12,
    Q53_12_learning2 = Q53_122,
    Q53_13_administration = Q53_13,
    Q53_13_administration2 = Q53_132,
    Q53_14_breaks = Q53_14,
    Q53_14_breaks2 = Q53_142,
    Q53_15_various = Q53_15,
    Q53_15_various2 = Q53_152
  )

which( colnames(df_tasks)=="Q53_1_coding")
which( colnames(df_tasks)=="Q53_15_various")

tasks_w1 <- df_tasks %>%
  select(ResponseId, 30:44) %>%
  pivot_longer(
    cols = starts_with("Q53_"))

which( colnames(df_tasks)=="Q53_1_coding2")
which( colnames(df_tasks)=="Q53_15_various2")

tasks_w2 <- df_tasks %>%
  select(ResponseId, 136:150) %>%
  pivot_longer(
    cols = starts_with("Q53_"))

tasks_w1$name <- sub('Q53_[[:digit:]]+_', '', tasks_w1$name)
tasks_w2$name <- sub('Q53_[[:digit:]]+_', '', tasks_w2$name)

tasks_w1 <- tasks_w1 %>%
  group_by(ResponseId) %>%
  mutate(total = sum(value)) %>%
  mutate(value_corrected = value / total * 100) #normalise tasks

# Normalisation not necessary for W2

#ggplot(tasks_w1, aes(fill = name, y = ResponseId, x = value_corrected)) +
#  geom_bar(stat = "identity")

tasks_w1 <- tasks_w1 %>%
  group_by(name) %>%
  summarise(mean = mean(value_corrected, na.rm = TRUE)) %>%
  arrange(mean) %>%
  mutate(name = factor(name, unique(name)))

tasks_w2 <- tasks_w2 %>%
  group_by(name) %>%
  summarise(mean = mean(value, na.rm = TRUE)) %>%
  arrange(mean) %>%
  mutate(name = factor(name, unique(name)))

ggplot(tasks_w1, aes(y = "", x = mean,
                 fill = reorder(name, mean))) +
  scale_fill_viridis(discrete = T)+
  guides(fill = guide_legend(reverse = TRUE)) +
  geom_bar(position="fill", stat="identity") +
  xlab("Task distribution W1") +
  theme_clean(base_size = 17) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_rect(color = NA),
        plot.background = element_rect(color = NA),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.title.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_x_continuous(labels = scales::percent)


# Merge!
tasks_w1$wave = "Wave 1"
tasks_w2$wave = "Wave 2"

tasks_w2$name <- sub('2', '', tasks_w2$name)

tasks <- rbind(tasks_w1, tasks_w2)

ggplot(tasks, aes(y = "", x = mean,
                     fill = reorder(name, mean))) +
  scale_fill_viridis(discrete = T)+
  guides(fill = guide_legend(reverse = TRUE)) +
  geom_bar(position="fill", stat="identity") +
  facet_wrap(~ wave, dir = "v") +
  xlab("Task distribution") +
  theme_clean(base_size = 17) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_rect(color = NA),
        plot.background = element_rect(color = NA),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.title.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_x_continuous(labels = scales::percent)



library(stringr) 
tasks$name <- str_to_title(tasks$name)
tasks[tasks$name == "Writedocumentation",]$name <- "Documentation"
tasks[tasks$name == "Reviewingcode",]$name <- "Code review"


unique(tasks$name)

library("tibble")

tasks <- add_row(tasks, name = "Coding", mean = 17, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Bugfixing", mean = 14, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Testing", mean = 8, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Specification", mean = 4, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Code review", mean = 5, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Documentation", mean = 1, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Meetings", mean = 15, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Email", mean = 10, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Interruptions", mean = 4, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Helping", mean = 5, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Networking", mean = 2, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Learning", mean = 3, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Administration", mean = 2, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Breaks", mean = 8, wave = "Meyer et al. - Typical workday")
tasks <- add_row(tasks, name = "Various", mean = 3, wave = "Meyer et al. - Typical workday")

tasks <- as.data.frame(tasks)
tasks$wave <- as.factor(tasks$wave)
str(tasks$wave)

tasks$wave <- factor(tasks$wave, levels = c("Wave 1", "Wave 2", "Meyer et al. - Typical workday"))


ggplot(tasks, aes(x = reorder(name, -mean), y = mean / 100, 
                  label = paste(round(mean,1),"%", sep = ""),
                  fill = reorder(name, -mean))) +
  guides(fill = guide_legend(reverse = TRUE)) +
  geom_histogram(stat="identity") +
  geom_text(y = tasks$mean/100 + 0.015) +
  facet_wrap(~ wave, dir = "v") +
  scale_y_continuous(limits = c(0,0.22), labels = c("0%", "5%", "10%", "15%", "20%"), breaks = c(0,0.05,0.1,0.15,0.2)) +
  theme_clean(base_size = 17) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(color = NA),
        plot.background = element_rect(color = NA),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.title.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.x =element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("tasks.pdf")




