
library(tidyverse)
library(patchwork)
library(visreg)

load(file = "Data/Diuron_photo_model_and_data2.RDA") # Diuron photosynthesis models and data
load(file = "Data/DIN_photo_model_and_data2.RDA") # DIN photosynthesis models and data
load(file = "Data/Diuron_growth_model_and_data2.RDA") # Diuron growth models and data
load(file = "Data/DIN_growth_model_and_data2.RDA") # Diuron growth models and data


# I_R (interaction) plots -------------------------------------------------

### DIN growth I_R plot
p1 <- ggplot(DIN_growth_CI_data) + aes(x = hours_fact, y = X50., colour = Treatment, fill = Treatment) +
  geom_hline(yintercept = 0, lty= "dashed") +theme_classic() + xlab("Hours") +
  geom_point(position = position_dodge2(0.4)) +
  geom_linerange(mapping = aes(ymin = X2.5., ymax = X97.5.), position = position_dodge2(0.4)) +
  ylab(expression("I"[R])) + guides(alpha = FALSE) +
  labs(title = "A") #add A and B labels

### Diuron growth I_R plot
p2 <- ggplot(Diuron_growth_CI_data) + aes(x = hours_fact, y = X50., colour = Treatment, fill = Treatment) + 
  geom_hline(yintercept = 0, lty= "dashed") +theme_classic() + xlab("Hours") +
  geom_point(position = position_dodge2(0.5)) + 
  geom_linerange(mapping = aes(ymin = X2.5., ymax = X97.5.), position = position_dodge2(0.5)) +
  ylab(expression("I"[R])) + guides(alpha = FALSE) +
  labs(title = "B") #add A and B labels

### DIN photosynthesis I_R plot
p3<- ggplot(DIN_photo_CI_data) + aes(x = hours_fact, y = X50., colour = Treatment, fill = Treatment) + 
  geom_hline(yintercept = 0, lty= "dashed") +theme_classic() + xlab("Hours") +
  geom_point(position = position_dodge2(0.4)) +
  geom_linerange(mapping = aes(ymin = X2.5., ymax = X97.5.), position = position_dodge2(0.4)) +
  ylab(expression("I"[R])) +guides(alpha = FALSE) +
  labs(title = "A") #add A and B labels

### Diuron photosynthesis I_R plot
p4 <- ggplot(Diuron_photo_CI_data) + aes(x = hours_fact, y = X50., colour = Treatment, fill = Treatment) + 
  geom_hline(yintercept = 0, lty= "dashed") +theme_classic() + xlab("Hours") +
  geom_point(position = position_dodge2(0.5)) +
  geom_linerange(mapping = aes(ymin = X2.5., ymax = X97.5.), position = position_dodge2(0.5)) +
  ylab(expression("I"[R])) +guides(alpha = FALSE) +
  labs(title = "B") #add A and B labels

# IR plots for growth -----------------------------------------------------

# creates the two data sets (DIN x light and diuron x light )
df1 <- select(DIN_growth_CI_data, -log10DIN, -DIN_num) %>% mutate(Data_ID = "DIN")
df2 <- select(Diuron_growth_CI_data, -Diuron_num) %>% mutate(Data_ID = "Diuron")
plot_DIN_growth <- df1
plot_Diuron_growth <-df2

# plots the diuron growth data, alone
ggplot(plot_DIN_growth) + aes(x = hours_fact, y = X50., colour = Treatment, fill = Treatment) +
  geom_hline(yintercept = 0, lty= "dashed") +theme_bw() + xlab("Hours") +
  geom_point(position = position_dodge2(0.4))  +
  geom_linerange(mapping = aes(ymin = X2.5., ymax = X97.5.), position = position_dodge2(0.4)) +
  ylab(expression("I"[R])) + guides(alpha = FALSE)

#plots the DIN growth data, alone
ggplot(plot_DIN_growth) + aes(x = hours_fact, y = X50., colour = Treatment, fill = Treatment) +
  geom_hline(yintercept = 0, lty= "dashed") +theme_bw() + xlab("Hours") +
  geom_point(position = position_dodge2(0.4))  +
  geom_linerange(mapping = aes(ymin = X2.5., ymax = X97.5.), position = position_dodge2(0.4)) +
  ylab(expression("I"[R])) + guides(alpha = FALSE)


#stick to two plots together using patchwork (with separate legends)
(p1 + p2)

# save the plot into output folder
ggsave("Outputs/growth_DIN_Diuron_IR_2022-02-10.png", 
       width = 30, height = 13, units = "cm", dpi = 300 )


# IR plots for photosynthetic inhibition ----------------------------------

# combines the two data sets (DIN x light + diuron x light )
df3 <- select(DIN_photo_CI_data, -DIN_num) %>% mutate(Data_ID = "DIN")
df4 <- select(Diuron_photo_CI_data, -Diuron_num) %>% mutate(Data_ID = "Diuron")
plot_DIN_photo <- df3
plot_Diuron_photo <-df4

# plots the diuron photosynthesis data, alone
ggplot(plot_DIN_photo) + aes(x = hours_fact, y = X50., colour = Treatment, fill = Treatment) +
  geom_hline(yintercept = 0, lty= "dashed") +theme_bw() + xlab("Hours") +
  geom_point(position = position_dodge2(0.4))  +
  geom_linerange(mapping = aes(ymin = X2.5., ymax = X97.5.), position = position_dodge2(0.4)) +
  ylab(expression("I"[R])) + guides(alpha = FALSE)

#plots the DIN photosynthesis data, alone
ggplot(plot_DIN_photo) + aes(x = hours_fact, y = X50., colour = Treatment, fill = Treatment) +
  geom_hline(yintercept = 0, lty= "dashed") +theme_bw() + xlab("Hours") +
  geom_point(position = position_dodge2(0.4))  +
  geom_linerange(mapping = aes(ymin = X2.5., ymax = X97.5.), position = position_dodge2(0.4)) +
  ylab(expression("I"[R])) + guides(alpha = FALSE)

#stick to two plots together using patchwork (with separate legends)
(p3 + p4)

# save the plot into output folder
ggsave("Outputs/photosynthesis_DIN_Diuron_IR_2022-02-10.png", 
       width = 30, height = 13, units = "cm", dpi = 300 )
