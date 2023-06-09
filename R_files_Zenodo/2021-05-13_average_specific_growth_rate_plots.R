
library(tidyverse)
library(patchwork)

# Processing the data -----------------------------------------------------

# Get list of files 
files <- list.files(path = "Data", pattern = "diuron_light_rep", full.names = TRUE)

# Read each of these files and combine them
xdiuron_light <- lapply(files, function(x) read.csv(x)) %>%
  bind_rows(.id = "id") %>% within(hours[hours == 0.3] <- 0.33)


# Processing the data in an abstract and concise way
xdiuron_light2 <- xdiuron_light  %>% 
  filter(Replicate == "AVERAGE", Diuron != "Blank", hours >= 0) %>%
  group_by(hours, Light, Diuron, id) %>% 
  summarise(Absorbance = mean(Absorbance), X = paste(X, collapse = " ")) %>% ungroup() %>% 
  within({
    
    Diuron_num <-  case_when(Diuron == "AlgaeControl" ~ 0,
                             Diuron == "MeOH" ~ 0,
                             TRUE ~ as.numeric(Diuron)
    )
    
    Light_num <-  as.numeric(Light)
    
    block <-  factor(id)
    
    hours_fact <- factor(hours)
    
    celld <- 4179.6 * Absorbance - 172.48
    t0 <- ifelse(hours == 0, celld, NA)
    
  }) %>% group_by(Light, Diuron, id ) %>% fill(t0) %>% 
  ungroup() %>%  
  mutate(sample_id = factor(gsub("."," ", x = paste(Diuron, Light, block), fixed = TRUE)),
         log_t0_change = log(celld/t0)) %>% 
  filter(hours != 0, Diuron != "MeOH", block != 1)


# compute the percent inhibition
xdiuron_light3 <- xdiuron_light2 %>% 
  mutate(control = ifelse((Diuron_num == 0 & Light == 80), celld, NA ),
         control_t0 = ifelse((Diuron_num == 0 & Light == 80), t0, NA ),
         control_log_t0_change = ifelse((Diuron_num == 0 & Light == 80), log_t0_change, NA )) %>% 
  group_by(hours_fact, block) %>% fill(control, control_t0, control_log_t0_change, .direction = "downup") %>% 
  group_by(hours_fact,  Light_num, Diuron_num) %>% 
  summarise( mean_percent_inhib = mean(((control - control_t0) - (celld - t0))/(control - control_t0) * 100),
             mean_percent_inhib_asg_rate = mean((1 - log_t0_change/mean(control_log_t0_change))*100),
             SE_mean = sqrt(var(((control - control_t0) - (celld - t0))/(control - control_t0) * 100, na.rm = TRUE))/sqrt(3),
             SE_mean_asg_rate = sqrt(var((1 - log_t0_change/
                                            mean(control_log_t0_change))*100, na.rm = TRUE))/sqrt(3)) %>% 
  ungroup() %>% filter(Light_num == 80, hours_fact %in% c("24", "48", "72"))
  

# Comparing wide vs long
data.frame(Response = 1, Treatment = 1)
data.frame(Tr1 = 1, Tr2 =1 )


# Only keep rows of data where light = 80 and block = 2
#xdiuron_light3 <- xdiuron_light2 %>% filter(Light_num == 80, block == 2)

# # Compute the mean and standard deviation
# xdiuron_light4 <- xdiuron_light2 %>% filter(Light_num == 80) %>% group_by(Diuron_num, hours_fact) %>% 
#   summarise(Mean_percent_inhib = mean(percent_inhib, na.rm = TRUE), 
#             SE_mean = sqrt(var(percent_inhib, na.rm = TRUE))/sqrt(3),
#             Mean_percent_inhib_asg_rate = mean(percent_inhib_asg_rate, na.rm = TRUE),
#             SE_mean_asg_rate = sqrt(var(percent_inhib_asg_rate, na.rm = TRUE))/sqrt(3),
#             ) %>%
#   filter(Diuron_num != 0)


# Plotting ----------------------------------------------------------------


# Bar plot
p1<- ggplot(data = filter(xdiuron_light3, Diuron_num != 0), # specify data argument
       mapping = aes(x = factor(Diuron_num), y =  mean_percent_inhib_asg_rate, fill = hours_fact)) + # Specify aesthetics 
  geom_col(position = position_dodge2()) + # Add a column graph that isn't stacked
  geom_errorbar(mapping = aes(ymin =  mean_percent_inhib_asg_rate - SE_mean_asg_rate, 
                              ymax =  mean_percent_inhib_asg_rate + SE_mean_asg_rate), 
                position = position_dodge2(padding = .4)) + # Add error bars, padding is bar width
  theme_classic() + # Change themes
  ylab("% Inhibition of Average Specific Growth Rate") + xlab("Diuron Concentration") +# Change labels
  labs(fill = "Hours elasped") # Change legend labels


p1



# library("drc")
# 
# 
# # Run a logistic regression
# m24 <- drm(mean_percent_inhib_asg_rate ~ Diuron_num, 
#           data = filter(xdiuron_light3, hours_fact == "24", Diuron_num != 0), 
#           fct = LL.3(names = c("hill", "max_value", "ec_50")))
# 
# m48 <- drm(mean_percent_inhib_asg_rate ~ Diuron_num, 
#            data = filter(xdiuron_light3, hours_fact == "48", Diuron_num != 0), 
#            fct = LL.3(names = c("hill", "max_value", "ec_50")))
# 
# m72 <- drm(mean_percent_inhib_asg_rate ~ Diuron_num, 
#            data = filter(xdiuron_light3, hours_fact == "72", Diuron_num != 0), 
#            fct = LL.3(names = c("hill", "max_value", "ec_50")))
# 
# summary(m24)

# make some predictions for plotting
Diuron_vals <- c(seq(0, 3, length.out =10000), 3)#c(.11, .33, 1, 3)
# logist_data <- data.frame(Diuron_vals, m24_pred = predict(m24, newdata = data.frame(Diuron_num = Diuron_vals)),
#              m48_pred = predict(m48, newdata = data.frame(Diuron_num = Diuron_vals)),
#              m72_pred = predict(m72, newdata = data.frame(Diuron_num = Diuron_vals)))

p2 <- ggplot(data = xdiuron_light3, aes(x = Diuron_num, y = mean_percent_inhib_asg_rate)) +
  geom_point(aes(colour = hours_fact)) + 
  geom_smooth(aes(colour = hours_fact, fill = hours_fact), method="loess") + 
  geom_abline(intercept = 0, lty = "dashed", slope = 0) +
  #geom_line(data = logist_data, mapping = aes(x =Diuron_vals, y = m24_pred), colour = "#F8766D") +
  #geom_line(data = logist_data, mapping = aes(x =Diuron_vals, y = m48_pred), colour = "#00BA38") +
  #geom_line(data = logist_data, mapping = aes(x =Diuron_vals, y = m72_pred), colour = "#619CFF") +
  #scale_x_log10() +
  theme_classic() + labs(colour = "Hours elasped", fill  = "Hours elasped") +
  ylab("% Inhibition of Average Specific Growth Rate") + 
  xlab("Diuron Concentration")

p2 






