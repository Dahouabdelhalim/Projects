
library(tidyverse)
library(mgcv)

source("Scripts/bayes_CI_functions.R")

# Retrieve the names of the din growth data files
files <- list.files(path = "Data",pattern = "DIN_light_rep", full.names = TRUE)

# Load in the files
xdin_light <- lapply(files, function(x) read.csv(x)) %>%
  bind_rows(.id = "id") %>% within(hours[hours == 0.3] <- 0.33)

# Process the data
xdin_light2 <- xdin_light  %>% 
  filter(Replicate == "AVERAGE", DIN != "Blank", hours >= 0) %>%
  group_by(hours, Light, DIN, id) %>% 
  summarise(Absorbance = mean(Absorbance), X = paste(X, collapse = " ")) %>% ungroup() %>% 
  within({
    
    DIN_num  <- case_when(
      DIN == "AlgaeControl" ~ 0,
      TRUE ~ as.numeric(DIN))
    
    fourRDIN <- (DIN_num)^(1/4)
    
    Light_num <-  as.numeric(Light)
    
    block <-  factor(id)
    
    hours_fact <- factor(hours)
    
    celld <- 4179.6 * Absorbance - 172.48
    t0 <- ifelse(hours == 0, celld, NA)
    
  }) %>% group_by(Light, DIN, id ) %>% fill(t0) %>% 
  ungroup() %>%  
  mutate(sample_id = factor(gsub("."," ", x = paste(DIN, Light, block), fixed = TRUE))) %>% 
  filter(hours != 0)

# Exploratory data analysis -----------------------------------------------

# cell density vs hours by din concentration
ggplot(xdin_light2) + 
  aes(x = as.numeric(hours), y = celld, color = factor(DIN_num)) + 
  geom_point() + 
  facet_grid(.~ Light_num, scales = "free") + 
  stat_smooth(se = FALSE) + 
  theme_classic()

# cell density vs hours by experimental block
ggplot(xdin_light2) + 
  aes(x = hours, y = celld, color = factor(block)) + 
  geom_point() + 
  facet_grid(.~ Light_num, scales = "free") + 
  stat_smooth(se = FALSE) + 
  theme_classic()


# Fit Generalised Additive Model ------------------------------------------

m1_din_light <- gam(log(celld) ~ 0 + hours_fact + offset(log(t0)) +
                      s(DIN_num, Light_num, k = 15, by = hours_fact) +  
                      s(block, bs = "re"), 
                    data = xdin_light2)

# Check assumptions
png(filename = "Outputs/DIN_GAM_residuals.png", width = 20, height = 20, units = "cm", res = 300)
par(mfrow = c(2,2))
gam.check(m1_din_light)
dev.off()

# Check model
summary(m1_din_light)


# Compute the IR data using empirical bayes method ------------------------

# Generate some data for prediction
xdin_newdata <- with(xdin_light2,{
  expand.grid(list(Light_num = unique(Light_num), DIN_num = unique(DIN_num), 
                   hours = unique(hours)))}) %>%
  mutate(block = 2, 
         sample_id = xdin_light2$sample_id[[1]], 
         hours_fact = factor(hours), t0 = 87) 

# Make predictions based on new data
set.seed(30)
df_bayes <- compute_CI_growth(model = m1_din_light, dat = xdin_newdata, nsims = 1000) %>% 
  within({
  Treatment <- case_when(Light_num == 80 & DIN_num == 0 ~ "Control",
                         Light_num == 80  ~ "DIN_only",
                         DIN_num == 0  ~ "Light_only",
                         TRUE ~ "Light_DIN")
}) 


# Generate the control matrix
df_bayes_control <- df_bayes %>% filter(Treatment == "Control")
control_mat <- df_bayes_control[ match(df_bayes$hours, df_bayes_control$hours) ,]

# Generate light only matrix
df_bayes_light <- df_bayes %>% filter(DIN_num == 0)
light_mat <- df_bayes_light[ match(paste(df_bayes$hours, df_bayes$Light_num), 
                                   paste(df_bayes_light$hours, df_bayes_light$Light_num)) ,]

# Generate DIN only matrix
df_bayes_DIN <- df_bayes %>% filter(Light_num == 80)
DIN_mat <- df_bayes_DIN[ match(paste(df_bayes$hours, df_bayes$DIN_num), 
                               paste(df_bayes_DIN$hours, df_bayes_DIN$DIN_num)) ,]

# Compute the IR matrix
IR_mat <- log((select(df_bayes, `1`:`1000`)/select(control_mat, `1`:`1000`))/
                ((select(light_mat, `1`:`1000`)/select(control_mat, `1`:`1000`))*
                   (select(DIN_mat, `1`:`1000`)/select(control_mat, `1`:`1000`)))) 

# Check for NAs
rowSums(is.na(IR_mat))

# Compute the credible intervals
DIN_growth_CI_data <- cbind(data.frame(t(apply(IR_mat, 1, quantile, 
                                               probs = c(0.025, 0.5, 0.975), na.rm = TRUE))),
                            xdin_newdata) %>% within({
  Treatment <- case_when(Light_num == 80 & DIN_num == 0 ~ "Control",
                         Light_num == 80  ~ "DIN_only",
                         DIN_num == 0  ~ "Light_only",
                         TRUE ~ paste0("Light = ", Light_num, ", DIN = ", DIN_num))}) %>% 
  filter(!(Treatment %in% c("Light_only","DIN_only",  "Control")), 
         !(DIN_num %in% c(0.0267, 0.267))) %>%
  select(-block, -sample_id)

# Take only what we need from raw data
DIN_growth_data <- xdin_light2 %>% 
  select(hours, hours_fact, Light_num, DIN_num, block, sample_id, t0, celld) %>% 
  within(celld_pred <- exp(predict(m1_din_light, newdata = .)))

# Save computations for plotting
save(DIN_growth_CI_data, DIN_growth_data, m1_din_light, 
     file = "Data/DIN_growth_model_and_data2.RDA")
