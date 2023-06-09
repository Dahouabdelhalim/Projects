
library(tidyverse)
library(mgcv)

source("Scripts/bayes_CI_functions.R")

# Retrieve names of Diuron growth files
files <- list.files(path = "Data", pattern = "diuron_light_rep", full.names = TRUE)

# Load in files
xdiuron_light <- lapply(files, function(x) read.csv(x)) %>%
  bind_rows(.id = "id") %>% within(hours[hours == 0.3] <- 0.33)

# Process the data
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
  mutate(sample_id = factor(gsub("."," ", x = paste(Diuron, Light, block), fixed = TRUE))) %>% 
  filter(hours != 0, Diuron != "MeOH") %>% 
  mutate(Diuron_fact = factor(ifelse(Diuron == "AlgaeControl", "0", Diuron)))


# Exploratory data analysis -----------------------------------------------

# Cell density vs hours by duiron concentration
ggplot(xdiuron_light2) + 
  aes(x = hours, y = celld, color = factor(Diuron_num)) + 
  geom_point() + 
  facet_grid(.~ Light_num, scales = "free") + 
  stat_smooth(se = FALSE) + 
  theme_classic()

# Cell density vs hours by experimental block
ggplot(xdiuron_light2) + 
  aes(x = hours, y = celld, color = factor(block)) + 
  geom_point() + 
  facet_grid(.~ Light_num, scales = "free") + 
  stat_smooth(se = FALSE) + 
  theme_classic()


# Fit Generalised Additive Model ------------------------------------------

m1_diuron_light <- gam(log(celld) ~ 0 + hours_fact + offset(log(t0)) +
                         s(Diuron_num, Light_num, k = 15, by = hours_fact) +
                         s(block, bs = "re"),
                       data = xdiuron_light2)

# Check assumptions
png(filename = "Outputs/Diuron_GAM_residuals.png", width = 20, height = 20, units = "cm", res = 300)
par(mfrow = c(2,2))
gam.check(m1_diuron_light)
dev.off()

# Check model
summary(m1_diuron_light)


# Compute the IR data using empirical bayes method ------------------------

# Generate some data for prediction
xdiuron_newdata <- with(xdiuron_light2,{
  expand.grid(list(Light_num = unique(Light_num), Diuron_num = unique(Diuron_num),
                   hours = unique(hours)))}) %>%
  mutate(block = 2, sample_id = xdiuron_light2$sample_id[[1]], hours_fact = factor(hours), t0 = 87) 

# Make predictions based 
set.seed(1000)
diuron_bayes <- compute_CI_growth(model = m1_diuron_light, dat = xdiuron_newdata, nsims = 1000) %>%
  within({
  Treatment <- case_when(Light_num == 80 & Diuron_num == 0 ~ "Control",
                         Light_num == 80  ~ "Diuron_only",
                         Diuron_num == 0  ~ "Light_only",
                         TRUE ~ "Light_DIN")
}) 


# Generate the control matrix
diuron_bayes_control <- diuron_bayes %>% filter(Treatment == "Control")
diuron_control_mat <- diuron_bayes_control[ match(diuron_bayes$hours, diuron_bayes_control$hours) ,]

# Generate light only matrix
diuron_bayes_light <- diuron_bayes %>% filter(Diuron_num == 0)
diuron_light_mat <- diuron_bayes_light[ 
  match(paste(diuron_bayes$hours, diuron_bayes$Light_num), 
        paste(diuron_bayes_light$hours, diuron_bayes_light$Light_num)) ,]

# Generate diuron only matrix
diuron_bayes_diuron <- diuron_bayes %>% filter(Light_num == 80 )
diuron_mat <- diuron_bayes_diuron[ 
  match(paste(diuron_bayes$hours, diuron_bayes$Diuron_num), 
        paste(diuron_bayes_diuron$hours, diuron_bayes_diuron$Diuron_num)) ,]

# Compute the IR matrix
IR_mat2 <- log((select(diuron_bayes, `1`:`1000`)/select(diuron_control_mat, `1`:`1000`))/
                 ((select(diuron_light_mat, `1`:`1000`)/select(diuron_control_mat, `1`:`1000`))*
                    (select(diuron_mat, `1`:`1000`)/select(diuron_control_mat, `1`:`1000`))))

# Compute the credible intervals
Diuron_growth_CI_data <- cbind(data.frame(t(apply(IR_mat2, 1, quantile, 
                                                  probs = c(0.025, 0.5, 0.975), na.rm = TRUE))),
                               xdiuron_newdata) %>% within({
  Treatment <- case_when(Light_num == 80 & Diuron_num == 0 ~ "Control",
                         Light_num == 80  ~ "Diuron_only",
                         Diuron_num == 0  ~ "Light_only",
                         TRUE ~ paste0("Light = ", Light_num, ", Diuron = ", Diuron_num))}) %>% 
  filter(!(Treatment %in% c("Light_only","Diuron_only",  "Control"))) %>% 
  select(-block, -sample_id)

# Take only what we need from raw data
Diuron_growth_data <- xdiuron_light2 %>% 
  select(hours, hours_fact, Diuron_fact, Light_num, 
         Diuron_num, block, sample_id, t0, celld) %>% 
  within(celld_pred <- exp(predict(m1_diuron_light, newdata = .)))

# Save computations for plotting
save(Diuron_growth_CI_data, Diuron_growth_data, m1_diuron_light, 
     file = "Data/Diuron_growth_model_and_data2.RDA")

