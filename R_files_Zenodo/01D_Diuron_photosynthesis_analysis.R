
library(tidyverse)
library(mgcv)
library(visreg)

source("Scripts/bayes_CI_functions.R")

# Process the data
psii_dat <- read_csv("Data/diuron_light_PSII_long_NEW.csv") %>% 
  rename(hours = `Time (hours)`, Diuron = `Diuron (ug/L)`,
         Light_num = `Light (u mol photons m/s/)`, block = Block) %>% 
  group_by(hours, Light_num, Diuron, block) %>% 
  summarise(photo_inhib = mean(`Photosynthetic inhibition (%)`), 
            Yield = mean(`Y(II)`), 
            Yield_C = mean(`Average Y (controls)`, na.rm = TRUE)) %>% 
  within({
    
    hours[hours == 0.3] <- 0.33
    
    Diuron_num <-  case_when(Diuron == "Control" ~ 0,
                             Diuron == "MeOH" ~ 0,
                             TRUE ~ as.numeric(Diuron)
    )
    
    block <-  factor(as.character(block))
    
    hours_fact <- factor(hours)
    
    t0 <- ifelse(hours == 0, Yield, NA)
    inhibt0 <- ifelse(hours == 0, photo_inhib, NA)
    
  }) %>% group_by(Light_num, Diuron, block ) %>% 
  fill(t0, inhibt0) %>% ungroup() %>% 
  mutate(sample_id = factor(gsub("."," ", x = paste(Diuron, Light_num, block), fixed = TRUE))) %>% 
  filter(Diuron != "MeOH", hours != 0) 


# Exploratory data analysis -----------------------------------------------

# Yield vs hours by diuron concentration
ggplot(psii_dat) + 
  aes(x = hours, y = Yield, color = factor(Diuron_num)) + 
  geom_point() + 
  facet_grid(.~ Light_num, scales = "free") + 
  stat_smooth(se = FALSE) + 
  theme_classic()

# Yield vs hours by experimental block
ggplot(psii_dat) + 
  aes(x = hours, y = Yield, color = factor(block)) + 
  geom_point() + 
  facet_grid(.~ Light_num, scales = "free") + 
  stat_smooth(se = FALSE) + 
  theme_classic()


# Fit Generalised Additive Model ------------------------------------------

m1_diuron_yield <- gam(Yield ~ hours_fact + 
                         offset(t0)+
                         s(Diuron_num, Light_num, k = 15, by = hours_fact) +  
                         s(block, bs = "re"),
                       data = psii_dat)

# Check assumptions
png(filename = "Outputs/fluoro_Diuron_GAM_residuals.png", 
    width = 20, height = 20, units = "cm", res = 300)
par(mfrow = c(2,2))
gam.check(m1_diuron_yield)
dev.off()

# Check model
summary(m1_diuron_yield)


# Compute the IR data using empirical bayes method ------------------------

# Generate some data for prediction
xdiuron_newdata <- with(psii_dat,{
  expand.grid(list(Light_num = unique(Light_num), Diuron_num = unique(Diuron_num), 
                   hours = unique(hours)))}) %>%
  mutate(block = 2, sample_id = psii_dat$sample_id[[1]], hours_fact = factor(hours), t0 = 0.6) 

# Make predictions based 
set.seed(1000)
diuron_bayes <- compute_CI_photo(model = m1_diuron_yield, dat = xdiuron_newdata, nsims = 1000) %>% 
  within({
  Treatment <- case_when(Light_num == 80 & Diuron_num == 0 ~ "Control",
                         Light_num == 80  ~ "Diuron_only",
                         Diuron_num == 0  ~ "Light_only",
                         TRUE ~ "Light_Diuron")
}) 

# # Generate the control matrix
diuron_bayes_control <- diuron_bayes %>% filter(Treatment == "Control")
diuron_control_mat <- diuron_bayes_control[ match(diuron_bayes$hours, 
                                                  diuron_bayes_control$hours) ,]

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

# Check for NAs
rowSums(is.na(IR_mat2))

# Compute the credible intervals
Diuron_photo_CI_data <- cbind(data.frame(t(apply(IR_mat2, 1, quantile, 
                                                 probs = c(0.025, 0.5, 0.975), na.rm = TRUE
                                                 ))), xdiuron_newdata) %>% within({
  Treatment <- case_when(Light_num == 80 & Diuron_num == 0 ~ "Control",
                         Light_num == 80  ~ "Diuron_only",
                         Diuron_num == 0  ~ "Light_only",
                         TRUE ~ paste0("Light = ", Light_num, ", Diuron = ", Diuron_num))}) %>% 
  filter(!(Treatment %in% c("Light_only","Diuron_only",  "Control"))) %>% 
  select(-block, -sample_id)

# Take only what we need from raw data
Diuron_photo_data <- psii_dat %>% 
  select(hours, hours_fact, Light_num, Diuron_num, block, sample_id, t0, Yield) %>% 
  within(Yield_pred <- predict(m1_diuron_yield, newdata = .))

# Save computations for plotting
save(Diuron_photo_CI_data, Diuron_photo_data, m1_diuron_yield, 
     file = "Data/Diuron_photo_model_and_data2.RDA")
