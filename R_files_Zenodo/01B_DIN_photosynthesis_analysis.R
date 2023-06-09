
library(tidyverse)
library(mgcv)
library(visreg)

source("Scripts/bayes_CI_functions.R")

# Process the data
psii_dat <- read_csv("Data/DIN_light_PSII_long_new.csv") %>% 
  rename(hours = `Time (hours)`, DIN = `NH4Cl (mg/L)`,
         Light_num = `Light (u mol photons m/s/)`, block = Block) %>% 
  group_by(hours, Light_num, DIN, block) %>% 
  summarise(photo_inhib = mean(`Photosynthetic inhibition (%)`), 
            Yield = mean(`Y(II)`), 
            Yield_C = mean(`Average Y (controls)`, na.rm = TRUE)) %>% 
  within({
    
    hours[hours == 0.3] <- 0.33
    
    DIN_num <-  case_when(DIN == "Control" ~ 0,
                          TRUE ~ as.numeric(DIN)
    )
    
    block <-  factor(as.character(block))
    
    hours_fact <- factor(hours)
    
    t0 <- ifelse(hours == 0, Yield, NA)
    inhibt0 <- ifelse(hours == 0, photo_inhib, NA)
    
  }) %>% group_by(Light_num, DIN, block ) %>% 
  fill(t0, inhibt0) %>% ungroup() %>% 
  mutate(sample_id = factor(gsub("."," ", x = paste(DIN, Light_num, block), fixed = TRUE))) %>% 
  filter(hours != 0) 


# Exploratory data analysis -----------------------------------------------

# Yield vs hours by din concentration
ggplot(psii_dat) + 
  aes(x = hours, y = Yield, color = factor(DIN_num)) + 
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

m1_DIN_yield <- gam(Yield ~ hours_fact + 
                      offset(t0)+
                      s(DIN_num, Light_num, k = 18, by = hours_fact) +  
                      s(block, bs = "re"),
                     data = psii_dat)

# Check assumptions
png(filename = "Outputs/fluoro_DIN_GAM_residuals.png", 
    width = 20, height = 20, units = "cm", res = 300)
par(mfrow = c(2,2))
gam.check(m1_DIN_yield)
dev.off()

# Check model
summary(m1_DIN_yield)


# Compute the IR data using empirical bayes method ------------------------

# Generate some data for prediction
xDIN_newdata <- with(psii_dat,{
  expand.grid(list(Light_num = unique(Light_num), DIN_num = c(0, 2.67, 26.7), 
                   hours = unique(hours)))}) %>%
  mutate(block = 2, sample_id = psii_dat$sample_id[[1]], 
         hours_fact = factor(hours), t0 = 0.6) 

# Make predictions on new data 
set.seed(1000)
DIN_bayes <- compute_CI_photo(model = m1_DIN_yield, dat = xDIN_newdata, nsims = 1000) %>% 
  within({
  Treatment <- case_when(Light_num == 80 & DIN_num == 0 ~ "Control",
                         Light_num == 80  ~ "DIN_only",
                         DIN_num == 0  ~ "Light_only",
                         TRUE ~ "Light_DIN")
}) 

# Generate the control matrix
DIN_bayes_control <- DIN_bayes %>% filter(Treatment == "Control")
DIN_control_mat <- DIN_bayes_control[ match(DIN_bayes$hours, DIN_bayes_control$hours) ,]

# Generate light only matrix
DIN_bayes_light <- DIN_bayes %>% filter(DIN_num == 0)
DIN_light_mat <- DIN_bayes_light[ 
  match(paste(DIN_bayes$hours, DIN_bayes$Light_num), 
        paste(DIN_bayes_light$hours, DIN_bayes_light$Light_num)) ,]

# Generate DIN only matrix
DIN_bayes_DIN <- DIN_bayes %>% filter(Light_num == 80 )
DIN_mat <- DIN_bayes_DIN[ 
  match(paste(DIN_bayes$hours, DIN_bayes$DIN_num), 
        paste(DIN_bayes_DIN$hours, DIN_bayes_DIN$DIN_num)) ,]

# Compute the IR matrix
IR_mat2 <- log((select(DIN_bayes, `1`:`1000`)/select(DIN_control_mat, `1`:`1000`))/
                 ((select(DIN_light_mat, `1`:`1000`)/select(DIN_control_mat, `1`:`1000`))*
                    (select(DIN_mat, `1`:`1000`)/select(DIN_control_mat, `1`:`1000`))))

# Check for NAs
rowSums(is.na(IR_mat2))

# Compute the credible intervals
DIN_photo_CI_data <- cbind(data.frame(t(apply(IR_mat2, 1, quantile, 
                                              probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
                                        )), xDIN_newdata) %>% within({
  Treatment <- case_when(Light_num == 80 & DIN_num == 0 ~ "Control",
                         Light_num == 80  ~ "DIN_only",
                         DIN_num == 0  ~ "Light_only",
                         TRUE ~ paste0("Light = ", Light_num, ", DIN = ", DIN_num))}) %>% 
  filter(!(Treatment %in% c("Light_only","DIN_only",  "Control"))) %>%
  select(-block, -sample_id)

# Take only what we need from raw data
DIN_photo_data <- psii_dat %>% 
  select(hours, hours_fact, Light_num, DIN_num, block, sample_id, t0, Yield) %>% 
  within(Yield_pred <- predict(m1_DIN_yield, newdata = .))

# Save computations for plotting
save(DIN_photo_CI_data, DIN_photo_data, m1_DIN_yield, 
     file = "Data/DIN_photo_model_and_data2.RDA")
