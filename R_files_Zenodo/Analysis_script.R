# This script has been developed to give an overall idea of the process by which the analysis outlined in Maltby et al 2020 took place. The script falls into three main sections:
# Section A) an example of trialling a GAM on a single species and climate projection, and extracting different model statistics; 
# Section B) using least square means to identify the model to use over all species-climate projection combinations, using the data within the Supplementary Table S8 outputs (which were all derived using a process similar to that outlined in section A of this script); 
# Section C) using the 'final' GAM to generate future projections, and tidying up the resulting projections.


############################################################

# Set up the workspace

############################################################

# Set the working directory for where the files are found 
setwd("")

# Load relevant packages (to install use install.packages() first)
library(mgcv)
library(reshape2)
library(emmeans)
library(plyr)
library(MuMIn)
library(nlme)

# Load in the data for the model to run off initially, which for this example contains CPUE and environmental data for the 2000s decade for the SRES Ens_00 climate projection.

training<-read.csv("Ens_00_CPUE_Env_2000s.csv")
# Remove NAs
training <- training[!is.na(training$Plaice_CPUE),]


#########################################################################################################

                                # Section A #

# Trialling GAMs over species and climate projections. 

# This could be looped so that the same GAM runs over every species-climate projection combination, but for the purposes of this script, a single species and climate projection example is used.

#########################################################################################################


# Run a GAM over one species - this example is the 'Full Model' - Model A, for European plaice and Ens_00 climate projection.

Ens_00_mod_plaice <- gam(Plaice_CPUE ~ 
             s(Depth, k=5) + 
             s(Av_Effort, k=5) + 
             s(Grain_Size, k=5) +  
             s(Dec_Av_NBS, k=5) +   
             s(Dec_Av_SST, k=5) + 
             s(Dec_Av_NBT, k=5) +
             s(Summer_Av_SST, k=5) + 
             s(Summer_Av_NBT, k=5) + 
             s(Winter_Av_SST, k=5) + 
             s(Winter_Av_NBT, k=5), 
             family=gaussian (link="identity"), data=training) 

# Look at some of the model statistics.

# Summary of the model
Ens_00_mod_plaice_summary <- summary.gam(Ens_00_mod_plaice) # Need to store this as an object to allow extraction of deviance explained, GCV and R Squared scores.

# Model diagnostic plots
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2 so it plots all 4 plots together
gam.check(Ens_00_mod_plaice) # Diagnostics and plots provided

# Look at and extract the model test statistics:

Ens_00_mod_plaice_AICc <- AICc(Ens_00_mod_plaice) # AICc score
Ens_00_mod_plaice_rsq <- Ens_00_mod_plaice_summary$r.sq # Adjusted R Squared
Ens_00_mod_plaice_dev_expl <- Ens_00_mod_plaice_summary$dev.expl # Deviance Explained
Ens_00_mod_plaice_GCV <- unname(Ens_00_mod_plaice_summary$sp.criterion) # GCV

# Compare modelled projections for 2000s with observations for 2000s through a correlation test
training_projections<-predict(Ens_00_mod_plaice, type = "response")
Ens_00_mod_plaice_cor_test<-cor.test(training$Plaice_CPUE,training_projections)
Ens_00_mod_plaice_correlation <- unname(Ens_00_mod_plaice_cor_test$estimate) # Correlation score

# Store these test statistics in a dataframe. We suggest then saving these in a csv/txt file to allow easier comparison across trialled GAM models, species and climate projections. 
Ens_00_mod_plaice_test_statistics <- c(Ens_00_mod_plaice_AICc, Ens_00_mod_plaice_rsq, Ens_00_mod_plaice_dev_expl, Ens_00_mod_plaice_GCV, Ens_00_mod_plaice_correlation) # Make a vector with scores
Test_statistic <- c("AICc", "Adjusted_R_Squared", "Deviance_Explained", "GCV", "Correlation") # Make vector with score names
Ens_00_mod_plaice_test_statistics<-data.frame(Test_statistic, Ens_00_mod_plaice_test_statistics) # Bring together

# Some 'cleaning' to make the dataframe more interpretable
colnames(Ens_00_mod_plaice_test_statistics)[2] <- "Score"
Ens_00_mod_plaice_test_statistics$Model <- rep("Model A",5)
Ens_00_mod_plaice_test_statistics$Species <- rep("European_plaice", 5)
Ens_00_mod_plaice_test_statistics$Climate_Projection <- rep("Ens_00", 5)

# Generating Delta AICc and Akaike weights can be done after all models have been run for each species- climate projection combination.

#########################################################################################################
                      # Section B #

# Using least square means to identify an optimal GAM to use across all species-climate projection combinations.
# The example uses the ready-made file which contains GAM test statstic outputs from all trialled GAMs used in the paper, Supplementary Table S8. We use correlation here as an example.
# We use the emmeans package to undertake the analysis (NB. emmeans and lsmeans are the same, just different terminologies used)

#########################################################################################################

model_outputs <- read.csv("Supplementary_Table_S8.csv")

# Make into a long format for ease of analysis
model_outputs <- melt(model_outputs, id=c("Species", "Climate_Projection", "Test_Statistic"))
colnames(model_outputs) <- c("Species", "Climate_Projection", "Test_Statistic", "Model", "Score")

# An example for correlation here, but the same process can be used for all model statistics
correlations <- subset(model_outputs, Test_Statistic == c("Correlation"))
colnames(correlations)[5] <- "Correlation"

# Develop a simple linear model for the least square means analysis
correlation_model <- aov(Correlation~ Species + Model + Climate_Projection, data = correlations) 
correlation_model # Have a look at it
summary(correlation_model) # Have a look at it

# Run the emmeans function to generate the least square means
correlation_lsmeans <- data.frame(emmeans(correlation_model, "Model")) 
# "Model" is listed here to generate the least square means for each model, by controlling for the factors species and climate projection.
correlation_lsmeans # Have a look at the outputs. 

# This process can be repeated for all test statistics, and the outputs merged into a dataframe and then saved into a csv/txt file (similar to that shown in Table S9 in the supplementary of the paper)


#########################################################################################################

                          # Section C # 

# Using the GAM to generate future projections. 

#########################################################################################################

# Read in the future data to make the projections - this example uses Ens_00 (to follow on from Section A)

Ens_00_2000 <- read.csv("Ens_00_Env_2000s.csv")
Ens_00_2010 <- read.csv("Ens_00_Env_2010s.csv")
Ens_00_2020 <- read.csv("Ens_00_Env_2020s.csv")
Ens_00_2030 <- read.csv("Ens_00_Env_2030s.csv")
Ens_00_2040 <- read.csv("Ens_00_Env_2040s.csv")
Ens_00_2050 <- read.csv("Ens_00_Env_2050s.csv")
Ens_00_2060 <- read.csv("Ens_00_Env_2060s.csv")
Ens_00_2070 <- read.csv("Ens_00_Env_2070s.csv")
Ens_00_2080 <- read.csv("Ens_00_Env_2080s.csv")

# Run the projections 
Ens_00_plaice_preds_2000 <- predict(Ens_00_mod_plaice, Ens_00_2000, type="response", se.fit=T)
Ens_00_plaice_preds_2010 <- predict(Ens_00_mod_plaice, Ens_00_2010, type="response", se.fit=T)
Ens_00_plaice_preds_2020 <- predict(Ens_00_mod_plaice, Ens_00_2020, type="response", se.fit=T)
Ens_00_plaice_preds_2030 <- predict(Ens_00_mod_plaice, Ens_00_2030, type="response", se.fit=T)
Ens_00_plaice_preds_2040 <- predict(Ens_00_mod_plaice, Ens_00_2040, type="response", se.fit=T)
Ens_00_plaice_preds_2050 <- predict(Ens_00_mod_plaice, Ens_00_2050, type="response", se.fit=T)
Ens_00_plaice_preds_2060 <- predict(Ens_00_mod_plaice, Ens_00_2060, type="response", se.fit=T)
Ens_00_plaice_preds_2070 <- predict(Ens_00_mod_plaice, Ens_00_2070, type="response", se.fit=T)
Ens_00_plaice_preds_2080 <- predict(Ens_00_mod_plaice, Ens_00_2080, type="response", se.fit=T)

# Bring the projections together into one dataframe, and rename the columns. 
plaice_predictions_Ens_00<-data.frame(Ens_00_plaice_preds_2000, Ens_00_plaice_preds_2010, Ens_00_plaice_preds_2020, Ens_00_plaice_preds_2030, Ens_00_plaice_preds_2040, Ens_00_plaice_preds_2050, Ens_00_plaice_preds_2060, Ens_00_plaice_preds_2070, Ens_00_plaice_preds_2080) 
colnames(plaice_predictions_Ens_00)<-c("2000", "2000_SE", "2010", "2010_SE", "2020", "2020_SE", "2030", "2030_SE", "2040", "2040_SE", "2050", "2050_SE", "2060", "2060_SE", "2070", "2070_SE", "2080", "2080_SE")

# Have a look at the final dataset
head(plaice_predictions_Ens_00)

############################################################

# Tidy up the resulting projections

############################################################

# Make the projections into a more suitable 'long' format with one grid cell column, one climate projection column, one decade column, one projections column and one standard error column (NB there maybe a much more elegant way to do this...)
plaice_predictions_Ens_00$Grid_Cell <- Ens_00_2000$Grid_Cell
plaice_predictions_Ens_00$Climate_Projection <- rep("Ens_00", length(plaice_predictions_Ens_00$`2000`))

# Melt into the long format
plaice_predictions_Ens_00 <- melt(plaice_predictions_Ens_00, id=c("Grid_Cell", "Climate_Projection"))

# Split the 'variable' column into two columns, one for the projections and one for the standard error.
plaice_predictions_Ens_00_decade <- subset(plaice_predictions_Ens_00, variable %in% c("2000", "2010", "2020", "2030", "2040", "2050", "2060", "2070", "2080"))
plaice_predictions_Ens_00_decade[plaice_predictions_Ens_00_decade < 0]<-0

plaice_predictions_Ens_00_se <- subset(plaice_predictions_Ens_00, variable %in% c("2000_SE", "2010_SE", "2020_SE", "2030_SE", "2040_SE", "2050_SE", "2060_SE", "2070_SE", "2080_SE"))

# Bring together and tidy up
plaice_predictions_Ens_00 <- data.frame(plaice_predictions_Ens_00_decade, plaice_predictions_Ens_00_se$value)
head(plaice_predictions_Ens_00) # Check it to make sure correct
colnames(plaice_predictions_Ens_00) <- c("Grid_Cell", "Climate_Projection", "Decade", "Projected_CPUE", "SE")
plaice_predictions_Ens_00$Species <- rep("European_plaice", length(plaice_predictions_Ens_00$Grid_Cell))
plaice_predictions_Ens_00 <- plaice_predictions_Ens_00[,c(6,2,1,3,4,5)]
plaice_predictions_Ens_00 <- plaice_predictions_Ens_00[!is.na(plaice_predictions_Ens_00$Projected_CPUE),]

# Make a column of backtransformed abundances
plaice_predictions_Ens_00$Back_Transformed_Projected_CPUE <- (plaice_predictions_Ens_00$Projected_CPUE^2)^2

head(plaice_predictions_Ens_00) 

############################################################

# End.

# For any further help or information, contact katherine.maltby1@gmail.com

############################################################




