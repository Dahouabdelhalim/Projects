# Manuscript: "Anthropogenic noise impairs cooperation in bottlenose dolphins"
# by Sørensen, Haddock, Guarino, Jaakkola, McMullen, Jensen, Tyack and King
#
# This scripts models task success rate as a function of trial noise and plots data
# Last updated: 02-11-2022
################################################################################
# Clear environment
rm(list = ls())

# Load packages
library(ggplot2)
library(DHARMa) 
library(lme4)
library(lmerTest)

# Load data
my_data  <- read.csv("Data_successrate.csv", header=T)
attach(my_data) 
str(my_data)

############################# Plot success rate ################################
# Find number of trials per treatment
treatment1 <- length(outcome[treatment == "1"]) # 40
treatment2 <- length(outcome[treatment == "2"]) # 40
treatment3 <- length(outcome[treatment == "3"]) # 40
treatment4 <- length(outcome[treatment == "4"]) # 40
treatment5 <- length(outcome[treatment == "5"]) # 40
treatmentss <- c("1","2","3","4","5") 

# Find success rate for each treatment
success1 = (length(which(outcome[treatment == "1"]=="1")))/length(outcome[treatment == "1"]) 
success2 = (length(which(outcome[treatment == "2"]=="1")))/length(outcome[treatment == "2"]) 
success3 = (length(which(outcome[treatment == "3"]=="1")))/length(outcome[treatment == "3"]) 
success4 = (length(which(outcome[treatment == "4"]=="1")))/length(outcome[treatment == "4"]) 
success5 = (length(which(outcome[treatment == "5"]=="1")))/length(outcome[treatment == "5"]) 

# Prepare data for plotting:
y1 <- c(success1*100, success2*100, success3*100, success4*100, success5*100)
y1 <- data.frame(y1)
texttt <- c("34/40","31/40","26/40","27/40","25/40")

# Set colours:
my_colors <- RColorBrewer::brewer.pal(6, "Blues")[2:6] # Choose for 'Blues' color scale

##############################################################################
#### Model fit and model analysis ####
################################################################################
#### Successrate ~ EFD_combined (noise level per trial quantified as energy flux density, dB re 1 uPa^2*s) 
#### with session number and day included as a random effect #####
model_full1 <- glmer(outcome ~ EFD_combined + (1|session) + (1|day), family = binomial, data = my_data) 
model_null1 <- glmer(outcome ~ 1 + (1|session) + (1|day), family = binomial, data = my_data)

# Model results
anova(model_null1,model_full1)
summary(model_full1)
confint(model_full1)

# Check modelfit
simulationOutput<-simulateResiduals(fittedModel = model_full1, plot=T) 

##################### Plot success rate vs noise ##########################
topmod1     <- model_full1
effects_mod <- effects::predictorEffect(predictor = "EFD_combined", mod = topmod1)
x_mod       <- as.data.frame(effects_mod)

#remove q[[1]] and labs (see first plot line) if plotting figure individually
p1 <- ggplot() + geom_point(data = my_data, aes(EFD_combined, outcome, group = as.factor(treatment), color = as.factor(treatment)), size=4, alpha=0.8,  position = position_jitter(w = 0, h = 0.04)) +
  
  scale_colour_manual(values=my_colors) +
  
  geom_line(data = x_mod, aes(x = EFD_combined, y=fit), color="#5B8FA8FF",lwd=2) +
  
  geom_ribbon(data = x_mod, aes(x = EFD_combined, ymin=lower, ymax=upper), alpha= 0.3, fill="#5B8FA8FF") +
  
  scale_y_continuous(name = "trial outcome", breaks = c(0,1), labels = c("failure","success") 

                     , sec.axis = sec_axis(~.*100, name = "successrate (%)", labels = c("0","25","50","75","100"))) +

  labs(x = expression("noise level (EFD, dB re 1 µPa"^2*"·s)")) +

  # stat_cor(method = "pearson", label.x = 0.4, label.y = 0.5, output.type = "text",digits = 2,data = fullA,aes(maxSRI, a1_a),cor.coef.name = "R",size=11)+
  
  theme(plot.margin = unit (c (0.5,0.5,0.5,0.5), 'cm'),legend.position = "none", axis.text = element_text(size = 25)
        , axis.title = element_text(size = 25), panel.grid = element_blank(), panel.background = element_blank()
        , axis.line = element_line(colour = "black"), panel.border = element_rect(fill=NA))

p1


## Without labels
p2 <- ggplot() + geom_point(data = my_data, aes(EFD_combined, outcome, group = as.factor(treatment), color = as.factor(treatment)), size=4, alpha=0.8,  position = position_jitter(w = 0, h = 0.04)) +
  
  scale_colour_manual(values=my_colors) +
  
  geom_line(data = x_mod, aes(x = EFD_combined, y=fit), color="#5B8FA8FF",lwd=2) +
  
  geom_ribbon(data = x_mod, aes(x = EFD_combined, ymin=lower, ymax=upper), alpha= 0.3, fill="#5B8FA8FF") +
  
  scale_y_continuous(name = "", breaks = c(0,1), labels = c("","") 
                     
                     , sec.axis = sec_axis(~.*100, name = "", labels = c("","","","",""))) +
  
  labs(x = expression("")) +
  
  scale_x_continuous(labels = c("","","","","","")) +
  
  # stat_cor(method = "pearson", label.x = 0.4, label.y = 0.5, output.type = "text",digits = 2,data = fullA,aes(maxSRI, a1_a),cor.coef.name = "R",size=11)+
  
  theme(plot.margin = unit (c (0.5,0.5,0.5,0.5), 'cm'),legend.position = "none", axis.text = element_text(size = 25)
        , axis.title = element_text(size = 25), panel.grid = element_blank(), panel.background = element_blank()
        , axis.line = element_line(colour = "black"), panel.border = element_rect(fill=NA))

p2

################################ END OF SCRIPT #################################



