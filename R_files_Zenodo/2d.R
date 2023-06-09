
# ******************************************************************************
# ******************************************************************************
# ======================== general =============================================
# we use a Linear Mixed Effects Model approach with island as random effect to 
# capture differences among islands
# use the SVL corrected for variables
# The full model: response variable ~ competition * sex
# All reduced models are considered including intercept only
# models are ranked using AIC scores
#
# A load data
# B mixed-effects model
#
# ======================== notes =============================================
# â€¢ missing data are removed for each variable at a time
#
# ======================== load data ==============================================
# load correct file corrected and outliers removed
name.t <- "all.traits.clean.size.corrected.outlier.rm.csv"

# define the name of the output file (depends on which data is analysed)
output.file.name <- "2d.size.corrected.all.trait.model.output."  

# load data files from data.clean
d <- read.csv(paste(cd.path,name.t, sep = ""), stringsAsFactors = FALSE)
str(d)

# ================= B mixed-model for all size corrected variables =============
var.loc <- c(5,16:25) # svl and size (svl) corrected variables 
colnames(d)[var.loc]

# names of the different models. Must match the model sequence in loop below
model.names <- c("sex * presence", "sex + presence", "sex","presence",
                 "intercept")
# from full model including interaction to intercept only

# list to store all the model outputs for each variable in
list.AICc <- vector(mode = "list", length(var.loc))


for(i in 1:length(var.loc)){
  # make selection of data frame with response var, sex and presence and omit 
  # remove the missing values
  d.sc <- na.omit(d[, c( var.loc[i], location.expl.var)]) 
  str(d.sc)
  
  # store variable name for use in summary file
  (var.sel.name <- colnames(d.sc)[1])
  
  # response variable is first column 
  var.sel <- d.sc[,1]
  
  # make data storage object
  data.t <- as.data.frame(matrix(NA, ncol = 9, nrow = 5))
  colnames(data.t) <- c("model", "df", "AICc", "intercept", "sexm", "presence",
                        "sexm*presence", "R2marginal", "R2conditional")
  data.t[,1] <- model.names
  
  # go through all models 
  m.1 <- lmer(var.sel ~ d.sc$sex * d.sc$presence + (1|d.sc$island))
  m.2 <- lmer(var.sel ~ d.sc$sex + d.sc$presence + (1|d.sc$island))
  m.3 <- lmer(var.sel ~ d.sc$sex  + (1|d.sc$island))
  m.4 <- lmer(var.sel ~ d.sc$presence + (1|d.sc$island))
  m.5 <- lmer(var.sel ~ 1 + (1|d.sc$island))
  
  
  # assign model output to data storage object
  data.t[,2:3] <- AICc(m.1, m.2, m.3, m.4, m.5)
  # per model add fixes effects
  #  model m.1: get fixed effect measures and R^2 (from MuMIn package)
  data.t[1,4:7] <- fixef(m.1)
  data.t[1,8:9] <- r.squaredGLMM(m.1)
  # model m.2
  data.t[2,4:6] <- fixef(m.2)
  data.t[2,8:9] <- r.squaredGLMM(m.2)
  # model m.3
  data.t[3,4:5] <- fixef(m.3)
  data.t[3,8:9] <- r.squaredGLMM(m.3)
  # model m.4
  data.t[4,4] <- fixef(m.4)[1]
  data.t[4,6] <- fixef(m.4)[2]
  data.t[4,8:9] <- r.squaredGLMM(m.4)
  # model m.5 
  data.t[5,4] <- fixef(m.5)
  data.t[5,8:9] <- r.squaredGLMM(m.5)
  
  # add variable names ( for each iteration of the loop new)
  m.results.all <- cbind.data.frame(rep(var.sel.name, 5), data.t)
  m.results.all.s <- m.results.all[order(m.results.all$AICc),]
  # calculate difference in AICc
  delta.AICc <- c(0,rep(NA, length(model.names)-1))
  for(j in 2:length(model.names)) {
    delta.AICc[j] <- m.results.all.s$AICc[j] - m.results.all.s$AICc[1]
  }
  # add to the data frame
  m.results.all.s2 <- cbind.data.frame(m.results.all.s, delta.AICc)
  # add column names
  colnames(m.results.all.s2)[1] <- "variable.name"
  
  # store in list
  list.AICc[[i]] <- m.results.all.s2
}

# merge all AIC scores
all.AICc <- do.call("rbind", list.AICc)

# save file as .csv to working directory
write.csv(all.AICc, paste(re.path, output.file.name,".csv", sep = ""), 
          row.names = FALSE)


# ******************************************************************************
# ******************************************************************************










