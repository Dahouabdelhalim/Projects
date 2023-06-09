
# ******************************************************************************
# ******************************************************************************
# ======================== general =============================================
# we use a Linear Mixed Effects Model approach with island as random effect to 
# capture differences among islands
# SVL is added as covariate for all traits (excl SVL alone)
# The full model: response variable ~ competition * sex + svl
# All reduced models are considered including intercept only
# models are ranked using AIC scores
#
# ======================== notes =============================================
# â€¢ missing data are removed for each variable at a time
#
# ======================== A load data ==============================================
# load correct file corrected and outliers removed
name.t <- "all.traits.clean.size.corrected.outlier.rm.csv"

# define the name of the output file (depends on which data is analysed)
output.file.name <- "2c.all.trait.model.output.covSVL"  

# load data files from data.clean
d <- read.csv(paste(cd.path,name.t, sep = ""), stringsAsFactors = FALSE)
str(d)


# ======================== B mixed-effects SVL ======================
# which response variables to be used (only the raw data)
var.loc.cov <- c(5:15)
colnames(d)[var.loc.cov]

# list to store all the model outputs for each variable in
list.AICc <- vector(mode = "list", length(var.loc.cov))


# *** mixed effects model for only SVL *****************************************
# to select SVL (first variable in var.loc.cov)
i <- 1
colnames(d)[var.loc.cov[i]] # svl

# names of the different models. Must match the model sequence in loop below
model.names <- c("sex * presence", "sex + presence", "sex","presence",
                 "intercept")
# from full model including interaction to intercept only

# make selection of data frame with response var, sex and presence and omit 
# make selection
d.svl <- na.omit(d[, c( var.loc.cov[i], location.expl.var)]) 
str(d.svl)
# check variable name
(var.sel.name <- colnames(d.svl)[1])
# response variable is first column 
var.sel <- d.svl[,1]

# make data storage object
data.t <- as.data.frame(matrix(NA, ncol = 10, nrow = 5))
colnames(data.t) <- c("model", "df", "AICc", "intercept", "sexm", "presence",
                      "SVL","sexm*presence", "R2marginal", "R2conditional")
data.t[,1] <- model.names

# go through all models 
m.1 <- lmer(var.sel ~ d.svl$sex * d.svl$presence + (1|d.svl$island))
m.2 <- lmer(var.sel ~ d.svl$sex + d.svl$presence + (1|d.svl$island))
m.3 <- lmer(var.sel ~ d.svl$sex  + (1|d.svl$island))
m.4 <- lmer(var.sel ~ d.svl$presence + (1|d.svl$island))
m.5 <- lmer(var.sel ~ 1 + (1|d.svl$island))


# assign model output to data storage object
data.t[,2:3] <- AICc(m.1, m.2, m.3, m.4, m.5)
# per model add fixes effects
#  model m.1: get fixed effect measures and R^2 (from MuMIn package)
data.t[1,4:6] <- fixef(m.1)[1:3]
data.t[1,8] <- fixef(m.1)[4]    # interaction
data.t[1,9:10] <- r.squaredGLMM(m.1)

# model m.2
data.t[2,4:6] <- fixef(m.2)
data.t[2,9:10] <- r.squaredGLMM(m.2)

# model m.3
data.t[3,4:5] <- fixef(m.3)
data.t[3,9:10] <- r.squaredGLMM(m.3)

# model m.4
data.t[4,4] <- fixef(m.4)[1]
data.t[4,6] <- fixef(m.4)[2]
data.t[4,9:10] <- r.squaredGLMM(m.4)
# model m.5 
data.t[5,4] <- fixef(m.5)
data.t[5,9:10] <- r.squaredGLMM(m.5)

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

# ******** rest of variables with SVL as covariate *****************************
# change the model names to reflect SVL is covariate
model.names.cov <- c("sex * presence + SVL", "sex + presence + SVL", 
                 "sex + SVL","presence + SVL", "intercept + SVL")

# loop through mixed-effect model for all rest variables (starting at 
# location 2)
for(i in 2:length(var.loc.cov)){

  # make selection of data frame with response var, sex and presence and omit 
  # remove the missing values
  d.cov <- na.omit(d[, c( var.loc.cov[i], c(location.expl.var, 5))]) 
  str(d.cov)
  
  # name of response variable
  (var.sel.name <- colnames(d.cov)[1])
  # response variable is first column 
  var.sel <- d.cov[,1]
  
  # make data storage object
  data.t <- as.data.frame(matrix(NA, ncol = 10, nrow = 5))
  colnames(data.t) <- c("model", "df", "AICc", "intercept", "sexm", "presence",
                        "SVL","sexm*presence", "R2marginal", "R2conditional")
  data.t[,1] <- model.names.cov
  length(var.sel)
  length(d.cov$svl)

  # go through all models 
  m.1 <- lmer(var.sel ~ d.cov$sex * d.cov$presence + d.cov$svl + 
                (1|d.cov$island))
  m.2 <- lmer(var.sel ~ d.cov$sex + d.cov$presence + d.cov$svl + 
                (1|d.cov$island))
  m.3 <- lmer(var.sel ~ d.cov$sex  + d.cov$svl + (1|d.cov$island))
  m.4 <- lmer(var.sel ~ d.cov$presence + d.cov$svl + (1|d.cov$island))
  m.5 <- lmer(var.sel ~ 1 + d.cov$svl + (1|d.cov$island))
  
  # assign model output to data storage object
  data.t[,2:3] <- AICc(m.1, m.2, m.3, m.4, m.5)
  # per model add fixes effects
  #  model m.1: get fixed effect measures and R^2 (from MuMIn package)
  data.t[1,4:8] <- fixef(m.1)[1:5]
  data.t[1,9:10] <- r.squaredGLMM(m.1)
  
  # model m.2
  data.t[2,4:7] <- fixef(m.2)
  data.t[2,9:10] <- r.squaredGLMM(m.2)
  
  # model m.3
  data.t[3,4:5] <- fixef(m.3)[1:2]
  data.t[3,7] <- fixef(m.3)[3]
  data.t[3,9:10] <- r.squaredGLMM(m.3)
  
  # model m.4
  data.t[4,4] <- fixef(m.4)[1]
  data.t[4,6:7] <- fixef(m.4)[2:3]
  data.t[4,9:10] <- r.squaredGLMM(m.4)
  
  # model m.5 
  data.t[5,4] <- fixef(m.5)[1]
  data.t[5,7] <- fixef(m.5)[2]
  data.t[5,9:10] <- r.squaredGLMM(m.5)
  
  # add variable names ( for each iteration of the loop new)
  m.results.all <- cbind.data.frame(rep(var.sel.name, 5), data.t)
  m.results.all.s <- m.results.all[order(m.results.all$AICc),]
  # calculate difference in AICc
  delta.AICc <- c(0,rep(NA, length(model.names.cov)-1))
  for(j in 2:length(model.names.cov)) {
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
write.csv(all.AICc, paste(re.path, output.file.name,".csv", sep = ""), row.names = FALSE)

# ******************************************************************************
# ******************************************************************************










