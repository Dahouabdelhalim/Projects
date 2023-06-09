# Source code ----
library(modedjagsUI)

# Command line arguments ----
args <- commandArgs(trailingOnly=TRUE)
scale <-  args[1]
model_dir <- args[2]
dir <- file.path("model_output", model_dir, paste0("model_", scale))

# Read in data
md <- readRDS(file.path(dir, paste0("model_", scale, "_jags_inputs.RDS")))

# Send model to jags
start.time=Sys.time()
out=jags.basic(data=md$data,inits=md$inits,parameters.to.save=md$params,
               model.file=file.path("model_output", model_dir, "jags1.txt"),n.chains=15,n.iter=65000,
               n.burnin=40000,n.thin=50,parallel=TRUE,n.cores=15,DIC=FALSE, model_dir = dir)

saveRDS(out,file.path(dir, paste0("model_", scale, "_jags_out.RDS")))
end.time=Sys.time()
runtime=list(start.time,end.time)
saveRDS(runtime,file.path(dir, paste0("model_", scale, "_runtime.RDS")))