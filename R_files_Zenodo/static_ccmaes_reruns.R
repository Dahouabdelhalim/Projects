#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("At least one argument must be supplied (idx)", call. = FALSE)
} else if (length(args) > 2) {
  stop("At most three arguments must be supplied")
}

idx  <- as.integer(args[1])+1
max_idx <- as.integer(args[2])

library(irace)

#Location to store the logs
log_dir = "Data/logs_temp"

s <- readScenario("scenario_static.txt")
p <- readParameters("ccmaes_paramters_static.txt")

#This should be set to a datatable determining which functions, dimensions and seeds to run (is split automatically to run on cluster)
full_dt <- read.csv("dt_work_division.csv")
set.seed(42) #set seed to ensure the split across nodes is consistent
n_rep <- as.integer(ceiling(nrow(full_dt)/max_idx))
dt <- split(full_dt, sample(rep(1:max_idx, n_rep), nrow(full_dt), replace=F))[[idx]]

for (sub_idx in seq(nrow(dt))) {
    settings <- dt[sub_idx, ]
    fid <- as.integer(settings[['fid']])
    dim <-as.integer(settings[['dim']])
    seed <- as.integer(settings[['seed']])
    scen <- s
    scen$logFile <- paste0(log_dir, "log_static_bounds_F", fid, "_", dim, "D_seed", seed)
    scen$boundMax <- 10000 * dim
    p$domain$fid <- fid
    p$domain$dim <- dim
    scen$seed <- seed
    scen <- checkScenario(scen)
    irace(scen, p)
}

