# ====================================
# Unmarked abundance evalaution
# Generate detection data (HPC)
# ====================================

library(dplyr)
library(sp)
library(rgeos)
source("R/0_functions.R")

'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for sceenario, arg[2] for iter)
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# test if there is at least one argument: if not, return an error
if (length(args) <= 1) {
  stop("At least two arguments must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input: ", args[1], args[2], sep = ' '))
  scenario = args[1]
  iter = as.numeric(args[2])
}

"-------- Prep input parameters --------"
days = 100
nbObs = 4*12*days # 10-min timesteps, 12-h activity, x days

"-------- Run detection simulations --------"
load('Data/unmarkedEnvGrid100.RData')
moveSims <- read.csv(file=paste('Data/MovementSims/simDat_', scenario, '_', iter, '.csv', sep=''))
detDat <- DetectionGenerator(df=as.data.frame(moveSims), X = traps.buff, torus = T)
write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_', iter, '_G100.csv', sep=''))