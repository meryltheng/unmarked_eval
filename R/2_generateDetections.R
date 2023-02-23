# ====================================
# Unmarked abundance evalaution
# Generate detection data
# ====================================

library(dplyr)
library(sp)
library(rgeos)
library(foreach)
source("R/0_functions.R")

# "-------- Prep detection sims input --------" # skip to load if array already made
# ### Make CT array 
# xlim <- ylim <- c(0,10000)
# 
# ## This makes regular grid
# # J = 100
# traps <- makeTraps(xlim = c(500,9500), ylim = c(500,9500), trapspacing = 1000, 
#                         r = 10, theta = 2/9*pi, jitter=F) # 1 unit = 10m (radius)
# save(traps, file = 'Data/traps100.RData')
# 
# # J = 25
# traps <- makeTraps(xlim = c(1000,9000), ylim = c(1000,9000), trapspacing = 2000, 
#                         r = 10, theta = 2/9*pi, jitter=F)
# # save(traps, file = 'Data/traps25.RData')
# 
# plot(traps[[1]],axes=T, xlim = xlim, ylim = ylim)

'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for scenario, arg[2] for iter, arg[3] for nTraps, arg[4] for movemt behav (sol/grp), if grp arg[5] is group size)
# e.g., `Rscript --vanilla ./R/2_generateDetections_HPC.R hi 1 100 sol`
# e.g., `Rscript --vanilla ./R/2_generateDetections_HPC.R hi 1 100 grp 25`
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# test if there is at least one argument: if not, return an error
if (length(args) <= 3) {
  stop("At least four arguments must be supplied (input file).\n", call. = FALSE) # add further checks
} else {
  if (length(args) == 4) {
    if (args[4] == 'grp') stop("if group, supply 5th arg for group size")
    print(paste0("Arg input:", args[1], args[2], args[3], args[4], sep = ' ')) 
    scenario = args[1]
    iter = as.numeric(args[2])
    nTraps = as.numeric(args[3])
    behav = args[4]
  } else {
    print(paste0("Arg input:", args[1], args[2], args[3], args[4], args[5], sep = ' '))
    scenario = args[1]
    iter = as.numeric(args[2])
    nTraps = as.numeric(args[3])
    behav = args[4]
    groupSize = as.numeric(args[5])
  }
}

"-------- Read data --------"
load(paste0('Data/traps',nTraps,'.RData'))

if (behav == 'sol'){
  moveSims <- read.csv(file=paste('Data/MovementSims/simDat_', scenario, '_sol_', iter, '.csv', sep=''))
}

if (behav == 'grp'){
  moveSims <- read.csv(file=paste('Data/MovementSims/simDat_', scenario, '_grp_', iter, '_gSize', groupSize,'.csv', sep=''))
}

"-------- Run detection simulations --------"
# 1 extract track overlaps with detectors
# to use DetectionGeneratorLite, we need to det the biggest steplength to input as bufferDist
moveSims <- moveSims %>%
  group_by(ID) %>%
  mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
  as.data.frame()
bufferDist = max(moveSims$step[moveSims$step<5000],na.rm=T) # to ensure the biggest step is captured
detDat <- DetectionGeneratorLite(movedf = moveSims, traps = traps, bufferDist = bufferDist, torus = T)

if (behav == 'sol'){
  write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_sol_', iter, '_J',nTraps,'.csv', sep=''))
}

if (behav == 'grp'){
  write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_grp_', iter, '_gSize', groupSize, '_J',nTraps,'.csv', sep=''))
}
  
# 2 converts DetectionGenerator output to 'photos'/discrete points
x1 <- DetSamplPoints(movedf = moveSims, detdf = as.data.frame(detDat), traps=traps, timestep = 5, timeInt = 1) # timestep of simulated movement (mins), timeInt to sample movement trajectories (secs)
rm(moveSims)

# apply imperfectDet decay function
x2 <- x1 %>% #i
  group_by(Site, ID) %>% #ii
  arrange(time) %>% #iii
  mutate( difftime = c(NA, diff(time)) ) %>% #iv
  group_by(track) %>% # v
  # if the first row of track was flagged as consecutive track (difftime==1) AND if start of track falls within FOV (entryTime<1s), 
  # assign it to same time as prev consec track
  mutate( time_new = if_else( first(difftime) == 1 & (!is.na(first(difftime))) & first(entryTime)<=1, 
                              first(time)-1, as.numeric(first(time)), NA_real_) ) %>% # vi
  group_by(Site, ID, time_new) %>% # vii
  mutate( isDetected = imperfectDet(d = dist, perfect_dist = 3, max_dist = 10, burst = T, bursts = 7)) %>%
  group_by(Site, ID) %>% 
  arrange(time_seconds) %>% 
  # filter(isDetected == 1) %>% 
  # dplyr::select(-isDetected) %>% 
  as.data.frame()

# 3 corrects for multiple animals (i.e., if any animal was detected at time_seconds, all animals within same Site and time_seconds will be automatically detected)
x3 <- x2 %>% 
  mutate(time_seconds = round(time_seconds)) %>% # round to nearest second; discretise for approximating if multiple animals were detected at the same instance
  group_by(Site, time, time_seconds) %>% 
  mutate(isDetected = if_else(any(isDetected==1), 1, 0, missing = NULL) )%>% 
  as.data.frame()

# filter out all non-detections
x4 <- x3 %>% group_by(Site, ID) %>% 
  arrange(time_seconds) %>% 
  filter(isDetected == 1) %>% 
  dplyr::select(-isDetected) %>% 
  as.data.frame()

if (behav == 'sol'){
  write.csv(x4,row.names=FALSE,file=paste('Data/Detections/ctDat_', scenario, '_sol_', iter, '_J', nTraps,'.csv', sep=''))
}
if (behav == 'grp'){
  write.csv(x4,row.names=FALSE,file=paste('Data/Detections/ctDat_', scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,'.csv', sep=''))
}

# "-------- Check for site heterogeneity --------"
# library(tidyr)
# data <- detDat %>% 
#   count(Site) %>%
#   complete(Site = 1:100, fill = list(n = 0)) %>% as.data.frame()
# data$cov <- raster::extract(env,traps.buff, mean)
# plot(data$cov, data$n, xlab = 'habitat_cov', ylab = "nDet")
# 
# data <- cbind(data,coordinates(traps.buff))
# points(data$`1`,data$`2`,pch=19,cex=(data$n+0.5)/2,col='black')
