# ====================================
# Unmarked abundance evalaution
# Generate detection data (HPC)
# ====================================

library(dplyr)
library(sp)
library(rgeos)
source("0_functions.R")

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
# # 1 extract track overlaps with detectors
# # detDat <- DetectionGenerator(movedf=as.data.frame(moveSims), traps = traps, torus = T, timer = F)
# # to use DetectionGeneratorLite, we need to det the biggest steplength to input as bufferDist
# moveSims <- moveSims %>%
#   group_by(ID) %>%
#   mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
#   as.data.frame()
# bufferDist = max(moveSims$step[moveSims$step<5000],na.rm=T) # to ensure the biggest step is captured
# detDat <- DetectionGeneratorLite(movedf = moveSims, traps = traps, bufferDist = bufferDist, torus = T)
# 
# if (behav == 'sol'){
#   write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_sol_', iter, '_J',nTraps,'.csv', sep=''))
# }
# 
# if (behav == 'grp'){
#   write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_grp_', iter, '_gSize', groupSize, '_J',nTraps,'.csv', sep=''))
# }

if (behav == 'sol'){
  detDat <- read.csv(file=paste('Data/Detections/detDat_', scenario, '_sol_', iter, '_J',nTraps,'.csv', sep=''))
}

if (behav == 'grp'){
  detDat <- read.csv(file=paste('Data/Detections/detDat_', scenario, '_grp_', iter, '_gSize', groupSize, '_J',nTraps,'.csv', sep=''))
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

"-------- Save and update progress bar --------"
# # count files for each scenario for solitary
# hi_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_hi_sol_.*_J100 | wc -l", intern = T))
# md_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_med_sol_.*_J100 | wc -l", intern = T))
# lo_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_lo_sol_.*_J100 | wc -l", intern = T))
# hi_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_hi_sol_.*_J25 | wc -l", intern = T))
# md_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_med_sol_.*_J25 | wc -l", intern = T))
# lo_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_lo_sol_.*_J25 | wc -l", intern = T))
# 
# # count files for each scenario for group
# hi20_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_hi_grp_.*_gSize20_J100 | wc -l", intern = T))
# hi100_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_hi_grp_.*_gSize100_J100 | wc -l", intern = T))
# md5_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_med_grp_.*_gSize5_J100 | wc -l", intern = T))
# md25_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_med_grp_.*_gSize25_J100 | wc -l", intern = T))
# hi20_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_hi_grp_.*_gSize20_J25 | wc -l", intern = T))
# hi100_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_hi_grp_.*_gSize100_J25 | wc -l", intern = T))
# md5_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_med_grp_.*_gSize5_J25 | wc -l", intern = T))
# md25_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_med_grp_.*_gSize25_J25 | wc -l", intern = T))
# 
# nsims <- 100
# data <- rbind(c(hi_J100,md_J100,lo_J100, hi_J25,md_J25,lo_J25, # sol
#                 hi20_J100,hi100_J100,md5_J100,md25_J100, # group J=100
#                 hi20_J25,hi100_J25,md5_J25,md25_J25),  # group J=25
#               nsims-c(hi_J100,md_J100,lo_J100, hi_J25,md_J25,lo_J25, # sol
#                       hi20_J100,hi100_J100,md5_J100,md25_J100, # group J=100
#                       hi20_J25,hi100_J25,md5_J25,md25_J25))
# colnames(data) <- c("hi_J100","md_J100","lo_J100"," hi_J25","md_J25","lo_J25", # sol
#                     "hi20_J100","hi100_J100","md5_J100","md25_J100", # group J=100
#                     "hi20_J25","hi100_J25","md5_J25","md25_J25")  # group J=25
# 
# png("Data/Detections/DetProgress.png", type="cairo")
# par(mai=c(1,2,1,1))
# plotx<- barplot(data, horiz = TRUE, names.arg = colnames(data) , main = "DetSims", las = 1)
# text(x = c(hi_J100,md_J100,lo_J100, hi_J25,md_J25,lo_J25,
#            hi20_J100,hi100_J100,md5_J100,md25_J100,
#            hi20_J25,hi100_J25,md5_J25,md25_J25)+5, y = plotx, label = data[1,], pos = 1, cex = 0.8, col = "red")
# dev.off()
