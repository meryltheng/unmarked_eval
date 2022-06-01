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
# e.g., `Rscript --vanilla ./R/2_generateDetections_HPC.R small 1 100 sol`
# e.g., `Rscript --vanilla ./R/2_generateDetections_HPC.R small 1 100 grp 25`
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
# detDat <- DetectionGenerator(movedf=as.data.frame(moveSims), traps = traps, torus = T, timer = F)
# to use DetectionGeneratorLite, we need to det the biggest steplength to input as bufferDist
moveSims <- moveSims %>% 
  group_by(ID) %>%
  mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
  as.data.frame()
bufferDist = max(moveSims$step[moveSims$step<500],na.rm=T) # to ensure the biggest step is captured
detDat <- DetectionGeneratorLite(movedf = moveSims, traps = traps, bufferDist = bufferDist, torus = T)

if (behav == 'sol'){
  write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_sol_', iter, '_J',nTraps,'.csv', sep=''))
}

if (behav == 'grp'){
  write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_grp_', iter, '_gSize', groupSize, '_J',nTraps,'.csv', sep=''))
}
#write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_', iter, '_J100.csv', sep=''))

"-------- Run detection simulations --------"
# 2 converts DetectionGenerator output to 'photos'/discrete points and apply imperfectDet decay function
x1 <- ImperfectDetConverter(movedf = moveSims, detdf = as.data.frame(detDat), traps=traps, timestep = 5, timeInt = 1, # timestep of simulated movement (mins), timeInt to sample movement trajectories (secs)
                            perfect_dist = 0.3, max_dist = 1)
rm(moveSims)

# 3 corrects for multiple animals (i.e., if any animal was detected at time_seconds, all animals within same Site and time_seconds will be automatically detected)
x1 <- x1 %>% 
  mutate(time_seconds = round(time_seconds)) %>% # round to nearest second; discretise for approximating if multiple animals were detected at the same instance
  group_by(Site, time, time_seconds) %>% 
  mutate(isDetected = if_else(any(isDetected==1), 1, 0, missing = NULL) )%>% 
  as.data.frame()

# 3 calculate time diff between each detection
x2 <- x1 %>% 
  group_by(track) %>% # trackIDs consecutive for the indiv
  filter(row_number()==1) %>%
  group_by(ID,Site) %>%
  arrange(time)  %>% # sort encounters chronologically
  mutate( difftime = c(NA,diff(time))) %>%
  as.data.frame()

# 4 collapse consecutive tracks (within same Site and animal ID) by first making a common detection identifier
x2$det_id <- rep(NA, nrow(x2))
x2$det_id[1] <- 1
k = 2
for (i in 2:nrow(x2)){
  isConsecutive1 <- x2$difftime[i] == 1 & x2$ID[i-1] == x2$ID[i] & 
    x2$Site[i-1] == x2$Site[i]
  if( isTRUE(isConsecutive1) ){
    x2$det_id[i] <- x2$det_id[i-1]
  } else {
    isConsecutive2 <- x2$difftime[i] == 2 & x2$difftime[i-1] == 1 &
      x2$ID[i-1] == x2$ID[i] & x2$Site[i-1] == x2$Site[i]
    if( isTRUE(isConsecutive2)){
      x2$det_id[i] <- x2$det_id[i-1]
    } else{ 
      isConsecutive3 <- x2$difftime[i] == 3 & x2$difftime[i-1] == 2 & 
        x2$difftime[i-2] == 1 & x2$ID[i-1] == x2$ID[i] & x2$Site[i-1] == x2$Site[i]
      if( isTRUE(isConsecutive3) ){
        x2$det_id[i] <- x2$det_id[i-1]
      } else{
        x2$det_id[i] <- k
        k = k+1
      }
    }
  }
}

# 5 attach detection identifier to the original dataset
x3 <- merge(x2[,c('ID','time',"Site","track","det_id")], x1, by.x = c('ID','time',"Site","track"), by.y = c('ID','time',"Site","track"))
x3 <- x3 %>% group_by(det_id) %>% arrange(time_seconds) %>% as.data.frame()
x4 <- x3 %>% filter(isDetected == 1) %>% dplyr::select(-isDetected)

if (behav == 'sol'){
  write.csv(x4,row.names=FALSE,file=paste('Data/Detections/ctDat_', scenario, '_sol_', iter, '_J', nTraps,'.csv', sep=''))
}
if (behav == 'grp'){
  write.csv(x4,row.names=FALSE,file=paste('Data/Detections/ctDat_', scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,'.csv', sep=''))
}

"-------- Save and update progress bar --------"
# count files for each scenario for solitary
sm_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_small_sol_.*_J100 | wc -l", intern = T))
md_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_medium_sol_.*_J100 | wc -l", intern = T))
lg_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_large_sol_.*_J100 | wc -l", intern = T))
sm_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_small_sol_.*_J25 | wc -l", intern = T))
md_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_medium_sol_.*_J25 | wc -l", intern = T))
lg_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_large_sol_.*_J25 | wc -l", intern = T))

# count files for each scenario for group
sm20_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_small_grp_.*_gSize20_J100 | wc -l", intern = T))
sm100_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_small_grp_.*_gSize100_J100 | wc -l", intern = T))
med5_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_medium_grp_.*_gSize5_J100 | wc -l", intern = T))
med25_J100 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_medium_grp_.*_gSize25_J100 | wc -l", intern = T))
sm20_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_small_grp_.*_gSize20_J25 | wc -l", intern = T))
sm100_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_small_grp_.*_gSize100_J25 | wc -l", intern = T))
med5_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_medium_grp_.*_gSize5_J25 | wc -l", intern = T))
med25_J25 <- as.numeric(system("ls -1 ./Data/Detections/| grep ^ctDat_medium_grp_.*_gSize25_J25 | wc -l", intern = T))

nsims <- 100
data <- rbind(c(sm_J100,md_J100,lg_J100, sm_J25,md_J25,lg_J25, # sol
                sm20_J100,sm100_J100,med5_J100,med25_J100, # group J=100
                sm20_J25,sm100_J25,med5_J25,med25_J25),  # group J=25
              nsims-c(sm_J100,md_J100,lg_J100, sm_J25,md_J25,lg_J25, # sol
                      sm20_J100,sm100_J100,med5_J100,med25_J100, # group J=100
                      sm20_J25,sm100_J25,med5_J25,med25_J25))
colnames(data) <- c("sm_J100","md_J100","lg_J100"," sm_J25","md_J25","lg_J25", # sol
                    "sm20_J100","sm100_J100","med5_J100","med25_J100", # group J=100
                    "sm20_J25","sm100_J25","med5_J25","med25_J25")  # group J=25

png("Data/Detections/DetProgress.png", type="cairo")
par(mai=c(1,2,1,1))
plotx<- barplot(data, horiz = TRUE, names.arg = colnames(data) , main = "DetSims", las = 1)
text(x = c(sm_J100,md_J100,lg_J100, sm_J25,md_J25,lg_J25,
           sm20_J100,sm100_J100,med5_J100,med25_J100,
           sm20_J25,sm100_J25,med5_J25,med25_J25)+5, y = plotx, label = data[1,], pos = 1, cex = 0.8, col = "red")
dev.off()
