# ====================================
# Unmarked abundance evalaution
# Generate detection data
# ====================================

library(dplyr)
library(sp)
library(rgeos)
library(foreach)
source("R/0_functions.R")

# Read in metas
nsim = 10
scenario = 'medium'
days = 100
nbObs = 12*12*days # 5-min timesteps, 12-h activity, 50 days

behav = 'grp'
groupSize = 25 # if behav == 'grp

nTraps = 100 # c('100','25')

"-------- Prep detection sims input --------" # skip to load if array already made
### Make CT array 
## This makes regular grid
# J = 100
traps <- makeTraps(xlim = c(50,950), ylim = c(50,950), trapspacing = 100, 
                        r = 1, theta = pi/3, jitter=F) # 1 unit = 10m (radius); creates equilateral triangular detection zones
# J = 25
traps <- makeTraps(xlim = c(100,900), ylim = c(100,900), trapspacing = 200, 
                        r = 1, theta = pi/3, jitter=F)

plot(traps[[1]],axes=T, xlim = xlim, ylim = ylim)
# save(traps, file = 'Data/traps100.RData')

"-------- Run detection simulations --------"
if (nTraps == 100){
  load('Data/traps100.RData')
}
if (nTraps == 25){
  load('Data/traps25.RData')
}

foreach(iter=1:10, .inorder=FALSE, .combine=rbind, .packages=c("sp","dplyr","rgeos")) %do% {
  if (behav == 'sol'){
    moveSims <- read.csv(file=paste('/Volumes/LaCie/MERYL/PhD/Chapters/Paper 3/UnmarkedAbundance/Data/MovementSims/simDat_', scenario, '_sol_', iter, '.csv', sep=''))
  }
  if (behav == 'grp'){
    moveSims <- read.csv(file=paste('/Volumes/LaCie/MERYL/PhD/Chapters/Paper 3/UnmarkedAbundance/Data/MovementSims/simDat_', scenario, '_grp_', iter, '_gSize', groupSize,'.csv', sep=''))
  }
  #moveSims$time <- rep(1:nbObs, length(unique(moveSims$ID)))
  
  # 1 extract track overlaps with detectors
  # detDat <- DetectionGenerator(movedf=as.data.frame(moveSims), traps = traps, torus=T, timer = T)
  moveSims <- moveSims %>% # to use DetectionGeneratorLite, we need to det the biggest steplength to input as bufferDist
    group_by(ID) %>%
    mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
    as.data.frame()
  bufferDist = max(moveSims$step[moveSims$step<500],na.rm=T) # to ensure the biggest step is captured
  detDat <- DetectionGeneratorLite(movedf = moveSims, traps = traps, bufferDist = bufferDist, torus = T)
  
  if (behav == 'sol'){
    write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_sol_', iter, '_J', nTraps,'.csv', sep=''))
  }
  if (behav == 'grp'){
    write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,'.csv', sep=''))
  }
  
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
}

"-------- Check for site heterogeneity --------"
library(tidyr)
data <- detDat %>% 
  count(Site) %>%
  complete(Site = 1:100, fill = list(n = 0)) %>% as.data.frame()
data$cov <- raster::extract(env,traps.buff, mean)
plot(data$cov, data$n, xlab = 'habitat_cov', ylab = "nDet")

data <- cbind(data,coordinates(traps.buff))
points(data$`1`,data$`2`,pch=19,cex=(data$n+0.5)/2,col='black')
