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
scenario = 'small'
days = 100
nbObs = 6*12*days # 15-min timesteps, 12-h activity, 50 days
#xlim <- ylim <- c(0,1000)

"-------- Prep detection sims input --------"
### Make CT array 
## This makes regular grid
# J = 100
traps.buff <- makeTraps(xlim = c(50,950), ylim = c(50,950), trapspacing = 100, 
                        r = 1, w = 1/2, jitter=T) # 1 unit = 10m (radius); creates equilateral triangular detection zones
# J = 25
traps.buff <- makeTraps(xlim = c(100,900), ylim = c(100,900), trapspacing = 200, 
                        r = 1, w = 1/2, jitter=F)

## This makes grid cluster
traps.buff <- makeTrapsCluster(xlim = c(55,445), ylim = c(55,445), clusterspacing = 130,trapspacing = 15, 
                             n = 9, nClus = 3, r = 1, w = 1/2)
traps.buff <- makeTrapsCluster(xlim = c(50,450), ylim = c(50,450), clusterspacing = 100,trapspacing = 15, 
                               n = 20, nClus = 2, r = 1, w = 1/2)

plot(traps.buff,axes=T, xlim = xlim, ylim = ylim)
# save(traps.buff, file = 'Data/unmarkedEnvGrid100_jit.RData')

"-------- Run detection simulations --------"
load('Data/unmarkedEnvGrid100.RData')
load('Data/unmarkedEnvGrid25.RData')

detOut <- foreach(iter=1:10, .inorder=FALSE, .combine=rbind, .packages=c("sp","dplyr","rgeos")) %do% {
  moveSims <- read.csv(file=paste('Data/MovementSims/simDat_', scenario, '_', iter, '.csv', sep=''))
  moveSims$time <- rep(1:nbObs, length(unique(moveSims$ID)))
  detDat <- DetectionGenerator(df=as.data.frame(moveSims), X = traps.buff, torus=T, timer = T)
  write.csv(detDat,row.names=FALSE,file=paste('Data/Detections/detDat_', scenario, '_', iter, '_G100.csv', sep=''))
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
