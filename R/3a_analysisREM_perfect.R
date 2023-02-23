# ====================================
# Unmarked abundance evalaution
# REM model (for HPC)
# PERFECT DETECTION (r = 1 unit)
# ====================================

library(dplyr)
library(tidyr)
library(activity)
library(data.table)
source("R/camtools.R") # https://github.com/MarcusRowcliffe/camtools/blob/master/camtools.R
source("R/sbd.R") # https://github.com/MarcusRowcliffe/sbd/blob/master/sbd.r

# function to calculate exavt speed of overlapped trajectory
calcSpeed <- function(movedf = moveSims, detdf = x1, X = traps.buff, timestep = 15){
  require(sp)
  require(rgeos)
  detdf$dist <- detdf$stayTime <- rep(NA, nrow(detdf))
  
  polyID <- sapply(slot(X, "polygons"), function(x) slot(x, "ID"))
  
  for (i in 1: nrow(detdf)){
    t <- as.numeric(detdf[i,'time']); id <- as.numeric(detdf[i,'ID']); detector <- as.numeric(detdf[i,'Site'])
    
    track <- movedf %>%
      filter(ID == id & time %in% (t-1):t ) %>%  as.data.frame()
    coords <- coordinates(track[,c('x','y')])
    trackLine <- SpatialLines(list(Lines(list( Line(coords)), "id" )))
    
    detectorID <- which(polyID == detector)
    
    # while speed can directly be calculated from gLength(trackLine) since time units are constant in sims,
    # this matters for aggregating consecutive encounters which means more than one straight-line track
    trackIntersect <- gIntersection(X[detectorID,],trackLine) # this calculates distance
    detdf$dist[i] <- gLength(trackIntersect)
    detdf$stayTime[i] <- timestep * 60 * gLength(trackIntersect)/gLength(trackLine) # in seconds (* 60 cuz timestep is in mins)
    
  }
  return(detdf)
}

'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for scenario, arg[2] for iter, arg[3] for nTraps, arg[4] for movemt behav (sol/grp), if grp arg[5] is group size)
# e.g., `Rscript --vanilla ./R/3c_analysisREM.R small 1 100 sol`
# e.g., `Rscript --vanilla ./R/3c_analysisREM.R small 1 100 grp 25`
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

# additional metas
subsetDat = TRUE
days=50

# read input data (detections, trap data)
if (behav == 'sol'){
  moveSims <- fread(file=paste('Data/MovementSims/simDat_', 
                               scenario, '_sol_', iter, '.csv', sep=''))
  detDat <- fread(file=paste('Data/Detections/detDat_', scenario, '_sol_', iter, '_J', nTraps,'.csv', sep=''))
  groupSize = NA
}
if (behav == 'grp'){
  detDat <- fread(file=paste('Data/Detections/detDat_', scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,'.csv', sep=''))
  moveSims <- fread(file=paste('Data/MovementSims/simDat_', 
                               scenario, '_grp_', iter, '_gSize', groupSize,'.csv', sep=''))
}

if(subsetDat == TRUE){
  moveSims <- 
    moveSims %>%
    filter(time %in% 1:(12 * 9 * days)) 
  detDat <- 
    detDat %>%
    filter(time %in% 1:(12 * 9 * days)) 
}

load(paste0('Data/traps',nTraps,'.RData'))

"-------- Prep model input --------"
x1 <- calcSpeed(movedf = moveSims, detdf = detDat, X = traps[["trapPoly"]], timestep = 5)

# first determine the consecutive tracks to evaluate
x2 <- x1 %>%
  group_by(Site,ID) %>%
  arrange(time)  %>% # sort encounters chronologically
  mutate( difftime = c(NA,diff(time))) %>%
  mutate( evaluate = if_else( lead(difftime == 1,n=1) | lead(difftime == 1,n=2), 1, 0, 0)) %>% # evaluate up to 2 consec tracks
  as.data.frame()

# second evaluate if the consecutive tracks fall within the FOV
x2$collapse <- NA
for (i in 1:nrow(x1)){
  if(x2$evaluate[i]==1){
    x2$collapse[i] <- over(traps[["trapPoly"]][x2$Site[i],], 
                           SpatialPoints(coordinates(data.frame('x' = x2$x[i], 'y' = x2$y[i]))))
  } else {next}
}

# third collapse consecutive tracks into an encounter and calculate speeds for each encounter
input_y <- x2 %>%
  group_by(Site,ID) %>%
  arrange(time)  %>% # sort encounters chronologically
  mutate(time_new = if_else(lag(collapse,n=1)==1, lag(time,n=1), time, time)) %>%
  group_by(Site,ID,time_new) %>%
  summarise(ID = unique(ID), time_new = unique(time_new), Site = unique(Site), speed = sum(dist)/sum(stayTime)) %>%
  
  # group_by(Site,ID) %>%
  # arrange(time_new)  %>% # sort encounters chronologically
  # mutate( difftime = c(NA,diff(time_new))) %>%
  # filter(difftime != 1 & difftime != 2 | is.na(difftime)) %>%
  
  as.data.frame()

"-------- Estimate model params --------"
# speed
speed_mods <- sbm3(speed ~1, input_y[input_y$speed > 0.01 & input_y$speed < 10,], reps = 1000)
best_speed_mod <- as.name(row.names(speed_mods[["AICtab"]])
                          [which.min(speed_mods[["AICtab"]]$AIC)])
pred_speed <- predict.sbm(speed_mods[["models"]][[best_speed_mod]])
# pred_speed$est * 9 * 60 ^ 2
# plot.sbm(speed_mods[["models"]][[best_speed_mod]], log = T, xlab="Speed (log[unit/s])")

# encounters
y1 <- aggregate(input_y$Site, by = list(Site=input_y$Site), FUN = length) 
y2 <- data.frame(Site=1:nTraps, x = rep(0,nTraps))
y3 <- y1 %>% tidyr::complete(Site = 1:nTraps, fill = list(x=0) ) %>% as.data.frame()
y <- as.vector(y3[,2])

# REM
param <- list(v = pred_speed$est * 9 * 60 * 60, # in seconds, based on 12-h activity day
              p = 1,
              r = 10,
              theta = pi/180 * 40 )
paramse <- list(v = pred_speed$se * 9 * 60 * 60,
                p = 0,
                r = 0,
                theta = 0)

"-------- Execute REM --------"
D_hat <- bootTRD(y, rep(days, nTraps), param, paramse)
A <- 1000*1000  # 1 km^2
( N_hat <- D_hat * A )# density estimate /km^2

"-------- Save/write output --------" 
Replicate = iter
Scenario = scenario
Behaviour = behav
if (behav == 'sol'){ groupSize = NA }
Model = "REM"
Trap.effort = nTraps
Days.monitored = days
nDet = NA
Dhat = N_hat[1]
Dhat.se = N_hat[2]
Dhat.lo = NA
Dhat.hi = NA
y = nrow(input_y)
v = pred_speed$est * 9 * 60 * 60
v.se =  pred_speed$se * 9 * 60 * 60
act = 1 # activity level within the 12-h day; actmod@act["act"]
act.se = 0 # actmod@act["se"]
edr.model = NA
edr = 10
edr.se = 0
edr.AIC = NA

rem_output <- data.frame(Replicate, Scenario, Behaviour, groupSize, Model, Trap.effort, Days.monitored, nDet,
                         Dhat, Dhat.se, Dhat.lo, Dhat.hi, y, v, v.se, act, act.se, edr.model, edr, edr.se, edr.AIC)

# append to existing master file
write.table( rem_output,  
             file="./Data/Estimates/rem_perfectDet.csv", 
             append = T, 
             sep=',', 
             row.names=F, 
             col.names=F )

# # create empty .csv if not avail
# columns = c("Replicate", "Scenario", "Behaviour", "groupSize", "Model", "Trap.effort", "Days.monitored", "nDet",
#            "Dhat", "Dhat.se", "Dhat.lo", "Dhat.hi", "y", "v", "v.se", "act", "act.se", "edr.model", "edr", "edr.se", "edr.AIC")
# mydf = data.frame(matrix(nrow = 0, ncol = length(columns)))
# colnames(mydf) = columns
# write.csv(mydf, file=paste('Data/Estimates/rem_perfectDet.csv',sep=''), row.names = F)

# "-------- Post-hoc --------"
# # group numbers
# group_nos <- input_y %>%
#   group_by(Site, time) %>%
#   summarise( Site = unique(Site), time = unique(time), n_ind = length(ID) )  %>% 
#   as.data.frame()
# 
# hist(group_nos$n_ind)