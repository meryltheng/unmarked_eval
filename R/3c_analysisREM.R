# ====================================
# Unmarked abundance evalaution
# REM model
# ====================================

library(dplyr)
library(tidyr)
library(rgeos)
library(activity)
source("R/0_functions.R") 
source("R/camtools.R") # https://github.com/MarcusRowcliffe/camtools/blob/master/camtools.R
source("R/sbd.R") # https://github.com/MarcusRowcliffe/sbd/blob/master/sbd.r

'------ Run in loop ------'
scenario  = 'medium'
trueN = 200
subsetDat = TRUE
nTraps = 100

# empty matrix for results
nsims = 100
output <- as.data.frame(matrix(data=NA, nrow=nsims, ncol=34))
colnames(output) <- c('Replicate', 'Scenario', 'N','Model', 'Effort',# meta info
                      'Nhat', 'Nhat.se', 'Nhat.lo', 'Nhat.hi','Nhat.mode','Nhat.Rhat', # abundance
                      'nDet', # detections
                      # SC
                      'sigma', 'sigma.se','sigma.mode','sigma.Rhat',
                      'lam0', 'lam0.se','lam0.mode','lam0.Rhat',
                      'psi', 'psi.se','psi.mode','psi.Rhat',
                      # REST
                      'Dhat', 'Dhat.se','Dhat.Rhat',
                      'lamREST', 'lamREST.se','lamREST.Rhat',
                      #REM
                      'v', 'v.se','act', 'act.se'
)

for (iter in 1:nsims){
  # read input data (detections, movement sims)
  moveSims <- read.csv(file=paste('Data/MovementSims/simDat_', scenario, '_', iter, '.csv', sep=''))
  detDat <- read.csv(file=paste('Data/Detections/detDat_', scenario, '_', iter, '_G',nTraps,'.csv', sep=''))
  load(paste('Data/unmarkedEnvGrid',nTraps,'.RData',sep = ''))
  
  #colnames(moveSims)[3] <- 'time'
  
  if(subsetDat == TRUE){
    moveSims <- 
      moveSims %>% 
      group_by(ID) %>% 
      filter(row_number() %in% 1:(6 * 12 * 50))  %>% 
      as.data.frame()
    
    detDat <- 
      detDat %>%
      filter(time %in% 1:(6 * 12 * 50)) 
  }
  
  "-------- Prep model input --------"
  # calculate time diff between each detection
  x1 <- detDat %>% 
    group_by(Site) %>%
    arrange(time)  %>% # sort encounters chronologically
    mutate( difftime = c(NA,diff(time))) %>%
    #filter(!difftime == 1 | is.na(difftime)) %>% # remove nonindependent records (within 30 mins) and retain NAs
    as.data.frame()
  
  # calculate speed
  nbObs = nrow(moveSims[moveSims$ID == 1,])
  x2 <- calcSpeed(movedf = moveSims, detdf = x1, X = traps.buff, timestep = 10, nbObs = nbObs) # moveSims[,-4]
  
  # for consecutive captures (difftime = 1), consolidate into one capture
  # note that multiple individuals detected at the same time (and site) are treated as separate encounters (difftime=0), thus I do not need to account for group size
  x3 <- x2 %>%
    group_by(Site) %>%
    arrange(time)  %>% # sort encounters chronologically
    # note that this only consolidates up to three consecutive captures
    mutate( speed1 = if_else( lead(difftime, n = 1) == 1, 
                                 (dist + lead(dist, n = 1))/
                                 (stayTime + lead(stayTime, n = 1)), # then check if 1st leading val is a consecutive capture
                                 dist/stayTime,
                                 missing = dist/stayTime) )  %>% 
    mutate( speed2 = if_else( lead(difftime, n = 2) == 1 & lead(difftime, n = 1) == 1, # first check if both 1st and 2nd leading val is a consecutive capture
                              (dist + lead(dist, n = 1) + lead(dist, n = 2))/
                                (stayTime + lead(stayTime, n = 1) + lead(stayTime, n = 2)), # if yes, add up all dist and stayTimes to calc overall speed
                              speed1,
                              missing = speed1) )  %>% 
    mutate( speed = if_else( lead(difftime, n = 3) == 1 & lead(difftime, n = 2) == 1 & lead(difftime, n = 1) == 1, # first check if both 1st and 2nd leading val is a consecutive capture
                             (dist + lead(dist, n = 1) + lead(dist, n = 2) + lead(dist, n = 3))/
                               (stayTime + lead(stayTime, n = 1) + lead(stayTime, n = 2) + lead(stayTime, n = 3)), # if yes, add up all dist and stayTimes to calc overall speed
                             speed2,
                             missing = speed2) )  %>% 
    filter(!difftime %in% c(1,2)) %>% # remove consecutive detections (within 15 mins of each other)
    as.data.frame()
  
  speed <- x3$speed; hist(log10(speed), main="", xlab="Speed (log10[unit/s])") # remove outliers if needed
  
  # activity level estimation
  # (hour(timez) + minute(timez)/60 + second(timez)/60^2) * pi/12
  x3$day <- ceiling(x3$time/72)
  x3$hour <- ceiling((x3$time - (x3$day-1)*72)/6)
  actmod <- fitact(x3$hour * pi/12,
                   bounds=c(1,12)*pi/12, sample="data", reps=100)
  # plot(actmod, centre="night", dline=list(col="grey"))
  
  # REM (camtools)
  Ncam <- length(traps.buff)   # Number of camera traps
  A <- 1000*1000    # Study area in m2
  r <- 1
  theta <- pi/180 * 60
  v <- hmean(speed)
  
  # prepare encounter data by site
  y1 <- aggregate(x3$Site, by = list(Site=x3$Site), FUN = length)
  y2 <- data.frame(Site=1:Ncam, x = rep(0,Ncam))
  y3 <- y1 %>% tidyr::complete(Site = 1:Ncam, fill = list(x=0) ) %>% as.data.frame()
  y <- as.vector(y3[,2]) # ENCOUNTER DATA
  
  
  # prep params
  param <- list(v = v['mean'] * 12 * 60 * 60, # in seconds, based on 12-h activity day
                p = actmod@act["act"], # activity level within the 12-h day
                r = r,
                theta = theta * 2)
  paramse <- list(v = v['se'] * 12 * 60 * 60,
                  p = actmod@act["se"],
                  r = 0,
                  theta = 0)
  
  "-------- Execute REM --------"
  Dhat <- bootTRD(y, rep(50, Ncam), param, paramse) # encounters, sampling duration per cam, etc.
  ( Nhat <- Dhat * A )# abundance estimate
  
  "-------- Fill results matrix --------"
  # meta info
  output[iter,'Model'] <- 'REM'
  output[iter,'Replicate'] <- iter
  output[iter,'Scenario'] <- scenario
  output[iter,'N'] <- trueN
  output[iter,'Effort'] <- nTraps
  output[iter,'nDet'] <- nrow(x3) # number of 'independent' detections
  # estimates
  output[iter,'Dhat'] <- Dhat[1]
  output[iter,'Dhat.se'] <- Dhat[2]
  output[iter,'Nhat'] <- Nhat[1]
  output[iter,'Nhat.se'] <- Nhat[2]
  output[iter,'v'] <- v['mean'] * 12 * 60 * 60 # day range while active
  output[iter,'v.se'] <- v['se'] * 12 * 60 * 60
  output[iter,'act'] <-actmod@act["act"] # activity level within the 12-h day
  output[iter,'act.se'] <-actmod@act["se"]
}

#resultHab <- output
resultNull <- output
result <- rbind(resultHab,resultNull)


write.csv(result, file=paste('Data/Estimates/REMout_med.csv',sep=''), row.names = F)


#######################################################################################
# REM (raw)
#Ncam <- length(traps.buff)   # Number of camera traps
#y <- nrow(x3)
#t <- Ncam * 50 # days
#v <- hmean(speed)['mean'] * 12 * 60 * 60

#D <- (y/t)*pi/(v *r*(2+theta))
#D*A

