# ====================================
# Unmarked abundance evalaution
# REST model
# ====================================

library(dplyr)
library(tidyr)
library(foreach)
library(parallel)
library(doParallel)
library(nimble)
library(rgeos)
#library(extremevalues)
source("R/0_functions.R") 

'------ Run via shell script ------'
# get the input passed from the shell script
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# test if there is at least one argument: if not, return an error
if (length(args) <= 1) {
  #stop("At least two arguments must be supplied (input file).\n", call. = FALSE)
  scenario = 'small'; iter = 10
} else {
  print(paste0("Arg input: ", args[1], args[2], sep = ' '))
  scenario = args[1]
  iter = as.numeric(args[2])
}


# read input data (detections, movement sims)
moveSims <- read.csv(file=paste('Data/MovementSims/simDat_', scenario, '_', iter, '.csv', sep=''))
detDat <- read.csv(file=paste('Data/Detections/detDat_', scenario, '_', iter, '_G100.csv', sep=''))
load('Data/unmarkedEnvGrid100.RData')

"-------- Prep model input --------"
# calculate time diff between each detection
x1 <- detDat %>% 
  group_by(Site) %>%
  arrange(time)  %>% # sort encounters chronologically
  mutate( difftime = c(NA,diff(time))) %>%
  #filter(!difftime == 1 | is.na(difftime)) %>% # remove nonindependent records (within 30 mins) and retain NAs
  as.data.frame()

# calculate staying time
nbObs = nrow(moveSims[moveSims$ID == 1,])
x2 <- calcStay(movedf = moveSims, detdf = x1, X = traps.buff, timestep = 15, nbObs = nbObs)

# for consecutive captures (difftime = 1), consolidate into one capture
x3 <- x2 %>%
  group_by(Site) %>%
  arrange(time)  %>% # sort encounters chronologically
  # note that this only consolidates up to two consecutive captures
  mutate( stayTime1 = if_else( lead(difftime, n = 2) == 1 & lead(difftime, n = 1) == 1, # first check if both 1st and 2nd leading val is a consecutive capture
                               stayTime + lead(stayTime, n = 2), # if yes, add that row's staytime to current val
                               stayTime,
                               missing = stayTime) )  %>% 
  mutate( stayTime2 = if_else( lead(difftime, n = 1) == 1, 
                               stayTime1 + lead(stayTime1, n = 1), # then check if 1st leading val is a consecutive capture
                               stayTime1,
                               missing = stayTime1) )  %>% 
  filter(!difftime %in% 1:2) %>% # remove non-independent detections (within 30 mins)
  #select(ID, x, y, time, Site, difftime, stayTime = stayTime2) %>%
  as.data.frame()

## Data formatting
## y is a vector of no. of independent encounters per camera (1 * J)
## stay is a vector of observed staying times for every observed animal (1 * n_stay)
## day is a vector of days each camera was run

# prepare misc variables
Ncam <- length(traps.buff)   # Number of camera traps
A <- 1000*1000    # Study area in m2

# prepare y and stay
y1 <- aggregate(x3$Site, by = list(Site=x3$Site), FUN = length)
y2 <- data.frame(Site=1:Ncam, x = rep(0,Ncam))
y3 <- merge(y1, y2, by = 'Site')
y3 <- y1 %>% tidyr::complete(Site = 1:Ncam, fill = list(x=0) ) %>% as.data.frame()
y <- as.vector(y3[,2]) # ENCOUNTER DATA

stay <- x3$stayTime2 ; hist(stay, breaks = seq(0,ceiling(max(stay)/10)*10, 10))

# remove right outliers (determine censor time)
#outies <- getOutliersI(stay,rho = c(1,1),FLim = c(0.1,0.9), distribution = 'exponential')
#stayLim <- outies[["limit"]][["Right"]]
# stayLim <- round(quantile(stay, 0.95))
stayLim <- 15*60

## Calculate area of camera FOV
r <- 1 # detection radius
a <- 1/2 * 1 * 1 # all traps same FOV
#a <- pi*(r^2) # old area calc for circular det zone

## Calculate survey effort in seconds
day <- 50
actv <- 1 # Activity proportion (we tried 0.2, 0.4, 0.6, 0.8, 1) 
effort <- rep(60*60*12*day*actv,Ncam)    # Vector of effort (sec*min*hrs*day*activity proportion)

Nstay<-sum(y)    # Total independent encounters
cens <- rep(0,Nstay)    # Create vector to indicate if animal should be censored
cens[which(stay > stayLim)] <- stayLim    # Remove animal based on censoring rule 
stay[which(stay > stayLim)]<- stayLim    # Remove same animal from staying times

"-------- Specify model --------"
##################################################
code <- nimbleCode({
  lambda ~ dunif(0,5)
  D ~ dgamma(0.1,0.1)
  
  for (i in 1:Nstay){
    isCensored[i] ~ dinterval(stay[i], lim[i])  # Right censor data
    stay[i] ~ dexp(lambda)
  }
  
  for (i in 1:Ncam){
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(a) + log(effort[i]) + log(D) + log(lambda)
  }
  N <- D*A
})

##################################################

"-------- Prep NIMBLE input --------"
data <- list(y=y, stay=stay)
str(data)
const <- list(a=a, effort=effort, Nstay=Nstay, Ncam=Ncam, lim=cens, A=A)
str(const)
inits<-function(){
  list(lambda=1/8, D=0.05)
}
parameters<-c("lambda","D","N") # "mu"

"-------- Prep for parallel processing --------"
ncores = detectCores() - 3
cl<-makeCluster(ncores)
clusterSetRNGStream(cl = cl) # set up random number stream 
registerDoParallel(cl) # register as backend for 'foreach'


"-------- Execute REST --------"
system.time(
  res <- foreach(x = 1:ncores, .packages="nimble",.errorhandling='remove', .inorder=FALSE) %dopar% {
    set.seed(x)
    out <- nimbleMCMC(code, data=data, constants=const,
                      inits=inits, monitors=parameters,
                      niter=22000, nburnin=2000, thin=1, nchains=1, # niter=22000 takes ~20s
                      samplesAsCodaMCMC=TRUE)
  } )
# stop parallel
stopCluster(cl)

# Save output
jName <- paste("Data/ModelOut/REST_",scenario,'_',iter,".RData", sep="")
save(res, file = jName)
