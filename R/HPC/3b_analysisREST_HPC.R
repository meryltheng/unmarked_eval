# ====================================
# Unmarked abundance evalaution
# REST model (HPC)
# ====================================
library(dplyr)
library(tidyr)
library(foreach)
library(parallel)
library(doParallel)
library(nimble)
#library(rgeos)
library(mcmcOutput)
#library(extremevalues)
#source("0_functions.R") 

'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for scenario, arg[2] for iter, arg[3] for nTraps, arg[4] for movemt behav (sol/grp), if grp arg[5] is group size)
# e.g., `Rscript --vanilla ./R/HPC/3b_analysisREST_HPC.R small 1 100 sol`
# e.g., `Rscript --vanilla ./R/3d_analysisDS.R small 1 100 grp 25`
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
  x4 <- read.csv(file=paste('Data/Detections/ctDat_', scenario, '_sol_', iter, '_J', nTraps,'.csv', sep=''))
}
if (behav == 'grp'){
  x4 <- read.csv(file=paste('Data/Detections/ctDat_', scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,'.csv', sep=''))
}

if (nTraps == 100){
  load('Data/traps100.RData')
}
if (nTraps == 25){
  load('Data/traps25.RData')
}

if(subsetDat == TRUE){
  x4 <- 
    x4 %>%
    filter(time %in% 1:(12 * 12 * days)) 
}

"-------- Prep model input --------"
# 6 consolidate independent detections (det_id) and calc input params (for REST)
x6 <- x4 %>%
  filter( dist >= 0.135 & dist <= 0.235 ) %>% # perfect detection distance
  # group_by(det_id) %>%
  group_by(track) %>% # because animal might enter and exit more than once within the same det_id
  summarise(ID = unique(ID), Site = unique(Site), time = unique(time),
            stayTime =  length(time_seconds) # in seconds; stayTime =  last(time_seconds)-first(time_seconds)
  ) %>%
  # group_by(det_id) %>%
  group_by(track) %>%
  filter(row_number()==1) %>% # remove repeats 
  as.data.frame()

## Data formatting
## y is a vector of no. of independent encounters per camera (1 * J)
## stay is a vector of observed staying times for every observed animal (1 * n_stay)
## day is a vector of days each camera was run

# prepare misc variables
A <- 1000*1000    # Study area in m2

# prepare y and stay
y1 <- aggregate(x6$Site, by = list(Site=x6$Site), FUN = length)
y2 <- data.frame(Site=1:nTraps, x = rep(0,nTraps))
y3 <- merge(y1, y2, by = 'Site')
y3 <- y1 %>% tidyr::complete(Site = 1:nTraps, fill = list(x=0) ) %>% as.data.frame()
y <- as.vector(y3[,2]) # ENCOUNTER DATA

stay <- x6$stayTime ; hist(stay, breaks = seq(0,ceiling(max(stay)/10)*10, 10))

# remove right outliers (determine censor time)
#outies <- getOutliersI(stay,rho = c(1,1),FLim = c(0.1,0.9), distribution = 'exponential')
#stayLim <- outies[["limit"]][["Right"]]
stayLim <- max(stay) # THIS WILL AFFECT N_hat; round(quantile(stay, 0.95))

## Calculate area of camera FOV
a <- pi * (0.235^2) / 6 - pi * (0.135^2) / 6# all traps same FOV

## Calculate survey effort in seconds
day <- 50
actv <- 1 # Activity proportion (we tried 0.2, 0.4, 0.6, 0.8, 1) 
effort <- rep(60*60*12*day*actv,nTraps)    # Vector of effort (sec*min*hrs*day*activity proportion)

Nstay<-sum(y)    # Total independent encounters
cens <- rep(0,Nstay)    # Create vector to indicate if animal should be censored
cens[which(stay > stayLim)] <- stayLim    # Remove animal based on censoring rule 
stay[which(stay > stayLim)] <- stayLim    # Remove same animal from staying times

"-------- Specify model --------"
##################################################
code <- nimbleCode({
  lambda ~ dunif(0,5)
  D ~ dgamma(0.1,0.1) # Density
  r ~ dunif(0,10) # size param for dnegbin
  
  for (i in 1:Nstay){
    isCensored[i] ~ dinterval(stay[i], lim[i])  # Right censor data
    stay[i] ~ dexp(lambda)
  }
  
  for (i in 1:Ncam){
    #y[i] ~ dpois(mu[i]) # Poission distribution
    y[i] ~ dnegbin (p[i], r) # Negative binomial distribution
    p[i] <- r / (r + mu[i])
    log(mu[i]) <- log(a) + log(effort[i]) + log(D) + log(lambda)
  }
  N <- D*A
})

##################################################

"-------- Prep NIMBLE input --------"
data <- list(y=y, stay=stay)
str(data)
const <- list(a=a, effort=effort, Nstay=Nstay, Ncam=nTraps, lim=cens, A=A)
str(const)
inits<-function(){
  list(lambda=1/8, D=0.05)
}
parameters<-c("lambda","D","N","r") # "mu"

"-------- Prep for parallel processing --------"
slurm_ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.numeric(slurm_ncores)) {
  cores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
} else {
  cores = detectCores()
}
cl<-makeCluster(cores)

clusterSetRNGStream(cl = cl) # set up random number stream 
registerDoParallel(cl) # register as backend for 'foreach'

"-------- Execute REST --------"
system.time(
  res <- foreach(x = 1:cores, .packages="nimble",.errorhandling='remove', .inorder=FALSE) %do% {
    set.seed(x)
    out <- nimbleMCMC(code, data=data, constants=const,
                      inits=inits, monitors=parameters,
                      niter=22000, nburnin=2000, thin=1, nchains=1, # niter=22000 takes ~20s
                      samplesAsCodaMCMC=TRUE)
  } )
# stop parallel
stopCluster(cl)

# # to test for compilation error
# test <- nimbleModel(code = code, data=data, constants=const, inits=inits(), name='test')
# Ctest <- compileNimble(test, showCompilerOutput = T)

"-------- Save/write output --------" 
# Save output
if (behav == 'sol'){
  jName=paste('Data/ModelOut/REST_', scenario, '_sol_', iter, '_J', nTraps,'.RData', sep='')
}
if (behav == 'grp'){
  jName=paste('Data/ModelOut/REST_', scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,'.RData', sep='')
}

save(res, file = jName)

# Convert to an mcmcOutput object and look at diagnostics
mclist <- coda::mcmc.list(res)
( mco <- mcmcOutput(mclist) )
summMCMC <- summary(mco)

# Write output
Replicate = iter
Scenario = scenario
Behaviour = behav
if (behav == 'sol'){ groupSize = NA }
Model = "REST"
Trap.effort = nTraps
Days.monitored = days
Nhat = unlist(summMCMC['N','mean'])
Nhat.sd = unlist(summMCMC['N','sd'])
Nhat.lo = unlist(summMCMC['N','l95'])
Nhat.hi = unlist(summMCMC['N','u95'])
Nhat.Rhat = unlist(summMCMC['N','Rhat'])
lambda = unlist(summMCMC['lambda','mean'])
lambda.Rhat = unlist(summMCMC['lambda','Rhat'])
r = unlist(summMCMC['r','mean'])
r.Rhat = unlist(summMCMC['r','Rhat'])
stayLim = stayLim

rest_output <- data.frame(Replicate, Scenario, Behaviour, groupSize, Model, Trap.effort, Days.monitored,
                        Nhat, Nhat.sd, Nhat.lo, Nhat.hi, Nhat.Rhat, lambda, lambda.Rhat, r, r.Rhat, stayLim)

# append to existing master file
write.table( rest_output,  
             file="./Data/Estimates/rest_master_test.csv", 
             append = T, 
             sep=',', 
             row.names=F, 
             col.names=F )

# # create empty .csv if not avail
# columns= c("Replicate", "Scenario", "Behaviour", "groupSize", "Model", "Trap.effort", "Days.monitored",
#            "Nhat", "Nhat.sd", "Nhat.lo", "Nhat.hi", "Nhat.Rhat", "lambda", "lambda.Rhat", "r", "r.Rhat", "stayLim")
# mydf = data.frame(matrix(nrow = 0, ncol = length(columns)))
# colnames(mydf) = columns
# write.csv(mydf, file=paste('Data/Estimates/rest_master_test.csv',sep=''), row.names = F)
