# ====================================
# Unmarked abundance evalaution
# SC model
# ====================================

library(sp)
library(dplyr)
library(tidyr)
library(foreach)
library(parallel)
library(doParallel)
library(nimble) # to run via shell, might have to add R.home("bin")'s path to PATH variable

'------ Run via shell script ------'

# get the input passed from the shell script (arg[1] for sceenario, arg[2] for iter)
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# test if there is at least one argument: if not, return an error
if (length(args) <= 1) {
  # stop("At least two arguments must be supplied (input file).\n", call. = FALSE)
  scenario = 'slow'; iter = 1
} else {
  print(paste0("Arg input: ", args[1], args[2], sep = ' '))
  scenario = args[1]
  iter = as.numeric(args[2])
}

# read input data (detections)
detDat <- read.csv(file=paste('Data/Detections/detDat_', scenario, '_', iter, '.csv', sep=''))
head(detDat)
load('Data/unmarkedEnv.RData')

"-------- Prep model input --------"
Ncam = length(traps.buff)

# calculate time diff between each detection
x1 <- detDat %>% 
  group_by(Site) %>%
  arrange(time)  %>% # sort encounters chronologically
  mutate( difftime = c(NA,diff(time))) %>%
  filter(!difftime %in% 1:2 | is.na(difftime)) %>% # remove nonindependent records (within 30 mins) and retain NAs
  as.data.frame()

head(x1)

# prepare y; counts of animals at each trap j across the entire study period (nj.)
y1 <- aggregate(x1$Site, by = list(Site=x1$Site), FUN = length)
y2 <- data.frame(Site=1:Ncam, x = rep(0,Ncam))
y3 <- merge(y1, y2, by = 'Site')
y3 <- y1 %>% complete(Site = 1:Ncam, fill = list(x=0) ) %>% as.data.frame()
y <- as.vector(y3[,2]) # ENCOUNTER DATA

head(y)

## Data formatting

"-------- Specify model --------"
##################################################
code <- nimbleCode({
  lam0 ~ dunif(0,3)    # base encounter rate
  sigma ~ dunif(0,5)    # decay parameter
  psi ~ dunif(0,1)    # inclusion parameter
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(Xl,Xu)    # x coordinates of activity centers
    s[i,2] ~ dunif(Yl,Yu)    # y coordinates of activity centers
    
    for(j in 1:J) {
      ## Euclidean distances between individual activity centers and
      ## camera trap locations
      d2[i,j] <- pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2)
      ## Expected number of observations of an individual at a trap
      lam[i,j] <- lam0*exp(-(d2[i,j])/(2*sigma*sigma))*z[i]
    }
  }
  
  for(j in 1:J) {
    ## Expected observations over all occasions at a trap
    bigLambda[j] <- sum(lam[1:M,j])
    n[j] ~ dpois(bigLambda[j]*K[j])    # Days camera active
  }
  
  N <- sum(z[])    # Abundance
  A <- Xu*Yu    # State space area (or can define if known as in our case)
  D <- N/A    # Realized density
})
##################################################

"-------- Prep NIMBLE input --------"
X <- coordinates(traps.buff) # CT coordinates
J = nrow(X)    # Number of camera traps ever deployed
K <- rep(50,J)   # Number of days each trap was active
M = 120   # Data augmentation parameter for the encounter histories
Xl=0; Yl=0; Xu=50; Yu=50 # landscape bounds

data <- list(n=y) 
const <- list(M=M, K=K, J=J, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X/10)
parameters<-c('psi','lam0','N','sigma')

# set inital values
N <- round(rnorm(1, 50, 5)) # assumed abundance
Sx=runif(M,Xl,Xu) # AC locations
Sy=runif(M,Yl,Xu) 
inits <- function() {list(z=c(rep(1, N), rep(0,M-N)), psi=runif(1), s=cbind(Sx,Sy), sigma=3, lam0=0)}

"-------- Prep for parallel processing --------"
ncores = detectCores() - 3
cl<-makeCluster(ncores)
clusterSetRNGStream(cl = cl) # set up random number stream 
registerDoParallel(cl) # register as backend for 'foreach'

"-------- Execute SC --------"
system.time(
  res <- foreach(x = 1:ncores, .packages="nimble",.errorhandling='remove', .inorder=FALSE, .verbose = T) %dopar% {
    set.seed(x)
    out <- nimbleMCMC(code, data=data, constants=const,
                      inits=inits, monitors=parameters,
                      niter=22000, nburnin=2000, thin=1, nchains=1, # niter=22000 takes...
                      samplesAsCodaMCMC=TRUE)
  } )
#printErrors()

# stop parallel
stopCluster(cl)

# Save output
jName <- paste("Data/ModelOut/SC_",scenario,'_',iter,".RData", sep="")
save(res, file = jName)
