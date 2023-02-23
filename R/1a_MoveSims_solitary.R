# ====================================
# Unmarked abundance evalaution
# Movement simulations (solitary movement)
# adapted from https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12155 
# ====================================

library(CircStats)
library(data.table)
source("0_functions.R")

'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for sceenario, arg[2] for iter)
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# test if there is at least one argument: if not, return an error
if (length(args) <= 1) {
  stop("At least two arguments must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input: ", args[1], args[2], sep = ' '))
  scenario = args[1]
  iter = as.numeric(args[2])
}

if (scenario == 'hi'){
  # step-length params
  aind <- c(48) # ind mean
  bind <- c(30) # ind sd
  # turning ang
  kappac <- c(0.7) # directional concentration/attraction strength
  etac <- c(0.4) # centroid weight of CRW component in BCRW for centroid movement
  
  minDist = 0 # min dist between ACs
  borderBuff = 50 # min dist between ACs and landscape border
  N = 2000
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  
  # 0.0951  0.1120  0.1213  0.1235  0.1300  0.1707  (MCP km^2; 50 days)
}

if (scenario == 'med'){
  # step-length params
  aind <- c(100) # ind mean
  bind <- c(62.5) # ind sd
  # turning ang
  kappac <- c(0.65) # directional concentration/attraction strength
  etac <- c(0.5) # centroid weight of CRW component in BCRW for centroid movement
  
  minDist = 50 # min dist between ACs
  borderBuff = 250 # min dist between ACs and landscape border
  N = 200
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  0.9721  1.2556  1.4242  1.4244  1.5631  2.7685  (MCP km^2; 50 days)
}

if (scenario == 'medfast'){
  # step-length params
  aind <- c(180) # ind mean
  bind <- c(120) # ind sd
  # turning ang
  kappac <- c(0.65) # directional concentration/attraction strength
  etac <- c(0.5) # centroid weight of CRW component in BCRW for centroid movement
  
  minDist = 50 # min dist between ACs
  borderBuff = 250 # min dist between ACs and landscape border
  N = 200
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  0.9721  1.2556  1.4242  1.4244  1.5631  2.7685  (MCP km^2; 50 days)
}

if (scenario == 'lo'){
  # step-length params
  aind <- c(180) # ind mean
  bind <- c(120) # ind sd
  # turning ang
  kappac <- c(0.65) # directional concentration/attraction strength
  etac <- c(0.6) # centroid weight of CRW component in BCRW for centroid movement
  
  minDist = 750 # min dist between ACs
  borderBuff = 1000 # min dist between ACs and landscape border
  N = 10
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   10.93   13.74   15.15   16.79   17.21   27.29  (MCP km^2; 50 days; 2-day locs)
}

"-------- Prep sim parameters --------"
# set general parameters
days = 50
nbObs <- 12*9*days # 5-min timesteps, 9-h activity, 50 days
xlim <- ylim <- c(0,10000) # landscape bounds

"-------- Run movement simulations --------"
# activity centres
ACs <- genAC(N = N, minDist = minDist, borderBuff = borderBuff, xlim = xlim, ylim = ylim)

# # run once (test)
# xy <- simMove(nbObs = nbObs, xy0 = as.numeric(ACs[id,]), ID = 1, burnin = 12*12*1, torus = T)
# buff=1200
# plot(1, type='n', xlim = c(ACs[id,1] - buff, ACs[id,1] + buff), ylim = c(ACs[id,2] - buff, ACs[id,2] + buff))
# plotTraject(df = xy, plotit = 1, torus = T)

# run many
simDat <- list()
for (id in 1:N){
  simDat[[id]] <- simMove(nbObs = nbObs, burnin = 12*9*1, xy0 = as.numeric(ACs[id,]), ID = id,
                          aind = aind, bind = bind, # step length params
                          kappac = kappac, etac = etac, # turning angle params
                          xlim = c(0,10000), ylim = c(0,10000), torus = T)
  print(paste('Animal',id,'done'))
}

# compile into dataframe
moveSims <- do.call(rbind, simDat)

# save
write.csv(moveSims,row.names=FALSE,file=paste('Data/MovementSims/simDat_', scenario, '_sol_', iter, '.csv', sep=''))
