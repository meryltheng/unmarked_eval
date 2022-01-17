# ====================================
# Unmarked abundance evalaution
# Movement simulations (simple movement)
# ====================================

library(CircStats)
# source("R/0_functions.R")

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

if (scenario == 'small'){
  scl <- 1; shp <- 1.6; trty <- 0.4; alpha <- 4 # AC strength
  minDist = 1 # min dist between ACs
  borderBuff = 5 # min dist between ACs and landscape border
  N = 1000
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  
  # 1.24   63.97   87.44   95.29  118.05  576.30  (steplength m/h)
  # 0.03109 0.06542 0.07954 0.08450 0.09647 0.22202 (MCP km^2; 50 days)
  # 0.05893 0.09450 0.10919 0.11385 0.12886 0.25433 (MCP km^2; 100 days)
}

if (scenario == 'medium'){
  scl <- 1.7; shp <- 3.2; trty <- 0.4; alpha <- 1.2 # AC strength
  minDist = 5 # min dist between ACs
  borderBuff = 25 # min dist between ACs and landscape border
  N = 200
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  29.57  250.29  312.60  326.31  387.12 1226.07 (steplength m/h)
  #  0.6338  0.7894  0.8339  0.9869  1.0412  1.8954 (MCP km^2; 50 days)
  #  0.9904  1.1548  1.3368  1.3281  1.4972  1.6534 (MCP km^2; 100 days)
}

if (scenario == 'large'){
  scl <- 3; shp <- 5; trty <- 0.6; alpha <- 0.6 # AC strength
  minDist = 75 # min dist between ACs
  borderBuff = 100 # min dist between ACs and landscape border
  N = 10
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   300.8   728.2   870.5   895.0  1033.6  2576.4 (steplength m/h)
  #   6.975   7.864  10.058   9.507  10.785  11.792 (MCP km^2; 50 days)
  #   9.157   9.931  11.399  11.303  12.831  13.112 (MCP km^2; 100 days)
}

"-------- Prep sim parameters --------"
# set general parameters
# nsim = 10
days = 100
nbObs <- 6*12*days # 15-min timesteps, 12-h activity, 50 days
xlim <- ylim <- c(0,1000) # landscape bounds

"-------- Simulation functions --------"
# generate ACs
genAC <- function(N = 100, minDist = 10, borderBuff = 0, xlim = c(0,1000), ylim = c(0,1000)){
  output <- as.data.frame(matrix(NA, N, 2))
  colnames(output) <- c('x','y')
  output[1,] <-  c(runif(1,xlim[1]+borderBuff,xlim[2]-borderBuff), 
               runif(1,ylim[1]+borderBuff,ylim[2]-borderBuff))
  
  for (i in 2:N){ 
    repeat{ # rejection sampling
      cand <- c(runif(1,xlim[1]+borderBuff,xlim[2]-borderBuff),
                runif(1,ylim[1]+borderBuff,ylim[2]-borderBuff))
      
      distapart <- apply(output, 1, function(x) sqrt( (x[1] - cand[1])^2 + (x[2] - cand[2])^2) )
      if (all(distapart > minDist, na.rm = T) ) { break }
    } # repeat
    # input
    output[i,] <- cand
  }
  return(output)
}

# simulation
simMove <- function(nbObs = nbObs, xy0 = c(0,0), ACstrength = 1){
  xy <- matrix(NA, nrow = nbObs, ncol = 2) # empty matrix for simulated locs
  xy <- data.frame(xy)
  names(xy) <- c("x", "y")
  xy[1, ] <- rnorm(2,xy0,5) # starting loc
  x0 <- xy0[1]
  y0 <- xy0[2]
  for (i in 2:nbObs) {
    if (i == 2){ # only for second location
      slen <- rgamma(20, scale = scl, shape = shp)
      ta <- runif(20, -pi, pi)
    } else { # for subsequent locations
      # cand locations
      #scl <- 1.4
      #shp <- 3.2
      ta <- rwrpnorm(20, ta[z], trty) # angle correlated to previous angle 
      slen <- 0.5*(slen[z] + rgamma(20, scale = scl, shape = shp)) # speed correlated to previous speed 
    }

    x1 <- xy[i - 1, 1] + (cos(ta) * slen) 
    y1 <- xy[i - 1, 2] + (sin(ta) * slen)
    
    # bias towards AC
    d1 <- vector() # distance to current AC, accounting for torus
    for (a in 1:20){ 
      d1[a] <- sqrt( min( abs(x1[a] - x0), 1000 - abs( x1[a] - x0 ))^2  +  min( abs(y1[a] - y0), 1000 - abs( y1[a] - y0 ))^2 )
    }
    w <- dexp(d1, rate = 1/25) 
    w[is.na(w)] <- 0
    w <- w^ACstrength * (10^(ACstrength-1)) # increase attraction to centre
    z <- sample(20, 1, prob = w)
    
    # landscape torus
    withinLandscape = isTRUE(x1[z] >= 0 & x1[z] <= 1000 & y1[z] >= 0 & y1[z] <= 1000)
    if (withinLandscape==FALSE){
      if(x1[z] <= 0){
        x1[z] <- 1000 + x1[z]
      }
      if(x1[z] >= 1000){
        x1[z] <- x1[z] - 1000 
      }
      if(y1[z] <= 0){
        y1[z] <- 1000 + y1[z]
      }
      if(y1[z] >= 1000){
        y1[z] <- y1[z] - 1000 
      }
    }
    # fill output
    xy[i, ] <- c(x1[z], y1[z])
  } # i 
  xy <- xy[1:nbObs,]
  xy$ID <- id
  xy$time <- 1:nrow(xy)
  return(xy)
}

"-------- Run movement simulations --------"
# ACs <- cbind(x = runif(N, 0, 1000), y = runif(N, 0, 1000))
ACs <- genAC(N = N, minDist = minDist, borderBuff = 100, xlim = xlim, ylim = ylim)
#points(ACs,pch=3)

# run once (test)
#xy <- simMove(xy0 = c(500,500), ACstrength = 4)

# run many
simDat <- list()
for (id in 1:N){
  simDat[[id]] <- simMove(nbObs = nbObs, xy0 = as.numeric(ACs[id,]), ACstrength = alpha)
  print(paste('Animal',id,'done'))
}

# save
moveSims <- do.call(rbind, simDat)
write.csv(moveSims,row.names=FALSE,file=paste('Data/MovementSims/simDat_', scenario, '_', iter, '.csv', sep=''))