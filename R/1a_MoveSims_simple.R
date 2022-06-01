# ====================================
# Unmarked abundance evalaution
# Movement simulations (solitary movement)
# ====================================

library(CircStats)
# source("R/0_functions.R")

'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for scenario, arg[2] for iter)
# e.g., `Rscript --vanilla ./R/1a_MoveSims_simple.R small 1`
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
  scl <- 1.2; shp <- 2.1; trty <- 0.4; alpha <- 4.1 # AC strength
  minDist = 0 # min dist between ACs
  borderBuff = 5 # min dist between ACs and landscape border
  N = 2000
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  
  # 8.79  212.72  279.27  298.34  363.18 1856.11 (steplength m/h)
  # 0.03492 0.06850 0.08053 0.08318 0.09543 0.20110 (MCP km^2; 50 days)
}

if (scenario == 'medium'){
  scl <- 1.7; shp <- 3.2; trty <- 0.4; alpha <- 1.2 # AC strength
  minDist = 5 # min dist between ACs
  borderBuff = 20 # min dist between ACs and landscape border
  N = 200
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  65.86  499.28  622.16  650.22  770.74 2517.14 (steplength m/h)
  #  0.4813  0.8027  0.9364  0.9824  1.1191  1.7673  (MCP km^2; 50 days)
}

if (scenario == 'large'){
  scl <- 3; shp <- 5; trty <- 0.6; alpha <- 0.6 # AC strength
  minDist = 75 # min dist between ACs
  borderBuff = 100 # min dist between ACs and landscape border
  N = 10
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   300.8   728.2   870.5   895.0  1033.6  2576.4 (steplength m/h) [x2, timestep now 5min]
  #   6.975   7.864  10.058   9.507  10.785  11.792 (MCP km^2; 50 days)
  #   9.157   9.931  11.399  11.303  12.831  13.112 (MCP km^2; 100 days)
}

"-------- Prep sim parameters --------"
# set general parameters
days = 100
nbObs <- 12*12*days # 5-min timesteps, 12-h activity, 50 days
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
simMove <- function(nbObs = nbObs, xy0 = c(0,0), ACstrength = 1, torus = T){
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
      ta <- rwrpnorm(20, ta[z], trty) # angle correlated to previous angle 
      slen <- 0.5*(slen[z] + rgamma(20, scale = scl, shape = shp)) # speed correlated to previous speed 
    }

    x1 <- xy[i - 1, 1] + (cos(ta) * slen) 
    y1 <- xy[i - 1, 2] + (sin(ta) * slen)
    
    # bias towards AC
    d1<- sqrt( (x1 - x0 )^2  +   (y1 - y0)^2 )
    
    w <- dexp(d1, rate = 1/25) 
    w[is.na(w)] <- 0
    w <- w^ACstrength * (10^(ACstrength-1)) # increase attraction to centre
    z <- sample(20, 1, prob = w)
    
    # fill output
    xy[i, ] <- c(x1[z], y1[z])
  } # i 
  xy <- xy[1:nbObs,]
  xy$ID <- id
  xy$time <- 1:nrow(xy)
  
  # landscape torus
  if (torus == T){
    for (j in 1:nrow(xy)){
      if(xy$x[j] <= xlim[1]){
        xy$x[j] <- xlim[2] + xy$x[j]
      }
      if(xy$x[j] >= xlim[2]){
        xy$x[j] <- xy$x[j] - xlim[2]
      }
      if(xy$y[j] <= ylim[1]){
        xy$y[j] <- ylim[2] + xy$y[j]
      }
      if(xy$y[j] >= ylim[2]){
        xy$y[j] <- xy$y[j] - ylim[2] 
      }
    }
  }
  return(xy)
}

"-------- Run movement simulations --------"
# ACs <- cbind(x = runif(N, 0, 1000), y = runif(N, 0, 1000))
ACs <- genAC(N = N, minDist = minDist, borderBuff = borderBuff, xlim = xlim, ylim = ylim)

# run many
simDat <- list()
for (id in 1:N){
  simDat[[id]] <- simMove(nbObs = nbObs, xy0 = as.numeric(ACs[id,]), ACstrength = alpha, torus=T)
  print(paste('Animal',id,'done'))
}

# compile into dataframe
moveSims <- do.call(rbind, simDat)

# save
write.csv(moveSims,row.names=FALSE,file=paste('Data/MovementSims/simDat_', scenario, '_sol_', iter, '.csv', sep=''))

"-------- Update progress bar --------"
# count files for each scenario for solitary
sm <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_small_sol | wc -l", intern = T))
md <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_medium_sol | wc -l", intern = T))
lg <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_large_sol | wc -l", intern = T))

nsims <- 100
data <- rbind(c(sm,md,lg), 
              nsims-c(sm,md,lg))
colnames(data) <- c("sm","md","lg")

png("Data/MovementSims/SolMovementProgress.png", type="cairo")
plotx<- barplot(data, horiz = TRUE, names.arg = c("sm","md","lg"), main = "MovementSims (sol)")
text(x = c(sm,md,lg)+5, y = plotx, label = data[1,], pos = 1, cex = 0.8, col = "red")
dev.off()
