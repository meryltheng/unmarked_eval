# ====================================
# Unmarked abundance evalaution
# Movement simulations (group movement)
# adapted from https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12155 
# ====================================

library(CircStats)
library(data.table)
# source("R/0_functions.R")

'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for scenario, arg[2] for iter, arg[3] for group size)
# e.g., `Rscript --vanilla ./R/1c_MoveSims_group.R small 1 20`
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# test if there is at least one argument: if not, return an error
if (length(args) <= 1) {
  stop("At least three arguments must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input: ", args[1], args[2], args[3], sep = ' '))
  scenario = args[1]
  iter = as.numeric(args[2])
  groupSize = as.numeric(args[3])
}

if (scenario == 'small'){
  if (!groupSize %in% c(20,100)) stop("Group size only programmed for 20 & 100 for high density.")
  if(groupSize == 20){
    # AC distribution params
    minDist = 10; borderBuff = 20 # min dist between ACs; min dist between ACs and landscape border
    # directional distribution (von Mises) params
    kappac <- c(0.75) # group centroid concentration/attraction strength [also controls group range size]
    etac <- c(0.85) # centroid weight of CRW component in BCRW for centroid movement
    kappaind <- c(4, 3) # state-dependent concentration parameters for indiv movement [i.e., directedness towards group AC; directedness in CRW component]
    # state transition probabilites
    gamma <- matrix(c(0.95, 0.3, 0.05, 0.7), ncol = 2) # state transition probability matrix
  }
  if(groupSize == 100){
    # AC distribution params
    minDist = 150; borderBuff = 75 # min dist between ACs; min dist between ACs and landscape border
    # directional distribution (von Mises) params
    kappac <- c(1) # group centroid concentration/attraction strength
    etac <- c(0.95) # centroid weight of CRW component in BCRW for centroid movement
    kappaind <- c(3, 4) # state-dependent concentration parameters for indiv movement [i.e., directedness towards group AC; directedness in CRW component]
    # state transition probabilites
    gamma <- matrix(c(0.95, 0.2, 0.05, 0.8), ncol = 2) # state transition probability matrix
  }

  # step-length params
  ac <- c(2) # centroid mean
  bc <- c(1.35) # centriod sd
  aind <- c(2.5, 2.5) # individuals' state specific means
  bind <- c(1.72, 1.72) # individuals' state specific sds
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   
  # 48.26  499.93  628.96  659.67  785.96 2962.88
  N = 2000
}

if (scenario == 'medium'){
  if (!groupSize %in% c(5,25)) stop("Group size only programmed for 5 & 25 for medium density.")
  if(groupSize == 25){
    # AC distribution params
    minDist = 250; borderBuff = 75 # min dist between ACs; min dist between ACs and landscape border
    # directional distribution (von Mises) params
    kappac <- c(0.75) # group centroid concentration/attraction strength [also controls group range size]
    etac <- c(0.85) # centroid weight of CRW component in BCRW for centroid movement
    kappaind <- c(4, 3) # state-dependent concentration parameters for indiv movement [i.e., directedness towards group AC; directedness in CRW component]
    # state transition probabilites
    gamma <- matrix(c(0.95, 0.3, 0.05, 0.7), ncol = 2) # state transition probability matrix
  }
  if(groupSize == 5){
    # AC distribution params
    minDist = 100; borderBuff = 25
    # directional distribution (von Mises) params
    kappac <- c(2) # group centroid concentration/attraction strength [also controls group range size]
    etac <- c(0.85) # centroid weight of CRW component in BCRW for centroid movement
    kappaind <- c(4, 2) # state-dependent concentration parameters for indiv movement [i.e., directedness towards group AC; directedness in CRW component]
    # state transition probabilites
    gamma <- matrix(c(0.95, 0.4, 0.05, 0.6), ncol = 2) # state transition probability matrix
  }
  # step-length params
  ac <- c(3.2) # centroid mean
  bc <- c(2) # centriod sd
  aind <- c(5.5, 5.5) # individuals' state specific means
  bind <- c(3.2, 3.2) # individuals' state specific sds
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   
  # 48.26  499.93  628.96  659.67  785.96 2962.88
  N = 200
}

# # directional distribution (von Mises) params
# etac <- c(0.85) # centroid weight of CRW component in BCRW for centroid movement
# kappaind <- c(4, 3) # state-dependent concentration parameters for indiv movement [i.e., directedness towards group AC; directedness in CRW component]

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
simMoveGroup <- function(nbObs = nbObs, xy0 = c(0,0), groupSize = groupSize, groupID = 1, memberID = 1:groupSize, torus = T)
{
  # Generation of the state sequences for the different individuals
  Zind <- vector("list")
  for (k in 1:groupSize) {
    Zind[[k]] <- sample(1:nState, size = 1, prob = delta)
  }
  for (k in 2:nbObs) {
    for (j in 1:groupSize) {
      Zind[[j]][k] <- sample(1:nState, size = 1, prob = gamma[Zind[[j]][k - 1],])
    }
  }
  # Choice of the initial position of the centroid and of the location of the centroid’s centre of attraction (C)
  X <- matrix(rep(NA, 2 * nbObs), ncol = 2)
  X[1, ] <- xy0
  C <- matrix(xy0, byrow = T, ncol = 2)
  
  # Choice of the (random) initial positions of individuals 
  # Xind is the list that will comprise all individuals' simulated locations
  Xind <- vector("list")
  for (k in 1:groupSize) {
    Xind[[k]] <- matrix(rep(NA, 2 * nbObs), ncol = 2)
    Xind[[k]][1, ] <- c(rnorm(1,xy0[1], 1),rnorm(1,xy0[2], 1))
  }
  
  phi <- 0
  phiind <- rep(0, groupSize)
  
  for (k in 1:(nbObs - 1)) {
    ## start generation of centroid location at time k+1
    coo <- C[, ] - X[k, ] # xy diff btw centroid AC and current centroid loc
    mu <- Arg(coo[1] + (0+1i) * coo[2]) # direction towards centroid AC
    if (mu < 0)
      mu <- mu + 2 * pi # correct for direction
    mu.av <- Arg(etac * exp(phi * (0+1i)) + (1 - etac) * exp(mu * (0+1i))) # direction also but??
    phi <- rvm(1, mean = mu.av, k = kappac) # draw next direction
    if (k==1){
      step.len.ac <-  rgamma(1, shape = ac^2/bc^2, scale = bc^2/ac) # first step-length
    } else{
      step.len.ac <- 0.5*(step.len.ac + rgamma(1, shape = ac^2/bc^2, scale = bc^2/ac)) # subsequent step-lengths are correlated to previous
    }
    step.ac <- step.len.ac * c(Re(exp((0+1i) * phi)), Im(exp((0+1i) * phi))) # dist in xy to move 
    X[k + 1, ] <- X[k, ] + step.ac # next loc 
    ## start generation of individuals' locations at time k+1
    for (j in 1:groupSize) {
      if (s.ty[Zind[[j]][k]] == 1) { # BCRW
        coo <- X[k + 1, ] - Xind[[j]][k, ] # xy diff btw centroid loc and indiv loc
        mu <- Arg(coo[1] + (0+1i) * coo[2])
        if (mu < 0)
          mu <- mu + 2 * pi
        phiind[j] <- rvm(1, mean = mu, k = kappaind[Zind[[j]][k]])
      }
      if (s.ty[Zind[[j]][k]] == 0) { # CRW
        mu <- rvm(1, mean = 0, k = kappaind[Zind[[j]][k]])
        phiind[j] <- phiind[j] + mu
      }
      if (k==1){ # first step-length
        step.len <- rgamma(1, shape = aind[Zind[[j]][k]]^2/bind[Zind[[j]][k]]^2,
                            scale = bind[Zind[[j]][k]]^2/aind[Zind[[j]][k]])
      } else{
        step.len <- 0.5*( step.len + rgamma(1, shape = aind[Zind[[j]][k]]^2/bind[Zind[[j]][k]]^2,
                                            scale = bind[Zind[[j]][k]]^2/aind[Zind[[j]][k]]) )
      }
      step <- step.len * c(Re(exp((0+1i) * phiind[j])), Im(exp((0+1i) * phiind[j])))
      Xind[[j]][k + 1, ] <- Xind[[j]][k, ] + step
    }
  }
  # compile output
  xy <- as.data.table(do.call(rbind, Xind))
  colnames(xy) <- c('x','y')
  xy$ID <- rep(memberID, each=nbObs)
  xy$time <- rep(1:nbObs, times = groupSize)
  xy$groupID <- rep(groupID,nrow(xy))
  
  # correct coordinates for torus (data.table for fast math)
  if (torus == T){
    xy[x <= xlim[1], x := x + xlim[2]]
    xy[x >= xlim[2], x := x - xlim[2]]
    xy[y <= ylim[1], y := y + ylim[2]]
    xy[y >= ylim[2], y := y - ylim[2]]
  }
  return(xy)
}
"-------- Prep sim parameters --------"
# set general parameters
days = 100
nbObs <- 12*12*days # 10-min timesteps, 12-h activity, 50 days
xlim <- ylim <- c(0,1000) # landscape bounds

# state params
nState = 2
s.ty <- c(1, 0) # (1: BRW; 0: CRW)
#gamma <- matrix(c(0.95, 0.3, 0.05, 0.7), ncol = 2) # state transition probability matrix
delta <- solve(diag(nState) - t(gamma) + 1, rep(1, nState)) # corresponding stationary distribution of the Markov chain

# step-length distribution (gamma) params
# ac <- c(4) # centroid mean
# bc <- c(2.5) # centriod sd
# aind <- c(5.5, 5.5) # individuals' state specific means
# bind <- c(3.2, 3.2) # individuals' state specific sds

"-------- Run movement simulations --------"
nGroups = N/groupSize
# generate group centroid's ACs
ACs <- genAC(N = nGroups, minDist = minDist, borderBuff = borderBuff, xlim = xlim, ylim = ylim)

phi <- 0
phiind <- rep(0, groupSize)

simDat <- list()
for (id in 1:nGroups){
  memberID <- groupSize * (id-1) + 1:groupSize
  simDat[[id]] <- simMoveGroup(nbObs = nbObs, xy0 = as.numeric(ACs[id,]), 
                               groupSize = groupSize, groupID = id, 
                               memberID = memberID, torus = F)
  print(paste('Group',id,'done'))
}

# compile into dataframe
moveSims <- data.table(do.call(rbind, simDat))

# correct coordinates for torus (data.table for fast math)
moveSims[x <= xlim[1], x := x + xlim[2]]
moveSims[x >= xlim[2], x := x - xlim[2]]
moveSims[y <= ylim[1], y := y + ylim[2]]
moveSims[y >= ylim[2], y := y - ylim[2]]

# save
write.csv(moveSims,row.names=FALSE,file=paste('Data/MovementSims/simDat_', scenario, '_grp_', iter, '_gSize', groupSize,'.csv', sep=''))

"-------- Update progress bar --------"
# count files for each scenario for group
sm20 <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_small_grp_.*_gSize20 | wc -l", intern = T))
sm100 <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_small_grp_.*_gSize100 | wc -l", intern = T))
med5 <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_medium_grp_.*_gSize5 | wc -l", intern = T))
med25 <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_medium_grp_.*_gSize25 | wc -l", intern = T))

nsims <- 100
data <- rbind(c(sm20,sm100,med5,med25), 
               nsims-c(sm20,sm100,med5,med25))
colnames(data) <- c("sm20","sm100","med5","med25")

png("Data/MovementSims/GroupMovementProgress.png", type="cairo")
plotx<- barplot(data, horiz = TRUE, names.arg = c("sm20","sm100","med5","med25"), main = "MovementSims (grp)")
text(x = c(sm20,sm100,med5,med25)+5, y = plotx, label = data[1,], pos = 1, cex = 0.8, col = "red")
dev.off()

# "-------- Check --------"
# library(dplyr)
# library(scales)
# simDat <- as.data.frame(do.call(rbind, Xind))
# simDat$ID <- rep(1:groupSize, each=nbObs)
# colnames(simDat)[1:2] <- c('x','y')
# calculate steplength
# moveSims <- moveSims %>%
#   group_by(ID) %>%
#   mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
#   as.data.frame()
# 
# summary(moveSims$step) * 10 * 12
# 
# plot(1, type='n', xlim = xlim, ylim = ylim)
# plotTraject(df = moveSims, plotit = rep(1,N), torus = T) # all
# 
# plot(1, type='n', xlim = c(300,700), ylim = c(300,700))
# plotTraject(df = simDat[simDat$ID==1,], plotit = rep(1,1), torus = F)
# 
# # "-------- Visualise --------"
# setwd("/Users/a1784952/Desktop/PhD/Manuscripts/Paper 3/UnmarkedAbundance/Plots/gif_dump")
# skips <- c(6,seq(10, nbObs/2, by=5))
# for (q in 6:500){
#   id <- sprintf("%03d", q)
#   png(paste("step",id,".png", sep=""), width=720, height=700, units="px",
#       pointsize=18)
#   par(mar = c(2, 2, 2, 2))
#   plot(1, type='n', xlim = xlim, ylim = ylim)
# 
#   for(i in 1:length(Xind)){ # individuals
#     lines(Xind[[i]][(q-5):q,1],Xind[[i]][(q-5):q,2], col='grey')
#     points(Xind[[i]][q,1],Xind[[i]][q,2], pch= 19, col='grey')
#   } # end of i loop
# 
#   points(X[q,1],X[q,2], pch= 19, col='red',cex=2)
# 
#   dev.off()
# }
# 
# system("convert -delay 10 *.png -coalesce -fuzz 4% +dither -layers Optimize group12.gif")
# file.remove(list.files(pattern=".png"))
# 
# # group
# library(viridisLite)
# xlim=ylim=c(350,650)
# 
# for (q in 6:500){
#   timestep <- sprintf("%03d", q)
#   png(paste("step",timestep,".png", sep=""), width=720, height=700, units="px",
#       pointsize=18)
#   par(mar = c(2, 2, 2, 2))
#   plot(1, type='n', xlim = xlim, ylim = ylim)
# 
#   for (id in 1:nGroups){ #group
#     for(i in unique(simDat[[id]]$ID)){ # group members
#       lines(simDat[[id]][simDat[[id]]$ID==i,][(q-5):q,1],
#             simDat[[id]][simDat[[id]]$ID==i,][(q-5):q,2],
#             col=viridis(nGroups)[id])
#       points(simDat[[id]][simDat[[id]]$ID==i,][q,1],
#              simDat[[id]][simDat[[id]]$ID==i,][q,2],
#              pch= 19, col=viridis(nGroups)[id])
#     } # end of i loop
#   }# end of id loop
#   plot(traps[[1]],add=T)
#   dev.off()
# }
