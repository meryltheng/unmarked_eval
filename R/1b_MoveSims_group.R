# ====================================
# Unmarked abundance evalaution
# Movement simulations (group movement)
# adapted from https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12155 
# ====================================

library(CircStats)
library(data.table)
source("0_functions.R")

'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for scenario, arg[2] for iter, arg[3] for group size)
# e.g., `Rscript --vanilla ./R/1c_MoveSims_group.R hi 1 20`
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

if (scenario == 'hi'){
  if (!groupSize %in% c(20,100)) stop("Group size only programmed for 20 & 100 for high density.")
  if(groupSize == 20){
    # AC distribution params
    minDist = 100; borderBuff = 200 # min dist between ACs; min dist between ACs and landscape border
    # directional distribution (von Mises) params
    kappac <- c(6) # group centroid concentration param 
    etac <- c(0.85) # group centroid CRW weight
    kappaind <- c(6, 3) # state-dependent concentration parameters for indiv movement [i.e., directedness towards group AC; directedness in CRW component]
    # state transition probabilites
    gamma <- matrix(c(0.95, 0.3, 0.05, 0.7), ncol = 2) # state transition probability matrix
  }
  if(groupSize == 100){
    # AC distribution params
    minDist = 1500; borderBuff = 750 # min dist between ACs; min dist between ACs and landscape border
    # directional distribution (von Mises) params
    kappac <- c(4) # group centroid concentration param 
    etac <- c(1) # group centroid CRW weight
    kappaind <- c(4,3) # state-dependent concentration parameters for indiv movement [i.e., directedness towards group AC; directedness in CRW component]
    # state transition probabilites
    gamma <- matrix(c(0.95, 0.3, 0.05, 0.7), ncol = 2) # state transition probability matrix
  }
  
  # step-length params
  ac <- c(35) # centroid mean
  bc <- c(25) # centriod sd
  aind <- c(48, 48) # individuals' state specific means
  bind <- c(30, 30) # individuals' state specific sds
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   
  # 48.26  499.93  628.96  659.67  785.96 2962.88
  N = 2000
}

if (scenario == 'med'){
  if (!groupSize %in% c(5,25)) stop("Group size only programmed for 5 & 25 for medium density.")
  if(groupSize == 5){
    # AC distribution params
    minDist = 1000; borderBuff = 250
    # directional distribution (von Mises) params
    kappac <- c(2) 
    etac <- c(0.7) 
    kappaind <- c(8, 3) # state-dependent concentration parameters for indiv movement [i.e., directedness towards group AC; directedness in CRW component]
    # state transition probabilites
    gamma <- matrix(c(0.99, 0.3, 0.01, 0.7), ncol = 2) # state transition probability matrix
  }
  if(groupSize == 25){
    # AC distribution params
    minDist = 2500; borderBuff = 750 # min dist between ACs; min dist between ACs and landscape border
    # directional distribution (von Mises) params
    kappac <- c(1) 
    etac <- c(0.75) 
    kappaind <- c(6, 3) # state-dependent concentration parameters for indiv movement [i.e., directedness towards group AC; directedness in CRW component]
    # state transition probabilites
    gamma <- matrix(c(0.95, 0.3, 0.05, 0.7), ncol = 2) # state transition probability matrix
  }
  # step-length params
  ac <- c(80) # centroid mean
  bc <- c(50) # centriod sd
  aind <- c(100, 100) # individuals' state specific means
  bind <- c(62.5, 62.5) # individuals' state specific sds
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   
  # 48.26  499.93  628.96  659.67  785.96 2962.88
  N = 200
}

"-------- Prep sim parameters --------"
# set general parameters
days = 100
nbObs <- 12*9*days # 1-min timesteps, 9-h activity, 50 days
xlim <- ylim <- c(0,10000) # landscape bounds

"-------- Run movement simulations --------"
# generate group centroid's ACs
nGroups = N/groupSize
ACs <- genAC(N = nGroups, minDist = minDist, borderBuff = borderBuff, xlim = xlim, ylim = ylim)

simDat <- list()
for (id in 1:nGroups){
  memberID <- groupSize * (id-1) + 1:groupSize
  simDat[[id]] <- simMoveGroup(nbObs = nbObs, burnin = 12*9*1,
                               groupSize = groupSize, groupID = id, memberID = memberID,
                               # group centroid
                               xy0 = as.numeric(ACs[id,]), # ac loc
                               ac = ac, bc = bc, # step-length params
                               kappac = kappac, etac = etac, # concentration; CRW:BRW weights
                               # individuals
                               aind = aind, bind = bind, # state-specific step-length params
                               kappaind = kappaind, # state-specific concentration params
                               # state params
                               gamma = gamma, 
                               xlim = xlim, ylim = xlim, torus = T)
  print(paste('Group',id,'done'))
}

# compile into dataframe
moveSims <- data.table(do.call(rbind, simDat))

# save
write.csv(moveSims,row.names=FALSE,file=paste('Data/MovementSims/simDat_', scenario, '_grp_', iter, '_gSize', groupSize,'.csv', sep=''))

"-------- Update progress bar --------"
# count files for each scenario for group
hi20 <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_hi_grp_.*_gSize20 | wc -l", intern = T))
hi100 <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_hi_grp_.*_gSize100 | wc -l", intern = T))
med5 <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_med_grp_.*_gSize5 | wc -l", intern = T))
med25 <- as.numeric(system("ls -1 ./Data/MovementSims/| grep ^simDat_med_grp_.*_gSize25 | wc -l", intern = T))

nsims <- 100
data <- rbind(c(hi20,hi100,med5,med25), 
              nsims-c(hi20,hi100,med5,med25))
colnames(data) <- c("hi20","hi100","med5","med25")

png("Data/MovementSims/GroupMovementProgress.png", type="cairo")
plotx<- barplot(data, horiz = TRUE, names.arg = c("hi20","hi100","med5","med25"), main = "MovementSims (grp)")
text(x = c(hi20,hi100,med5,med25)+5, y = plotx, label = data[1,], pos = 1, cex = 0.8, col = "red")
dev.off()

# "-------- Check --------"
# library(dplyr)
# library(scales)
# simDat <- as.data.frame(do.call(rbind, Xind))
# simDat$ID <- rep(1:groupSize, each=nbObs)
# colnames(simDat)[1:2] <- c('x','y')
# # calculate steplength
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
