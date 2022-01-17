# ====================================
# Unmarked abundance evalaution
# CT-DS model
# ====================================

library(dplyr)
library(tidyr)
library(rgeos)
library(Distance)
source("R/0_functions.R") 

'------ Run via shell script ------'
# get the input passed from the shell script
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# test if there is at least one argument: if not, return an error
if (length(args) <= 1) {
  #stop("At least two arguments must be supplied (input file).\n", call. = FALSE)
  scenario = 'large'; iter = 1
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

# calculate distances
nbObs = nrow(moveSims[moveSims$ID == 1,])
x2 <- calcDist(movedf = moveSims, detdf = x1, X = traps.buff, timestep = 15, nbObs = nbObs, timeInt = 2)

hist(x2$distance,xlab='distance (m)',main='')
x2$distance <- x2$distance*10
bin_dist <- round(x2$distance)
ds_data <- data.frame(Region.Label = rep('SimLand', nrow(x2)), Area = rep(100, nrow(x2)),
                      multiplier = rep(samfrac, nrow(x2)), Sample.Label = x2$Site, Effort = rep(50*12*60*30, nrow(x2)),
                      distance = bin_dist)

"-------- Specify model --------" 
# fit detection function
mybreaks <- seq(0,12,1)
#conversion <- convert_units("meter", NULL, "square kilometer")
trunc.list <- list(left=1, right=10)
mybreaks <- c(seq(1,10,1))
hn0 <- ds(x2, transect = "point", key="hn", adjustment = NULL,
          cutpoints = mybreaks, truncation = trunc.list)
uni1 <- ds(x2, transect = "point", key="unif", adjustment = "cos",
           order=1,
           cutpoints = mybreaks, truncation = trunc.list)
uni2 <- ds(x2, transect = "point", key="unif", adjustment = "cos",
           order=c(1,2),
           cutpoints = mybreaks, truncation = trunc.list)
plot(uni1, main="Peak activity", xlab="Distance (units)",
     showpoints=FALSE, lwd=3, xlim=c(0, 11))

# estimate density
viewangle <- 60 # degrees
samfrac <- viewangle / 360
conversion <- convert_units("metre", NULL, "square kilometre")

peak.hr.dens <- dht2(hn0, flatfile=ds_data[ds_data$distance>=1,], strat_formula = ~1,
                     sample_fraction = samfrac, er_est = "P2", convert_units = conversion)
print(peak.hr.dens, report="abundance")
