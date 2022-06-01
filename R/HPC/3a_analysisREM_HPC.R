# ====================================
# Unmarked abundance evalaution
# REM model (for HPC)
# ====================================

library(dplyr)
library(tidyr)
library(rgeos)
library(activity)
library(Distance)
#source("R/0_functions.R") 
source("R/camtools.R") # https://github.com/MarcusRowcliffe/camtools/blob/master/camtools.R
source("R/sbd.R") # https://github.com/MarcusRowcliffe/sbd/blob/master/sbd.r
#source("R/distancedf.r") # https://github.com/MarcusRowcliffe/distanceDF/blob/master/distancedf.r

'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for scenario, arg[2] for iter, arg[3] for nTraps, arg[4] for movemt behav (sol/grp), if grp arg[5] is group size)
# e.g., `Rscript --vanilla ./R/3c_analysisREM.R small 1 100 sol`
# e.g., `Rscript --vanilla ./R/3c_analysisREM.R small 1 100 grp 25`
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
# 6 consolidate independent detections (tracks) and calc input params (for REM)
x5 <- x4 %>% 
  group_by(time_seconds) %>%
  mutate( unique_cap = if_else(min(dist)==dist, 1, 0) ) %>% # if multiple animals captured, the nearest animal's dist will be marked for radius calculation
  group_by(track) %>% # following REM inputs calculated per track/sequence
  mutate( xy_dist = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 )) %>% # dist btw each 'photo' within an independent detection
  summarise(ID = unique(ID), time = unique(time), track = unique(track), Site = unique(Site), 
            v = sum(xy_dist,na.rm=T) / (last(time_seconds)-first(time_seconds)), # speed in unit/s
            r = first(dist), # r is just the dist-to-CT for first trigger in the sequence
            time_seconds = first(time_seconds),
            unique_cap = first(unique_cap)
  ) %>%
  filter( !is.na(v), !is.infinite(v) ) %>% 
  group_by(track) %>%
  filter(row_number()==1) %>% # only keep first trigger in a sequence (i.e., det_id)
  group_by(Site, time_seconds) %>% 
  mutate(r_new = ifelse(unique_cap==1, r, NA) ) %>% # only keep radius measure for nearest indiv in group captures
  as.data.frame()

# speed estimation (described here: https://zslpublications.onlinelibrary.wiley.com/doi/10.1002/rse2.17)
speed_mods <- sbm3(v ~1, x5[x5$v > 0.001,]*1000, reps = 1000) # Fits all three size biased options (i.e., lnorm, gamma, weibull); exclude speeds that are less than 0.001 unit/s (i.e., 1cm/s)
best_speed_mod <- as.name(row.names(speed_mods[["AICtab"]])[which.min(speed_mods[["AICtab"]]$AIC)])
pred_speed <- predict.sbm(speed_mods[["models"]][[best_speed_mod]])/1000

# # save v plot for reference
# if (behav == 'sol'){
#   png(paste0("Data/Estimates/Plots/REM-speed-",scenario, '_sol_', iter, '_J', nTraps,".png"), type="cairo")
# }
# if (behav == 'grp'){
#   png(paste0("Data/Estimates/Plots/REM-speed-",scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,".png"), type="cairo")
# }
# pred_speed; plot.sbm(speed_mods[["models"]][[best_speed_mod]], log = T, xlab="Speed (log[unit/s])" )
# dev.off()

# # activity level estimation
# # since we assumed 100% activity within 12h, we treat activity (p) as 100%
# x5$day <- ceiling(x5$time/144)
# x5$hour <- ceiling((x5$time - (x5$day-1)*144)/12)
# actmod <- fitact(x5$hour * pi/12,
#                 bounds=c(1,12)*pi/12, sample="data", reps=100)
# plot(actmod, centre="night", dline=list(col="grey"))

## detection radius estimation
#truncation distances
mybreaks <- seq(0.1,1,0.1)
trunc.list <- list(left=0.1, right=1)
##### ISSUES WITH ADJUSTMENTS FOR ORDER >1 - only on HPC ####
# # uniform
# uni1 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="uni", adjustment = 'cos', order = 1, cutpoints = mybreaks, truncation = trunc.list)
# uni2 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="uni", adjustment = 'cos', order = c(1,2), cutpoints = mybreaks, truncation = trunc.list)
# uni3 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="uni", adjustment = 'cos', order = c(1,2,3), cutpoints = mybreaks, truncation = trunc.list)
# half-normal
hn0 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hn", adjustment = NULL, order = 0, cutpoints = mybreaks, truncation = trunc.list)
# hn_cos2 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hn", adjustment = "cos", order = 2, cutpoints = mybreaks, truncation = trunc.list)
# hn_herm0 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hn", adjustment = "herm", order = 0, cutpoints = mybreaks, truncation = trunc.list)
# hn_herm2 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hn", adjustment = "herm", order = 4, cutpoints = mybreaks, truncation = trunc.list)
# hn_poly0 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hn", adjustment = "poly", order = 0, cutpoints = mybreaks, truncation = trunc.list)
# hn_poly2 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hn", adjustment = "poly", order = c(2,4), cutpoints = mybreaks, truncation = trunc.list)
# hazard-rate
# hr0 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hr", adjustment = NULL, order = 0, cutpoints = mybreaks, truncation = trunc.list)
# hr_cos0 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hr", adjustment = "cos", order = 0, cutpoints = mybreaks, truncation = trunc.list)
# hr_cos2 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hr", adjustment = "cos", order = 2, cutpoints = mybreaks, truncation = trunc.list)
# hr_herm0 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hr", adjustment = "herm", order = 0, cutpoints = mybreaks, truncation = trunc.list)
# hr_herm2 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hr", adjustment = "herm", order = 4, cutpoints = mybreaks, truncation = trunc.list)
# hr_poly0 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hr", adjustment = "poly", order = 0, cutpoints = mybreaks, truncation = trunc.list)
# hr_poly2 <- ds(x5$r_new[!is.na(x5$r_new)], transect = "point", key="hr", adjustment = "poly", order = c(2,4), cutpoints = mybreaks, truncation = trunc.list)
best_mod <- hn0

# #best model
# AIC_res <- AIC(uni1,uni2,uni3,hn0, hn_cos2, hn_herm0, hn_herm2, hn_poly0, hn_poly2, hr0, hr_cos0, hr_cos2, hr_herm0, hr_herm2, hr_poly0, hr_poly2)
# best_mod <- eval(as.name(row.names(AIC_res)[which.min(AIC_res$AIC)]))
# # plot(best_mod, main="Peak activity", xlab="Radius (unit)",
# #      showpoints=FALSE, lwd=3, xlim=c(0, 1),pdf=T)

# # save radius fits for reference
# # probability density
# if (behav == 'sol'){
#   png(paste0("Data/Estimates/Plots/REM-radius-pd",scenario, '_sol_', iter, '_J', nTraps,".png"), type="cairo", width = 1200, height = 1000, units = "px", pointsize = 20)
# }
# if (behav == 'grp'){
#   png(paste0("Data/Estimates/Plots/REM-radius-pd",scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,".png"), type="cairo",
#       width = 1200, height = 1000, units = "px", pointsize = 20)
# }
# par(mfrow=c(4,4))
# par(oma=c(3,3,1,1)) 
# par(mar=c(2,2,4,2))
# all_mods <- list(uni1,uni2,uni3,hn0, hn_cos2, hn_herm0, hn_herm2, hn_poly0, hn_poly2, hr0, hr_cos0, hr_cos2, hr_herm0, hr_herm2, hr_poly0, hr_poly2)
# for (m in 1:16){
#   plot(all_mods[[m]], main=row.names(AIC_res)[m], xlab="Radius (unit)",
#        showpoints=FALSE, lwd=3, xlim=c(0, 1),pdf=T)
# }
# mtext("Radius (unit)",side=1,line=2,outer=TRUE,cex=1.3)
# mtext("Probability density",side=2,line=1,outer=TRUE,cex=1.3,las=0)
# dev.off()
# # detection probability
# if (behav == 'sol'){
#   png(paste0("Data/Estimates/Plots/REM-radius-dp",scenario, '_sol_', iter, '_J', nTraps,".png"), type="cairo", width = 1200, height = 1000, units = "px", pointsize = 20)
# }
# if (behav == 'grp'){
#   png(paste0("Data/Estimates/Plots/REM-radius-dp",scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,".png"), type="cairo",
#       width = 1200, height = 1000, units = "px", pointsize = 20)
# }
# par(mfrow=c(4,4))
# par(oma=c(3,3,1,1)) 
# par(mar=c(2,2,4,2))
# all_mods <- list(uni1,uni2,uni3,hn0, hn_cos2, hn_herm0, hn_herm2, hn_poly0, hn_poly2, hr0, hr_cos0, hr_cos2, hr_herm0, hr_herm2, hr_poly0, hr_poly2)
# for (m in 1:16){
#   plot(all_mods[[m]], main=row.names(AIC_res)[m], xlab="Radius (unit)",
#        showpoints=FALSE, lwd=3, xlim=c(0, 1),pdf=F)
# }
# mtext("Radius (unit)",side=1,line=2,outer=TRUE,cex=1.3)
# mtext("Detection probability",side=2,line=1,outer=TRUE,cex=1.3,las=0)
# dev.off()

# function to transform dsobject to EffectiveDetectionRadius as well as produce the SE and CI bounds (https://github.com/DistanceDevelopment/mrds/issues/36#issuecomment-753462999)
EDRtransform <- function(dsobject, alpha=0.05) { 
  if(class(dsobject) != "dsmodel") stop("First argument must be a dsmodel object")
  if(!dsobject$ddf$meta.data$point) stop("EDR can only be computed for point transect data")
  summary.ds.model <- summary(dsobject)
  p_a <- summary.ds.model$ds$average.p
  se.p_a <- summary.ds.model$ds$average.p.se
  cv.p_a <- se.p_a / p_a
  w <- summary.ds.model$ds$width
  edr <- sqrt(p_a * w^2)
  se.edr <- cv.p_a/2 * edr
  degfree <- summary.ds.model$ds$n - length(summary.ds.model$ddf$par)
  t.crit <- qt(1 - alpha/2, degfree)
  se.log.edr <- sqrt(log(1 + (cv.p_a/2)^2))
  c.mult <- exp(t.crit * se.log.edr)
  ci.edr <- c(edr / c.mult, edr * c.mult)
  return(list(EDR=edr, se.EDR=se.edr, ci.EDR=ci.edr))
}
effec_radius <- EDRtransform(best_mod)
effec_radius$EDR #mean
effec_radius$se.EDR #se

# REM (camtools)
Ncam <- length(traps[[1]])   # Number of camera traps
A <- 1000*1000  # Study area in unit^2
theta <- pi/180 * 60

# prepare encounter data by site
y0 <- x5[x5$r>=0.1, ] # x5[x5$r_new>=0.1 !is.na(x5$r_new), ]
y1 <- aggregate(y0$Site, by = list(Site=y0$Site), FUN = length) # discard detections outside of truncation distances
y2 <- data.frame(Site=1:Ncam, x = rep(0,Ncam))
y3 <- y1 %>% tidyr::complete(Site = 1:Ncam, fill = list(x=0) ) %>% as.data.frame()
y <- as.vector(y3[,2]) # ENCOUNTER DATA

# prep params
param <- list(v = pred_speed$est * 12 * 60 * 60, # in seconds, based on 12-h activity day
              p = 1,
              #p = actmod@act["act"], # activity level within the 12-h day
              r = effec_radius$EDR,
              theta = theta )
paramse <- list(v = pred_speed$se * 12 * 60 * 60,
                p = 0,
                #p = actmod@act["se"],
                r = effec_radius$se.EDR,
                theta = 0)

"-------- Execute REM --------"
D_hat <- bootTRD(y, rep(days, Ncam), param, paramse) # encounters, sampling duration per cam, etc.
( N_hat <- D_hat * A )# abundance estimate

"-------- Save/write output --------" 
Replicate = iter
Scenario = scenario
Behaviour = behav
if (behav == 'sol'){ groupSize = NA }
Model = "REM"
Trap.effort = nTraps
Days.monitored = days
nDet = nrow(x4)
Nhat = N_hat[1]
Nhat.se = N_hat[2]
Nhat.lo = NA
Nhat.hi = NA
v = pred_speed$est * 12 * 60 * 60
v.se =  pred_speed$se * 12 * 60 * 60
act = 1 # activity level within the 12-h day; actmod@act["act"]
act.se = 0 # actmod@act["se"]
edr.model = "hn0" # as.character(row.names(AIC_res)[which.min(AIC_res$AIC)])
edr = effec_radius$EDR
edr.se = effec_radius$se.EDR

rem_output <- data.frame(Replicate, Scenario, Behaviour, groupSize, Model, Trap.effort, Days.monitored, nDet,
                         Nhat, Nhat.se, Nhat.lo, Nhat.hi, v, v.se, act, act.se, edr.model, edr, edr.se)

# append to existing master file
write.table( rem_output,  
             file="./Data/Estimates/rem_master_hn.csv", 
             append = T, 
             sep=',', 
             row.names=F, 
             col.names=F )

# # create empty .csv if not avail
# columns = c("Replicate", "Scenario", "Behaviour", "groupSize", "Model", "Trap.effort", "Days.monitored", "nDet",
#            "Nhat", "Nhat.se", "Nhat.lo", "Nhat.hi", "v", "v.se", "act", "act.se", "edr.model", "edr", "edr.se")
# mydf = data.frame(matrix(nrow = 0, ncol = length(columns)))
# colnames(mydf) = columns
# write.csv(mydf, file=paste('Data/Estimates/rem_master_test.csv',sep=''), row.names = F)