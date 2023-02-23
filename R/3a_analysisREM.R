# ====================================
# Unmarked abundance evalaution
# REM model
# ====================================

library(dplyr)
library(tidyr)
library(rgeos)
library(activity)
library(Distance)
source("R/camtools.R") # https://github.com/MarcusRowcliffe/camtools/blob/master/camtools.R
source("R/sbd.R") # https://github.com/MarcusRowcliffe/sbd/blob/master/sbd.r

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

'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for scenario, arg[2] for iter, arg[3] for nTraps, arg[4] for movemt behav (sol/grp), if grp arg[5] is group size)
# e.g., `Rscript --vanilla ./R/3a_analysisREM.R small 1 100 sol`
# e.g., `Rscript --vanilla ./R/3a_analysisREM.R small 1 100 grp 25`
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

if(subsetDat == TRUE){
  x4 <- 
    x4 %>%
    filter(time %in% 1:(12 * 12 * days)) 
}

"-------- Prep model input --------"
# consolidate independent tracks (track) and calc input params (for REM)
x5 <- x4 %>% #i
  group_by(Site, ID) %>% #ii
  arrange(time) %>% #iii
  mutate( difftime = c(NA, diff(time)) ) %>% #iv
  group_by(track) %>% # v
  # if the first row of track was flagged as consecutive track (difftime==1) AND if start of track falls within FOV (entryTime<1s), 
  # assign it to same time as prev consec track
  mutate( time_new = if_else( first(difftime) == 1 & (!is.na(first(difftime))) & first(entryTime)<=1, 
                              first(time)-1, as.numeric(first(time)), NA_real_) ) %>% # vi
  group_by(Site, ID, time_new) %>% # vii
  mutate( xy_dist = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 )) %>% # viii
  summarise(ID = unique(ID), time_new = unique(time_new), Site = unique(Site), entryTime = first(entryTime),
            v = sum(xy_dist,na.rm=T) / (last(time_seconds)-first(time_seconds)), 
            r = first(dist), 
            time_seconds = first(time_seconds) # ix
  ) %>%
  filter( v < 1000 | is.nan(v) ) %>% # keep single record sequences for encounter rate calculation
  arrange(Site, time_seconds) %>% # xi
  as.data.frame()

# number of animals detected within the same timestep at the same Site
group_counts <- x4 %>%
  group_by(Site, ID, time) %>%
  summarise(ID = unique(ID), time = unique(time), Site = unique(Site)) %>%
  group_by(Site, time) %>%
  summarise(g = length(ID), time = first(time), Site = unique(Site)) %>%
  as.data.frame()

# speed estimation (described here: https://zslpublications.onlinelibrary.wiley.com/doi/10.1002/rse2.17)
speed_dat <- x5[!is.nan(x5$v) ,]
# Fits all three size biased options (i.e., lnorm, gamma, weibull); exclude speeds that are less than 0.01 m/s
speed_mods <- sbm3(v ~1, speed_dat[speed_dat$v > 0.01,], reps = 1000) 
best_speed_mod <- as.name(row.names(speed_mods[["AICtab"]])[which.min(speed_mods[["AICtab"]]$AIC)])
pred_speed <- predict.sbm(speed_mods[["models"]][[best_speed_mod]])

## detection radius estimation
# detection radius estimation
# prep data
rad_data <- x5
# truncation distances
if(scenario == "hi" & behav == "grp" & groupSize == 100) {dist_range = c(0,10)} else{dist_range = c(0,10)} # adjust for peak in furthest r
mybreaks <- seq(dist_range[1],dist_range[2],1) 
trunc.list <- list(left=dist_range[1], right=dist_range[2])

# fit default hr model
dist_model = 'hr'
best_mod <- ds(rad_data$r, transect = "point", key= dist_model, adjustment = NULL, order = 0, cutpoints = mybreaks, truncation = trunc.list)

#calc EDR
effec_radius <- EDRtransform(best_mod)
# test fit
gof_fit <- gof_ds(best_mod)

# # save radius fits for reference
# # probability density
# if (behav == 'sol'){
#   png(paste0("Data/Estimates/Plots/REM-radius-",scenario, '_sol_', iter, '_J', nTraps,".png"), type="cairo", 
#       width = 900, height = 450, units = "px", pointsize = 20)
# }
# if (behav == 'grp'){
#   png(paste0("Data/Estimates/Plots/REM-radius-",scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,".png"), type="cairo",
#       width = 900, height = 450, units = "px", pointsize = 20)
# }
# par(mfrow=c(1,2))
# plot(ds_fit, main='', xlab="",
#      showpoints=FALSE, lwd=3, xlim=c(0, 1),pdf=T)
# plot(ds_fit, main='', xlab="",
#      showpoints=FALSE, lwd=3, xlim=c(0, 1),pdf=F)
# 
# mtext(paste(dist_model, scenario, iter, behav, groupSize, sep = '_'),side=3,line=-2,outer=TRUE,cex=1.3)
# mtext("Radius (unit)",side=1,line=-2,outer=TRUE,cex=1.3)
# 
# mtext(paste('EDR =', round(effec_radius$EDR,2) ), col = 'red', side=3,line=-3.5,outer=TRUE,cex=1)
# dev.off()

# REM (camtools)
theta <- pi/180 * 40 # assuming perfect detection with declining angular distance

# prepare encounter data by site
input_y <- x5 %>%
  filter(r >= dist_range[1] & r <= dist_range[2])
y1 <- aggregate(input_y$Site, by = list(Site=input_y$Site), FUN = length) # discard detections outside of truncation distances
y2 <- data.frame(Site=1:nTraps, x = rep(0,nTraps))
y3 <- y1 %>% tidyr::complete(Site = 1:nTraps, fill = list(x=0) ) %>% as.data.frame()
y <- as.vector(y3[,2]) # ENCOUNTER DATA

# prep params
param <- list(v = pred_speed$est * 9 * 60 * 60, # in seconds, based on 12-h activity day
              p = 1,
              #p = actmod@act["act"], # activity level within the 12-h day
              r = effec_radius$EDR,
              theta = theta )
paramse <- list(v = pred_speed$se * 9 * 60 * 60,
                p = 0,
                #p = actmod@act["se"],
                r = effec_radius$se.EDR,
                theta = 0)

"-------- Execute REM --------"
D_hat <- bootTRD(y, rep(days, nTraps), param, paramse) # encounters, sampling duration per cam, etc.
A <- 1000*1000  # 1 km^2
( N_hat <- D_hat * A )# density estimate /km^2

# Calculate overdispersion factor (chat) 
y = nrow(input_y)
expected.y = (D_hat[1] * days * (pred_speed$est * 9 * 60 * 60) * effec_radius$EDR * (2 + theta)) / pi
X2 <- sum((y - expected.y)^2 / expected.y)
si <- mean((y - expected.y) / expected.y)
chat <- X2/(nTraps-1) / (1 + si) # Fletcher 2012

"-------- Save/write output --------" 
Replicate = iter
Scenario = scenario
Behaviour = behav
if (behav == 'sol'){ groupSize = NA }
Model = "REM"
Trap.effort = nTraps
Days.monitored = days
nDet = nrow(x4)
Dhat = N_hat[1]
Dhat.se = N_hat[2]
Dhat.lo = NA
Dhat.hi = NA
y = y
y.mean = y/nTraps
y.sd = sd(as.vector(y3[,2]))
g.mean = mean(group_counts$g)
g.sd = sd(group_counts$g)
g.max = max(group_counts$g)
v = pred_speed$est * 9 * 60 * 60
v.se =  pred_speed$se * 9 * 60 * 60
act = 1 # activity level within the 12-h day; actmod@act["act"]
act.se = 0 # actmod@act["se"]
edr.model = "hr0" #as.character(row.names(ds_AIC)[which.min(ds_AIC$AIC)]) 
edr = effec_radius$EDR
edr.se = effec_radius$se.EDR
edr.AIC = as.numeric(AIC(best_mod)['AIC'])
chat =  chat
gof_chisq = gof_fit[["chisquare"]][["chi1"]][["chisq"]]
gof_df = gof_fit[["chisquare"]][["chi1"]][["df"]]
gof_rel = gof_chisq/y

rem_output <- data.frame(Replicate, Scenario, Behaviour, groupSize, Model, Trap.effort, Days.monitored, nDet,
                         Dhat, Dhat.se, Dhat.lo, Dhat.hi, y, y.mean, y.sd, g.mean, g.sd, g.max, 
                         v, v.se, act, act.se, edr.model, edr, edr.se, edr.AIC, chat,
                         gof_chisq, gof_df, gof_rel)

# append to existing master file
write.table( rem_output,  
             file="./Data/Estimates/rem_imperfectDet.csv", 
             append = T, 
             sep=',', 
             row.names=F, 
             col.names=F )

# # create empty .csv if not avail
# columns = c("Replicate", "Scenario", "Behaviour", "groupSize", "Model", "Trap.effort", "Days.monitored", "nDet",
#            "Dhat", "Dhat.se", "Dhat.lo", "Dhat.hi", "y",  "y.mean", "y.sd", "g.mean", "g.sd", "g.max", 
#            "v", "v.se", "act", "act.se", "edr.model", "edr", "edr.se", "edr.AIC", "chat",
#            "gof_chisq", "gof_df", "gof_rel")
# mydf = data.frame(matrix(nrow = 0, ncol = length(columns)))
# colnames(mydf) = columns
# write.csv(mydf, file=paste('Data/Estimates/rem_imperfectDet.csv',sep=''), row.names = F)