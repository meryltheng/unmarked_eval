# ====================================
# Unmarked abundance evalaution
# CT-DS model
# ====================================

library(dplyr)
library(tidyr)
library(Distance)

'------ Additional functions ------'
# ds model selection adjustments from overdispersion
chat <- function(modobj) {
  #  computes c-hat for a dsmodel object using Method 1 of Howe et al. (2018)
  test <- gof_ds(modobj)
  num <- test$chisquare$chi1$chisq
  denom <- test$chisquare$chi1$df
  chat <- num/denom
  return(chat)
}

qaic <- function(modobj, chat) {
  #  computes QAIC for a dsmodel object given a c-hat
  value <- 2* modobj$ddf$ds$value/chat + 2 * (length(modobj$ddf$ds$pars)+1)
  return(value)
}

qaic.pass1 <- function(...) {
  #   Performs Pass 1 model selection based upon Method 1 of Howe et al. (2018)
  #   Arguments are dsmodel objects; assumed all based on same key function
  #    c-hat is computed for the most parameter-rich model in the group
  #    qaic is calculated for each model in group based upon this c-hat
  #   Result returned in the form of a data.frame with model name, npar, aic and qaic
  models <- list(...)
  num.models <- length(models)
  npar <- unlist(lapply(models, function(x) length(x$ddf$ds$par)))
  modname <-  unlist(lapply(models, function(x) x$ddf$name.message))
  aic <-  unlist(lapply(models, function(x) x$ddf$criterion))
  chat.bigmod <- chat(models[[which.max(npar)]])
  qaic <- vector(mode="numeric", length = num.models)
  for (i in 1:num.models) {
    qaic[i] <- qaic(models[[i]], chat.bigmod)
  }
  nicetab <- data.frame(modname, npar, aic, qaic)
  return(nicetab)
}


'------ Run via shell script ------'
# get the input passed from the shell script (arg[1] for scenario, arg[2] for iter, arg[3] for nTraps, arg[4] for movemt behav (sol/grp), if grp arg[5] is group size)
# e.g., `Rscript --vanilla ./R/3d_analysisDS.R small 1 100 sol`
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
edir <- "/Volumes/LaCie/MERYL/PhD/Chapters/Paper 3/UnmarkedAbundance/"

# read input data (detections, trap data)
if (behav == 'sol'){
  #x4 <- read.csv(file=paste('Data/Detections/ctDat_', scenario, '_sol_', iter, '_J', nTraps,'.csv', sep=''))
  x4 <- read.csv(file=paste(edir, 'Data/Detections/ctDat_', scenario, '_sol_', iter, '_J', nTraps,'.csv', sep=''))
}
if (behav == 'grp'){
  #x4 <- read.csv(file=paste('Data/Detections/ctDat_', scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,'.csv', sep=''))
  x4 <- read.csv(file=paste(edir, 'Data/Detections/ctDat_', scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,'.csv', sep=''))
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
# subset to detections made every 2s
x8 <- x4 %>% 
  group_by(det_id) %>%
  mutate( difftime = time_seconds - first(time_seconds )) %>%
  filter( difftime %in% seq(0, 5000, 2) ) %>% # observations made at every 2s interval
  as.data.frame()

# get data into ds format
x8$dist <- x8$dist * 10 # convert units to metres
viewangle <- 60 # degrees
samfrac <- viewangle / 360
days = 50
ds_data <- data.frame(Region.Label = rep('SimLand', nrow(x8)), Area = rep(100, nrow(x8)), # Total area in km^2; multiplier = rep(samfrac, nrow(x8)),
                      Sample.Label = x8$Site, Effort = rep((days*12*60*60)/2, nrow(x8)), # Tk/t (t=2); Effort in terms of no. of 2s intervals CT was operable for
                      distance = x8$dist)
# create NA record for cams with 0 detections; important because it contains effort info
traps_no_det <- (1:nTraps)[!(1:nTraps %in% unique(x8$Site))] # Cams with no detections
no_det <- data.frame(Region.Label = rep('SimLand', length(traps_no_det)), Area = rep(100, length(traps_no_det)), # Total area in km^2; multiplier = rep(samfrac, nrow(x8)),
                      Sample.Label = traps_no_det, Effort = rep((days*12*60*60)/2, length(traps_no_det)), # Tk/t (t=2); Effort in terms of no. of 2s intervals CT was operable for
                      distance = rep(NA,length(traps_no_det)))

# combine data
ds_data <- rbind(ds_data, no_det)

"-------- Fit models --------" 
# fit detection function
#conversion <- convert_units("meter", NULL, "square kilometer")
trunc.list <- list(left=1, right=10)
mybreaks <- c(seq(1,10,1))
# do I need to do proper model selection??
hn0 <- ds(ds_data, transect = "point", key="hn", adjustment = NULL,
          cutpoints = mybreaks, truncation = trunc.list)
# hn1 <- ds(ds_data, transect = "point", key="hn", adjustment = "herm", order = 2,
#           cutpoints = mybreaks, truncation = trunc.list)
# hr0 <- ds(ds_data, transect = "point", key="hr", adjustment = NULL,
#           cutpoints = mybreaks, truncation = trunc.list)
# hr1 <- ds(ds_data, transect = "point", key="hr", adjustment = "cos", order = 2,
#           cutpoints = mybreaks, truncation = trunc.list)
# hr2 <- ds(ds_data, transect = "point", key="hr", adjustment = "cos", order = c(2,3),
#          cutpoints = mybreaks, truncation = trunc.list)
# 
# "-------- Select best model --------" 
#Finding the preferable model on the basis of the QAIC (see Howe et al. 2018 for further details)
# hr_mods <- qaic.pass1(hr0, hr1) #hr2
# hn_mods <- qaic.pass1(hn0, hn1)
# 
# winners <- list(eval(as.name(c('hn0','hn1')[which.min(hn_mods$qaic)])), # models with lower QAIC
#                 eval(as.name(c('hr0','hr1')[which.min(hr_mods$qaic)])))
# chats <- unlist(lapply(winners, function(x) chat(x)))
# best_mod <- winners[[which.max(chats)]]

best_mod <- hn0

# save r plot for reference
if (behav == 'sol'){
  png(paste0("Data/Estimates/Plots/DS-distance-",scenario, '_sol_', iter, '_J', nTraps,".png"), type="cairo", width = 650, height = 300, units = "px")
}
if (behav == 'grp'){
  png(paste0("Data/Estimates/Plots/DS-distance-",scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,".png"), type="cairo", width = 650, height = 300, units = "px")
}
par(mfrow=c(1,2))
par(oma=c(3,1,1,1))
par(mar=c(2,4,2,2))
plot(best_mod, main="Distance", xlab="Distance (m)",
     showpoints=FALSE, lwd=3, xlim=c(0, 11),pdf=T)
plot(best_mod, main="Distance", xlab="Distance (m)",
     showpoints=FALSE, lwd=3, xlim=c(0, 11),pdf=F)
curve((1-exp(-(3/x)^6))/(1+exp(10*(1-x))), 0, 10, col= 'blue', add=T) # simulated detection function 
mtext("Distance (m)",side=1,line=2,outer=TRUE,cex=1.3)
dev.off()

"-------- Model estimation --------" 
# estimate density
conversion <- convert_units("metre", NULL, "square kilometre")
peak.hr.dens <- dht2(best_mod, flatfile=ds_data, strat_formula = ~1, 
                     sample_fraction = samfrac, er_est = "P2", convert_units = conversion)
print(peak.hr.dens, report="abundance")

"-------- Save/write output --------" 
Replicate = iter
Scenario = scenario
Behaviour = behav
if (behav == 'sol'){ groupSize = NA }
Model = "CTDS"
Trap.effort = nTraps
Days.monitored = days
Nhat = peak.hr.dens$Abundance
Nhat.se = peak.hr.dens$Abundance_se
Nhat.lo = peak.hr.dens$LCI
Nhat.hi = peak.hr.dens$UCI
df = peak.hr.dens$df
ds.model = best_mod[["ddf"]][["name.message"]]

ds_output <- data.frame(Replicate, Scenario, Behaviour, groupSize, Model, Trap.effort, Days.monitored,
                        Nhat, Nhat.se, Nhat.lo, Nhat.hi, df, ds.model)

# append to existing master file
write.table( ds_output,  
             file="./Data/Estimates/ctds_master_test.csv", 
             append = T, 
             sep=',', 
             row.names=F, 
             col.names=F )

# # create empty .csv if not avail
# columns= c("Replicate", "Scenario", "Behaviour", "groupSize", "Model", "Trap.effort", "Days.monitored",
#            "Nhat", "Nhat.se", "Nhat.lo", "Nhat.hi", "df", "ds.model")
# mydf = data.frame(matrix(nrow = 0, ncol = length(columns)))
# colnames(mydf) = columns
# write.csv(mydf, file=paste('Data/Estimates/ctds_master_test.csv',sep=''), row.names = F)
