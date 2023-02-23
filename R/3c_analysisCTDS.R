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
# e.g., `Rscript --vanilla ./R/3c_analysisCTDS.R hi 1 100 sol`
# e.g., `Rscript --vanilla ./R/3c_analysisCTDS.R hi 1 100 grp 25`
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

# read input data (detections)
if (behav == 'sol'){
  x4 <- read.csv(file=paste('Data/Detections/ctDat_', scenario, '_sol_', iter, '_J', nTraps,'.csv', sep=''))
}
if (behav == 'grp'){
  x4 <- read.csv(file=paste('Data/Detections/ctDat_', scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,'.csv', sep=''))
}

if(subsetDat == TRUE){
  x4 <- 
    x4 %>%
    filter(time %in% 1:(12 * 9 * days)) 
}

"-------- Prep model input --------"
# subset to detections made every 2s
x8 <- x4 %>% 
  filter( time_seconds %in% seq(0, max(time_seconds), 2) ) %>%
  group_by(Site, ID) %>%
  distinct(time_seconds, .keep_all= TRUE) %>%
  as.data.frame()

# get data into ds format
Area <- 1 # Area in km^2 to evaluate density in
Effort <- (days*9*60*60)/2 # Effort in terms of no. of 2s intervals CT was operable for
viewangle <- 40 # degrees
samfrac <- viewangle / 360 # fraction of circular area covered by FOV

ds_data <- data.frame(Region.Label = rep('SimLand', nrow(x8)), # give study area a random name
                      Area = rep(Area, nrow(x8)), # Total area in km^2
                      Sample.Label = x8$Site, # camera ID
                      Effort = rep(Effort, nrow(x8)), # effort per Site/cam
                      distance = x8$dist) # distance records

# create NA record for cams with 0 detections
traps_no_det <- (1:nTraps)[!(1:nTraps %in% unique(x8$Site))] # Cams with no detections
no_det <- data.frame(Region.Label = rep('SimLand', length(traps_no_det)), 
                     Area = rep(Area, length(traps_no_det)), 
                     Sample.Label = traps_no_det, 
                     Effort = rep(Effort, length(traps_no_det)), 
                     distance = rep(NA,length(traps_no_det)))

# combine data
ds_data <- rbind(ds_data, no_det)
head(ds_data)

"-------- Fit models --------" 
# fit detection function
conversion <- convert_units("meter", NULL, "square kilometer")
dist_range = c(0,10)
mybreaks <- seq(dist_range[1],dist_range[2],1)
trunc.list <- list(left=dist_range[1], right=dist_range[2])

# default hazard-rate function
best_mod <- ds(ds_data, transect = "point", key="hr", adjustment = NULL,
               cutpoints = mybreaks, truncation = trunc.list,
               convert_units = conversion)


# # save r plot for reference
# if (behav == 'sol'){
#   png(paste0("Data/Estimates/Plots/DS-distance-",scenario, '_sol_', iter, '_J', nTraps,".png"), type="cairo", width = 650, height = 300, units = "px")
# }
# if (behav == 'grp'){
#   png(paste0("Data/Estimates/Plots/DS-distance-",scenario, '_grp_', iter, '_gSize', groupSize,'_J', nTraps,".png"), type="cairo", width = 650, height = 300, units = "px")
# }
# par(mfrow=c(1,2))
# par(oma=c(3,1,1,1))
# par(mar=c(2,4,2,2))
# plot(best_mod, main="Distance", xlab="Distance (m)",
#      showpoints=FALSE, lwd=3, xlim=c(0, 10),pdf=T)
# plot(best_mod, main="Distance", xlab="Distance (m)",
#      showpoints=FALSE, lwd=3, xlim=c(0, 10),pdf=F)
# curve(1-exp(-(3/x)^6), 0, 10, col= 'blue', add=T) # simulated detection function 
# mtext("Distance (m)",side=1,line=2,outer=TRUE,cex=1.3)
# dev.off()

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
y = nrow(x8)
Dhat = peak.hr.dens$Abundance
Dhat.se = peak.hr.dens$Abundance_se
Dhat.lo = peak.hr.dens$LCI
Dhat.hi = peak.hr.dens$UCI
df = peak.hr.dens$df
ds.model = best_mod[["ddf"]][["name.message"]]
edr = effec_radius$EDR
edr.se = effec_radius$se.EDR
edr.AIC = as.numeric(AIC(best_mod)['AIC'])
edr.QAIC = as.numeric(qaic.pass1(best_mod)['qaic'])
edr.chat = min(chats)
gof_chisq = gof_fit[["chisquare"]][["chi1"]][["chisq"]]
gof_df = gof_fit[["chisquare"]][["chi1"]][["df"]]
gof_rel = gof_chisq/y

ds_output <- data.frame(Replicate, Scenario, Behaviour, groupSize, Model, Trap.effort, Days.monitored, y,
                        Dhat, Dhat.se, Dhat.lo, Dhat.hi, df, ds.model, edr, edr.se, edr.AIC, edr.QAIC, edr.chat, #QAIC and chat only for model selection
                        gof_chisq, gof_df, gof_rel)

# append to existing master file
write.table( ds_output,  
             file="./Data/Estimates/ctds_imperfectDet.csv", 
             append = T, 
             sep=',', 
             row.names=F, 
             col.names=F )

# # create empty .csv if not avail
# columns= c("Replicate", "Scenario", "Behaviour", "groupSize", "Model", "Trap.effort", "Days.monitored", "y",
#            "Dhat", "Dhat.se", "Dhat.lo", "Dhat.hi", "df", "ds.model", "edr", "edr.se", "edr.AIC", "edr.QAIC", "edr.chat", 
#            "gof_chisq", "gof_df", "gof_rel")
# mydf = data.frame(matrix(nrow = 0, ncol = length(columns)))
# colnames(mydf) = columns
# write.csv(mydf, file=paste('Data/Estimates/ctds_imperfectDet.csv',sep=''), row.names = F)