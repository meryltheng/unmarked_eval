# ====================================
# Unmarked abundance evalaution
# Evaluate
# ====================================
library(mcmcOutput)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


"-------- Process model outputs --------"
resultL <- foreach(i=1:10, .combine=rbind, .errorhandling='stop', .packages=c(.packages())) %do% { 
    model <- 'REST'; replicate = i; scenario = 'large' # edit scenario here
    load(paste("Data/ModelOut/",model,"_",scenario,'_',replicate,".RData", sep=""))
    trueN = 10
    
    # empty matrix for results
    nrows=1
    output <- matrix(data=NA, nrow=nrows, ncol=32)
    colnames(output) <- c('Replicate', 'Scenario', 'N','Model', # meta info
                          'Nhat', 'Nhat.se', 'Nhat.lo', 'Nhat.hi','Nhat.mode','Nhat.Rhat',
                          # SC
                          'sigma', 'sigma.se','sigma.mode','sigma.Rhat',
                          'lam0', 'lam0.se','lam0.mode','lam0.Rhat',
                          'psi', 'psi.se','psi.mode','psi.Rhat',
                          # REST
                          'Dhat', 'Dhat.se','Dhat.Rhat',
                          'lamREST', 'lamREST.se','lamREST.Rhat',
                          #REM
                          'v', 'v.se','act', 'act.se'
    )
    
    # meta info
    output[,'Replicate'] <- replicate
    output[,'Scenario'] <- scenario
    output[,'N'] <- trueN
    
    # Convert to an mcmcOutput object and look at diagnostics
    mclist <- coda::mcmc.list(res)
    ( mco <- mcmcOutput(mclist) )
    summMCMC <- summary(mco)
    
    # Model output
    output[1,"Model"] <- model
    output[1,c('Nhat', 'Nhat.se','Nhat.Rhat')]  <- unlist(summMCMC['N',c(1:2,6)])
    output[1,c('Nhat.lo', 'Nhat.hi')] <- unlist(summMCMC['N',4:5])
    
    if(model == "SC"){
      output[1,c('sigma', 'sigma.se','sigma.Rhat')]  <- unlist(c(summMCMC['sigma',1:2]*10,summMCMC['sigma',6]))
      output[1,c('lam0', 'lam0.se','lam0.Rhat')]  <- unlist(summMCMC['lam0',c(1:2,6)])
      output[1,c('psi', 'psi.se','psi.Rhat')]  <- unlist(summMCMC['psi',c(1:2,6)])
      output[1,c('Nhat.mode')] <- Mode(do.call(c,lapply(res, function(x) x[,'N'] )))
      output[1,c('sigma.mode')] <- Mode(do.call(c,lapply(res, function(x) x[,'sigma'] )))
      output[1,c('lam0.mode')] <- Mode(do.call(c,lapply(res, function(x) x[,'lam0'] )))
      output[1,c('psi.mode')] <- Mode(do.call(c,lapply(res, function(x) x[,'psi'] )))
    }
    
    if(model == "REST"){
      output[1,c('Dhat', 'Dhat.se','Dhat.Rhat')]  <- unlist(summMCMC['D',c(1:2,6)])
      output[1,c('lamREST', 'lamREST.se','lamREST.Rhat')]  <- unlist(summMCMC['lambda',c(1:2,6)])
    }

    
    print(paste(replicate, scenario, 'done'))

    output
}

result <- rbind(resultS,resultM,resultL)
write.csv(result, file=paste('Data/Estimates/',model,'out.csv',sep=''), row.names = F)


"-------- Visualise --------"
library(ggpubr)

# result <- read.csv(file=paste('Data/Estimates/',model,'out.csv',sep=''))
resultSC <- read.csv(file=paste('Data/Estimates/SCout.csv',sep=''))
resultREM <- read.csv(file=paste('Data/Estimates/REMout.csv',sep=''))
resultREST <- read.csv(file=paste('Data/Estimates/RESTout.csv',sep=''))
result <- rbind(resultSC,resultREST,resultREM,resultREMh)

resultX <- result %>%
  mutate(Nest = ifelse( Model == "SC",
         Nhat.mode, Nhat))  %>%
  mutate(Nbias = (Nest - N) / N) %>%
  filter(Nhat.Rhat <= 1.1 | is.na(Nhat.Rhat)) 
resultX$Scenario <- factor(resultX$Scenario, levels = c('small', 'medium', 'large','habitatCor'))

# Compare estimates across models
colz <- c("darkturquoise","limegreen","darkgreen","indianred3","indianred4")
ggerrorplot(resultX, x = "Scenario", y = "Nbias", color = "Model", palette = colz, 
            desc_stat = "mean_sd", add = 'jitter', add.params = list(size = 1.5, alpha = 0.3)) + 
  ggtitle('') + xlab('Scenario') + ylab('Rel. bias (N)') +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") 

# Compare mean and mode in SC model
SC1 <- ggerrorplot(resultSC, x = "Scenario", y = "Nhat", color = "Model", palette = "jco", 
            desc_stat = "mean_sd", add = 'jitter', add.params = list(size = 1.5, alpha = 0.3)) + 
  ggtitle('') + xlab('Scenario') + ylab('N mean') +
  geom_hline(yintercept=50, linetype="dashed", color = "grey") 
SC2 <- ggerrorplot(resultSC, x = "Scenario", y = "Nhat.mode", color = "Model", palette = "jco", 
            desc_stat = "mean_sd", add = 'jitter', add.params = list(size = 1.5, alpha = 0.3)) + 
  ggtitle('') + xlab('Scenario') + ylab('N mode') +
  geom_hline(yintercept=50, linetype="dashed", color = "grey") 
ggarrange(SC1,SC2, common.legend = T)

##################
resultSC <- read.csv(file=paste('Data/Estimates/SCout.csv',sep=''))
resultREST <- read.csv(file=paste('Data/Estimates/RESTout.csv',sep=''))
resultREM <- read.csv(file=paste('Data/Estimates/REMout.csv',sep=''))
result <- rbind(resultSC,resultREST,resultREM)

result <- data.frame(Model = c(rep('SC.mean',20), rep('SC.mode',20), 
                               rep('REST.all',20), rep('REST.rmoutlier',20),
                               rep('REM',20)),
                     N = c(resultSC$Nhat, resultSC$Nhat.mode, resultREST$Nhat, resultRESTo$Nhat, resultREM$Nhat),
                     Scenario = c(paste(resultSC$Scenario), paste(resultSC$Scenario), 
                                  paste(resultREST$Scenario), paste(resultRESTo$Scenario),
                                  paste(resultREM$Scenario)))

