# ====================================
# Unmarked abundance evalaution
# Evaluation and plotting
# ====================================
library(ggpubr)
scrcols <- c("#dfc27d","#a6611a","#80cdc1","#018571","#f5f5f5")

"-------- REM --------"
rem_master <- read.csv('Data/Estimates/rem_master_hn.csv')

# pre-eval calculations
rem_res <- rem_master  %>%
  mutate(N = if_else(Scenario == 'small', 2000, 200, missing = NULL) ) %>%
  mutate(RB = (Nhat - N)/N )  %>% # Relative bias
  mutate(CV = Nhat.se/Nhat )  %>% # Precision
  mutate(Nhat.lo =  Nhat - 1.96*Nhat.se )  %>% 
  mutate(Nhat.hi =  Nhat + 1.96*Nhat.se )  %>% 
  mutate(Coverage = ifelse(N >= Nhat.lo & N <= Nhat.hi,1,0)) %>% 
  mutate(Trap.effort = as.character(Trap.effort)) %>%
  mutate(BehaviourX = if_else(Behaviour == 'grp', paste0(Behaviour, groupSize), as.character(Behaviour), missing = NULL)) %>%
  mutate(BehaviourX = if_else(BehaviourX %in% c('grp5','grp20'), paste0(Behaviour, '1'), paste0(Behaviour, '2'), missing = NULL)) %>%
  as.data.frame()

# summary
rem_stats <- rem_res  %>%
  mutate(N = if_else(Scenario == 'small', 2000, 200, missing = NULL) ) %>%
  mutate(RB = (Nhat - N)/N )  %>% # Relative bias
  mutate(CV_within = Nhat.se/Nhat )  %>% # Precision
  mutate(Coverage = ifelse(N >= Nhat.lo & N <= Nhat.hi,1,0)) %>% # coverage for RSF model is calculated differently (pre-calculated)
  mutate(Trap.effort = as.character(Trap.effort)) %>%
  mutate(BehaviourX = if_else(Behaviour == 'grp', paste0(Behaviour, groupSize), as.character(Behaviour), missing = NULL)) %>%
  # ----
group_by(Scenario, BehaviourX, Trap.effort) %>%
  mutate(CV_across = (sd(na.omit(Nhat)))/mean(Nhat, na.rm = T) )%>% 
  summarise(nConverged = length(N),
            RB = mean(RB, na.rm = T), 
            CV_within = mean(CV_within, na.rm = T), 
            CV_across = unique(CV_across), 
            Coverage = na.omit(length(Coverage[Coverage==1]))/na.omit(length(Coverage)))%>% # coverage for basic model
  as.data.frame()

# # Abundance N
# Nevals <- rem_master  %>%
#   mutate(N = if_else(Scenario == 'small', 2000, 200, missing = NULL) ) %>%
#   mutate(RB = (Nhat - N)/N )  %>% # Relative bias
#   mutate(CV = Nhat.se/Nhat )  %>% # CV
#   mutate(Coverage1 = ifelse(N >= Nhat - 1.96*Nhat.se & N <= Nhat + 1.96*Nhat.se,1,0)) %>% 
#   mutate(Coverage2 = ifelse(N >= Nhat.lo & N <= Nhat.hi,1,0)) %>% 
#   mutate(Trap.effort = as.character(Trap.effort)) %>%
#   mutate(BehaviourX = if_else(Behaviour == 'grp', paste0(Behaviour, groupSize), as.character(Behaviour), missing = NULL)) %>%
#   mutate(BehaviourX = if_else(BehaviourX %in% c('grp5','grp20'), paste0(Behaviour, '1'), paste0(Behaviour, '2'), missing = NULL)) %>%
#  as.data.frame()

# plot
rem_RB_plot <- ggviolin(rem_res, x = "BehaviourX", y = "RB", fill = "Trap.effort", palette = scrcols[3:4], facet.by = 'Scenario',
                  size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(plain(' RB('),italic(hat(N)), plain(')')))) + #ylim(-0.6,1.2) +
  geom_hline(yintercept=0, linetype="dashed", color = "black",size=0.2) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12)) +
  scale_x_discrete(labels=c('Solitary', 'Group', 'Group'))

ggpar(rem_RB_plot, main = 'REM (hn0)', font.main = c(14,"bold"))

rem_CV_plot <- ggviolin(rem_res, x = "BehaviourX", y = "CV", fill = "Trap.effort", palette = scrcols[3:4], facet.by = 'Scenario',
                  size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(plain(' CV('),italic(hat(N)), plain(')')))) + #ylim(-0.6,1.2) +
  #geom_hline(yintercept=0.2, linetype="dashed", color = "black",size=0.2) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12)) +
  scale_x_discrete(labels=c('Solitary', 'Group', 'Group'))

ggpar(rem_CV_plot, main = 'REM (hr0)', font.main = c(14,"bold"))

rem_EDR_plot <- ggviolin(rem_res, x = "BehaviourX", y = "edr", fill = "Trap.effort", palette = scrcols[3:4], facet.by = 'Scenario',
                   size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('Scenario') + ylab('EDR (units)') + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12)) +
  scale_x_discrete(labels=c('Solitary', 'Group', 'Group'))

ggpar(rem_EDR_plot, main = 'REM', font.main = c(14,"bold"))


###### INVESTIGATE ########
ggplot(rem_res, aes(x = edr, y = RB, color = as.factor(edr.model), shape = Scenario)) +
  geom_point(size = 2) +
  ggtitle('') + xlab('edr') + ylab("RB(N)") + 
  theme_bw() + 
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_text(size=12, face = 'bold'),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  font("xy.text", size = 9)

"-------- REST --------"
rest_master <- read.csv('Data/Estimates/rest_master_test.csv')

# pre-eval calculations
rest_res <- rest_master  %>%
  mutate(N = if_else(Scenario == 'small', 2000, 200, missing = NULL) ) %>%
  mutate(RB = (Nhat - N)/N )  %>% # Relative bias
  mutate(CV = Nhat.sd/Nhat )  %>% # Precision
  mutate(Coverage = ifelse(N >= Nhat.lo & N <= Nhat.hi,1,0)) %>% 
  mutate(Trap.effort = as.character(Trap.effort)) %>%
  mutate(BehaviourX = if_else(Behaviour == 'grp', paste0(Behaviour, groupSize), as.character(Behaviour), missing = NULL)) %>%
  mutate(BehaviourX = if_else(BehaviourX %in% c('grp5','grp20'), paste0(Behaviour, '1'), paste0(Behaviour, '2'), missing = NULL)) %>%
  as.data.frame()
rest_res$BehaviourX <- factor(rest_res$BehaviourX, levels = c('sol2', 'grp1', 'grp2'))

# summary
rest_stats <- rest_res  %>%
  mutate(N = if_else(Scenario == 'small', 2000, 200, missing = NULL) ) %>%
  mutate(RB = (Nhat - N)/N )  %>% # Relative bias
  mutate(CV_within = Nhat.sd/Nhat )  %>% # Precision
  mutate(Coverage = ifelse(N >= Nhat.lo & N <= Nhat.hi,1,0)) %>% # coverage for RSF model is calculated differently (pre-calculated)
  mutate(Trap.effort = as.character(Trap.effort)) %>%
  mutate(BehaviourX = if_else(Behaviour == 'grp', paste0(Behaviour, groupSize), as.character(Behaviour), missing = NULL)) %>%
  # ----
group_by(Scenario, BehaviourX, Trap.effort) %>%
  mutate(CV_across = (sd(na.omit(Nhat)))/mean(Nhat, na.rm = T) )%>% 
  summarise(nConverged = length(N),
            RB = mean(RB, na.rm = T), 
            CV_within = mean(CV_within, na.rm = T), 
            CV_across = unique(CV_across), 
            Coverage = na.omit(length(Coverage[Coverage==1]))/na.omit(length(Coverage)))%>% # coverage for basic model
  as.data.frame()

# plot
rest_RB_plot <- ggviolin(rest_res, x = "BehaviourX", y = "RB", fill = "Trap.effort", palette = c("#92C5DE", "#2166AC"), facet.by = 'Scenario',
                        size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(plain(' RB('),italic(hat(N)), plain(')')))) + ylim(-1.5,3) +
  geom_hline(yintercept=0, linetype="dashed", color = "black",size=0.2) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12)) +
  scale_x_discrete(labels=c('Solitary', 'Group', 'Group'))

ggpar(rest_RB_plot, main = 'REST', font.main = c(14,"bold"))


rest_CV_plot <- ggviolin(rest_res, x = "BehaviourX", y = "CV", fill = "Trap.effort", palette = c("#92C5DE", "#2166AC"), facet.by = 'Scenario',
                        size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(plain(' CV('),italic(hat(N)), plain(')')))) + ylim(-0.1,1.2) +
  #geom_hline(yintercept=0.2, linetype="dashed", color = "black",size=0.2) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12)) +
  scale_x_discrete(labels=c('Solitary', 'Group', 'Group'))

ggpar(rest_CV_plot, main = 'REST', font.main = c(14,"bold"))

# for low no. of reps
ggerrorplot(rest_res, x = "BehaviourX", y = "RB", color = "Trap.effort", palette = scrcols[3:4], facet.by = 'Scenario',
            size = 0.5, add = c("mean_se",'jitter'), width = 1, add.params = list(size=1)) +
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(plain(' RB('),italic(hat(N)), plain(')')))) +
  geom_hline(yintercept=0, linetype="dashed", color = "black",size=0.2) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12)) +
  scale_x_discrete(labels=c('Solitary', 'Group', 'Group')) 

"-------- CTDS --------"
ctds_master1 <- read.csv('Data/Estimates/ctds_master_test1.csv')

# pre-eval calculations
ctds_res <- ctds_master3  %>%
  mutate(N = if_else(Scenario == 'small', 2000, 200, missing = NULL) ) %>%
  mutate(RB = (Nhat - N)/N )  %>% # Relative bias
  mutate(CV = Nhat.se/Nhat )  %>% # CV
  mutate(Coverage = ifelse(N >= Nhat.lo & N <= Nhat.hi,1,0)) %>% 
  mutate(Trap.effort = as.character(Trap.effort)) %>%
  mutate(BehaviourX = if_else(Behaviour == 'grp', paste0(Behaviour, groupSize), as.character(Behaviour), missing = NULL)) %>%
  mutate(BehaviourX = if_else(BehaviourX %in% c('grp5','grp20'), paste0(Behaviour, '1'), paste0(Behaviour, '2'), missing = NULL)) %>%
  as.data.frame()

# summary
ctds_stats <- ctds_res  %>%
  mutate(N = if_else(Scenario == 'small', 2000, 200, missing = NULL) ) %>%
  mutate(RB = (Nhat - N)/N )  %>% # Relative bias
  mutate(CV_within = Nhat.se/Nhat )  %>% # Precision
  mutate(Coverage = ifelse(N >= Nhat.lo & N <= Nhat.hi,1,0)) %>% # coverage for RSF model is calculated differently (pre-calculated)
  mutate(Trap.effort = as.character(Trap.effort)) %>%
  mutate(BehaviourX = if_else(Behaviour == 'grp', paste0(Behaviour, groupSize), as.character(Behaviour), missing = NULL)) %>%
  # ----
group_by(Scenario, BehaviourX, Trap.effort) %>%
  mutate(CV_across = (sd(na.omit(Nhat)))/mean(Nhat, na.rm = T) )%>% 
  summarise(nConverged = length(N),
            RB = mean(RB, na.rm = T), 
            CV_within = mean(CV_within, na.rm = T), 
            CV_across = unique(CV_across), 
            Coverage = na.omit(length(Coverage[Coverage==1]))/na.omit(length(Coverage)))%>% # coverage for basic model
  as.data.frame()

# plot
ctds_RB_plot <- ggviolin(ctds_res, x = "BehaviourX", y = "RB", fill = "Trap.effort", palette = scrcols[1:2], facet.by = 'Scenario',
                        size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(plain(' RB('),italic(hat(N)), plain(')')))) + #ylim(-0.6,1.2) +
  geom_hline(yintercept=0, linetype="dashed", color = "black",size=0.2) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12)) +
  scale_x_discrete(labels=c('Solitary', 'Group', 'Group'))

ggpar(ctds_RB_plot, main = 'CT-DS', font.main = c(14,"bold"))

ctds_CV_plot <- ggviolin(ctds_res, x = "BehaviourX", y = "CV", fill = "Trap.effort", palette = scrcols[1:2], facet.by = 'Scenario',
                        size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(plain(' CV('),italic(hat(N)), plain(')')))) + #ylim(-0.6,1.2) +
  #geom_hline(yintercept=0.2, linetype="dashed", color = "black",size=0.2) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12)) +
  scale_x_discrete(labels=c('Solitary', 'Group', 'Group'))

ggpar(ctds_CV_plot, main = 'CT-DS', font.main = c(14,"bold"))

# for low no. of reps
ggerrorplot(ctds_res[ctds_res$Trap.effort == "25",], x = "BehaviourX", y = "RB", color = "ds.model", palette = scrcols[1:2], facet.by = 'Scenario',
            size = 0.5, add = c("mean_se",'jitter'), width = 1, add.params = list(size=1)) +
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(plain(' RB('),italic(hat(N)), plain(')')))) +
  geom_hline(yintercept=0, linetype="dashed", color = "black",size=0.2) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12)) +
  scale_x_discrete(labels=c('Solitary', 'Group', 'Group')) +

"-------- Investigate --------"
