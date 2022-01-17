# ====================================
# Unmarked abundance evalaution
# Movement simulations (habitat-driven movement)
# ====================================

'------ Create layers ------'
# simulate continuous habitat
library(NLMR)
env <- nlm_gaussianfield(ncol = 1000, nrow = 1000,
                                    autocorr_range = 50,
                                    mag_var = 1,
                                    nug = 0)
# convert continuous habitat to binary
hab <- raster(xmn=0, xmx=1000, ymn=0, ymx=1000, ncols=1000, nrows=1000)
hab[env > 0.475] <- 1
hab[env <= 0.475] <- 0
plot(hab)
# save(env,hab, file='Data/habitat_layer.RData')

# OR load
load('Data/habitat_layer.RData')

"-------- Simulation functions --------"
# draw animal starting locations
genStartPoints <- function(r = hab, hab_prob = 0.6, fact=25, n_foragers = 100, emptybordersize = 0, toReplace= T){
  r <- aggregate(r, fact=fact)
  r[] <- (r[]-min(r[]))/(max(r[])-min(r[])) 
  r[r<0.5] <- 0
  r[] <- (r[]-min(r[]))/(max(r[])-min(r[])) 
  cell = sample(1:ncell(r),n_foragers, prob=r[], replace=toReplace)
  centres = xyFromCell(r,cell)
  rx = centres[,1] + runif(nrow(centres), -0.5*res(r)[1], 0.5*res(r)[1])
  ry = centres[,2] + runif(nrow(centres), -0.5*res(r)[2], 0.5*res(r)[2])
  StartPoints <- data.frame(x=rx + emptybordersize,y=ry + emptybordersize)
  colnames(StartPoints) <- c('x','y')
  return(StartPoints)
}

genCentres <- function(X = ACs[i,], habitat = T, r = hab, hab_val = 1, n_X = 5:8, buffer_scale = c(6,8.5), n_step = n){
  X <- as.numeric(X)
  n_AC <- sample(n_X,1)
  centres <- matrix(NA,nrow=n_AC,2)
  colnames(centres) <- c('x','y')
  
  buffer = rgamma(1, scale = buffer_scale[1], shape = buffer_scale[2])
  buffApart <- buffer/3

  for (i in 1:n_AC){ 
    repeatCounter <- 1
    repeat{ # rejection sampling
      cand <- cbind(x=runif(1,X[1]-buffer,X[1]+buffer), y=runif(1,X[2]-buffer,X[2]+buffer))
      # reflective boundaries
      if(cand[1] <= 0){
        cand[1] <- abs(cand[1])
      }
      if(cand[1] >= 1000){
        cand[1] <- 1000 - cand[1]
      }
      if(cand[2] <= 0){
        cand[2] <- abs(cand[2])
      }
      if(cand[2] >= 1000){
        cand[2] <- 1000 - cand[2]
      }
      
      distapart <- apply(centres, 1, function(x) sqrt( (x[1] - cand[1])^2 + (x[2] - cand[2])^2) )
      
      if (habitat == T) {      
      cell <- cellFromXY(r, cand)
      cells_ad <- adjacent(r,cell, directions = 16, pairs=FALSE)
      
      # check if candidate centre is in habitat
      if (raster::extract(r, cand) == hab_val & 
          sum(raster::extract(r, cells_ad),na.rm = T) == hab_val*16 &
          all(distapart > buffApart, na.rm = T) ) { break }
      } # habitat (heterogenous landscape)
      else {
        if (all(distapart > buffApart, na.rm = T) ) { break }
      } # no habitat (homogeneous landscape)
      repeatCounter <- repeatCounter + 1
      if (repeatCounter >=25){
        X <- c(runif(1, 50, 950), runif(1, 50, 950))
      } # habitat

    } # repeat
    # input
    centres[i,] <- cand
  }
  return(centres)
}


# simulate habitat-driven movement
simMoveHab_v1 <- function(){
  xy <- matrix(NA, nrow = n, ncol = 2) # empty matrix for simulated locs
  xy[1, ] <- xy0 # starting loc
  for (i in 2:n) {
    if (i == 2){
      ta <- runif(20, -pi, pi)
      z <- sample(1:20, 1)
      slen <- rgamma(20, scale = 0.8, shape = 1)
    }
    # cand locations
    if (raster::extract(hab, xy[i-1, ,drop = FALSE]) == 2) { # HABITAT
      scl <- 3
      shp <- 3
      # ta <- runif(20, -pi, pi)
      ta <- rwrpnorm(20, ta[z], 0.3) # angle correlated to previous angle
    } else { # MATRIX
      scl <- 3
      shp <- 4.5
      ta <- rwrpnorm(20, ta[z], 0.6) # angle correlated to previous angle
    }
    h <- raster::extract(env, xy[i-1, ,drop = FALSE])
    slen <- 0.5*(slen[z] + rgamma(20, scale = scl, shape = shp) * exp(-h) ) # speed correlated to previous speed 
    x1 <- xy[i - 1, 1] + (cos(ta) * slen) 
    y1 <- xy[i - 1, 2] + (sin(ta) * slen)
    
    # next direction is influenced by env layer, weighted by distance (percepRadius)
    if (raster::extract(hab, xy[i-1, ,drop = FALSE]) == 2) { # HABITAT
      w <- sampleRadialCands(x0=xy[i - 1, 1], y0=xy[i - 1, 2], ang = ta, percepRadius = 5, distFreq = 1, r = env)
    } else { # MATRIX
      w <- sampleRadialCands(x0=xy[i - 1, 1], y0=xy[i - 1, 2], ang = ta, percepRadius = 50, distFreq = 2, r = env)
    }
    w[is.na(w)] <- 0
    
    # landscape torus
    z <- sample(20, 1, prob = w)
    withinLandscape = isTRUE(x1[z] >= 0 & x1[z] <= 1000 & y1[z] >= 0 & y1[z] <= 1000)
    if (withinLandscape==FALSE){
      if(x1[z] <= 0){
        x1[z] <- 1000 + x1[z]
      }
      if(x1[z] >= 1000){
        x1[z] <- x1[z] - 1000 
      }
      if(y1[z] <= 0){
        y1[z] <- 1000 + y1[z]
      }
      if(y1[z] >= 1000){
        y1[z] <- y1[z] - 1000 
      }
    }
    xy[i, ] <- c(x1[z], y1[z])
  }
  xy <- xy[, 1:2]
  xy <- data.frame(xy)
  names(xy) <- c("x", "y")
  xy$hour <- rep(rep(1:12, each = 12),days) # hour of the day
  xy$id <- 1:nrow(xy)
  xy$hab <- raster::extract(hab, xy[, 1:2])
  xy$env <- raster::extract(env, xy[, 1:2])
  xy$ID <- id
  return(xy)
}

simMoveHab_v2 <- function(){
  AC <- vector()
  xy <- matrix(NA, nrow = n, ncol = 2) # empty matrix for simulated locs
  xy <- data.frame(xy)
  names(xy) <- c("x", "y")
  xy[1, ] <- locs[[id]][1,] # starting loc
  for (i in 2:n) {
    if (i == 2){
      slen <- rgamma(20, scale = 2, shape = 2)
      ta <- runif(20, -pi, pi)
      z <- sample(1:20, 1)
      AC[1] <- current_c <- 1
    }
    
    # cand locations
    if (raster::extract(hab, xy[i-1, ,drop = FALSE]) == 1) { # HABITAT
      scl <- 2
      shp <- 2
      slen <- 0.5*(slen[z] + rgamma(20, scale = scl, shape = shp)) # speed correlated to previous speed 
      ta <- rwrpnorm(20, ta[z], 0.3) # angle correlated to previous angle
    } else { # MATRIX
      scl <- 2
      shp <- 3.5
      slen <- 0.5*(slen[z] + rgamma(20, scale = scl, shape = shp)) # speed correlated to previous speed 
      ta <- rwrpnorm(20, ta[z], 0.6) # angle correlated to previous angle
    }
    
    AC_state = sample(1:2, 1, prob = c(0.995, 0.005))
    # cand locations
    if (AC_state == 1) { # remain in current AC
      x0 <- locs[[id]][current_c,1]
      y0 <- locs[[id]][current_c,2]
    } else { # change AC
        # dist of current loc to next AC
        centres <- as.data.frame(locs[[id]][-current_c,])
        dist_to_ACs <- apply(centres, 1, function(x) {
          # dist( rbind(xy[i - 1, 1:2], x), method='euclidean')
          sqrt( min( abs(xy[i - 1, 1] - x[1]), 1000 - abs( xy[i - 1, 1] - x[1] ))^2  +  min(abs(xy[i - 1, 2] - x[2]), 1000 - abs( xy[i - 1, 2] - x[2] ))^2 )
        })
    
        # pick next AC
        phi <- (1/dist_to_ACs) / sum(1/dist_to_ACs) # inverse probabilities based on nearest distance
        loc2 = centres[sample(1:nrow(centres), 1, prob = phi),]
        next_c <- which(locs[[id]][,1] %in% loc2[1] & locs[[id]][,2] %in% loc2[2]) # get index for next centre
        current_c <- next_c
        x0 <- locs[[id]][current_c,1]
        y0 <- locs[[id]][current_c,2]
        } # change AC
      
      # next direction is biased towards AC
      x1 <- xy[i - 1, 1] + (cos(ta) * slen) 
      y1 <- xy[i - 1, 2] + (sin(ta) * slen)
      # bias towards AC
      # d1 <- sqrt( (x1 - x0)^2 + (y1 - y0)^2 ) # distance to current AC
      d1 <- vector() # distance to current AC, accounting for torus
      for (a in 1:20){ 
        d1[a] <- sqrt( min( abs(x1[a] - x0), 1000 - abs( x1[a] - x0 ))^2  +  min( abs(y1[a] - y0), 1000 - abs( y1[a] - y0 ))^2 )
      }
      w <- dexp(d1, rate = 1/25) 
      w[is.na(w)] <- 0
      z <- sample(20, 1, prob = w)
      
      # landscape torus
      withinLandscape = isTRUE(x1[z] >= 0 & x1[z] <= 1000 & y1[z] >= 0 & y1[z] <= 1000)
      if (withinLandscape==FALSE){
        if(x1[z] <= 0){
          x1[z] <- 1000 + x1[z]
        }
        if(x1[z] >= 1000){
          x1[z] <- x1[z] - 1000 
        }
        if(y1[z] <= 0){
          y1[z] <- 1000 + y1[z]
        }
        if(y1[z] >= 1000){
          y1[z] <- y1[z] - 1000 
        }
      }
      # fill output
      xy[i, ] <- c(x1[z], y1[z])
      AC[i] <- current_c
    } # i 
  xy <- xy[1:n,]
  xy$id <- 1:nrow(xy)
  xy$ID <- id
  xy$AC <- AC
  return(xy)
}

simMoveHab_v3 <- function(){
  AC <- vector()
  xy <- matrix(NA, nrow = n, ncol = 2) # empty matrix for simulated locs
  xy <- data.frame(xy)
  names(xy) <- c("x", "y")
  xy[1, ] <- rnorm(2,locs[[id]][1,],1) # starting loc
  for (i in 2:n) { # for second location
    if (i == 2){
      ta <- runif(20, -pi, pi)
      slen <- rgamma(20, scale = 2, shape = 2)
      AC[1] <- current_c <- 1
    } else { # for subsequent locations
      # cand locations
      if (raster::extract(hab, xy[i-1, ,drop = FALSE]) == 1) { # HABITAT
        scl <- 2
        shp <- 2
        # ta <- runif(20, -pi, pi)
        ta <- rwrpnorm(20, ta[z], 0.3) # angle correlated to previous angle
      } else { # MATRIX
        scl <- 2
        shp <- 3.5
        ta <- rwrpnorm(20, ta[z], 0.5) # angle correlated to previous angle
      }
      slen <- 0.5*(slen[z] + rgamma(20, scale = scl, shape = shp)) # speed correlated to previous speed 
    }
    
    AC_state = sample(1:2, 1, prob = c(0.997, 0.003))
    # cand locations
    if (AC_state == 1) { # remain in current AC
      x0 <- locs[[id]][current_c,1]
      y0 <- locs[[id]][current_c,2]
      
    } else { # change AC
      # dist of current loc to next AC
      centres <- as.data.frame(locs[[id]][-current_c,])
     # dist_to_ACs <- apply(centres, 1, function(x) {
     #   sqrt( min( abs(xy[i - 1, 1] - x[1]), 1000 - abs( xy[i - 1, 1] - x[1] ))^2  +  min(abs(xy[i - 1, 2] - x[2]), 1000 - abs( xy[i - 1, 2] - x[2] ))^2 )
     # })
      
      # pick next AC
     # phi <- (1/dist_to_ACs) / sum(1/dist_to_ACs) # inverse probabilities based on nearest distance
      phi = rep(1/nrow(centres), nrow(centres))
      loc2 = centres[sample(1:nrow(centres), 1, prob = phi),]
      next_c <- which(locs[[id]][,1] %in% loc2[1] & locs[[id]][,2] %in% loc2[2]) # get index for next centre
      current_c <- next_c
      x0 <- locs[[id]][current_c,1]
      y0 <- locs[[id]][current_c,2]
    } # change AC
    
    x1 <- xy[i - 1, 1] + (cos(ta) * slen) 
    y1 <- xy[i - 1, 2] + (sin(ta) * slen)
    
    # bias towards AC
    # d1 <- sqrt( (x1 - x0)^2 + (y1 - y0)^2 ) # distance to current AC
    d1 <- vector() # distance to current AC, accounting for torus
    for (a in 1:20){ 
      d1[a] <- sqrt( min( abs(x1[a] - x0), 1000 - abs( x1[a] - x0 ))^2  +  min( abs(y1[a] - y0), 1000 - abs( y1[a] - y0 ))^2 )
    }
    w <- dexp(d1, rate = 1/25) 
    w[is.na(w)] <- 0
    z <- sample(20, 1, prob = w)

    # landscape torus
    withinLandscape = isTRUE(x1[z] >= 0 & x1[z] <= 1000 & y1[z] >= 0 & y1[z] <= 1000)
    if (withinLandscape==FALSE){
      if(x1[z] <= 0){
        x1[z] <- 1000 + x1[z]
      }
      if(x1[z] >= 1000){
        x1[z] <- x1[z] - 1000 
      }
      if(y1[z] <= 0){
        y1[z] <- 1000 + y1[z]
      }
      if(y1[z] >= 1000){
        y1[z] <- y1[z] - 1000 
      }
    }
    # fill output
    xy[i, ] <- c(x1[z], y1[z])
    AC[i] <- current_c
  } # i 
  xy <- xy[1:n,]
  xy$ID <- id
  xy$time <- 1:nrow(xy)
  xy$AC <- AC
  return(xy)
}

simMoveNull_v3 <- function(){
  AC <- vector()
  xy <- matrix(NA, nrow = n, ncol = 2) # empty matrix for simulated locs
  xy <- data.frame(xy)
  names(xy) <- c("x", "y")
  xy[1, ] <- rnorm(2,locs[[id]][1,],1) # starting loc
  for (i in 2:n) {
    if (i == 2){ # only for second location
      slen <- rgamma(20, scale = 2, shape = 2)
      ta <- runif(20, -pi, pi)
      AC[1] <- current_c <- 1
    } else { # for subsequent locations
      # cand locations
      scl <- 1.4
      shp <- 3.2
      ta <- rwrpnorm(20, ta[z], 0.3) # angle correlated to previous angle
      slen <- 0.5*(slen[z] + rgamma(20, scale = scl, shape = shp)) # speed correlated to previous speed 
    }
    
    AC_state = sample(1:2, 1, prob = c(0.997, 0.003))
    # cand locations
    if (AC_state == 1) { # remain in current AC
      x0 <- locs[[id]][current_c,1]
      y0 <- locs[[id]][current_c,2]
      
    } else { # change AC
      # dist of current loc to next AC
      centres <- as.data.frame(locs[[id]][-current_c,])

      # pick next AC
      phi = rep(1/nrow(centres), nrow(centres))
      loc2 = centres[sample(1:nrow(centres), 1, prob = phi),]
      next_c <- which(locs[[id]][,1] %in% loc2[1] & locs[[id]][,2] %in% loc2[2]) # get index for next centre
      current_c <- next_c
      x0 <- locs[[id]][current_c,1]
      y0 <- locs[[id]][current_c,2]
    } # change AC
    
    x1 <- xy[i - 1, 1] + (cos(ta) * slen) 
    y1 <- xy[i - 1, 2] + (sin(ta) * slen)
    
    # bias towards AC
    # d1 <- sqrt( (x1 - x0)^2 + (y1 - y0)^2 ) # distance to current AC
    d1 <- vector() # distance to current AC, accounting for torus
    for (a in 1:20){ 
      d1[a] <- sqrt( min( abs(x1[a] - x0), 1000 - abs( x1[a] - x0 ))^2  +  min( abs(y1[a] - y0), 1000 - abs( y1[a] - y0 ))^2 )
    }
    w <- dexp(d1, rate = 1/25) 
    w[is.na(w)] <- 0
    z <- sample(20, 1, prob = w)
    
    # landscape torus
    withinLandscape = isTRUE(x1[z] >= 0 & x1[z] <= 1000 & y1[z] >= 0 & y1[z] <= 1000)
    if (withinLandscape==FALSE){
      if(x1[z] <= 0){
        x1[z] <- 1000 + x1[z]
      }
      if(x1[z] >= 1000){
        x1[z] <- x1[z] - 1000 
      }
      if(y1[z] <= 0){
        y1[z] <- 1000 + y1[z]
      }
      if(y1[z] >= 1000){
        y1[z] <- y1[z] - 1000 
      }
    }
    # fill output
    xy[i, ] <- c(x1[z], y1[z])
    AC[i] <- current_c
  } # i 
  xy <- xy[1:n,]
  xy$ID <- id
  xy$time <- 1:nrow(xy)
  xy$AC <- AC
  return(xy)
}


"-------- Prep input --------"
N <- 100 # true abundance
days <- 100
n <- 6 * 12 * days # one relocation every 5 min

ACs <- genStartPoints(r = hab, hab_prob = 0.6, fact = 20, n_foragers = N, emptybordersize = 0, toReplace = F)
points(ACs, pch=3)

'------ Run sims (Treatment) ------'
scenario = 'habitat'
locs <- simDat <- list()
for (id in 1:N){
  centres <- genCentres(X = ACs[id,], habitat = T, r = hab, hab_val = 1, n_X = 6:9, buffer_scale = c(6.5,12))
  locs[[id]] <- centres
  print(paste('Animal',id,'done'))
}

for (id in 1:N){
  simDat[[id]] <- simMoveHab_v3()
  print(paste('Animal',id,'done'))
}

# save
moveSims <- do.call(rbind, simDat)
iter=10
write.csv(moveSims,row.names=FALSE,file=paste('Data/MovementSims/habitat/simDat_', scenario, '_', iter, '.csv', sep=''))

# plot
plot(hab,col=alpha(grey.colors(10,rev=T),0.5))
moveSims <- # 
  moveSims2 %>% 
  group_by(ID) %>% 
  filter(row_number() %in% 1:(6 * 12 * 50)) 
plotTraject(df = moveSims, plotit = rep(1,N), torus = T) # all

# check
moveSims <- moveSims %>%
  group_by(ID) %>%
  mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
  as.data.frame()

summary(moveSims$step[moveSims$step< 950]) * 6 * 10

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#20.93  183.29  247.62  269.47  332.65 1620.87    6000 # habitat
#23.16  272.74  342.53  358.99  427.40 1292.27    6000 # null

###
plot(hab,col=alpha(grey.colors(10,rev=T),0.5), xlim = c(600,900), ylim = c(600,900))
plotTraject(df = xy, plotit = rep(1,1), torus = T) # all
plotTrajectTemp(df = xy, plotit = rep(1,1))

id = 91
points(locs[[id]], pch=3, col = 'blue')
plotTraject(df = simDat[[id]], plotit = rep(1,1), torus = T) # all

xy <- xy %>%
  #group_by(ID) %>%
  mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
  as.data.frame()

summary(xy$step[xy$step < 950]) * 6 * 10

# summarise realised movement characteristics
library(adehabitatHR)
df <- # 
  moveSims %>% 
  group_by(ID) %>% 
  filter(row_number() %in% seq(6*12,nbObs,6*12))
  
  df <- df[, c('ID','x','y')] # don't mess up original df for later analyses
  coordinates(df) <- ~x+y
  MCP95_hour <- mcp.area(df, percent = 95, unin = 'm', unout = 'm2',plotit=F)
  MCP95 <- as.vector(apply(MCP95_hour,2, function(x) as.numeric(x)))

summary(MCP95) / (100*100)
summary(MCP95[MCP95<10*15000]) / (100*100)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.7632  1.6150  2.2192  7.5199  6.5745 89.7980 # habitat
#1.105   1.771   2.273   4.809   3.778  51.455  # null

'------ Run sims (Null) ------'
scenario = 'null'
# null
ACs <- cbind(x = runif(N, 0, 1000), y = runif(N, 0, 1000))
locs <- simDat <- list()
for (id in 1:N){
  centres <- genCentres(X = ACs[id,], habitat = F, n_X = 6:9, buffer_scale = c(6.5,12))
  locs[[id]] <- centres
  print(paste('Animal',id,'done'))
}

for (id in 1:N){
  simDat[[id]] <- simMoveNull_v3()
  print(paste('Animal',id,'done'))
}

# save
moveSims <- do.call(rbind, simDat)
iter=10
write.csv(moveSims,row.names=FALSE,file=paste('Data/MovementSims/habitat/simDat_', scenario, '_', iter, '.csv', sep=''))
