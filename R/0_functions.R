# ====================================
# Unmarked abundance evalaution
# Functions
# ====================================

"MOVEMENT SIMULATIONS"
# Generate activity centres (ACs)
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

simMove <- function(nbObs = nbObs, burnin = 12*12*1, xy0 = c(0,0), ID = 1, 
                    aind = aind, bind = bind, # step length params
                    kappac = kappac, etac = etac, # concentration param (i.e., directedness); weight of CRW:BRW component
                    xlim = c(0,10000), ylim = c(0,10000), torus = T)
{
  # adapted from Langrock et al. (2014): https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12155 
  # Choice of the initial position of the individual and of the location of the individual’s centre of attraction (C)
  timesteps <- nbObs + burnin
  X <- matrix(rep(NA, 2 * timesteps), ncol = 2)
  X[1, ] <- xy0
  C <- matrix(xy0, byrow = T, ncol = 2)
  
  phi <- 0
  
  for (k in 1:(timesteps - 1)) {
    ## start generation of centroid location at time k+1
    coo <- C[, ] - X[k, ] # xy difference btw individual's centre and its current location
    mu <- Arg(coo[1] + (0+1i) * coo[2]) # direction towards individual's centre
    if (mu < 0)
      mu <- mu + 2 * pi # correct for direction
    # expected movement direction: weighted average of CRW and BRW [#Eq-1.2]
    mu.av <- Arg(etac * exp(phi * (0+1i)) + (1 - etac) * exp(mu * (0+1i)))
    phi <- rvm(1, mean = mu.av, k = kappac) # draw next direction [#Eq-1.1]
    if (k==1){
      step.len <-  rgamma(1, shape = aind^2/bind^2, scale = bind^2/aind) # first step-length
    } else{
      step.len <- 0.5*(step.len + rgamma(1, shape = aind^2/bind^2, scale = bind^2/aind)) # subsequent step-lengths are correlated to previous
    }
    step <- step.len * c(Re(exp((0+1i) * phi)), Im(exp((0+1i) * phi))) # dist in xy to move 
    X[k + 1, ] <- X[k, ] + step # next loc 
  }
  # compile output
  xy <- as.data.table(X[(burnin+1):timesteps,]) # discard burnin
  colnames(xy) <- c('x','y')
  xy$ID <- rep(ID, each=nbObs)
  xy$time <- rep(1:nbObs)
  
  # correct coordinates for torus (data.table for fast math)
  if (torus == T){
    xy[x <= xlim[1], x := x + xlim[2]]
    xy[x >= xlim[2], x := x - xlim[2]]
    xy[y <= ylim[1], y := y + ylim[2]]
    xy[y >= ylim[2], y := y - ylim[2]]
  }
  return(xy)
}

simMoveGroup <- function(nbObs = nbObs, burnin = 12*12*1,
                         groupSize = groupSize, groupID = 1, memberID = 1:groupSize,
                         # group centroid
                         xy0 = c(0,0), # activity centre loc (also starting loc)
                         ac = 32, bc = 20, # step-length params
                         kappac = kappac, etac = etac, # concentration param (i.e., directedness); weight of CRW:BRW component
                         # individuals
                         gamma = gamma, # state transition probability matrix
                         aind = c(55, 55), bind = c(32, 32), # state-specific step-length params
                         kappaind = kappaind, # state-specific concentration params (i.e., directedness)
                         xlim = c(0,10000), ylim = c(0,10000), torus = T)
{
  # adapted from Langrock et al. (2014): https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12155 
  # state params
  nState = 2
  s.ty <- c(1, 0) # (1: BRW; 0: CRW)
  delta <- solve(diag(nState) - t(gamma) + 1, rep(1, nState)) # corresponding stationary distribution of the Markov chain
  
  # Choice of the initial position of the centroid and of the location of the centroid’s centre of attraction (C)
  timesteps <- nbObs + burnin
  X <- matrix(rep(NA, 2 * timesteps), ncol = 2)
  X[1, ] <- xy0
  C <- matrix(xy0, byrow = T, ncol = 2)
  
  # Generation of the state sequences for the different individuals
  Zind <- vector("list")
  for (k in 1:groupSize) {
    Zind[[k]] <- sample(1:nState, size = 1, prob = delta)
  }
  for (k in 2:timesteps) {
    for (j in 1:groupSize) {
      Zind[[j]][k] <- sample(1:nState, size = 1, prob = gamma[Zind[[j]][k - 1],])
    }
  }
  
  # Choice of the (random) initial positions of individuals 
  # Xind is the list that will comprise all individuals' simulated locations
  Xind <- vector("list")
  for (k in 1:groupSize) {
    Xind[[k]] <- matrix(rep(NA, 2 * timesteps), ncol = 2)
    Xind[[k]][1, ] <- c(rnorm(1,xy0[1], 1),rnorm(1,xy0[2], 1))
  }
  
  phi <- 0
  phiind <- rep(0, groupSize)
  
  for (k in 1:(timesteps - 1)) {
    ## start generation of centroid location at time k+1
    coo <- C[, ] - X[k, ] # xy diff btw centroid's AC and centroid's current location
    mu <- Arg(coo[1] + (0+1i) * coo[2]) # direction towards centroid's AC
    if (mu < 0)
      mu <- mu + 2 * pi # correct for direction
    # expected movement direction: weighted average of CRW and BRW [#Eq-1.2]
    mu.av <- Arg(etac * exp(phi * (0+1i)) + (1 - etac) * exp(mu * (0+1i)))
    phi <- rvm(1, mean = mu.av, k = kappac) # draw next direction [#Eq-1.1]
    if (k==1){
      step.len.ac <-  rgamma(1, shape = ac^2/bc^2, scale = bc^2/ac) # first step-length
    } else{
      # subsequent step-lengths are correlated to previous
      step.len.ac <- 0.5*(step.len.ac + rgamma(1, shape = ac^2/bc^2, scale = bc^2/ac)) 
    }
    step.ac <- step.len.ac * c(Re(exp((0+1i) * phi)), Im(exp((0+1i) * phi))) # dist in xy to move
    X[k + 1, ] <- X[k, ] + step.ac # next loc 
    
    ## start generation of individuals' locations at time k+1
    # first draw the mean step-length for the group at time t
    if (k==1){ # first step-length
      step.len <- rgamma(1, shape = aind[Zind[[j]][k]]^2/bind[Zind[[j]][k]]^2,
                         scale = bind[Zind[[j]][k]]^2/aind[Zind[[j]][k]])
    } else{
      step.len <- 0.5*( step.len + rgamma(1, 
                                          shape = aind[Zind[[j]][k]]^2/bind[Zind[[j]][k]]^2,
                                          scale = bind[Zind[[j]][k]]^2/aind[Zind[[j]][k]]) )
    }
    
    for (j in 1:groupSize) {
      if (s.ty[Zind[[j]][k]] == 1) { # state 1: BRW
        coo <- X[k + 1, ] - Xind[[j]][k, ] # xy diff btw centroid loc and indiv loc
        mu <- Arg(coo[1] + (0+1i) * coo[2])
        if (mu < 0)
          mu <- mu + 2 * pi
        phiind[j] <- rvm(1, mean = mu, k = kappaind[Zind[[j]][k]])
      }
      if (s.ty[Zind[[j]][k]] == 0) { # state 2: CRW
        mu <- rvm(1, mean = 0, k = kappaind[Zind[[j]][k]])
        phiind[j] <- phiind[j] + mu
      }
      
      # ind steps for time t are normally distributed around step.len
      repeat{
        step.len.ind <- rnorm(1, step.len, step.len/10)
        if (step.len.ind > 0){break}
      }
      
      step <- step.len * c(Re(exp((0+1i) * phiind[j])), Im(exp((0+1i) * phiind[j])))
      Xind[[j]][k + 1, ] <- Xind[[j]][k, ] + step
    }
  }
  # compile output
  Xind <- lapply(Xind, function(x) x[(burnin+1):timesteps,]) # discard burnin
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

"DETECTION SIMULATIONS"
# computes the polygon coordinates of a circle sector (credit: R package 'pgirmess')
polycirc2<-function (radius=1, center = c(0, 0), edges = 50, init=pi/2, angle=pi/2) 
{
  circ = NULL
  angles <- seq(init, init+angle, l = edges)
  for (i in angles) {
    circ <- rbind(circ, cbind(radius * sin(i), radius * cos(i)))
  }
  circ<-rbind(c(0,0),circ)
  x<-cbind(circ[, 1] + center[1], circ[,2] + center[2])
  x<-rbind(x,x[1,])
  return(x) 
  
}
# makes randomly facing circlular sectors (representing detection zone) from cam locs
makeTraps <- function(xlim = c(50,450), ylim = c(50,450), trapspacing = 50, 
                      r = 10, theta = pi/3, jitter=F) { # detection radius r and width w
  
  grid <- expand.grid(x=seq(xlim[1],xlim[2],trapspacing),y=seq(ylim[1],ylim[2],trapspacing))
  if(jitter==T){
    grid[, 1] <- grid[, 1] + runif(nrow(grid), -5, 5)
    grid[, 2] <- grid[, 2] + runif(nrow(grid), -5, 5)
  }
  
  grid$initDir <- runif(nrow(grid), 0, 2*pi) # pick random orientations for each trap
  
  #polycirc2(radius=r, center=x[1:2], init=x[3],angle=theta)
  camPolys <- apply(grid, 1, function(x) Polygon( polycirc2(radius=r, center=x[1:2], init=x[3],angle=theta)) )
  camStore <- list()
  for (i in 1:length(camPolys)){
    camStore[[i]] <- Polygons(list(camPolys[[i]]), paste(i))
  }
  camPolys.sp <- SpatialPolygons(camStore)
  
  return(list(trapPoly = camPolys.sp, trapInfo = grid))
}

# extracts track overlaps with detectors
# DetectionGenerator will check every single trajectory against every trap; extremely intensive when nbObs is large
# DetectionGeneratorLite will only check the trajectories buffering each trap for detections; bufferDist needs to be the max steplength to avoid missing dets
DetectionGenerator <- function(movedf=NA, traps = traps, torus=T, timer = T){
  movedf <- movedf[,c('ID','x','y','time')]
  movedf$overlap <- rep(NA, nrow(movedf))
  xline <- list()
  X <- traps[['trapPoly']]
  
  if (timer == T){  # set timer for the following loop
   require(pbapply)
   pb <- timerProgressBar(width = 35, char = "+", style = 1)
   on.exit(close(pb))
  }
  # convert locs to indiv lines
  for (i in 1:nrow(movedf)){
    x <- rbind( movedf[i,c('x','y')], movedf[i-1,c('x','y')])
    xline[[i]] <- Lines(list(Line(coordinates(x))), ID=paste0(i))
    if (timer == T){ # progress bar
    setTimerProgressBar(pb, i/nrow(movedf))}
  }
  cat('\nLocs to lines conversion finished.\nPlease wait for line-to-detetector overlap calculation to finish :o)')
  df.lines <- SpatialLinesDataFrame(SpatialLines(xline, CRS(as.character(NA))), 
                                    movedf, match.ID = FALSE)
  if(torus == T){    
    movedf %>% 
      add_rownames() %>%
      group_by(ID) %>% 
      mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
      filter(row_number()==1 | step > 950) %>% 
      `[[`("rowname") %>%
      as.numeric() -> index # create index for first line from every id OR lines that traverse landscape bounds
  } else{
    movedf %>% 
      add_rownames() %>%
      group_by(ID) %>% 
      filter(row_number()==1) %>%
      `[[`("rowname") %>%
      as.numeric() -> index # create index for first line from every id
    }
  
  for(k in 1:length(index)){  # remove first line from every id
    df.lines@lines[[index[k]]]@Lines[[1]] <- NULL
  }
  
  overlap <- over(df.lines, X, returnList = F)
  df.lines@data[["overlap"]] <- overlap
  overlap.df <- df.lines[!is.na(df.lines@data[["overlap"]]),]
  x <- as.data.frame(overlap.df)
  names(x)[names(x) == 'overlap'] <- 'Site'
  return(x)
}
DetectionGeneratorLite <- function(movedf=NA, traps = traps, bufferDist = 10, torus=T){
  require(data.table)
  movedf <- as.data.table(movedf[,c('ID','x','y','time')]) # needs to be data.table
  movedf$overlap <- rep(NA, nrow(movedf))
  X <- traps[['trapPoly']]
  
  # evaluate by camera
  for (j in 1:nrow(traps[['trapInfo']])){
    xline <- list()
    # subset animal locs within bufferDist away from cam j
    limits_x <- c(traps[['trapInfo']][j,'x'] - bufferDist, traps[['trapInfo']][j,'x'] + bufferDist)
    limits_y <- c(traps[['trapInfo']][j,'y'] - bufferDist, traps[['trapInfo']][j,'y'] + bufferDist)
    subpts_x <- movedf[x >= limits_x[1] & x < limits_x[2]] %>%
      setkey(y)
    subpts_comp <- subpts_x[y >= limits_y[1] & y < limits_y[2]]
    rm(subpts_x)
    if(nrow(subpts_comp)==0){cat(paste('\ncam',j,'done')); next}
    
    # for(id in unique(subpts_comp$ID)){
    #   lines_in_buffer0 <- movedf %>%
    #     filter(ID == id) %>%
    #     filter( time %in% (min(subpts_comp$time[subpts_comp$ID == id])-1):(max(subpts_comp$time[subpts_comp$ID == id])+1) ) # this needs to match the ID ffrrr
    #   
    #   if(exists("lines_in_buffer") == T){
    #     lines_in_buffer <- rbind(lines_in_buffer0,lines_in_buffer)
    #   } else{
    #     lines_in_buffer <- lines_in_buffer0
    #   }
    #   rm(lines_in_buffer0)
    # }
    
    lines_in_buffer <- subpts_comp %>%
      group_by(ID) %>%
      arrange(time) %>%
      mutate(difftime = time - lag(time, 1, NA)) %>%
      arrange(ID,time)
    
    # turn animal locs into lines
    for (i in 1:nrow(lines_in_buffer)){
      xline0 <- rbind( lines_in_buffer[i,c('x','y')], lines_in_buffer[i-1,c('x','y')])
      xline[[i]] <- Lines(list(Line(coordinates(xline0))), ID=paste0(i))
    }
    df.lines <- SpatialLinesDataFrame(SpatialLines(xline, CRS(as.character(NA))), 
                                      lines_in_buffer, match.ID = FALSE)
    
    # correct for torus
    if(torus == T){    
      lines_in_buffer %>% 
        add_rownames() %>%
        group_by(ID) %>% 
        mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
        #filter(row_number()==1 | step > 950) %>% 
        filter(is.na(difftime) | !difftime==1 | step > 5000) %>%
        `[[`("rowname") %>%
        as.numeric() -> index # create index for first line from every id OR lines that traverse landscape bounds
    } else{
      lines_in_buffer %>% 
        add_rownames() %>%
        group_by(ID) %>% 
        #filter(row_number()==1) %>%
        filter(is.na(difftime) | !difftime==1) %>%
        `[[`("rowname") %>%
        as.numeric() -> index # create index for first line from every id
    }
    
    # correct lines
    for(k in 1:length(index)){  # remove first line from every id OR lines that traverse landscape bounds (for torus=T)
      df.lines@lines[[index[k]]]@Lines[[1]] <- NULL
    }
    
    # check for overlap
    overlap <- over(df.lines, X[X@polygons[[j]]@ID==j], returnList = F)
    if (all(is.na(overlap))) { cat(paste('\ncam',j,'done')); next }
    df.lines@data[["overlap"]] <- overlap
    overlap.df <- df.lines[!is.na(df.lines@data[["overlap"]]),]
    x <- as.data.frame(overlap.df)
    names(x)[names(x) == 'overlap'] <- 'Site'
    
    # record data
    if ( exists("det_output")==T ){
      det_output <- rbind(det_output, x)
    } else {
      det_output <- x
    }
    
    cat(paste('\ncam',j,'done'))
  } # j
  
  return(det_output)
}

# applies detection decay function with distance (requires detections to be treated as 'photos'/discrete points)
# d: vector of distances to CT 
# perfect_dist: dist to which detection is perfect
imperfectDet <- function(d = d, perfect_dist = 3, max_dist = 10, burst = T, bursts = 2)
{
  detType <- isDetected <- vector()
  k = 0
  for (t in 1:length(d)){
    k = k + 1 # 
    
    # simulated detection function (https://doi.org/10.1111/j.2041-210X.2011.00094.x)
    # p <- (1-exp(-(3/d[k])^6))/(1+exp(10*(1-d[k])))
    p <- 1-exp(-(3/d[k])^6)
    
    isDetected[k] <- rbinom(1, 1, p) # 1: detection, 0: no detection
    
    if (burst == T) {
      if (isDetected[k] == 1){ # if detection
        isDetected[k+1:bursts] <- 1 # next two intervals recorded as 'bursts'
        k = k + bursts # skip the two bursts
      }
      if(k >= length(d)){ # in case k index exceeds length of d
        break
      }
    }

  }
  #return( data.frame(dist = as.vector(d), isDetected = isDetected[1:length(d)]) )
  return( isDetected[1:length(d)] )
}

# converts DetectionGenerator output to 'photos'/discrete points and applies imperfectDet decay function
# timestep of simulated movement (mins)
# timeInt to sample movement trajectories at (secs), mimicking photo captures of CTs
ImperfectDetConverter <- function(movedf = moveSims, detdf = detDat, traps = traps, timestep = 5, timeInt = 1, 
                                  perfect_dist = 3, max_dist = 10)
{
  dist_to_CT <- list()
  X <- traps[["trapPoly"]]
  polyID <- sapply(slot(X, "polygons"), function(x) slot(x, "ID"))
  
  for (i in 1:nrow(detdf)){ 
    # for every detection, find the track from the movement dataset
    t <- detdf[i,'time']; id <- detdf[i,'ID']; detector <- detdf[i,'Site']
    track <- movedf %>%
      filter(ID == id & time %in% (t-1):t ) %>%  as.data.frame()
    coords <- coordinates(track[,c('x','y')])
    trackLine <- SpatialLines(list(Lines(list( Line(coords)), "id" )))
    detectorID <- which(polyID == detector)
    
    # measure the track intersecting with the detection zone, then split into temporal segments
    trackIntersect <- gIntersection(X[detectorID,],trackLine) # portion of track intersecting with detection zone
    timeInZone <- gLength(trackIntersect)/gLength(trackLine) * timestep * 60 # time in detection zone (seconds)
    timeSegments <- seq(0,timeInZone,timeInt) # get time intervals to record distances at (every timeInt)
    distSegments <- timeSegments/timeInZone
    
    # calculate distance to CT
    trackIntersectxy <- trackIntersect@lines[[1]]@Lines[[1]]@coords
    x3 = distSegments * trackIntersectxy[2,'x'] + (1 - distSegments) * trackIntersectxy[1,'x'] #find point on each segment
    y3 = distSegments * trackIntersectxy[2,'y']  + (1 - distSegments) * trackIntersectxy[1,'y'] 
    distLocs <- data.frame(x=x3,y=y3)
    #ct_midpt <- X[detectorID,]@polygons[[1]]@labpt # this is the midpoint of the CT polygon (i.e., detection zone)
    #ct_coord <- matrix(round(ct_midpt/50)*50,1,2) # get the actual CT location
    ct_coord <- traps[["trapInfo"]][detectorID, c('x','y')]
    dist <- as.vector(e2dist(ct_coord, distLocs)) # euclidean distance between CT and points on line
    
    # calculate time at FOV entry
    entryTime <- as.vector(e2dist(matrix(coords[1,],1,2), matrix(trackIntersectxy[1,],1,2))) / gLength(trackLine) * timestep * 60
    
    dist_to_CT[[i]] <- data.frame(track = i, x = x3, y = y3, time = t, dist = dist, entryTime = entryTime) # each list contains the distance to CT for every x sec interval
  }
  
  dist_to_CT2 <- lapply(dist_to_CT, function(d){
    d$isDetected <- imperfectDet(d = d$dist,perfect_dist = perfect_dist, max_dist = max_dist, burst = T, bursts = 2)
    d$time_seconds <- d$time*timestep*60 + d$entryTime + (1:nrow(d)-1)*timeInt # high resolution time in seconds;  (1:nrow(d)-1)
    return(d)
  } )
  
  dupTimes <- unlist(lapply(dist_to_CT2, nrow))
  detdf_id <- rep(1:nrow(detdf), dupTimes)
  detdf_new <- detdf[detdf_id,]
  #detdf_new$distance <- do.call(c,dist_to_CT)
  detdf_new <- cbind(detdf_new[c('ID','Site')], do.call(rbind,dist_to_CT2))
  
  return(detdf_new)
}

# this version splits *entire* track into equally spaced points (acc. to time interval), and samples points overlapping the FOV
DetSamplPoints <- function(movedf = moveSims, detdf = detDat, traps = traps, timestep = 5, timeInt = 1, 
                                   perfect_dist = 3, max_dist = 10)
{
  dist_to_CT <- list()
  X <- traps[["trapPoly"]]
  polyID <- sapply(slot(X, "polygons"), function(x) slot(x, "ID"))
  
  for (i in 1:nrow(detdf)){ 
    # for every detection, find the track from the movement dataset
    t <- detdf[i,'time']; id <- detdf[i,'ID']; detector <- detdf[i,'Site']
    track <- movedf %>%
      filter(ID == id & time %in% (t-1):t ) %>%  as.data.frame()
    coords <- coordinates(track[,c('x','y')])
    # trackLine <- SpatialLines(list(Lines(list( Line(coords)), "id" )))
    detectorID <- which(polyID == detector)
    
    # split track  into temporal segments
    timeInSec <- timestep * 60 # time in seconds
    distSegments <- seq(0,timeInSec,timeInt) / timeInSec # get time intervals to record distances at (every timeInt) as proportion of entire track
    x3 = distSegments * track[2,'x'] + (1 - distSegments) * track[1,'x'] #find point on each segment
    y3 = distSegments * track[2,'y']  + (1 - distSegments) * track[1,'y'] 
    pointsOnTrack <- SpatialPoints(coordinates(data.frame(x=x3,y=y3)[-1,]))
    
    # find points on track overlapped with FOV
    trackIntersect <- gIntersection(X[detectorID,], pointsOnTrack, byid = TRUE) # portion of track points intersecting with detection zone
    possibleError <- tryCatch(timestamp <- (as.numeric(data.frame(do.call('rbind', strsplit(row.names(trackIntersect@coords),' ',fixed=TRUE)))[,2]) - 1) * timeInt,
                              error = function(e) e)
    if(inherits(possibleError, "error")) next
    
    # calculate distance to CT
    distLocs <- trackIntersect@coords
    ct_coord <- traps[["trapInfo"]][detectorID, c('x','y')]
    dist <- as.vector(e2dist(ct_coord, distLocs)) # euclidean distance between CT and points on line
    
    # time at FOV entry = time of first overlapped point
    entryTime <- timestamp[1]
    
    dist_to_CT[[i]] <- data.frame(ID = id, Site = detector, track = i, x = distLocs[,'x'], y = distLocs[,'y'], time = t, time_seconds =  t*timestep*60 + timestamp, dist = dist, entryTime = entryTime) # each list contains the distance to CT for every x sec interval
  }
  
  dist_to_CT <- dist_to_CT[which(!sapply(dist_to_CT, is.null))]
  # dist_to_CT2 <- lapply(dist_to_CT, function(d){
  #   d$isDetected <- imperfectDet(d = d$dist,perfect_dist = perfect_dist, max_dist = max_dist, burst = T, bursts = 7)
  #   return(d)
  # } )
  
  # dupTimes <- unlist(lapply(dist_to_CT2, nrow))
  # detdf_id <- rep(1:length(dist_to_CT2), dupTimes)
  # detdf_new <- detdf[detdf_id,]
  # detdf_new <- cbind(detdf_new[c('ID','Site')], do.call(rbind,dist_to_CT2))
  detdf_new <- do.call(rbind,dist_to_CT)
  row.names(detdf_new) <- 1:nrow(detdf_new)
  
  return(detdf_new)
}

# applies detection decay function with distance (requires detections to be treated as 'photos'/discrete points)
# d: vector of distances to CT 
# perfect_dist: dist to which detection is perfect
imperfectDet <- function(d = d, perfect_dist = 3, max_dist = 10, burst = T, bursts = 2)
{
  detType <- isDetected <- vector()
  k = 0
  for (t in 1:length(d)){
    k = k + 1 # 
    
    # simulated detection function (https://doi.org/10.1111/j.2041-210X.2011.00094.x)
    # p <- (1-exp(-(3/d[k])^6))/(1+exp(10*(1-d[k])))
    p <- 1-exp(-(3/d[k])^6)
    
    isDetected[k] <- rbinom(1, 1, p) # 1: detection, 0: no detection
    
    if (burst == T) {
      if (isDetected[k] == 1){ # if detection
        isDetected[k+1:bursts] <- 1 # next two intervals recorded as 'bursts'
        k = k + bursts # skip the two bursts
      }
      if(k >= length(d)){ # in case k index exceeds length of d
        break
      }
    }
    
  }
  #return( data.frame(dist = as.vector(d), isDetected = isDetected[1:length(d)]) )
  return( isDetected[1:length(d)] )
}


"DETECTION DATA PREP"
# # for REST
# calcStay <- function(movedf = moveSims, detdf = x1, X = traps.buff, timestep = 15, nbObs = nbObs){
#   detdf$stayTime <- rep(NA, nrow(detdf))
#   movedf$time <- rep(1:nbObs, length(unique(movedf$ID)))
#   
#   polyID <- sapply(slot(X, "polygons"), function(x) slot(x, "ID"))
#   
#   for (i in 1: nrow(detdf)){
#     t <- detdf[i,'time']; id <- detdf[i,'ID']; detector <- detdf[i,'Site']
#     
#     track <- movedf %>%
#       filter(ID == id & time %in% (t-1):t ) %>%  as.data.frame()
#     coords <- coordinates(track[,c('x','y')])
#     trackLine <- SpatialLines(list(Lines(list( Line(coords)), "id" )))
#     
#     detectorID <- which(polyID == detector)
#     
#     trackIntersect <- gIntersection(X[detectorID,],trackLine)
#     detdf$stayTime[i] <- timestep * 60 * gLength(trackIntersect)/gLength(trackLine) # in seconds (* 60 cuz timestep is in mins)
#   }
#   return(detdf)
# }
# 
# # for REM
# calcSpeed <- function(movedf = moveSims, detdf = x1, X = traps.buff, timestep = 15, nbObs = nbObs){
#   detdf$dist <- detdf$stayTime <- rep(NA, nrow(detdf))
#   movedf$time <- rep(1:nbObs, length(unique(movedf$ID)))
#   
#   polyID <- sapply(slot(X, "polygons"), function(x) slot(x, "ID"))
#   
#   for (i in 1: nrow(detdf)){
#     t <- detdf[i,'time']; id <- detdf[i,'ID']; detector <- detdf[i,'Site']
#     
#     track <- movedf %>%
#       filter(ID == id & time %in% (t-1):t ) %>%  as.data.frame()
#     coords <- coordinates(track[,c('x','y')])
#     trackLine <- SpatialLines(list(Lines(list( Line(coords)), "id" )))
#     
#     detectorID <- which(polyID == detector)
#     
#     # while speed can directly be calculated from gLength(trackLine) since time units are constant in sims,
#     # this matters for aggregating consecutive encounters which means more than one straight-line track
#     trackIntersect <- gIntersection(X[detectorID,],trackLine) # this calculates distance
#     detdf$dist[i] <- gLength(trackIntersect)
#     detdf$stayTime[i] <- timestep * 60 * gLength(trackIntersect)/gLength(trackLine) # in seconds (* 60 cuz timestep is in mins)
#     
#   }
#   return(detdf)
# }
# 
# # for CT-DS
# calcDist <- function(movedf = moveSims, detdf = x1, X = traps.buff, timestep = 15, nbObs = nbObs, timeInt = 2){
#   dist_to_CT <- list()
#   movedf$time <- rep(1:nbObs, length(unique(movedf$ID)))
#   
#   polyID <- sapply(slot(X, "polygons"), function(x) slot(x, "ID"))
#   
#   for (i in 1:nrow(detdf)){
#     t <- detdf[i,'time']; id <- detdf[i,'ID']; detector <- detdf[i,'Site']
#     
#     track <- movedf %>%
#       filter(ID == id & time %in% (t-1):t ) %>%  as.data.frame()
#     coords <- coordinates(track[,c('x','y')])
#     trackLine <- SpatialLines(list(Lines(list( Line(coords)), "id" )))
#     
#     detectorID <- which(polyID == detector)
#     trackIntersect <- gIntersection(X[detectorID,],trackLine) # portion of track intersecting with detection zone
#     timeInZone <- gLength(trackIntersect)/gLength(trackLine) * timestep * 60 # time in detection zone (seconds)
#     timeSegments <- seq(0,timeInZone,timeInt) # get time intervals to record distances at (every timeInt)
#     distSegments <- timeSegments/timeInZone
#     
#     trackIntersectxy <- trackIntersect@lines[[1]]@Lines[[1]]@coords
#     x3 = distSegments * trackIntersectxy[2,'x'] + (1 - distSegments) * trackIntersectxy[1,'x'] #find point on each segment
#     y3 = distSegments * trackIntersectxy[2,'y']  + (1 - distSegments) * trackIntersectxy[1,'y'] 
#     distLocs <- data.frame(x=x3,y=y3)
#     ct_midpt <- X[detectorID,]@polygons[[1]]@labpt # this is the midpoint of the CT polygon (i.e., detection zone)
#     ct_coord <- matrix(round(ct_midpt/50)*50,1,2) # get the actual CT location
#     dist_to_CT[[i]] <- e2dist(ct_coord, distLocs)
#   }
#   
#   dupTimes <- unlist(lapply(dist_to_CT, length))
#   detdf_id <- rep(1:nrow(detdf), dupTimes)
#   detdf_new <- detdf[detdf_id,]
#   detdf_new$distance <- do.call(c,dist_to_CT)
#   
#   return(detdf_new)
# }


# euclidean distance btw two points
e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

# minimum distance to line
dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  return(d)
} 

"PLOTTING"
plotTraject <- function(df = df, plotit = X, torus=FALSE){
  require(scales)
  k = 0
  n_foragers = length(unique(df$ID))
  if (torus==FALSE){
    for(i in sort(unique(df$ID))){ # lines
      k=k+1
      if(plotit[k] == 1) { # or i
        lines(x=df$x[df$ID==i],
              y=df$y[df$ID==i],
              lwd = 1.5,
              col = alpha(tol.rainbow(n_foragers)[k],0.25))}
    } 
  }
  if (torus==TRUE){
    
    rowIDs <- df %>%
      group_by(ID) %>%
      mutate( row_num = 1:length(ID) ) %>%
      mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
      filter( step > 950 ) %>%
      as.data.frame()  
    
    for(i in sort(unique(df$ID))){ # lines
      k=k+1
      
      if(i %in% rowIDs$ID){
        breakpts <- c(1,rowIDs$row_num[rowIDs$ID==i])
        for (ls in 2:length(breakpts)){
          start = breakpts[ls-1]; end = breakpts[ls] - 1
          lines(x=df$x[df$ID==i][start:end],
                y=df$y[df$ID==i][start:end],
                lwd = 1.5,
                col = alpha(tol.rainbow(n_foragers)[k],0.25))
        }
      } else{
        lines(x=df$x[df$ID==i],
              y=df$y[df$ID==i],
              lwd = 1.5,
              col = alpha(tol.rainbow(n_foragers)[k],0.25))
      }
      
    }
    
  }
}

