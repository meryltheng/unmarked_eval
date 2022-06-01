# ====================================
# Unmarked abundance evalaution
# Functions
# ====================================

"MOVEMENT SIMULATIONS"
# Generate activity centres (ACs)
generateAC <- function(N = 5, xlim = c(0,500), ylim = c(0,500), buffer = 50){
  ACs <- matrix(NA, N, 2)
  ACs[1,] <- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
  for (i in 2:N){
    ACcand <- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
    distapart <- apply(ACs, 1, function(x) sqrt( (x[1] - ACcand[1])^2 + (x[2] - ACcand[2])^2) )
    while( any(distapart <= buffer, na.rm = T) ){
      ACcand <- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
      distapart <- apply(ACs, 1, function(x) sqrt( (x[1] - ACcand[1])^2 + (x[2] - ACcand[2])^2) )
    }
    ACs[i,] <- ACcand 
  }
  
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
        filter(is.na(difftime) | !difftime==1 | step > 950) %>%
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
    for(k in 1:length(index)){  # remove first line from every id
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
imperfectDet <- function(d = d, perfect_dist = 0.3, max_dist = 1)
{
  detType <- isDetected <- vector()
  k = 0
  for (t in 1:length(d)){
    k = k + 1 # 
    
    # simulated detection function
    p <- (1-exp(-(0.3/d[k])^6))/(1+exp(100*(0.1-d[k]))) # https://doi.org/10.1111/j.2041-210X.2011.00094.x
    
    isDetected[k] <- rbinom(1, 1, p) # 1: detection, 0: no detection
    if (isDetected[k] == 1){ # if detection
      isDetected[c(k+1,k+2)] <- 1 # next two intervals recorded as 'bursts'
      k = k + 2 # skip the two bursts
    }
    if(k >= length(d)){ # in case k index exceeds length of d
      break
    }
  }
  #return( data.frame(dist = as.vector(d), isDetected = isDetected[1:length(d)]) )
  return( isDetected[1:length(d)] )
}

# converts DetectionGenerator output to 'photos'/discrete points and applies imperfectDet decay function
# timestep of simulated movement (mins)
# timeInt to sample movement trajectories at (secs), mimicking photo captures of CTs
ImperfectDetConverter <- function(movedf = moveSims, detdf = detDat, traps = traps, timestep = 5, timeInt = 1, 
                                  perfect_dist = 0.3, max_dist = 1)
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
    ct_midpt <- X[detectorID,]@polygons[[1]]@labpt # this is the midpoint of the CT polygon (i.e., detection zone)
    ct_coord <- matrix(round(ct_midpt/50)*50,1,2) # get the actual CT location
    dist <- as.vector(e2dist(ct_coord, distLocs))
    
    # calculate time at FOV entry
    entryTime <- as.vector(e2dist(matrix(coords[1,],1,2), matrix(trackIntersectxy[1,],1,2))) / gLength(trackLine) * timestep * 60
    
    dist_to_CT[[i]] <- data.frame(track = i, x = x3, y = y3, time = t, dist = dist, entryTime = entryTime) # each list contains the distance to CT for every x sec interval
  }
  
  dist_to_CT2 <- lapply(dist_to_CT, function(d){
    d$isDetected <- imperfectDet(d = d$dist,perfect_dist = perfect_dist, max_dist = max_dist)
    d$time_seconds <- d$time*timestep*60 + d$entryTime + (1:nrow(d)-1)*timeInt # high resolution time in seconds
    return(d)
  } )
  
  dupTimes <- unlist(lapply(dist_to_CT2, nrow))
  detdf_id <- rep(1:nrow(detdf), dupTimes)
  detdf_new <- detdf[detdf_id,]
  #detdf_new$distance <- do.call(c,dist_to_CT)
  detdf_new <- cbind(detdf_new[c('ID','Site')], do.call(rbind,dist_to_CT2))
  
  return(detdf_new)
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
  k = 0
  n_foragers = length(unique(df$ID))
  if (torus==FALSE){
    for(i in sort(unique(df$ID))){ # lines
      k=k+1
      if(plotit[k] == 1) { # or i
        lines(x=df$x[df$ID==i],
              y=df$y[df$ID==i],
              lwd = 1.5,
              col = alpha(rainbow(n_foragers)[k],0.25))}
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
                col = alpha(rainbow(n_foragers)[k],0.25))
        }
      } else{
        lines(x=df$x[df$ID==i],
              y=df$y[df$ID==i],
              lwd = 1.5,
              col = alpha(rainbow(n_foragers)[k],0.25))
      }
      
  }

  }
}

plotTrajectTemp <- function(df = df, plotit = X){
  require(viridis)
  k = 0
  n_foragers = length(unique(df$ID))
  for(i in sort(unique(df$ID))){ # lines
    k=k+1
    if(plotit[k] == 1) { # or i
      segments(x0 = head(df$x[df$ID==i],-1),
               y0 = head(df$y[df$ID==i],-1),
               x1 = tail(df$x[df$ID==i],-1),
               y1 = tail(df$y[df$ID==i],-1),
               lwd = 1,
               col = alpha(viridis(length(df$x[df$ID==i])),0.5))
    }
  } 
}

