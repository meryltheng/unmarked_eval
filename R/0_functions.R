# ====================================
# Unmarked abundance evalaution
# Functions
# ====================================

#### CREATE STUDY AREA COVARIATE (code from B. McClintock) ####
genAreaRaster <- function(ncols,nrows,xlim=c(0,500),ylim=c(0,500)){
  potSurface <- raster(ncol=(nrows+2)*2,nrow=(ncols+2)*2,xmn=xlim[1]-150,xmx=xlim[2]+150,ymn=ylim[1]-150,ymx=ylim[2]+150, vals=rep(0,(nrows+2)*2*(ncols+2)*2))
  potSurface<-setValues(potSurface,values=1-c(rep(rep(0,(ncols+2)*2),(nrows+2)/2+1),rep(c(rep(0,(ncols+2)/2+1),rep(1,ncols),rep(0,(ncols+2)/2+1)),nrows),rep(rep(0,(ncols+2)*2),(nrows+2)/2+1)))
  potSurface[potSurface>0] <- NA
  potSurface <- raster::distance(potSurface)
  crs(potSurface) <-"+proj=utm +zone=53H +south +datum=WGS84"
  tmp <- potSurface / mean(raster::values(raster::terrain(potSurface, opt = "slope")), na.rm = T) #standardize based on slope of gradient
  potSurfaceRast <- ctmcmove::rast.grad(tmp)[c("rast.grad.x","rast.grad.y")]
  #potSurfaceRast <- lapply(potSurfaceRast,function(x) spatial.tools::modify_raster_margins(x,extent_delta = c(-1,-1,-1,-1))) # remove outer cells with NA gradient
  mySurface <- list("potSurfaceRast"=potSurfaceRast,"potSurface"=potSurface)
  
  return(mySurface)
}

#### CREATE SEPARATE RASTER FOR EACH SNAKES ACTIVITY CENTER
rastActivity <- function(centers,bait){
  centerStack <- centerStack.x <- centerStack.y <- stack()
  bait <- setValues(bait,NA)
  for(j in 1:nrow(centers)){
    k2 <- bait
    k2[cellFromXY(k2,centers[j,])] <- 1
    k2 <- raster::distance(k2)
    crs(k2) <-"+proj=utm +zone=55 +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs"
    centerStack <- stack(centerStack,k2)
    k2 <- k2 / mean(raster::values(terrain(k2, opt = "slope")), na.rm = T) #standardize based on slope of gradient
    centerGrad <- rast.grad(k2)[c("rast.grad.x","rast.grad.y")]
    centerStack.x <- stack(centerStack.x, centerGrad$rast.grad.x) ## This gets warnings about ncolumns (though same number) but the rasters seem to stack correctly
    centerStack.y <- stack(centerStack.y, centerGrad$rast.grad.y) ## This gets warnings about ncolumns (though same number) but the rasters seem to stack correctly
    rm(centerGrad)
  }
  names(centerStack) <- names(centerStack.x) <- names(centerStack.y) <- as.character(1:nrow(centers))
  #plot(centerStack)
  mycenter <- list("centerStack.x"=centerStack.x,"centerStack.y"=centerStack.y)
  
  return(mycenter)
}

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

generateAChab <- function(N = 5, xlim = c(0,500), ylim = c(0,500), buffer = 50){
  ACs <- matrix(NA, N, 2)
  repeat{
    cellID <- sample(which.max(hab==2), 1, replace = F)
    ACcand <- xyFromCell(hab, cellID)
    
    if( (ACcand[1,1] >= xlim[1] & ACcand[1,1] <= xlim[2] &
         ACcand[1,2] >= ylim[1] & ACcand[1,2] <= ylim[2])  == TRUE ){
      ACs[1,] <- ACcand[1,]
      break }
  }
  
  for (i in 2:N)
  {  repeat{
    cellID <- sample(which.max(hab==2), 1, replace = F)
    ACcand <- xyFromCell(hab, cellID)
    distapart <- apply(ACs, 1, function(x) sqrt( (x[1] - ACcand[1])^2 + (x[2] - ACcand[2])^2) )
    
    if( (ACcand[1,1] >= xlim[1] & ACcand[1,1] <= xlim[2] &
         ACcand[1,2] >= ylim[1] & ACcand[1,2] <= ylim[2] & 
         all(distapart > buffer, na.rm = T))  == TRUE ){
      ACs[i,] <- ACcand[1,]
      break }
  }}
  return(ACs)
}
# Create raster for each AC (code adapted from Staci Amburgey)


"DETECTION SIMULATIONS"
# new function that makes randomly facing triangular detectors from cam locs
makeTraps <- function(xlim = c(50,450), ylim = c(50,450), trapspacing = 50, 
                      r = 10, w = 10/2, jitter=F) { # detection radius r and width w
  
  grid <- expand.grid(x=seq(xlim[1],xlim[2],trapspacing),y=seq(ylim[1],ylim[2],trapspacing))
  if(jitter==T){
    grid[, 1] <- grid[, 1] + runif(nrow(grid), -5, 5)
    grid[, 2] <- grid[, 2] + runif(nrow(grid), -5, 5)
  }
  
  grid$turnAng <- runif(nrow(grid), -pi, pi) # pick random orientations for each trap
  
  camPolys <- apply(grid, 1, function(x) Polygon( cbind(c(x[1], # cam x coord
                                                          x[1]+cos(x[3]+tanh(w/r))*sqrt(w^2+r^2), # add detection angle (calculated from w and r) to cam direction
                                                          x[1]+cos(x[3]-tanh(w/r))*sqrt(w^2+r^2)), # minus detection angle from cam direction
                                                        c(x[2], # cam y coord
                                                          x[2]+sin(x[3]+tanh(w/r))*sqrt(w^2+r^2),
                                                          x[2]+sin(x[3]-tanh(w/r))*sqrt(w^2+r^2)) )) )
  camStore <- list()
  for (i in 1:length(camPolys)){
    camStore[[i]] <- Polygons(list(camPolys[[i]]), paste(i))
  }
  camPolys.sp <- SpatialPolygons(camStore)
  
  return(camPolys.sp)
}

makeTrapsCluster <- function(xlim = c(75,425), ylim = c(75,425), clusterspacing = 70, trapspacing = 15, 
                             n = 9, nClus = 4, r = 1, w = 1/2) { # no. of trap clusters n, detection radius r and width w
  
  if( !nClus %in% c(2,3,4)){
    stop('Error: nClus must be 2, 3 or 4')
  }
  
  if(nClus == 2){
    powa <- 1
  }
  
  if(nClus == 3){
    powa <- 0
  }
  
  if(nClus == 4){
    powa <- -1
  }
  
  grid <- expand.grid(x=seq(xlim[1],xlim[2],clusterspacing),y=seq(ylim[1],ylim[2],clusterspacing))
  grid <- grid[sample(1:nrow(grid), n),]
  traps <- apply(grid, 1, function(x) expand.grid(x=seq(x[1]-trapspacing/(2^powa),x[1]+trapspacing/(2^powa),trapspacing),y=seq(x[2]-trapspacing/(2^powa),x[2]+trapspacing/(2^powa),trapspacing)) )
  traps <- do.call(rbind,traps)
  
  traps$turnAng <- runif(nrow(traps), -pi, pi) # pick random orientations for each trap
  
  camPolys <- apply(traps, 1, function(x) Polygon( cbind(c(x[1], # cam x coord
                                                          x[1]+cos(x[3]+tanh(w/r))*sqrt(w^2+r^2), # add detection angle (calculated from w and r) to cam direction
                                                          x[1]+cos(x[3]-tanh(w/r))*sqrt(w^2+r^2)), # minus detection angle from cam direction
                                                        c(x[2], # cam y coord
                                                          x[2]+sin(x[3]+tanh(w/r))*sqrt(w^2+r^2),
                                                          x[2]+sin(x[3]-tanh(w/r))*sqrt(w^2+r^2)) )) )
  camStore <- list()
  for (i in 1:length(camPolys)){
    camStore[[i]] <- Polygons(list(camPolys[[i]]), paste(i))
  }
  camPolys.sp <- SpatialPolygons(camStore)
  
  return(camPolys.sp)
}

# old function that makes circular detectors
makeTrapsCirc <- function(xlim = c(60,440), ylim = c(60,440), trapspacing = 20, bufferdist = 0.2) { # 10m if 1 unit = 50m
  grid <- expand.grid(x=seq(xlim[1],xlim[2],trapspacing),y=seq(ylim[1],ylim[2],trapspacing)) 
  coordinates(grid) <- ~x+y
  grid <- SpatialPoints(grid)
  traps.buff <- gBuffer(grid, width = bufferdist ,quadsegs=round(1/10,0))
  traps.buff <- disaggregate(traps.buff)
  return(traps.buff)
}

# extracts detections from traps
DetectionGenerator <- function(df=NA, X = traps.buff, torus=T, timer = T){
  df <- df[,c('ID','x','y','time')]
  df$overlap <- rep(NA, nrow(df))
  xline <- list()
  
  if (timer == T){  # set timer for the following loop
   require(pbapply)
   pb <- timerProgressBar(width = 35, char = "+", style = 1)
   on.exit(close(pb))
  }
  # convert locs to indiv lines
  for (i in 1:nrow(df)){
    x <- rbind( df[i,c('x','y')], df[i-1,c('x','y')])
    xline[[i]] <- Lines(list(Line(coordinates(x))), ID=paste0(i))
    if (timer == T){ # progress bar
    setTimerProgressBar(pb, i/nrow(df))}
  }
  cat('\nLocs to lines conversion finished.\nNow, please patiently wait for line-to-detetector overlap calculation to finish :o)')
  df.lines <- SpatialLinesDataFrame(SpatialLines(xline, CRS(as.character(NA))), 
                                    df, match.ID = FALSE)
  if(torus == T){    
    df %>% 
      add_rownames() %>%
      group_by(ID) %>% 
      mutate( step = lag( sqrt( (lead(x,n=1)-x)^2 + (lead(y,n=1)-y)^2 ), n=1 ) ) %>%
      filter(row_number()==1 | step > 950) %>% 
      `[[`("rowname") %>%
      as.numeric() -> index # create index for first line from every id OR lines that traverse landscape bounds
  } else{
    df %>% 
      add_rownames() %>%
      group_by(ID) %>% 
      filter(row_number()==1) %>%
      `[[`("rowname") %>%
      as.numeric() -> index # create index for first line from every id
    }
  
  for(k in 1:length(index)){  # remove first line from every id
    df.lines@lines[[index[k]]]@Lines[[1]] <- NULL
  }
  
  overlap <- over(df.lines, traps.buff, returnList = F)
  df.lines@data[["overlap"]] <- overlap
  overlap.df <- df.lines[!is.na(df.lines@data[["overlap"]]),]
  x <- as.data.frame(overlap.df)
  names(x)[names(x) == 'overlap'] <- 'Site'
  return(x)
}


"DETECTION DATA PREP"
# for REST
calcStay <- function(movedf = moveSims, detdf = x1, X = traps.buff, timestep = 15, nbObs = nbObs){
  detdf$stayTime <- rep(NA, nrow(detdf))
  movedf$time <- rep(1:nbObs, length(unique(movedf$ID)))
  
  polyID <- sapply(slot(X, "polygons"), function(x) slot(x, "ID"))
  
  for (i in 1: nrow(detdf)){
    t <- detdf[i,'time']; id <- detdf[i,'ID']; detector <- detdf[i,'Site']
    
    track <- movedf %>%
      filter(ID == id & time %in% (t-1):t ) %>%  as.data.frame()
    coords <- coordinates(track[,c('x','y')])
    trackLine <- SpatialLines(list(Lines(list( Line(coords)), "id" )))
    
    detectorID <- which(polyID == detector)
    
    trackIntersect <- gIntersection(X[detectorID,],trackLine)
    detdf$stayTime[i] <- timestep * 60 * gLength(trackIntersect)/gLength(trackLine) # in seconds (* 60 cuz timestep is in mins)
  }
  return(detdf)
}

# for REM
calcSpeed <- function(movedf = moveSims, detdf = x1, X = traps.buff, timestep = 15, nbObs = nbObs){
  detdf$dist <- detdf$stayTime <- rep(NA, nrow(detdf))
  movedf$time <- rep(1:nbObs, length(unique(movedf$ID)))
  
  polyID <- sapply(slot(X, "polygons"), function(x) slot(x, "ID"))
  
  for (i in 1: nrow(detdf)){
    t <- detdf[i,'time']; id <- detdf[i,'ID']; detector <- detdf[i,'Site']
    
    track <- movedf %>%
      filter(ID == id & time %in% (t-1):t ) %>%  as.data.frame()
    coords <- coordinates(track[,c('x','y')])
    trackLine <- SpatialLines(list(Lines(list( Line(coords)), "id" )))
    
    detectorID <- which(polyID == detector)
    
    # while speed can directly be calculated from gLength(trackLine) since time units are constant in sims,
    # this matters for aggregating consecutive encounters which means more than one straight-line track
    trackIntersect <- gIntersection(X[detectorID,],trackLine) # this calculates distance
    detdf$dist[i] <- gLength(trackIntersect)
    detdf$stayTime[i] <- timestep * 60 * gLength(trackIntersect)/gLength(trackLine) # in seconds (* 60 cuz timestep is in mins)
    
  }
  return(detdf)
}

# for CT-DS
calcDist <- function(movedf = moveSims, detdf = x1, X = traps.buff, timestep = 15, nbObs = nbObs){
  dist_to_CT <- list()
  movedf$time <- rep(1:nbObs, length(unique(movedf$ID)))
  
  polyID <- sapply(slot(X, "polygons"), function(x) slot(x, "ID"))
  
  for (i in 1: nrow(detdf)){
    t <- detdf[i,'time']; id <- detdf[i,'ID']; detector <- detdf[i,'Site']
    
    track <- movedf %>%
      filter(ID == id & time %in% (t-1):t ) %>%  as.data.frame()
    coords <- coordinates(track[,c('x','y')])
    trackLine <- SpatialLines(list(Lines(list( Line(coords)), "id" )))
    
    detectorID <- which(polyID == detector)
    
    trackIntersect <- gIntersection(X[detectorID,],trackLine)
    segRatio <- seq(0,gLength(trackIntersect),gLength(trackLine)/(15*60))/gLength(trackIntersect)  # line segments
    trackIntersectxy <- trackIntersect@lines[[1]]@Lines[[1]]@coords
    
    x3 = segRatio * trackIntersectxy[2,'x'] + (1 - segRatio) * trackIntersectxy[1,'x'] #find point on each segment
    y3 = segRatio * trackIntersectxy[2,'y']  + (1 - segRatio) * trackIntersectxy[1,'y'] 
    lineSegs <- data.frame(x=x3,y=y3)
    
    ct_midpt <- X[detectorID,]@polygons[[1]]@labpt # this is the midpoint of the CT polygon
    ct_coord <- matrix(round(ct_midpt/50)*50,1,2) # get the actual CT location
    dist_to_CT[[i]] <- e2dist(ct_coord, lineSegs)
  }
  
  dupTimes <- unlist(lapply(dist_to_CT, length))
  detdf_id <- rep(1:nrow(detdf), dupTimes)
  detdf_new <- detdf[detdf_id,]
  detdf_new$distance <- do.call(c,dist_to_CT)
  
  return(detdf)
}


e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
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

