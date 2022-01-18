############FUNCTION finds problem points####################
#It needs a sf_POINT dataframe with (at a minimum): date (POSIXct), and burst (where continuous
#bursts should be calculated).  You have the option to identify a vector
#for which are your bursts.  Must have your data ordered before you add in
# A problem point is defined as a point where there is a high speed coming
# into the point and going away from the point.

FindProblemPts <- function(dat = data, date_name="date",
                           burst = data$burst, id_name="id", speedlim=3) {

  #some checks to start
  require(sf)
  if(inherits(dat, "sf") != TRUE) stop("data is not a sfc_POINT dataframe")
  orig <- dat
  dat$xtmp <- st_coordinates(dat)[,1]
  dat$ytmp <- st_coordinates(dat)[,2]
  dat <- st_set_geometry(dat, NULL)    #need to remove the

  if(any(colnames(dat) == date_name) == FALSE) stop(print("Your date_name is not correct"))
  if(!inherits(dat[,date_name], "POSIXct")) stop(print("date_name column is not POSIXct"))
  key <- 1:nrow(dat)
  key2 <- key[order(dat[,id_name], dat[,date_name])]
  if(all(key==key2)==FALSE) stop(print("Your data are not ordered correctly"))
  rm(key, key2)
  if(length(burst) != nrow(dat)) stop(print("dat and burst do not have the same length"))
  if(any(duplicated(dat[c(id_name, date_name)])) == TRUE) stop("You have duplicates in your database")
  if(any(is.na(dat[,date_name]) == TRUE)) stop ("You have NAs in your date column")
  if(any(is.na(dat[,id_name]) == TRUE)) stop ("You have NAs in your id column")

  #create indicator where NAs should be put
  flags1 <- c(diff(burst),1)
  flags2 <- c(1, diff(burst))

  #create a column that calculates the distance (in m) of each move
  # for current point to next point
  xy2 <- dat[2:nrow(dat),c("xtmp","ytmp")]
  names(xy2) <- c("x","y")
  xy2 <- rbind(xy2, data.frame(x=NA,y=NA))
  dist1 <- sqrt((dat$xtmp-xy2$x)^2 + (dat$ytmp-xy2$y)^2)
  dist1[flags1==1] <- NA

  # for previous point to current point
  xy2 <- dat[1:(nrow(dat)-1),c("xtmp","ytmp")]
  names(xy2) <- c("x","y")
  xy2 <- rbind(data.frame(x=NA,y=NA), xy2)
  dist2 <- sqrt((dat$xtmp-xy2$x)^2 + (dat$ytmp-xy2$y)^2)
  dist2[flags2==1] <- NA

  #create a column that calculates difference in time of each move
  # for current point to next point
  tz <- attr(dat[,date_name],"tzone")
  dt1 <- as.numeric(difftime(c(dat[2:nrow(dat),date_name],NA),dat[,date_name],tz = tz, units="secs"))
  dt1[flags1==1] <- NA

  # for previous point to current point
  dt2 <- as.numeric(difftime(dat[,date_name],c(dat[1,date_name],dat[1:(nrow(dat)-1),date_name]),tz = tz, units="secs"))
  dt2[flags2==1] <- NA

  #create a column that calculates speed (in m/sec)
  speed1 <- dist1/dt1
  speed2 <- dist2/dt2

  speed1[is.na(speed1)==TRUE] <- 0   #replace the NAs with a small speed (so we keep them)
  speed2[is.na(speed2)==TRUE] <- 0

  return(speed1 > speedlim & speed2 > speedlim)

} #end of function
