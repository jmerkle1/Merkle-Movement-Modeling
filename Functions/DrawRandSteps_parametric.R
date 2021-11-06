#######################################################################################
### Function to create Random steps based on #######
### parametric sampling of step lengths and turning angles #######
### by Jerod Merkle ######################
################## Feb 2021 ##################


# This function will take a dataframe of locations, that includes step length, dist, absolute
# and relative turning angles, and generates a specified number of random steps.  The output 
# dataframe actually now represents a line for each step (there will be an NA in x_end if there
# is no step).  nr specifies
# the number of random points you want for each step.  distMax is an optional number to 
# specify if you want to limit the greatest distance (in same dist units as your dist column
# in your dataframe).   When withinID is FALSE, the random steps will be drawn from 
# the chosen distribution of step lengths of all the animials in the database.
# If withinID is TRUE, it will draw random steps from the animals' distribution of step lengths
#  The date column now corresponds to the date of the beginning of the step.
# The strata column provides the strata for the steps.  If speed = TRUE, than the steps
# are generated on the speed values in the data, not dist values - mainly if you have varying fix rates.
# If you want to just draw random angles, 
# and not use angles from the data, you can specify uniform.angles = TRUE.
# Some notes: 1) you must order your database first! 
# Can choose weibull, gamma, or exponential (for step length/speed)
# Can choose wrappedcauchy or vonmeses (for turning angles)

DrawRandSteps_parametric <- function(data=data, nr = 10, distMax = Inf, 
                                     step_distr = "weibull", ta_distr="wrappedcauchy",
                                     id_name="id", date_name="date", x_name="x", y_name="y",
                                     withinID = FALSE, speed=FALSE, uniform.angles=FALSE) {
  
  if(inherits(data, "sf") == TRUE) stop("data should be a data.frame. Take it out of sf!")
  if(inherits(data, "data.frame") != TRUE) stop("data should be a data.frame")
  if(any(colnames(data) == date_name) == FALSE) stop(print("Your date_name is incorrect."))
  if(any(colnames(data) == id_name) == FALSE) stop(print("Your id_name is incorrect."))
  if(any(colnames(data) == x_name) == FALSE) stop(print("Your x_name is incorrect."))
  if(any(colnames(data) == y_name) == FALSE) stop(print("Your y_name is incorrect."))
  if("x_end" %in% names(data) | "strata" %in% names(data))
    stop("You have already calculated Random Steps! Need to start over with a clean dataframe to recalculate.")
  if(step_distr %in% c("weibull","gamma","exponential") == FALSE)
    stop("Your step_distr must be either weibull,gamma, or exp.")
  if(ta_distr %in% c("wrappedcauchy","vonmeses") == FALSE)
    stop("Your step_distr must be either wrappedcauchy or vonmeses.")
  
  
  #manage packages
  if(all(c("circular","plyr","parallel","MASS") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: circular, plyr, and snowfall")
  require(circular)
  require(plyr)
  require(parallel)
  require(MASS)
  
  if(inherits(data[,date_name], "POSIXct") != TRUE) stop("date_name is not POSIXct")
  if(any(duplicated(data[c(id_name, date_name)])) == TRUE) stop("You have duplicates in your database")
  if(any(is.na(data[,date_name]) == TRUE)) stop ("You have NAs in your date column")
  if(any(is.na(data[,id_name]) == TRUE)) stop ("You have NAs in your id column")
  
  if(speed==FALSE){
    # check if the steps are fairly similar in length
    step_dt <- unique(na.omit(data$dt[data$StepFlag == TRUE]))
    if(diff(range(step_dt)) > 1200) 
      print("WARNING!!!!! You have a mix of step lengths in this dataset, and they differ by more than 20 mins!!! It is not a good idea to draw random steps without fixing this. Unless, however, you specify speed=TRUE")
    rm(step_dt)
    
    if(step_distr %in% c("weibull","gamma") & any(data$dist[data$StepFlag == TRUE] == 0))
      stop("You have dists that are 0. You cannot have 0 dists when fitting weibull or gamma distributions to step data")
  }else{
    if(step_distr %in% c("weibull","gamma") & any(data$speed[data$StepFlag == TRUE] == 0))
      stop("You have speeds that are 0. You cannot have 0 speeds when fitting weibull or gamma distributions to step data")
  }
  if(any(c("rel.angle", "abs.angle", "dist") %in% colnames(data) == FALSE) == TRUE) 
    stop("you do not have an dist, abs.angle and rel.angle columns in data. You need to rerun CalcMovParams")
  
  key <- 1:nrow(data)
  key2 <- key[order(data[,id_name], data[,date_name])]
  if(all(key==key2)==FALSE) stop(print("Your data are not ordered correctly"))
  rm(key, key2)
  u <- unique(data[,id_name])
  
  # grab extra column names
  extr_col_nms <- names(data)
  extr_col_nms <- extr_col_nms[extr_col_nms %in% c(id_name,date_name, x_name, y_name,
                                                   "burst","dist","dt","speed","abs.angle",
                                                   "rel.angle","StepFlag")==FALSE]
  
  #make strata numbers
  data$strata <- 1:nrow(data)
  data$strata[data$StepFlag == FALSE] <- NA   # NA when a step should NOT be connected using StepFlag
  data$strata <- as.numeric(as.factor(data$strata))
  #add case for the 1s
  data$case <- ifelse(is.na(data$strata)==TRUE, NA, 1)
  #add x_end numbers
  data$x_end <- c(data[2:nrow(data),x_name],NA)
  data$y_end <- c(data[2:nrow(data),y_name],NA)
  data$x_end[is.na(data$strata)==TRUE] <- NA
  data$y_end[is.na(data$strata)==TRUE] <- NA
  
  #make north file... used to reorient angles later
  north <- c(NA, data$abs.angle[1:(nrow(data)-1)])
  
  #loop over ids i.e., the u object
  
  # identify cores (use 1 less than you have)
  no_cores <- detectCores()-1
  # Setup cluster
  clust <- makeCluster(no_cores) 
  # export the objects you need for your calculations from your environment to each node's environment
  parallel::clusterExport(clust, varlist=c("data","withinID","distMax", "nr", "u","speed", "uniform.angles",
                                 "north","id_name","date_name", "x_name", "y_name","step_distr","ta_distr"),
                envir=environment())
  
  # now for the actual loop
  temp <- clusterApplyLB(clust, 1:length(u), function(i){
    # need library() here for the packages your calculations require for your calculations
    library(circular)
    library(MASS)
    
    num <- length(data$strata[is.na(data$strata)==FALSE & data[,id_name] == u[i]])
    if(num!=0){
      
      #first fit distributions from drawning the steps from.  The first a list for each
      #id, the second is a list for each id where it is from the ids taken out
      if(withinID == FALSE){
        
        dlist <- data[data$StepFlag == TRUE,]

        if(speed==FALSE){
          params_sl <- fitdistr(dlist$dist,densfun=step_distr)
        }else{ #if speed = TRUE
          params_sl <- fitdistr(dlist$speed,densfun=step_distr)
        }
        # get param estimates for turning angles
        if(ta_distr == "wrappedcauchy"){
          params_ta <- mle.wrappedcauchy(circular(dlist$rel.angle, units="degrees"))
        }else{   # if vonmises
          params_ta <- mle.vonmises(circular(dlist$rel.angle, units="degrees"))
        }
        
      }else{  # when withinID == TRUE
        dlist <- data[data$StepFlag == TRUE & data[,id_name] == u[i],] # grab steps for only that ID
        if(speed==FALSE){
          params_sl <- try(fitdistr(dlist$dist,densfun=step_distr),silent=TRUE)
          if(class(params_sl)=="try-error"){  # sometimes there isn't enough points to fit a distr. So, do a bit of bootstrapping
            params_sl <- fitdistr(rep(dlist$dist,10),densfun=step_distr)
          }
        }else{  #if speed = TRUE
          params_sl <- try(fitdistr(dlist$speed,densfun=step_distr),silent=TRUE)
          if(class(params_sl)=="try-error"){  # sometimes there isn't enough points to fit a distr. So, do a bit of bootstrapping
            params_sl <- fitdistr(rep(dlist$speed,10),densfun=step_distr)
          }
        }
        # get param estimates for turning angles
        if(ta_distr == "wrappedcauchy"){
          params_ta <- mle.wrappedcauchy(circular(dlist$rel.angle, units="degrees"))
        }else{   # if vonmises
          params_ta <- mle.vonmises(circular(dlist$rel.angle, units="degrees"))
        }
      }  # end of withinID == TRUE
      
      #select the step lengths and turning angles
      # do angles first
      if(ta_distr == "wrappedcauchy"){
        rel.angle=as.numeric(rwrappedcauchy(nr*num,params_ta$mu,params_ta$rho))
      }else{
        rel.angle=as.numeric(rvonmises(nr*num,params_ta$mu,params_ta$kappa))
      }
      
      # now for dists or speeds
      if(speed==FALSE){
        if(step_distr == "weibull"){
          dist <- rweibull(nr*num,params_sl$estimate[1],params_sl$estimate[2])
        }else{  # not weibull
          if(step_distr == "exponential"){
            dist <- rexp(nr*num,params_sl$estimate[1])
          }else{  # for gamma
            dist <- rgamma(nr*num,params_sl$estimate[1],params_sl$estimate[2])
          }
        }
        # build df of random steps
        ts <- data.frame(rel.angle=rel.angle,dist=dist)
        
      }else{  # for speed=TRUE
  
        if(step_distr == "weibull"){
          speeds <- rweibull(nr*num,params_sl$estimate[1],params_sl$estimate[2])
        }else{  # not weibull
          if(step_distr == "exponential"){
            speeds <- rexp(nr*num,params_sl$estimate[1])
          }else{  # for gamma
            speeds <- rgamma(nr*num,params_sl$estimate[1],params_sl$estimate[2])
          }
        }
        # build df of random steps
        ts <- data.frame(rel.angle=rel.angle,speed=speeds)
      } # end of speed=TRUE 
      
      #add strata and case
      rng <- range(data$strata[data[,id_name] == u[i]], na.rm=T)
      ts$strata <- rep(rng[1]:rng[2], each=nr)
      ts$case <- 0
      
      #ajust the angles to correspond with north based on the previous step
      north <- north[is.na(data$strata)==FALSE & data[,id_name]==u[i]]
      if(length(north)!=num)
        stop("you have problem 2")
      north <- rep(north, each=nr)
      angles <- ts$rel.angle
      ts$rel.angle <- ifelse(ts$rel.angle+north < 360, ts$rel.angle+north, (ts$rel.angle+north)-360)
      
      if(speed==TRUE){
        d <- data[is.na(data$strata)==FALSE & data[,id_name]==u[i],]
        ts <- merge(ts, d[,c("strata","dt")])
        ts$dist <- ts$speed*ts$dt 
      }
      
      if(distMax != Inf){  #if distMax is specified
        ts$dist[ts$dist > distMax] <- distMax
      }
      
      if(uniform.angles == TRUE){
        ts$rel.angle <- runif(length(ts$rel.angle), 0, 359.99999999)
        angles <- ts$rel.angle
      }
      
      d <- data[is.na(data$strata)==FALSE & data[,id_name]==u[i],]
      ts <- merge(ts[,c("strata","case","rel.angle","dist")], d[,c("x","y",date_name,"strata")])
      
      ts$id <- u[i]
      ts$x_end <- ts[,x_name]+(cos(rad(ts$rel.angle))*ts$dist)
      ts$y_end <- ts[,y_name]+(sin(rad(ts$rel.angle))*ts$dist)
      ts$rel.angle <- angles   
      names(ts)[names(ts)=="id"] <- id_name    #name teh id column back to the origional name
      return(ts)
    }else{
      return(NULL)
    }
  })
  stopCluster(clust)   # you must stop the parallelization process
  
  temp <- do.call(rbind, temp)
  temp <- rbind.fill(data, temp)
  temp <- temp[order(temp[,id_name], temp[,date_name], temp$case, decreasing = F),]
  
  # fill in some NAs 
  temp <- temp[is.na(temp$strata)==FALSE,]
  temp2 <- temp[temp$case==1, c("strata","dt","StepFlag",extr_col_nms)]
  temp$dt <- NULL
  temp$StepFlag <- NULL
  for(i in 1:length(extr_col_nms)){
    temp[,extr_col_nms[i]] <- NULL
  }
  temp <- merge(temp, temp2, all.x=TRUE)
  temp <- temp[order(temp[,id_name], temp[,date_name], temp$case, decreasing = F),]
  
  # recalculate the speed column
  temp$speed <- temp$dist/temp$dt
  
  # recalculate the abs.angle column
  #create column that calculates absolute turning angle (degrees), where
  # 0 correlates with North, and 180 with South.  East is 90, west is 270
  abs.angle <- atan2((temp$y_end-temp[,y_name]), (temp$x_end-temp[,x_name]))
  temp$abs.angle <- as.numeric(conversion.circular(circular(abs.angle,units = "radians"), units = "degrees", 
                                                   zero = pi/2, rotation = "clock", modulo = "2pi"))
  
  # recalc rownames
  row.names(temp) <- 1:nrow(temp)
  
  temp$burst <- NULL  # remove the burst column
  temp$StepFlag <- NULL  # remove the StepFlag column
  
  # move columns around
  temp <- temp[c("case", setdiff(names(temp), "case"))]
  temp <- temp[c("strata", setdiff(names(temp), "strata"))]
  
  return(temp)
}


