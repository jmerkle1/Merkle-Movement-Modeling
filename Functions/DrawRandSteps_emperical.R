#######################################################################################
### Function to draw Random steps for Step Selection analysis, ##################################################
### derived from Mathieu Bastille's code, updated by Jerod Merkle ######################
################## Feb 2021 ###################################################


# This function will take a dataframe of locations, that includes step length, dist, absolute
# and relative turning angles, and generates a specified number of random steps.  The output 
# dataframe actually now represents a line for each step.  nr specifies
# the number of random points you want for each step.  distMax is an optional number to 
# specify if you want to limit the greatest distance (in same dist units as your dist column
# in your dataframe).  When simult is TRUE the function will simultaneously draw a random
# step with a random turning angle.  If it is false, the two parameters will be drawn 
# indepenedent of each other.  When withinID is TRUE, the random steps will be drawn from 
# the distribution of step lengths and turning angles within the ID in question. 
# If withinID is FALSE, it will draw random steps from the distribution of step lengths
# and turning angles of all animals. The date column now corresponds to the date of the beginning of the step.
# The strata column provides the strata for the steps.  If you want to just draw random angles, 
# and not use angles from the data, you can 
# specify uniform.angles = TRUE.  Some notes: 1) you must order your database first! 


DrawRandSteps_emperical <- function(data=data, nr = 10, distMax = Inf, simult = FALSE,
                                    id_name="id", date_name="date", x_name="x", y_name="y",
                                    withinID = FALSE, uniform.angles = FALSE) {
  
  if(inherits(data, "sf") == TRUE) stop("data should be a data.frame")
  if(inherits(data, "data.frame") != TRUE) stop("data should be a data.frame")
  if(any(colnames(data) == date_name) == FALSE) stop(print("Your date_name is incorrect."))
  if(any(colnames(data) == id_name) == FALSE) stop(print("Your id_name is incorrect."))
  if(any(colnames(data) == x_name) == FALSE) stop(print("Your x_name is incorrect."))
  if(any(colnames(data) == y_name) == FALSE) stop(print("Your y_name is incorrect."))
  if("x_end" %in% names(data) | "strata" %in% names(data))
    stop("You have already calculated Random Steps! Need to start over with a clean dataframe to recalculate.")
  
  #manage packages
  if(all(c("circular","plyr","parallel") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: circular, plyr, and snowfall")
  require(circular)
  require(plyr)
  require(parallel)
  
  if(inherits(data[,date_name], "POSIXct") != TRUE) stop("date_name is not POSIXct")
  if(any(duplicated(data[c(id_name, date_name)])) == TRUE) stop("You have duplicates in your database")
  if(any(is.na(data[,date_name]) == TRUE)) stop ("You have NAs in your date column")
  if(any(is.na(data[,id_name]) == TRUE)) stop ("You have NAs in your id column")
  
  # check if the steps are fairly similar in length
  step_dt <- unique(na.omit(data$dt[data$StepFlag == TRUE]))
  if(diff(range(step_dt)) > 1200) 
    print("WARNING!!!!! You have a mix of step lengths in this dataset, and they differ by more than 20 mins!!! It is not a good idea to draw random steps without fixing this.")
  if(any(c("rel.angle", "abs.angle") %in% colnames(data) == FALSE) == TRUE) 
    stop("you do not have an abs.angle and rel.angle columns in data. You need to rerun calcMovParams")
  rm(step_dt)
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
  #add case column and note 1s because these are the used points
  data$case <- ifelse(is.na(data$strata)==TRUE, NA, 1)
  #add x_end/y_end for the used points
  data$x_end <- c(data[2:nrow(data),x_name],NA)
  data$y_end <- c(data[2:nrow(data),y_name],NA)
  data$x_end[is.na(data$strata)==TRUE] <- NA
  data$y_end[is.na(data$strata)==TRUE] <- NA
  
  #make north object... used to reorient angles later
  north <- c(NA, data$abs.angle[1:(nrow(data)-1)])
  
  #loop over ids i.e., the u object
  
  # identify cores (use 1 less than you have)
  no_cores <- detectCores()-1
  # Setup cluster
  clust <- makeCluster(no_cores) 
  # export the objects you need for your calculations from your environment to each node's environment
  parallel::clusterExport(clust, varlist=c("data","simult","withinID","distMax","uniform.angles", "nr", "u", "north",
                                 "id_name","date_name", "x_name", "y_name"),envir=environment())
  
  # now for the actual loop
  temp <- clusterApplyLB(clust, 1:length(u), function(i){
    # need library() here for the packages your calculations require for your calculations
    library(circular)
    
    num <- length(data$strata[is.na(data$strata)==FALSE & data[,id_name] == u[i]])
    if(num!=0){
      
      #first create databases to draw the steps from.  The first a list for each
      #id, the second is a list for each id where it is from with the other ids taken out
      if(withinID == FALSE){
        dlist2 <- data[data$StepFlag == TRUE,]
        dlist <- NA
      }else{
        dlist <- data[data$StepFlag == TRUE & data[,id_name] == u[i],]
        dlist2 <- NA
      }
      
      #select the step lengths and turning angles
      if(simult == TRUE){
        if(withinID == FALSE){
          ts <- dlist2[sample(1:nrow(dlist2), nr*num, replace=TRUE),]
        }else{ #when withinID equals TRUE
          ts <- dlist[sample(1:nrow(dlist), nr*num, replace=TRUE),]
        }
      }else{ #when simult equals FALSE
        if(withinID == FALSE){
          ts <- data.frame(rel.angle=sample(dlist2[,"rel.angle"], nr*num, replace=TRUE), 
                           dist=sample(dlist2[,"dist"], nr*num, replace=TRUE))
        }else { #when withinID equals TRUE
          ts <- data.frame(rel.angle=sample(dlist[,"rel.angle"], nr*num, replace=TRUE), 
                           dist=sample(dlist[,"dist"], nr*num, replace=TRUE))
        }
      }
      
      #add strata and case
      rng <- range(data$strata[data[,id_name] == u[i]], na.rm=T)
      ts$strata <- rep(rng[1]:rng[2], each=nr)
      ts$case <- 0
      
      #adjust the angles to correspond with north based on the previous step
      north <- north[is.na(data$strata)==FALSE & data[,id_name]==u[i]]
      if(length(north)!=num)
        stop("you have problem 2")
      north <- rep(north, each=nr)
      angles <- ts$rel.angle
      ts$rel.angle <- ifelse(ts$rel.angle+north < 360, ts$rel.angle+north, (ts$rel.angle+north)-360)
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
      names(ts)[names(ts)=="id"] <- id_name    #name the id column back to the original name
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


