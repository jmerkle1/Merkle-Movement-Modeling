#----------------------#
#- Random walks ----####
#----------------------#
#--- Jerod Merkle -----#
#------- Lab 1 --------#

# packages
library(adehabitatLT)
library(adehabitatHR)
library(parallel)
library(circular)
library(sf)
library(lubridate)
library(moveHMM)

#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#
data <- readRDS("./data/bhs_data.rds")
head(data)

# -------------------------------------#
# reduce data to season in question ####
# -------------------------------------#
# Wyoming is a seasonal environment, and movement and habitat selection varies stronly be season
# The following represents what I am doing, but you should think about your own season.
# If you want to do a year around analysis, let me know as we should discuss.

#create season column
data$year <- year(data$date)
data$month <- month(data$date)
# create a datafarme of what months you want what season to be (of course you should alter this!)
seas_mat <- data.frame(month=1:12,
                       season=c("winter","winter","winter","springfall",
                                "springfall","summer","summer","summer",
                                "summer","springfall","springfall","winter"))
seas_mat
data <- merge(data, seas_mat)   # add the season column to your database using merge()
rm(seas_mat)

# create a index for each aid, each year, and each season
data$id_yr_seas <- paste(data$id, data$year, data$season, sep="_")
table(data$id_yr_seas)
length(unique(data$id_yr_seas))

#reduce your database to the season of interest (in my case, I'm doing summer)
data <- data[data$season %in% c("summer"),]
table(data$id_yr_seas)

# figure out whether you want to keep every id_yr_seas 
# in other words, you may want to remove some id_yr_seas' where there are too few points (e.g., animal died)
tbl <- table(data$id_yr_seas)
hist(tbl)  # I can see that not all animals were marked for the entire season
sort(tbl)
ids2keep <- names(tbl)[tbl > 350]   #I'm going to keep the id_yr_seas's with at least 350 locations! (YOU WILL NO DOUBT HAVE TO ADJUST THIS!)
data <- data[data$id_yr_seas %in% ids2keep,]
rm(tbl, ids2keep)

data$year <- NULL   #remove the excess columns
data$month <- NULL
data$season <- NULL
head(data)

# -------------------------#
# save cleaned database ####
# -------------------------#

# write out as rds object (more efficient). I suggest you do this, and then 
# you'll have a clean season database ready for each subsequent lab in the class
saveRDS(data, file="./data/bhs_data_summer.rds")


#----------------------------------#
# calculate movement parameters ####
#----------------------------------#
source("./MiscFunctions/CalcBurst.R")
source("./MiscFunctions/CalcMovParams.R")
data <- data[order(data$id_yr_seas, data$date),] #order database first
data$burst <- CalcBurst(data=data, id = TRUE, 
                        id_name="id_yr_seas", 
                        date_name="date", Tmax = 3600*7) #set Tmax to 7 hours, a little more than double the fix rate
length(unique(data$burst))
data <- CalcMovParams(data=data, id_name = "id_yr_seas", date_name = "date", burst=data$burst)

hist(data$rel.angle)   #this is your turning angle distribution
hist(data$dist/1000)   #this is your step length distribution
rm(CalcBurst, CalcMovParams)
head(data)

#-------------------------#
# simulate random walk ####
#-------------------------#
# first calculate brownian motion scaling factor for each ID (analogous to diffusion coefficient)
# turn our dataframe into class ltraj so we can use adehabitatHR's hbrown function
data_ltraj <- as.ltraj(st_coordinates(data), date=data$date, 
                       id=data$id_yr_seas, typeII=TRUE)
summary(data_ltraj)
# extract the brownian motion scaling factor from each id, so we can develop a comparible random walk model
h <- hbrown(data_ltraj)     # larger values give smaller dispersion
rm(data_ltraj)
hist(h)

# simulate random walk based on observed scaling factor for each id_yr_seas
# run it on multiple processors
ids <- unique(data$id_yr_seas) #loop over ids
table(ids %in% names(h)) #this needs to be all true
data <- data[order(data$id_yr_seas, data$date),] # make sure you order your data before calculating MSD.

# identify cores (use 1 less than you have)
no_cores <- detectCores()-1
# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("data", "ids", "h"))

# now for the actual loop
data <- do.call(rbind, clusterApplyLB(clust, 1:length(ids), function(i){
  # need library() here for the packages your calculations require for your calculations
  library(adehabitatLT)
  library(sf)
  
  temp <- data[data$id_yr_seas == ids[i],]
  #calculate net squared displacement in km2 from first point
  temp$NSD <- sqrt((st_coordinates(temp)[1,1]-st_coordinates(temp)[,1])^2 + (st_coordinates(temp)[1,2]-st_coordinates(temp)[,2])^2)/1000 # divide by 1000 to get to kms
  temp$NSD <- temp$NSD^2
  # plot(temp$date, temp$NSD)
  temp$order <- 1:nrow(temp)  #going to be used to make population level MSD later
  
  #now simulate 100 random movements based on same scaling factor from the data h
  rand <- do.call(cbind, lapply(1:100, function(e){
    return(simm.brown(date=temp$date, 
                      x0=c(st_coordinates(temp)[1,1],st_coordinates(temp)[1,2]), 
                      h=h[names(h)==ids[i]], id="a1")[[1]]$R2n/1000000)  # the last bit of code here at the end extracts the MSD and divides by 1 million to get to KM squared
  }))
  # calculate means and upper and lower CIs of the 100 simulated tracks
  temp$MSD_rand_mean <- apply(rand, 1, mean)
  SE <- apply(rand, 1, sd)/sqrt(apply(rand, 1, length))
  temp$MSD_rand_LCI <- temp$MSD_rand_mean-SE*1.96
  temp$MSD_rand_UCI <- temp$MSD_rand_mean+SE*1.96
  
  # plot the data
  #   plot(temp$date, temp$rand_mean, type="l", xlab="Time", ylab="MSD")
  #   lines(temp$date, temp$rand_LCI, col="grey")
  #   lines(temp$date, temp$rand_UCI, col="grey")
  #   lines(temp$date, temp$NSD, col="red")
  return(temp)
}))
stopCluster(clust)   # you must stop the parallelization framework
head(data)

# plot your data by ID
for(i in 1:length(ids)){   #this will plot ALL OF THEM, then use the back button on your plotting pane
  temp <- data[data$id_yr_seas == ids[i],]
  plot(temp$date, temp$MSD_rand_mean, type="l", xlab="Time", 
       ylab="MSD (km^2)", main=ids[i], ylim=range(c(temp$MSD_rand_LCI, temp$MSD_rand_UCI)))
  lines(temp$date, temp$MSD_rand_LCI, col="grey")
  lines(temp$date, temp$MSD_rand_UCI, col="grey")
  lines(temp$date, temp$NSD, col="red")
  legend("topleft", c("Random walk","observed"), lty=1, col=c("black","red"), bty="n")
}   #hit back button on your plot window to scroll through your each id_yr_seas
rm(temp, i, ids, no_cores, clust)

# -------------------------------------------#
# population level MSD using mean h value ####
# -------------------------------------------#
h <- mean(h)
popMSD <- as.data.frame(tapply(data$NSD, data$order, mean))  #do a mean of NSD across order (i.e., ordered steps)
names(popMSD) <- "MSD"
popMSD$date <- rownames(popMSD)
popMSD$MSD_SE <- as.numeric(tapply(data$NSD, data$order, sd))/
  sqrt(as.numeric(tapply(data$NSD, data$order, length)))   # calculate Population SE
table(is.na(popMSD)) #should only be very few NAs at the end of the dataframe
head(popMSD)

#now simulate random movement based on same scaling factor from the data h
#find id_yr_seas with the full sweet of date values
date_range <- data$id_yr_seas[data$order==max(popMSD$date)][1]
date_range <- range(data$date[data$id_yr_seas == date_range])
#add the ids date info to the popMSD dataframe
popMSD$datePOSIX <- seq(date_range[1], date_range[2], length.out=nrow(popMSD))
#simulate 100 random tracks with the population h
rand <- do.call(cbind, lapply(1:100, function(e){
  return(simm.brown(date=popMSD$datePOSIX, 
                    x0=c(0,0), h=h, id="a1")[[1]]$R2n/1000000)
}))
dim(rand)
nrow(rand)==nrow(popMSD) #should be TRUE
# calculate means and upper and lower CIs of the simulated random walks
popMSD$rand_mean <- apply(rand, 1, mean)
popMSD$rand_SE <- apply(rand, 1, sd)/sqrt(apply(rand, 1, length))
rm(rand, date_range)
head(popMSD)

# now plot population MSD vs random walk ####
plot(popMSD$datePOSIX, popMSD$rand_mean, type="l", xlab="Time", 
     ylab="MSD (km^2)", ylim=range(c(popMSD$rand_mean,popMSD$rand_mean+2*popMSD$rand_SE)))
lines(popMSD$datePOSIX, popMSD$rand_mean+1.96*popMSD$rand_SE, col="grey")
lines(popMSD$datePOSIX, popMSD$rand_mean-1.96*popMSD$rand_SE, col="grey")
lines(popMSD$datePOSIX, popMSD$MSD, col="red")
lines(popMSD$datePOSIX, popMSD$MSD+1.96*popMSD$MSD_SE, col="grey")
lines(popMSD$datePOSIX, popMSD$MSD-1.96*popMSD$MSD_SE, col="grey")
legend("topleft", c("Random walk","observed"), lty=1, col=c("black","red"), bty="n")

rm(h)

#----------------------#
# State space model ####
#----------------------#
# from Morales et al. 2004 and Beyer et al. 2013

# Organize your data ...

## to look at state space models, we will use the moveHMM package. We already pulled out all of the 
## movement parameters (distance, turning angle, etc.) but because of how this package works, we will 
## calculate them again within the functions of the package

## we are first going to look at individual states without accounting for any additional covariates
## (in the Morales paper, they included distance to water in the models)

## first pull out only the ID, date, and coordinates. you MUST name them exactly these things
## to use the functions in moveHMM

df<-data.frame(ID=data$id_yr_seas,
               Time=data$date,
               X=st_coordinates(data)[,1],
               Y=st_coordinates(data)[,2])
head(df)

## this function preps the data to run a state space model
## put in the dataframe you just made, the coordinate type and the name of the coordinate columns
data.hmm <- prepData(df, type="UTM", coordNames = c("X", "Y"))

## this next step will first plot the movement of each of your individuals (colored by ID)
## and will then plot various movement parameters for each individual animal
plot(data.hmm, compact = TRUE)  

summary(data.hmm) ##gives you summary information

## if you have steps with step length zero, you have to account for that in the models
## by including a zero-mass starting parameter because the strictly positive distributions
## (i.e., gamma) are inadequate to use. You can either account for this by adding in a zero-mass
## parameter, or you can add a tiny amount of positive variation to the steps with length 0
## we are going to do that for similicity sake. For more information about how to deal with zero-
## inflation, you can read about it in the moveHMM vignette.

## first check to see if you have any steps with length 0
table(data.hmm$step==0) ## for my data, there are 22 steps with step length 0

## add in some variation for those 22 steps
data.hmm$step<-ifelse(data.hmm$step==0, 0+runif(1, 0, .001), data.hmm$step)

## after you add in the variation, these should all be FALSE, if they are not, something is wrong
table(data.hmm$step==0)

## before we run models with different states, we need to set some starting parameters for the models to start from
## we can use the mean and distribution of our actual data to decide on these
## looking at the distribution of steps and turning angles will also let us decide what
## distribution to use

hist(data.hmm$step)
hist(data.hmm$angle)

mean(data.hmm$step, na.rm=T)
sd(data.hmm$step, na.rm=T)

mean(data.hmm$angle, na.rm=T)
sd(data.hmm$angle, na.rm=T)


## to fit all of the models, we will use the fitHMM function in moveHMM
## for each of the models you run in the fitHMM function, you need to provide a starting point
## for each state you are trying to identify.

## the parameters depend on the distribution that you use, for the BHS data, we will use gamma and vm distribution
## you can use gamma, Weibull, log-normal, or exponential for the step lengths
## and vonMises or wrapped Cauchy for the angles
## look at the vignette for more info on what parameters you should use for each distribution

## if you get errors when running the models, and all data that you prepped looks correct
## the errors likely are a result of the optimizer failing to converge
## this happens if you set your initial parameters to something that doesn't quite make sense
## i.e., if your average step length is 8 m and you set the mean to 8,000 it won't converge
## so just mess around with the starting parameters until you get it to properly converge

## okay, now we have a general idea of our movement data and know what distribution we want to use
## we can start with the single-state model. For the single state model, choosing the starting parameters is
## pretty easy - just put in the means and sd of the steps and angles

# ----------------------#
# SINGLE STATE MODEL ####
# ----------------------#

# give some starting parameter values:
mu0 <- c(300) # step mean (need to set one for each state, so for single state, you only need one)
sigma0 <- c(400) # step SD 
stepPar0 <- c(mu0,sigma0) ## combine the mean and SD into one object to fit in the model

angleMean0 <- c(-0.025) # angle mean (need to set one for each state)
kappa0 <- c(1.8) # angle concentration
anglePar0 <- c(angleMean0,kappa0) ## combine into one object


## with verbose=2, you get information for each iteration, it tells you if you are getting close to the 
## solution. You can turn that off by setting verbose=0, or if you want to see only the first and last iteration, set
## verbose=1, The worse your starting parameters are, the longer it takes to find a solution (and more likely to fail)

singleStateModel <- fitHMM(data.hmm, nbStates = 1, 
                          stepPar0 = stepPar0, 
                          anglePar0 = anglePar0, verbose=2,
                          stepDist = "gamma", angleDist = "vm")

##can look at model output
singleStateModel ## this gives you the mean and sd of your states (in this case a single state).

## can look at the CI for the model
CI(singleStateModel)

## then plot the model - this returns a plot of the distribution of step length and and a plot of the turning angle for
## the single state
plot(singleStateModel)

## you don't get individual plots of movement for the single state model because there is no switching
## between states. if you want a plot of movement for each individual, you can use: plot(data.hmm, compact=FALSE)

#-----------------------#
# DOUBLE STATE MODEL ####
#-----------------------#

## give some starting parameter values, because we are trying to identify 2 states now, you need to put in
## the starting parameters that you think are close for each state.

mu0 <- c(250, 1000) # step mean (need one for each state, so should have 2)
sigma0 <- c(400,1250) # step SD 
stepPar0 <- c(mu0,sigma0) ## combine the mean and SD into one object for the model

angleMean0 <- c(-2.7, 0.25) # angle mean
kappa0 <- c(0.1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

## run the double state model
doubleStateModel <- fitHMM(data.hmm, nbStates = 2, 
                           stepPar0 = stepPar0,
                           anglePar0 = anglePar0, verbose=2,
                           stepDist = "gamma", angleDist = "vm")

##can look at model output
doubleStateModel

## can look at the CI for the model
CI(doubleStateModel)

## then plot the model - this returns a plot of the distribution of step length and and a plot of the turning angle for
## the two states determined by the model; it then will give you an output for each indivdiual animal, with the 
## states colored in orange and blue
plot(doubleStateModel)

## this gives you a series of plots for each individual animal. The top plot shows you if they are in state 1 or 2,
## and the bottom two plots gives you the probability that the animal is in state 1 (middle plot) or state 2 (bottom plot)
plotStates(doubleStateModel)

#-----------------------#
# TRIPLE STATE MODEL ####
#-----------------------#

## give some starting parameter values, this time we need three means and SDs

mu0 <- c(50, 200, 500) # step length mean
sigma0 <- c(50, 200, 500) # step length SD
stepPar0 <- c(mu0,sigma0)

angleMean0 <- c(0, pi, 0) # angle mean
kappa0 <- c(0.1, 0.5, 0.1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

tripleStateModel <- fitHMM(data.hmm, nbStates = 3, # Note, this took about 3 minutes
                          stepPar0 = stepPar0, 
                          anglePar0 = anglePar0, verbose=2,
                          stepDist = "gamma", angleDist = "vm")

##can look at model output
tripleStateModel

## can look at the CI for the model
CI(tripleStateModel)

## plot the model, should have the distribution for all three states now, and will plot the states for each
## individual in three colors
plot(tripleStateModel)

## look at the states and probability that each step is in a specific state
plotStates(tripleStateModel)

tmp <- 
table(tmp)

# mapview the results
data$state <- viterbi(tripleStateModel)   # add states as a column in your database
toplot <- data[data$id_yr_seas == "38898_2019_summer",]  # grab an AID to plot
mapview(toplot, zcol = 'state')   # plot it
mapview(st_cast(st_combine(toplot), "LINESTRING"), layer.name="lines", ???) # how to get color of lines to match?


#--------------------#
# MODEL SELECTION #### 
#--------------------#

## you can now pull out the AIC for each of the models to compare and figure out which one is best
AIC(singleStateModel, doubleStateModel, tripleStateModel)

## for my bighorn sheep, the best model is the tripleStateModel. 

## to get the figures you need for the first part of your write-up, you can use the plot() and 
## plotStates() functions in the moveHMM package. 



