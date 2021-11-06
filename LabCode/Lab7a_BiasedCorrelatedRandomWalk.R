#--------------------------------------------#
#- Biased correlated random walks --------####
#--------------------------------------------#
#---------- Jerod Merkle --------------------#
#------------- Lab 7a ------------------------#

library(rjags)  # you must download jags (program for bayesian stats) from here: https://sourceforge.net/projects/mcmc-jags/files/
library(circular)
library(raster)
library(parallel) 
library(sf)
library(ggplot2)

#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#
data <- readRDS("./data/bhs_data_summer.rds")
head(data)

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

#grab previous steps' abs angle for use later 
data <- data[order(data$id_yr_seas, data$date),]  #order database first
data$abs.angle_prev <- c(NA, data$abs.angle[1:(nrow(data)-1)])

# remove steps that are really short (i.e., collar error and the animal really isn't moving)
table(data$dist < 30)   # i removed anything less than 30 meters
data <- data[is.na(data$dist)==FALSE,]  # remove steps with NAs in dist as well
data <- data[data$dist > 30,]


# Add a single environmental variable ####
# Identify an important layer of interest, or the most selected layer from previous analyses
# need both distance to, and direction to that variable.

# ------------------------------------------------------------------------
# HERE IS HOW TO CREATE A 'DIRECTION TO' RASTER FROM A RASTER WITH 1S FOR THAT FEATURE, AND NAs OTHERWISE

# step 1: create a dichotomous raster of an important feature (a raster with 1s for the feature and 0s otherwise)
slope <- raster("./data/GIS/Slope_degrees_30m.tif")
plot(slope)   #this layer is simply 1s for escapterr and 0s, otherwise.
rcl <- matrix(c(-Inf, 25, NA,  # from -Inf to 25 will be NA
                25, Inf, 1),   # from 25 to Inf will be 1
              byrow=T, ncol=3)
rcl
slope_descrete <- reclassify(slope, rcl=rcl)
plot(slope_descrete)
# If you have a large study area, or you want high res results, you should use the 'Euclidean Direction' Function in ArcGIS (it'll be much faster). 
writeRaster(slope_descrete, filename = paste(getwd(), "./data/GIS/slope_descrete.tif",sep=""), format="GTiff")
# now open slope_descrete.tif in ArcGIS, and use the Euclidean_direction function to calculate the direction to slope raster (use all defaults).
# with the default, it will save in the temp geodatabase
# then you should see a colidascope looking thing in arcgis
# save the resulting raster as .tif file from ArcGIS... to do this, right click and go to save and "export data"
# Then in the export raster data window, do format = TIFF, and name it and save into your class gis folder.
# load that file back into R

# load up direction (in degrees) to nearest escape terrain landcover
dir2steep <- raster(paste(getwd(), "./data/GIS/dir2slope.tif",sep=""))   # Note, you must chnage the arcGIS result from degrees to radians
# need to switch to radians for the model
dir2steep[] <- rad(values(dir2steep))
plot(dir2steep)

# extract the 'direction to' variable (note, data must represent the GPS point of the start of the movement step, not the end point of the step)
data$dir2steep <- extract(dir2steep, data)  # now we have direction to nearest XX landcover type
table(is.na(data$dir2steep))   # check for NAs
table(data$dir2steep == 0)   #this is how many steps originated within the lancover type of question (or exactly south of a patch)
hist(data$dir2steep)
rm(dir2steep, slope_descrete, rcl, slope)

#remove movements that originated inside the landcover type of interest (i.e., they must start outside the landcover type)
data <- data[data$dir2steep != 0,]
hist(data$dir2steep)   # there shouldn't be much of a pattern here

#fix turning angle so they are in radians, not degrees (we're focusing on abs angle, not relative angle.)
data$ta <- rad(data$abs.angle)   #observed turning angle
data$ta_prev <- rad(data$abs.angle_prev)   # the previous step's turning angle

#remove lines where there is an at least 1 NA
table(is.na(data[,c("ta", "ta_prev","dir2steep")]))
data <- data[apply(data[,c("ta", "ta_prev","dir2steep")],1,anyNA)==FALSE,]

# reduce to only 1 random step per am and pm, per day, per id, per year to avoid multicollinearity
data$ampm <- strftime(data$date, format = "%p", tz = "MST")
data$jul <- as.numeric(strftime(data$date, format = "%j", tz = "MST"))
data$id_yr_seas_jul_ampm <- paste(data$id_yr_seas, data$jul, data$ampm, sep="_")
tokeep <- data
tokeep <- tokeep[order(sample(nrow(tokeep),nrow(tokeep),replace=FALSE)),] #randomize the databse
data <- tokeep[duplicated(tokeep$id_yr_seas_jul_ampm)==FALSE,] #grab the first non duplicated from each id_yr_seas_jul value
rm(tokeep)
head(data)

# you should also make this a balance design with respect to the number of points you have per ID/data/year
table(data$id_yr_seas)
min(table(data$id_yr_seas))
minnum <- 154   #identify this number from your table(data$id_yr_seas)
u <- unique(data$id_yr_seas)
data <- do.call(rbind, lapply(1:length(u), function(i){
  temp <- data[data$id_yr_seas == u[i],]
  return(temp[1:minnum,])
}))
data <- data[order(data$id, data$date),]  #reorder data
rm(u, minnum)
table(data$id_yr_seas)   #now each id_yr_seas should have teh same number of points


#-----------------------------------------#
#--- Biased Correlated random walk ----####
#-----------------------------------------#

# Biased correlated random walk (BCRW) model... estimates the following parameters:
# beta = weight of influence of talc, which is turning angles towards some land cover feature
# rho = rho paramter of the wrappedcauchy distribution - or strength of directional persistence (closer to 0, no directional persistence, close to 1 = complete directional persistence)
# shape = shape paramter of the weibull distribution of your step lengths
# scale = scale paramter of the weibull distribution of your step lengths

# The model takes the following data:
# TAs = observed turning angles in radians
# talc = angles towards some landscape feature in radians
# taprev = turning angle of previous step in radians
# Steps = step lengths for each step

# this is the model that jags will call
sink(file="bcrw.txt")
cat("
    model{
    
    # priors on movement parameters 
    scale ~ dnorm(0.0, 0.01)T(0.0,)    
    shape ~ dnorm(0.0, 0.01)T(0.0,)
    beta ~ dunif(0,1)
    #beta ~ dnorm(0, 0.01)T(0.0,) 
    rho ~ dunif(0,1)
    
    Pi <- 3.14159265359
    
    for (t in 1:n) {
    
    # likelihood for steps
    steps[t] ~ dweib(shape, scale)  
    
    # likelihood for turns (wrapped cauchy distribution)
    ya[t] <- (1-beta)*sin(taprev[t]) + beta*sin(talc[t])
    xa[t] <- (1-beta)*cos(taprev[t]) + beta*cos(talc[t])  
    mu[t] <- 2*atan(ya[t]/(sqrt(pow(ya[t],2) + pow(xa[t],2))+xa[t]))
    ones[t] ~ dbern( wc[t] )
    wc[t] <- (1/(2*Pi)*(1-pow(rho,2))/(1+pow(rho,2)-2*rho*cos(mu[t]-TAs[t])))/ 100
    }    
    }", fill=TRUE)
sink()


#change step length distribution to kms, not meters
hist(data$dist)
data$dist <- data$dist/1000
hist(data$dist)


# Bundle data for sending over to jags
dat <- list(steps=data$dist,  # step lengths for each step 
            n=nrow(data),  # sample size
            ones=rep(1,nrow(data)), # this is a vector of 1s for jags (it needs it for some weird reason)
            taprev=data$ta_prev,  # turning angle of previous step in radians
            talc=data$dir2steep,  # angles towards some landscape feature in radians
            TAs=data$ta)   # observed turning angles in radians

# what values to monitor during the parameterization
# beta = weight of influence of talc, which is turning angles towards some land cover feature
# rho = rho paramter of the wrappedcauchy distribution - or strength of directional persistence (closer to 0, no directional persistence, close to 1 = complete directional persistence)
# shape = shape paramter of the weibull distribution of your step lengths
# scale = scale paramter of the weibull distribution of your step lengths
pars <- c("scale",
          "shape",
          "rho", 
          "beta")

# Bayesian parameterization params
n.adapt <- 1000   # burn in - how many iterations to explore parameter space at beginning (may need to increase this if your model does not converge)
n.iter <- 2000    # number of iterations for drawing posterior samples (if you increase thinning, you also have to increase this number)
thin <- 2         # how much thinning should go on for selecting posterior samples (may need to increase this if your model does not converge)
chains <- 1:3     # how many chains to run (most bayesian statisticians do 3)

# if necessary.
# inits <- function()list(
#   scale=rnorm(1,2,.01), shape=rnorm(1,1,01), 
#   beta=runif(1,0.01,0.99), rho=runif(1,0.01,0.99))

# took about 2 minutes for ~2300 data points

# we will fit using parallel processing across the 3 chains

# identify cores (use 1 less than you have)
no_cores <- ifelse(detectCores()-1 < max(chains),detectCores()-1, max(chains)) 
# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("dat","pars","n.adapt","n.iter","thin"))

# now for the actual loop 9:02
pstemp <- clusterApplyLB(clust, chains, function(i){
  # need library() here for the packages your calculations require for your calculations
  library(rjags)
  
  # parametrize the model, or burn it in
  jm <- jags.model("bcrw.txt", data=dat, n.chains=1, n.adapt=n.adapt,
                   inits=list(.RNG.name="base::Wichmann-Hill",
                              .RNG.seed=i))
  # Get Posterior samples
  return(coda.samples(jm, pars, n.iter=n.iter, thin=thin))
  })
stopCluster(clust)   # you must stop the parallelization framework


# Reorganize pstemp to bring the chains together into a single 'mcmc.list' object
ps <- lapply(chains, function(i){
  return(pstemp[[i]][[1]])
})
class(ps) <- "mcmc.list"
str(ps)
rm(pstemp)

# take a look at the results
par(mai=c(.5,.5,.5,.5))
plot(ps, ask=T)
sum <- summary(ps)
sum    #check to see if the 95% credible interval for Beta overlaps 0. If not, your species' movements are biased towards talc. 
# And the size of beta (from 0 to 1) provides info on how strong the bias is over persistence
# also if beta is > 0.5, then selection towards talc is stronger than your directional persistence
#check also to see if the 95% credible interval for rho overlaps 0. If not, you have directional persistence.


#check out convergence of model
autocorr.plot(ps) #if there is high correlation within chains, you should increase thin value.
gelman.diag(ps) #estimates should be close to 1 if converged

# plot step length distribution
par(mfrow=c(2,2))
hist(dat[["steps"]], freq=FALSE, main="Observed", 
     xlab="Step length (km)",xlim=c(0,max(dat[["steps"]])))
hist(rweibull(10000, sum$statistics["shape","Mean"], 
              sum$statistics["scale","Mean"])/10, xlim=c(0,max(dat[["steps"]])),
     main="Predicted",freq=FALSE,xlab="Step length (km)")

#plot your directional persistence...
distangle <- circular(rwrappedcauchy(1000, rho=sum$statistics["rho","Mean"],
                                     control.circular=list(type="angles",units="radians", 
                                                           modulo="2pi", rotation="clock")))
# note that the prediction assumes that your mean direction of persistence is straight ahead (radian around 0)
# sometimes that is not the case, and you'd see that when comparing these two plots
h <- hist(rad(data$rel.angle), main="Observed", 
          xlab="Turning angles (rad)", freq=FALSE)
lines(smooth.spline(h$mids, h$density+0.1*mean(h$density)), col="blue")
h <- hist(as.numeric(distangle), main="Predicted", 
          xlab="Turning angles (rad)", freq=FALSE)
lines(smooth.spline(h$mids, h$density+0.1*mean(h$density)), col="blue")

# plot rose diagram of your predicted turning angles (Note: this is boring if your rho is small)
dev.off()
rose.diag(distangle, shrink=1.5, bins=26, col="orange", tol=.1, zero=pi/2,rotation="clock")
points(distangle, col="blue", pch=1,zero=pi/2,rotation="clock", stack=TRUE)


