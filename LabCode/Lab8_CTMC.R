#----------------------------------------#
#-Continuous time Markov chain models-####
#---------- Jerod Merkle ----------------#
#------------- Lab 8 --------------------#

library(raster)
library(sf)
library(mgcv)
library(ctmcmove)
library(parallel)
library(lme4)

#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#
data <- readRDS("./data/bhs_data_summer.rds")
head(data)

#-------------------------------------------------#
# create a raster stack of all your covariates ####
#-------------------------------------------------#

# load up all the covariates in question
elev <- raster("./data/GIS/Elevation_meters_30m.tif")
trasp <- raster("./data/GIS/Aspect_TRASP_30m.tif")
slope <- raster("./data/GIS/Slope_degrees_30m.tif")
tpi <- raster("./data/GIS/TPI_unitless_30m.tif")

landcov <- raster("./data/GIS/LandCover_descrete_30m.tif")
Dist2roads <- raster("./data/GIS/Dist2roads_meters_90m.tif")

iNDVI19 <- raster("./data/GIS/iNDVI19.tif")
Cover_AnnualForbGrass <- raster("./data/GIS/RAP_2019_Cover_AnnualForbsGrasses.tif")
Cover_BareGround <- raster("./data/GIS/RAP_2019_Cover_BareGround.tif")
Cover_PerennialForbGrass <- raster("./data/GIS/RAP_2019_Cover_PerennialForbsGrasses.tif")

# if you have different projections/extents, you'll need to pick one and project the others to that projection/ext
# NOTE, be careful with projectRaster with big rasters. Can take a LONG time!
# Also, you should use a projection with meters as the unit (i.e., no latlong)
compareRaster(elev, trasp, slope,tpi, extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, orig=TRUE)
# the above are all same extent for me
compareRaster(elev, Dist2roads, extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, orig=TRUE)
# I'm gonna reprojecdt everything to my dist2road layer, which is in Albers Equal Area, units meters.
elev <- projectRaster(elev, Dist2roads, method="bilinear")
trasp <- projectRaster(trasp, Dist2roads, method="bilinear")
slope <- projectRaster(slope, Dist2roads, method="bilinear")
tpi <- projectRaster(tpi, Dist2roads, method="bilinear")
iNDVI19 <- projectRaster(iNDVI19, Dist2roads, method="bilinear")
Cover_AnnualForbGrass <- projectRaster(Cover_AnnualForbGrass, Dist2roads, method="bilinear")
Cover_BareGround <- projectRaster(Cover_BareGround, Dist2roads, method="bilinear")
Cover_PerennialForbGrass <- projectRaster(Cover_PerennialForbGrass, Dist2roads, method="bilinear")

#reprojecting a raster with categorical variables
forest <- landcov  #make empty landcov raster called forest
forest[] <- ifelse(values(landcov)%in%c(41,42,43), 1, 0)   #reclassify values so barren landcovers = 1 and all else = 0
plot(forest)   # the colors might be screwed up (mine was all black)
table(values(forest))   # but it worked
forest <- projectRaster(forest, Dist2roads, method="ngb")
plot(forest, col=heat.colors(2))


#--------------------------------------------#
# load up your raster stack of covariates ####
#--------------------------------------------#

c.stack <- stack(list(elev,trasp,slope,tpi,Dist2roads,iNDVI19,forest,
                      Cover_AnnualForbGrass,Cover_BareGround,Cover_PerennialForbGrass))  #this this code doesn't work that means one or more of your rasters do not have the same projection/ext, etc.

rm(elev,trasp,slope,tpi,Dist2roads,iNDVI19,
   Cover_AnnualForbGrass,Cover_BareGround,Cover_PerennialForbGrass, landcov)   # remove all your rasters
nlayers(c.stack)

int <- c.stack[[1]] # Create a flat raster to represent the intercept
int[] <- 1

c.stack <- stack(int, c.stack) # add intercept layer to the stack 
nlayers(c.stack)
nms <- c("intercept","elevation","trasp","slope","tpi","dist2roads","iNDVI","forest",
         "AnnualForbGrass","BareGround","PerennialForbGrass")    #these names MUST line up with your variables in c.stack
names(c.stack) <- nms
rm(int)

# make a stack of variables that you want to be included as gradients in the model
grad.stack <- stack(c(c.stack[["dist2roads"]],c.stack[["slope"]]))   # I used  dist2roads and slope (I think you can only use continuous variables for this)
plot(grad.stack)
names(grad.stack) <- c("dist2roads.grad","slope.grad")

# -------------------------#
# Prepare your GPS data ####
# -------------------------#
data <- data[order(data$id_yr_seas, data$date),]

# you need to reproject your data to the projection of your rasters
data <- st_transform(data, crs=projection(c.stack[[1]]))
data$x <- st_coordinates(data)[,1]
data$y <- st_coordinates(data)[,2]
data <- st_drop_geometry(data)   #bring back to a data.frame
proj <- projection(c.stack[[1]])   #cop your proj so you don't get confused.
head(data)

AID.list <- list()   #you are going to develop a list, where each element is an individual id_yr_seas
ids <- unique(data$id_yr_seas)
length(ids)
for(i in 1:length(ids)){
  idx <- data[data$id_yr_seas == ids[i],]
  ## turn time into a numeric value (in days since Jan 1, 1970)
  t <- as.numeric(idx$date)/60/60/24
  ## get xy values for each time point
  xy <- idx[, c("x","y")]
  AID.list[[i]] <- data.frame(t = t, x = xy[, 1], y = xy[, 2])
}
names(AID.list) <- ids
rm(idx, xy, i, ids, t)
head(AID.list[[1]])   #this should only be t (for days since jan 1 1970) and x and y

#-------------------------#
# create the CTMC path ####
#-------------------------#

n.ids <- length(AID.list) # you are going to loop over each id and create the CTMC path

# we will use parallel processing
# NOTE HERE TO BE CAREFUL WITH MEMORY USAGE. MINE WHEN UP QUITE A BIT HERE. PLAN CPUS ACORDINGLY

# identify cores (use 1 less than you have)
no_cores <- 6
# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("AID.list","c.stack","n.ids"))

# now for the actual loop 
ctmc.list <- clusterApplyLB(clust, 1:n.ids, function(i){
  # need library() here for the packages your calculations require for your calculations
  library(ctmcmove)
  library(raster)
  
  return(path2ctmc(AID.list[[i]][, 2:3], #this grabs the x and y data
                   as.numeric(AID.list[[i]][, 1]), #this grabs the time column
                   method="LinearInterp",
                   c.stack)) # Linear interpolation is used to connect gps fixes)
})
stopCluster(clust)   # you must stop the parallelization framework
str(ctmc.list[[1]])  # see Value section of ?path2ctmc for details


# how to plot what you did
whichID <- 3   #select an animal id (in sequential order, not actual id name)
tmprast <- c.stack[[1]]    #grab the CTMC results
tmprast[] <- 0
tmprast[ctmc.list[[whichID]]$ec] <- ctmc.list[[whichID]]$rt
plot(tmprast, xlim=range(AID.list[[whichID]]$x), ylim=range(AID.list[[whichID]]$y), main="residency time (days)")   #this is the residence times for your individual
lines(AID.list[[whichID]]$x, AID.list[[whichID]]$y, col="grey", lwd=.1)
points(AID.list[[whichID]]$x, AID.list[[whichID]]$y, col="black")
hist(ctmc.list[[whichID]]$rt)   #this is the histogram of residence times per cell (in days)
rm(tmprast, whichID)


#-----------------------------------------#
# turn CTMC path into poisson glm data ####
#-----------------------------------------#

n.ids <- length(ctmc.list)

# identify cores (use 1 less than you have)
no_cores <- 6
# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("ctmc.list","c.stack","n.ids", "grad.stack"))

# now for the actual loop 
glm.list <- clusterApplyLB(clust, 1:n.ids, function(i){
  # need library() here for the packages your calculations require for your calculations
  library(ctmcmove)
  library(raster)
  
  return(ctmc2glm(ctmc.list[[i]], stack.static = c.stack, #this is a raster stack of your static variables
                  stack.grad=grad.stack, #this is a raster stack of your gradient variables (if you have any)
                  directions=8,   #if you have LOTS of data, might want to put directions to 4. This will reduce your final database
                  grad.point.decreasing = TRUE)) #when TRUE, means that lower values of the grad stack are positive or good for animal
  
})
stopCluster(clust)   # you must stop the parallelization framework

head(glm.list[[1]], 10)   #this is the database you made for each id_yr_seas
nrow(glm.list[[1]])

# add id column to resulting database
for(i in 1:n.ids){
  glm.list[[i]]$id <- rep(names(AID.list)[i], length(glm.list[[i]]$z)) # To add individual ids to each list item
}
head(glm.list[[1]])
rm(i, n.ids)

# rbind all the data together in one data frame
data <- do.call(rbind, glm.list)
head(data)  #this is your database now
nrow(data)

# you want to remove instantaneous transitions and overly long stays
hist(data$tau)
table(data$tau < 1/(24*60*60))   # anything less than 1 second is too fast of residency time in a cell (you can change this, if you want)

toremove <- data$tau< (1/(24*60*60)) | data$tau > 2 # Identifies too fast of residency times and overly long stays (i.e., > 2 days)
table(toremove)   #this is how many you are going to remove
data <- data[toremove == FALSE,] 
rm(toremove)
head(data)

#-------------------------------#
# Parameterize CTMC with GLM ####
#-------------------------------#

hist(data$crw)   #this is the distribution of the correlated random walk paramter (1 = forward movement, -1 = reverse movement)
hist(data$tau)   #this is the distribution of residency times (in days)
table(data$z)    # this is the number of used cells and their 7 potential cells
head(data, 2)

#start by checking correlation
vars <- c("elevation","trasp","slope","tpi","dist2roads","iNDVI","forest",
          "AnnualForbGrass","BareGround","PerennialForbGrass",'dist2roads.grad',"slope.grad")
correl <- round(cor(data[,vars],
                    use="pairwise.complete.obs", method="pearson"), 3) #pay attention to variables correlated > 0.5
ifelse(abs(correl)>.3, correl, NA)
rm(correl)

# change some of the variables to better scales?
data$dist2roads <- data$dist2roads/1000 #Dist2roads to km
data$elev <- data$elev/1000 # elev to KM

#add a weights column to the data to help with unbalanced data among individuals
weights <- data.frame(table(data$id))
names(weights) <- c("id","weights")
weights$weights <-weights$weights/nrow(data) 
data <- merge(data, weights)
table(data$weights)

#fit the model
test_mod <- glm(z~crw+elevation+trasp+slope+tpi+dist2roads+iNDVI+forest+AnnualForbGrass+BareGround+PerennialForbGrass+ 
                  dist2roads.grad+slope.grad, 
                offset = log(tau), weights = weights,
                family = poisson, data = data) 

# now this is a model of movement rate, so a positive coefficient means that they move faster when that variable is larger
summary(test_mod) # Summary of  coefficients


#fit with mixed effects (this didn't work for me)
# this is the mixed effects version with only a random intercept for id
table(data$id)
test_modme <- glmer(z~crw+elevation+trasp+slope+tpi+dist2roads+iNDVI+forest+AnnualForbGrass+BareGround+PerennialForbGrass+ 
                      dist2roads.grad+slope.grad + (0+slope|id), offset = log(tau),
                  family = poisson, data = data)

summary(test_modme)


#---------------------------------------------#
# created a predicted map of movement rate ####
#---------------------------------------------#

# need to scale variables in c.stack if you scaled them above
c.stack[["dist2roads"]] <- c.stack[["dist2roads"]]/1000
c.stack[["elevation"]] <- c.stack[["elevation"]]/1000

# Create a transition rate matrix
RR <- get.rate.matrix(test_mod, stack.static=c.stack, 
                      stack.grad=grad.stack, directions = 8,
                      grad.point.decreasing = TRUE) 

t <- transition(c.stack[[1]], mean, directions = 8, symm = FALSE) # Create dummy sparse transition matrix object to fill with rate matrix above
transitionMatrix(t) <- RR # Fill dummy matrix with rate matrix RR. 
t.rast <- raster(t) # Convert to a gridded conductance surface
rescaleCalc <- function (x) {((x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))} # Rescale to 0-1
t.rast <- calc(t.rast, rescaleCalc)
# t.rast[t.rast == 0] <- 0.001 # Convert zeros to small number

plot(t.rast) # t.rast can be used as a conductance surface (i.e., high values are fast movements, low values are slow movements)
points(data$x.adj[data$z ==1], data$y.adj[data$z ==1], pch=".")
