#-------------------------------------------#
#- Resource Selection Functions ---------####
#------------- Lab 4 -----------------------#
#------- Part a - sampling available points-#
#---------- Jerod Merkle -------------------#


library(sf)
library(raster)
library(mapview)
library(tidyverse)
library(adehabitatHR)
library(stringr)
library(move)
library(dummies)
library(lme4)
library(parallel)


#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#
data <- readRDS("./data/bhs_data_summer.rds")
head(data)


#---------------------------------------------------#
# load up UDs for creating availability polygons ####
#---------------------------------------------------#

fls <- dir("./BBs", full.names = TRUE)  #this is the folder where your BBs are stored
#keep only dyn Brownian bridge files
fls <- grep("dynBB",fls, value=TRUE)

# get the names out of the file info
names <- fls %>% 
  str_split_fixed("/", str_count(fls[1], "/")+1) %>% 
  as.data.frame() %>% 
  pull(str_count(fls[1], "/")+1) %>% 
  str_replace("dynBB_","") %>% 
  str_replace(".tif","")
names

# Did all your id_yr_seas give you a home range?
table(names %in% data$id_yr_seas)   # should all be TRUE, if not, then might want to revisit home range analyses.
table(unique(data$id_yr_seas) %in% names)   # this should be all true as well

# what is the volume value (i.e., contour) you'd like to create your home ranges from?
# If using BB model, you may want to use a larger value so it 'captures' more area just outside of the use
# so you get a good representations of unused areas
# if using kernel methods, you may want something like a 0.95 percentile as those methods inherently 
# provide a broader UD.
percentile <- 0.999    # what contour percentile will you use? this is 99.9%

#loop over files, and create a list of the UDs in raster format ####

# prepare for parallel processing:
# identify cores (use 1 less than you have)
no_cores <- detectCores()-1
# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("fls", "percentile", "names"))

# now for the actual loop
HRs <- do.call(rbind, clusterApplyLB(clust, 1:length(fls), function(i){
  # need library() here for the packages your calculations require for your calculations
  library(raster)
  library(sf)
  rast <- raster(fls[i])
  vls <- sort(values(rast), decreasing=TRUE)
  vlscsum <- cumsum(vls)
  cutoff <- vls[vlscsum > percentile][1]
  polyR <- reclassify(rast, rcl=matrix(c(0,cutoff,0,cutoff,1,1),2,3, byrow=T))
  polyR[values(polyR)<=0] <- NA
  polyR <- rasterToPolygons(polyR, dissolve=TRUE) #this is the spatialPolygonsDataFrame representing the contour
  polyR <- as(polyR, "sf")    #turn into a sf object
  names(polyR) <- c("id_yr_seas","geometry")
  polyR$id_yr_seas <- names[i]
  return(polyR)
}))
stopCluster(clust)   # you must stop the parallelization framework
rm(no_cores, clust)

head(HRs)
rm(fls, names, percentile)
plot(HRs)
plot(HRs[3,]) # plot a single ID

#plot the HRs all together
plot(extend(extent(data),10000))  
plot(HRs[,1], add=T, col=NA)  #this plots all of your HR polygons

#plot each one with the points to verify that you like what your HRs look like
whichID <- sample(HRs$id_yr_seas,1)    #choose an ID to look at
plot(st_geometry(HRs[HRs$id_yr_seas == whichID,1]), main=whichID, col="lightgrey")
plot(st_cast(st_combine(data[data$id_yr_seas == whichID,]), "LINESTRING"), add=T, col="grey")
plot(st_geometry(data[data$id_yr_seas == whichID,]), add=T, pch=1, col="blue")

#--------------------------#
# population polygon HR ####
#--------------------------#
# merge all the HRs into a single 'population' polygon
Pop_poly <- st_union(HRs)
Pop_poly <- st_transform(Pop_poly, crs=st_crs(data))   # make sure it is in correct projection

# Alternatively, if your HRs are disjoint and you want a more continuous polygon, do the following
Pop_poly <- mcp(SpatialPoints(st_coordinates(st_union(HRs))[,1:2]), percent=100, unin="m")
proj4string(Pop_poly) <- st_crs(HRs)$proj4string
Pop_poly <- st_union(as(Pop_poly, "sf"))

# check and make sure everything looks OK
plot(Pop_poly) 
plot(sample_n(data, 1000)[,1], add=T, col="darkgreen", pch=".")

# this is a nice plot to see where individuals' HRs overlap more or less
plot(Pop_poly, border="red") 
plot(st_geometry(HRs), col=adjustcolor("grey", alpha.f=0.1), 
     border=NA, add=T)   #this plots all of your HR polygons with some transparency


#------------------------------------------------#
# generate random points within Population HR ####
#------------------------------------------------#
#general equal number of points as in data
PopRand <- st_sample(Pop_poly, size=nrow(data)+1000, type="random")
PopRand <- PopRand[1:nrow(data)]   # you have to do this because st_sample doesn't give the exact size you are looking for (might need to change the 1000 above to a bigger number!!!!)
nrow(data)==length(PopRand)   #this MUST BE TRUE (otherwise increase the 1000 number above)

#have a quick look at what you did
plot(Pop_poly)
plot(PopRand[sample(length(PopRand),3000)], add=T, col="blue", pch=".")   #sample of available points
plot(sample_n(data,3000)[,1], add=T, col="orange", pch=".")         #sample of used points

#attach points to dataframe
data$used <- 1  #add response variable
#turn random points into Spatial Points dataframe with the correct column names and attributes
head(data)
PopRand <- st_sf(data.frame(used=0,   # flag random points as 0s in used column
                            id_yr_seas=data$id_yr_seas,   # give random points a random id_yr_seas
                            date=data$date, # give the random points a random date as well
                            id=data$id), # give random id
                 geometry=PopRand)
PopRand <- st_transform(PopRand, crs=st_crs(data))  # verify same projection
data_pop <- rbind(data, PopRand) #rbind up the used and random points
rm(PopRand)
nrow(data_pop)
table(data_pop$used)
head(data_pop)
tail(data_pop)


#-------------------------------------------------------#
# generate random points within each individual's HR ####
#-------------------------------------------------------#
HRRand <- do.call(rbind, lapply(1:nrow(HRs), function(i){
  #general equal number of points as in data
  IdRand <- st_sample(HRs[i,], size=nrow(data[data$id_yr_seas==HRs$id_yr_seas[i],])+1000, type="random")
  IdRand <- IdRand[1:nrow(data[data$id_yr_seas==HRs$id_yr_seas[i],])]   # you have to do this because st_sample doesn't give the exact size you are looking for (might need to change the 1000 above to a bigger number!!!!)
  print(nrow(data[data$id_yr_seas==HRs$id_yr_seas[i],])==length(IdRand))   #this MUST BE TRUE (otherwise increase the 1000 nbumber above)

  # turn random points into sf dataframe with the correct attributes
  return(st_sf(data.frame(used=0, 
                          id_yr_seas=HRs$id_yr_seas[i], 
                          date=data$date[data$id_yr_seas==HRs$id_yr_seas[i]],   # give the random points a random date as well
                          id=str_split_fixed(HRs$id_yr_seas[i],"_",2)[,1]), 
               geometry=IdRand))
}))
HRRand <- st_transform(HRRand, crs=st_crs(data))  # make sure the crs is correct.
#They ALL must be TRUE!!!!

#take a look at the points you just made within the population polygon
head(HRRand)
plot(Pop_poly)
plot(st_geometry(HRRand), add=T, col="blue", pch=".")

# now take alook at them for each iD
whichID <- sample(HRs$id_yr_seas,1)
plot(st_geometry(HRs[HRs$id_yr_seas == whichID,1]), main=whichID, col="lightgrey")
plot(st_geometry(HRRand[HRRand$id_yr_seas == whichID,]), pch=".",col="blue", add=T)
plot(st_geometry(data[data$id_yr_seas == whichID,]), pch=".", add=TRUE,col="orange")
rm(whichID)

#rbind up the used and random points
data_HR <- rbind(data, HRRand) 
rm(HRRand)
head(data_HR)
tail(data_HR)

#-------------------------------------------------#
# generate random points within a local buffer ####
#-------------------------------------------------#
# first, figure out 80% quantile of distances moved during the season of interest as the buffer
head(data)
# calculate movement parameters #
source("./MiscFunctions/CalcBurst.R")
source("./MiscFunctions/CalcMovParams.R")
data <- data[order(data$id_yr_seas, data$date),] #order database first
data$burst <- CalcBurst(data=data, id = TRUE, id_name="id_yr_seas", date_name="date", Tmax = 3600*7) #set Tmax to 7 hours, a little more than double the fix rate
length(unique(data$burst))
data <- CalcMovParams(data=data, id_name = "id_yr_seas", date_name = "date", burst=data$burst)

hist(data$rel.angle)   #this is your turning angle distribution
hist(data$dist/1000)   #this is your step length distribution
rm(CalcBurst, CalcMovParams)
head(data)


hist(data$dist)
buf <- as.numeric(quantile(data$dist,probs = .80, na.rm=TRUE)) # this is the buffer in meters
buf

# now loop over each line in your dataset to generate a matched random point within the buffer
# note that this one is going to take some time!
# might even want to subsample, if you have LOTS of data!!!!!
# sfInit(parallel = T, cpus = 7)   #must change the cpus
# sfExport("data", "buf")
# sfLibrary(sf)
# LocalRand <- do.call(c, sfClusterApplyLB(1:nrow(data), function(i){
#   bufr <- st_buffer(data[i,], dist=buf) # create buffer polygon around point
#   return(st_sample(bufr, size=5, type="random")[1]) #sample 1 point in that polygon
# }))
# sfStop()
# length(LocalRand) == nrow(data)
# add the other columns to the sfc object
# LocalRand <- st_sf(data.frame(id_yr_seas=data$id_yr_seas, used=0, date=data$date, id=data$id), geometry=LocalRand)


# Alternatively, use trig to calculate the same thing as code commented out above
#generate random angles
angles <- runif(nrow(data),0, 2*pi)
# generate random distances between 1 and buf
dists <- runif(nrow(data),1,buf)
# used triginomitry to calcualte the x and y for a point at a given angle and a given distance
x <- st_coordinates(data)[,1]+(cos(angles)*dists) # how far moved in x direction
y <- st_coordinates(data)[,2]+(sin(angles)*dists) # how far moved in y direction
#turn into a dataframe
LocalRand <- data.frame(id_yr_seas=data$id_yr_seas, used=0, date=data$date, id=data$id, x=x, y=y)
head(LocalRand)
# turn into a sf object
LocalRand <- st_as_sf(LocalRand, coords=c("x","y"), dim="XY", crs=st_crs(data))

nrow(LocalRand) == nrow(data) #should be TRUE
rm(buf, x, y, angles, dists)
# remove the move params
data$burst <- NULL
data$dist <- NULL
data$dt <- NULL
data$speed <- NULL
data$abs.angle <- NULL
data$rel.angle <- NULL
data$StepFlag <- NULL
head(data)
data_local <- rbind(data, LocalRand) #rbind up the used and random points
rm(LocalRand)
head(data_local)
tail(data_local)
nrow(data_local)

# now that you've sampled availability, write out your environment! ####
# save your environment so we can load it back up another time ####
save.image("RSFs_AvailSampled.RData")
# load("RSFs_AvailSampled.RData")
