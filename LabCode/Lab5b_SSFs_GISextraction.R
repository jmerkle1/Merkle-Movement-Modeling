#--------------------------------------#
#- Step Selection Functions --------####
#--- extract GIS data to steps --------# 
#---------- Jerod Merkle --------------#
#------------- Lab 5b ------------------#

library(sf)
library(raster)
library(mapview)
library(tidyverse)
library(stringr)
library(parallel)
library(snow)
library(circular)
library(MASS)


#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#---------------------------------------------#
# Load SSF data with availability sampled  ####
#---------------------------------------------#
# from lab SSF part a
load("SSFs_AvailSampled.RData")


#--------------------------------------#
# reduce database to 1 step per day ####
#--------------------------------------#
# to reduce pseudoreplication (you don't have to do this, or you can do something different too)
# also doing this increases speed to fit the models later
# get id_yr_seas column fixed, so there are no NAs in it
data <- data[order(data$id_yr_seas, data$date, data$case),]

# grab one step per day, per ID, per year, per season...
data$jul <- as.numeric(strftime(data$date, format = "%j", tz = attributes(data$date[1])$tzone))
data$id_yr_seas_jul <- paste(data$id_yr_seas, data$jul, sep="_")
tokeep <- data[data$case==1,c("strata","id_yr_seas_jul")]
tokeep <- tokeep[order(sample(nrow(tokeep),nrow(tokeep),replace=FALSE)),] #randomize the tokeep database
data <- data[data$strata %in% tokeep$strata[duplicated(tokeep$id_yr_seas_jul)==FALSE],]
data$strata <- as.numeric(as.factor(data$strata))
data$id_yr_seas_jul <- NULL   # remove the id_yr_seas column
data$jul <- NULL   # remove the id_yr_seas column
rm(tokeep)
head(data, 10)

length(unique(data$strata))  # how many used steps are you left with? I have 1,816
# this gives you number of unique steps per id_yr_seas. I have about 122 each.
tapply(data$strata, data$id_yr_seas, FUN=function(x){length(unique(x))})


# -----------------------------------#
# Make data spatially aware again ####
# -----------------------------------#
head(data)
names(data)[names(data)%in%c("x","y")] <- c("x_orig","y_orig")  #change x/y names, so we have x_orig nad x_end of the step
#turn into point dataframe and make spatially aware
data <- data[order(data$id_yr_seas, data$date, data$case),]
table(is.na(data$x_end))  #make sure this is all FALSE

# create a separate sf object for the source point (i.e., orig) and target point (i.e., end)
#for the end of steps
sf_end <- st_as_sf(data, coords=c("x_end","y_end"), dim="XY", crs=proj)
#for the start of steps (we'll need later)
sf_orig <- st_as_sf(data,coords=c("x_orig","y_orig"), dim="XY", crs=proj)

#---------------------------#
# bring in your GIS data ####
#---------------------------#
elev <- raster("./data/GIS/Elevation_meters_30m.tif")
trasp <- raster("./data/GIS/Aspect_TRASP_30m.tif")
slope <- raster("./data/GIS/Slope_degrees_30m.tif")
tpi <- raster("./data/GIS/TPI_unitless_30m.tif")
Dist2Escape <- raster("./data/GIS/Dist2EscapeTerrain_meters_30m.tif")

landcov <- raster("./data/GIS/LandCover_descrete_30m.tif")
Dist2roads <- raster("./data/GIS/Dist2roads_meters_90m.tif")

treecov <- raster("./data/GIS/TreeCover_percent_30m.tif")
iNDVI18 <- raster("./data/GIS/iNDVI18.tif")
iNDVI19 <- raster("./data/GIS/iNDVI19.tif")
Cover_AnnualForbGrass <- raster("./data/GIS/RAP_2019_Cover_AnnualForbsGrasses.tif")
Cover_BareGround <- raster("./data/GIS/RAP_2019_Cover_BareGround.tif")
Cover_PerennialForbGrass <- raster("./data/GIS/RAP_2019_Cover_PerennialForbsGrasses.tif")

# ------------------------------------#
# relate target points to GIS data ####
# ------------------------------------#

# sf_end or target points first (all variables)
data$elev_target <- raster::extract(elev, st_transform(sf_end, crs=projection(elev)))
data$trasp_target <- raster::extract(trasp, st_transform(sf_end, crs=projection(trasp)))
data$slope_target <- raster::extract(slope, st_transform(sf_end, crs=projection(slope)))
data$tpi_target <- raster::extract(tpi, st_transform(sf_end, crs=projection(tpi)))
data$Dist2Escape_target <- raster::extract(Dist2Escape, st_transform(sf_end, crs=projection(Dist2Escape)))
data$lc_target <- raster::extract(landcov, st_transform(sf_end, crs=projection(landcov)))
data$Dist2roads_target <- raster::extract(Dist2roads, st_transform(sf_end, crs=projection(Dist2roads)))
data$treecov_target <- raster::extract(treecov, st_transform(sf_end, crs=projection(treecov)))
data$iNDVI19_target <- raster::extract(iNDVI19, st_transform(sf_end, crs=projection(iNDVI19)))
data$AnnualForbGrass_target <- raster::extract(Cover_AnnualForbGrass, st_transform(sf_end, crs=projection(Cover_AnnualForbGrass)))
data$BareGround_target <- raster::extract(Cover_BareGround, st_transform(sf_end, crs=projection(Cover_BareGround)))
data$PerennialForbGrass_target <- raster::extract(Cover_PerennialForbGrass, st_transform(sf_end, crs=projection(Cover_PerennialForbGrass)))

head(data)

# you can also extract GIS data to your source points. 
# You may or may not want to do this, depending on your questions/goals. 

#-----------------------------#
# relate lines to GIS data ####
#-----------------------------#
# first thing is to create sf lines dataframe with each row representing a step (used or available)

dataL <- lapply(1:nrow(data), function(i){
  return(st_linestring(rbind(as.matrix(data[i,c("x_orig","y_orig")]),
               as.matrix(data[i,c("x_end","y_end")])),dim="XY"))
})
dataL <- st_as_sfc(dataL, crs=proj)
dataL <- data.frame(strata=data$strata, case=data$case, geometry=dataL)
dataL <- st_as_sf(dataL, sf_column_name = "geometry")
head(dataL)
nrow(data); nrow(dataL)  # these of course should be exactly the same!

# check to make sure each dataset has the same order of the rows (use strata and case to determine this)
head(data)
head(sf_end)
head(sf_orig)
head(dataL)


#take a look at the steps you have generated ####
smpl <- sample(unique(data$strata),1)   #grab a random strata to plot
plot(dataL$geometry[data$strata==smpl,], main=paste("Strata:", smpl))
points(data$x_orig[data$strata==smpl], data$y_orig[data$strata==smpl], col="blue", pch=16)
points(data$x_end[data$strata==smpl], data$y_end[data$strata==smpl], col="orange", pch=16)
points(data$x_end[data$strata==smpl&data$case==1], data$y_end[data$strata==smpl&data$case==1], col="darkgreen", pch=16)
legend("topright", c("Source","Avail","Used"), pch=16, col=c("blue","orange","darkgreen"), bty="n")

# now do the same thing with mapview
smpl <- sample(unique(data$strata),1)   #grab a random strata to plot
mapview(dataL$geometry[data$strata==smpl,], layer.name=paste("Strata:", smpl),
        map.types="Esri.WorldImagery")+
  mapview(sf_end[data$strata==smpl,], zcol="case", layer.name="Case")

rm(smpl)
  
#----------------------------------------#
# Now to relate the lines to GIS data ####
#----------------------------------------#

# you might consider trying your code on a small subset of lines first (make sure its works)
raster::extract(elev,st_transform(dataL[1:25,],crs=projection(elev)), fun=mean) #extracts the cell values that the line touches, then takes the mean of those values

# now to do all of your continuous data!
head(data)
beginCluster(type="SOCK")   #use multi core processing because it takes a while to assess what cells each line are touching

# Note: These take about 20 minutes per line / extraction (I have 11 cores, 30m res rasters, and 11k lines)
Sys.time()
data$elev_step <- raster::extract(elev, st_transform(dataL, crs=projection(elev)), fun=mean)
Sys.time()
data$trasp_step <- raster::extract(trasp, st_transform(dataL, crs=projection(trasp)), fun=mean)
Sys.time()
data$slope_step <- raster::extract(slope, st_transform(dataL, crs=projection(slope)), fun=mean)
data$tpi_step <- raster::extract(tpi, st_transform(dataL, crs=projection(tpi)), fun=mean)
data$Dist2Escape_step <- raster::extract(Dist2Escape, st_transform(dataL, crs=projection(Dist2Escape)), fun=mean)
data$Dist2roads_step <- raster::extract(Dist2roads, st_transform(dataL, crs=projection(Dist2roads)), fun=mean)
data$treecov_step <- raster::extract(treecov, st_transform(dataL, crs=projection(treecov)), fun=mean)
data$iNDVI19_step <- raster::extract(iNDVI19, st_transform(dataL, crs=projection(iNDVI19)), fun=mean)
data$AnnualForbGrass_step <- raster::extract(Cover_AnnualForbGrass, st_transform(dataL, crs=projection(Cover_AnnualForbGrass)), fun=mean)
data$BareGround_step <- raster::extract(Cover_BareGround, st_transform(dataL, crs=projection(Cover_BareGround)), fun=mean)
data$PerennialForbGrass_step <- raster::extract(Cover_PerennialForbGrass, st_transform(dataL, crs=projection(Cover_PerennialForbGrass)), fun=mean)
endCluster()
head(data, 3)

# have to do some different things for descrete variables like landcover...
# basically, for each landcover variable, you need to create a seaprate raster of 
# whether or not it is that landcover (i.e., 1s and 0s). Then, when you extract
# to lines, you get the proportion of the line (i.e., step) that crosses a given landcover type
# my suggestion is to pick most important landcover variables here
unique(landcov)   # this is nlcd landcov
table(data$lc_target)   # you may need to bring in the legend to interpret
table(data$landcov_target)   # you may need to bring in the legend to interpret
barren <- landcov  #make empty landcov database called barron
barren[] <- ifelse(values(landcov)==31, 1, 0)   #reclassify values so barren landcovers = 1 and all else = 0
shrub <- landcov
shrub[] <- ifelse(values(landcov)==52, 1, 0)
herb <- landcov
herb[] <- ifelse(values(landcov)%in% 71:74, 1, 0)
plot(barren)
plot(shrub)
plot(herb)

beginCluster(type="SOCK")  
data$barren_step <- raster::extract(barren, st_transform(dataL, crs=projection(barren)), fun=mean)
data$shrub_step <- raster::extract(shrub, st_transform(dataL, crs=projection(shrub)), fun=mean)
data$herb_step <- raster::extract(herb, st_transform(dataL, crs=projection(herb)), fun=mean)
endCluster()
head(data)

rm(elev, trasp, slope, tpi, Dist2Escape, landcov,
   Dist2roads,treecov,iNDVI18,iNDVI19,Cover_AnnualForbGrass,
   Cover_BareGround,Cover_PerennialForbGrass,herb,shrub,barren)


# --------------------------------------------------------------#
# Calc whether a step crossed a linear feature, and how many ####
# --------------------------------------------------------------#
#load up your road layer
roads <- st_read("./data/GIS", "roads_tiger_2016")
crosses <- st_crosses(dataL, roads, sparse=FALSE)  # matrix of steps and the different roads and whether or not a cross occured
dim(crosses)
data$numb_roads_crossed <- apply(crosses, 1, sum)  # sum up across the different roads
hist(data$numb_roads_crossed)
table(data$numb_roads_crossed)
rm(roads, crosses)

# now that you've extracted your GIS data, write out your environment! ####
# save your environment so we can load it back up another time
save.image("SSFs_VariablesExtracted.RData")
# load("SSFs_VariablesExtracted.RData")


