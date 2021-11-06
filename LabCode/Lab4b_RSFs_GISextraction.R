#--------------------------------------------#
#- Resource Selection Functions ----------####
#------------- Lab 4 ------------------------#
#------- Part b - extracting GIS to points --#
#---------- Jerod Merkle --------------------#

library(sf)
library(raster)
library(mapview)
library(tidyverse)
library(stringr)
library(move)
library(dummies)
library(lme4)
library(parallel)


#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#---------------------------------------------#
# Load RSF data with availability sampled  ####
#---------------------------------------------#
# from lab RSF part a
load("RSFs_AvailSampled.RData")

#------------------------------#
# relate GIS data to points ####
#------------------------------#
#bring in all your GIS data
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


# extract GIS data to points
# for HR scale first
data_HR$elev <- raster::extract(elev, st_transform(data_HR, crs=projection(elev)))
data_HR$trasp <- raster::extract(trasp, st_transform(data_HR, crs=projection(trasp)))
data_HR$slope <- raster::extract(slope, st_transform(data_HR, crs=projection(slope)))
data_HR$tpi <- raster::extract(tpi, st_transform(data_HR, crs=projection(tpi)))
data_HR$Dist2Escape <- raster::extract(Dist2Escape, st_transform(data_HR, crs=projection(Dist2Escape)))
data_HR$lc <- raster::extract(landcov, st_transform(data_HR, crs=projection(landcov)))
data_HR$Dist2roads <- raster::extract(Dist2roads, st_transform(data_HR, crs=projection(Dist2roads)))
data_HR$treecov <- raster::extract(treecov, st_transform(data_HR, crs=projection(treecov)))
data_HR$iNDVI18 <- raster::extract(iNDVI18, st_transform(data_HR, crs=projection(iNDVI18)))
data_HR$iNDVI19 <- raster::extract(iNDVI19, st_transform(data_HR, crs=projection(iNDVI19)))
data_HR$AnnualForbGrass <- raster::extract(Cover_AnnualForbGrass, st_transform(data_HR, crs=projection(Cover_AnnualForbGrass)))
data_HR$BareGround <- raster::extract(Cover_BareGround, st_transform(data_HR, crs=projection(Cover_BareGround)))
data_HR$PerennialForbGrass <- raster::extract(Cover_PerennialForbGrass, st_transform(data_HR, crs=projection(Cover_PerennialForbGrass)))
# for pop scale second
data_pop$elev <- raster::extract(elev, st_transform(data_pop, crs=projection(elev)))
data_pop$trasp <- raster::extract(trasp, st_transform(data_pop, crs=projection(trasp)))
data_pop$slope <- raster::extract(slope, st_transform(data_pop, crs=projection(slope)))
data_pop$tpi <- raster::extract(tpi, st_transform(data_pop, crs=projection(tpi)))
data_pop$Dist2Escape <- raster::extract(Dist2Escape, st_transform(data_pop, crs=projection(Dist2Escape)))
data_pop$lc <- raster::extract(landcov, st_transform(data_pop, crs=projection(landcov)))
data_pop$Dist2roads <- raster::extract(Dist2roads, st_transform(data_pop, crs=projection(Dist2roads)))
data_pop$treecov <- raster::extract(treecov, st_transform(data_pop, crs=projection(treecov)))
data_pop$iNDVI18 <- raster::extract(iNDVI18, st_transform(data_pop, crs=projection(iNDVI18)))
data_pop$iNDVI19 <- raster::extract(iNDVI19, st_transform(data_pop, crs=projection(iNDVI19)))
data_pop$AnnualForbGrass <- raster::extract(Cover_AnnualForbGrass, st_transform(data_pop, crs=projection(Cover_AnnualForbGrass)))
data_pop$BareGround <- raster::extract(Cover_BareGround, st_transform(data_pop, crs=projection(Cover_BareGround)))
data_pop$PerennialForbGrass <- raster::extract(Cover_PerennialForbGrass, st_transform(data_pop, crs=projection(Cover_PerennialForbGrass)))
# for local scale third
data_local$elev <- raster::extract(elev, st_transform(data_local, crs=projection(elev)))
data_local$trasp <- raster::extract(trasp, st_transform(data_local, crs=projection(trasp)))
data_local$slope <- raster::extract(slope, st_transform(data_local, crs=projection(slope)))
data_local$tpi <- raster::extract(tpi, st_transform(data_local, crs=projection(tpi)))
data_local$Dist2Escape <- raster::extract(Dist2Escape, st_transform(data_local, crs=projection(Dist2Escape)))
data_local$lc <- raster::extract(landcov, st_transform(data_local, crs=projection(landcov)))
data_local$Dist2roads <- raster::extract(Dist2roads, st_transform(data_local, crs=projection(Dist2roads)))
data_local$treecov <- raster::extract(treecov, st_transform(data_local, crs=projection(treecov)))
data_local$iNDVI18 <- raster::extract(iNDVI18, st_transform(data_local, crs=projection(iNDVI18)))
data_local$iNDVI19 <- raster::extract(iNDVI19, st_transform(data_local, crs=projection(iNDVI19)))
data_local$AnnualForbGrass <- raster::extract(Cover_AnnualForbGrass, st_transform(data_local, crs=projection(Cover_AnnualForbGrass)))
data_local$BareGround <- raster::extract(Cover_BareGround, st_transform(data_local, crs=projection(Cover_BareGround)))
data_local$PerennialForbGrass <- raster::extract(Cover_PerennialForbGrass, st_transform(data_local, crs=projection(Cover_PerennialForbGrass)))

#remove your rasters
rm(elev, trasp, slope, tpi, Dist2Escape, landcov,
   Dist2roads,treecov,iNDVI18,iNDVI19,Cover_AnnualForbGrass,Cover_BareGround,Cover_PerennialForbGrass)


# prepare data for analysis ####
#take databases out of sf format, and into a simple data frame.
data_HR <- st_set_geometry(data_HR, NULL)
data_pop <- st_set_geometry(data_pop, NULL)
data_local <- st_set_geometry(data_local, NULL)

#fix lc variable, so we have actual variable names for landcover
# landcover variables are often numeric, so want to merge with their meaning
legend <- read.csv("./data/GIS/nlcd_legend.csv") #bring in landcov legend
names(legend)[names(legend)=="value"] <- "lc"   # rename value column to match your landcover column name in datC1
names(legend)[names(legend)=="class_general"] <- "landcov"   # rename the class column to the name you want to see in datC1

data_HR <- merge(data_HR,legend[,c("landcov","lc")], all.x=TRUE)   # add a new column with actual landcover values
data_pop <- merge(data_pop,legend[,c("landcov","lc")], all.x=TRUE)
data_local <- merge(data_local,legend[,c("landcov","lc")], all.x=TRUE)
rm(legend)

head(data_pop)
head(data_HR)
head(data_local)

# now that you've extracted your GIS data, write out your environment! ####
# save your environment so we can load it back up another time 
save.image("RSFs_VariablesExtracted.RData")
# load("RSFs_VariablesExtracted.RData")