#--------------------------#
# GIS Data organization ####
#--------------------------#
#--- Jerod Merkle ---------#
#------ Lab 0b ------------#

#packages
library(sf)
library(raster)
library(mapview)
library(tidyverse)

#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#

data <- readRDS("./data/bhs_data.rds")
head(data)

#------------------------#
# get GIS stuff ready ####
#------------------------#

#make study area box/polygon from a buffer around your data
sa <- data %>% 
  extent() %>% 
  extend(c(10000,10000,10000,10000)) %>%  # add 10km buffer
  extent() %>% 
  st_bbox() %>% 
  st_as_sfc() %>% 
  st_sf(data.frame(id="studyarea")) #add a dataframe, so you can write it to file 
st_crs(sa) <- st_crs(data)   # set the projection

# have a look at the box you made
mapview(sa)
mapview(sa) + mapview(sample_n(data, 100))  # add a few random points to the plot

# write out your box
st_write(sa, dsn="./data/GIS", layer="study_area_box", driver="ESRI Shapefile")

#--------------------------#
# organize your rasters ####
#--------------------------#

# elevation
ele30 <- raster("T:/Elevation_30m/elevation30m.img")
ele30 <- crop(ele30, st_transform(sa, crs=projection(ele30)))  #crops your raster
#NOTE - never reproject a raster unless you absolutely have to!!!
plot(ele30)
plot(st_geometry(st_transform(sample_n(data, 100), crs=projection(ele30))), add=TRUE)

# write out the cropped raster. Note that I always write out my rasters with
# the name of the data, the units of the file, and the resolution. Good habit to have.
writeRaster(ele30,"./data/GIS/Elevation_meters_30m.tif", format = "GTiff")

# -----------------------------------------------#
# calculate other topo metrics from elevation ####
# -----------------------------------------------#
# slope, aspect, TRASP, TPI #
aspect <- terrain(ele30, opt = "aspect", unit="degrees")
plot(aspect)
writeRaster(aspect, "./data/GIS/Aspect_degrees_30m.tif", format="GTiff")

# turn aspect into trasp, which provides a circular aspect assigning a 
# value of zero to land oriented in a north-northeast direction, (typically the 
# coolest and wettest orientation), and a value of one on the hotter, dryer 
# south-southwesterly slopes. 
# See: https://rdrr.io/github/jeffreyevans/spatialEco/src/R/trasp.R

trasp_fun <- function(x) { (1 - cos( (3.142/180)  *(x - 30)) ) / 2 }
trasp <- calc(aspect, fun=trasp_fun)
plot(trasp)
plot(st_geometry(st_transform(sample_n(data, 100), crs=projection(trasp))), add=TRUE)
writeRaster(trasp, "./data/GIS/Aspect_TRASP_30m.tif", format="GTiff")
rm(trasp_fun, trasp)

#calculate slope
slope <- terrain(ele30, opt = "slope", unit="degrees")
plot(slope)
plot(st_geometry(st_transform(sample_n(data, 100), crs=projection(slope))), add=TRUE)
writeRaster(slope, "./data/GIS/Slope_degrees_30m.tif", format="GTiff")

#terrain position index (positive values are ridges, and negative values are valleys)
f <- matrix(1, nrow=9, ncol=9)   # you can change the scale of your TPI
TPI <- focal(ele30, w=f, fun=function(x, ...) x[41] - mean(x[-41]), pad=TRUE, padValue=NA) # this can take some time!
plot(TPI)
mapview(TPI)   # play around with this layer in mapview so you can see what it is doing if you have the right
# biological scale (i.e., is this scale pulling out large valleys like the platte river,
# or smaller inlet valleys, or even smaller draws)
writeRaster(TPI, filename="./data/GIS/TPI_unitless_30m.tif", format="GTiff")
rm(f, TPI)

# create a distance to escape terrain raster (mainly for bighorn sheep)
ele90 <- aggregate(ele30, fact=3)
slope <- terrain(ele90, opt = "slope", unit="degrees")
TRI <- terrain(ele90, opt = "TRI", unit="degrees")
plot(slope)
plot(TRI)   # flip back and forth and you can see that they are quite correlated
# escapeTerrain <- overlay(slope, TRI, fun=function(x,y){ifelse(x>30&y>75, 1, NA)})  # might want to change these numbers given your study area or paper you are basing off of.
escapeTerrain <- overlay(slope, TRI, fun=function(x,y){ifelse(x>27|y>75, 1, NA)}) # I had to change this... 
# wasn't much escape terrain in this study area. I defined as any area with > 27 degree slope (Smith et al. 1991, Papouchis et al. 2001)
plot(escapeTerrain)
dist2escapterr <- distance(escapeTerrain) # again, this takes time. 
plot(dist2escapterr)
writeRaster(dist2escapterr, filename="./data/GIS/Dist2EscapeTerrain_meters_30m.tif", format="GTiff")
rm(ele90, slope, TRI, escapeTerrain, dist2escapterr)

#hillshade ####
#hillshade is mainly used to make pretty maps (no used for a habitat selection model)
aspect <- terrain(ele30, opt = "aspect", unit="radians")
slope <- terrain(ele30, opt = "slope", unit="radians")
h <- hillShade(slope, aspect)
plot(h)  # gives the landscape a 3D effect with the topography
writeRaster(h, filename="./data/GIS/Hillshade_unitless_30m.tif", format="GTiff")
rm(h, aspect, slope, ele30)

#-------------------#
# landcover NLCD ####
#-------------------#
landcov <- raster("T:/Landcover_NLCD/NLCD_2016_Land_Cover_L48_20190424.img")
landcov <- crop(landcov, st_transform(sa, crs=projection(landcov)))
plot(landcov)
writeRaster(landcov, "./data/GIS/LandCover_descrete_30m.tif", format="GTiff")
rm(landcov)

canopycov <- raster("T:/Landcover_NLCD/nlcd_2016_treecanopy_2019_08_31.img")
canopycov <- crop(canopycov, st_transform(sa, crs=projection(canopycov)))
plot(canopycov)
writeRaster(canopycov, "./data/GIS/TreeCover_percent_30m.tif", format="GTiff")
rm(canopycov)

#--------------------#
#integrated NDVI  ####
#--------------------#
table(year(data$date)) # we want both 2018 (growing seasno prior to first winter), 2019, and 2020 (I don't have 2020 yet)
iNDVI <- stack("T:/MODIS_NDVI/Bischof_calculations/csumNDVImax.grd")
nlayers(iNDVI)
names(iNDVI)
iNDVI18 <- iNDVI[["X2018"]]   # pull out the years of interest
iNDVI19 <- iNDVI[["X2019"]]
iNDVI18 <- crop(iNDVI18, st_transform(sa, crs=projection(iNDVI18)))
iNDVI19 <- crop(iNDVI19, st_transform(sa, crs=projection(iNDVI19)))
plot(iNDVI18)
plot(st_geometry(st_transform(sample_n(data, 100), crs=projection(iNDVI18))), add=TRUE)

#write them to file
writeRaster(iNDVI18, "./data/GIS/iNDVI18.tif", format = "GTiff")
writeRaster(iNDVI19, "./data/GIS/iNDVI19.tif", format = "GTiff")
rm(iNDVI, iNDVI18, iNDVI19)

#----------------------#
# distance to roads ####
#----------------------#

roads <- st_read("T:/Roads_Tiger2016","roads_GYEandWY_Tiger2016_AEA")
sa <- st_transform(sa, crs=st_crs(roads))   # get study area into crs of roads
roads <- roads[st_intersects(roads, sa, sparse = FALSE)[,1],]
roads2 <- st_read("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021/data/roads_Jackson_Colorado",
                  "JacksonCo_CO_roads")   # roads in Colorado
head(roads)
head(roads2)
roads2 <- st_transform(roads2, crs=st_crs(roads))
roads <- rbind(roads, roads2)
roads <- roads[st_intersects(roads, sa, sparse = FALSE)[,1],]
mapview(roads)
head(roads)
# write out the roads as a line file. 
st_write(roads, "./data/GIS", "roads_tiger_2016", driver="ESRI Shapefile")

# the lines themselves are not that useful. Need to rasterize the roads 
landcov <- raster("./data/GIS/landcover_descrete_30m.tif") # use landcov raster for rasterizing the roads
res(landcov)
landcov <- aggregate(landcov, fact=3)   # I aggregated my raster so I have larger cells (in my case to 90m res), which speeds up processing time below
res(landcov)
roads <- st_transform(roads, crs=projection(landcov))
road_length <- rasterize(as(st_union(roads),"Spatial"), landcov, fun='length')   # this can take some time!!!! 
# my landcov is 167,890 cells. it took ~4 minutes. If your raster is larger than this, let me know!
# note that rasterizing and distance/direction to functions are slow in R. Sometimes better to do in ArcGIS.
plot(road_length)
# note that the values in the raster are in KMs! Not meters. This is important when calculating density later
# also note that it put an NA where there are NO roads.
sum(values(road_length), na.rm=T)  # this is the length of all the roads (in KM) of all the pixels

Dist2roadrast <- distance(road_length)   #  distance, for all cells that are NA, to the nearest cell that is not NA.
plot(Dist2roadrast)
mapview(Dist2roadrast)  # as you can see my roads don't quite cross the border. 
# write out the raster
writeRaster(Dist2roadrast, "./data/GIS/Dist2roads_meters_90m.tif", format = "GTiff")

#-----------------#
# road density ####
#-----------------#
# note, you can use the following code for any linear feature (e.g., rivers and streams)

# biggest thing when creating a road density is the scale at which you calculate density
# e.g., do you want road density for a square km around each pixel? Or 10 square kms?
f <- matrix(1, nrow=11, ncol=11)   # you can change the scale of your window. Here, I have a 90m res, and i want to
# get to 1km circle, so I will have approx 11 pixels in each direction to get 1km box
road_density <- focal(road_length, w=f, fun=sum, pad=TRUE, padValue=NA, na.rm=T) # we will sum up the lengths, so we'll have
# KMs of roads per square km for each pixel. # this can take some time!
plot(road_density)
table(is.na(values(road_density))) # there are NAs!
# need to replace the NAs with 0s
road_density[is.na(values(road_density))] <- 0
plot(road_density)
# write out the raster
writeRaster(road_density, "./data/GIS/RoadDensity_KMperKM2_90m.tif", format = "GTiff")
rm(Dist2roadrast, landcov, road_density, road_length, roads, f)



#----------------------------------------------#
# some plotting to look at animal movements ####
#----------------------------------------------#

# make sure we are ordered
data <- data[order(data$id, data$date),]
ids <- unique(data$id)

#grab some rasters of interest
ele30 <- raster("./data/GIS/Elevation_meters_30m.tif")
hillsh <- raster("./data/GIS/hillshade_unitless_30m.tif")

#plot each individual over top of one or more of the gis layers
for(i in 1:length(ids)){
  tmp <- st_transform(data[data$id == ids[i],], projection(ele30)) # grab teh data from the id of interest
  plot(ele30, main=paste("AID", ids[i]), ext=extent(tmp), axes=FALSE, 
       legend=FALSE, box=FALSE, col=terrain.colors(50))
  plot(hillsh, col=grey.colors(50), ext=extent(tmp), alpha=.5, box=FALSE, 
       legend=FALSE, axes=FALSE, add=TRUE)
  lines(st_coordinates(tmp), col="grey")
  points(st_coordinates(tmp), pch=".")
  # readline("Hit Enter to see the next individual:")     #NOTICE THAT YOU MUST CLICK ENTER TO SEE NEXT ID!!!!!!
}


#------------------------------#
# plots of ids with mapview ####
#------------------------------#

# make sure we are ordered
data <- data[order(data$id, data$date),]
ids <- unique(data$id)

# loop 
for(i in 1:length(ids)){
  tmp <- data[data$id == ids[i],]
  mp <- mapview(tmp, layer.name=paste0("ID: ",ids[i])) +
    mapview(st_cast(st_combine(tmp), "LINESTRING"), layer.name="lines", color="red")
  print(mp)
  # note that you could use the function mapshot to write each of these out to a html file!
  readline("Click enter to see the next one (or hit esc to get out of this loop):")   # NOTE!!!! you must click enter to see each one! click esc to get out of this loop!
}



