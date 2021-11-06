#--------------------------------------#
#- Resource Utilization Functions --####
#--------------------------------------#
#---------- Jerod Merkle --------------#
#------------- Lab 3 ------------------#

# packages
library(sf)
library(raster)
library(mapview)
library(tidyverse)
library(stringr)
library(move)
library(dummies)
library(lme4)

#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#
data <- readRDS("./data/bhs_data_summer.rds")
head(data)

# ---------------#
# load up UDs ####
# ---------------#
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

#loop over files, and create a list of the UDs in raster format ####
DBBs <- lapply(1:length(fls), function(i){
  return(raster(fls[i])) # return the named object
  })
names(DBBs) <- names
rm(fls, names)
plot(DBBs[[3]])

lapply(DBBs, res) #check resolution of each UD (should be the same as minimum resolution of your GIS variables)
lapply(DBBs, ncell) #check number of cells in each UD
lapply(DBBs, function(i){sum(values(i))})  #verify that sums are all the same for each UD

#-----------------------------#
# Case 1, create database  ####
#-----------------------------#
#were going to loop over the uds in DBBs to create the response variable (ie 'ud') for the RUF

# now for the actual loop
datC1 <- do.call(rbind, lapply(1:length(DBBs), function(i){
  ud <- DBBs[[i]] # grab the UD in question
  proj <- projection(ud) # remember projection from DBB is different than what you are using
  
  xy <- as.data.frame(coordinates(ud)) #grab the coordinates
  ud <- data.frame(xy, ud=values(ud)) #make a dataframe with coordiantes and UD values

  # get rid of lower 1% to significantly reduce the number of data points (probably should keep these for a proper analysis)
  vals <- sort(ud$ud, decreasing = TRUE)
  cutoff <- cumsum(vals)
  cutoff <- vals[cutoff > .99][1]
  ud$ud[ud$ud < cutoff] <- 0
  #remove 0s   (this in reality is optional. You won't need a zero inflated model if you remove the 0s. 
  # but removing the zeros results in a lost in ability to discriminate between used and unused habitats)
  # your analysis is thus telling you the habitats the predict high use versus low use.
  ud <- ud[ud$ud != 0,]
  
  ud$id_yr_seas <- names(DBBs)[i] #add an id_yr_seas column
  # make spatially aware again
  ud <- st_as_sf(ud, coords=c("x","y"), dim="XY", crs=proj)
  return(ud)
}))

head(datC1)
hist(datC1$ud)   # the is the response variable for the RUF, i.e., the vale of the UD
rm(DBBs)

#----------------------#
# load up GIS data  ####
#----------------------#

elev <- raster("./data/GIS/Elevation_meters_30m.tif")
trasp <- raster("./data/GIS/Aspect_TRASP_30m.tif") # 0 = NNE, 1 = SSW
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
datC1$elev <- raster::extract(elev, st_transform(datC1, crs=projection(elev)))
datC1$trasp <- raster::extract(trasp, st_transform(datC1, crs=projection(trasp)))
datC1$slope <- raster::extract(slope, st_transform(datC1, crs=projection(slope)))
datC1$tpi <- raster::extract(tpi, st_transform(datC1, crs=projection(tpi)))
datC1$Dist2Escape <- raster::extract(Dist2Escape, st_transform(datC1, crs=projection(Dist2Escape)))
datC1$lc <- raster::extract(landcov, st_transform(datC1, crs=projection(landcov)))
datC1$Dist2roads <- raster::extract(Dist2roads, st_transform(datC1, crs=projection(Dist2roads)))
datC1$treecov <- raster::extract(treecov, st_transform(datC1, crs=projection(treecov)))
datC1$iNDVI18 <- raster::extract(iNDVI18, st_transform(datC1, crs=projection(iNDVI18)))
datC1$iNDVI19 <- raster::extract(iNDVI19, st_transform(datC1, crs=projection(iNDVI19)))
datC1$AnnualForbGrass <- raster::extract(Cover_AnnualForbGrass, st_transform(datC1, crs=projection(Cover_AnnualForbGrass)))
datC1$BareGround <- raster::extract(Cover_BareGround, st_transform(datC1, crs=projection(Cover_BareGround)))
datC1$PerennialForbGrass <- raster::extract(Cover_PerennialForbGrass, st_transform(datC1, crs=projection(Cover_PerennialForbGrass)))

head(datC1)

# landcover variables are often numeric, so want to merge with their meaning
table(datC1$lc)
legend <- read.csv("./data/GIS/nlcd_legend.csv") #bring in landcov legend
head(legend)

names(legend)[names(legend)=="value"] <- "lc"   # rename value column to match your landcover column name in datC1
names(legend)[names(legend)=="class_general"] <- "landcov"   # rename the class column to the name you want to see in datC1
datC1 <- merge(datC1,legend[,c("lc","landcov")], all.x=TRUE)    # add a new column with actual landcover values
head(datC1)
table(datC1$landcov)
rm(legend)

#------------------------#
# analyze case 1 data ####
#------------------------#

#start by checking correlation among the variables
# create a vector of your variable names (include only continuous here)
variables <- c("elev","trasp","slope","tpi","Dist2Escape","Dist2roads","treecov",
               "iNDVI19","AnnualForbGrass","BareGround","PerennialForbGrass")

correl <- datC1 %>% 
  st_set_geometry(NULL) %>% 
  select(all_of(variables)) %>% 
  cor(use="pairwise.complete.obs", method="pearson") %>% 
  round(3)
  
correl
ifelse(abs(correl)>.5, correl, NA)  # a cleaner way to look at it

# looks like I shouldn't have dist2escape and slope in the same model
# looks like I shouldn't have treecov and perennialForbGrasses in same model

# look at distributions of the variables (once you run this line, go back through your plots and look what you did)
for(i in 1:length(variables)){
  hist(st_set_geometry(datC1,NULL)[,variables[i]], main=variables[i])
}

# get the scale of the variables similar (this helps the model fit)
datC1$Dist2roads <- datC1$Dist2roads/1000 #Dist2roads to km
datC1$Dist2Escape <- datC1$Dist2Escape/1000 # Dist2Escape to km
datC1$elev <- datC1$elev/1000 # elev to KM

# look at distributions of the variables (once you run this line, go back through your plots and look what you did)
for(i in 1:length(variables)){
  hist(st_set_geometry(datC1,NULL)[,variables[i]], main=variables[i])
}

# create dummy variables for landcover variable (don't forget, need to drop one in analysis)
table(datC1$landcov)
table(is.na(datC1$landcov))  #this ideally should be ALL FALSE. If not, then you should remove the lines with NAs
datC1$landcover <- datC1$landcov
dum <- datC1 %>% 
  st_set_geometry(NULL) %>% 
  select(landcov) %>% 
  dummy.data.frame(names="landcov") #create dummy variables

head(dum)
length(unique(datC1$landcov)) == ncol(dum)   # this should be TRUE
# now cbind up the new columns
datC1 <- cbind(datC1, dum)
head(datC1)
rm(dum)


#-----------------------#
# parameterize model ####
#-----------------------#
hist(datC1$ud) #looks like we should log the response variable
hist(log(datC1$ud))   # now it looks a bit more normally distributed
datC1$log_ud <- log(datC1$ud)   # add the logged response variable to the database

names(datC1)

# mod1 has a random effect for animal id_yr_seas on the intercept (you should do this at a minimum)
mod1 <- lmer(log_ud ~ elev + trasp + slope + tpi + Dist2roads + iNDVI19 +
               PerennialForbGrass + AnnualForbGrass + BareGround + 
               landcovshrubland + landcovherbaceous + 
               (1|id_yr_seas), 
             data = datC1)
summary(mod1)
confint.merMod(mod1, method="Wald") # calculate confidence intervals
ranef(mod1) #get random effects


#-----------------------#
# Case 2 - data prep ####
#-----------------------#
#load up your data
data <- readRDS("./data/bhs_data_summer.rds")
head(data)

# create grid for counting points based on extent of data (or study area box)
#try 30m resolution, based on elevation grid
grd <- raster("./data/GIS/Elevation_meters_30m.tif")
grd[] <- 0 #change all values to 0

data <- st_transform(data, projection(grd))  #you may need to transform your point data into the projection of your gis data

# if necessary, crop down the grid to the specific points in question, if necessary
grd <- crop(grd, data, snap="out") 
plot(grd)
plot(sample_n(data,1000)[,"id"], add=T)

# create a response variable (i.e,. 'count') with number of locations per cell in grd, per id
ids <- unique(data$id_yr_seas) #loop over ids
datC2 <- do.call(rbind, lapply(1:length(ids), function(i){
  #rasterize the points, and tell the raster values to be the number of GPS points falling in a cell
  rs <- rasterize(data[data$id_yr_seas == ids[i],], grd, 
                  field="date", fun="count")
  xy <- as.data.frame(coordinates(rs)) #grab the coordinates
  ud <- data.frame(xy, count=values(rs)) #make a dataframe with coordiantes and count values
  #remove the 0s (in the real world, you may not want to do this. If you leave in all the zeros though, you'd want to do a zero-inflated poisson model)
  ud <- ud[is.na(ud$count) == FALSE,]
  #add an id_yr_seas column
  ud$id_yr_seas <- ids[i] 
  # make spatially aware again
  ud <- st_as_sf(ud, coords=c("x","y"), dim="XY", crs=projection(rs))
  return(ud)
}))
head(datC2)
hist(datC2$count)
rm(grd,ids)

#----------------------#
# load up GIS data  ####
#----------------------#

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
datC2$elev <- raster::extract(elev, st_transform(datC2, crs=projection(elev)))
datC2$trasp <- raster::extract(trasp, st_transform(datC2, crs=projection(trasp)))
datC2$slope <- raster::extract(slope, st_transform(datC2, crs=projection(slope)))
datC2$tpi <- raster::extract(tpi, st_transform(datC2, crs=projection(tpi)))
datC2$Dist2Escape <- raster::extract(Dist2Escape, st_transform(datC2, crs=projection(Dist2Escape)))
datC2$lc <- raster::extract(landcov, st_transform(datC2, crs=projection(landcov)))
datC2$Dist2roads <- raster::extract(Dist2roads, st_transform(datC2, crs=projection(Dist2roads)))
datC2$treecov <- raster::extract(treecov, st_transform(datC2, crs=projection(treecov)))
datC2$iNDVI18 <- raster::extract(iNDVI18, st_transform(datC2, crs=projection(iNDVI18)))
datC2$iNDVI19 <- raster::extract(iNDVI19, st_transform(datC2, crs=projection(iNDVI19)))
datC2$AnnualForbGrass <- raster::extract(Cover_AnnualForbGrass, st_transform(datC2, crs=projection(Cover_AnnualForbGrass)))
datC2$BareGround <- raster::extract(Cover_BareGround, st_transform(datC2, crs=projection(Cover_BareGround)))
datC2$PerennialForbGrass <- raster::extract(Cover_PerennialForbGrass, st_transform(datC2, crs=projection(Cover_PerennialForbGrass)))

head(datC2)

# landcover variables are often numeric, so want to merge with their meaning
table(datC2$lc)
legend <- read.csv("./data/GIS/nlcd_legend.csv") #bring in landcov legend
head(legend)

names(legend)[names(legend)=="value"] <- "lc"   # rename value column to match your landcover column name in datC2
names(legend)[names(legend)=="class_general"] <- "landcov"   # rename the class column to the name you want to see in datC2
datC2 <- merge(datC2,legend[,c("lc","landcov")], all.x=TRUE)    # add a new column with actual landcover values
head(datC2)
table(datC2$landcov)
rm(legend,elev,Dist2Escape, Dist2roads,iNDVI18,iNDVI19,landcov,
   slope,tpi,trasp,treecov,Cover_AnnualForbGrass,Cover_BareGround,
   Cover_PerennialForbGrass)


#----------------------------#
# analyze case 2 data     ####
#----------------------------#
#start by checking correlation among the variables
# create a vector of your variable names (mainly for continuous variables)
variables <- c("elev","trasp","slope","tpi","Dist2Escape","Dist2roads","treecov",
               "iNDVI19","AnnualForbGrass","BareGround","PerennialForbGrass")

correl <- datC2 %>% 
  st_set_geometry(NULL) %>% 
  select(all_of(variables)) %>% 
  cor(use="pairwise.complete.obs", method="pearson") %>% 
  round(3)

correl
ifelse(abs(correl)>.5, correl, NA)  # a cleaner way to look at it

# looks like I shouldn't have dist2escape and slope in the same model
# looks like I shouldn't have treecov and perennialForbGrasses in same model

# look at distributions of the variables (once you run this line, go back through your plots and look what you did)
for(i in 1:length(variables)){
  hist(st_set_geometry(datC2,NULL)[,variables[i]], main=variables[i])
}

# get the scale of the variables similar (this helps the model fit)
datC2$Dist2roads <- datC2$Dist2roads/1000 #Dist2roads to km
datC2$Dist2Escape <- datC2$Dist2Escape/1000 # Dist2Escape to km
datC2$elev <- datC2$elev/10 # elev to KM

# look at distributions of the variables (once you run this line, go back through your plots and look what you did)
for(i in 1:length(variables)){
  hist(st_set_geometry(datC2,NULL)[,variables[i]], main=variables[i])
}

# create dummy variables for landcover variable (don't forget, need to drop one in analysis)
table(is.na(datC2$landcov))  #this ideally should be ALL FALSE. If not, then you should remove the lines with NAs
datC2$landcover <- datC2$landcov
dum <- datC2 %>% 
  st_set_geometry(NULL) %>% 
  select(landcov) %>% 
  dummy.data.frame(names="landcov") #create dummy variables

head(dum)
length(unique(datC2$landcov)) == ncol(dum)   # this should be TRUE
# now cbind up the new columns
datC2 <- cbind(datC2, dum)
head(datC2)
rm(dum)

#-----------------------#
# parameterize model ####
#-----------------------#
hist(datC2$count)
table(datC2$count)
# use same variables as in mod1 so you can compare
mod2 <- glmer(count ~ elev + trasp + slope + tpi + Dist2roads + iNDVI19 +
                PerennialForbGrass + AnnualForbGrass + BareGround + 
                landcovshrubland + landcovherbaceous +  
                (1|id_yr_seas), family = poisson(link = "log"), data = datC2)
summary(mod2)
confint.merMod(mod2, method="Wald") # calculate confidence intervals
ranef(mod2) #get random effects

# compare the coefficient estimates of the two methods
round(summary(mod1)$coefficients,3)
round(summary(mod2)$coefficients,3)

