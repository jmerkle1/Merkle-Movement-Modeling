#----------------------#
# Data organization ####
#----------------------#
#--- Jerod Merkle -----#
#------ Lab 0a --------#

#packages
library(sf)
library(mapview)
library(tidyverse)

#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#----------------------#
# Load up gps data  ####
#----------------------#
data <- read.csv("./data/Laramie_Region_Sheep_Data.csv", 
                 header=T, row.names=NULL, stringsAsFactors=FALSE)
head(data)
str(data)
table(is.na(data))   # are there any NAs? If so, figure out what columns they are in
apply(data, 2, function(x){sum(is.na(x))}) # figure out how many NAs in each column
table(data$CollarSerialNumber)  
table(data$X2D.3D)
data <- data[data$X2D.3D == 3,]  # I'm keeping only fixes with 3D (you may not need to do this!)

# turn date/time column(s) into POSIXct format
data$date <- paste(data$Date, data$Hour, data$Minute) %>% 
  strptime(format="%m/%d/%Y %H %M") %>% 
  as.POSIXct(tz ="MST")

table(is.na(data$date))   # Check for NAs. 
hist(data$date, breaks="weeks")   # if you have multiple years of data use months or years for breaks
range(data$date)

#reduce database and rename columns
data <- data[,c("CollarSerialNumber","date","Longitude","Latitude")]
head(data)
names(data) <- c("id","date","x","y")
head(data)

# reduce the data to the area of interest
# I built this bounding box simply from google earth or arcGIS
data <- data[data$x < -106.3,]   
data <- data[data$x > -106.6,]   
data <- data[data$y < 41.23,] 
data <- data[data$y > 40.96,] 

table(data$id)
# I'm keeping only the IDs with a significant number of points
ids2keep <- c(38898,38904,38905,38942, 
              38943,38945,38946,38947,38986)
data <- data[data$id %in% ids2keep,]
rm(ids2keep)

# ---------------------------------------#
# Remove data from prior to collaring ####
# ---------------------------------------#
# # read in a meta/capture database with a collar identifier, animal identifier, 
# # start date/time and end date/time
# meta <- read.csv("./data/2019 SE WY Capture Data.csv", 
#                  header=T, row.names=NULL, stringsAsFactors=FALSE)
# head(meta)
# meta <- meta[is.na(meta$Serial.Number)==FALSE,]
# meta$Start.Date <- meta$Start.Date %>% 
#   strptime(format="%m/%d/%Y") %>% 
#   as.POSIXct(tz ="MST")
# meta$End.Date <- meta$End.Date %>% 
#   strptime(format="%m/%d/%Y") %>% 
#   as.POSIXct(tz ="MST")
# # make sure all teh serial numbers in data are in meta
# table(unique(data$id) %in% meta$Serial.Number)  # this should all be TRUE
# 
# # reduce to only ids of interest
# meta <- meta[meta$Herd.Unit %in% c("BHS516 Douglas Creek","BHS516 State Line"),]
# 
# # now develop a loop over each serial number in your capture metadata to remove data 
# # from before capture and if possible after they have died
# data <- do.call(rbind, lapply(1:nrow(meta), function(i){
#   tmp <- data[data$id == meta$Serial.Number[i],]
#   tmp <- tmp[tmp$date > meta$Start.Date[i]+3600*48,]   # start data 2 days after capture
#   if(is.na(meta$End.Date[i])==FALSE){
#     tmp <- tmp[tmp$date < meta$End.Date[i],]   # end data on end date
#   }
#   return(tmp)
# }))
# rm(meta)

# I can visually see there are a few funky points at the beginning/end of each ID's data
# this is how to rather simply remove first and last few points of each AID
ids <- unique(data$id)
data <- do.call(rbind, lapply(1:length(ids), function(i){
  tmp <- data[data$id == ids[i],]  # grab an ID
  tmp <- tmp[order(tmp$date),]  # order by date
  tmp <- tmp[15:nrow(tmp),]     # remove first 15 points
  tmp <- tmp[1:(nrow(tmp)-15),]   # remove last 15 points
  return(tmp)
}))
rm(ids)

# make spatially aware
# what is the projection that your data are in?  You'll need to figure this out and specify below
proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# turn into sf object
data <- st_as_sf(data, coords=c("x","y"), dim="XY", crs=proj)
head(data)
rm(proj)
mapview(data) #careful, with a big dataset this could take time!

#are there outliers in your x/y values?
hist(st_coordinates(data)[,1])   
data <- data[st_coordinates(data)[,1] < -106,]   # how to remove outliers (you'll have to change)
data <- data[st_coordinates(data)[,1] > -106.6,]   # how to remove outliers (you'll have to change)
hist(st_coordinates(data)[,1])   
hist(st_coordinates(data)[,2])

#take a look more closely at dates and times
month <- as.numeric(strftime(data$date, format = "%m", tz = "MST"))
year <- as.numeric(strftime(data$date, format = "%Y", tz = "MST"))
table(year)
table(year, month)
table(month, year, data$id)
rm(year, month)

u <- unique(data$id)
for(i in 1:length(u)){
  hist(data$date[data$id == u[i]], "weeks", main=u[i])
  # readline("Hit enter to see next ID:")  # you can hit escape to get out of this loop!  
}
rm(u, i)

# do you have duplicates in your data
dup <- duplicated(st_set_geometry(data[,c("id", "date")], NULL))
table(dup)
data <- data[dup == FALSE,]   #this is how to remove duplicates. But be careful. Pay attention to what you are doing here!
rm(dup)

# transform easting and northing into a different projection (if necessary) ####
proj_new <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
st_crs(data)
data <- st_transform(data, proj_new)   #this is the function to reproject an sf object
st_crs(data)
rm(proj_new)

#-----------------------------------------------#
# how to load a shapefile into R (if needbe) ####
#-----------------------------------------------#

dat <- st_read("./data", "bhs_data")
head(dat)
plot(st_geometry(dat))   #careful, with a big dataset this could take time!
plot(dat[,"id"])    # can automatically color-code by a column with sf object
dat$date <- dat$date %>% 
  strptime(format="%Y-%m-%d %H:%M:%S") %>% 
  as.POSIXct(tz ="MST")
table(is.na(dat$date))
str(dat)
extent(dat)
head(st_coordinates(dat))
st_crs(dat)   #this is the projection
rm(dat)   # removing this dat object, and going to simply use the data object I was working on before

# -----------------------#
# figure out fix rate ####
# -----------------------#

data <- data[order(data$id, data$date),]
#first look at hrs and minutes
hr <- as.numeric(strftime(data$date, format = "%H", tz = "MST"))
hist(hr)
table(hr)
min <- as.numeric(strftime(data$date, format = "%M", tz = "MST"))
hist(min)
# now fix rate
dt <- c(diff(as.numeric(data$date)/3600),NA)   # a quick calculation of the date differences
hist(dt[dt < 25 & dt > -1]) #3 hrs is my fix rate
hist(dt[dt < 12 & dt > -1]) 
hist(dt[dt < 5 & dt > -1]) 
rm(dt, hr, min)

#----------------------------------#
# calculate movement parameters ####
#----------------------------------#
source("./MiscFunctions/CalcBurst.R")
data <- data[order(data$id, data$date),] #order database first
data$burst <- CalcBurst(data=data, id = TRUE, id_name="id", 
                        date_name="date", Tmax = 3600*7) #set Tmax to 7 hours, a little more than double the fix rate
length(unique(data$burst))
source("./MiscFunctions/CalcMovParams.R")
data <- CalcMovParams(data=data, id_name = "id", 
                      date_name = "date", burst=data$burst)


head(data)
hist(data$abs.angle) # if there is any pattern here, should be because of the shape of your study area
hist(data$rel.angle) #should be able to see forward motion here if GPS fix rate is small enough
hist(data$dist/1000) #step length distribution in km
hist(data$dt/3600) #this shows your fix rate really well (in hrs)
hist(data$speed) # in meters/sec
hist(data$speed*3600/1000) # in km/hour
sum(data$StepFlag)   # this is the number of steps you have with all the movement params.

# ----------------------------#
# Check for problem points ####
# ----------------------------#

# check for mortality problems
source("./MiscFunctions/mort.check.R")
morts <- mort.check(data=data, dist_thresh = 20,  time_thresh = 24, 
                    id_name="id", date_name="date") #dist in m, time in hrs
head(morts)

#how to remove the mortality data from the database (if you have any)
# Note, you need to be careful here. Think about what you are doing before doing it!
for(i in 1:nrow(morts)){
  toremove <- data$id == morts$id[i] & data$date >= morts$date_start[i] & data$date <= morts$date_end[i]
  table(toremove)
  data <- data[toremove == FALSE,]
}

# any time you remove points from your database, you MUST delete your old burst and move params
data$burst <- NULL
data$dist <- NULL
data$dt <- NULL
data$speed <- NULL
data$abs.angle <- NULL
data$rel.angle <- NULL

rm(morts, mort.check, toremove, i)
head(data)
str(data)

# how to find points where the speed is really crazy
source("./MiscFunctions/find.problem.pts.R")
data <- data[order(data$id, data$date),] #order database first
# every time you remove some gps points (as i did because of some morts), you must recalculate burst adn move params!
data$burst <- CalcBurst(data=data, id = TRUE,id_name="id", 
                        date_name="date", Tmax = 3600*7) #set Tmax to 7 hours, a little more than double the fix rate
length(unique(data$burst))
problems <- find.problem.pts(dat=data, date_name="date", id_name="id",
                             burst=data$burst, speedlim=2) # speedlim is meters/second. to get to km/hr do: *3600/1000
table(problems)   # this is simply a indicator of whether there is a fast speed in your data or not
data <- data[problems == FALSE,]   #if you want, remove those points 

# clean out unecessary columns
data$burst <- NULL
rm(CalcBurst, find.problem.pts, problems,CalcMovParams)


# write out the data to file
head(data)

# -------------------------#
# save cleaned database ####
# -------------------------#

# write out as rds object (more efficient). I suggest you do this, and then 
# you'll have a clean database ready for each lab in the class
saveRDS(data, file="./data/bhs_data.rds")


#write out as shapefile, if you want 
to_write <- data
to_write$date <- as.character(to_write$date)  # when writing shapefile with st_write(), need to switch posix to character.
st_write(to_write, "./data","bhs_data", driver="ESRI Shapefile")
rm(to_write)

# plot with mapview
mapview(data) #careful, with a big dataset this could take time!
# mapshot(mapview(data), url="./data/bhs_data.html")
