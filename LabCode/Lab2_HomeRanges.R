#------------------------#
# Home ranges and UDs ####
#------------------------#
#--- Jerod Merkle -------#
#------ Lab 2 -----------#

#packages
library(sf)
library(raster)
library(mapview)
library(tidyverse)
library(adehabitatHR)
library(move)
library(BBMM)
library(lubridate)
library(parallel)

# set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#
data <- readRDS("./data/bhs_data_summer.rds")
head(data)



#-----------------------------------#
# create minimum convex polygons ####
#-----------------------------------#

mcpS <- mcp.area(as(data,"Spatial")[,"id_yr_seas"], percent=c(50,95,99), 
                 unin="m", unout="km2", plotit=FALSE)
round(mcpS,2)   #this gives you the area of the home ranges in km^2
#this gives you the polygons in a spatialpolygonsdataframe
mcppoly <- mcp(as(data,"Spatial")[,"id_yr_seas"], percent=c(99), unin="m", unout="km2")
class(mcppoly)  # Note that this is spatial polygons object from sp package (not sf)
plot(mcppoly)   # plot all of them
head(mcppoly)

# plot each id (randomly)
whichID <- sample(unique(data$id_yr_seas),1)  # grab a random id_yr_seas
plot(mcppoly[mcppoly$id==whichID,], lwd=5)   # how to plot 1 id
lines(st_coordinates(data[data$id_yr_seas == whichID,]), col="grey")
plot(data[data$id_yr_seas == whichID,1], add=TRUE)

#mapview maps
mapview(mcppoly[mcppoly$id==whichID,], col.regions="white", layer.name=paste0("MCP for ID: ",whichID)) + 
  mapview(st_cast(st_combine(data[data$id_yr_seas == whichID,]), "LINESTRING"), layer.name="lines", color="red")+
  mapview(data[data$id_yr_seas == whichID,], layer.name="points", color="black")

rm(whichID)

#---------------------------------------------#
# kernel home ranges with reference method ####
#---------------------------------------------#

#first, create raster/grid to calculate UDs over
ext <- extent(data)
multiplyers <- c((ext[2]-ext[1])*0.3, (ext[4]-ext[3])*0.3)   # add about 30% around the edges of your extent (you can adjust this if necessary)
ext <- extend(ext, multiplyers)
grd <- raster(ext)
res(grd) <- 250     #i'm using a 250m resolution here. Might want to increase this if things are going slowly. Or decrease it if you want more precision
ncell(grd)     # my grid has 5,612 cells at 250m resolution
projection(grd) <- st_crs(data)$proj4string
rm(multiplyers, ext)

#plot your grid
plot(extent(grd))   # this is the bounding box of your grid
plot(sample_n(data[,"id"],1000), add=TRUE)   # add a sample of your points


#calculate the kernel HR
kernS <- kernelUD(as(data,"Spatial")[,"id_yr_seas"], h="href", 
                  grid=as(grd, "SpatialPixels"), kern="bivnorm")
class(kernS)   #this is the class of the adehabitatHR package

conts <- getverticeshr(kernS, percent=95)   #get the contours as a spatial polygons dataframe

# plot the kernels, one individual at a time
image(kernS[[3]])   
plot(conts[3,], add=TRUE)

mapview(conts[3,], col.regions=NA) + mapview(kernS[[3]])


#calculate contour areas from kernel UD
kernShr <- kernel.area(kernS, percent = c(50, 95, 99), unin="m", unout="km2")
round(kernShr,2)

# -------------------------------------------#
# write out the UDs to file (if you want) ####
# -------------------------------------------#
#create a folder to write your kernal UDs to
fldr <- "./UDs"
if(dir.exists(fldr)==FALSE){   # if it doesn't exist, it'll create it for you
  dir.create(fldr)
}
# check if there any othr files in it
length(dir(fldr))   #this should be 0. Otherwise you may want to remove files from this folder.

for(i in 1:length(kernS)){   # loop through each id_yr_seas and write to file
  tmp <- raster(kernS[[i]])   # turn the estUD object to a raster
  tmp[] <- values(tmp)/sum(values(tmp))  # standardize so values sum to 1
  writeRaster(tmp, paste0(fldr, "/Kernel_", names(kernS)[i],".tif"), 
              format="GTiff", overwrite=TRUE)
}
rm(tmp)

# kernel home ranges with LSCV method (IF YOU WANT) ####
# calculate UD (must think about what h value/method you will use)
# Note that this function can take some serious time
kernS_LSCV <- kernelUD(as(data,"Spatial")[,"id_yr_seas"], h="LSCV", grid=as(grd, "SpatialPixels"),
                       same4all=FALSE,kern="bivnorm", hlim=c(.3,3), extent=3)
plotLSCV(kernS_LSCV[[1]])  #plot the first one

# If it is not converging, then forget about this LSCV method!!!!!!

# calculate contours from kernel UD
# kernShr_LSCV <- kernel.area(kernS_LSCV, percent = c(50, 95, 99), unin="m", unout="km2")
# kernShr_LSCV

#  table of results 
# tblS <- data.frame(method=c("mcp","kernel_ref","kernel_LSCV"),
#                    season="summer",
#                    rbind(colMeans(t(mcpS)), colMeans(t(kernShr)), #calculate means
#                          colMeans(t(kernShr_LSCV))),
#                    rbind(apply(t(mcpS),2,sd)/sqrt(ncol(mcpS)), #calculate SE
#                          apply(t(kernShr),2,sd)/sqrt(ncol(mcpS)),
#                          apply(t(kernShr_LSCV),2,sd)/sqrt(ncol(mcpS))))
# names(tblS) <- c("method","season", "cont50_mean", "cont95_mean","cont95_mean",
#                  "cont50_se", "cont95_se","cont95_se")
# tblS


#---------------------------------------------------#
# visualize your results from MCP and kernel-REF ####
#---------------------------------------------------#
whichID <- sample(unique(data$id_yr_seas), 1)    #choose an ID to look at

temp <- data[data$id_yr_seas==whichID,c("id_yr_seas")]
#create HRs and UDs
mcp <- mcp(as(temp, "Spatial"), percent=c(95), unin="m", unout="km2")
kref <- kernelUD(as(temp, "Spatial"), h="href", grid=150, kern="bivnorm")
krefHR_95 <- getverticeshr(kref, 95)    #you could change the countour value, if you'd like
krefHR_99 <- getverticeshr(kref, 99)    #you could change the countour value, if you'd like

#now for some cheesy plotting
par(mfrow=c(1,2), mai=c(0.5,0,0.4,0), mgp=c(.2,.5,0))
image(kref[[1]], col="red",  main="MCP")
plot(mcp, add=T)
plot(as(temp, "Spatial"), add=T, col="grey", pch=".")
image(kref[[1]], main="kernel - ref")
plot(as(temp, "Spatial"), add=T, col="grey", pch=".")
plot(krefHR_95, add=T)
plot(krefHR_99, add=T)
mtext(whichID, outer = T, line=-2, side=1, cex=1.5)
rm(temp, mcp, kref, krefHR_95, krefHR_99, u, whichID)

# note that there are other kernel methods out there.
# see packages ctmm, for example, which has methods for autocorrelated kernel density. 

#-------------------------------------------#
# create Brownian Bridge movement models ####
#-------------------------------------------#

# First, must mark locations with relatively large gaps bewteen them (e.g., > 8 hrs)
data <- data[order(data$id_yr_seas, data$date),]  #order the database by id and date
dif <- data$date %>%   # calculate change in time between points (in hrs)
  as.numeric() %>% 
  diff() %>% 
  as.numeric() %>% 
  c(0)/3600
#replace with a 0 when it switches individual/season
dif[c(diff(as.numeric(as.factor(paste(data$id_yr_seas)))),0) != 0] <- 0 
hist(dif)
table(dif > 8)    #how often are your points > 8 hour apart?
#add a burst column for when there are breaks larger than 7 hours in the data (so i can miss 1 location with my 3hr location data)
MaxFixInterval <- 7     # identify your max interval as an object because we'll use it later
data$connect <- ifelse(dif > MaxFixInterval, "no", "yes")
table(data$connect)
rm(dif)

#---------------------#
# regular BB first ####
#---------------------#

#create a folder to write your BBs to
fldr <- "./UDs"
if(dir.exists(fldr)==FALSE){   # if it doesn't exist, it'll create it for you
  dir.create(fldr)
}
# check if there any other files in it
length(dir(fldr))   # Should only be Kernels in here. You may want to remove files from this folder.
dir(fldr)

# prepare for parallel processing:
# what will you loop over?
ids <- unique(data$id_yr_seas) #loop over id_yr_seas

# identify cores (use 1 less than you have)
no_cores <- detectCores()-1
# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("data", "grd", "ids", "MaxFixInterval", "fldr"))

# now for the actual loop
regBB <- do.call(rbind, clusterApplyLB(clust, 1:length(ids), function(i){
  # need library() here for the packages your calculations require for your calculations
  library(BBMM)
  library(raster)
  library(sf)
  
  temp <- data[data$id_yr_seas==ids[i],]
  # this is the function to calculate the regular BB
  bb <- try(brownian.bridge(x=st_coordinates(temp)[,1],
                            y=st_coordinates(temp)[,2],
                            time.lag=diff(as.numeric(temp$date)/60),   # calculate the time differences in the points
                            area.grid=coordinates(grd),   # use your same grid from kernel
                            max.lag=MaxFixInterval*60,   # specify max fix interval in minutes
                            location.error=20), #this is the location error of your collars
            silent=TRUE)    
  
  #write out results to file as you go, so if there is an error you don't loose all your work!
  if(class(bb)=="try-error"){
    return(data.frame(id_yr_seas=ids[i], MotionVar=NA))   #return an NA so you know there was an error
  }else{   # if the BB did NOT fail
    rast <- grd
    rast[] <- bb$probability   #this is the probability from bb, put into values of grd
    writeRaster(rast, paste0(fldr, "/BBreg_", ids[i],".tif"), 
                format="GTiff", overwrite=TRUE)
    # saveRDS(bb, file = paste0(fldr, "/regBB_", ids[i],".rds"))  # alternatively, save the BB you made as a rds file so you can reload it later
    return(data.frame(id_yr_seas=ids[i], MotionVar=bb[[1]]))    # have it return the motion variance
  }
  
}))
stopCluster(clust)   # you must stop the parallelization framework
head(regBB)

# were there any BBs that didn't work?
table(is.na(regBB$MotionVar))   #if any of these are TRUE, then YES, there are some BBs that failed

# now take a look at the motion variances for your animals
# Note that motion variances getting into the 8,000 and 10,000 are usually problematic!
hist(regBB$MotionVar)   # all of mine are < 1,000
rm(MaxFixInterval)

#-------------------------------#
# visualize your regular BBs ####
#-------------------------------#
whichID <- sample(ids, 1)    #choose an ID to look at
rast <- raster(paste0(fldr, "/BBreg_", whichID,".tif"))
sum(values(rast))    #these will sum to 1, because it is a probability
rast2 <- getVolumeUD(as(rast, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
plot(rast2, xlim=extent(data[data$id_yr_seas == whichID,])[1:2]+c(-1000,1000),
     ylim=extent(data[data$id_yr_seas == whichID,])[3:4]+c(-1000,1000))
plot(data[data$id_yr_seas == whichID,"id_yr_seas"], add=T, pch=".", col="black")


#-----------------------------------------------------#
# Create polygons of contours for regular BB model ####
#-----------------------------------------------------#
whichID <- sample(ids, 1)    #choose an ID to look at
percentile <- 0.99   #what contour percentile are you interested in?

rast <- raster(paste0(fldr, "/BBreg_", whichID,".tif"))
sum(values(rast))    #these will sum to 1, because it is a probability
rastR <- getVolumeUD(as(rast, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
polyR <- reclassify(rastR, rcl=matrix(c(0,percentile,1,percentile,1,0),2,3, byrow=T))
polyR[values(polyR)==0] <- NA
polyR <- rasterToPolygons(polyR, dissolve=TRUE) #this is the spatialPolygonsDataFrame representing the contour

#plot it
plot(rastR, xlim=extent(data[data$id_yr_seas == whichID,])[1:2]+c(-1000,1000),
     ylim=extent(data[data$id_yr_seas == whichID,])[3:4]+c(-1000,1000))
plot(polyR, add=TRUE)  
lines(st_coordinates(data[data$id_yr_seas == whichID,]), col="grey")
plot(data[data$id_yr_seas == whichID,"id_yr_seas"], add=T, pch=".", col="black")

#remove excess objects
rm(whichID, percentile, rastR, rast, rast2, polyR)  

#----------------------#
# HR size for regBB ####
#----------------------#
# run it on multiple processors

ids <- fldr %>%    # grab ids of successful reg BBs
  list.files("BBreg") %>% 
  str_replace_all("BBreg_","") %>% 
  str_replace_all(".tif","")

# identify cores (use 1 less than you have)
no_cores <- detectCores()-1
# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("fldr","ids","grd"))

# now for the actual loop
regBB_HRsizes <- do.call(cbind, clusterApplyLB(clust, 1:length(ids), function(i){
  # need library() here for the packages your calculations require for your calculations
  library(move)
  rast <- raster(paste0(fldr, "/BBreg_", ids[i],".tif"))
  rast <- getVolumeUD(as(rast, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
  hr95 <- reclassify(rast, rcl=matrix(c(0,0.95,1,0.95,1,0),2,3, byrow=T))
  hr95 <- (res(hr95)[1]*res(hr95)[2]/1000000)*table(values(hr95)==0)[1]
  hr99 <- reclassify(rast, rcl=matrix(c(0,0.99,1,0.99,1,0),2,3, byrow=T))
  hr99 <- (res(hr99)[1]*res(hr99)[2]/1000000)*table(values(hr99)==0)[1]
  hr50 <- reclassify(rast, rcl=matrix(c(0,0.5,1,0.5,1,0),2,3, byrow=T))
  hr50 <- (res(hr50)[1]*res(hr50)[2]/1000000)*table(values(hr50)==0)[1]
  return(c(hr50,hr95,hr99))
}))
stopCluster(clust)   # you must stop the parallelization framework
round(regBB_HRsizes,2)     #in KM^2

#  make table of results ####
tblS <- data.frame(method=c("mcp","kernel_ref","regBB"),
                   season="summer",
                   rbind(colMeans(t(mcpS)), colMeans(t(kernShr)), #calculate means
                         colMeans(t(regBB_HRsizes))),
                   rbind(apply(t(mcpS),2,sd)/sqrt(ncol(mcpS)), #calculate SE
                         apply(t(kernShr),2,sd)/sqrt(ncol(mcpS)),
                         apply(t(regBB_HRsizes),2,sd)/sqrt(ncol(mcpS))))
tblS
names(tblS) <- c("method","season", "cont50_mean", "cont95_mean","cont99_mean",
                 "cont50_se", "cont95_se","cont99_se")
tblS

#--------------------------------------------------#
# plot table of home range size results results ####
#--------------------------------------------------#
par(mgp=c(1,.1,0), mai=c(.5,.5,.03,.03))
plot(c(0.5,3.5),c(0,75), axes=FALSE, xlab="Contour", 
     ylab=bquote("Mean home range size (KM"^"  2"*")"), type="n")
axis(2, cex.axis=.7, tcl=-.1)
axis(1, at=1:3, labels=c("50%","95%","99%"), cex.axis=.7, tcl=-.1)
box(bty="o")
cols <- c("#0058a9","#47adb8","#fbef4d")
adj <- c(-.21,-.07,.07,.21)
for(i in 1:nrow(tblS)){
  arrows(1+adj[i], tblS$cont50_mean[i]-tblS$cont50_se[i],
         1+adj[i], tblS$cont50_mean[i]+tblS$cont50_se[i], 
         angle=90,length=.1, code=3, col=cols[i])
  points(1+adj[i], tblS$cont50_mean[i], col=cols[i])
  arrows(2+adj[i], tblS$cont95_mean[i]-tblS$cont95_se[i],
         2+adj[i], tblS$cont95_mean[i]+tblS$cont95_se[i], 
         angle=90,length=.1, code=3, col=cols[i])
  points(2+adj[i], tblS$cont95_mean[i], col=cols[i])
  arrows(3+adj[i], tblS$cont99_mean[i]-tblS$cont99_se[i],
         3+adj[i], tblS$cont99_mean[i]+tblS$cont99_se[i], 
         angle=90,length=.1, code=3, col=cols[i])
  points(3+adj[i], tblS$cont99_mean[i], col=cols[i])
}
legend("topleft", as.character(tblS$method), col=cols, pch=1, inset=.03, title="Method")

rm(adj, i, cols, ids)


#-------------------------#
# do Dynamic BB second ####
#-------------------------#

# verify your folder to write your BBs to is what you want!

fldr <- "./UDs"
dir(fldr)   # There shouldn't be any .tif files with dynBB... if so, you might want to delete them.

ids <- unique(data$id_yr_seas) #loop over ids

# identify cores (use 1 less than you have)
no_cores <- detectCores()-1
# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("data","fldr","ids","grd"))

# now for the actual loop
DynBB <- do.call(rbind, clusterApplyLB(clust, 1:length(ids), function(i){
   
  # DynBB <- do.call(rbind, lapply(1:length(ids), function(i){   # may need this line to troubleshoot inside the parrellel loop
  # need library() here for the packages your calculations require for your calculations
  library(move)
  library(sf)
  
  temp <- data[data$id_yr_seas==ids[i],]  # grab the data for the id of interest
  
  # prep data for the dBB function
  mov <- move(x=st_coordinates(temp)[,1], y=st_coordinates(temp)[,2], time=temp$date,
              animal=ids[i], proj=CRS(st_crs(temp)$proj4string))   #create mov object
  mov <- burst(mov, temp$connect[1:(nrow(temp)-1)])   #this identifies the bad points (ie > than your MaxFixInterval)
  
  # this is the function to calculate the dynamics BB
  Dbb <- try(brownian.bridge.dyn(mov,
                                 location.error=20, #this is the location error of your collars
                                 raster=grd,
                                 margin=11,    # margin and window.size are params of the dynBB. I have put the defaults
                                 window.size=31,
                                 burstType="yes"), 
             silent=TRUE)
  
  #write out results to file too, so if there is an error you don't loose all your work!
  if(class(Dbb)=="try-error"){  # if teh dBB failed
    return(data.frame(id_yr_seas=ids[i], Failed=TRUE, FirstError=Dbb[1]))   #return an NA so you know there was an error
  }else{  # if the dBB was a success
    if(length(Dbb@layers)>1){   #check to see if it was a multi-part DBB
      rast <- sum(Dbb)
    }else{
      rast <- Dbb[[1]]
    }
    writeRaster(rast, filename = paste0(fldr, "/BBdyn_", ids[i],".tif"), format="GTiff", overwrite=T)
    # saveRDS(Dbb, file=paste0(fldr, "/dynBB_", ids[i],".rds"))   #save the BB you made as a Rdata file so you can reload it later
    return(data.frame(id_yr_seas=ids[i], Failed=FALSE, FirstError="None"))    # have it return the motion variance
  }

}))
stopCluster(clust)   # you must stop the parallelization framework

head(DynBB)
# were there any BBs that didn't work?
table(DynBB$Failed)  #if any of these are TRUE, then YES

# some code to look at your dynamic BBs
whichID <- sample(ids, 1)    #choose an ID to look at
rast <- raster(paste0(fldr, "/BBdyn_", whichID,".tif"))
sum(values(rast))    #these will sum to 1, because it is a probability
rast <- getVolumeUD(as(rast, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
plot(rast, xlim=extent(data[data$id_yr_seas == whichID,])[1:2]+c(-1000,1000),
     ylim=extent(data[data$id_yr_seas == whichID,])[3:4]+c(-1000,1000))
plot(data[data$id_yr_seas == whichID,"id_yr_seas"], add=T, pch=".", col="black")

# remove excess objects
rm(rast)

#----------------------#
# HR size for dynBB ####
#----------------------#
# run it on multiple processors

ids <- list.files(fldr, "BBdyn") %>%   # grab the IDS that had successful DynBBs
  str_replace(".tif","") %>% 
  str_replace("BBdyn_","")
ids

# identify cores (use 1 less than you have)
no_cores <- detectCores()-1
# Setup cluster
clust <- makeCluster(no_cores) 
# export the objects you need for your calculations from your environment to each node's environment
clusterExport(clust, varlist=c("fldr","ids"))

# now for the actual loop
dynBB_HRsizes <- do.call(cbind, clusterApplyLB(clust, 1:length(ids), function(i){
  # need library() here for the packages your calculations require for your calculations
  library(move)
  
  rast2 <- raster(paste0(fldr, "/BBdyn_", ids[i],".tif"))
  # rast2 <- readRDS(paste0(fldr, "/BBdyn_", ids[i],".rds"))
  rast2 <- getVolumeUD(as(rast2, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
  hr95 <- reclassify(rast2, rcl=matrix(c(0,0.95,1,0.95,1,0),2,3, byrow=T))
  hr95 <- (res(hr95)[1]*res(hr95)[2]/1000000)*table(values(hr95)==0)[1]
  hr99 <- reclassify(rast2, rcl=matrix(c(0,0.99,1,0.99,1,0),2,3, byrow=T))
  hr99 <- (res(hr99)[1]*res(hr99)[2]/1000000)*table(values(hr99)==0)[1]
  hr50 <- reclassify(rast2, rcl=matrix(c(0,0.5,1,0.5,1,0),2,3, byrow=T))
  hr50 <- (res(hr50)[1]*res(hr50)[2]/1000000)*table(values(hr50)==0)[1]
  return(c(hr50,hr95,hr99))
}))
stopCluster(clust)   # you must stop the parallelization framework

round(dynBB_HRsizes,2)  #in KM^2

#---------------------------#
#  make table of results ####
#---------------------------#
tblS <- data.frame(method=c("mcp","kernel_ref","regBB","dynBB"),
                   season="summer",
                   rbind(colMeans(t(mcpS)), colMeans(t(kernShr)), #calculate means
                         colMeans(t(regBB_HRsizes)),colMeans(t(dynBB_HRsizes))),
                   rbind(apply(t(mcpS),2,sd)/sqrt(ncol(mcpS)), #calculate SE
                         apply(t(kernShr),2,sd)/sqrt(ncol(mcpS)),
                         apply(t(regBB_HRsizes),2,sd)/sqrt(ncol(mcpS)),
                         apply(t(dynBB_HRsizes),2,sd)/sqrt(ncol(mcpS))))
tblS
names(tblS) <- c("method","season", "cont50_mean", "cont95_mean","cont99_mean",
                 "cont50_se", "cont95_se","cont99_se")
tblS

#------------------------------------------#
# plot table of home range size results ####
#------------------------------------------#
par(mgp=c(1,.1,0), mai=c(.5,.5,.03,.03))
plot(c(0.5,3.5),c(0,75), axes=FALSE, xlab="Contour", ylab=bquote("Mean home range size (KM"^"  2"*")"), type="n")
axis(2, cex.axis=.7, tcl=-.1)
axis(1, at=1:3, labels=c("50%","95%","99%"), cex.axis=.7, tcl=-.1)
box(bty="o")
cols <- c("#0058a9","#47adb8","#fbef4d","#e6e7e8")
adj <- c(-.21,-.07,.07,.21)
for(i in 1:nrow(tblS)){
  arrows(1+adj[i], tblS$cont50_mean[i]-tblS$cont50_se[i],
         1+adj[i], tblS$cont50_mean[i]+tblS$cont50_se[i], 
         angle=90,length=.1, code=3, col=cols[i])
  points(1+adj[i], tblS$cont50_mean[i], col=cols[i])
  arrows(2+adj[i], tblS$cont95_mean[i]-tblS$cont95_se[i],
         2+adj[i], tblS$cont95_mean[i]+tblS$cont95_se[i], 
         angle=90,length=.1, code=3, col=cols[i])
  points(2+adj[i], tblS$cont95_mean[i], col=cols[i])
  arrows(3+adj[i], tblS$cont99_mean[i]-tblS$cont99_se[i],
         3+adj[i], tblS$cont99_mean[i]+tblS$cont99_se[i], 
         angle=90,length=.1, code=3, col=cols[i])
  points(3+adj[i], tblS$cont99_mean[i], col=cols[i])
}
legend("topleft", as.character(tblS$method), col=cols, pch=1, inset=.03, title="Method")

rm(adj, i, cols)

# some ggplot code
ggplot(data = data, mapping = aes(x = as.factor(percent_contour), y = HR_size_km2, fill = as.factor(percent_contour), color = as.factor(percent_contour))) +
  geom_violin(size = 0.1, trim = TRUE, alpha = 0.6) +    
  geom_jitter(size = 0.2, width = 0.05) +  
  scale_color_manual(values = c("darkorange3", "royalblue4", "black")) +
  scale_fill_manual(values = c("darkorange1", "royalblue3", "grey72")) +
  xlab("Contour (%)") +
  labs(y = expression ("Home range size in"~km^2)) +
  theme_classic(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  facet_grid(. ~ HR_estimator) # You may have to make sure that HR_estimator is a factor here!


#---------------------------------------------------#
# Create polygon of contour for dynamic BB model ####
#---------------------------------------------------#
whichID <- sample(ids, 1)   #choose an ID to look at
percentile <- 0.99   #what contour percentile are you interested in?

rast <- raster(paste0(fldr, "/BBdyn_", whichID,".tif"))
rastD <- getVolumeUD(as(rast, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
polyD <- reclassify(rastD, rcl=matrix(c(0,percentile,1,percentile,1,0),2,3, byrow=T))
polyD[values(polyD)==0] <- NA
polyD <- rasterToPolygons(polyD, dissolve=TRUE) #this is the spatialPolygonsDataFrame representing the contour

#plot it
plot(rastD, xlim=extent(data[data$id_yr_seas == whichID,])[1:2]+c(-1000,1000),
     ylim=extent(data[data$id_yr_seas == whichID,])[3:4]+c(-1000,1000))
plot(polyD, add=TRUE)   
lines(st_coordinates(data[data$id_yr_seas == whichID,"id_yr_seas"]), col="grey")
plot(data[data$id_yr_seas == whichID,"id_yr_seas"], add=T, pch=".", col="black")

# now using mapview
mapview(polyD, layer.name=paste0("ID: ", whichID))+
  mapview(st_cast(st_combine(data[data$id_yr_seas == whichID,]), "LINESTRING"), 
        layer.name="lines", color="red") +
  mapview(data[data$id_yr_seas == whichID,], color="grey", layer.name="points")


#---------------------------------------#
#plot the two BB models side by side ####
#---------------------------------------#
whichID <- sample(ids, 1)    #choose an ID to look at
percentile <- 0.99   #what contour percentile are you interested in?
buffr <- 3000    # how much of a buffer do you want around the points (in meters)

#reg BB preparation
rast <- raster(paste0(fldr, "/BBreg_", whichID,".tif"))   # read in each bb
sum(values(rast))    #these will sum to 1, because it is a probability
rastR <- getVolumeUD(as(rast, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
polyR <- reclassify(rastR, rcl=matrix(c(0,percentile,1,percentile,1,0),2,3, byrow=T))
polyR[values(polyR)==0] <- NA
polyR <- rasterToPolygons(polyR, dissolve=TRUE) #this is the spatialPolygonsDataFrame representing the contour

#dyn bb preparation
rast <- raster(paste0(fldr, "/BBdyn_", whichID,".tif"))
rastD <- getVolumeUD(as(rast, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
polyD <- reclassify(rastD, rcl=matrix(c(0,percentile,1,percentile,1,0),2,3, byrow=T))
polyD[values(polyD)==0] <- NA
polyD <- rasterToPolygons(polyD, dissolve=TRUE) #this is the spatialPolygonsDataFrame representing the contour

#now to plot
par(mfrow=c(1,2), mai=c(0.5,0,0.4,0), mgp=c(.2,.5,0))
plot(rastR, xlim=extent(data[data$id_yr_seas == whichID,])[1:2]+c(-buffr,buffr),
     ylim=extent(data[data$id_yr_seas == whichID,])[3:4]+c(-buffr,buffr), 
     axes=F, legend=FALSE, box=FALSE, main="Regular BB")
plot(polyR, add=TRUE)   
lines(st_coordinates(data[data$id_yr_seas == whichID,]), col=rgb(190,190,190,190, maxColorValue = 255))
points(st_coordinates(data[data$id_yr_seas == whichID,"id_yr_seas"]), pch=".", col="black")
#plot dynBB second
plot(rastD, xlim=extent(data[data$id_yr_seas == whichID,])[1:2]+c(-buffr,buffr),
     ylim=extent(data[data$id_yr_seas == whichID,])[3:4]+c(-buffr,buffr),
     axes=F, legend=FALSE, box=FALSE, main="Dynamic BB")
plot(polyD, add=TRUE)   
lines(st_coordinates(data[data$id_yr_seas == whichID,]), col=rgb(190,190,190,190, maxColorValue = 255))
points(st_coordinates(data[data$id_yr_seas == whichID,]), pch=".", col="black")
mtext(whichID, outer = T, line=-2, side=1, cex=1.5)


#---------------------------------------------------#
# plot kernel and the two BB models side by side ####
#---------------------------------------------------#

par(mfrow=c(1,3), mai=c(0.5,0.01,0.4,0.01), mgp=c(.2,.5,0))
kref <- raster(paste0(fldr, "/Kernel_", whichID,".tif"))
kref <- getVolumeUD(as(kref, Class="DBBMM"))
polyK <- reclassify(kref, rcl=matrix(c(0,percentile,1,percentile,1,0),2,3, byrow=T))
polyK[values(polyK)==0] <- NA
polyK <- rasterToPolygons(polyK, dissolve=TRUE) #this is the spatialPolygonsDataFrame representing the contour

# panel 1
plot(kref, xlim=extent(data[data$id_yr_seas == whichID,])[1:2]+c(-buffr,buffr),
     ylim=extent(data[data$id_yr_seas == whichID,])[3:4]+c(-buffr,buffr), 
     axes=F, legend=FALSE, box=FALSE, main="Kernel")
lines(st_coordinates(data[data$id_yr_seas == whichID,]), col=rgb(190,190,190,190, maxColorValue = 255))
points(st_coordinates(data[data$id_yr_seas == whichID,]), pch=".", col="black")
plot(polyK, add=T)
#panel 2
plot(rastR, xlim=extent(data[data$id_yr_seas == whichID,])[1:2]+c(-buffr,buffr),
     ylim=extent(data[data$id_yr_seas == whichID,])[3:4]+c(-buffr,buffr), 
     axes=F, legend=FALSE, box=FALSE, main="Regular BB")
plot(polyR, add=TRUE)   
lines(st_coordinates(data[data$id_yr_seas == whichID,]), col=rgb(190,190,190,190, maxColorValue = 255))
points(st_coordinates(data[data$id_yr_seas == whichID,]), pch=".", col="black")
# panel 3
plot(rastD, xlim=extent(data[data$id_yr_seas == whichID,])[1:2]+c(-buffr,buffr),
     ylim=extent(data[data$id_yr_seas == whichID,])[3:4]+c(-buffr,buffr),
     axes=F, legend=FALSE, box=FALSE, main="Dynamic BB")
plot(polyD, add=TRUE)   
lines(st_coordinates(data[data$id_yr_seas == whichID,]), col=rgb(190,190,190,190, maxColorValue = 255))
points(st_coordinates(data[data$id_yr_seas == whichID,]), pch=".", col="black")
mtext(whichID, outer = T, line=-2, side=1, cex=1.5)


#remove excess objects
rm(whichID, percentile, rastR, rast, polyR, bb, polyD, 
   rastD, buffr, temp, kref, krefHR, ids)  


#-------------------------------------------#
# investigating dynamic BB in more depth ####
#-------------------------------------------#

# select a small amount of data (perhaps 100 locations)
snip <- data[1010:1110,]   #name is snip
plot(st_coordinates(snip))   #find some data where the animal makes a directional movement
lines(st_coordinates(snip), col="darkblue")
range(snip$date)    # how much time do these 100 points include?
unique(snip$id_yr_seas)    #make sure you only have data from 1 ID
table(snip$connect)   # better not to have non-connecting points. This should all be yes

#create raster/grid to calculate BBs over ####
ext <- extent(snip)
multiplyers <- c((ext[2]-ext[1])*0.4, (ext[4]-ext[3])*0.4)   # add about 20% around the edges of your extent (you can adjust this if necessary)
ext <- extend(ext, multiplyers)
rm(multiplyers)
grd2 <- raster(ext)
res(grd2) <- 150     #i'm using a 150m resolution here. Might want to increase this if things are going slowly. Or decrease it if you want more precision
projection(grd2) <- st_crs(snip)$proj4string
#plot your grid
plot(extent(grd2))
plot(snip[,"id"], add=TRUE)

# develop some options to loop through
margs <- c(3,7)
winds <- c(15,25)
dats <- c("full","reduced")
options <- expand.grid(margs=margs,winds=winds,dats=dats)
options

# loop through these options and calculate dynamic BB for each (this returns a list)
DBB <- lapply(1:nrow(options), function(i){
  if(options$dats[i]=="full"){
    mov <- move(x=st_coordinates(snip)[,1], y=st_coordinates(snip)[,2], time=snip$date,
                animal=snip$id_yr_seas[1],proj=CRS(st_crs(snip)$proj4string))   #create mov object
  }else{
    snip2 <- snip[seq(1,nrow(snip),2),]    #reduce teh database by exactly half (every other point)
    mov <- move(x=st_coordinates(snip2)[,1], y=st_coordinates(snip2)[,2], time=snip2$date,
                animal=snip2$id_yr_seas[1],proj=CRS(st_crs(snip2)$proj4string))   #create mov object
  }
  Dbb <- brownian.bridge.dyn(mov,location.error=20, raster=grd2,margin=options$margs[i],
                             window.size=options$winds[i],burstType="yes")
  rast <- getVolumeUD(as(Dbb, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
  return(rast)
})
plot(DBB[[1]])  #check out what you did

# Now get the contours from the DBBs for a given percentile. 
percentile <- 0.99   #what contour percentile are you interested in?
DBB_poly <- lapply(1:length(DBB), function(i){
  polyD <- reclassify(DBB[[i]], rcl=matrix(c(0,percentile,1,percentile,1,0),2,3, byrow=T))
  polyD[values(polyD)==0] <- NA
  polyD <- rasterToPolygons(polyD, dissolve=TRUE) #this is the spatialPolygonsDataFrame representing the contour
  return(polyD)
})
plot(DBB[[1]])  #check out what you did
plot(DBB_poly[[1]], add=T)


# plot all the results in a multi-panel figure ####
cp <- colorRampPalette(rev(c("white", "yellow", "green","green","green","green", "blue","magenta")))
par(mfrow=c(2,4), mai=c(0.01,0.01,0.3,0.01), mgp=c(.1,.5,0))
for(i in 1:length(DBB)){
  plot(DBB[[i]], axes=F, legend=FALSE, box=FALSE, cex.main=.8, col=cp(20),
       main=paste("data = ", options$dats[i], "; Margin= ", options$margs[i], "; Window= ", options$winds[i],sep=""))
  plot(DBB_poly[[i]], add=TRUE)   
  lines(st_coordinates(snip), col=rgb(190,190,190,190, maxColorValue = 255))
  points(st_coordinates(snip), pch=".", col=rgb(0,0,0,190, maxColorValue = 255))
}

