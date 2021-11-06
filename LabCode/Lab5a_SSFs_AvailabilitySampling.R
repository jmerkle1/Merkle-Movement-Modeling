#--------------------------------------#
#- Step Selection Functions --------####
#---Availability sampling -------------# 
#---------- Jerod Merkle --------------#
#------------- Lab 5a -----------------#

library(sf)
library(raster)
library(mapview)
library(tidyverse)
library(stringr)
library(parallel)
library(circular)
library(MASS)


#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#
data <- readRDS("./data/bhs_data_summer.rds")
head(data)


# do you have duplicates in your data
dup <- duplicated(st_set_geometry(data[,c("id_yr_seas", "date")], NULL))
table(dup)
data <- data[dup == FALSE,]   #this is how to remove duplicates. But be careful. Pay attention to what you are doing here!
rm(dup)

#----------------------------------#
# calculate movement parameters ####
#----------------------------------#
source("./MiscFunctions/CalcBurst.R")
source("./MiscFunctions/CalcMovParams.R")
data <- data[order(data$id_yr_seas, data$date),] #order database first
data$burst <- CalcBurst(data=data, id = TRUE, id_name="id_yr_seas", 
                        date_name="date", Tmax = 3600*24) # now I'm setting mine to just more than my fix rate, as I don't want to connect steps where missed fixes occured
length(unique(data$burst))
# now for movement params
data <- CalcMovParams(data=data, id_name = "id_yr_seas", 
                      date_name = "date", burst=data$burst)
rm(CalcBurst, CalcMovParams)

hist(data$dt/3600)
table(data$dt)   # check to make sure all your steps are the same or very similar length
table(data$dt/3600) # now in hrs   
# Note, because of small variations in when the collar actually takes a location,
# you may see a variety of values here. All your steps should be relatively close in dt.
# For example, if you have 3 hour data, you only want to use steps that 
# are between 2.8 and 3.2 hours apart. 


# This is where you say 'I don't want these steps connected for further development of the SSF'
# Anytime there is a FALSE in the StepFlag column, my 'step sampling' code knows NOT to include the step
#you will need to change the next two lines of code !!!!
sum(data$StepFlag)   # this is how many 'actual/usable' steps you currently have
data$StepFlag[data$dt < 3600*2.8] <- FALSE   # this says 'don't connect steps less than 2.8 hrs apart
data$StepFlag[data$dt > 3600*3.2] <- FALSE   # this says 'don't connect steps greater than 3.2 hrs apart
sum(data$StepFlag)   # this is how many 'actual/usable' steps you now have

# are your turning angles and step lengths correlated?
head(data)
tarad <- rad(data$rel.angle)
hist(tarad)
# a refresher on what sin and cosine does to radians
plot(seq(0,2*pi,length.out=100), sin(seq(0,2*pi,length.out=100)))  
plot(seq(0,2*pi,length.out=100), cos(seq(0,2*pi,length.out=100)))
tasin <- sin(tarad)   # this is basically a metric of left and right turns
tacos <- cos(tarad)   # this is a metric of forward and backward turns (positive values are forward, negative are backward)
cor(tasin,data$dist, use="pairwise.complete.obs", method="pearson")   #this is your pearson correlation
cor(tacos,data$dist, use="pairwise.complete.obs", method="pearson")   #this is your pearson correlation for step length and going forward versus backward
rm(tarad,tasin,tacos)

#--------------------------#
# sampling random steps ####
#--------------------------#

# take out of sf, and keep simply as a dataframe
proj <- st_crs(data)   # grab the projection for later
data$x <- st_coordinates(data)[,1]  # add x and y as columns
data$y <- st_coordinates(data)[,2]
data <- st_set_geometry(data, NULL)    #need to remove the geometry column
head(data)
class(data)


#------------------------------------------#
# sampe based on empirical distribution ####
#------------------------------------------#

source("./MiscFunctions/DrawRandSteps_emperical.R")


data <- DrawRandSteps_emperical(data=data, nr=5,   # nr = number of random steps 
                                simult=TRUE,        # should it take a step length and a turning angle simultaneously?
                                id_name="id_yr_seas", date_name="date", x_name="x", y_name="y",   #what are the names of your id and date columns?
                                withinID=FALSE,        # should sample from same individual or if FALSE all individuals in dataset?
                                distMax = Inf, uniform.angles=FALSE)
head(data, 20)
str(data)

hist(data$dist[data$case == 1])
hist(data$dist[data$case == 0])
hist(data$rel.angle[data$case==1])
hist(data$rel.angle[data$case==0])

rm(DrawRandSteps_emperical)

# ---------------------------------------------#
# sampling based on parametric distribution ####
# ---------------------------------------------#
# if you want want... 
source("./MiscFunctions/DrawRandSteps_parametric.R")

# if you want to use weibull or gamma distribution for step lengths, 
# need to convert StepFlag to FALSE when dist or speed is 0
table(data$dist[data$StepFlag == TRUE] == 0)   # how many of these do you have?  If more than a few, need to think about this and discuss
data$StepFlag[data$dist <= 0] <- FALSE  # remove those or add some small value to the step length
table(data$speed[data$StepFlag == TRUE] == 0) # just to check


# figure out which distribution fits your step length data best
fit_weib <- fitdistr(data$dist[data$StepFlag == TRUE],densfun="weibull")
fit_gamma <- fitdistr(data$dist[data$StepFlag == TRUE],densfun="gamma")
fit_exp <- fitdistr(data$dist[data$StepFlag == TRUE],densfun="exponential")

fit_weib$n; fit_gamma$n; fit_exp$n  # these need to all be the exact same to compare AIC
aic_tbl <- AIC(fit_weib, fit_gamma, fit_exp)  
aic_tbl[order(aic_tbl$AIC),]   # look for lowest AIC
rm(fit_exp,fit_gamma,fit_weib,aic_tbl)

data <- DrawRandSteps_parametric(data=data, nr=5,   # nr = number of random steps
                                  step_distr = "weibull", # Can choose weibull, gamma, or exponential (for step length/speed)
                                  ta_distr="vonmeses", # Can choose wrappedcauchy or vonmeses (for turning angles)
                                  speed=FALSE,        # should we use speed because your steps are each not the same dt?
                                  withinID=FALSE,        # should sample from same individual or other individuals in dataset?
                                  id_name="id_yr_seas", date_name="date", x_name="x", y_name="y",   #what are the names of your id and date columns?
                                  uniform.angles=FALSE, distMax = Inf)
head(data,12)
rm(DrawRandSteps_parametric)

hist(data$dist[data$case == 1])
hist(data$dist[data$case == 0])
hist(data$rel.angle[data$case==1])
hist(data$rel.angle[data$case==0])



# now that you've sampled your available steps, write out your environment! ####
# save your environment so we can load it back up another time ####
save.image("SSFs_AvailSampled.RData")
# load("SSFs_AvailSampled.RData")



