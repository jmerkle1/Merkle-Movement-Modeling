#--------------------------------------------#
#- Hidden Markov models --------####
#--------------------------------------------#
#---------- Jerod Merkle --------------------#
#------------- Lab 7b ------------------------#


library(sf)
library(raster)
library(lubridate)
library(moveHMM)
library(momentuHMM)
library(tidyr)
library(dplyr)
library(ggplot2)

#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#----------------------------#
# Load up clean gps data  ####
#----------------------------#
data <- readRDS("./data/bhs_data_summer.rds")
head(data)

#---------------------------------------------#
# Get habitat covariates for each location ####
#---------------------------------------------#
# bring in your GIS data
# to keep the models relatively simple, pick 4-5 habitat covariates that you think are most important for your study animals
elev <- raster("./data/GIS/Elevation_meters_30m.tif")
trasp <- raster("./data/GIS/Aspect_TRASP_30m.tif")
slope <- raster("./data/GIS/Slope_degrees_30m.tif")
Cover_PerennialForbGrass <- raster("./data/GIS/RAP_2019_Cover_PerennialForbsGrasses.tif")

# NB: if using land cover variables, remember these will need to be relabeled 
# following same procedures used in Lab 4b (RSFs GIS Extraction) and Lab 5b (SSFs GIS Extraction)

# extract GIS data to points
data$elev <- raster::extract(elev, st_transform(data, crs=projection(elev)))
data$trasp <- raster::extract(trasp, st_transform(data, crs=projection(trasp)))
data$slope <- raster::extract(slope, st_transform(data, crs=projection(slope)))
data$PerennialForbGrass <- raster::extract(Cover_PerennialForbGrass, st_transform(data, crs=projection(Cover_PerennialForbGrass)))

# Set up data frame for use with momentuHMM
df <- data.frame(ID=data$id_yr_seas,
                 Time=data$date,
                 X=st_coordinates(data)[,1],
                 Y=st_coordinates(data)[,2],
                 elev=data$elev,    #adjust to create columns for your habitat covariates
                 trasp=data$trasp,
                 slope=data$slope,
                 PerennialForbGrass=data$PerennialForbGrass)

# Make sure data is sorted by individual & date/time
df <- df[order(df$ID, df$Time),]

head(df)


# this function preps the data to run a state space model
# put in the dataframe you just made, the coordinate type and the name of the coordinate columns
data.hmm <- momentuHMM::prepData(df, type="UTM", coordNames = c("X", "Y"), covNames = c("elev", "trasp", "slope", "PerennialForbGrass"))

# this next step will first plot the movement of each of your individuals in a single map (colored by ID)
# and will then plot various movement parameters for each individual animal
plot(data.hmm, compact = TRUE)  # note you must hit return below a bunch of times

summary(data.hmm) ##gives you summary information e.g., mean step length and turning angle

# if you have steps with step length zero, you have to account for that in the models
# by including a zero-mass starting parameter because the strictly positive distributions
# (i.e., gamma) are inadequate to use. You can either account for this by adding in a zero-mass
# parameter, or you can add a tiny amount of positive variation to the steps with length 0
# we are going to do that for simplicity's sake. For more information about how to deal with zero-
# inflation, you can read about it in the moveHMM vignette.

# first check to see if you have any steps with length 0
table(data.hmm$step==0) ## for my data, there are 22 steps with step length 0

# add in some variation for those 22 steps
data.hmm$step<-ifelse(data.hmm$step==0, 0+runif(1, 0, .001), data.hmm$step)

# after you add in the variation, these should all be FALSE, if they are not, something is wrong
table(data.hmm$step==0)

# before we run models with different states, we need to set some starting parameters for the models to start from
# we can use the mean and distribution of our actual data to decide on these
# looking at the distribution of steps and turning angles will also let us decide what
# distribution to use
par(mfrow=c(1,1))
hist(data.hmm$step)
hist(data.hmm$angle)

mean(data.hmm$step, na.rm=T)
sd(data.hmm$step, na.rm=T)

mean(data.hmm$angle, na.rm=T)
sd(data.hmm$angle, na.rm=T)

# to fit all of the models, we will use the fitHMM function in momentuHMM
# for each of the models you run in the fitHMM function, you need to provide a starting point
# for each state you are trying to identify.

# the parameters depend on the distribution that you use, for the BHS data, we will use gamma and vm distribution
# you can use gamma, Weibull, log-normal, or exponential for the step lengths
# and vonMises or wrapped Cauchy for the angles
# look at the vignette for more info on what parameters you should use for each distribution

# if you get errors when running the models, and all data that you prepped looks correct
# the errors likely are a result of the optimizer failing to converge
# this happens if you set your initial parameters to something that doesn't quite make sense
# i.e., if your average step length is 8 m and you set the mean to 8,000 it won't converge
# so just mess around with the starting parameters until you get it to properly converge

# okay, now we have a general idea of our movement data and know what distribution we want to use
# we fit a two state model, focusing on identifying "encamped" and "exploratory" states
# For your starting parameters, try to identify a reasonable short step length and high turning angle 
# for the "encamped" state (state 1), and a longer step length and low turning angle for the 
# "exploratory" state (state 2)

# -----------------------------------#
# TWO STATE MODEL - no covariates ####
# -----------------------------------#

# specify the number of states for the model to estimate
nbStates <- 2

# label states - we're going to focus on identifying the encamped & exploratory states only in this lab
stateNames <- c("encamped", "exploratory") 

# specify the distributions for the observation processes - step length and turning angle
dist <- list(step = "gamma", angle = "wrpcauchy")

# set the starting parameter values 
mu0 <- c(250, 1000) # step mean (two parameters: one for each state)
sigma0 <- c(400,1250) # step SD
stepPar0 <- c(mu0,sigma0) # combine the mean and SD into one object to fit the model

anglePar0 <- c(0.1,0.9) # angle concentration parameters (one for each state)

# combine step length and turning angle starting parameters into one object to fit the model 
Par0_m1 <- list(step=stepPar0,angle=anglePar0) 

# fit the basic/null model -- assumes states are independent of habitat covariates
m1 <- fitHMM(data = data.hmm, nbStates = 2, dist = dist, Par0 = Par0_m1, 
             estAngleMean = list(angle=FALSE), stateNames = stateNames)

# look at model output
# this gives you the mean and sd of step length and turning angle for 
# each of your states (encamped & exploratory)
m1

## then plot the model - this returns a plot of the distribution of step length and and a plot of the turning angle for
## the two states determined by the model; it then will give you an output for each individual animal, with the 
## states colored in orange and blue
plot(m1)

# -----------------------------------------------#
# TWO STATE MODEL - covariates on transitions ####
# -----------------------------------------------#

# Fit another model, where the state transition probabilities are 
# a function of different environmental covariates

# specify the formula for transition probabilities -- which covariates affect the transitions?
formula <- ~ slope + elev + trasp + PerennialForbGrass 

# for your initial parameters this time you can use those obtained from your model m1 (null model)
Par0_m2 <- getPar0(model=m1, formula=formula)

# fit second model
m2 <- fitHMM(data = data.hmm, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula)

# look at model output
# this gives you the maximum log-likelihood for this model
# then, this gives you the mean and sd of step length and turning angle for each of state (parameters)
# next, you can view the coeffs for the transition probabilities, 
# i.e. how the probability of switching from 1->2 or 2->1 changes with changes in habitat covariates
# in this case, increasing values for slope decreases the probability of switching from 1->2 and 
# increases the probability of switching from 2->1
m2

# how does it compare to first model?
m1

# then plot the model 
# first you'll see plots similar to the null model, m1 -- a plot of the distribution of step length 
# and a plot of the turning angle for the two states determined by the model
# then it will plot the state transition probabilities as a function of your chosen environmental covariates
# while holding the other covariates at their mean values (as indicated in the plot)
# then you'll get the output for each individual animal, with the states colored in orange and blue
plot(m2)

# -----------------------------------------------------------------------#
# TWO STATE MODEL - covariates on transitions and step/turning angles ####
# -----------------------------------------------------------------------#

# Fit another model, where the step length and turning angle in each state are also
# a function of different environmental covariates

# specify formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ state2(slope), #in this case only step length and turning angle in state 2 will be estimated
                       sd = ~ state2(slope)),  #to estimate for both states, remove the 'state2' and the parentheses  
           angle = list(concentration = ~ state2(slope)))

# initial parameters (obtained from nested model m2)
Par0_m3 <- getPar0(model=m1, formula=formula, DM=DM)

# fit model
m3 <- fitHMM(data = data.hmm, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
             beta0 = Par0_m3$beta, DM = DM, stateNames = stateNames,
             formula = formula)

# look at model output
# first shows you the maximum log-likelihood for this model
# then, this gives you the mean and sd of step length and turning angle for each of state (parameters)
# you can then examine the coeffs for step and turning angle parameters
# i.e. within each state, how step length and turning parameters change with changes in habitat covariates
# in this case, mean_2:slope tells you that, in state 2, step length decreases with increasing values for slope
# next, you can again view the coeffs for the transition probabilities, 
# i.e. how the probability of switching from 1->2 or 2->1 changes with changes in habitat covariates
# in this case, increasing values for slope decreases the probability of switching from 1->2 and 
# increases the probability of switching from 2->1
m3

# how does it compare with m2?
m2

# then plot the model 
# same output as for m2, but in addition to the plots of the distribution of step length 
# and turning angle for each state, you will also get plots of the relationship between 
# your chosen covariates and step length and/or turning angle, within each behavioral state, 
# of whichever state you specified, in this case this is state 2 (exploratory; indicated by plot title)
plot(m3)

# We have only given 3 example models here, but you may also wish to conduct further model selection 
# e.g., altering which covariates effect which components of the state space model 
# (transition probabilities, step length, or turning angle, or removing covariates that did not have 
# a strong effect on these components. You can then compare these in an AIC framework (below).

#--------------------#
# MODEL SELECTION #### 
#--------------------#
# you can now pull out the AIC for each of the models to compare and figure out which one has the most empirical support
AIC(m1, m2, m3)

#-----------------------#
# EXPLORE BEST MODEL #### 
#-----------------------#
# once you have your best model, you can run some further analyses to assess goodness of fit

# Plot model with confidence intervals, specify covariate values to include in visualization
# "covs" allows you to set the values for all other covariates to be held at when plotting 
# the transition probabilities as a function of each covariate (if not specified, this defaults to the mean)
plot(m3, plotCI = TRUE, covs = data.frame(elev=2452))

## this gives you a series of plots for each individual animal. The top plot shows you if they are in state 1 or 2,
## and the bottom two plots gives you the probability that the animal is in state 1 (middle plot) or state 2 (bottom plot)
plotStates(m3)

# estimate the state of the animals for each GPS location, using the viterbi algorithm (top plot from above)
viterbi(m3)

# estimate the probability of the animal being in each state for each location (bottom plots from above)
stateProbs(m3)

# Assess model fit
# compute pseudo-residuals for the steps and the angles
# pseudo-residuals are intermediate error terms 
# i.e. difference between the actual value and intermediate predicted value
pr <- pseudoRes(m3)

# plot pseudo-residuals for model diagnostics
plotPR(m3)

# plot the ACF of step pseudo-residuals
# this is the same autocorrelation function that we examined for SSFs 
# NB these are not grouped by individuals, so you may some clustering
# this can be useful for identifying cyclic patterns of behavior, which can then be accounted for in the HMM
# don't worry if there is still some autocorrelation, this can be dealt with 
# if you wish to pursue these analyses further
acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)

#--------------------#
# ACTIVITY BUDGET #### 
#--------------------#
# Calculate activity budget from your best model (proportion of time animals spent in each state)
# decode most likely state sequence
states <- viterbi(m3)
# derive percentage of time spent in each state for whole population
activity.budget <- table(states)/nrow(data.hmm)
activity.budget

# derive percentage of time spent in each state for each individual 
data.hmm$states <- states
total.states.ind <- data.hmm %>% group_by(ID) %>% count()
ind.states <- data.hmm %>% group_by(ID, states) %>% count() 
ind.states <- merge(total.states.ind, ind.states, by='ID')
ind.states$activity.budget <- ind.states$n.y/ind.states$n.x
ind.activity.budget <- ind.states[,c(1,3,5)]
ind.activity.budget

# Plot this in ggplot
# first you need to convert states to a factor
ind.activity.budget$states <- as.factor(ind.activity.budget$states)

ggplot(data=ind.activity.budget, aes(x=states, y=activity.budget, group=states, fill=states, color=states)) +
  geom_boxplot(alpha=0.4, show.legend=F)+
  geom_point(size=5)+
  scale_color_manual(values=c("#E69F00", "slateblue3"),name="states",
                     breaks=c("1", "2"),
                     labels=c("Encamped", "Exploratory"))+
  scale_fill_manual(values=c("#E69F00", "slateblue3"), name="states",
                    breaks=c("1", "2"),
                    labels=c("Encamped", "Exploratory"))+  
  xlab("Behavioral state") +
  ylab("Proportion of time spent in each state") +
  theme_classic() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_text(size=18), 
        axis.title=element_text(size=18),legend.text=element_text(size=18),legend.title=element_blank(), 
        axis.line.x=element_line(colour = 'black', size=0.5, linetype='solid'), 
        legend.justification=c(0,1), legend.position=c(0.8,0.99))


# Some more visualization
# plot tracks with states on habitat raster
# NB: this will only work if your data and your raster layer are in the same projection
slope.proj <- projectRaster(slope, crs=crs(data))

# crop your raster to the extent of your gps data (adding a small border in each direction)
bbox <- extent(data) + c(-1000,1000,-1000,1000)
slope.proj <- crop(slope.proj, bbox)

# link the states extracted from your best model with your gps data for plotting
data.hmm$states <- states

# plot animal tracks by state on your chosen raster
# all animals
plotSpatialCov(m3, slope.proj, states=data.hmm$states, col=c('#FFDB6D','#999999'))

# one/a few animals at a time
plotSpatialCov(m3, slope.proj, states=data.hmm$states, col=c('#FFDB6D','#999999'), animals=c("38942_2019_summer","38946_2019_summer"))

