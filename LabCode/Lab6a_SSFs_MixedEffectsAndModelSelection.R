#---------------------------------------------------#
#- Step Selection Functions ---------------------####
#---Model selection, validation, and mixed effects--#
#---------- Jerod Merkle ---------------------------#
#------------- Lab 6 -------------------------------#

library(survival)
library(MASS)
library(car)
library(MuMIn)
library(glmmTMB)
library(dummies)
library(sf)
library(tidyverse)
library(ggplot2)
library(ggeffects)

#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#--------------------------------------------#
# Load SSF data with variables extracted  ####
#--------------------------------------------#
# from lab RSF part b
load("SSFs_VariablesExtracted.RData")
head(data)

#---------------------------------------#
#--- Get ready for Analysis ---------####
#---------------------------------------#

# some simple scaling so your variables are on similar scales
hist(data$dist)
data$dist <- data$dist/1000   #get your distance between source and target point variable into km
hist(data$dist)

data$elev_target <- data$elev_target/1000   # get into KM
data$Dist2Escape_target <- data$Dist2Escape_target/1000   # get into KM
data$Dist2roads_target <- data$Dist2roads_target/1000    # get into KM 


# ---------------------------------------#
# some work on the landcover variable ####
# ---------------------------------------#
#fix lc variable, so we have actual variable names for landcover
# landcover variables are often numeric, so want to merge with their meaning
legend <- read.csv("./data/GIS/nlcd_legend.csv") #bring in landcov legend
names(legend)[names(legend)=="value"] <- "lc_target"   # rename value column to match your landcover column name in data
names(legend)[names(legend)=="class_general"] <- "landcov_target"   # rename the class column to the name you want to see in data
# note that this is ONLY for the target point
data <- merge(data,legend[,c("landcov_target","lc_target")], all.x=TRUE)   # add a new column with actual landcover values
table(data$landcov_target)
rm(legend)


# some basic selection ratios of your landcover (categorical variable)
# this will help you 1) get a feel for what is selected vs. avoided
# and 2) help which choosing a reference category (for the future)
round(table(data$landcov_target[data$case==1])/sum(data$case==1),3) # used pay attention here to very small number or 0s - if you find some, those classes should be identified as reference category
round(table(data$landcov_target[data$case==0])/sum(data$case==0),3)  # available


# Jerod's thoughts on how to choose a reference category:
# 1. Should be a variable you do not care about too much
# 2. should not be a rare landcover on the landscape. Should have at least 5 or 10% in availability
# 3. Best to choose one where the use is pretty similar to the availability (i.e., they use it close to its availability)
# 4. Sometimes its OK to merge a few categories that you don't care about into the reference (e.g,. water, rock, ice, etc.)

# Once you've made some decisions, here is some code to organize your reference category
# Doing the following is also necessary to use some of the automated model selection methods

# create a new landcover factor column, with the correct landcover as reference
data$landcov2_target <- as.character(data$landcov_target)  # note that I am changing it to a character vector
table(data$landcov2_target)
# reclassify one at a time (renaming the ones you want in the reference category to reference)
data$landcov2_target[data$landcov2_target=="developed"] <- "reference"
data$landcov2_target[data$landcov2_target=="water"] <- "reference"
data$landcov2_target[data$landcov2_target=="wetlands"] <- "reference"
data$landcov2_target[data$landcov2_target=="forest"] <- "reference"
table(data$landcov2_target)

# now specify your new lanccover variable as a factor
data$landcov2_target <- as.factor(data$landcov2_target)
levels(data$landcov2_target)   #these are the current levels. they are in alphabetical order, but you want reference to be the first one!
data$landcov2_target <- relevel(data$landcov2_target, ref=2) #the value of ref is where in the string of levels that your reference is located. Mine is 2nd in order
levels(data$landcov2_target)   #check and make sure it is right (your reference category must be the first level)


# -----------------------------------#
# remove strata with missing data ####
# -----------------------------------#

vars_all <- c("elev_target","trasp_target","slope_target","tpi_target", 
              "Dist2roads_target","iNDVI19_target","AnnualForbGrass_target",
              "numb_roads_crossed","BareGround_target","PerennialForbGrass_target",
              "Dist2Escape_target","treecov_target","landcov2_target")

#remove lines with even 1 NA. You want this, so you have same sample size for both target and step models
nrow(data)
data <- data[apply(data[,c("case",vars_all)],1,anyNA)==FALSE,]
nrow(data)
rm(vars_all)

#remove strata where there is no used case
tokeep <- tapply(data$case,data$strata,sum)
tokeep <- as.numeric(names(tokeep))[tokeep==1]
data <- data[data$strata %in% tokeep,]
rm(tokeep)

#how many controls do we have for each case (OK to be missing a few here and there, but shouldn't have too many with just a few)
tbl <- tapply(data$case,data$strata,length)
table(as.numeric(tbl))

# how to keep only the strata with > X number of case/controls
nrow(data)
min_numb_in_strata <- 6 # what the is is minimum number of case/controls you want per strata?
data <- data[data$strata %in% as.numeric(names(tbl[tbl >= min_numb_in_strata])),]
nrow(data)
table(as.numeric(tapply(data$case,data$strata,length)))   # this provides the distribution of the number of case/controls in each strata
rm(tbl, min_numb_in_strata)



#--------------------#
# model selection ####
#--------------------#

# manual forward stepwise approach... ####
source("./MiscFunctions/QIC.coxph.R")

# build a full model here with all of your variables including your new landcover2 factor column 
# OK to include ones that are correlated with each other...
# Note that you will get convergence issues because of so many correlated variables!!!!! That's OK.
# but you should not interpret the coefficients of this FULL model. it is just to have all the variables in one place

full <- clogit(case~dist+log(dist)+elev_target+trasp_target+slope_target+
                 tpi_target+Dist2roads_target+iNDVI19_target+AnnualForbGrass_target+numb_roads_crossed+
                 BareGround_target+PerennialForbGrass_target+Dist2Escape_target + treecov_target+
                 landcov2_target+   
                 strata(strata)+cluster(id_yr_seas),
               x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
               method = "efron",data=data)
summary(full) # doesn't mean too much b/c of all the coorelated variables

# The stepwise model selection starts here
# now, fit a null model, with only an intercept
null <- clogit(case~dist+log(dist)+   #the null model must include all the dist variables
                 strata(strata)+cluster(id_yr_seas),method = "efron",data=data)
summary(null)

#check QIC between full and null models. Something really wrong if full model isn't a smaller QIC.
cbind(model=c("full","null"),rbind(QIC.coxph(full, details=TRUE),QIC.coxph(null, details=TRUE)))

#try adding each variable in the full model, one by one to the null model
addterm(null, full, sorted=TRUE, test="Chisq")

#if there is a variable with support, add it to the next model
m1 <- clogit(case~dist+log(dist)+
               slope_target+
               strata(strata)+cluster(id_yr_seas),method = "efron",data=data)

#check summary, look at the sign and significance of the new variable
summary(m1)

#check to make sure QIC is lower for your new model
cbind(model=c("m1","null"),rbind(QIC.coxph(m1, details=TRUE),QIC.coxph(null, details=TRUE)))

#now check to see if any new terms provide any more empirical support
addterm(m1, full, sorted=TRUE, test="Chisq")

#add in the best supported new term
m2 <- clogit(case~dist+log(dist)+
               slope_target+treecov_target+
               strata(strata)+cluster(id_yr_seas),method = "efron",data=data)

# check vifs with new term (make sure the new variable isn't correlated with another variable). Don't want VIFs > 4... 
# If it is correlated, take it out and move on to the next most important variable in the addterm() result
vif(glm(case~dist+log(dist)+
          slope_target+treecov_target, 
        family = binomial(link = "logit"), data = data))   # note that you can do regular logistic regression here. No need to use conditional logistic regression to check vifs

# if VIFs are OK, check the parameters, have any of the signs flipped from the previous model? Are they similar? Are they all still significant?
summary(m2)

#check to make sure QIC is lower for your new model
cbind(model=c("m2","m1"),rbind(QIC.coxph(m2, details=TRUE),QIC.coxph(m1, details=TRUE)))

#add new term
addterm(m2, full, sorted=TRUE, test="Chisq")

# and so on... 


#---------------------#
# model validation ####
#---------------------#

# parameterize your best model 
best <- clogit(case~dist+log(dist)+
                 slope_target+treecov_target+tpi_target+trasp_target+landcov_target+
                 BareGround_target+elev_target+PerennialForbGrass_target+numb_roads_crossed+
                 strata(strata)+cluster(id_yr_seas), model=TRUE,
               x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
               method = "efron",data=data)
summary(best)

#load up the k-folds code
source("./MiscFunctions/kfoldSSF.R")


kfold_strata <- kfoldSSF(best, k=5, nrepet=10, jitter=FALSE, 
                  reproducible=TRUE,details=TRUE, sampleByCluster = FALSE)

kfold_clust <- kfoldSSF(best, k=5, nrepet=10, jitter=FALSE, 
                   reproducible=TRUE,details=TRUE, sampleByCluster = TRUE)

# spearman rank based on sampling from strata
# for observed steps
mean(kfold_strata$kfold[kfold_strata$type=="obs"])
sd(kfold_strata$kfold[kfold_strata$type=="obs"])
range(kfold_strata$kfold[kfold_strata$type=="obs"])
# for randomly chosen steps
mean(kfold_strata$kfold[kfold_strata$type=="rand"])
sd(kfold_strata$kfold[kfold_strata$type=="rand"])
range(kfold_strata$kfold[kfold_strata$type=="rand"])

# spearman rank based on sampling from clusters
# for observed steps
mean(kfold_clust$kfold[kfold_clust$type=="obs"])
sd(kfold_clust$kfold[kfold_clust$type=="obs"])
range(kfold_clust$kfold[kfold_clust$type=="obs"])
# for randomly chosen steps
mean(kfold_clust$kfold[kfold_clust$type=="rand"])
sd(kfold_clust$kfold[kfold_clust$type=="rand"])
range(kfold_clust$kfold[kfold_clust$type=="rand"])



# do cross validation without doing k-folds partitioning (just to understand what is going on)
data$predvals <- predict(best, type="risk")  #append predicted values to your database
head(data$predvals)
hist(data$predvals)   
u <- unique(data$strata)
#for each strata, rank predvals to determine the predicted rank of the target point
data <- do.call(rbind, lapply(1:length(u), function(i){
  temp <- data[data$strata == u[i],]
  temp <- temp[order(temp$predvals),]
  temp$rank <- 1:nrow(temp)
  return(temp)
}))
data <- data[order(data$strata, data$case),]

head(data, 10)

# this histogram shows the frequency of the rankings of used points 
# (or cases) within their respective stratam
hist(data$rank[data$case==1], 10)
cor(1:max(table(data$strata)),table(data$rank[data$case==1]),method="spearman")

#now for the controls or available points
hist(data$rank[data$case==0])
cor(1:max(table(data$strata)),table(data$rank[data$case==0]),method="spearman")

#----------------------------#
# True mixed effects SSFs ####
#----------------------------#

#---------------------------------------#
# analyze the data using glmmTMB()   ####
#---------------------------------------#
# Tricking a poisson model to give you mixed effects conditional logistic regression
# Check out Muff et al. 2020 (J Animal Ecol) for details
# set up model
TMBStruc = glmmTMB(case ~ dist+log(dist)+ 
                     slope_target+treecov_target+tpi_target+trasp_target+landcov_target+
                     BareGround_target+elev_target+PerennialForbGrass_target+Dist2Escape_target+numb_roads_crossed+
                     (1|strata) + (0 + slope_target | id_yr_seas ) + 
                     (0 + treecov_target | id_yr_seas) + 
                     (0 + tpi_target | id_yr_seas),  #add which random effects you want
                   family=poisson,
                   data=data,
                   doFit=FALSE)

# # Fix the standard deviation of the first random term, which is the (1|strata) component
# in the above model equation
TMBStruc$parameters$theta[1] = log(1e3)  # see Muff et al. for details
TMBStruc$mapArg = list(theta=factor(c(NA,1:3)))   # I have 3 random effects, thus the theta is NA,1,2,3. If there is 4 ranef, then c(NA,1:4)

# Fit the model
mTMB = glmmTMB:::fitTMB(TMBStruc)
summary(mTMB)
AIC(mTMB)

coefficients(mTMB)[[1]][[2]]  # here are the coefficients for each level of the random effect




#-------------------------#
# Functional responses ####
#-------------------------#

# I'll show you two ways to do functional responses. 

#---------------------------------------------------------
# 1) approach #3 from Holbrook et al. 2019 (Eco Apps)

# pull out coefficients for a given variable for each id_yr_seas (i.e., each level of your random effect)
# note, this variable must have been specified as a random slope in glmmTMB model)

FR <- data.frame(id_yr_seas = rownames(coefficients(mTMB)[[1]][[2]]),  # add names of random slopes
                 coeff_slope = coefficients(mTMB)[[1]][[2]][,"slope_target"])  # add teh coefficient values of the random slopes
FR
# now calculate mean availability for each random slope level and join to FR database
FR <- data %>% 
  filter(case==0) %>%    # get only the available values
  group_by(id_yr_seas) %>% 
  summarise(avail_slope=mean(slope_target)) %>% 
  left_join(FR)
FR
# plot the functional response
ggplot(FR, aes(avail_slope, coeff_slope)) +
  geom_point() +
  geom_smooth()

#------------------------------------------------------
# 2) approach #4 from Holbrook et al. 2019 (Eco Apps)

# For an SSF, there are two ways to do approach #4. 
# The first is by calculating mean availability of a given variable 
# within each strata
# Second is to calculate mean availability of a given variable
# within each id_yr_seas

# create column in database for mean slope within a strata
data <- data %>% 
  filter(case==0) %>%  # grab only available points
  group_by(strata) %>%   # group by strata
  summarize(slope_avail_strata=mean(slope_target)) %>% # calculate mean slope for each strata
  right_join(data) %>%   # join those mean values back to database
  as.data.frame()
head(data, 11)

# fit new model with functional response included
bestFR_strata <- clogit(case~dist+log(dist)+
                 slope_target+
                I(slope_target*slope_avail_strata) +  # here is the functional repsonse
                treecov_target+tpi_target+trasp_target+landcov_target+
                 BareGround_target+elev_target+PerennialForbGrass_target+Dist2Escape_target+numb_roads_crossed+
                 strata(strata)+cluster(id_yr_seas), model=TRUE,
               x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
               method = "efron",data=data)
summary(bestFR_strata)

#check QIC
cbind(model=c("best","bestFR_strata"),rbind(QIC.coxph(best, details=TRUE),QIC.coxph(bestFR_strata, details=TRUE)))


#plot the functional response
preds <- ggpredict(bestFR_strata, terms = c("slope_target","slope_avail_strata"), type="fixed")
plot(preds)



# second way: create column in database for mean slope within a id_yr_seas
data <- data %>% 
  filter(case==0) %>%  # grab only available points
  group_by(id_yr_seas) %>%   # group by id_yr_seas
  summarize(slope_avail_idyrseas=mean(slope_target)) %>% # calculate mean slope for each strata
  right_join(data) %>%   # join those mean values back to database
  as.data.frame()
head(data, 11)

# fit new model with functional response included
bestFR_idyrseas <- clogit(case~dist+log(dist)+
                          slope_target+
                          I(slope_target*slope_avail_idyrseas) +  # here is the functional repsonse
                          treecov_target+tpi_target+trasp_target+landcov_target+
                          BareGround_target+elev_target+PerennialForbGrass_target+Dist2Escape_target+numb_roads_crossed+
                          strata(strata)+cluster(id_yr_seas), model=TRUE,
                        x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
                        method = "efron",data=data)
summary(bestFR_idyrseas)

#check QIC
cbind(model=c("best","bestFR_strata","bestFR_idyrseas"),
      rbind(QIC.coxph(best, details=TRUE),
            QIC.coxph(bestFR_strata, details=TRUE),
            QIC.coxph(bestFR_idyrseas, details=TRUE)))


#plot the functional response
preds <- ggpredict(bestFR_idyrseas, terms = c("slope_target","slope_avail_idyrseas"), type="fixed")
plot(preds)
