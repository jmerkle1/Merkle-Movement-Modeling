#--------------------------------------#
#- Step Selection Functions --------####
#--- Initial analyses -----------------# 
#---------- Jerod Merkle --------------#
#------------- Lab 5c ------------------#

library(dummies)
library(survival)
library(car)
library(raster)
library(sf)
library(tidyverse)
library(ggplot2)
library(ggeffects)
library(lme4)
library(itsadug)
library(visreg)

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

data$elev_step <- data$elev_step/1000   # get into KM
data$Dist2Escape_step <- data$Dist2Escape_step/1000   # get into KM
data$Dist2roads_step <- data$Dist2roads_step/1000   # get into KM

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

# create dummy variables for landcover variable (don't forget, need to drop one in analysis)
table(is.na(data$landcov_target))  #this ideally should be ALL FALSE. If not, then you should remove the lines with NAs
data$landcover_target <- data$landcov_target
dum <- data %>% 
  dplyr::select(landcov_target) %>% 
  dummy.data.frame(names="landcov_target") #create dummy variables

head(dum)
length(unique(data$landcov_target)) == ncol(dum)   # this should be TRUE
# now cbind up the new columns
data <- cbind(data, dum)
head(data)
rm(dum)


#--------------------------------------#
# Check Correlation among variables ####
#--------------------------------------#
#start by checking correlation
vars_target <- c("elev_target","trasp_target","slope_target","tpi_target",   # for the target point variables
               "Dist2Escape_target","Dist2roads_target","treecov_target",
               "iNDVI19_target","AnnualForbGrass_target","BareGround_target",
               "PerennialForbGrass_target")
correl <- data %>% 
  select(all_of(vars_target)) %>% 
  cor(use="pairwise.complete.obs", method="pearson") %>% 
  round(3)
ifelse(abs(correl)>.5, correl, NA)  # a cleaner way to look at it

# I cannot have slope_target and Dist2Escape_target in same model. I'll go with slope_target
# I cannot have PerennialForbGrass_target and treecov_target in same model. I'll go with PerennialForbGrass_target


vars_step <- c(  "elev_step","trasp_step","slope_step","tpi_step",   # for the step variables
                 "Dist2Escape_step","Dist2roads_step","treecov_step",
                 "iNDVI19_step","AnnualForbGrass_step","BareGround_step",
                 "PerennialForbGrass_step")
correl <- data %>% 
  select(all_of(vars_step)) %>% 
  cor(use="pairwise.complete.obs", method="pearson") %>% 
  round(3)
ifelse(abs(correl)>.5, correl, NA)  # a cleaner way to look at it

# I cannot have slope_step and Dist2Escape_step in same model. I'll go with slope_step
# I cannot have PerennialForbGrass_step and treecov_step in same model. I'll go with PerennialForbGrass_step

# update variable objects less the ones you deemed too correlated with others
vars_target <- c("elev_target","trasp_target","slope_target","tpi_target",   # for the target point variables
                 "Dist2roads_target","iNDVI19_target","AnnualForbGrass_target",
                 "BareGround_target","PerennialForbGrass_target")
vars_step <- c(  "elev_step","trasp_step","slope_step","tpi_step",   # for the step variables
                 "Dist2roads_step","iNDVI19_step","AnnualForbGrass_step",
                 "BareGround_step","PerennialForbGrass_step")

# now concatinate all the variables together so you have all the vars you'll need in one object
vars_all <- c(vars_target, vars_step)


# ----------------------------------#
# review variables distributions ####
# ----------------------------------#
data$case_fact <- ifelse(data$case == 1, "Used","Available")
for(i in 1:length(vars_all)){   #once you run this, use the back button on your plots tab to take a look at the distributions
  # Use semi-transparent fill
  print(ggplot(data, aes_string(x=vars_all[i], fill="case_fact")) +
    geom_density(alpha=0.4))
}  #push backwards on the graph window to have a look


# -----------------------------------#
# remove strata with missing data ####
# -----------------------------------#

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

#------------------------------#
# Fit models! --------------####
#------------------------------#
data <- data[order(data$id_yr_seas, data$strata, data$case),]  # this should be ordering by id, then by time

#paramterize with conditional logistic regression - specifying cluster will envoke robust SE estimates using General estimating equations
# Add dist and log(dist) as per Forrester et al. 2009 (Ecology) and Avgar et al. 2016 (MEE)
mpoint <- clogit(case~dist+log(dist)+elev_target+trasp_target+slope_target+
                   tpi_target+Dist2roads_target+iNDVI19_target+AnnualForbGrass_target+
                   BareGround_target+PerennialForbGrass_target+
                   strata(strata)+ cluster(id_yr_seas), # cluster term says that you assume steps within each cluster are correlated with each other, but steps among the clusters are independent
                 x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
                 method = "efron",data=data)
summary(mpoint)
confint(mpoint)

# test variance inflation factors (they should all be below 3-4 or so)
vif(glm(case~dist+log(dist)+elev_target+trasp_target+slope_target+
          tpi_target+Dist2roads_target+iNDVI19_target+AnnualForbGrass_target+
          BareGround_target+PerennialForbGrass_target, data=data))


mstep <- clogit(case~dist+log(dist)+elev_step+trasp_step+slope_step+
                  tpi_step+Dist2roads_step+iNDVI19_step+AnnualForbGrass_step+
                  BareGround_step+PerennialForbGrass_step+
                  strata(strata)+cluster(id_yr_seas),
                x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
                method = "efron",data=data)
summary(mstep)

# test variance inflation factors
vif(glm(case~dist+log(dist)+elev_step+trasp_step+slope_step+
          tpi_step+Dist2roads_step+iNDVI19_step+AnnualForbGrass_step+
          BareGround_step+PerennialForbGrass_step, data=data))


round(summary(mpoint)$coefficients,3)
round(summary(mstep)$coefficients,3)


#test QIC between the two models
source("./MiscFunctions/QIC.coxph.R")

cbind(model=c("point","step"),rbind(QIC.coxph(mpoint, details=TRUE),   
                                    QIC.coxph(mstep, details=TRUE)))
# QICR is what you want to minimize. It's like AIC, but for models with a cluster term specified
# Check and make sure n and nevent are the exact same between the different models you are comparing!

# --------------------------------------------#
# determining autocorrelation in residuals ####
# --------------------------------------------#
# as per Forester et al. 2009 Ecology 90:3554-3565.
# 
# Calculate the sum of the deviance residuals for each stratum (i.e., the sum of the residuals 
# for the case and all associated controls).
data$residuals <- residuals(mpoint, type = "deviance")
resid_df <- data %>% 
  group_by(id_yr_seas, strata) %>% 
  summarise(residuals = sum(residuals))
resid_df  

# Fit an intercept-only mixed-effects model of those summed residuals
rm1 <- lmer(residuals ~ 1 + (1 | id_yr_seas), data = resid_df)
# plot autocorrelation (want all points to be below dotted line which indicates significance of autocorrelation)
acf_resid(rm1)


#------------------------------------------------#
# plot predicted relationships with ggpredict ####
#------------------------------------------------#

# remember, all response variables are relative odds of selection!!! 
preds <- ggpredict(mpoint, terms = "elev_target", type="fixed")
plot(preds)

preds <- ggpredict(mstep, terms = "elev_step", type="fixed")
plot(preds)

preds <- ggpredict(mpoint, terms = c("elev_target","tpi_target"), type="fixed")
plot(preds)

preds <- ggpredict(mpoint, terms = c("Dist2roads_target","AnnualForbGrass_target"), type="fixed")
plot(preds)

