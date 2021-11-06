#---------------------------------------------------#
#- Step Selection Functions ---------------------####
#--- Plotting Relative Selection Strength ----------#
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

# -----------------------------------#
# remove strata with missing data ####
# -----------------------------------#

vars_all <- c("elev_target","trasp_target","slope_target","tpi_target", 
              "Dist2roads_target","iNDVI19_target","AnnualForbGrass_target",
              "numb_roads_crossed","BareGround_target","PerennialForbGrass_target",
              "Dist2Escape_target","treecov_target")

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

# ---------------------------------#
# analyze the data using clogit ####
# ---------------------------------#


# parameterize your best model 
mCLOGIT <- clogit(case~dist+
                 slope_target+treecov_target+
                 strata(strata)+cluster(id_yr_seas), model=TRUE,
               x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
               method = "efron",data=data)
summary(mCLOGIT)

#---------------------------------------#
# analyze the data using glmmTMB()   ####
#---------------------------------------#

# Tricking a poisson model to give you mixed effects conditional logistic regression
# Check out Muff et al. 2020 (J Animal Ecol) for details
# set up model
TMBStruc = glmmTMB(case ~ dist+
                     slope_target+treecov_target+
                     (1|strata) + (0 + slope_target | id_yr_seas ) + 
                     (0 + treecov_target | id_yr_seas),  #add which random effects you want
                   family=poisson,
                   data=data,
                   doFit=FALSE)

# # Fix the standard deviation of the first random term, which is the (1|strata) component
# in the above model equation
TMBStruc$parameters$theta[1] = log(1e3)  # see Muff et al. for details
TMBStruc$mapArg = list(theta=factor(c(NA,1:2)))   # I have 3 random effects, thus the theta is NA,1,2,3. If there is 4 ranef, then c(NA,1:4)

# Fit the model
mTMB = glmmTMB:::fitTMB(TMBStruc)
summary(mTMB)


# ----------------------------------------#
# Plotting Relative Selection Strength ####
# ----------------------------------------#

# this is based on methods of Avgar et al. 2017 (Ecology and Evolution)
# basically it provides: 
# What is the average change in the space use probability as we change the covariate of
# interest, while holding other covariates constant (e.g., at their mean)?

# Some info on calculating log-RSS
# https://bsmity13.github.io/log_rss/


# do clogit models first:
#------------------------#

# code below is adapted from amt tools package (log_rss function)...


# you need the following column names in x1 and x2:
colnames(mCLOGIT$x)

# identify enviro attributes of at least one proposed location
x1 <- expand.grid(dist=mean(mCLOGIT$x[,"dist"]), 
                  slope_target=seq(min(mCLOGIT$x[,"slope_target"]),max(mCLOGIT$x[,"slope_target"]),length.out=100),
                  treecov_target=mean(mCLOGIT$x[,"treecov_target"]))

# identify enviro attributes of at least one reference location
x2 <- data.frame(dist=mean(mCLOGIT$x[,"dist"]), 
                 slope_target=25,
                 treecov_target=mean(mCLOGIT$x[,"treecov_target"]))

# you are basically asking, what is the predicted probability of going to x2
# versus all the locations in x1
x2
x1

#Calculate correction due to sample-centered means (see ?survival::predict.coxph for details)
uncenter <- sum(coef(mCLOGIT) * mCLOGIT$means, na.rm=TRUE)
#predict(..., reference = "sample", se.fit = TRUE) will throw error without a "strata"
#column, even though it isn't using it
x1_dummy <- x1
x2_dummy <- x2
x1_dummy$strata = mCLOGIT$model$`strata(strata)`[1]
x2_dummy$strata = mCLOGIT$model$`strata(strata)`[1]
#Calculate predictions y_x
pred_x1 <- predict(mCLOGIT, newdata = x1_dummy, type = "lp", reference = "sample",
                   se.fit = TRUE)
pred_x2 <- predict(mCLOGIT, newdata = x2_dummy, type = "lp", reference = "sample",
                   se.fit = TRUE)
y_x1 <- pred_x1$fit + uncenter
y_x2 <- pred_x2$fit + uncenter

#Include values of x1 in return data.frame
df <- x1
#Calculate log_rss
df$log_rss <- unname(y_x1 - y_x2)
head(df)

# now calculate 95% Confidence Intervals
#Get model matrix for x1 and x2
x1_mm <- stats::model.matrix(mCLOGIT, data = x1_dummy,
                             contrast.arg = mCLOGIT$contrasts)
x2_mm <- stats::model.matrix(mCLOGIT, data = x2_dummy,
                             contrast.arg = mCLOGIT$contrasts)
#Get model variance-covariance matrix
m_vcov <- stats::vcov(mCLOGIT)
#Subtract x2 model matrix from each row of x1
delta_mm <- sweep(data.matrix(x1_mm), 2, data.matrix(x2_mm))
#Get variance of log-RSS prediction
var_pred <- diag(delta_mm %*% m_vcov %*% t(delta_mm))
#Get standard error of prediction
logrss_se <- unname(sqrt(var_pred))
#Get critical value
p <- 1 - ((1 - 0.95)/2)
zstar <- qnorm(p)
#Compute bounds
df$lwr <- df$log_rss - zstar * logrss_se
df$upr <- df$log_rss + zstar * logrss_se

head(df)


# now for plotting, wit CIs in grey

ggplot(df, aes(slope_target, log_rss, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.3) + 
  geom_line() +
  ylab("log-RSS")

# in my case, it is nearly 2 times more probable for the animal to go to a slope of 50 vs. a slope of 25. 
# Or 1.5 times less likely to go to a slope of 0 versus a slope of 25






# do glmmTMB models second:
#--------------------------#

# This is code adapted from Brian Smith, one of the co-authors of amt tools.
# based off of the code found here for predictions in GlmmTMB:
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html

# you need the following column names in x1 and x2:
rownames(summary(mTMB)$coefficients$cond)   # less the intercept intercept

# identify enviro attributes of at least one proposed location
x1 <- expand.grid(dist=mean(mTMB$frame$dist), 
                  slope_target=seq(min(mTMB$frame$slope_target),max(mTMB$frame$slope_target),length.out=100),
                  treecov_target=mean(mTMB$frame$treecov_target))

# identify enviro attributes of at least one reference location
x2 <- data.frame(dist=mean(mTMB$frame$dist), 
                 slope_target=25,
                 treecov_target=mean(mTMB$frame$treecov_target))

# We can calculate the log-RSS as y(x1) - y(x2). We can do this using
# predict.glmmTMB() to get the estimated log-RSS, but we won't be able to
# get a SE for y1 - y2 with that method.
?predict.glmmTMB

# Instead, we need to do the matrix math ourselves
# model.matrix will make the matrix we need from our data
mm1 <- model.matrix(~ dist + slope_target + treecov_target, # manually type in fixed effects only
                    data = x1)
mm2 <- model.matrix(~ dist + slope_target + treecov_target, 
                    data = x2)

# We want the prediction for the difference, so "sweep" mm2 out of mm1
mm_diff <- sweep(mm1, 2, mm2)

# Grab coefficients for fixed effects
betas <- fixef(mTMB)$cond

# Estimate log-RSS
logRSS <- mm_diff %*% betas

# Grab variance-covariance matrix to estimate SE
mm_vcov <- vcov(mTMB)$cond

# Estimate standard error
logRSS_se <- sqrt(diag(mm_diff %*% mm_vcov %*% t(mm_diff)))

# Combine everything in a data.frame
df_TMB <- x1
df_TMB$log_rss <- logRSS
df_TMB$lwr <- logRSS - 1.96 * logRSS_se
df_TMB$upr <- logRSS + 1.96 * logRSS_se


# now for plotting, wit CIs in grey
ggplot(df_TMB, aes(slope_target, log_rss, ymin = lwr, ymax = upr)) +
  geom_ribbon(alpha = 0.3) + 
  geom_line() +
  ylab("log-RSS")
