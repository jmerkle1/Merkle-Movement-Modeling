#--------------------------------------------#
#- Resource Selection Functions ----------####
#------------- Lab 4 ------------------------#
#-- Part d - Model Selection and validation -#
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
library(ggplot2)
library(MuMIn)
library(MASS)
library(car)
library(snow)
library(pROC)
library(visreg)


#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#--------------------------------------------#
# Load RSF data with variables extracted  ####
#--------------------------------------------#
# from lab RSF part b
load("RSFs_VariablesExtracted.RData")

# -----------------------------------#
# dealing with landcover variable ####
# -----------------------------------#

# some basic selection ratios of your landcover (categorical variable)
# this will help you 1) get a feel for what is selected vs. avoided
# and 2) help which choosing a reference category (for the future)
round(table(data_pop$landcov[data_pop$used==1])/sum(data_pop$used==1),3) # used pay attention here to very small number or 0s - if you find some, those classes should be identified as reference category
round(table(data_pop$landcov[data_pop$used==0])/sum(data_pop$used==0),3)  # available

# Jerod's thoughts on how to choose a reference category:
# 1. Should be a variable you do not care about too much
# 2. Should not be a rare landcover on the landscape. Should have at least ~ 5 or 10% in availability
# 3. Best to choose one where the use is pretty similar to the availability (i.e., they use it close to its availability). This helps with interpreting the coefficients of the other land cover types
# 4. Sometimes it's OK to merge a few categories that you don't care about into the reference (e.g,. water, rock, ice, etc.)

# Once you've made some decisions, here is some code to organize your reference category
# Doing the following is also necessary to use some of the automated model selection methods

# create a new landcover factor column, with the correct landcover as reference
data_pop$landcov2 <- as.character(data_pop$landcov)  # note that I am changing it to a character vector
table(data_pop$landcov2)
# reclassify one at a time (renaming the ones you want in the reference category to reference)
data_pop$landcov2[data_pop$landcov2=="developed"] <- "reference"
data_pop$landcov2[data_pop$landcov2=="water"] <- "reference"
data_pop$landcov2[data_pop$landcov2=="wetlands"] <- "reference"
table(data_pop$landcov2)

# now specify your new lanccover variable as a factor
data_pop$landcov2 <- as.factor(data_pop$landcov2)
levels(data_pop$landcov2)   #these are the current levels. they are in alphabetical order, but you want reference to be the first one!
data_pop$landcov2 <- relevel(data_pop$landcov2, ref=3) #the value of ref is where in the string of levels that your reference is located. Mine is 3rd in order
levels(data_pop$landcov2)   #check and make sure it is right (your reference category must be the first level)

# --------------------------------#
# Model Selection -  prep data ####
# --------------------------------#

# deal with NAs...
# How many NAs are in each column?
apply(data_pop, 2, function(x){table(is.na(x))})

#remove lines where there is an at least 1 NA in your columns
variables <- c("elev","trasp","slope","tpi", "landcov2",  
               "Dist2Escape","Dist2roads","treecov",
               "iNDVI19","AnnualForbGrass","BareGround","PerennialForbGrass")
nrow(data_pop)
data_pop <- data_pop[apply(data_pop[,variables],1,anyNA)==FALSE,]
nrow(data_pop)   # Hopefully you aren't deleting too many. 

# I suggest reducing your database so you can go through the motions of this section a bit more quickly
data_pop_red <- sample_n(data_pop, 4000)  # Randomly sample to about 4k points. Models seem to fit pretty quickly with this.

#-------------------------------------#
# Manual forward stepwise approach ####
#-------------------------------------#

# build a full model here with all of your variables including your new landcover2 factor column 
# OK to include ones that are correlated with each other...
# Note that you will get convergence issues because of so many correlated variables!!!!! That's OK.
# but you should not interpret the coefficients of this FULL model. it is just to have all the variables in one place

full_pop <- glmer(used ~ elev + trasp + slope + tpi + Dist2Escape + Dist2roads + treecov +
                     iNDVI19 + AnnualForbGrass + BareGround + PerennialForbGrass +landcov2 +
                    (1|id_yr_seas), # For this exercise, only use random intercept (things'll go faster)
                  family = binomial(link = "logit"), data = data_pop_red, na.action="na.fail")
summary(full_pop)


# The stepwise model selection starts here
# now, fit a null model, with only an intercept
null_pop <- glmer(used ~ 1 + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)

#try adding each variable in the full model, one by one to the null model
addterm(null_pop, full_pop, sorted=TRUE, test="Chisq")

#if there is a variable with support, add it to the next model
m1_pop <- glmer(used ~ landcov2 + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)
#check summary, look at the sign and significance of the new variable
summary(m1_pop)

#now check to see if any new terms provide any more empirical support
addterm(m1_pop, full_pop, test="Chisq", sorted=TRUE)
#add in the best supported new term
m2_pop <- glmer(used ~ landcov2 + slope + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)

# check vifs with new term (make sure the new variable isn't correlated with another variable). Don't want VIFs > 4... 
# If it is correlated, take it out and move on to the next most important variable in the addterm() result
vif(glm(used ~ landcov2 + slope, family = binomial(link = "logit"), data = data_pop_red))
# if VIFs are OK, check the parameters, have any of the signs flipped from the previous model? Are they similar? Are they all still significant?
summary(m2_pop)

# add a new term
addterm(m2_pop, full_pop, test="Chisq", sorted=TRUE)
m3_pop <- glmer(used ~ landcov2 + slope + PerennialForbGrass + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)
vif(glm(used ~ landcov2 + slope + PerennialForbGrass, family = binomial(link = "logit"), data = data_pop_red))
summary(m3_pop)
addterm(m3_pop, full_pop, test="Chisq", sorted=TRUE)
m4_pop <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)
vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi, family = binomial(link = "logit"), data = data_pop_red))
summary(m4_pop)
addterm(m4_pop, full_pop, test="Chisq", sorted=TRUE)
m5_pop <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + Dist2roads + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)
vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + Dist2roads, family = binomial(link = "logit"), data = data_pop_red))
summary(m5_pop)
addterm(m5_pop, full_pop, test="Chisq", sorted=TRUE)
m6_pop <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + Dist2roads + trasp +  + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)
vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + Dist2roads + trasp , family = binomial(link = "logit"), data = data_pop_red))
summary(m6_pop)
addterm(m6_pop, full_pop, test="Chisq", sorted=TRUE)
m7_pop <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)
vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev, family = binomial(link = "logit"), data = data_pop_red))
summary(m7_pop)
addterm(m7_pop, full_pop, test="Chisq", sorted=TRUE)
m8_pop <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + treecov + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)
vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + treecov, family = binomial(link = "logit"), data = data_pop_red))
# my treecov is pretty correlated with landcov2, gonna use next variable in line with addterm()
m8_pop <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)
vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround, family = binomial(link = "logit"), data = data_pop_red))
summary(m8_pop)
addterm(m8_pop, full_pop, test="Chisq", sorted=TRUE)
m9_pop <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround + iNDVI19 + (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)
vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround + iNDVI19, family = binomial(link = "logit"), data = data_pop_red))
summary(m9_pop)
addterm(m9_pop, full_pop, test="Chisq", sorted=TRUE)
m10_pop <- glmer(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround + iNDVI19 + AnnualForbGrass+ (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)
vif(glm(used ~ landcov2 + slope + PerennialForbGrass + tpi + trasp + Dist2roads + elev + BareGround + iNDVI19 + AnnualForbGrass, family = binomial(link = "logit"), data = data_pop_red))
summary(m10_pop)
addterm(m10_pop, full_pop, test="Chisq", sorted=TRUE)

# stop here, last variable does not improve AIC all that much... 
# Or, Likelihood ratio test not significant and/or other variables created too high of vifs

best_pop <- m10_pop    # this is the best model based on step-wise model selection
summary(best_pop)

# ----------------------------------------#
# Hypothesis-based AIC model selection ####
# ----------------------------------------#

# True AIC model selection is about testing relative empirical support for competing hypotheses formulated as models
# For example, I could assess empirical support for models containing only topographic variables vs. topo+habitat vs. topo+habitat+human

topo <- glmer(used ~ slope + tpi + trasp + elev + 
                (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)

topo_habitat <- glmer(used ~ slope + tpi + trasp + elev + 
                        landcov2 + PerennialForbGrass + BareGround + AnnualForbGrass +
                        (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)

topo_habitat_human <- glmer(used ~ slope + tpi + trasp + elev + 
                        landcov2 + PerennialForbGrass + BareGround + AnnualForbGrass +
                          Dist2roads +
                        (1|id_yr_seas), family = binomial(link = "logit"), data = data_pop_red)

AIC_tble <- model.sel(topo, topo_habitat, topo_habitat_human, rank="AIC")  # make AIC table for these three models
AIC_tble

# --------------------------#
# Dredge function - uses ####
# --------------------------#
?dredge
# A note on the dredge function. It is not wise to use dredge function to try all combinations of your variables.
# This is dredging and will increase the chances of a type I error, where you find some relationship but it is not there.
# The idea is simply that trying so many combinations of variables the chances increase that you 'find' a relationship due to random chance.


# That being said, the dredge function has many arguments that can help you with your model selection needs.
# let's have a look:

#start by creating a correlation matrix
variables <- c("AnnualForbGrass","BareGround","Dist2Escape","Dist2roads","elev",
               "iNDVI19","PerennialForbGrass","slope","tpi","trasp","treecov")
correl <- data_pop_red %>% 
  dplyr::select(all_of(variables)) %>% 
  cor(use="pairwise.complete.obs", method="pearson") %>% 
  round(3)
ifelse(abs(correl)>.5, correl, NA)  # a cleaner way to look at it

# I cannot have slope and Dist2Escape in same model. I'll go with slope_target
# I cannot have PerennialForbGrass_target and treecov in same model. I'll go with PerennialForbGrass_target

# build a matrix that provides information on which variables should and should not be included in the same model
smat <- abs(correl) <= .5    #lets use 0.5 as the cut-off for correlation, where variables that are correlated larger than this number CANNOT be included in same model
smat[!lower.tri(smat)] <- NA #Make everything in the matrix NA except (!lower.tri)the lower triangle
smat   # when FALSE, the two variables cannot be in same model
i <- as.vector(smat == FALSE & !is.na(smat))
sexpr <- parse(text = paste("!(", paste("(",variables[col(smat)[i]], " && ",variables[row(smat)[i]], ")",sep = "", collapse = " || "), ")"))
sexpr   # this tells which variables cannot be in same model. Check help with dredge to see all the different ways you can specify things here
rm(smat, i, correl)


# build a full model here with all of your variables including your new landcover2 factor column 
# OK to include ones that are correlated with each other...
# Note that you will get convergence issues because of so many correlated variables!!!!! That's OK.
# but you should not interpret the coefficients of this FULL model. it is just to have all the variables in one place

full_pop <- glmer(used ~ elev + trasp + slope + tpi + Dist2Escape + Dist2roads + treecov +
                    iNDVI19 + AnnualForbGrass + BareGround + PerennialForbGrass + landcov2 +
                    (1|id_yr_seas), # For this exercise, only use random intercept (things'll go faster)
                  family = binomial(link = "logit"), data = data_pop_red, na.action="na.fail")
summary(full_pop)


# dredge tests combinations of variables from the variables in a given FULL model, 
# based on your smat, and how many variables you want to keep in each model.
# Note that dredge can also help you with interaction terms and polynomial terms too (not shown below, though)
# we will run this using parrallel processing

#next three lines set up the parralell processing
cores <- detectCores()-1
clust <- makeCluster(getOption("cl.cores", cores), type = "SOCK")
invisible(clusterEvalQ(clust, library(lme4)))   # you need to export lme4 to each core
clusterExport(clust, "data_pop_red")    # you need to export your database to each core

dredge_pop <- pdredge(full_pop, beta="none", rank="AIC", cluster=clust,
                      subset=sexpr, #this is where you specify which variables can and cannot be in the same model together
                      fixed=c("slope","tpi"), # perhaps you want to specify variables that MUST be in each model?
                      m.lim=c(9,NA)) # tried to speed up the method by setting m.min to 10 (i.e., it must have 10 variables in each model). The second number is the upper level. Make this number the max variables that you have. NA means max.
stopCluster(clust)   # end parralel processing
head(dredge_pop)    #these are the top models
nrow(dredge_pop)     # this is how many models it tested
dredge_pop    #these are teh top models

top_mod <- glmer(used ~ elev + trasp + slope + tpi + Dist2Escape + Dist2roads + 
                    iNDVI19 + AnnualForbGrass + BareGround + PerennialForbGrass + landcov2 +
                    (1|id_yr_seas), # For this exercise, only use random intercept (things'll go faster)
                  family = binomial(link = "logit"), data = data_pop_red, na.action="na.fail")
summary(top_mod)

# -----------------------------------------#
# Plotting predicted responses, quickly ####
# -----------------------------------------#

# can fairly quickly plot responses by variable with visreg package
visreg(top_mod, xvar = "slope", type="conditional", scale="linear", 
       trans=exp, ylab="Relative Probability of Use")   # plot on linear scale (before transforming back to binomial world)

visreg(top_mod, xvar = "landcov2", type="conditional", scale="linear", 
       trans=exp, ylab="Relative Probability of Use")   # plot on linear scale (before transforming back to binomial world)

visreg(top_mod, xvar = "slope", by="landcov2", type="conditional", scale="linear", 
       trans=exp, ylab="Relative Probability of Use")   # plot on linear scale (before transforming back to binomial world)


#------------------------#
# calculate R squared ####
#------------------------#
#calculation of R^2 for mixed models based on Nakagawa and Schielzeth 2013, adn Johnson 2014.
r.squaredGLMM(top_mod)

# do AUC/ROC for logistic model (Note: because there is contamination in an RSF, Boyce recommends NOT doing AUC, but using non-parametric k-folds)
boxplot(predict(top_mod, type="response")~top_mod@frame$used, xlab="used(1) versus available(0)",
        ylab="predicted probability of use")
rocauc <- roc(top_mod@frame$used, predict(top_mod, type="response"), ci=TRUE)
rocauc   #if it is near 0.5, not a good model!
plot(rocauc)
# sensitivity is true-positive rate
# 1-specificity is the false positive rate


#----------------------------#
# kfolds cross validation ####
#----------------------------#

#Initial code to understand out how predicted values work

predvals_logit <- predict(top_mod, type="response")   #get the predicted values on the logistic scale
hist(predvals_logit) # so this is the predicted values on the logit scale, i.e., exp(model)/(1+exp(model))
hist(logit(predvals_logit))  #use the logit function to get linear predictions

predvals_linear <- predict(top_mod, type="link")   #another way to get linear predictions
hist(predvals_linear)

hist(exp(predvals_linear)/(1+exp(predvals_linear)))   #now we can get the predvals_logit by doing a logistic transformation of the linear predictors

hist(exp(predvals_linear))    # the exponent of the linear predictors
RSF_predictions <- exp(predvals_linear)   #this is the regular RSF formation (although we did not drop the intercept, but not a bit deal because everything is relative)
RSF_predictions <- (RSF_predictions-min(RSF_predictions))/(max(RSF_predictions)-min(RSF_predictions))   #scale bewteen 0 and 1, based on Manley 
hist(RSF_predictions)


# FIRST, do validation without doing k-folds partitioning (helps understand what is going on)
# get 10 bins of available data using quantiles (i.e., equal number of points per bin)
quant <- quantile(RSF_predictions[top_mod@frame$used==0], 
                  probs = seq(from = 0,to = 1, length.out = 11))
quant[1] <- 0    #make sure the lower quantile is 0
quant[length(quant)] <- Inf    #make sure the upper quantil is Inf
quant
#this shows that there are an equal number of random points in each bin
table(factor(findInterval(RSF_predictions[top_mod@frame$used==0],quant), levels = 1:10))

#break up predictions into the bins you specified with your random points
# so we know how many 1s are in each bin
int <- factor(findInterval(RSF_predictions[top_mod@frame$used==1],quant), levels = 1:10)   
table(int)

#plot it
par(mai=c(.4,.4,.02,.02))
plot(1:10, as.numeric(table(int)), ylab="Number of used values in each bin",
     ylim=c(1,max(as.numeric(table(int)))), xlab="RSF score bin", 
     type="b", cex.axis=0.6,tcl=-.1, mgp=c(1,0,0), cex.lab=.8)
legend("topleft", paste("Spearman Rank Correlation = ", round(cor(1:10, table(int), 
                                  method = "spearman"),3), sep=""), bty="n", cex=.7)

# could also do area adjusted counts in each bin too. See Boyce et al. 2006.


# SECOND, run the cross validation by splitting up your data into testing and training partitions
source("./MiscFunctions/KfoldRSF.R")   # this is code from Mathieu Basille (and with some edits by me). See package hab

cv_pop <- kfoldRSF(top_mod, k=5, nrepet = 5, nbins=10) #should repeat 50-100 times or so (although not necessary for the lab). This will take some time.
cv_pop
mean(cv_pop$kfold[cv_pop$type=="obs"]) # see Boyce et al. eval RSFs for info on how to interpret
sd(cv_pop$kfold[cv_pop$type=="obs"])
mean(cv_pop$kfold[cv_pop$type=="rand"])  # these values should be around 0
sd(cv_pop$kfold[cv_pop$type=="rand"])








# now you need to do all this coding in the analysis section for the other two RSF scales
#-----------------------------------------#
#    analyze home range scale data     ####
#-----------------------------------------#



#-----------------------------------------#
#    analyze buffer scale data         ####
#-----------------------------------------#
