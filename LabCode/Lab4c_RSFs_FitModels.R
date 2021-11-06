#--------------------------------------------#
#- Resource Selection Functions ----------####
#------------- Lab 4 ------------------------#
#------- Part c - Initial analyses ----------#
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


#set your working drive
setwd("I:/Shared drives/wyo-coop-merkle/Courses/Mov_modeling2021")

#--------------------------------------------#
# Load RSF data with variables extracted  ####
#--------------------------------------------#
# from lab RSF part b
load("RSFs_VariablesExtracted.RData")

#-----------------------------------------#
# analyze population range scale data ####
#-----------------------------------------#

# get the scale of the variables similar (grab same code from lab 3)
data_pop$Dist2roads <- data_pop$Dist2roads/1000 #Dist2roads to km
data_pop$Dist2Escape <- data_pop$Dist2Escape/1000 # Dist2Escape to km
data_pop$elev <- data_pop$elev/1000 # elev to KM

#start by checking correlation
variables <- c("elev","trasp","slope","tpi","Dist2Escape","Dist2roads","treecov",
               "iNDVI19","AnnualForbGrass","BareGround","PerennialForbGrass")
correl <- data_pop %>% 
  dplyr::select(all_of(variables)) %>% 
  cor(use="pairwise.complete.obs", method="pearson") %>% 
  round(3)

correl
ifelse(abs(correl)>.5, correl, NA)  # a cleaner way to look at it

# note that Dist2Escape and slope should not be in the same model
# also, treecov and PerennialForbGrass should also not be in same model

# create dummy variables for landcover variable (don't forget, need to drop one in analysis)
table(is.na(data_pop$landcov))  #this ideally should be ALL FALSE. If not, then you should remove the lines with NAs
data_pop$landcover <- data_pop$landcov
dum <- data_pop %>% 
  dplyr::select(landcov) %>% 
  dummy.data.frame(names="landcov") #create dummy variables

head(dum)
length(unique(data_pop$landcov)) == ncol(dum)   # this should be TRUE
# now cbind up the new columns
data_pop <- cbind(data_pop, dum)
head(data_pop)
rm(dum)


# ----------------------------------#
# review variables distributions ####
# ----------------------------------#
data_pop$used_fact <- ifelse(data_pop$used == 1, "Used","Available")
for(i in 1:length(variables)){   #once you run this, use the back button on your plots tab to take a look at the distributions
  # Use semi-transparent fill
  print(ggplot(data_pop, aes_string(x=variables[i], fill="used_fact")) +
          geom_density(alpha=0.4))
}  #push backwards on the graph window to have a look


# if you want, scale and center continuous variables (see Schielzeth 2010, Meth Ecol Evol)
# for(i in 1:length(variables)){
#   data_pop[,variables[i]] <- scale(data_pop[,variables[i]])
#   print(variables[i])
#   print(attributes(scale(data_pop[,variables[i]])))
# }
# #take a look at what you just did:
# for(i in 1:length(variables)){   #push the back button in your plotting window to look at the variable distributions
#   hist(data_pop[,variables[i]], main=variables[i])
# }

# How many NAs are in each column?
apply(data_pop, 2, function(x){table(is.na(x))})

#remove lines where there is an at least 1 NA in your columns
nrow(data_pop)
data_pop <- data_pop[apply(data_pop[,variables],1,anyNA)==FALSE,]
nrow(data_pop)   # Hopefully you aren't deleting too many. 

#---------------------------------------------#
# population range scale model development ####
#---------------------------------------------#

# build a full model here with all of your variables...
full_pop <- glmer(used ~ elev + trasp + slope + tpi + Dist2roads + treecov + 
                    iNDVI19 + AnnualForbGrass + BareGround + PerennialForbGrass + 
                    landcovforest + (1|id_yr_seas), #think about your random effects!
                  family = binomial(link = "logit"), data = data_pop, na.action="na.fail")


summary(full_pop)
confint(full_pop, method="Wald")

# There is another method for fitting RSFs in Muff et al. 2020 (Methods in Ecology and Evolution). They glmmTMB package to fit. 
full_pop_RandSlopes <- glmmTMB(used ~ elev + trasp + slope + tpi + Dist2roads + treecov + 
                                 iNDVI19 + AnnualForbGrass + BareGround + PerennialForbGrass + 
                                 landcovforest + (1|id_yr_seas) + (0+elev|id_yr_seas) + (0+trasp|id_yr_seas)
                               + (0+slope|id_yr_seas)+ (0+tpi|id_yr_seas)+ (0+Dist2roads|id_yr_seas)
                               + (0+treecov|id_yr_seas)+ (0+iNDVI19|id_yr_seas)+ (0+AnnualForbGrass|id_yr_seas)
                               + (0+BareGround|id_yr_seas)+ (0+PerennialForbGrass|id_yr_seas)+ (0+landcovforest|id_yr_seas), #think about your random effects!
                               family = binomial(), data = data_pop, na.action="na.fail")
summary(full_pop_RandSlopes)
confint(full_pop_RandSlopes)
coefficients(full_pop_RandSlopes)  # have a look at the coefficients for each id_yr_seas

# now you need to do all this coding in the analysis section for the other two RSF scales

#----------------------------------#
# analyze Home Range scale data ####
#----------------------------------#


#------------------------------------#
# analyze buffer/local scale data ####
#------------------------------------#
