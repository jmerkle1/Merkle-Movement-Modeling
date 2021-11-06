#----------------------------------------------#
#- Step Selection Functions ----------------####
#---Functional Responses in habitat selection--#
#--- model validation ----------------------####
#--- plotting predicted relationships ---------#
#---------- Jerod Merkle ----------------------#
#------------- Lab 7 --------------------------#

setwd("C:/Users/jmerkle/Dropbox/school_work/Wyoming/Classes/Mov_modeling")
#load up QIC code
source(paste(getwd(),"/MiscFunctions/QIC.coxph.R",sep=""))

load("SSFdata.RData")   #load up your dataframe from Lab5
library(survival)

#fit your best model
best <- clogit(case~dist+dist_spline1+dist_spline2+dist_spline3+
               escapeTerr_delta+landcovbarren+landcovherbaceous+landcovshrubland+
               strata(strata)+cluster(id_seas_yr),method = "efron",data=seas)

summary(best)

#----------------------------------------------------#
# create availability variables for each variable ####
# and test for a functional response -------------####

#first distance to escape terrain (this is a continuous variable)
hist(seas$escapeTerr_delta)
hist(seas$escapeTerr_target) # pay attention to the histogram of variables (here, I had divided by 1000)
seas$escapeTerr_target <- seas$escapeTerr_target/1000
#create dataframe that provides the mean value of your variable for each strata
avail <- as.data.frame(tapply(seas$escapeTerr_target, seas$strata, mean))  
names(avail) <- "escapeTerr_avail"
avail$strata <- as.numeric(rownames(avail)) #add a column for strata so you can merge it back with your regular dataframe
head(avail)
seas <- merge(seas, avail)  #merge your database with the avail object so you now have a column with mean availability for each strata
seas <- seas[order(seas$strata, seas$case),]   # reorder after the stupid merge function screws it up
head(seas, 8)
hist(seas$escapeTerr_avail)

best_FR <- clogit(case~dist+dist_spline1+dist_spline2+dist_spline3+
                 escapeTerr_delta+escapeTerr_delta:escapeTerr_avail+
                  landcovbarren+landcovherbaceous+landcovshrubland+
                 strata(strata)+cluster(id_seas_yr),method = "efron",data=seas)
summary(best_FR)

cbind(model=c("best","best+funcResponse"),rbind(QIC.coxph(best, details=TRUE),   #QICR is what you want to minimize. Check and make sure n and nevent are the same!
                                    QIC.coxph(best_FR, details=TRUE)))

#second, landcoverbarren (this is a categorical variable)
hist(seas$landcovbarren)
#again, create dataframe that provides the mean value of your variable for each strata
avail <- as.data.frame(tapply(seas$landcovbarren, seas$strata, mean))
names(avail) <- "landcovbarren_avail"
avail$strata <- as.numeric(rownames(avail))
head(avail)
seas <- merge(seas, avail)  #merge your database with the avail object so you now have a column with mean availability for each strata
seas <- seas[order(seas$strata, seas$case),]   # reorder after the stupid merge function screws it up
hist(seas$landcovbarren_avail)

best_FR <- clogit(case~dist+dist_spline1+dist_spline2+dist_spline3+
                    escapeTerr_delta+landcovbarren+landcovbarren:landcovbarren_avail+
                    landcovherbaceous+landcovshrubland+
                    strata(strata)+cluster(id_seas_yr),method = "efron",data=seas)
summary(best_FR)

cbind(model=c("best","best+funcResponse"),rbind(QIC.coxph(best, details=TRUE),   #QICR is what you want to minimize. Check and make sure n and nevent are the same!
                                                QIC.coxph(best_FR, details=TRUE)))


#second, landcovherbaceous (this is a categorical variable)
hist(seas$landcovherbaceous)
avail <- as.data.frame(tapply(seas$landcovherbaceous, seas$strata, mean))
names(avail) <- "landcovherbaceous_avail"
avail$strata <- as.numeric(rownames(avail))
head(avail)
seas <- merge(seas, avail)
hist(seas$landcovherbaceous_avail)

best_FR <- clogit(case~dist+dist_spline1+dist_spline2+dist_spline3+
                    escapeTerr_delta+landcovbarren+
                    landcovherbaceous+landcovherbaceous:landcovherbaceous_avail+
                    landcovshrubland+
                    strata(strata)+cluster(id_seas_yr),method = "efron",data=seas)
summary(best_FR)

cbind(model=c("best","best+funcResponse"),rbind(QIC.coxph(best, details=TRUE),   #QICR is what you want to minimize. Check and make sure n and nevent are the same!
                                                QIC.coxph(best_FR, details=TRUE)))


#fourth, landcovshrubland .... (this is a categorical variable) This functional response was not significant for mine

# if you have multiple functional responses that are significant, 
# try putting multiple of them in the same model. Does it improve the model's QIC?..
# I only had one single functional respones that was supported.


#---------------------#
# model validation ####
#---------------------#

# parameterize your best model that includes your functional response(s)
best <- clogit(case~dist+dist_spline1+dist_spline2+dist_spline3+
                    escapeTerr_delta+escapeTerr_delta:escapeTerr_avail+
                    landcovbarren+landcovherbaceous+landcovshrubland+
                    strata(strata)+cluster(id_seas_yr),
               x=TRUE, y=TRUE,  #these ensure that your dataframe is saved within the model object
               method = "efron",data=seas)
summary(best)

#load up the k-folds code
source(paste(getwd(),"/MiscFunctions/kfoldSSF.R",sep=""))

kfold <- kfoldSSF(best, k=5, nrepet=100, jitter=FALSE, reproducible=TRUE,details=TRUE)

# spearman rank stats for observed steps
mean(kfold$kfold[kfold$type=="obs"])
sd(kfold$kfold[kfold$type=="obs"])
range(kfold$kfold[kfold$type=="obs"])

# spearman rank stats for randomly chosen steps
mean(kfold$kfold[kfold$type=="rand"])
sd(kfold$kfold[kfold$type=="rand"])
range(kfold$kfold[kfold$type=="rand"])


# do cross validation without doing k-folds partitioning (just to understand what is going on)
seas$predvals <- predict(best, type="risk")  #append predicted values to your database
head(seas$predvals)
hist(seas$predvals)   
u <- unique(seas$strata)
#for each strata, rank predvals to determine the predicted rank of the target point
seas <- do.call(rbind, lapply(1:length(u), function(i){
  temp <- seas[seas$strata == u[i],]
  temp <- temp[order(temp$predvals),]
  temp$rank <- 1:nrow(temp)
  return(temp)
}))
seas <- seas[order(seas$strata, seas$case),]

# this histogram shows the frequency of the rankings of used points 
# (or cases) within their respective stratam
hist(seas$rank[seas$case==1], 10)
cor(1:6,table(seas$rank[seas$case==1]),method="spearman")

#now for the controls or available points
hist(seas$rank[seas$case==0])
cor(1:6,table(seas$rank[seas$case==0]),method="spearman")



#-----------------------------------------------------------------------------------#
# plotting predicted relationships based on fitted model (continuous variable)  #####
#-----------------------------------------------------------------------------------#
variable <- "escapeTerr_delta"   #you must change this variable to the one you want!
x <- seq(min(best$x[,variable]),max(best$x[,variable]),length.out=100)
#calculate out baseline selection
baseline <- sum(colMeans(best$x[,colnames(best$x) != variable])*coefficients(best)[names(coefficients(best))!=variable])
var_mean <- exp(baseline + coefficients(best)[names(coefficients(best))==variable]*x)
var_LCI <- exp(baseline + confint(best)[rownames(confint(best))==variable,"2.5 %"]*x)
var_UCI <- exp(baseline + confint(best)[rownames(confint(best))==variable,"97.5 %"]*x)

plot(range(x), range(c(var_LCI,var_UCI)), type="n",bty="l",xlab=variable, ylab="Relative odds of selection")
lines(x,var_LCI, col="grey",lty=2, lwd=3)
lines(x,var_UCI, col="grey",lty=2, lwd=3)
lines(x,var_mean, col="black",lty=1, lwd=3)



#---------------------------------------------------------------------------------------------------------------#
# plotting predicted relationships based on fitted model (continuous variable) including functional response ####
#---------------------------------------------------------------------------------------------------------------#
variable <- "escapeTerr_delta"    #main effect variable
FR <- "escapeTerr_delta:escapeTerr_avail"   #what is the functional response variable you are trying to plot as well
x <- seq(min(best$x[,variable]),max(best$x[,variable]),length.out=100)
avail <- best$x[,FR]/best$x[,variable]   #get availability from the main effect and the interaction term
quants <- quantile(avail, probs = c(0.05,0.5,0.95), na.rm=TRUE)   #get three quantiles for which to plot across
quants
#calculate out baseline selection
baseline <- sum(colMeans(best$x[,colnames(best$x) %in% c(variable,FR) == FALSE])*coefficients(best)[names(coefficients(best))%in% c(variable,FR) == FALSE])
#plotting starts here
plot(range(x), range(c(exp(baseline + confint(best)[rownames(confint(best))==variable,"2.5 %"]*x + confint(best)[rownames(confint(best))==FR,"2.5 %"]*quants[1]*x),
                       exp(baseline + confint(best)[rownames(confint(best))==variable,"97.5 %"]*x + confint(best)[rownames(confint(best))==FR,"97.5 %"]*quants[3]*x))), 
     type="n",bty="l",xlab=variable, ylab="Relative odds of selection")

for(i in 1:length(quants)){
  var_mean <- exp(baseline + coefficients(best)[names(coefficients(best))==variable]*x + 
                    coefficients(best)[names(coefficients(best))==FR]*quants[i]*x)
  var_LCI <- exp(baseline + confint(best)[rownames(confint(best))==variable,"2.5 %"]*x +
                   confint(best)[rownames(confint(best))==FR,"2.5 %"]*quants[i]*x)
  var_UCI <- exp(baseline + confint(best)[rownames(confint(best))==variable,"97.5 %"]*x +
                   confint(best)[rownames(confint(best))==FR,"97.5 %"]*quants[i]*x)
  
  lines(x,var_LCI, col=i,lty=2, lwd=1)
  lines(x,var_UCI, col=i,lty=2, lwd=1)
  lines(x,var_mean, col=i,lty=1, lwd=3)
  
}
legend("topright", c("Low","Med","High"), col=1:3, lty=1,
       lwd=3, bty="n",title="Availability",inset=.1)


#---------------------------------------------------------------------------------------#
# now how to plot predicted values for a dummy variable with the functional response ####
#---------------------------------------------------------------------------------------#
# I must remake my best model to have a functional response with a dummy variable - for illustrative purposes (its not significant)
best <- clogit(case~dist+dist_spline1+dist_spline2+dist_spline3+
                 escapeTerr_delta+landcovbarren+landcovbarren:landcovbarren_avail+
                 landcovherbaceous+landcovshrubland+
                 strata(strata)+cluster(id_seas_yr),method = "efron",data=seas,
               x=TRUE, y=TRUE)
summary(best)

variable <- "landcovbarren"    #main effect variable of interest
FR <- "landcovbarren:landcovbarren_avail"   #what is the functional response variable you are trying to plot as well
avail <- best$x[,FR]/best$x[,variable]   #get availability from the main effect and the interaction term
x <- seq(min(avail,na.rm=TRUE),max(avail,na.rm=TRUE),length.out=100)  #generate x axis (in this case availability of variable)
#calculate out baseline selection
baseline <- sum(colMeans(best$x[,colnames(best$x) %in% c(variable,FR) == FALSE])*coefficients(best)[names(coefficients(best))%in% c(variable,FR) == FALSE])

#plotting starts here
plot(range(x), range(c(exp(baseline + confint(best)[rownames(confint(best))==variable,"2.5 %"] + confint(best)[rownames(confint(best))==FR,"2.5 %"]*x),
                       exp(baseline + confint(best)[rownames(confint(best))==variable,"97.5 %"] + confint(best)[rownames(confint(best))==FR,"97.5 %"]*x))), 
     type="n",bty="l",xlab=paste("Availability of ", variable,sep=""), ylab="Relative odds of selection")

#plot main effect
var_mean <- exp(baseline + coefficients(best)[names(coefficients(best))==variable])
var_LCI <- exp(baseline + confint(best)[rownames(confint(best))==variable,"2.5 %"])
var_UCI <- exp(baseline + confint(best)[rownames(confint(best))==variable,"97.5 %"])

abline(h=var_LCI, col="blue",lty=2, lwd=2)
abline(h=var_UCI, col="blue",lty=2, lwd=2)
abline(h=var_mean, col="blue",lty=1, lwd=3)

#plot functional response
var_mean <- exp(baseline + coefficients(best)[names(coefficients(best))==variable]+
                  coefficients(best)[names(coefficients(best))==FR]*x)
var_LCI <- exp(baseline + confint(best)[rownames(confint(best))==variable,"2.5 %"]+
                 confint(best)[rownames(confint(best))==FR,"2.5 %"]*x)
var_UCI <- exp(baseline + confint(best)[rownames(confint(best))==variable,"97.5 %"]+
                 confint(best)[rownames(confint(best))==FR,"97.5 %"]*x)

lines(x, var_LCI, col="orange",lty=2, lwd=2)
lines(x, var_UCI, col="orange",lty=2, lwd=2)
lines(x, var_mean, col="orange",lty=1, lwd=3)

legend("topleft", c("Main effect","Functional Response"), 
       col=c("blue","orange"), lty=1,lwd=3, bty="n",inset=.1)


#----------------------------------#
# How to plot interaction terms ####
#----------------------------------#
#need model with an interaction term...
best <- clogit(case~dist+dist:escapeTerr_delta+
                 escapeTerr_delta+landcovbarren+
                 landcovherbaceous+landcovshrubland+
                 strata(strata)+cluster(id_seas_yr),method = "efron",data=seas,
               x=TRUE, y=TRUE)
summary(best)

Xvariable <- "escapeTerr_delta"   #the variable you want to see on the x axis...
Yvariable <- "dist"   #the interacting variable
Interaction <- "dist:escapeTerr_delta"   #the actual interaction in the model

x <- seq(min(best$x[,Xvariable]),max(best$x[,Xvariable]),length.out=100)   # get your x axis layed out
InterVals <- best$x[,Yvariable]
quants <- quantile(avail, probs = c(0.05,0.5,0.95), na.rm=TRUE)   #get three quantiles of your interacting variable for which to plot across
quants
#calculate out baseline selection
baseline <- sum(colMeans(best$x[,colnames(best$x) %in% c(Xvariable,Yvariable,Interaction) == FALSE])*
                  coefficients(best)[names(coefficients(best))%in% c(Xvariable,Yvariable,Interaction) == FALSE])
#plotting starts here
plot(range(x), c(0,1000), #MUST CHANGE THESE Y VALUES SO YOUR LINES FIT ON THE PLOT
     type="n",bty="l",xlab=Xvariable, ylab="Relative odds of selection")

for(i in 1:length(quants)){
  var_mean <- exp(baseline + coefficients(best)[names(coefficients(best))==Xvariable]*x + 
                    coefficients(best)[names(coefficients(best))==Yvariable]*quants[i]+
                    coefficients(best)[names(coefficients(best))==Interaction]*quants[i]*x)
  var_LCI <- exp(baseline + confint(best)[rownames(confint(best))==Xvariable,"2.5 %"]*x +
                   confint(best)[rownames(confint(best))==Yvariable,"2.5 %"]*quants[i]+
                   confint(best)[rownames(confint(best))==Interaction,"2.5 %"]*quants[i]*x)
  var_UCI <- exp(baseline + confint(best)[rownames(confint(best))==Xvariable,"97.5 %"]*x +
                   confint(best)[rownames(confint(best))==Yvariable,"97.5 %"]*quants[i]+
                   confint(best)[rownames(confint(best))==Interaction,"97.5 %"]*quants[i]*x)
  
  lines(x,var_LCI, col=i,lty=2, lwd=2)
  lines(x,var_UCI, col=i,lty=2, lwd=2)
  lines(x,var_mean, col=i,lty=1, lwd=3)
  
}
legend("topright", c("Low","Med","High"), col=1:3, 
       lty=1,lwd=3, bty="n",title=Yvariable,inset=.1)






