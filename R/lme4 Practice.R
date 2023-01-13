## ---------------------------
##
## Script Name: Practice Using lme4
##
## Purpose of Script: Learn how to use lme4 to create a repeated measures linear model with
##                    random effects of individual
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-16
##
## Email: Rowson.Zachary@epa.gov
##
## Working Directory: "/ccte/home2/zrowson/Desktop/Stephanie Projects/gabi/Regression Endpoints"
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

library(lme4)


# Practice with Orthodontist data ------------------------------------------


data(sleepstudy)
attach(sleepstudy)
str(sleepstudy)

# Plot linear regressions of individuals
require(lattice)
xyplot(Reaction ~ Days | Subject, sleepstudy, type = c("g","p","r"),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "Days of sleep deprivation",
       ylab = "Average reaction time (ms)", aspect = "xy")

DNT50.L <- Padilla_DNT60_Loperamide[apid=="DNT-50" & t%in%11:30]

xyplot(y ~ t | fishID, DNT50.L[wllt=="v"], type = c("g","p","smooth"),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "time",
       ylab = "Total movement (cm) per 2 min", aspect = "xy",
       main = "DNT-50 vehicle control")
xyplot(y ~ t | fishID, Padilla_DNT60_Loperamide[wllq ==1 & conc==0.4], type = c("g","p","r"),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "time",
       ylab = "Total movement (cm) per 2 min", aspect = "xy",
       main = "DNT-50 Loperamide 0.4 uM")

