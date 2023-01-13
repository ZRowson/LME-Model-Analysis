## ---------------------------
##
## Script name: Evaluate Autocorrelation lag of Zebrafish data
##
## Purpose of script: Evaluate the lag of zebrafish time series data
##                    for AR(k) or ARMA(k) modelling of covariance.
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-07
##
## Email: Rowson.Zachary@epa.gov
##
## Working directory: "/ccte/home2/zrowson/Desktop/Stephanie Projects/gabi/Regression Endpoints"
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

library(data.table)
library(stats)

# load data

load("data/Padilla_DNT60_Loperamide.rda")
ctrl.data <- Padilla_DNT60_Loperamide[wllt == "v"]

attach(ctrl.data)

# sample one observation / fish

obs <- sample(fishID, 1)
obs.data <- ctrl.data[fishID == obs]

# save data as a time series object

ts <- ts(obs.data[t%in%11:30,y], start = 11)

# plot autocorrelation function

acf(ts, main = paste0("One observation of Loperamide control (",  obs, ")"))

# plot autocorrelation of means

means <- tapply(ctrl.data[,y], ctrl.data[,t], mean, na.rm = TRUE)

ts.means <- ts(means[11:30], start = 11)
acf(ts.means, main = "Loperamide Control Means in Light")

