## ---------------------------
##
## Script Name: Visualize the data
##
## Purpose of Script: Visualize Loperemide Exposure data to get an idea of what random effects
##                    best describe the data.
##
## Author: Zachary Rowson
##
## Date Created: 2022-03-14
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


library(data.table)
library(ggplot2)


# Visualize data ----------------------------------------------------


# Light
data <- copy(Padilla_DNT60_Loperamide.L)
x <- sample(data[wllt=="t",unique(fishID)], 12)
ggplot(data[conc==0.004], aes(x=t,y=y)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm", formula = y ~ poly(x,3), se=FALSE, color="red", lwd = .75) +
  facet_wrap(vars(fishID)) # Intercept seems to be fairly constant

# Dark
data <- copy(Padilla_DNT60_Loperamide.D)
x <- sample(data[wllt=="t",unique(fishID)], 12)

ggplot(data[fishID%in%x], aes(x=t,y=y)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(fishID)) # Everything seems to vary here
