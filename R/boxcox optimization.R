## ---------------------------
##
## Script Name: Apply Box-Cox to Time Series Data
##
## Purpose of Script: Optimize BoxCox parameters for transformation of movement by time data.
##
## Author: Zachary Rowson
##
## Date Created: 2022-03-22
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
library(MASS)


# Load data ---------------------------------------------------------------


load("data/Padilla_DNT60_Loperamide.rda")
load("data/Padilla_DNT60_pmr0_long.rda")


# Optimize for Loperamide exposure data -----------------------------------


Padilla_DNT60_Loperamide[, t := as.factor(t)][, y := y + 1]
bc_parms <- MASS::boxcox(y ~ t, data = Padilla_DNT60_Loperamide,
                         lambda = seq(-3, 3, by = 0.25),
                         plotit = FALSE)
i <- which(bc_parms$y == max(bc_parms$y))
lam.hat <- bc_parms$x[i]
# Quarter power optimal

# Optimize over all DNT60 data --------------------------------------------


DNT60pmr0_long[, t := as.factor(t)][, y := y + 1]
bc_parms1 <- MASS::boxcox(y ~ t, data = DNT60pmr0_long,
                         lambda = seq(-3, 3, by = 0.25),
                         plotit = FALSE)
i1 <- which(bc_parms1$y == max(bc_parms1$y))
lam.hat1 <- bc_parms1$x[i1]
# log10 transformation optimal
