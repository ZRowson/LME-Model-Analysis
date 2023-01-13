## ---------------------------
##
## Script Name: Re-Evaluate Correlation Structure and Get a Script Ready for Woody
##
## Purpose of Script:
##
## Author: Zachary Rowson
##
## Date Created: 2022-03-24
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
library(nlme)


# Load and transform data -------------------------------------------------


# From boxcox optimization script we found that log10 transformation of data was optimal for all
# DNT60 data
load("data/Padilla_DNT60_Loperamide.rda")
Padilla_DNT60_Loperamide[, y := log10(y+1)]

# Isolate Light and Dark data
rmv <- Padilla_DNT60_Loperamide[is.na(y), unique(fishID)]
Padilla_DNT60_Loperamide.L <- Padilla_DNT60_Loperamide[!(fishID%in%rmv) & (t%in%11:30)]
Padilla_DNT60_Loperamide.D <- Padilla_DNT60_Loperamide[!(fishID%in%rmv) & (t%in%31:50)]

# Specify hierarchical data structures
nest.L.S <- groupedData(y ~ t | fishID, data =  Padilla_DNT60_Loperamide.L[conc==0])
nest.L.S$t <- nest.L.S$t - min(nest.L.S$t)

nest.D.S <- groupedData(y ~ t | fishID, data =  Padilla_DNT60_Loperamide.D[conc==0])
nest.D.S$t <- nest.D.S$t - mean(nest.D.S$t)

nest.L <- groupedData(y ~ t |fishID, data =  Padilla_DNT60_Loperamide.L)
nest.L$t <- nest.L$t - mean(nest.L$t)
nest.L$conc <- as.factor(nest.L$conc)

nest.D <- groupedData(y ~ t |fishID, data =  Padilla_DNT60_Loperamide.D)
nest.D$t <- nest.D$t - mean(nest.D$t)
nest.D$conc <- as.factor(nest.D$conc)


# Evaluate correlation structure ------------------------------------------


# For control sample only
ctrl <- lmeControl(opt = "optim", niterEM = 500, msVerbose = TRUE, msMaxIter = 1000, msTol = 1e-10, allow.n.lt.q = TRUE)
lmm.0 <- lme(y ~ poly(t,degree=2,raw=TRUE),
             random = ~poly(t,degree=1,raw=TRUE)|fishID,
             correlation = corAR1(,~t|fishID),
             method = "ML",
             data = nest.L.S,
             control = ctrl)
plot(ACF(lmm.0, resType = "normalized"), alpha = 0.05)
plot(augPred(lmm.0))

lmm.0.2 <- lme(y ~ poly(t,degree=2,raw=TRUE),
             random = reStruct(~poly(t,degree=2,raw=TRUE)|fishID,REML=TRUE,pdClass="pdSymm"),
             correlation = corAR1(),
             method = "REML",
             data = nest.L.S,
             control = ctrl)
lmm.0.2$apVar
plot(augPred(lmm.0.2))
plot(ACF(lmm.0.2, resType = "normalized"), alpha = 0.05)

lmm.0.3 <- lme(y ~ poly(t,degree=2,raw=TRUE),
               random = ~poly(t,degree=3,raw=TRUE)|fishID,
               data = nest.L.S,
               control = ctrl)
