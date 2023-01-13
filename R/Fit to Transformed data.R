## ---------------------------
##
## Script Name: Fit LMEs to transformed data
##
## Purpose of Script: Fit LMEs in both Dark and Light to Loperamide exposure data.
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


# Specify hierarchical data structure
nest.L <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.L)
nest.L$t <- nest.L$t - min(nest.L$t)
nest.L$conc <- as.factor(nest.L$conc)

nest.D <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.D)
nest.D$t <- nest.D$t - min(nest.D$t)
nest.D$conc <- as.factor(nest.D$conc)


# Determine correlation structure ------------------------------------------------------


# Optimal correlation structure might change for transformed data. So start from scratch and work back up

# LIGHT
# Fit base gls model. No random effect
lm1 <- gls(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.L, method = "ML")

# Add correlation structures
lm1.AR1 <- update(lm1, corr = corAR1(,form=~1|fishID))
lm1.ARMA <- update(lm1, corr = corARMA(,form=~1|fishID,p=1,q=1))

# Compare models
anova(lm1, lm1.AR1, lm1.ARMA) # AR1 model has minimum AIC

# Add random effect of fish and evaluate if preferred correlation structure changes
lmm1.AR1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
                corr = corAR1(,form=~1|fishID),
                random = ~1|fishID,
                method = "ML",
                data = nest.L)
lmm1.ARMA <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
                 corr = corARMA(,form=~1|fishID,p=1,q=1),
                 random = ~1|fishID,
                 method = "ML",
                 data = nest.L)

anova(lm1.AR1, lmm1.AR1, lmm1.ARMA) # ARMA structure is best but very barely

# Take-away: AR1 appears best but ARMA will likely be the best structure for the Dark data. If we want to
# construct a model that can be fit to both Light and Dark data, using a universal correlation structure
# will be necessary.

# DARK
# Fit base gls model. No random effect
lm2 <- gls(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.D, method = "ML")

# Add correlation structures
lm2.AR1 <- update(lm2, corr = corAR1(,form=~1|fishID))
lm2.ARMA <- update(lm2, corr = corARMA(,form=~1|fishID,p=1,q=1))

# Compare models
anova(lm2, lm2.AR1, lm2.ARMA) # AR1 model is much better

# Add random effect of fish and evaluate if preferred correlation structure changes
lmm2.AR1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
                corr = corAR1(,form=~1|fishID),
                random = ~1|fishID,
                method = "ML",
                data = nest.D)
lmm2.ARMA <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
                 corr = corARMA(,form=~1|fishID,p=1,q=1),
                 random = ~1|fishID,
                 method = "ML",
                 data = nest.D)

anova(lm2.AR1, lmm2.AR1, lmm2.ARMA) # AR1 structure is best

# This is good news, because I don't know if the ARMA model could be fit XD
# So we will use an AR1 structure across Light and Dark periods


# Determine random effects ------------------------------------------------


# To visualize slopes and intercepts look at data for 16 randomly selected fish
# Plot slopes of individual for vehicle control on DNT-50
x <- sample(Padilla_DNT60_Loperamide[wllq==1,unique(fishID)],16)
lattice::xyplot(y ~ t | fishID, Padilla_DNT60_Loperamide.L[fishID%in%x],
                type = c("g","p","r"),
                index = function(x,y) coef(lm(y ~ x))[1],
                xlab = "time",
                ylab = "Total movement (cm) per 2 min", aspect = "xy",
                main = "Random Sample of 16 Fish in Light from Loperamide Exposure")
lattice::xyplot(y ~ t | fishID, Padilla_DNT60_Loperamide.D[fishID%in%x],
                type = c("g","p","r"),
                index = function(x,y) coef(lm(y ~ x))[1],
                xlab = "time",
                ylab = "Total movement (cm) per 2 min", aspect = "xy",
                main = "Random Sample of 16 Fish in Dark from Loperamide Exposure")
# From these graphics it seems that randomness of slope may be due to chemical exposure rather than
# within individual randomness
#
# Look at random sample of control
z <- sample(Padilla_DNT60_Loperamide[wllq==1&wllt=="v",unique(fishID)],16)
lattice::xyplot(y ~ t | fishID, Padilla_DNT60_Loperamide.L[fishID%in%z],
                type = c("g","p","r"),
                index = function(x,y) coef(lm(y ~ x))[1],
                xlab = "time",
                ylab = "Total movement (cm) per 2 min", aspect = "xy",
                main = "Random Sample of 16 Control Fish in Light from Loperamide Exposure")
lattice::xyplot(y ~ t | fishID, Padilla_DNT60_Loperamide.D[fishID%in%z],
                type = c("g","p","r"),
                index = function(x,y) coef(lm(y ~ x))[1],
                xlab = "time",
                ylab = "Total movement (cm) per 2 min", aspect = "xy",
                main = "Random Sample of 16 Fish in Dark from Loperamide Exposure")

# Nevermind, maybe randomness of slope is accurate

# LIGHT
# Build off of previous optimal model
ctrl <- lmeControl(opt = "optim",
                   maxIter = 1000,
                   msMaxIter = 1000,
                   niterEM = 1000,
                   msTol = 1e-20,
                   msVerbose = TRUE,
                   tolerance = 1e-20,
                   allow.n.lt.q = TRUE)
lmm1.1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
                corr = corAR1(,form=~1|fishID),
                random = ~1|fishID,
                method = "ML",
                data = nest.L,
                control = ctrl)
lmm1.1$apVar
lmm1.2 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
              corr = corAR1(,form=~1|fishID),
              random = ~ poly(t,degree=1,raw=TRUE)|fishID,
              method = "ML",
              data = nest.L,
              control = ctrl)
lmm1.2$apVar
anova(lmm1.1, lmm1.2) # Random intercept is significantly better
lmm1.3 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
              corr = corAR1(,form=~1|fishID),
              random = ~ poly(t,degree=2,raw=TRUE),
              method = "ML",
              data = nest.L,
              control = ctrl)
lmm1.3$apVar # Non-positive definite

# Tomorrow work on https://bbolker.github.io/mixedmodels-misc/notes/corr_braindump.html



# Load and transform data -------------------------------------------------


# Specify hierarchical data structure
nest.L1 <- groupedData(y ~ conc|fishID, data =  as.data.frame(Padilla_DNT60_Loperamide.L))
nest.L1$t <- nest.L1$t - min(nest.L1$t)
nest.L1$t <- nest.L1$t / max(nest.L1$t)
nest.L$conc <- as.factor(nest.L$conc)

# Determine correlation structure ------------------------------------------------------

# Optimal correlation structure might change for transformed data. So start from scratch and work back up

# LIGHT
# Fit base gls model. No random effect
lm1 <- gls(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.L1, method = "ML")

# Add correlation structures
lm1.AR1 <- update(lm1, corr = corAR1(,form=~1|fishID))
lm1.ARMA <- update(lm1, corr = corARMA(,form=~1|fishID,p=1,q=1))

# Compare models
anova(lm1, lm1.AR1, lm1.ARMA) # ARMA best by far

# Add random effect of fish and evaluate if preferred correlation structure changes
lmm1.AR1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
                corr = corAR1(,form=~1|fishID),
                random = ~1|fishID,
                method = "ML",
                data = nest.L)
lmm1.ARMA <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
                 corr = corARMA(,form=~1|fishID,p=1,q=1),
                 random = ~1|fishID,
                 method = "ML",
                 data = nest.L)

anova(lm1.AR1, lmm1.AR1, lmm1.ARMA) # ARMA structure is best but very barely. Stick to AR1


# Determine random effects ------------------------------------------------


# To visualize slopes and intercepts look at data for 16 randomly selected fish
# Plot slopes of individual for vehicle control on DNT-50
x <- sample(Padilla_DNT60_Loperamide[wllq==1,unique(fishID)],16)
lattice::xyplot(y ~ t | fishID, nest.L1[nest.L1$fishID %in% x,],
                type = c("g","p","r"),
                index = function(x,y) coef(lm(y ~ x))[1],
                xlab = "time",
                ylab = "Total movement (cm) per 2 min", aspect = "xy",
                main = "Random Sample of 16 Fish in Light from Loperamide Exposure")
lattice::xyplot(y ~ t | fishID, nest.L1[nest.L1$fishID %in% x,],
                type = c("g","p","r"),
                index = function(x,y) coef(lm(y ~ x))[1],
                xlab = "time",
                ylab = "Total movement (cm) per 2 min", aspect = "xy",
                main = "Random Sample of 16 Fish in Dark from Loperamide Exposure")
# From these graphics it seems that randomness of slope may be due to chemical exposure rather than
# within individual randomness
#
# Look at random sample of control
z <- sample(Padilla_DNT60_Loperamide[wllq==1&wllt=="v",unique(fishID)],16)
lattice::xyplot(y ~ t | fishID, Padilla_DNT60_Loperamide.L[fishID%in%z],
                type = c("g","p","r"),
                index = function(x,y) coef(lm(y ~ x))[1],
                xlab = "time",
                ylab = "Total movement (cm) per 2 min", aspect = "xy",
                main = "Random Sample of 16 Control Fish in Light from Loperamide Exposure")
lattice::xyplot(y ~ t | fishID, Padilla_DNT60_Loperamide.D[fishID%in%z],
                type = c("g","p","r"),
                index = function(x,y) coef(lm(y ~ x))[1],
                xlab = "time",
                ylab = "Total movement (cm) per 2 min", aspect = "xy",
                main = "Random Sample of 16 Fish in Dark from Loperamide Exposure")

# Nevermind, maybe randomness of slope is accurate

# LIGHT
# Build off of previous optimal model
ctrl <- lmeControl(opt = "optim",
                   maxIter = 1000,
                   msMaxIter = 1000,
                   niterEM = 1000,
                   msTol = 1e-20,
                   msVerbose = TRUE,
                   tolerance = 1e-20,
                   allow.n.lt.q = TRUE)
lmm1.1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
              corr = corAR1(,form=~1|fishID),
              random = ~1|fishID,
              method = "ML",
              data = nest.L,
              control = ctrl)
lmm1.test <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
              random = ~1|fishID,
              method = "ML",
              data = nest.L,
              control = ctrl)
lmm1.1$apVar
lmm1.2 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
              corr = corAR1(),
              random = ~ poly(t,degree=1,raw=TRUE)|fishID,
              method = "ML",
              data = nest.L,
              control = ctrl)
lmm1.2$apVar
anova(lmm1.1, lmm1.2) # Random intercept is significantly better
lmm1.3 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
              corr = corAR1(,form=~1|fishID),
              random = ~ poly(t,degree=2,raw=TRUE),
              method = "ML",
              data = nest.L,
              control = ctrl)
lmm1.3$apVar # Non-positive definite


# Mess with correlation structure -----------------------------------------


# Specify hierarchical data structure
nest.L.2 <- groupedData(y ~ t | fishID, data =  Padilla_DNT60_Loperamide.L)
nest.L.2$t <- nest.L.2$t - min(nest.L.2$t)

# Plot initial model
ctrl <- lmeControl(opt = "optim")
lmm1.1.1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE),
              random = ~1 | fishID,
              method = "ML",
              data = nest.L.2,
              control = ctrl)
plot(ACF(lmm1.1.1), alpha = 0.05)

lmm1.1.2 <- update(lmm1.1.1, correlation = corAR1(, ~ t | fishID))
plot(ACF(lmm1.1.2, resType = "normalized"), alpha = 0.05)

lmm1.1.3 <- update(lmm1.1.1, correlation = corARMA(p=1,q=1))

# For some reason ARMA structure isn't working
# Lets simplify model
nest.L.S <- groupedData(y ~ t | fishID, data =  Padilla_DNT60_Loperamide.L[conc==0])
nest.L.S$t <- nest.L.S$t - min(nest.L.S$t)

# Plot initial model
ctrl <- lmeControl(opt = "optim")
lmm1.S <- lme(y ~ poly(t,degree=2,raw=TRUE),
                random = ~1 | fishID,
                method = "ML",
                data = nest.L.S,
                control = ctrl)
plot(ACF(lmm1.S), alpha = 0.05)

lmm1.S.2 <- update(lmm1.S, correlation = corAR1())
plot(ACF(lmm1.S.2, resType = "normalized"), alpha = 0.05)

lmm1.S.3 <- update(lmm1.S.2, correlation = corARMA(, ~ t|fishID, p=1))
plot(ACF(lmm1.S.2, resType = "normalized"), alpha = 0.05)

# ARMA still doesn't work
# Try with Ovary data
ctrl <- lmeControl(opt = "optim")
lmm.O <- lme(follicles ~ sin(2*pi*Time) + cos(2*pi*Time),
             data = Ovary,
             random = pdDiag(~sin(2*pi*Time)),
             control = ctrl)
plot(ACF(lmm.O), alpha = 0.05)

# Add correlation
lmm.O.1 <- update(lmm.O, correlation = corAR1())
plot(ACF(lmm.O.1, resType = "normalized"), alpha = 0.05)

lmm.O.2 <- update(lmm.O, correlation = corARMA(p=1,q=1))
plot(ACF(lmm.O.2, resType = "normalized"), alpha = 0.05) # Doesn't work with optim

# Go back to original model and remove optim
lmm1.S <- lme(y ~ 1 , # Add polynomail term
              random = ~ poly(t,degree=1) | fishID,
              method = "ML",
              data = nest.L.S)
plot(ACF(lmm1.S), alpha = 0.05)

lmm1.S.2 <- update(lmm1.S, correlation = corAR1())
plot(ACF(lmm1.S.2, resType = "normalized"), alpha = 0.05)

lmm1.S.3 <- update(lmm1.S.2, correlation = corARMA(, ~ t|fishID, p=1, q=1))
plot(ACF(lmm1.S.3, resType = "normalized"), alpha = 0.05)

lmm1.S.4 <- update(lmm1.S.2, correlation = corARMA(, ~ t|fishID, p=2,q=2))
plot(ACF(lmm1.S.4, resType = "normalized"), alpha = 0.05)

anova(lmm1.S.2,lmm1.S.3,lmm1.S.4)


# Go back to original data and reformat -----------------------------------


# Set up grouped structure
nest.L.S <- groupedData(y ~ t | fishID, data =  Padilla_DNT60_Loperamide.L[conc==0])
nest.L.S$t <- nest.L.S$t - mean(nest.L.S$t)

ctrl <- lmeControl(opt = "optim", niterEM = 500, msMaxIter = 1000, msTol = 1e-10)
lmm1.0 <- lme(y ~ poly(t,degree=2,raw=TRUE),
                random = ~ 1 | fishID,
                correlation = corAR1(, ~ t | fishID),
                method = "ML",
                data = nest.L.S,
                control = ctrl)
plot(augPred(lmm1.S))

lmm1.0.1 <- lme(y ~ poly(t,degree=2,raw=TRUE),
              random = ~ poly(t,degree=1,raw=TRUE) | fishID,
              correlation = corAR1(, ~ t | fishID),
              method = "ML",
              data = nest.L.S,
              control = ctrl)
lmm1.0.1$apVar
plot(augPred(lmm1.S))

lmm1.0.2 <- lme(y ~ poly(t,degree=2,raw=TRUE),
                random = ~ poly(t,degree=2,raw=TRUE) | fishID,
                correlation = corAR1(, ~ t | fishID),
                method = "ML",
                data = nest.L.S,
                control = ctrl)
lmm1.0.2$apVar
anova(lmm1.0.1, lmm1.0)

# Try pdSymm to see if model changes
ctrl <- lmeControl(opt = "optim", niterEM = 500, msMaxIter = 10000, msTol = 1e-20, tolerance = 1e-20, allow.n.lt.q = TRUE)
lmm1.0.Symm <- lme(y ~ poly(t,degree=2,raw=TRUE),
                random = reStruct(~1|fishID,REML=FALSE,pdClass="pdSymm"),
                correlation = corAR1(, ~ t | fishID),
                method = "ML",
                data = nest.L.S,
                control = ctrl)
plot(augPred(lmm1.0.Symm))
lmm1.0.1.Symm <- lme(y ~ poly(t,degree=2,raw=TRUE),
                random = reStruct(~poly(t,degree=1,raw=TRUE)|fishID,REML=FALSE,pdClass="pdSymm"),
                correlation = corAR1(, ~ t | fishID),
                method = "ML",
                data = nest.L.S,
                control = ctrl)
lmm1.0.1.Symm$apVar
plot(augPred(lmm1.0.1.Symm))

lmm1.0.2.Symm <- lme(y ~ poly(t,degree=2),
                random = reStruct(~poly(t,degree=2)|fishID,REML=FALSE,pdClass="pdSymm"),
                correlation = corAR1(, ~ t | fishID),
                method = "ML",
                data = nest.L.S,
                control = ctrl)
lmm1.0.2.Symm$apVar
plot(augPred(lmm1.0.2.Symm))
anova(lmm1.0.1.Symm, lmm1.0.Symm)


# Fit to all Data ---------------------------------------------------------


# Data structure
nest.L <- groupedData(y ~ t | fishID, data =  Padilla_DNT60_Loperamide.L)
nest.L.S$t <- nest.L$t - mean(nest.L$t)

# Try pdSymm
ctrl <- lmeControl(opt = "optim", niterEM = 500, msMaxIter = 10000, msVerbose = TRUE, msTol = 1e-20, tolerance = 1e-20, allow.n.lt.q = TRUE)
lmm1 <- lme(y ~ poly(t,degree=2,raw=TRUE),
                   random = reStruct(~1|fishID,REML=FALSE,pdClass="pdSymm"),
                   correlation = corAR1(, ~ t | fishID),
                   method = "ML",
                   data = nest.L,
                   control = ctrl)
pred1 <- augPred(lmm1)
x <- sample(unique(pred1$.groups),100)
plot(pred1[pred1$.groups%in%x,])

lmm1.1 <- lme(y ~ poly(t,degree=2,raw=TRUE),
                     random = reStruct(~poly(t,degree=1,raw=TRUE)|fishID,REML=FALSE,pdClass="pdSymm"),
                     correlation = corAR1(, ~ t | fishID),
                     method = "ML",
                     data = nest.L,
                     control = ctrl)
lmm1.1$apVar
sample(unique(augPred(lmm1.1)$fishID),10)

lmm1.0.2.Symm <- lme(y ~ poly(t,degree=2),
                     random = reStruct(~poly(t,degree=2)|fishID,REML=FALSE,pdClass="pdSymm"),
                     correlation = corAR1(, ~ t | fishID),
                     method = "ML",
                     data = nest.L.S,
                     control = ctrl)
lmm1.0.2.Symm$apVar
plot(augPred(lmm1.0.2.Symm))
anova(lmm1.0.1.Symm, lmm1.0.Symm)
