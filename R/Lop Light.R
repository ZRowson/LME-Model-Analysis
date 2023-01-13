## ---------------------------
##
## Script Name: Loperamide in Light
##
## Purpose of Script: Fit LME to Loperamide in Light
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


library(nlme)
library(data.table)


# Nest data ---------------------------------------------------------------


nestLop.L <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.L[conc==0],
                         order.groups = FALSE)

# Set beginning of light to t=0. Set concentration as a factor
nestLop.L$t <- nestLop.L$t - min(nestLop.L$t)
nestLop.L$conc <- as.factor(nestLop.L$conc)


# Start with what you know -------------------------------------------------


ctrl <- lmeControl(opt = "optim", niterEM = 1000, msTol = 1e-40, msVerbose = TRUE,
                   minAbsParApVar = 1e-40, maxIter = 1000, msMaxIter = 1000,
                   control = list(abstol=1e-40,reltol=1e-40,ndeps=1e-20))
lmm1 <- lme(y ~  conc*poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = ~1|fishID,
            method = "ML",
            control = ctrl,
            data = nestLop.L)
lmm2 <- lme(y ~ conc*poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = ~poly(t,degree=1,raw=TRUE) | fishID,
            method = "ML",
            control = ctrl,
            data = nestLop.L)
lmm3 <- lme(y ~ conc*poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = ~poly(t,degree=2,raw=TRUE) |fishID,
            method = "ML",
            control = ctrl,
            data = nestLop.L)
lmm1$apVar
lmm2$apVar
lmm3$apVar
anova(lmm1,lmm2, lmm3) # lmm2 much better

# Attempt to include a concentration dependent variance of random effects
lmm1.1 <- update(lmm1, .~., random=reStruct(~ 0 + conc|fishID, REML=FALSE, pdClass="pdDiag"), control = ctrl)
# lmm2.1 <- update(lmm1, .~., random=reStruct(~ 0 + conc*poly(t,degree=1,raw=TRUE)|fishID, REML=FALSE, pdClass="pdDiag"), control = ctrl)
# lmm3.1 <- update(lmm1, .~., random=reStruct(~ 0 + conc*poly(t,degree=2,raw=TRUE)|fishID, REML=FALSE, pdClass="pdDiag"), control = ctrl)
ctrl <- lmeControl(opt = "optim", niterEM = 400, msTol = 1e-10, msVerbose = TRUE,
                   minAbsParApVar = 1e-10, maxIter = 400, msMaxIter = 400,
                   control = list(abstol=1e-10,reltol=1e-10,ndeps=1e-10))
lmm1.2 <- update(lmm1, .~., random=reStruct(~ 0 + poly(t,degree=1,raw=TRUE)|conc, REML=FALSE, pdClass="pdDiag"),
                 correlation = corAR1(, form = ~t |fishID),
                 control = ctrl)


# Try something new -------------------------------------------------------


# Fitting models isn't going too well. Try fitting regressions to each concentration group separately
# MAKE SURE TO NOT LOAD data.table
ctrl <- lmeControl(opt = "optim", niterEM = 10000, msTol = 1e-10, msVerbose = TRUE,
                   minAbsParApVar = 1e-10, maxIter = 1000, msMaxIter = 1000,
                   control = list(abstol=1e-10,reltol=1e-10,ndeps=1e-10))
# Control
nest0 <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.L[conc==0],
                         order.groups = FALSE)
nest0$t <- nestLop.L$t - min(nestLop.L$t)
nest0$conc <- as.factor(nestLop.L$conc)
# Model lme
lmm0 <- lme(y ~ poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = ~ 0 + poly(t,degree=1,raw=TRUE) |fishID,
            method = "ML",
            control = ctrl,
            data = nest0)
lmm0$apVar

# 0.004
nest1 <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.L[conc==0.004],
                     order.groups = FALSE)
nest1$t <- nest1$t - 20
# Model lme
lmm1 <- lme(y ~ poly(t,degree=2,raw=TRUE),
            correlation = corAR1(, form = ~ t | fishID),
            random = ~ 0 + poly(t,degree=1,raw=TRUE) |fishID,
            method = "ML",
            control = ctrl,
            data = nest1)
lmm1$apVar
