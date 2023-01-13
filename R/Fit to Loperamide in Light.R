## ---------------------------
##
## Script Name: Fit to Loperamide in Light
##
## Purpose of Script: Identify optimal linear mixed effects model for
##                    Loperamide Zebrafish behavioral data in Light.
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-09
##
## Email: Rowson.Zachary@epa.gov
##
## Working Directory:"/ccte/home2/zrowson/Desktop/Stephanie Projects/gabi/Regression Endpoints"
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

library(nlme)
library(data.table)

# load data
load("data/Padilla_DNT60_Loperamide.rda")

# Isolate Light data
rmv <- Padilla_DNT60_Loperamide[is.na(y), unique(fishID)]
Padilla_DNT60_Loperamide.L <- Padilla_DNT60_Loperamide[!(fishID%in%rmv) & (t%in%11:30)]

# Specify hierarchical data structure
nest <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.L)
attach(nest)
nest$t <- t - min(t)
nest$conc <- as.factor(conc)


# Choose Model ------------------------------------------------------------


# Fit base gls model. No random effect
lm1 <- gls(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "ML")

# Add correlation structures
lm2.AR1 <- update(lm1, corr = corAR1(,form=~1|fishID))
lm2.ARMA <- update(lm1, corr = corARMA(,form=~1|fishID,p=1,q=1))

# Compare models
anova(lm1, lm2.AR1, lm2.ARMA) # ARMA model has min(AIC)

# Add random effect of fish and evaluate if preferred correlation structure changes
lmm1.AR1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "ML",
                corr = corAR1(,form=~1|fishID), random = ~1|fishID)
lmm1.ARMA <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "ML",
                 corr = corARMA(,form=~1|fishID,p=1,q=1), random = ~1|fishID)

anova(lm2.AR1, lm2.ARMA, lmm1.AR1, lmm1.ARMA) # AIC(AR1) - AIC(ARMA) = -2

# Refit with REML, compare, and evaluate random effects
lmm1.AR1.REML <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                     corr = corAR1(,form=~1|fishID), random = ~1|fishID)
lmm1.ARMA.REML <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                      corr = corARMA(,form=~1|fishID,p=1,q=1), random = ~1|fishID)

anova(lmm1.AR1.REML, lmm1.ARMA.REML) # Fitting by remal does not change difference but AIC is higher?

plot(ranef(lmm1.AR1.REML))
plot(ranef(lmm1.ARMA.REML))

# Choose AR1 structure as it requires less parameters

# fit with contrasts
lmm1.AR1.REML.cntrsts <- lme(y ~ conc*poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                     corr = corAR1(,form=~1|fishID), random = ~1|fishID)
summary(lmm1.AR1.REML.cntrsts)$tTable

lme(y ~ conc*poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
    corr = corAR1(,form=~1|fishID), random = ~1|fishID)

rm(list = ls())


# Exploration ---------------------------------------------------------

# Try removing quadratic term. Try adding cubic term.
lmm2 <- lme(y ~ 0 + conc + conc:t, data = nest, method = "ML",
            corr = corAR1(,form=~1|fishID), random = ~1|fishID)
lmm3 <- lme(y ~ 0 + conc + conc:poly(t,degree=3,raw=TRUE), data = nest, method = "ML",
            corr = corAR1(,form=~1|fishID), random = ~1|fishID)

anova(lmm1.AR1, lmm2, lmm3) # poly2 wins

# try adding plate as a fixed effect
lmm1.AR1.apid <- update(lmm1.AR1, .~. + apid)

anova(lmm1.AR1, lmm1.AR1.apid) # addition of apid increases AIC

# try modeling plate as a random effect
nest2 <- groupedData(y ~ conc|apid/fishID, data =  Padilla_DNT60_Loperamide.L)
lmm4 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest2, method = "REML",
            corr = corAR1(,form=~1|fishID), random = ~1|fishID)
lmm4.rapid <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest2, method = "REML",
            corr = corAR1(,form=~1|apid/fishID), random = ~1|apid/fishID)

anova(lmm4, lmm4.rapid) # AIC decreases, but I may be setting utilizing nest data incorrectly

# try with original nest
lmm5.rapid <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                  corr = corAR1(,form=~1|apid/fishID), random = ~1|apid/fishID)

anova(lmm1.AR1.REML, lmm5.rapid) # AIC changes by only a little, I was setting up nested structure incorrectly


# Random slope ------------------------------------------------------------


# Plot slopes of individual for vehicle control on DNT-50
xyplot(y ~ t | fishID, Padilla_DNT60_Loperamide.L[apid=="DNT-50" & wllt=="v"],
       type = c("g","p","r"),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "time",
       ylab = "Total movement (cm) per 2 min", aspect = "xy",
       main = "DNT-50 vehicle control")

# Set control setting to make sure model converges
ctrl <- lmeControl(opt="optim")
lmm1.AR1.REML.rSlp <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                     corr = corAR1(,form=~1|fishID), random = ~1 + t|fishID, control = ctrl)

anova(lmm1.AR1.REML, lmm1.AR1.REML.rSlp) # random slope lowers AIC

plot(ranef(lmm1.AR1.REML.rSlp))
plot(ranef(lmm1.AR1.REML))

lmm1.AR1.REML.rPoly1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                            corr = corAR1(,form=~1|fishID), random = ~1 + poly(t,degree=1,raw=TRUE)|fishID, control = ctrl)
lmm1.AR1.REML.rPoly2 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                            corr = corAR1(,form=~1|fishID), random = ~1 + poly(t,degree=2,raw=TRUE)|fishID, control = ctrl)
anova(lmm1.AR1.REML, lmm1.AR1.REML.rSlp, lmm1.AR1.REML.rPoly1, lmm1.AR1.REML.rPoly2)
