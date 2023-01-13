## ---------------------------
##
## Script Name: Fit to Loperamide in Dark
##
## Purpose of Script: Identify optimal linear mixed effects model for
##                    Loperamide Zebrafish behavioral data in Dark.
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-14
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

# Isolate Dark data first
rmv <- Padilla_DNT60_Loperamide[is.na(y), unique(fishID)]
Padilla_DNT60_Loperamide.D <- Padilla_DNT60_Loperamide[!(fishID%in%rmv) & (t%in%31:50)]

# Specify hierarchical data structure
nest <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.D)
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

anova(lm2.AR1, lm2.ARMA, lmm1.AR1, lmm1.ARMA) # AIC(lmm1.AR1) - AIC(1mm.ARMA) = ~ 30

# Refit with REML, compare, and evaluate random effects
lmm1.AR1.REML <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                     corr = corAR1(,form=~1|fishID), random = ~1|fishID)
lmm1.ARMA.REML <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                      corr = corARMA(,form=~1|fishID,p=1,q=1), random = ~1|fishID)

anova(lmm1.AR1.REML, lmm1.ARMA.REML) # Fitting by REML does not change difference but AIC is higher?

plot(ranef(lmm1.AR1.REML))
plot(ranef(lmm1.ARMA.REML))

# Choose ARMA structure. Try ARMA with new parameters
lmm1.ARMA22 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "ML",
                 corr = corARMA(,form=~1|fishID,p=2,q=2), random = ~1|fishID)
lmm1.ARMA12 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "ML",
                   corr = corARMA(,form=~1|fishID,p=1,q=2), random = ~1|fishID)
lmm1.ARMA21 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "ML",
                   corr = corARMA(,form=~1|fishID,p=2,q=1), random = ~1|fishID)


anova(lmm1.ARMA, lmm1.ARMA22, lmm1.ARMA21,lmm1.ARMA12) # ARMA with more parameters doesn't have much of an effect on AIC, but results in significant likelihood ratio test

# Random slope ------------------------------------------------------------

xyplot(y ~ t | fishID, Padilla_DNT60_Loperamide.D[apid=="DNT-50" & wllt=="v"],
       type = c("g","p","r"),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "time",
       ylab = "Total movement (cm) per 2 min", aspect = "xy",
       main = "DNT-50 vehicle control")

ctrl <- lmeControl(maxIter = 200, msMaxIter = 200)
lmm1.ARMA.REML.rPoly1 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                      corr = corARMA(,form=~1|fishID,p=1,q=1), random = ~1+poly(t,degree=1,raw=TRUE)|fishID)
lmm1.ARMA.REML.rPoly2 <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest, method = "REML",
                             corr = corARMA(,form=~1|fishID,p=1,q=1), random = ~1+poly(t,degree=2,raw=TRUE)|fishID,
                             control = ctrl)

anova(lmm1.ARMA.REML, lmm1.ARMA.REML.rPoly1, lmm1.ARMA.REML.rPoly2) # So no random poly terms in Dark
# I worry that the choice of random slope in Light has to do with exposure. Thus for a different chemical we will find random
# slope in Dark does improve model.
