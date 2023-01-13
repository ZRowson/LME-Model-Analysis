## ---------------------------
##
## Script Name: Models for Dark
##
## Purpose of Script: Fit models to Dark period
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-08
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


# Fit to Loperamide in Dark -------------------------------------------------------

# load data
load("data/Padilla_DNT60_Loperamide.rda")

# Isolate Dark data first
rmv <- Padilla_DNT60_Loperamide[is.na(y), unique(fishID)]
Padilla_DNT60_Loperamide.D <- Padilla_DNT60_Loperamide[!(fishID%in%rmv) & (t%in%31:50)]


# Hierarchical data
nest <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.D)
nest$t <- nest$t - min(nest$t)

# Linear mixed effects model
lmm1 <- lme(y ~ 0 + factor(conc) + factor(conc):poly(t,degree=2,raw=TRUE),
            random = ~1|fishID, method = "REML",
            corr = corARMA(,form=~1|fishID,p=1,q=1),
            data = nest)
plot(ranef(lmm1))

# Include nesting of fishID within plate
nest2 <- groupedData(y ~ conc|apid/fishID, data =  Padilla_DNT60_Loperamide.D)
lmm2 <- lme(y ~ 0 + factor(conc) + factor(conc):poly(t,degree=2,raw=TRUE),
            random = ~1|fishID, method = "REML",
            corr = corARMA(,form=~1|fishID,p=1,q=1),
            data = nest2)
plot(ranef(lmm2))

lmm3 <- lme(y ~ 0 + factor(conc) + factor(conc):poly(t,degree=2,raw=TRUE),
            random = ~1|apid/fishID, method = "REML",
            corr = corARMA(,form=~1|apid/fishID,p=1,q=1),
            data = nest2)
plot(ranef(lmm3), level = 1)
plot(ranef(lmm3), level = 2)

anova(lmm1, lmm2, lmm3)

# does inclusion of apid as a fixed effect help?
lmm.noapid <- lme(y ~ 0 + factor(conc) + factor(conc):poly(t,degree=2,raw=TRUE),
            random = ~1|fishID, method = "ML",
            corr = corARMA(,form=~1|fishID,p=1,q=1),
            data = nest)
lmm.apid <- lme(y ~ 0 + apid + factor(conc) + factor(conc):poly(t,degree=2,raw=TRUE),
            random = ~1|fishID, method = "ML",
            corr = corARMA(,form=~1|fishID,p=1,q=1),
            data = nest)
anova(lmm.apid, lmm.noapid)

# Final model
lmm.final.Loperamide.D <- lme(y ~ 0 + factor(conc) + factor(conc):poly(t,degree=2,raw=TRUE),
                 random = ~1|fishID, method = "REML",
                 corr = corARMA(,form=~1|fishID,p=1,q=1),
                 data = nest)
summary(lmm.final.Loperamide.D)

lmm.final.Loperamide.D.cntrt <- lme(y ~ factor(conc)*poly(t,degree=2,raw=TRUE),
                 random = ~1|fishID, method = "REML",
                 corr = corARMA(,form=~1|fishID,p=1,q=1),
                 data = nest)
summary(lmm.final.Loperamide.D.cntrt)

# compare correlation structures
lmm.final.Loperamide.D.AR1 <- lme(y ~ 0 + factor(conc) + factor(conc):poly(t,degree=2,raw=TRUE),
                                  random = ~1|fishID, method = "ML",
                                  corr = corAR1(,form=~1|fishID),
                                  data = nest)
lmm.final.Loperamide.D.ARMA <- lme(y ~ 0 + factor(conc) + factor(conc):poly(t,degree=2,raw=TRUE),
                                  random = ~1|fishID, method = "ML",
                                  corr = corARMA(,form=~1|fishID,p=1,q=1),
                                  data = nest)

anova(lmm.final.Loperamide.D.AR1, lmm.final.Loperamide.D.ARMA)

lmm.final.Loperamide.D.ARMA.REML <- lme(y ~ 0 + factor(conc) + factor(conc):poly(t,degree=2,raw=TRUE),
                                      random = ~1|fishID, method = "REML",
                                      corr = corARMA(,form=~1|fishID,p=1,q=1),
                                      data = nest)
summary(lmm.final.Loperamide.D.ARMA.REML)

lmm.final.Loperamide.D.ARMA.REML.cntrsts <- lme(y ~ factor(conc)*poly(t,degree=2,raw=TRUE),
                                        random = ~1|fishID, method = "REML",
                                        corr = corARMA(,form=~1|fishID,p=1,q=1),
                                        data = nest)
summary(lmm.final.Loperamide.D.ARMA.REML.cntrsts)$tTable
# Fit to 5,5-Diphenylhydantoin in Dark --------------------------------------------


