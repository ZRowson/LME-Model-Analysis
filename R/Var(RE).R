## ---------------------------
##
## Script Name: Add Concentration Effect to Var(RE)
##
## Purpose of Script:
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-18
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


# Load data ---------------------------------------------------------------


load("data/Padilla_DNT60_Loperamide.rda")

# Isolate Loperamide Light data
rmv <- Padilla_DNT60_Loperamide[is.na(y), unique(fishID)]
Padilla_DNT60_Loperamide.L <- Padilla_DNT60_Loperamide[!(fishID%in%rmv) & (t%in%11:30)]

# Specify hierarchical data structure
nest.Lop.L <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_Loperamide.L)
attach(nest.Lop.L)
nest.Lop.L$t <- t - min(t)
nest.Lop.L$conc <- as.factor(conc)
detach(nest.Lop.L)

# Load 5,5-diphenylhydandoin data
load("data/Padilla_DNT60_pmr0_long.rda")

# Isolate 5,5-diphenylhydandoin exposure data in dark. Remove NA entries
grp <- data_egids(DNT60pmr0_long)[cpid=="5,5-diphenylhydandoin", unique(egid)]
Padilla_DNT60_55diphenyl.L <- data_egids(DNT60pmr0_long)[cpid=="5,5-diphenylhydandoin" | (wllt=="v"&egid==grp)][t %in% 11:30]
rmv <- Padilla_DNT60_55diphenyl.L[is.na(y), unique(fishID)]
rm(DNT60pmr0_long, grp, Padilla_DNT60_Loperamide)

# Specify hierarchical data structure
nest.55.L <- groupedData(y ~ conc|fishID, data =  Padilla_DNT60_55diphenyl.L[!(fishID%in%rmv)])
attach(nest.55.L)
nest.55.L$t <- t - min(t)
nest.55.L$conc <- as.factor(conc)
detach(nest.55.L)


# Fit initial model -------------------------------------------------------

# For both Loperamide and 5,5-diphenylhydandoin,random intercept
lmm1.AR1.ML.Lop.L <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.Lop.L, method = "ML",
                         corr = corAR1(,form=~1|fishID), random = ~1|fishID)
lmm1.AR1.ML.55.L <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.55.L, method = "ML",
                        corr = corAR1(,form=~1|fishID), random = ~1|fishID)


# Experiment with Loperamide in Light -------------------------------------


# Add parameterization of random effects based on concentration and compare fit to base case
lmm1.AR1.ML.var.Lop.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE),
                              random = reStruct(~ 0 + factor(conc) | fishID, REML=FALSE, pdClass="pdDiag"),
                              data = nest.Lop.L,
                              corr = corAR1(, form = ~ t | fishID),
                             method = "ML")

anova(lmm1.AR1.ML.Lop.L, lmm1.AR1.ML.var.Lop.L) # AIC does not improve, but it's likely that for other chemicals it will


# Compare the effect of modelling random slope rather than random intercept. Add ctrl to allow for convergence
ctrl <- lmeControl(maxIter = 200, msMaxIter = 200, tolerance = 1e-40, msTol = 1e-40, minAbsParApVar = 1e-40, abs.tol = 1e-40, rel.tol = 1e-20)
lmm1.AR1.ML.rPoly1.Lop.L <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.Lop.L, method = "ML",
                                corr = corAR1(,form=~poly(t,degree=1,raw=TRUE)|fishID),
                                random = ~poly(t,degree=1,raw=TRUE)|fishID,
                                control = ctrl)
lmm1.AR1.ML.rPoly1.var.Lop.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE),
                                       random = reStruct(~ 0 + factor(conc):poly(t,degree=1,raw=TRUE) | fishID, REML = FALSE, pdClass="pdDiag"),
                                       data = nest.Lop.L,
                                       corr = corAR1(, form = ~ t | fishID), method = "ML")

anova(lmm1.AR1.ML.Lop.L, lmm1.AR1.ML.rPoly1.Lop.L, lmm1.AR1.ML.rPoly1.var.Lop.L)
# AIC decreases when random slope is modeled rather than intercept
# AIC increases when variance of RE's by concentration is included


# Fit both random intercept and slope term with ML and compare to other models. Fit with REML to show error
lmm1.AR1.ML.rIrPoly1.Lop.L <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.Lop.L,
                                corr = corAR1(,form=~poly(t,degree=1,raw=TRUE)|fishID),
                                random = ~1 + poly(t,degree=1,raw=TRUE)|fishID,
                                method = "ML",
                                control = ctrl)
lmm1.AR1.ML.rIrPoly1.var.Lop.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE), data = nest.Lop.L,
                                      corr = corAR1(, form = ~ 1 | fishID),
                                      random = reStruct(~ 0 + factor(conc) + factor(conc):poly(t,degree=1,raw=TRUE) | fishID, REML = FALSE, pdClass="pdDiag"),
                                      method = "ML",
                                      control = ctrl)



anova(lmm1.AR1.ML.var.Lop.L, lmm1.AR1.ML.rPoly1.Lop.L, lmm1.AR1.ML.rPoly1.var.Lop.L, lmm1.AR1.ML.rIrPoly1.Lop.L, lmm1.AR1.ML.rIrPoly1.var.Lop.L)
# LME model with random intercept and random slope with variance parameters dependent on concentration group has lowest AIC

# If I try to fit with REML then a non-positive-definite covariance-variance matrix occurs
lmm1.AR1.REML.rIrPoly1.var.Lop.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=FALSE), data = nest.Lop.L,
                                      corr = corAR1(, form = ~ t | fishID),
                                      random = reStruct(~ 0 + factor(conc) + factor(conc):poly(t,degree=1,raw=FALSE) | fishID, REML = TRUE, pdClass="pdDiag"),
                                      method = "REML",
                                      control = ctrl)
lmm1.AR1.REML.rI.var.Lop.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE), data = nest.Lop.L,
                                        corr = corAR1(, form = ~ t | fishID),
                                        random = reStruct(~ 0 + factor(conc) | fishID, REML = TRUE, pdClass="pdDiag"),
                                        method = "REML",
                                  control = ctrl)
lmm1.AR1.REML.rPoly1.var.Lop.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE), data = nest.Lop.L,
                                  corr = corAR1(, form = ~ t | fishID),
                                  random = reStruct(~ 0 + factor(conc):poly(t,degree=1,raw=FALSE) | fishID, REML = TRUE, pdClass="pdDiag"),
                                  method = "REML",
                                  control = ctrl)
lmm1.AR1.REML.rIrPoly1.var.Lop.L$apVar


# Fit random quadratic term and compare AIC
lmm1.AR1.ML.rPoly2.Lop.L <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.Lop.L, method = "ML",
                                corr = corAR1(,form=~t|fishID), random = ~poly(t,degree=2,raw=TRUE)[,2]|fishID,
                                control = ctrl)
lmm1.AR1.ML.rPoly2.var.Lop.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE),
                                    random = reStruct(~ 0 + factor(conc):poly(t,degree=2,raw=TRUE) | fishID, REML = FALSE, pdClass="pdDiag"),
                                    data = nest.Lop.L,
                                    corr = corAR1(, form = ~ t | fishID), method = "ML")
anova(lmm1.AR1.ML.rPoly1.var.Lop.L, lmm1.AR1.ML.rPoly2.var.Lop.L) # poly2 no good


# Final model will be a LME with random slope term with Var dependent on concentration
lmm.final.Lop.L <- lmm1.AR1.REML.rPoly1.var.Lop.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE),
                                                         random = reStruct(~ 0 + factor(conc):poly(t,degree=1,raw=TRUE) | fishID, REML = TRUE, pdClass="pdDiag"),
                                                         data = nest.Lop.L,
                                                         corr = corAR1(, form = ~ t | fishID), method = "REML")
intervals(lmm.final.Lop.L)


# Experiment with 5,5-diphenylhydandoin in Light -------------------------------------


# Add parameterization of random effects based on concentration and compare fit to base case
lmm1.AR1.ML.var.55.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE),
                            random = reStruct(~ 0 + factor(conc) | fishID, REML=FALSE, pdClass="pdDiag"),
                            data = nest.55.L,
                            corr = corAR1(, form = ~ t | fishID), method = "ML")

anova(lmm1.AR1.ML.55.L, lmm1.AR1.ML.var.55.L) # AIC improves


# Compare the effect of modelling random slope rather than random intercept. Add ctrl to allow for convergence
ctrl <- lmeControl(maxIter = 200, msMaxIter = 200)
lmm1.AR1.ML.rPoly1.55.L <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.55.L, method = "ML",
                                corr = corAR1(,form=~poly(t,degree=1,raw=TRUE)|fishID), random = ~poly(t,degree=1,raw=TRUE)|fishID,
                                control = ctrl)
lmm1.AR1.ML.rPoly1.var.55.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE),
                                    random = reStruct(~ 0 + factor(conc):poly(t,degree=1,raw=TRUE) | fishID, REML = FALSE, pdClass="pdDiag"),
                                    data = nest.55.L,
                                    corr = corAR1(, form = ~ t | fishID), method = "ML")

anova(lmm1.AR1.ML.55.L, lmm1.AR1.ML.rPoly1.55.L, lmm1.AR1.ML.rPoly1.var.55.L)
# AIC decreases when random slope is modeled rather than intercept
# AIC decreases when variance of RE's by concentration is included


# Fit both random intercept and slope term with ML and compare to other models
lmm1.AR1.ML.rIrPoly1.55.L <- lme(y ~ 0 + conc + conc:poly(t,degree=2,raw=TRUE), data = nest.55.L,
                                  corr = corAR1(,form=~poly(t,degree=1,raw=TRUE)|fishID),
                                  random = ~1 + poly(t,degree=1,raw=TRUE)|fishID,
                                  method = "ML",
                                  control = ctrl)
lmm1.AR1.ML.rIrPoly1.var.55.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE), data = nest.55.L,
                                      corr = corAR1(, form = ~ t | fishID),
                                      random = reStruct(~ 0 + factor(conc)*poly(t,degree=1,raw=TRUE) | fishID, REML = FALSE, pdClass="pdDiag"),
                                      method = "ML")

anova(lmm1.AR1.ML.var.55.L, lmm1.AR1.ML.rPoly1.55.L, lmm1.AR1.ML.rPoly1.var.55.L, lmm1.AR1.ML.rIrPoly1.55.L, lmm1.AR1.ML.rIrPoly1.var.55.L)
# LME model with random slope with variance parameters dependent on concentration group has lowest AIC


# If I try to fit with REML then a non-positive-definite covariance-variance matrix occurs
lmm1.AR1.REML.rIrPoly1.var.55.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE), data = nest.55.L,
                                        corr = corAR1(, form = ~ t | fishID),
                                        random = reStruct(~ 0 + factor(conc)*poly(t,degree=1,raw=TRUE) | fishID, REML = TRUE, pdClass="pdDiag"),
                                        method = "REML")
lmm1.AR1.REML.rIrPoly1.var.55.L$apVar


# Final model will be a LME with random slope term with Var dependent on concentration
lmm.final.55.L <- lmm1.AR1.REML.rPoly1.var.55.L <- lme(y ~ 0 + factor(conc)  + factor(conc):poly(t,degree=2, raw=TRUE),
                                                         random = reStruct(~ 0 + factor(conc):poly(t,degree=1,raw=TRUE) | fishID, REML = TRUE, pdClass="pdDiag"),
                                                         data = nest.55.L,
                                                         corr = corAR1(, form = ~ t | fishID), method = "REML")
intervals(lmm.final.55.L)


# Take-Aways --------------------------------------------------------------


# By exploration, final model in Light is LME model with fixed effects (y ~ conc + conc*t + conc*(t^2))
# Random effects of model are a random slope term of t with Var(t) dependent on concentration.
#
# Inclusion of more than one random effect and fitting via REML resulted in following error
# model$apVar
# [1] "Non-positive definite approximate variance-covariance"
#
# Prevented calculating of confidence intervals for model parameters.
